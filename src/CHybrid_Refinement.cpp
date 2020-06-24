

#include "CHybrid.h"

#include "utils.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>






//!-------------------------------------------------------------//
//! refine_Mesh: - call by reference function			//
//!-------------------------------------------------------------//
void CHybrid::refine_Mesh(INT32 TL)
{


	//! static refinement
	bool refine_cuboid = false;
	bool refine_sphere = false;

	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		if(TL_REFINE_STATIC_CUBOID[level] == TL)
		refine_cuboid = true;

		if(TL_REFINE_STATIC_SPHERE[level] == TL)
		refine_sphere = true;
	}

	if(refine_cuboid || refine_sphere)
	{

		//! remove gather blocks & mpi memory
		pre_mesh_refinement();

		log_file << endl;
		log_file <<  " /---------------------------------------------/" << endl;
		log_file <<  " /       REFINING STATIC MESH                 /" << endl;
		log_file <<  " /--------------------------------------------/" << endl;

		if(refine_cuboid)
		static_refinement_cuboid();
	
		if(refine_sphere)
		static_refinement_sphere();

		//! create buffer blocks for moment gathering etc
		post_mesh_refinement();

	}




	//! dynamic refinement
	if(!TL) return;
	if(!TL_REFINE_MESH) return;
	if( !(TL % TL_REFINE_MESH ==0) )
	{
		log_file << " Next refine mesh in "<< TL_REFINE_MESH -TL%TL_REFINE_MESH <<" TL." << endl;
		return;
	}





  	log_file << endl;
  	log_file <<  " /---------------------------------------------/" << endl;
  	log_file <<  " /       REFINING DYNAMIC MESH                /" << endl;
  	log_file <<  " /--------------------------------------------/" << endl;

	clock_t start,finish;
	double time;
	start = clock();






	//! set all oct to status unflagged
	reset_flags();



	for(INT32 criteria=0; criteria<num_refcrit; criteria++)
	{


		//! check whether derived fields have to be
		//! pre-calculated
		if(refcrit_field_IDs[criteria] == id_rotB)
		{
			//! BEven is used to calc rotB
			calc_rot(id_BEven, id_rotB);
			FULL_GN_UPDATE(id_rotB);
		}

		if(refcrit_field_IDs[criteria] == id_PISpecies1)
		{

			INT32 species = 0;
			//! provide rho
			//! gather Ui 
			collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);
	
			//! provide rho
			//! provide Ui
			//! gather vth2
			collect_Ui_vth2_of_species(species, id_UI_Species1, id_PISpecies1 ,getVTH2);
		}

		if(refcrit_field_IDs[criteria] == id_PISpecies1     +0)
		{
			//! BEven is used to calc rotB
			collect_Ui_vth2_of_species(0, id_UI_Species1, id_PISpecies1 ,getVTH2);
			FULL_GN_UPDATE(id_rotB);
		}


		//! use respective criteria to flag octs
		set_average_ref_value(criteria);
	
		//! distribute values to all procs and build global min/max
		set_global_MinMax_of_refValues();
	
		//! now that values for each oct are set and global min/max set,
		//! flag octs
		flag_octs(criteria);

	}


	//! flag environment for refinement to avoid
	//! level boundaries on discontinuities
 	flag_environment();



	//! force_refine_entire_block resulting in standard block AMR
// 	flag_full_block_Array();

	//! remove gather blocks & mpi memory
	pre_mesh_refinement();


	//! add blocks corresponding to curlB
	add_Blocks();

	//! remove blocks corresponding to curlB
	remove_Blocks();

	//! check refinement efficiency (ratio of flagged and refined blocks)
	//! NOTE:
	//! FOR SOME REASON THE INTEL COMPILER CRASHES IN LOOP VECTORIZATION WHEN THIS
	//! FUNCTION IS USED
	//! THIS OCCOURS IN VERSION 1.8 BUT NOT IN 1.7, EVEN THOUGH THE FUNCTION
	//! HAS NOT BEEN CHANGED SINCE THEN !!!
	//! ALSO KDEVELOP CRASHES OCCASIONALLY
// 	estimate_refinement_efficiency(id_rotB);




	//! create buffer blocks for moment gathering etc
	post_mesh_refinement();



	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file <<" refining time: " << time << "s." << endl << endl;

	//! the Moments have to be recollected in order to
	//! get a correct J2u convertion after collecting 
	//! u_minus next TL
	collect_RHOnp1_UIplus_LAM_GAM();
	copy_Field(id_UI_minus, id_UI_plus);
	copy_Field(id_rho_n, id_rho_np1);

#ifdef use_neutral_species_as_field	
	init_Neutral_Profile();
#endif
	
#ifdef use_ion_production_as_field	
	init_IonProduction_Profile();
#endif

#ifdef use_dust_species_as_field
	init_Dust_Profile();
#endif	

	init_Eta_Profile();
	
	build_B_total();
	calc_first_E();


}


//!-------------------------------------------------------------//
//! pre_mesh_refinement:
//! - always apply this function before either
//!   -> mesh has changed
//!   -> block have been redistributed (mpi)
//!-------------------------------------------------------------//
void CHybrid::pre_mesh_refinement(void)
{

	//! remove buffer blocks for moment gathering
	remove_gather_Blks();
	manage_parent_recv_memory(false);

}


//!-------------------------------------------------------------//
//! post_mesh_refinement:
//! - always apply this function after either
//!   -> mesh has changed
//!   -> block have been redistributed (mpi)
//!-------------------------------------------------------------//
void CHybrid::post_mesh_refinement(void)
{


	log_file << endl;
	log_file << " POST MESH REFINEMENT:" << endl;

	//! build BTotal rather than restoring it
	log_file << " building BTotal ...";
	build_B_total();
	log_file << " done" << endl;

	//! allocate some aditional memory
	create_gather_Blocks();
	manage_parent_recv_memory(true);

	//! set oPiC / flags
        set_Blk_optimal_MPiC();
	set_block_tags();


	//! reset:
	//! - do_send_particle_to_childArray
	//! - do_receive_particle_from_parent
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			memset(temp_Block->do_send_particle_to_childArray, 0, 8*sizeof(bool));
		 	temp_Block->do_receive_particle_from_parent = false;

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}


	//! mark blocks that must send/recv particle
	for(INT32 level=1; level<=MAX_LEVEL; level++)
	{
		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			temp_Block->check_if_receive_particle_from_parent_is_required();
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}


	log_file << " done." << endl << endl;
        
        //! Gather densities to enable the calculation of RecombinationAlpha after a refinement
        collect_RHOnp1_UIplus_LAM_GAM();
        
}




//!---------------------------------------------------------------//
//! manage_parent_recv_memory: 						  //
//!---------------------------------------------------------------//
void CHybrid::manage_parent_recv_memory(bool allocate)
{

	if(allocate)
	log_file << "  Allocating parent recv memory ..." << endl;
	else
	log_file << "  Deleting parent recv memory ..." << endl;
	
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			bool memory_managed = false;
			INT32 parent_process = temp_Block->responsible_mpi_process;
	
			//! check for every block:
			//! -> has different responsible process?
			//! -> if yes, has children
			//! -> if yes, allocate memory for parent receive
			//! -> make sure that field is only allocated once from
			//!    every process, since all 8 octs share
			//!    the same parent field
			if(mpi_myRank!=parent_process)
			 for(INT32 oct=0; oct<8; oct++)
			  if(temp_Block->child_array[oct])
			  {
	
				INT32 child_process = temp_Block->child_array[oct]->responsible_mpi_process;
				
				if(child_process==mpi_myRank && !memory_managed)
				{
	
					temp_Block->manage_parent_recv_memory(allocate);
					memory_managed = true;
				}
	
			  }
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
	
		}

	}

	log_file << "  done." << endl;

}



//!-------------------------------------------------------------//
//! flag_octs:
//!-------------------------------------------------------------//
void CHybrid::reset_flags(void)
{


	log_file << " RESETTING FLAGS ..." << endl;

	//! Climb up the tree:
	//! Refinement rating of new born Blocks is set to zero
	//! therefore new born blocks won't be refined in this loop
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			//! reset flags of block even in case there are
			//! not processed by my_rank
			memset(temp_Block->child_flag_refinement, 0, 8*sizeof(bool));
			memset(temp_Block->child_flag_removal, 0, 8*sizeof(bool));
			memset(temp_Block->flagged_by_neighbour,  0, 8*sizeof(bool));
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	log_file << " done." <<endl <<endl;


}

//!-------------------------------------------------------------//
//! flag_octs:
//!-------------------------------------------------------------//
void CHybrid::flag_octs(INT32 criteria)
{


	log_file << " FLAGGING OCTS ..." << endl;
	log_file << " (including already refined blocks)" << endl;

	log_file << "  Using Field ID:" << refcrit_field_IDs[criteria] << endl;
	log_file << "  -> related Field Name: '" << Field_Name[refcrit_field_IDs[criteria]] <<"'"<< endl;
	log_file << "  -> component option:" << refcrit_field_comp[criteria] << endl;

	D_REAL rating = 0.;


	INT32 octs_flagged = 0;

	INT32* num_octs_flagged_refinement = new INT32[MAX_LEVEL+1];
	memset(num_octs_flagged_refinement,0, (MAX_LEVEL+1)*sizeof(INT32));

	INT32* num_octs_flagged_removal = new INT32[MAX_LEVEL+1];
	memset(num_octs_flagged_removal,0, (MAX_LEVEL+1)*sizeof(INT32));



	//! if values shall not be normalized, set
	//! global extremums to 1
	if(!refcrit_normalize_to_global_extremum[criteria])
	 for(INT32 level=0; level<MAX_LEVEL; level++)
	 {
		global_min_refValue[level] = 1.;
		global_max_refValue[level] = 1.;
	 }
	

	//! Climb up the tree:
	//! Refinement rating of new born Blocks is set to zero
	//! therefore new born blocks won't be refined in this loop
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{

	   CBlock *temp_Block = BlockList_of_Lev[level];
	   while(temp_Block)
	   {


		//! loop accross 8 octs
		 for(INT32 oct=0; oct<8; oct++)
	         {


			//!---------------------------------------------------------------------
			//! MAXIMUM CRITERIA
			//!---------------------------------------------------------------------
			//! check whether this oct should be refined
			//! (independent on whether it already is)
			if(refcrit_maximum_based[criteria])
			{
				if(temp_Block->average_ref_value[oct]/global_max_refValue[level] >= refine_threshold[level][criteria])
				{
					temp_Block->child_flag_refinement[oct] = true;
					num_octs_flagged_refinement[level]++;
				}
	
				//! check whether this oct should be removed
				//! (independent on whether it already is)
				if(     temp_Block->child_array[oct]
				    && !temp_Block->child_flag_refinement[oct]
				    && temp_Block->average_ref_value[oct]/global_max_refValue[level] < remove_threshold[level][criteria])
				{
					temp_Block->child_flag_removal[oct] = true;
					num_octs_flagged_removal[level]++;
				}
			}

			//!---------------------------------------------------------------------
			//! MINIMUM CRITERIA
			//!---------------------------------------------------------------------
			//! check whether this oct should be refined
			//! (independent on whether it already is)
			if(!refcrit_maximum_based[criteria])
			{
				if(temp_Block->average_ref_value[oct]/global_min_refValue[level] <= refine_threshold[level][criteria])
				{
					temp_Block->child_flag_refinement[oct] = true;
					num_octs_flagged_refinement[level]++;
				}
	
				//! check whether this oct should be removed
				//! (independent on whether it already is)
				if(     temp_Block->child_array[oct]
				    && !temp_Block->child_flag_refinement[oct]
				    && temp_Block->average_ref_value[oct]/global_min_refValue[level] > remove_threshold[level][criteria])
				{
					temp_Block->child_flag_removal[oct] = true;
					num_octs_flagged_removal[level]++;
				}
			}

		 }//! end for oct


		temp_Block = temp_Block->next_Blk_of_BlockList;

	   }//! end while temp_Block
	}//! end for level




	for(INT32 level=0; level<MAX_LEVEL; level++)
	log_file << "  -> num_octs_flagged_refinement L" << level <<" :  " << num_octs_flagged_refinement[level] << endl;
	log_file << endl;

	for(INT32 level=0; level<MAX_LEVEL; level++)
	log_file << "  -> num_octs_flagged_removal L" << level <<" :  " << num_octs_flagged_removal[level] << endl;
	log_file << endl;

	log_file << " done." <<endl <<endl;


	delete[] num_octs_flagged_refinement;
	delete[] num_octs_flagged_removal;

}


//!-------------------------------------------------------------//
//! flag_environment:
//!-------------------------------------------------------------//
void CHybrid::flag_environment(void)
{


	if(!force_refine_environment)
	return;

	log_file << " FLAGGING ENVIRONMENT OF BLOCKS ..." << endl;

	D_REAL rating = 0.;


	INT32 blocks_flagged = 0;

	INT32* num_flagged_neighbours = new INT32[2* (MAX_LEVEL+1)];
	memset(num_flagged_neighbours,0,(2* (MAX_LEVEL+1))*sizeof(INT32));



	//! Climb up the tree:
	//! Refinement rating of new born Blocks is set to zero
	//! therefore new born blocks won't be refined in this loop
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
	
			//! check whether this oct should be refined
			//! (independent on whether it already is)
			//! -> if yes, also flag neighbours for refinement
			//! BUT: only do so if block was flagged due to refinement criteria
			//!      and not by its neighbour
			//!      (flagged_by_neighbour set to sero in set_average_ref_value function)
			for(INT32 child=0; child<8; child++)
			if(temp_Block->child_flag_refinement[child] && !temp_Block->flagged_by_neighbour[child])
			temp_Block->flag_full_environment(child, num_flagged_neighbours);
	
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
	
		}//! end while temp_Block
	}//! end for level




	for(INT32 level=0; level<MAX_LEVEL; level++)
	log_file << "  -> num neighbours flagged in estimated L" << level <<" :  " << num_flagged_neighbours[level] << endl;
	log_file << endl;

	for(INT32 level=0; level<MAX_LEVEL; level++)
	log_file << "  -> num neighbours flagged in lower than estimated L" << level <<" :  " << num_flagged_neighbours[level +(MAX_LEVEL+1)] << endl;


	delete[] num_flagged_neighbours;

       log_file << " done." <<endl <<endl <<endl;
	

}

//!-------------------------------------------------------------//
//! flag_environment:
//!-------------------------------------------------------------//
void CHybrid::flag_full_block_Array(void)
{


	if(!force_refine_entire_block)
	return;

	log_file << " FLAGGING FULL BLOCK ARRAYS ..." << endl;

	D_REAL rating = 0.;


	INT32 blocks_flagged = 0;

	INT32* num_flagged_octs = new INT32[(MAX_LEVEL+1)];
	memset(num_flagged_octs,0,((MAX_LEVEL+1))*sizeof(INT32));



	//! Climb up the tree:
	//! Refinement rating of new born Blocks is set to zero
	//! therefore new born blocks won't be refined in this loop
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{

	   CBlock *temp_Block = BlockList_of_Lev[level];
	   while(temp_Block)
	   {


		//! check whether any child of block is flagged
		bool any_child_flagged = false;
		for(INT32 child=0; child<8; child++)
	        if(temp_Block->child_flag_refinement[child])
	        any_child_flagged = true;

		//! if any child is flagged, flag entire block
		if(any_child_flagged)
		 for(INT32 child=0; child<8; child++)
		  if(!temp_Block->child_flag_refinement[child])
		  {
	        	temp_Block->child_flag_refinement[child] = true;
			num_flagged_octs[level]++;
		  }

		temp_Block = temp_Block->next_Blk_of_BlockList;

	   }//! end while temp_Block
	}//! end for level




	for(INT32 level=0; level<MAX_LEVEL; level++)
	log_file << "  -> num octs in equal child array flagged in L" << level <<" :  " << num_flagged_octs[level] << endl;
	log_file << endl;


	delete[] num_flagged_octs;

       log_file << " done." <<endl <<endl <<endl;
	

}


//!-------------------------------------------------------------//
//! add_Blocks:
//!-------------------------------------------------------------//
void CHybrid::add_Blocks(void)
{


	log_file << " REFINIG BLOCKS ..." << endl;

	D_REAL rating = 0.;
	INT32 total_Blocks = 0;

	INT32 blocks_created = 0;

	INT32* num_created_Blocks = new INT32[MAX_LEVEL+1];
	memset(num_created_Blocks,0,(MAX_LEVEL+1)*sizeof(INT32));

	INT32* request_rejected = new INT32[MAX_LEVEL+1];
	memset(request_rejected,0,(MAX_LEVEL+1)*sizeof(INT32));

	INT32* num_created_env_Blocks = new INT32[MAX_LEVEL+1];
	memset(num_created_env_Blocks,0,(MAX_LEVEL+1)*sizeof(INT32));

	INT32* num_removed_Blocks = new INT32[MAX_LEVEL+1];
	memset(num_removed_Blocks,0,(MAX_LEVEL+1)*sizeof(INT32));



	//! now at each Block Min_Value & Max_Value are set
	//! Also min_Value[level] and max_Value[level] for 
	//! the entire computational Domain are set
	//! rotB is not usesd any more
       log_file << "  -> Number of Blocks before Refinement:" << endl; 
       for(INT32 level=0; level<=MAX_LEVEL; level++)
       log_file << "  -> Blocks in L" << level <<":  " << total_Blocks_L[level] << endl;

       total_Blocks = 0;
       for(INT32 level=0; level<=MAX_LEVEL; level++) total_Blocks += total_Blocks_L[level];

       log_file << "     --------------------------" << endl;
       log_file << "     Total Blocks: " << total_Blocks << endl << endl;


	//! Climb up the tree:
	//! Refinement rating of new born Blocks is set to zero
	//! therefore new born blocks won't be refined in this loop
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{

	   CBlock *temp_Block = BlockList_of_Lev[level];
	   while(temp_Block)
	   {

		for(INT32 child=0; child<8; child++)
		 //! check whether this oct shall be refined
		 if(!temp_Block->child_array[child] && temp_Block->child_flag_refinement[child])
		 {


			if(temp_Block->is_refinement_permitted(child))
			{
				//! CHILD Refinement
				temp_Block->refine_Oct(child);
				num_created_Blocks[level+1]++;
		
	
			}
			//! Request rejected so first refine Enviroment,
			//! then try again
			else 
			{
				//! ENVIROMENT Refinement in ACTIVE LEVEL
				//! refined the blocks around temp_Block if permitted
				//! NOTE:
				//! If environmental refinement is not allowed, way to few
				//! blks will be refined which will cause current sheet on
				//! level boundaries which in turn will cause strong artefacts
				num_refined_Octs=0;
				temp_Block->refine_oct_environment(child);
	
				num_created_env_Blocks[level]+=num_refined_Octs;
	
				//! Env. Arrays are created in equal level (not +1)
				num_created_Blocks[level]+=blocks_created;
	
	
				//! TRY AGAIN TO REFINE
				if(temp_Block->is_refinement_permitted(child))
				{
					//! CHILD Refinement in ACTIVE+1 LEVEL
					temp_Block->refine_Oct(child);
					num_created_Blocks[level+1]++;
	
				}
				else request_rejected[level]++;
			}
		
		 }


	     temp_Block = temp_Block->next_Blk_of_BlockList;

	   }//! end while temp_Block
	}//! end for level




	for(INT32 level=1; level<=MAX_LEVEL; level++)
	log_file << "  -> created Block in L" << level <<" :  " << num_created_Blocks[level] << endl;
	log_file << endl;

	for(INT32 level=1; level<=MAX_LEVEL; level++)
	log_file << "  -> created Block due to env. request in L" 
		<< level <<" :  " << num_created_env_Blocks[level] << endl;
	log_file << endl;

	
	for(INT32 level=1; level<=MAX_LEVEL; level++)
	log_file << "  -> refine request rejected (after env. refine) in L" 
		<< level <<" :  " << request_rejected[level] << endl;
	log_file << endl;
	
	
	
	log_file << "  -> Number of Blocks after Refinement:" << endl; 
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	log_file << "  -> Blocks in L" << level <<":  " << total_Blocks_L[level] << endl;
	
	total_Blocks = 0;
	for(INT32 level=0; level<=MAX_LEVEL; level++) total_Blocks += total_Blocks_L[level];
	
	log_file << "     --------------------------" << endl;
	log_file << "     Total Blocks: " << total_Blocks << endl;


	delete[] num_created_Blocks ;
	delete[] request_rejected;
	delete[] num_created_env_Blocks;
	delete[] num_removed_Blocks;



       log_file << " done." <<endl <<endl <<endl;
	

}


//!-------------------------------------------------------------//
//! refine_Mesh: - call by reference function			//
//!-------------------------------------------------------------//
void CHybrid::remove_Blocks()
{


        log_file << " REMOVING BLOCKS ..." << endl;

	D_REAL rating = 0.;

	num_rem_reject_child_is_refined = 0;
	num_rem_reject_initial_refined = 0;
	num_rem_reject_refined_neighbour = 0;

	INT32* num_removed_Blocks = new INT32[MAX_LEVEL+1];
	memset(num_removed_Blocks,0,(MAX_LEVEL+1)*sizeof(INT32));

	INT32* num_removed_rejected = new INT32[MAX_LEVEL+1];
	memset(num_removed_rejected,0,(MAX_LEVEL+1)*sizeof(INT32));

	//! Climb up the tree:
	//! Refinement rating of new born Blocks is set to zero
	//! therefor new born blocks won't be refined in this loop
	for(INT32 level=MAX_LEVEL-1; level>=0; level--)
	{

	   CBlock *temp_Block = BlockList_of_Lev[level];
	   while(temp_Block)
	   {

		//! child might be flagged for both refinement and removal,
		//! since it has been flaged my neighbour, environment etc...
		//! in this case DO NOT remove child
		for(INT32 child=0; child<8; child++)
		 if(    temp_Block->child_array[child]
		     &&  temp_Block->child_flag_removal[child]
		     && !temp_Block->child_flag_refinement[child])
		    if(temp_Block->is_removing_permitted(child))
		    {

			//! REMOVE
			temp_Block->remove_Child(child);
			num_removed_Blocks[level+1]++;

		    }
		    else
		    num_removed_rejected[level+1]++;




	       temp_Block = temp_Block->next_Blk_of_BlockList;

	   }//! end while temp_Block
	}//! end for level

       for(INT32 level=1; level<=MAX_LEVEL; level++)
       log_file << "  -> removed Blocks in L" << level <<":  " << num_removed_Blocks[level] << endl;
       log_file << endl;


       for(INT32 level=1; level<=MAX_LEVEL; level++)
       log_file << "  -> remove request rejected in L"
            << level <<" :  " << num_removed_rejected[level] << endl;
       log_file << endl;

       log_file << "  reason for rejection:" << endl;
       log_file << "  -> num_rem_reject_initial_refined:   " << num_rem_reject_initial_refined << endl;
       log_file << "  -> num_rem_reject_child_is_refined:  " << num_rem_reject_child_is_refined << endl;
       log_file << "  -> num_rem_reject_refined_neighbour: " << num_rem_reject_refined_neighbour << endl;
       log_file << endl;


       log_file << "  Number of Blocks after Removal:" << endl; 
       for(INT32 level=0; level<=MAX_LEVEL; level++)
       log_file << "  -> Blocks in L" << level <<":  " << total_Blocks_L[level] << endl;

       INT32 total_Blocks = 0;
       for(INT32 level=0; level<=MAX_LEVEL; level++) total_Blocks += total_Blocks_L[level];

       log_file << "     --------------------------" << endl;
       log_file << "     Total Blocks: " << total_Blocks << endl;

       delete[] num_removed_Blocks;
       delete[] num_removed_rejected;

       log_file << " done." <<endl <<endl <<endl;


}


//!--------------------------------------------------------
//!- create_gather_Blocks:
//!--------------------------------------------------------
void CHybrid::create_gather_Blocks(void)
{

	num_refined_Octs = 0;

	
	log_file << "  Creating Gather-Blocks: " << endl;


	//! no gather Blocks possible at Level 0, so start at Level 1
	if(use_gather_blocks)
	for(INT32 level=1; level<=MAX_LEVEL; level++)
	{
		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			temp_Block->refine_gatherEnvironment();
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}



	log_file 	<< "  -> refined gather Octs: " << num_refined_Octs << endl;
        log_file << "  done." <<endl;


}

//!--------------------------------------------------------
//!- remove_gather_Blks:
//!--------------------------------------------------------
void CHybrid::remove_gather_Blks(void)
{

	if(!use_gather_blocks)
	return;

	num_removed_gatherBlks = 0;

	log_file 	<< " Removing gather blocks... " << endl;
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock *temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{

		     for(INT32 gChild=0; gChild<8; gChild++)
		      if(temp_Block->gather_child_array[gChild])
		      temp_Block->remove_gather_Blk(gChild);

		    temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}
	}
	log_file 	<< "  -> removed gather Blks: " << num_removed_gatherBlks << endl;

       log_file << " done." <<endl <<endl;


}

/*
//!------------------------------------------------------------------------
//!- estimate_refinement_efficiency:
//!  global_max_refValue must be set now
//!  (is done in "set_average_ref_value")
//!------------------------------------------------------------------------
void CHybrid::estimate_refinement_efficiency(INT32 src_type)
{


   INT32* num_refined_nodes = new INT32[MAX_LEVEL];
   INT32* num_flagged_nodes = new INT32[MAX_LEVEL];


   memset(num_refined_nodes, 0, MAX_LEVEL *sizeof(INT32));
   memset(num_flagged_nodes, 0, MAX_LEVEL *sizeof(INT32));


    INT32 num_phys_nodes_in_block = (BlkNds_X-2) *(BlkNds_Y-2) *(BlkNds_Z-2);
   //!-----------------------------------------//
   //! 	Loop over all Blks		     //
   //!----------------------------------------//
   for(INT32 level=0; level<MAX_LEVEL; level++)
   {
	CBlock *temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{


		//! COUNT NUMBER OFF REFINED NODES (that are refined)
		//! in 3D each node is refined in 8 nodes
		//! in case alchildren exist
		//! -> 8 additional num_phys_nodes_in_block exist in L+1
		//! -> num_phys_nodes_in_block in L are refined
		for(INT32 child=0; child<8; child++)
		if(temp_Block->child_array[child])
		num_refined_nodes[level] += int(1./8.* num_phys_nodes_in_block);

		
		//! COUNT NUMBER OFF FLAGGED NODES (that should be refined)
		temp_Block->count_flagged_nodes(src_type, global_max_refValue, num_flagged_nodes);
	

		temp_Block = temp_Block->next_Blk_of_BlockList;


	}
   }

   log_file << " ESTIMATING REFINEMENT EFFICIENCY: " << endl;
   for(INT32 level=0; level<MAX_LEVEL; level++)
   {

   	log_file << " Lev"  <<level<< " -> flagged: " << num_flagged_nodes[level] << endl;
   	log_file << "    " << " -> refined: " << num_refined_nodes[level] << endl;

   	log_file << "    " << " -> ratio flaged/refined nodes: " << (1.*num_flagged_nodes[level])/num_refined_nodes[level] << endl;
   }

       log_file << " done." <<endl <<endl <<endl;


    delete[] num_refined_nodes;
    delete[] num_flagged_nodes;

}*/


//!------------------------------------------------------------------------
//!- set_average_ref_value:
//!------------------------------------------------------------------------
void CHybrid::set_average_ref_value(INT32 id_criteria)
{

	for(INT32 level=0; level<MAX_LEVEL; level++)
	{


		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(temp_Block->responsible_mpi_process == mpi_myRank)
			temp_Block->set_average_ref_value(id_criteria);
		
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}

}

//!------------------------------------------------------------------------
//!- calc_rotB:
//!------------------------------------------------------------------------
void CHybrid::set_global_MinMax_of_refValues(void)
{



	log_file << " SETTING GLOBAL MINIMA/MAXIMA ...       " << endl;
	


	CBlock** BlkList = NULL;
	INT32  num_global_blocks = 0;
	INT32* num_blocks_at_process = new INT32[mpi_num_processes];
	memset(num_blocks_at_process, 0, mpi_num_processes *sizeof(INT32));
	
	
	
	for(INT32 level=0; level<=MAX_LEVEL;level++)
	{
		global_min_refValue[level] =  1.e18;
		global_max_refValue[level] = -1.e18;
	
	}


	F_REAL* global_average_Values = NULL;



	//!--------------------------------------------//
	//! 	Loop over all Blks			     //
	//! (MAX_LEVEL requied for removing procedure) //
	//!--------------------------------------------//
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{

		if(!BlockList_of_Lev[level])
		break;


		//! count blocks of respective level and
		//! copy block pointer in array
		count_Blks_create_List(level,
				       BlkList,
				       num_global_blocks,
				       num_blocks_at_process);

		//! allocate array to receive values from all processes
		global_average_Values = new F_REAL[8*num_global_blocks];
		memset(global_average_Values, 0, 8*num_global_blocks *sizeof(F_REAL));

		//! receive values from all processes
		distribute_values_all_procs(BlkList,
					    num_global_blocks,
					    num_blocks_at_process,
					    global_average_Values,
					    AVERAGE_REF_VALUES);

		//! - set all block properties of this process
		//! - estimate global maximum
		for(INT32 blk=0; blk<num_global_blocks; blk++)
		 for(INT32 oct=0; oct<8; oct++)
		 {

			F_REAL rev_value = global_average_Values[blk*8 +oct];

			BlkList[blk]->average_ref_value[oct] = rev_value;

			if(global_min_refValue[level] > rev_value) global_min_refValue[level] = rev_value;
			if(global_max_refValue[level] < rev_value) global_max_refValue[level] = rev_value;
		 }
		 
		//! release Memory that was allocated for this level
		delete[] BlkList;
		delete[] global_average_Values;

	}

	delete[] num_blocks_at_process;
	
	for(INT32 level=0; level<= MAX_LEVEL; level++)
	log_file << "  -> (Min,Max) at Level["<<level<<"]: ("<< global_min_refValue[level] << " , "
							     << global_max_refValue[level] <<")" << endl;
	log_file << " done." <<endl <<endl <<endl;



}



//!------------------------------------------------------------------------
//!  set_block_tags:
//!  - Set individual integer value fo each block
//!  - required for MPI send / recv
//!    always call this fucntion when
//!    mesh was changed
//!------------------------------------------------------------------------
void CHybrid::set_block_tags(void)
{



	log_file << "  Setting block tags...       " << endl;
	
	INT32 blk_tag = 0;
	//!--------------------------------------------//
	//! 	Loop over all Blks			      //
	//!--------------------------------------------//
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock *temp_Block =  GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			//! set individual tag to each block for Send/Recv
			temp_Block->mpi_tag = blk_tag;
			blk_tag++;
	
			temp_Block = temp_Block->next_Blk_of_GatherBlockList;
		}
	}


	total_num_mpi_tags = blk_tag;
	log_file << "  -> total_num_mpi_tags:  " << total_num_mpi_tags << endl;
	log_file << "  -> total_active_Blocks: " << total_active_Blocks << endl;
	log_file << "  done."  << endl;
	
	//! does not match up in case gather block are activated
	if(!use_gather_blocks   &&     blk_tag!=total_active_Blocks)
	{

		log_file << endl;
		log_file << " ERROR: " << total_active_Blocks << endl;
		log_file << " maximal blk_tag:     " << blk_tag << endl;
		log_file << " total_active_Blocks: " << total_active_Blocks << endl;
		log_file << " -> Both values must match up !!!" << endl;
		log_file << "    Else it may cause serious errors in mpi send / recv." << endl;
		log_file << " Exiting ..." << endl;
		exit(1);
	}





}



