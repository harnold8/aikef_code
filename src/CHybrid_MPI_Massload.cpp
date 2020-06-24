
#include "utils.h"
#include "CHybrid.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <math.h>

// #include <stdio.h>
#include <mpi.h>


INT32 num_values_each_block[NUM_REDIST_OPTIONS] = {8,1,1,1,1,1,1,1,1,1,1};






//!--------------------------------------------------------
//!- count_Blks_create_List
//!--------------------------------------------------------
void CHybrid::count_Blks_create_List(INT32 level,
				     CBlock** &BlkList,
				     INT32  &num_global_blocks,
				     INT32* &num_blocks_at_process)
{


	//! reset values to be safe
	num_global_blocks = 0;
	memset(num_blocks_at_process, 0, mpi_num_processes *sizeof(INT32));


	//! count how many blocks are inside this level and protocol
	//! responsible process
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{
		num_global_blocks++;
		num_blocks_at_process[temp_Block->responsible_mpi_process]++;

		temp_Block = temp_Block->next_Blk_of_BlockList;
	}



	BlkList = new CBlock*[num_global_blocks];
	INT32 blk=0;
	temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		BlkList[blk] = temp_Block;
		blk++;


		temp_Block = temp_Block->next_Blk_of_BlockList;
	}

}



//!--------------------------------------------------------------------
//! - distribute_values_all_procs:
//! - this function scatters the mass load values of all blocks
//!   of the respective level. The blocks of one process do not have to be within
//!   a closed interval along a SFC but can be distributed arbitrarily along
//!   any curve.
//!
//!   The following variables have to be set:
//!   - time_process_fields_incChilds
//!   - time_process_particle_incChilds
//!  these variable of all processes are then scattered into the array:
//!  - global_block_Values
//! 
//!--------------------------------------------------------------------
void CHybrid::distribute_values_all_procs(CBlock** BlkList,
					  INT32 num_global_blocks,
					  INT32* num_blocks_at_process,
					  F_REAL* global_block_Values,
					  INT32 value_option)
{



	log_file << "  distributing values via MPI ...          " << endl;
	clock_t start,finish;
	double time;
	start = clock();


	if(TL<=0 && value_option==PARTICLE_NUMBER)
	log_file << "  ->using estimated particle number ...          " << endl;



	//! just for abreviation and clarification:
	INT32 my_num_blocks = num_blocks_at_process[mpi_myRank];

	//! now the number of root blocks at this process is known,
	//! alloc memory to store massload value for each block of
	//! this process
	F_REAL* my_blocks_Values = new F_REAL[my_num_blocks *num_values_each_block[value_option]];
	memset( my_blocks_Values, 0, my_num_blocks *sizeof(F_REAL));


	//! write values of this process into array
	INT32 my_block_counter=0;
	for(INT32 blk=0; blk < num_global_blocks; blk++)
	{


// 		log_file << BlkList[blk]->Blk_optimal_MPiC[0] << "	";

		//! gather my block massload values into array
		//! -> decomposition along BlkList does not need
		//!    to be connected for one process
		if(BlkList[blk]->responsible_mpi_process == mpi_myRank)
		{

			switch(value_option)
			{
				case AVERAGE_REF_VALUES:
					for(INT32 oct=0; oct<8; oct++)
					my_blocks_Values[my_block_counter*8 +oct] = BlkList[blk]->average_ref_value[oct];
					break;


				case TOTAL_TIME:
					my_blocks_Values[my_block_counter] =  BlkList[blk]->time_process_fields_incChilds
									     +BlkList[blk]->time_process_particle_incChilds;
					break;


				case FIELD_TIME:
					my_blocks_Values[my_block_counter] = BlkList[blk]->time_process_fields;
					break;

				case BLOCK_NUMBER:
					my_blocks_Values[my_block_counter] = 1;
					break;


				case PARTICLE_TIME:
					my_blocks_Values[my_block_counter] = BlkList[blk]->time_process_particle;
					break;

				case PARTICLE_NUMBER:

					if(TL>0)
					my_blocks_Values[my_block_counter] = BlkList[blk]->num_particle;
					else
					{
					//! In TL zero when particle are not filled, precalculated
					//! how many particle will be inserted
						my_blocks_Values[my_block_counter] = 0;

						for(INT32 id_oct=0; id_oct<8; id_oct++)
						 if(!BlkList[blk]->child_array[id_oct])
						 {

							//! get index of oct (0 or 1 for each dimension)
							INT32 a = id_oct/4;
							INT32 b = (id_oct -a*4)/2;
							INT32 c = (id_oct -a*4 -b*2);
							
						
							//! see move method for detailed comments
							INT32 start[3] = {!a *!BlkList[blk]->is_box_boundary[0] +a*BlkNds_X/2,
									 !b *!BlkList[blk]->is_box_boundary[2] +b*BlkNds_Y/2,
									 !c *!BlkList[blk]->is_box_boundary[4] +c*BlkNds_Z/2};
						
						
						
							INT32 end[3]={BlkNds_X -a*!BlkList[blk]->is_box_boundary[1] -!a*(BlkNds_X/2),
								     BlkNds_Y -b*!BlkList[blk]->is_box_boundary[3] -!b*(BlkNds_Y/2),
								     BlkNds_Z -c*!BlkList[blk]->is_box_boundary[5] -!c*(BlkNds_Z/2)};


							for(short species=0; species<num_Charged_Species; species++)
							 for(INT32 i = start[0]; i < end[0]; i++)
							  for(INT32 j = start[1]; j < end[1]; j++)
							   for(INT32 k = start[2]; k < end[2]; k++)
							    if(!BlkList[blk]->Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Y +k])
							    my_blocks_Values[my_block_counter] += BlkList[blk]->Blk_optimal_MPiC[species];

						 }



					}
					break;


				//! -------------- INC_CHILDREN ----------------------------------
				case TOTAL_TIME_INC_CHILDREN:
					my_blocks_Values[my_block_counter] = BlkList[blk]->time_process_fields_incChilds
									    +BlkList[blk]->time_process_particle_incChilds;
					break;


				case FIELD_TIME_INC_CHILDREN:
					my_blocks_Values[my_block_counter] = BlkList[blk]->time_process_fields_incChilds;
					break;


				case PARTICLE_TIME_INC_CHILDREN:
					my_blocks_Values[my_block_counter] = BlkList[blk]->time_process_particle_incChilds;
					break;

				case PARTICLE_NUMBER_INC_CHILDREN:
					my_blocks_Values[my_block_counter] = BlkList[blk]->num_particle_incChilds;
					break;

				default: log_file << "  -> Incorrect values: " << value_option << endl << " Exiting ... " << endl;
					   exit(1);

			}

			my_block_counter++;

		}



	}

	//! specify displacements for receive array
	INT32* displacements = new INT32[mpi_num_processes];
	memset( displacements, 0, mpi_num_processes *sizeof(INT32));

	//! set displacement from:
	//! displacement of process 0 is 0
	for(INT32 process=1; process<mpi_num_processes; process++)
	displacements[process] = displacements[process-1] +(num_blocks_at_process[process-1] *num_values_each_block[value_option]);


	INT32 my_num_elmts = my_num_blocks *num_values_each_block[value_option];


	INT32* num_elmts_at_process = new INT32[mpi_num_processes];
	for(INT32 proc=0; proc<mpi_num_processes; proc++)
	num_elmts_at_process[proc] = num_blocks_at_process[proc] *num_values_each_block[value_option];

	//! 	Allgatherv(const void* sendbuf,
	//! 		   int sendcount,
	//!		   const MPI::Datatype& sendtype,
	//!		   void* recvbuf,
	//!         	   const int recvcounts[],
	//!         	   const int displs[],
	//!         	   const MPI::Datatype& recvtype);
	
	//! MPI_ALLGATHERV(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm)
	//! IN sendbuf:		starting address of send buffer (choice)
	//! IN sendcount:		number of elements in send buffer (integer)
	//! IN sendtype:		data type of send buffer elements (handle)
	//! OUT recvbuf:		address of receive buffer (choice)
	//! IN recvcounts:	integer array (of length group size) containing the number
	//! 				of elements that are received from each process
	//! IN displs:		integer array (of length group size). Entry i specifies
	//! 				the displacement (relative to recvbuf) at which to place
	//! 				the incoming data from process i
	//! IN recvtype:		data type of receive buffer elements (handle)
	//! IN comm:		communicator (handle)


// 	log_file << "  -> Allgatherv ..." << endl;
	mpi_myComm.Allgatherv(my_blocks_Values,
			      my_num_elmts,
			      MPI_F_REAL,
			      global_block_Values,
			      num_elmts_at_process,
			      displacements,
			      MPI_F_REAL);
// 	log_file << "     done." << endl;





	delete[] displacements;
	delete[] my_blocks_Values;
	delete[] num_elmts_at_process;


	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "  -> took: " << time << "s." << endl;

}


//!--------------------------------------------------------
//!- estimate_new_block_process_dependency:
//!  compute which block should be processed by
//!  which process
//!--------------------------------------------------------
bool CHybrid::estimate_new_block_process_dependency(INT32 num_global_blocks,
						    INT32* num_blocks_at_process,
						    F_REAL* global_block_Values,
						    INT32* active_responsible_mpi_process,
						    INT32* new_responsible_mpi_process)
{


	log_file << "  estimating new block <-> process dependency ..." << endl;
	clock_t start,finish;
	double time;
	start = clock();

	//! protocol active mass load for each block
	F_REAL* active_massload_of_process = new F_REAL[mpi_num_processes];
	memset( active_massload_of_process, 0, mpi_num_processes *sizeof(F_REAL));

	//! protocol newly estimated mass load for each block
	F_REAL* new_massload_of_process = new F_REAL[mpi_num_processes];
	memset( new_massload_of_process, 0, mpi_num_processes *sizeof(F_REAL));


	//! sum up total massload (all processes)
	D_REAL global_sum_massload_Values = 0;
	for(INT32 block=0; block<num_global_blocks; block++)
	{



		global_sum_massload_Values += global_block_Values[block];

		INT32 process = active_responsible_mpi_process[block];
		active_massload_of_process[process] += global_block_Values[block];
	}

// 	log_file << "  -> total massload: " << global_sum_massload_Values << endl;


	//! estimate massload / process in case of uniform distribution:
	D_REAL massload_each_process = global_sum_massload_Values/mpi_num_processes;


// 	log_file << "  -> aspired massload each process: " << massload_each_process << endl << endl;



	//! estimate new responsible processes
	const INT32 blks_each_process = ceil( 1.*num_global_blocks/mpi_num_processes);


	if(distribution_criteria == BLOCK_NUMBER)
	{


		//! calculate number of processes that may process blks_each_process-1
		INT32 num_underfull_procs = blks_each_process * mpi_num_processes -num_global_blocks;
		INT32 num_underfull_blks = num_underfull_procs *(blks_each_process-1);

		INT32 active_process = 0;
		INT32 active_process_blks = 0;

        	log_file << "  -> num underfull procs: " << num_underfull_procs << endl;

		
		for(INT32 block=0; block<num_underfull_blks; block++)
		{

			new_responsible_mpi_process[block] = active_process;
			new_massload_of_process[active_process] += global_block_Values[block];

			active_process_blks++;

			if(active_process_blks==(blks_each_process-1))
			{

				active_process++;
				active_process_blks=0;

			}

		}

		for(INT32 block=num_underfull_blks; block<num_global_blocks; block++)
		{

			new_responsible_mpi_process[block] = active_process;
			new_massload_of_process[active_process] += global_block_Values[block];

			active_process_blks++;

			if(active_process_blks==blks_each_process)
			{

				active_process++;
				active_process_blks=0;

			}

		}




	}


	log_file << "  -> max iteration: " << max_iteration_redistribute << endl;
	
	D_REAL inc_fac = 4./max_iteration_redistribute;
	D_REAL fac =inc_fac;

	INT32 Try=0;


	//! for counting blocks at each process
	INT32* new_num_blocks_at_process = new INT32[mpi_num_processes];
	memset(new_num_blocks_at_process, 0, mpi_num_processes*sizeof(INT32));

	if(distribution_criteria != BLOCK_NUMBER)
	 while(Try < max_iteration_redistribute && new_massload_of_process[mpi_num_processes-1] <= massload_each_process)
	 {


		memset(new_num_blocks_at_process, 0, mpi_num_processes*sizeof(INT32));
		memset(new_massload_of_process,   0, mpi_num_processes *sizeof(F_REAL));



		//! estimate new responsible processes
		INT32 active_process = 0;
		D_REAL total_massload_distributed = 0.;
		for(INT32 block=0; block<num_global_blocks; block++)
		{
	

	
			//! NOTE:
			//! use < fac *global_block_Values[block] to allow for some more weight,
			//! else final blk will be loaded with all the rest

			//! NOTE:
			//! use slightly higher value up to 1% above massload_each_process for comparison
			//! -> otherwise it can happen that the first time that process N-1 is assigned a
			//!    huge massload
	
			//! in case this process massload is smaller than average massload
			//! attach block to this process
			if(new_massload_of_process[active_process] +fac*global_block_Values[block]
				<= (1. +0.02*gsl_rng_uniform(randGen_general_synchronized))*massload_each_process
				&& new_num_blocks_at_process[active_process] < multiple_of_average*blks_each_process)
			{
	
				new_num_blocks_at_process[active_process]++;
				new_responsible_mpi_process[block] = active_process;
				new_massload_of_process[active_process] += global_block_Values[block];
	
			}
			else
			{
	
				//! increase active block except last process
				//! is reached
				if(active_process<mpi_num_processes-1)
				active_process++;
	
				new_num_blocks_at_process[active_process]++;
				new_responsible_mpi_process[block] = active_process;
				new_massload_of_process[active_process] += global_block_Values[block];
			}
			
		}

// 		log_file << "  -> massload process[0]: " << new_massload_of_process[0] << endl;
// 		log_file << "  -> massload process[N-1]: " << new_massload_of_process[mpi_num_processes-1] << endl;

		fac += inc_fac;
		Try ++;

	 }

        log_file << "  -> Try: " << Try << endl;
        log_file << "  -> massload each process: " << massload_each_process << endl;	
        log_file << "  -> massload process[N-1]: " << new_massload_of_process[mpi_num_processes-1] << endl;

	F_REAL total_new_deviation = 0.;
	F_REAL total_active_deviation = 0.;

	for(INT32 process=0; process<mpi_num_processes; process++)
	{

		F_REAL new_deviation    = fabs(   new_massload_of_process[process]/massload_each_process -1.);
		F_REAL active_deviation = fabs(active_massload_of_process[process]/massload_each_process -1.);

		total_new_deviation    += new_deviation;
		total_active_deviation += active_deviation;

// 		log_file << "  -> active massload p" << process <<": "
// 			 << active_massload_of_process[process] << " (= "
// 			 << active_deviation
// 			 << " deviation) " << endl;
// 
// 		log_file << "  -> new massload p" << process <<"   : "
// 			 << new_massload_of_process[process] << " (= "
// 			 << new_deviation
// 			 << " deviation) " << endl << endl;
	}

	log_file << "  -> active deviation: " << total_active_deviation << endl;
	log_file << "  -> new deviation:    " << total_new_deviation    << endl;
	log_file << "  -> global massload:  " << global_sum_massload_Values    << endl;


	if(total_active_deviation<=total_new_deviation)
	{


		log_file << "  -> estimate yields no improvement -> using active dependency. "  << endl;
		memcpy(new_responsible_mpi_process, active_responsible_mpi_process, num_global_blocks*sizeof(INT32));

		delete[] new_num_blocks_at_process;
		delete[] new_massload_of_process;
		delete[] active_massload_of_process;

		return false;
	}




	memset(new_num_blocks_at_process, 0, mpi_num_processes*sizeof(INT32));

	for(INT32 block=0; block<num_global_blocks; block++){
//             log_file << block << " --- " << " new: " << new_responsible_mpi_process[block] << " old: " << active_responsible_mpi_process[block] << endl;
            new_num_blocks_at_process[new_responsible_mpi_process[block]]++;
        }

	//! - compare old / new number
	log_file << "  Number of blocks at process: (total "<<num_global_blocks<<"->"<<blks_each_process<<" average)" << endl;
	log_file << "       ( NEW  /  OLD )" << endl; 
	for(INT32 proc=0; proc<mpi_num_processes; proc++)
	{


		log_file << "   p" << proc 
			 << ": ( "  << new_num_blocks_at_process[proc] 
			 << "  /  " << num_blocks_at_process[proc] <<" )"<< endl;
	}

	delete[] new_num_blocks_at_process;
	delete[] new_massload_of_process;
	delete[] active_massload_of_process;


	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "  -> took: " << time << "s." << endl;

	return true;
}





//!--------------------------------------------------------
//!- exchange_blocks_MPI
//!  based on new block <-> process responsibility, send blocks
//!  to processes that take over those of my_rank and receive
//!  those that are no at my_ranks responsibility
//!--------------------------------------------------------
void CHybrid::exchange_blocks_MPI(CBlock** BlkList,
				  INT32 num_global_blocks,
				  INT32* active_responsible_mpi_process,
				  INT32* new_responsible_mpi_process)
{


	log_file << "  exchanging blocks inbetween procs ..." << endl;
	clock_t start,finish;
	double time;
	start = clock();


	//! always use id_Package=0 since only full blocks are send
	//! (not only fractions of the blocks such as octs/edges)
	INT32 id_Package = 0;

	INT32 blks_send=0;
	INT32 blks_received=0;


	//! choose individual request id for each send
	INT32 id_request_Fields =   0;
	INT32 id_request_MPiCInfo = 1;
	INT32 id_request_Particle = 2;

	//! send particle of entire physical Volume to account for
	//! newly injected particles in ghost cells 
	INT32 start_ijk[3] = {       0,        0,        0};
	INT32   end_ijk[3] = {BlkNds_X, BlkNds_Y, BlkNds_Z};
	INT32   entries_in_MPiC_Info_package = BlkNds_X*BlkNds_Y*BlkNds_Z *num_Charged_Species +1;


	//! 1)
	//! send blocks in case my_rank is not responsible any maore
// 	log_file  << "  -> sending Blks ... ";
	for(INT32 blk=0; blk < num_global_blocks; blk++)
	 if(    active_responsible_mpi_process[blk] == mpi_myRank
	     &&  new_responsible_mpi_process[blk] != mpi_myRank)
	 {


		//! set tags
		INT32 Field_Package_Tag    = BlkList[blk]->mpi_tag   +id_request_Fields*total_num_mpi_tags;
		INT32 MPiCInfo_Package_Tag = BlkList[blk]->mpi_tag +id_request_MPiCInfo*total_num_mpi_tags;
		INT32 Particle_Package_Tag = BlkList[blk]->mpi_tag +id_request_Particle*total_num_mpi_tags;

		//! send blk field memory
		//! remove blk memory AFTER receive since
		//! send might be non buffered
		send_MPI_Field(id_ALL_FIELDS,
			       id_request_Fields,
			       new_responsible_mpi_process[blk],
			       Field_Package_Tag,
			       BlkList[blk]);
                
                log_file << "send block: " << blk << "by proc:" << mpi_myRank << endl;

		//! delete blocks info package memory and particle
		//! package memory from last move
		BlkList[blk]->delete_package_memory();


		//! prepare particle package
		BlkList[blk]->prepare_pPackage(start_ijk,
					       end_ijk,
					       id_Package);
		

		//! send Info
		BlkList[blk]->send_MPiC_Package_MPI(id_request_MPiCInfo,
						    MPiCInfo_Package_Tag,
						    new_responsible_mpi_process[blk],
						    entries_in_MPiC_Info_package,
						    BlkList[blk]->MPiC_Info_package[id_Package]);

		
		//! send paticle package
		BlkList[blk]->send_partPackage_MPI(id_request_Particle,
						   Particle_Package_Tag,
						   new_responsible_mpi_process[blk],
						   entries_in_MPiC_Info_package,
						   BlkList[blk]->particle_package[id_Package],
						   BlkList[blk]->MPiC_Info_package[id_Package]);

		blks_send++;

	 }
// 	log_file  << " done. " << endl;

        synchronize_allProcesses_MPI();
        
	//! 2)
	//! receive in case my_rank shall take over
// 	log_file <<"  -> receiving Blks ... ";
	for(INT32 blk=0; blk < num_global_blocks; blk++)
	 if(   active_responsible_mpi_process[blk] != mpi_myRank
	     && new_responsible_mpi_process[blk] == mpi_myRank)
	 {


		//! alloc memory
		BlkList[blk]->alloc_process_specific_Memory();
		
		//! set all fields like Flag eta Flag
		//! (all others will be overwritten just)
		//! (simplicities sake set all here)
		BlkList[blk]->set_Fields();

		//! set tags
		INT32 Field_Package_Tag    = BlkList[blk]->mpi_tag +id_request_Fields*total_num_mpi_tags;
		INT32 MPiCInfo_Package_Tag = BlkList[blk]->mpi_tag +id_request_MPiCInfo*total_num_mpi_tags;
		INT32 Particle_Package_Tag = BlkList[blk]->mpi_tag +id_request_Particle*total_num_mpi_tags;


		//! receive blk field
		recv_MPI_Field(id_ALL_FIELDS,
			       active_responsible_mpi_process[blk],
			       Field_Package_Tag,
			       BlkList[blk]->Field_Type[id_ALL_FIELDS]);
                
                log_file << "receive block: " << blk << "by proc:" << mpi_myRank << endl;

		//! receive MPiCInfo
		BlkList[blk]->receive_infoPackage_MPI(MPiCInfo_Package_Tag,
						      active_responsible_mpi_process[blk],
						      entries_in_MPiC_Info_package,
						      BlkList[blk]->MPiC_Info_package[id_Package]);


		//! receive particle pakackage
		BlkList[blk]->receive_partPackage_MPI(Particle_Package_Tag,
						      active_responsible_mpi_process[blk],
						      entries_in_MPiC_Info_package,
						      BlkList[blk]->particle_package[id_Package],
						      BlkList[blk]->MPiC_Info_package[id_Package]);


		//! insert particle package
		//! NOTE:
		//! if this works use faster method:
		//! insert_pPackage_emptyBlock
		//! which has not been tested yet
		BlkList[blk]->insert_pPackage_emptyBlock(start_ijk,
							 end_ijk,
							 BlkList[blk]->particle_package[id_Package],
							 BlkList[blk]->MPiC_Info_package[id_Package]);

		blks_received++;

	 }
// 	log_file  << " done. " << endl;



	//! 3)
	//! - clean MPI memory by waiting for
	//!   completed send
	//! - delete particle & field memory
	//!   from Blocks that have been send
	//!   (not before since it might be non buffered)
	//! - update "responsible_mpi_process" of EVERY block
// 	log_file  << "  -> cleaning / updating ... ";
	for(INT32 blk=0; blk < num_global_blocks; blk++)
	{
		if(    active_responsible_mpi_process[blk] == mpi_myRank
		    && new_responsible_mpi_process[blk] != mpi_myRank)
		{
	
			//! if new responsibility is at different blk,
			//! send has been called from this blk
			if(sync_mpi_send_rec) {
			  BlkList[blk]->req_is_MPI_send_completed[id_request_Fields].Wait();
			  BlkList[blk]->req_is_MPI_send_completed[id_request_MPiCInfo].Wait();
			  BlkList[blk]->req_is_MPI_send_completed[id_request_Particle].Wait();
			}
			else {
			  BlkList[blk]->wait_all_MPI();
			}

			//! WAIT with delete delete of info package memory
			//! and particle package memory until send complete
			BlkList[blk]->delete_package_memory();
	
			//! WAIT with delete of until send complete
			BlkList[blk]->delete_process_specific_Memory();
			
		}

		//! Finally update responsible_mpi_process:
		//! Update responsible_mpi_process at EVERY process !!!
		//! otw. it may happen that p1 and p3 exchange blocks
		//! which is not kept update at p0.
		//! -> p0 may try to receive from p1 while
		//!    p3 is the new neighbour (sender)
		BlkList[blk]->responsible_mpi_process = new_responsible_mpi_process[blk];
		
	}
// 	log_file  << " done. " << endl;

// 	log_file << "  blks_send:      " << blks_send << endl;
// 	log_file << "  blks_received:  " << blks_received << endl;

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "  -> took: " << time << "s." << endl;


}




//!--------------------------------------------------------
//!- send_child_blocks_to_parent_process
//!--------------------------------------------------------
void CHybrid::send_child_blocks_to_parent_process(INT32 level)
{


	clock_t start,finish;
	double time;
	start = clock();

	log_file << " SEND CHILD BLOCKS TO PARENT PROCESS OF LEVEL " << level << " ..." <<endl;

	CBlock** BlkList = NULL;
	INT32  num_global_blocks = 0;
	INT32* num_blocks_at_process = new INT32[mpi_num_processes];
	memset(num_blocks_at_process, 0, mpi_num_processes *sizeof(INT32));


	//! 1)
	//! count blocks of respective level and
	//! copy block pointer in array for
	//! improved overview
	count_Blks_create_List(level,
			       BlkList,
			       num_global_blocks,
			       num_blocks_at_process);



	//! now the number of blocks at this level is known,
	//! alloc array to store active block process dependency
	INT32* active_responsible_mpi_process = new INT32[num_global_blocks];
	memset( active_responsible_mpi_process, 0, num_global_blocks *sizeof(INT32));

	//! alloc array to store new block process dependency
	INT32* new_responsible_mpi_process = new INT32[num_global_blocks];
	memset( new_responsible_mpi_process, 0, num_global_blocks *sizeof(INT32));

	//! 2)
	//! set new_responsible_mpi_process to parent
	for(INT32 blk=0; blk < num_global_blocks; blk++)
	{
		active_responsible_mpi_process[blk] = BlkList[blk]->responsible_mpi_process;
		new_responsible_mpi_process[blk] = BlkList[blk]->parent->responsible_mpi_process;

	}

	//! 2)
	//! based on new block <-> process responsibility, send blocks
	//! to processes that take over those of my_rank and receive
	//! those that are no at my_ranks responsibility
	exchange_blocks_MPI(BlkList,
			    num_global_blocks,
			    active_responsible_mpi_process,
			    new_responsible_mpi_process);
 


	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " send to parent process time: " << time << "s." << endl << endl;



	delete[] BlkList;
	delete[] num_blocks_at_process;
	delete[] active_responsible_mpi_process;
	delete[] new_responsible_mpi_process;

}


//!--------------------------------------------------------
//!- redistribute_blocks_RB_independent
//!--------------------------------------------------------
void CHybrid::redistribute_blocks(INT32 TL)
{

	if(mpi_num_processes<2)
	return;

	//! force method by passing -1
	if(TL==-1)
	{

		//! only continue in case redistribution is not time based
		//! (which has not been set yet)
		if(   distribution_criteria != PARTICLE_NUMBER
		    && distribution_criteria != PARTICLE_NUMBER_INC_CHILDREN
		    && distribution_criteria != BLOCK_NUMBER)
		return;

	}
	else
	{
	
		if(!TL || !mpi_num_processes) return;
		if(!TL_REDISTRIBUTE_BLOCKS) return;
		if(!(TL % TL_REDISTRIBUTE_BLOCKS==OFFSET_REDISTRIBUTE_BLOCKS) )
		{
			log_file << " Next redistribute blocks in "<< TL_REDISTRIBUTE_BLOCKS -TL%TL_REDISTRIBUTE_BLOCKS +OFFSET_REDISTRIBUTE_BLOCKS<<" TL." << endl;
			return;
		}

	}



	//! only count particles when required
	if(   distribution_criteria == PARTICLE_NUMBER
	    || distribution_criteria == PARTICLE_NUMBER_INC_CHILDREN)
	count_particle_each_block();


	//! only sum up time when required
	if(   distribution_criteria == TOTAL_TIME_INC_CHILDREN
	    || distribution_criteria == FIELD_TIME_INC_CHILDREN
	    || distribution_criteria == PARTICLE_TIME_INC_CHILDREN)
	sum_up_block_timing();



	if(distribute_RB_based)
	redistribute_blocks_RB_based();
	else
	redistribute_blocks_RB_independent();

	//! should be zero anyway ...
	//! just to be safe 
	set_zero_field_inCore(id_PhiDC);
	set_zero_field_inCore(id_BEven);


}
//!--------------------------------------------------------
//!- redistribute_blocks_RB_based
//!--------------------------------------------------------
void CHybrid::redistribute_blocks_RB_based(void)
{


	INT32 level = 0;

	clock_t start,finish;
	double time;
	start = clock();

	log_file << " REDISTRIBUTING BLOCKS OF LEVEL " << level << " RB BASED..." <<endl;

	CBlock** BlkList = NULL;
	INT32  num_global_blocks = 0;
	INT32* num_blocks_at_process = new INT32[mpi_num_processes];
	memset(num_blocks_at_process, 0, mpi_num_processes *sizeof(INT32));


	pre_mesh_refinement();


	//! 1)
	//! count blocks of respective level and
	//! copy block pointer in array
	count_Blks_create_List(level,
			       BlkList,
			       num_global_blocks,
			       num_blocks_at_process);


	//! mass_load values for all blocks that can be exchanged inbetween processes
	//! of all processes 
	F_REAL* global_block_Values = new F_REAL[num_global_blocks];
	memset( global_block_Values, 0, num_global_blocks *sizeof(F_REAL));




	//! 3)
	//! distrubute values to each process via mpi
	distribute_values_all_procs(BlkList,
				    num_global_blocks,
				    num_blocks_at_process,
				    global_block_Values,
				    distribution_criteria);



	//! now the number of blocks at this level is known,
	//! alloc array to store active block process dependency
	INT32* active_responsible_mpi_process = new INT32[num_global_blocks];
	memset( active_responsible_mpi_process, 0, num_global_blocks *sizeof(INT32));

	//! alloc array to store new block process dependency
	INT32* new_responsible_mpi_process = new INT32[num_global_blocks];
	memset( new_responsible_mpi_process, 0, num_global_blocks *sizeof(INT32));


	//! initialize active_responsible_mpi_process with current block distribution 
	for(INT32 blk=0; blk < num_global_blocks; blk++)
	active_responsible_mpi_process[blk] = BlkList[blk]->responsible_mpi_process;


	bool do_exchange;
	//! 4)
	//! based on massload values create new block <-> process responsibility
	do_exchange = estimate_new_block_process_dependency(num_global_blocks,
							    num_blocks_at_process,
							    global_block_Values,
							    active_responsible_mpi_process,
							    new_responsible_mpi_process);



	//! 5)
	//! based on new block <-> process responsibility, send blocks
	//! to processes that take over those of my_rank and receive
	//! those that are no at my_ranks responsibility'
	if(do_exchange)
	exchange_blocks_MPI(BlkList,
			    num_global_blocks,
			    active_responsible_mpi_process,
			    new_responsible_mpi_process);


	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " redistribute time: " << time << "s." << endl << endl;



	delete[] BlkList;

	delete[] num_blocks_at_process;
	delete[] global_block_Values;

	delete[] active_responsible_mpi_process;
	delete[] new_responsible_mpi_process;


	for(INT32 child_level=1; child_level<=MAX_LEVEL; child_level++)
	send_child_blocks_to_parent_process(child_level);


	post_mesh_refinement();


}



//!--------------------------------------------------------
//!- redistribute_blocks_RB_independent
//!--------------------------------------------------------
void CHybrid::redistribute_blocks_RB_independent(void)
{





	log_file << " REDISTRIBUTING BLOCKS RB INDEPENDENT ..." <<endl;
	synchronize_allProcesses_MPI();

	clock_t start,finish;
	double time;
	start = clock();


	CBlock** BlkList = NULL;
	INT32  num_global_blocks = 0;
	INT32* num_blocks_at_process = new INT32[mpi_num_processes];
	memset(num_blocks_at_process, 0, mpi_num_processes *sizeof(INT32));


	pre_mesh_refinement();



	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		log_file << "  redistributing level " << level << " ..." <<endl;
// 		synchronize_allProcesses_MPI();


		//! 1)
		//! count blocks of respective level and
		//! copy block pointer in array
		count_Blks_create_List(level,
				       BlkList,
				       num_global_blocks,
				       num_blocks_at_process);
// 		synchronize_allProcesses_MPI();
	
		//! mass_load values for all blocks that can be exchanged inbetween processes
		//! of all processes 
		F_REAL* global_block_Values = new F_REAL[num_global_blocks];
		memset( global_block_Values, 0, num_global_blocks *sizeof(F_REAL));
	

	
		//! 2)
		//! distrubute values to each process via mpi
		distribute_values_all_procs(BlkList,
					    num_global_blocks,
					    num_blocks_at_process,
					    global_block_Values,
					    distribution_criteria);

// 		synchronize_allProcesses_MPI();
	
	
		//! now the number of blocks at this level is known,
		//! alloc array to store active block process dependency
		INT32* active_responsible_mpi_process = new INT32[num_global_blocks];
		memset( active_responsible_mpi_process, 0, num_global_blocks *sizeof(INT32));
	
		//! alloc array to store new block process dependency
		INT32* new_responsible_mpi_process = new INT32[num_global_blocks];
		memset( new_responsible_mpi_process, 0, num_global_blocks *sizeof(INT32));
	
	
		//! initialize active_responsible_mpi_process with current block distribution 
		for(INT32 blk=0; blk < num_global_blocks; blk++)
		active_responsible_mpi_process[blk] = BlkList[blk]->responsible_mpi_process;
	
	

		bool do_exchange;

		//! 3)
		//! based on massload values create new block <-> process responsibility
		do_exchange = estimate_new_block_process_dependency(num_global_blocks,
								    num_blocks_at_process,
								    global_block_Values,
								    active_responsible_mpi_process,
								    new_responsible_mpi_process);
	
// 		synchronize_allProcesses_MPI();
	
		//! 4)
		//! based on new block <-> process responsibility, send blocks
		//! to processes that take over those of my_rank and receive
		//! those that are no at my_ranks responsibility
		if(do_exchange)
		exchange_blocks_MPI(BlkList,
				    num_global_blocks,
				    active_responsible_mpi_process,
				    new_responsible_mpi_process);

// 		synchronize_allProcesses_MPI();
		
		
		log_file << endl;
	

	
		delete[] BlkList;
		delete[] global_block_Values;
	
		delete[] active_responsible_mpi_process;
		delete[] new_responsible_mpi_process;

	}


	delete[] num_blocks_at_process;

	post_mesh_refinement();

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " redistribute time: " << time << "s." << endl << endl;


}
