



#include "CHybrid.h"

#include "utils.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <math.h>

#define set_InnerRho_after_convert 	true
#define do_not_set_MCD		false

#define EVEN_ITERATION		true
#define ODD_ITERATION		false




D_REAL initial_sum_magnetic_kinetic_Erg;


//!--------------------------------------------------------------------
//! - count_particle_each_block
//!--------------------------------------------------------------------
void CHybrid::count_particle_each_block(void)
{


	//! set incChild to zero
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->count_particle();

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}


	//! climb down from MAX_LEVEL to level 0 in order to sum up time
	//! to root blocks
	for(INT32 level=MAX_LEVEL; level>=1; level--)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			//! add total sum of this block to parents total sum
			temp_Block->parent->num_particle_incChilds += temp_Block->num_particle_incChilds;
		

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}



}

//!--------------------------------------------------------------------
//! - sum_up_block_timing
//!--------------------------------------------------------------------
void CHybrid::sum_up_block_timing(void)
{

// 	if(!TL_reset_block_timing) return;
// 	if( !(TL % TL_reset_block_timing ==0) ) return;

	log_file << "  summing up block timing ...          " << endl;
	clock_t start,finish;
	double time;
	start = clock();


	//! set incChild to zero
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			temp_Block->time_process_fields_incChilds   = temp_Block->time_process_fields;
			temp_Block->time_process_particle_incChilds = temp_Block->time_process_particle;

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}



	//! climb down from MAX_LEVEL to level 0 in order to sum up time
	//! to root blocks
	for(INT32 level=MAX_LEVEL; level>=1; level--)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{

				//! add total sum of this block to parents total sum
				temp_Block->parent->time_process_fields_incChilds +=
				 temp_Block->time_process_fields_incChilds;
				temp_Block->parent->time_process_particle_incChilds +=
				 temp_Block->time_process_particle_incChilds;
		
			}
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}

	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "  -> took: " << time << "s." << endl << endl;

}

//!--------------------------------------------------------------------
//! - reset_block_timing
//!--------------------------------------------------------------------
void CHybrid::reset_block_timing(void)
{

	if(!TL_reset_block_timing) return;
	if( !(TL % TL_reset_block_timing ==0) ) return;

	log_file << " resetting block timing ...          " << endl;
	clock_t start,finish;
	double time;
	start = clock();

	INT32 my_rootBlock_counter=0;

	//! climb down from MAX_LEVEL to level 0 in order to sum up time
	//! to root blocks
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{

				temp_Block->time_process_fields   = 0.;
				temp_Block->time_process_particle = 0.;
				
			}
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}

	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "took: " << time << "s." << endl << endl;

}



//!--------------------------------------------------------------------
//!- calc_first_E:
//!  input Fields:
//!  - Ui_plus: collected last TL
//!  - rho_n:   collected last TL, should now be identically rho_np1
//!  - B_EVEN:  averaged solution of last TL
//!--------------------------------------------------------------------
void CHybrid::calc_first_E(void)
{


	if(TestParticle_Simulation)
	return;

	log_file << " CALC FIRST E...          " <<endl;
	clock_t start,finish;
	double time;
	start = clock();


	//! id_BEven GN must be update now
	INT32 counter=0;

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->calc_E(id_BEven, id_UI_plus, id_rho_n);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	FULL_GN_UPDATE(id_EField);


	smooth_Field(id_EField, smooth_E);


	//! field_to_parent_smooth to get correct forces ?!?
	//! NOT REQUIRED WHEN MOMENTS & BFIELD ARE SMOOTHED TO PARENT

	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file <<" first EField time: " << time << "s." << endl << endl;



}

//!--------------------------------------------------------------------
//!- CAM:
//!  input Fields:
//!  - Ui_plus:    collected last TL
//!  - rho_n:      collected last TL, should now be identically to rho_np1
//!  - B_EVEN:     averaged solution of last TL
//!  - id_rotB:     used as temp scratch for Ji_plus
//!  - id_UI_minus: used as temp scratch for Ji_np1
//!--------------------------------------------------------------------
void CHybrid::CAM(void)
{


	if(TestParticle_Simulation)
	return;

	log_file << " CAM ...          ";
	clock_t start,finish;
	double time;
	start = clock();

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->CAM();
		
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << endl <<" CAM time: " << time << "s." << endl << endl;

}

//!---------------------------------------------------------------------
//!- calc_second_E:
//!---------------------------------------------------------------------
void CHybrid::calc_second_E(void)
{

	if(TestParticle_Simulation)
	return;


	log_file << " CALC SECOND E...          "<<endl;
	clock_t start,finish;
	double time;
	start = clock();

	//! 1) id_BEven GN must be update now


	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->calc_E(id_BEven, id_UI_plus, id_rho_n);
	
		temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	FULL_GN_UPDATE(id_EField);

	smooth_Field(id_EField,smooth_E);
// 	FULL_GN_UPDATE(id_EField);

	//! to get correct forces, smooth moments to parent ?!?
	//! NOT REQUIRED WHEN MOMENTS & BFIELD ARE SMOOTHED TO PARENT
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file <<" second EField time: " << time << "s." << endl << endl;

}

//!---------------------------------------------------------------------
//!- calc_second_E:
//!---------------------------------------------------------------------
void CHybrid::add_force(INT32 field_id)
{

	if(TestParticle_Simulation)
	return;


	log_file << " ADDING FORCE TO FIELD "<< Field_Name[field_id] << " ...          "<<endl;
	clock_t start,finish;
	double time;
	start = clock();

	//! 1) id_BEven GN must be update now


	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->add_force(field_id);
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}



	//! to get correct forces, smooth moments to parent ?!?
	//! NOT REQUIRED WHEN MOMENTS & BFIELD ARE SMOOTHED TO PARENT
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file <<" adding force time: " << time << "s." << endl << endl;

}



//!--------------------------------------------------------
//!- accelerate_Particle:
//!--------------------------------------------------------
void CHybrid::accelerate_Particle(void)
{


	log_file 	<< " ACCELERATING PARTICLE ..." << endl;
	num_accelerated = 0;
	clock_t start,finish;
	double time;
	start = clock();
	
	//!----- ACCELERATE -----------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			for(INT32 oct=0; oct<8; oct++)
			if(!temp_Block->child_array[oct])
			temp_Block->accelerate_particle(oct);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	INT64 local_info_values[2] = {num_accelerated,
				      num_total_particles};


	stringstream info_names[2];
	info_names[0] << "   ->Accelerated Particles: ";
	info_names[1] << "   ->total Particles: ";

	show_information(local_info_values,
			 info_names,
			 2, BUILD_SUM);

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " accelerate time: " << time << "s." << endl << endl;


	//! no syncronization rewuired since no data is exchanged
// 	synchronize_allProcesses_MPI();

}


//!--------------------------------------------------------
//! collect_Ui_vth2_of_species:
//! - this function is exclusively called for output purpouse
//! - collect_vth2=true:
//!    -> collect vth2 of species in id_vth2
//!    -> Ui has to be provided by id_Ji
//! - collect_vth2=false:
//!    -> collect Ui of species in id_Ji
//!--------------------------------------------------------
void CHybrid::collect_Ui_vth2_of_species(INT32 species,
					 INT32 id_Ji,
					 INT32 id_vth2,
					 bool collect_vth2)
{


	//! use id_to_collect to avoid if interceptions and to make code more readable
	INT32 id_to_collect = 0;

	//! Either collect Ji^2 or Ji
	if(collect_vth2)
	{
		id_to_collect = id_vth2;
		log_file 	<< " GATHERING " << Field_Name[id_PISpecies1 +species] << ": "  << endl;
	}
	else
	{
		id_to_collect = id_Ji;
		log_file 	<< " GATHERING " << Field_Name[id_UI_Species1 +species] << ": "  << endl;
	}
	num_total_particles_collected = 0;
	num_particles_collected_from_parent = 0;

	clock_t start,finish;
	double time;
	start = clock();

	//! Collecting are to loop over all levels:
	//! 
	//! 1) collect in every Block which has no children
	//!    and inject collected values down to parent
	//!
	//! 2) 


	//! First add, then devide !!!

	//! 1)  collect all rho_species & Ji_total & Gam total
	//! 2)      add all rho_species & Ji_total & Gam total to neighbors
	//! 3) exchange all rho_species & Ji_total & Gam total Ghost Cells
	//! 4) in one loop over all blocks:
	//!	 a) build rho total
	//!	 b) buils Lam
	//!	 c) calc Ji2Ui


	set_zero_field_incGather(id_to_collect);


	//!------------------------------------------------------------
	//!-- 1) GATEHER LOOP -----------------------------------------
	//!------------------------------------------------------------
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{
		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{
			

				//! gather block
				temp_Block->collect_Ji(species,
						       id_Ji,
						       id_Ji,
						       id_vth2,
						       collect_vth2);



				//! simple inject to parent in case gather blocks are
				//! not used
				if(!use_gather_blocks && level && !temp_Block->is_gatherBlk)
				temp_Block->add_field_to_parent(id_to_collect, id_to_collect);

			}


			temp_Block = temp_Block->next_Blk_of_GatherBlockList;
		}


		//! send injected field to parent via MPI
		if(!use_gather_blocks)
		zero_parent_add_children_field_MPI(level, id_to_collect);


		//! now gathering in level is finished, add boundary moments
		Add_Boundary_Moments_of_L(level,id_to_collect);

	}


	//!------------------------------------------------------------
	//!-- 2) UPDATE GN --------------------------------------------
	//!------------------------------------------------------------
	if(use_gather_blocks)
	moments_to_parent_smooth(id_to_collect);
	else
	FULL_MOMENT_UPDATE(id_to_collect);


	//!------------------------------------------------------------
	//!-- 3) CONVERT CURRENTS AND DENSITIES TO VELOCITIES       ---
	//!--    and inject to parent				             ---
	//!--    leave out this step in case RhoV2mean is collected ---
	//!------------------------------------------------------------


	for(INT32 level=MAX_LEVEL; level>= 0; level--)
	{

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			//! leaving out this conversion basically 
			//! results in pressure p=rho_m*v_th^2
			//! This step must not left out in case 
			//! vp² is calculated above instead of (vp-vm)²
			//! NOTE:
			//! Do this is only true for q=m=1
			//! Do not confuse mass and charge density
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{

				//! decide whether vth2 or common J2U conversion
				if(collect_vth2)
				temp_Block->convert_vth2rho_to_vth2(species,
								    id_allRhoSpecies,
								    id_vth2);
				else
				temp_Block->convert_Ji2Ui(species,
							  id_allRhoSpecies,
							  id_Ji,
							  id_Ji);
			}

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	//! calculate vth² = <vp²> - <vp>² and multiply with Ion_Masses[species],
	//! to get a meassure for thermal energy ~ temperature
	//! the factor 2./3. is du to normilization of ideal gas equation
	//! and Definition of temperature:
	//! T = (m*vth²)/(3*kb)
	//! p = n*kb*T
	//! do to pressure normilization a factor 2 has to be considered:
	//! p0 = (B0²)/(2*mu0)


	//! TODO:
	//! Einzelne quadrierte komponenten <vx²>, <vy²> und <vz²> sammeln und dann quadrierte
	//! komponenten <vx>², <vy>² und <vz>² abziehen. Das als standard Vektorfeld ausgeben!

	if(collect_vth2)
	sub_squaredField_Multiply(id_vth2, id_Ji, 1./*2./3.*Ion_Masses[species]*/);
	
	//! use 2./3. when pressure shall be calculated
//  	sub_squaredField_Multiply(id_vth2, id_Ji, 2./3.*Ion_Masses[species]);
//      sub_squaredField_Multiply(id_vth2, id_Ji, Ion_Masses[species]*16.*amu/elementary_charge*SI_v0*SI_v0);
	//! only do below for caculating pressure
// 	if(collect_vth2)
// 	Multiply_fields(id_vth2,id_vth2, id_rhoSpecies1 +species);

	//! Show information:
	INT64 local_info_values[3] = {num_total_particles,
				      num_total_particles_collected,
				      num_particles_collected_from_parent};

	stringstream info_names[3];
	info_names[0] << "   num_total_particles: ";
	info_names[1] << "   particles collected: ";
	info_names[2] << "   particles collected from parent:";


	show_information(local_info_values,
			 info_names,
			 3, BUILD_SUM);


	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " gather time: " << time << "s." << endl;

}



//!--------------------------------------------------------
//! collect_Ui_minus:
//! Ui_minus denotes the update velocity of this TL 
//! (as particles are already accelerated) but are at 
//! old positions (of Last TL).
//!--------------------------------------------------------
void CHybrid::collect_Ui_minus(void)
{

	if(TestParticle_Simulation)
	return;


	log_file 	<< " GATHERING UI_minus: "  << endl;

	num_total_particles_collected = 0;
	num_particles_collected_from_parent = 0;

	clock_t start,finish;
	double time;
	start = clock();

	//! Collecting are to loop over all levels:
	//! 
	//! 1) collect in every Block which has no children 
	//!    and inject collected values down to parent
	//!
	//! 2) 


	//! First add, then devide !!!

	//! 1)  collect all rho_species & Ji_total & Gam total
	//! 2)      add all rho_species & Ji_total & Gam total to neighbors
	//! 3) exchange all rho_species & Ji_total & Gam total Ghost Cells
	//! 4) in one loop over all blocks:
	//!	 a) build rho total
	//!	 b) buils Lam
	//!	 c) calc Ji2Ui



	set_zero_field_incGather(id_UI_minus);


	//!------------------------------------------------------------
	//!-- 1) GATHER LOOP ------------------------------------------
	//!------------------------------------------------------------
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{
		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{


			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{
			
				//! gather moments
				temp_Block->collect_Ji(ALL_SPECIES,
						       id_UI_minus,
						       id_rotB,
						       id_notDefined,
						       noVTH2);

	
				//! simple inject to parent in case gather blocks are
				//! not used
				if(!use_gather_blocks && level && !temp_Block->is_gatherBlk)
				temp_Block->add_field_to_parent(id_UI_minus, id_UI_minus);
				
				
			}

			temp_Block = temp_Block->next_Blk_of_GatherBlockList;
		}


		//! send injected field to parent
		if(!use_gather_blocks)
		zero_parent_add_children_field_MPI(level, id_UI_minus);

		//! now gathering in level is finished, add boundary moments
		Add_Boundary_Moments_of_L(level, id_UI_minus);

	}

#if defined(use_dust_species_as_field)
	add_dust_current_to_field(id_UI_minus);
#endif
	
	//!------------------------------------------------------------
	//!-- 2) UPDATE GN --------------------------------------------
	//!------------------------------------------------------------
	if(use_gather_blocks)
	moments_to_parent_smooth(id_UI_minus);
	else
	FULL_MOMENT_UPDATE(id_UI_minus);

	//!------------------------------------------------------------
	//!-- 3) CONVERT CURRENTS AND DENSITIES TO VELOCITIES ---------
	//!--    and inject to parent				    ---------
	//!------------------------------------------------------------
	for(INT32 level=MAX_LEVEL; level>= 0; level--)
	{

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{
				
			temp_Block->convert_Ji2Ui(0,
						  id_rho_np1,
						  id_UI_minus,
						  id_UI_minus);
			}
			
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	//! Show information:
	INT64 local_info_values[3] = {num_total_particles,
				      num_total_particles_collected,
				      num_particles_collected_from_parent};

	stringstream info_names[3];
	info_names[0] << "   num_total_particles: ";
	info_names[1] << "   particles collected: ";
	info_names[2] << "   particles collected from parent:";

	show_information(local_info_values,
			 info_names,
			 3, BUILD_SUM);

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " gather time: " << time << "s." << endl << endl;


}






//!--------------------------------------------------------
//!- move_Particle:
//!--------------------------------------------------------
void CHybrid::move_Particle(void)
{


	log_file 	<< " MOVING PARTICLE..." << endl;

#ifdef TL_PROTOCOL_ANY_PARTICLE
if(TL_PROTOCOL_ANY_PARTICLE && TL%TL_PROTOCOL_ANY_PARTICLE==0)
{
	char filename[100];
	sprintf(filename,"%s/phase_func_prc%d_TL%d.txt",data_output_path, mpi_myRank, TL);
	AnyParticle_FILE.open(filename);
}
#endif


	max_PiC = 0;
	num_reset = 0;
	num_moved = 0;
	num_unhooked = 0;
	num_injected_particles = 0;
	num_deleted_too_fast_particle = 0;
	num_pArrayReallocs = 0;
        for(int level=0;level<MAX_LEVEL+1;level++)
            fastest_particle_v2[level] = 0.;


	memset(num_in_obst, 0,num_Charged_Species *sizeof(INT64));
	memset(num_out_of_Box, 0,num_Charged_Species *sizeof(INT64));

	clock_t start,finish;
	clock_t start_wait;
	D_REAL move_time, wait_time, reassign_time;
	start = clock();
	
	//! After particles are moved some have to be attached to a new 
	//! root Block. This should be fairly simple:
	//! just take the first of the respective cells linear list 
	//! and attach it to next Block. Of course all particles 
	//! connected to the first follow automatically. 
	//! However, this won't work with MPI.

	//! NaN occurs in case particle shall be attached to not existing node
	bool NaN_in_Move = false;

	//!----- MOVE -----------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			if(!temp_Block->block_intern_move_particle())
			NaN_in_Move = true;

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	//! check for NaN Before Inter Block communication takes place
	check_for_NaN_MPI(NaN_in_Move);

// 	start_wait = clock();
// 	synchronize_allProcesses_MPI();

	finish = clock();
	move_time = (double(finish)-double(start))/CLOCKS_PER_SEC;

	//! wait time meassures how long the porcess has to wait
	//! inbetween end of movement and start of reassignement
// 	wait_time = (double(finish)-double(start_wait))/CLOCKS_PER_SEC;


	//! NOTE:
	//! just one reassignement is required. Meanwhile ALL particle are
	//! gathered in case of periodic boundaires with and without
	//! ion temperature (exact number, even after several 100 TL) !!!
	log_file 	<< " reassigning particle..." << endl;
	reassign_particle(reassign_time);



	INT32 num_total_in_obst=0;
	INT32 num_total_out_of_Box=0;
	for(INT32 species=0; species<num_Charged_Species; species++)
	{
		num_total_in_obst     += num_in_obst[species];
		num_total_out_of_Box  += num_out_of_Box[species];

	}
	
	
	log_file 	<< "   finished." << endl;

	INT64 local_info_values[8 +2*num_Charged_Species +(MAX_LEVEL+1)];
	local_info_values[0] = num_moved;
	local_info_values[1] =     (num_total_particles
					+num_deleted_too_fast_particle
					+num_total_in_obst
					+num_total_out_of_Box
					-num_injected_particles);
	local_info_values[2] = num_total_particles;
	local_info_values[3] = num_injected_particles;
	local_info_values[4] = num_unhooked;
	local_info_values[5] = num_pArrayReallocs;
	local_info_values[6] = num_particle_total_storage;
	local_info_values[7] = num_deleted_too_fast_particle;

	stringstream info_names[INFO_ARRAY_SIZE];
	info_names[0] << "   ->Moved Particles:     ";
	info_names[1] << "   ->new +obs +out +fast -inj:  ";
	info_names[2] << "   ->new total Particles: ";
	info_names[3] << "   ->Injected Particles:  ";
	info_names[4] << "   ->Unhooked Particles:  ";
	info_names[5] << "   ->num_pArrayReallocs:  ";
	info_names[6] << "   ->num_particle_total_storage: ";

#ifdef DELETE_PARTICLE_IF_V_EXCEEDS
	info_names[7] << "    ->deleted (v>"<<DELETE_PARTICLE_IF_V_EXCEEDS<<"): ";
#else
	info_names[7] << "    ->deleted (v>inf): ";
#endif
	


	//! Deleted Particles in Obstacle:
	for(INT32 species=0; species<num_Charged_Species; species++)
	{

		info_names[8 +species] << "    ->deleted (in Obs) " << Field_Name[id_rhoSpecies1 +species] << ": ";
		local_info_values[8 +species] = num_in_obst[species];
	}

	//! Deleted Particles out of Box:
	for(INT32 species=0; species<num_Charged_Species; species++)
	{
	
		info_names[8 +num_Charged_Species +species] << "    ->deleted (out Box) " 
					<<Field_Name[id_rhoSpecies1 +species]<<": ";
		local_info_values[8 +num_Charged_Species +species] = num_out_of_Box[species];
	}


	//! particle in respective level:
	for(int z=0; z<= MAX_LEVEL; z++)
	{

		info_names[8 +2*num_Charged_Species +z] << "     -> L[" << z << "]:  ";
		local_info_values[8 +2*num_Charged_Species +z] = num_total_particles_in_L[z];
	}

	show_information(local_info_values,
			 info_names,
			 8 +2*num_Charged_Species +(MAX_LEVEL+1),
			 BUILD_SUM);




	//! check whether all particle have been moved 
// 	if(    TL_LOGFILE_globalINFO  && TL%TL_LOGFILE_globalINFO==0
// 	    && local_info_values[1]!=local_info_values[0])
// 	log_file << " WARNING: num particle that have NOT been moved: "
// 		 << local_info_values[1]-local_info_values[0] << endl;


	log_file 	<< "   ->PARTICLE MEMORY used/available:  "
			<< 100.*local_info_values[2]/local_info_values[6]
			<< "%" << endl;

	local_info_values[0] = max_PiC;
        for(int level=0;level<MAX_LEVEL+1;level++)
        {
            local_info_values[1+level] = int(sqrt(fastest_particle_v2[level]));
        }


	stringstream info_names2[1+MAX_LEVEL+1];

	info_names2[0] << "   ->new max_PiC: ";
        for(int level=0;level<MAX_LEVEL+1;level++)
        {
            info_names2[1+level] << "   ->Fastest particle (int) in Level " << level << ":  v=";
        }
	

#ifdef ESTIMATE_FASTEST_PARTICLE
	show_information(local_info_values,
			 info_names2,
			 1+MAX_LEVEL+1,
			 BUILD_MAX);

	//! if velocity is <= than global int(velocity),
	//! set local to global for statistics
        for(int level=0;level<MAX_LEVEL+1;level++)
        {
            if(fastest_particle_v2[level] < local_info_values[1+level]*local_info_values[1+level])
                fastest_particle_v2[level] = local_info_values[1+level]*local_info_values[1+level];
        }
	
#else
	show_information(local_info_values,
			 info_names2,
			 1,
			 BUILD_MAX);
#endif

	if(local_info_values[0]>12000)
	{
		log_file << " WARNING: max_PiC by far to high !!" << endl;
		log_file << "          Inefficient computation !" << endl;
		log_file << " -> increase 'num_merge_each_cell'!" << endl;
	
	}

// 	log_file << " wait time:     " << wait_time << "s. ("
// 		 << int(100.*wait_time/(reassign_time+move_time))<<"%)" << endl;
	log_file << " reassign time: " << reassign_time << "s. ("
		 << int(100.*reassign_time/(reassign_time+move_time))<<"%)" << endl;
	log_file << " move time:     " << move_time << "s." << endl << endl;


#ifdef TL_PROTOCOL_ANY_PARTICLE
if(TL_PROTOCOL_ANY_PARTICLE && TL%TL_PROTOCOL_ANY_PARTICLE==0)
AnyParticle_FILE.close();
#endif


}




//!--------------------------------------------------------
//!- collect_RHOnp1_UIplus_LAM_GAM:
//! UIplus denotes the update velocity of this TL 
//! (as particles are already accelerated) and also
//! update positions (of this TL).
//!TODO:
//! Vielleicht macht es anders herum mehr Sinn:
//! 1) level 0 sammeln bound values to level 1
//! 2) level 1 sammeln (ist nun vollständig)
//! 3) smooth to L0
//! anstelle
//! 1) sammel L1 und bound values to level 0
//! 2) sammel L0 (ist nun vollständig)
//! 3) update penultimate line
//! 4) smooth to LO
//! -> geht nicht weil einige Zellen z.B. 0.5*(0.25 +0.5)=0.375 Dichte
//!    bekommen. Beim Sammeln werden jedoch nur dir Werte:
//!    0.25, 0.5, 0.75 und 1.0 angenommen. Keiner dieser Werte
//!    ergänzt 0.375 zu 1.0. Damit ist offensichtlich, dass sobald
//!    Interpolation bei der Ergänzung der Dichte notwendig ist,
//!    das Schema nicht mehr funktioniert. Es muß damit auf dem
//!    feinen level gesammelt und in den groben injeziert werden !!!
//!    Damit kann der grobe voständig gesammelt werden und kritische
//!    Cellen in den feinen interpolieren
//!    (am besten mal Skizze machen !!!)
//!--------------------------------------------------------
void CHybrid::collect_RHOnp1_UIplus_LAM_GAM(void)
{


	log_file 	<< " GATHERING RHO_np1, UI_plus, LAM, GAM: "  << endl;

	num_total_particles_collected = 0;
	num_particles_collected_from_parent = 0;

	clock_t start,finish;
	double time;
	start = clock();

	//! Collecting are two loop over all levels:
	//! 
	//! 1) collect in every Block which has no children 
	//!    and inject collected values down to parent
	//!
	//! 2) 


	//! First add, then devide !!!

	//! 1)  collect all rho_species & Ji_total & Gam total
	//! 2)      add all rho_species & Ji_total & Gam total to neighbors
	//! 3) exchange all rho_species & Ji_total & Gam total Ghost Cells
	//! 4) in one loop over all blocks:
	//!	 a) build rho total
	//!	 b) buils Lam
	//!	 c) calc Ji2Ui


	//!NOTE:
	//! Do not delete Moments of entire Block directly before it is processed, 
	//! else moments that have been injected by children would be deleted
	//! Rather set zero ALL blocks BEFORE gather procedure
	set_zero_field_incGather(id_UIp_Gam_allRhoSpec);


	//!------------------------------------------------------------
	//!-- 1) GATHER LOOP ------------------------------------------
	//!------------------------------------------------------------
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{

		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{


				//! gather Block
				temp_Block->collect_rho_Ji_LAM_GAM(id_UI_plus,
								   id_Gam,
								   id_rotB);



				//! simple inject to parent in case gather blocks are
				//! not used
				if(!use_gather_blocks && level && !temp_Block->is_gatherBlk)
				temp_Block->add_field_to_parent(id_UIp_Gam_allRhoSpec,
							         id_UIp_Gam_allRhoSpec);

			}
			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}


		//! send injected field to parent
		if(!use_gather_blocks)
		zero_parent_add_children_field_MPI(level, id_UIp_Gam_allRhoSpec);

		//! now gathering in level is finished, add boundary moments
		Add_Boundary_Moments_of_L(level, id_UIp_Gam_allRhoSpec);
	}

#if defined(use_dust_species_as_field)

	add_dust_to_field();
	add_dust_current_to_field(id_UI_plus);
	add_dust_current_to_field(id_Gam);

#endif
	

	//!------------------------------------------------------------
	//!-- 2) UPDATE GN --------------------------------------------
	//!------------------------------------------------------------
	if(use_gather_blocks)
	moments_to_parent_smooth(id_UIp_Gam_allRhoSpec);
	else
	FULL_MOMENT_UPDATE(id_UIp_Gam_allRhoSpec);


	//!------------------------------------------------------------
	//!-- 3) CONVERT CURRENTS AND DENSITIES TO VELOCITIES ---------
	//!--    and inject to parent			      ---------
	//!------------------------------------------------------------
	for(INT32 level=MAX_LEVEL; level>= 0; level--)
	{

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{		
					
				temp_Block->sum_up_rhoSpecies(id_rho_np1, id_Lam);
				
				temp_Block->convert_Ji2Ui(0,
							  id_rho_np1,
							  id_UI_plus,
							  id_UI_plus);
			}

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}




	INT64 local_info_values[3] = {num_total_particles,
				      num_total_particles_collected,
				      num_particles_collected_from_parent};


	stringstream info_names[3];
	info_names[0] << "   num_total_particles: ";
	info_names[1] << "   particles collected: ";
	info_names[2] << "   particles collected from parent:";

	show_information(local_info_values,
			 info_names,
			 3, BUILD_SUM);


	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " gather time: " << time << "s." << endl << endl;



}




//!--------------------------------------------------------
//!- average_Ui_rho_setREZrho:
//!--------------------------------------------------------
void CHybrid::average_Ui_rho_setREZrho(void)
{

	if(TestParticle_Simulation)
	return;

	log_file 	<< " AVERAGE RHO & UI, set rezRHO: "  << endl;


	clock_t start,finish;
	double time;
	start = clock();

	num_total_particles_collected = 0;
	
	total_particle_mass   = 0;
	total_particle_energy = 0;



	//! Prepare Moments that are collected on FIMesh
	//! in case a intermediate mesh is used
	if(mesh_type==STAGGERED)
	{

		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{
			//! store HI velocity in EField
			FIMesh_to_HIMesh(level, id_rho_np1_HIMesh, id_rho_np1);
			FIMesh_to_HIMesh(level,   id_rho_n_HIMesh,   id_rho_n);
		}
	}

	memset(total_particle_momentum, 0, 3*sizeof(D_REAL));
	//!--------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->average_Ui_rho_setREZrho();

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	//! Prepare Moments that are collected on FIMesh
	//! in case a intermediate mesh is used
	if(mesh_type==STAGGERED)
	{
		//! temporary copy id_UI_minus to id_rotB
		//! in order to shift id_UI_minus on HIMesh
		//! afterwards
		//! (id_rotB is a scratch field)
		copy_Field(id_rotB, id_UI_minus);


		for(INT32 level=0; level<=MAX_LEVEL; level++)
		FIMesh_to_HIMesh(level, id_UI_minus, id_rotB);

	}

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " average time: " << time << "s." << endl << endl;



}



//!--------------------------------------------------------
//!- apply_new_BField_Boundaries:
//!--------------------------------------------------------
void CHybrid::apply_new_BField_Boundaries(void)
{

	log_file << endl;
	log_file << "   Setting new BField Boundaries." << endl;
	log_file << "   Step: " << TL -TL_new_Bsw +1<< "/" << num_TL_adapt << endl;

	int TL_rotate =  TL -TL_new_Bsw;

	//! rotate from 0 tp pi
	D_REAL angle = 0.5 *M_PI *TL_rotate/(num_TL_adapt-1);

	D_REAL BX_Bound = cos(angle)*B_sw[0] +sin(angle)*new_B_sw[0];
	D_REAL BY_Bound = cos(angle)*B_sw[1] +sin(angle)*new_B_sw[1];
	D_REAL BZ_Bound = cos(angle)*B_sw[2] +sin(angle)*new_B_sw[2];

	log_file << "   BX_Bound = " << BX_Bound << endl;
	log_file << "   BY_Bound = " << BY_Bound << endl;
	log_file << "   BZ_Bound = " << BZ_Bound << endl;;


	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->apply_new_BField_Boundaries(id_BEven);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}


//!--------------------------------------------------------
//!- FIMesh_to_HIMesh:
//!--------------------------------------------------------
void CHybrid::FIMesh_to_HIMesh(INT32 level, INT32 id_dest, INT32 id_src)
{


	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->FIMesh_to_HIMesh(id_dest, id_src);

		temp_Block = temp_Block->next_Blk_of_BlockList;
	}

	update_GN_of_level(level, id_dest, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);

}

//!--------------------------------------------------------
//!- HIMesh_to_FIMesh:
//!--------------------------------------------------------
void CHybrid::HIMesh_to_FIMesh(INT32 level, INT32 id_dest, INT32 id_src)
{


	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->HIMesh_to_FIMesh(id_dest, id_src);

		temp_Block = temp_Block->next_Blk_of_BlockList;
	}

	update_GN_of_level(level, id_dest, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);

}


//!--------------------------------------------------------
//!- LF_at_cycle:
//!--------------------------------------------------------
bool CHybrid::LF_at_cycle(INT32 level, INT32 CYCLE)
{

	//! mark NaN at averageing procedure
	bool NaN_in_BField = false;


	//! put either BField of HI or FI Mesh in rotB
	if(mesh_type==STAGGERED)
	FIMesh_to_HIMesh(level, id_B_HI, GN_type_after[CYCLE-1]);


	
	//! either compute: 
	//! one step:
	//! 1) dtB = rot(uxB-...)
	//! two step
	//! 1)   E = -uxB
	//! 2) dtB = -rotE
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->LF_at_cycle(CYCLE);

		temp_Block = temp_Block->next_Blk_of_BlockList;
	}


	//! in case two steps LF is used:
	if(!LF_one_step)
	{

		//! E* is stored in id_BDerivative, so update boundaries
		update_GN_of_level(level, id_BDerivative, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);


		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			//! calculate dtB = -rotE
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->dtB_Faraday_at_cycle(CYCLE);
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}

	//! cpy GC before field_to_parent_smooth, as in the
	//! averaging procedure GC are needed
	update_GN_of_level(level, GN_type_after[CYCLE], SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
	
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	update_GN_of_level(level, GN_type_Pe_after[CYCLE], SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
#endif

	//! at final cycle BOdd and BEven are at equal points
	//! in time and need to be averaged into field BEven
	if(CYCLE==NUM_SUB_CYCLE+1)
	{
		//! Now average even & odd solution
		temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			if(temp_Block->LF_average_B_Solution())
			NaN_in_BField = 1;
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}

	//! smooth field to lower level
	//! GN_type_after[odd  number] = BOdd
	//! GN_type_after[even number] = Beven
	smooth_field_to_parent(level, GN_type_after[CYCLE]);

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	smooth_field_to_parent(level, GN_type_Pe_after[CYCLE]);
#endif




	return NaN_in_BField;

}

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
//!--------------------------------------------------------
//!- explicit_midpoint_method_for_Pe:
//!--------------------------------------------------------
bool CHybrid::explicit_midpoint_method_for_Pe(INT32 level)
{

	//! mark NaN at averageing procedure
	bool NaN_in_BField = false;

	log_file << " ADVANCING PLASME PE ..." << endl;
	
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		{
		  temp_Block->explicit_midpoint_method_for_Pe_1();
		}

		temp_Block = temp_Block->next_Blk_of_BlockList;
	}


	//! cpy GC before field_to_parent_smooth, as in the
	//! averaging procedure GC are needed
	update_GN_of_level(level, id_PE_odd, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
	
	

	//! smooth field to lower level
	smooth_field_to_parent(level, id_PE_odd);
	
	
	
	
	temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		{
		  if(temp_Block->explicit_midpoint_method_for_Pe_2())
		  {
		    NaN_in_BField = true;
		  }
		}

		temp_Block = temp_Block->next_Blk_of_BlockList;
	}


	//! cpy GC before field_to_parent_smooth, as in the
	//! averaging procedure GC are needed
	update_GN_of_level(level, id_PE_even, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
	
	

	//! smooth field to lower level
	smooth_field_to_parent(level, id_PE_even);
	
	
	
	return NaN_in_BField;
}
#endif


//!--------------------------------------------------------
//!- smooth_field_to_parent:
//!--------------------------------------------------------
void CHybrid::smooth_field_to_parent(INT32 level, INT32 field_type)
{

		//! this function must not be called in level 0
		if(!level) return;

		//! smooth field to parent
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->field_to_parent_smooth(field_type,
							   field_type);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}


		//! send smoothed to parent field to parent process
		copy_children_field_to_parent_MPI(level, field_type);

}




//!--------------------------------------------------------
//!- advanceB_Staggered:
//!--------------------------------------------------------
void CHybrid::advanceB_LF_globalDT(void)
{

	if(TestParticle_Simulation)
	return;

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	log_file << " ADVANCING PLASMA BFIELD AND PE TOGETHER..." << endl;
#elif defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	log_file << " ADVANCING PLASMA BFIELD AND PE CONSECUTIVELY..." << endl;
#else
	log_file << " ADVANCING PLASMA BFIELD ..." << endl;
#endif

	if(mesh_type==STAGGERED)
	log_file << " (staggered, ";
	else
	log_file << " (un-staggered, ";

	if(LF_one_step==true)
	log_file << "single step)" << endl;
	else
	log_file << "double step)" << endl << "  ";
	
	
	clock_t start,finish;
	double time;
	start = clock();
	
	//! to check whether NaN at any block
	INT32 NaN_in_BField = 0;
	
	
	//! check whether new boundaries shall be applied to box inflow boundaries
	if(TL_new_Bsw && TL>=TL_new_Bsw && TL < TL_new_Bsw+num_TL_adapt)
	apply_new_BField_Boundaries();





	//!---------------------------------------------
	//! Leap Frog BField Method:
	//!---------------------------------------------
	for(INT32 CYCLE=1; CYCLE<=NUM_SUB_CYCLE+1; CYCLE++)
	{

		if(CYCLE<=NUM_SUB_CYCLE)
		log_file <<"[" << CYCLE <<"]";
		else
		log_file <<"[" << NUM_SUB_CYCLE <<"*]";
		
		//! advance either BEven or Bodd
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		NaN_in_BField = LF_at_cycle(level, CYCLE);
	
	}
	
	
	check_for_NaN_MPI(NaN_in_BField);
	FULL_GN_UPDATE(id_BEven);
	
	log_file << endl;
	
	smooth_Field(id_BEven,smooth_B);
	FULL_GN_UPDATE(id_BEven);
	
	//! full gc update inside smooth function
	//! (leave it in, else reults change slightly)
	
	//! now id_BEven is update and smoothed,
	//! build total BField from:
	//! B_total = B_even + B_cfbg needed in:
	//! - EField Calculation
	//! - CAM
	//! - particle Acceleration
	build_B_total();
	
	
#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	//! Explicit Midpoint method for Pe
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	  NaN_in_BField = explicit_midpoint_method_for_Pe(level);
	}
	check_for_NaN_MPI(NaN_in_BField);
#endif
	
#if defined(nonadiabatic_gradPE_TERM)
	
#if defined(nonadiabatic_gradPE_TERM_smooth_PE)
        FULL_GN_UPDATE(id_PE_even);
	smooth_Field(id_PE_even,smooth_PE);
#endif
	FULL_GN_UPDATE(id_PE_even);
	
#endif
	
	
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file <<" advancing time: " << time << "s." << endl << endl;
}

//!--------------------------------------------------------
//!- advance_LfB_globalDT:
//!--------------------------------------------------------
void CHybrid::advanceB_unStaggered(void)
{

	log_file << " ADVANCING PLASMA BFIELD (LF un-staggered) ...";
	
	clock_t start,finish;
	double time;
	start = clock();
	
	
	INT32 NaN_in_BField = 0;
	
	
	//! check whether new boundaries shall be applied to box inflow boundaries
	if(TL_new_Bsw && TL>=TL_new_Bsw && TL < TL_new_Bsw+num_TL_adapt)
	apply_new_BField_Boundaries();

	log_file << endl << "  ";


	//!---------------------------------------------
	//! Leap Frog BField Method:
	//!---------------------------------------------
	for(INT32 CYCLE=1; CYCLE<=NUM_SUB_CYCLE+1; CYCLE++)
	{

		if(CYCLE<=NUM_SUB_CYCLE)
		log_file <<"[" << CYCLE <<"]";
		else
		log_file <<"[" << NUM_SUB_CYCLE <<"*]";
		
		//! advance either BEven or Bodd
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		NaN_in_BField = LF_at_cycle(level, CYCLE);
	
	}

	
	check_for_NaN_MPI(NaN_in_BField);
	FULL_GN_UPDATE(id_BEven);
	
	log_file << endl;
	
	smooth_Field(id_BEven,smooth_B);
	FULL_GN_UPDATE(id_BEven);

	//! full gc update inside smooth function
	//! (leave it in, else reults change slightly)
	
	//! now id_BEven is update and smoothed,
	//! build total BField from:
	//! B_total = B_even + B_cfbg needed in:
	//! - EField Calculation
	//! - CAM
	//! - particle Acceleration
	build_B_total();
	
	
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file <<" advancing time: " << time << "s." << endl << endl;



}


//!-------------------------------------------------------------
//!- RK_Step:
//!--------------------------------------------------------------
#ifndef nonadiabatic_gradPE_TERM
bool CHybrid::RK_Step(INT32 level, INT32 field_type, bool calc_B_dev, D_REAL c1, D_REAL c2)
{

	bool NaN_in_BField = false;

	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		NaN_in_BField = temp_Block->RK_advance_B(calc_B_dev, c1, c2);

		temp_Block = temp_Block->next_Blk_of_BlockList;
	}


	//! As B_dev = B_old + c * K1,
	//! Ghost Nodes of K1 must be updated now !
	update_GN_of_level(level, field_type, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);


	return NaN_in_BField;

}
#endif /* nonadiabatic_gradPE_TERM */



//!--------------------------------------------------------
//!- advance_LfB_globalDT:
//!--------------------------------------------------------
#ifndef nonadiabatic_gradPE_TERM
void CHybrid::advanceB_RK_globalDT(void)
{

	
	log_file << " ADVANCING PLASMA BFIELD (RK) ...";
	clock_t start,finish;
	double time;
	start = clock();

	//! check whether new boundaries shall be applied to box inflow boundaries
	if(TL_new_Bsw && TL>=TL_new_Bsw && TL < TL_new_Bsw+num_TL_adapt)
	apply_new_BField_Boundaries();

	log_file << endl << "  ";
	
	INT32 NaN_in_BField = 0;
	
	for(INT32 CYCLE=1; CYCLE<=NUM_SUB_CYCLE; CYCLE++)
	{
	
		log_file <<"[" << CYCLE <<"]";
		if(CYCLE%10==0) log_file << endl << " ";
	
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{
			

			//! When using larger dt it is NOT sufficient to perform
			//! the 4 RK steps without GC update and finally only
			//! update B_new (as in old Version -> RK_sequence()).
			//! The id_KX GC have to be updated after each computation,
			//! otherwise errors are introduced at Blk Boundaries
			//! and the code crashes.
			//! (this was tested with dx=5, dt=0.1, mercury dipol)
	
			//! Ghost Nodes of B_old must be update !
			//!-------------------------------------------------
			//!------------ K1 = f(B_old) ----------------------
			//!-------------------------------------------------
			CBlock* temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{
	
				if(mpi_myRank == temp_Block->responsible_mpi_process)
				temp_Block->RK_advance_B(false,0.,(1./6.));
	
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}
	
	
			//! As B_dev = B_old + c * K1,
			//! Ghost Nodes of K1 must be updated now !
			update_GN_of_level(level, id_KX, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);

// 			RK_Step(level, id_KX, false, 0., (1./6.));
	
			//!-------------------------------------------------
			//!--------- B_dev = B_old + c * K1 ----------------
			//!------------ K2 = f(B_dev) ----------------------
			//!-------------------------------------------------
			temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{
	
				if(mpi_myRank == temp_Block->responsible_mpi_process)
				temp_Block->RK_advance_B(true,0.5,(2./6.));
	
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}
	
			//! As B_dev = B_old + c * K2,
			//! Ghost Nodes of K2 must be updated now !
			update_GN_of_level(level, id_KX, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
	
			//!-------------------------------------------------
			//!--------- B_dev = B_old + c * K2 ----------------
			//!------------ K3 = f(B_dev) ----------------------
			//!-------------------------------------------------
			temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{
	
				if(mpi_myRank == temp_Block->responsible_mpi_process)
				temp_Block->RK_advance_B(true,0.5,(2./6.));
	
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}
	
			//! As B_dev = B_old + c * K3,
			//! Ghost Nodes of K3 must be updated now !
			update_GN_of_level(level, id_KX, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
	
			//!-------------------------------------------------
			//!--------- B_dev = B_old + c * K3 ----------------
			//!------------ K4 = f(B_dev) ----------------------
			//!-------------------------------------------------
			temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{
	
	
				//! check final advancement (-> c1 = 1.) for NaN in BField
				if(mpi_myRank == temp_Block->responsible_mpi_process)
				if(temp_Block->RK_advance_B(true,1.,(1./6.)))
				NaN_in_BField = 1;
	
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}
	
			//! Cp GC before field_to_parent_smooth, as in the
			//! averaging procedure GC are needed
			update_GN_of_level(level, id_BNew, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
	
			//! smooth final field to parent
// 			smooth_field_to_parent(level, id_BNew);


		}

	}
	
	
	check_for_NaN_MPI(NaN_in_BField);
	
	
	FULL_GN_UPDATE(id_BNew);
	
	log_file << endl;
	
	smooth_Field(id_BNew,smooth_B);
	FULL_GN_UPDATE(id_BNew);

	//! full gc update inside smooth function
	//! (leave it in, else reults change slightly)

	
	//! now id_BNew is update and smoothed,
	//! build total BField from:
	//! B_total = B_even + B_cfbg needed in:
	//! - EField Calculation
	//! - CAM
	//! - particle Acceleration
	build_B_total();
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file <<" advancing time: " << time << "s." << endl << endl;
	

}
#endif /* nonadiabatic_gradPE_TERM */




//!--------------------------------------------------------
//!- advanceB_Obstacle:
//!--------------------------------------------------------
void CHybrid::advanceB_Obstacle(void)
{


	if(TestParticle_Simulation)
	return;


	if(advance_obstacle_B)
	solve_Obstacle_BField();
	
	
	
	if(div_cleaner)
	clean_divergence();
	



}


//!--------------------------------------------------------
//!- Split_Particle:
//!--------------------------------------------------------
void CHybrid::Split_Particle(void)
{

	if(TestParticle_Simulation)
	return;

	if(!TL_SPLIT || !(TL % TL_SPLIT==0))
	return;


	clock_t start,finish;
	double time;
	start = clock();

	log_file 	<< " SPLITTING PARTICLE..." << endl;

// 	if(split_heaviest_particle)
	log_file 	<< "  (heaviest)" << endl;
// 	else
// 	log_file 	<< "  (most cell centred)" << endl;



	num_pArrayReallocs = 0;
	num_split_canceled_low_weight = 0;
	num_split_canceled_out_of_cell = 0;
	
  	memset(num_split_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));

	//!----- SPLIT -----------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];

		if(do_split_in_L[level])
		 while(temp_Block)
		 {


			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{
				if(split_heaviest_particle)
				temp_Block->split_heaviest_particle();
				else
				temp_Block->split_most_centred_particle();
			}

			//! to use merge and split at the same time for some reason does not work well
			//! -> 1) it looks incorrect 
			//! -> 2) conservation quantities (Erg etc.) are 
			//! ->    conserved quite badly

			//! Independent of that it seems that the intervals where 
			//! to merge or to split must not overlap, otherwise 
			//! conservation of Erg etc. does not work well
			temp_Block = temp_Block->next_Blk_of_BlockList;

		 }
	}



	//! show information
	long num_split = 0;
	for(int z=0; z<= MAX_LEVEL; z++)
	num_split += num_split_in_L[z];

	INT64 local_info_values[6+ (MAX_LEVEL+1)];
	local_info_values[0] = num_split;
	local_info_values[1] = num_split_canceled_low_weight;
	local_info_values[2] = num_split_canceled_out_of_cell;
	local_info_values[3] = num_total_particles;
	local_info_values[4] = num_pArrayReallocs;



	stringstream info_names[INFO_ARRAY_SIZE];
	info_names[0] << "   ->splitted Particle:    ";
	info_names[1] << "   ->Canceled low weight:  ";
	info_names[2] << "   ->Canceled out of cell: ";
	info_names[3] << "   ->total Particles:      ";
	info_names[4] << "   ->num_pArrayReallocs:   ";


	//! particle in respective level:
	for(int z=0; z<= MAX_LEVEL; z++)
	{

		info_names[5 +z] << "     -> L[" << z  << "]:  ";
		local_info_values[5 +z] = num_split_in_L[z];
	}


	show_information(local_info_values,
			 info_names,
			 5+ (MAX_LEVEL+1), BUILD_SUM);

	log_file 	<< "   ->ratio canceld/splited:  " <<
			 1.*(local_info_values[1]+local_info_values[2])/local_info_values[0]
					 << endl;

	num_split_since_start += local_info_values[0];



	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " split time: " << time << "s." << endl << endl;

}


//!--------------------------------------------------------
//!- Merge_Particle:
//!--------------------------------------------------------
void CHybrid::Merge_Particle(void)
{

	if(TestParticle_Simulation)
	return;

	if(!TL_MERGE || !(TL % TL_MERGE==0))
	return;

	clock_t start,finish;
	double time;
	start = clock();

	//!----- MERGE -----------------------------------------
	log_file << " MERGING PARTICLE..." << endl;
  	memset(num_merge_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));
	num_merge_canceled_out_of_cell = 0;



	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];

		if(do_merge_in_L[level])
		 while(temp_Block)
		 {

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->merge_particle();


			temp_Block = temp_Block->next_Blk_of_BlockList;
		 }
	}



	//! show information
	long num_merge = 0;
	for(int z=0; z<= MAX_LEVEL; z++)
	num_merge += num_merge_in_L[z];

	INT64 local_info_values[3+ (MAX_LEVEL+1)];
	local_info_values[0] = num_merge;
	local_info_values[1] = num_merge_canceled_out_of_cell;
	local_info_values[2] = num_total_particles;


	stringstream info_names[INFO_ARRAY_SIZE];
	info_names[0] << "   ->merged Particle:      ";
	info_names[1] << "   ->Canceled out of cell: ";
	info_names[2] << "   ->total Particles:      ";


	//! merges in respective level:
	for(int z=0; z<= MAX_LEVEL; z++)
	{
		info_names[3 +z] << "     -> L[" << z <<"]:  ";
		local_info_values[3 +z] = num_merge_in_L[z];
	}


	show_information(local_info_values,
			 info_names,
			 3 +(MAX_LEVEL+1), BUILD_SUM);

	log_file 	<< "   ->ratio canceld/merged:  " <<
			 1.*(local_info_values[1])/local_info_values[0] << endl;


	num_merge_since_start += local_info_values[0];


	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " merge time: " << time << "s." << endl << endl;



}

//!--------------------------------------------------------
//!- build_B_total:
//!- basically add CFBF to BEven
//!--------------------------------------------------------
void CHybrid::build_B_total(void)
{

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			//! also necessary in case has child
			//! as EField boundaries at level
			//! borders come from those Blocks
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->build_B_total();
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
}

//!--------------------------------------------------------
//!- show_MPiC:
//!--------------------------------------------------------
void CHybrid::check_weight_sorting(void)
{

	if(TestParticle_Simulation)
	return;

#ifdef SORT_BY_WEIGHT
	if(   (!TL_SPLIT || !(TL % TL_SPLIT==0))
	    &&(!TL_MERGE || !(TL % TL_MERGE==0)))
	return;

	clock_t start,finish;
	double time;
	start = clock();
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->check_weight_sorting();


			temp_Block = temp_Block->next_Blk_of_BlockList;


		}
	}
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " check weight Sorting time: " << time << "s." << endl << endl;

#endif
}

//!--------------------------------------------------------
//!- show_MPiC:
//!--------------------------------------------------------
void CHybrid::resize_pArrays(void)
{


	if(!TL_RESIZE_pARRAYS || !(TL % TL_RESIZE_pARRAYS==0))
	return;


	log_file << " RESIZING PARTICLE ARRAYS ..." << endl;

	//! show information
	//! (has to be 3 to fit all values in, see below)
	INT64 local_info_values[3];
	local_info_values[0] = num_total_particles;
	local_info_values[1] = num_particle_total_storage;



	stringstream info_names[3];
	info_names[0] << "   ->num_total_particles:        ";
	info_names[1] << "   ->num_particle_total_storage: ";

	show_information(local_info_values,
			 info_names,
			 2, BUILD_SUM);

	log_file << "   ->used/available before realloc:  " 
		 << 100.*local_info_values[0]/local_info_values[1] << "%" << endl;


	num_pArrayReallocs=0;

	clock_t start,finish;
	double time;
	start = clock();
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->resize_all_pArrays_of_Block();


			temp_Block = temp_Block->next_Blk_of_BlockList;


		}
	}
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;

	//! show information
	local_info_values[0] = num_total_particles;
	local_info_values[1] = num_particle_total_storage;
	local_info_values[2] = num_pArrayReallocs;


	stringstream info_names2[3];
	info_names2[0] << "   ->num_total_particles:        ";
	info_names2[1] << "   ->num_particle_total_storage: ";
	info_names2[2] << "   ->num_pArrayReallocs:         ";

	show_information(local_info_values,
			 info_names2,
			 3, BUILD_SUM);

	log_file << "   ->used/available after realloc:  " 
		 << 100.*local_info_values[0]/local_info_values[1] << "%" << endl;

	log_file << " resize particle array time: " << time << "s." << endl << endl;

}


//!--------------------------------------------------------
//!- show_PE:
//!--------------------------------------------------------
void CHybrid::show_PE(INT32 species)
{

#if defined nonadiabatic_gradPE_TERM
	//  nothing to do because PE is already computed
	return;
#else

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			//! in case Block has no children it 
			//! is filled with particle
			if(temp_Block->responsible_mpi_process==mpi_myRank)
			temp_Block->show_PE(species);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

#endif


}

//!--------------------------------------------------------
//!- show_responsible_Proc:
//!--------------------------------------------------------
void CHybrid::calc_macro_Force(INT32 MacroForce_id, INT32 EField_id, INT32 UField_id,  INT32 BField_id)
{



	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			//! in case Block has no children it 
			//! is filled with particle
			if(temp_Block->responsible_mpi_process==mpi_myRank)
			temp_Block->calc_macro_Force( MacroForce_id, EField_id, UField_id, BField_id);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}


}


//!--------------------------------------------------------
//!- show_responsible_Proc:
//!--------------------------------------------------------
void CHybrid::show_responsible_Proc(INT32 field_id)
{

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(temp_Block->responsible_mpi_process==mpi_myRank)
			for(INT32 b=0; b<num_nodes_in_block; b++)
			temp_Block->Field_Type[field_id][b] = mpi_myRank;
// 			temp_Block->Field_Type[field_id][b] = temp_Block->Block_Nr;


			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------
//!- show_MPiC:
//!--------------------------------------------------------
void CHybrid::show_MPiC(INT32 field_id)
{

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(temp_Block->responsible_mpi_process==mpi_myRank)
			{
				memset(temp_Block->Field_Type[field_id], 0, num_nodes_in_block*sizeof(D_REAL));
                                
				for(INT32 b=0; b<num_nodes_in_block; b++)
                                {
                                    for(INT32 species=0; species<num_Charged_Species; species++)
                                    {
//                                         log_file << "temp_Block->num_MPiC[species][b] : " << temp_Block->num_MPiC[species][b] << endl;
                                        temp_Block->Field_Type[field_id][b]+= 1.* temp_Block->num_MPiC[species][b];
                                    }
                                }
			}
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
}

//!--------------------------------------------------------
//!- show_MPiC:
//!--------------------------------------------------------
void CHybrid::show_time(INT32 opt, INT32 field_id)
{


	sum_up_block_timing();

	//! field time
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(temp_Block->responsible_mpi_process==mpi_myRank)
			{

				//! field time
				if(opt == 0)
				for(INT32 b=0; b<num_nodes_in_block; b++)
				temp_Block->Field_Type[field_id][b]
						= ( temp_Block->time_process_fields_incChilds );
	
				//! part time
				if(opt == 1)
				for(INT32 b=0; b<num_nodes_in_block; b++)
				temp_Block->Field_Type[field_id][b]
						= ( temp_Block->time_process_particle_incChilds );
	
				//! field time + part time
				if(opt == 2)
				for(INT32 b=0; b<num_nodes_in_block; b++)
				temp_Block->Field_Type[field_id][b]
						= ( temp_Block->time_process_fields_incChilds
						   +temp_Block->time_process_particle_incChilds );

			}
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}


}

//!--------------------------------------------------------
//!- show_Block_Level:
//!--------------------------------------------------------
void CHybrid::show_Block_Level(INT32 field_id)
{

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			//! in case Block has no children it 
			//! is filled with particle
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{
				for(INT32 b=0; b<num_nodes_in_block; b++)
				temp_Block->Field_Type[field_id][b]
					= 1. *level;

			}
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}








//!--------------------------------------------------------
//!- fill_Paticle_in_Box:
//!--------------------------------------------------------
void CHybrid::collect_Energy_Momentum_Mass(void)
{

	


	total_particle_mass     = 0;
	total_particle_energy   = 0;
	total_collected_density = 0;

	memset(total_magnetic_field,        0, 3*sizeof(D_REAL));
	memset(total_particle_momentum,     0, 3*sizeof(D_REAL));
	memset(total_collected_momentum,    0, 3*sizeof(D_REAL));
	memset(total_collected_temperature, 0, 3*sizeof(D_REAL));





	//! provide rho
	//! gather Ui 
	collect_Ui_vth2_of_species(0, id_UI_Species1, 0, noVTH2);

	//! provide rho
	//! provide Ui
	//! gather vth2
	collect_Ui_vth2_of_species(0, id_UI_Species1, id_PISpecies1 ,getVTH2);



	//!--------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->calc_Energy_Momentum_Mass();

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}





	D_REAL global_values[16];
	D_REAL local_values[16] = {total_particle_mass,
				   total_particle_energy,
				   total_collected_density,
				   total_magnetic_field[0],
				   total_magnetic_field[1],
				   total_magnetic_field[2],
				   total_particle_momentum[0],
				   total_particle_momentum[1],
				   total_particle_momentum[2],
				   total_collected_momentum[0],
				   total_collected_momentum[1],
				   total_collected_momentum[2],
				   total_collected_temperature[0],
				   total_collected_temperature[1],
				   total_collected_temperature[2],
				   num_total_particles};
// 
// 
// 
	mpi_build_sum(local_values, global_values, 16);

	total_particle_mass	= global_values[0];
	total_particle_energy   = global_values[1];
	total_collected_density = global_values[2];

	total_magnetic_field[0] = global_values[3];
	total_magnetic_field[1] = global_values[4];
	total_magnetic_field[2] = global_values[5];

	total_particle_momentum[0] = global_values[6];
	total_particle_momentum[1] = global_values[7];
	total_particle_momentum[2] = global_values[8];

	total_collected_momentum[0] = global_values[9];
	total_collected_momentum[1] = global_values[10];
	total_collected_momentum[2] = global_values[11];

	total_collected_temperature[0] = global_values[12];
	total_collected_temperature[1] = global_values[13];
	total_collected_temperature[2] = global_values[14];

	INT64 global_num_total_particles   = global_values[15];



	if(!mpi_myRank)
	{

		ofstream Particle_ErgMomMass_FILE;

		ofstream OUT_FILE;


		char filename[200];

		//! Particle quantities
		sprintf(filename,"%s/MASS_ERG_TEMP/MassErgMom.txt",data_output_path);
		OUT_FILE.open(filename, ios_base::app);

		OUT_FILE << TL			       << "   ";
		OUT_FILE << total_particle_mass        << "   ";
		OUT_FILE << total_particle_energy      << "   ";
		OUT_FILE << total_particle_momentum[0] << "   ";
		OUT_FILE << total_particle_momentum[1] << "   ";
		OUT_FILE << total_particle_momentum[2] << "   ";
		OUT_FILE << endl;
		OUT_FILE.close();

		//! BField
		sprintf(filename,"%s/MASS_ERG_TEMP/BField.txt",data_output_path);
		OUT_FILE.open(filename, ios_base::app);

		OUT_FILE << TL				  << "   ";
		OUT_FILE << total_magnetic_field[0]       << "   ";
		OUT_FILE << total_magnetic_field[1]       << "   ";
		OUT_FILE << total_magnetic_field[2]       << "   ";
		OUT_FILE << vec_len(total_magnetic_field) << "   ";
		OUT_FILE << endl;
		OUT_FILE.close();

		//! UField
		sprintf(filename,"%s/MASS_ERG_TEMP/UField.txt",data_output_path);
		OUT_FILE.open(filename, ios_base::app);

		OUT_FILE << TL				      << "   ";
		OUT_FILE << total_collected_momentum[0]       << "   ";
		OUT_FILE << total_collected_momentum[1]       << "   ";
		OUT_FILE << total_collected_momentum[2]       << "   ";
		OUT_FILE << vec_len(total_collected_momentum) << "   ";
		OUT_FILE << endl;
		OUT_FILE.close();


		//! RhoField
		sprintf(filename,"%s/MASS_ERG_TEMP/Rho.txt",data_output_path);
		OUT_FILE.open(filename, ios_base::app);

		OUT_FILE << TL			    << "   ";
		OUT_FILE << total_collected_density << "   ";
		OUT_FILE << endl;
		OUT_FILE.close();

		//! TField
		sprintf(filename,"%s/MASS_ERG_TEMP/TemperatureField.txt",data_output_path);
		OUT_FILE.open(filename, ios_base::app);

		OUT_FILE << TL				   << "   ";
		OUT_FILE << total_collected_temperature[0] << "   ";
		OUT_FILE << total_collected_temperature[1] << "   ";
		OUT_FILE << total_collected_temperature[2] << "   ";
		OUT_FILE << total_collected_temperature[0] +total_collected_temperature[1] +total_collected_temperature[2] << "   ";
		OUT_FILE << endl;
		OUT_FILE.close();

		//! Split Merge
		sprintf(filename,"%s/MASS_ERG_TEMP/SplitMerge_sinceStart.txt",data_output_path);
		OUT_FILE.open(filename, ios_base::app);

		OUT_FILE << TL			  << "   ";

		OUT_FILE << num_split_since_start << "   ";
		OUT_FILE << num_merge_since_start << "   ";
		OUT_FILE << global_num_total_particles   << "   ";
		OUT_FILE << global_num_total_particles -num_split_since_start +num_merge_since_start << "   ";

		OUT_FILE << endl;
		OUT_FILE.close();



	}

}


//!--------------------------------------------------------
//! solve_Obstacle_BField:
//! - id_BEven is advanced 
//! - since dtB = cur(curlB) computed, this method should 
//!   work for CFBG-Field Scenarios as well.
//! - core is NOT taken into account 
//! - no domain decomposition support yet
//!--------------------------------------------------------
void CHybrid::solve_Obstacle_BField(void)
{

	log_file << " ADVANCING OBSTACLE BFIELD (CN) ..." << endl;

	clock_t start,finish;
	double time;
	start = clock();


	//! alloc array to store error in phi of each block
	//! LOCAL
	D_REAL *local_error_phi_L = new D_REAL[MAX_LEVEL +1];
	memset(local_error_phi_L, 0, (MAX_LEVEL +1)*sizeof(D_REAL));

	//! GLOBAL (all processes)
	D_REAL *global_error_phi_L = new D_REAL[MAX_LEVEL +1];
	memset(global_error_phi_L, 0, (MAX_LEVEL +1)*sizeof(D_REAL));

	//!--------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->calc_RHS_CN(id_rotB, id_BEven);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	FULL_GN_UPDATE(id_rotB);
	
	//! TODO:
	//! cancelling when error has been minimized does not work with
	//! MPI unless info is broadcasted to every process

	//!--------------------------------------------------
	for(INT32 counter=0; counter<=B_SOR_max_Iteration; counter++)
	{


		//! - itrrate in each level, get boundary values from lower levels
		//! - do not inject solution (may not match up with divB in respective level)
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{

			//! reset local error of respective level
			local_error_phi_L[level] = 0.;


			//! - iterate more often in higher level since it takes more steps on finer meshes
			//! - using powers of two seems to converge phi of different levels equally fast
		      #if defined(use_vectorclass)
			for(INT32 z=0; z<INT32(pow(double(B_SOR_num_cycle_base),double(level))); z++)
		      #else
			for(INT32 z=0; z<INT32(pow(B_SOR_num_cycle_base,level)); z++)
		      #endif
			{

				CBlock* temp_Block = BlockList_of_Lev[level];
				while(temp_Block)
				{

					//! iterate BEven at EVEN cell indices and store error of respective block
					if(mpi_myRank == temp_Block->responsible_mpi_process)
					local_error_phi_L[level] += temp_Block->Obstacle_BField_SOR_Step(counter, id_BEven, id_rotB, EVEN_ITERATION);

		
					temp_Block = temp_Block->next_Blk_of_BlockList;
				}


				update_GN_of_level(level, id_BEven, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);

				temp_Block = BlockList_of_Lev[level];
				while(temp_Block)
				{
	
					//! iterate BEven at ODD cell indices and store error of respective block
					if(mpi_myRank == temp_Block->responsible_mpi_process)
					local_error_phi_L[level] += temp_Block->Obstacle_BField_SOR_Step(counter, id_BEven, id_rotB, ODD_ITERATION);

		
					temp_Block = temp_Block->next_Blk_of_BlockList;
				}


				update_GN_of_level(level, id_BEven, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);

			}
		}

		if(counter%B_SOR_calc_error_step==0)
		{

			D_REAL error_all_levels = 0.;

			mpi_build_sum(local_error_phi_L,
				      global_error_phi_L,
				      MAX_LEVEL +1);

			if(num_optimize_SOR_steps)
			{	
				if(TL==TLstart_optimize_B_SOR_omega)
				{
					BSOR_opti_array = new D_REAL*[MAX_LEVEL+1];
					for(INT32 level=0; level<MAX_LEVEL+1; level++)
					BSOR_opti_array[level] = new D_REAL[num_optimize_SOR_steps];
					
					for(INT32 level=0; level<MAX_LEVEL+1; level++)
					 for(INT32 steps=0; steps<num_optimize_SOR_steps; steps++)
					  BSOR_opti_array[level][steps]=0;					
				}

				//! free memory here in next time steps and not at end of function
				//! because of multiple BField loops per TL
/*				if(TL==TLstart_optimize_B_SOR_omega+num_optimize_SOR_steps+1)
				{	
					for(INT32 level=0; level<MAX_LEVEL+1; level++)
					delete[] BSOR_opti_array[level];

					delete[] BSOR_opti_array;	
				}*/	
				
				if(counter==B_SOR_calc_error_step)
				 if(TL>=TLstart_optimize_B_SOR_omega && TL<TLstart_optimize_B_SOR_omega + num_optimize_SOR_steps)	 
				  for(INT32 level=0; level<=MAX_LEVEL; level++)	
				  BSOR_opti_array[level][TL-TLstart_optimize_B_SOR_omega] = global_error_phi_L[level];	   
			}
		
			///! protocol errors and sum up errors of level
			log_file << "  iteration step: " << counter << endl;
			for(INT32 level=0; level<=MAX_LEVEL; level++)
			{
				log_file << "  Level:  " << level << "	";
				log_file << "  error (global): " << global_error_phi_L[level] << endl;
				error_all_levels += global_error_phi_L[level];
			}

			log_file << "  error_all_levels: " << error_all_levels << endl;
			log_file  << endl;


			//! check whether error is small enough
			if(error_all_levels < B_SOR_max_error)
			{
				log_file << "  error < " << B_SOR_max_error << endl;
				log_file << "  exiting iteration ..." << endl;
				counter = B_SOR_max_Iteration+1;
				
			}



		}


	}


	//! now id_BEven is update,
	//! build total BField from:
	//! B_total = B_even + B_cfbg needed in:
	//! - EField Calculation
	//! - CAM
	//! - particle Acceleration
	build_B_total();

	delete[] local_error_phi_L;
	delete[] global_error_phi_L;

	
	if(num_optimize_SOR_steps)
	{	
		if(TL==TLstart_optimize_B_SOR_omega+num_optimize_SOR_steps)
		{
			D_REAL best_B_SOR_error[MAX_LEVEL+1];
			D_REAL best_B_SOR_value[MAX_LEVEL+1];
			
			//! find best B_SOR
			for(INT32 level=0; level<=MAX_LEVEL; level++)	
			{	
				best_B_SOR_error[level] = 1.e+06; 
				
				for(INT32 step=0;step<num_optimize_SOR_steps;step++)
				if(BSOR_opti_array[level][step]<best_B_SOR_error[level] && BSOR_opti_array[level][step]>0)
				{	
					best_B_SOR_error[level]=BSOR_opti_array[level][step];
					best_B_SOR_value[level]=B_SOR_omega_L[level] - 0.01*num_optimize_SOR_steps/2+0.01*step;
				}	
			
				log_file << " OPTIMAL B_SOR[level="<<level<<"] = "<< best_B_SOR_value[level] <<"   (B_SOR_error = ";
				log_file << best_B_SOR_error[level]<<")"<<endl;
				
			}	
			
		}	 
		
		
	}
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " Obstacle BField time: " << time << "s." << endl << endl;



}

//!--------------------------------------------------------
//!- clean_divergence 
//! - id_BEven is cleaned
//! - core is NOT taken into account 
//! - no domain decomposition support yet
//!--------------------------------------------------------
void CHybrid::clean_divergence(void)
{

	log_file 	<< " CLEANING DIVERGENCE ..."  << endl;


	clock_t start,finish;
	double time;
	start = clock();



	calc_div(id_BEven, id_divB);
	FULL_GN_UPDATE(id_divB);



	//! alloc array to store error in phi of each block
	//! LOCAL
	D_REAL *local_error_phi_L = new D_REAL[MAX_LEVEL +1];
	memset(local_error_phi_L, 0, (MAX_LEVEL +1)*sizeof(D_REAL));

	//! GLOBAL (all processes)
	D_REAL *global_error_phi_L = new D_REAL[MAX_LEVEL +1];
	memset(global_error_phi_L, 0, (MAX_LEVEL +1)*sizeof(D_REAL));
	

	//!--------------------------------------------------
	for(INT32 counter=0; counter<=DC_max_Iteration; counter++)
	{




		//! - itrrate in each level, get boundary values from lower levels
		//! - do not inject solution (may not match up with divB in respective level)
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{

			//! - iterate more often in higher level since it takes more steps on finer meshes
			//! - using powers of two seems to converge phi of different levels equally fast
		      #if defined(use_vectorclass)
			//! Compute z_max = pow(DC_num_cycle_base,level)
			INT32 z_max;
			if(DC_num_cycle_base==2)
			{
			  z_max = 1 << level;
			}
			else
			{
			  z_max = 1;
			  for(INT32 i=0; i<level; i++)
			  {
			    z_max *= DC_num_cycle_base;
			  }
			}
			for(INT32 z=0; z<z_max; z++)
		      #else
			for(INT32 z=0; z<INT32(pow(DC_num_cycle_base,level)); z++)
		      #endif
			{

				//! reset local error of respective level
				//! -> only consider error of last z iterarion 
				local_error_phi_L[level] = 0.;


				//! EVEN ITERATION
				CBlock* temp_Block = BlockList_of_Lev[level];
				while(temp_Block)
				{
					
					//! iterate phi at EVEN cell indices and store error of respective block
					if(mpi_myRank == temp_Block->responsible_mpi_process)
					local_error_phi_L[level] += temp_Block->Psi_DC_SOR_Step(counter, id_divB, EVEN_ITERATION);


					//! since psi values are updated "in Place", need to
					//! copy values within loop
					temp_Block = temp_Block->next_Blk_of_BlockList;
				}
	
				if(z==0)
				update_GN_of_level(level, id_PhiDC, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
				else
				update_GN_of_level(level, id_PhiDC, SKIP_PHYS_FACE, NO_PARENT_BUFFER_UPDATE);

				//! ODD ITERATION
				temp_Block = BlockList_of_Lev[level];
				while(temp_Block)
				{
					
					//! iterate phi at ODD cell indices and store error of respective block
					if(mpi_myRank == temp_Block->responsible_mpi_process)
					local_error_phi_L[level] += temp_Block->Psi_DC_SOR_Step(counter, id_divB, ODD_ITERATION);


					//! since psi values are updated "in Place", need to
					//! copy values within loop
					temp_Block = temp_Block->next_Blk_of_BlockList;
				}
	
				if(z==0)
				update_GN_of_level(level, id_PhiDC, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
				else
				update_GN_of_level(level, id_PhiDC, SKIP_PHYS_FACE, NO_PARENT_BUFFER_UPDATE);

	
			}

		}



		if(counter%DC_calc_error_step==0)
		{


			D_REAL error_all_levels = 0.;

			mpi_build_sum(local_error_phi_L,
				      global_error_phi_L,
				      MAX_LEVEL +1);

			if(num_optimize_SOR_steps)
			{
				if(TL==TLstart_optimize_B_SOR_omega)
				{
					DC_omega_opti_array = new D_REAL*[MAX_LEVEL+1];
					for(INT32 level=0; level<MAX_LEVEL+1; level++)
					DC_omega_opti_array[level] = new D_REAL[num_optimize_SOR_steps];
					
					for(INT32 level=0; level<MAX_LEVEL+1; level++)
					 for(INT32 steps=0; steps<num_optimize_SOR_steps; steps++)
					  DC_omega_opti_array[level][steps]=0;					
				}
				
				//! free memory here in next time steps and not at end of function
				//! because of multiple BField loops per TL
/*				if(TL==TLstart_optimize_B_SOR_omega+num_optimize_SOR_steps+1)
				{	
					for(INT32 level=0; level<MAX_LEVEL+1; level++)
					delete[] DC_omega_opti_array[level];

					delete[] DC_omega_opti_array;	
				}*/	
				
				

				if(counter==DC_calc_error_step)
				{					
				  if(TL>=TLstart_optimize_B_SOR_omega && TL<TLstart_optimize_B_SOR_omega + num_optimize_SOR_steps)	 
				   for(INT32 level=0; level<=MAX_LEVEL; level++)	
				    DC_omega_opti_array[level][TL-TLstart_optimize_B_SOR_omega] = global_error_phi_L[level];	
	
				} 
			} 

			///! protocol errors and sum up errors of level
			log_file << "  iteration step: " << counter << endl;
			for(INT32 level=0; level<=MAX_LEVEL; level++)
			{
				log_file << "  Level:  " << level << "	";
				log_file << "  error (global): " << global_error_phi_L[level] << endl;
				error_all_levels += global_error_phi_L[level];
			}

			log_file << "  error_all_levels: " << error_all_levels << endl;
			log_file  << endl;


			//! check whether error is small enough
			if(error_all_levels < DC_max_error)
			{
				log_file << "  error < " << DC_max_error << endl;
				log_file << "  exiting iteration ..." << endl;
				counter = DC_max_Iteration+1;
				
			}
		}

	}



	//!--------------------------------------------------
	 for(INT32 level=0; level<=MAX_LEVEL; level++)
	 {
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->project_divB_out_of_BField(id_BEven);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}

	FULL_GN_UPDATE(id_BEven);

	//! now id_BEven is update,
	//! build total BField from:
	//! B_total = B_even + B_cfbg needed in:
	//! - EField Calculation
	//! - CAM
	//! - particle Acceleration
	build_B_total();


	delete[] local_error_phi_L;
	delete[] global_error_phi_L;

	if(num_optimize_SOR_steps)
	{	
		if(TL==TLstart_optimize_B_SOR_omega+num_optimize_SOR_steps)
		{
			D_REAL best_DC_omega_error[MAX_LEVEL+1];
			D_REAL best_DC_omega_value[MAX_LEVEL+1];
			
			//! find best B_SOR
			for(INT32 level=0; level<=MAX_LEVEL; level++)	
			{	
				best_DC_omega_error[level] = 1.e+06; 
				
				for(INT32 step=0;step<num_optimize_SOR_steps;step++)		
				if(DC_omega_opti_array[level][step]<best_DC_omega_error[level] && DC_omega_opti_array[level][step]>0)
				{	
					best_DC_omega_error[level]=DC_omega_opti_array[level][step];
					best_DC_omega_value[level]=DC_omega_L[level]-0.01*num_optimize_SOR_steps/2+0.01*step;
	
				}
				
				log_file << " OPTIMAL DC_omega[level="<<level<<"] = "<< best_DC_omega_value[level] <<"   (DC_omega_error = ";
				log_file << best_DC_omega_error[level]<<")"<<endl;
				
			}	
			
		}	
	}		
	
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " cleaning divergence time: " << time << "s." << endl << endl;

}



//!------------------------------------------------------------------------
//!- calc_rotB:
//!------------------------------------------------------------------------
void CHybrid::calc_rot(INT32 BField_type, INT32 rot_Field_type)
{


	log_file << "  calc rotB...       ";
	
	//! in order to advance odd BField derivatives of
	//! even BField are necessary. So first cp even GC.
	FULL_GN_UPDATE(id_BEven);
	
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->calc_rot(BField_type, rot_Field_type);

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	FULL_GN_UPDATE(rot_Field_type);




// 	for(INT32 level=MAX_LEVEL; level>0; level--)
// 	{
// 	
// 		CBlock* temp_Block = BlockList_of_Lev[level];
// 		while(temp_Block)
// 		{
// 
// 			if(mpi_myRank == temp_Block->responsible_mpi_process)
// 			temp_Block->field_to_parent_smooth(rot_Field_type,rot_Field_type);
// 
// 			temp_Block = temp_Block->next_Blk_of_BlockList;
// 		}
// 	}
// 
// 	FULL_GN_UPDATE(rot_Field_type);

	
	log_file << "done." << endl;

}


//!------------------------------------------------------------------------
//!- calc_rotB:
//!------------------------------------------------------------------------
void CHybrid::calc_ExB_Drift(INT32 id_Result, INT32 id_Field1, INT32 id_Field2)
{


	log_file << "  calc "<< Field_Name[id_Result]
		      <<" = "<< Field_Name[id_Field1]
		      <<" x "<< Field_Name[id_Field2];
	
	//! in order to advance odd BField derivatives of
	//! even BField are necessary. So first cp even GC.

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{

				D_REAL *R1 = temp_Block->Field_Type[id_Result] +0*num_nodes_in_block;
				D_REAL *R2 = temp_Block->Field_Type[id_Result] +1*num_nodes_in_block;
				D_REAL *R3 = temp_Block->Field_Type[id_Result] +2*num_nodes_in_block;

				D_REAL *A1 = temp_Block->Field_Type[id_Field1] +0*num_nodes_in_block;
				D_REAL *A2 = temp_Block->Field_Type[id_Field1] +1*num_nodes_in_block;
				D_REAL *A3 = temp_Block->Field_Type[id_Field1] +2*num_nodes_in_block;

				D_REAL *B1 = temp_Block->Field_Type[id_Field2] +0*num_nodes_in_block;
				D_REAL *B2 = temp_Block->Field_Type[id_Field2] +1*num_nodes_in_block;
				D_REAL *B3 = temp_Block->Field_Type[id_Field2] +2*num_nodes_in_block;


				for(INT32 elm=0; elm<num_nodes_in_block; elm++)
				{

					D_REAL R[3] = {     0.,      0.,      0.};
					D_REAL A[3] = {A1[elm], A2[elm], A3[elm]};
					D_REAL B[3] = {B1[elm], B2[elm], B3[elm]};

					vec_cross(R, A, B);

					D_REAL rez_B2 = 1./vec_len2(B);

					R1[elm] = rez_B2*R[0];
					R2[elm] = rez_B2*R[1];
					R3[elm] = rez_B2*R[2];
				}

			}

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	
	log_file << " done." << endl;

}

//!--------------------------------------------------------
//!- calc_local_gyro_radius:
//!--------------------------------------------------------
void CHybrid::calc_local_gyro_radius(INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField)
{

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->calc_local_gyro_radius(id_Gyro_Radius, species, id_BField, id_UField, SI_x0/1000);

	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------
//!- calc_local_electron_gyro_radius:
//!--------------------------------------------------------
void CHybrid::calc_local_electron_gyro_radius(INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField)
{
    
    for(INT32 level=0; level<=MAX_LEVEL; level++)
    {
        CBlock* temp_Block = BlockList_of_Lev[level];
        while(temp_Block)
        {
            
            
            if(mpi_myRank == temp_Block->responsible_mpi_process)
                temp_Block->calc_local_electron_gyro_radius(id_Gyro_Radius, species, id_BField, id_UField,SI_x0/1000);
            
            
            temp_Block = temp_Block->next_Blk_of_BlockList;
        }
    }
    
}


//!------------------------------------------------------------------------
//!- calc_grad:
//!------------------------------------------------------------------------
void CHybrid::calc_grad(INT32 in_Field_type, INT32 grad_Field_type)
{


	log_file << "  calc grad " << Field_Name[in_Field_type] << "...       ";
	
	
	//! in order to advance odd BField derivatives of
	//! even BField are necessary. So first cp even GC.
	FULL_GN_UPDATE(in_Field_type);
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->calc_grad(in_Field_type, grad_Field_type);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}
	


	log_file << "done." << endl;

}



//!------------------------------------------------------------------------
//!- calc_divB:
//!------------------------------------------------------------------------
void CHybrid::calc_div(INT32 Field_type, INT32 div_Field_type)
{


	log_file << "  calc " << Field_Name[div_Field_type] << "...       ";
	
	
	//! in order to advance odd BField derivatives of
	//! even BField are necessary. So first cp even GC.
	FULL_GN_UPDATE(id_BEven);
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->calc_div(Field_type, div_Field_type);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}
	
	log_file << "done." << endl;

}


//!--------------------------------------------------------------
//!- smooth_Field:
//!--------------------------------------------------------------
void CHybrid::smooth_Field(INT32 type, const D_REAL smooth_value[])
{

	//! return in case all smooth values are zero
	bool do_smooth = false;
	for(INT32 level=0; level<= MAX_LEVEL; level++)
	if(smooth_value[level] > 0.) do_smooth = true;

	if(!do_smooth)
	return;


	log_file <<"  smoothing Field "<< Field_Name[type] << " ...  " << endl  << "   ";

	clock_t start,finish;
	double time;
	start = clock();

	log_file <<"|";
	
	if(several_smooth_in_levels)
	{
#if defined(use_vectorclass)
		for(INT32 level=0; level<= MAX_LEVEL; level++)
		{
		  log_file << "L" << level << ": " << pow(2.,double(level)) << "x"<< smooth_value[level] <<"| |";
		}
		log_file << endl;
#else
		for(INT32 level=0; level<= MAX_LEVEL; level++)
		log_file << "L" << level << ": " << pow(2.,level) << "x"<< smooth_value[level] <<"| |";
		log_file << endl;
#endif
	}


     	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		INT32 num_smooth = 1;
#if defined(use_vectorclass)
		if(several_smooth_in_levels)
		{
		  num_smooth = int( pow(2.,double(level)) );
		}
#else
		if(several_smooth_in_levels)
		num_smooth = int(pow(2.,level));
#endif
	

		for(INT32 smooth_step=0; smooth_step<num_smooth; smooth_step++)
		{


			if(!several_smooth_in_levels)
			log_file << "L" << level << ": " << smooth_value[level] <<"| |";


			//! TODO:
			//! maybe update is required here as well in order to
			//! get values from lower level before this is smoothed.

			//! smooth field
			CBlock* temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{

				if(mpi_myRank == temp_Block->responsible_mpi_process)
				temp_Block->smooth_Field(type, smooth_value[level]);
			
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}



			//! update GN in level before smooth to parent
			//! (applying this after S2P does work equally well)
			update_GN_of_level(level, type, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);


			smooth_field_to_parent(level, type);


		}//! end for num smooth
	}//! end for level

  if(!several_smooth_in_levels)
  log_file << endl;


  finish = clock();
  time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  log_file << "  smooth time: " << time << "s." << endl;

}

//!------------------------------------------------------------
//! negative particles:
//!
//! - Check how many cells with rho < 0
//! - delete negative particle in those cells
//!------------------------------------------------------------
void CHybrid::negative_particles(void)
{
	//! check at TL 0 if there are negative charges 
	if(!TL || TL%(TL_SAVE_STATE+1)==0 )
	{  
	 for(INT32 species=0; species<num_Charged_Species; species++)
	  if(Ion_Charges[species]<0)
	    negative_particles_in_Simu = true;
	  
	    return;
	}  
  
	log_file << " NEGATIVE IONS/DUST: " <<endl;

	//! this parameter should NOT be neccessary
// 	stopp_inserting_negative_ions=false;

	clock_t start,finish;
	double time;
	start = clock();

	//! array for counting negative cells (index 0)
	//! and number of deleted particles (index 1)
	//! and deleted particles in respective level(index 2 + level)
	INT64* num_negCells_delPart = new INT64[2+MAX_LEVEL + 1]; 
	memset(num_negCells_delPart,0,(2 + MAX_LEVEL + 1)*sizeof(INT64));

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->negative_particles(num_negCells_delPart, level);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	//! Provide Information about how many cells with negative charge were found
	//! and how many negative particles have been deleted
	INT64 local_info_values[2];
	stringstream info_names[INFO_ARRAY_SIZE];
	INT32 counter=2;
	local_info_values[0] = num_negCells_delPart[0];
	local_info_values[1] = num_negCells_delPart[1];
	info_names[0] << "   -> Num cells with rho+ < rho- : ";
	info_names[1] << "   -> Negative Particles deleted : ";

	//! Reduce via MPI if required
	show_information(local_info_values,
		info_names,
		counter,
		BUILD_SUM);

// 	if(local_info_values[0]>max_number_of_cells_with_negative_rho)
// 	stopp_inserting_negative_ions=true;

	//! update statistics
	num_total_particles += -num_negCells_delPart[1];
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		num_total_particles_in_L[level] += - num_negCells_delPart[2 + level];
	}

	delete[] num_negCells_delPart;

	//! measure time
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " time for negative particles: " << time << "s." << endl << endl;

}


//!------------------------------------------------------------
//! collect rho / J :
//!
//! - decide if only marked particle shall be collected
//!------------------------------------------------------------
void CHybrid::collect_rho_Ji(INT32 species_to_collect,
			     INT32 id_rho_type,
			     INT32 id_Ji_type,
			     bool skip_unmarked_particle)
{


	log_file 	<< " GATHERING " << Field_Name[id_rho_type] << ": "  << endl;
	log_file 	<< " GATHERING " << Field_Name[id_Ji_type] << ": "  << endl;

	num_total_particles_collected = 0;
	num_particles_collected_from_parent = 0;

	clock_t start,finish;
	double time;
	start = clock();

	//! Collecting are two loop over all levels:
	//! 
	//! 1) collect in every Block which has no children 
	//!    and inject collected values down to parent
	//!
	//! 2) 


	//! First add, then devide !!!

	//! 1)  collect all rho_species & Ji_total & Gam total
	//! 2)      add all rho_species & Ji_total & Gam total to neighbors
	//! 3) exchange all rho_species & Ji_total & Gam total Ghost Cells
	//! 4) in one loop over all blocks:
	//!	 a) build rho total
	//!	 b) buils Lam
	//!	 c) calc Ji2Ui


	//!NOTE:
	//! Do not delete Moments of entire Block directly before it is processed, 
	//! else moments that have been injected by children would be deleted
	//! Rather set zero ALL blocks BEFORE gather procedure
	set_zero_field_incGather(id_rho_type);
	set_zero_field_incGather(id_Ji_type);


	//!------------------------------------------------------------
	//!-- 1) GATHER LOOP ------------------------------------------
	//!------------------------------------------------------------
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{

		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{


				//! gather Block
				temp_Block->collect_rho_Ji(species_to_collect,
							   id_rho_type,
							   id_Ji_type,
							   skip_unmarked_particle);



				//! simple inject to parent in case gather blocks are
				//! not used
				if(!use_gather_blocks && level && !temp_Block->is_gatherBlk)
				{
					temp_Block->add_field_to_parent(id_rho_type,
							         	id_rho_type);

					temp_Block->add_field_to_parent(id_Ji_type,
							         	id_Ji_type);
				}

			}
			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}


		//! send injected field to parent
		if(!use_gather_blocks)
		{
			zero_parent_add_children_field_MPI(level, id_rho_type);
			zero_parent_add_children_field_MPI(level, id_Ji_type);
		}

		//! now gathering in level is finished, add boundary moments
		Add_Boundary_Moments_of_L(level, id_rho_type);
		Add_Boundary_Moments_of_L(level, id_Ji_type);
	}


	//!------------------------------------------------------------
	//!-- 2) UPDATE GN --------------------------------------------
	//!------------------------------------------------------------
	if(use_gather_blocks)
	{
		moments_to_parent_smooth(id_rho_type);
		moments_to_parent_smooth(id_Ji_type);

	}
	else
	{
		FULL_MOMENT_UPDATE(id_rho_type);
		FULL_MOMENT_UPDATE(id_Ji_type);
	}

	//!------------------------------------------------------------
	//!-- 3) CONVERT CURRENTS AND DENSITIES TO VELOCITIES ---------
	//!--    and inject to parent			      ---------
	//!------------------------------------------------------------
	for(INT32 level=MAX_LEVEL; level>= 0; level--)
	{

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->convert_Ji2Ui(0,
						  id_rho_type,
						  id_Ji_type,
						  id_Ji_type);
			

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}




	INT64 local_info_values[3] = {num_total_particles,
				      num_total_particles_collected,
				      num_particles_collected_from_parent};


	stringstream info_names[3];
	info_names[0] << "   num_total_particles: ";
	info_names[1] << "   particles collected: ";
	info_names[2] << "   particles collected from parent:";

	show_information(local_info_values,
			 info_names,
			 3, BUILD_SUM);


	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " gather time: " << time << "s." << endl << endl;



}

//!
//! prepare_Recombination_Density()
//!

//! The function will copy the field 
//! id_rho_np1
//! to 
//! id_rho_np1_recombined

//! The function is needed by calc_Recombination_Density()!

void CHybrid::prepare_Recombination_Density()
{
        //!-----------------------------------------//
        //!     Loop over all Blks                       //
        //!----------------------------------------//
        for(INT32 level=0; level<=MAX_LEVEL; level++)
        {
                CBlock* temp_Block = BlockList_of_Lev[level];
                while(temp_Block)
                {

                if(mpi_myRank==temp_Block->responsible_mpi_process)
                        temp_Block->prepare_Recombination_Density();
                        temp_Block = temp_Block->next_Blk_of_BlockList;
                }
        }
}


//! 
//! calc_Recombination_Density()
//! 

//! The function will calculate the loss of density 
//! by subtracting the fields
//! id_rho_np1 and id_rho_np1_recombined

void CHybrid::calc_Recombination_Density()
{
        //!-----------------------------------------//
        //!     Loop over all Blks                       //
        //!----------------------------------------//
        for(INT32 level=0; level<=MAX_LEVEL; level++)
        {
                CBlock* temp_Block = BlockList_of_Lev[level];
                while(temp_Block)
                {

                if(mpi_myRank==temp_Block->responsible_mpi_process)
                        temp_Block->calc_Recombination_Density();
                        temp_Block = temp_Block->next_Blk_of_BlockList;
                }
        }
}


//!------------------------------------------------------------------------
//!- average reaction rate
//!------------------------------------------------------------------------
void CHybrid::calc_reaction_rate_field(INT32 dest, INT32 neutral_species)
{


	log_file << "  Calc mean reaction rate to field "<< Field_Name[dest];
		   
	D_REAL reaction_Rate=0;
	
	D_REAL *rho_species, *rho_species_array, *nu_mean, *rho, *u_speciesX; 	
	

	for(short species=0; species<num_Inflow_Species; species++)
	{	
		
		//! use id_ForceSpecies1 as temp scratch 
		collect_Ui_vth2_of_species(species, id_ForceSpecies1, 0, noVTH2);		
						
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{
		
			CBlock* temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{

				if(mpi_myRank == temp_Block->responsible_mpi_process)
				{	
		
					//! total ion density
					rho = temp_Block->Field_Type[id_rho_np1];					
						
					//! set pointer to respective species
					rho_species = temp_Block->Field_Type[id_rhoSpecies1 +species];					

					//! set pointers to ion velocity (in scratch field)
					u_speciesX = temp_Block->Field_Type[id_ForceSpecies1];
					D_REAL* u_speciesY = u_speciesX + 1*num_nodes_in_block;
					D_REAL* u_speciesZ = u_speciesX + 2*num_nodes_in_block;
					
					//! set pointer to nu_mean field
					nu_mean = temp_Block->Field_Type[dest];
					
					for(INT32 elm=0; elm<num_nodes_in_block; elm++)
					{
						PARTICLE_REAL vec_u[3] = {u_speciesX[elm],u_speciesY[elm],u_speciesZ[elm]};
						

						//! calc total rate for all reactions of species										
						reaction_Rate=0;
						
						for(INT32 destSpec=0;destSpec<num_Charged_Species;destSpec++)
						reaction_Rate += calc_reaction_rate(vec_u,species,destSpec,neutral_species);
						
						nu_mean[elm] += rho_species[elm]*reaction_Rate/rho[elm];

											
					}
						
						
				}
				
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}
		}
	}
	
	log_file << " done." << endl;

}

//!-----------------------------------------------------------
//!-----------------------------------------------------------
void CHybrid::show_gradPE(INT32 field_id)
{

	for(INT32 level=MAX_LEVEL; level>= 0; level--)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->show_gradPE(field_id);	
	
			temp_Block = temp_Block->next_Blk_of_BlockList;		
		}

	}
	
	//! UPDATE GN --------------------------------------------
	if(use_gather_blocks)
	moments_to_parent_smooth(field_id);
	else
	FULL_MOMENT_UPDATE(field_id);
	
}	
