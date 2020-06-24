

#include "defines.h"
#include "CHybrid.h"
#include "CHybrid.h"
#include "parameters.h"
#include "absolute_Globals.h"


#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>



using namespace std;

CHybrid* Hybrid;



//!-------------------------------------------------------------//
//! main:
//! - Acces routine of programm.
//! - Hybrid Class is instanced and Hybrid Cycle is performed
//!   within this routine
//!-------------------------------------------------------------//
int main(int argc, char *argv[])
{



	//! Allocate Instance of Hybrid Class
	Hybrid = new CHybrid;

	//! init MPI, get myRnak and num_processes
   	Hybrid->init_MPI(argc, argv);

	//! read parameter, alloc memory, restore state
   	Hybrid->init();
	log_file<<"init complete"<<endl;
	//! --------------------------------------------
	//! --- Initialize and collect particle --------
	//! - (some initialization steps are left out)--
	//! --------------------------------------------
	Hybrid->fill_Particle_in_Box();
// 	if(TL==0)
// 	for( int i=0;i<1000; i++ )
// 	Hybrid->inject_obstacle_ions();

	Hybrid->collect_RHOnp1_UIplus_LAM_GAM();
        
	Hybrid->copy_Field(id_rho_n, id_rho_np1);
        
	
#if defined nonadiabatic_gradPE_TERM
        if(TL==0)
        {
            Hybrid->init_PE_Profile();
        }
#endif
	
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
	Hybrid->init_Neutral_Profile();
#endif	

#if defined(use_ion_production_as_field)
	Hybrid->init_IonProduction_Profile();
#endif	
	
#if defined(use_dust_species_as_field)
	Hybrid->init_Dust_Profile();
#endif		
	
// 	Hybrid->add_force(id_EField);

	Hybrid->negative_particles();


	Hybrid->output_all(TL);


	//! syncronize all process before time is measured,
	//! else conflict when writing state file may appear.
	Hybrid->synchronize_allProcesses_MPI();
	time(&Hybrid->run_start_phys);
	Hybrid->run_start_comp = clock();

	Hybrid->synchronize_allProcesses_MPI();
	Hybrid->global_MPI_lastState_time = MPI_Wtime();
	Hybrid->run_start_mpi = Hybrid->global_MPI_lastState_time;




	//! Time Loop defined as in TB 2004 PhD thesis pp. 27
	for(; TL<TL_MAX;)
	{


		TL++;
		
		//ionization routine, for ionospheres that build slowly
		//if switch not activated calculate everything 
		if (!ionizationonly || TL>TLIO){
		  //! measure time for statistics
		  Hybrid->measure_time();

		  
		  //! --- MAIN HYBRID CYCLE ---
		  //! 1.)
		  if(run_calc_first_E) {
			  Hybrid->calc_first_E();
		  }
		  
		  //! 2.)
		  if(run_CAM) {
			  Hybrid->CAM();
		  }

		  //! 3.)
		  if(run_calc_second_E) {
			  Hybrid->calc_second_E();
		  }

  // 		Hybrid->add_force(id_EField);

		  //! 4.)
		  if(run_accelerate_Particle) {
			  Hybrid->accelerate_Particle();
		  }

		  //! 5.) 
		  if(run_collect_Ui_minus) {
			  Hybrid->collect_Ui_minus();
		  }

		  //! 6a.)
		  if(run_move_Particle){
			  Hybrid->move_Particle();
		  }

		  //! 6b)
		  if(run_Split_Merge_Particle){
			  Hybrid->Split_Particle();
			  Hybrid->Merge_Particle();
		  }
		  
		  //! scheint zu Abstürzen auf SWARM zu führen:
		  //! 20 < 20 ... exiting ... ?!?
  // 		Hybrid->check_weight_sorting();

		  //! 6b2)
		  if(run_negative_Particles){
			  if(negative_particles_in_Simu)
				  Hybrid->negative_particles();		
		  }
		  

		  
		  if(run_collect_Rho_prepare_Recombination_Density){
			  if(!TestParticle_Simulation){
				  Hybrid->collect_RHOnp1_UIplus_LAM_GAM();
				  Hybrid->prepare_Recombination_Density();
			  }
		  }
		  
		  //! 6d.)
		  if(run_chemical_Reactions){
			  if(!TestParticle_Simulation) {
				  if(TL!=0 && TL !=1) {                
					  Hybrid->chemical_Reactions();
				  }
			  }
		  }
		  
		  //! 6c.)
		  if(run_inject_obstacle_ions){
			  Hybrid->inject_obstacle_ions();
		  }	
		  
		  //! 6f.)
		  if(run_resize_pArrays){
			  Hybrid->resize_pArrays();
		  }
		  
		  //! 7.)
		  if(run_collect_Rho_calc_Recombination_Density){
			  if(!TestParticle_Simulation){
				  Hybrid->collect_RHOnp1_UIplus_LAM_GAM();
				  Hybrid->calc_Recombination_Density();
			  }
		  }
		  

		  //! 8.)
		  if(run_average_Ui_rho_setREZrho){
			  Hybrid->average_Ui_rho_setREZrho();
		  }

		  //! subcycle plasma and obstacle B if required
		  if(run_advanceB){
			  for(INT32 loop=0; loop<num_advance_B_loops; loop++)
			  {
				  //! 9.a)
				  Hybrid->advanceB_Plasma();

				  //! 9.b)
				  Hybrid->advanceB_Obstacle();
			  }
		  }


		  //! --- DATA ADMINISTRATION ---

		  //! A.) Redistribute Blocks to processes
		  Hybrid->redistribute_blocks(TL);

		  //! B.) Refine Mesh
		  Hybrid->refine_Mesh(TL);

		  //! C.) Output Data
		  Hybrid->output_all(TL);

		  //! D.) Reset Block timing
		  Hybrid->reset_block_timing();

		  //! E.) Build Statistic
		  Hybrid->statistic();

		  //! F.) Save State
		  if(TL_SAVE_STATE && TL_SAVE_STATE<=TL_MAX) 
		  {
			  if( !(TL % TL_SAVE_STATE ==0) )
			  log_file << " Next save state in "<< TL_SAVE_STATE -TL%TL_SAVE_STATE <<" TL." << endl;
			  else
			  Hybrid->save_state();
		  }				
		  
		  //! G.) Resub job
		  Hybrid->resubjob();
		}
		
		
		//if switch activated, calculate only ionospheric and statistic stuff
		else{
		  
		  Hybrid->measure_time();
		  if(run_collect_Rho_prepare_Recombination_Density){
			  if(!TestParticle_Simulation){
				  Hybrid->collect_RHOnp1_UIplus_LAM_GAM();
				  Hybrid->prepare_Recombination_Density();
			  }
		  }
		  
		  //! 6d.)
		  /*if(run_chemical_Reactions){
			  if(!TestParticle_Simulation) {
				  if(TL!=0 && TL !=1) {                
					  Hybrid->chemical_Reactions();
				  }
			  }
		  }*/
		  
		  //! 6c.)
		  if(run_inject_obstacle_ions){
			  Hybrid->inject_obstacle_ions();
		  }	
		  
		  //! 6f.)
		  if(run_resize_pArrays){
			  Hybrid->resize_pArrays();
		  }
		  //! --- DATA ADMINISTRATION ---

		  //! A.) Redistribute Blocks to processes
		  Hybrid->redistribute_blocks(TL);

		  //! B.) Refine Mesh
		  Hybrid->refine_Mesh(TL);

		  //! C.) Output Data
		  Hybrid->output_all(TL);

		  //! D.) Reset Block timing
		  Hybrid->reset_block_timing();

		  //! E.) Build Statistic
		  Hybrid->statistic();

		  //! F.) Save State
		  if(TL_SAVE_STATE && TL_SAVE_STATE<=TL_MAX) 
		  {
			  if( !(TL % TL_SAVE_STATE ==0) )
			  log_file << " Next save state in "<< TL_SAVE_STATE -TL%TL_SAVE_STATE <<" TL." << endl;
			  else
			  Hybrid->save_state();
		  }				
		  
		  //! G.) Resub job
		  Hybrid->resubjob();
		}

	}

	log_file << " Simulation finished at TL "<<TL<<endl;

	//! finish MPI process
	Hybrid->finalize_MPI();


}

