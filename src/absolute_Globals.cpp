

	#include <iostream>
	#include <fstream>
	
	#include "defines.h"
	#include "parameters.h"
	#include "gsl/gsl_rng.h"
	#include "absolute_Globals.h"



	using namespace std;
	
	
	//! MPI related globals
	MPI::Intracomm mpi_myComm;
	INT32 mpi_myRank, mpi_num_processes;
// 	INT32 *num_rootBlocks_at_process;
	ofstream log_file;
	
	
	//! Random Generators
	gsl_rng **randGen_of_species;
	gsl_rng *randGen_general_asynchronous;
	gsl_rng *randGen_general_synchronized;
	
	//! refinement variable
	INT64 num_rem_reject_child_is_refined, num_rem_reject_initial_refined, num_rem_reject_refined_neighbour;
	
	
	//! Field related gloabals
	INT32  *COMPs_FType;
	INT32  *VisuNr_FType;
	char **Field_Name;
	INT32   Field_Name_Size;


	bool collect_only_ion_species;
	
	//! Particle Collecting associated variables
	INT64  num_total_particles;
	INT64  num_injected_particles;
	INT64  num_total_particles_collected,  num_particles_collected_from_parent;
	INT64 *num_total_particles_in_L;
	INT64 *num_split_in_L, *num_merge_in_L;
	INT64  num_merge_canceled_out_of_cell, num_split_canceled_low_weight, num_split_canceled_out_of_cell;
	bool *is_inflow_species;

	INT64* num_Blocks_transfered_L;
	INT64* num_particle_processed_L;
	
	
	//! Particle Movement associated variables
	ofstream ObstacleParticle_FILE, AnyParticle_FILE;
	INT64 num_particle_total_storage, num_pArrayReallocs;
	INT64 num_unhooked, num_reset, num_moved, num_accelerated;
	INT64 num_deleted_too_fast_particle;
	INT64 max_PiC;
	ofstream TraceParticle_FILE;
	
	INT32* level_time, *step_time, *CYCLE_of_L;
	INT64* move_time;
	

	INT64 num_marked_particle;

	INT64 local_mark_particle_counter;
	
	WEIGHT_REAL weight_limit;
// 	PARTICLE_REAL *thermal_velocity_of;
	D_REAL *thermal_velocity_perp_of;
	D_REAL *thermal_velocity_para_of;
        PARTICLE_REAL fastest_particle_v2[MAX_LEVEL_DEF+1];
	
	INT64 *num_out_of_Box, *num_in_obst;
	
	//! charge exchange related
	D_REAL velocity_rate[30][NUM_PARTICLE_SPECIES];
		
        bool species_does_chemical_reactions[NUM_PARTICLE_SPECIES];
	
	INT64 num_newlyionised_particles;
	INT32 insert_every_x_TL = 10;
	
	D_REAL obstacle_MP_weight_each_t0[NUM_PARTICLE_SPECIES];
	D_REAL Neutral_vth[NUM_NEUTRAL_SPECIES]; 		
	D_REAL Neutral_Betas[NUM_NEUTRAL_SPECIES];
	
        //! Chemical Reaction Related
        D_REAL max_ChemicalReaction_Probability;
        D_REAL max_ChemicalReaction_Rates[NUM_PARTICLE_SPECIES*(NUM_PARTICLE_SPECIES+1)];
        
	//! Negative Particles related
	bool negative_particles_in_Simu = false;
	
	D_REAL norm_IonProduction_Rate_fromNeutSpec[NUM_PARTICLE_SPECIES];

	
	//! geometry variables
	D_REAL** Blk_Length_of;
	CBlock* Root_Block_Array;
	F_REAL *global_min_refValue, *global_max_refValue;
	INT64 num_refined_Octs, num_removed_gatherBlks;
	INT64  *total_Blocks_L;
	INT64 total_num_mpi_tags;

	INT32 SFC_RB_power;
	
	//! extern profiles
	D_REAL **extern_1D_Field;
	D_REAL  **values_EBFs;
	INT32 *num_rows_EBFs;
	
	
	//! for Cross Section
	INT32  blk_indedx_L0;
	INT32* cell_indedx_L;
	
	WEIGHT_REAL startWeight_of_smallestCell;
	

	D_REAL total_collected_density, total_particle_energy, total_particle_mass;
	D_REAL total_magnetic_field[3], total_collected_temperature[3], total_collected_momentum[3], total_particle_momentum[3];

	D_REAL initial_magnetic_energy, initial_particle_energy, initial_particle_number;
	D_REAL initial_particle_momentum[3];
	WEIGHT_REAL initial_particle_mass;


	INT64 num_split_since_start, num_merge_since_start;
	
	bool DIR_of_particleTracks2D_exist;
	bool DIR_of_particleTracks3D_exist;


	bool DIR_of_trajectories_exist;
	bool DIR_of_lines_exist;
	bool DIR_of_FieldID_exists[NUM_FIELDS];
	
	//!-------------------------------------------------------------//
	//!--- misc Variables -> set to extern in parameters.h ---------//
	//!-------------------------------------------------------------//
	CBlock**  BlockList_of_Lev;
	CBlock**  GATHER_BlockList_of_Lev;
	INT32 total_active_Blocks = 0;
	
	INT32  TL = 0;
	INT32 num_scalar_average_fields;
	//! in order to build an average time of time level, time since prog start
	//! has to be devided by (TL-TL_at_last_restore_state)
	INT32  TL_at_last_restore_state;

	FILE_REAL* XBlock_CrossSection;
	FILE_REAL* YBlock_CrossSection;
	FILE_REAL* ZBlock_CrossSection;
	
	D_REAL **BSOR_opti_array;
	D_REAL **DC_omega_opti_array;

