


#include <iostream>
#include <fstream>

#include "defines.h"
#include "parameters.h"
#include "gsl/gsl_rng.h"



	using namespace std;
	
	
	//! MPI related globals
	extern MPI::Intracomm mpi_myComm;
	extern INT32 mpi_myRank, mpi_num_processes;
// 	extern INT32 *num_rootBlocks_at_process;
	extern ofstream log_file;
	
	
	//! Random Generators
	extern gsl_rng **randGen_of_species;
	extern gsl_rng *randGen_general_asynchronous;
	extern gsl_rng *randGen_general_synchronized;
	
	//! refinement variable
	extern INT64 num_rem_reject_child_is_refined, num_rem_reject_initial_refined, num_rem_reject_refined_neighbour;
	
	
	//! Field related gloabals
	extern INT32  *COMPs_FType;
	extern INT32  *VisuNr_FType;
	extern char **Field_Name;
	extern INT32   Field_Name_Size;

	extern	bool collect_only_ion_species;
	
	//! Particle Collecting associated variables
	extern INT64  num_total_particles;
	extern INT64  num_injected_particles;
	extern INT64  num_total_particles_collected,  num_particles_collected_from_parent;
	extern INT64 *num_total_particles_in_L;
	extern INT64 *num_split_in_L, *num_merge_in_L;
	extern INT64  num_merge_canceled_out_of_cell, num_split_canceled_low_weight, num_split_canceled_out_of_cell;
	extern bool *is_inflow_species;

	extern INT64* num_Blocks_transfered_L;
	extern INT64* num_particle_processed_L;
	
	
	//! Particle Movement associated variables
	extern ofstream ObstacleParticle_FILE, AnyParticle_FILE;
	extern INT64 num_particle_total_storage, num_pArrayReallocs;
	extern INT64 num_unhooked, num_reset, num_moved, num_accelerated;
	extern INT64 num_deleted_too_fast_particle;
	extern INT64 max_PiC;
	extern ofstream TraceParticle_FILE;
	
	extern INT32* level_time, *step_time, *CYCLE_of_L;
	extern INT64* move_time;


	extern INT64 num_marked_particle;
	
	extern INT64 local_mark_particle_counter;
	
	extern WEIGHT_REAL weight_limit;
//         extern PARTICLE_REAL *thermal_velocity_of; 
	extern D_REAL *thermal_velocity_perp_of;
	extern D_REAL *thermal_velocity_para_of;
        extern PARTICLE_REAL fastest_particle_v2[MAX_LEVEL_DEF+1];
	
	extern INT64 *num_out_of_Box, *num_in_obst;
	
	//! charge exchange related
	extern D_REAL velocity_rate[30][NUM_PARTICLE_SPECIES];
	
	extern bool species_does_chemical_reactions[NUM_PARTICLE_SPECIES];
		
	extern INT64 num_newlyionised_particles;
	extern INT32 insert_every_x_TL;		
	
	extern D_REAL obstacle_MP_weight_each_t0[NUM_PARTICLE_SPECIES];
	extern D_REAL Neutral_vth[NUM_NEUTRAL_SPECIES]; 		
	extern D_REAL Neutral_Betas[NUM_NEUTRAL_SPECIES];

        //! chemiacl reaction related
        extern D_REAL max_ChemicalReaction_Probability;
        extern D_REAL max_ChemicalReaction_Rates[NUM_PARTICLE_SPECIES*(NUM_PARTICLE_SPECIES+1)];
        
	//! Negative Particles related
	extern bool negative_particles_in_Simu;
	
	extern D_REAL norm_IonProduction_Rate_fromNeutSpec[NUM_PARTICLE_SPECIES];
	
	//! geometry variables
	extern D_REAL** Blk_Length_of;
	extern CBlock* Root_Block_Array;
	extern F_REAL *global_min_refValue, *global_max_refValue;
	extern INT64 num_refined_Octs, num_removed_gatherBlks;
	extern INT64  *total_Blocks_L;
	extern INT64 total_num_mpi_tags;

	extern INT32 SFC_RB_power;
	
	//! extern profiles
	extern D_REAL **extern_1D_Field;
	extern D_REAL **values_EBFs;
	extern INT32 *num_rows_EBFs;
	
	//! for Cross Section
	extern INT32  blk_indedx_L0;
	extern INT32* cell_indedx_L;
	
	extern WEIGHT_REAL startWeight_of_smallestCell;
	
	extern D_REAL total_collected_density, total_particle_energy, total_particle_mass;
	extern D_REAL total_magnetic_field[3], total_collected_temperature[3], total_collected_momentum[3], total_particle_momentum[3];


	extern D_REAL initial_magnetic_energy, initial_particle_energy, initial_particle_number;
	extern D_REAL initial_particle_momentum[3];
	extern WEIGHT_REAL initial_particle_mass;

	extern INT64 num_split_since_start, num_merge_since_start;


	extern bool DIR_of_particleTracks2D_exist;
	extern bool DIR_of_particleTracks3D_exist; 

	extern bool DIR_of_trajectories_exist;
	extern bool DIR_of_lines_exist;
	extern bool DIR_of_FieldID_exists[NUM_FIELDS];
	
	
	//!-------------------------------------------------------------//
	//!--- misc Variables -> set to extern in parameters.h ---------//
	//!-------------------------------------------------------------//
	extern CBlock**  BlockList_of_Lev;
	extern CBlock**  GATHER_BlockList_of_Lev;
	extern INT32 total_active_Blocks;
	
	extern INT32  TL, num_scalar_average_fields;
	//! in order to build an average time of time level, time since prog start
	//! has to be devided by (TL-TL_at_last_restore_state)
	extern INT32  TL_at_last_restore_state;

	extern FILE_REAL* XBlock_CrossSection;
	extern FILE_REAL* YBlock_CrossSection;
	extern FILE_REAL* ZBlock_CrossSection;

	extern D_REAL **BSOR_opti_array;
	extern D_REAL **DC_omega_opti_array;

