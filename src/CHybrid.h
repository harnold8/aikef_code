

#ifndef CHYBRID_H
#define CHYBRID_H

#include "defines.h"
#include "CBlk.h"
#include "parameters.h"


#include <iostream>
#include <fstream>
#include <time.h>
#include <assert.h>



using namespace std;

//!-------------------------------------------------------------//
//!------------- some general Comments -------------------------//
//!-------------------------------------------------------------//
//! All Hybrid Methods are grouped into one class.


class CHybrid
{

public:
//!*************************************************************//
//!************** member Variables *****************************//
//!*************************************************************//


//!-------------------------------------------------------------//
//!-------------- File Communication Variables -----------------//
//!-------------------------------------------------------------//


double global_MPI_lastState_time;



double  run_start_mpi,  loop_start_mpi;
double  variation_time_mpi, variation_time_comp;
clock_t run_start_comp, loop_start_comp;
time_t  run_start_phys, loop_start_phys;
time_t last_save_state_time;

bool particle_restored;

ofstream Outfile2D;



// INT32* GN_type_after;
INT32 GN_type_after[NUM_SUB_CYCLE+3];
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
// INT32* GN_type_Pe_after;
INT32 GN_type_Pe_after[NUM_SUB_CYCLE+3];
#endif

double run_time_mpi;


//!*************************************************************//
//!************** member functions *****************************//
//!*************************************************************//


//!-------------------------------------------------------------//
//!-------------- MPI related Functions ------------------------//
//!-------------------------------------------------------------//

void init_MPI(int argc, char *argv[]);

void finalize_MPI(void);
void synchronize_allProcesses_MPI(void);

void show_information(INT64* local_info_values,
			stringstream* info_names,
			const INT32 num_values,
			INT32 option);

void check_for_NaN_MPI(bool is_NaN_this_process);
bool estimate_new_block_process_dependency(INT32 global_num_blocks,
				     INT32* num_blocks_at_process,
				     F_REAL* global_block_Values,
				     INT32* active_responsible_mpi_process,
				     INT32* new_responsible_mpi_process);


//!-------------------------------------------------------------//
//!-------------- allocation & initialization Functions --------//
//!-------------------------------------------------------------//
void init(void);
void write_Info_Structure(void);

void init_RootBlocks(void);
void init_RootBlocks_set_mpi_process(INT32* root_mpi_process);

void alloc_RootBlks(void);
void alloc_variables(void);
void read_EBF(void);
void static_refinement_cuboid(void);
void static_refinement_sphere(void);


void static_refinement_ZylinderX(void);
void static_refinement_ZylinderY(void);
void static_refinement_ZylinderZ(void);

void do_consistency_check_init_derived_variables(void);



void create_gather_Blocks(void);
void remove_gather_Blks(void);
void set_Blk_optimal_MPiC(void);

void init_Eta_Profile(void);


#if defined nonadiabatic_gradPE_TERM
void init_PE_Profile(void);
#endif
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
void set_analytical_neutral_profile(INT32 neutral_species); 
void init_Neutral_Profile(void);
#endif
#if defined(use_dust_species_as_field)
void init_Dust_Profile(void);
#endif

//!-------------------------------------------------------------//
//!-------------- Particle related Functions -------------------//
//!-------------------------------------------------------------//
void fill_Particle_in_Box(void);


void move_Particle(void);
void mark_particle(void);
void trace_particle_output(void);
void reassign_particle(D_REAL &reassign_time);



void reassign_particle_child_to_parent_LB(void);
void reassign_particle_parent_to_child_LB(void);


void reassign_particle_iDirection(void);
void reassign_particle_jDirection(void);
void reassign_particle_kDirection(void);

void reassign_particle_equal_level(INT32 _m1NB_, INT32 _p1NB_,
				   INT32* SENTm1_start_ijk,
				   INT32* SENTm1_end_ijk,
				   INT32* SENTp1_start_ijk,
				   INT32* SENTp1_end_ijk,
				   INT32* RECVm1_start_ijk,
				   INT32* RECVm1_end_ijk,
				   INT32* RECVp1_start_ijk,
				   INT32* RECVp1_end_ijk,
				   INT32  entries_in_MPiC_Info_package);



void accelerate_Particle(void);
void accelerate_Particle_of_L(INT32 level, bool do_show_info);


void reset_particle_time_of_L(INT32 level);
void move_particle_of_L(INT32 level, bool do_show_info);
void split_particle_of_L(INT32 level, bool do_show_info);
void merge_particle_of_L(INT32 level, bool do_show_info);

void inject_particle_to_track(INT32 species);

void charge_exchange(void);
void chemical_Reactions(void);
void prepare_Recombination_Density();
void calc_Recombination_Density();

void read_extern_IonProfile(INT32 rho_extern, char* filename);
void read_extern_IonProfile_uniform_grid(INT32 neutralSpec, INT32 id_densityfield, INT32 id_velocityfield);

void send_injected_Ions_to_Blks(INT32 species,
				PARTICLE_REAL start_weight_HI,
				PARTICLE_REAL** positions,
				PARTICLE_REAL** velocities,
				const INT32* num_particles_to_inject);
void send_injected_Ions_to_Blks_calc_weight(INT32            species,
						    PARTICLE_REAL** positions,
						    PARTICLE_REAL** velocities);
void inject_externRhoVelocity_ions(void);

void negative_particles(void);

void Split_Particle(void);
void Merge_Particle(void);
void check_weight_sorting(void);
void resize_pArrays(void);


void collect_rho_Ji(INT32 species_to_collect,
	            INT32 id_rho_type,
	            INT32 id_Ji_type,
	            bool skip_unmarked_particle);

void collect_Ui_vth2_of_species(INT32 species,
				INT32 Ji_Type,
				INT32 RhoV2mean_Type,
				bool collect_vth2);


void collect_Ui_minus(void);
void collect_RHOnp1_UIplus_LAM_GAM(void);



void average_Ui_rho_setREZrho(void);
void average_Ui_rho_setREZrho_of_L(INT32 level, bool do_show_info);

//!-------------------------------------------------------------//
//!-------------- Field Calculation Functions ------------------//
//!-------------------------------------------------------------//
// typedef void(CHybrid::*ptr2void)(void);

void Local_TL_Cycle_LF(void);
void full_hybrid_cycle_of_L(INT32 level, INT32 *moved_L);
#ifndef nonadiabatic_gradPE_TERM
void RK_advance_BField(void);

bool RK_Step(INT32 level, INT32 field_type, bool calc_B_dev, D_REAL c1, D_REAL c2);
#endif /* nonadiabatic_gradPE_TERM */

void advanceB_Obstacle(void);
void advanceB_LF_globalDT(void);
#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
bool explicit_midpoint_method_for_Pe(INT32 level);
#endif
#ifndef nonadiabatic_gradPE_TERM
void advanceB_RK_globalDT(void);
#endif /* nonadiabatic_gradPE_TERM */


void advanceB_Staggered(void);
void advanceB_unStaggered(void);

void clean_divergence(void);
void solve_Obstacle_BField(void);


bool LF_at_cycle(INT32 level, INT32 CYCLE);
void smooth_field_to_parent(INT32 level, INT32 field_type);


void FIMesh_to_HIMesh(INT32 level, INT32 id_dest, INT32 id_src);
void HIMesh_to_FIMesh(INT32 level, INT32 id_dest, INT32 id_src);

void apply_new_BField_Boundaries(void);


#if defined nonadiabatic_gradPE_TERM
void advanceB_Plasma(void){assert(advance_B_Algorithm==0); advanceB_LF_globalDT();}
#else
void advanceB_Plasma(void){if(advance_B_Algorithm)advanceB_RK_globalDT();
					     else advanceB_LF_globalDT();}
#endif /* nonadiabatic_gradPE_TERM */

void mpi_build_sum(D_REAL *local_values,
		   D_REAL *global_values,
		   INT32 num_values);

void reset_block_timing(void);
void sum_up_block_timing(void);

void count_particle_each_block(void);

void distribute_values_all_procs(CBlock** BlkList,
				 INT32 global_num_blocks,
				 INT32* num_blocks_at_process,
				 F_REAL* global_block_Values,
				 INT32 value_type);

void calc_first_E(void);
void CAM(void);
void calc_second_E(void);

void add_force(INT32 field_id);


void calc_first_E_of_L(INT32 level, bool do_show_info);
void CAM_of_L(INT32 level, bool do_show_info);
void calc_second_E_of_L(INT32 level, bool do_show_info);

void calc_scalar_fields_sum(INT32 id_Result, INT32 id_Field1, INT32 id_Field2);
void calc_vector_fields_sum(INT32 id_Result, INT32 id_Field1, INT32 id_Field2);
void multiply_field(INT32 id_Field1, D_REAL factor);

void calc_grad(INT32 in_Field_type, INT32 grad_Field_type);
void calc_rot(INT32 BField_type, INT32 rot_Field_type);
void calc_div(INT32 Field_type, INT32 div_Field_type);
void calc_ExB_Drift(INT32 id_Result, INT32 id_Field1, INT32 id_Field2);
void calc_local_gyro_radius(INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField);
void calc_local_electron_gyro_radius(INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField);
void calc_electron_temperature(INT32 Field_type);
void calc_RecombinationRate(void);

//!-------------------------------------------------------------//
//!-------------- Utility Functions ---------------------------//
//!-------------------------------------------------------------//
void measure_time(void);
void statistic(void);


void estimate_refinement_efficiency(INT32 src_type);
void set_global_MinMax_of_refValues(void);
void set_average_ref_value(INT32 id_criteria);

void smooth_Field(INT32 type, const D_REAL smooth_value[]);
void smooth_Field_of_L(INT32 level, INT32 type, D_REAL smooth_value);
void copy_Field(INT32 dest, INT32 src);
void copy_Field_L(INT32 level, INT32 dest, INT32 src);

void add_Vector2Field(INT32 dest, INT32 src, D_REAL* vec);
void add_multipliedField(INT32 dest, INT32 src, D_REAL factor);
void sub_squaredField_Multiply(INT32 dest, INT32 src, D_REAL factor);
void Multiply_fields(INT32 dest, INT32 src1, INT32 src2);
void Multiply_field(INT32 dest, INT32 src1, D_REAL src2);
void Add_fields(INT32 dest, INT32 src1, INT32 src2);
void Devide_fields(INT32 dest, INT32 src1, INT32 src2, D_REAL min);
void calc_RecombinationAlphaField(void);

void set_zero_field(INT32 dest);
void set_zero_field_inCore(INT32 dest);
void set_zero_field_incGather(INT32 dest);

void square_Field(INT32 dest, INT32 src);


//!-------------------------------------------------------------//
//!-------------- Block Communication Functions ----------------//
//!-------------------------------------------------------------//


//!-------------- copy GN related -----------------------------------
void FULL_GN_UPDATE(INT32 field);
void FULL_MOMENT_UPDATE(INT32 field);

void cp_GN_send_MPI(     INT32 level, INT32 _m1_, INT32 _p1_, INT32 field_type);
void cp_GN_recv_MPI(     INT32 level, INT32 _m1_, INT32 _p1_, INT32 field_type);
void cp_GN_equal_process(INT32 level, INT32 _m1_, INT32 _p1_, INT32 field_type);

void wait_all_MPI(INT32 level);
void wait_SendFields_MPI(INT32 level,  INT32 _m1_, INT32 _p1_);

void copy_children_field_to_parent_MPI(INT32 level, INT32 id_field);
void zero_parent_add_children_field_MPI(INT32 level, INT32 id_field);
void field_to_children_MPI(INT32 level, INT32 id_field);


void parent_send_recv_field_MPI(INT32 level, INT32 id_field,  bool do_send, bool compose);
void child_send_recv_field_MPI(INT32 level, INT32 id_field,  bool do_send);
void wait_child_send_recv_field_MPI(INT32 level, INT32 id_field);
void wait_parent_send_recv_field_MPI(INT32 level, INT32 id_field);

void manage_parent_recv_memory(bool allocate);

void update_GN_of_level(INT32 level, INT32 id_field, bool incl_first_phys_face, bool update_parent_buffer);
void update_GATHER_GC_of_level(INT32 level,INT32 field_type);

void get_GN_from_parent(INT32 level,INT32 _m1_, INT32 _p1_, INT32 field_type, INT32* end_ijk, bool incl_first_phys_face);

void moments_to_parent_smooth(INT32 field_type);


//!--------------- add GN related -----------------------------------
void Add_Boundary_Moments_of_L(INT32 level, INT32 field_type);
void add_GN_send_MPI(     INT32 level, INT32 _p1_, INT32 field_type);
void add_GN_recv_MPI(     INT32 level, INT32 _m1_, INT32 field_type);
void add_GN_equal_process(INT32 level, INT32 _m1_, INT32 field_type);





//!--------------- refinement related -----------------------------------

void pre_mesh_refinement(void);
void post_mesh_refinement(void);

void add_Blocks(void);
void remove_Blocks(void);
void set_block_tags(void);
void refine_Mesh(INT32 TL);

void redistribute_blocks(INT32 TL);
void redistribute_blocks_RB_based(void);
void redistribute_blocks_RB_independent(void);


void send_child_blocks_to_parent_process(INT32 level);


void count_Blks_create_List(INT32 level,
			    CBlock** &BlkList,
			    INT32  &num_global_blocks,
			    INT32* &num_blocks_at_process);

void exchange_blocks_MPI(CBlock** BlkList,
			 INT32 num_global_blocks,
			 INT32* active_responsible_mpi_process,
			 INT32* new_responsible_mpi_process);

void flag_octs(INT32 criteria);
void reset_flags(void);
void flag_environment(void);
void flag_full_block_Array(void);
void update_Neighbours_of_all_Blocks(void);



//!-------------------------------------------------------------//
//!-------- File Communication/Visualization Functions ---------//
//!-------------------------------------------------------------//


void show_PE(INT32 species);
void show_MPiC(INT32 field_id);
void calc_macro_Force(INT32 MacroForce_id, INT32 EField_id, INT32 UField_id,  INT32 BField_id);


void show_time(INT32 opt, INT32 field_id);
void show_responsible_Proc(INT32 field_id);

void handle_average_fields(INT32 TL);
void reset_average_fields(void);
void add_fields_to_average_fields(void);
void normalize_average_fields(void);


void show_Block_Level(INT32 field_id);
void collect_Energy_Momentum_Mass(void);

void backup_state_files(void);
void save_state(void);
bool restore_state(void);
void secure_state_files();

void build_B_total(void);
void build_B_total_of_L(INT32 level);
INT32 count_BlockTree_mpi_processes(void);

void reset_responsible_mpi_process(INT32 number_of_processes_in_BlockFile);
void reorganize_memory_for_assigned_blocks(INT32 restore_start_id, INT32 restore_end_id);

bool read_BlockFile(void);
void write_BlockFile(void);

void read_Field3D(INT32 type,  INT32 restorer_mpi_process);
void write_Field3D(INT32 type, INT32 TL, bool write_state);

void write_Particle();
void read_Particle(INT32 restorer_mpi_process);

void write_Particle_of_Block(INT32 species, CBlock *active_block);
bool read_Particle_of_Block(INT32 species, CBlock *active_block);
bool read_Particle_of_Block_convert_to_TRACK_particle(INT32 species, CBlock *active_block);
bool read_Particle_of_Block_convert_to_ORDINARY_particle(INT32 species, CBlock *active_block);



void set_Field_Comps(void);
void set_Field_Names_IDs(void);


void alloc_init_random_Generators(void);

// bool calc_cells_to_cut(float norm_pos, INT32 X);
bool is_Block_in_CrossSection(INT32 X, INT32 level, CBlock* test_Block);

void write_XBlock_Cut(INT32 Level, INT32 F_Type,CBlock* Block);
void write_YBlock_Cut(INT32 Level, INT32 F_Type,CBlock* Block);
void write_ZBlock_Cut(INT32 Level, INT32 F_Type,CBlock* Block);


void write_num_of_Blks_in_XS(void);

void write_all_native_CS(INT32 type);
void write_XCrossSection(INT32 type);
void write_YCrossSection(INT32 type);
void write_ZCrossSection(INT32 type);


void output_all(INT32 TL);
void gnuplot_output(INT32 TL);

void trajectory_output();
void trajectory_count_values(INT32 id_trajectory, INT32 &num_values, INT32 &num_positions);


void parallel_particle_track_to_single_file(void);

void parallel_output_to_single_file(INT32 id_trajectory, INT32 id_field);
void read_trajectory(INT32 id_trajectory, INT32& num_positions, D_REAL**& positions, char**& time_string_array, D_REAL**& SCField);
void trace_trajectory(INT32 id_trajectory, INT32 Field_id, INT32 num_positions, D_REAL** positions, char** time_string_array, D_REAL** SCField);
void trace_trajectory_energy_spectrum(INT32 id_trajectory, INT32 Field_id, INT32 num_positions, D_REAL** positions, char** time_string_array, D_REAL** SCField);
void parallel_output_to_single_file_energy_spectrum(INT32 id_trajectory, INT32 id_field);


void native_output_2D(INT32 TL);
void native_output_3D(INT32 TL);

void silo_output_2D(void);
void silo_output_3D(void);

void silo_2D_Trajectory(INT32 id_trajectory);
void silo_3D_Trajectory(INT32 id_trajectory);

//! -- Lineout related
void lineout_output(INT32 id_line);
void trace_line(INT32 id_line, INT32 Field_id, INT32 num_positions, D_REAL** positions);
void trace_line_energy_spectrum(INT32 id_line, INT32 Field_id, INT32 num_positions, D_REAL** positions);
void parallel_output_to_single_lineout_file_energy_spectrum(INT32 id_line, INT32 id_field);
void parallel_output_to_single_file_lineout(INT32 id_line, INT32 id_field);
void get_positions_from_analytical_formula(INT32 id_line, INT32 &num_positions, D_REAL** &positions);

void silo_2D_LineOut(INT32 id_line);
void silo_3D_LineOut(INT32 id_line);
//! -- Lineout related end

    void init_global_variables();

//! Get Mesh related
void get_Mesh(void);
//! Get Mesh related end

//! -- Particle Detector related
void particle_detector_output(INT32 id_detector);
void parallel_particle_detector_to_single_file(INT32 id_detector);
//! -- Particle Detector related end

//! -- Statefile related
void resubjob(void);
//! -- Statefile related end
void set_fields(void);

void multiply_scalar_vector_field(int arg1, int arg2, int arg3);

void calc_reaction_rate_field(INT32 dest, INT32 neutral_species);
void analyze_chemistry();
void init_velocity_dependent_rates();
void init_ion_neutral_variables(void);

void insert_enceladus_ions(INT32 species);
void inject_Moon_ions(INT32 species);
void insert_Titan_Ions2(INT32 species);
void insert_Mars_Ions(void);
void inject_obstacle_ions(void);
void inject_cometary_ions(INT32 species);
void inject_cometary_ions_perBox(INT32 species);
void inject_shpere_ions(INT32 species);
void inject_box_ions(INT32 species);
void inject_planetary_ions(INT32 species);
void inject_simplified_Titan_Ions(INT32 species);
void insert_testScenario_ions(int arg1);
void insert_sphere_cylinder_ions(INT32 species);
void insert_Titan_Ions();
void inject_callisto_particle_to_track(INT32 species);
void insert_ions_from_neutral_profile();

void show_gradPE(INT32 field_id);

void get_height_integrated_Field(int src_field, int dest_field);

    
#if defined(use_ion_production_as_field)     
    void init_IonProduction_Profile(void);
    void set_ion_production_profile(INT32 neutralSpec);
    void calculate_Chapman(void);
    void calc_lineOfSight_Positions(D_REAL**& positions, INT32& num_positions, D_REAL* coor_of_node, D_REAL& dx, D_REAL* vector_to_sun);
    void calc_lineOfSight_ColumnDens(D_REAL& ColumnDens, D_REAL**& positions, INT32& num_positions, INT32 neutral_spec, D_REAL& dx);
    void set_value_to_field(D_REAL* ColumnDens, D_REAL* coor_of_node, INT32 field_id);
#endif    
  
    bool read_BlockFile_Dust(void); 
    bool restore_dust(void);
    void read_DustParticle(INT32 arg1);
    bool read_DustParticle_of_Block(INT32 arg1, CBlock* arg2);
    
#if defined(use_dust_species_as_field)
    void read_extern_Dust(INT32 dustSpec, INT32 id_dustSpecDens, INT32 id_dustSpecVel);
    void set_analytical_dust_plume_extern(INT32 id_dustSpecDens, INT32 id_dustSpecVel);
    void add_dust_to_field();
    void add_dust_current_to_field(int arg1);
#endif     
    
    
};

#endif
