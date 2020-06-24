


#ifndef CBLOCK_H
#define CBLOCK_H

#include "defines.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <mpi.h>

using std::ofstream;

//! particle struct for ion description
struct particle
{

	//! relative cell coords = [0;1[
	PARTICLE_REAL rel_r[3];

	//! particle velocity
	PARTICLE_REAL v[3];

	//! weight indicates how many "real" ions are represented
	//! by this particular particle
	WEIGHT_REAL weight;


#ifdef TRACK_PARTICLE
	INT32 number;
	INT32 exchange;
#endif


	//!
// 	INT32 number;
// 
// 	INT32 moved;

	//! move time of the particle
// 	INT32 time;

	//!-------------!!!!!!!!!!!!!!!!!------------------------------
	//! NOTE:
	//! Either use double for WEIGHT_REAL instead of float
	//! or define "PARTICLE_REAL empty":
	//! Otherwise sizeof(particle) is 28 byte and 
	//! REALLY STRANGE/UNSYSTEMATIC ERRORS OCCUR !!!!
	//! So specify this variable to comple the 32 byte and
	//! leave it unused.
	//! The Erros only occur when "resize_pArrays" function is
	//! activated (see main loop). So it seems as std delete
	//! function fails.
	//!-------------!!!!!!!!!!!!!!!!!------------------------------
// 	PARTICLE_REAL empty;


};

//! particle struct for ion description
struct TRACK_particle
{

	//! relative cell coords = [0;1[
	PARTICLE_REAL rel_r[3];

	//! particle velocity
	PARTICLE_REAL v[3];

	//! weight indicates how many "real" ions are represented
	//! by this particular particle
	WEIGHT_REAL weight;


	INT32 number;
	INT32 exchange;

};

struct ORDINARY_particle
{

	//! relative cell coords = [0;1[
	PARTICLE_REAL rel_r[3];

	//! particle velocity
	PARTICLE_REAL v[3];

	//! weight indicates how many "real" ions are represented
	//! by this particular particle
	WEIGHT_REAL weight;

};




void init_Block_Globals(void);

class CBlock
{
	
public:

	
	CBlock *parent;
	CBlock *child_array[8];
	CBlock *gather_child_array[8];

	INT32 num_children;
	
// 	CBlock *Neighbour[6];
// 	CBlock *gather_Neighbour[6];

	CBlock  *Neighbour[12];
	CBlock **gather_Neighbour;


	//! one "next_Block_Array pointer" every eith block would do it,
	//! however the presented way is more straight
	//! forward and only negligible memory is wasted

	CBlock* prev_Blk_of_BlockList;
	CBlock* next_Blk_of_BlockList;

	CBlock* prev_Blk_of_GatherBlockList;
	CBlock* next_Blk_of_GatherBlockList;

	INT32  mpi_tag;

	//! Block_Nr can get huge really quickly in max level
	INT64 Block_Nr;
// 	INT32 associated_RB_Nr;
	//! 0: xMin, 1:xMax
	//! 2: yMin, 3:yMax
	//! 4: zMin, 5:zMax
	//! 6: is any ?
	//! 7: blk overlaps with obstalce
	bool is_box_boundary[8];
// 	bool package_received[2];

	bool is_gatherBlk;
	bool initial_refined;
	bool do_gather;
	bool do_receive_particle_from_parent;
	bool do_send_particle_to_childArray[8];

	//! flagged_by_neighbour is to distinguish
	//! whether block was flagged due to neighbour
	//! or refinement criteria
	bool child_flag_refinement[8];
	bool child_flag_removal[8];
	bool flagged_by_neighbour[8];


	//! MPI related
	INT32 responsible_mpi_process;


	//! Refinement Level of Block
	INT32 RLevel;

	//! Blk_Index in parent or root Array
	INT32 Blk_Index[3];

	//! origin meaning the position of the i,j,k={1,1,1} point
	//! of the current Block (Ghost cells are ignored)!!!
	F_REAL origin[3];

	F_REAL average_ref_value[8];


	//!------- Fields ----------------------
	//! ghost cells are included -----------
	//!-------------------------------------
	
	//! Field_Type[] is introduced in order to
	//! acces B,E,U etc. by numbers 3,4,5
	//! This allows several routines to be written for an
	//! arbitrary field.
	//! Check Block.cpp to see which number belongs to which field
	//! (Eg. Field_Type[3] = B1;)

	bool *Flag, *Core_Flag;

	D_REAL* Field_Type[NUM_FIELDS];

	//! Field Transmission Realated
	D_REAL *send_buffer[2];

	//! Particle Realated
	INT32** num_MPiC;
	INT32*  Blk_optimal_MPiC;
	particle*** pArray;
	INT32 **size_of_pArray;

	//! Particle Transmission Realated
	INT32*     MPiC_Info_package[8];
	particle* particle_package[8];



	//! Eight request each block are required since
	//! one parent need to sent information to all
	//! eight children that may be on different
	//! processes
	MPI::Request req_is_MPI_send_completed[NUM_REQUESTS];


	//! comp time that is required to process the field of this block
	F_REAL time_process_fields, time_process_fields_incChilds;

	//! comp time that is required to process the particle of this block
	F_REAL time_process_particle, time_process_particle_incChilds;
	
	
	//! count number of particle in this Block
	INT32 num_particle, num_particle_incChilds;


	//!------ Member Functions --------------------------------
	void get_a_b_c_from_abc(INT32 abc, INT32* a_b_c);
	void set_Fields(void);
	void copy_Field(INT32 dest, INT32 src);
	void calc_CFBG_Vector(D_REAL* Vector, D_REAL* position);

	void square_Field(INT32 dest, INT32 src);
	void Multiply_fields(INT32 dest, INT32 src1, INT32 src2);
        void Multiply_field(INT32 dest, INT32 src1, D_REAL src2);
        void Add_fields(INT32 dest, INT32 src1, INT32 src2);
        void Devide_fields(INT32 dest, INT32 src1, INT32 src2, D_REAL min);
	void add_Vector2Field(INT32 dest, INT32 src, D_REAL* vec);
	void add_multipliedField(INT32 dest, INT32 src, D_REAL factor);
	void sub_squaredField_Multiply(INT32 dest, INT32 src, D_REAL factor);
        void calc_RecombinationAlphaField();

	void init_blk_origin(void);
	void init_blk_boundaries(void);

	void init_Eta_Profile(void);

// added by HR
#if defined nonadiabatic_gradPE_TERM
	void init_PE_Profile(void);
#endif
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
  void which_neutral_profile(INT32 neutralSpecies, D_REAL *x, D_REAL *Nneutral, D_REAL *VNeutral);
  void set_analytical_neutral_profile(INT32 neutralSpecies);
  void barometric_neutral(INT32 neutralSpecies, D_REAL *x, D_REAL *Nneutral, D_REAL *VNeutral);
  void init_Cometary_Neutral_Profile(INT32 neutralSpecies,  D_REAL *x, D_REAL *Nneutral, D_REAL *Vneutral);
  void neutral_profile(INT32 neutralSpecies, D_REAL *x, D_REAL *Nneutral, D_REAL *Vneutral);
  void neutral_profile_Enceladus(INT32 neutralSpecies, D_REAL *x, D_REAL *Nneutral, D_REAL *Vneutral);
  void insert_ions_from_ionprod_file(INT32 insert_every_x_TL, INT32 id_oct, INT64* num_injected_particles_thisTL);
  void insert_ions_from_neutral_profile(INT32 insert_every_x_TL, INT32 id_oct, INT64* num_injected_particles_thisTL);
  void init_neutral_profile_test(INT32 neutralSpecies, D_REAL* x, D_REAL* Nneutral, D_REAL* VNeutral);
#endif


/*	void set_RhoUi_extern(INT32  rho_extern,
			     INT32  *num_Nodes,
			     D_REAL *Origin,
			     D_REAL *Length,
			     D_REAL *rho,
			     INT32 &num_values_not_in_extern_box);
			     */

	void set_RhoUi_extern(INT32  id_densityfield,
			      INT32 id_velocityfield,
			     INT32  *num_Nodes,
			     D_REAL *Origin,
			     D_REAL *Length,
			     D_REAL *rho,
			     INT32 &num_values_not_in_extern_box);

	

	
	void set_zero_Field(INT32 field_type);
	void set_zero_Field_inCore(INT32 field_type);


	void manage_parent_recv_memory(bool allocate);
	void alloc_GatherMemory(void);
	void alloc_process_specific_Memory(void);

	void delete_Memory(void);
	void delete_GatherMemory(void);
	void delete_process_specific_Memory(void);


	void init_Block(INT32 Blk_Nr, INT32 *ind, INT32 Level, CBlock *Parent);
	void smooth_Field(INT32 type, D_REAL smooth_value);


	void reset_average_fields(void);
	void add_fields_to_average_fields(void);
	void normalize_average_fields(void);

	//! ----------- General physical Functions -----------------
	void calc_E(INT32 BField_type, INT32 UField_type, INT32 rho_type);
	void show_PE(INT32 species);
	void calc_macro_Force(INT32 MacroForce_id, INT32 EField_id, INT32 UField_id,  INT32 BField_id);
	void CAM(void);

	void add_force(INT32 field_id);

	void calc_div(INT32 Field_type, INT32 div_Field_type);
	void calc_rot(INT32 Field_type, INT32 rot_Field_type);
	void calc_grad(INT32 in_Field_type, INT32 grad_Field_type);
	void calc_local_gyro_radius(INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField, D_REAL SI_x0);
        void calc_local_electron_gyro_radius(INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField, D_REAL SI_x0);


// 	void count_flagged_nodes(INT32 src_type, F_REAL *max_Value, INT32* num_flagged_nodes);
	void set_average_ref_value(INT32 id_criteria);

	void FIMesh_to_HIMesh(INT32 id_dest, INT32 id_src);
	void HIMesh_to_FIMesh(INT32 id_dest, INT32 id_src);


	//! ------------ Advance BField Functions ------------------



#if defined(nonadiabatic_gradPE_TERM)
	inline void set_pointers_to_calc_derivatives(INT32 B_FIELD_ID,
						       INT32 PE_FIELD_ID,
						       INT32 U_FIELD_ID,
						       INT32 rezRHO_FIELD_ID,
						       INT32 ETA_FIELD_ID);

	inline void set_pointers_to_calc_derivatives_CFBG_BField(INT32 B_SW_ID,
								  INT32 B_TOTAL_ID,
								  INT32 PE_FIELD_ID,
								  INT32 U_FIELD_ID,
								  INT32 rezRHO_FIELD_ID,
								  INT32 ETA_FIELD_ID);
#endif
	inline void set_pointers_to_calc_derivatives(INT32 B_FIELD_ID,
						       INT32 U_FIELD_ID,
						       INT32 rezRHO_FIELD_ID,
						       INT32 ETA_FIELD_ID);

	inline void set_pointers_to_calc_derivatives_CFBG_BField(INT32 B_SW_ID,
								  INT32 B_TOTAL_ID,
								  INT32 U_FIELD_ID,
								  INT32 rezRHO_FIELD_ID,
								  INT32 ETA_FIELD_ID);





	
	inline void calc_derivatives(void);
	inline void calc_derivatives_conservative(void);
	inline void calc_derivatives_conservative_CFBG_BField(void);
	
#if defined(nonadiabatic_gradPE_TERM) && !(defined nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	inline void calc_derivatives_for_Pe_conservative(void);
	inline void calc_derivatives_for_Pe_conservative_CFBG_BField(void);
#endif


	void (CBlock::*LF_cycle[NUM_SUB_CYCLE+2])(void);
	inline void LF_at_cycle(INT32 num_cycle){(this->*LF_cycle[num_cycle])();}

	void LF_odd_B0toB1(void);
	void LF_advance_B_even(void);
	void LF_advance_B_odd(void);
#if defined(nonadiabatic_gradPE_TERM) && !(defined nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	void explicit_midpoint_method_for_Pe_1(void);
	bool explicit_midpoint_method_for_Pe_2(void);
#endif
	void LF_even_B_Solution(void);
	bool LF_average_B_Solution(void);

	void build_B_total(void);


	void LF_B_Step(INT32 B_to_update,
		       INT32 B_derivatives,
		       D_REAL dt_step);


#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	void LF_B_Step_conservative(INT32 B_to_update,
				    INT32 PE_to_update,
				    INT32 B_derivatives,
				    INT32 PE_derivatives,
				    D_REAL dt_step);

	void LF_B_Step_conservative_CFBG_BField(INT32 B_to_update,
						INT32 PE_to_update,
						INT32 B_derivatives,
						INT32 PE_derivatives,
						D_REAL dt_step);
#else
	void LF_B_Step_conservative(INT32 B_to_update,
				    INT32 B_derivatives,
				    D_REAL dt_step);

	void LF_B_Step_conservative_CFBG_BField(INT32 B_to_update,
						INT32 B_derivatives,
						D_REAL dt_step);
#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	void LF_Pe_Step_conservative( INT32 PE_to_update,
				      INT32 PE_derivatives,
				      D_REAL dt_step);
	
	void LF_Pe_Step_conservative_CFBG_BField( INT32 PE_to_update,
						  INT32 PE_derivatives,
						  D_REAL dt_step);
#endif




	void calc_Faraday_EField(INT32 id_B_FIMesh);
	void calc_Faraday_EField_staggered(INT32 id_B_FIMesh, INT32 id_B_HIMesh);

	void dtB_Faraday_at_cycle(INT32 CYCLE);

#ifndef nonadiabatic_gradPE_TERM
        bool RK_advance_B(bool calc_B_dev, D_REAL c1, D_REAL c2);
	void RK_B_Step(INT32 KX_type, INT32 B_dev_type);
	void RK_B_Step_conservative_CFBG_BField();
#endif /* nonadiabatic_gradPE_TERM */

	void calc_RHS_CN(INT32 type_RHS, INT32 type_BField);

	D_REAL Psi_DC_SOR_Step(INT32 counter, INT32 div_type, bool Odd_Even_Step);
	void project_divB_out_of_BField(INT32 type_BField);

	D_REAL Obstacle_BField_SOR_Step(INT32 counter, INT32 B_type, INT32 type_RHS, bool Odd_Even_Step);

	void estimate_nodes_BBound_value(D_REAL *BBound, INT32* cell_indices);
	void apply_new_BField_Boundaries(INT32 B_Field_type);
	//! -------- Refinement Functions ------
	bool is_removing_permitted(INT32 id_oct);
	CBlock* has_refined_neighbour(void);
	bool is_refinement_permitted(INT32 id_oct);


	CBlock* seek_Block_at(INT32 req_level, D_REAL *position, INT32 &child_id_at_position);
	void refine_Oct_at(INT32 level, D_REAL *position, bool ref_in_gatherOct);
	void refine_gatherEnvironment(void);


	void flag_full_environment(INT32 child, INT32 *num_flagged_neighbours);
	void refine_oct_environment(INT32 id_child);
	bool remove_child_array_OLD(CBlock* &List_Root);


	void refine_Oct(INT32 id_oct);
	void refine_Oct_set_mpi_process(INT32 id_oct, INT32 mpi_process);

	void remove_Child(INT32 id_child);
	void alloc_set_OctFields(INT32 id_oct);

	void remove_Block_from_BlkList(INT32 id_child);
	void remove_GthrBlock_from_GthrBlkList(INT32 id_gChild);

	void remove_gather_Blk(INT32 id_gChild);
	void refine_gatherOct(INT32 id_oct);


// 	void update_BlockArray_List_L();
	void add_Block_to_BlkList(INT32 id_oct);
	void add_Block_to_GthrBlkList(INT32 id_oct);

	//! -------- Parent Child Neighbour Communication Functions ------
	void wait_SendFields_MPI(INT32 _m1_, INT32 _p1_);
	void wait_for_completedSendParticle_MPI(INT32 _m1_, INT32 _p1_);
	void wait_all_MPI();

	void cp_GN_send_MPI	(INT32 _m1_, INT32 _p1_, INT32 field_type);
	void cp_GN_recv_MPI	(INT32 _m1_, INT32 _p1_, INT32 field_type);
	void cp_GN_equal_process(INT32 _m1_, INT32 _p1_, INT32 field_type);



	void add_GN_send_MPI	 (INT32 _p1_, INT32 field_type);
	void add_GN_recv_MPI	 (INT32 _m1_, INT32 field_type);
	void add_GN_equal_process(INT32 _m1_, INT32 field_type);

	void compose_recvd_fields(INT32 id_field, INT32 num_distinct_mpi_procs, D_REAL** buffer_proc);

	void copy_recvd_fields(INT32 id_field, INT32* process_buffer_id, D_REAL** buffer_proc);


// 	void BV_from_parent_OddNP(INT32 parent_type, INT32 child_type);
// 	void BV_from_parent_EvenNP(INT32 parent_type, INT32 child_type);
//  	void BV_from_parent(INT32 parent_type,INT32 child_type);
// 	{(this->*ptr2BV_from_parent)(parent_type,child_type);};
//  	void (CBlock::*ptr2BV_from_parent)(INT32 parent_type,INT32 child_type);


	void cp_imin_Boundaries(INT32 field);
	void cp_imax_Boundaries(INT32 field);
	void cp_jmin_Boundaries(INT32 field);
	void cp_jmax_Boundaries(INT32 field);
	void cp_kmin_Boundaries(INT32 field);
	void cp_kmax_Boundaries(INT32 field);

 	void set_inhom_B_field(INT32 field, INT32 *start, INT32 *end, INT32 *src);


	void apply_BField_Boundaries(INT32 id_field);
#if defined nonadiabatic_gradPE_TERM
	void apply_PE_Boundaries(INT32 id_field);
#endif
	void cp_Boundaries(INT32 field, const bool *use_inflow);

	void find_im1_Neighbour(void);
	void find_ip1_Neighbour(void);
	void find_jm1_Neighbour(void);
	void find_jp1_Neighbour(void);
	void find_km1_Neighbour(void);
	void find_kp1_Neighbour(void);

	void find_gather_im1_Neighbour(void);
	void find_gather_ip1_Neighbour(void);
	void find_gather_jm1_Neighbour(void);
	void find_gather_jp1_Neighbour(void);
	void find_gather_km1_Neighbour(void);
	void find_gather_kp1_Neighbour(void);


	void set_Blk_optimal_MPiC(void);

	void find_all_ordanary_Neighbours(void);
	void find_all_gather_Neighbours(void);

	void link_all_ordanary_Neighbours(void);
	void link_all_gather_Neighbours(void);

	void unlink_all_Neighbours(void);

	void field_from_parent(INT32 parent_type, INT32 child_type, INT32* start_ijk, INT32* end_ijk);


	void add_field_to_parent(INT32 parent_type, INT32 child_type);
	void field_to_parent_smooth(INT32 parent_type, INT32 child_type);
	void moments_to_parent_smooth(INT32 field_type);

// 	void add_inject_GC_to_parent(INT32 parent_type, INT32 child_type);
#ifndef nonadiabatic_gradPE_TERM
	void RK_sequence(void);
#endif /* nonadiabatic_gradPE_TERM */

	//!-----------------------------------------------
	//! Particle related functoins
	//!-----------------------------------------------

	//! Particle Initialization
	void fill_Block(void);
	void fill_Oct(INT32 id_oct);
	void fill_empty_Cell(INT32 species, INT32 i_j_k, PARTICLE_REAL *v_init, PARTICLE_REAL rho_init);
	void fill_empty_Volume(INT32 species, const INT32* start, const INT32* end);

	void copy_empty_Cell(INT32 species, INT32 dest_i_j_k, INT32 src_i_j_k);
	void copy_empty_Volume(INT32 species, const INT32* start, const INT32* end, const INT32* offset);

	void inject_extern_profile_ions(INT32 species, INT32 insert_every_x_TL, INT32 id_oct);
        void inject_cometary_profile_ions(INT32 oct, INT32 species, INT32 insert_every_x_TL);
        D_REAL getNeutralDensityAtPoint(PARTICLE_REAL *r_intern,INT32 * cell_indices,INT32 ion_species);
        void getNeutralVelocityAtPoint(PARTICLE_REAL *r_intern,INT32 * cell_indices, INT32 ion_species,PARTICLE_REAL *vec_v);
	D_REAL getIonProdRateAtPoint(PARTICLE_REAL* r_intern, INT32* cell_indices, INT32 ion_species, INT32 neutral_species);
        void getIonProdVelocityAtPoint(PARTICLE_REAL* r_intern, INT32* cell_indices, INT32 ion_species, INT32 neutral_species, PARTICLE_REAL* vec_v);
	void negative_particles(INT64* num_negCells_delPart, INT32 level);

	void get_velocity(D_REAL**& velocity, INT32 species, INT32 i_j_k, D_REAL*& gewicht);
	
	void charge_exchange(INT64* num_exchanged_particle);
        inline bool calculate_charge_exchange(PARTICLE_REAL* x, PARTICLE_REAL* v, D_REAL dt_particle, D_REAL reaction_rate, D_REAL Neutral_dens, D_REAL* Neutral_vel, D_REAL RecombinationRate, D_REAL IonDensity, D_REAL species);
//         inline void calc_Reaction_CrossSection(PARTICLE_REAL *v, D_REAL *Reaction_CrossSection, INT32 species, D_REAL *total_crossSec, D_REAL RecombinationAlpha);
//         inline void calc_Reaction_CrossSection(PARTICLE_REAL *v, D_REAL *Reaction_CrossSection, INT32 species, D_REAL *total_crossSec);
        inline INT32 determine_dest_species(D_REAL *Reaction_CrossSection, D_REAL total_crossSec);
	inline INT32 determineNeutralSpecies(INT32 destinationSpecies, D_REAL *processRates, D_REAL *processRates_eachNeutralSpec);
	void estimate_extern_value(D_REAL &v, INT32 cord_id, INT32 field_id, PARTICLE_REAL *x);

        void chemical_Reactions(INT64* num_exchanged_particle);
	D_REAL pickScalarFieldValue(INT32 field_id, PARTICLE_REAL* vec_r, INT32* ind);
        void pickVectorFieldValue(INT32 field_id, PARTICLE_REAL* vec_r, INT32* ind, PARTICLE_REAL *vec_v);
        INT32 determineDestinationSpecies(D_REAL totalProcessRate, D_REAL* processRates);
        bool checkChemicalReaction(D_REAL totalProcessRate, D_REAL dt_particle);
        void setVelocityAfterReaction(PARTICLE_REAL *x,PARTICLE_REAL* v);
        void prepare_Recombination_Density();
        void calc_Recombination_Density();
        
	void cell_centre_normedCoords(double* x, INT32* cell_indices);
	void normed2internCoords(PARTICLE_REAL* intern_r, PARTICLE_REAL* x);
	void intern2normedCoords(double* x, double* cell_intern_r, INT32* cell_indices);
	void intern2normedCoords_HIMesh(double* x, double* cell_intern_r, INT32* cell_indices);
	

	//! Particle Moment functions


	void average_Ui_rho_setREZrho(void);



	void calc_Energy_Momentum_Mass(void);

	void collect_Ji(INT32 species_to_collect,
			INT32 total_Ji_type,
			INT32 species_Ji_type,
			INT32 species_vth2_type,
			bool collect_vth2);


	void collect_rho_Ji_LAM_GAM(INT32 Ji_type,
				    INT32 Gam_type,
				    INT32 species_Ji_type);

	void collect_rho_Ji(INT32 species_to_collect,
			    INT32 species_rho_type,
			    INT32 species_Ji_type,
			    bool skip_unmarked_particle);


	void sum_up_rhoSpecies(INT32 rho_type,
                           	INT32 Lam_type);

	void convert_Ji2Ui(INT32 species,
			   INT32 rho_type,
			   INT32 Ji_type,
			   INT32 Ui_type/*,
			   bool set_inner_Dens*/);

	void convert_vth2rho_to_vth2(INT32 species,
				     INT32 rho_type,
				     INT32 vth2rho_type);

	void gather_Oct(short species,
			  bool collect_rho,
			  bool collect_Ji,
			  bool collect_vth2,
			  bool skip_unmarked_particle,
			  INT32 id_oct,
			  D_REAL* rho,
			  D_REAL* JiX);



	void gather_from_parent_Block(short species,
			 	  bool collect_rho,
			  	  bool collect_Ji,
			  	  bool collect_vth2,
			  	  D_REAL* rho,
			  	  D_REAL* JiX);





	//! Particle Refinement functions
	void estimate_Cords_of_split_Particle_CoM_ON(PARTICLE_REAL* x_new1,
						     PARTICLE_REAL* x_new2,
						     PARTICLE_REAL* x_old,
						     PARTICLE_REAL* v_part);

	void estimate_Cords_of_split_Particle_CoM_OFF(PARTICLE_REAL* x_new1,
						      PARTICLE_REAL* x_new2,
						      PARTICLE_REAL* x_old,
						      PARTICLE_REAL* v_part);

	void set_Cords_perp_to_velocity_CoM_OFF(PARTICLE_REAL* x_part,
						       PARTICLE_REAL* v_part);

	void split_most_centred_particle(void);
	void split_heaviest_particle(void);
	void split_1to6particle(void);
	void check_weight_sorting(void);
	void merge_particle(void);

	//! Particle Movement functions

	void trace_particle(INT64* particle_counter_of_species);
	void mark_particle(INT64* particle_counter_of_species, INT64 num_partcle_in_list, INT64* Partilce_number_List);
	void write_particle_number(void);


	void delete_particle(INT32 species, INT32 i_j_k, INT32 &part_index);
	void move_particle(INT32 id_oct);
	void move_particle_noObsBoxBoundCheck(INT32 id_oct);
	bool block_intern_move_particle(void);
	bool cell_reassign_particle(void);
	void inject_boundary_particle(INT32 id_oct);
	void apply_inhom_particle_bounds(const INT32 comp, INT32 species, const INT32 *start, const INT32 *end, const INT32 *src);
	void reset_particle_time_of_L(void);


	void count_particle(void);

	void resize_pArray(INT32 species,
			   INT32 node,
			   INT32 new_size);

	void resize_all_pArrays_of_Block(void);

	void add_particle_to_pArray(INT32 species,
				    INT32 new_a_b_c,
				    particle* particle_to_insert);

	void swap_particles_pArray(INT32     species,
				   INT32     part_index,
				   INT32     a_b_c,
				   INT32 new_a_b_c,
				   particle* &active_particle);

	void accelerate_particle(INT32 id_oct);





	//!--------------------------------------
	//! pArray based Block Communication
	//!--------------------------------------

	void delete_package_memory(void);


	//! particle transfer for entire octs
	void prepare_pPackage_toParent(void);
	void prepare_pPackage_toChild(INT32 id_child);

	void insert_pPackage_fromParent(void);
	void insert_pPackage_fromChild(INT32 id_child);


	void send_particle_toChild(INT32 id_child);
	void recv_particle_fromParent(void);

	//! particle transfer at Level Boundaries (LB)
	void send_particle_toParent_LB(void);
	void send_particle_toChildren_LB(void);
	void wait_send_particle_toParent_LB(void);
	void wait_send_particle_toChildren_LB(void);

	void recv_particle_fromChildren_LB(void);
	void recv_particle_fromParent_LB(void);


	void send_pPackage_toParent_MPI(INT32 id_direc,
					INT32 entries_in_MPiC_Info_package);

	void recv_pPackage_fromChild_MPI(INT32 id_child,
					 INT32 id_direc,
					 INT32 entries_in_MPiC_Info_package);

	void send_pPackage_toChild_MPI(INT32 id_child, INT32 entries_in_MPiC_Info_package);
	void recv_pPackage_fromParent_MPI(INT32 entries_in_MPiC_Info_package);



	void check_if_receive_particle_from_parent_is_required(void);






	//! 1)
	void prepare_pPackage_toNeighbour(INT32* m1_start_ijk,
				      INT32* m1_end_ijk,
				      INT32* p1_start_ijk,
				      INT32* p1_end_ijk,
				      INT32 _m1_,
				      INT32 _p1_);

	//! 2)
	void recv_MPiC_Info_fromNeighbour(INT32 entries_in_MPiC_Info_package,
				        INT32 _m1_,
				        INT32 _p1_);


	//! 3)
	void send_MPiC_Info_toNeighbour(INT32 entries_in_MPiC_Info_package,
			    INT32 _m1_,
			    INT32 _p1_);

	//! 4)
	void send_particle_toNeighbour(INT32 entries_in_MPiC_Info_package,
			   INT32 _m1_,
			   INT32 _p1_);



	//! 5)
	void insert_pPackage_fromNeighbour(INT32* m1_start_ijk,
					   INT32* m1_end_ijk,
					   INT32* p1_start_ijk,
					   INT32* p1_end_ijk,
					   INT32 _m1_,
					   INT32 _p1_);

	//! 6)
	void recv_particle_fromNeighbour(INT32 entries_in_MPiC_Info_package,
				   INT32 _m1_,
				   INT32 _p1_);


	//! 7)
	void inserting_received_MPiC_Part(INT32* m1_start_ijk,
					  INT32* m1_end_ijk,
					  INT32* p1_start_ijk,
					  INT32* p1_end_ijk,
					  INT32 _m1_,
					  INT32 _p1_);





	void send_MPiC_Package_MPI(INT32  id_request,
				   INT32  id_package,
				   INT32  dest_process,
				   INT32  num_entries_in_InfoPackage,
				   INT32* &MPiC_Info_package);


	void send_partPackage_MPI(INT32  id_request,
				  INT32  id_package,
				  INT32  dest_process,
				  INT32  num_entries_in_InfoPackage,
				  particle* &particle_package,
				  INT32* &MPiC_Info_package);

	void receive_infoPackage_MPI(INT32  id_package,
				     INT32  src_process,
				     INT32  num_entries_in_InfoPackage,
				     INT32* &MPiC_Info_package);

	void receive_partPackage_MPI(INT32  id_package,
				     INT32  src_process,
				     INT32  num_entries_in_InfoPackage,
				     particle* &particle_package,
				     INT32* &MPiC_Info_package);




	void prepare_pPackage(INT32* start_ind,
			      INT32* end_ind,
			      INT32  pID);

	void insert_add_pPackage(INT32* start_ijk,
				 INT32* end_ijk,
				 particle* &particle_package,
				 INT32* &MPiC_Info_package);

	void insert_pPackage_emptyBlock(INT32* start_ijk,
					INT32* end_ijk,
					particle* &particle_package,
					INT32* &MPiC_Info_package);

	void insert_child_pPackage(INT32* start_ijk,
				   INT32* end_ijk,
				   INT32 child,
				   INT32 id_package);



	void insert_parent_pPackage(INT32* start_ijk,
				    INT32* end_ijk,
				    INT32 id_package);


	void insert_parent_pPackage_into_childArray(INT32* start_ijk,
						    INT32* end_ijk,
						    INT32 id_package);



//! -- Particle Detector related
void particle_detector_output(INT32 id_detector, std::ostream &particle_detector_file);
//! -- Particle Detector related end

//! -- GetMesh related
void getMesh(INT32 id_oct, ofstream& outfile);
//! -- GetMesh related end
   

    
    void show_gradPE(INT32 field_id);
    
    void analyze_chemistry(D_REAL*& density, INT32& nodes);

    
   void set_Hint_to_field(INT32 src_field_id, INT32 dest_field_id, short *num_Nodes, D_REAL *Origin,
			     D_REAL *Length,
			     D_REAL *fieldarray,
			     bool *integration_axis);

    
    
#if defined(use_ion_production_as_field)    
	void set_ion_production_profile(INT32 neutral_species);

	void set_ion_production_profile_Enceladus(INT32 neutral_species);
#endif

	
#if defined(use_dust_species_as_field)
	
	void set_analytical_dust_field(INT32 id_dustSpecDens, INT32 id_dustSpecVel);
	
	void add_dust_to_field(INT32 id_rho);
	void add_dust_velocity_to_field(INT32 Ioncurrent_id);
	
	void set_Dust_extern(INT32  id_densityfield,
			INT32 id_velocityfield,
			short  *num_Nodes,
			FILE_REAL *Origin,
			FILE_REAL *Length,
			FILE_REAL *rho,
			INT32 &num_values_not_in_extern_box);
#endif    
    
    bool insert_at_rnormed(D_REAL *r_normed);
	

};

#endif
