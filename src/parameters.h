


#include "defines.h"
#include "CBlk.h"

//! below has to be included on SuSE 11
#include <string.h>
#include <stdlib.h>


//ionization switch
extern const INT32 TLIO;
extern const bool ionizationonly;


//!-------------------------------------------------------------//
//!------------- 0) Normalisation Parameter --------------------//
//!-------------------------------------------------------------// 
extern const D_REAL pi;
extern const D_REAL mu_0;
extern const D_REAL m_p;
extern const D_REAL e;
extern const D_REAL kB;
extern const D_REAL mass_electron_norm;

extern const D_REAL COM_v0;
extern const D_REAL SI_x0;
extern const D_REAL SI_B0;
extern const D_REAL SI_n0;
extern const D_REAL SI_m0;
extern const D_REAL SI_v0;
extern const D_REAL SI_t0;
extern const D_REAL COM_Q[];
extern const D_REAL COM_Q_Total;

extern const D_REAL SSL;
extern const D_REAL SLT;

//! Angle of Magnetic Field to Solar Wind Velocity in deg
extern const D_REAL Bx;
extern const D_REAL By;
extern const D_REAL Bz;
extern const D_REAL B_angle;
//! Radius of Cometary Ion Injection
extern const D_REAL COM_a;
//! Characteristic Frequency of Photoionisation
//! Neutral Drag 
// extern const D_REAL ND_k;
//! Raw Neutral Density
extern const D_REAL COM_RawDensity;
extern const D_REAL calcBeta;
extern const D_REAL SW_v;
extern const D_REAL MA;
extern const D_REAL SI_M0;

extern const short flyby;

extern const D_REAL SI_v0;

extern const D_REAL SI_n0;
// extern const D_REAL nu_photo;

extern const D_REAL Te;
extern const D_REAL Te_ionosphere; 


extern const D_REAL base_density;


extern const bool analytical_dust_plume; //! set nano dust according to formula
extern const bool analytical_dust_plume_only; //! neutral dens only by formula

//! use realistic Enceladus plume model or simple test geometry
extern const bool use_test_scenario;




//! disable photoionization in Ence's and Saturn's shadow

extern const D_REAL theta0;	
extern const D_REAL phi0;


extern const D_REAL xi_nD_0;


extern const D_REAL H_d;
extern const D_REAL H_d_dust;
extern const D_REAL opening_angle_gas;
extern const D_REAL opening_angle_dust;


//!-------------------------------------------------------------//
//!------------- 1) Geometry Parameter -------------------------//
//!-------------------------------------------------------------//



//! Number of Nodes in each Block:
//! - only even values allowed (otw. false interpolation)
//! - Ghost Nodes are included
//!   -> only BlkNds_-2 Nodes are used for physical Field !
extern const INT32 BlkNds_X;
extern const INT32 BlkNds_Y;
extern const INT32 BlkNds_Z;


//! Number of Root Blocks in Szenario
extern const INT32 RB_X;
extern const INT32 RB_Y;
extern const INT32 RB_Z;

//! Radius of Obstacle
extern const D_REAL R_Moon; 
extern const D_REAL R_Obstacle;

//! Size of simulation box 
extern const D_REAL LX;
extern const D_REAL LY;
extern const D_REAL LZ;


//! Origin of simulaionbox (= Position of Obstacle)
extern const D_REAL Box_Origin[3];


//! use periodic Boundaries
extern const bool use_periodic_bounds[3];


//! specify number of cells respective boundary set are
//! set to eta value defined below
extern const INT32 resistive_bound_cells[];

//! specify value of eta at boundary
extern const F_REAL resistive_bound_eta[];

extern const bool eta_Alfven_Wing_boundary;

//! specify distance for eta boundary
extern const F_REAL resistive_bound_dist[];

//! specify how many file shall be loaded
extern const INT32 num_EBF;


//! specify how many coordinates values have to be read for each data point
//! (this indicates whether the dataset boundary values are a 1D, 2D or 3D function)
extern const INT32 num_coords_EBF[];

//! specify how many components the fields inside the files have
extern const INT32 num_comps_EBF[];

//! specify name of the respective file
extern const char file_names_EBF[][RUN_NAME_SIZE];


extern const bool use_hom_particle_bounds[7];

//! for each side of the simulation box has to be specified, which
//! type of boundary conditions shall be used for the particle:
//! true:  inject particle at respective side (inflow)
//! false: do not inject particle at respective side (outflow)

//! 1. entry: x_min-side	2. entry: x_max-side
//! 3. entry: y_min-side	4. entry: y_max-side
//! 5. entry: z_min-side	6. entry: z_max-side
//! in case "use_periodic_bounds" is set true, this parameter
//! will be ignored.
extern const bool use_particle_inflow_bounds[6];

extern const bool use_particle_copy_bounds[6];
//! for each side of the simulation box has to be specified, which
//! type of boundary conditions shall be used for the EM-Fields:
//! true:  set fields constant to background value (inflow)
//! false: set fields derivative to zero (outflow)
//! in case "use_periodic_bounds" is set true, this parameter
//! will be ignored.

extern const bool use_hom_B_bounds[7];

extern const bool set_BFieldInflow_at_VelocityInflow;

extern const bool use_B_inflow_bounds[6];
extern const bool use_E_inflow_bounds[6];
extern const bool use_U_inflow_bounds[6];

#if defined nonadiabatic_gradPE_TERM
extern const bool use_PE_inflow_bounds[6];
#endif



//! specify mesh type
extern const INT32 mesh_type;


extern const bool LF_one_step;





//!-------------------------------------------------------------//
//!------------- 2) Output Parameter ---------------------------//
//!-------------------------------------------------------------//

//! name of simulation run (used as a prefix for each file)
extern const char Run_Name[RUN_NAME_SIZE];

extern const char data_output_path[RUN_NAME_SIZE];


extern INT32  TL;
extern INT32  TL_at_last_restore_state;

//! End Run when TL_MAX is reached
extern const INT32 TL_MAX;

extern const INT32 TL_DUST_MAX;

//! output in units of TL
//! SET 0 TO SWITCH OFF RESPECTIVE OUTPUT

//! vtk    is fileformat for Paraview
//! silo   is fileformat for Visit
//! native is fileformat for own Visualization


//! zonal Data: Each Zone/Cell is rendered using one color value.
//! nodal Data: Each Node is asociated with one color value,
//!		    color of areas inbetween nodes are interpolated.
//! set true for zonal style.
//! set false for nodal style.
extern bool SILO_STYLE_ZONAL;

extern const INT32  TL_OUTPUT_2D_VTK;
extern const INT32  TL_OUTPUT_2D_SILO;
extern const INT32  TL_OUTPUT_2D_NATIVE;
extern const INT32  TL_OUTPUT_2D_GNUPLOT;

extern const INT32  TL_OUTPUT_3D_VTK;
extern const INT32  TL_OUTPUT_3D_SILO;
extern const INT32  TL_OUTPUT_3D_NATIVE;
extern const INT32  TL_OUTPUT_3D_uniform_grid;
extern const INT32  TL_OUTPUT_3D_ASCII;
extern const INT32  TL_OUTPUT_LINEOUT;
extern const INT32  TL_OUTPUT_PARTICLE_DETECTOR;
extern const INT32  TL_OUTPUT_GETMESH;

extern bool COMPRESS_SILO;


//! track path of trajectory and write 1dim field data
extern const INT32  TL_OUTPUT_TRAJECTORY;

//!-----------------------------------------------------------------
//!---------- statefile related -------------------------------
//!-----------------------------------------------------------------
 
//! resubing job if TL_Max is not reached at the end of the Walltime
extern bool resubmit_enable;
extern const INT32 resubmit_walltime;
extern const char resubmit_jobscript[100];
extern const INT32 resubmit_security_factor;

//!-----------------------------------------------------------------
//!---------- LineOut output related -------------------------------
//!-----------------------------------------------------------------

extern const INT32 num_lines;
extern const D_REAL lineout_pos_conversion_fact;

extern const char LineOut_FileName[][RUN_NAME_SIZE];

extern const INT32 num_fields_to_trace_lineout;

extern const INT32 ID_of_Fields_to_trace_lineout[];       

extern const bool pack_parallel_output_to_single_file_lineout;

extern const bool create_new_file_each_TL_lineout;

//!-----------------------------------------------------------------
//!---------- Average fields related -------------------------------
//!-----------------------------------------------------------------



//! specify when averaging shall be finished
//! (most likely equal to output)
//! set zero to turn off averaging
extern const INT32  TL_FINISH_AVERAGE_FIELDS;

//! specify how many time levels averaging shall
//! start before "TL_FINISH_AVERAGE_FIELDS"
extern const INT32  num_TL_average_Fields;



//! specify number of fields to average
extern const INT32 num_average_fields;


extern const char average_field_name_prefix[];

//! specify fields that shall be averaged via their ID
//! as specified in defines.h
extern const INT32 IDs_of_Fields_to_average[];


extern const bool pack_parallel_output_to_single_file;

//! decide whether new file shall be opened each time level
//! or output shall be appended to existing file
//! (e.g. for time series)
extern const bool create_new_file_each_TL;

extern const INT32 num_trajectories;

//! name of trajectory file
extern const char Trajectory_FileName[][RUN_NAME_SIZE];

//! decide whether to trace trajectory or only to display in 2D out
//! (latter one applyies especially to marks such as MP or BS)
extern const bool trajectory_2D_out_only[];

extern const bool do_read_timestring;
extern const bool do_read_SpaceCraftField;

extern const INT32 num_fields_to_trace;
extern const INT32 ID_of_Fields_to_trace[];

extern const D_REAL trajectory_pos_conversion_fact;

//! write file with total enegy, momentum and mass of ions.
extern const INT32  TL_OUTPUT_ErgMomMass;

//! show global information (particles moved, collected, moved ...)
//! since all information have to be gathered collecively
//! together which is expensive, this should not be done to frequently
extern const INT32  TL_LOGFILE_globalINFO;

extern const bool show_proc0_logfile_only;


extern const bool secure_state_file;

//! Only for 2D output
//! Do not cut through GN Planes in depth direction
//! GC surroundig the Blks are always printed !!!
extern const INT32  avoid_GhostCXS;

//! save state every xx seconds
extern const INT32 TL_SAVE_STATE;


extern const INT32 convert_TRACK_PARTICLE_State;


//! in general all processes should write their state file at the same time, however
//! some filesystem cannot handle that much data at once. In such case this parameter 
//! should be set true in order to force the processes to write their state file 
//! one after the other.
extern const bool serialize_writing_of_StateFile;

//! Cross Section to cut (in normed units with respect to origin=(0,0,0))
extern const F_REAL CROSS_SECTION[3];




//! use factor eg. to scale coordinated from normalized to SI units
//! in visit visualisation.
extern const F_REAL factor_scale_mesh;


//! specify label name for VisIt visualization
extern const char visit_xlabel[20];
extern const char visit_ylabel[20];
extern const char visit_zlabel[20];

//! specify length unit of coordinate axis
extern const char visit_length_unit[20];



//! some specifications for uniform_grid output
extern const F_REAL uniform_grid_B0SI;
extern const F_REAL uniform_grid_n0SI;
extern const F_REAL uniform_grid_x0SI;
extern const F_REAL uniform_grid_v0SI;
extern const INT32 numNds_uniform_gridMesh[3];
extern const F_REAL uniform_grid_shrink_domain[3];
extern const INT32 uniform_grid_num_fields_to_trace;
extern const INT32 uniform_grid_ID_of_Fields_to_trace[];
extern const INT32 uniform_grid_Field_type[];


//!--------------------------------------------------------------
//!------------ Track Particle (Test Particle Simulation) -------
//!--------------------------------------------------------------

//!--------------------------------------------
//! NOTE:                                     -
//! TRACK_PARTICLE has to be set in defines.h -
//!--------------------------------------------


//! In test particle simulations any function is deactivated except
//! for particle movement and acceleration
extern const bool TestParticle_Simulation;
//! In test particle simulations any function is deactivated except
//! for particle movement and acceleration

extern const bool INJECT_PARTICLE_TO_TRACK;

extern const INT32  num_mark_particle_in_species[];
extern const D_REAL mark_particle_volume[6];


extern const INT32 TL_INJECT_PARTICLE_TO_TRACK;

// extern const INT32  num_mark_particle_in_species[];

extern const INT32  TL_MARK_PARTICLE_TRACKS;
extern const INT32  TL_OUTPUT_PARTICLE_TRACKS;

extern const INT32 TL_CONVERT_PARTICLE_TRACKS_TO_SILO_MESH;
extern const bool binary_particle_tracks;


extern const INT32 num_tracks_each_group;

//!--------------------------------------------------------------
//!------------ Track Particle (Test Particle Simulation) -------
//!--------------------------------------------------------------

//! Number of detector boxes
extern const INT32 num_particle_detector;

extern const char Detector_FileName[][RUN_NAME_SIZE];

//! Detector Box Geometry
			
extern const D_REAL detector_box_xmin[];
extern const D_REAL detector_box_xmax[];
extern const D_REAL detector_box_ymin[];
extern const D_REAL detector_box_ymax[];
extern const D_REAL detector_box_zmin[];
extern const D_REAL detector_box_zmax[];

//!-------------------------------------------------------------//
//!-------------- 3a) Numerical Parameter -----------------------//
//!-------------------------------------------------------------//



//! Numerical time step at Level 0.
//! NOTE:
//! dt should be at least 5 times smaller than
//! Courant Criteria suggest. Otherwise density 
//! lacks at Blk Level borders occur in flow direction.
extern const D_REAL  dt;


//! Choose advance_B_Algorithm for Plasma Equations:
//! 0: Leap Frog
//! 1: Runge Kutta
extern const INT32 advance_B_Algorithm;


extern const INT32 num_advance_B_loops;


//! strength of EM-Field smoothing
extern const bool several_smooth_in_levels;
extern const D_REAL smooth_E[];
extern const D_REAL smooth_B[];
#if defined(nonadiabatic_gradPE_TERM_smooth_PE)
extern const D_REAL smooth_PE[];
#endif

//! Minimal Charge Density (MCD_J2U) is used in "convert_Ji2Ui" function
//! as an lower limit for rho. Applying MCD_BField results in strongly
//! underestimated velocities.
extern const D_REAL MCD_J2U;

//! Minimal Charge Density (MCD_BField) is used 
//! in "magnetice Field" function as an lower limit for rho.
extern const D_REAL MCD_BField;



//!-------------------------------------------------------------//
//!-------------- 3b) Refinement Parameter ---------------------//
//!-------------------------------------------------------------//

//! Number of maximal refinement-levels used
extern const INT32 MAX_LEVEL;

//! Apply particle splitting and merging every x TL
//! set 0 to switch off
extern const INT32 TL_SPLIT;
extern const INT32 TL_MERGE;

//! decide whether to split exclusively particle of certain
//! species, eg. do not split planetary Ion species
extern const bool do_split_in_species[];
extern const bool do_merge_in_species[];

//! define how many merge/split processes should be carried out in each cell
//! -> the higher oPiC the higer these parameters should be set
extern const int num_merge_each_cell;
extern const int num_split_each_cell;


//! ---- DEFINE HOW TO SPLIT/MERGE -------------------------------
//! decide wheter to split heaviest or most centred particle
//! -> in case most centred is chosen, CoM is conserved by
//!    default an parameters below do not affect splitting
extern const bool split_heaviest_particle;


extern const INT32 num_particle_in_MergeList;

extern const D_REAL merge_tail_distance;

//! to slightly violate centre of mass conservation drastically
//! reduces noise and avoids numerical artifacts
//! -> in summary this seems to reflect the physics much better
//! -> THIS SHOULD ALWAYS BE SET TO FALSE EXCEPT FOR TESTING !!!!

extern const bool conserve_CoM_at_merging;
extern const bool conserve_CoM_at_splitting;


//! new randomly generated positions may be out of cell
//! -> define how often to retry spliting a given particle
extern const INT32 num_randomize_position_split;
extern const INT32 num_randomize_position_merge;

//! level individual activation of splitting/merging
//! in general it should not be required to merge
//! or split in the highest level when "fac_oMPiC_at_LevBoundBlks"
//! is set to 8 or higher
extern const bool do_merge_in_L[];
extern const bool do_split_in_L[];
//! --------------------------------------------------------------

//! Apply Mesh Refinement Algorithm every x TL
//! set 0 to switch off
extern const INT32 TL_REFINE_MESH;




extern const INT32 num_refcrit;
extern const INT32 refcrit_field_IDs[];
extern const INT32 refcrit_field_comp[];
extern const bool refcrit_maximum_based[];

//! decide wether oct_average of respective field shall be normalized
//! by its simulation global extremum of wheter the absolute value
//! shall be used
extern const bool refcrit_normalize_to_global_extremum[];


extern const F_REAL refine_threshold[][NUM_CRITERIA];
extern const F_REAL remove_threshold[][NUM_CRITERIA];



//! increase OPiC in Block neighbours in order to get smoother transition
extern const F_REAL fac_oMPiC_at_LevBoundBlks;


//! decide wheter to smooth or to simply inject moments to parent:
//! Assume L0 is refined to L1 hence particle are in L1.

//! Using the smooth version gives exactly the same density profile for L0,
//! as if particle would have been collected in L0.
//! This only cannot be fullfilled at level boundaries:
//! At minus boundaries for first node inside refined block
//! At plus boundaries for first node outside refined block
//! (can be seen in visu, follow reason using sketch)
//! This in-accuricy appears with and without using gather blocks,
//! so there seems to be no real advantage in using gather blocks.
extern const bool smooth_moments_to_parent;

//! This in-accuricy described above appears with and without using gather blocks,
//! so there seems to be no real advantage in using gather blocks.
extern const bool use_gather_blocks;

//! in order to avoid discontonuities on level boundaries,
//! also environment of flagged block s refined
extern const bool force_refine_environment;

//! in case one Block is flagged, every Block of
//! the respective child_array will be flagged
//! resulting in standard block AMR
extern const bool force_refine_entire_block;

//! Initial Refinement
//! Specify minimal and maximal Coordinate relative to Box Origin
extern const bool never_remove_static_refinement;

extern const INT32 TL_REFINE_STATIC_CUBOID[];
extern const D_REAL minX_refBox_of_L[];
extern const D_REAL minY_refBox_of_L[];
extern const D_REAL minZ_refBox_of_L[];

extern const D_REAL maxX_refBox_of_L[];
extern const D_REAL maxY_refBox_of_L[];
extern const D_REAL maxZ_refBox_of_L[];

extern const INT32 TL_REFINE_STATIC_SPHERE[];
extern const D_REAL radius_refSphere_of_L[];
extern const D_REAL SPHERE_INNER_BOUNDARY[];

extern const INT32 TL_REFINE_STATIC_ZYLINDER_X[];

extern const D_REAL radiusX_refZylinder_of_L[];

extern const D_REAL minX_refZylinder_of_L[];
extern const D_REAL maxX_refZylinder_of_L[];



//! D) Zylinder along y
extern const INT32 TL_REFINE_STATIC_ZYLINDER_Y[];

extern const D_REAL radiusY_refZylinder_of_L[];

extern const D_REAL minY_refZylinder_of_L[];
extern const D_REAL maxY_refZylinder_of_L[];



//! E) Zylinder along z
extern const INT32 TL_REFINE_STATIC_ZYLINDER_Z[];

extern const D_REAL radiusZ_refZylinder_of_L[];

extern const D_REAL minZ_refZylinder_of_L[];
extern const D_REAL maxZ_refZylinder_of_L[];

//!-------------------------------------------------------------//
//!-------------- 3c) MPI Optimization Parameter ---------------//
//!-------------------------------------------------------------//

//! decide whether to order blocks along
//! Space Filling Curve (SFC) or to use the
//! common linear schema
extern bool use_SFC;

//! in order to distribute blocks to intervals of SFC,
//! processing time of each block is recorded.
extern const INT32 TL_reset_block_timing;

extern const INT32 TL_REDISTRIBUTE_BLOCKS;

extern const INT32 OFFSET_REDISTRIBUTE_BLOCKS;

extern const  bool distribute_RB_based;
extern const  bool redistribute_after_restore;

extern const F_REAL multiple_of_average;

extern INT32 max_iteration_redistribute;

//! choose redistibrution criteria:
//! available are:
//! TOTAL_TIME
//! FIELD_TIME
//! BLOCK_NUMBER
//! PARTICLE_TIME
//! PARTICLE_NUMBER

//! 
//! TOTAL_TIME_INC_CHILDREN
//! FIELD_TIME_INC_CHILDREN
//! BLOCK_NUMBER_INC_CHILDREN
//! PARTICLE_TIME_INC_CHILDREN
//! PARTICLE_NUMBER_INC_CHILDREN
extern const INT32 distribution_criteria;


//! Synchronisation of MPI Sends and Recieves
extern const bool sync_mpi_send_rec;


//!-------------------------------------------------------------//
//!-------------- 3d) Simulations Configuration ----------------//
//!-------------------------------------------------------------//

//! These variables allows a fast configuration of A.I.K.E.F. for simple test, e.g. test of the ionisation routine
extern const bool run_calc_first_E;
extern const bool run_CAM;
extern const bool run_calc_second_E;
extern const bool run_accelerate_Particle;
extern const bool run_collect_Ui_minus;
extern const bool run_move_Particle;
extern const bool run_Split_Merge_Particle;
extern const bool run_negative_Particles;
extern const bool run_inject_obstacle_ions;
extern const bool run_collect_Rho_prepare_Recombination_Density;
extern const bool run_chemical_Reactions;
extern const bool run_resize_pArrays;
extern const bool run_collect_Rho_calc_Recombination_Density;
extern const bool run_average_Ui_rho_setREZrho;
extern const bool run_advanceB;


//!-------------------------------------------------------------//
//!-------------- 4) Ion / Ionospheric Parameter ---------------//
//!-------------------------------------------------------------//


//! Number of ion species:
//! This has to be at least one
extern const INT32 num_Charged_Species;

extern const INT32 num_Particle_Species;

// added by HR
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
	//! Number of neutral species:
	//! used as parameters in nonadiabatic_gradPE_TERM computation
	extern const INT32 num_Neutral_Species;

	//! Particle-masses of neutral species (total number of species is defined above)
	extern const D_REAL Neutral_Masses[];

	//! Beta of newly ionised electrons
	extern const D_REAL NewElectron_Betas[];

	//! degrees of freedom of neutral species
	//! for example 3 for a pointmass
	extern const D_REAL Neutral_DOF[];
	
	extern const D_REAL Neutral_Temperature[];
	
	extern const bool analytical_neutral_profile; 	//! modification of neutral dens with formula
	
	extern const bool neutral_profile_from_file; //! neutral dens only by formula
	
	extern const D_REAL norm_externNeutralDensity[];

	extern const D_REAL norm_externNeutralVelocity[];

	//! file names for reading extern fields
	extern const char extern_NeutralField_name[][80];
		
#endif

	extern const D_REAL ReactionRate[NUM_PARTICLE_SPECIES][NUM_PARTICLE_SPECIES][NUM_NEUTRAL_SPECIES];
	
	extern const bool use_velocity_dependent_reaction_rates;
	
#if defined(use_dust_species_as_field)
	
	extern const INT32 num_Dust_Species;

	extern const D_REAL norm_externDustDensity[];

	extern const D_REAL norm_externDustVelocity[];

	//! file names for reading extern fields
	extern const char extern_DustField_name[][80];

#endif
	
	
#if defined(use_ion_production_as_field)
	
	extern const INT32 num_ion_prod_fields;
	
	//! use ion_prod density profile from file
	extern const bool ion_prod_profile_from_file; 
	
	//! file names for reading extern fields
	extern const char extern_IonProdField_name[][80];

	//! calculate the ion production profile of the neutral field
	//! that results from photoabsorption (Chapman profile)	
	extern const bool calc_ionProd_from_neutral_field;

	//! set analytical ion production profile to field
	//! e.g. analytical Chapman profile or modify neutral profile
	//! to consider shadow of moon or planet
	extern const bool set_analytical_ionProd_profile;
	
	extern const D_REAL norm_externIonProdRate[];

	extern const D_REAL norm_externIonProdVelocity[];
#endif	
	

#if defined nonadiabatic_gradPE_COLLISION_TERM

//! collision parameter for electron-neutral-interactions
//! the product of the electron density (ion charge density in normalized units)
//! and the neutral density
//! and the en_coll_param is the factor
//! n_e*m_e /(m_e+m_n) * nu_{e,n}
//! -> see Bachelorthesis (Hendrik Ranocha)
extern const D_REAL en_coll_param[];


#endif

//! Number of ion species emmitted from the box-boundary:
extern const INT32 num_Inflow_Species;

//! indeces of upstream species;
extern const INT32 index_Inflow_Species[];


//! define when start to inject respective species
//! - will be ignored in case is no inflow species
//! - set to zero for continues inflow
extern const D_REAL start_inject_species_each_xx_t0[];

//! define duration of injection
//! - will be ignored in case start is zero
extern const D_REAL duration_inject_species_each_xx_t0[];

//! switch ON/OFF force on ions by neutral Atmosphere
//! define in "NeutralDrag.h"
extern const bool NeutralDrag_Species[];

//! activate recombination for respective ion species
//! recombination rate is calculated in CBlk_Chemical_Reaction.cpp
//! in function void CBlock::calc_RecombinationAlphaField(void)
extern const bool Recombination_for_Species[];

extern const bool check_max_reaction_rates;
extern const bool check_max_reaction_probability;

//! CONCERNING BRACKETS BELOW:
//! Only the first entry up to num_Charged_Species will be used:

//! Aspired number of particle in each cell. Each cell is initialized with this number.
//! During the simulation the number of particles in a cell will changed, so splitting
//! and merging procedures try to adjust the active number to optimal_MPiC.
extern const INT32 optimal_MPiC[];

//! Special Velocity Distribution
extern const bool special_Velocity_Distribution[];

//! Specify Mass, Charge, beta and Obstacle rho for all species:
//!   Mass of 1. represents one upstream ion Mass   (eg. 1 proton mass)
//! Charge of 1. represents one upstream ion Charge (eg. 1 proton Charge)
extern const D_REAL  Ion_Masses[];
extern const D_REAL Ion_Charges[];



//! Ion plasma betas for each Species:
//! THIS VARIABLE MAY ONLY BE USED, IN CASE "rho_sw" AND "Ion_Masses" ARE
//! SPECIFIED FOR THE RESPECTIVE SPECIES (as mass-density is needed to derive
//! thermal velocity from plasma beta).
//! In case "rho_sw" or "Ion_Masses" are not set for respective
//! species (e.g. in case of obstacle ions), directly specify
//! temperature by using "Ion_Thermal_Velocities" variable.
extern const D_REAL Ion_Betas[];

//! Ion-Thermal-Velocities (ITV) for each Species:
//! Instead of specifying a plasma beta, the ITV can be directly specified.
//! - If ITV shall be derived from palsma beta, set Ti to zero.
//! - If Ti shall be zero, set both Ti and plasma beta to zero.
extern const D_REAL Ti_para[];
extern const D_REAL Ti_perp[];



extern const D_REAL Tneutral;

//! cut off velocity for maxwellian distribution used for each species
//! (in thermal velocities vth = 3kT/m of respetive species)
extern const D_REAL max_vth;

//! Define Electron Betas for each Ion species
extern const D_REAL Electron_Betas[];

//! Adiabatic exponent for equation of state
extern const D_REAL kappa_electron;

extern const D_REAL PhotoionisationRate[][NUM_NEUTRAL_SPECIES];

extern const bool ElectronionisationRate_Te_dependent;
extern const D_REAL ElectronionisationRate[][NUM_NEUTRAL_SPECIES];

extern const D_REAL Global_IonProduction_Rate[][NUM_NEUTRAL_SPECIES];

//! num macro Particle each TL
//! representing the heavy ions
extern const INT32 obstacle_MP_num_each_TL[];

//! specify number of "lookup" density fields for Ion injection
extern const INT32 num_externRhoVelocityFields;
 
extern const INT32 index_externRho_Species[];

//! normalization of extern Density:
extern const D_REAL norm_externRho[];

//! normalization of extern Velocity:
extern const D_REAL norm_externUi[];

extern const bool phi_cloud_from_extern;

//! file names for reading extern fields
extern const char extern_Field_name[][RUN_NAME_SIZE];


extern const bool read_extern_field_from_state;

extern const bool serialize_reading_of_extern_field;

extern const INT32 num_smooth_extern;
extern const D_REAL smooth_extern[];

extern const INT32 TL_read_extern;

//!-------------------------------------------------------------//
//!-------------- 5) Background Plasma Parameter ---------------//
//!-------------------------------------------------------------//

//! IMF (Magnetic Background BField)
extern const D_REAL B_sw[3];

//! set new IMF at inflow box boundaries when TL_new_Bsw is reached
//! (e.g. sector border, magnetopause crossing)
extern const D_REAL new_B_sw[3];

//! specify time level for new IMF
extern const INT32 TL_new_Bsw;

//! in order to avois sharp transition rotate old to
//! new field within finite time
extern const INT32 num_TL_adapt;

//! Solar Wind Plasma Velocity
extern const D_REAL V_sw[3];

//! Plasma Densities for each inflow species.
//! (values for at index of obstacle species will be ignored)
extern const D_REAL rho_sw[];

//! plasma resistivity
extern const D_REAL Eta_sw;


//!-------------------------------------------------------------//
//!-------------- 6) Obstacle Parameter ------------------------//
//!-------------------------------------------------------------//

extern const bool use_resistive_obstacle;


//! for two body simulations: activate second obstacle
extern const bool use_second_obstacle;
extern const D_REAL R_SecondObstacle;
extern const D_REAL Position_SecondObstacle[3]; 

//! number & strength of smoothing the Resistivity profile
extern const INT32  num_smooth_eta;
extern const D_REAL smooth_eta[];

extern const D_REAL omega_rotating_obs[3];

//! Magnetic Field calculation inside the Obstacle (BO).
//! OM1) Switch on/off calculation inside the Obstacle.
extern const bool advance_obstacle_B;

//! OM2) choose relaxations parameter for SOR
extern const D_REAL B_SOR_omega_L[];

//! OM3) choose break criteria for SOR
extern const D_REAL B_SOR_max_error;

//! DC4) B_SOR_num_cycle_base^L cycles will be performed
//!      in level L each iterations
extern const INT32 B_SOR_num_cycle_base;

//! OM5) choose maximal number of iterations for SOR
extern const INT32 B_SOR_max_Iteration;

//! OM6) calc error every xx TL
extern const INT32 B_SOR_calc_error_step;

//! Divergence Cleaner (DC)
//! DC1) Switch on/off BField-divergence cleaning inside the Obstacle.
extern const bool div_cleaner;

//! DC2) choose relaxations parameter for SOR
extern const D_REAL DC_omega_L[];

//! DC3) choose break criteria for SOR
extern const D_REAL DC_max_error;

//! DC4) DC_num_cycle_base^L cycles will be performed
//!      in level L each iterations
extern const INT32 DC_num_cycle_base;

//! DC5) choose maximal number of iterations for SOR
extern const INT32 DC_max_Iteration;

//! DC6) calc error every xx TL
extern const INT32 DC_calc_error_step;


//! check if different values converge faster
//! TL starting the optimization
//! (testing  values from B_SOR_omega_L-steps/2 in steps of 0.01)
extern const INT32 TLstart_optimize_B_SOR_omega;
extern const INT32 num_optimize_SOR_steps;





//! Switch on/off electric field calculation inside the Obstacle.
extern const bool advance_obstacle_E;


//! within the core not magnetic field is advanced nor smoothed
//! (-> it will remain the initial field forever)
//! in % of obstacle radius
extern const D_REAL obstacle_core_fraction;

//! Obstacle's intrinsic Dipole Field
//! (so far it points in z direction)
extern const D_REAL Magnetic_Moment[3];
extern const D_REAL MM[3];

extern const bool is_intern_Magnetic_Moment;

extern const D_REAL Magnetic_Moment_Angle[2];

extern const D_REAL Magnetic_Moment_Offset[3];


//! densities for each ion species to eliminate strong electron pessure 
//! Gradients at Obstacle surface.
extern const D_REAL obstacle_rho[];




//!-------------------------------------------------------------//
//!--- misc Variable -> declarations in absolute_Globals.h -----//
//!-------------------------------------------------------------//
//! For horizontal search
extern CBlock**  BlockList_of_Lev;
extern CBlock**  GATHER_BlockList_of_Lev;


extern INT32 total_active_Blocks;

extern const D_REAL min_pArray_Size_factor;
extern const D_REAL max_pArray_Size_factor;
extern const INT32 TL_RESIZE_pARRAYS;

extern const D_REAL ram;
extern const D_REAL thermal;
extern const D_REAL mag;



extern const INT32 num_nodes_in_block, num_root_blocks;



