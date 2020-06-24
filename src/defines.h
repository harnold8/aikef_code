//when defined, the dipole and/or the atmosphere is off
//#define dipoleoff
//#define atmosphereoff

//on
#ifndef atmosphereoff
  //! number of species that are treated as particles                                                                                                                                                                                       
  #define NUM_PARTICLE_SPECIES 2
  //! for chemical reactions array                                                                                                                                                                                                          
  #define NUM_NEUTRAL_SPECIES 1
//off
#else
  #define NUM_PARTICLE_SPECIES 1                                                                                                                                                                                                            
  #define NUM_NEUTRAL_SPECIES 1
#endif


#define use_neutral_species_as_field
#define use_ion_production_as_field
// #define use_dust_species_as_field


//! -------------- FIELD RELATED -----------------------------------------

//dipole=on
#ifndef dipoleoff
  //! define this to use a dipole field as background
  #define USE_CFBG_BFIELD
#endif

//! define this to use an inhomogeneous background magnetic field
//#define use_inhom_B

//! Define this to use a better smooth_Field method
#define use_new_smooth_Field_method


#define NUM_SUB_CYCLE  7


//! defines to switch off/on physical parts
//! of the E/BField Equation
#define CONVEC_TERM
#define HALL_TERM
//#define gradPE_TERM

//! switch on/off use eta in Leap Frog method.
//! Must be off if advance_obstacle_B
//! in parameter.cpp is set true
// #define ETA_TERM_BField

//! switch on/off use eta in EField method.
// #define ETA_TERM_EField


// #define nonadiabatic_gradPE_TERM
// #define nonadiabatic_gradPE_COLLISION_TERM
// #define nonadiabatic_gradPE_SOURCE_TERM
// #define nonadiabatic_gradPE_DRAIN_TERM

//! Set nonadiabatic_gradPE_TERM_derivation_form to
//! 0 : non-conservative form
//! 1 : conservative form
//! 2 : conservative form using extra field (id_scratch_scalar) for pow(PE, 1./kappa_electron)
#define nonadiabatic_gradPE_TERM_derivation_form 2
// #define nonadiabatic_gradPE_TERM_set_divUe_to_zero_DEBUG

//! Set nonadiabatic_gradPE_TERM_curl_derivation_form to
//! 0 : Compute rot( 1/rho * grad(Pe) )
//! 1 : Compute - rot( Pe * grad(1/rho) )
//! 2 : Compute grad(1/rho) x grad(Pe)
//! If rot(grad(...)) = 0, all terms are equivalent
#define nonadiabatic_gradPE_TERM_curl_derivation_form 1


//! define this if you want to smooth PE similar to E and B
#define nonadiabatic_gradPE_TERM_smooth_PE

//! define this to ensure PE >= 0
#define nonadiabatic_gradPE_TERM_ensure_PE_not_negative

//! Only use in case of Debug
//! If you use conservative form from above
//! define nonadiabatic_gradPE_TERM_adiabatic_term_only
//! to use only the adiabatic term in the time advancement
// #define nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
// #define nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only

//! Set nonadiabatic_gradPE_initial_condition to
//! 1 if you want to have homogenous electron pressure as initial value
//!   the initial value for the electron pressure will be set via Electron_Betas[0]
//! 0 if you want to have homogenous electron temperature as initial value
//!   the initial value for the electron pressure will be computed
//!   using the ion density and Electron_Betas[0]
//!   PE[u_v_w] = RHO[u_v_w] * Electron_Betas[0]
//! NOTE: good old AIKEF with adiabatic electron pressure assumes
//!       homogeous electron pressure defined by Electron_Betas
//!       as initial value
#define nonadiabatic_gradPE_initial_condition 0

//! If you use nonadiabatic_gradPE_TERM
//! you have to decide, whether Pe shall be advanced
//! together with B or in an extra loop using 
//! the explicit midpoint method.
//! In the first case, define nonadiabatic_gradPE_TERM_advance_Pe_with_B
#define nonadiabatic_gradPE_TERM_advance_Pe_with_B



//! If kappa_electron is 2.0, you should define this
//! to get better performance, especially with vectorclass
//! defined below
#define kappa_electron_is_2


//! Define use_vectorclass if you want to use the optimized
//! vectorclass from http://www.agner.org/optimize/#vectorclass
//!
//! If use_vectorclass is defined,
//! NOTE: num_nodes_in_block has to be a multiple of 4
//! because of CBlock::sub_squaredField_Multiply
//#define use_vectorclass


#if defined(use_vectorclass)

//! Use Intels Short Vector Math library (SVML)
//! see vectorclass/VectorClass.pdf, page 32 etc.
//! for details about it
//! TODO: Does not work correctly up to now...
#define VECTORMATH 2

#include "vectorclass/vectorclass.h"
#include "vectorclass/special/vectormath.h"

//! NOTE: If you change something in the other typedefs,
//! you have to make the same changes here!
typedef Vec4d VEC4_PARTICLE_REAL;
typedef Vec4f VEC4_F_REAL;
typedef Vec4d VEC4_D_REAL;
typedef Vec2d VEC2_D_REAL;

#endif


//! ----------------------------------------------------------------------


//! to run on NESS and HECToR below has to be defined
#define MPICH_IGNORE_CXX_SEEK

//! choose between buffered and synchronized send
//! - set to "Isend" for performance code
//! - set to "Issend" or "Ibsend" for debugging
#define IB_OR_IS_SEND	Isend




//! -------------- PARTICLE RELATED --------------------------------------
//! JUST FOR STATISTICS:
//! FOR SOME REASON CAN SLOW SIGNIFICANTLY DOWN THE CODE !!!
#define  ESTIMATE_FASTEST_PARTICLE

//! MAX_LEVEL_DEF
//! to determine the fastest particle per level set min. MAX LEVEL 
#define MAX_LEVEL_DEF 7

//! Delete particle at very high velocities which
//! do not fulfill the courant criteria and thereby
//! will make the move method crash
//! -> this should rather be used for debugging
//!    since ANY particle should fulfill the CC
// #define DELETE_PARTICLE_IF_V_EXCEEDS 60


//! define if particle shall be sorted by weight into cell in order to
//! accelerate merging procedure. Since code does barely decelerate
//! with sorting, this should be always switched on.
#define SORT_BY_WEIGHT

//! mark and trace particles
#define TRACK_PARTICLE

//! ----------------------------------------------------------------------

#define DO_PARENT_BUFFER_UPDATE		true
#define NO_PARENT_BUFFER_UPDATE		false

#define INCL_PHYS_FACE			true
#define SKIP_PHYS_FACE			false

//! The more Ion Species are allocated, the more
//! Fields must be allowed

#if (defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)) && defined(use_dust_species_as_field) && defined(use_ion_production_as_field)

#define NUM_FIELDS  	114
#define RUN_NAME_SIZE   40
#define NUM_CRITERIA  	5
#define NUM_SCALAR_FIELDS_EACH_BLK  (49 +13*num_Charged_Species + NUM_PARTICLE_SPECIES+ 6*num_Neutral_Species+ 4*num_Dust_Species+5*num_ion_prod_fields+4*num_externRhoVelocityFields +num_scalar_average_fields)

#elif (defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)) && defined(use_dust_species_as_field) && !defined(use_ion_production_as_field)

#define NUM_FIELDS  	114
#define RUN_NAME_SIZE   40
#define NUM_CRITERIA  	5
//! TODO NOTE I have absolutely no clue why it has to be 5*num_Dust_Species instead of 4*num_Dust_Species, but 4 produces a seg fault
#define NUM_SCALAR_FIELDS_EACH_BLK  (49 +13*num_Charged_Species + NUM_PARTICLE_SPECIES+6*num_Neutral_Species+ 5*num_Dust_Species+4*num_externRhoVelocityFields +num_scalar_average_fields)

#elif (defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)) && defined(use_ion_production_as_field)  && !defined(use_dust_species_as_field)

#define NUM_FIELDS  	200
#define RUN_NAME_SIZE   40
#define NUM_CRITERIA  	5
//! TODO NOTE I have absolutely no clue why it has to be 5*num_ion_prod_fields instead of 4*num_ion_prod_fields, but 4 produces a seg fault
#define NUM_SCALAR_FIELDS_EACH_BLK  (100 +13*num_Charged_Species + 6*num_Neutral_Species+ 5*num_ion_prod_fields+4*num_externRhoVelocityFields +num_scalar_average_fields)

#elif (defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)) && !defined(use_dust_species_as_field) && !defined(use_ion_production_as_field)

#define NUM_FIELDS  	114
#define RUN_NAME_SIZE   40
#define NUM_CRITERIA  	5
#define NUM_SCALAR_FIELDS_EACH_BLK  (49 +13*num_Charged_Species + NUM_PARTICLE_SPECIES+ 6*num_Neutral_Species+4*num_externRhoVelocityFields +num_scalar_average_fields)

#else

#define NUM_FIELDS  	114
#define RUN_NAME_SIZE   40
#define NUM_CRITERIA  	5
#define NUM_SCALAR_FIELDS_EACH_BLK  (49 +13*num_Charged_Species+ NUM_PARTICLE_SPECIES +4*num_externRhoVelocityFields +num_scalar_average_fields)

#endif




//! unfortunately the IBM compiler does not
//! support dynamic stack allocation of strings
//! so leave this value constant
#define INFO_ARRAY_SIZE 30

//! perform aditional interrogations to check for rarely occuring errors
// #define DEBUG_MODUS

//! protocol particle that impact on the obstacle
// #define startTL_PROTOCOL_OBSTACLE_PARTICLE	3
// #define TL_PROTOCOL_ANY_PARTICLE	2


//! real Inprecission interception:
//! used in move and insertion procedure
//! should be 1.e-6 in case PARTICLE_REAL is float
//! should be 1.e-16 in case PARTICLE_REAL is double
//! (see comments in CBlock::move_particle for details)
#define PART_REAL_PRECISION		1.e-16

#define BUILD_SUM	0
#define BUILD_MAX	1


#define UNIFORM		0
#define STAGGERED	1



//! define SORTED_PARTICLE_LIST for faster merge procedure:
//! if list are weight sorted it is sufficiant to compare
//! first ~40 particles.
//! As this procedure does not need significant
//! time, it should always be defined
#define SORTED_PARTICLE_LIST



//! In case Neumann (outflow) Boundaries are used:
//! MARS  (fast plasma): set NM_BOUND_NODE to 1
//! DIONE (slow plasma): set NM_BOUND_NODE to 2
//! For some reason, Mars-Run got instable (NaN in B)
//! in case NM_BOUND_NODE was set to 2.
//! In case of slow plasma both cases do work, but in 
//! case of 1 strong wave reflection occured at the 
//! outflow boundary. This could be strongly decreased
//! by setting source to 2.
#define NM_BOUND_NODE	2



//! As "int" is 4 BYTE on all todays architechturs (which
//! results in a Range of [-2e9;+2e.9]), use int and avoid
//! using long.



//! All variable types are replaced with own types for 
//! the following reasons:
//! 1) If an error accurs, by setting all types to double precision
//!    an memory overflow can be easyly excluded
//! 2) Different systems use different precison for the same type 
//! 	(eg. int 2byte or 4 byte).
//! 3) Plotting File Stream is set to single Precesion by default.
//!    This can be easily replaced if necessary.


//! integer datatypes
#define MPI_INT64  MPI_LONG_LONG
typedef long long  INT64;

#define MPI_INT32  MPI_INT
typedef int  INT32;



//! floating point variables:
#define MPI_D_REAL  MPI_DOUBLE
typedef double  D_REAL;

#define MPI_F_REAL  MPI_FLOAT
typedef float  F_REAL;

#define MPI_PARTICLE_REAL  MPI_DOUBLE
typedef double  PARTICLE_REAL;

#define MPI_WEIGHT_REAL  MPI_DOUBLE
typedef double  WEIGHT_REAL;

//! always use float for FILE_REAL, else
//! silo functions will return errors
#define MPI_FILE_REAL  MPI_FLOAT
typedef float  FILE_REAL;



#define NUM_REQUESTS	16

#define gathNEIB	6

//! Box - Boundarie Abbriviations
#define Xmin_BB	0
#define Xmax_BB	1
#define Ymin_BB	2
#define Ymax_BB	3
#define Zmin_BB	4
#define Zmax_BB	5

#define _Im1_	0
#define _Ip1_	1
#define _Jm1_	2
#define _Jp1_	3
#define _Km1_	4
#define _Kp1_	5
#define _ANY_	6
#define _INIT_	6
#define _OBS_	7

//! RESTORE STATE CONVERTION:
#define LEAVE_UNCHANGED     0
#define CONVERT_TO_TRACK    1
#define CONVERT_TO_ORDINARY 2


//! For gathering ... 
#define ALL_SPECIES -1

#define noVTH2  false
#define getVTH2  true

//! Block redistribution
#define NUM_REDIST_OPTIONS		11
#define AVERAGE_REF_VALUES		 0
#define TOTAL_TIME			 1
#define FIELD_TIME			 2
#define BLOCK_NUMBER			 3
#define PARTICLE_TIME			 4
#define PARTICLE_NUMBER			 5

#define TOTAL_TIME_INC_CHILDREN		 6
#define FIELD_TIME_INC_CHILDREN		 7
#define BLOCK_NUMBER_INC_CHILDREN	 8
#define PARTICLE_TIME_INC_CHILDREN	 9
#define PARTICLE_NUMBER_INC_CHILDREN	10



//! uniform_grid output
#define uniform_grid_TYPE_BFIELD        0
#define uniform_grid_TYPE_EFIELD        1
#define uniform_grid_TYPE_DENSITY_FIELD    2
#define uniform_grid_TYPE_VELOCITY_FIELD    3



//! using short for INT32 does not speed up the code
//! (maybe 1% at most)


//!------------- Geometry Field Flags -------------------//
#define FILE_TYPE_INFO  	200
#define FILE_TYPE_BLOCKS 	201



//!------------- Field IDs -------------------//

#define id_notDefined	        -1


//!-------------------------------------
//!-- PHYSICAL FIELDS FOR CALCULATION --
//!-- (-> allocate memory for each)   --
//!-------------------------------------
//! EM Field ids
#define id_BTotal		 0
#define id_BEven 		 1
#define id_BOdd			 2
#define id_Bcfbg		 3
#define id_EField		 4


//! density Field ids
#define id_Lam			 5
#define id_rho_n		 6
#define id_rho_np1 		 7
#define id_rho_rez		 8
// #define id_allRhoSpecies 	 9

//! ion mean velocity Field ids
#define id_Gam			9
#define id_UI_plus		10
#define id_UI_minus		11

//! misc field id's
#define id_Eta 			12
#define id_PhiDC		13
#define id_BDerivative		14

//! Field Groups
#define id_ALL_FIELDS		15
#define id_UIp_Gam_allRhoSpec 	16



//!------------------------------------------
//!-- FIELDS EXCLUSIVELY FOR VISUALIZATION --
//!--   (-> share one "scratch" field )    --
//!------------------------------------------

//! derived physical fields
#define id_gradPI0		17
#define id_rotB 		18
#define id_divB 		19
#define id_divE 		20
#define id_PMagnetic		21
#define id_PTotal		22

//! time 
#define id_FieldTime		23
#define id_ParticleTime		24
#define id_scratch_vector	25
#define id_scratch_scalar	26

#define id_ExB			25 //just used as scratchpad by HK-Code

//! refinement
#define id_gradB		27
#define id_Refine_Rating	28
#define id_ElectronTemperature  29
#define id_rho_np1_recombined   30

#define id_Bcfbg_HIMesh         31


//! species related
#define id_allRhoSpecies        32
#define id_allPISpecies         33
#define id_allGyroSpecies       34
#define id_allGyroELSpecies     35
#define id_allPESpecies         36
#define id_allUISpecies         37
#define id_allForcesSpecies     38


//! NOTE:
//! all numbers from 

//! id_rhoSpecies1 to id_rhoSpecies1 +5*num_Charged_Species-1
//! are used for each species density, velocity and temperature !!!

//! E.g. in case of 2 ion species:
//! id 28: rho of species 1
//! id 29: rho of species 2
//! id 30: velocity of species 1
//! id 31: velocity of species 2
//! id 32: gyroradius of species 1
//! id 33: gyroradius of species 2
//! id 34: gyroradius of electrons species 1
//! id 35: gyroradius of electrons species 2
//! id 36: temperature of species 1
//! id 37: temperature of species 2

//!----------------------------------------
//!- DO NOT USE THESE NUMBERS ELSEWHERE!  -
//!----------------------------------------
#define id_rhoSpecies1          39
#define id_UI_Species1		(id_rhoSpecies1   +num_Charged_Species)
#define id_gyro_Species1        (id_UI_Species1   +num_Charged_Species)
#define id_gyro_el_Species1     (id_gyro_Species1   +num_Charged_Species)
#define id_recomb_Species1	(id_gyro_el_Species1   +num_Charged_Species)
#define id_ForceSpecies1	(id_recomb_Species1   +num_Charged_Species)
#define id_PISpecies1		(id_ForceSpecies1 +num_Charged_Species)
#define id_PESpecies1		(id_PISpecies1    +num_Charged_Species)

#define id_PEtotal		(id_PESpecies1  +num_Charged_Species)
//! neutral particles are used as parameters in
//! nonadiabatic_gradPE_TERM computations
//! NOTE: You may have to adjust NUM_FIELDS defined above
//!       if you use nonadiabatic_gradPE_TERM
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)

	//#define id_allRhoNeutralSpecies                 (id_PEtotal + 1)
	#define id_allRhoNeutralSpecies                 (id_PESpecies1  +num_Charged_Species +1)
	//#define id_numberdensity_neutralSpecies1	(id_allRhoNeutralSpecies  + 1)
	#define id_numberdensity_neutralSpecies1        (id_PESpecies1  +num_Charged_Species +2)
	#define id_allUNeutralSpecies                   (id_numberdensity_neutralSpecies1 + num_Neutral_Species)
	#define id_velocity_neutralSpecies1		(id_numberdensity_neutralSpecies1 + num_Neutral_Species + 1)
	#define id_allPNeutralSpecies                   (id_velocity_neutralSpecies1 + num_Neutral_Species)
	#define id_pressure_neutralSpecies1		(id_velocity_neutralSpecies1 + num_Neutral_Species + 1)
	#define id_allnewBetaNeutralSpecies             (id_pressure_neutralSpecies1 + num_Neutral_Species)        
	#define id_new_electron_beta_neutralSpecies1	(id_pressure_neutralSpecies1 + num_Neutral_Species + 1)

	#define id_externRho1		(id_new_electron_beta_neutralSpecies1 + num_Neutral_Species)

	//! Use id_PESpecies1 as temporary scratch for PE_odd
	#define id_PE_even	id_PEtotal
	#define id_PE_odd	id_PESpecies1

#else

	#define id_externRho1		(id_PESpecies1  +num_Charged_Species +1)

#endif


#define id_extern_Ui1		(id_externRho1  +num_externRhoVelocityFields)


#if defined(use_dust_species_as_field) && defined(use_ion_production_as_field)

	#define id_density_ionProdSpecies1	(id_extern_Ui1  +num_externRhoVelocityFields)
	#define id_velocity_ionProdSpecies1	(id_density_ionProdSpecies1  +num_ion_prod_fields)

	#define id_density_dustSpecies1        (id_velocity_ionProdSpecies1  +num_ion_prod_fields)
	#define id_velocity_dustSpecies1	(id_density_dustSpecies1 + num_Dust_Species)

	#define id_average_Field1 	(id_velocity_dustSpecies1 + num_Dust_Species)

#elif !defined(use_dust_species_as_field) && defined(use_ion_production_as_field)

	#define id_density_ionProdSpecies1	(id_extern_Ui1  +num_externRhoVelocityFields)
	#define id_velocity_ionProdSpecies1	(id_density_ionProdSpecies1  +num_ion_prod_fields)

	#define id_average_Field1 	(id_velocity_ionProdSpecies1 + num_ion_prod_fields)

#elif defined(use_dust_species_as_field) && !defined(use_ion_production_as_field)

	#define id_density_dustSpecies1        (id_extern_Ui1  +num_externRhoVelocityFields)
	#define id_velocity_dustSpecies1	(id_density_dustSpecies1 + num_Dust_Species)

	#define id_average_Field1 	(id_velocity_dustSpecies1 + num_Dust_Species)
#else

	#define id_average_Field1	(id_extern_Ui1  +num_externRhoVelocityFields)

#endif 



#define id_divrotB		(id_average_Field1 + num_average_fields +1)
#define id_divU			(id_divrotB + 1)
#define id_gradPE		(id_divU + 1)

#define id_BxgradB		(id_gradPE+1)
#define id_RespProc		(id_BxgradB+1)

//! refinement
#define id_Refine_Status	(id_RespProc+1)

#define id_extern_PhiC	(id_Refine_Status+1)

//! Note:
//! Dot not use numnbers higher 
//! than NUM_FIELDS defined above !!!
#define id_BOld			id_BOdd
#define id_BNew			id_BEven
#define id_KX			id_rotB
#define id_BMidStep		id_rotB
#define id_B_HI			id_rotB

#define id_rho_np1_HIMesh	id_PESpecies1
#define id_rho_n_HIMesh 	id_PISpecies1

//! define id for output of energy spectrum
//! (has to be larger than highest field id )
#define ENERGY_SPECTRUM (id_extern_PhiC+1)



