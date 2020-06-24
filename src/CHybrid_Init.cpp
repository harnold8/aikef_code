


#include <math.h>
#include "utils.h"
#include "CHybrid.h"
#include "parameters.h"
#include "hilbert.h"
#include "gsl/gsl_rng.h"
#include "absolute_Globals.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctime>



extern D_REAL *dt_particle_of_L, **delta_of_L;


//! extern array of function pointers for
//! cp GN procedure (initialize in this file)
extern void (*send_GN_MPI[3])(bool direc,
				CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32 id_package,
				INT32 field_type);

extern void (*recv_GN_MPI[3])(bool direc,
				CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32 id_package,
				INT32 field_type);


extern void (*get_GN_equal_process[3])(bool direc,
					CBlock* dest_Blk,
					 CBlock* src_Blk,
					 INT32    field_type);

//! extern array of function pointers for
//! add GN procedure
extern void (*send_add_GN_MPI[3])(CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32 id_package,
				INT32 field_type);

extern void (*recv_add_GN_MPI[3])(CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32 id_package,
				INT32 field_type);


extern void (*get_add_GN_equal_process[3])(CBlock* dest_Blk,
					 CBlock* src_Blk,
					 INT32    field_type);


extern void (*BV_from_parent[6])(CBlock* dest_Blk,
				  INT32    field_type);





//!--------------------------------------------------------
//!- init
//!  Main Initialization Function
//!--------------------------------------------------------
void CHybrid::init(void)
{


	log_file<<"id_extern_PhiC : "<<id_extern_PhiC<<endl;
	//! check whether user specified variables conflict
	log_file<<"entering check sucessful"<<endl;
	do_consistency_check_init_derived_variables();
	log_file<<"consistency check sucessful"<<endl;
	//! Call this functions at the very Beginning
   	set_Field_Comps();
	log_file<<"set field components sucessful"<<endl;
   	set_Field_Names_IDs();
	log_file<<"set field names sucessful"<<endl;
	//! write infostructure for native file format
	//! (only master process)
	log_file<<"attempting to write info structure"<<endl;
	write_Info_Structure();
	log_file<<"successful to write info structure"<<endl;
	//! init random generators
	alloc_init_random_Generators();

	//! Read parameters first !
	log_file<<"attempting to alloc variables"<<endl;
	alloc_variables();
	log_file<<"successful to alloc variables"<<endl;
	

	//! Read Extern Boundary Fields if specified
	read_EBF();

	//! set Field IDs before RootBlks are created
	alloc_RootBlks();

	//! restore Blks before children are created
    	if(!restore_state())
	{
	   //! initialize position in box, memory and set fields
	  //! only in case they are not restored
	  init_RootBlocks();
	
	}
	else
	{
	  if(secure_state_file)
	    secure_state_files();
	}




	//! create initially refined Sphere built out of blocks
	//! or add these blocks in case of restore if not yet refined
	//! -> also mark respective Blocks in case of restore
	//!    as initial refined
	static_refinement_ZylinderX();
	static_refinement_ZylinderY();


        static_refinement_sphere();

	static_refinement_ZylinderZ();

	//! create rectangular initially refined rect built out of blocks
	//! or add these blocks in case of restore if not yet refined
	//! -> also mark respective Blocks in case of restore
	//!    as initial refined
        static_refinement_cuboid();


	//! required for MPI send / recv
	//! always call this fucntion when
	//! mesh was changed
	post_mesh_refinement();


	//!NOTE:
	//! Aways redistribute BEFORE box is filled, otherwise particle might not fit
	//! on every processor -> CRASH / SWAP 
	//! redistribute by means of optimal_PiC in case box ist not
	//! filled but PARTICLE_NUMBER speciefied
	//! (optimal_PiC is automatically set in method)
	if((TL==0 || redistribute_after_restore) && (distribution_criteria == BLOCK_NUMBER || distribution_criteria == PARTICLE_NUMBER))
	redistribute_blocks(-1);

	//! Init Eta Profile
	init_Eta_Profile();


	//! Smooth Eta Field after children are generated either by means of
	//! restore or initial generation
 	for(INT32 a=0; a<num_smooth_eta; a++)
 	smooth_Field(id_Eta,smooth_eta);


	if(TL_INJECT_PARTICLE_TO_TRACK==TL)
	//inject_particle_to_track(0);
	//inject_particle_to_track(0);
	//inject_particle_to_track(1);
	//inject_particle_to_track(2);
	//inject_particle_to_track(3);
	  inject_callisto_particle_to_track(3);

	 //! initialize global variables such as negative_particles_in_simu and reaction_rate[species]
	//! this should be done here, to avoid problems when simu is started from state file
	init_global_variables();

}



void CHybrid::set_fields(void)
{

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->set_Fields();
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}
//!--------------------------------------------------------
//!- init_Eta_Profile:
//!--------------------------------------------------------
void CHybrid::init_Eta_Profile(void)
{


	log_file <<  endl;
	log_file << " SETTING ETA PROFILE ... " << endl;


	//! set Eta Profile
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->init_Eta_Profile();

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}

	log_file << " done. " << endl << endl;


}

// added by HR
#if defined nonadiabatic_gradPE_TERM
//!--------------------------------------------------------
//!- init_PE_Profile:
//!--------------------------------------------------------
void CHybrid::init_PE_Profile(void)
{


	log_file <<  endl;
	log_file << " INITIALIZING PE PROFILE ... " << endl;


	//! set PE Profile
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->init_PE_Profile();

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}

	log_file << " done. " << endl << endl;


}
#endif


//!--------------------------------------------------------
//!- read_EBF:
//!--------------------------------------------------------
void CHybrid::read_EBF(void)
{

	if(!num_EBF)
	return;

	log_file <<  endl;
	log_file << " READING EXTERN BOUNDARY FIELDS (EBF) ... " << endl;
	

	values_EBFs = new D_REAL*[num_EBF];
	num_rows_EBFs = new INT32[num_EBF];

	D_REAL buffer;
	char axis[3][3] = {"X","Y","Z"};

	//! loop across all fields to read
	for(INT32 field=0; field<num_EBF; field++)
	{

		INT32 num_columns = num_coords_EBF[field]+num_comps_EBF[field];

		log_file << "  Reading " << file_names_EBF[field] << " ..." << endl;
		log_file << "   -> num_coords: " << num_coords_EBF[field]  << endl;
		log_file << "   -> num_comps:  " << num_comps_EBF[field]  << endl;
		log_file << "   -> file is assumed to consist of " << num_columns <<" columns " << endl;



	
		ifstream infile;
		infile.open(file_names_EBF[field]);

		if(!infile)
		{
			log_file << "   no " << file_names_EBF[field] <<" File found." << endl;
			log_file << "   Exiting ..." << endl;
			exit(1);
			
		}

		//! It is implicitely assumed that extern boundary conditions are 1 dimensional
		//! which means that they are either a function of X or Y or Z, but not a
		//! function of eg. X and Y at the same time.

		//! -> hence each file most contain of 1+num_comps_EBF columns, where the first column
		//!    is the coordinate in either X,Y or Z direction


		//! read file to buffer for counting
		num_rows_EBFs[field] = 0;
		for(;;)
		{

			for(INT32 coord=0; coord<num_coords_EBF[field]; coord++)
			infile >> buffer;

			if(infile.eof()) break;

			for(INT32 comp=0; comp<num_comps_EBF[field]; comp++)
			infile >> buffer;

			num_rows_EBFs[field]++;

		}
		infile.close();

		log_file << "   -> rows in file estimated: " << num_rows_EBFs[field] << endl;


		//! alloc 1 column for coordinate plus num_comps_EBF[field]
		values_EBFs[field] = new D_REAL[num_rows_EBFs[field]*num_columns];

		log_file << "   -> writing file to array ..." << endl;
		infile.open(file_names_EBF[field]);
		for(INT32 row=0; row<num_rows_EBFs[field]; row++)
		{
		


			for(INT32 coord=0; coord<num_coords_EBF[field]; coord++)
			infile >> values_EBFs[field][row +coord*num_rows_EBFs[field]];

			if(infile.eof())
			{
				log_file << "   -> error while reading file at row: " << row << endl;
				log_file << "   -> Exiting ... " << endl;
				exit(1);
			}
	
			for(INT32 comp=0; comp<num_comps_EBF[field]; comp++)
			infile >> values_EBFs[field][row +(num_coords_EBF[field]+comp)*num_rows_EBFs[field]];
	

		}


		infile.close();
		log_file << "   done" << endl;
		


	}

	log_file <<  endl;

}





//!--------------------------------------------------------
//!- set_Field_Comps_Names_IDs: 
//!--------------------------------------------------------
void CHybrid::set_Field_Comps(void)
{

      //!--------------------------------------------
      //!------ COMPS -------------------------------
      //!--------------------------------------------
      //! COMPs_FType defined in absolut_Globals.h

	COMPs_FType = new INT32[NUM_FIELDS];

	COMPs_FType[id_BEven ] = 3;
	COMPs_FType[id_BOdd  ] = 3;
	COMPs_FType[id_EField] = 3;

	COMPs_FType[id_UI_plus] = 3;
	COMPs_FType[id_UI_minus] = 3;

	COMPs_FType[id_rho_n] = 1;
	COMPs_FType[id_rho_np1] = 1;
	COMPs_FType[id_rho_rez] = 1;

	COMPs_FType[id_Eta   ] = 1;
	COMPs_FType[id_rotB  ] = 3;
	COMPs_FType[id_ExB  ] = 3;
	COMPs_FType[id_gradPI0  ] = 3;
	COMPs_FType[id_divB  ] = 1;
	COMPs_FType[id_divE  ] = 1;
	COMPs_FType[id_divrotB] = 1;
	COMPs_FType[id_divU] = 1;
	COMPs_FType[id_gradPE] = 3;
	COMPs_FType[id_PhiDC ] = 1;
	COMPs_FType[id_PMagnetic] = 1;
	COMPs_FType[id_PTotal] = 1;

	COMPs_FType[id_FieldTime   ] = 1;
	COMPs_FType[id_ParticleTime] = 1;
	COMPs_FType[id_BxgradB   ] = 1;
	COMPs_FType[id_RespProc    ] = 1;
	
	COMPs_FType[id_scratch_scalar] = 1;
	COMPs_FType[id_scratch_vector] = 3;
        
        COMPs_FType[id_ElectronTemperature] = 1;
        COMPs_FType[id_rho_np1_recombined] = 1;
	
	COMPs_FType[id_Refine_Status] = 1;
	COMPs_FType[id_gradB] = 3;
	COMPs_FType[id_Refine_Rating] = 1;
	COMPs_FType[id_Lam    ] = 1;
	COMPs_FType[id_Gam    ] = 3;
	//! sequence of all ion densities
	COMPs_FType[id_allRhoSpecies  ] = num_Charged_Species;
	//! sequence of U, Gam and all ion densities
	COMPs_FType[id_UIp_Gam_allRhoSpec] =  3 +3 +num_Charged_Species;
	COMPs_FType[id_Bcfbg] = 3;
	COMPs_FType[id_Bcfbg_HIMesh] = 3;

	COMPs_FType[id_BTotal] = 3;
	COMPs_FType[id_BDerivative] = 3;

	//! handle to entire block field memory
	//! (including scratch)
	COMPs_FType[id_ALL_FIELDS] = NUM_SCALAR_FIELDS_EACH_BLK;


	for(INT32 species=0; species<num_Charged_Species; species++)
	{
		//! species density, velocity and temperature
		COMPs_FType[id_rhoSpecies1    +species] = 1;
		COMPs_FType[id_UI_Species1    +species] = 3;
                COMPs_FType[id_gyro_Species1  +species] = 1;
                COMPs_FType[id_gyro_el_Species1  +species] = 1;
                COMPs_FType[id_recomb_Species1  +species] = 1;
		COMPs_FType[id_ForceSpecies1  +species] = 3;
		COMPs_FType[id_PISpecies1     +species] = 3;
		COMPs_FType[id_PESpecies1     +species] = 1;
	}

	COMPs_FType[id_PEtotal] = 1;
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
	for(INT32 species=0; species<num_Neutral_Species; species++)
	{
		//! species density, velocity and temperature
		COMPs_FType[id_numberdensity_neutralSpecies1     + species] = 1;
		COMPs_FType[id_velocity_neutralSpecies1          + species] = 3;
		COMPs_FType[id_pressure_neutralSpecies1          + species] = 1;
		COMPs_FType[id_new_electron_beta_neutralSpecies1 + species] = 1;
	}
#endif

#if defined(use_ion_production_as_field)
	for(INT32 species=0; species<num_ion_prod_fields; species++)
	{
		//! species density, velocity and temperature
		COMPs_FType[id_density_ionProdSpecies1     + species] = 1;
		COMPs_FType[id_velocity_ionProdSpecies1    + species] = 3;
	}	
#endif

#if defined(use_dust_species_as_field)
	for(INT32 species=0; species<num_Dust_Species; species++)
	{
		//! species density, velocity and temperature
		COMPs_FType[id_density_dustSpecies1     + species] = 1;
		COMPs_FType[id_velocity_dustSpecies1    + species] = 3;
	}
#endif


	for(INT32 rho_extern=0; rho_extern<num_externRhoVelocityFields; rho_extern++)
	{
		//! extern species density and velocity
		COMPs_FType[id_externRho1 +rho_extern] = 1;
		COMPs_FType[id_extern_Ui1 +rho_extern] = 3;
	}

	COMPs_FType[id_extern_PhiC] = 1;
	
	num_scalar_average_fields = 0;

	//! average fields
	for(INT32 average_field=0; average_field<num_average_fields; average_field++)
	{
		COMPs_FType[id_average_Field1 +average_field] = COMPs_FType[IDs_of_Fields_to_average[average_field]];

		num_scalar_average_fields += COMPs_FType[IDs_of_Fields_to_average[average_field]];

	}





}

//!--------------------------------------------------------
//!- alloc_init_random_Generators: 
//!--------------------------------------------------------
void CHybrid::alloc_init_random_Generators(void)
{

	log_file << "Initializing random generators ..." << endl;
	unsigned long int seed;

	//! Initialize random generator INDEPENDENT of MPI proc
	//! -> equal numbers are generated on each process
	//! -> this is used to generate equal heavy ions
	//!    on each process
	//! -> last rand_gen (at number num_Charged_Species) is
	//!    for massload redistribution used in "CHybrid_MPI_Massload.cpp"
	randGen_of_species = new gsl_rng*[num_Charged_Species];


	for(INT32 species=0; species<num_Charged_Species; species++)
	{


		//!------------------------------------------------------------------------
		//!- GSL random Generators:
		//!------------------------------------------------------------------------
		//!This function returns a pointer to a newly-created instance of a random
		//!number generator of type gsl_rng_type
		randGen_of_species[species] = gsl_rng_alloc(gsl_rng_taus);


		
		//! Use:
		//! 1) Different seeds each Process for inflow   ions
		//! 2) Equal     seeds each Process for obstacle ions
		if(is_inflow_species[species])
		seed = 1000 *mpi_myRank  *(TL+1);
		else
		seed = 1000 *species;



		//!This function initializes (or `seeds') the random number generator.
		//!If the generator is seeded with the same value of 'seed' on two different runs,
		//!the same stream of random numbers will be generated by successive calls
		//!to the routines below. If different values of seed >= 1 are supplied, then
		//!the generated streams of random numbers should be completely different.
		//!If the seed s is zero then the standard seed from the original implementation
		//!is used instead.
		gsl_rng_set(randGen_of_species[species], seed);

	}


	//! init random generator for several functions that
	//! doesn't have to be synchronized across mpi processes
	seed = 123*mpi_myRank+1;

	randGen_general_asynchronous = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(randGen_general_asynchronous, seed);

	//! init random generator for several functions that
	//! has to be synchronized across mpi processes
	seed = 500;

	randGen_general_synchronized = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(randGen_general_synchronized, seed);
	

	log_file << "done." << endl;

}





//!--------------------------------------------------------
//!- alloc_RootBlks: 
//!--------------------------------------------------------
void CHybrid::alloc_RootBlks(void)
{

	log_file << "Allocating root block memory..." << endl;
	
	//! allocate rootblocks on every process
	Root_Block_Array = new CBlock[num_root_blocks];
	total_active_Blocks = num_root_blocks;

	//! Alloc memory for Block Lists
	BlockList_of_Lev = new CBlock*[MAX_LEVEL+1];
	memset(BlockList_of_Lev, 0,(MAX_LEVEL+1)*sizeof(CBlock*));

	GATHER_BlockList_of_Lev = new CBlock*[MAX_LEVEL+1];
	memset(GATHER_BlockList_of_Lev, 0,(MAX_LEVEL+1)*sizeof(CBlock*));

	//! Initialize Block lists
	BlockList_of_Lev[0] = Root_Block_Array;
	GATHER_BlockList_of_Lev[0] = Root_Block_Array;

	Root_Block_Array[num_root_blocks-1].next_Blk_of_BlockList = 0;
	Root_Block_Array[num_root_blocks-1].next_Blk_of_GatherBlockList = 0;


	log_file << "done." << endl;


}


//!--------------------------------------------------------
//!- init_RootBlocks
//! using default mpi process
//!--------------------------------------------------------
void CHybrid::init_RootBlocks(void)
{


	INT32* empty_array = NULL;
	init_RootBlocks_set_mpi_process(empty_array);



}


//!--------------------------------------------------------
//!- init_RootBlocks_SFC
//!  initialize position in box, memory and set fields
//!--------------------------------------------------------
void CHybrid::init_RootBlocks_set_mpi_process(INT32* root_mpi_process)
{

	if(use_SFC)
	log_file << "  Initializing Root Blocks along SFC..." << endl;
	else
	log_file << "  Initializing Root Blocks lineary..." << endl;


	if(root_mpi_process!=NULL)
	log_file << "   -> using restored MPI processes..." << endl;

	//! to store indices of space filling curve
	INT32 neighbour_BlkNr;
	INT32 blk_indices[3], neighbour_indices[3];
	INT32 RB[3] = {RB_X, RB_Y, RB_Z};

	//!--------------------------------------------------
	//! Apply SFC Indices to Root_Block_Array
	//!--------------------------------------------------
	for(INT32 blk_nr=0; blk_nr<num_root_blocks; blk_nr++)
	{

		//! get indices of root block
		//! -> always level 0
		if(use_SFC)
		SFC_BlkNr_to_Indices(blk_nr, blk_indices, 0);
		else
		LINEAR_BlkNr_to_Indices(blk_nr, blk_indices);


		//! initialize Block
		//! e.g. - set responsible_mpi_process
		Root_Block_Array[blk_nr].init_Block(blk_nr, blk_indices, 0, 0);

		//! set next_Blk_of_BlockList & next_Blk_of_GatherBlockList to next Array Element
		if(blk_nr<num_root_blocks-1)
		{
 			Root_Block_Array[blk_nr].next_Blk_of_BlockList = Root_Block_Array +(blk_nr +1);
 			Root_Block_Array[blk_nr].next_Blk_of_GatherBlockList = Root_Block_Array +(blk_nr +1);
 		}


		//! set MPI process if restored
		if(root_mpi_process!=NULL)
		Root_Block_Array[blk_nr].responsible_mpi_process = root_mpi_process[blk_nr];
// 		Root_Block_Array[blk_nr].responsible_mpi_process = 0;

		//! Only allocate memory and set fields in case myRank
		//! equals mpi_responsible_process
		if(mpi_myRank == Root_Block_Array[blk_nr].responsible_mpi_process)
		{
			Root_Block_Array[blk_nr].alloc_process_specific_Memory();
			Root_Block_Array[blk_nr].set_Fields();
		}

		//!-----------------------------------------------
		//! SET RESPECTIVE NEIGHBOURS
		//!-----------------------------------------------
		//!  Estimate _m1_ Neighbour
		for(INT32 direc=0; direc<3; direc++)
		 if(blk_indices[direc]>0)
		 {

			//! copy indices of this block to temporary buffer
			memcpy(neighbour_indices, blk_indices, 3*sizeof(INT32));

			//! decrease ind of respective direction
			neighbour_indices[direc]--;

			//! get index of root neighbour_BlkNr
			//! -> always level 0
			if(use_SFC)
			neighbour_BlkNr = SFC_Indices_to_BlkNr(neighbour_indices, 0);
			else
			neighbour_BlkNr = LINEAR_Indices_to_BlkNr(neighbour_indices);

			Root_Block_Array[blk_nr].Neighbour[2*direc] = Root_Block_Array +neighbour_BlkNr;
		 }

		//!  Estimate _p1_ Neighbour
		for(INT32 direc=0; direc<3; direc++)
		 if(blk_indices[direc] < RB[direc]-1 )
		 {
	
			//! copy indices of this block to temporary buffer
			memcpy(neighbour_indices, blk_indices, 3*sizeof(INT32));

			//! increase ind of respective direction
			neighbour_indices[direc]++;

			//! get index of root neighbour_BlkNr
			//! -> always level 0
			if(use_SFC)
			neighbour_BlkNr = SFC_Indices_to_BlkNr(neighbour_indices, 0);
			else
			neighbour_BlkNr = LINEAR_Indices_to_BlkNr(neighbour_indices);

			Root_Block_Array[blk_nr].Neighbour[2*direc +1] = Root_Block_Array +neighbour_BlkNr;
		 }

		//!-----------------------------------------------
		//! SET RESPECTIVE NEIGHBOURS in PERIODIC BOX
		//!-----------------------------------------------

		//! PERIODIC _m1_ Neighbour
		for(INT32 direc=0; direc<3; direc++)
		 if(use_periodic_bounds[direc] && blk_indices[direc]==0)
		 {

			//! copy indices of this block to temporary buffer
			memcpy(neighbour_indices, blk_indices, 3*sizeof(INT32));

			//! set index of respective direction to opposit Box Bondary
			neighbour_indices[direc] = RB[direc] -1;

			//! get index of root neighbour_BlkNr
			//! -> always level 0
			if(use_SFC)
			neighbour_BlkNr = SFC_Indices_to_BlkNr(neighbour_indices, 0);
			else
			neighbour_BlkNr = LINEAR_Indices_to_BlkNr(neighbour_indices);

			Root_Block_Array[blk_nr].Neighbour[2*direc] = Root_Block_Array +neighbour_BlkNr;
		 }


		//! PERIODIC _P1_ Neighbour
		for(INT32 direc=0; direc<3; direc++)
		if(use_periodic_bounds[direc] && blk_indices[direc]==RB[direc]-1)
		 {

			//! copy indices of this block to temporary buffer
			memcpy(neighbour_indices, blk_indices, 3*sizeof(INT32));

			//! set index of respective direction to opposit Box Bondary
			neighbour_indices[direc] = 0;

			//! get index of root neighbour_BlkNr
			//! -> always level 0
			if(use_SFC)
			neighbour_BlkNr = SFC_Indices_to_BlkNr(neighbour_indices, 0);
			else
			neighbour_BlkNr = LINEAR_Indices_to_BlkNr(neighbour_indices);

			Root_Block_Array[blk_nr].Neighbour[2*direc +1] = Root_Block_Array +neighbour_BlkNr;
		 }

		//! COPY ordanary neighbours to gather neighbours
		memcpy(Root_Block_Array[blk_nr].gather_Neighbour,
		       Root_Block_Array[blk_nr].Neighbour,
		       6*sizeof(CBlock*));
	}






	log_file << " done." << endl;


}


//!--------------------------------------------------------
//!- write_Info_Structure: Box Infos in are collected 
//!   in a structure wich is necessary for file output and
//!   eg. addres BlkNds_[i] (with varying "i" )
//!--------------------------------------------------------
void CHybrid::write_Info_Structure(void)
{


	if(mpi_myRank) return;
	
	SBox_Info Box_Info;
	
	Box_Info.max_level = MAX_LEVEL;
	Box_Info.R_obstacle = R_Obstacle;
	Box_Info.dt = dt;
	
	Box_Info.BoxOrigin[0] =  Box_Origin[0];
	Box_Info.BoxOrigin[1] =  Box_Origin[1];
	Box_Info.BoxOrigin[2] =  Box_Origin[2];
	
	Box_Info.BoxLength[0] =  LX;
	Box_Info.BoxLength[1] =  LY;
	Box_Info.BoxLength[2] =  LZ;
	
	Box_Info.RB_[0] = RB_X;
	Box_Info.RB_[1] = RB_Y;
	Box_Info.RB_[2] = RB_Z;
	
	Box_Info.BlkNds_[0] = BlkNds_X;
	Box_Info.BlkNds_[1] = BlkNds_Y;
	Box_Info.BlkNds_[2] = BlkNds_Z;
	
	Box_Info.Output2D = TL_OUTPUT_2D_NATIVE;
	Box_Info.Output3D = TL_OUTPUT_3D_NATIVE;
	
	char filename[200];
	sprintf(filename,"%s/native/Hybrid_Info_%s.info", data_output_path, Run_Name);
	
	ofstream HybridInfoFile;
	HybridInfoFile.open(filename, ofstream::binary);
	
	INT32 file_type = FILE_TYPE_INFO;
	INT32 st_length = RUN_NAME_SIZE;
	
	//! Szenario/Root Values
	HybridInfoFile.write(reinterpret_cast<char*> (&file_type),sizeof(INT32));
	HybridInfoFile.write(reinterpret_cast<char*> (&st_length),sizeof(INT32));
	HybridInfoFile.write( (Run_Name),RUN_NAME_SIZE*sizeof(char));
	
	//! this straight way maybe does not work when 
	//! visualisation and hybrid code run on different architectures 
	//! (eg. 32 vs. 64 Bit)
	HybridInfoFile.write( reinterpret_cast<char*> (&Box_Info),sizeof(SBox_Info));
	
	HybridInfoFile.flush();
	HybridInfoFile.close();



}






//!--------------------------------------------------------
//!- fill_Paticle_in_Box:
//!--------------------------------------------------------
void CHybrid::fill_Particle_in_Box(void)
{

	if(TestParticle_Simulation)
	return;

	log_file << endl;
	if(particle_restored)
	{
		log_file << " Particles initialization skipped." << endl << endl;
		return;
	}
	
	log_file << " Initializing particles ... " << endl;


	INT32 fill_counter = 0;
	num_total_particles = 0;
	num_injected_particles = 0;

	clock_t start,finish;
	D_REAL fill_time;
	start = clock();
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		log_file  << " -> Filling Level " << level
		          << " Blocks ... ("<<total_active_Blocks<<" in total)" << endl;

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank==temp_Block->responsible_mpi_process)
			 for(INT32 oct=0; oct<8; oct++)
			 {
	
				if(!temp_Block->child_array[oct])
				temp_Block->fill_Oct(oct);
			 }


			fill_counter++;
			if(fill_counter%100==0)
			log_file << "    [" << fill_counter << "] top level Blocks filled" << endl;

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}



	INT64 local_info_values[2] = {num_injected_particles,
				      num_total_particles};


	stringstream info_names[2];

	info_names[0] << " -> num_injected_particles:  ";
	info_names[1] << " -> num_total_particles: ";

	show_information(local_info_values,
			 info_names,
			 2, BUILD_SUM);
	
	

	finish = clock();
	fill_time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " fill time:     " << fill_time << "s." << endl << endl;



}


//!--------------------------------------------------------
//!- do_consistency_check_init_derived_variables
//!--------------------------------------------------------
void CHybrid::do_consistency_check_init_derived_variables(void)
{



	local_mark_particle_counter = 100000000 *mpi_myRank;

	//! syncronize here, else error messages are not displayed at every processors logfile
	synchronize_allProcesses_MPI();

	//! check whether size of info Array is large enuogh
	if(7 +2*num_Charged_Species +(MAX_LEVEL+1) > INFO_ARRAY_SIZE)
	{
		log_file << endl << endl << endl;
		log_file << " INFO_ARRAY_SIZE to small: " << INFO_ARRAY_SIZE<< endl;
		log_file << " please increase INFO_ARRAY_SIZE to at least "
			 << 7 +2*num_Charged_Species +(MAX_LEVEL+1) << " ." <<endl;
		log_file << " see defines.h for further details. " <<endl;
		log_file << " Exiting ... " << endl;
		exit(1);


	}

	//! check whether num Blk nodes are even
	if(BlkNds_X%2 != 0 || BlkNds_Y%2 != 0 || BlkNds_Z%2 != 0)
	{
		log_file << endl << endl << endl;
		log_file << " BlkNds must be EVEN: " << endl;
		log_file << " BlkNds_X: " << BlkNds_X << endl;
		log_file << " BlkNds_Y: " << BlkNds_Y << endl;
		log_file << " BlkNds_Z: " << BlkNds_Z << endl;
		log_file << " Exiting ... " << endl;
		exit(1);

	}

	log_file << endl;
	log_file << "   Fields in use: " << id_average_Field1 +num_average_fields-1 << endl;
	log_file << endl;

	//! has to be larger (equal is not sufficient, since index=[0,NUM_FIELDS-1])
	if( id_average_Field1 +(num_average_fields-1) >= NUM_FIELDS)
	{

		//! NOTE:
		//! Even though precalcRho Fiels is not use, id_externRho1 will be used to 
		//! access  COMPs_FType, Field_Type, etc ...
		log_file << endl << endl << endl;
		log_file << "   id_average_Field1"		<< endl
		     << "  +num_average_fields" << endl
		     << "  exceeds/equals NUM_FIELDS: " << endl
		     << "  " << id_average_Field1 + (num_average_fields-1)
		     << "  >=" << NUM_FIELDS << endl << endl
		     //! plus 1 below since index=[0,NUM_FIELDS-1]
		     << " increase NUM_FIELDS in defines.h to at least "
		     << id_extern_Ui1 +(num_externRhoVelocityFields-1) +1 << endl
		     << endl
		     << " num_Charged_Species: " << num_Charged_Species << endl
		     << " id_rhoSpecies1: " << id_rhoSpecies1 << endl
		     << " id_UI_Species1: " << id_UI_Species1 << endl
		     << " id_PISpecies1: " << id_PISpecies1 << endl
		     << " id_PESpecies1: " << id_PESpecies1 << endl << endl
		     << " num_externRhoVelocityFields: " << num_externRhoVelocityFields << endl
		     << " id_externRho1: " << id_externRho1 << endl
		     << " id_extern_Ui1: " << id_extern_Ui1 << endl
		     << " id_average_Field1: " << id_average_Field1 << endl
		     << " num_average_fields: " << num_average_fields << endl;
		exit(1);
	}



	if(mpi_num_processes > 1)
	{


		if(TL_OUTPUT_2D_NATIVE || TL_OUTPUT_2D_GNUPLOT)
		{

			log_file << "ERROR:" << endl;
			log_file << "Only Silo2D and Silo3D supported for parallel runs" << endl;
			log_file << "Please adjust output parameter in parameter.cpp." << endl;
			log_file << "Exiting ..." << endl;
			exit(1);
		}
	}

	INT32 fields_in_use = id_UI_Species1 +(num_Charged_Species-1) +num_externRhoVelocityFields;

	INT32 field_memory = fields_in_use *num_nodes_in_block *sizeof(D_REAL);

	//! do not calc first UI twice ->(num_Charged_Species-1)
	log_file << "Fields in use: " << fields_in_use << endl;
	log_file << "resulting in " << field_memory/1000 <<"kB field memory per Block." <<  endl;


	log_file << endl << endl;
	log_file << "This System uses the denoted number of bytes " << endl
	         << "for the respective data type:" << endl;
	log_file << "sizeof(INT32):      " << sizeof(INT32) << endl;
	log_file << "sizeof(INT32):      " << sizeof(INT32) << endl;
	log_file << "sizeof(INT32):      " << sizeof(INT32) << endl;
	
	log_file << "sizeof(D_REAL):     " << sizeof(D_REAL) << endl;
	log_file << "sizeof(F_REAL):     " << sizeof(F_REAL) << endl;
	log_file << "sizeof(u long long):     " << sizeof(unsigned long long) << endl;
	
	log_file << "sizeof(CHybrid): " << sizeof(CHybrid) << endl;

	log_file << "sizeof(CBlock):    "  << sizeof(CBlock)  << endl;
	log_file << "sizeof(CBlock*):   "  << sizeof(CBlock*) << endl;

	stringstream info_names[INFO_ARRAY_SIZE];
	log_file << "sizeof(info_names[INFO_ARRAY_SIZE]):   " 
	<< sizeof(info_names) << endl;

	log_file << "sizeof(particle):  " << sizeof(particle) << endl;
	log_file << "sizeof(particle*): " << sizeof(particle*) << endl;

	
	log_file << "sizeof(SBox_Info): " << sizeof(SBox_Info) << endl;

	log_file << endl << "Run_Name: " << Run_Name << "." << endl << endl;

#ifdef USE_CFBG_BFIELD
	if( !vec_len(Magnetic_Moment))
	{
		log_file << endl << endl << endl;
		log_file << "   no magnetic Moment set." << endl;
		log_file << "   Either set Magnetic_Moment != 0 or " << endl;
		log_file << "   remove 'USE_CFBG_BFIELD' statement in " << endl;
		log_file << "   file 'defines.h'. " << endl;
		log_file << "   EXITING ... " << endl;
		exit(1);
		

	}
#endif

#ifndef USE_CFBG_BFIELD
	if(vec_len(Magnetic_Moment))
	{
		log_file << endl << endl << endl;
		log_file << "   Magnetic Moment set but required field method is turned off." << endl;
		log_file << "   Either set Magnetic_Moment == 0 or " << endl;
		log_file << "   incomment 'USE_CFBG_BFIELD' statement in " << endl;
		log_file << "   file 'defines.h'. " << endl;
		log_file << "   EXITING ... " << endl;
		exit(1);
		

	}
#endif

	//! Stand off distance is
	//! x = pow((Magneti_Moment/(v_sw*sqrt(2.))),(1./3.))
	//! no, i'm not sure it is. this may apply for mercury, but not for Callisto --LL
	if(vec_len(Magnetic_Moment))
	  {
	    log_file << endl << endl;
	    log_file << "-------------------------Dipole Parameters-------------------------" << endl << endl;
	    log_file << "   Normalized magnetic moment magnitude: " << vec_len(Magnetic_Moment) << endl;
	    //! NOTE: think this stand off distance is only good for mercury? commented out. see below --LL
	    //log_file << "   Stand off distance: " << pow((vec_len(Magnetic_Moment)/(vec_len(V_sw)*sqrt(2.))),(1./3.)) << endl;
	    log_file << "   Dipole field at equator = "<<vec_len(MM)/(R_Moon*R_Moon*R_Moon*SI_x0*SI_x0*SI_x0*1e7*SI_B0)<<" B0" << endl;
	    log_file << "   Dipole field at poles   = "<<2*vec_len(MM)/(R_Moon*R_Moon*R_Moon*SI_x0*SI_x0*SI_x0*1e7*SI_B0)<<" B0" << endl << endl;

	    log_file << "   Upstream ram pressure           = " << ram << " Pa" << endl;
	    log_file << "   Upstream thermal pressure       = " << thermal << " Pa" << endl;
	    log_file << "   Upstream mag pressure           = " << mag << " Pa" << endl;
	    
	    log_file << "   Standoff distance from Callisto = " << pow(mu_0 * pow(vec_len(MM),2)/(8*pi*(ram + thermal + mag)),1./6.)/1560800 << " Rc" << endl << endl;
	    log_file << "-------------------------------------------------------------------" << endl << endl;
	  }


	//! check whether equal distributuin of blocks among np is possible
// 	if(num_root_blocks%mpi_num_processes)
// 	{
// 
// 		log_file << endl;
// 		log_file << "ERROR:" << endl;
// 		log_file << "num_root_blocks%mpi_num_processes: " << num_root_blocks%mpi_num_processes << endl;
// 		log_file << "please improve code." << endl;
// 		exit(0);
// 	}

#ifdef DUST_FIELD

	if(num_externRhoVelocityFields<2)
	{
		log_file << " DUST_FIELD is set in defines.h but num_externRhoVelocityFields is < 2 ... " <<endl;
		log_file << " EXITING ... " << endl;
		finalize_MPI();
		exit(1);
	}
#endif

#ifdef DUST_PARTICLE

	if(num_externRhoVelocityFields<2)
	{
		log_file << " DUST_PARTICLE is set in defines.h but num_externRhoVelocityFields is < 2 ... " <<endl;
		log_file << " EXITING ... " << endl;
		finalize_MPI();
		exit(1);
	}
#endif

#ifdef DUST_FIELD
	
	#ifdef DUST_PARTICLE
		log_file << " DUST_FIELD and DUST_PARTICLE are set in defines.h ... " <<endl;
		log_file << " EXITING ... " << endl;
		finalize_MPI();
		exit(1);
	#endif
#endif


#ifdef DUST_FIELD
	if(num_Charged_Species==num_Particle_Species)
	{
		log_file << " num_Charged_Species = num_Particle_Species but DUST_FIELD is set in defines.h ... " <<endl;
		log_file << " EXITING ... " << endl;
		finalize_MPI();
		exit(1);
	}	
#endif

#ifdef DUST_PARTICLE
	if(num_Inflow_Species==num_Particle_Species)
	{
		log_file << " num_Inflow_Species = num_Particle_Species but DUST_PARTICLE is set in defines.h ... " <<endl;
		log_file << " EXITING ... " << endl;
		finalize_MPI();
		exit(1);
	}
#endif


#ifndef TRACK_PARTICLE
	for(INT32 species=0; species<num_Charged_Species; species++)
	{


		if(TL_INJECT_PARTICLE_TO_TRACK >=0 || TL_OUTPUT_PARTICLE_TRACKS > 0)
		{
			log_file << "ERROR: "  << endl;
			log_file << " but TRACK_PARTICLE is not defined in defines.h" << endl;
	
			log_file << "EXITING ... "  << endl;
			finalize_MPI();
			exit(1);
		}
	}
#endif

#ifdef TRACK_PARTICLE

// 	INT32 num_mark_total_particle = 0;
// 	for(INT32 species=0; species<num_Charged_Species; species++)
// 	num_mark_total_particle += num_mark_particle_in_species[species];


/*
	if(num_mark_total_particle<=0)
	{
		log_file << "WARNING: "  << endl;
		log_file << "  num_mark_total_particle = "<< num_mark_total_particle << endl
			 << " but TRACK_PARTICLE is defined in defines.h" << endl;

// 		log_file << "EXITING ... "  << endl;
// 		finalize_MPI();
// 		exit(1);
	}*/

#endif


	if( sizeof(particle)%8!=0 )
	{

		log_file << "WARNING: "  << endl;
		log_file << "sizeof(particle) =  " << sizeof(particle) << " is not multiple of 8" << endl;
		log_file << "!!!--------------------------------------------!!!" << endl;
		log_file << "!!!     THIS WILL RESULT IN SERIOUS ERRORS     !!!" << endl;
		log_file << "!!! INCREASE sizeof(particle) TO MULTIPLE OF 8 !!!" << endl;
		log_file << "!!!   (see CBlk.h for declaratin of particle)  !!!" << endl;
		log_file << "!!!--------------------------------------------!!!" << endl;
		log_file << "Exiting ... " << endl;
		exit(0);
	
	}

	if(num_Inflow_Species>num_Charged_Species)
	{
		log_file << "Error:" << endl
			 << "num_Charged_Species:     " << num_Charged_Species << endl
			 << "num_Inflow_Species: " << num_Inflow_Species << endl
			 << "To few ionSpecies defined." << endl;
		log_file << "Exiting ... " << endl;
		exit(0);
	}

	if(num_Neutral_Species>NUM_NEUTRAL_SPECIES)
	{
		log_file << "Error:" << endl
			<< "num_Neutral_Species (parameters.cpp):     " << num_Neutral_Species << endl
			<< "NUM_NEUTRAL_SPECIES (defines.h): " << NUM_NEUTRAL_SPECIES << endl
			<< "Increase NUM_NEUTRAL_SPECIES!" << endl;
		log_file << "Exiting ... " << endl;
		finalize_MPI();
		exit(1);
	}	
	
	if(num_Particle_Species>NUM_PARTICLE_SPECIES)
	{
		log_file << "Error:" << endl
			<< "num_Particle_Species (parameters.cpp):     " << num_Particle_Species << endl
			<< "NUM_PARTICLE_SPECIES (defines.h): " << NUM_PARTICLE_SPECIES << endl
			<< "Increase NUM_PARTICLE_SPECIES!" << endl;
		log_file << "Exiting ... " << endl;
		finalize_MPI();
		exit(1);
	}
	
	//! give information on simulation
	log_file << endl<<" dt = "<<dt<<endl;	
	log_file <<  " SI_v0 = "<<SI_v0*1.e-03<<" km/s"<<endl;
	log_file <<  " SI_x0 = "<<SI_x0*1.e-03<<" km"<<endl;
	log_file <<  " SI_t0 = "<<SI_t0<<" s"<<endl;
	log_file <<  " SI_B0 = "<<SI_B0*1.e+09<<" nT"<<endl<<endl;

	log_file <<  " M_A   = "<<MA<<endl;
	log_file <<  " M_MS  = "<<MA/sqrt((Ion_Betas[0]+Electron_Betas[0]+1)*kappa_electron/2)<<endl;
	log_file <<  " M_S   = "<<SW_v/sqrt(kappa_electron*kB*(Ion_Betas[0]+Electron_Betas[0])/(calcBeta*SI_m0))<<endl<<endl;
	//! check for sub magnetosonic
	if(MA/sqrt((Ion_Betas[0]+Electron_Betas[0]+1)*kappa_electron/2) >= 1)
	  {
	    log_file << "  !!!!!!!!" << endl
		     << "  WARNING:" << endl
		     << "  !!!!!!!!" << endl
		     << "  M_MS is larger than 1! " << endl
		     << "  This may not be intended" << endl
		     << "  Try increasing Beta_i or reducing SW_v."  << endl << endl;
	  }
	log_file <<  " calcBeta = "<<calcBeta<<endl;
	log_file << " Ion temperature is: "<<kB*Ion_Betas[0]/(e*calcBeta)<<" eV"<<endl<<endl;

	//! special che in case SFC is used
	if(use_SFC)
	{

		//! check whether equal number R_B nodes are
		if(RB_X!=RB_Y || RB_Y!=RB_Z )
		{
			log_file << " RB_X: " << RB_X << endl;
			log_file << " RB_Y: " << RB_Y << endl;
			log_file << " RB_Z: " << RB_Z << endl;
	
			log_file << " Equal number of root blocks \
				each direction required to build \
				hilbert SFC! " << endl;
	
			log_file << " Exiting ..." << endl;
			exit(1);
		}
	
		//! estimate power
		D_REAL power = log(1. *RB_X) / log(2.);
	
		//! - check whether root blocks are power of two
		//! - add small number to avoid round errors
		//!   (case RB_X=8  power may equal 2.999999999)
		D_REAL offset = 1.e-2;
		if( power < 1.*int(power+offset)-offset  ||  power > 1.*int(power+offset) +offset )
		{
			log_file << " RB_X no power of 2: " << power << endl;
			log_file << "         int(power): " << int(power+1.e-2) << endl;
			log_file << " Exiting ..."          << endl;
			exit(1);
		}

		SFC_RB_power = int(power+offset);

		INT32 SFC_ML_power = SFC_RB_power +MAX_LEVEL;

		log_file << endl << endl;
		log_file << "SFC in use:" << endl;
		log_file << "MAX_LEVEL is set to: " << MAX_LEVEL << endl;
		log_file << "-> number of root blocks equates 2^" << SFC_RB_power << endl;
		log_file << "-> number of max level blocks equates 2^" << SFC_ML_power << endl;

		log_file << "-> SFC bimask requires nDim*SFC_ML_power bits = " << 3*SFC_ML_power << endl;
		log_file << "-> sizeof active bitmask (bits): " << 8*sizeof(bitmask_t) << endl;

		if(8*sizeof(bitmask_t) < 3*SFC_ML_power)
		{

			log_file << " ERROR:" << endl;
			log_file << " Too few bits in bitmask !!!" << endl;
			log_file << " Use integer of higher precision:" << endl;
			log_file << " file:    'hilbert.h'." << endl;
			log_file << " typedef: bitmask_t " << endl;
			log_file << " Exiting ..." << endl;
			exit(1);
		}
		
		
		log_file << endl << endl;

	}

	//! ------------- Set/derive thermal velocities ----------------------------
	PARTICLE_REAL norm_value, v2mean;
	
	thermal_velocity_para_of = new PARTICLE_REAL[num_Charged_Species];
	memset(thermal_velocity_para_of, 0, num_Charged_Species*sizeof(PARTICLE_REAL));
	thermal_velocity_perp_of = new PARTICLE_REAL[num_Charged_Species];
	memset(thermal_velocity_perp_of, 0, num_Charged_Species*sizeof(PARTICLE_REAL));
	
  	is_inflow_species = new bool[num_Charged_Species];
	memset(is_inflow_species, 0, num_Charged_Species*sizeof(bool));

	for(INT32 species=0; species<num_Charged_Species; species++)
	 for(INT32 inflow_sp=0; inflow_sp<num_Inflow_Species; inflow_sp++)
	  if(index_Inflow_Species[inflow_sp]==species)
	  is_inflow_species[species] = true;

	for(INT32 species=0; species<num_Charged_Species; species++)
	  if(!Ion_Charges[species] || !Ion_Masses[species])
	  {
		  log_file<<endl<<endl;
		  log_file << " ERROR:"<<endl
			   << " Ion_Charges and Ion_Masses have to be set for all species!"<<endl
			   << " (num_Charged_Species = "<<num_Charged_Species<<" )"<<endl 
		  	   << "Exiting ..." <<endl;
		  exit(1);
	  }
	 
	 
	for(INT32 extern_rho=0; extern_rho< num_externRhoVelocityFields; extern_rho++)
	 if(index_externRho_Species[extern_rho] >= num_Charged_Species)
	 {
	    log_file << endl << endl;
	    log_file << " ERROR:" << endl
		     << " index_externRho_Species["<<extern_rho<<"] = " 
		     << index_externRho_Species[extern_rho] << endl
		     << " maximal allowed value is num_Charged_Species-1 = "<< num_Charged_Species-1 <<" ." << endl
		     << " Exiting ..." << endl << endl;
		exit(0);
	 }

	for(INT32  species=0; species < num_Charged_Species; species++)
	{

		if(Ion_Betas[species] !=0.)
		{
			if(Ion_Masses[species]==0.)
			{
			    log_file << " ERROR:" << endl
				     << " Plasma Beta of species " << species<< " may not used, "   << endl
				     << " as no 'Ion_Masses' is speciefied for respective species." << endl
				     << " Use 'Ion_Thermal_Velocities' instead." << endl << endl;
				exit(0);
			}

			if(rho_sw[species]==0.)
			{
			    log_file << " ERROR:" << endl
				     << " Plasma Beta of species " << species<< " may not be used, " << endl
				     << " as no 'rho_sw' is speciefied for respective species. " << endl
				     << " Use 'Ion_Thermal_Velocities["<<species<<"]'"
				     << " instead to specify temperature." << endl << endl;
				exit(0);
			}

			//! see S.Simon's PhD Thesis pp. 261
			norm_value = 1. /(Ion_Masses[species]*rho_sw[species]);
				
		
			v2mean = 1.5 * Ion_Betas[species] *norm_value;
			thermal_velocity_para_of[species] = sqrt(v2mean)/sqrt(3.);
			thermal_velocity_perp_of[species] = sqrt(v2mean)/sqrt(3.);

			//! division by sqrt(3) is necessary since thermal_velocity_para_of 
			//! is used as the standard deviation sigma of a gaussian
			//! and this is given by sqrt(kT/m) (without 3)
		}
		else
		{
			thermal_velocity_para_of[species] = sqrt(Ti_para[species]/(Ion_Masses[species]*SI_m0))/SI_v0;
			thermal_velocity_perp_of[species] = sqrt(Ti_perp[species]/(Ion_Masses[species]*SI_m0))/SI_v0;
		}	
			

		log_file << " Thermal Velocity (parallel to B0) of species" << species<<": "
		         << thermal_velocity_para_of[species]<< " = " << thermal_velocity_para_of[species]*SI_v0*1.e-3<<" km/s"<<endl;
		log_file << " Thermal Velocity (perpendicular to B0) of species" << species<<": "
		         << thermal_velocity_perp_of[species]<< " = " << thermal_velocity_perp_of[species]*SI_v0*1.e-3<<" km/s"<<endl<<endl;
	}



	//! ------------- Choose BField_Method ----------------------------
	if(advance_B_Algorithm!=0 && advance_B_Algorithm!=1)
	{
	    log_file << "Error:" << endl
		     << "advance_B_Algorithm:     " << advance_B_Algorithm << endl
		     << "Set advance_B_Algorithm to either 0 or 1." << endl;
		exit(0);
	}





#ifdef ETA_TERM_BField
	if(advance_obstacle_B)
	{
		log_file << "Inconsistency:" <<endl
		<< "Both BField Methods (LF and CN) are switched on"<< endl
		<< "Switch of at least one of them." <<endl
		<< "CN in parameter.cpp" <<endl
		<< "LF in defines.h" <<endl
		<< "Exiting ..." <<endl;
		exit(1);
	}

	if(!use_resistive_obstacle && Eta_sw==0.)
	{
		log_file << endl;
		log_file << "WARNING:" <<endl
		<< "No Resistivity specified but resistivity Term activated!"<< endl
		<< "Switch off:" <<endl
		<< "LF in defines.h" <<endl;
		log_file << endl;
// 		<< "Exiting ..." <<endl;
// 		exit(1);
	}

#endif


#ifndef ETA_TERM_BField
	if(!advance_obstacle_B && (use_resistive_obstacle || Eta_sw!=0.))
	{
		log_file << "Inconsistency:" <<endl
		<< "Resistivity specified but resistivity Term deactivated!"<< endl
		<< "Switch on at least one of either:" <<endl
		<< "CN in parameter.cpp" <<endl
		<< "LF in defines.h" <<endl
		<< "Exiting ..." <<endl;
		exit(1);
	}
#endif

	if(advance_obstacle_B && !use_resistive_obstacle && Eta_sw==0.)
	{
		log_file << endl;
		log_file << "WARNING:" <<endl
		<< "No Resistivity specified but advance_obstacle_B activated!"<< endl
		<< "Switch off:" <<endl
		<< "CN in parameter.cpp" <<endl;
		log_file << endl;

	}
	log_file<<"paramter checks"<<endl;


}

//!--------------------------------------------------------
//!- initialize some global variables at simulation begin
//!--------------------------------------------------------
void CHybrid::init_global_variables(void)
{
	log_file<<"before ion neutral init"<<endl;
	init_ion_neutral_variables();
	log_file<<"after ion neutral init"<<endl;
	negative_particles_in_Simu = false;
	
	//! check if negative particles in simu
	 for(INT32 species=0; species<num_Particle_Species; species++)
	  if(Ion_Charges[species]<0)
	    negative_particles_in_Simu = true; 

	
}	



//!--------------------------------------------------------
//!- alloc_variables
//!--------------------------------------------------------
void CHybrid::alloc_variables(void)
{


	//!***************************************//
	//!***** COPY GN RELATED *****************//
	//!***************************************//

	log_file << "Setting function pointers..." << endl;
	//! set some utils.h function pointer:
	send_GN_MPI[0] = &cp_I_GN_send_MPI;
	send_GN_MPI[1] = &cp_J_GN_send_MPI;
	send_GN_MPI[2] = &cp_K_GN_send_MPI;

	//! set some utils.h function pointer:
	recv_GN_MPI[0] = &cp_I_GN_recv_MPI;
	recv_GN_MPI[1] = &cp_J_GN_recv_MPI;
	recv_GN_MPI[2] = &cp_K_GN_recv_MPI;

	//! set some utils.h function pointer:
	get_GN_equal_process[0] = &get_I_GN_equal_process;
	get_GN_equal_process[1] = &get_J_GN_equal_process;
	get_GN_equal_process[2] = &get_K_GN_equal_process;

	//!***************************************//
	//!***** ADD GN RELATED ******************//
	//!***************************************//

	//! set some utils.h function pointer:

	send_add_GN_MPI[0] = &add_I_GN_send_MPI;
	send_add_GN_MPI[1] = &add_J_GN_send_MPI;
	send_add_GN_MPI[2] = &add_K_GN_send_MPI;

	//! set some utils.h function pointer:
	recv_add_GN_MPI[0] = &add_I_GN_receive_MPI;
	recv_add_GN_MPI[1] = &add_J_GN_receive_MPI;
	recv_add_GN_MPI[2] = &add_K_GN_receive_MPI;

	//! set some utils.h function pointer:
	get_add_GN_equal_process[0] = &add_I_GN_same_process;
	get_add_GN_equal_process[1] = &add_J_GN_same_process;
	get_add_GN_equal_process[2] = &add_K_GN_same_process;


	//! set some utils.h function pointer:
	BV_from_parent[0] = &im1_BV_from_parent;
	BV_from_parent[1] = &ip1_BV_from_parent;
	BV_from_parent[2] = &jm1_BV_from_parent;
	BV_from_parent[3] = &jp1_BV_from_parent;
	BV_from_parent[4] = &km1_BV_from_parent;
	BV_from_parent[5] = &kp1_BV_from_parent;



	log_file << "Allocating variables..." << endl;
	//! array to indicate which field (even vs. odd)
	//! has to be updated
// 	GN_type_after = new INT32[NUM_SUB_CYCLE+3];


	for(int cycle=0; cycle<NUM_SUB_CYCLE+2; cycle+=2)
	{
		GN_type_after[cycle] = id_BEven;
		GN_type_after[cycle+1] = id_BOdd;
	}
	log_file<<"this was ok1"<<endl;
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
// 	GN_type_Pe_after = new INT32[NUM_SUB_CYCLE+3];

	log_file<<"this was ok2.0"<<endl;
	for(int cycle=0; cycle<NUM_SUB_CYCLE+2; cycle+=2)
	{
	  		log_file<<"cycle test "<<cycle<<endl;
		GN_type_Pe_after[cycle] = id_PE_even;
		GN_type_Pe_after[cycle+1] = id_PE_odd;
	}
	log_file<<"this was ok2"<<endl;
#endif

	log_file<<"this was ok3"<<endl;
	num_particle_total_storage = 0;

  	level_time = new INT32[MAX_LEVEL+1];
   	memset(level_time,0,(MAX_LEVEL+1)*sizeof(INT32));

   	step_time = new INT32[MAX_LEVEL+1];
   	memset(step_time,0,(MAX_LEVEL+1)*sizeof(INT32));

   	CYCLE_of_L = new INT32[MAX_LEVEL+1];
   	memset(CYCLE_of_L,0,(MAX_LEVEL+1)*sizeof(INT32));

	//! --------- Refinement --------------
	global_min_refValue = new F_REAL[MAX_LEVEL+1];
   	memset(global_min_refValue,0,(MAX_LEVEL+1)*sizeof(F_REAL));

	global_max_refValue = new F_REAL[MAX_LEVEL+1];
   	memset(global_max_refValue,0,(MAX_LEVEL+1)*sizeof(F_REAL));

	total_Blocks_L = new INT64[MAX_LEVEL+1];
   	memset(total_Blocks_L,0,(MAX_LEVEL+1)*sizeof(INT64));

	total_Blocks_L[0] = num_root_blocks;
	log_file<<"this was ok4"<<endl;

	num_total_particles_in_L = new INT64[MAX_LEVEL+1];
  	memset(num_total_particles_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));

	num_split_in_L = new INT64[MAX_LEVEL+1];
  	memset(num_split_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));

	num_merge_in_L = new INT64[MAX_LEVEL+1];
  	memset(num_merge_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));

	num_out_of_Box = new INT64[num_Charged_Species+1];
  	memset(num_out_of_Box,0,(num_Charged_Species+1)*sizeof(INT64));

	num_in_obst = new INT64[num_Charged_Species+1];
  	memset(num_in_obst,0,(num_Charged_Species+1)*sizeof(INT64));

	move_time = new INT64[MAX_LEVEL+1];
  	memset(move_time,0,(MAX_LEVEL+1)*sizeof(INT64));



	log_file<<"this was ok5"<<endl;
   	for(INT32 level=0; level<=MAX_LEVEL; level++)
   	{

		step_time[level] = 1;
		INT32 potenz = (MAX_LEVEL-level);

		for(INT32 a=0; a<potenz;a++)
		step_time[level]*=2;

		//! multiply with to in order to allow
		//! half time steps in Leap Frog Method
		step_time[level]*=2;

// 		log_file << "step_time["<< level<<"]: " << step_time[level] << endl;

// 		if(advance_B_Algorithm==0)
// 		step_time[level]*=7;

   	}
	log_file<<"this was ok6"<<endl;
	//! init general Block variables (these are no members)
	init_Block_Globals();
	log_file<<"this was ok7"<<endl;
	//! 2D-Cuts
	num_Blocks_transfered_L = new INT64[MAX_LEVEL+2];
	num_particle_processed_L = new INT64[MAX_LEVEL+2];

	cell_indedx_L = new INT32[MAX_LEVEL+1];
	XBlock_CrossSection = new FILE_REAL[3*BlkNds_Y*BlkNds_Z];
	YBlock_CrossSection = new FILE_REAL[3*BlkNds_X*BlkNds_Z];
	ZBlock_CrossSection = new FILE_REAL[3*BlkNds_X*BlkNds_Y];


	//! reset initial values
	initial_particle_mass   = 0.;
	initial_particle_energy = 0.;
	initial_particle_number = 0.;
	memset(initial_particle_momentum, 0, 3*sizeof(D_REAL));

	num_split_since_start = 0;
	num_merge_since_start = 0;


	TL_at_last_restore_state = 0;
	variation_time_mpi  = 0;
	variation_time_comp = 0;

#ifdef startTL_PROTOCOL_OBSTACLE_PARTICLE
	sprintf(filename,"%s/ObstacleParticle_prc%d.txt",data_output_path, mpi_myRank);
	ObstacleParticle_FILE.open(filename);
#endif

	log_file << "done." << endl;
}

//!--------------------------------------------------------
//!- build_statistic
//!--------------------------------------------------------
void CHybrid::measure_time(void)
{

	log_file << endl;
  	log_file <<  "/******************************************/" << endl;
  	log_file <<  "/*           Time Level: " << TL <<"              */" << endl;
  	log_file <<  "/******************************************/" << endl;

	loop_start_comp = clock();
	time (&loop_start_phys);
	loop_start_mpi =  MPI_Wtime();


}

//!--------------------------------------------------------
//!- build_statistic
//!--------------------------------------------------------
void CHybrid::statistic(void)
{


// 	time_t active_time_phys;

// 	double loop_time_phys, run_time_phys;
	double loop_time_mpi,  run_time_mpi,  active_time_mpi, average_time_mpi;
	double loop_time_comp, run_time_comp, active_time_comp, average_time_comp;

	INT32 TL_processed = (TL -TL_at_last_restore_state);

	log_file << " STATISTICS OF TL "<<TL<<":"  << endl;


	//! measure mpi time
	active_time_mpi = MPI_Wtime();
	run_time_mpi  = active_time_mpi -run_start_mpi;
	loop_time_mpi = active_time_mpi -loop_start_mpi;


	//! measure comp time
	active_time_comp = clock();
	run_time_comp = (double(active_time_comp)-double(run_start_comp))/CLOCKS_PER_SEC;
	loop_time_comp = (double(active_time_comp)-double(loop_start_comp))/CLOCKS_PER_SEC;

	//! measure physical time
// 	time (&active_time_phys);
// 	run_time_phys 	= difftime (active_time_phys, run_start_phys);
// 	loop_time_phys  = difftime (active_time_phys, loop_start_phys);


	average_time_mpi  = run_time_mpi / TL_processed;
	average_time_comp = run_time_comp / TL_processed;

	variation_time_mpi  += (loop_time_mpi -average_time_mpi)
			      *(loop_time_mpi -average_time_mpi)/ TL_processed;

	variation_time_comp += (loop_time_comp -average_time_comp)
			      *(loop_time_comp -average_time_comp)/ TL_processed;

	log_file << " ---------------------------------------" << endl;
	log_file << " - Run time  (mpi   time): " << run_time_mpi << "s." << endl;
	log_file << " - Run time  (comp. time): " << run_time_comp << "s." << endl;
// 	log_file << " - Run time  (phys. time): " << run_time_phys << "s." << endl;
	log_file << " ---------------------------------------" << endl;
	log_file << " - Loop time  (mpi   time): " << loop_time_mpi  << "s." << endl;
	log_file << " - Loop time  (comp. time): " << loop_time_comp << "s." << endl;
// 	log_file << " - Loop time  (phys. time): " << loop_time_phys << "s." << endl;
	log_file << " ---------------------------------------" << endl;
	log_file << "   -> Avrg. s/TL (mpi   time): " << average_time_mpi  << "s." << endl;
	log_file << "   -> Avrg. s/TL (comp. time): " << average_time_comp << "s." << endl;
// 	log_file << "   -> Avrg. s/TL (phys. time): " << run_time_phys / TL_processed << "s." << endl;
	log_file << " ---------------------------------------" << endl;
	log_file << "   -> Variation s/TL (mpi   time): " << sqrt(variation_time_mpi) << "s." << endl;
	log_file << "   -> Variation s/TL (comp. time): " << sqrt(variation_time_comp) << "s." << endl;
	log_file << " ---------------------------------------" << endl;




	F_REAL TL_cell_crossing = delta_of_L[MAX_LEVEL][0]/( vec_len(V_sw) * dt);

	F_REAL min_sw_crosses_cell = (TL_cell_crossing*loop_time_mpi)/60.;
#if defined(use_vectorclass)
	F_REAL min_sw_crosses_Box  = RB_X*(BlkNds_X-2)*min_sw_crosses_cell *pow(double(2),double(MAX_LEVEL));
	F_REAL BoxCrossings_since_start = (1.*TL) / (pow(2.,double(MAX_LEVEL))*RB_X*(BlkNds_X-2)*TL_cell_crossing);
#else
	F_REAL min_sw_crosses_Box  = RB_X*(BlkNds_X-2)*min_sw_crosses_cell *pow(2,MAX_LEVEL);
	F_REAL BoxCrossings_since_start = (1.*TL) / (pow(2,MAX_LEVEL)*RB_X*(BlkNds_X-2)*TL_cell_crossing);
#endif

	char crossing_value[10];
	sprintf(crossing_value,"%2.3f",BoxCrossings_since_start);

	F_REAL min_Simu_finished = (loop_time_comp*(TL_MAX-TL))/60.;

#ifdef ESTIMATE_FASTEST_PARTICLE
        log_file <<" - Fastest particle: " << endl;
        bool tofast = false;
        for(int level=0;level<MAX_LEVEL+1;level++)
        {
            log_file << "   -> Level " << level << ": v="<< sqrt(fastest_particle_v2[level]) << " -> "
                     << delta_of_L[level][0]/(sqrt(fastest_particle_v2[level])*dt_particle_of_L[level]) 
                     << " TL to cross cell " << endl;
            if(delta_of_L[level][0]/(sqrt(fastest_particle_v2[level])*dt_particle_of_L[level]) < 1.)
                tofast = true;
        }
        
        if(tofast)
        {

                log_file <<" !!!-----------------------------------------------!!!" << endl;
                log_file <<" !!! WARNING: time step to large for highest level !!!" << endl;
                log_file <<" !!!     THIS MAY RESULT IN SERIOUS ERRORS         !!!" << endl;
                log_file <<" !!!              IN MOVE METHODS                  !!!" << endl;
                log_file <<" !!!          (set dt to lower value)              !!!" << endl;
                log_file <<" !!!                   (or)                        !!!" << endl;
                log_file <<" !!!     (limit maximal v_part in defines.h)       !!!" << endl;
                log_file <<" !!!-----------------------------------------------!!!" << endl;
        }
#endif
	log_file <<" - At active computational speed, solar wind takes" << endl;
	log_file <<"   -> to cross max level cell:    " << int(TL_cell_crossing)       << "TL."  << endl;
	log_file <<"   -> to cross max level cell:    " << int(min_sw_crosses_cell)    << "min"  << endl;
	log_file <<"   -> to cross the entire Box:    " << int(min_sw_crosses_Box/60.) << "h " << int(min_sw_crosses_Box)%60<< "min"  << endl;
	log_file <<"   -> to finish Simulation Run:   " << int(min_Simu_finished/(60.*24.)) << "days " 
						  << int(min_Simu_finished/60.)%24 << "h " 
						  << int(min_Simu_finished)%60<< "min"  << endl;
	log_file <<" - Box Crossings since start:   " << crossing_value << endl;

	

}




