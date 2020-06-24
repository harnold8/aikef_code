


#include "CBlk.h"
#include "parameters.h"
#include "utils.h"
#include "CBlk_Globals.h"
#include "absolute_Globals.h"
#include <assert.h>

#include "defines.h"
#include <iostream>
#include <fstream>
#include <math.h>


using namespace std;


extern D_REAL **Blk_Length_of;
extern WEIGHT_REAL startWeight_of_smallestCell;



extern INT32 *COMPs_FType;
extern D_REAL *q_of, *q2m_of, *q_m_ratio;



//!-------------------------------------------------------------//
//! init_Block_Globals: -					//
//!	call from arbitrary block at appplication begin		//
//!-------------------------------------------------------------//
void init_Block_Globals(void)
{





	//! Alloc free store memory for smooth_coeff
	smooth_coeff = new D_REAL**[3];
	for (short u=0; u<3; u++)
	{
	     smooth_coeff[u] = new D_REAL*[3];

	     for (short v=0; v<3; v++)
	     smooth_coeff[u][v] = new D_REAL[3];

	}

#if defined(use_vectorclass)
	//! Initialize smooth_coeff
	for (short a=-1; a<=1; a++)
	{
	  for (short b=-1; b<=1; b++)
	  {
	    for (short c=-1; c<=1; c++)
	    {
	      smooth_coeff[a+1][b+1][c+1] = pow( 2.,double( -(a*a+b*b+c*c+3) ) );
	    }
	  }
	}
#else
	//! Initialize smooth_coeff
	for (short a=-1; a<=1; a++)
	 for (short b=-1; b<=1; b++)
          for (short c=-1; c<=1; c++)
            smooth_coeff[a+1][b+1][c+1] = pow(2.,-(a*a+b*b+c*c+3));
#endif
	temp_i = new INT32[MAX_LEVEL+1];
	temp_j = new INT32[MAX_LEVEL+1];
	temp_k = new INT32[MAX_LEVEL+1];

	q_of = new D_REAL[num_Charged_Species];
	q2m_of = new D_REAL[num_Charged_Species];
	q_m_ratio = new D_REAL[num_Charged_Species];




	memset(q_of, 0,num_Charged_Species*sizeof(D_REAL));
	memset(q2m_of, 0,num_Charged_Species*sizeof(D_REAL));
	memset(q_m_ratio, 0,num_Charged_Species*sizeof(D_REAL));


	for(INT32 species=0; species<num_Charged_Species; species++)
	{
	  q_of[species]      =  Ion_Charges[species];
	  q2m_of[species]    = (Ion_Charges[species]*Ion_Charges[species])  / Ion_Masses[species];
	  q_m_ratio[species] =  Ion_Charges[species]			    / Ion_Masses[species];

	}

	dt_field_of_L = new D_REAL[MAX_LEVEL+1];
	rdt_field_of_L = new D_REAL[MAX_LEVEL+1];
	dt_particle_of_L = new D_REAL[MAX_LEVEL+1];

	dt_particle_of_L[0] = dt;
	dt_field_of_L[0]    = dt/(num_advance_B_loops*NUM_SUB_CYCLE);
	rdt_field_of_L[0] = 1./dt_field_of_L[0];

	rd_of_L = new D_REAL*[MAX_LEVEL+1];
	rd2_of_L = new D_REAL*[MAX_LEVEL+1];
	delta_of_L = new D_REAL*[MAX_LEVEL+1];
	CellVol_of_L = new D_REAL[MAX_LEVEL+1];
	Blk_Length_of = new D_REAL*[MAX_LEVEL+1];

	//! allocate Level 0
	rd_of_L[0] = new D_REAL[3];
	rd2_of_L[0] = new D_REAL[3];
	delta_of_L[0] = new D_REAL[3];
	Blk_Length_of[0] = new D_REAL[3];



	//! Size of Block in level 0
	Blk_Length_of[0][0] = LX/(RB_X);
	Blk_Length_of[0][1] = LY/(RB_Y);
	Blk_Length_of[0][2] = LZ/(RB_Z);


	//! Size of delta in level 0
	//! distance of two Gridnodes
	//! usually this would be (NP-1) but 2 cells are ghost cells
	//! Tis can be easily understood with sketch:
	//! "BlockGeometry_even_cells.fig"
	delta_of_L[0][0] = Blk_Length_of[0][0]/(BlkNds_X-2);
	delta_of_L[0][1] = Blk_Length_of[0][1]/(BlkNds_Y-2);
	delta_of_L[0][2] = Blk_Length_of[0][2]/(BlkNds_Z-2);

	CellVol_of_L[0] = delta_of_L[0][0] 
			* delta_of_L[0][1]
			* delta_of_L[0][2];

	//! 2. always appears in derivatives
	rd_of_L[0][0] = 1./(2. *delta_of_L[0][0]);
	rd_of_L[0][1] = 1./(2. *delta_of_L[0][1]);
	rd_of_L[0][2] = 1./(2. *delta_of_L[0][2]);

// 	cout << delta_of_L[0][0] << endl;
// 	cout << delta_of_L[0][1] << endl;
// 	cout << delta_of_L[0][2] << endl;

	rd2_of_L[0][0] = 1./(delta_of_L[0][0]*delta_of_L[0][0]);
	rd2_of_L[0][1] = 1./(delta_of_L[0][1]*delta_of_L[0][1]);
	rd2_of_L[0][2] = 1./(delta_of_L[0][2]*delta_of_L[0][2]);

	//! if global time stepping is used, every level has the same dt

 	for(int lev=1; lev <= MAX_LEVEL; lev++)
	{

		//! allocate Level 1 to MAX_LEVEL
		rd_of_L[lev] = new D_REAL[3];
		rd2_of_L[lev] = new D_REAL[3];
		delta_of_L[lev] = new D_REAL[3];
		Blk_Length_of[lev] = new D_REAL[3];


		//! spatial differences (independent of local vs. global stepping)
		for(int comp=0; comp<3; comp++)
		{
			rd_of_L[lev][comp]   = rd_of_L[lev-1][comp] *2.;
			rd2_of_L[lev][comp]   = rd2_of_L[lev-1][comp] *4.;
			delta_of_L[lev][comp]   = delta_of_L[lev-1][comp] *0.5;
			Blk_Length_of[lev][comp] = Blk_Length_of[lev-1][comp]*0.5;
		}
	
		CellVol_of_L[lev] = CellVol_of_L[lev-1]  *0.125;


		//! time differences (dependent of local vs. global stepping)
// 		if(global_time_stepping)
// 		{
		dt_particle_of_L[lev] = dt;
		dt_field_of_L[lev]    = dt/(num_advance_B_loops*NUM_SUB_CYCLE);
		rdt_field_of_L[lev]   = 1./dt_field_of_L[lev];


// 		}
// 		else
// 		{
// 
// 			dt_particle_of_L[lev] = dt_particle_of_L[lev-1]*0.5;
// 			dt_field_of_L[lev]    = dt_field_of_L[lev-1]*0.5;
// 			rdt_field_of_L[lev]   = 1./dt_field_of_L[lev];
// 
// 		}

		log_file << "dt_field_of_L["<<lev<< "]  "
			 << dt_field_of_L[lev]<< endl;
		log_file << "dt_particle_of_L["<<lev<< "]  "
			 << dt_particle_of_L[lev]<< endl;
	}

	//! estimate particle start of smallest cell
	//! (take MPiC od species 1 as reference)
	startWeight_of_smallestCell = CellVol_of_L[MAX_LEVEL]/(1.*optimal_MPiC[0]);

}

//!-------------------------------------------------------------//
//! init_Block: 						//
//!-------------------------------------------------------------//
void CBlock::init_Block(INT32 Blk_Nr, INT32 *blk_ind, INT32 Level, CBlock *Parent)
{


	//! - set "gather_Neighbour" to sixth element of Neighbour
	//!   array do have both right behind each other
	//! ->simplifies many functions 
	gather_Neighbour = Neighbour +6;

	//! store local copy of Level in member variable
	RLevel = Level;

	//! store local copy of Parent in member variable
	parent = Parent;

	//! store local copy of Blk_Index in member variable
	Blk_Index[0] = blk_ind[0];
	Blk_Index[1] = blk_ind[1];
	Blk_Index[2] = blk_ind[2];

	//! set origin of block
	init_blk_origin();

	//! set Blk number for sorting into linked list
	if(!RLevel)
	Block_Nr = Blk_Nr;
	else
	{



		//! In order to use a SFC estimate indices in Box with respective resolution
		//! -> centre position in middle of cell to avois round off errors !!!
		INT32 box_index[3] = {int( (Box_Origin[0] +origin[0])/ Blk_Length_of[RLevel][0] +0.5),
				     int( (Box_Origin[1] +origin[1])/ Blk_Length_of[RLevel][1] +0.5),
				     int( (Box_Origin[2] +origin[2])/ Blk_Length_of[RLevel][2] +0.5)};




		if(use_SFC)
		Block_Nr = SFC_Indices_to_BlkNr(box_index, RLevel);
		else
		Block_Nr = parent->Block_Nr *8
				  +blk_ind[0]*2*2
				  +blk_ind[1]*2
				  +blk_ind[2];




	}

	//! it is assumed that num_root_blocks is a multiple of mpi_num_processes
	INT32 num_RB_each_process;

	if(num_root_blocks%mpi_num_processes==0)
	num_RB_each_process = (num_root_blocks/mpi_num_processes);
	else
	num_RB_each_process = (num_root_blocks/mpi_num_processes) +1;

	//! set responsible_mpi_process
	if(!RLevel)
	{
		//! select responsible process
		//! eg.:
		//! - num_root_blocks = 100;
		//! - mpi_num_processes = 10;
		//! -> num_RB_each_process = 10;
		//! -> Block_Nr 00-09  -> p0
		//! -> Block_Nr 10-19  -> p1
		//! -> Block_Nr 20-29  -> p2
		//! -> .............   -> .............
		//! -> Block_Nr 90-99  -> p9

		//! TODO: no integer at division (num_root_blocks/mpi_num_processes)
		//! eg.:
		//! - num_root_blocks = 10;
		//! - mpi_num_processes = 4;
		//! -> num_RB_each_process = 3;
		//! -> Block_Nr 0-2  -> p0
		//! -> Block_Nr 3-5  -> p1
		//! -> Block_Nr 6-8  -> p2
		//! -> Block_Nr 9-10 -> p3
		responsible_mpi_process = Block_Nr/ num_RB_each_process;
		
		//! count how many blocks are assigned to respective process
		//! UPDATE:
		//! Recount when required, eg. at redistribute blocks
// 		num_rootBlocks_at_process[responsible_mpi_process]++;
	}
	else
	responsible_mpi_process = parent->responsible_mpi_process;





	is_gatherBlk = false;
	initial_refined = false;
	do_receive_particle_from_parent = false;



	num_children = 0;

	memset(do_send_particle_to_childArray, 0, 8*sizeof(bool));

	memset(child_array,0,8*sizeof(CBlock*));
	memset(gather_child_array,0,8*sizeof(CBlock*));

	memset(Neighbour,        0, 6*sizeof(CBlock*));
	memset(gather_Neighbour, 0, 6*sizeof(CBlock*));

	memset(child_flag_refinement, 0, 8*sizeof(bool));
	memset(child_flag_removal, 0, 8*sizeof(bool));

	memset(average_ref_value, 0, 8*sizeof(F_REAL));

	//! this has to be allocated at EVERY process
	Blk_optimal_MPiC = new INT32[num_Charged_Species];

	//! allocate buffer with any number to ensure
	//! there is no conflict in first delete[]
	send_buffer[0] = new D_REAL[5];
	send_buffer[1] = new D_REAL[5];

// 	req_is_MPI_send_completed = new MPI::Request[NUM_REQUESTS];

	next_Blk_of_BlockList = NULL;
	next_Blk_of_GatherBlockList = NULL;



	//! mark block as box boundary if required
	init_blk_boundaries();

	time_process_fields   = 0.;
	time_process_particle = 0.;

	time_process_fields_incChilds  = 0.;
	time_process_particle_incChilds = 0.;



	for(INT32 req=0; req<8; req++)
	req_is_MPI_send_completed[req] = MPI_REQUEST_NULL;

	memset(Field_Type, 0, NUM_FIELDS*sizeof(D_REAL*));



	//! ------------FIELD SOLVER ---------------------------------------
	//! --- set function pointers---------------------------------------
	//! ----------------------------------------------------------------
// 	if(BlkNds_X%2==1) {log_file << "no odd Node Number supported"; exit(1);}
//  	else		 ptr2BV_from_parent = &CBlock::BV_from_parent_EvenNP;

	//! LF_cycle[0] remains unused
	LF_cycle[1] = &CBlock::LF_odd_B0toB1;

	for(int cycle=2; cycle<NUM_SUB_CYCLE; cycle+=2)
	{
	 LF_cycle[cycle] = &CBlock::LF_advance_B_even;
	 LF_cycle[cycle+1] = &CBlock::LF_advance_B_odd;
// 	 log_file << cycle+1 << endl;
	}

	LF_cycle[NUM_SUB_CYCLE+1] = &CBlock::LF_even_B_Solution;


	//! Particle Transmission Realated
	memset(MPiC_Info_package, 0, 8*sizeof(INT32*));
	memset(particle_package , 0, 8*sizeof(particle*));




}

//!-------------------------------------------------------------//
//! set_box_boundaries: 					//
//!-------------------------------------------------------------//
void CBlock::init_blk_origin(void)
{

	//! distinguish between root level and L>0
	if(!RLevel)
	{
		//! meaning the left lower cell vertex
		//! which is no ghost cell ->  Cell(1,1,1)
		//! NOTE:
		//! The geometry is NOT perfectly symmetric
		//! in this formulation. Howerver, in any refinement
		//! level exist a node at the origins position
		//! (see XFig images for details)
		for(int comp=0; comp<3; comp++)
		origin[comp] = Blk_Index[comp] * Blk_Length_of[RLevel][comp] - Box_Origin[comp];

	}	
	else for(int comp=0; comp<3; comp++)
	origin[comp] = parent->origin[comp] + Blk_Index[comp] * Blk_Length_of[RLevel][comp];

}


//!-------------------------------------------------------------//
//! init_blk_boundaries: 					//
//!-------------------------------------------------------------//
void CBlock::init_blk_boundaries(void)
{

	//! Marc as Boundariy Block if necessary
	memset(is_box_boundary,0,8*sizeof(bool));

	//! distinguish between root level and L>0
	if(!RLevel)
	{
		//! 0: xMin, 1:xMax
		if(!use_periodic_bounds[0] && Blk_Index[0] ==      0) is_box_boundary[0] = true;
		if(!use_periodic_bounds[0] && Blk_Index[0] == RB_X-1) is_box_boundary[1] = true;
	
		//! 2: yMin, 3:yMax
		if(!use_periodic_bounds[1] && Blk_Index[1] ==      0) is_box_boundary[2] = true;
		if(!use_periodic_bounds[1] && Blk_Index[1] == RB_Y-1) is_box_boundary[3] = true;
	
		//! 4: zMin, 5:zMax
		if(!use_periodic_bounds[2] && Blk_Index[2] ==      0) is_box_boundary[4] = true;
		if(!use_periodic_bounds[2] && Blk_Index[2] == RB_Z-1) is_box_boundary[5] = true;

	}
	else
	{
		//! 0: xMin, 1:xMax
		if(parent->is_box_boundary[0] && !Blk_Index[0]) is_box_boundary[0] = true;
		if(parent->is_box_boundary[1] &&  Blk_Index[0]) is_box_boundary[1] = true;
	
		//! 2: yMin, 3:yMax
		if(parent->is_box_boundary[2] && !Blk_Index[1]) is_box_boundary[2] = true;
		if(parent->is_box_boundary[3] &&  Blk_Index[1]) is_box_boundary[3] = true;
	
		//! 4: zMin, 5:zMax
		if(parent->is_box_boundary[4] && !Blk_Index[2]) is_box_boundary[4] = true;
		if(parent->is_box_boundary[5] &&  Blk_Index[2]) is_box_boundary[5] = true;
	
	}

		if(   is_box_boundary[0] || is_box_boundary[1] || is_box_boundary[2]
		   || is_box_boundary[3] || is_box_boundary[4] || is_box_boundary[5])
		is_box_boundary[6] = true;


	//! check whether obstacle overlaps with this block
	//! (to avoid some calculations in particle move)

	//!---------- X Direction ----------------------------------
	if(     origin[0]			     >= +R_Obstacle
	     || origin[0] + Blk_Length_of[RLevel][0] <= -R_Obstacle)
	is_box_boundary[7] = false;
	else
	{

		//!---------- Y Direction ----------------------------------
		if(     origin[1]			     >= +R_Obstacle
		     || origin[1] + Blk_Length_of[RLevel][1] <= -R_Obstacle)
		is_box_boundary[7] = false;
		else
		{
	
			//!---------- Z Direction ----------------------------------
			if(   origin[2]			            >= +R_Obstacle
			    || origin[2] + Blk_Length_of[RLevel][2] <= -R_Obstacle)
			is_box_boundary[7] = false;
			else
			is_box_boundary[7] = true;
			
	
		}

	}

	if(use_second_obstacle)
	{	
		//! TODO not sure if it has to be  "-" or "+" Position_SecondObstacle[0]
		//!---------- X Direction ----------------------------------
		if(     origin[0]			     >= +R_SecondObstacle - Position_SecondObstacle[0]
		|| origin[0] + Blk_Length_of[RLevel][0] <= -R_SecondObstacle - Position_SecondObstacle[0])
		is_box_boundary[7] = false;
		else
		{

			//!---------- Y Direction ----------------------------------
			if(     origin[1]			     >= +R_SecondObstacle - Position_SecondObstacle[1]
			|| origin[1] + Blk_Length_of[RLevel][1] <= -R_SecondObstacle - Position_SecondObstacle[1])
			is_box_boundary[7] = false;
			else
			{
		
				//!---------- Z Direction ----------------------------------
				if(   origin[2]			            >= +R_SecondObstacle - Position_SecondObstacle[2]
				|| origin[2] + Blk_Length_of[RLevel][2] <= -R_SecondObstacle - Position_SecondObstacle[2])
				is_box_boundary[7] = false;
				else
				is_box_boundary[7] = true;
				
		
			}

		}
	}
	
}



//!-------------------------------------------------------------//
//! delete_Memory: 								//
//!-------------------------------------------------------------//
void CBlock::delete_Memory(void)
{



	//! - Some memory is allocated at every process,
	//!   indepentend of responsible_mpi_process
	//! - This memory was allocated in the "init_Block",
	//!   while process specific memory is allocated in
	//!   the "alloc_process_specific_Memory" function.
	delete[] Blk_optimal_MPiC;


	delete[] send_buffer[0];
	delete[] send_buffer[1];

	Blk_optimal_MPiC = NULL;
	send_buffer[0] = NULL;
	send_buffer[1] = NULL;

	//! memory below is only allocated at respective process
	if(responsible_mpi_process==mpi_myRank)
	delete_process_specific_Memory();



}


//!-------------------------------------------------------------//
//! alloc_process_specific_Memory: 					//
//!-------------------------------------------------------------//
void CBlock::alloc_GatherMemory(void)
{


	//! Just allocate Memory relevant for Gather procedure
	Field_Type[id_rho_np1] =  new D_REAL[(10+num_Charged_Species)*num_nodes_in_block];
	memset(Field_Type[id_rho_np1],0,
		(10+num_Charged_Species)*num_nodes_in_block*sizeof(D_REAL));

	//! (Array 1)
	Field_Type[id_UIp_Gam_allRhoSpec] = Field_Type[id_rho_np1] +num_nodes_in_block;

	//! (Array 2-4)
	//!NOTE:
	//! as moment memory is only used by gather blocks which
	//! do not need the data after gather function is finished,
	//! both Ui_plus as well as Ui_minus may use the same memory.
	//! this does not work for common blocks wgere Ui_plus
	//! must not be ovwerwritten !
	//! The same is true when vth2 is gathered, so use memory Block 
	//! for three fields
	//! 
	//! Maybe it would be sufficient to allocate a single vector field for 
	//! the entire gather memory ?!?
	Field_Type[id_UI_plus]      = Field_Type[id_UIp_Gam_allRhoSpec];
	Field_Type[id_UI_minus]    = Field_Type[id_UIp_Gam_allRhoSpec];
	Field_Type[id_PISpecies1] = Field_Type[id_UIp_Gam_allRhoSpec];
	//! (Array 5-7)
	Field_Type[id_Gam]   = Field_Type[id_UI_plus] +3*num_nodes_in_block;

	//! (num_Charged_Species Arrays)
	Field_Type[id_allRhoSpecies] = Field_Type[id_Gam] +3*num_nodes_in_block;

	for(int species=0; species<num_Charged_Species; species++)
	Field_Type[id_rhoSpecies1 +species] = Field_Type[id_allRhoSpecies] + species *num_nodes_in_block;

	//! rotB is used as a temp field for Ji_species
	//! (Array 8-10)
	Field_Type[id_rotB]  = Field_Type[id_allRhoSpecies] +num_Charged_Species*num_nodes_in_block;

	for(int species=0; species<num_Charged_Species; species++)
	Field_Type[id_UI_Species1 +species] = Field_Type[id_rotB];
}


//!-------------------------------------------------------------//
//! delete_GatherMemory: 							//
//!-------------------------------------------------------------//
void CBlock::delete_GatherMemory(void)
{

	delete[] Field_Type[id_rho_np1];


}


//!-------------------------------------------------------------//
//! alloc_process_specific_Memory: 								//
//!-------------------------------------------------------------//
void CBlock::alloc_process_specific_Memory(void)
{

	//!------- Memory Allocation -----------------------------
	//! The following allocation turned out 
	//! to allow the fastest memory acces
	//!-------------------------------------------------------

	//! allocate Flag (inside or outside planet)
	Flag = new bool[num_nodes_in_block];
	memset(Flag,0,num_nodes_in_block*sizeof(bool));

	//! allocate Flag (inside or outside planetary core)
	Core_Flag = new bool[num_nodes_in_block];
	memset(Core_Flag,0,num_nodes_in_block*sizeof(bool));




	//! Allocate Memory for all Block Fields at once
	//! to achieve improved Data locality

	//! Order of Fields:
	//! BEven->BOdd->BTotal->Bcfbg->EField-> etc...


	//! - 'num_Charged_Species' for each density field.
	//! - 4*num_externRhoVelocityFields for ExDens, ExV1, ExV2, ExV3;
	//! - Even Magnetic Field (Array 1-3)
	Field_Type[id_BEven] = new D_REAL[NUM_SCALAR_FIELDS_EACH_BLK*num_nodes_in_block];
	memset(Field_Type[id_BEven],0,NUM_SCALAR_FIELDS_EACH_BLK*num_nodes_in_block*sizeof(D_REAL));
        
	//! set id to start of entire block memory
	Field_Type[id_ALL_FIELDS] = Field_Type[id_BEven];

	//! Odd Magnetic Field  (Array 4-6)
	Field_Type[id_BOdd] = Field_Type[id_BEven] +3*num_nodes_in_block;

	//! Odd Magnetic Field  (Array 7-9)
	Field_Type[id_BDerivative] = Field_Type[id_BOdd] +3*num_nodes_in_block;

	//! NOTE: use same memory for id_BDerivative & id_gradB 
	Field_Type[id_gradB] = Field_Type[id_BDerivative];

	//! Total Magnetic Field (Array 10-12)
	Field_Type[id_BTotal] = Field_Type[id_BDerivative] +3*num_nodes_in_block;

	//! Curl Free Magnetic Backround Field  (Array 13-15)
	Field_Type[id_Bcfbg] = Field_Type[id_BTotal] +3*num_nodes_in_block; 

	//! Curl Free Magnetic Backround Field  (Array 13-15)
	Field_Type[id_Bcfbg_HIMesh] = Field_Type[id_Bcfbg] +3*num_nodes_in_block; 

	//! Electric Field (Array 16-18)
	Field_Type[id_EField]  = Field_Type[id_Bcfbg_HIMesh]  +3*num_nodes_in_block;


	//! density (Array 19)
	Field_Type[id_rho_np1] = Field_Type[id_EField] + 3*num_nodes_in_block;

	//!-----------------------------------------------------------------
	//!---- U GAM and rho species are sequently arranged in Memory -----
	//!---- and builad a new "vector field with 3 +3 + num_Charged_Species --
	//!---- components.								
	//!-----------------------------------------------------------------
	//! Mean Velocity Field (Array 20-22)
	Field_Type[id_UI_plus]       = Field_Type[id_rho_np1] + num_nodes_in_block;
	Field_Type[id_UIp_Gam_allRhoSpec] = Field_Type[id_UI_plus];
	//! Gam Field (Array 23-25)
	Field_Type[id_Gam]     = Field_Type[id_UI_plus] +3*num_nodes_in_block;

	//! allocate mem for all ionspecies density, velocity, ....
	//! (needed for grad pE)
	//! NOTE:
	//! here are num_Charged_Species adational D_REAL
	//! which are left out in the following counting
        
        //! Ion Density for all species
        //! Scalar
	Field_Type[id_allRhoSpecies]    = Field_Type[id_Gam] +3*num_nodes_in_block;
	for(int species=0; species<num_Charged_Species; species++) 
            Field_Type[id_rhoSpecies1 +species] = Field_Type[id_allRhoSpecies] + species *num_nodes_in_block;          
        
        //! PI for all species
        //! Vector (vth2x, vth2y. vth2z)
        Field_Type[id_allPISpecies]    = Field_Type[id_allRhoSpecies] +num_Charged_Species*num_nodes_in_block;
        for(int species=0; species<num_Charged_Species; species++) 
            Field_Type[id_PISpecies1 +species] = Field_Type[id_allPISpecies] + species *3*num_nodes_in_block;
        
        //! Gyroradius for all species
        //! Scalar
        Field_Type[id_allGyroSpecies]    = Field_Type[id_allPISpecies] +num_Charged_Species*3*num_nodes_in_block;
        for(int species=0; species<num_Charged_Species; species++) 
            Field_Type[id_gyro_Species1   +species] = Field_Type[id_allGyroSpecies] + species *num_nodes_in_block;
        
        //! Electron Gyroradius for all species
        //! Scalar
        Field_Type[id_allGyroELSpecies]    = Field_Type[id_allGyroSpecies] +num_Charged_Species*num_nodes_in_block;
        for(int species=0; species<num_Charged_Species; species++) 
            Field_Type[id_gyro_el_Species1   +species] = Field_Type[id_allGyroELSpecies] + species *num_nodes_in_block;
	
	//! recombination field for all species
        //! Scalar
        Field_Type[id_recomb_Species1]    = Field_Type[id_allGyroELSpecies] +num_Charged_Species*num_nodes_in_block;
        for(int species=1; species<num_Charged_Species; species++) 
            Field_Type[id_recomb_Species1   +species] = Field_Type[id_recomb_Species1] + species *num_nodes_in_block;
        
        //! PE for all species
        //! Scalar
        Field_Type[id_allPESpecies]    = Field_Type[id_recomb_Species1] +num_Charged_Species*num_nodes_in_block;
        for(int species=0; species<num_Charged_Species; species++) 
            Field_Type[id_PESpecies1 +species] = Field_Type[id_allPESpecies] + species *num_nodes_in_block;
        
        //! UI for all species
        //! Vector    
        Field_Type[id_allUISpecies]    = Field_Type[id_allPESpecies] +num_Charged_Species*num_nodes_in_block;
        for(int species=0; species<num_Charged_Species; species++) 
            Field_Type[id_UI_Species1   +species] = Field_Type[id_allUISpecies] + species *3*num_nodes_in_block;
        
        //! Forces for all species
        //! vector    
        Field_Type[id_allForcesSpecies]    = Field_Type[id_allUISpecies] +num_Charged_Species*3*num_nodes_in_block;
        for(int species=0; species<num_Charged_Species; species++)  
            Field_Type[id_ForceSpecies1 +species] = Field_Type[id_allForcesSpecies] + species *3*num_nodes_in_block;
        
	//!------------------------------------------------------------------
	//!---- END: U, GAM ,rho species ------------------------------------
	//!------------------------------------------------------------------

	//! Total electron pressure
        Field_Type[id_PEtotal] = Field_Type[id_allForcesSpecies] + num_Charged_Species*3*num_nodes_in_block;

	//! LAM     (Array 26)
        Field_Type[id_Lam] = Field_Type[id_PEtotal] + num_nodes_in_block;

	//! rho_n     (Array 27) 
	Field_Type[id_rho_n] = Field_Type[id_Lam] +num_nodes_in_block;

	//! id_rho_rez (Array 28)
	Field_Type[id_rho_rez] = Field_Type[id_rho_n] +num_nodes_in_block;

	//! id_UI_minus (Array 29-31)
	Field_Type[id_UI_minus] = Field_Type[id_rho_rez] +num_nodes_in_block;


	//! id_Eta (Array 32)
	Field_Type[id_Eta]  = Field_Type[id_UI_minus] + 3*num_nodes_in_block;

	//! id_PhiDC (Array 33)
	Field_Type[id_PhiDC]  = Field_Type[id_Eta] + num_nodes_in_block;
        
        Field_Type[id_ElectronTemperature]  = Field_Type[id_PhiDC] + num_nodes_in_block;
        Field_Type[id_rho_np1_recombined] =  Field_Type[id_ElectronTemperature] + num_nodes_in_block;


	//!---------------SRC1-FIELDS------------------------------------
	//! Use (Array 34-36) to calculate temp fields like 
	//! div & rot & rho, so several Field_Type's
	//! will point to start of this memory Block
        Field_Type[id_divB] = Field_Type[id_rho_np1_recombined] +num_nodes_in_block;
	Field_Type[id_rotB] = Field_Type[id_divB] +num_nodes_in_block;
	Field_Type[id_gradPI0 ] = Field_Type[id_rotB];
	Field_Type[id_divE]     = Field_Type[id_rotB];
	Field_Type[id_divrotB]  = Field_Type[id_divB];
	Field_Type[id_divU]     = Field_Type[id_divB];
	Field_Type[id_gradPE]   = Field_Type[id_rotB];
        Field_Type[id_PMagnetic] = Field_Type[id_rotB];
        Field_Type[id_FieldTime   ] = Field_Type[id_rotB];
        Field_Type[id_ParticleTime] = Field_Type[id_rotB];


	Field_Type[id_Refine_Rating] = Field_Type[id_rotB] +3*num_nodes_in_block;

        Field_Type[id_scratch_scalar] = Field_Type[id_Refine_Rating] + num_nodes_in_block;
	Field_Type[id_scratch_vector] = Field_Type[id_scratch_scalar] + num_nodes_in_block;

	
#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
        Field_Type[id_allRhoNeutralSpecies]    = Field_Type[id_scratch_vector] + 3*num_nodes_in_block;
        for(INT32 species=0; species<num_Neutral_Species; ++species)
            Field_Type[id_numberdensity_neutralSpecies1 + species] = Field_Type[id_allRhoNeutralSpecies] + species *num_nodes_in_block;
        
        Field_Type[id_allUNeutralSpecies]    = Field_Type[id_allRhoNeutralSpecies] + num_Neutral_Species*num_nodes_in_block;
        for(INT32 species=0; species<num_Neutral_Species; ++species)
            Field_Type[id_velocity_neutralSpecies1 + species] = Field_Type[id_allUNeutralSpecies] + species *3*num_nodes_in_block;
        
        Field_Type[id_allPNeutralSpecies]    = Field_Type[id_allUNeutralSpecies] + num_Neutral_Species*3*num_nodes_in_block;
        for(INT32 species=0; species<num_Neutral_Species; ++species)
            Field_Type[id_pressure_neutralSpecies1 + species] = Field_Type[id_allPNeutralSpecies] + species *num_nodes_in_block;
        
        Field_Type[id_allnewBetaNeutralSpecies]    = Field_Type[id_allPNeutralSpecies] + num_Neutral_Species*num_nodes_in_block;
        for(INT32 species=0; species<num_Neutral_Species; ++species)
            Field_Type[id_new_electron_beta_neutralSpecies1 + species] = Field_Type[id_allnewBetaNeutralSpecies] + species*num_nodes_in_block;
        
        //! PTotal-Field and Pressure- /Velocity Fields must not overlap
        Field_Type[id_PTotal] = Field_Type[id_allnewBetaNeutralSpecies] + num_Neutral_Species*num_nodes_in_block;
#else
        Field_Type[id_PTotal] = Field_Type[id_scratch_vector] + 3*num_nodes_in_block;
        
#endif
	

	
#if defined(use_dust_species_as_field) && defined(use_ion_production_as_field)
	
	Field_Type[id_density_ionProdSpecies1]    = Field_Type[id_PTotal] + num_nodes_in_block;
	for(INT32 species=1; species<num_ion_prod_fields; ++species)
            Field_Type[id_density_ionProdSpecies1 + species] = Field_Type[id_density_ionProdSpecies1] + species *num_nodes_in_block;
	
	Field_Type[id_velocity_ionProdSpecies1]    = Field_Type[id_density_ionProdSpecies1] + num_ion_prod_fields*num_nodes_in_block;
        for(INT32 species=1; species<num_ion_prod_fields; ++species)
            Field_Type[id_velocity_ionProdSpecies1 + species] = Field_Type[id_velocity_ionProdSpecies1] + species *3*num_nodes_in_block;
	
	Field_Type[id_density_dustSpecies1]    = Field_Type[id_velocity_ionProdSpecies1] + 3*num_ion_prod_fields*num_nodes_in_block;
        for(INT32 species=1; species<num_Dust_Species; ++species)
            Field_Type[id_density_dustSpecies1 + species] = Field_Type[id_density_dustSpecies1] + species *num_nodes_in_block;
	
	Field_Type[id_velocity_dustSpecies1]    = Field_Type[id_density_dustSpecies1] + num_Dust_Species*num_nodes_in_block;
        for(INT32 species=1; species<num_Dust_Species; ++species)
            Field_Type[id_velocity_dustSpecies1 + species] = Field_Type[id_velocity_dustSpecies1] + species *3*num_nodes_in_block;
	
	Field_Type[id_externRho1] = Field_Type[id_velocity_dustSpecies1] + 3*num_Dust_Species*num_nodes_in_block;
	
#elif !(defined(use_dust_species_as_field)) && defined(use_ion_production_as_field)
	
	Field_Type[id_density_ionProdSpecies1]    = Field_Type[id_PTotal] + num_nodes_in_block;
	for(INT32 species=1; species<num_ion_prod_fields; ++species)
            Field_Type[id_density_ionProdSpecies1 + species] = Field_Type[id_density_ionProdSpecies1] + species *num_nodes_in_block;
	
	Field_Type[id_velocity_ionProdSpecies1]    = Field_Type[id_density_ionProdSpecies1] + num_ion_prod_fields*num_nodes_in_block;
        for(INT32 species=1; species<num_ion_prod_fields; ++species)
            Field_Type[id_velocity_ionProdSpecies1 + species] = Field_Type[id_velocity_ionProdSpecies1] + species *3*num_nodes_in_block;

	Field_Type[id_externRho1] = Field_Type[id_velocity_ionProdSpecies1] + 3*num_ion_prod_fields*num_nodes_in_block;
	
#elif defined(use_dust_species_as_field) && !(defined(use_ion_production_as_field))

	Field_Type[id_density_dustSpecies1]    = Field_Type[id_PTotal] + num_nodes_in_block;
        for(INT32 species=1; species<num_Dust_Species; ++species)
            Field_Type[id_density_dustSpecies1 + species] = Field_Type[id_density_dustSpecies1] + species *num_nodes_in_block;
	
	Field_Type[id_velocity_dustSpecies1]    = Field_Type[id_density_dustSpecies1] + num_Dust_Species*num_nodes_in_block;
        for(INT32 species=1; species<num_Dust_Species; ++species)
            Field_Type[id_velocity_dustSpecies1 + species] = Field_Type[id_velocity_dustSpecies1] + species *3*num_nodes_in_block;
	
	Field_Type[id_externRho1] = Field_Type[id_velocity_dustSpecies1] + 3*num_Dust_Species*num_nodes_in_block;
	
#else 	
	      
	Field_Type[id_externRho1] = Field_Type[id_PTotal] + num_nodes_in_block;
	
#endif	
	
	//!-------------------------------------------------------------
	//!--------------------------------------------------------------

	for(INT32 rho_extern=1; rho_extern<num_externRhoVelocityFields; rho_extern++)
            Field_Type[id_externRho1 +rho_extern] = Field_Type[id_externRho1] + rho_extern*num_nodes_in_block;

	//! put id_extern_Ui1 right behind externRho1 field
	Field_Type[id_extern_Ui1] = Field_Type[id_externRho1] +num_externRhoVelocityFields*num_nodes_in_block;
	for(INT32 rho_extern=1; rho_extern<num_externRhoVelocityFields; rho_extern++)
            Field_Type[id_extern_Ui1 +rho_extern] = Field_Type[id_extern_Ui1] + rho_extern*(3*num_nodes_in_block);

	//! average fields
	Field_Type[id_average_Field1] = Field_Type[id_extern_Ui1] +num_externRhoVelocityFields*(3*num_nodes_in_block);

	for(INT32 avarage_field=1; avarage_field<num_average_fields; avarage_field++)
            Field_Type[id_average_Field1 +avarage_field] = Field_Type[id_average_Field1 +avarage_field-1] +COMPs_FType[id_average_Field1 +avarage_field-1]*num_nodes_in_block;



	//!----------------------------------------------------------------
	//! --- Particle related memory allocation ------------------------
	//!----------------------------------------------------------------
	pArray = new particle**[num_Charged_Species];
	num_MPiC = new INT32*[num_Charged_Species];
	size_of_pArray = new INT32*[num_Charged_Species];
	
	for(short species=0; species<num_Charged_Species; species++)
	{

		
		pArray[species] = new particle*[num_nodes_in_block];
		memset(pArray[species],0,num_nodes_in_block*sizeof(particle*));

		for(INT32 node=0; node<num_nodes_in_block; node++)
		pArray[species][node] = NULL;

		//! storage for number of particle for each species and node
		//! (=num_Charged_Species *num_nodes_in_block entries)
		num_MPiC[species] = new INT32[num_nodes_in_block];
		memset(num_MPiC[species],0,num_nodes_in_block*sizeof(INT32));
                
		//! storage for size of PArray for each species and node
		//! (=num_Charged_Species *num_nodes_in_block entries)
		size_of_pArray[species] = new INT32[num_nodes_in_block];
		memset(size_of_pArray[species],0 ,num_nodes_in_block*sizeof(INT32));

	}

	//! !!! DO NOT ALLOCATE PATICLE MEMORY IN ADVANCE !!!
	//! In case particle memory is assigned to each block in advance, many parents whose
	//! cell are empty will get memory as well as rarely used ghost cells. This can
	//! easily result in twice as much memory consumption
/*
	INT32 u_v_w;

	for(short species=0; species<num_Charged_Species; species++)
	 for(u_v_w=0; u_v_w<num_nodes_in_block; u_v_w++)
// 	for(INT32 u=1; u< BlkNds_X-1; u++)
// 	  for(INT32 v=1; v< BlkNds_Y-1; v++)
// 	    for(INT32 w=1; w< BlkNds_Z-1; w++)
	    {

// 		u_v_w  = u*BlkNds_Y*BlkNds_Z 
// 			+v*BlkNds_Z 
// 			+w;

		pArray[species][u_v_w] = new particle[initial_pArray_Size[species]];
		//! allocate memory for particle
// 		pArray[species][u_v_w] = new particle[initial_pArray_Size[species]];

		//! protocol allocated memory for each cell
// 		size_of_pArray[species][u_v_w] = initial_pArray_Size[species];


		//! protocol total allocated memory
// 		num_particle_total_storage +=initial_pArray_Size[species];
	 }

*/


// 	static bool ids_written = false;
// 	
// 	if( ! ids_written)
// 	{
// 	  ids_written = true;
// 	  
// 	  log_file << "   id_PESpecies1                        = " << id_PESpecies1 << endl;
// 	  log_file << "   id_PEtotal                           = " << id_PEtotal << endl;
//   #if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
// 	  for(INT32 species=0; species<num_Neutral_Species; ++species)
// 	  {
// 	    log_file << "   id_numberdensity_neutralSpecies1 + " << species << " = " << id_numberdensity_neutralSpecies1+species << endl;
// 	    log_file << "   id_velocity_neutralSpecies1      + " << species << " = " << id_velocity_neutralSpecies1+species << endl;
// 	    log_file << "   id_pressure_neutralSpecies1      + " << species << " = " << id_pressure_neutralSpecies1+species << endl;
// 	  }
// 	  log_file << "   id_PE_even                           = " << id_PE_even << endl;
// 	  log_file << "   id_PE_odd                            = " << id_PE_odd << endl;
//   #endif
// 	  for(INT32 extern_field=0; extern_field<num_externRhoVelocityFields; ++extern_field)
// 	  {
// 	    log_file << "   id_externRho1 + " << extern_field << "                   = " << id_externRho1+extern_field << endl;
// 	    log_file << "   id_extern_Ui1 + " << extern_field << "                   = " << id_extern_Ui1+extern_field << endl;
// 	  }
// 	  for(INT32 average_field=0; average_field<num_externRhoVelocityFields; ++average_field)
// 	  {
// 	    log_file << "   id_average_Field1 + " << average_field << "               = " << id_average_Field1+average_field << endl;
// 	  }
// 	  log_file << "   id_divrotB                           = " << id_divrotB << endl;
// 	  log_file << "   id_divU                              = " << id_divU << endl;
// 	  log_file << "   id_gradPE                            = " << id_gradPE << endl;
// 	  
// 	  
// 	  
// 	  log_file << "   Field_Type[id_PESpecies1]                        = " << Field_Type[id_PESpecies1] << endl;
// 	  log_file << "   Field_Type[id_PEtotal]                           = " << Field_Type[id_PEtotal] << endl;
//   #if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
// 	  for(INT32 species=0; species<num_Neutral_Species; ++species)
// 	  {
// 	    log_file << "   Field_Type[id_numberdensity_neutralSpecies1 + " << species << "] = "
// 		      << Field_Type[id_numberdensity_neutralSpecies1+species] << endl;
// 	    log_file << "   Field_Type[id_velocity_neutralSpecies1      + " << species << "] = "
// 		      << Field_Type[id_velocity_neutralSpecies1+species] << endl;
// 	    log_file << "   Field_Type[id_pressure_neutralSpecies1      + " << species << "] = "
// 		      << Field_Type[id_pressure_neutralSpecies1+species] << endl;
// 	  }
// 	  log_file << "   Field_Type[id_PE_even]                           = " << Field_Type[id_PE_even] << endl;
// 	  log_file << "   Field_Type[id_PE_odd]                            = " << Field_Type[id_PE_odd] << endl;
//   #endif
// 	  for(INT32 extern_field=0; extern_field<num_externRhoVelocityFields; ++extern_field)
// 	  {
// 	    log_file << "   Field_Type[id_externRho1 + " << extern_field << "]                   = "
// 		    << Field_Type[id_externRho1+extern_field] << endl;
// 	    log_file << "   Field_Type[id_extern_Ui1 + " << extern_field << "]                   = "
// 		    << Field_Type[id_extern_Ui1+extern_field] << endl;
// 	  }
// 	  for(INT32 average_field=0; average_field<num_externRhoVelocityFields; ++average_field)
// 	  {
// 	    log_file << "   Field_Type[id_average_Field1 + " << average_field << "]               = "
// 		    << Field_Type[id_average_Field1+average_field] << endl;
// 	  }
// 	  log_file << "   Field_Type[id_divrotB]                           = " << Field_Type[id_divrotB] << endl;
// 	  log_file << "   Field_Type[id_divU]                              = " << Field_Type[id_divU] << endl;
// 	  log_file << "   Field_Type[id_gradPE]                            = " << Field_Type[id_gradPE] << endl;
// 	  
// 	}


}



//!-------------------------------------------------------------//
//! delete_Memory: 						//
//! - memory below is exclusively allocated				//
//!   at respective process							//
//!-------------------------------------------------------------//
void CBlock::delete_process_specific_Memory(void)
{


	delete[] Flag;
	delete[] Core_Flag;
	delete[] Field_Type[id_BEven];


	Flag      = NULL;
	Core_Flag = NULL;

	memset(Field_Type, 0, NUM_FIELDS*sizeof(D_REAL*));

	//! delete pArrays of each node
	for(short species=0; species<num_Charged_Species; species++)
	 for(INT32 node=0; node<num_nodes_in_block; node++)
	 {

		if(pArray[species][node]!=NULL)
		{
			delete[] pArray[species][node];
	
			//! update memory information
			num_particle_total_storage -= size_of_pArray[species][node];
		}
	 }


	
	for(short species=0; species<num_Charged_Species; species++)
	{
		delete[] num_MPiC[species];
		delete[] pArray[species];
		delete[] size_of_pArray[species];
	}


	delete[] num_MPiC;
	delete[] pArray;
	delete[] size_of_pArray;

	memset(num_MPiC,        0, num_Charged_Species*sizeof(INT32*));
	memset(pArray,          0, num_Charged_Species*sizeof(particle*));
	memset(size_of_pArray , 0, num_Charged_Species*sizeof(INT32*));


}

//!-------------------------------------------------------------//
//! manage_parent_recv_memory							//
//!-------------------------------------------------------------//
void CBlock::manage_parent_recv_memory(bool allocate)
{

	if(allocate)
	{

		//! alloc field of size "id_UIp_Gam_allRhoSpec"
		//! which is the largest array
		Field_Type[0] = new D_REAL[COMPs_FType[id_UIp_Gam_allRhoSpec] *num_nodes_in_block];
	
	
		//!assign pointer to each field
		for(INT32 field=1; field<NUM_FIELDS; field++)
		Field_Type[field] = Field_Type[0];
// 		alloc_process_specific_Memory();
// 		set_Fields();

	}
	else
	{


		//! delete
		delete[] Field_Type[0];
		memset(Field_Type, 0, NUM_FIELDS*sizeof(D_REAL*));

	}


}



/* old VErsion by CK - replaced by new Version(HK) in CBlk_Plume
//!-------------------------------------------------------------//
//! set_RhoUi_extern: 								//
//!-------------------------------------------------------------//
void CBlock::set_RhoUi_extern(INT32 nr_extern_rho,
			      INT32  *num_extern_Nodes,
			      D_REAL *extern_Origin,
			      D_REAL *extern_Length,
			      D_REAL *extern_rho,
			      INT32 &num_values_not_in_extern_box)
{



	INT32  ind[3];
	INT32  num_extern_box_nodes;
	INT32 a,b,c, a_b_c;
	INT32 ap1_b_c, a_bp1_c, a_b_cp1;
	INT32 ap1_bp1_c, ap1_b_cp1, a_bp1_cp1, ap1_bp1_cp1;

	D_REAL *RHO, *UiX, *UiY, *UiZ;
	D_REAL *extern_UiX, *extern_UiY, *extern_UiZ;
	D_REAL r[3], extern_delta[3], shape_func[8];

	PARTICLE_REAL x_BlockNode[3], x_extern[3], cell_intern_r[3];
	memset(cell_intern_r,0,3*sizeof(PARTICLE_REAL));


	for(INT32 comp=0; comp<3; comp++)
	extern_delta[comp] = extern_Length[comp]/(num_extern_Nodes[comp]-1);


	//! Set pointer to extern field
	num_extern_box_nodes =   num_extern_Nodes[0]
				*num_extern_Nodes[1]
				*num_extern_Nodes[2];

	extern_UiX = extern_rho +num_extern_box_nodes;
	extern_UiY = extern_UiX +num_extern_box_nodes;
	extern_UiZ = extern_UiY +num_extern_box_nodes;

	//! Set pointer to intern (Block) field
	RHO = Field_Type[id_externRho1 +nr_extern_rho];

	UiX = Field_Type[id_extern_Ui1 +nr_extern_rho];
	UiY = UiX +num_nodes_in_block;
	UiZ = UiY +num_nodes_in_block;

	//! use a,b,c for extern mesh 
	//! use i,j,k for intern mesh
	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
	
	
		i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
		//! get Coordinate of intern Block Node 
		intern2normedCoords(x_BlockNode, cell_intern_r, ind);
		
		//! find lower next neighbour node in read Mesh
		x_extern[0] = (x_BlockNode[0] +extern_Origin[0])/extern_delta[0];
		x_extern[1] = (x_BlockNode[1] +extern_Origin[1])/extern_delta[1];
		x_extern[2] = (x_BlockNode[2] +extern_Origin[2])/extern_delta[2];

		 a  = int(x_extern[0]);
		 b  = int(x_extern[1]);
		 c  = int(x_extern[2]);


		if(   (x_extern[0]<0) || a>num_extern_Nodes[0]-2
		    ||  x_extern[1]<0 || b>num_extern_Nodes[1]-2
		    ||  x_extern[2]<0 || c>num_extern_Nodes[2]-2)
		  {

			num_values_not_in_extern_box++;


		  }
		  else
		  {
		  
			r[0] = x_extern[0]-a;
			r[1] = x_extern[1]-b;
			r[2] = x_extern[2]-c;
			
			shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2]);
			shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2]);
			
			shape_func[4] = (   r[0])*(   r[1])*(1.-r[2]);
			shape_func[5] = (   r[0])*(1.-r[1])*(   r[2]);
			shape_func[6] = (1.-r[0])*(   r[1])*(   r[2]);
			shape_func[7] = (   r[0])*(   r[1])*(   r[2]);
			
			//! -----------------------------------------------
			a_b_c   =      a*num_extern_Nodes[1]*num_extern_Nodes[2]    +b*num_extern_Nodes[2]      +c;
			
			//! -----------------------------------------------
			ap1_b_c = (a+1)*num_extern_Nodes[1]*num_extern_Nodes[2]     +b*num_extern_Nodes[2]     +c;
			a_bp1_c =     a*num_extern_Nodes[1]*num_extern_Nodes[2] +(b+1)*num_extern_Nodes[2]     +c;
			a_b_cp1 =     a*num_extern_Nodes[1]*num_extern_Nodes[2]     +b*num_extern_Nodes[2] +(c+1);
			
			//! ------------------------------------------------
			ap1_bp1_c = (a+1)*num_extern_Nodes[1]*num_extern_Nodes[2] +(b+1)*num_extern_Nodes[2]     +c;
			ap1_b_cp1 = (a+1)*num_extern_Nodes[1]*num_extern_Nodes[2]     +b*num_extern_Nodes[2] +(c+1);
			a_bp1_cp1 =     a*num_extern_Nodes[1]*num_extern_Nodes[2] +(b+1)*num_extern_Nodes[2] +(c+1);
			
			ap1_bp1_cp1 = (a+1)*num_extern_Nodes[1]*num_extern_Nodes[2] +(b+1)*num_extern_Nodes[2] +(c+1);

			RHO[i_j_k]  =  extern_rho[  a_b_c] * shape_func[0]
				      +extern_rho[ap1_b_c] * shape_func[1]
				      +extern_rho[a_bp1_c] * shape_func[2]
				      +extern_rho[a_b_cp1] * shape_func[3]
			
				      +extern_rho[  ap1_bp1_c] * shape_func[4]
				      +extern_rho[  ap1_b_cp1] * shape_func[5]
				      +extern_rho[  a_bp1_cp1] * shape_func[6]
				      +extern_rho[ap1_bp1_cp1] * shape_func[7];

				if(RHO[i_j_k]  < 0.)
				{
		
					log_file << endl << "ERROR:" << endl
					<< " Negative Values not allowed." << endl;

					log_file <<  "ind[0]: " << ind[0] << endl;
					log_file <<  "ind[1]: " << ind[1] << endl;
					log_file <<  "ind[2]: " << ind[2] << endl;

					log_file <<  " x_extern[0]: " << x_extern[0] << endl;
					log_file <<  " x_extern[1]: " << x_extern[1] << endl;
					log_file <<  " x_extern[2]: " << x_extern[2] << endl;

					log_file <<  " a: " << a << endl;
					log_file <<  " b: " << b << endl;
					log_file <<  " c: " << c << endl;

					for(int f=0; f<8; f++)
					log_file <<  "shape_func["<<f<<"]: " << shape_func[f] << endl;


					log_file <<  "RHO[i_j_k]: " << RHO[i_j_k]<< endl;
					
						exit(1);
				

				}


			UiX[i_j_k]  =  extern_UiX[  a_b_c] * shape_func[0]
				      +extern_UiX[ap1_b_c] * shape_func[1]
				      +extern_UiX[a_bp1_c] * shape_func[2]
				      +extern_UiX[a_b_cp1] * shape_func[3]
			
				      +extern_UiX[  ap1_bp1_c] * shape_func[4]
				      +extern_UiX[  ap1_b_cp1] * shape_func[5]
				      +extern_UiX[  a_bp1_cp1] * shape_func[6]
				      +extern_UiX[ap1_bp1_cp1] * shape_func[7];

			UiY[i_j_k]  =  extern_UiY[  a_b_c] * shape_func[0]
				      +extern_UiY[ap1_b_c] * shape_func[1]
				      +extern_UiY[a_bp1_c] * shape_func[2]
				      +extern_UiY[a_b_cp1] * shape_func[3]
			
				      +extern_UiY[  ap1_bp1_c] * shape_func[4]
				      +extern_UiY[  ap1_b_cp1] * shape_func[5]
				      +extern_UiY[  a_bp1_cp1] * shape_func[6]
				      +extern_UiY[ap1_bp1_cp1] * shape_func[7];

			UiZ[i_j_k]  =  extern_UiZ[  a_b_c] * shape_func[0]
				      +extern_UiZ[ap1_b_c] * shape_func[1]
				      +extern_UiZ[a_bp1_c] * shape_func[2]
				      +extern_UiZ[a_b_cp1] * shape_func[3]
			
				      +extern_UiZ[  ap1_bp1_c] * shape_func[4]
				      +extern_UiZ[  ap1_b_cp1] * shape_func[5]
				      +extern_UiZ[  a_bp1_cp1] * shape_func[6]
				      +extern_UiZ[ap1_bp1_cp1] * shape_func[7];

		}
	  }

}
*/

//!-------------------------------------------------------------//
//! set_zero_Field: 								//
//!-------------------------------------------------------------//
void CBlock::set_zero_Field_inCore(INT32 field_type)
{

	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	{

		for(INT32 node=0; node<num_nodes_in_block; node++)
		if(Core_Flag[node])
		Field_Type[field_type][node +comp*num_nodes_in_block] = 0.;
		

	}


}


//!-------------------------------------------------------------//
//! set_zero_Field: 								//
//!-------------------------------------------------------------//
void CBlock::set_zero_Field(INT32 field_type)
{
	memset(Field_Type[field_type],0,
	       COMPs_FType[field_type] *num_nodes_in_block*sizeof(D_REAL));

}

//!-------------------------------------------------------------//
//! init_Block: 								//
//!-------------------------------------------------------------//
void CBlock::copy_Field(INT32 dest, INT32 src)
{

	memcpy(Field_Type[dest], Field_Type[src],
		COMPs_FType[src] *num_nodes_in_block*sizeof(D_REAL));


}

//!-------------------------------------------------------------//
//! add_Vector2Field: 							//
//!-------------------------------------------------------------//
void CBlock::add_Vector2Field(INT32 dest, INT32 src, D_REAL* vec)
{

	D_REAL* Dest1 = Field_Type[dest];
	D_REAL* Dest2 = Dest1 +num_nodes_in_block;
	D_REAL* Dest3 = Dest2 +num_nodes_in_block;

	D_REAL* Src1 = Field_Type[src];
	D_REAL* Src2 = Src1 +num_nodes_in_block;
	D_REAL* Src3 = Src2 +num_nodes_in_block;

	for (INT32 i=1; i < BlkNds_X-1; i++)
	 for (INT32 j=1; j < BlkNds_Y-1; j++)
	  for (INT32 k=1; k < BlkNds_Z-1; k++)
	  {
	
	
		i_j_k = i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k;

		Dest1[i_j_k] = Src1[i_j_k] +vec[0];
		Dest2[i_j_k] = Src2[i_j_k] +vec[1];
		Dest3[i_j_k] = Src3[i_j_k] +vec[2];


	  }


}

//!-------------------------------------------------------------//
//! add_Vector2Field: 							//
//!-------------------------------------------------------------//
void CBlock::add_multipliedField(INT32 dest, INT32 src, D_REAL factor)
{
	D_REAL* Src1 = Field_Type[src];
	D_REAL* Dest1 = Field_Type[dest];
	
	for(INT32 node=0; node<COMPs_FType[dest] *num_nodes_in_block; node++) {
	  Dest1[node] += Src1[node] *factor;
	}
}

//!-------------------------------------------------------------//
//! sub_squaredField: 							//
//!-------------------------------------------------------------//
#if defined(use_vectorclass)
void CBlock::sub_squaredField_Multiply(INT32 dest, INT32 src, D_REAL factor)
{


	D_REAL* Dest1 = Field_Type[dest] +0*num_nodes_in_block;
	D_REAL* Dest2 = Field_Type[dest] +1*num_nodes_in_block;
	D_REAL* Dest3 = Field_Type[dest] +2*num_nodes_in_block;


	D_REAL* Src1 = Field_Type[src] +0*num_nodes_in_block;
	D_REAL* Src2 = Field_Type[src] +1*num_nodes_in_block;
	D_REAL* Src3 = Field_Type[src] +2*num_nodes_in_block;
	
	VEC4_D_REAL vec_dest1;
	VEC4_D_REAL vec_dest2;
	VEC4_D_REAL vec_dest3;
	
	VEC4_D_REAL vec_src1;
	VEC4_D_REAL vec_src2;
	VEC4_D_REAL vec_src3;


	//! NOTE: num_nodes_in_block has to be a multiple
	//! of 4 if you use this code!
	for (INT32 node=0; node < num_nodes_in_block; node+=4)
	{
		vec_dest1.load( Dest1+node );
		vec_dest2.load( Dest2+node );
		vec_dest3.load( Dest3+node );
		
		vec_src1.load( Src1+node );
		vec_src2.load( Src2+node );
		vec_src3.load( Src3+node );
		
		vec_dest1 -= vec_src1*vec_src1;
		vec_dest2 -= vec_src2*vec_src2;
		vec_dest3 -= vec_src1*vec_src3;
		
		vec_dest1 *= factor;
		vec_dest2 *= factor;
		vec_dest3 *= factor;
		
		vec_dest1.store( Dest1+node );
		vec_dest2.store( Dest2+node );
		vec_dest3.store( Dest3+node );
	}

}
#else
void CBlock::sub_squaredField_Multiply(INT32 dest, INT32 src, D_REAL factor)
{


	D_REAL* Dest1 = Field_Type[dest] +0*num_nodes_in_block;
	D_REAL* Dest2 = Field_Type[dest] +1*num_nodes_in_block;
	D_REAL* Dest3 = Field_Type[dest] +2*num_nodes_in_block;


	D_REAL* Src1 = Field_Type[src] +0*num_nodes_in_block;
	D_REAL* Src2 = Field_Type[src] +1*num_nodes_in_block;
	D_REAL* Src3 = Field_Type[src] +2*num_nodes_in_block;



	for (INT32 node=0; node < num_nodes_in_block; node++)
	{

		Dest1[node] -=  (Src1[node]*Src1[node]);
		Dest2[node] -=  (Src2[node]*Src2[node]);
		Dest3[node] -=  (Src3[node]*Src3[node]);
		


		Dest1[node] *= factor;
		Dest2[node] *= factor;
                Dest3[node] *= factor;

	}

}
#endif

//!-------------------------------------------------------------//
//! sub_squaredField: 							//
//!-------------------------------------------------------------//
void CBlock::Multiply_fields(INT32 dest, INT32 src1, INT32 src2)
{

	D_REAL* Dest = Field_Type[dest];
	D_REAL* Src1 = Field_Type[src1];
	D_REAL* Src2 = Field_Type[src2];

        //! If src1=scalar and src2=scalar
        if(COMPs_FType[src1]==1&&COMPs_FType[src2]==3) {
            for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
            Dest[node] = Src1[node] * Src2[node];
        }
        //! If src1=scalar and src2=vector
        if(COMPs_FType[src1]==1&&COMPs_FType[src2]==3) {
            for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++) {
                Dest[node] = Src1[node%num_nodes_in_block] * Src2[node];
            }
        }
        //! If src1=vector and src2=scalar
        if(COMPs_FType[src1]==3&&COMPs_FType[src2]==1) {
            for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
                Dest[node] = Src1[node] * Src2[node%num_nodes_in_block];
        }
        //! If src1=vector and src2=vector
        if(COMPs_FType[src1]==3&&COMPs_FType[src2]==3) {
            for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
                Dest[node] = Src1[node] * Src2[node];
        }
        

}

//!-------------------------------------------------------------//
//! sub_squaredField:                                                   //
//!-------------------------------------------------------------//
void CBlock::Multiply_field(INT32 dest, INT32 src1, D_REAL src2)
{
     
    D_REAL* Dest = Field_Type[dest];
    D_REAL* Src1 = Field_Type[src1];
    
    //! If src1=scalar and src2=scalar
    if(COMPs_FType[src1]==1) {
        for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
            Dest[node] = Src1[node] * src2;
    }
    //! If src1=vector and src2=scalar
    if(COMPs_FType[src1]==3) {
        for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
             Dest[node] = Src1[node] * src2;
    }

}

//!-------------------------------------------------------------//
//! Add_fields:                                                   //
//!-------------------------------------------------------------//
void CBlock::Add_fields(INT32 dest, INT32 src1, INT32 src2)
{
    
    D_REAL* Dest = Field_Type[dest];
    D_REAL* Src1 = Field_Type[src1];
    D_REAL* Src2 = Field_Type[src2];
    
    //! If src1=scalar and src2=scalar
    if(COMPs_FType[src1]==1&&COMPs_FType[src2]==1) {
        for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
            Dest[node] = Src1[node] + Src2[node];
    }
    //     //! If src1=scalar and src2=vector
    //     if(COMPs_FType[src1]==1&&COMPs_FType[src2]==3) {
        //         for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++) {
            //             Dest[node] = Src1[node%num_nodes_in_block] * Src2[node];
            //         }
            //     }
            //     //! If src1=vector and src2=scalar
            //     if(COMPs_FType[src1]==3&&COMPs_FType[src2]==1) {
                //         for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
                //             Dest[node] = Src1[node] * Src2[node%num_nodes_in_block];
                //     }
    //! If src1=vector and src2=vector
    if(COMPs_FType[src1]==3&&COMPs_FType[src2]==3) {
        for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
            Dest[node] = Src1[node] + Src2[node];
    }
                
                
}

//!-------------------------------------------------------------//
//! Devide_fields:                                                   //
//! Devide scalar fields  src1/src2
//! as long as src>min otherwise src1/minF
//!-------------------------------------------------------------//
void CBlock::Devide_fields(INT32 dest, INT32 src1, INT32 src2, D_REAL min)
{
    
    D_REAL* Dest = Field_Type[dest];
    D_REAL* Src1 = Field_Type[src1];
    D_REAL* Src2 = Field_Type[src2];
    
    //! If src1=scalar and src2=scalar
    if(COMPs_FType[src1]==1&&COMPs_FType[src2]==1) {
        for (INT32 node=0; node < num_nodes_in_block*COMPs_FType[dest]; node++)
            if(Src2[node]>min)
                Dest[node] = Src1[node] / Src2[node];
            else
                Dest[node] = Src1[node] / min;
    }
}



//!-------------------------------------------------------------//
//! sub_squaredField: 							//
//!-------------------------------------------------------------//
void CBlock::square_Field(INT32 dest, INT32 src)
{


	//! scalar Source
	if(COMPs_FType[src]==1)
	{

		D_REAL* Src  = Field_Type[src];
		D_REAL* Dest = Field_Type[dest];


		for (INT32 node=0; node < num_nodes_in_block; node++)
		Dest[node] = Src[node] * Src[node];

	}

	//! vector Source
	if(COMPs_FType[src]==3)
	{

		D_REAL* Src1  = Field_Type[src];
		D_REAL* Src2  = Src1 +num_nodes_in_block;
		D_REAL* Src3  = Src2 +num_nodes_in_block;

		D_REAL* Dest = Field_Type[dest];


		for (INT32 node=0; node < num_nodes_in_block; node++)
		Dest[node] =  Src1[node] * Src1[node]
			     +Src2[node] * Src2[node]
			     +Src3[node] * Src3[node];

	}




}

//!-------------------------------------------------------------//
//! init_Block: 								//
//!-------------------------------------------------------------//
void CBlock::calc_CFBG_Vector(D_REAL* Vector, D_REAL* position)
{


		//! convert degree to rad
		D_REAL x_angle = 3.14159 / 180. *Magnetic_Moment_Angle[0];
		D_REAL z_angle = 3.14159 / 180. *Magnetic_Moment_Angle[1];


		D_REAL pos_RotX[3];
		D_REAL pos_RotZ[3];

		//! rotate around X axis
		pos_RotX[0] = 1.*position[0]           +0.*position[1]           +0.*position[2];
		pos_RotX[1] = 0.*position[0] +cos(x_angle)*position[1] +sin(x_angle)*position[2];
		pos_RotX[2] = 0.*position[0] -sin(x_angle)*position[1] +cos(x_angle)*position[2];


		//! rotate around Z axis
		pos_RotZ[0] = +cos(z_angle)*pos_RotX[0]  +sin(z_angle)*pos_RotX[1] +0.*pos_RotX[2];
		pos_RotZ[1] = -sin(z_angle)*pos_RotX[0]  +cos(z_angle)*pos_RotX[1] +0.*pos_RotX[2];
		pos_RotZ[2] =            0.*pos_RotX[0]            +0.*pos_RotX[1] +1.*pos_RotX[2];


		position[0] = pos_RotZ[0] -Magnetic_Moment_Offset[0];
		position[1] = pos_RotZ[1] -Magnetic_Moment_Offset[1];
		position[2] = pos_RotZ[2] -Magnetic_Moment_Offset[2];


		//! --get-r--------------
		D_REAL r = vec_len(position);



		if(r<=0.)
		{
			log_file << " WARNING: in init CFBG:" <<endl<<" r="<<r<<endl<<endl;
			log_file << "          Better choose Mesh without Node at {0,0,0}"<<endl<<endl;

		}

// 		M_r = vec_scalar(Magnetic_Moment, position);

		D_REAL M_r = (   position[0] *Magnetic_Moment[0]
				+position[1] *Magnetic_Moment[1]
				+position[2] *Magnetic_Moment[2]);

		D_REAL r_3 = r*r*r;
		D_REAL r_5 = r*r*r*r*r;




		Vector[0] = (3.*(M_r)*position[0])/r_5 -Magnetic_Moment[0]/r_3;
		Vector[1] = (3.*(M_r)*position[1])/r_5 -Magnetic_Moment[1]/r_3;
		Vector[2] = (3.*(M_r)*position[2])/r_5 -Magnetic_Moment[2]/r_3;


}


//!-------------------------------------------------------------//
//! init_Block: 								//
//!-------------------------------------------------------------//
void CBlock::set_Fields(void)
{

	//! ATTENTION:
	//! Indices i_j_k  are alredy used in outer loop
	//! across all Root Blocks. So use u_v_w here to 
	//! not confuse



	D_REAL null_vec[3] = {0,0,0};
	D_REAL r_vec[3], r_vec_HIMesh[3];//, eR_vec[3], eTheta_vec[3];
	
	D_REAL* B1_even = Field_Type[id_BEven];
	D_REAL* B2_even = B1_even +num_nodes_in_block;
	D_REAL* B3_even = B2_even +num_nodes_in_block;

	D_REAL* B1_total = Field_Type[id_BTotal];
	D_REAL* B2_total = B1_total +num_nodes_in_block;
	D_REAL* B3_total = B2_total +num_nodes_in_block;

	D_REAL* B1_cfbg = Field_Type[id_Bcfbg];
	D_REAL* B2_cfbg = B1_cfbg +num_nodes_in_block;
	D_REAL* B3_cfbg = B2_cfbg +num_nodes_in_block;

	D_REAL* B1_cfbg_HIMesh = Field_Type[id_Bcfbg_HIMesh];
	D_REAL* B2_cfbg_HIMesh = B1_cfbg_HIMesh +num_nodes_in_block;
	D_REAL* B3_cfbg_HIMesh = B2_cfbg_HIMesh +num_nodes_in_block;


	D_REAL* U1 = Field_Type[id_UI_plus];
	D_REAL* U2 = U1 +num_nodes_in_block;
	D_REAL* U3 = U2 +num_nodes_in_block;

	D_REAL* E1 = Field_Type[id_EField];
	D_REAL* E2 = E1 +num_nodes_in_block;
	D_REAL* E3 = E2 +num_nodes_in_block;


	//! First set id_Bcfbg field:
	//! - Do not set in case Box Boundary:
	
	
	INT32 start[3] = {is_box_boundary[0],
			 is_box_boundary[2],
			 is_box_boundary[4]};
	
	INT32 end[3]   = {is_box_boundary[1],
			 is_box_boundary[3],
			 is_box_boundary[5]};

	if(!is_intern_Magnetic_Moment)
	{
		for(INT32 i=0; i<3; i++)
		{
			start[i]=0;
			end[i]=0;
		}
	}


	if(vec_len(Magnetic_Moment)!=0.)
	for(INT32 u=start[0]; u< BlkNds_X-end[0];u++)
	  for(INT32 v=start[1]; v< BlkNds_Y-end[1];v++)
	    for(INT32 w=start[2]; w< BlkNds_Z-end[2];w++)
	    {

		INT32 u_v_w  = u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;


		INT32 cell_indices[3] = {u,v,w};

		D_REAL CFBG_Vector[3];


		//! set CFBG Field at full integer mesh
		intern2normedCoords(r_vec, null_vec, cell_indices);
		calc_CFBG_Vector(CFBG_Vector, r_vec);

		B1_cfbg[u_v_w] = CFBG_Vector[0];
		B2_cfbg[u_v_w] = CFBG_Vector[1];
		B3_cfbg[u_v_w] = CFBG_Vector[2];


		if(mesh_type==STAGGERED)
		{
			//! set CFBG Field at half integer mesh
			intern2normedCoords_HIMesh(r_vec_HIMesh, null_vec, cell_indices);
			calc_CFBG_Vector(CFBG_Vector, r_vec_HIMesh);
	
			B1_cfbg_HIMesh[u_v_w] = CFBG_Vector[0];
			B2_cfbg_HIMesh[u_v_w] = CFBG_Vector[1];
			B3_cfbg_HIMesh[u_v_w] = CFBG_Vector[2];
		}


	    }
		

		//! Only for dispersion test:
// 		D_REAL A_Wave = 0.3;
// 		D_REAL num_Max_in_Box = 11.;
// 
// 		D_REAL l_Wave = LZ/num_Max_in_Box;
// 		D_REAL k_Wave = 2.*M_PI/l_Wave;
// 
// 		cout << "k_Wave: " << k_Wave << endl;
// 		cout << "l_Wave: " << l_Wave << endl << endl;

	//! Now initialize all Nodes (including GN). 
	//! Boundaries values are never changed 
	//! in Hybrid Cycle, so what is initialized 
	//! here will never change.
	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		INT32 u_v_w  = u*BlkNds_Y*BlkNds_Z 
			     +v*BlkNds_Z 
			     +w;
		

		INT32 cell_indices[3] = {u,v,w};
		intern2normedCoords(r_vec, null_vec, cell_indices);

		

		Field_Type[id_rho_np1][u_v_w] = 1.;
		Field_Type[id_rho_rez][u_v_w] = 1.;
		

		//!--------------------------------------------------------------
		//!---------- FIRST OBSTACLE ------------------------------------
		//!--------------------------------------------------------------
		D_REAL x_shift = 0.;
		D_REAL y_shift = 0.;



		D_REAL dist_to_orig = sqrt( (r_vec[0]-x_shift)*(r_vec[0]-x_shift)
					   +(r_vec[1]-y_shift)*(r_vec[1]-y_shift) 
					   + r_vec[2]*r_vec[2]
					 );
	

		//! ---------------------------------------------------------
	

                U1[u_v_w] =  V_sw[0];
                U2[u_v_w] =  V_sw[1];
                U3[u_v_w] =  V_sw[2];
                

		//! set zero field inside core
		//! -> only cfbg is fixed inside core
		if(dist_to_orig < obstacle_core_fraction *R_Obstacle)
		{
			Core_Flag[u_v_w] = true;

			//! eta should be zero inside obstacle core
			Field_Type[id_Eta][u_v_w] = 0.;
		
			B1_even[u_v_w] = B_sw[0];
			B2_even[u_v_w] = B_sw[1];
			B3_even[u_v_w] = B_sw[2];


	
		}
		else
		{
			//! initialize Block with solar wind BField
			if(!TL_new_Bsw || TL<TL_new_Bsw)
			{


				if(use_hom_B_bounds[_INIT_])
				      {
					//! if using an inhomogeneous background field
					//! define here the analytical equation for the field
					//! TODO: only works with homogeneous B field boundaries so far
					//! i.e. in/outflow B field at boundaries
                                        #if defined(use_inhom_B)
					      B1_even[u_v_w] = 0.;
				              B2_even[u_v_w] = ((r_vec[1]/R_Moon  + 28.74) * 1e-9)/SI_B0;
					      B3_even[u_v_w] = ((-r_vec[2]/R_Moon - 11.00) * 1e-9)/SI_B0;
					 //! otherwise use background field as set in params.cpp
                                         #else
				              B1_even[u_v_w] = B_sw[0];
					      B2_even[u_v_w] = B_sw[1];
					      B3_even[u_v_w] = B_sw[2];
                                         #endif
				}
				else
				{

					D_REAL B_init[3];
					PARTICLE_REAL x[3],v[3];
					cell_centre_normedCoords(x, cell_indices);

					set_inflow_BField(B_init, x);
                                        set_inflow_velocity(v,x,0);

					B1_even[u_v_w] = B_init[0];
					B2_even[u_v_w] = B_init[1];
					B3_even[u_v_w] = B_init[2];
                                        
                                        U1[u_v_w] =  v[0];
                                        U2[u_v_w] =  v[1];
                                        U3[u_v_w] =  v[2];

				}

// 				if(r_vec[0]>0.) B3_even[u_v_w] = -B_sw[2];

// 				if(r_vec[2]<-2.)
// 				{
// 					B1_even[u_v_w] += +A_Wave *cos(k_Wave*r_vec[2]);
// 					B2_even[u_v_w] += -A_Wave *sin(k_Wave*r_vec[2]);
// 				}


			}
			else
			{
				B1_even[u_v_w] = new_B_sw[0];
				B2_even[u_v_w] = new_B_sw[1];
				B3_even[u_v_w] = new_B_sw[2];
	
			}
		}









		if(dist_to_orig < R_Obstacle)
		{

			D_REAL omega[3] = {omega_rotating_obs[0],
					   omega_rotating_obs[1],
					   omega_rotating_obs[2]};

			if(vec_len2(omega)>0.)
			{
				D_REAL v_tang[3];
				vec_cross(v_tang, omega, r_vec);
	
				U1[u_v_w] = v_tang[0];
				U2[u_v_w] = v_tang[1];
				U3[u_v_w] = v_tang[2];
			}
			else
			{

				U1[u_v_w] =  0.;
				U2[u_v_w] =  0.;
				U3[u_v_w] =  0.;

			}

		}

		//!--------------------------------------------------------------
		//!---------- SECOND OBSTACLE ------------------------------------
		//!--------------------------------------------------------------
		if(use_second_obstacle)
		{	
			D_REAL x_shift = Position_SecondObstacle[0];
			D_REAL y_shift = Position_SecondObstacle[1];
			D_REAL z_shift = Position_SecondObstacle[2];

			D_REAL dist_to_orig_obs2 = sqrt( (r_vec[0]-x_shift)*(r_vec[0]-x_shift)
						+(r_vec[1]-y_shift)*(r_vec[1]-y_shift) 
						+ (r_vec[2]-z_shift)*(r_vec[2]-z_shift)
						);
			
			if(dist_to_orig_obs2 < R_SecondObstacle)
			{
			
				B1_even[u_v_w] = B_sw[0];
				B2_even[u_v_w] = B_sw[1];
				B3_even[u_v_w] = B_sw[2];
				
				U1[u_v_w] =  0.;
				U2[u_v_w] =  0.;
				U3[u_v_w] =  0.;
							
				Flag[u_v_w] = true;
			
				Field_Type[id_rho_np1][u_v_w] = MCD_BField;
				Field_Type[id_rho_rez][u_v_w] = 1./MCD_BField;
				
			}
			
		}

		//! Set B_Total to calc EBounds and acceleration in first cycle
		B1_total[u_v_w] = B1_even[u_v_w] + B1_cfbg[u_v_w];
		B2_total[u_v_w] = B2_even[u_v_w] + B2_cfbg[u_v_w];
		B3_total[u_v_w] = B3_even[u_v_w] + B3_cfbg[u_v_w];




		E1[u_v_w] = -(
				    U2[u_v_w] *B3_total[u_v_w]
				  - U3[u_v_w] *B2_total[u_v_w]
			  );

		E2[u_v_w] = -(
				    U3[u_v_w] *B1_total[u_v_w]
				  - U1[u_v_w] *B3_total[u_v_w]
			    );

		E3[u_v_w] = -(
				    U1[u_v_w] *B2_total[u_v_w]
				  - U2[u_v_w] *B1_total[u_v_w]
			  );


	}


	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		INT32 u_v_w  = u*BlkNds_Y*BlkNds_Z 
			     +v*BlkNds_Z 
			     +w;
		

		INT32 cell_indices[3] = {u,v,w};

		if(mesh_type == STAGGERED)
		intern2normedCoords_HIMesh(r_vec, null_vec, cell_indices);
		else
		intern2normedCoords(r_vec, null_vec, cell_indices);
	


		D_REAL dist_to_orig = vec_len(r_vec);

	

		//! set values inside obstacle in obstacle to zero
		if(dist_to_orig < R_Obstacle)
		{

		
			Field_Type[id_rho_rez][u_v_w] = 0.;
		
			Flag[u_v_w] = true;
		
			Field_Type[id_rho_np1][u_v_w] = MCD_BField;
			Field_Type[id_rho_rez][u_v_w] = 1./MCD_BField;
	
		}

	}


	memcpy(Field_Type[id_rho_n], Field_Type[id_rho_np1],num_nodes_in_block *sizeof(D_REAL));
	memcpy(Field_Type[id_UI_minus], Field_Type[id_UI_plus],3 *num_nodes_in_block *sizeof(D_REAL));
	memcpy(Field_Type[id_BOdd], Field_Type[id_BEven],3 *num_nodes_in_block *sizeof(D_REAL));



}



//!-------------------------------------------------------------//
//! reset_average_fields: 						//
//!-------------------------------------------------------------//
void CBlock::reset_average_fields(void)
{



	for(INT32 average_field=0; average_field < num_average_fields; average_field++)
	memset(Field_Type[id_average_Field1 +average_field],
	       0,
	       COMPs_FType[id_average_Field1 +average_field] *num_nodes_in_block*sizeof(D_REAL));



}


//!-------------------------------------------------------------//
//! add_fields_to_average_fields: 					//
//!-------------------------------------------------------------//
void CBlock::add_fields_to_average_fields(void)
{

#if defined(use_vectorclass)
	INT32 max_node;
	VEC4_D_REAL avrg_vec;
	VEC4_D_REAL field_vec;
	D_REAL* avrg;
	
	for(INT32 average_field=0; average_field < num_average_fields; average_field++)
	{
	  max_node = num_nodes_in_block *COMPs_FType[id_average_Field1 +average_field];
	  for(INT32 node=0; node<max_node; node+=4)
	  {
	    avrg = Field_Type[id_average_Field1 +average_field] + node;
	   
	    field_vec.load( Field_Type[IDs_of_Fields_to_average[average_field]] + node );
	    avrg_vec.load(  avrg );
	   
	    avrg_vec += field_vec;
	   
	    avrg_vec.store( avrg );
	  }
	}
#else
	for(INT32 average_field=0; average_field < num_average_fields; average_field++)
	 for(INT32 node=0; node<num_nodes_in_block *COMPs_FType[id_average_Field1 +average_field]; node++)
	  Field_Type[id_average_Field1 +average_field][node]
	   += Field_Type[IDs_of_Fields_to_average[average_field]][node];
#endif


}


//!-------------------------------------------------------------//
//! normalize_average_fields: 						//
//!-------------------------------------------------------------//
void CBlock::normalize_average_fields(void)
{
  
D_REAL rez_num_average_TL = 1./num_TL_average_Fields;
  
#if defined(use_vectorclass)
	INT32 max_node;
	VEC4_D_REAL avrg_vec;
	D_REAL* avrg;
	
	for(INT32 average_field=0; average_field < num_average_fields; average_field++)
	{
	  max_node = num_nodes_in_block *COMPs_FType[id_average_Field1 +average_field];
	  for(INT32 node=0; node<max_node; node+=4)
	  {
	    avrg = Field_Type[id_average_Field1 +average_field] + node;
	    
	    avrg_vec.load(  avrg );
	   
	    avrg_vec *= rez_num_average_TL;
	   
	    avrg_vec.store( avrg );
	  }
	}
#else
	for(INT32 average_field=0; average_field < num_average_fields; average_field++)
	 for(INT32 node=0; node<num_nodes_in_block *COMPs_FType[id_average_Field1 +average_field]; node++)
	  Field_Type[id_average_Field1 +average_field][node] *= rez_num_average_TL;
#endif
}





