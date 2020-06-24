
#include "CHybrid.h"
#include "output_silo.h"
#include "utils.h"
#include "parameters.h"
#include "absolute_Globals.h"
#include "unistd.h"

#include <cmath>

//! time to wait in micro-seconds
#define TIME_TO_WAIT 200000

#if defined(use_ion_production_as_field)     

//!------------------------------------------------------------------------
//! initialize chapman profile from neutral profile
//!------------------------------------------------------------------------
void CHybrid::calculate_Chapman(void)
{

	log_file << "  Calculating Chapman profile for all neutral species  ...       ";

	
	//! calculate vector towards sun
	D_REAL vector_to_sun[3];	
	D_REAL LT=3.;
	D_REAL SSLChapman=-(-19.06*M_PI)/180.+M_PI/2.;
	D_REAL SLTChapman=-M_PI/12.*LT+M_PI*0.5;
	vector_to_sun[0]=1./*sin(SSLChapman)*cos(SLTChapman)*/;
	vector_to_sun[1]=0/*sin(SSLChapman)*sin(SLTChapman)*/;
	vector_to_sun[2]=0/*cos(SSLChapman)*/;	
		
	
	D_REAL** positions;
	INT32 num_positions;
	
	D_REAL* ion_prod;
	D_REAL* ion_prod_allSpec;
	
	D_REAL* neutral_dens;
	D_REAL* neutral_dens_allSpec;
	
	D_REAL ColumnDens;
	
	INT32  ind[3];	

	D_REAL coor_of_node[3], dx;
	
	 for(INT32 neutral_spec=0; neutral_spec<num_Neutral_Species; neutral_spec++)
	  for(INT32 species=0; species<num_Particle_Species; species++)  
	  set_zero_field(id_density_ionProdSpecies1+neutral_spec*num_Particle_Species + species);
	
	//! constants
	//! sigma_i * 10^-22 for m^2	
	//! Source: Schunk and Nagy, 2004
	const D_REAL neutralSpec_abs[3/*NUM_NEUTRAL_SPECIES*/][37] =
	{ //! H2
		{0.0108,0.0798,0.2085,0.4333,0.6037,0.8388,0.7296,1.0180,1.0220,1.4170,1.9420,1.9010,3.0250,3.8700,4.5020,5.3560,6.1680,7.0210,6.8640,7.8110,8.4640,8.4450,
		9.9000,10.7310,11.3720,10.7550,8.6400,7.3390,8.7480,8.2530,0.4763,0.1853,0.,0.0456,0.,0.,0.} ,
	  //! N2
		{0.72,2.261,4.958,8.392,10.210,10.900,10.493,11.670,11.700,13.857,16.910,16.395,21.675,23.160,23.471,24.501,24.130,22.400,22.787,22.790,23.370,23.339,
		31.755,26.540,24.662,120.490,14.180,16.487,33.578,16.992,20.249,9.680,2.240,50.988,0.,0.,0.},
	  //! CH4
		{0.204,0.593,1.496,2.794,3.857,5.053,4.360,6.033,6.059,7.829,10.165,9.776,14.701,18.770,21.449,24.644,27.924,31.052,30.697,33.178,35.276,34.990,39.280,
		41.069,42.927,45.458,45.716,46.472,45.921,48.327,48.968,48.001,41.154,38.192,32.700,30.121,29.108}/*,
	  //! H2O
		{0.699,1.971,4.069,6.121,7.520,8.934,8.113,9.907,9.930,11.350,13.004,12.734,16.032,18.083,18.897,20.047,21.159,21.908,21.857,22.446,22.487,22.502,22.852,
		22.498,22.118,19.384,20.992,16.975,18.151,16.632,19.837,20.512,15.072,15.176,18.069,15.271,8.001}*/
	};


	//! solar flux at particular day
	//! T9 on 26.12.2005, Wakeflyby	
	D_REAL F107=89.5; //! daily flux
	D_REAL F107A=84.88; //! averaged flux +-1 month
	D_REAL dist_to_sun=9.106;  //!AE
		
	//! constants
	//! sigma_i * 10^-22 for m^2	
	//! consider both, multiple neutral species and
	//! multiple ion species created from one neutral species
	//! (e.g. H2O -> H2O+, OH+, O+, H+)
	//! array entries are addressed as
	//! [neutral_spec*num_Particle_Species + species]
	const D_REAL neutralSpec_ion[3/*NUM_NEUTRAL_SPECIES*NUM_PARTICLE_SPECIES*/][37] =
	{ //! H2
		{0.0097,0.0758,0.2009,0.4028,0.5509,0.7454,0.6538,0.8999,0.9041,1.2960,1.7840,1.7420,2.8900,3.7780,4.0470,5.2540,6.0500,6.9000,6.7410,7.6680,8.2990,8.2880,
		9.7020,10.7310,9.7610,8.6240,7.0710,5.0720,6.6290,0.0899,0.,0.,0.,0.,0.,0.,0.} ,
	//! N2
		{0.443,1.479,3.153,5.226,6.781,8.100,7.347,9.180,9.210,11.600,15.350,14.669,20.692,22.100,22.772,24.468,24.130,22.400,22.787,22.790,23.370,23.339,29.235,
		25.480,15.060,65.800,8.500,8.860,14.274,0.,0.,0.,0.,0.,0.,0.,0.},
	//! CH4
		{0.051,0.147,0.387,0.839,1.192,1.681,1.398,2.095,2.103,2.957,3.972,3.820,6.255,8.442,9.837,11.432,13.398,14.801,14.640,15.734,17.102,16.883,19.261,20.222,
		21.314,22.599,22.763,23.198,22.886,25.607,24.233,13.863,0.136,0.475,0.,0.,0.}
	};
	
	//! constants for solar flux
	//! length of wave lengths intervall (1 means single line)
	double dlambda[37]={50, 50, 50, 50, 1, 1, 50, 1, 1, 50, 1, 50, 50, 1, 50, 50, 1, 1, 50,1, 1, 50, 50, 1, 50, 1, 1, 1, 50, 50, 50, 50, 1, 50, 1, 1, 50};
	//! I_Inf_i * 10^13 for SI  ##(F74133 reference flux)
	double I_inf_i[37]={1.2,0.45,4.8,3.1,0.46,0.21,1.679,0.8,6.9,0.965,0.650,0.314,0.383,0.29,0.285,0.452,0.72,1.27,0.357,0.53,1.59,0.342,0.32,0.36,0.141,0.17,0.26,0.702,0.758,1.625,3.537
		,3,4.4,1.475,3.5,2.1,2.467};
	double A[37]={1.0017e-2,7.1250e-3,1.3375e-2,1.9450e-2,2.7750e-3,1.3768e-1,2.6467e-2,2.5000e-2,3.3333e-3,2.2450e-2,6.5917e-3,3.6542e-2,7.4082e-3,7.4917e-3,2.0225e-2,8.7583e-3,3.2667e-3,
		5.1583e-3,3.6583e-3,1.6175e-2,3.3250e-3,1.18e-2,4.2667e-3,3.0417e-3,4.7500e-3,3.8500e-3,1.2808e-2,3.2750e-3,4.7667e-3,4.8167e-3,5.6750e-3,4.9833e-3,3.9417e-3,4.4168e-3,
		5.1833e-3,5.2833e-3,4.3750e-3};	

		
	//! calculate solar flux at particular day and distance to sun
	D_REAL P=(F107+F107A)/2.;
	for(int i=0;i<37;i++)
	I_inf_i[i]*=(1+A[i]*(P-80))/(dist_to_sun*dist_to_sun);   
	
	INT32 num_Nds_x = pow(2.,(double)MAX_LEVEL)*(BlkNds_X-2)*RB_X;
	INT32 num_Nds_y = pow(2.,(double)MAX_LEVEL)*(BlkNds_Y-2)*RB_Y;
	INT32 num_Nds_z = pow(2.,(double)MAX_LEVEL)*(BlkNds_Z-2)*RB_Z;
	
	
	//! loop over all coordinates
	for (ind[0]=0; ind[0] < num_Nds_x; ind[0]++)
	 for (ind[1]=0; ind[1] < num_Nds_y; ind[1]++)
	  for (ind[2]=0; ind[2] < num_Nds_z; ind[2]++)
	  {
		  
		coor_of_node[0] = - Box_Origin[0]  + ind[0] * float(1.*LX/(1.*num_Nds_x)) ;
		coor_of_node[1] = - Box_Origin[1]  + ind[1] * float(1.*LY/(1.*num_Nds_y)) ;
		coor_of_node[2] = - Box_Origin[2]  + ind[2] * float(1.*LZ/(1.*num_Nds_z)) ;
		 
		//! TODO theoretically, one could at first check if node exists at this position...			
		
		//! get positions along line of sight path from coor_of_node
		calc_lineOfSight_Positions(positions, num_positions, coor_of_node,dx, vector_to_sun);  	
	
		
		for(INT32 neutral_spec=0; neutral_spec<num_Neutral_Species; neutral_spec++)
		{	
		
			//! integrate (add) neutral density along line of sight
			ColumnDens = 0;
			calc_lineOfSight_ColumnDens(ColumnDens, positions, num_positions, neutral_spec,dx);


			//! sum up value from all processes
			INT32 num_values = 1;
			D_REAL ColumnDens_local[1] = {ColumnDens};
			D_REAL ColumnDens_global[1] = {0};				
			mpi_build_sum(ColumnDens_local,ColumnDens_global,num_values);					
			
			//! set integrated neutral density to ion production field
			set_value_to_field(ColumnDens_global,coor_of_node,id_density_ionProdSpecies1+neutral_spec*num_Particle_Species);
		
		}
	
	}
	
	//! for each neutral species: copy integrated neutral density to ion production fields of the different ion species
	 for(INT32 neutral_spec=0; neutral_spec<num_Neutral_Species; neutral_spec++)
	  for(INT32 species=0; species<num_Particle_Species; species++)  
	  copy_Field(id_density_ionProdSpecies1+neutral_spec*num_Particle_Species+species,id_density_ionProdSpecies1+neutral_spec*num_Particle_Species);
	 
	
	 //! ghost nodes
	 for(INT32 neutral_spec=0; neutral_spec<num_Neutral_Species; neutral_spec++)
	  for(INT32 species=0; species<num_Particle_Species; species++)  	 
	   FULL_GN_UPDATE(id_density_ionProdSpecies1 +neutral_spec*num_Particle_Species+species);
	 
	synchronize_allProcesses_MPI();
	
	ColumnDens = 0;
	
	//! now integrated neutral density along line of sight is known at all nodes
	//! -> calculate photoabsorption	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank==temp_Block->responsible_mpi_process)
			{
				ion_prod_allSpec = temp_Block->Field_Type[id_density_ionProdSpecies1];
				neutral_dens_allSpec = temp_Block->Field_Type[id_numberdensity_neutralSpecies1];		
				
				
				for(INT32 node=0; node<num_nodes_in_block; node++)
				{
					
					//! sum up line of sight absorption of all neutral species
					//! (separately for each wavelength)
					D_REAL Absorption_all_Spec[37];			
					memset(Absorption_all_Spec,0,37*sizeof(D_REAL));
					
					D_REAL Ion_Production[NUM_NEUTRAL_SPECIES];
					memset(Ion_Production,0,NUM_NEUTRAL_SPECIES*sizeof(D_REAL));

					
					for(INT32 i=0;i<37;i++)
					{
						
						//! sum over all neutral species since absorption at lambda given by
						//! exp(- sum_neutral sigma_abs(lambda) * NeutralColumnDens )
						for(INT32 neutral_spec=0; neutral_spec<num_Neutral_Species; neutral_spec++)
						{ 
							//! ion production field is still neutral column density (at ion species 0)
							ion_prod = ion_prod_allSpec +(neutral_spec*num_Particle_Species)*num_nodes_in_block;
							
							//! store neutral column density in temporary variable
							//! for i=0 to allow overwriting it but use it for all other bins
							if(i==0)
							Ion_Production[neutral_spec] = ion_prod[node];
							
							Absorption_all_Spec[i] += 1.e-22*SI_n0*neutralSpec_abs[neutral_spec][i]*Ion_Production[neutral_spec];
						}	
												
						//! now set ion production for every ion species from every neutral species	
						for(INT32 neutral_spec=0; neutral_spec<num_Neutral_Species; neutral_spec++)	
						{	
							neutral_dens = neutral_dens_allSpec +neutral_spec*num_nodes_in_block;
							
							for(INT32 species=0; species<num_Particle_Species; species++)  
							{
								ion_prod = ion_prod_allSpec +(neutral_spec*num_Particle_Species + species)*num_nodes_in_block;
											
								//! ion production at node is stored in Ion_Production[neutral_spec] for i=0
								//! thus, it could be set to zero to allow summation over wavelengths
								if(i==0)
								ion_prod[node] = 0;
								
								ion_prod[node] += dlambda[i] * neutralSpec_ion[neutral_spec*num_Particle_Species + species][i]*1.e-22
										*I_inf_i[i]*1.e+13 * neutral_dens[node]*exp(-Absorption_all_Spec[i]);
								
							}	 
						}

					
					} //! end wavelength bins

				} //! end nodes	

			} //! end responsible_mpi_process	
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}
	

	log_file << "done." << endl;

}




void CHybrid::calc_lineOfSight_Positions(D_REAL** &positions, INT32 &num_positions, D_REAL *coor_of_node, D_REAL &dx, D_REAL *vector_to_sun)
{

	//! pointing AWAY from box
	D_REAL lineOfSight_Vec[3] = { vector_to_sun[0], vector_to_sun[1] , vector_to_sun[2]};

	
	//! Stepsize
	//! 1/2 * Smallest Box Size 
	D_REAL delta_r[3];	
	delta_r[0] = LX/(2.*pow(2.,(double)MAX_LEVEL)*(BlkNds_X-2)*RB_X);
 	delta_r[1] = LY/(2.*pow(2.,(double)MAX_LEVEL)*(BlkNds_Y-2)*RB_Y);
 	delta_r[2] = LZ/(2.*pow(2.,(double)MAX_LEVEL)*(BlkNds_Z-2)*RB_Z); 
	
	//! set dx to min of delta_r
	if(delta_r[0]<delta_r[1])
	{
		if(delta_r[0]<delta_r[2])
		dx=delta_r[0];
		else
		dx=delta_r[2];	

	}
	else
	{
		if(delta_r[1]<delta_r[2])
		dx=delta_r[1];
		else
		dx=delta_r[2];
		
	}	
	
	
	//! estimate number of positions
	bool in_box = true;
	D_REAL tempx[3];
	INT32 i=0;
	while(in_box==true)
	{
		tempx[0] = coor_of_node[0] + lineOfSight_Vec[0]*dx*i;
		tempx[1] = coor_of_node[1] + lineOfSight_Vec[1]*dx*i;
		tempx[2] = coor_of_node[2] + lineOfSight_Vec[2]*dx*i;	
	
		//! normalize coordinate to box length
		//! get normalized position in Simu-Box
		tempx[0] = (tempx[0]+Box_Origin[0])/LX;
		tempx[1] = (tempx[1]+Box_Origin[1])/LY;
		tempx[2] = (tempx[2]+Box_Origin[2])/LZ;
	
		i++;
		
		if(    tempx[0]<0. || tempx[0] >= 1.
		    || tempx[1]<0. || tempx[1] >= 1.
		    || tempx[2]<0. || tempx[2] >= 1.)
		{
			in_box = false;
		}		
	}	
	

	
	num_positions = i;
	
	//! alloc coordinates memory
	positions = new D_REAL*[3];
	positions[0] = new D_REAL[num_positions];
	positions[1] = new D_REAL[num_positions]; 
	positions[2] = new D_REAL[num_positions];
	

	INT32 counter = 0;
	
	while(counter<num_positions)
	{
	
		positions[0][counter] = coor_of_node[0] + lineOfSight_Vec[0]*dx*counter;
		positions[1][counter] = coor_of_node[1] + lineOfSight_Vec[1]*dx*counter;
		positions[2][counter] = coor_of_node[2] + lineOfSight_Vec[2]*dx*counter;
		
		//! normalize coordinate to box length
		//! get normalized position in Simu-Box
		positions[0][counter] = (positions[0][counter]+Box_Origin[0])/LX;
		positions[1][counter] = (positions[1][counter]+Box_Origin[1])/LY;
		positions[2][counter] = (positions[2][counter]+Box_Origin[2])/LZ;
		
		counter++;		
	
	}
	
}	

//!------------------------------------------------------------------------
//! add neutral density from each position to ColumnDens
//! function is a slightly adapted version of trace_trajectory
//!------------------------------------------------------------------------
void CHybrid::calc_lineOfSight_ColumnDens(D_REAL &ColumnDens, D_REAL** &positions, INT32 &num_positions, INT32 neutral_spec, D_REAL &dx)
{
	
	//! indices for interpolation
	INT32  i,j,k;
	INT32  i_j_k, ip1_jp1_kp1;
	INT32  ip1_j_k, i_jp1_k,  i_j_kp1;
	INT32  ip1_jp1_k, ip1_j_kp1,  i_jp1_kp1;

	ofstream trajectory_outfile;
	INT32  *index_blk[3], *index_cell[3], blk_index;
	D_REAL *x_cell[3], r_vec[3],  shape_func[8], *BlkField[3];

	CBlock *temp_Block, *top_level_Block;



	//!-----------------------------------
	//! 1) alloc variable memory
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		 index_blk[comp] = new INT32[MAX_LEVEL+1];
		index_cell[comp] = new INT32[MAX_LEVEL+1];
		    x_cell[comp] = new D_REAL[MAX_LEVEL+1];
	}


	//!-----------------------------------
	//! 2) trace
	//!-----------------------------------

	for(INT32 counter=0; counter<num_positions; counter++)
	{
	
		//!--------------------------------------------------
		//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
		//!     (even though target blk might not be refined)
		//!--------------------------------------------------
		r_vec[0] = positions[0][counter];
		r_vec[1] = positions[1][counter];
		r_vec[2] = positions[2][counter];


		//! calculate indices and check whether position is in Simu Box
		if(calc_BlockNodeXCell_of_Pos(index_blk, index_cell, x_cell, r_vec))
		{


			//!--------------------------------------------------
			//! 2b) climb to highest possible level at respective
			//!    position.
			//!--------------------------------------------------
		
	
			//! copy indices to temp array
			INT32 blk_indices[3] = {index_blk[0][0],
					       index_blk[1][0],
					       index_blk[2][0]};

			
			//! first calculate root block
			//! -> always level 0
			if(use_SFC)
			blk_index = SFC_Indices_to_BlkNr(blk_indices, 0);
			else
			blk_index = LINEAR_Indices_to_BlkNr(blk_indices);
		
		
			//! initialize to respective root block
			temp_Block      = Root_Block_Array +blk_index;
			top_level_Block = Root_Block_Array +blk_index;;
		
			//! descent to highest existing level
			for(INT32 level=1; level<=MAX_LEVEL; level++)
			{
		
				//! set block index to precalculated indices
				blk_index = index_blk[0][level] *2 *2
					   +index_blk[1][level] *2
					   +index_blk[2][level];
				
		
				//! Decent by ordanary child array
				temp_Block = temp_Block->child_array[blk_index];
		
				//! cancel in case block does not exist
				if(!temp_Block)
				break;
		
				//! Set top_level_Block to temp_Block
				top_level_Block = temp_Block;
		
				
			}
	

			//!-----------------------------------------------------
			//! 2c) Now target Block and Cell are known, interpolate
			//!     Field at respective position.
			//!-----------------------------------------------------
	
			if(mpi_myRank == top_level_Block->responsible_mpi_process)
			{
				//! use pointers to block's field 
// 				for(INT32 comp=0; comp<COMPs_FType[Field_id]; comp++)
// 				BlkField[comp] = top_level_Block->Field_Type[Field_id] +comp *num_nodes_in_block;
				
				BlkField[0] = top_level_Block->Field_Type[id_numberdensity_neutralSpecies1+neutral_spec];

		
				INT32 level=top_level_Block->RLevel;
		
				//! use common i,j,k notation:
				i = index_cell[0][level];
				j = index_cell[1][level];
				k = index_cell[2][level];
		
				//! ------------------------------------------
				i_j_k   =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
				
				//! ------------------------------------------
				ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
				i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
				i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
				
				//! -------------------------------------------
				ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
				ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
				i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
				
				ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
		
				
				//! use now r_vec for psition in respective cell
				r_vec[0] = x_cell[0][level];
				r_vec[1] = x_cell[1][level];
				r_vec[2] = x_cell[2][level];
		
				shape_func[0] = (1.-r_vec[0])*(1.-r_vec[1])*(1.-r_vec[2]);
				shape_func[1] = (   r_vec[0])*(1.-r_vec[1])*(1.-r_vec[2]);
				shape_func[2] = (1.-r_vec[0])*(   r_vec[1])*(1.-r_vec[2]);
				shape_func[3] = (1.-r_vec[0])*(1.-r_vec[1])*(   r_vec[2]);
				
				shape_func[4] = (   r_vec[0])*(   r_vec[1])*(1.-r_vec[2]);
				shape_func[5] = (   r_vec[0])*(1.-r_vec[1])*(   r_vec[2]);
				shape_func[6] = (1.-r_vec[0])*(   r_vec[1])*(   r_vec[2]);
				shape_func[7] = (   r_vec[0])*(   r_vec[1])*(   r_vec[2]);
	
		
				//! interpolate and store field in variable A
				for(INT32 comp=0; comp<1; comp++)
				{		
				      ColumnDens += dx *
						 ( BlkField[comp][  i_j_k] * shape_func[0]
						+BlkField[comp][ip1_j_k] * shape_func[1]
						+BlkField[comp][i_jp1_k] * shape_func[2]
						+BlkField[comp][i_j_kp1] * shape_func[3]
				
						+BlkField[comp][  ip1_jp1_k] * shape_func[4]
						+BlkField[comp][  ip1_j_kp1] * shape_func[5]
						+BlkField[comp][  i_jp1_kp1] * shape_func[6]
						+BlkField[comp][ip1_jp1_kp1] * shape_func[7] );
		
						 
				}				
			
			}//! end if responsible_mpi_process
		}

	}//! end for positions
	//! NOTE: remove after DFT


	//!-----------------------------------
	//! ) clean up
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		delete index_blk[comp];
		delete index_cell[comp];
		delete x_cell[comp] ;		
	}

	//! delete coordinates memory
	delete[] positions[0];
	delete[] positions[1];
	delete[] positions[2];

	delete[] positions;


		
	
}	




//!------------------------------------------------------------------------
//! add neutral density from each position to ColumnDens
//! function is a slightly adapted version of trace_trajectory
//!------------------------------------------------------------------------
void CHybrid::set_value_to_field(D_REAL* ColumnDens, D_REAL* coor_of_node, INT32 field_id)
{
	
	//! indices for interpolation
	INT32  i,j,k;
	INT32  i_j_k;



	INT32  *index_blk[3], *index_cell[3], blk_index;
	D_REAL *x_cell[3], r_vec[3],  *BlkField[3];

	CBlock *temp_Block, *top_level_Block;



	//!-----------------------------------
	//! 1) alloc variable memory
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		 index_blk[comp] = new INT32[MAX_LEVEL+1];
		index_cell[comp] = new INT32[MAX_LEVEL+1];
		    x_cell[comp] = new D_REAL[MAX_LEVEL+1];
	}


	//!-----------------------------------
	//! 2) trace
	//!-----------------------------------


	//!--------------------------------------------------
	//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
	//!     (even though target blk might not be refined)
	//!--------------------------------------------------

	//! normalize coordinate to box length
	//! get normalized position in Simu-Box
	r_vec[0] = (coor_of_node[0]+Box_Origin[0])/LX;
	r_vec[1] = (coor_of_node[1]+Box_Origin[1])/LY;
	r_vec[2] = (coor_of_node[2]+Box_Origin[2])/LZ;
	
	
	//! calculate indices and check whether position is in Simu Box
	if(calc_BlockNodeXCell_of_Pos(index_blk, index_cell, x_cell, r_vec))
	{

		//!--------------------------------------------------
		//! 2b) climb to highest possible level at respective
		//!    position.
		//!--------------------------------------------------
	

		//! copy indices to temp array
		INT32 blk_indices[3] = {index_blk[0][0],
					index_blk[1][0],
					index_blk[2][0]};

		
		//! first calculate root block
		//! -> always level 0
		if(use_SFC)
		blk_index = SFC_Indices_to_BlkNr(blk_indices, 0);
		else
		blk_index = LINEAR_Indices_to_BlkNr(blk_indices);
	
	
		//! initialize to respective root block
		temp_Block      = Root_Block_Array +blk_index;
		top_level_Block = Root_Block_Array +blk_index;;
	
		//! descent to highest existing level
		for(INT32 level=1; level<=MAX_LEVEL; level++)
		{
	
			//! set block index to precalculated indices
			blk_index = index_blk[0][level] *2 *2
					+index_blk[1][level] *2
					+index_blk[2][level];
			
	
			//! Decent by ordanary child array
			temp_Block = temp_Block->child_array[blk_index];
	
			//! cancel in case block does not exist
			if(!temp_Block)
			break;
	
			//! Set top_level_Block to temp_Block
			top_level_Block = temp_Block;
	
			
		}


		//!-----------------------------------------------------
		//! 2c) Now target Block and Cell are known, interpolate
		//!     Field at respective position.
		//!-----------------------------------------------------

		if(mpi_myRank == top_level_Block->responsible_mpi_process)
		{
			//! use pointers to block's field 
// 				for(INT32 comp=0; comp<COMPs_FType[Field_id]; comp++)
// 				BlkField[comp] = top_level_Block->Field_Type[Field_id] +comp *num_nodes_in_block;
			
			BlkField[0] = top_level_Block->Field_Type[field_id];

	
			INT32 level=top_level_Block->RLevel;
	
			//! use common i,j,k notation:
			i = index_cell[0][level];
			j = index_cell[1][level];
			k = index_cell[2][level];
	
			//! ------------------------------------------
			i_j_k   =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;			
			
	
			BlkField[0][i_j_k] = ColumnDens[0];
			

			
		}//! end if responsible_mpi_process


	}//! end for positions
	//! NOTE: remove after DFT

	//!-----------------------------------
	//! ) clean up
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		delete index_blk[comp];
		delete index_cell[comp];
		delete x_cell[comp] ;
	}

	
}	

#endif