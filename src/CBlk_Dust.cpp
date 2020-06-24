


#include "CBlk.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"



#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>


using namespace std;

#ifdef use_dust_species_as_field



extern D_REAL **delta_of_L;

extern D_REAL *CellVol_of_L;
extern D_REAL **Blk_Length_of;

extern INT32 i_j_k;

extern D_REAL *q_of, *q2m_of, *q_m_ratio;

//!-------------------------------------------------------------//
//! set_RhoUi_extern: 								//
//!-------------------------------------------------------------//
void CBlock::set_Dust_extern(INT32 id_densityfield, INT32 id_velocityfield, short int* num_Nodes, FILE_REAL* Origin, FILE_REAL* Length, FILE_REAL* rho, INT32& num_values_not_in_extern_box)
{



	INT32  ind[3];
	INT32  num_extern_box_nodes;
	INT32 a,b,c, a_b_c;
	INT32 ap1_b_c, a_bp1_c, a_b_cp1;
	INT32 ap1_bp1_c, ap1_b_cp1, a_bp1_cp1, ap1_bp1_cp1;

	D_REAL *RHO, *UiX, *UiY, *UiZ;
	FILE_REAL *extern_UiX, *extern_UiY, *extern_UiZ;
	D_REAL r[3], extern_delta[3], shape_func[8];

	PARTICLE_REAL x_BlockNode[3], x_extern[3], cell_intern_r[3];
	memset(cell_intern_r,0,3*sizeof(PARTICLE_REAL));


	for(INT32 comp=0; comp<3; comp++)
	extern_delta[comp] = Length[comp]/(num_Nodes[comp]-1);


	//! Set pointer to extern field
	num_extern_box_nodes =   num_Nodes[0]
				*num_Nodes[1]
				*num_Nodes[2];

	extern_UiX = rho +num_extern_box_nodes;
	extern_UiY = extern_UiX +num_extern_box_nodes;
	extern_UiZ = extern_UiY +num_extern_box_nodes;
	
	//! Set pointer to intern (Block) field
	RHO = Field_Type[id_densityfield];

	UiX = Field_Type[id_velocityfield];
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
		x_extern[0] = (x_BlockNode[0] +Origin[0])/extern_delta[0];
		x_extern[1] = (x_BlockNode[1] +Origin[1])/extern_delta[1];
		x_extern[2] = (x_BlockNode[2] +Origin[2])/extern_delta[2];

		 a  = int(x_extern[0]);
		 b  = int(x_extern[1]);
		 c  = int(x_extern[2]);


		if(   (x_extern[0]<0) || a>num_Nodes[0]-2
		    ||  x_extern[1]<0 || b>num_Nodes[1]-2
		    ||  x_extern[2]<0 || c>num_Nodes[2]-2)
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
			a_b_c   =      a*num_Nodes[1]*num_Nodes[2]    +b*num_Nodes[2]      +c;
			
			//! -----------------------------------------------
			ap1_b_c = (a+1)*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2]     +c;
			a_bp1_c =     a*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2]     +c;
			a_b_cp1 =     a*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2] +(c+1);
			
			//! ------------------------------------------------
			ap1_bp1_c = (a+1)*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2]     +c;
			ap1_b_cp1 = (a+1)*num_Nodes[1]*num_Nodes[2]     +b*num_Nodes[2] +(c+1);
			a_bp1_cp1 =     a*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2] +(c+1);
			
			ap1_bp1_cp1 = (a+1)*num_Nodes[1]*num_Nodes[2] +(b+1)*num_Nodes[2] +(c+1);

			RHO[i_j_k]  =  rho[  a_b_c] * shape_func[0]
				      +rho[ap1_b_c] * shape_func[1]
				      +rho[a_bp1_c] * shape_func[2]
				      +rho[a_b_cp1] * shape_func[3]
			
				      +rho[  ap1_bp1_c] * shape_func[4]
				      +rho[  ap1_b_cp1] * shape_func[5]
				      +rho[  a_bp1_cp1] * shape_func[6]
				      +rho[ap1_bp1_cp1] * shape_func[7];

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
			
				      
			//! catch errors in dust velocity	      
			if(fabs(UiX[i_j_k])>V_sw[0]*3.)
			{	
				if(UiX[i_j_k]>0)
				UiX[i_j_k] = V_sw[0];
				else
				UiX[i_j_k] = -V_sw[0];
			}
			
			if(fabs(UiY[i_j_k])>V_sw[0]*3.)
			{	
				if(UiY[i_j_k]>0)
				UiY[i_j_k] = V_sw[0];
				else
				UiY[i_j_k] = -V_sw[0];
			}
				      
			if(fabs(UiZ[i_j_k])>V_sw[0]*3.)
			{	
				if(UiZ[i_j_k]>0)
				UiZ[i_j_k] = V_sw[0];
				else
				UiZ[i_j_k] = -V_sw[0];
			}
     
 
			if(Flag[i_j_k])
			{
				RHO[i_j_k] = 0;
				UiX[i_j_k] = 0;
				UiY[i_j_k] = 0;
				UiZ[i_j_k] = 0;
			}	

		}
	  }

}
 
void CBlock::add_dust_to_field(INT32 id_rho)
{
	D_REAL *dust_rho_extern, *rho; 

	for(short dustSpec=0; dustSpec<num_Dust_Species; dustSpec++)
	{	
			
		//! total ion density
		dust_rho_extern = Field_Type[id_density_dustSpecies1+dustSpec];					

		//! set pointer to respective species
		rho = Field_Type[id_rho];					

	
		for(INT32 node=0; node<num_nodes_in_block; node++)
		{
			if(TL<TL_DUST_MAX)
				rho[node] += dust_rho_extern[node]*TL/TL_DUST_MAX;
			else
				rho[node] += dust_rho_extern[node];
	
		}
			
	}
}	


void CBlock::add_dust_velocity_to_field(INT32 Ioncurrent_id)
{
	D_REAL *dust_rho_extern, *dust_u_extern_x, *Ji_x; 

	for(short dustSpec=0; dustSpec<num_Dust_Species; dustSpec++)
	{	
			
		//! total ion density
		dust_rho_extern = Field_Type[id_density_dustSpecies1+dustSpec];					

		dust_u_extern_x = Field_Type[id_velocity_dustSpecies1+dustSpec];
		D_REAL* dust_u_extern_y = dust_u_extern_x + 1*num_nodes_in_block;
		D_REAL* dust_u_extern_z = dust_u_extern_x + 2*num_nodes_in_block;
		
		//! set pointers to Ji
		Ji_x = Field_Type[Ioncurrent_id];
		D_REAL* Ji_y = Ji_x + 1*num_nodes_in_block;
		D_REAL* Ji_z = Ji_x + 2*num_nodes_in_block;					

	
		for(INT32 node=0; node<num_nodes_in_block; node++)
		{
			if(TL<TL_DUST_MAX)
			{	
				Ji_x[node] += dust_rho_extern[node]*dust_u_extern_x[node]*TL/TL_DUST_MAX;
				Ji_y[node] += dust_rho_extern[node]*dust_u_extern_y[node]*TL/TL_DUST_MAX;
				Ji_z[node] += dust_rho_extern[node]*dust_u_extern_z[node]*TL/TL_DUST_MAX;
			}	
			else
			{	
				Ji_x[node] += dust_rho_extern[node]*dust_u_extern_x[node];
				Ji_y[node] += dust_rho_extern[node]*dust_u_extern_y[node];
				Ji_z[node] += dust_rho_extern[node]*dust_u_extern_z[node];
			}	
	
		}
			
	}
}	



//!-------------------------------------------------------------//
//! for neutral density					//
//!-------------------------------------------------------------//
void CBlock::set_analytical_dust_field(INT32 id_dustSpecDens, INT32 id_dustSpecVel)
{

	const D_REAL H_theta_dust = 1./((opening_angle_dust*M_PI/180.)*(opening_angle_dust*M_PI/180.));

	INT32  ind[3];

	D_REAL *RHO, *UiX, *UiY, *UiZ;

	D_REAL x_BlockNode[3], cell_intern_r[3];
	memset(cell_intern_r,0,3*sizeof(D_REAL));

	//! Set pointer to intern (Block) field
	RHO = Field_Type[id_dustSpecDens];

	UiX = Field_Type[id_dustSpecVel];
	UiY = UiX +num_nodes_in_block;
	UiZ = UiY +num_nodes_in_block;

	//! use a,b,c for extern mesh 
	//! use i,j,k for intern mesh
	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
	
	
		INT32 i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
		//! get Coordinate of intern Block Node 
		intern2normedCoords(x_BlockNode, cell_intern_r, ind);		
	
		if(!use_test_scenario)
		{	
			
			D_REAL absr = vec_len(x_BlockNode);

			D_REAL temp_X[3];
			intern2normedCoords(temp_X, cell_intern_r, ind);

			//! NOTE MODIFICATION FOR ASYMMETRY
			
	// 		if(x_BlockNode[0]>0 && x_BlockNode[0] < x_down1 && x_BlockNode[2]<-1*RE)
	// 		x_BlockNode[0] *= fac_downstream;
	// 		
	// 		if(x_BlockNode[0]>0 && x_BlockNode[0] < x_down1 && x_BlockNode[2]>=-1*RE && x_BlockNode[2] < 1*RE)
	// 		x_BlockNode[0] *= (1.-fac_downstream)/(2.*RE)* x_BlockNode[2] + (1.+fac_downstream)/2.;  
	// 
	// 		D_REAL	fac_downstream_x = (1.-fac_downstream)/(x_down2 - x_down1)*temp_X[0] + 1. - (1.-fac_downstream)/(x_down2 - x_down1) * x_down2;
	// 
	// 		if(x_BlockNode[0]> x_down1 && x_BlockNode[0] < x_down2 && x_BlockNode[2]<-1*RE)
	// 		x_BlockNode[0] *= fac_downstream_x;
	// 		
	// 		if(x_BlockNode[0]> x_down1 && x_BlockNode[0] < x_down2 && x_BlockNode[2]>=-1*RE && x_BlockNode[2] < 1*RE)
	// 		x_BlockNode[0] *= (1.-fac_downstream_x)/(2.*RE)* x_BlockNode[2] + (1.+fac_downstream_x)/2.;  		
	// 
	//                if(x_BlockNode[0]<0 && x_BlockNode[2]<-1*RE)
	//                 x_BlockNode[0] *= fac_upstream;
	// 
	//                 if(x_BlockNode[0]<0 && x_BlockNode[2]>=-1*RE && x_BlockNode[2] < 1*RE)
	//                 x_BlockNode[0] *= (1.-fac_upstream)/(2.*RE)* x_BlockNode[2] + (1.+fac_upstream)/2.;

			if(absr*absr > R_Moon*R_Moon )
			{
					
				D_REAL nneutral=0;	    
			
	// 			D_REAL axis[3] = {0,0,1};
	// 		       
				D_REAL asym = 0;
	// 			
	// 			//! NOTE: BE CAREFUL ABOUT LENGTHS IN RE AND IN l0 !!!!!!!!!!!!!!!!!!!!!!!!!!!
	// 			D_REAL pos_wrt_source[3] ={x_BlockNode[0]-RE*axis[0]+plume_verbeulung*(RE + x_BlockNode[2]*RE),
	// 						x_BlockNode[1]-RE*axis[1],x_BlockNode[2]-RE*axis[2]};
	// 			
	// 			D_REAL betrag = vec_len(pos_wrt_source);
	// 		
	// 			D_REAL theta = acos( (axis[0]*pos_wrt_source[0] + axis[1]*pos_wrt_source[1] + axis[2]*pos_wrt_source[2])/(1.*betrag));
	// 			
	// 			nneutral = xi_nD_0/(1.*SI_n0)*RE*RE/(betrag*betrag)*exp(-(theta)*(theta)*H_theta_dust)*exp(-(betrag-RE)/H_d_dust);
				
			
	//          
		
	// 			for(short source=0;source<8;source++)
	// 			{
				
				
				
					//! positions of jets from spitale & porco in spherical coordinates
	// 				double jets[9][2]={	  { 172*2*M_PI/360 ,  57*2*M_PI/360 },
	// 							{ 169*2*M_PI/360 , 135*2*M_PI/360 },
	// 							{ 171*2*M_PI/360 , 157*2*M_PI/360 },
	// 							{ 163*2*M_PI/360 , 301*2*M_PI/360 },
	// 							{ 170*2*M_PI/360 , 290*2*M_PI/360 }, //! alternative source IV from Dong 2011
	// 							{ 169*2*M_PI/360 ,  18*2*M_PI/360 },
	// 							{ 177*2*M_PI/360 , 219*2*M_PI/360 },
	// 							{ 165*2*M_PI/360 ,  60*2*M_PI/360 },
	// 							{ 172*2*M_PI/360 , 334*2*M_PI/360 } };


	// 				D_REAL axis[3] = {cos(jets[source][1])*sin(jets[source][0]),sin(jets[source][1])*sin(jets[source][0]),cos(jets[source][0])};
					D_REAL axis[3] = {0,0,-1};
		
					//! NOTE: BE CAREFUL ABOUT LENGTHS IN RE AND IN l0 !!!!!!!!!!!!!!!!!!!!!!!!!!!
					D_REAL pos_wrt_source[3] ={x_BlockNode[0]-R_Moon*axis[0]+ asym*(x_BlockNode[2]-R_Moon*axis[2]),x_BlockNode[1]-R_Moon*axis[1],x_BlockNode[2]-R_Moon*axis[2]};
					
					D_REAL betrag = vec_len(pos_wrt_source);
				
					D_REAL theta = acos( (axis[0]*pos_wrt_source[0] + axis[1]*pos_wrt_source[1] + axis[2]*pos_wrt_source[2])/(1.*betrag));
		
			
					//! neutral plume from Saur 2008 
					//! ----------------------------------------------------------------------------------------------------------------
					nneutral += xi_nD_0/(1.*SI_n0)*R_Moon/(absr)*exp(-(theta)*(theta)*H_theta_dust)*exp(-(absr-R_Moon)/H_d_dust);
					//! ----------------------------------------------------------------------------------------------------------------

					//! pickup modification
	
					D_REAL phi = atan2(x_BlockNode[1],x_BlockNode[0]) + M_PI;

					D_REAL pickupmod = 3*cos(0.5*(phi+M_PI/4.))*cos(0.5*(phi+M_PI/4.));   
					 
// 					log_file << " phi = " << phi << "    pickup: "<<pickupmod<<endl;  
					
					nneutral *= 1 + pickupmod;


					//! neutral plume from Dong 2011 
					//! ----------------------------------------------------------------------------------------------------------------	
				
					//! for E3
	// 				D_REAL Mach = 1.6;
	// 				D_REAL SourceStrength_SI[9] = {1.3e+22,1.3e+22,1.3e+22,0,0,0,1.3e+22,0,0} ;
					//!NOTE: in 1/cm

					//! for E5
	// 				D_REAL Mach = 1.4;
	// 				D_REAL SourceStrength_SI[9] = {3.2e+22,3.2e+22,3.2e+22,0,0,0,3.2e+22,0,0} ;
					//!NOTE: in 1/cm
					
					//! for E7 
	// 				D_REAL Mach = 3;
	// 				D_REAL SourceStrength_SI[9] = {1,1,1,3.5,0.e+22,0,0,0,0} ;
	// //				D_REAL SourceStrength_SI[9] = {1.e+22,1.e+22,1.e+22,3.5e+22,0.e+22,0,0,0,0} ;  //! Dong 2011
	// //				D_REAL SourceStrength_SI[9] = {1.1e+22,1.1e+22,1.1e+22,0.e+22,1.e+22,0,0,0,0} ;//! NOTE: Best E7 from Dong2011
	// 				//!NOTE: in 1/cm
	// 		 				 				 		
	//  				nneutral += xi_nD_0*(15.e+5*15.e+5)/(M_PI*betrag*betrag*l0_SI*l0_SI*1.e+4)*(
	//  								    2.*Mach*cos(theta)/(1.*sqrt(M_PI))*exp(-(Mach*Mach)) 
	//  								     + exp(-(Mach*Mach*sin(theta)*sin(theta)))*(1.+2.*Mach*Mach*cos(theta)*cos(theta))
	//  									*(1.+erf(Mach*cos(theta))))*exp(-(betrag-RE)/H_d_dust);
	//				nneutral += xi_nD_0*(15.e+5*15.e+5)/(M_PI*betrag*betrag*l0_SI*l0_SI*1.e+4)*betrag/RE*(
	//								    2.*Mach*cos(theta)/(1.*sqrt(M_PI))*exp(-(Mach*Mach)) 
	//								     + exp(-(Mach*Mach*sin(theta)*sin(theta)))*(1.+2.*Mach*Mach*cos(theta)*cos(theta))
	//									*(1.+erf(Mach*cos(theta))))*exp(-(betrag-RE)/H_d_dust);

					//! ----------------------------------------------------------------------------------------------------------------
									
	// 			}
					

	//                 	 D_REAL n0 = 7.4e+05/70.5;
	//         	         D_REAL nbg = 3.3e+04/70.5;
		
				RHO[i_j_k] += nneutral ;//+  RE*RE/(absr*absr)*n0 + nbg;

				//! velocity
				
	// 			D_REAL vmax = 1./3.*V_sw[0];
	// 			D_REAL zmin = -7.5*RE;
				
				UiX[i_j_k] = 0;
				
	// 			if(x_BlockNode[2]>-RE)
				UiY[i_j_k] = 0;
				
	// 			if(x_BlockNode[2]<zmin)
	// 			UiY[i_j_k] = vmax;
					
	// 			if(x_BlockNode[2]>zmin && x_BlockNode[2]<-RE)
	// 			UiY[i_j_k] = vmax/(zmin + RE)*x_BlockNode[2] + vmax*RE/(zmin + RE);
				
				UiZ[i_j_k] = 0;
			}
	// 		else
	// 		RHO[i_j_k] = 0;
		}
		else{ //! test scenario
			
			
			if(x_BlockNode[0]*x_BlockNode[0]+x_BlockNode[1]*x_BlockNode[1]<4*R_Moon*4*R_Moon && fabs(x_BlockNode[2])<2*R_Moon)
			{ 
		
				if(!Flag[i_j_k])
				{	
					RHO[i_j_k] = 2./3.;
					UiX[i_j_k] = 0;
					UiY[i_j_k] = 0;
					UiZ[i_j_k] = 0;
				}
			}
			else
			RHO[i_j_k] = 0;

		}	

		
	  }
	  
}	

#endif
