



#include "CBlk.h"
#include "utils.h"
#include "parameters.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>






//! Some general comments to field solver:
//! Vorschlag Schuele:
//! Neben einer Refinement grenze eine Puffer Zone einführen,
//! das ist ein ebenfalls refinter Block der die Werte von
//! seinem parent bekommt aber keine Werte auf diesen injeziert


using namespace std;

extern INT32  *COMPs_FType;
extern D_REAL ***smooth_coeff;
extern D_REAL *q_of;



extern D_REAL *dt_field_of_L, *dt_particle_of_L, **rd_of_L, **rd2_of_L;

//! pointers for EM Fields
extern D_REAL *BX, *BY, *BZ;
extern D_REAL *EX, *EY, *EZ;
extern D_REAL *UX, *UY, *UZ;
extern D_REAL *eta;
//! pointers to calc derivatives
extern D_REAL *rd;


//! 
D_REAL *rho_n, *rho_total;
	
//! use id_rotB as temp scratch for Ji_plus
D_REAL* JiX_plus, *JiY_plus,  *JiZ_plus;

D_REAL* UiX_plus;
D_REAL* UiY_plus;
D_REAL* UiZ_plus = UiY_plus +num_nodes_in_block;

//! CAM pointers:
//! use id_UI_minus as temp scratch for Ji_np1
D_REAL *JiX_np1, *JiY_np1, *JiZ_np1;
D_REAL *Lam, *GamX, *GamY, *GamZ;

D_REAL *rotBX, *rotBY, *rotBZ;
D_REAL *BX_total, *BY_total, *BZ_total;

//! MPI related variables
extern ofstream log_file;






//!---------------------------------------------------------------------
//!------ Calculation of BField PDE, EField, rotB, divB, etc...---------
//!---------------------------------------------------------------------


//!-------------------------------------------------------------//
//! show_PE: -							//
//!-------------------------------------------------------------//
void CBlock::show_PE(INT32 species)
{


	D_REAL* PE  = Field_Type[id_PESpecies1  + species];
	D_REAL* rho = Field_Type[id_rhoSpecies1 + species];
// 	D_REAL* rho = Field_Type[id_allRhoSpecies +species];
	
#if defined nonadiabatic_gradPE_TERM
	//  nothing to do because PE is already computed
	return;
#else
        for(INT32 node=0; node<num_nodes_in_block; node++)
	PE[node] = Electron_Betas[species] *pow(rho[node],kappa_electron);
#endif
	
}

void CBlock::show_gradPE(INT32 field_id)
{
	
	D_REAL* gradPE_X = Field_Type[field_id];
	D_REAL* gradPE_Y = gradPE_X +num_nodes_in_block;
	D_REAL* gradPE_Z = gradPE_Y +num_nodes_in_block;
	
	memset(gradPE_X,0,3*num_nodes_in_block*sizeof(D_REAL));
	
	D_REAL* rho_species;
	D_REAL* temp_rho_species = Field_Type[id_rho_rez];
	
	
      for(INT32 species=0; species<num_Charged_Species; species++)
      {

	 //! copy rho_species to temporary rho
	 rho_species = Field_Type[id_allRhoSpecies] +species*num_nodes_in_block;
         memcpy(temp_rho_species,rho_species, num_nodes_in_block*sizeof(D_REAL));

//!NOTE: modification of grad PE:
//!NOTE: component of gradPE due to obstacle is set to zero

	 //! apply inner density to temporary rho and use number density !!!
         for(INT32 node=0; node<num_nodes_in_block; node++)
         {
		//! apply inner density in case of obstacle node
		if(Flag[node])
		temp_rho_species[node] = obstacle_rho[species];
	
		//! Always use ion charge density but not number density!!!
		//! See Normalization for further details

		temp_rho_species[node] = Electron_Betas[species]*pow(temp_rho_species[node],kappa_electron);
	 }

         for (INT32 i=1; i < BlkNds_X-1; i++)
          for (INT32 j=1; j < BlkNds_Y-1; j++)
           for (INT32 k=1; k < BlkNds_Z-1; k++)
	   {


		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
		//! --- indices for firs derivatives ---------
		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

		//! two times faster than gradPE calculation above since
		//! several pow() calls are avoided
                
                if(Flag[ip1_j_k] || Flag[im1_j_k])
		gradPE_X[i_j_k] += 0;
		else                
                gradPE_X[i_j_k] +=
                (
			 temp_rho_species[ip1_j_k]
			-temp_rho_species[im1_j_k]
		)*rd[0];

		if(Flag[i_jp1_k] || Flag[i_jm1_k])
		gradPE_Y[i_j_k] += 0;
		else      
                gradPE_Y[i_j_k] +=
                (
			 temp_rho_species[i_jp1_k] 
			-temp_rho_species[i_jm1_k]
		)*rd[1];

		if(Flag[i_j_kp1] || Flag[i_j_km1])
		gradPE_Z[i_j_k] += 0;
		else      
                gradPE_Z[i_j_k] +=
                (
			 temp_rho_species[i_j_kp1] 
			-temp_rho_species[i_j_km1]
		)*rd[2];


	   }
	}
}



//!-------------------------------------------------------------//
//! calc_macro_Force: -							//
//!-------------------------------------------------------------//
void CBlock::calc_macro_Force(INT32 MacroForce_id, INT32 EField_id, INT32 UField_id,  INT32 BField_id)
{


	D_REAL* F1 = Field_Type[MacroForce_id];
	D_REAL* F2 = F1 +num_nodes_in_block;
	D_REAL* F3 = F2 +num_nodes_in_block;

	D_REAL* E1 = Field_Type[EField_id];
	D_REAL* E2 = E1 +num_nodes_in_block;
	D_REAL* E3 = E2 +num_nodes_in_block;

	D_REAL* U1 = Field_Type[UField_id];
	D_REAL* U2 = U1 +num_nodes_in_block;
	D_REAL* U3 = U2 +num_nodes_in_block;

	D_REAL* B1 = Field_Type[BField_id];
	D_REAL* B2 = B1 +num_nodes_in_block;
	D_REAL* B3 = B2 +num_nodes_in_block;


        for(INT32 node=0; node<num_nodes_in_block; node++)
        {
	
		F1[node] = E1[node] + (U2[node]*B3[node] - U3[node]*B2[node]);
		F2[node] = E2[node] + (U3[node]*B1[node] - U1[node]*B3[node]);
		F3[node] = E3[node] + (U1[node]*B2[node] - U2[node]*B1[node]);
		

	}
}


//!-------------------------------------------------------------//
//! estimate_nodes_BBound_value:						//
//!-------------------------------------------------------------//
void CBlock::estimate_nodes_BBound_value(D_REAL *BBound, INT32* cell_indices)
{


/*
	//! TWO STREAM BOUNDARY
	D_REAL width = 3.;

	PARTICLE_REAL x[3], abs_x1;
	PARTICLE_REAL offset[3] = {0., 0., 0.};
	intern2normedCoords(x, offset, cell_indices);
	memset(BBound, 0, 3 *sizeof(D_REAL));


	abs_x1 = fabs(x[1]);


	if(abs_x1 < width)
	BBound[2] =  -sin(x[1]/width *0.5 *M_PI);
	else
	BBound[2] =  -x[1]/abs_x1;
*/



	int TL_rotate =  TL -TL_new_Bsw;

	//! rotate from 0 tp pi
	D_REAL angle = 0.5 *M_PI *TL_rotate/(num_TL_adapt-1);

	BBound[0] = cos(angle)*B_sw[0] +sin(angle)*new_B_sw[0];
	BBound[1] = cos(angle)*B_sw[1] +sin(angle)*new_B_sw[1];
	BBound[2] = cos(angle)*B_sw[2] +sin(angle)*new_B_sw[2];


}


//!-------------------------------------------------------------//
//! apply_new_BField_Boundaries:						//
//!-------------------------------------------------------------//
void CBlock::apply_new_BField_Boundaries(INT32 B_Field_type)
{

	INT32 i,j,k,i_j_k;
	D_REAL BBound[3];

	D_REAL* BX = Field_Type[B_Field_type];
	D_REAL* BY = BX +num_nodes_in_block;
	D_REAL* BZ = BY +num_nodes_in_block;


	//! left box boundary
	i=0;
	if(is_box_boundary[0] && use_B_inflow_bounds[0])
	 for (j=0; j < BlkNds_Y; j++)
	  for (k=0; k < BlkNds_Z; k++)
	  {
	
		i_j_k =  i*BlkNds_Y*BlkNds_Z
			+j*BlkNds_Z
			+k;

		INT32 cell_indices[3] = {i,j,k};
		estimate_nodes_BBound_value(BBound, cell_indices);

		BX[i_j_k] = BBound[0];
		BY[i_j_k] = BBound[1];
		BZ[i_j_k] = BBound[2];

	  }

	//! right box boundary
	i=BlkNds_X-1;
	if(is_box_boundary[1] && use_B_inflow_bounds[1])
	 for (j=0; j < BlkNds_Y; j++)
	  for (k=0; k < BlkNds_Z; k++)
	  {
	
		i_j_k =  i*BlkNds_Y*BlkNds_Z
			+j*BlkNds_Z
			+k;

		INT32 cell_indices[3] = {i,j,k};
		estimate_nodes_BBound_value(BBound, cell_indices);

		BX[i_j_k] = BBound[0];
		BY[i_j_k] = BBound[1];
		BZ[i_j_k] = BBound[2];

	  }

	//! front box boundary
	j=0;
	if(is_box_boundary[2] && use_B_inflow_bounds[2])
	 for (i=0; i < BlkNds_X; i++)
	  for (k=0; k < BlkNds_Z; k++)
	  {
	
		i_j_k =  i*BlkNds_Y*BlkNds_Z
			+j*BlkNds_Z
			+k;

		INT32 cell_indices[3] = {i,j,k};
		estimate_nodes_BBound_value(BBound, cell_indices);

		BX[i_j_k] = BBound[0];
		BY[i_j_k] = BBound[1];
		BZ[i_j_k] = BBound[2];

	  }

	//! back box boundary
	j=BlkNds_Y-1;
	if(is_box_boundary[3] && use_B_inflow_bounds[3])
	 for (i=0; i < BlkNds_X; i++)
	  for (k=0; k < BlkNds_Z; k++)
	  {
	
		i_j_k =  i*BlkNds_Y*BlkNds_Z
			+j*BlkNds_Z
			+k;

		INT32 cell_indices[3] = {i,j,k};
		estimate_nodes_BBound_value(BBound, cell_indices);

		BX[i_j_k] = BBound[0];
		BY[i_j_k] = BBound[1];
		BZ[i_j_k] = BBound[2];

	  }

	//! top box boundary
	k=0;
	if(is_box_boundary[4] && use_B_inflow_bounds[4])
	 for (i=0; i < BlkNds_X; i++)
	  for (j=0; j < BlkNds_Y; j++)
	  {
	
		i_j_k =  i*BlkNds_Y*BlkNds_Z
			+j*BlkNds_Z
			+k;

		INT32 cell_indices[3] = {i,j,k};
		estimate_nodes_BBound_value(BBound, cell_indices);

		BX[i_j_k] = BBound[0];
		BY[i_j_k] = BBound[1];
		BZ[i_j_k] = BBound[2];

	  }

	//! bottom box boundary
	k=BlkNds_Z-1;
	if(is_box_boundary[5] && use_B_inflow_bounds[5])
	 for (i=0; i < BlkNds_X; i++)
	  for (j=0; j < BlkNds_Y; j++)
	  {
	
		i_j_k =  i*BlkNds_Y*BlkNds_Z
			+j*BlkNds_Z
			+k;

		INT32 cell_indices[3] = {i,j,k};
		estimate_nodes_BBound_value(BBound, cell_indices);

		BX[i_j_k] = BBound[0];
		BY[i_j_k] = BBound[1];
		BZ[i_j_k] = BBound[2];

	  }




}



//!-------------------------------------------------------------//
//! FIMesh_to_HIMesh: -								//
//!-------------------------------------------------------------//
void CBlock::FIMesh_to_HIMesh(INT32 id_dest, INT32 id_src)
{


      D_REAL* Dest = Field_Type[id_dest];
      D_REAL* Src = Field_Type[id_src];


	//! loop over Full Integer Mesh
	for(INT32 i=1; i < BlkNds_X-1; i++)
	 for(INT32 j=1; j < BlkNds_Y-1; j++)
	  for(INT32 k=1; k < BlkNds_Z-1; k++)
	  {

		//!INTERPOLATION (full to half):
		//!ijk  0   1   2   3   (FULL)
		//!    (*)  *   *  (*)
		//!          \ / \ /
		//!      (-)  -   -  (-)
		//!uvw    0   1   2   3 (HALF)
		//! (do not interpolate to ghost nodes)

		//! ------------------------------------------
		INT32 i_j_k   =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		//! here u_v_w = i_j_k (see sketch)

		//! ------------------------------------------
		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		
		//! -------------------------------------------
		INT32 ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
		
		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);


		//! loop over all components
		for(INT32 comp=0; comp<COMPs_FType[id_dest]; comp++)
		{

			INT32 offset = comp*num_nodes_in_block;

			Dest[i_j_k +offset] =

			0.125*(Src[i_j_k +offset]

			      +Src[ip1_j_k +offset]
			      +Src[i_jp1_k +offset]
			      +Src[i_j_kp1 +offset]

			      +Src[ip1_jp1_k +offset]
			      +Src[ip1_j_kp1 +offset]
			      +Src[i_jp1_kp1 +offset]

			      +Src[ip1_jp1_kp1 +offset]);
		}

	  }



}


//!-------------------------------------------------------------//
//! HIMesh_to_FIMesh: -								//
//!-------------------------------------------------------------//
void CBlock::HIMesh_to_FIMesh(INT32 id_dest, INT32 id_src)
{


	
 
      D_REAL* Dest = Field_Type[id_dest];
      D_REAL* Src = Field_Type[id_src];


	//! loop over Half Integer Mesh
	for(INT32 u=0; u < BlkNds_X-2; u++)
	 for(INT32 v=0; v < BlkNds_Y-2; v++)
	  for(INT32 w=0; w < BlkNds_Z-2; w++)
	  {


		//!INTERPOLATION (half to full):
		//!ijk  0   1   2   3   (FULL)
		//!    (*)  *   *  (*)
		//!        / \ / \
		//!      (-)  -   -  (-)
		//!uvw    0   1   2   3 (HALF)
		//! (do not interpolate to ghost nodes)


		INT32 u_v_w   =      u*BlkNds_Y*BlkNds_Z    +v*BlkNds_Z      +w;

		//! ------------------------------------------
		INT32 up1_v_w = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z     +w;
		INT32 u_vp1_w =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z     +w;
		INT32 u_v_wp1 =     u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +(w+1);
		
		//! -------------------------------------------
		INT32 up1_vp1_w = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z     +w;
		INT32 up1_v_wp1 = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +(w+1);
		INT32 u_vp1_wp1 =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +(w+1);
		
		INT32 up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +(w+1);
		
		
		//! ------------------------------------------
		INT32 i = u+1;
		INT32 j = v+1;
		INT32 k = w+1;

		INT32 i_j_k   =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;


		//! loop over all components
		for(INT32 comp=0; comp<COMPs_FType[id_dest]; comp++)
		{

			INT32 offset = comp*num_nodes_in_block;

			Dest[i_j_k +offset] =

			0.125*(  Src[u_v_w +offset]

				+Src[up1_v_w +offset]
				+Src[u_vp1_w +offset]
				+Src[u_v_wp1 +offset]

				+Src[up1_vp1_w +offset]
				+Src[up1_v_wp1 +offset]
				+Src[u_vp1_wp1 +offset]

				+Src[up1_vp1_wp1 +offset]);
		}

	}
	
}






//!-------------------------------------------------------------//
//! calc_grad: -									//
//!-------------------------------------------------------------//
void CBlock::calc_grad(INT32 in_Field_type, INT32 grad_Field_type)
{

	
	D_REAL *A, *AX, *AY, *AZ, *gradAX, *gradAY, *gradAZ;
	
	if(COMPs_FType[in_Field_type]==1)
	D_REAL* A = Field_Type[in_Field_type];	
	
	if(COMPs_FType[in_Field_type]==3)
	{
		AX = Field_Type[in_Field_type];	
		AY = AX +num_nodes_in_block;
		AZ = AY +num_nodes_in_block;
	}

	gradAX = Field_Type[grad_Field_type];	
	gradAY = gradAX +num_nodes_in_block;
	gradAZ = gradAY +num_nodes_in_block;
	
	rd = rd_of_L[RLevel];
	

	//! set zero
	memset(gradAX, 0, 3*num_nodes_in_block*sizeof(D_REAL));
	
	
	for (INT32 i=1; i < BlkNds_X-1; i++)
	 for (INT32 j=1; j < BlkNds_Y-1; j++)
	  for (INT32 k=1; k < BlkNds_Z-1; k++)
	  {

		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);


		if(COMPs_FType[in_Field_type]==1)
		{

			gradAX[i_j_k] = rd[0]*(A[ip1_j_k]-A[im1_j_k]);
			gradAY[i_j_k] = rd[1]*(A[i_jp1_k]-A[i_jm1_k]);
			gradAZ[i_j_k] = rd[2]*(A[i_j_kp1]-A[i_j_km1]);
		}

		if(COMPs_FType[in_Field_type]==3)
		{

			gradAX[i_j_k] = rd[0]*( sqrt(AX[ip1_j_k]*AX[ip1_j_k] +AY[ip1_j_k]*AY[ip1_j_k] +AZ[ip1_j_k]*AZ[ip1_j_k])
						-sqrt(AX[im1_j_k]*AX[im1_j_k] +AY[im1_j_k]*AY[im1_j_k] +AZ[im1_j_k]*AZ[im1_j_k]));

			gradAY[i_j_k] = rd[1]*( sqrt(AX[i_jp1_k]*AX[i_jp1_k] +AY[i_jp1_k]*AY[i_jp1_k] +AZ[i_jp1_k]*AZ[i_jp1_k])
						-sqrt(AX[i_jm1_k]*AX[i_jm1_k] +AY[i_jm1_k]*AY[i_jm1_k] +AZ[i_jm1_k]*AZ[i_jm1_k]));

			gradAZ[i_j_k] = rd[2]*( sqrt(AX[i_j_kp1]*AX[i_j_kp1] +AY[i_j_kp1]*AY[i_j_kp1] +AZ[i_j_kp1]*AZ[i_j_kp1])
						-sqrt(AX[i_j_km1]*AX[i_j_km1] +AY[i_j_km1]*AY[i_j_km1] +AZ[i_j_km1]*AZ[i_j_km1]));

		}
	  }
	
}



//!-------------------------------------------------------------//
//! calc_divB: -									//
//!-------------------------------------------------------------//
void CBlock::calc_div(INT32 B_Field_type, INT32 div_Field_type)
{

      D_REAL* AX = Field_Type[B_Field_type];	
      D_REAL* AY = AX +num_nodes_in_block;
      D_REAL* AZ = AY +num_nodes_in_block;

      rd = rd_of_L[RLevel];

      D_REAL* divA = Field_Type[div_Field_type];
      //! set zero for visu auto color
      memset(divA, 0, num_nodes_in_block*sizeof(D_REAL));


      for (INT32 i=1; i < BlkNds_X-1; i++)
       for (INT32 j=1; j < BlkNds_Y-1; j++)
        for (INT32 k=1; k < BlkNds_Z-1; k++)
	{

      		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

	    	divA[i_j_k] = rd[0]*(AX[ip1_j_k]-AX[im1_j_k])
	    		     +rd[1]*(AY[i_jp1_k]-AY[i_jm1_k])
	    		     +rd[2]*(AZ[i_j_kp1]-AZ[i_j_km1]);
	}
	
}

//!-------------------------------------------------------------//
//! calc_rotB: -							//
//!-------------------------------------------------------------//
void CBlock::calc_rot(INT32 B_Field_type, INT32 rot_Field_type)
{


      rotBX = Field_Type[rot_Field_type];
      rotBY = rotBX +num_nodes_in_block;
      rotBZ = rotBY +num_nodes_in_block;

      //! set zero for visu auto color
      memset(rotBX, 0, 3*num_nodes_in_block* sizeof(D_REAL));

      BX = Field_Type[B_Field_type];		
      BY = BX +num_nodes_in_block;
      BZ = BY +num_nodes_in_block;

      rd = rd_of_L[RLevel];



      for (INT32 i=1; i < BlkNds_X-1; i++)
       for (INT32 j=1; j < BlkNds_Y-1; j++)
        for (INT32 k=1; k < BlkNds_Z-1; k++)
// 	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{

      		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);



		rotBX[i_j_k] = (  rd[1]*(BZ[i_jp1_k]-BZ[i_jm1_k]) 
				- rd[2]*(BY[i_j_kp1]-BY[i_j_km1]));

		rotBY[i_j_k] = (  rd[2]*(BX[i_j_kp1]-BX[i_j_km1]) 
				- rd[0]*(BZ[ip1_j_k]-BZ[im1_j_k]));

		rotBZ[i_j_k] = (  rd[0]*(BY[ip1_j_k]-BY[im1_j_k])
		 		- rd[1]*(BX[i_jp1_k]-BX[i_jm1_k]));

	}
	
}



//!-------------------------------------------------------------//
//! calc_divB: -								//
//!-------------------------------------------------------------//
void CBlock::calc_local_gyro_radius(INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField, D_REAL SI_x0)
{


	D_REAL  abs_B, abs_U;

	D_REAL* BX = Field_Type[id_BField];		
	D_REAL* BY = BX +num_nodes_in_block;
	D_REAL* BZ = BY +num_nodes_in_block;
	
	D_REAL* UX = Field_Type[id_UField];		
	D_REAL* UY = UX +num_nodes_in_block;
	D_REAL* UZ = UY +num_nodes_in_block;
	
	D_REAL* Gyro_Radius = Field_Type[id_Gyro_Radius];


	for(INT32 elm=0; elm<num_nodes_in_block; elm++)
	{

		abs_B = sqrt(BX[elm]*BX[elm] +BY[elm]*BY[elm] +BZ[elm]*BZ[elm]);
		abs_U = sqrt(UX[elm]*UX[elm] +UY[elm]*UY[elm] +UZ[elm]*UZ[elm]);

		Gyro_Radius[elm] = (abs_U * Ion_Masses[species])*SI_x0 / (abs_B * Ion_Charges[species]);


	}
}

//!-------------------------------------------------------------//
//! calc_local_electron_gyro_radius: -                                                           //
//!-------------------------------------------------------------//
void CBlock::calc_local_electron_gyro_radius (INT32 id_Gyro_Radius, INT32 species, INT32 id_BField, INT32 id_UField, D_REAL SI_x0)
{
    
    
    D_REAL  abs_B, abs_U;
    
    D_REAL* BX = Field_Type[id_BField];             
    D_REAL* BY = BX +num_nodes_in_block;
    D_REAL* BZ = BY +num_nodes_in_block;
    
    D_REAL* UX = Field_Type[id_UField];             
    D_REAL* UY = UX +num_nodes_in_block;
    D_REAL* UZ = UY +num_nodes_in_block;
    
    D_REAL* Gyro_Radius = Field_Type[id_Gyro_Radius];
    
    //! Electronen Masse
    //! 1800*m_e = m_p
    for(INT32 elm=0; elm<num_nodes_in_block; elm++)
    {
        
        abs_B = sqrt(BX[elm]*BX[elm] +BY[elm]*BY[elm] +BZ[elm]*BZ[elm]);
        abs_U = sqrt(UX[elm]*UX[elm] +UY[elm]*UY[elm] +UZ[elm]*UZ[elm]);
        
        Gyro_Radius[elm] = (abs_U)*SI_x0 / (abs_B * 1800 * Ion_Charges[species]);
        
        
    }
}



//!-------------------------------------------------------------//
//! calc_divB: -								//
//!-------------------------------------------------------------//
void CBlock::add_force(INT32 field_id)
{


	D_REAL  r_vec[3], null_vec[3]={0,0,0};

	D_REAL* AX = Field_Type[field_id];		
	D_REAL* AY = AX +num_nodes_in_block;
	D_REAL* AZ = AY +num_nodes_in_block;
	

	D_REAL region = 2.*R_Obstacle;
	D_REAL force_mag = -0.1;

	for (INT32 i=0; i < BlkNds_X; i++)
	 for (INT32 j=0; j < BlkNds_Y; j++)
	  for (INT32 k=0; k < BlkNds_Z; k++)
	  {
	


			INT32 cell_indices[3]={i,j,k};
			intern2normedCoords(r_vec, null_vec, cell_indices);

			INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

			D_REAL dist = vec_len(r_vec);

			r_vec[0] /= dist;
			r_vec[1] /= dist;
			r_vec[2] /= dist;

			D_REAL zAxis[3] = {0., 0., 1.};
			D_REAL force[3] = {r_vec[0], r_vec[1], r_vec[2]};

// 			vec_cross(force, r_vec, zAxis);

			if(dist > region)
			{

				AX[i_j_k] += force_mag * (dist-region) *force[0];
				AY[i_j_k] += force_mag * (dist-region) *force[1];
				AZ[i_j_k] += force_mag * (dist-region) *force[2];

			}


	  }
	
}


//!-------------------------------------------------------------//
//! calc_E: -							//
//!-------------------------------------------------------------//
void CBlock::calc_E(INT32 BField_type, INT32 UField_type, INT32 rho_type)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

	//! POINTER to B & rotB are set in calc_rotB
	//! B_Total has to be set extra !
	calc_rot(BField_type, id_rotB);
	
	
	EX = Field_Type[id_EField];
	EY = EX +num_nodes_in_block;
	EZ = EY +num_nodes_in_block;
	
	BX_total = Field_Type[id_BTotal];
	BY_total = BX_total +num_nodes_in_block;
	BZ_total = BY_total +num_nodes_in_block;
	
	UX = Field_Type[UField_type];		
	UY = UX +num_nodes_in_block;
	UZ = UY +num_nodes_in_block;
	
	eta = Field_Type[id_Eta];
	rho_total = Field_Type[rho_type];
	
#if defined nonadiabatic_gradPE_TERM
	D_REAL* PE = Field_Type[id_PEtotal];
#endif
	
	
	//! In hybrid cycle the "id_UI_minus" is only used
	//! in between the following functions: 
	//! 5.) Ui_minus is collected
	//! 	and 
	//! 9.) Magnetic field is advance 
	
	//! So at calc first and second E (1. & 3.) it
	//! can be used as a schratch to add all "id_allRhoSpecies"
	//! pressure gradients:
	
	
#ifndef nonadiabatic_gradPE_TERM
	//! rho_rez only required in and directly before BField method, so use it
	//! as a sratch here
	D_REAL* temp_rho_species = Field_Type[id_rho_rez];
#endif
	
	
	D_REAL* gradPE_X = Field_Type[id_UI_minus];
	D_REAL* gradPE_Y = gradPE_X +num_nodes_in_block;
	D_REAL* gradPE_Z = gradPE_Y +num_nodes_in_block;
	
	memset(gradPE_X,0,3*num_nodes_in_block*sizeof(D_REAL));
	
	D_REAL* rho_species = NULL;
	
	
	rd = rd_of_L[RLevel];

#if (defined gradPE_TERM) && !(defined  nonadiabatic_gradPE_TERM)
      for(INT32 species=0; species<num_Charged_Species; species++)
      {

	 //! copy rho_species to temporary rho
	 rho_species = Field_Type[id_allRhoSpecies] +species*num_nodes_in_block;
         memcpy(temp_rho_species,rho_species, num_nodes_in_block*sizeof(D_REAL));

	 //! apply inner density to temporary rho and use number density !!!
         for(INT32 node=0; node<num_nodes_in_block; node++)
         {
		//! apply inner density in case of obstacle node
		if(Flag[node])
		temp_rho_species[node] = obstacle_rho[species];
	
		//! Always use ion charge density but not number density!!!
		//! See Normalization for further details

		temp_rho_species[node] = Electron_Betas[species]*pow(temp_rho_species[node],kappa_electron);
	 }

         for (INT32 i=1; i < BlkNds_X-1; i++)
          for (INT32 j=1; j < BlkNds_Y-1; j++)
           for (INT32 k=1; k < BlkNds_Z-1; k++)
	   {


		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
		//! --- indices for firs derivatives ---------
		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

		//! two times faster than gradPE calculation above since
		//! several pow() calls are avoided
                gradPE_X[i_j_k] +=
                (
			 temp_rho_species[ip1_j_k]
			-temp_rho_species[im1_j_k]
		)*rd[0];

                gradPE_Y[i_j_k] +=
                (
			 temp_rho_species[i_jp1_k] 
			-temp_rho_species[i_jm1_k]
		)*rd[1];


                gradPE_Z[i_j_k] +=
                (
			 temp_rho_species[i_j_kp1] 
			-temp_rho_species[i_j_km1]
		)*rd[2];


	   }
	}
#endif
#if defined nonadiabatic_gradPE_TERM
	for (INT32 i=1; i < BlkNds_X-1; i++)
	  for (INT32 j=1; j < BlkNds_Y-1; j++)
	    for (INT32 k=1; k < BlkNds_Z-1; k++)
	    {
		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
		//! --- indices for first derivatives ---------
		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
		
		
		gradPE_X[i_j_k] = ( PE[ip1_j_k] - PE[im1_j_k] ) * rd[0];
		
		gradPE_Y[i_j_k] = ( PE[i_jp1_k] - PE[i_jm1_k] ) * rd[1];
		
		gradPE_Z[i_j_k] = ( PE[i_j_kp1] - PE[i_j_km1] ) * rd[2];
	    }
#endif

	//! As EField is used as a scratch in RK Method,
	//! also boundaries have to be recalculated. 
	//! This can be easily done within the the loop that
	//! includes the boundary Nodes: 
	//! - E_Convec  is the usual Boundary Condition 
	//! - E_Hall    is zero as rotB zero on Boundaries 
	//! - E_gradPE  is zero (not calculated on Bounds)
	//! - E_Eta	   is zero as rotB zero on Boundaries 
	//! No loss in computation speed is visible when 
	//! recalculating the EBound Values


	D_REAL rez_rho_value = 0;

//       for (i=1-is_box_boundary[0]; i < BlkNds_X-1+is_box_boundary[1]; i++)
//        for (j=1-is_box_boundary[2]; j < BlkNds_Y-1+is_box_boundary[3]; j++)
//         for (k=1-is_box_boundary[4]; k < BlkNds_Z-1+is_box_boundary[5]; k++)

      for (INT32 i=1; i < BlkNds_X-1; i++)
       for (INT32 j=1; j < BlkNds_Y-1; j++)
        for (INT32 k=1; k < BlkNds_Z-1; k++)
	if(!Flag[i*BlkNds_Y*BlkNds_Z+j*BlkNds_Z+k] || advance_obstacle_E)
	{


      	    INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;


	    //! MCD check can be removed, when already tested in
	    //! J2u. However, as is does not slow down the Code,
	    //! it can be left here.
	    if(rho_total[i_j_k]<MCD_BField)
	    rez_rho_value = 1./MCD_BField;
	    else
	    rez_rho_value = 1./rho_total[i_j_k];


//!--------------- E1 ----------------------------------
	    EX[i_j_k] =
#ifdef CONVEC_TERM
  	    //! ----- -(uxB) -----------------------
		- (   UY[i_j_k] * BZ_total[i_j_k] 
	             -UZ[i_j_k] * BY_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rez_rho_value *( rotBY[i_j_k] * BZ_total[i_j_k] 
						-  rotBZ[i_j_k] * BY_total[i_j_k])
#endif
#ifdef gradPE_TERM
	    //!------ grad PE --------------------
                - 0.5 * rez_rho_value  *gradPE_X[i_j_k]
#endif

#ifdef ETA_TERM_EField
	    //!------ eta * rotB--------------------
	       + eta[i_j_k] * rotBX[i_j_k]
#endif

	    ;

//!--------------- E2 ----------------------------------
	    EY[i_j_k] =
#ifdef CONVEC_TERM
	    //! ----- -(uxB) -----------------------
		-(   UZ[i_j_k] * BX_total[i_j_k] 
		    -UX[i_j_k] * BZ_total[i_j_k])
#endif

#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rez_rho_value *( rotBZ[i_j_k] * BX_total[i_j_k] 
		      				-  rotBX[i_j_k] * BZ_total[i_j_k])
#endif

#ifdef gradPE_TERM
	    //!------ grad PE --------------------
                - 0.5 * rez_rho_value *gradPE_Y[i_j_k]
#endif

#ifdef ETA_TERM_EField
	    //!------ eta * rotB--------------------
	       + eta[i_j_k] * rotBY[i_j_k];
#endif

	    ;

//!--------------- E3 ----------------------------------
	    EZ[i_j_k] =
#ifdef CONVEC_TERM
	    //! ----- -(uxB) -----------------------
		-(   UX[i_j_k] * BY_total[i_j_k] 
		    -UY[i_j_k] * BX_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rez_rho_value *( rotBX[i_j_k] * BY_total[i_j_k] 
						-  rotBY[i_j_k] * BX_total[i_j_k])
#endif
#ifdef gradPE_TERM
	    //!------ grad PE --------------------
                - 0.5 * rez_rho_value *gradPE_Z[i_j_k]

#endif
#ifdef ETA_TERM_EField
	    //!------ eta * rotB--------------------
	       + eta[i_j_k] * rotBZ[i_j_k];
#endif

	    ;

	
	}

	//! cpy boundaries in case field outflow condition
	//! in case inflow correct BV are calculated, see above
	//! for details.
	cp_Boundaries(id_EField, use_E_inflow_bounds);

// 	if(is_box_boundary[Xmin_BB] && origin[1]>=0)
// 	{
// 		bool cp_bound[6] = {0,1,1,1,1,1};
// 		cp_Boundaries(id_EField, cp_bound);
// 
// 	}
// // 
// 	if(is_box_boundary[Xmax_BB] && origin[1]<0)
// 	{
// 		bool cp_bound[6] = {1,0,1,1,1,1};
// 		cp_Boundaries(id_EField, cp_bound);
// 
// 	}

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

	
}

//!-------------------------------------------------------------//
//! CAM: -							//
//!-------------------------------------------------------------//
void CBlock::CAM(void)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	//! id_rho_np1 and id_rho_n should be equal here
	rho_n = Field_Type[id_rho_n];
	
	//! use id_rotB as temp scratch for Ji_plus
	JiX_plus = Field_Type[id_rotB];
	JiY_plus = JiX_plus +num_nodes_in_block;
	JiZ_plus = JiY_plus +num_nodes_in_block;

	UiX_plus = Field_Type[id_UI_plus];
	UiY_plus = UiX_plus +num_nodes_in_block;
	UiZ_plus = UiY_plus +num_nodes_in_block;

	//! use id_UI_minus as temp scratch for Ji_np1
	JiX_np1 = Field_Type[id_UI_minus];
	JiY_np1 = JiX_np1 +num_nodes_in_block;
	JiZ_np1 = JiY_np1 +num_nodes_in_block;


	Lam = Field_Type[id_Lam];

	GamX = Field_Type[id_Gam];
	GamY = GamX +num_nodes_in_block;
	GamZ = GamY +num_nodes_in_block;


	BX = Field_Type[id_BTotal];
	BY = BX +num_nodes_in_block;
	BZ = BY +num_nodes_in_block;

	EX = Field_Type[id_EField];
	EY = EX +num_nodes_in_block;
	EZ = EY +num_nodes_in_block;


	//! To be 100% precise here also could be a MCD_J2U check.
	//! But assume rho to be zero, than Lam also is zero.

	//! Hence Ji_plus is zero (due to rho=0),
	//! and Ji_np1 also remains zero (due to Ji_plus=0, LAM=0).
	//! This means UiX_plus is zero also, even though  
	//! Ji_plus its multiplied by 1/MCD_J2U.

	//! As UiX_plus would be expected to be zero, this yields 
	//! the correct expression.

	//! CAM must not be performed on Box Boundary, as no boundary values
	//! are provided for Lam & Gam
	//! -> UI_Plus remains unchanged
	INT32 start[3] = {is_box_boundary[0],
		         is_box_boundary[2],
		         is_box_boundary[4]};

	INT32 end[3] =   {is_box_boundary[1],
		         is_box_boundary[3],
		         is_box_boundary[5]};

      for(INT32 i=start[0]; i < BlkNds_X-end[0]; i++)
       for(INT32 j=start[1]; j < BlkNds_Y-end[1]; j++)
        for(INT32 k=start[2]; k < BlkNds_Z-end[2]; k++)
	{

		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		JiX_plus[i_j_k] = rho_n[i_j_k] *UiX_plus[i_j_k];
		JiY_plus[i_j_k] = rho_n[i_j_k] *UiY_plus[i_j_k];
		JiZ_plus[i_j_k] = rho_n[i_j_k] *UiZ_plus[i_j_k];
	}


      for (INT32 i=start[0]; i < BlkNds_X-end[0]; i++)
       for (INT32 j=start[1]; j < BlkNds_Y-end[1]; j++)
        for (INT32 k=start[2]; k < BlkNds_Z-end[2]; k++)
	{

		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;


		//! Check why half time step ???
		JiX_np1[i_j_k] =  JiX_plus[i_j_k] +0.5 *dt_particle_of_L[RLevel]

			     *(
				 Lam[i_j_k] * EX[i_j_k]
				+(
				   GamY[i_j_k] * BZ[i_j_k]
				  -GamZ[i_j_k] * BY[i_j_k]
				 )
			     );

		JiY_np1[i_j_k] =  JiY_plus[i_j_k] +0.5 *dt_particle_of_L[RLevel]

			     *(
				 Lam[i_j_k] * EY[i_j_k]
				+(
				   GamZ[i_j_k] * BX[i_j_k]
				  -GamX[i_j_k] * BZ[i_j_k]
				 )
			     );

		JiZ_np1[i_j_k] =  JiZ_plus[i_j_k] +0.5 *dt_particle_of_L[RLevel]

			     *(
				 Lam[i_j_k] * EZ[i_j_k]
				+(
				   GamX[i_j_k] * BY[i_j_k]
				  -GamY[i_j_k] * BX[i_j_k]
				 )
			     );
	}


      D_REAL rez_rho_value = 0.;
	
      for (INT32 i=start[0]; i < BlkNds_X-end[0]; i++)
       for (INT32 j=start[1]; j < BlkNds_Y-end[1]; j++)
        for (INT32 k=start[2]; k < BlkNds_Z-end[2]; k++)
	{

		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		//! Since MCD is only set internally it has to
		//! be set here again.
		if(rho_n[i_j_k]<MCD_J2U)
		rez_rho_value = 1./MCD_J2U;
		else
		rez_rho_value = 1./rho_n[i_j_k];

		UiX_plus[i_j_k] = JiX_np1[i_j_k] *rez_rho_value;
		UiY_plus[i_j_k] = JiY_np1[i_j_k] *rez_rho_value;
		UiZ_plus[i_j_k] = JiZ_np1[i_j_k] *rez_rho_value;



		//! set tangential velocity inside obstacle if required
		if(Flag[i_j_k])
		{

			D_REAL omega[3] = {omega_rotating_obs[0],
					   omega_rotating_obs[1],
					   omega_rotating_obs[2]};

			if(vec_len2(omega)>0.)
			{
				D_REAL r_vec[3], v_tang[3];
				D_REAL null_vec[3] = {0,0,0};
				INT32 cell_indices[3] = {i,j,k};


				intern2normedCoords(r_vec, null_vec, cell_indices);
				vec_cross(v_tang, omega, r_vec);
	
				UiX_plus[i_j_k] = v_tang[0];
				UiY_plus[i_j_k] = v_tang[1];
				UiZ_plus[i_j_k] = v_tang[2];
			}

		}



	}


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}





//!-------------------------------------------------------------//
//! smooth_Field: -						//
//!-------------------------------------------------------------//
#if defined(use_new_smooth_Field_method)
//! new version
void CBlock::smooth_Field(INT32 type, D_REAL smooth_value)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	INT32 a, b, c;
	INT32 i, j, k, i_j_k;
	INT32 comp, temp1;
	
	D_REAL* field = Field_Type[type];
	
	//! use divB as scratch
	D_REAL* scratch = Field_Type[id_divB];
	
#if defined(use_vectorclass)
	VEC2_D_REAL vec_field;
	VEC2_D_REAL vec_scratch;
#endif
	
	for(comp=0; comp < COMPs_FType[type]; ++comp)
	{
		memset(scratch,0,num_nodes_in_block*sizeof(D_REAL));
		
		for (i=1; i < BlkNds_X-1; ++i)
		{
		  for (j=1; j < BlkNds_Y-1; ++j)
		  {
		    for (k=1; k < BlkNds_Z-1; ++k)
		    {
			  i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		  
			  if(!Core_Flag[i_j_k])
			  {
			    for (a=-1; a<=1; ++a)
			    {
			      for (b=-1; b<=1; ++b)
			      {
				for (c=-1; c<=1; ++c)
				{
				  temp1 = comp*num_nodes_in_block + (i+a)*BlkNds_Y*BlkNds_Z + (j+b)*BlkNds_Z + (k+c);
				  scratch[i_j_k] += smooth_coeff[a+1][b+1][c+1] * field[temp1];
				}
			      }
			    }
			  }
			  else
			  {
				  scratch[i_j_k] = field[i_j_k];
			  }
		    }
		  }
		}
		
		for (i=1; i < BlkNds_X-1; ++i)
		{
		  for (j=1; j < BlkNds_Y-1; ++j)
		  {
#if defined(use_vectorclass)
		    for (k=1; k < BlkNds_Z-1; k+=2)
		    {
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
			
			vec_field.load( field + (comp*num_nodes_in_block + i_j_k) );
			vec_scratch.load( scratch + i_j_k );
			
			vec_field = (1.-smooth_value) * vec_field + smooth_value * vec_scratch;
			
			vec_field.store( field + (comp*num_nodes_in_block + i_j_k) );
		    }
#else
		    for (k=1; k < BlkNds_Z-1; ++k)
		    {
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
			
			field[comp*num_nodes_in_block + i_j_k] = 
			  (1.-smooth_value) * field[comp*num_nodes_in_block + i_j_k] + smooth_value * scratch[i_j_k];
		    }
#endif
		  }
		}
	}

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;
}
#else
//! old version 
void CBlock::smooth_Field(INT32 type,D_REAL smooth_value)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	INT32 a,b,c;
	
	//! use divB as scratch
	D_REAL* scr1 = Field_Type[id_divB];
	
	for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
	{
	
	memset(scr1,0,num_nodes_in_block*sizeof(D_REAL));
	
	for (INT32 i=1; i < BlkNds_X-1; i++)
	 for (INT32 j=1; j < BlkNds_Y-1; j++)
	  for (INT32 k=1; k < BlkNds_Z-1; k++)
		{
	
		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
	
		if(!Core_Flag[i_j_k])
			for (a=-1; a<=1; a++)
			for (b=-1; b<=1; b++)
			for (c=-1; c<=1; c++)
			{
	
			INT32 temp1 = comp*num_nodes_in_block +(i+a)*BlkNds_Y*BlkNds_Z +(j+b)*BlkNds_Z +(k+c);
			scr1[i_j_k] += smooth_coeff[a+1][b+1][c+1] * Field_Type[type][temp1];
			}
	
		}
	
	for(INT32 i=1; i < BlkNds_X-1; i++)
	 for(INT32 j=1; j < BlkNds_Y-1; j++)
	  for(INT32 k=1; k < BlkNds_Z-1; k++)
	  {
	
		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
	
		if(!Core_Flag[i_j_k])
		Field_Type[type][comp*num_nodes_in_block +i_j_k] = (1.-smooth_value) *
			Field_Type[type][comp*num_nodes_in_block +i_j_k]+ smooth_value * scr1[i_j_k];
	
	
	  }
	}

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;
}
#endif




//!-------------------------------------------------------------//
//! LF_B0toB1: -						//
//!-------------------------------------------------------------//
void CBlock::apply_BField_Boundaries(INT32 id_field)
{


	//! return when block is not simulation box boundary
	if(!is_box_boundary[_ANY_])
	return;



	//! ------------------------------------------
	//! --- Xmin Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Xmin_BB])
	 if(use_hom_B_bounds[Xmin_BB])
	 {

		if(!use_B_inflow_bounds[Xmin_BB])
		cp_imin_Boundaries(id_field);

	 }
	 else
	 {
		//! specify values to set
		INT32 start[3]  = {             0,        0,        0};
		INT32   end[3]  = {             1, BlkNds_Y, BlkNds_Z};
		INT32   src[3]  = {+NM_BOUND_NODE,        0,        0};

		set_inhom_B_field(id_field, start, end, src);

	 }



	//! ------------------------------------------
	//! --- Xmax Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Xmax_BB])
	 if(use_hom_B_bounds[Xmax_BB])
	 {

		if(!use_B_inflow_bounds[Xmax_BB])
		cp_imax_Boundaries(id_field);

	 }
	 else
	 {
		//! specify values to set
		INT32 start[3]  = {    BlkNds_X-1,        0,        0};
		INT32   end[3]  = {    BlkNds_X-0, BlkNds_Y, BlkNds_Z};
		INT32   src[3]  = {-NM_BOUND_NODE,        0,        0};

		set_inhom_B_field(id_field, start, end, src);

	 }


	//! ------------------------------------------
	//! --- Ymin Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Ymin_BB])
	 if(use_hom_B_bounds[Ymin_BB])
	 {

		if(!use_B_inflow_bounds[Ymin_BB])
		cp_jmin_Boundaries(id_field);

	 }
	 else
	 {

		//! specify values to set
		INT32 start[3]  = {       0,              0,        0};
		INT32   end[3]  = {BlkNds_X,              1, BlkNds_Z};
		INT32   src[3]  = {       0, +NM_BOUND_NODE,         0};

		set_inhom_B_field(id_field, start, end, src);
	 }

	//! ------------------------------------------
	//! --- Ymax Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Ymax_BB])
	 if(use_hom_B_bounds[Ymax_BB])
	 {

		if(!use_B_inflow_bounds[Ymax_BB])
		cp_jmax_Boundaries(id_field);

	 }
	 else
	 {

		//! specify values to set
		INT32 start[3]  = {       0,     BlkNds_Y-1,        0};
		INT32   end[3]  = {BlkNds_X,     BlkNds_Y-0, BlkNds_Z};
		INT32   src[3]  = {       0, -NM_BOUND_NODE,        0};

		set_inhom_B_field(id_field, start, end, src);

	 }


	//! ------------------------------------------
	//! --- Zmin Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Zmin_BB])
	 if(use_hom_B_bounds[Zmin_BB])
	 {

		if(!use_B_inflow_bounds[Zmin_BB])
		cp_kmin_Boundaries(id_field);

	 }
	 else
	 {

		//! specify values to set
		INT32 start[3]  = {       0,        0,              0};
		INT32   end[3]  = {BlkNds_X, BlkNds_Y,              1};
		INT32   src[3]  = {       0,        0, +NM_BOUND_NODE};

		set_inhom_B_field(id_field, start, end, src);

	 }

	//! ------------------------------------------
	//! --- Zmax Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Zmax_BB])
	 if(use_hom_B_bounds[Zmax_BB])
	 {

		if(!use_B_inflow_bounds[Zmax_BB])
		cp_kmax_Boundaries(id_field);

	 }
	 else
	 {

		//! specify values to set
		INT32 start[3]  = {       0,        0,     BlkNds_Z-1};
		INT32   end[3]  = {BlkNds_X, BlkNds_Y,     BlkNds_Z-0};
		INT32   src[3]  = {       0,        0, -NM_BOUND_NODE};

		set_inhom_B_field(id_field, start, end, src);

	 }



}


#if defined nonadiabatic_gradPE_TERM
//!-------------------------------------------------------------//
//! apply_PE_Boundaries: -						//
//!-------------------------------------------------------------//
void CBlock::apply_PE_Boundaries(INT32 id_field)
{


	//! return when block is not simulation box boundary
	if(!is_box_boundary[_ANY_])
	{
	  return;
	}


	//! ------------------------------------------
	//! --- Xmin Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Xmin_BB])
	{
	  if(!use_PE_inflow_bounds[Xmin_BB])
	  {
	    cp_imin_Boundaries(id_field);
	  }
	}
	
	
	//! ------------------------------------------
	//! --- Xmax Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Xmax_BB])
	{
	  if(!use_PE_inflow_bounds[Xmax_BB])
	  {
	    cp_imax_Boundaries(id_field);
	  }
	}
	
	
	//! ------------------------------------------
	//! --- Ymin Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Ymin_BB])
	{
	  if(!use_PE_inflow_bounds[Ymin_BB])
	  {
	    cp_jmin_Boundaries(id_field);
	  }
	}
	
	
	//! ------------------------------------------
	//! --- Ymax Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Ymax_BB])
	{
	  if(!use_PE_inflow_bounds[Ymax_BB])
	  {
	    cp_jmax_Boundaries(id_field);
	  }
	}
	
	
	//! ------------------------------------------
	//! --- Zmin Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Zmin_BB])
	{
	  if(!use_PE_inflow_bounds[Zmin_BB])
	  {
	    cp_kmin_Boundaries(id_field);
	  }
	}
	
	
	//! ------------------------------------------
	//! --- Zmax Boundary ------------------------
	//! ------------------------------------------
	if(is_box_boundary[Zmax_BB])
	{
	  if(!use_PE_inflow_bounds[Zmax_BB])
	  {
	    cp_kmax_Boundaries(id_field);
	  }
	}

}
#endif


//!------------------------------------------------------------------//
//!-------------- Leap Frog related Functions -----------------------//
//!------------------------------------------------------------------//

//!-------------------------------------------------------------//
//! LF_B0toB1: -						//
//!-------------------------------------------------------------//
void CBlock::LF_odd_B0toB1(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();



	D_REAL* B1_odd = Field_Type[id_BOdd];
	D_REAL* B1_even = Field_Type[id_BEven];
	
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	//! use id_rotB as temporary scratch
	D_REAL* PE_odd = Field_Type[id_PE_odd];
	D_REAL* PE_even = Field_Type[id_PE_even];
	
	//! equivalent comment as below
	memcpy(PE_odd, PE_even, num_nodes_in_block*sizeof(D_REAL));
#endif
	

	 //! at the beginnig B_Even and B_Odd have to be
   	 //! the average solution of last TL which is stored
   	 //! in B_Even. So first cp B_Even to B_Odd
  	 memcpy(B1_odd, B1_even,3*num_nodes_in_block*sizeof(D_REAL));

	//! id_BEven = BEven at FI Mesh
	apply_BField_Boundaries(id_BEven);
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	apply_PE_Boundaries(id_PE_even);
#endif

	//! apply boundaries to half integer mesh
	if(mesh_type==STAGGERED)
	apply_BField_Boundaries(id_B_HI);



#ifdef USE_CFBG_BFIELD


	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative_CFBG_BField(id_BOdd,
						   id_PE_odd,
						   id_BEven,
						   id_PE_even,
						   1.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative_CFBG_BField(id_BOdd,
						   id_BEven,
						   1.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BEven, id_B_HI);
		else
		calc_Faraday_EField(id_BEven);
	}

#else

	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative(id_BOdd,
					id_PE_odd,
					id_BEven,
					id_PE_even,
					1.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative(id_BOdd,
					id_BEven,
					1.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BEven, id_B_HI);
		else
		calc_Faraday_EField(id_BEven);
	}
#endif

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}

//!-------------------------------------------------------------//
//! advance_B_even: -							//
//!-------------------------------------------------------------//
void CBlock::LF_advance_B_even(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	//! apply boundaries NM condiotn if necessary
	apply_BField_Boundaries(id_BOdd);
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	apply_PE_Boundaries(id_PE_odd);
#endif

	//! apply boundaries to half integer mesh
	if(mesh_type==STAGGERED)
	apply_BField_Boundaries(id_B_HI);



#ifdef USE_CFBG_BFIELD

	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative_CFBG_BField(id_BEven,
						   id_PE_even,
						   id_BOdd,
						   id_PE_odd,
						   2.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative_CFBG_BField(id_BEven,
						   id_BOdd,
						   2.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BOdd, id_B_HI);
		else
		calc_Faraday_EField(id_BOdd);
	}

#else

	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative(id_BEven,
				       id_PE_even,
				       id_BOdd,
				       id_PE_odd,
				       2.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative(id_BEven,
				       id_BOdd,
				       2.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BOdd, id_B_HI);
		else
		calc_Faraday_EField(id_BOdd);
	}


#endif


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}

//!-------------------------------------------------------------//
//! advance_B_odd: -							//
//!-------------------------------------------------------------//
void CBlock::LF_advance_B_odd(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

	//! apply boundaries NM condiotn if necessary
	apply_BField_Boundaries(id_BEven);
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	apply_PE_Boundaries(id_PE_even);
#endif

	//! apply boundaries to half integer mesh
	if(mesh_type==STAGGERED)
	apply_BField_Boundaries(id_B_HI);





#ifdef USE_CFBG_BFIELD

	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative_CFBG_BField(id_BOdd,
						   id_PE_odd,
						   id_BEven,
						   id_PE_even,
						   2.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative_CFBG_BField(id_BOdd,
						   id_BEven,
						   2.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BEven, id_B_HI);
		else
		calc_Faraday_EField(id_BEven);
	}

#else

	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative(id_BOdd,
				       id_PE_odd,
				       id_BEven,
				       id_PE_even,
				       2.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative(id_BOdd,
				       id_BEven,
				       2.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BEven, id_B_HI);
		else
		calc_Faraday_EField(id_BEven);
	}

#endif



	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}

//!-------------------------------------------------------------//
//! LF_even_B_Solution: -							//
//!-------------------------------------------------------------//
void CBlock::LF_even_B_Solution(void)
{
	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	//! apply boundaries NM condiotn if necessary
	apply_BField_Boundaries(id_BOdd);
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	apply_PE_Boundaries(id_PE_odd);
#endif

	//! apply boundaries to half integer mesh
	if(mesh_type==STAGGERED)
	apply_BField_Boundaries(id_B_HI);


	//! advance B_Even(2N) to  B_(2N+1) "Even Solution (at odd TL)"


#ifdef USE_CFBG_BFIELD


	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative_CFBG_BField(id_BEven,
						   id_PE_even,
						   id_BOdd,
						   id_PE_odd,
						   1.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative_CFBG_BField(id_BEven,
						   id_BOdd,
						   1.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BOdd, id_B_HI);
		else
		calc_Faraday_EField(id_BEven);
	}


#else

	if(LF_one_step)
	{

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
		LF_B_Step_conservative(id_BEven,
				       id_PE_even,
				       id_BOdd,
				       id_PE_odd,
				       1.*dt_field_of_L[RLevel]);
#else
		LF_B_Step_conservative(id_BEven,
				       id_BOdd,
				       1.*dt_field_of_L[RLevel]);
#endif
	}
	else
	{
		if(mesh_type==STAGGERED)
		calc_Faraday_EField_staggered(id_BOdd, id_B_HI);
		else
		calc_Faraday_EField(id_BEven);
	}

#endif


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}




//!-------------------------------------------------------------//
//! LF_average_B_Solution: -							//
//!-------------------------------------------------------------//
bool CBlock::LF_average_B_Solution(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	bool NaN_in_BField = false;

	D_REAL* B1_odd   = Field_Type[id_BOdd];
	D_REAL* B1_even  = Field_Type[id_BEven];
	

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	//! id_PE_even should be id_PESpecies1
	D_REAL* PE_even = Field_Type[id_PE_even];
	D_REAL* PE_odd  = Field_Type[id_PE_odd ];
#endif
	
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	for (INT32 elm=0; elm < num_nodes_in_block; ++elm)
	{
		//! Average solution to PE
		PE_even[elm] = 0.5 * ( PE_even[elm] + PE_odd[elm] );
		
		//! Check whether solution diverges
		if( 0.*PE_even[elm]!=0.)
		{
			NaN_in_BField = true;
		}
	}
#endif
	
	for (INT32 elm=0; elm < 3*num_nodes_in_block; elm++)
	{
	
		//! Average fields to B_Even
		B1_even[elm]  = 0.5 *(B1_even[elm] + B1_odd[elm]);
		
	
		//! Check whether solution diverges
		if( fabs(B1_even[elm]) > 1.e8 || 0.*B1_even[elm]!=0.)
		{
// 			log_file << "NaN in B at:" << endl
// 				 << "B: " << B1_even[elm] << endl
// 				 << "RLevel: " <<RLevel<<" " << endl
// 				 << "Cell: (" <<elm<<") " << endl
// 				 << "Block: (" <<Blk_Index[0]<<","
// 						<<Blk_Index[1]<<","
// 						<<Blk_Index[2]<<")" << endl;

			NaN_in_BField = true;
	
		}
	
	}


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

	return NaN_in_BField;
}

//!-------------------------------------------------------------//
//! build_B_total: -							  //
//!-------------------------------------------------------------//
void CBlock::build_B_total(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


#ifdef USE_CFBG_BFIELD

      D_REAL* B1_even  = Field_Type[id_BEven];
      D_REAL* B1_cfbg  = Field_Type[id_Bcfbg];
      D_REAL* B1_total = Field_Type[id_BTotal];

      
#if defined(use_vectorclass)
	VEC4_D_REAL vec_B1_even;
	VEC4_D_REAL vec_B1_cfbg;
	for (INT32 elm=0; elm < 3*num_nodes_in_block; elm+=4)
	{
	  vec_B1_even.load(B1_even+elm);
	  vec_B1_cfbg.load(B1_cfbg+elm);
	  vec_B1_even += vec_B1_cfbg;
	  vec_B1_even.store(B1_total+elm);
	}
#else
      //! Add Backround Field and solution to B_Total
      for (INT32 elm=0; elm < 3*num_nodes_in_block; elm++)
      B1_total[elm] = B1_even[elm] + B1_cfbg[elm];
#endif
      
      

#else

      D_REAL* B1_even  = Field_Type[id_BEven];
      D_REAL* B1_total = Field_Type[id_BTotal];

      memcpy(B1_total, B1_even, 3*num_nodes_in_block*sizeof(D_REAL));

#endif


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}




#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
void CBlock::explicit_midpoint_method_for_Pe_1(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

	
	D_REAL* PE_odd = Field_Type[id_PE_odd];
	D_REAL* PE_even = Field_Type[id_PE_even];
	
	//! at the beginnig PE_even and PE_odd have to be
   	 //! the solution of last TL which is stored
   	 //! in PE_even. So first cp PE_even to PE_odd
	memcpy(PE_odd, PE_even, num_nodes_in_block*sizeof(D_REAL));
	

	apply_PE_Boundaries(id_PE_even);

#ifdef USE_CFBG_BFIELD
	LF_Pe_Step_conservative_CFBG_BField(id_PE_odd,
					    id_PE_even,
					    0.5*static_cast<D_REAL>(NUM_SUB_CYCLE)*dt_field_of_L[RLevel]);
#else
	LF_Pe_Step_conservative(id_PE_odd,
				id_PE_even,
				0.5*static_cast<D_REAL>(NUM_SUB_CYCLE)*dt_field_of_L[RLevel]);
#endif

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;
}

bool CBlock::explicit_midpoint_method_for_Pe_2(void)
{
	bool NaN_in_BField = false;

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

	
	D_REAL* PE_odd = Field_Type[id_PE_odd];
	D_REAL* PE_even = Field_Type[id_PE_even];
	
	
	
	apply_PE_Boundaries(id_PE_odd);

#ifdef USE_CFBG_BFIELD
	LF_Pe_Step_conservative_CFBG_BField(id_PE_even,
					    id_PE_odd,
					    1.0*static_cast<D_REAL>(NUM_SUB_CYCLE)*dt_field_of_L[RLevel]);
#else
	LF_Pe_Step_conservative(id_PE_even,
				id_PE_odd,
				1.0*static_cast<D_REAL>(NUM_SUB_CYCLE)*dt_field_of_L[RLevel]);
#endif
	
	//! Check whether solution diverges
	for (INT32 elm=0; elm < num_nodes_in_block; ++elm)
	{
		if( 0.*PE_even[elm]!=0.)
		{
			NaN_in_BField = true;
		}
	}

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;
	
	return NaN_in_BField;
}
#endif



//!------------------------------------------------------------------//
//!-------------- Runge Kutta related Functions ---------------------//
//!------------------------------------------------------------------//

//!-------------------------------------------------------------//
//! RK_sequence: -							//
//!-------------------------------------------------------------//
#ifndef nonadiabatic_gradPE_TERM
void CBlock::RK_sequence(void)
{

	//! K1
	RK_advance_B(false,0.,(1./6.));

	//! K2
	RK_advance_B(true,0.5,(2./6.));

	//! K3
	RK_advance_B(true,0.5,(2./6.));
	 
	//! K4
	RK_advance_B(true,1.,(1./6.));



}
#endif /* nonadiabatic_gradPE_TERM */
//!-------------------------------------------------------------//
//! RK_advance_B: -							//
//!-------------------------------------------------------------//
#ifndef nonadiabatic_gradPE_TERM
bool CBlock::RK_advance_B(bool calc_B_dev, D_REAL c1, D_REAL c2)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();



	bool NaN_in_BField = false;
        INT32 elmt;
	D_REAL dt_block = dt_field_of_L[RLevel];

	D_REAL *B_old = Field_Type[id_BOld];
	D_REAL *B_new = Field_Type[id_BNew];

	D_REAL *B_dev = Field_Type[id_BDerivative];
	D_REAL *K_X   = Field_Type[id_KX];

	//! ---- Calc B_dev = B+dt*KX_1_dest) if necessary ---------
	if(!calc_B_dev)
	{

		//! if function is called first time the last "new" solution
		//! to be set to old
		memset(K_X,0,3*num_nodes_in_block*sizeof(D_REAL));
		memcpy(B_old, B_new, 3*num_nodes_in_block*sizeof(D_REAL));
		memcpy(B_dev, B_new, 3*num_nodes_in_block*sizeof(D_REAL));


	}
	else
	{
		//! K_X GN must be update now
		for(elmt=0; elmt<3*num_nodes_in_block;elmt++)
        	B_dev[elmt] = B_old[elmt] + c1 *dt_block *K_X[elmt];

	}

	//! ---- Calc KX_dest = f(B_dev) ---------------------------
	RK_B_Step_conservative_CFBG_BField();


	//! ---- Calc id_BNew = c2 *dt KX_dest------------------
	if(c1 != 1.)
	for(elmt=0; elmt<3*num_nodes_in_block;elmt++)
	B_new[elmt] += c2 *dt_block *K_X[elmt];
	//! just check in last "RK_advance_B" for NaN values
	else
	for(elmt=0; elmt<3*num_nodes_in_block;elmt++)
	{
		B_new[elmt] += c2 *dt_block *K_X[elmt];
		//! Check whether solution diverges
		if( fabs(B_new[elmt]) > 1.e8 || 0.*B_new[elmt]!=0.)
		{
// 			log_file << "NaN in B at:" << endl
// 				 << "B: " << B_new[elmt] << endl
// 				 << "RLevel: " <<RLevel<<" " << endl
// 				 << "Cell: (" <<elmt<<") " << endl
// 				 << "Block: (" <<Blk_Index[0]<<","
// 						<<Blk_Index[1]<<","
// 						<<Blk_Index[2]<<")" << endl;

			NaN_in_BField = true;
	
		}
	}


	//! apply boundaries NM condiotn if necessary
	apply_BField_Boundaries(id_BNew);



	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


	return NaN_in_BField;


}
#endif /* nonadiabatic_gradPE_TERM */





//!-------------------------------------------------------------//
//! cp_imin_Boundaries: -							//
//!-------------------------------------------------------------//
 void CBlock::cp_imin_Boundaries(INT32 field)
{

	//! source is ultimate(penultimate) CS with physical fields:
	const INT32 i_source = NM_BOUND_NODE * BlkNds_Y * BlkNds_Z;

	//! cp to last CS
	const INT32 i_dest   = 0;

	//! size of X Cross Section 
	const INT32 XCS_size =  BlkNds_Y *BlkNds_Z *sizeof(D_REAL);


	for(INT32 comp=0; comp<COMPs_FType[field]; comp++)
	memcpy(Field_Type[field] +comp*num_nodes_in_block + i_dest,
	       Field_Type[field] +comp*num_nodes_in_block + i_source,
	       XCS_size);


}


//!-------------------------------------------------------------//
//! cp_imax_Boundaries: -							//
//!-------------------------------------------------------------//
void CBlock::cp_imax_Boundaries(INT32 field)
{



	//! source is ultimate(penultimate) CS with physical fields:
	const INT32 i_source = (BlkNds_X-1 -NM_BOUND_NODE) * BlkNds_Y * BlkNds_Z;

	//! cp to last CS
	const INT32 i_dest   = (BlkNds_X-1) * BlkNds_Y * BlkNds_Z;

	//! size of X Cross Section 
	const INT32 XCS_size =  BlkNds_Y *BlkNds_Z *sizeof(D_REAL);


	for(INT32 comp=0; comp<COMPs_FType[field]; comp++)
	memcpy(Field_Type[field] +comp*num_nodes_in_block + i_dest,
	       Field_Type[field] +comp*num_nodes_in_block + i_source,
	       XCS_size);


}


//!-------------------------------------------------------------//
//! cp_jmin_Boundaries: -							//
//!-------------------------------------------------------------//
void CBlock::cp_jmin_Boundaries(INT32 field)
{


	//! source is ultimate(penultimate) CS with physical fields:
	const INT32 k_row_source =  NM_BOUND_NODE *BlkNds_Z;

	const INT32 k_row_dest = 0;

	//! size of a k-row
	const INT32 k_row_size =  BlkNds_Z *sizeof(D_REAL);

	//! this loop copies the k-rows
	for(INT32 comp=0; comp<COMPs_FType[field]; comp++)
	 for(INT32 a=0; a<BlkNds_X; a++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		INT32 XCS_pos = a *BlkNds_Y  *BlkNds_Z;
	
		memcpy(Field_Type[field] +comp*num_nodes_in_block +XCS_pos +k_row_dest,
		       Field_Type[field] +comp*num_nodes_in_block +XCS_pos +k_row_source,
		       k_row_size);
	 }	

}

//!-------------------------------------------------------------//
//! cp_jmax_Boundaries: -							//
//!-------------------------------------------------------------//
void CBlock::cp_jmax_Boundaries(INT32 field)
{


	//! source is ultimate(penultimate) CS with physical fields:
	const INT32 k_row_source =  (BlkNds_Y-1 -NM_BOUND_NODE) *BlkNds_Z;

	const INT32 k_row_dest = (BlkNds_Y-1)*BlkNds_Z;

	//! size of a k-row
	const INT32 k_row_size =  BlkNds_Z *sizeof(D_REAL);

	//! this loop copies the k-rows
	for(INT32 comp=0; comp<COMPs_FType[field]; comp++)
	 for(INT32 a=0; a<BlkNds_X; a++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		INT32 XCS_pos = a *BlkNds_Y  *BlkNds_Z;
	
		memcpy(Field_Type[field] +comp*num_nodes_in_block +XCS_pos +k_row_dest,
		       Field_Type[field] +comp*num_nodes_in_block +XCS_pos +k_row_source,
		       k_row_size);
	 }	

}

//!-------------------------------------------------------------//
//! cp_kmin_Boundaries: -							//
//!-------------------------------------------------------------//
void CBlock::cp_kmin_Boundaries(INT32 field)
{



	D_REAL *Field_Dest = Field_Type[field];
	D_REAL *Field_Src = Field_Type[field];


	//! source is ultimate(penultimate) CS with physical fields:
	const INT32 z_node_src = NM_BOUND_NODE;

	const INT32 z_node_dest = 0;

	//! this loop copies the cells

	for(INT32 comp=0; comp<COMPs_FType[field]; comp++)
	 for(INT32 a=0; a<BlkNds_X-0; a++)
	 {

		INT32 XCS_pos = a *BlkNds_Y *BlkNds_Z;
		for(INT32 b=0; b<BlkNds_Y-0; b++)
		{
			INT32 Y_row = b *BlkNds_Z;

			  Field_Dest[comp*num_nodes_in_block +XCS_pos + Y_row + z_node_dest] = 
			   Field_Src[comp*num_nodes_in_block +XCS_pos + Y_row + z_node_src];
		}
	 }
}

//!-------------------------------------------------------------//
//! cp_kmax_Boundaries: -							//
//!-------------------------------------------------------------//
void CBlock::cp_kmax_Boundaries(INT32 field)
{



	D_REAL *Field_Dest = Field_Type[field];
	D_REAL *Field_Src = Field_Type[field];


	//! source is ultimate(penultimate) CS with physical fields:
	const INT32 z_node_src = (BlkNds_Z-1 -NM_BOUND_NODE);

	const INT32 z_node_dest = (BlkNds_Z-1);

	//! this loop copies the cells

	for(INT32 comp=0; comp<COMPs_FType[field]; comp++)
	 for(INT32 a=0; a<BlkNds_X-0; a++)
	 {

		INT32 XCS_pos = a *BlkNds_Y *BlkNds_Z;
		for(INT32 b=0; b<BlkNds_Y-0; b++)
		{
			INT32 Y_row = b *BlkNds_Z;

			  Field_Dest[comp*num_nodes_in_block +XCS_pos + Y_row + z_node_dest] = 
			   Field_Src[comp*num_nodes_in_block +XCS_pos + Y_row + z_node_src];
		}
	 }
}


//!-------------------------------------------------------------//
//! cp_Boundaries: -							//
//!-------------------------------------------------------------//
void CBlock::cp_Boundaries(INT32 field, const bool *use_inflow)
{

	//! return in case this block is no box boundary
	if(!is_box_boundary[_ANY_])
	return;

	if(is_box_boundary[0] && !use_inflow[0])
	cp_imin_Boundaries(field);

	if(is_box_boundary[1] && !use_inflow[1])
	cp_imax_Boundaries(field);

	if(is_box_boundary[2] && !use_inflow[2])
	cp_jmin_Boundaries(field);

	if(is_box_boundary[3] && !use_inflow[3])
	cp_jmax_Boundaries(field);

	if(is_box_boundary[4] && !use_inflow[4])
	cp_kmin_Boundaries(field);

	if(is_box_boundary[5] && !use_inflow[5])
	cp_kmax_Boundaries(field);

}



//!-------------------------------------------------------------//
//! set_inhom_B_field: -							//
//!-------------------------------------------------------------//
 void CBlock::set_inhom_B_field(INT32 id_field, INT32 *start, INT32 *end, INT32 *src)
{

		D_REAL *BX = Field_Type[id_field] +0*num_nodes_in_block;
		D_REAL *BY = Field_Type[id_field] +1*num_nodes_in_block;
		D_REAL *BZ = Field_Type[id_field] +2*num_nodes_in_block;

		INT32 ind[3];

		PARTICLE_REAL x[3];
		D_REAL B_init[3] = {0., 0., 0.};
		PARTICLE_REAL v_init[3] = {0., 0., 0.};



		for(ind[0] = start[0]; ind[0] < end[0]; ind[0]++)
		 for(ind[1] = start[1]; ind[1] < end[1]; ind[1]++)
		  for(ind[2] = start[2]; ind[2] < end[2]; ind[2]++)
		  {
			
			INT32 dest_i_j_k  =   ind[0]*BlkNds_Y*BlkNds_Z 
					    +ind[1]*BlkNds_Z 
					    +ind[2];

			INT32 src_i_j_k  =   (ind[0] +src[0])*BlkNds_Y*BlkNds_Z 
					   +(ind[1] +src[1])*BlkNds_Z 
					   +(ind[2] +src[2]);


			//! estimate coordinate of cell's centre
			cell_centre_normedCoords(x, ind);

			//! choose velocity depending on coordinate
			if(set_BFieldInflow_at_VelocityInflow)
			set_inflow_velocity(v_init, x, 0);

		
			//! - copy cell when velocity is directed towards BB
			//! - in case of min BB src is -1
			//! - in case of max BB src is +1
			//!   ->   equal sign:  inflow -> set
			//!   -> oposite sign: outflow -> copy
			if(set_BFieldInflow_at_VelocityInflow && src[0]!=0 && v_init[0] *src[0]<0.)
			{
				BX[dest_i_j_k] = BX[src_i_j_k];
				BY[dest_i_j_k] = BY[src_i_j_k];
				BZ[dest_i_j_k] = BZ[src_i_j_k];
			}
			else
			{

				set_inflow_BField(B_init, x);

				BX[dest_i_j_k] = B_init[0];
				BY[dest_i_j_k] = B_init[1];
				BZ[dest_i_j_k] = B_init[2];
			}
		  }

}





