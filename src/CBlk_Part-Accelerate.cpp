
#include "CBlk.h"
#include "parameters.h"
#include "utils.h"
#include "NeutralDrag.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>




extern D_REAL *q_m_ratio, *dt_particle_of_L;


//! pointers for EM Fields
extern D_REAL *BX, *BY, *BZ;
extern D_REAL *EX, *EY, *EZ;

// inline PARTICLE_REAL vec_scalar_inl(PARTICLE_REAL *v1, PARTICLE_REAL *v2)
// { 
//   return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
// }

// inline D_REAL vec_scalar_inl(PARTICLE_REAL *v1, D_REAL *v2)
// { 
//   return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
// }

inline D_REAL vec_scalar_inl(D_REAL *v1, D_REAL *v2)
{ 
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

//!-----------------------------------------------------------
//! accelerate_particle
//!-----------------------------------------------------------
void CBlock::accelerate_particle(INT32 id_oct)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	INT32 comp, ind[3];
	particle *active_particle;


	//! double precision is used for to reasons:
	//! - it is more precise
	//! - it is faster to typecast the pos and velocities into
	//!   double once and do the remaining calculations in
	//!   double precision
	D_REAL r[3], v[3], coord[3], E[3], B[3];
	D_REAL h,K, sq_B, E_B, v_B, E_x_B[3], v_x_B[3];
	D_REAL shape_func[8], E_corner[24], B_corner[24];
	
	BX = Field_Type[id_BTotal];
	BY = BX +num_nodes_in_block;
	BZ = BY +num_nodes_in_block;
	
	EX = Field_Type[id_EField];
	EY = EX +num_nodes_in_block;
  	EZ = EY +num_nodes_in_block;


	//! get index of oct (0 or 1 for each dimension)
	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);
	

	//! START:
	//! eg: BlkNds=6
	//! 1)
	//! a=0
	//! -> start = 1
	//! 2)
	//! a=1
	//! -> start = 3
	const INT32 oct_start[3] = { 1 +a*(BlkNds_X/2-1),
				    1 +b*(BlkNds_Y/2-1),
				    1 +c*(BlkNds_Z/2-1)};


	//! END:
	//! eg: BlkNds=6
	//! a=0
	//! -> end = 3
	//! a=1
	//! -> end = 5
	const INT32 oct_end[3]={ BlkNds_X/2 +a*(BlkNds_X/2-1),
				BlkNds_Y/2 +b*(BlkNds_Y/2-1),
				BlkNds_Z/2 +c*(BlkNds_Z/2-1)};

	for(INT32 i = oct_start[0]; i < oct_end[0]; i++)
	 for(INT32 j = oct_start[1]; j < oct_end[1]; j++)
	  for(INT32 k = oct_start[2]; k < oct_end[2]; k++)
	  {
	

		
		//! ------------------------------------------
		INT32 i_j_k   =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		
		//! ------------------------------------------
		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		
		//! -------------------------------------------
		INT32 ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
		
		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);



		//! store EM Fields in local Arrays for faster cache acces

		//! ----- EX ----------------
		E_corner[0 +0*8] = EX[  i_j_k];
		E_corner[1 +0*8] = EX[ip1_j_k];
		E_corner[2 +0*8] = EX[i_jp1_k];
		E_corner[3 +0*8] = EX[i_j_kp1];
		E_corner[4 +0*8] = EX[  ip1_jp1_k];
		E_corner[5 +0*8] = EX[  ip1_j_kp1];
		E_corner[6 +0*8] = EX[  i_jp1_kp1];
		E_corner[7 +0*8] = EX[ip1_jp1_kp1];

		//! ----- EY ----------------
		E_corner[0 +1*8] = EY[  i_j_k];
		E_corner[1 +1*8] = EY[ip1_j_k];
		E_corner[2 +1*8] = EY[i_jp1_k];
		E_corner[3 +1*8] = EY[i_j_kp1];
		E_corner[4 +1*8] = EY[  ip1_jp1_k];
		E_corner[5 +1*8] = EY[  ip1_j_kp1];
		E_corner[6 +1*8] = EY[  i_jp1_kp1];
		E_corner[7 +1*8] = EY[ip1_jp1_kp1];

		//! ----- EZ ----------------
		E_corner[0 +2*8] = EZ[  i_j_k];
		E_corner[1 +2*8] = EZ[ip1_j_k];
		E_corner[2 +2*8] = EZ[i_jp1_k];
		E_corner[3 +2*8] = EZ[i_j_kp1];
		E_corner[4 +2*8] = EZ[  ip1_jp1_k];
		E_corner[5 +2*8] = EZ[  ip1_j_kp1];
		E_corner[6 +2*8] = EZ[  i_jp1_kp1];
		E_corner[7 +2*8] = EZ[ip1_jp1_kp1];

		//! ----- BX ----------------
		B_corner[0 +0*8] = BX[  i_j_k];
		B_corner[1 +0*8] = BX[ip1_j_k];
		B_corner[2 +0*8] = BX[i_jp1_k];
		B_corner[3 +0*8] = BX[i_j_kp1];
		B_corner[4 +0*8] = BX[  ip1_jp1_k];
		B_corner[5 +0*8] = BX[  ip1_j_kp1];
		B_corner[6 +0*8] = BX[  i_jp1_kp1];
		B_corner[7 +0*8] = BX[ip1_jp1_kp1];

		//! ----- BY ----------------
		B_corner[0 +1*8] = BY[  i_j_k];
		B_corner[1 +1*8] = BY[ip1_j_k];
		B_corner[2 +1*8] = BY[i_jp1_k];
		B_corner[3 +1*8] = BY[i_j_kp1];
		B_corner[4 +1*8] = BY[  ip1_jp1_k];
		B_corner[5 +1*8] = BY[  ip1_j_kp1];
		B_corner[6 +1*8] = BY[  i_jp1_kp1];
		B_corner[7 +1*8] = BY[ip1_jp1_kp1];

		//! ----- BZ ----------------
		B_corner[0 +2*8] = BZ[  i_j_k];
		B_corner[1 +2*8] = BZ[ip1_j_k];
		B_corner[2 +2*8] = BZ[i_jp1_k];
		B_corner[3 +2*8] = BZ[i_j_kp1];
		B_corner[4 +2*8] = BZ[  ip1_jp1_k];
		B_corner[5 +2*8] = BZ[  ip1_j_kp1];
		B_corner[6 +2*8] = BZ[  i_jp1_kp1];
		B_corner[7 +2*8] = BZ[ip1_jp1_kp1];
	
		

		//! even though using species as inner loop is not perfectly
		//! cache efficient, just do it here to save all the field
		//! interpolation for each cell
		for(INT32 species=0; species<num_Charged_Species; species++)
		 for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
		 {

			active_particle = pArray[species][i_j_k] +part_index;
		
		
			//! copy ponter only -> few percent faster
			v[0] = active_particle->v[0];
			v[1] = active_particle->v[1];
			v[2] = active_particle->v[2];


			r[0] = active_particle->rel_r[0];
			r[1] = active_particle->rel_r[1];
			r[2] = active_particle->rel_r[2];
		

			

			//! --- now get shape function for interpolation ---
			shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2]);
			shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2]);
			shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2]);
			
			shape_func[4] = (   r[0])*(   r[1])*(1.-r[2]);
			shape_func[5] = (   r[0])*(1.-r[1])*(   r[2]);
			shape_func[6] = (1.-r[0])*(   r[1])*(   r[2]);
			shape_func[7] = (   r[0])*(   r[1])*(   r[2]);
		
		
			//! --- INTERPOLATION ----------------------------------------
			//! -- Electrical Field -------------------
			E[0] =   E_corner[0 +0*8] * shape_func[0]
				+E_corner[1 +0*8] * shape_func[1]
				+E_corner[2 +0*8] * shape_func[2]
				+E_corner[3 +0*8] * shape_func[3]
		
				+E_corner[4 +0*8] * shape_func[4]
				+E_corner[5 +0*8] * shape_func[5]
				+E_corner[6 +0*8] * shape_func[6]
				+E_corner[7 +0*8] * shape_func[7];
		
			E[1] =   E_corner[0 +1*8] * shape_func[0]
				+E_corner[1 +1*8] * shape_func[1]
				+E_corner[2 +1*8] * shape_func[2]
				+E_corner[3 +1*8] * shape_func[3]
		
				+E_corner[4 +1*8] * shape_func[4]
				+E_corner[5 +1*8] * shape_func[5]
				+E_corner[6 +1*8] * shape_func[6]
				+E_corner[7 +1*8] * shape_func[7];
		
			E[2] =   E_corner[0 +2*8] * shape_func[0]
				+E_corner[1 +2*8] * shape_func[1]
				+E_corner[2 +2*8] * shape_func[2]
				+E_corner[3 +2*8] * shape_func[3]
	
				+E_corner[4 +2*8] * shape_func[4]
				+E_corner[5 +2*8] * shape_func[5]
				+E_corner[6 +2*8] * shape_func[6]
				+E_corner[7 +2*8] * shape_func[7];
		
			//! -- Magnetic Field --------------------
			B[0] =   B_corner[0 +0*8] * shape_func[0]
				+B_corner[1 +0*8] * shape_func[1]
				+B_corner[2 +0*8] * shape_func[2]
				+B_corner[3 +0*8] * shape_func[3]
		
				+B_corner[4 +0*8] * shape_func[4]
				+B_corner[5 +0*8] * shape_func[5]
				+B_corner[6 +0*8] * shape_func[6]
				+B_corner[7 +0*8] * shape_func[7];
		
			B[1] =   B_corner[0 +1*8] * shape_func[0]
				+B_corner[1 +1*8] * shape_func[1]
				+B_corner[2 +1*8] * shape_func[2]
				+B_corner[3 +1*8] * shape_func[3]
		
				+B_corner[4 +1*8] * shape_func[4]
				+B_corner[5 +1*8] * shape_func[5]
				+B_corner[6 +1*8] * shape_func[6]
				+B_corner[7 +1*8] * shape_func[7];
		
			B[2] =   B_corner[0 +2*8] * shape_func[0]
				+B_corner[1 +2*8] * shape_func[1]
				+B_corner[2 +2*8] * shape_func[2]
				+B_corner[3 +2*8] * shape_func[3]
	
				+B_corner[4 +2*8] * shape_func[4]
				+B_corner[5 +2*8] * shape_func[5]
				+B_corner[6 +2*8] * shape_func[6]
				+B_corner[7 +2*8] * shape_func[7];


			
			/*
			//!----------------------------------------
			//! ANALYTIC DIPOL FIELD
			D_REAL r_vec[3];


			INT32 cell_indices[3] = {i,j,k};
			intern2normedCoords(r_vec, r, cell_indices);

			D_REAL u[3] = {1. ,0. ,0.};



			D_REAL M_r = (   r_vec[0] *Magnetic_Moment[0]
					+r_vec[1] *Magnetic_Moment[1]
					+r_vec[2] *Magnetic_Moment[2]);

			D_REAL r = vec_len(r_vec);
	
			D_REAL r_3 = r*r*r;
			D_REAL r_5 = r*r*r*r*r;
	

			B[0] = (3.*(M_r)*r_vec[0])/r_5 -Magnetic_Moment[0]/r_3;
			B[1] = (3.*(M_r)*r_vec[1])/r_5 -Magnetic_Moment[1]/r_3;
			B[2] = (3.*(M_r)*r_vec[2])/r_5 -Magnetic_Moment[2]/r_3;

			vec_cross(E,B,u);*/
			//!----------------------------------------



			//! now the new velocities are calculated
			//! first introduce some abbreviations
			E_B =  vec_scalar_inl(E,B);
			v_B =  vec_scalar_inl(v,B);
			sq_B = vec_scalar_inl(B,B);

			//! inlining the code here make code run slower for some reason
			//! (ATHLON 3000)
			vec_cross(E_x_B,E,B);
			vec_cross(v_x_B,v,B);
		
		
			h = dt_particle_of_L[RLevel] * q_m_ratio[species];
			K = 1./(1.+0.25*h*h*sq_B);


			v[0]=  K * (
					(1.-0.25*h*h*sq_B)*v[0] + h*(v_x_B[0]+E[0])
					+0.5*h*h*(E_x_B[0] + v_B * B[0]) + 0.25*h*h*h*E_B*B[0]
					);

			v[1]=  K * (
					(1.-0.25*h*h*sq_B)*v[1] + h*(v_x_B[1]+E[1])
					+0.5*h*h*(E_x_B[1] + v_B * B[1]) + 0.25*h*h*h*E_B*B[1]
					);

			v[2]=  K * (
					(1.-0.25*h*h*sq_B)*v[2] + h*(v_x_B[2]+E[2])
					+0.5*h*h*(E_x_B[2] + v_B * B[2]) + 0.25*h*h*h*E_B*B[2]
					);
		


#ifdef NEUTRAL_DRAG
			//! calc drag force
			if(NeutralDrag_Species[species])
			{
		
				ind[0]=i; ind[1]=j; ind[2]=k;
		
				intern2normedCoords(coord, r, ind);
				apply_NeutralDrag(coord, v, dt_particle_of_L[RLevel]);
		
			}
#endif

			//! copy new velocities back to particle
			active_particle->v[0] = v[0];
			active_particle->v[1] = v[1];
			active_particle->v[2] = v[2];

			num_accelerated++;


		 }//! end for species /cell
	
	
	}//! end for ijk


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;



}
