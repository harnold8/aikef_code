
#ifndef CBLOCK_GLOBALS_H
#define CBLOCK_GLOBALS_H


// #include <stdio.h>
#include <mpi.h>

//!------------------------------------------------------------
//!---------- Variables which have to be initialized ----------
//!------------------------------------------------------------

D_REAL ***smooth_coeff;

//!------------------------------------------------------------
//!---------- indices for derivatives and interpolation -------
//!------------------------------------------------------------
//INT32 a, b, c; //! also negetive numbers are used for a,b,c !!!
INT32 temp1, temp2;

INT32 *temp_i, *temp_j, *temp_k;

INT32 i,j,k, ip1_jp1_kp1;
INT32 i_j_k, ip1_j_k, im1_j_k, i_jp1_k , i_jm1_k ,i_j_kp1 ,i_j_km1;

//! --- indices for mixed derivatives ---------------------
INT32 ip1_j_kp1, ip1_j_km1, im1_j_kp1, im1_j_km1;
INT32 ip1_jp1_k, ip1_jm1_k, im1_jp1_k, im1_jm1_k;
INT32 i_jp1_kp1, i_jp1_km1, i_jm1_kp1, i_jm1_km1;

//! particle related
INT32 *statistics_of;



D_REAL *q_of, *q2m_of, *q_m_ratio;
D_REAL *dt_field_of_L, *rdt_field_of_L, *dt_particle_of_L;

D_REAL **rd_of_L, **rd2_of_L, **delta_of_L, *CellVol_of_L;
D_REAL *rd,*rd2;



//! pointers for EM Fields
D_REAL *BX, *BY, *BZ;
D_REAL *EX, *EY, *EZ;
D_REAL *UX, *UY, *UZ;
D_REAL *eta;

//! MPI related variables
MPI::Datatype PARTICLE_MPI;

#endif
