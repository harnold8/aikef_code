





#include "CBlk.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>


using namespace std;
extern ofstream log_file;

extern INT32  *COMPs_FType;
extern D_REAL ***smooth_coeff;

extern INT32 a, b, c;
extern INT32 i,j,k, ip1_jp1_kp1;
extern INT32 i_j_k, ip1_j_k, im1_j_k, i_jp1_k , i_jm1_k ,i_j_kp1 ,i_j_km1;

//! --- indices for mixed derivatives & loop ----------------
extern INT32 temp1, temp2;
extern INT32 ip1_j_kp1, ip1_j_km1, im1_j_kp1, im1_j_km1;
extern INT32 ip1_jp1_k, ip1_jm1_k, im1_jp1_k, im1_jm1_k;
extern INT32 i_jp1_kp1, i_jp1_km1, i_jm1_kp1, i_jm1_km1;

extern D_REAL *dt_field_of_L, *dt_particle_of_L, **rd_of_L, **rd2_of_L;

extern D_REAL *rd,*rd2;



//!------------------------------------------------------------
//!---------- CBlock: Global Field Solver Variables -----------
//!------------------------------------------------------------

//!----------------- rho ------------------------------
D_REAL rP,dxP,dyP,dzP;

//!----------------- Eta ------------------------------
D_REAL Eta;
D_REAL dxEta, dyEta, dzEta;

//!----------------- U --------------------------------
D_REAL ux,uy,uz;
D_REAL dxUX,dxUY,dxUZ,dyUX,dyUY,dyUZ,dzUX,dzUY,dzUZ;

//!----------------- B --------------------------------
D_REAL bx,by,bz;

//!------ first derivatives B -------------------------
D_REAL dxBX,dxBY,dxBZ,dyBX,dyBY,dyBZ,dzBX,dzBY,dzBZ;

//!------ second derivatives B --(just mixed occure)---
D_REAL d2xBY,d2xBZ,d2yBX,d2yBZ,d2zBX,d2zBY;

//!------ mixed derivatives B ------------------------------
D_REAL dxdzBX,dxdyBX,dydzBX,dydxBY,dzdxBY,dydzBY,dzdxBZ,dydxBZ,dzdyBZ;

//! pointers to calc derivatives

D_REAL *E1,*E2,*E3;
D_REAL *B1,*B2,*B3;
D_REAL *U1,*U2,*U3;
D_REAL *ETA, *rRHO;


D_REAL *B1_cfbg;
D_REAL *B1_sw,    *B2_sw,    *B3_sw;
D_REAL *B1_total, *B2_total, *B3_total;




//!--------------------------------------------------------------
//!----------- derivatives for CONSERVATIVE FORM ----------------
//!-------------- (= using no product rule) ---------------------
//!--------------------------------------------------------------

//! In order to distinguish between conservative and not conservative
//! form, the conservative derivatives are indexed with 1,2,3, whereas
//! the alternative formulation uses X,Y,Z

//! The conservative form better conserves the divergence free condition
//! and produces smoother results. The solution is more stable.

//! Also:
//! - It needs about 50% longer time/step
//! - What is more physical correct I dont know
//! - The original code used the conservative form

//!---------- Convective Deriivatives CONSERVATIVE FORM ---------
D_REAL d1_U1B2, d1_U2B1, d1_U1B3, d1_U3B1; 
D_REAL d2_U1B2, d2_U2B1, d2_U2B3, d2_U3B2; 
D_REAL d3_U1B3, d3_U3B1, d3_U2B3, d3_U3B2;

//!---------- Resistive Deriivatives CONSERVATIVE FORM ---------
D_REAL d2_Eta_d1_b2, d3_Eta_d1_b3, d1_Eta_d2_b1;
D_REAL d3_Eta_d2_b3, d1_Eta_d3_b1, d2_Eta_d3_b2;



//!---------- Resivetive Deriivatives CONSERVATIVE FORM ---------
D_REAL d1_Eta, d2_Eta, d3_Eta;

D_REAL d1_B2, d1_B3;
D_REAL d2_B1, d2_B3;
D_REAL d3_B1, d3_B2;

D_REAL d11_B2, d11_B3;
D_REAL d22_B1, d22_B3;
D_REAL d33_B1, d33_B2;

D_REAL d1d3_B1, d1d2_B1;
D_REAL d2d1_B2, d2d3_B2;
D_REAL d3d1_B3, d3d2_B3;


//! HALL Derivatives 
D_REAL B1P, B2P, B3P;
D_REAL d1_B1P, d2_B2P, d3_B3P;

D_REAL d1_B1P_d2_b1, d1_B1P_d3_b1;
D_REAL d1_B2P_d2_b3, d1_B2P_d3_b2;
D_REAL d1_B3P_d2_b3, d1_B3P_d3_b2;

D_REAL d2_B1P_d1_b3, d2_B1P_d3_b1;
D_REAL d2_B2P_d1_b2, d2_B2P_d3_b2;
D_REAL d2_B3P_d1_b3, d2_B3P_d3_b1;


D_REAL d3_B1P_d1_b2, d3_B1P_d2_b1;
D_REAL d3_B2P_d1_b2, d3_B2P_d2_b1;
D_REAL d3_B3P_d1_b3, d3_B3P_d2_b3;


//! PE Derivatives
#if defined nonadiabatic_gradPE_TERM
//! NOTE: Only implemented for conservative derivatives
//!       and UNIFORM grid
//!       and LF (not RK!)

D_REAL* PE;
D_REAL* neutral_n;
D_REAL* neutral_u1;
D_REAL* neutral_u2;
D_REAL* neutral_u3;
D_REAL* neutral_p;
D_REAL* neutral_beta;

//! P means 1/rho, PE means electron pressure

#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only

#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
D_REAL d2_P_d3_PE, d3_P_d2_PE;
D_REAL d3_P_d1_PE, d1_P_d3_PE;
D_REAL d1_P_d2_PE, d2_P_d1_PE;
#endif

#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
D_REAL d2_PE_d3_P, d3_PE_d2_P;
D_REAL d3_PE_d1_P, d1_PE_d3_P;
D_REAL d1_PE_d2_P, d2_PE_d1_P;
#endif

#endif /* #ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only */

//! non-conservative form
#if (nonadiabatic_gradPE_TERM_derivation_form == 0) || (nonadiabatic_gradPE_TERM_curl_derivation_form == 2)
D_REAL d1_PE, d2_PE, d3_PE;
D_REAL d1_P, d2_P, d3_P;
#endif

//! conservative form
#if (nonadiabatic_gradPE_TERM_derivation_form == 1) || (nonadiabatic_gradPE_TERM_derivation_form == 2)

D_REAL d1_PEkappa_U1, d2_PEkappa_U2, d3_PEkappa_U3;
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
D_REAL d1_PEkappa_P, d2_PEkappa_P, d3_PEkappa_P;
#endif

#if nonadiabatic_gradPE_TERM_derivation_form == 2
D_REAL* PEkappa;
#endif

#endif




#endif



using namespace std;

//!-------------------------------------------------------------//
//! set_pointers_to_calc_derivatives: //
//!-------------------------------------------------------------//
#if defined(nonadiabatic_gradPE_TERM)
inline void CBlock::set_pointers_to_calc_derivatives(INT32 B_FIELD_ID,
						       INT32 U_FIELD_ID,
					              INT32 rezRHO_FIELD_ID,
					              INT32 ETA_FIELD_ID,
						      INT32 PE_FIELD_ID	)
{
	//! first rd of respective level
	 rd = rd_of_L[RLevel];
	rd2 = rd2_of_L[RLevel];

	B1 = Field_Type[B_FIELD_ID];
	B2 = B1 +num_nodes_in_block;
	B3 = B2 +num_nodes_in_block;

	U1 = Field_Type[U_FIELD_ID];
	U2 = U1 +num_nodes_in_block;
	U3 = U2 +num_nodes_in_block;

	ETA = Field_Type[ETA_FIELD_ID];
	rRHO = Field_Type[rezRHO_FIELD_ID];
	
	PE = Field_Type[PE_FIELD_ID];
#if nonadiabatic_gradPE_TERM_derivation_form == 2
	PEkappa = Field_Type[id_scratch_scalar];
#endif
	
}
#endif


inline void CBlock::set_pointers_to_calc_derivatives(INT32 B_FIELD_ID,
					       	      INT32 U_FIELD_ID,
					              INT32 rezRHO_FIELD_ID,
					              INT32 ETA_FIELD_ID)
{
	//! first rd of respective level
	 rd = rd_of_L[RLevel];
	rd2 = rd2_of_L[RLevel];

	B1 = Field_Type[B_FIELD_ID];
	B2 = B1 +num_nodes_in_block;
	B3 = B2 +num_nodes_in_block;

	U1 = Field_Type[U_FIELD_ID];
	U2 = U1 +num_nodes_in_block;
	U3 = U2 +num_nodes_in_block;

	ETA = Field_Type[ETA_FIELD_ID];
	rRHO = Field_Type[rezRHO_FIELD_ID];
	
#if defined(nonadiabatic_gradPE_TERM)
	PE = Field_Type[id_PEtotal];
#endif
}

//!-------------------------------------------------------------//
//! set_pointers_to_calc_derivatives: //
//!-------------------------------------------------------------//
#if defined(nonadiabatic_gradPE_TERM)
inline void CBlock::set_pointers_to_calc_derivatives_CFBG_BField(INT32 B_SW_ID,
								   INT32 B_TOTAL_ID,
								   INT32 PE_FIELD_ID,
								   INT32 U_FIELD_ID,
								   INT32 rezRHO_FIELD_ID,
								   INT32 ETA_FIELD_ID)
{
	//! first rd of respective level
	 rd = rd_of_L[RLevel];
	rd2 = rd2_of_L[RLevel];

	B1_sw = Field_Type[B_SW_ID];
	B2_sw = B1_sw +num_nodes_in_block;
	B3_sw = B2_sw +num_nodes_in_block;

	B1_total = Field_Type[B_TOTAL_ID];
	B2_total = B1_total +num_nodes_in_block;
	B3_total = B2_total +num_nodes_in_block;

	//! Only first comp of id_Bcfbg needed
	B1_cfbg = Field_Type[id_Bcfbg];

	U1 = Field_Type[U_FIELD_ID];
	U2 = U1 +num_nodes_in_block;
	U3 = U2 +num_nodes_in_block;

	ETA = Field_Type[ETA_FIELD_ID];
	rRHO = Field_Type[rezRHO_FIELD_ID];
	
	PE = Field_Type[PE_FIELD_ID];
#if nonadiabatic_gradPE_TERM_derivation_form == 2
	PEkappa = Field_Type[id_scratch_scalar];
#endif
	
}
#endif

inline void CBlock::set_pointers_to_calc_derivatives_CFBG_BField(INT32 B_SW_ID,
								   INT32 B_TOTAL_ID,
								   INT32 U_FIELD_ID,
								   INT32 rezRHO_FIELD_ID,
								   INT32 ETA_FIELD_ID)
{
	//! first rd of respective level
	 rd = rd_of_L[RLevel];
	rd2 = rd2_of_L[RLevel];

	B1_sw = Field_Type[B_SW_ID];
	B2_sw = B1_sw +num_nodes_in_block;
	B3_sw = B2_sw +num_nodes_in_block;

	B1_total = Field_Type[B_TOTAL_ID];
	B2_total = B1_total +num_nodes_in_block;
	B3_total = B2_total +num_nodes_in_block;

	//! Only first comp of id_Bcfbg needed
	B1_cfbg = Field_Type[id_Bcfbg];

	U1 = Field_Type[U_FIELD_ID];
	U2 = U1 +num_nodes_in_block;
	U3 = U2 +num_nodes_in_block;

	ETA = Field_Type[ETA_FIELD_ID];
	rRHO = Field_Type[rezRHO_FIELD_ID];
	
#if defined(nonadiabatic_gradPE_TERM)
	PE = Field_Type[id_PEtotal];
#endif
}


//!------------------------------------------------------------------//
//! calc_derivatives:
//! NOTE: before calling this function, the function 
//! 	   "set_pointers_to_calc_derivatives" has to 
//! 	   be called, otherwise pointers B1,.. are undefined !!!
//!------------------------------------------------------------------//
inline void CBlock::calc_derivatives()
{

      i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

      //! --- indices for firs derivatives ----------------------
      ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
      im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;

      i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      //! --- indices for mixed derivatives ---------------------
      ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
      im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
      im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
      i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
      i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
      i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);


      //!----------- rho - derivatives ---------------------------
      rP    = rRHO[i_j_k];

      dxP = rd[0]*(rRHO[ip1_j_k]-rRHO[im1_j_k]);
      dyP = rd[1]*(rRHO[i_jp1_k]-rRHO[i_jm1_k]);
      dzP = rd[2]*(rRHO[i_j_kp1]-rRHO[i_j_km1]);

      //!----------- Eta - derivatives ---------------------------
      Eta    = ETA[i_j_k];

      dxEta = rd[0]*(ETA[ip1_j_k]-ETA[im1_j_k]);
      dyEta = rd[1]*(ETA[i_jp1_k]-ETA[i_jm1_k]);
      dzEta = rd[2]*(ETA[i_j_kp1]-ETA[i_j_km1]);

      //!----------- U - derivatives -----------------------------
      ux   = U1[i_j_k];
      uy   = U2[i_j_k];
      uz   = U3[i_j_k];

      dxUX = rd[0]*(U1[ip1_j_k]-U1[im1_j_k]);
      dxUY = rd[0]*(U2[ip1_j_k]-U2[im1_j_k]);
      dxUZ = rd[0]*(U3[ip1_j_k]-U3[im1_j_k]);


      dyUX = rd[1]*(U1[i_jp1_k]-U1[i_jm1_k]);
      dyUY = rd[1]*(U2[i_jp1_k]-U2[i_jm1_k]);
      dyUZ = rd[1]*(U3[i_jp1_k]-U3[i_jm1_k]);


      dzUX = rd[2]*(U1[i_j_kp1]-U1[i_j_km1]);
      dzUY = rd[2]*(U2[i_j_kp1]-U2[i_j_km1]);
      dzUZ = rd[2]*(U3[i_j_kp1]-U3[i_j_km1]);

      //!----------- B - derivatives ------------------------------
      bx   = B1[i_j_k];
      by   = B2[i_j_k];
      bz   = B3[i_j_k];

      //!------ first derivatives ---------------------------------
      dxBX = rd[0]*(B1[ip1_j_k]-B1[im1_j_k]);
      dxBY = rd[0]*(B2[ip1_j_k]-B2[im1_j_k]);
      dxBZ = rd[0]*(B3[ip1_j_k]-B3[im1_j_k]);

      dyBX = rd[1]*(B1[i_jp1_k]-B1[i_jm1_k]);
      dyBY = rd[1]*(B2[i_jp1_k]-B2[i_jm1_k]);
      dyBZ = rd[1]*(B3[i_jp1_k]-B3[i_jm1_k]);

      dzBX = rd[2]*(B1[i_j_kp1]-B1[i_j_km1]);
      dzBY = rd[2]*(B2[i_j_kp1]-B2[i_j_km1]);
      dzBZ = rd[2]*(B3[i_j_kp1]-B3[i_j_km1]);


      //!------ second derivatives --(just mixed occure) ----------
      d2xBY = rd2[0]*(B2[ip1_j_k]-2.*B2[i_j_k]+B2[im1_j_k]);
      d2xBZ = rd2[0]*(B3[ip1_j_k]-2.*B3[i_j_k]+B3[im1_j_k]);

      d2yBX = rd2[1]*(B1[i_jp1_k]-2.*B1[i_j_k]+B1[i_jm1_k]);
      d2yBZ = rd2[1]*(B3[i_jp1_k]-2.*B3[i_j_k]+B3[i_jm1_k]);

      d2zBX = rd2[2]*(B1[i_j_kp1]-2.*B1[i_j_k]+B1[i_j_km1]);
      d2zBY = rd2[2]*(B2[i_j_kp1]-2.*B2[i_j_k]+B2[i_j_km1]);



      //!------ mixed derivatives ---------------------------------
      dxdzBX = rd[0]*rd[2] *(B1[ip1_j_kp1]-B1[ip1_j_km1]
                            -B1[im1_j_kp1]+B1[im1_j_km1]);

      dxdyBX = rd[0]*rd[1] *(B1[ip1_jp1_k]-B1[ip1_jm1_k]
                            -B1[im1_jp1_k]+B1[im1_jm1_k]);

      dydzBX = rd[1]*rd[2] *(B1[i_jp1_kp1]-B1[i_jp1_km1]
                            -B1[i_jm1_kp1]+B1[i_jm1_km1]);


      dydxBY = rd[1]*rd[0] *(B2[ip1_jp1_k]-B2[ip1_jm1_k]
                            -B2[im1_jp1_k]+B2[im1_jm1_k]);

      dzdxBY = rd[2]*rd[0] *(B2[ip1_j_kp1]-B2[ip1_j_km1]
                            -B2[im1_j_kp1]+B2[im1_j_km1]);

      dydzBY = rd[1]*rd[2] *(B2[i_jp1_kp1]-B2[i_jp1_km1]
                            -B2[i_jm1_kp1]+B2[i_jm1_km1]);



      dzdxBZ = rd[2]*rd[0] *(B3[ip1_j_kp1]-B3[ip1_j_km1]
                            -B3[im1_j_kp1]+B3[im1_j_km1]);

      dydxBZ = rd[1]*rd[0] *(B3[ip1_jp1_k]-B3[ip1_jm1_k]
                            -B3[im1_jp1_k]+B3[im1_jm1_k]);


      dzdyBZ = rd[2]*rd[1] *(B3[i_jp1_kp1]-B3[i_jm1_kp1]
                            -B3[i_jp1_km1]+B3[i_jm1_km1]);

	
}


//!------------------------------------------------------------------//
//! calc_derivatives_conservative:
//! basically means not using product rule
//! (unless 2. derivatives occur)
//!
//! NOTE: before calling this function, the function 
//! 	   "set_pointers_to_calc_derivatives" has to 
//! 	    be called, otherwise pointers B1,.. are undefined !!!	
//!------------------------------------------------------------------//
inline void CBlock::calc_derivatives_conservative()
{


	i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

      //! --- indices for firs derivatives ----------------------
      ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
      im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;

      i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      //! --- indices for mixed derivatives ---------------------
      ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
      im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
      im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
      i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
      i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
      i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);

	
	//!----- Konvective Derivatives ----------------------------------------
	//!--------- d1 --------------------------------------------------------
	d1_U1B2 =  ( U1[ip1_j_k] * B2[ip1_j_k]
		    -U1[im1_j_k] * B2[im1_j_k])*rd[0];
	
	d1_U2B1 =  ( U2[ip1_j_k] * B1[ip1_j_k]
		    -U2[im1_j_k] * B1[im1_j_k])*rd[0];
	
	d1_U1B3 =  ( U1[ip1_j_k] * B3[ip1_j_k]
		    -U1[im1_j_k] * B3[im1_j_k])*rd[0];
	
	d1_U3B1 =  ( U3[ip1_j_k] * B1[ip1_j_k]
		    -U3[im1_j_k] * B1[im1_j_k])*rd[0];
	
	
	//!--------- d2 --------------------------------------------------------
	d2_U1B2 = ( U1[i_jp1_k] * B2[i_jp1_k]
		   -U1[i_jm1_k] * B2[i_jm1_k])*rd[1];
	
	d2_U2B1 = ( U2[i_jp1_k] * B1[i_jp1_k]
		   -U2[i_jm1_k] * B1[i_jm1_k])*rd[1];
	
	
	d2_U2B3 = ( U2[i_jp1_k] * B3[i_jp1_k]
		   -U2[i_jm1_k] * B3[i_jm1_k])*rd[1];
	
	d2_U3B2 = ( U3[i_jp1_k] * B2[i_jp1_k]
		   -U3[i_jm1_k] * B2[i_jm1_k])*rd[1];
	
	
	//!--------- d3 --------------------------------------------------------
	d3_U1B3 = ( U1[i_j_kp1] * B3[i_j_kp1]
		   -U1[i_j_km1] * B3[i_j_km1])*rd[2];
	
	d3_U3B1 = ( U3[i_j_kp1] * B1[i_j_kp1]
		   -U3[i_j_km1] * B1[i_j_km1])*rd[2];
	
	d3_U2B3 = ( U2[i_j_kp1] * B3[i_j_kp1]
		   -U2[i_j_km1] * B3[i_j_km1])*rd[2];
	
	
	d3_U3B2 = ( U3[i_j_kp1] * B2[i_j_kp1]
		   -U3[i_j_km1] * B2[i_j_km1])*rd[2];
	
	
	
	
	//!---------- Resivetive Deriivatives CONSERVATIVE FORM ---------
	
	
	//!------ first derivatives --(just mixed occure)-----------
	d1_B2 = rd[0]*(B2[ip1_j_k]-B2[im1_j_k]);
	d1_B3 = rd[0]*(B3[ip1_j_k]-B3[im1_j_k]);
	
	d2_B1 = rd[1]*(B1[i_jp1_k]-B1[i_jm1_k]);
	d2_B3 = rd[1]*(B3[i_jp1_k]-B3[i_jm1_k]);
	
	d3_B1 = rd[2]*(B1[i_j_kp1]-B1[i_j_km1]);
	d3_B2 = rd[2]*(B2[i_j_kp1]-B2[i_j_km1]);
	
	//!------ second derivatives --(just mixed occure)-----------
	d11_B2 = rd2[0]*(B2[ip1_j_k]-2.*B2[i_j_k]+B2[im1_j_k]);
	d11_B3 = rd2[0]*(B3[ip1_j_k]-2.*B3[i_j_k]+B3[im1_j_k]);
	
	d22_B1 = rd2[1]*(B1[i_jp1_k]-2.*B1[i_j_k]+B1[i_jm1_k]);
	d22_B3 = rd2[1]*(B3[i_jp1_k]-2.*B3[i_j_k]+B3[i_jm1_k]);
	
	d33_B1 = rd2[2]*(B1[i_j_kp1]-2.*B1[i_j_k]+B1[i_j_km1]);
	d33_B2 = rd2[2]*(B2[i_j_kp1]-2.*B2[i_j_k]+B2[i_j_km1]);

	
	//!------ Eta derivatives -----------------------------------
#ifdef ETA_TERM_BField


	Eta = ETA[i_j_k];

	d1_Eta = rd[0]*(ETA[ip1_j_k]-ETA[im1_j_k]);
	d2_Eta = rd[1]*(ETA[i_jp1_k]-ETA[i_jm1_k]);
	d3_Eta = rd[2]*(ETA[i_j_kp1]-ETA[i_j_km1]);

	d2_Eta_d1_b2 = ( ETA[i_jp1_k]*(B2[ip1_jp1_k] - B2[im1_jp1_k])
			-ETA[i_jm1_k]*(B2[ip1_jm1_k] - B2[im1_jm1_k]))*rd[1]*rd[0];

	d3_Eta_d1_b3 = ( ETA[i_j_kp1]*(B3[ip1_j_kp1] - B3[im1_j_kp1])
			-ETA[i_j_km1]*(B3[ip1_j_km1] - B3[im1_j_km1]))*rd[2]*rd[0];


	d1_Eta_d2_b1 = ( ETA[ip1_j_k]*(B1[ip1_jp1_k] - B1[ip1_jm1_k])
			-ETA[im1_j_k]*(B1[im1_jp1_k] - B1[im1_jm1_k]))*rd[0]*rd[1];

	d3_Eta_d2_b3 = ( ETA[i_j_kp1]*(B3[i_jp1_kp1] - B3[i_jm1_kp1])
			-ETA[i_j_km1]*(B3[i_jp1_km1] - B3[i_jm1_km1]))*rd[2]*rd[1];

	d1_Eta_d3_b1 = ( ETA[ip1_j_k]*(B1[ip1_j_kp1] - B1[ip1_j_km1])
			-ETA[im1_j_k]*(B1[im1_j_kp1] - B1[im1_j_km1]))*rd[0]*rd[2];

	d2_Eta_d3_b2 = ( ETA[i_jp1_k]*(B2[i_jp1_kp1] - B2[i_jp1_km1])
			-ETA[i_jm1_k]*(B2[i_jm1_kp1] - B2[i_jm1_km1]))*rd[1]*rd[2];
			
#endif



  	//!------------------------------------------------------
	//!--------------- Hall Derivatives ---------------------
  	//!------------------------------------------------------

	//!-----------------------------
	B1P = B1[i_j_k]*rRHO[i_j_k];
        B2P = B2[i_j_k]*rRHO[i_j_k];
        B3P = B3[i_j_k]*rRHO[i_j_k];	


	//!------------------------------------------------------
	d1_B1P = (B1[ip1_j_k]*rRHO[ip1_j_k] - B1[im1_j_k]*rRHO[im1_j_k])*rd[0];
	d2_B2P = (B2[i_jp1_k]*rRHO[i_jp1_k] - B2[i_jm1_k]*rRHO[i_jm1_k])*rd[1];
	d3_B3P = (B3[i_j_kp1]*rRHO[i_j_kp1] - B3[i_j_km1]*rRHO[i_j_km1])*rd[2];

	

	//!-------------------------------------------------------
	d1_B1P_d2_b1 = ( B1[ip1_j_k]*rRHO[ip1_j_k]*(B1[ip1_jp1_k] - B1[ip1_jm1_k])
			-B1[im1_j_k]*rRHO[im1_j_k]*(B1[im1_jp1_k] - B1[im1_jm1_k]))*rd[0]*rd[1];

	d1_B1P_d3_b1 = ( B1[ip1_j_k]*rRHO[ip1_j_k]*(B1[ip1_j_kp1] - B1[ip1_j_km1])
			-B1[im1_j_k]*rRHO[im1_j_k]*(B1[im1_j_kp1] - B1[im1_j_km1]))*rd[0]*rd[2];

	d1_B2P_d2_b3 = ( B2[ip1_j_k]*rRHO[ip1_j_k]*(B3[ip1_jp1_k] - B3[ip1_jm1_k])
			-B2[im1_j_k]*rRHO[im1_j_k]*(B3[im1_jp1_k] - B3[im1_jm1_k]))*rd[0]*rd[1];

	d1_B2P_d3_b2 = ( B2[ip1_j_k]*rRHO[ip1_j_k]*(B2[ip1_j_kp1] - B2[ip1_j_km1])
			-B2[im1_j_k]*rRHO[im1_j_k]*(B2[im1_j_kp1] - B2[im1_j_km1]))*rd[0]*rd[2];

	d1_B3P_d2_b3 = ( B3[ip1_j_k]*rRHO[ip1_j_k]*(B3[ip1_jp1_k] - B3[ip1_jm1_k])
			-B3[im1_j_k]*rRHO[im1_j_k]*(B3[im1_jp1_k] - B3[im1_jm1_k]))*rd[0]*rd[1]; 

	d1_B3P_d3_b2 = ( B3[ip1_j_k]*rRHO[ip1_j_k]*(B2[ip1_j_kp1] - B2[ip1_j_km1])
			-B3[im1_j_k]*rRHO[im1_j_k]*(B2[im1_j_kp1] - B2[im1_j_km1]))*rd[0]*rd[2];


	//!------------------------------------------------------------------------------
	d2_B1P_d1_b3 = ( B1[i_jp1_k]*rRHO[i_jp1_k]*(B3[ip1_jp1_k] - B3[im1_jp1_k])
			-B1[i_jm1_k]*rRHO[i_jm1_k]*(B3[ip1_jm1_k] - B3[im1_jm1_k]))*rd[1]*rd[0];

	d2_B1P_d3_b1 = ( B1[i_jp1_k]*rRHO[i_jp1_k]*(B1[i_jp1_kp1] - B1[i_jp1_km1])
			-B1[i_jm1_k]*rRHO[i_jm1_k]*(B1[i_jm1_kp1] - B1[i_jm1_km1]))*rd[1]*rd[2]; 

	d2_B2P_d1_b2 = ( B2[i_jp1_k]*rRHO[i_jp1_k]*(B2[ip1_jp1_k] - B2[im1_jp1_k])
			-B2[i_jm1_k]*rRHO[i_jm1_k]*(B2[ip1_jm1_k] - B2[im1_jm1_k]))*rd[1]*rd[0];

	d2_B2P_d3_b2 = ( B2[i_jp1_k]*rRHO[i_jp1_k]*(B2[i_jp1_kp1] - B2[i_jp1_km1])
			-B2[i_jm1_k]*rRHO[i_jm1_k]*(B2[i_jm1_kp1] - B2[i_jm1_km1]))*rd[1]*rd[2];

	d2_B3P_d1_b3 = ( B3[i_jp1_k]*rRHO[i_jp1_k]*(B3[ip1_jp1_k] - B3[im1_jp1_k])
			-B3[i_jm1_k]*rRHO[i_jm1_k]*(B3[ip1_jm1_k] - B3[im1_jm1_k]))*rd[1]*rd[0];

	d2_B3P_d3_b1 = ( B3[i_jp1_k]*rRHO[i_jp1_k]*(B1[i_jp1_kp1] - B1[i_jp1_km1])
			-B3[i_jm1_k]*rRHO[i_jm1_k]*(B1[i_jm1_kp1] - B1[i_jm1_km1]))*rd[1]*rd[2];

	//!----------------------------------------------------------------------------
	d3_B1P_d1_b2 = ( B1[i_j_kp1]*rRHO[i_j_kp1]*(B2[ip1_j_kp1] - B2[im1_j_kp1])
			-B1[i_j_km1]*rRHO[i_j_km1]*(B2[ip1_j_km1] - B2[im1_j_km1]))*rd[2]*rd[0];

	d3_B1P_d2_b1 = ( B1[i_j_kp1]*rRHO[i_j_kp1]*(B1[i_jp1_kp1] - B1[i_jm1_kp1])
			-B1[i_j_km1]*rRHO[i_j_km1]*(B1[i_jp1_km1] - B1[i_jm1_km1]))*rd[2]*rd[1];

	d3_B2P_d1_b2 = ( B2[i_j_kp1]*rRHO[i_j_kp1]*(B2[ip1_j_kp1] - B2[im1_j_kp1])
			-B2[i_j_km1]*rRHO[i_j_km1]*(B2[ip1_j_km1] - B2[im1_j_km1]))*rd[2]*rd[0];

	d3_B2P_d2_b1 = ( B2[i_j_kp1]*rRHO[i_j_kp1]*(B1[i_jp1_kp1] - B1[i_jm1_kp1])
			-B2[i_j_km1]*rRHO[i_j_km1]*(B1[i_jp1_km1] - B1[i_jm1_km1]))*rd[2]*rd[1];

	d3_B3P_d1_b3 = ( B3[i_j_kp1]*rRHO[i_j_kp1]*(B3[ip1_j_kp1] - B3[im1_j_kp1])
			-B3[i_j_km1]*rRHO[i_j_km1]*(B3[ip1_j_km1] - B3[im1_j_km1]))*rd[2]*rd[0];

	d3_B3P_d2_b3 = ( B3[i_j_kp1]*rRHO[i_j_kp1]*(B3[i_jp1_kp1] - B3[i_jm1_kp1])
			-B3[i_j_km1]*rRHO[i_j_km1]*(B3[i_jp1_km1] - B3[i_jm1_km1]))*rd[2]*rd[1];
	
	
	
	
	//!------------------------------------------------------
	//!---------------- PE Derivatives ----------------------
  	//!------------------------------------------------------
#if defined(nonadiabatic_gradPE_TERM)
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
	d2_P_d3_PE = (    rRHO[i_jp1_k] * ( PE[i_jp1_kp1] - PE[i_jp1_km1] )
			- rRHO[i_jm1_k] * ( PE[i_jm1_kp1] - PE[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_P_d2_PE = (    rRHO[i_j_kp1] * ( PE[i_jp1_kp1] - PE[i_jm1_kp1] )
			- rRHO[i_j_km1] * ( PE[i_jp1_km1] - PE[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_P_d1_PE = (    rRHO[i_j_kp1] * ( PE[ip1_j_kp1] - PE[im1_j_kp1] )
			- rRHO[i_j_km1] * ( PE[ip1_j_km1] - PE[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_P_d3_PE = (    rRHO[ip1_j_k] * ( PE[ip1_j_kp1] - PE[ip1_j_km1] )
			- rRHO[im1_j_k] * ( PE[im1_j_kp1] - PE[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_P_d2_PE = (    rRHO[ip1_j_k] * ( PE[ip1_jp1_k] - PE[ip1_jm1_k] )
			- rRHO[im1_j_k] * ( PE[im1_jp1_k] - PE[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_P_d1_PE = (    rRHO[i_jp1_k] * ( PE[ip1_jp1_k] - PE[im1_jp1_k] )
			- rRHO[i_jm1_k] * ( PE[ip1_jm1_k] - PE[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
	d2_PE_d3_P = (    PE[i_jp1_k] * ( rRHO[i_jp1_kp1] - rRHO[i_jp1_km1] )
			- PE[i_jm1_k] * ( rRHO[i_jm1_kp1] - rRHO[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_PE_d2_P = (    PE[i_j_kp1] * ( rRHO[i_jp1_kp1] - rRHO[i_jm1_kp1] )
			- PE[i_j_km1] * ( rRHO[i_jp1_km1] - rRHO[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_PE_d1_P = (    PE[i_j_kp1] * ( rRHO[ip1_j_kp1] - rRHO[im1_j_kp1] )
			- PE[i_j_km1] * ( rRHO[ip1_j_km1] - rRHO[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_PE_d3_P = (    PE[ip1_j_k] * ( rRHO[ip1_j_kp1] - rRHO[ip1_j_km1] )
			- PE[im1_j_k] * ( rRHO[im1_j_kp1] - rRHO[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_PE_d2_P = (    PE[ip1_j_k] * ( rRHO[ip1_jp1_k] - rRHO[ip1_jm1_k] )
			- PE[im1_j_k] * ( rRHO[im1_jp1_k] - rRHO[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_PE_d1_P = (    PE[i_jp1_k] * ( rRHO[ip1_jp1_k] - rRHO[im1_jp1_k] )
			- PE[i_jm1_k] * ( rRHO[ip1_jm1_k] - rRHO[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#endif /* #ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only */
	
#if (nonadiabatic_gradPE_TERM_derivation_form == 0) || (nonadiabatic_gradPE_TERM_curl_derivation_form == 2)
	d1_PE = rd[0] * ( PE[ip1_j_k] - PE[im1_j_k] );
	d2_PE = rd[1] * ( PE[i_jp1_k] - PE[i_jm1_k] );
	d3_PE = rd[2] * ( PE[i_j_kp1] - PE[i_j_km1] );
	
	d1_P = rd[0] * ( rRHO[ip1_j_k] - rRHO[im1_j_k] );
	d2_P = rd[1] * ( rRHO[i_jp1_k] - rRHO[i_jm1_k] );
	d3_P = rd[2] * ( rRHO[i_j_kp1] - rRHO[i_j_km1] );
#endif
	
#endif /* #if defined(nonadiabatic_gradPE_TERM) */
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	
#if (nonadiabatic_gradPE_TERM_derivation_form == 0)
	//! Also needed in electron pressure equation
	dxUX = rd[0] * ( U1[ip1_j_k] - U1[im1_j_k] );
	dyUY = rd[1] * ( U2[i_jp1_k] - U2[i_jm1_k] );
	dzUZ = rd[2] * ( U3[i_j_kp1] - U3[i_j_km1] );
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 1
	d1_PEkappa_U1 = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * U1[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * U2[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * U3[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * rRHO[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * rRHO[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * rRHO[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * rRHO[i_j_km1]
				);
#endif
	
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 2
	d1_PEkappa_U1 = rd[0] * (
				    PEkappa[ip1_j_k] * U1[ip1_j_k]
				  - PEkappa[im1_j_k] * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    PEkappa[i_jp1_k] * U2[i_jp1_k]
				  - PEkappa[i_jm1_k] * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    PEkappa[i_j_kp1] * U3[i_j_kp1]
				  - PEkappa[i_j_km1] * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    PEkappa[ip1_j_k] * rRHO[ip1_j_k]
				  - PEkappa[im1_j_k] * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    PEkappa[i_jp1_k] * rRHO[i_jp1_k]
				  - PEkappa[i_jm1_k] * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    PEkappa[i_j_kp1] * rRHO[i_j_kp1]
				  - PEkappa[i_j_km1] * rRHO[i_j_km1]
				);
#endif
	
#endif
	
	
	
#endif /* #if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B) */
}


#if defined(nonadiabatic_gradPE_TERM) && !(defined nonadiabatic_gradPE_TERM_advance_Pe_with_B)
//!------------------------------------------------------------------//
//! calc_derivatives_for_Pe_conservative:
//! basically means not using product rule
//! (unless 2. derivatives occur)
//!
//! NOTE: before calling this function, the function 
//! 	   "set_pointers_to_calc_derivatives" has to 
//! 	    be called, otherwise pointers B1,.. are undefined !!!	
//!------------------------------------------------------------------//
inline void CBlock::calc_derivatives_for_Pe_conservative()
{


	i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

      //! --- indices for firs derivatives ----------------------
      ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
      im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;

      i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      //! --- indices for mixed derivatives ---------------------
      ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
      im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
      im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
      i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
      i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
      i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);


	//!---------- Resivetive Deriivatives CONSERVATIVE FORM ---------
	
	
	//!------ first derivatives --(just mixed occure)-----------
	d1_B2 = rd[0]*(B2[ip1_j_k]-B2[im1_j_k]);
	d1_B3 = rd[0]*(B3[ip1_j_k]-B3[im1_j_k]);
	
	d2_B1 = rd[1]*(B1[i_jp1_k]-B1[i_jm1_k]);
	d2_B3 = rd[1]*(B3[i_jp1_k]-B3[i_jm1_k]);
	
	d3_B1 = rd[2]*(B1[i_j_kp1]-B1[i_j_km1]);
	d3_B2 = rd[2]*(B2[i_j_kp1]-B2[i_j_km1]);
	
	
	//!------------------------------------------------------
	//!---------------- PE Derivatives ----------------------
  	//!------------------------------------------------------
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
	d2_P_d3_PE = (    rRHO[i_jp1_k] * ( PE[i_jp1_kp1] - PE[i_jp1_km1] )
			- rRHO[i_jm1_k] * ( PE[i_jm1_kp1] - PE[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_P_d2_PE = (    rRHO[i_j_kp1] * ( PE[i_jp1_kp1] - PE[i_jm1_kp1] )
			- rRHO[i_j_km1] * ( PE[i_jp1_km1] - PE[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_P_d1_PE = (    rRHO[i_j_kp1] * ( PE[ip1_j_kp1] - PE[im1_j_kp1] )
			- rRHO[i_j_km1] * ( PE[ip1_j_km1] - PE[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_P_d3_PE = (    rRHO[ip1_j_k] * ( PE[ip1_j_kp1] - PE[ip1_j_km1] )
			- rRHO[im1_j_k] * ( PE[im1_j_kp1] - PE[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_P_d2_PE = (    rRHO[ip1_j_k] * ( PE[ip1_jp1_k] - PE[ip1_jm1_k] )
			- rRHO[im1_j_k] * ( PE[im1_jp1_k] - PE[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_P_d1_PE = (    rRHO[i_jp1_k] * ( PE[ip1_jp1_k] - PE[im1_jp1_k] )
			- rRHO[i_jm1_k] * ( PE[ip1_jm1_k] - PE[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
	d2_PE_d3_P = (    PE[i_jp1_k] * ( rRHO[i_jp1_kp1] - rRHO[i_jp1_km1] )
			- PE[i_jm1_k] * ( rRHO[i_jm1_kp1] - rRHO[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_PE_d2_P = (    PE[i_j_kp1] * ( rRHO[i_jp1_kp1] - rRHO[i_jm1_kp1] )
			- PE[i_j_km1] * ( rRHO[i_jp1_km1] - rRHO[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_PE_d1_P = (    PE[i_j_kp1] * ( rRHO[ip1_j_kp1] - rRHO[im1_j_kp1] )
			- PE[i_j_km1] * ( rRHO[ip1_j_km1] - rRHO[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_PE_d3_P = (    PE[ip1_j_k] * ( rRHO[ip1_j_kp1] - rRHO[ip1_j_km1] )
			- PE[im1_j_k] * ( rRHO[im1_j_kp1] - rRHO[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_PE_d2_P = (    PE[ip1_j_k] * ( rRHO[ip1_jp1_k] - rRHO[ip1_jm1_k] )
			- PE[im1_j_k] * ( rRHO[im1_jp1_k] - rRHO[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_PE_d1_P = (    PE[i_jp1_k] * ( rRHO[ip1_jp1_k] - rRHO[im1_jp1_k] )
			- PE[i_jm1_k] * ( rRHO[ip1_jm1_k] - rRHO[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#endif
	
#if (nonadiabatic_gradPE_TERM_derivation_form == 0) || (nonadiabatic_gradPE_TERM_curl_derivation_form == 2)
	d1_PE = rd[0] * ( PE[ip1_j_k] - PE[im1_j_k] );
	d2_PE = rd[1] * ( PE[i_jp1_k] - PE[i_jm1_k] );
	d3_PE = rd[2] * ( PE[i_j_kp1] - PE[i_j_km1] );
	
	d1_P = rd[0] * ( rRHO[ip1_j_k] - rRHO[im1_j_k] );
	d2_P = rd[1] * ( rRHO[i_jp1_k] - rRHO[i_jm1_k] );
	d3_P = rd[2] * ( rRHO[i_j_kp1] - rRHO[i_j_km1] );
#endif
#if (nonadiabatic_gradPE_TERM_derivation_form == 0)
	//! Also needed in electron pressure equation
	dxUX = rd[0] * ( U1[ip1_j_k] - U1[im1_j_k] );
	dyUY = rd[1] * ( U2[i_jp1_k] - U2[i_jm1_k] );
	dzUZ = rd[2] * ( U3[i_j_kp1] - U3[i_j_km1] );
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 1
	d1_PEkappa_U1 = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * U1[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * U2[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * U3[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * rRHO[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * rRHO[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * rRHO[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * rRHO[i_j_km1]
				);
#endif
	
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 2
	d1_PEkappa_U1 = rd[0] * (
				    PEkappa[ip1_j_k] * U1[ip1_j_k]
				  - PEkappa[im1_j_k] * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    PEkappa[i_jp1_k] * U2[i_jp1_k]
				  - PEkappa[i_jm1_k] * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    PEkappa[i_j_kp1] * U3[i_j_kp1]
				  - PEkappa[i_j_km1] * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    PEkappa[ip1_j_k] * rRHO[ip1_j_k]
				  - PEkappa[im1_j_k] * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    PEkappa[i_jp1_k] * rRHO[i_jp1_k]
				  - PEkappa[i_jm1_k] * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    PEkappa[i_j_kp1] * rRHO[i_j_kp1]
				  - PEkappa[i_j_km1] * rRHO[i_j_km1]
				);
#endif
	
#endif
	
}
#endif



//!------------------------------------------------------------------//
//! calc_derivatives_conservative_CFBG_BField:
//! - "conservative" basically means not using product rule
//!    (unless 2. derivatives occur)
//! - id_Bcfbg denotes a Background Field which analytically
//!   is Curl Free (e.g. Dipol Field). Of course the numerical curl
//!   does not vanish identically. In order to correct for this effect,
//!   it is not considered when curl is calculated. 
//!   The total magnetic field is:
//!   B = B_sw + B_cfbg
//!
//! NOTE: before calling this function, the function 
//! 	   "set_pointers_to_calc_derivatives_CFBG_BField" has to 
//! 	    be called, otherwise pointers B1_sw,.. are undefined !!!	
//!------------------------------------------------------------------//
inline void CBlock::calc_derivatives_conservative_CFBG_BField()
{


	i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

      //! --- indices for firs derivatives ----------------------
      ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
      im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;

      i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      //! --- indices for mixed derivatives ---------------------
      ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
      im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
      im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
      i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
      i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
      i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);



//!---------------------------------------------------------------------
//!----- Convective Derivatives ----------------------------------------
//!---------------------------------------------------------------------

	//!--------- d1 --------------------------------------------------------
	d1_U1B2 =  ( U1[ip1_j_k] * B2_total[ip1_j_k]
		    -U1[im1_j_k] * B2_total[im1_j_k])*rd[0];
	
	d1_U2B1 =  ( U2[ip1_j_k] * B1_total[ip1_j_k]
		    -U2[im1_j_k] * B1_total[im1_j_k])*rd[0];
	
	d1_U1B3 =  ( U1[ip1_j_k] * B3_total[ip1_j_k]
		    -U1[im1_j_k] * B3_total[im1_j_k])*rd[0];
	
	d1_U3B1 =  ( U3[ip1_j_k] * B1_total[ip1_j_k]
		    -U3[im1_j_k] * B1_total[im1_j_k])*rd[0];
	
	
	//!--------- d2 --------------------------------------------------------
	d2_U1B2 = ( U1[i_jp1_k] * B2_total[i_jp1_k]
		   -U1[i_jm1_k] * B2_total[i_jm1_k])*rd[1];
	
	d2_U2B1 = ( U2[i_jp1_k] * B1_total[i_jp1_k]
		   -U2[i_jm1_k] * B1_total[i_jm1_k])*rd[1];
	
	
	d2_U2B3 = ( U2[i_jp1_k] * B3_total[i_jp1_k]
		   -U2[i_jm1_k] * B3_total[i_jm1_k])*rd[1];
	
	d2_U3B2 = ( U3[i_jp1_k] * B2_total[i_jp1_k]
		   -U3[i_jm1_k] * B2_total[i_jm1_k])*rd[1];
	
	
	//!--------- d3 --------------------------------------------------------
	d3_U1B3 = ( U1[i_j_kp1] * B3_total[i_j_kp1]
		   -U1[i_j_km1] * B3_total[i_j_km1])*rd[2];
	
	d3_U3B1 = ( U3[i_j_kp1] * B1_total[i_j_kp1]
		   -U3[i_j_km1] * B1_total[i_j_km1])*rd[2];
	
	d3_U2B3 = ( U2[i_j_kp1] * B3_total[i_j_kp1]
		   -U2[i_j_km1] * B3_total[i_j_km1])*rd[2];
	
	
	d3_U3B2 = ( U3[i_j_kp1] * B2_total[i_j_kp1]
		   -U3[i_j_km1] * B2_total[i_j_km1])*rd[2];
	
	
//!---------------------------------------------------------------------
//!---- FIRST & SECOND DERIVATIVES FOR HALL & ETA TERM -----------------
//!-------  (in both cases B_sw has to be derived)  --------------------
//!---------------------------------------------------------------------

	//!------ first derivatives --(just mixed occure)-----------
	//! these are used for resistive term where B_SW has to be considered
	d1_B2 = rd[0]*(B2_sw[ip1_j_k]-B2_sw[im1_j_k]);
	d1_B3 = rd[0]*(B3_sw[ip1_j_k]-B3_sw[im1_j_k]);
	
	d2_B1 = rd[1]*(B1_sw[i_jp1_k]-B1_sw[i_jm1_k]);
	d2_B3 = rd[1]*(B3_sw[i_jp1_k]-B3_sw[i_jm1_k]);
	
	d3_B1 = rd[2]*(B1_sw[i_j_kp1]-B1_sw[i_j_km1]);
	d3_B2 = rd[2]*(B2_sw[i_j_kp1]-B2_sw[i_j_km1]);
	
	//!------ second derivatives --(just mixed occure)-----------
	d11_B2 = rd2[0]*(B2_sw[ip1_j_k]-2.*B2_sw[i_j_k]+B2_sw[im1_j_k]);
	d11_B3 = rd2[0]*(B3_sw[ip1_j_k]-2.*B3_sw[i_j_k]+B3_sw[im1_j_k]);
	
	d22_B1 = rd2[1]*(B1_sw[i_jp1_k]-2.*B1_sw[i_j_k]+B1_sw[i_jm1_k]);
	d22_B3 = rd2[1]*(B3_sw[i_jp1_k]-2.*B3_sw[i_j_k]+B3_sw[i_jm1_k]);
	
	d33_B1 = rd2[2]*(B1_sw[i_j_kp1]-2.*B1_sw[i_j_k]+B1_sw[i_j_km1]);
	d33_B2 = rd2[2]*(B2_sw[i_j_kp1]-2.*B2_sw[i_j_k]+B2_sw[i_j_km1]);


	//!---------- Resivetive Deriivatives CONSERVATIVE FORM ---------	
#ifdef ETA_TERM_BField

	Eta = ETA[i_j_k];

	d1_Eta = rd[0]*(ETA[ip1_j_k]-ETA[im1_j_k]);
	d2_Eta = rd[1]*(ETA[i_jp1_k]-ETA[i_jm1_k]);
	d3_Eta = rd[2]*(ETA[i_j_kp1]-ETA[i_j_km1]);


	//! in Eta term ALWAYS B_sw has to be used
	d2_Eta_d1_b2 = ( ETA[i_jp1_k]*(B2_sw[ip1_jp1_k] - B2_sw[im1_jp1_k])
			-ETA[i_jm1_k]*(B2_sw[ip1_jm1_k] - B2_sw[im1_jm1_k]))*rd[1]*rd[0];

	d3_Eta_d1_b3 = ( ETA[i_j_kp1]*(B3_sw[ip1_j_kp1] - B3_sw[im1_j_kp1])
			-ETA[i_j_km1]*(B3_sw[ip1_j_km1] - B3_sw[im1_j_km1]))*rd[2]*rd[0];

	d1_Eta_d2_b1 = ( ETA[ip1_j_k]*(B1_sw[ip1_jp1_k] - B1_sw[ip1_jm1_k])
			-ETA[im1_j_k]*(B1_sw[im1_jp1_k] - B1_sw[im1_jm1_k]))*rd[0]*rd[1];

	d3_Eta_d2_b3 = ( ETA[i_j_kp1]*(B3_sw[i_jp1_kp1] - B3_sw[i_jm1_kp1])
			-ETA[i_j_km1]*(B3_sw[i_jp1_km1] - B3_sw[i_jm1_km1]))*rd[2]*rd[1];

	d1_Eta_d3_b1 = ( ETA[ip1_j_k]*(B1_sw[ip1_j_kp1] - B1_sw[ip1_j_km1])
			-ETA[im1_j_k]*(B1_sw[im1_j_kp1] - B1_sw[im1_j_km1]))*rd[0]*rd[2];

	d2_Eta_d3_b2 = ( ETA[i_jp1_k]*(B2_sw[i_jp1_kp1] - B2_sw[i_jp1_km1])
			-ETA[i_jm1_k]*(B2_sw[i_jm1_kp1] - B2_sw[i_jm1_km1]))*rd[1]*rd[2];
#endif



  	//!------------------------------------------------------
	//!--------------- Hall Derivatives ---------------------
  	//!------------------------------------------------------


	//!--- USE B_TOTAL FIELD: -----------------------
	B1P = B1_total[i_j_k]*rRHO[i_j_k];
        B2P = B2_total[i_j_k]*rRHO[i_j_k];
        B3P = B3_total[i_j_k]*rRHO[i_j_k];	


	//!--- USE B_TOTAL FIELD: -----------------------
	d1_B1P = (B1_total[ip1_j_k]*rRHO[ip1_j_k] - B1_total[im1_j_k]*rRHO[im1_j_k])*rd[0];
	d2_B2P = (B2_total[i_jp1_k]*rRHO[i_jp1_k] - B2_total[i_jm1_k]*rRHO[i_jm1_k])*rd[1];
	d3_B3P = (B3_total[i_j_kp1]*rRHO[i_j_kp1] - B3_total[i_j_km1]*rRHO[i_j_km1])*rd[2];

	

	//!--- USE B_TOTAL for First Derivatives,
	//!--- USE B_SW    for Second Derivatives
	d1_B1P_d2_b1 = ( B1_total[ip1_j_k]*rRHO[ip1_j_k]*(B1_sw[ip1_jp1_k] - B1_sw[ip1_jm1_k])
			-B1_total[im1_j_k]*rRHO[im1_j_k]*(B1_sw[im1_jp1_k] - B1_sw[im1_jm1_k]))*rd[0]*rd[1];

	d1_B1P_d3_b1 = ( B1_total[ip1_j_k]*rRHO[ip1_j_k]*(B1_sw[ip1_j_kp1] - B1_sw[ip1_j_km1])
			-B1_total[im1_j_k]*rRHO[im1_j_k]*(B1_sw[im1_j_kp1] - B1_sw[im1_j_km1]))*rd[0]*rd[2];

	d1_B2P_d2_b3 = ( B2_total[ip1_j_k]*rRHO[ip1_j_k]*(B3_sw[ip1_jp1_k] - B3_sw[ip1_jm1_k])
			-B2_total[im1_j_k]*rRHO[im1_j_k]*(B3_sw[im1_jp1_k] - B3_sw[im1_jm1_k]))*rd[0]*rd[1];

	d1_B2P_d3_b2 = ( B2_total[ip1_j_k]*rRHO[ip1_j_k]*(B2_sw[ip1_j_kp1] - B2_sw[ip1_j_km1])
			-B2_total[im1_j_k]*rRHO[im1_j_k]*(B2_sw[im1_j_kp1] - B2_sw[im1_j_km1]))*rd[0]*rd[2];

	d1_B3P_d2_b3 = ( B3_total[ip1_j_k]*rRHO[ip1_j_k]*(B3_sw[ip1_jp1_k] - B3_sw[ip1_jm1_k])
			-B3_total[im1_j_k]*rRHO[im1_j_k]*(B3_sw[im1_jp1_k] - B3_sw[im1_jm1_k]))*rd[0]*rd[1]; 

	d1_B3P_d3_b2 = ( B3_total[ip1_j_k]*rRHO[ip1_j_k]*(B2_sw[ip1_j_kp1] - B2_sw[ip1_j_km1])
			-B3_total[im1_j_k]*rRHO[im1_j_k]*(B2_sw[im1_j_kp1] - B2_sw[im1_j_km1]))*rd[0]*rd[2];


	//!------------------------------------------------------------------------------
	d2_B1P_d1_b3 = ( B1_total[i_jp1_k]*rRHO[i_jp1_k]*(B3_sw[ip1_jp1_k] - B3_sw[im1_jp1_k])
			-B1_total[i_jm1_k]*rRHO[i_jm1_k]*(B3_sw[ip1_jm1_k] - B3_sw[im1_jm1_k]))*rd[1]*rd[0];

	d2_B1P_d3_b1 = ( B1_total[i_jp1_k]*rRHO[i_jp1_k]*(B1_sw[i_jp1_kp1] - B1_sw[i_jp1_km1])
			-B1_total[i_jm1_k]*rRHO[i_jm1_k]*(B1_sw[i_jm1_kp1] - B1_sw[i_jm1_km1]))*rd[1]*rd[2]; 

	d2_B2P_d1_b2 = ( B2_total[i_jp1_k]*rRHO[i_jp1_k]*(B2_sw[ip1_jp1_k] - B2_sw[im1_jp1_k])
			-B2_total[i_jm1_k]*rRHO[i_jm1_k]*(B2_sw[ip1_jm1_k] - B2_sw[im1_jm1_k]))*rd[1]*rd[0];

	d2_B2P_d3_b2 = ( B2_total[i_jp1_k]*rRHO[i_jp1_k]*(B2_sw[i_jp1_kp1] - B2_sw[i_jp1_km1])
			-B2_total[i_jm1_k]*rRHO[i_jm1_k]*(B2_sw[i_jm1_kp1] - B2_sw[i_jm1_km1]))*rd[1]*rd[2];

	d2_B3P_d1_b3 = ( B3_total[i_jp1_k]*rRHO[i_jp1_k]*(B3_sw[ip1_jp1_k] - B3_sw[im1_jp1_k])
			-B3_total[i_jm1_k]*rRHO[i_jm1_k]*(B3_sw[ip1_jm1_k] - B3_sw[im1_jm1_k]))*rd[1]*rd[0];

	d2_B3P_d3_b1 = ( B3_total[i_jp1_k]*rRHO[i_jp1_k]*(B1_sw[i_jp1_kp1] - B1_sw[i_jp1_km1])
			-B3_total[i_jm1_k]*rRHO[i_jm1_k]*(B1_sw[i_jm1_kp1] - B1_sw[i_jm1_km1]))*rd[1]*rd[2];

	//!----------------------------------------------------------------------------
	d3_B1P_d1_b2 = ( B1_total[i_j_kp1]*rRHO[i_j_kp1]*(B2_sw[ip1_j_kp1] - B2_sw[im1_j_kp1])
			-B1_total[i_j_km1]*rRHO[i_j_km1]*(B2_sw[ip1_j_km1] - B2_sw[im1_j_km1]))*rd[2]*rd[0];

	d3_B1P_d2_b1 = ( B1_total[i_j_kp1]*rRHO[i_j_kp1]*(B1_sw[i_jp1_kp1] - B1_sw[i_jm1_kp1])
			-B1_total[i_j_km1]*rRHO[i_j_km1]*(B1_sw[i_jp1_km1] - B1_sw[i_jm1_km1]))*rd[2]*rd[1];

	d3_B2P_d1_b2 = ( B2_total[i_j_kp1]*rRHO[i_j_kp1]*(B2_sw[ip1_j_kp1] - B2_sw[im1_j_kp1])
			-B2_total[i_j_km1]*rRHO[i_j_km1]*(B2_sw[ip1_j_km1] - B2_sw[im1_j_km1]))*rd[2]*rd[0];

	d3_B2P_d2_b1 = ( B2_total[i_j_kp1]*rRHO[i_j_kp1]*(B1_sw[i_jp1_kp1] - B1_sw[i_jm1_kp1])
			-B2_total[i_j_km1]*rRHO[i_j_km1]*(B1_sw[i_jp1_km1] - B1_sw[i_jm1_km1]))*rd[2]*rd[1];

	d3_B3P_d1_b3 = ( B3_total[i_j_kp1]*rRHO[i_j_kp1]*(B3_sw[ip1_j_kp1] - B3_sw[im1_j_kp1])
			-B3_total[i_j_km1]*rRHO[i_j_km1]*(B3_sw[ip1_j_km1] - B3_sw[im1_j_km1]))*rd[2]*rd[0];

	d3_B3P_d2_b3 = ( B3_total[i_j_kp1]*rRHO[i_j_kp1]*(B3_sw[i_jp1_kp1] - B3_sw[i_jm1_kp1])
			-B3_total[i_j_km1]*rRHO[i_j_km1]*(B3_sw[i_jp1_km1] - B3_sw[i_jm1_km1]))*rd[2]*rd[1];

	
	
	
	//!------------------------------------------------------
	//!---------------- PE Derivatives ----------------------
  	//!------------------------------------------------------
#if defined(nonadiabatic_gradPE_TERM)
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
	d2_P_d3_PE = (    rRHO[i_jp1_k] * ( PE[i_jp1_kp1] - PE[i_jp1_km1] )
			- rRHO[i_jm1_k] * ( PE[i_jm1_kp1] - PE[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_P_d2_PE = (    rRHO[i_j_kp1] * ( PE[i_jp1_kp1] - PE[i_jm1_kp1] )
			- rRHO[i_j_km1] * ( PE[i_jp1_km1] - PE[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_P_d1_PE = (    rRHO[i_j_kp1] * ( PE[ip1_j_kp1] - PE[im1_j_kp1] )
			- rRHO[i_j_km1] * ( PE[ip1_j_km1] - PE[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_P_d3_PE = (    rRHO[ip1_j_k] * ( PE[ip1_j_kp1] - PE[ip1_j_km1] )
			- rRHO[im1_j_k] * ( PE[im1_j_kp1] - PE[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_P_d2_PE = (    rRHO[ip1_j_k] * ( PE[ip1_jp1_k] - PE[ip1_jm1_k] )
			- rRHO[im1_j_k] * ( PE[im1_jp1_k] - PE[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_P_d1_PE = (    rRHO[i_jp1_k] * ( PE[ip1_jp1_k] - PE[im1_jp1_k] )
			- rRHO[i_jm1_k] * ( PE[ip1_jm1_k] - PE[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
	d2_PE_d3_P = (    PE[i_jp1_k] * ( rRHO[i_jp1_kp1] - rRHO[i_jp1_km1] )
			- PE[i_jm1_k] * ( rRHO[i_jm1_kp1] - rRHO[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_PE_d2_P = (    PE[i_j_kp1] * ( rRHO[i_jp1_kp1] - rRHO[i_jm1_kp1] )
			- PE[i_j_km1] * ( rRHO[i_jp1_km1] - rRHO[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_PE_d1_P = (    PE[i_j_kp1] * ( rRHO[ip1_j_kp1] - rRHO[im1_j_kp1] )
			- PE[i_j_km1] * ( rRHO[ip1_j_km1] - rRHO[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_PE_d3_P = (    PE[ip1_j_k] * ( rRHO[ip1_j_kp1] - rRHO[ip1_j_km1] )
			- PE[im1_j_k] * ( rRHO[im1_j_kp1] - rRHO[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_PE_d2_P = (    PE[ip1_j_k] * ( rRHO[ip1_jp1_k] - rRHO[ip1_jm1_k] )
			- PE[im1_j_k] * ( rRHO[im1_jp1_k] - rRHO[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_PE_d1_P = (    PE[i_jp1_k] * ( rRHO[ip1_jp1_k] - rRHO[im1_jp1_k] )
			- PE[i_jm1_k] * ( rRHO[ip1_jm1_k] - rRHO[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#endif /* #ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only */
	
#if (nonadiabatic_gradPE_TERM_derivation_form == 0) || (nonadiabatic_gradPE_TERM_curl_derivation_form == 2)
	d1_PE = rd[0] * ( PE[ip1_j_k] - PE[im1_j_k] );
	d2_PE = rd[1] * ( PE[i_jp1_k] - PE[i_jm1_k] );
	d3_PE = rd[2] * ( PE[i_j_kp1] - PE[i_j_km1] );
	
	d1_P = rd[0] * ( rRHO[ip1_j_k] - rRHO[im1_j_k] );
	d2_P = rd[1] * ( rRHO[i_jp1_k] - rRHO[i_jm1_k] );
	d3_P = rd[2] * ( rRHO[i_j_kp1] - rRHO[i_j_km1] );
#endif
	
#endif /* #if defined(nonadiabatic_gradPE_TERM) */
	
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	
#if (nonadiabatic_gradPE_TERM_derivation_form == 0)
	//! Also needed in electron pressure equation
	dxUX = rd[0] * ( U1[ip1_j_k] - U1[im1_j_k] );
	dyUY = rd[1] * ( U2[i_jp1_k] - U2[i_jm1_k] );
	dzUZ = rd[2] * ( U3[i_j_kp1] - U3[i_j_km1] );
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 1
	d1_PEkappa_U1 = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * U1[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * U2[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * U3[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * rRHO[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * rRHO[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * rRHO[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * rRHO[i_j_km1]
				);
#endif
	
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 2
	d1_PEkappa_U1 = rd[0] * (
				    PEkappa[ip1_j_k] * U1[ip1_j_k]
				  - PEkappa[im1_j_k] * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    PEkappa[i_jp1_k] * U2[i_jp1_k]
				  - PEkappa[i_jm1_k] * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    PEkappa[i_j_kp1] * U3[i_j_kp1]
				  - PEkappa[i_j_km1] * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    PEkappa[ip1_j_k] * rRHO[ip1_j_k]
				  - PEkappa[im1_j_k] * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    PEkappa[i_jp1_k] * rRHO[i_jp1_k]
				  - PEkappa[i_jm1_k] * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    PEkappa[i_j_kp1] * rRHO[i_j_kp1]
				  - PEkappa[i_j_km1] * rRHO[i_j_km1]
				);
#endif
	
#endif
	
	
#endif /* #if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B) */
}


#if defined(nonadiabatic_gradPE_TERM) && !(defined nonadiabatic_gradPE_TERM_advance_Pe_with_B)
//!------------------------------------------------------------------//
//! calc_derivatives_for_Pe_conservative_CFBG_BField:
//! - "conservative" basically means not using product rule
//!    (unless 2. derivatives occur)
//! - id_Bcfbg denotes a Background Field which analytically
//!   is Curl Free (e.g. Dipol Field). Of course the numerical curl
//!   does not vanish identically. In order to correct for this effect,
//!   it is not considered when curl is calculated. 
//!   The total magnetic field is:
//!   B = B_sw + B_cfbg
//!
//! NOTE: before calling this function, the function 
//! 	   "set_pointers_to_calc_derivatives_CFBG_BField" has to 
//! 	    be called, otherwise pointers B1_sw,.. are undefined !!!	
//!------------------------------------------------------------------//
inline void CBlock::calc_derivatives_for_Pe_conservative_CFBG_BField()
{


	i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

      //! --- indices for firs derivatives ----------------------
      ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
      im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;

      i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      //! --- indices for mixed derivatives ---------------------
      ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
      im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;

      ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
      im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);

      i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
      i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
      i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
      i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);


//!---------------------------------------------------------------------
//!---- FIRST & SECOND DERIVATIVES FOR HALL & ETA TERM -----------------
//!-------  (in both cases B_sw has to be derived)  --------------------
//!---------------------------------------------------------------------

	//!------ first derivatives --(just mixed occure)-----------
	//! these are used for resistive term where B_SW has to be considered
	d1_B2 = rd[0]*(B2_sw[ip1_j_k]-B2_sw[im1_j_k]);
	d1_B3 = rd[0]*(B3_sw[ip1_j_k]-B3_sw[im1_j_k]);
	
	d2_B1 = rd[1]*(B1_sw[i_jp1_k]-B1_sw[i_jm1_k]);
	d2_B3 = rd[1]*(B3_sw[i_jp1_k]-B3_sw[i_jm1_k]);
	
	d3_B1 = rd[2]*(B1_sw[i_j_kp1]-B1_sw[i_j_km1]);
	d3_B2 = rd[2]*(B2_sw[i_j_kp1]-B2_sw[i_j_km1]);
	
	
	//!------------------------------------------------------
	//!---------------- PE Derivatives ----------------------
  	//!------------------------------------------------------
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
	d2_P_d3_PE = (    rRHO[i_jp1_k] * ( PE[i_jp1_kp1] - PE[i_jp1_km1] )
			- rRHO[i_jm1_k] * ( PE[i_jm1_kp1] - PE[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_P_d2_PE = (    rRHO[i_j_kp1] * ( PE[i_jp1_kp1] - PE[i_jm1_kp1] )
			- rRHO[i_j_km1] * ( PE[i_jp1_km1] - PE[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_P_d1_PE = (    rRHO[i_j_kp1] * ( PE[ip1_j_kp1] - PE[im1_j_kp1] )
			- rRHO[i_j_km1] * ( PE[ip1_j_km1] - PE[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_P_d3_PE = (    rRHO[ip1_j_k] * ( PE[ip1_j_kp1] - PE[ip1_j_km1] )
			- rRHO[im1_j_k] * ( PE[im1_j_kp1] - PE[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_P_d2_PE = (    rRHO[ip1_j_k] * ( PE[ip1_jp1_k] - PE[ip1_jm1_k] )
			- rRHO[im1_j_k] * ( PE[im1_jp1_k] - PE[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_P_d1_PE = (    rRHO[i_jp1_k] * ( PE[ip1_jp1_k] - PE[im1_jp1_k] )
			- rRHO[i_jm1_k] * ( PE[ip1_jm1_k] - PE[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
	d2_PE_d3_P = (    PE[i_jp1_k] * ( rRHO[i_jp1_kp1] - rRHO[i_jp1_km1] )
			- PE[i_jm1_k] * ( rRHO[i_jm1_kp1] - rRHO[i_jm1_km1] ) )*rd[1]*rd[2];
	
	d3_PE_d2_P = (    PE[i_j_kp1] * ( rRHO[i_jp1_kp1] - rRHO[i_jm1_kp1] )
			- PE[i_j_km1] * ( rRHO[i_jp1_km1] - rRHO[i_jm1_km1] ) )*rd[2]*rd[1];
	
	d3_PE_d1_P = (    PE[i_j_kp1] * ( rRHO[ip1_j_kp1] - rRHO[im1_j_kp1] )
			- PE[i_j_km1] * ( rRHO[ip1_j_km1] - rRHO[im1_j_km1] ) )*rd[2]*rd[0];
	
	d1_PE_d3_P = (    PE[ip1_j_k] * ( rRHO[ip1_j_kp1] - rRHO[ip1_j_km1] )
			- PE[im1_j_k] * ( rRHO[im1_j_kp1] - rRHO[im1_j_km1] ) )*rd[0]*rd[2];
	
	d1_PE_d2_P = (    PE[ip1_j_k] * ( rRHO[ip1_jp1_k] - rRHO[ip1_jm1_k] )
			- PE[im1_j_k] * ( rRHO[im1_jp1_k] - rRHO[im1_jm1_k] ) )*rd[0]*rd[1];
	
	d2_PE_d1_P = (    PE[i_jp1_k] * ( rRHO[ip1_jp1_k] - rRHO[im1_jp1_k] )
			- PE[i_jm1_k] * ( rRHO[ip1_jm1_k] - rRHO[im1_jm1_k] ) )*rd[1]*rd[0];
#endif
	
#endif
	
#if (nonadiabatic_gradPE_TERM_derivation_form == 0) || (nonadiabatic_gradPE_TERM_curl_derivation_form == 2)
	d1_PE = rd[0] * ( PE[ip1_j_k] - PE[im1_j_k] );
	d2_PE = rd[1] * ( PE[i_jp1_k] - PE[i_jm1_k] );
	d3_PE = rd[2] * ( PE[i_j_kp1] - PE[i_j_km1] );
	
	d1_P = rd[0] * ( rRHO[ip1_j_k] - rRHO[im1_j_k] );
	d2_P = rd[1] * ( rRHO[i_jp1_k] - rRHO[i_jm1_k] );
	d3_P = rd[2] * ( rRHO[i_j_kp1] - rRHO[i_j_km1] );
#endif
#if (nonadiabatic_gradPE_TERM_derivation_form == 0)
	//! Also needed in electron pressure equation
	dxUX = rd[0] * ( U1[ip1_j_k] - U1[im1_j_k] );
	dyUY = rd[1] * ( U2[i_jp1_k] - U2[i_jm1_k] );
	dzUZ = rd[2] * ( U3[i_j_kp1] - U3[i_j_km1] );
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 1
	d1_PEkappa_U1 = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * U1[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * U2[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * U3[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    pow(PE[ip1_j_k],1./kappa_electron) * rRHO[ip1_j_k]
				  - pow(PE[im1_j_k],1./kappa_electron) * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    pow(PE[i_jp1_k],1./kappa_electron) * rRHO[i_jp1_k]
				  - pow(PE[i_jm1_k],1./kappa_electron) * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    pow(PE[i_j_kp1],1./kappa_electron) * rRHO[i_j_kp1]
				  - pow(PE[i_j_km1],1./kappa_electron) * rRHO[i_j_km1]
				);
#endif
	
#endif
	
#if nonadiabatic_gradPE_TERM_derivation_form == 2
	d1_PEkappa_U1 = rd[0] * (
				    PEkappa[ip1_j_k] * U1[ip1_j_k]
				  - PEkappa[im1_j_k] * U1[im1_j_k]
				);
	d2_PEkappa_U2 = rd[1] * (
				    PEkappa[i_jp1_k] * U2[i_jp1_k]
				  - PEkappa[i_jm1_k] * U2[i_jm1_k]
				);
	d3_PEkappa_U3 = rd[2] * (
				    PEkappa[i_j_kp1] * U3[i_j_kp1]
				  - PEkappa[i_j_km1] * U3[i_j_km1]
				);
	
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
	d1_PEkappa_P  = rd[0] * (
				    PEkappa[ip1_j_k] * rRHO[ip1_j_k]
				  - PEkappa[im1_j_k] * rRHO[im1_j_k]
				);
	d2_PEkappa_P  = rd[1] * (
				    PEkappa[i_jp1_k] * rRHO[i_jp1_k]
				  - PEkappa[i_jm1_k] * rRHO[i_jm1_k]
				);
	d3_PEkappa_P  = rd[2] * (
				    PEkappa[i_j_kp1] * rRHO[i_j_kp1]
				  - PEkappa[i_j_km1] * rRHO[i_j_km1]
				);
#endif
	
#endif
}
#endif



//!-------------------------------------------------------------//
//! advance_B: -							//
//!-------------------------------------------------------------//
#ifndef nonadiabatic_gradPE_TERM
void CBlock::LF_B_Step(INT32 B_to_update,
		       INT32 B_derivatives,
		       D_REAL dt_step)
{

      D_REAL* B1_new = Field_Type[B_to_update];
      D_REAL* B2_new = B1_new +num_nodes_in_block;
      D_REAL* B3_new = B2_new +num_nodes_in_block;



      set_pointers_to_calc_derivatives(B_derivatives,
				       id_UI_minus,
				       id_rho_rez,
				       id_Eta);

      
      
      //! indices have to be i,j,k (as in calc_derivatives !!!)
      for (i=1; i < BlkNds_X-1; i++)
       for (j=1; j < BlkNds_Y-1; j++)
        for (k=1; k < BlkNds_Z-1; k++)
	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{


	    //! pointers to calc derivatives
	    //! have to be set !!!
	    calc_derivatives();


	    //!--------------------------------------------------------------------
            //!----------- B1 -----------------------------------------------------
	    //!--------------------------------------------------------------------

            B1_new[i_j_k] +=  dt_step *
            (

#ifdef CONVEC_TERM
		//!------ V x (u x B) ----------------------------------------------
		ux*dyBY + by*dyUX - uy*dyBX - bx*dyUY - uz*dzBX - bx*dzUZ + ux*dzBZ + bz*dzUX

#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		+bx*dyP*(dzBX-dxBZ) + rP*dyBX*(dzBX-dxBZ) - bx*rP*dydxBZ

		-(by*dyP+rP*dyBY)*(dyBZ-dzBY) - rP*by*(d2yBZ-dydzBY)

		-(bz*dzP + rP*dzBZ)*(dyBZ-dzBY)-rP*bz*(dzdyBZ-d2zBY)

		+bx*(dzP)*(dxBY-dyBX) + rP*(dzBX)*(dxBY-dyBX) + rP*bx*dzdxBY
#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------
		+ Eta*
		(
			d2yBX +d2zBX 
		      - dydxBY -dzdxBZ
		)
		
		-(
			dyEta * (dxBY - dyBX)
			-dzEta * (dzBX - dxBZ)
		 )


#endif
	   );


	    //!--------------------------------------------------------------------
            //!----------- B2 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            B2_new[i_j_k] +=  dt_step *
            (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		uy*dzBZ + bz*dzUY - uz*dzBY - by*dzUZ - ux*dxBY - by*dxUX + uy*dxBX + bx*dxUY

#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		+by*dzP*(dxBY-dyBX) + rP*dzBY*(dxBY-dyBX) - by*rP*dydzBX

		-(bz*dzP+rP*dzBZ)*(dzBX-dxBZ) - rP*bz*(d2zBX-dzdxBZ)

		-(bx*dxP + rP*dxBX)*(dzBX-dxBZ)-rP*bx*(dxdzBX-d2xBZ)

		+by*(dxP)*(dyBZ-dzBY) + rP*(dxBY)*(dyBZ-dzBY) + rP*by*dydxBZ

#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------
		+Eta*
		(
	        	d2xBY + d2zBY 
			- dxdyBX - dzdyBZ
		)

		-(
			dzEta * (dyBZ - dzBY)
			-dxEta * (dxBY - dyBX)
		)



#endif

	    );

	    //!--------------------------------------------------------------------
            //!----------- B3 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            B3_new[i_j_k] +=  dt_step *
	    (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		uz*dxBX + bx*dxUZ - ux*dxBZ - bz*dxUX - uy*dyBZ - bz*dyUY + uz*dyBY + by*dyUZ
#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		+bz*dxP*(dyBZ-dzBY) + rP*dxBZ*(dyBZ-dzBY) - bz*rP*dzdxBY

		-(bx*dxP+rP*dxBX)*(dxBY-dyBX) - rP*bx*(d2xBY-dxdyBX)

		-(by*dyP + rP*dyBY)*(dxBY-dyBX)-rP*by*(dydxBY-d2yBX)

		+bz*(dyP)*(dzBX-dxBZ) + rP*(dyBZ)*(dzBX-dxBZ) + rP*bz*dydzBX
#endif
#ifdef ETA_TERM_BField

		//!--------- - V x (Eta V X B)--------------------------------------
		+Eta*
		(
			   d2xBZ +d2yBZ
			 - dxdzBX - dydzBY
		)

		-(
			dxEta * (dzBX - dxBZ)
			-dyEta * (dyBZ - dzBY)
		)
#endif

	    );
	}

}
#endif /* nonadiabatic_gradPE_TERM */

//!-------------------------------------------------------------//
//! LF_B_Step_conservative: -							//
//!-------------------------------------------------------------//
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
void CBlock::LF_B_Step_conservative(INT32 B_to_update,
				    INT32 PE_to_update,
				    INT32 B_derivatives,
				    INT32 PE_derivatives,
				    D_REAL dt_step)
#else
void CBlock::LF_B_Step_conservative(INT32 B_to_update,
				    INT32 B_derivatives,
				    D_REAL dt_step)
#endif
{

      D_REAL* B1_new = Field_Type[B_to_update];
      D_REAL* B2_new = B1_new +num_nodes_in_block;
      D_REAL* B3_new = B2_new +num_nodes_in_block;
      
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
      D_REAL* PE_new = Field_Type[PE_to_update];
      D_REAL* DRAIN_e = Field_Type[id_rho_np1_recombined];
//       const D_REAL drain_e = 0.;
      
#if defined(nonadiabatic_gradPE_TERM_ensure_PE_not_negative)
      D_REAL* PE_old = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	if(PE_old[node] < 0.)
	{
	    PE_old[node] = 0.;
	}
      }
#endif

      
#endif



#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
      set_pointers_to_calc_derivatives(B_derivatives,
				       id_UI_minus,
				       id_rho_rez,
				       id_Eta,
				       PE_derivatives);
      
#if nonadiabatic_gradPE_TERM_derivation_form == 2
      
#if defined(use_vectorclass) && defined(kappa_electron_is_2)
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      VEC2_D_REAL vec_PE_old1;
      for(INT32 node=0; node<num_nodes_in_block; node+=2)
      {
	vec_PE_old1.load(PE_old1+node);
	vec_PE_old1 = sqrt(vec_PE_old1);
	vec_PE_old1.store(PEkappa+node);
      }
#elif defined(kappa_electron_is_2)
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	PEkappa[node] = sqrt( PE_old1[node] );
      }
#else
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	PEkappa[node] = pow( PE_old1[node] , 1./kappa_electron );
      }
#endif

#endif



      
#else
      set_pointers_to_calc_derivatives(B_derivatives,
				       id_UI_minus,
				       id_rho_rez,
				       id_Eta);
#endif
      
      

      //! indices have to be i,j,k (as in calc_derivatives !!!)
      for (i=1; i < BlkNds_X-1; i++)
       for (j=1; j < BlkNds_Y-1; j++)
        for (k=1; k < BlkNds_Z-1; k++)
	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{


	    //! pointers to calc derivatives
	    //! have to be set !!!
	    calc_derivatives_conservative();

	    //!--------------------------------------------------------------------
            //!----------- B1_odd -------------------------------------------------
	    //!--------------------------------------------------------------------

            B1_new[i_j_k] +=  dt_step *
            (

#ifdef CONVEC_TERM
		//!------ V x (u x B) ----------------------------------------------
		+( d2_U1B2 - d2_U2B1 - d3_U3B1 + d3_U1B3)

#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------

		-(	
	
			+ B2P*d22_B3   - B3P*d33_B2
			+ d2_B2P*d2_B3 - d3_B3P*d3_B2
			+ d2_B1P_d1_b3 - d2_B1P_d3_b1
			+ d3_B1P_d2_b1 - d3_B1P_d1_b2
			+ d3_B3P_d2_b3 - d2_B2P_d3_b2
	
		)


#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------

		-(
			 d2_Eta_d1_b2 - (d2_Eta*d2_B1) - Eta*d22_B1

			+d3_Eta_d1_b3 - (d3_Eta*d3_B1) - Eta*d33_B1
		 )

#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
		
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
		//!------ 1/2 * V x (rP grad(PE) ) ---------------------------------
		+ 0.5 * ( d2_P_d3_PE - d3_P_d2_PE )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
		//!------ 1/2 * V x (PE grad(rP) ) ---------------------------------
		- 0.5 * ( d2_PE_d3_P - d3_PE_d2_P )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 2
		//!------ 1/2 * grad(PE) * grad(rP) --------------------------------
		+ 0.5 * ( d2_P * d3_PE - d3_P * d2_PE )
	#endif

#endif
	   );


	    //!--------------------------------------------------------------------
            //!----------- B2 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            B2_new[i_j_k] +=  dt_step *
            (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		+(d3_U2B3 - d3_U3B2 - d1_U1B2 + d1_U2B1)



#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------

		-(
	
			+ B3P*d33_B1   - B1P*d11_B3
			+ d3_B3P*d3_B1 - d1_B1P*d1_B3

			+ d3_B2P_d2_b1 - d3_B2P_d1_b2 
			+ d1_B2P_d3_b2 - d1_B2P_d2_b3 
			+ d1_B1P_d3_b1 - d3_B3P_d1_b3
	
		)




#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------

		-(
			d1_Eta_d2_b1 - (d1_Eta*d1_B2) - Eta*d11_B2
	
			+d3_Eta_d2_b3 - (d3_Eta*d3_B2) - Eta*d33_B2
		)




#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
		
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
		//!------ 1/2 * V x (rP grad(PE) ) ---------------------------------
		+ 0.5 * ( d3_P_d1_PE - d1_P_d3_PE )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
		//!------ 1/2 * V x (PE grad(rP) ) ---------------------------------
		- 0.5 * ( d3_PE_d1_P - d1_PE_d3_P )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 2
		//!------ 1/2 * grad(PE) * grad(rP) --------------------------------
		+ 0.5 * ( d3_P * d1_PE - d1_P * d3_PE )
	#endif

#endif

	    );

	    //!--------------------------------------------------------------------
            //!----------- B3 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            B3_new[i_j_k] +=  dt_step *
	    (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		+(d1_U3B1 - d1_U1B3 - d2_U2B3 + d2_U3B2)



#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		-(	
	
			+ B1P*d11_B2   - B2P*d22_B1
			+ d1_B1P*d1_B2 - d2_B2P*d2_B1

			+ d1_B3P_d3_b2 - d1_B3P_d2_b3
			+ d2_B3P_d1_b3 - d2_B3P_d3_b1
			+ d2_B2P_d1_b2 - d1_B1P_d2_b1
	
		)



#endif
#ifdef ETA_TERM_BField

		//!--------- - V x (Eta V X B)--------------------------------------
	
		-(
			d1_Eta_d3_b1 - (d1_Eta*d1_B3) - Eta*d11_B3
	
			+d2_Eta_d3_b2 - (d2_Eta*d2_B3) - Eta*d22_B3
		)

#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
		
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
		//!------ 1/2 * V x (rP grad(PE) ) ---------------------------------
		+ 0.5 * ( d1_P_d2_PE - d2_P_d1_PE )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
		//!------ 1/2 * V x (PE grad(rP) ) ---------------------------------
		- 0.5 * ( d1_PE_d2_P - d2_PE_d1_P )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 2
		//!------ 1/2 * grad(PE) * grad(rP) --------------------------------
		+ 0.5 * ( d1_P * d2_PE - d2_P * d1_PE )
	#endif

#endif

	    );
	    
	    
	    
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)	    
	     //!--------------------------------------------------------------------
            //!----------- PE -----------------------------------------------------
	    //!--------------------------------------------------------------------
            PE_new[i_j_k] +=  dt_step *
	    (

#if nonadiabatic_gradPE_TERM_derivation_form == 0
		//!------ ( -u_i + rot(B)/rho ) * grad(PE) -----------------------------------
		+ ( -U1[i_j_k] + rRHO[i_j_k] * (d2_B3 - d3_B2) ) * d1_PE
		+ ( -U2[i_j_k] + rRHO[i_j_k] * (d3_B1 - d1_B3) ) * d2_PE
		+ ( -U3[i_j_k] + rRHO[i_j_k] * (d1_B2 - d2_B1) ) * d3_PE
		
		
		#ifndef nonadiabatic_gradPE_TERM_set_divUe_to_zero_DEBUG
		//!------ -kappa* PE *( div(U_i) - rot(B)*grad(1/rho_i) ) --------------------
		- kappa_electron * PE[i_j_k] *
			(
				dxUX + dyUY + dzUZ
				- (d2_B3 - d3_B2) * d1_P
				- (d3_B1 - d1_B3) * d2_P
				- (d1_B2 - d2_B1) * d3_P
			)
		#endif
		
#endif
			
#if (nonadiabatic_gradPE_TERM_derivation_form == 1) || (nonadiabatic_gradPE_TERM_derivation_form == 2)
	    #if defined(kappa_electron_is_2)
		+ 2.0 * sqrt(PE[i_j_k]) *
	    #else
		+ kappa_electron * pow(PE[i_j_k],1.-1./kappa_electron) *
	    #endif
		(
		  //!------   - div( PE^(1/kappa) * U )    -----------------------------------
		  - d1_PEkappa_U1
		  - d2_PEkappa_U2
		  - d3_PEkappa_U3
		  
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
		  //!------   + rot(B) * grad( pe^(1/kappa) / rho_i )    ---------------------
		  + (d2_B3 - d3_B2) * d1_PEkappa_P
		  + (d3_B1 - d1_B3) * d2_PEkappa_P
		  + (d1_B2 - d2_B1) * d3_PEkappa_P
#endif
		  
		)
			
#endif
		
			
		//!------ -drain_e * PE / rho_i -----------------------------------
#if defined nonadiabatic_gradPE_DRAIN_TERM
		- DRAIN_e[i_j_k] * PE[i_j_k] * rRHO[i_j_k]
#endif

	    );
	    
	    for(INT32 species=0; species<num_Neutral_Species; ++species)
	    {
	      neutral_n = Field_Type[id_numberdensity_neutralSpecies1 + species];
	      neutral_u1 = Field_Type[id_velocity_neutralSpecies1 + species];
	      neutral_u2 = neutral_u1 + num_nodes_in_block;
	      neutral_u3 = neutral_u2 + num_nodes_in_block;
	      neutral_p = Field_Type[id_pressure_neutralSpecies1 + species];
	      neutral_beta = Field_Type[id_new_electron_beta_neutralSpecies1 + species];
	      
	      if(neutral_n[i_j_k] < 10.*PART_REAL_PRECISION)
	      {
		continue;
	      }
	      

#if defined nonadiabatic_gradPE_COLLISION_TERM
	      if( rRHO[i_j_k] >= 10.*PART_REAL_PRECISION )
	      {
		PE_new[i_j_k] +=  dt_step *
		  (
		    //!------ 2/f_e * n_e*n_n*en_coll_param * (f_n p_n/n_n - f_e p_e/n_e + m_n (u_n-u_e)^2) ---
		    + (kappa_electron-1.) * neutral_n[i_j_k]/rRHO[i_j_k]*en_coll_param[species]
		    *(
			Neutral_DOF[species]*neutral_p[i_j_k]/neutral_n[i_j_k]
			- 2./(kappa_electron-1.) * PE[i_j_k] * rRHO[i_j_k]
			+ 2. * Neutral_Masses[species]
			  *(
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			   )
		    )
		  );
	      }
 #endif
	      
	      PE_new[i_j_k] +=  dt_step *
		(
		  0.
#if defined nonadiabatic_gradPE_SOURCE_TERM
		    //!------ source * m_e * 2/f_e * (-rot(B)/rho_i+u_i)^2 --------
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.) * mass_electron_norm
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			)
		     
		    //!------ -source * 4*m_e/f_e * (-rot(B)/rho_i+u_i)*u_n ---
		    - norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)*mass_electron_norm*2.
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *neutral_u1[i_j_k]
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *neutral_u2[i_j_k]
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *neutral_u3[i_j_k]
			)
		      
// 		    //!------ source * 2/f_e ( m_n*u_n^2 + f_n*p_n/(2*n_n) ) ---
// 		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)
// 		      * (
// 			      Neutral_Masses[species]
// 			      * (
// 				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
// 				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
// 				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
// 				)
// 			      + Neutral_DOF[species] * neutral_p[i_j_k] / (2. * neutral_n[i_j_k])
// 			)
// 		    
		    //!------ source * ( 2/f_e m_n*u_n^2 + T_{e,n} ) ---
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k]
		      * (
			      (kappa_electron-1.) * Neutral_Masses[species]
			      * (
				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
				)
			      + neutral_beta[i_j_k]
			)
		    
#endif
		);
	    }
#endif /* nonadiabatic_gradPE_TERM */
	    
	    
	}
}



//!-------------------------------------------------------------//
//! LF_B_Step_conservative: -							//
//!-------------------------------------------------------------//
#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
void CBlock::LF_Pe_Step_conservative( INT32 PE_to_update,
				      INT32 PE_derivatives,
				      D_REAL dt_step)
{
      D_REAL* PE_new = Field_Type[PE_to_update];
      D_REAL* DRAIN_e = Field_Type[id_rho_np1_recombined];
//       const D_REAL drain_e = 0.;
      
      
#if defined(nonadiabatic_gradPE_TERM_ensure_PE_not_negative)
      D_REAL* PE_old = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	if(PE_old[node] < 0.)
	{
	    PE_old[node] = 0.;
	}
      }
#endif


        set_pointers_to_calc_derivatives(id_BEven,
				       id_UI_minus,
				       id_rho_rez,
				       id_Eta,
				       PE_derivatives);
      
#if nonadiabatic_gradPE_TERM_derivation_form == 2

#if defined(use_vectorclass) && defined(kappa_electron_is_2)
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      VEC2_D_REAL vec_PE_old1;
      for(INT32 node=0; node<num_nodes_in_block; node+=2)
      {
	vec_PE_old1.load(PE_old1+node);
	vec_PE_old1 = sqrt(vec_PE_old1);
	vec_PE_old1.store(PEkappa+node);
      }
#elif defined(kappa_electron_is_2)
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	PEkappa[node] = sqrt( PE_old1[node] );
      }
#else
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	PEkappa[node] = pow( PE_old1[node] , 1./kappa_electron );
      }
#endif

#endif
      

      //! indices have to be i,j,k (as in calc_derivatives !!!)
      for (i=1; i < BlkNds_X-1; i++)
       for (j=1; j < BlkNds_Y-1; j++)
        for (k=1; k < BlkNds_Z-1; k++)
	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{


	    //! pointers to calc derivatives
	    //! have to be set !!!
	    calc_derivatives_for_Pe_conservative();
	    
	    
	     //!--------------------------------------------------------------------
            //!----------- PE -----------------------------------------------------
	    //!--------------------------------------------------------------------
            PE_new[i_j_k] +=  dt_step *
	    (

#if nonadiabatic_gradPE_TERM_derivation_form == 0
		//!------ ( -u_i + rot(B)/rho ) * grad(PE) -----------------------------------
		+ ( -U1[i_j_k] + rRHO[i_j_k] * (d2_B3 - d3_B2) ) * d1_PE
		+ ( -U2[i_j_k] + rRHO[i_j_k] * (d3_B1 - d1_B3) ) * d2_PE
		+ ( -U3[i_j_k] + rRHO[i_j_k] * (d1_B2 - d2_B1) ) * d3_PE
		
		#ifndef nonadiabatic_gradPE_TERM_set_divUe_to_zero_DEBUG
		//!------ -kappa* PE *( div(U_i) - rot(B)*grad(1/rho_i) ) --------------------
		- kappa_electron * PE[i_j_k] *
			(
				dxUX + dyUY + dzUZ
				- (d2_B3 - d3_B2) * d1_P
				- (d3_B1 - d1_B3) * d2_P
				- (d1_B2 - d2_B1) * d3_P
			)
		#endif
		
#endif
			
#if (nonadiabatic_gradPE_TERM_derivation_form == 1) || (nonadiabatic_gradPE_TERM_derivation_form == 2)
	    #if defined(kappa_electron_is_2)
		+ 2.0 * sqrt(PE[i_j_k]) *
	    #else
		+ kappa_electron * pow(PE[i_j_k],1.-1./kappa_electron) *
	    #endif
		(
		  //!------   - div( PE^(1/kappa) * U )    -----------------------------------
		  - d1_PEkappa_U1
		  - d2_PEkappa_U2
		  - d3_PEkappa_U3
		  
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
		  //!------   + rot(B) * grad( pe^(1/kappa) / rho_i )    ---------------------
		  + (d2_B3 - d3_B2) * d1_PEkappa_P
		  + (d3_B1 - d1_B3) * d2_PEkappa_P
		  + (d1_B2 - d2_B1) * d3_PEkappa_P
#endif
		  
		)
			
#endif
		
			
		//!------ -drain_e * PE / rho_i -----------------------------------
#if defined nonadiabatic_gradPE_DRAIN_TERM
		- DRAIN_e[i_j_k] * PE[i_j_k] * rRHO[i_j_k]
#endif

	    );
	    
	    for(INT32 species=0; species<num_Neutral_Species; ++species)
	    {
	      neutral_n = Field_Type[id_numberdensity_neutralSpecies1 + species];
	      neutral_u1 = Field_Type[id_velocity_neutralSpecies1 + species];
	      neutral_u2 = neutral_u1 + num_nodes_in_block;
	      neutral_u3 = neutral_u2 + num_nodes_in_block;
	      neutral_p = Field_Type[id_pressure_neutralSpecies1 + species];
	      neutral_beta = Field_Type[id_new_electron_beta_neutralSpecies1 + species];
	      
	      if(neutral_n[i_j_k] < 10.*PART_REAL_PRECISION)
	      {
		continue;
	      }
	      

#if defined nonadiabatic_gradPE_COLLISION_TERM
	      if( rRHO[i_j_k] >= 10.*PART_REAL_PRECISION )
	      {
		PE_new[i_j_k] +=  dt_step *
		  (
		    //!------ 2/f_e * n_e*n_n*en_coll_param * (f_n p_n/n_n - f_e p_e/n_e + m_n (u_n-u_e)^2) ---
		    + (kappa_electron-1.) * neutral_n[i_j_k]/rRHO[i_j_k]*en_coll_param[species]
		    *(
			Neutral_DOF[species]*neutral_p[i_j_k]/neutral_n[i_j_k]
			- 2./(kappa_electron-1.) * PE[i_j_k] * rRHO[i_j_k]
			+ 2. * Neutral_Masses[species]
			  *(
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			   )
		    )
		  );
	      }
 #endif
	      
	      PE_new[i_j_k] +=  dt_step *
		(
		  0.
#if defined nonadiabatic_gradPE_SOURCE_TERM
		    //!------ source * m_e * 2/f_e * (-rot(B)/rho_i+u_i)^2 --------
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.) * mass_electron_norm
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			)
		     
		    //!------ -source * 4*m_e/f_e * (-rot(B)/rho_i+u_i)*u_n ---
		    - norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)*mass_electron_norm*2.
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *neutral_u1[i_j_k]
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *neutral_u2[i_j_k]
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *neutral_u3[i_j_k]
			)
		      
		    //!------ source * 2/f_e ( m_n*u_n^2 + f_n*p_n/(2*n_n) ) ---
// 		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)
// 		      * (
// 			      Neutral_Masses[species]
// 			      * (
// 				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
// 				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
// 				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
// 				)
// 			      + Neutral_DOF[species] * neutral_p[i_j_k] / (2. * neutral_n[i_j_k])
// 			)
// 		    
		    //!------ source * ( 2/f_e m_n*u_n^2 + T_{e,n} ) ---
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k]
		      * (
			      (kappa_electron-1.) * Neutral_Masses[species]
			      * (
				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
				)
			      + neutral_beta[i_j_k]
			)
		    
#endif
		);
	    }
	}
}
#endif



//!-------------------------------------------------------------//
//! dtB_Faraday_at_cycle: -							//
//!-------------------------------------------------------------//
void CBlock::dtB_Faraday_at_cycle(INT32 CYCLE)
{

	INT32 newB, oldB;

	D_REAL dt_B;
	D_REAL *E1, *E2, *E3;
	D_REAL *B1_new, *B2_new, *B3_new;
	D_REAL *B1_old, *B2_old, *B3_old;


	E1 = Field_Type[id_BDerivative];
	E2 = E1 +num_nodes_in_block;
	E3 = E2 +num_nodes_in_block;

	bool cp_allBounds[6] = {0,0,0,0,0,0};
	cp_Boundaries(id_BDerivative, cp_allBounds);

	//! first rd of respective level
	 rd = rd_of_L[RLevel];


	if(CYCLE==1)
	{
		newB = id_BOdd;
		oldB = id_BEven;

		dt_B = 1.*dt_field_of_L[RLevel];
	}


	if(CYCLE==NUM_SUB_CYCLE+1)
	{

		newB = id_BEven;
		oldB = id_BOdd;

		dt_B = 1.*dt_field_of_L[RLevel];

	}

	if(CYCLE>1 && CYCLE<=NUM_SUB_CYCLE)
	{

		dt_B = 2.*dt_field_of_L[RLevel];

		if(CYCLE%2==0)
		{
	
			//! same fields for old & new
			newB = id_BEven;
			oldB = id_BEven;


		}
		else
		{
			//! same fields for old & new
			newB = id_BOdd;
			oldB = id_BOdd;

		}

	}

	B1_new = Field_Type[newB];
	B2_new = B1_new +num_nodes_in_block;
	B3_new = B2_new +num_nodes_in_block;

	B1_old = Field_Type[oldB];
	B2_old = B1_old +num_nodes_in_block;
	B3_old = B2_old +num_nodes_in_block;


	if(mesh_type==UNIFORM)
	//! indices have to be i,j,k (as in calc_derivatives !!!)
	for (INT32 i=1; i < BlkNds_X-1; i++)
	 for (INT32 j=1; j < BlkNds_Y-1; j++)
	  for (INT32 k=1; k < BlkNds_Z-1; k++)
	   if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	   {


      		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);


		//! dtB = -rotE
		B1_new[i_j_k] = B1_old[i_j_k] - dt_B*(    rd[1]*(E3[i_jp1_k]-E3[i_jm1_k])
							- rd[2]*(E2[i_j_kp1]-E2[i_j_km1]));

		B2_new[i_j_k] = B2_old[i_j_k] - dt_B*(    rd[2]*(E1[i_j_kp1]-E1[i_j_km1]) 
							- rd[0]*(E3[ip1_j_k]-E3[im1_j_k]));

		B3_new[i_j_k] = B3_old[i_j_k] - dt_B*(    rd[0]*(E2[ip1_j_k]-E2[im1_j_k])
							- rd[1]*(E1[i_jp1_k]-E1[i_jm1_k]));


	}


	if(mesh_type==STAGGERED)
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



		//!  dyE3
		D_REAL E3_uH_vp1_wH = 0.25 * ( E3[  u_vp1_w]
					      +E3[up1_vp1_w]
					      +E3[  u_vp1_wp1]
					      +E3[up1_vp1_wp1]);

		D_REAL E3_uH_v_wH =   0.25 * ( E3[  u_v_w]
					      +E3[up1_v_w]
					      +E3[  u_v_wp1]
					      +E3[up1_v_wp1]);


		//!  dzE2
		D_REAL E2_uH_vH_wp1 = 0.25 * ( E2[    u_v_wp1]
					      +E2[  u_vp1_wp1]
					      +E2[  up1_v_wp1]
					      +E2[up1_vp1_wp1]);

		D_REAL E2_uH_vH_w   = 0.25 * ( E2[    u_v_w]
					      +E2[  u_vp1_w]
					      +E2[  up1_v_w]
					      +E2[up1_vp1_w]);


		//! dtB1 = -rotE1
		B1_new[i_j_k] = B1_old[i_j_k] - dt_B*(    2.*rd[1]*(E3_uH_vp1_wH - E3_uH_v_wH)
							- 2.*rd[2]*(E2_uH_vH_wp1 - E2_uH_vH_w));



		//!  dzE1
		D_REAL E1_uH_vH_wp1 = 0.25 * ( E1[    u_v_wp1]
					      +E1[  u_vp1_wp1]
					      +E1[  up1_v_wp1]
					      +E1[up1_vp1_wp1]);

		D_REAL E1_uH_vH_w   = 0.25 * ( E1[    u_v_w]
					      +E1[  u_vp1_w]
					      +E1[  up1_v_w]
					      +E1[up1_vp1_w]);


		//!  dxE3
		D_REAL E3_up1_vH_wH = 0.25 * ( E3[up1_v_w]
					      +E3[up1_vp1_w]
					      +E3[up1_v_wp1]
					      +E3[up1_vp1_wp1]);

		D_REAL E3_u_vH_wH =   0.25 * ( E3[u_v_w]
					      +E3[u_vp1_w]
					      +E3[u_v_wp1]
					      +E3[u_vp1_wp1]);


		//! dtB2 = -rotE2
		B2_new[i_j_k] = B2_old[i_j_k] - dt_B*(    2.*rd[2]*(E1_uH_vH_wp1 - E1_uH_vH_w)
							- 2.*rd[0]*(E3_up1_vH_wH - E3_u_vH_wH));


		//!  dxE2
		D_REAL E2_up1_vH_wH = 0.25 * ( E2[up1_v_w]
					      +E2[up1_vp1_w]
					      +E2[up1_v_wp1]
					      +E2[up1_vp1_wp1]);

		D_REAL E2_u_vH_wH =   0.25 * ( E2[u_v_w]
					      +E2[u_vp1_w]
					      +E2[u_v_wp1]
					      +E2[u_vp1_wp1]);

		//!  dyE1
		D_REAL E1_uH_vp1_wH = 0.25 * ( E1[  u_vp1_w]
					      +E1[up1_vp1_w]
					      +E1[  u_vp1_wp1]
					      +E1[up1_vp1_wp1]);

		D_REAL E1_uH_v_wH =   0.25 * ( E1[  u_v_w]
					      +E1[up1_v_w]
					      +E1[  u_v_wp1]
					      +E1[up1_v_wp1]);



		//! dtB2 = -rotE2
		B3_new[i_j_k] = B3_old[i_j_k] - dt_B*(    2.*rd[0]*(E2_up1_vH_wH - E2_u_vH_wH)
							- 2.*rd[1]*(E1_uH_vp1_wH - E1_uH_v_wH));


	}





}


//!-------------------------------------------------------------//
//! calc_Faraday_CFBG_EField: -							//
//!-------------------------------------------------------------//
void CBlock::calc_Faraday_EField(INT32 id_B_FIMesh)
{

	//! first rd of respective level
	 rd = rd_of_L[RLevel];


	//! use id_BDerivative field for temporary EField
	//! and calculate
	//! dtB = -rotE afterwards
	E1 = Field_Type[id_BDerivative];
	E2 = E1 +num_nodes_in_block;
	E3 = E2 +num_nodes_in_block;

	//! field that is used to calc rotB
	D_REAL* B1_FI = Field_Type[id_B_FIMesh];
	D_REAL* B2_FI = B1_FI +num_nodes_in_block;
	D_REAL* B3_FI = B2_FI +num_nodes_in_block;


	U1 = Field_Type[id_UI_minus];
	U2 = U1 +num_nodes_in_block;
	U3 = U2 +num_nodes_in_block;

	ETA = Field_Type[id_Eta];
	rRHO = Field_Type[id_rho_rez];
	

#if defined nonadiabatic_gradPE_TERM
	PE = Field_Type[id_PEtotal];
#endif


#ifdef USE_CFBG_BFIELD
	B1_total = Field_Type[id_BTotal];
	B2_total = B1_total +num_nodes_in_block;
	B3_total = B2_total +num_nodes_in_block;

	//! Only first comp of id_Bcfbg needed
	//! either use mesh at full interger or
	//! half interger nodes
	B1_cfbg = Field_Type[id_Bcfbg];


#if defined(use_vectorclass)
	VEC4_D_REAL vec_B1_FI;
	VEC4_D_REAL vec_B1_cfbg;
	for (i=0; i < 3*num_nodes_in_block; i+=4)
	{
	  vec_B1_FI.load(B1_FI+i);
	  vec_B1_cfbg.load(B1_cfbg+i);
	  vec_B1_FI += vec_B1_cfbg;
	  vec_B1_FI.store(B1_total+i);
	}
#else
	//! loop over 3 components and all elements
	for (i=0; i < 3*num_nodes_in_block; i++)
	B1_total[i] = B1_FI[i] +  B1_cfbg[i];
#endif


#else

	B1_total = Field_Type[id_B_FIMesh];
	B2_total = B1_total +num_nodes_in_block;
	B3_total = B2_total +num_nodes_in_block;
#endif
	
	//! derivatives are calculated from 
	//! 1) B_FI
	//! 2) B_total
	
	//! indices have to be i,j,k (as in calc_derivatives !!!)
	for (INT32 i=1; i < BlkNds_X-1; i++)
	 for (INT32 j=1; j < BlkNds_Y-1; j++)
	  for (INT32 k=1; k < BlkNds_Z-1; k++)
	  if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
		{


      		INT32 i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

		INT32 ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		INT32 im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		INT32 i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		INT32 i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		INT32 i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		INT32 i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);



		D_REAL rotB1 = (  rd[1]*(B3_FI[i_jp1_k]-B3_FI[i_jm1_k])
			        - rd[2]*(B2_FI[i_j_kp1]-B2_FI[i_j_km1]));

		D_REAL rotB2 = (  rd[2]*(B1_FI[i_j_kp1]-B1_FI[i_j_km1]) 
				- rd[0]*(B3_FI[ip1_j_k]-B3_FI[im1_j_k]));

		D_REAL rotB3 = (  rd[0]*(B2_FI[ip1_j_k]-B2_FI[im1_j_k])
		 		- rd[1]*(B1_FI[i_jp1_k]-B1_FI[i_jm1_k]));


		E1[i_j_k] =

#ifdef CONVEC_TERM
  	    //! ----- -(uxB) -----------------------
		- (   U2[i_j_k] * B3_total[i_j_k]
	             -U3[i_j_k] * B2_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rRHO[i_j_k] * (   rotB2 * B3_total[i_j_k]
						 -  rotB3 * B2_total[i_j_k])
#endif


#ifdef ETA_TERM_BField
	    //!------ eta * rotB--------------------
	       + ETA[i_j_k] * rotB1
#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
	  //! ----- -0.5*(rP * gradPE) -----------------------
		- 0.5 * rRHO[i_j_k] * ( PE[ip1_j_k] - PE[im1_j_k] ) * rd[0]
#endif
	;

		E2[i_j_k] =

#ifdef CONVEC_TERM
  	    //! ----- -(uxB) -----------------------
		- (   U3[i_j_k] * B1_total[i_j_k]
	             -U1[i_j_k] * B3_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rRHO[i_j_k] * (   rotB3 * B1_total[i_j_k]
						 -  rotB1 * B3_total[i_j_k])
#endif


#ifdef ETA_TERM_BField
	    //!------ eta * rotB--------------------
	       + ETA[i_j_k] * rotB2
#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
	  //! ----- -0.5*(rP * gradPE) -----------------------
		- 0.5 * rRHO[i_j_k] * ( PE[i_jp1_k] - PE[i_jm1_k] ) * rd[1]
#endif
	;


		E3[i_j_k] =

#ifdef CONVEC_TERM
  	    //! ----- -(uxB) -----------------------
		- (   U1[i_j_k] * B2_total[i_j_k]
	             -U2[i_j_k] * B1_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rRHO[i_j_k] * (   rotB1 * B2_total[i_j_k]
						 -  rotB2 * B1_total[i_j_k])
#endif


#ifdef ETA_TERM_BField
	    //!------ eta * rotB--------------------
	       + ETA[i_j_k] * rotB3
#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
	  //! ----- -0.5*(rP * gradPE) -----------------------
		- 0.5 * rRHO[i_j_k] * ( PE[i_j_kp1] - PE[i_j_km1] ) * rd[2]
#endif
	;

	

	}

}

//!-------------------------------------------------------------//
//! calc_Faraday_CFBG_EField: -							//
//!-------------------------------------------------------------//
void CBlock::calc_Faraday_EField_staggered(INT32 id_B_FIMesh, INT32 id_B_HIMesh)
{

	//! first rd of respective level
	 rd = rd_of_L[RLevel];

	//! use id_BDerivative field for temporary EField
	//! and calculate
	//! dtB = -rotE afterwards
	E1 = Field_Type[id_BDerivative];
	E2 = E1 +num_nodes_in_block;
	E3 = E2 +num_nodes_in_block;

	//! field that is used to calc rotB
	D_REAL* B1_FI = Field_Type[id_B_FIMesh];
	D_REAL* B2_FI = B1_FI +num_nodes_in_block;
	D_REAL* B3_FI = B2_FI +num_nodes_in_block;


	U1 = Field_Type[id_UI_minus];
	U2 = U1 +num_nodes_in_block;
	U3 = U2 +num_nodes_in_block;

	ETA = Field_Type[id_Eta];
	rRHO = Field_Type[id_rho_rez];


#ifdef USE_CFBG_BFIELD
	B1_total = Field_Type[id_BTotal];
	B2_total = B1_total +num_nodes_in_block;
	B3_total = B2_total +num_nodes_in_block;


	//! Only first comp of id_Bcfbg needed
	//! either use half interger nodes
	B1_cfbg = Field_Type[id_Bcfbg_HIMesh];
	D_REAL* B1_HIMesh = Field_Type[id_B_HIMesh];

	//! loop over 3 components and all elements
	for (i=0; i < 3*num_nodes_in_block; i++)
	B1_total[i] = B1_HIMesh[i] +  B1_cfbg[i];
#else
	B1_total = Field_Type[id_B_HIMesh];
	B2_total = B1_total +num_nodes_in_block;
	B3_total = B2_total +num_nodes_in_block;
#endif
	
	//! derivatives are calculated from 
	//! 1) B_FI
	//! 2) B_total
	
	//! indices have to be i,j,k (as in calc_derivatives !!!)
	for (i=1; i < BlkNds_X-1; i++)
	 for (j=1; j < BlkNds_Y-1; j++)
	  for (k=1; k < BlkNds_Z-1; k++)
	  if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
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


		//!  dyB3
		D_REAL B3_iH_jp1_kH = 0.25 * ( B3_FI[  i_jp1_k]
					      +B3_FI[ip1_jp1_k]
					      +B3_FI[  i_jp1_kp1]
					      +B3_FI[ip1_jp1_kp1]);

		D_REAL B3_iH_j_kH =   0.25 * ( B3_FI[  i_j_k]
					      +B3_FI[ip1_j_k]
					      +B3_FI[  i_j_kp1]
					      +B3_FI[ip1_j_kp1]);


		//!  dzB2
		D_REAL B2_iH_jH_kp1 = 0.25 * ( B2_FI[    i_j_kp1]
					      +B2_FI[  i_jp1_kp1]
					      +B2_FI[  ip1_j_kp1]
					      +B2_FI[ip1_jp1_kp1]);

		D_REAL B2_iH_jH_k   = 0.25 * ( B2_FI[    i_j_k]
					      +B2_FI[  i_jp1_k]
					      +B2_FI[  ip1_j_k]
					      +B2_FI[ip1_jp1_k]);


		D_REAL rotB1 = (  2.*rd[1]*(B3_iH_jp1_kH - B3_iH_j_kH)
			        - 2.*rd[2]*(B2_iH_jH_kp1 - B2_iH_jH_k));



		//!  dzB1
		D_REAL B1_iH_jH_kp1 = 0.25 * ( B1_FI[    i_j_kp1]
					      +B1_FI[  i_jp1_kp1]
					      +B1_FI[  ip1_j_kp1]
					      +B1_FI[ip1_jp1_kp1]);

		D_REAL B1_iH_jH_k   = 0.25 * ( B1_FI[    i_j_k]
					      +B1_FI[  i_jp1_k]
					      +B1_FI[  ip1_j_k]
					      +B1_FI[ip1_jp1_k]);


		//!  dxB3
		D_REAL B3_ip1_jH_kH = 0.25 * ( B3_FI[ip1_j_k]
					      +B3_FI[ip1_jp1_k]
					      +B3_FI[ip1_j_kp1]
					      +B3_FI[ip1_jp1_kp1]);

		D_REAL B3_i_jH_kH =   0.25 * ( B3_FI[i_j_k]
					      +B3_FI[i_jp1_k]
					      +B3_FI[i_j_kp1]
					      +B3_FI[i_jp1_kp1]);


		D_REAL rotB2 = (  2.*rd[2]*(B1_iH_jH_kp1 - B1_iH_jH_k) 
				- 2.*rd[0]*(B3_ip1_jH_kH - B3_i_jH_kH));


		//!  dxB2
		D_REAL B2_ip1_jH_kH = 0.25 * ( B2_FI[ip1_j_k]
					      +B2_FI[ip1_jp1_k]
					      +B2_FI[ip1_j_kp1]
					      +B2_FI[ip1_jp1_kp1]);

		D_REAL B2_i_jH_kH =   0.25 * ( B2_FI[i_j_k]
					      +B2_FI[i_jp1_k]
					      +B2_FI[i_j_kp1]
					      +B2_FI[i_jp1_kp1]);

		//!  dyB1
		D_REAL B1_iH_jp1_kH = 0.25 * ( B1_FI[  i_jp1_k]
					      +B1_FI[ip1_jp1_k]
					      +B1_FI[  i_jp1_kp1]
					      +B1_FI[ip1_jp1_kp1]);

		D_REAL B1_iH_j_kH =   0.25 * ( B1_FI[  i_j_k]
					      +B1_FI[ip1_j_k]
					      +B1_FI[  i_j_kp1]
					      +B1_FI[ip1_j_kp1]);


		D_REAL rotB3 = (  2.*rd[0]*(B2_ip1_jH_kH - B2_i_jH_kH)
		 		- 2.*rd[1]*(B1_iH_jp1_kH - B1_iH_j_kH));


		E1[i_j_k] =

#ifdef CONVEC_TERM
  	    //! ----- -(uxB) -----------------------
		- (   U2[i_j_k] * B3_total[i_j_k]
	             -U3[i_j_k] * B2_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rRHO[i_j_k] * (   rotB2 * B3_total[i_j_k]
						 -  rotB3 * B2_total[i_j_k])
#endif


#ifdef ETA_TERM_BField
	    //!------ eta * rotB--------------------
	       + ETA[i_j_k] * rotB1
#endif
	;

		E2[i_j_k] =

#ifdef CONVEC_TERM
  	    //! ----- -(uxB) -----------------------
		- (   U3[i_j_k] * B1_total[i_j_k]
	             -U1[i_j_k] * B3_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rRHO[i_j_k] * (   rotB3 * B1_total[i_j_k]
						 -  rotB1 * B3_total[i_j_k])
#endif


#ifdef ETA_TERM_BField
	    //!------ eta * rotB--------------------
	       + ETA[i_j_k] * rotB2
#endif
	;


		E3[i_j_k] =

#ifdef CONVEC_TERM
  	    //! ----- -(uxB) -----------------------
		- (   U1[i_j_k] * B2_total[i_j_k]
	             -U2[i_j_k] * B1_total[i_j_k])
#endif
#ifdef HALL_TERM
	    //! eliminate Hall term inside Obstacle (by !Flag)
	    //! ----- -1/p*(rotBxB) -----------------
		+ (!Flag[i_j_k]) *rRHO[i_j_k] * (   rotB1 * B2_total[i_j_k]
						 -  rotB2 * B1_total[i_j_k])
#endif


#ifdef ETA_TERM_BField
	    //!------ eta * rotB--------------------
	       + ETA[i_j_k] * rotB3
#endif
	;

	

	}

}




//!-------------------------------------------------------------//
//! LF_B_Step_conservative_CFBG_BField: -							//
//!-------------------------------------------------------------//
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
void CBlock::LF_B_Step_conservative_CFBG_BField(INT32 B_to_update,
						INT32 PE_to_update,
				    		INT32 B_derivatives,
						INT32 PE_derivatives,
				    		D_REAL dt_step)
#else
void CBlock::LF_B_Step_conservative_CFBG_BField(INT32 B_to_update,
				    		INT32 B_derivatives,
				    		D_REAL dt_step)
#endif
{



      D_REAL *B1_new = Field_Type[B_to_update];
      D_REAL *B2_new = B1_new +num_nodes_in_block;
      D_REAL *B3_new = B2_new +num_nodes_in_block;
      

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
        D_REAL* PE_new = Field_Type[PE_to_update];
        D_REAL* DRAIN_e = Field_Type[id_rho_np1_recombined];
//       const D_REAL drain_e = 0.;
      
#if defined(nonadiabatic_gradPE_TERM_ensure_PE_not_negative)
      D_REAL* PE_old = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	if(PE_old[node] < 0.)
	{
	    PE_old[node] = 0.;
	}
      }
#endif
      
#endif


      //! - B_derivatives is equal to B_sw

#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
	set_pointers_to_calc_derivatives_CFBG_BField(B_derivatives,
						   id_BTotal,
						   PE_derivatives,
						   id_UI_minus,
						   id_rho_rez,
						   id_Eta);
#if nonadiabatic_gradPE_TERM_derivation_form == 2
      
#if defined(use_vectorclass) && defined(kappa_electron_is_2)
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      VEC2_D_REAL vec_PE_old1;
      for(INT32 node=0; node<num_nodes_in_block; node+=2)
      {
	vec_PE_old1.load(PE_old1+node);
	vec_PE_old1 = sqrt(vec_PE_old1);
	vec_PE_old1.store(PEkappa+node);
      }
#elif defined(kappa_electron_is_2)
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	PEkappa[node] = sqrt( PE_old1[node] );
      }
#else
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	PEkappa[node] = pow( PE_old1[node] , 1./kappa_electron );
      }
#endif

#endif

#else
	set_pointers_to_calc_derivatives_CFBG_BField(B_derivatives,
						   id_BTotal,
						   id_UI_minus,
						   id_rho_rez,
						   id_Eta);
#endif

      
      //! loop over 3 components and all elements
      //! ("set_pointers... " has to be called before this loop !!!)
#if defined(use_vectorclass)
	VEC4_D_REAL vec_B1_sw;
	VEC4_D_REAL vec_B1_cfbg;
	for (i=0; i < 3*num_nodes_in_block; i+=4)
	{
	  vec_B1_sw.load(B1_sw+i);
	  vec_B1_cfbg.load(B1_cfbg+i);
	  vec_B1_sw += vec_B1_cfbg;
	  vec_B1_sw.store(B1_total+i);
	}
#else
      for (i=0; i < 3*num_nodes_in_block; i++)
      B1_total[i] = B1_sw[i] +  B1_cfbg[i];
#endif

      //! derivatives are calculated from 
      //! 1) B_sw
      //! 2) B_total

      //! indices have to be i,j,k (as in calc_derivatives !!!)
      for (i=1; i < BlkNds_X-1; i++)
       for (j=1; j < BlkNds_Y-1; j++)
        for (k=1; k < BlkNds_Z-1; k++)
	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{


	    //! pointers to calc derivatives
	    //! have to be set !!!
	    calc_derivatives_conservative_CFBG_BField();

	    //!--------------------------------------------------------------------
            //!----------- B1 -----------------------------------------------------
	    //!--------------------------------------------------------------------

            B1_new[i_j_k] +=  dt_step *
            (

#ifdef CONVEC_TERM
		//!------ V x (u x B) ----------------------------------------------
		+( d2_U1B2 - d2_U2B1 - d3_U3B1 + d3_U1B3)

#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------

		-(	
	
			+   B2P*d22_B3 - B3P*d33_B2
			+ d2_B2P*d2_B3 - d3_B3P*d3_B2

			+ d2_B1P_d1_b3 - d2_B1P_d3_b1
			+ d3_B1P_d2_b1 - d3_B1P_d1_b2
			+ d3_B3P_d2_b3 - d2_B2P_d3_b2
	
		)


#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------

		-(
			 d2_Eta_d1_b2 - (d2_Eta*d2_B1) - Eta*d22_B1

			+d3_Eta_d1_b3 - (d3_Eta*d3_B1) - Eta*d33_B1
		 )

#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
		
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
		//!------ 1/2 * V x (rP grad(PE) ) ---------------------------------
		+ 0.5 * ( d2_P_d3_PE - d3_P_d2_PE )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
		//!------ 1/2 * V x (PE grad(rP) ) ---------------------------------
		- 0.5 * ( d2_PE_d3_P - d3_PE_d2_P )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 2
		//!------ 1/2 * grad(PE) * grad(rP) --------------------------------
		+ 0.5 * ( d2_P * d3_PE - d3_P * d2_PE )
	#endif

#endif
	   );



	    //!--------------------------------------------------------------------
            //!----------- B2 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            B2_new[i_j_k] +=  dt_step *
            (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		+(d3_U2B3 - d3_U3B2 - d1_U1B2 + d1_U2B1)



#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------

		-(

			+ B3P*d33_B1   - B1P*d11_B3
			+ d3_B3P*d3_B1 - d1_B1P*d1_B3

			+ d3_B2P_d2_b1 - d3_B2P_d1_b2 
			+ d1_B2P_d3_b2 - d1_B2P_d2_b3 
			+ d1_B1P_d3_b1 - d3_B3P_d1_b3
	
		)

#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------

		-(
			 d1_Eta_d2_b1 - (d1_Eta*d1_B2) - Eta*d11_B2
	
			+d3_Eta_d2_b3 - (d3_Eta*d3_B2) - Eta*d33_B2
		)

#endif

#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
		
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
		//!------ 1/2 * V x (rP grad(PE) ) ---------------------------------
		+ 0.5 * ( d3_P_d1_PE - d1_P_d3_PE )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
		//!------ 1/2 * V x (PE grad(rP) ) ---------------------------------
		- 0.5 * ( d3_PE_d1_P - d1_PE_d3_P )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 2
		//!------ 1/2 * grad(PE) * grad(rP) --------------------------------
		+ 0.5 * ( d3_P * d1_PE - d1_P * d3_PE )
	#endif

#endif

	    );

	    //!--------------------------------------------------------------------
            //!----------- B3 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            B3_new[i_j_k] +=  dt_step *
	    (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		+(d1_U3B1 - d1_U1B3 - d2_U2B3 + d2_U3B2)



#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		-(	
	
			+ B1P*d11_B2   - B2P*d22_B1
			+ d1_B1P*d1_B2 - d2_B2P*d2_B1

			+ d1_B3P_d3_b2 - d1_B3P_d2_b3
			+ d2_B3P_d1_b3 - d2_B3P_d3_b1
			+ d2_B2P_d1_b2 - d1_B1P_d2_b1
	
		)



#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------
	
		-(
			 d1_Eta_d3_b1 - (d1_Eta*d1_B3) - Eta*d11_B3
	
			+d2_Eta_d3_b2 - (d2_Eta*d2_B3) - Eta*d22_B3
		)
#endif
		
#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_adiabatic_term_in_B_only)
		
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 0
		//!------ 1/2 * V x (rP grad(PE) ) ---------------------------------
		+ 0.5 * ( d1_P_d2_PE - d2_P_d1_PE )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 1
		//!------ 1/2 * V x (PE grad(rP) ) ---------------------------------
		- 0.5 * ( d1_PE_d2_P - d2_PE_d1_P )
	#endif
	
	#if nonadiabatic_gradPE_TERM_curl_derivation_form == 2
		//!------ 1/2 * grad(PE) * grad(rP) --------------------------------
		+ 0.5 * ( d1_P * d2_PE - d2_P * d1_PE )
	#endif

#endif

	    );
	    
	    
	    
#if defined(nonadiabatic_gradPE_TERM) && defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)   
	    //!--------------------------------------------------------------------
            //!----------- PE -----------------------------------------------------
	    //!--------------------------------------------------------------------
            PE_new[i_j_k] +=  dt_step *
	    (

#if nonadiabatic_gradPE_TERM_derivation_form == 0
		//!------ ( -u_i + rot(B)/rho ) * grad(PE) -----------------------------------
		+ ( -U1[i_j_k] + rRHO[i_j_k] * (d2_B3 - d3_B2) ) * d1_PE
		+ ( -U2[i_j_k] + rRHO[i_j_k] * (d3_B1 - d1_B3) ) * d2_PE
		+ ( -U3[i_j_k] + rRHO[i_j_k] * (d1_B2 - d2_B1) ) * d3_PE
		
		#ifndef nonadiabatic_gradPE_TERM_set_divUe_to_zero_DEBUG
		//!------ -kappa* PE *( div(U_i) - rot(B)*grad(1/rho_i) ) --------------------
		- kappa_electron * PE[i_j_k] *
			(
				dxUX + dyUY + dzUZ
				- (d2_B3 - d3_B2) * d1_P
				- (d3_B1 - d1_B3) * d2_P
				- (d1_B2 - d2_B1) * d3_P
			)
		#endif
		
#endif
			
#if (nonadiabatic_gradPE_TERM_derivation_form == 1) || (nonadiabatic_gradPE_TERM_derivation_form == 2)
	    #if defined(kappa_electron_is_2)
		+ 2.0 * sqrt(PE[i_j_k]) *
	    #else
		+ kappa_electron * pow(PE[i_j_k],1.-1./kappa_electron) *
	    #endif
		(
		  //!------   - div( PE^(1/kappa) * U )    -----------------------------------
		  - d1_PEkappa_U1
		  - d2_PEkappa_U2
		  - d3_PEkappa_U3
		  
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
		  //!------   + rot(B) * grad( pe^(1/kappa) / rho_i )    ---------------------
		  + (d2_B3 - d3_B2) * d1_PEkappa_P
		  + (d3_B1 - d1_B3) * d2_PEkappa_P
		  + (d1_B2 - d2_B1) * d3_PEkappa_P
#endif
		  
		)
			
#endif
		
			
		//!------ -drain_e * PE / rho_i -----------------------------------
#if defined nonadiabatic_gradPE_DRAIN_TERM
		- DRAIN_e[i_j_k] * PE[i_j_k] * rRHO[i_j_k]
#endif

	    );
	    
	    for(INT32 species=0; species<num_Neutral_Species; ++species)
	    {
	      neutral_n = Field_Type[id_numberdensity_neutralSpecies1 + species];
	      neutral_u1 = Field_Type[id_velocity_neutralSpecies1 + species];
	      neutral_u2 = neutral_u1 + num_nodes_in_block;
	      neutral_u3 = neutral_u2 + num_nodes_in_block;
	      neutral_p = Field_Type[id_pressure_neutralSpecies1 + species];
	      neutral_beta = Field_Type[id_new_electron_beta_neutralSpecies1 + species];
	      
	      if(neutral_n[i_j_k] < 10.*PART_REAL_PRECISION)
	      {
		continue;
	      }
	      

#if defined nonadiabatic_gradPE_COLLISION_TERM
	      if( rRHO[i_j_k] >= 10.*PART_REAL_PRECISION )
	      {
		PE_new[i_j_k] +=  dt_step *
		  (
		    //!------ 2/f_e * n_e*n_n*en_coll_param * (f_n p_n/n_n - f_e p_e/n_e + m_n (u_n-u_e)^2) ---
		    + (kappa_electron-1.) * neutral_n[i_j_k]/rRHO[i_j_k]*en_coll_param[species]
		    *(
			Neutral_DOF[species]*neutral_p[i_j_k]/neutral_n[i_j_k]
			- 2./(kappa_electron-1.) * PE[i_j_k] * rRHO[i_j_k]
			+ 2. * Neutral_Masses[species]
			  *(
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			   )
		    )
		  );
	      }
 #endif
	      
	      PE_new[i_j_k] +=  dt_step *
		(
		  0.
#if defined nonadiabatic_gradPE_SOURCE_TERM
		    //!------ source * m_e * 2/f_e * (-rot(B)/rho_i+u_i)^2 --------
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.) * mass_electron_norm
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			)
		     
		    //!------ -source * 4*m_e/f_e * (-rot(B)/rho_i+u_i)*u_n ---
		    - norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)*mass_electron_norm*2.
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *neutral_u1[i_j_k]
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *neutral_u2[i_j_k]
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *neutral_u3[i_j_k]
			)
		      
		    //!------ source * 2/f_e ( m_n*u_n^2 + f_n*p_n/(2*n_n) ) ---
// 		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)
// 		      * (
// 			      Neutral_Masses[species]
// 			      * (
// 				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
// 				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
// 				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
// 				)
// 			      + Neutral_DOF[species] * neutral_p[i_j_k] / (2. * neutral_n[i_j_k])
// 			)
// 		    
		    //!------ source * ( 2/f_e m_n*u_n^2 + T_{e,n} ) ---
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k]
		      * (
			      (kappa_electron-1.) * Neutral_Masses[species]
			      * (
				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
				)
			      + neutral_beta[i_j_k]
			)
		    
#endif
		);
	    }
#endif /* nonadiabatic_gradPE_TERM */
	    
	}
}



//!-------------------------------------------------------------//
//! LF_B_Step_conservative_CFBG_BField: -							//
//!-------------------------------------------------------------//
#if defined(nonadiabatic_gradPE_TERM) && !defined(nonadiabatic_gradPE_TERM_advance_Pe_with_B)
void CBlock::LF_Pe_Step_conservative_CFBG_BField( INT32 PE_to_update,
						  INT32 PE_derivatives,
						  D_REAL dt_step)
{
      D_REAL* PE_new = Field_Type[PE_to_update];
      D_REAL* DRAIN_e = Field_Type[id_rho_np1_recombined];
//       const D_REAL drain_e = 0.;

      
#if defined(nonadiabatic_gradPE_TERM_ensure_PE_not_negative)
      D_REAL* PE_old = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	if(PE_old[node] < 0.)
	{
	    PE_old[node] = 0.;
	}
      }
#endif


      //! - B_derivatives is equal to B_sw

	set_pointers_to_calc_derivatives_CFBG_BField(id_BEven,
						   id_BTotal,
						   PE_derivatives,
						   id_UI_minus,
						   id_rho_rez,
						   id_Eta);
#if nonadiabatic_gradPE_TERM_derivation_form == 2
      D_REAL* PE_old1 = Field_Type[PE_derivatives];
      for(INT32 node=0; node<num_nodes_in_block; ++node)
      {
	PEkappa[node] = pow( PE_old1[node] , 1./kappa_electron );
      }
#endif



      
      //! loop over 3 components and all elements
      //! ("set_pointers... " has to be called before this loop !!!)
#if defined(use_vectorclass)
	VEC4_D_REAL vec_B1_sw;
	VEC4_D_REAL vec_B1_cfbg;
	for (i=0; i < 3*num_nodes_in_block; i+=4)
	{
	  vec_B1_sw.load(B1_sw+i);
	  vec_B1_cfbg.load(B1_cfbg+i);
	  vec_B1_sw += vec_B1_cfbg;
	  vec_B1_sw.store(B1_total+i);
	}
#else
      for (i=0; i < 3*num_nodes_in_block; i++)
      B1_total[i] = B1_sw[i] +  B1_cfbg[i];
#endif

      //! derivatives are calculated from 
      //! 1) B_sw
      //! 2) B_total

      //! indices have to be i,j,k (as in calc_derivatives !!!)
      for (i=1; i < BlkNds_X-1; i++)
       for (j=1; j < BlkNds_Y-1; j++)
        for (k=1; k < BlkNds_Z-1; k++)
	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{


	    //! pointers to calc derivatives
	    //! have to be set !!!
	    calc_derivatives_for_Pe_conservative_CFBG_BField();

	    //!--------------------------------------------------------------------
            //!----------- PE -----------------------------------------------------
	    //!--------------------------------------------------------------------
            PE_new[i_j_k] +=  dt_step *
	    (

#if nonadiabatic_gradPE_TERM_derivation_form == 0
		//!------ ( -u_i + rot(B)/rho ) * grad(PE) -----------------------------------
		+ ( -U1[i_j_k] + rRHO[i_j_k] * (d2_B3 - d3_B2) ) * d1_PE
		+ ( -U2[i_j_k] + rRHO[i_j_k] * (d3_B1 - d1_B3) ) * d2_PE
		+ ( -U3[i_j_k] + rRHO[i_j_k] * (d1_B2 - d2_B1) ) * d3_PE
		
		#ifndef nonadiabatic_gradPE_TERM_set_divUe_to_zero_DEBUG
		//!------ -kappa* PE *( div(U_i) - rot(B)*grad(1/rho_i) ) --------------------
		- kappa_electron * PE[i_j_k] *
			(
				dxUX + dyUY + dzUZ
				- (d2_B3 - d3_B2) * d1_P
				- (d3_B1 - d1_B3) * d2_P
				- (d1_B2 - d2_B1) * d3_P
			)
		#endif
		
#endif
			
#if (nonadiabatic_gradPE_TERM_derivation_form == 1) || (nonadiabatic_gradPE_TERM_derivation_form == 2)
	    #if defined(kappa_electron_is_2)
		+ 2.0 * sqrt(PE[i_j_k]) *
	    #else
		+ kappa_electron * pow(PE[i_j_k],1.-1./kappa_electron) *
	    #endif
		(
		  //!------   - div( PE^(1/kappa) * U )    -----------------------------------
		  - d1_PEkappa_U1
		  - d2_PEkappa_U2
		  - d3_PEkappa_U3
		  
#ifndef nonadiabatic_gradPE_TERM_adiabatic_term_in_PE_only
		  //!------   + rot(B) * grad( pe^(1/kappa) / rho_i )    ---------------------
		  + (d2_B3 - d3_B2) * d1_PEkappa_P
		  + (d3_B1 - d1_B3) * d2_PEkappa_P
		  + (d1_B2 - d2_B1) * d3_PEkappa_P
#endif
		  
		)
			
#endif
		
			
		//!------ -drain_e * PE / rho_i -----------------------------------
#if defined nonadiabatic_gradPE_DRAIN_TERM
		- DRAIN_e[i_j_k] * PE[i_j_k] * rRHO[i_j_k]
#endif

	    );
	    
	    for(INT32 species=0; species<num_Neutral_Species; ++species)
	    {
	      neutral_n = Field_Type[id_numberdensity_neutralSpecies1 + species];
	      neutral_u1 = Field_Type[id_velocity_neutralSpecies1 + species];
	      neutral_u2 = neutral_u1 + num_nodes_in_block;
	      neutral_u3 = neutral_u2 + num_nodes_in_block;
	      neutral_p = Field_Type[id_pressure_neutralSpecies1 + species];
	      neutral_beta = Field_Type[id_new_electron_beta_neutralSpecies1 + species];
	      
	      if(neutral_n[i_j_k] < 10.*PART_REAL_PRECISION)
	      {
		continue;
	      }
	      

#if defined nonadiabatic_gradPE_COLLISION_TERM
	      if( rRHO[i_j_k] >= 10.*PART_REAL_PRECISION )
	      {
		PE_new[i_j_k] +=  dt_step *
		  (
		    //!------ 2/f_e * n_e*n_n*en_coll_param * (f_n p_n/n_n - f_e p_e/n_e + m_n (u_n-u_e)^2) ---
		    + (kappa_electron-1.) * neutral_n[i_j_k]/rRHO[i_j_k]*en_coll_param[species]
		    *(
			Neutral_DOF[species]*neutral_p[i_j_k]/neutral_n[i_j_k]
			- 2./(kappa_electron-1.) * PE[i_j_k] * rRHO[i_j_k]
			+ 2. * Neutral_Masses[species]
			  *(
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) - neutral_u1[i_j_k] )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) - neutral_u2[i_j_k] )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) - neutral_u3[i_j_k] )
			   )
		    )
		  );
	      }
 #endif
	      
	      PE_new[i_j_k] +=  dt_step *
		(
		  0.
#if defined nonadiabatic_gradPE_SOURCE_TERM
		    //!------ source * m_e * 2/f_e * (-rot(B)/rho_i+u_i)^2 --------
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.) * mass_electron_norm
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			)
		     
		    //!------ -source * 4*m_e/f_e * (-rot(B)/rho_i+u_i)*u_n ---
		    - norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)*mass_electron_norm*2.
		      * (
			      ( U1[i_j_k] - rRHO[i_j_k] * (d2_B3 - d3_B2) )
			      *neutral_u1[i_j_k]
			      +
			      ( U2[i_j_k] - rRHO[i_j_k] * (d3_B1 - d1_B3) )
			      *neutral_u2[i_j_k]
			      +
			      ( U3[i_j_k] - rRHO[i_j_k] * (d1_B2 - d2_B1) )
			      *neutral_u3[i_j_k]
			)
		      
		    //!------ source * 2/f_e ( m_n*u_n^2 + f_n*p_n/(2*n_n) ) ---
// 		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k] * (kappa_electron-1.)
// 		      * (
// 			      Neutral_Masses[species]
// 			      * (
// 				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
// 				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
// 				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
// 				)
// 			      + Neutral_DOF[species] * neutral_p[i_j_k] / (2. * neutral_n[i_j_k])
// 			)
// 		    
		    //!------ source * ( 2/f_e m_n*u_n^2 + T_{e,n} ) ---
		    + norm_IonProduction_Rate_fromNeutSpec[species]*neutral_n[i_j_k]
		      * (
			      (kappa_electron-1.) * Neutral_Masses[species]
			      * (
				  neutral_u1[i_j_k]*neutral_u1[i_j_k]
				  + neutral_u2[i_j_k]*neutral_u2[i_j_k]
				  + neutral_u3[i_j_k]*neutral_u3[i_j_k]
				)
			      + neutral_beta[i_j_k]
			)
		    
#endif
		);
	    }
	}
}
#endif



//!-------------------------------------------------------------//
//! RK_B_Step: -							//
//!-------------------------------------------------------------//
#ifndef nonadiabatic_gradPE_TERM
void CBlock::RK_B_Step_conservative_CFBG_BField()
{

      D_REAL *KX_1 = Field_Type[id_KX];
      D_REAL *KX_2 = KX_1 +num_nodes_in_block;
      D_REAL *KX_3 = KX_2 +num_nodes_in_block;



      set_pointers_to_calc_derivatives_CFBG_BField(id_BDerivative,
						   id_BTotal,
						   id_UI_minus,
						   id_rho_rez,
						   id_Eta);



      //! loop over 3 components and all elements
      //! ("set_pointers... " has to be called before this loop !!!)
#if defined(use_vectorclass)
	VEC4_D_REAL vec_B1_sw;
	VEC4_D_REAL vec_B1_cfbg;
	for (i=0; i < 3*num_nodes_in_block; i+=4)
	{
	  vec_B1_sw.load(B1_sw+i);
	  vec_B1_cfbg.load(B1_cfbg+i);
	  vec_B1_sw += vec_B1_cfbg;
	  vec_B1_sw.store(B1_total+i);
	}
#else
      for (i=0; i < 3*num_nodes_in_block; i++)
      B1_total[i] = B1_sw[i] +  B1_cfbg[i];
#endif

      //! derivatives are calculated from 
      //! 1) B_sw
      //! 2) B_total

      //! indices have to be i,j,k (as in calc_derivatives !!!)
      for (i=1; i < BlkNds_X-1; i++)
       for (j=1; j < BlkNds_Y-1; j++)
        for (k=1; k < BlkNds_Z-1; k++)
	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{


	    //! pointers to calc derivatives
	    //! have to be set !!!
	    calc_derivatives_conservative_CFBG_BField();
	    //!--------------------------------------------------------------------
            //!----------- B1 -----------------------------------------------------
	    //!--------------------------------------------------------------------

            KX_1[i_j_k] =
            (

#ifdef CONVEC_TERM
		//!------ V x (u x B) ----------------------------------------------
		+( d2_U1B2 - d2_U2B1 - d3_U3B1 + d3_U1B3)

#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------

		-(	
	
			+   B2P*d22_B3 - B3P*d33_B2
			+ d2_B2P*d2_B3 - d3_B3P*d3_B2

			+ d2_B1P_d1_b3 - d2_B1P_d3_b1
			+ d3_B1P_d2_b1 - d3_B1P_d1_b2
			+ d3_B3P_d2_b3 - d2_B2P_d3_b2
	
		)


#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------

		-(
			 d2_Eta_d1_b2 - (d2_Eta*d2_B1) - Eta*d22_B1

			+d3_Eta_d1_b3 - (d3_Eta*d3_B1) - Eta*d33_B1
		 )

#endif
	   );


	    //!--------------------------------------------------------------------
            //!----------- B2 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            KX_2[i_j_k] =
            (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		+(d3_U2B3 - d3_U3B2 - d1_U1B2 + d1_U2B1)



#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------

		-(
	
			+ B3P*d33_B1   - B1P*d11_B3
			+ d3_B3P*d3_B1 - d1_B1P*d1_B3

			+ d3_B2P_d2_b1 - d3_B2P_d1_b2 
			+ d1_B2P_d3_b2 - d1_B2P_d2_b3 
			+ d1_B1P_d3_b1 - d3_B3P_d1_b3
	
		)




#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------

		-(
			 d1_Eta_d2_b1 - (d1_Eta*d1_B2) - Eta*d11_B2
	
			+d3_Eta_d2_b3 - (d3_Eta*d3_B2) - Eta*d33_B2
		)




#endif

	    );

	    //!--------------------------------------------------------------------
            //!----------- B3 -----------------------------------------------------
	    //!--------------------------------------------------------------------
            KX_3[i_j_k] =
	    (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		+(d1_U3B1 - d1_U1B3 - d2_U2B3 + d2_U3B2)



#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		-(	
	
			+ B1P*d11_B2   - B2P*d22_B1
			+ d1_B1P*d1_B2 - d2_B2P*d2_B1

			+ d1_B3P_d3_b2 - d1_B3P_d2_b3
			+ d2_B3P_d1_b3 - d2_B3P_d3_b1
			+ d2_B2P_d1_b2 - d1_B1P_d2_b1
	
		)

#endif
#ifdef ETA_TERM_BField

		//!--------- - V x (Eta V X B)--------------------------------------
		-(
			 d1_Eta_d3_b1 - (d1_Eta*d1_B3) - Eta*d11_B3
	
			+d2_Eta_d3_b2 - (d2_Eta*d2_B3) - Eta*d22_B3
		)
#endif

	    );

	}
}
#endif /* nonadiabatic_gradPE_TERM */


//!-------------------------------------------------------------//
//! RK_B_Step: -							//
//!-------------------------------------------------------------//
#ifndef nonadiabatic_gradPE_TERM
void CBlock::RK_B_Step(INT32 KX_type, INT32 B_dev_type)
{

      D_REAL *KX_1 = Field_Type[KX_type];
      D_REAL *KX_2 = KX_1 +num_nodes_in_block;
      D_REAL *KX_3 = KX_2 +num_nodes_in_block;



      set_pointers_to_calc_derivatives(B_dev_type,
				       id_UI_minus,
				       id_rho_rez,
				       id_Eta);

      //! indices have to be i,j,k (as in calc_derivatives !!!)
      for (i=1; i < BlkNds_X-1; i++)
       for (j=1; j < BlkNds_Y-1; j++)
        for (k=1; k < BlkNds_Z-1; k++)
	if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	{


	    //! pointers to calc derivatives
	    //! have to be set !!!
	    calc_derivatives();
	    //!--------------------------------------------------------------------
            //!----------- B1_odd -------------------------------------------------
	    //!--------------------------------------------------------------------

            KX_1[i_j_k] =
            (

#ifdef CONVEC_TERM
		//!------ V x (u x B) ----------------------------------------------
		ux*dyBY + by*dyUX - uy*dyBX - bx*dyUY - uz*dzBX - bx*dzUZ + ux*dzBZ + bz*dzUX

#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		+bx*dyP*(dzBX-dxBZ) + rP*dyBX*(dzBX-dxBZ) - bx*rP*dydxBZ

		-(by*dyP+rP*dyBY)*(dyBZ-dzBY) - rP*by*(d2yBZ-dydzBY)

		-(bz*dzP + rP*dzBZ)*(dyBZ-dzBY)-rP*bz*(dzdyBZ-d2zBY)

		+bx*(dzP)*(dxBY-dyBX) + rP*(dzBX)*(dxBY-dyBX) + rP*bx*dzdxBY
#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------
		+ Eta*
		(
			d2yBX +d2zBX 
		      - dydxBY -dzdxBZ
		)
		
		-(
			dyEta * (dxBY - dyBX)
			-dzEta * (dzBX - dxBZ)
		 )
#endif
	   );


	    //!--------------------------------------------------------------------
            //!----------- B2_odd -------------------------------------------------
	    //!--------------------------------------------------------------------
            KX_2[i_j_k] =
            (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		uy*dzBZ + bz*dzUY - uz*dzBY - by*dzUZ - ux*dxBY - by*dxUX + uy*dxBX + bx*dxUY

#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		+by*dzP*(dxBY-dyBX) + rP*dzBY*(dxBY-dyBX) - by*rP*dydzBX

		-(bz*dzP+rP*dzBZ)*(dzBX-dxBZ) - rP*bz*(d2zBX-dzdxBZ)

		-(bx*dxP + rP*dxBX)*(dzBX-dxBZ)-rP*bx*(dxdzBX-d2xBZ)

		+by*(dxP)*(dyBZ-dzBY) + rP*(dxBY)*(dyBZ-dzBY) + rP*by*dydxBZ

#endif
#ifdef ETA_TERM_BField
		//!--------- - V x (Eta V X B)--------------------------------------
		+Eta*
		(
	        	d2xBY + d2zBY 
			- dxdyBX - dzdyBZ
		)

		-(
			dzEta * (dyBZ - dzBY)
			-dxEta * (dxBY - dyBX)
		)
#endif

	    );

	    //!--------------------------------------------------------------------
            //!----------- B3_odd -------------------------------------------------
	    //!--------------------------------------------------------------------
            KX_3[i_j_k] =
	    (

#ifdef CONVEC_TERM

		//!-------- V x (u x B) --------------------------------------------
		uz*dxBX + bx*dxUZ - ux*dxBZ - bz*dxUX - uy*dyBZ - bz*dyUY + uz*dyBY + by*dyUZ
#endif
#ifdef HALL_TERM
		//!------- + V x (rP (B X (V X B)))---------------------------------
		+bz*dxP*(dyBZ-dzBY) + rP*dxBZ*(dyBZ-dzBY) - bz*rP*dzdxBY

		-(bx*dxP+rP*dxBX)*(dxBY-dyBX) - rP*bx*(d2xBY-dxdyBX)

		-(by*dyP + rP*dyBY)*(dxBY-dyBX)-rP*by*(dydxBY-d2yBX)

		+bz*(dyP)*(dzBX-dxBZ) + rP*(dyBZ)*(dzBX-dxBZ) + rP*bz*dydzBX
#endif
#ifdef ETA_TERM_BField

		//!--------- - V x (Eta V X B)--------------------------------------
		+Eta*
		(
			d2xBZ +d2yBZ
			 - dxdzBX - dydzBY
		)

		-(
			dxEta * (dzBX - dxBZ)
			-dyEta * (dyBZ - dzBY)
		)
#endif

	    );

	}
}
#endif /* nonadiabatic_gradPE_TERM */

//!-------------------------------------------------------------//
//! calc_RHS_CN: -						//
//!-------------------------------------------------------------//
void CBlock::calc_RHS_CN(INT32 type_RHS, INT32 type_BField)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	//! first rd of respective level
	 rd = rd_of_L[RLevel];
	rd2 = rd2_of_L[RLevel];

	D_REAL rdt_BLoops = (1. *num_advance_B_loops) / dt;

	ETA = Field_Type[id_Eta];

	B1 = Field_Type[type_BField];
	B2 = B1 +num_nodes_in_block;
	B3 = B2 +num_nodes_in_block;

	D_REAL* RHS1 = Field_Type[type_RHS];
	D_REAL* RHS2 = RHS1 +num_nodes_in_block;
	D_REAL* RHS3 = RHS2 +num_nodes_in_block;

	memset(RHS1, 0, 3*num_nodes_in_block*sizeof(D_REAL));

	//! indices have to be i,j,k (as in calc_derivatives !!!)
	for (i=1; i < BlkNds_X-1; i++)
	 for (j=1; j < BlkNds_Y-1; j++)
	  for (k=1; k < BlkNds_Z-1; k++)
	  if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	  {

		i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
		//! --- indices for firs derivatives ---------
		ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);


		//! --- indices for mixed derivatives ---------
		ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
		im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
		
		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
		i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
		i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
		i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);


		//!------ Eta - derivatives ------
		Eta    = ETA[i_j_k];
		
		dxEta = rd[0]*(ETA[ip1_j_k]-ETA[im1_j_k]);
		dyEta = rd[1]*(ETA[i_jp1_k]-ETA[i_jm1_k]);
		dzEta = rd[2]*(ETA[i_j_kp1]-ETA[i_j_km1]);

		//!------ first derivatives ---
		dxBY = rd[0]*(B2[ip1_j_k]-B2[im1_j_k]);
		dxBZ = rd[0]*(B3[ip1_j_k]-B3[im1_j_k]);
			
		dyBX = rd[1]*(B1[i_jp1_k]-B1[i_jm1_k]);
		dyBZ = rd[1]*(B3[i_jp1_k]-B3[i_jm1_k]);
			
		dzBX = rd[2]*(B1[i_j_kp1]-B1[i_j_km1]);
		dzBY = rd[2]*(B2[i_j_kp1]-B2[i_j_km1]);

		//!------ second derivatives --(just mixed occure) ----------
		d2xBY = rd2[0]*(B2[ip1_j_k]-2.*B2[i_j_k]+B2[im1_j_k]);
		d2xBZ = rd2[0]*(B3[ip1_j_k]-2.*B3[i_j_k]+B3[im1_j_k]);
		
		d2yBX = rd2[1]*(B1[i_jp1_k]-2.*B1[i_j_k]+B1[i_jm1_k]);
		d2yBZ = rd2[1]*(B3[i_jp1_k]-2.*B3[i_j_k]+B3[i_jm1_k]);
		
		d2zBX = rd2[2]*(B1[i_j_kp1]-2.*B1[i_j_k]+B1[i_j_km1]);
		d2zBY = rd2[2]*(B2[i_j_kp1]-2.*B2[i_j_k]+B2[i_j_km1]);


		//!------ mixed derivatives ---------------------------------
		dxdzBX = rd[0]*rd[2] *(B1[ip1_j_kp1]-B1[ip1_j_km1]
				      -B1[im1_j_kp1]+B1[im1_j_km1]);
		
		dxdyBX = rd[0]*rd[1] *(B1[ip1_jp1_k]-B1[ip1_jm1_k]
				      -B1[im1_jp1_k]+B1[im1_jm1_k]);
		

		
		
		dydxBY = rd[1]*rd[0] *(B2[ip1_jp1_k]-B2[ip1_jm1_k]
				      -B2[im1_jp1_k]+B2[im1_jm1_k]);
		
		
		dydzBY = rd[1]*rd[2] *(B2[i_jp1_kp1]-B2[i_jp1_km1]
				      -B2[i_jm1_kp1]+B2[i_jm1_km1]);
		
		
		dzdxBZ = rd[2]*rd[0] *(B3[ip1_j_kp1]-B3[ip1_j_km1]
				      -B3[im1_j_kp1]+B3[im1_j_km1]);
		
		
		dzdyBZ = rd[2]*rd[1] *(B3[i_jp1_kp1]-B3[i_jm1_kp1]
				      -B3[i_jp1_km1]+B3[i_jm1_km1]);



		//!--------- - V x (Eta V X B)--------------------------------------
		RHS1[i_j_k] =	B1[i_j_k]*rdt_BLoops
				+0.5*Eta*
				(
					 d2yBX +d2zBX 
					-dydxBY -dzdxBZ
				)
		
				-0.5*(
					 dyEta * (dxBY - dyBX)
					-dzEta * (dzBX - dxBZ)
				);

		//!--------- - V x (Eta V X B)--------------------------------------
		RHS2[i_j_k] =	B2[i_j_k]*rdt_BLoops
				+0.5*Eta*
				(
					  d2xBY + d2zBY 
					- dxdyBX - dzdyBZ
				)
		
				-0.5*(
					 dzEta * (dyBZ - dzBY)
					-dxEta * (dxBY - dyBX)
				);

		//!--------- - V x (Eta V X B)--------------------------------------
		RHS3[i_j_k] =	B3[i_j_k]*rdt_BLoops
				+0.5*Eta*
				(
					 d2xBZ +d2yBZ
					-dxdzBX - dydzBY
				)
		
				-0.5*(
					 dxEta * (dzBX - dxBZ)
					-dyEta * (dyBZ - dzBY)
				);


	  }

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;



}


//!-------------------------------------------------------------//
//! Obstacle_BField_SOR_Step: -						//
//!-------------------------------------------------------------//
D_REAL CBlock::Obstacle_BField_SOR_Step(INT32 counter, INT32 B_type, INT32 type_RHS, bool Odd_Even_Step)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	//! first rd of respective level
	 rd = rd_of_L[RLevel];
	rd2 = rd2_of_L[RLevel];

	D_REAL rdt_BLoops = (1. *num_advance_B_loops) / dt;

	ETA = Field_Type[id_Eta];

	B1 = Field_Type[B_type];
	B2 = B1 +num_nodes_in_block;
	B3 = B2 +num_nodes_in_block;


	D_REAL B_SOR_omega = B_SOR_omega_L[RLevel];

	if(num_optimize_SOR_steps)
	{	
		if(TL>=TLstart_optimize_B_SOR_omega && TL<TLstart_optimize_B_SOR_omega + num_optimize_SOR_steps)
		B_SOR_omega = 	B_SOR_omega_L[RLevel] - 0.01*num_optimize_SOR_steps/2 + 0.01*(TL - TLstart_optimize_B_SOR_omega);
		if(B_SOR_omega<1)
			B_SOR_omega = 1;
		if(B_SOR_omega>2)
			B_SOR_omega = 2;
	}

	D_REAL* RHS1 = Field_Type[type_RHS];
	D_REAL* RHS2 = RHS1 +num_nodes_in_block;
	D_REAL* RHS3 = RHS2 +num_nodes_in_block;


	D_REAL alphaX = 0.;
	D_REAL alphaY = 0.;
	D_REAL alphaZ = 0.;

	//! indices have to be i,j,k (as in calc_derivatives !!!)
	for (i=1; i < BlkNds_X-1; i++)
	 for (j=1; j < BlkNds_Y-1; j++)
	  for (k=1; k < BlkNds_Z-1; k++)
	  if((i+j+k)%2 == Odd_Even_Step)
	  if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	  {
	
		
		  i_j_k =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
		//! --- indices for firs derivatives ---------
		ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);


		//! --- indices for mixed derivatives ---------
		ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
		im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
		
		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
		i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
		i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
		i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);


		//!------ Eta - derivatives ------
		Eta    = ETA[i_j_k];
		
		dxEta = rd[0]*(ETA[ip1_j_k]-ETA[im1_j_k]);
		dyEta = rd[1]*(ETA[i_jp1_k]-ETA[i_jm1_k]);
		dzEta = rd[2]*(ETA[i_j_kp1]-ETA[i_j_km1]);

		//!------ first derivatives ---
		dxBY = rd[0]*(B2[ip1_j_k]-B2[im1_j_k]);
		dxBZ = rd[0]*(B3[ip1_j_k]-B3[im1_j_k]);
			
		dyBX = rd[1]*(B1[i_jp1_k]-B1[i_jm1_k]);
		dyBZ = rd[1]*(B3[i_jp1_k]-B3[i_jm1_k]);
			
		dzBX = rd[2]*(B1[i_j_kp1]-B1[i_j_km1]);
		dzBY = rd[2]*(B2[i_j_kp1]-B2[i_j_km1]);



		//!------ mixed derivatives ---------------------------------
		dxdzBX = rd[0]*rd[2] *(B1[ip1_j_kp1]-B1[ip1_j_km1]
				      -B1[im1_j_kp1]+B1[im1_j_km1]);
		
		dxdyBX = rd[0]*rd[1] *(B1[ip1_jp1_k]-B1[ip1_jm1_k]
				      -B1[im1_jp1_k]+B1[im1_jm1_k]);
		

		
		
		dydxBY = rd[1]*rd[0] *(B2[ip1_jp1_k]-B2[ip1_jm1_k]
				      -B2[im1_jp1_k]+B2[im1_jm1_k]);
		
		
		dydzBY = rd[1]*rd[2] *(B2[i_jp1_kp1]-B2[i_jp1_km1]
				      -B2[i_jm1_kp1]+B2[i_jm1_km1]);
		
		
		

		dzdxBZ = rd[2]*rd[0] *(B3[ip1_j_kp1]-B3[ip1_j_km1]
				      -B3[im1_j_kp1]+B3[im1_j_km1]);
		
		
		dzdyBZ = rd[2]*rd[1] *(B3[i_jp1_kp1]-B3[i_jm1_kp1]
				      -B3[i_jp1_km1]+B3[i_jm1_km1]);

		alphaX = B_SOR_omega / (rdt_BLoops + Eta *(rd2[1] +rd2[2]));
		alphaY = B_SOR_omega / (rdt_BLoops + Eta *(rd2[0] +rd2[2]));
		alphaZ = B_SOR_omega / (rdt_BLoops + Eta *(rd2[0] +rd2[1]));

	

		B1[i_j_k] = 	(1.-B_SOR_omega) * B1[i_j_k]

			       +alphaX *
			        (
				    //!-- 0.5 * V x (Eta V X B) ----
				    +0.5*
				    (
				     Eta*(
						  (B1[i_jp1_k] +B1[i_jm1_k]) *rd2[1]
						+ (B1[i_j_kp1] +B1[i_j_km1]) *rd2[2]
		
	
						- dydxBY -dzdxBZ
					 )
					-(
						 dyEta * (dxBY - dyBX)
						-dzEta * (dzBX - dxBZ)
					 )
				    )
				    +RHS1[i_j_k]

			        );

		B2[i_j_k] = 	(1.-B_SOR_omega) * B2[i_j_k]

			       +alphaY *
			        (
				    //!-- 0.5 * V x (Eta V X B) ----
				    +0.5*
				    (
				     Eta*(

						  (B2[ip1_j_k] +B2[im1_j_k]) *rd2[0]
						+ (B2[i_j_kp1] +B2[i_j_km1]) *rd2[2]
		
	
						- dxdyBX - dzdyBZ
					 )

					-(
						 dzEta * (dyBZ - dzBY)
						-dxEta * (dxBY - dyBX)
					 )
				     )
				     +RHS2[i_j_k]

			        );

		B3[i_j_k] = 	(1.-B_SOR_omega) * B3[i_j_k]

			       +alphaZ *
			        (
				    //!-- 0.5 * V x (Eta V X B) ----
				    +0.5*
				    (
				     Eta*(
	
						  (B3[ip1_j_k] +B3[im1_j_k]) *rd2[0]
						+ (B3[i_jp1_k] +B3[i_jm1_k]) *rd2[1]
		
	
						- dxdzBX - dydzBY
					 )

					-(
						 dxEta * (dzBX - dxBZ)
						-dyEta * (dyBZ - dzBY)
					 )
			
				     )
				     +RHS3[i_j_k]

			        );
	
	
	
	}


	D_REAL total_error_X = 0.;
	D_REAL total_error_Y = 0.;
	D_REAL total_error_Z = 0.;

	if(counter%B_SOR_calc_error_step==0)
	{

	      D_REAL error_X = 0.;
	      D_REAL error_Y = 0.;
	      D_REAL error_Z = 0.;



	      for (i=1; i < BlkNds_X-1; i++)
	       for (j=1; j < BlkNds_Y-1; j++)
		for (k=1; k < BlkNds_Z-1; k++)
		if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
		{
		

		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
			
			//! --- indices for first derivatives ---------
			ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
			im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
			
			i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
			i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
			
			i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
			i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
	
	
			//! --- indices for mixed derivatives ---------
			ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
			ip1_j_km1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
			im1_j_kp1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
			im1_j_km1 = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
			
			ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
			ip1_jm1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
			im1_jp1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
			im1_jm1_k = (i-1)*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
			
			i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
			i_jp1_km1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k-1);
			i_jm1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k+1);
			i_jm1_km1 =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z +(k-1);
	

			//!------ Eta - derivatives ------
			  Eta = ETA[i_j_k];
			
			dxEta = rd[0]*(ETA[ip1_j_k]-ETA[im1_j_k]);
			dyEta = rd[1]*(ETA[i_jp1_k]-ETA[i_jm1_k]);
			dzEta = rd[2]*(ETA[i_j_kp1]-ETA[i_j_km1]);
	
			//!------ first derivatives ---
			dxBY = rd[0]*(B2[ip1_j_k]-B2[im1_j_k]);
			dxBZ = rd[0]*(B3[ip1_j_k]-B3[im1_j_k]);
			
			dyBX = rd[1]*(B1[i_jp1_k]-B1[i_jm1_k]);
			dyBZ = rd[1]*(B3[i_jp1_k]-B3[i_jm1_k]);
			
			dzBX = rd[2]*(B1[i_j_kp1]-B1[i_j_km1]);
			dzBY = rd[2]*(B2[i_j_kp1]-B2[i_j_km1]);

			//!------ second derivatives --(just mixed occure) ----------
			d2xBY = rd2[0]*(B2[ip1_j_k]-2.*B2[i_j_k]+B2[im1_j_k]);
			d2xBZ = rd2[0]*(B3[ip1_j_k]-2.*B3[i_j_k]+B3[im1_j_k]);
			
			d2yBX = rd2[1]*(B1[i_jp1_k]-2.*B1[i_j_k]+B1[i_jm1_k]);
			d2yBZ = rd2[1]*(B3[i_jp1_k]-2.*B3[i_j_k]+B3[i_jm1_k]);
			
			d2zBX = rd2[2]*(B1[i_j_kp1]-2.*B1[i_j_k]+B1[i_j_km1]);
			d2zBY = rd2[2]*(B2[i_j_kp1]-2.*B2[i_j_k]+B2[i_j_km1]);

	
	
			//!------ mixed derivatives ---------------------------------
			dxdzBX = rd[0]*rd[2] *(B1[ip1_j_kp1]-B1[ip1_j_km1]
					      -B1[im1_j_kp1]+B1[im1_j_km1]);
			
			dxdyBX = rd[0]*rd[1] *(B1[ip1_jp1_k]-B1[ip1_jm1_k]
					      -B1[im1_jp1_k]+B1[im1_jm1_k]);
			
	
			
			
			dydxBY = rd[1]*rd[0] *(B2[ip1_jp1_k]-B2[ip1_jm1_k]
					      -B2[im1_jp1_k]+B2[im1_jm1_k]);
			
			
			dydzBY = rd[1]*rd[2] *(B2[i_jp1_kp1]-B2[i_jp1_km1]
					      -B2[i_jm1_kp1]+B2[i_jm1_km1]);
			
			
			
	
			dzdxBZ = rd[2]*rd[0] *(B3[ip1_j_kp1]-B3[ip1_j_km1]
					      -B3[im1_j_kp1]+B3[im1_j_km1]);
			
			
			dzdyBZ = rd[2]*rd[1] *(B3[i_jp1_kp1]-B3[i_jm1_kp1]
					      -B3[i_jp1_km1]+B3[i_jm1_km1]);
	
			//! ---------- Error X ------------------------
			error_X = B1[i_j_k]*rdt_BLoops
				  -0.5*Eta*
				  (
					 d2yBX +d2zBX 
					-dydxBY -dzdxBZ
				  )
		
				  +0.5*(
					 dyEta * (dxBY - dyBX)
					-dzEta * (dzBX - dxBZ)
				  )
				  -RHS1[i_j_k];

			total_error_X += error_X*error_X;
	
			//! ---------- Error Y ------------------------
			error_Y = B2[i_j_k]*rdt_BLoops
				  -0.5*Eta*
				  (
					  d2xBY + d2zBY 
					- dxdyBX - dzdyBZ
				  )
		
				  +0.5*(
					 dzEta * (dyBZ - dzBY)
					-dxEta * (dxBY - dyBX)
				  )
				  -RHS2[i_j_k];

			total_error_Y += error_Y*error_Y;
	
			//! ---------- Error Z ------------------------
			error_Z = B3[i_j_k]*rdt_BLoops
				  -0.5*Eta*
				  (
					 d2xBZ +d2yBZ
					-dxdzBX - dydzBY
				  )
		
				  +0.5*(
					 dxEta * (dzBX - dxBZ)
					-dyEta * (dyBZ - dzBY)
				  )
				  -RHS3[i_j_k];

			total_error_Z += error_Z*error_Z;
		
		
		
		}

		total_error_X = sqrt(total_error_X)/(BlkNds_X*BlkNds_Y*BlkNds_Z);
		total_error_Y = sqrt(total_error_Y)/(BlkNds_X*BlkNds_Y*BlkNds_Z);
		total_error_Z = sqrt(total_error_Z)/(BlkNds_X*BlkNds_Y*BlkNds_Z);


	}

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

	return (total_error_X +total_error_Y +total_error_Z);
}

//!-------------------------------------------------------------//
//! Obstacle_BField_SOR_Step: -						//
//!-------------------------------------------------------------//
D_REAL CBlock::Psi_DC_SOR_Step(INT32 counter, INT32 div_type, bool Odd_Even_Step)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	D_REAL DC_omega = DC_omega_L[RLevel];

	//! rd2 of respective level
	rd2 = rd2_of_L[RLevel];

	if(num_optimize_SOR_steps)
	{	
		if(TL>=TLstart_optimize_B_SOR_omega && TL<TLstart_optimize_B_SOR_omega + num_optimize_SOR_steps)
		DC_omega = DC_omega_L[RLevel] - 0.01*num_optimize_SOR_steps/2 + 0.01*(TL - TLstart_optimize_B_SOR_omega);
		if(DC_omega<1)
			DC_omega = 1;
		if(DC_omega>2)
			DC_omega = 2;
	}
	
	
	//! precal alpha
	D_REAL alpha_phi = DC_omega/(2.*(rd2[0]+rd2[1]+rd2[2]));


	D_REAL* Phi = Field_Type[id_PhiDC];
	D_REAL* divB = Field_Type[div_type];



	//! indices have to be i,j,k (as in calc_derivatives !!!)
	for (i=1; i < BlkNds_X-1; i++)
	 for (j=1; j < BlkNds_Y-1; j++)
	  for (k=1; k < BlkNds_Z-1; k++)
	  if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	  if((i+j+k)%2 == Odd_Even_Step)
	  {
	
		
		  i_j_k =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
		//! --- indices for firs derivatives ---------
		ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);


		Phi[i_j_k] = (1.-DC_omega) * Phi[i_j_k]
			
			    +alpha_phi *
			    (

				 (Phi[ip1_j_k] +Phi[im1_j_k])*rd2[0]
				+(Phi[i_jp1_k] +Phi[i_jm1_k])*rd2[1]
				+(Phi[i_j_kp1] +Phi[i_j_km1])*rd2[2]

				-divB[i_j_k]
			    );

	  }



	D_REAL total_error_Phi = 0.;

	if(counter%DC_calc_error_step==0)
	{


	    D_REAL error_Phi = 0.;

	    for (i=1; i < BlkNds_X-1; i++)
	     for (j=1; j < BlkNds_Y-1; j++)
	      for (k=1; k < BlkNds_Z-1; k++)
	      if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	      {
	
		
			 i_j_k =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
			
			//! --- indices for firs derivatives ---------
			ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
			im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
			
			i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
			i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
			
			i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
			i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);
	
	
			error_Phi =  (Phi[ip1_j_k] -2.*Phi[i_j_k] +Phi[im1_j_k])*rd2[0]
				+(Phi[i_jp1_k] -2.*Phi[i_j_k] +Phi[i_jm1_k])*rd2[1]
				+(Phi[i_j_kp1] -2.*Phi[i_j_k] +Phi[i_j_km1])*rd2[2]
	
				-divB[i_j_k];
	
			total_error_Phi += error_Phi*error_Phi;



	      }

	      total_error_Phi = sqrt(total_error_Phi)/(BlkNds_X*BlkNds_Y*BlkNds_Z);


	}


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


	return total_error_Phi;

}



//!-------------------------------------------------------------//
//! Obstacle_BField_SOR_Step: -						//
//!-------------------------------------------------------------//
void CBlock::project_divB_out_of_BField(INT32 type_BField)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	//! rd2 of respective level
	rd = rd_of_L[RLevel];


	D_REAL* Phi = Field_Type[id_PhiDC];

	B1 = Field_Type[type_BField];
	B2 = B1 +num_nodes_in_block;
	B3 = B2 +num_nodes_in_block;



	//! indices have to be i,j,k (as in calc_derivatives !!!)
	for (i=1; i < BlkNds_X-1; i++)
	 for (j=1; j < BlkNds_Y-1; j++)
	  for (k=1; k < BlkNds_Z-1; k++)
	  if(!Core_Flag[i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k])
	  {
	
		
		  i_j_k =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
		//! --- indices for firs derivatives ---------
		ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		im1_j_k = (i-1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
		
		i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
		i_jm1_k =     i*BlkNds_Y*BlkNds_Z +(j-1)*BlkNds_Z     +k;
		
		i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
		i_j_km1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k-1);


		B1[i_j_k] -= (Phi[ip1_j_k] - Phi[im1_j_k])*rd[0];
		B2[i_j_k] -= (Phi[i_jp1_k] - Phi[i_jm1_k])*rd[1];
		B3[i_j_k] -= (Phi[i_j_kp1] - Phi[i_j_km1])*rd[2];
			

	  }


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;



}
