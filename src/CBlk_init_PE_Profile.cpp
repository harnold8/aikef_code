




#include "CBlk.h"
#include "parameters.h"
#include "utils.h"



#include <iostream>
#include <fstream>
#include <math.h>
#include <cassert>


using namespace std;


// added by HR
#if defined nonadiabatic_gradPE_TERM
//!-------------------------------------------------------------//
//! init_PE_Profile: 									//
//!-------------------------------------------------------------//
void CBlock::init_PE_Profile(void)
{



	D_REAL r_vec[3];
	D_REAL null_vec[3] = {0.,0.,0.};



	D_REAL* PE = Field_Type[id_PEtotal];
	D_REAL* RHO = Field_Type[id_rho_n];

	
	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		INT32 u_v_w  = u*BlkNds_Y*BlkNds_Z 
			      +v*BlkNds_Z 
			      +w;
		

		INT32 cell_indices[3] = {u,v,w};

		assert(mesh_type == UNIFORM);
		intern2normedCoords(r_vec, null_vec, cell_indices);


		D_REAL dist_to_orig = vec_len(r_vec);

		//! see comment in defines.h for details about
		//! nonadiabatic_gradPE_initial_condition
		if(nonadiabatic_gradPE_initial_condition == 1)
		{
		  //! set PE value as p_e = \beta_e
		  PE[u_v_w] = Electron_Betas[0];
		}
		else
		{
		  //! set PE value as p_e = n_e k T_e = \rho_i \beta_e
		  PE[u_v_w] = RHO[u_v_w] * Electron_Betas[0];
		}

		// the rest of CBlk::init_Eta_Profile
		/*
		//! better use fermi function than define L0 & smooth:
		//! In case only defined on L0, Eta has to be interpolated
		//! to higher levels which lead to "rectangular shape".
	
		//! NOTE:
		//! In case "smooth ETA" version is used, do not forget to
		//! include field_from_parent(Eta) function in refine_Block !
		if(dist_to_orig < R_Eta_Obstacle && !use_eta_fermi)
		Eta[u_v_w] = Eta_Obstacle;

                if(dist_to_orig < 0.85*R_Eta_Obstacle && !use_eta_fermi)
                Eta[u_v_w] = 1.;
	
		if(use_eta_fermi)
		Eta[u_v_w]= Eta_Obstacle *1./(exp((dist_to_orig -R_Eta_Obstacle)*fermi_slope)+1.);
	

		//! eta should be zero inside obstacle core
		if(dist_to_orig < obstacle_core_fraction* R_Obstacle)
		Eta[u_v_w] = 0.;
		*/

	}

	// the rest of CBlk::init_Eta_Profile
	/*
	//! set Eta boundaries
	if(eta_Alfven_Wing_boundary)
	for(INT32 u=0; u< BlkNds_X; u++)
	  for(INT32 v=0; v< BlkNds_Y; v++)
	    for(INT32 w=0; w< BlkNds_Z; w++)
	    {

		INT32 u_v_w  = u*BlkNds_Y*BlkNds_Z 
			     +v*BlkNds_Z 
			     +w;
		


		INT32 cell_indices[3] = {u,v,w};
		intern2normedCoords(r_vec, null_vec, cell_indices);
		
		D_REAL eta ;
		
		//! upper boundary
		eta = resistive_bound_eta[4]/resistive_bound_dist[4] *( r_vec[2] - ( LZ - Box_Origin[2] - resistive_bound_dist[4]));
		
		if(eta>Eta_sw)
		Eta[u_v_w] = eta;
		
		//! lower boundary
		eta = resistive_bound_eta[5]/resistive_bound_dist[5] *( -r_vec[2] + ( -Box_Origin[2] + resistive_bound_dist[5]) );
		
		if(eta>Eta_sw)
		Eta[u_v_w] = eta;
		
	
	    }


	//! set Eta boundaries
	//!----------- -X boundary ------------------------
	if(is_box_boundary[_Im1_] && resistive_bound_cells[_Im1_])
	 for(INT32 u=0; u<resistive_bound_cells[_Im1_]; u++)
	  for(INT32 v=0; v<BlkNds_Y; v++)
	   for(INT32 w=0; w<BlkNds_Z; w++)
	   {
		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Im1_];

	    }



	//!----------- +X boundary ------------------------
	if(is_box_boundary[_Ip1_] && resistive_bound_cells[_Ip1_])
	 for(INT32 u=BlkNds_X-resistive_bound_cells[_Ip1_]; u<BlkNds_X; u++)
	  for(INT32 v=0;					     v<BlkNds_Y; v++)
	   for(INT32 w=0;				      w<BlkNds_Z; w++)
	   {
		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Ip1_];

	   }


	//!----------- -Y boundary ------------------------
	if(is_box_boundary[_Jm1_]  && resistive_bound_cells[_Jm1_])
	 for(INT32 u=0; u<BlkNds_X; u++)
	  for(INT32 v=0; v<resistive_bound_cells[_Jm1_]; v++)
	   for(INT32 w=0; w<BlkNds_Z; w++)
	   {
		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Jm1_];

	    }

	//!----------- +Y boundary ------------------------
	if(is_box_boundary[_Jp1_]  && resistive_bound_cells[_Jp1_])
	 for(INT32 u=0;					   u<BlkNds_X; u++)
	  for(INT32 v=BlkNds_Y-resistive_bound_cells[_Jp1_]; v<BlkNds_Y; v++)
	   for(INT32 w=0;				      w<BlkNds_Z; w++)
	   {

		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Jp1_];

	   }

	//!----------- -Z boundary ------------------------
	if(is_box_boundary[_Km1_]  && resistive_bound_cells[_Km1_])
	 for(INT32 u=0; u<BlkNds_X; u++)
	  for(INT32 v=0; v<BlkNds_Y; v++)
	   for(INT32 w=0; w<resistive_bound_cells[_Km1_]; w++)
	   {

		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Km1_];

	  }


	//!----------- +Z boundary ------------------------
	if(is_box_boundary[_Kp1_]  && resistive_bound_cells[_Kp1_])
	 for(INT32 u=0;					    u<BlkNds_X; u++)
	  for(INT32 v=0;					     v<BlkNds_Y; v++)
	   for(INT32 w=BlkNds_Z -resistive_bound_cells[_Kp1_]; w<BlkNds_Z; w++)
	   {

		INT32 u_v_w =  u*BlkNds_Y*BlkNds_Z
			     +v*BlkNds_Z
			     +w;

		Eta[u_v_w] = resistive_bound_eta[_Kp1_];

	   }

	*/



}
#endif


	