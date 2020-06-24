
//! some abriviations for particle gathering
#define noRHO   false
#define getRHO  true

#define noJI   false
#define getJI  true

#define noLAM_GAM  false
#define getLAM_GAM  true

#define noVTH2  false
#define getVTH2  true

#define EVERY_PARTICLE false


#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

#include "CBlk.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"



extern INT32 i,j,k;

extern INT32  i_j_k, ip1_jp1_kp1;
extern INT32  ip1_j_k, i_jp1_k,  i_j_kp1;
extern INT32  ip1_jp1_k, ip1_j_kp1,  i_jp1_kp1;

extern INT32 u,v,w, u_v_w, up1_v_w, u_vp1_w, u_v_wp1;
extern INT32 up1_vp1_w, up1_v_wp1, u_vp1_wp1, up1_vp1_wp1;

extern INT32 parent_Xstart, parent_Ystart, parent_Zstart;


extern D_REAL *q_of, *q2m_of, *q_m_ratio;
extern D_REAL *CellVol_of_L;

//! Global variables of this File


//!-------------------------------------------------------------//
//! collect_Ji:									//
//!-------------------------------------------------------------//
void CBlock::collect_Ji(INT32 species_to_collect,
			INT32 total_Ji_type,
			INT32 species_Ji_type,
			INT32 species_vth2_type,
			bool collect_vth2)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();




	D_REAL *Rho = NULL;
	D_REAL *total_JiX   = Field_Type[total_Ji_type];
	D_REAL *species_JiX = Field_Type[species_Ji_type];
	
	
	//! only set species_vth2 in case vth2rho shall be collected
	if(collect_vth2)
	species_JiX = Field_Type[species_vth2_type];
	
	

	//! ALL_SPECIES does not work with collect_vth2
	//! (all velocities would be required at the same time)
	//! -> 1) sum vth2rho of species up afterwards
	//! -> 2) use silo def for this
	if(species_to_collect == ALL_SPECIES)
	{
	
		
		for(short species=0; species<num_Particle_Species; species++)
		{
		
			memset(species_JiX,0,3*num_nodes_in_block*sizeof(D_REAL));
		
			//! Gather Block
			if(is_gatherBlk)
			gather_from_parent_Block(species, noRHO, getJI, noVTH2, NULL, species_JiX);
			else
			{
				for(INT32 oct=0; oct<8; oct++)
				if(!child_array[oct])
				gather_Oct(species, noRHO, getJI, noVTH2, EVERY_PARTICLE, oct, NULL, species_JiX);
			}
	
			for(i = 0; i < 3*num_nodes_in_block; i++)
			total_JiX[i] += q_of[species] * species_JiX[i];
		
		}
	}
	else
	{

	
		//! Gather Block
		if(is_gatherBlk)
		gather_from_parent_Block(species_to_collect, noRHO, getJI, collect_vth2, Rho, species_JiX);
		else
		{
			for(INT32 oct=0; oct<8; oct++)
			if(!child_array[oct])
			gather_Oct(species_to_collect, noRHO, getJI, collect_vth2, EVERY_PARTICLE, oct, Rho, species_JiX);
		}

		//! Do not use charge in case of temperature calculation.
		//! Never apply charge to rho in this function.
		//! It has to be applied to Ji in order to get correct J2u conversion
		if(!collect_vth2)
		for(i = 0; i < 3*num_nodes_in_block; i++)
		species_JiX[i] *= q_of[species_to_collect];
	
	
	}
	

	//! Also it is important to copy currents at PLUS boundaries, oth. GN at
	//! level boundaries will fail, since zero values in lower levels are used
	//! IT CAN BE SEEN AT BOX BOUNDARIES THAT ARE LEVEL BOUNDARIES AT THE SAME TIME !!!

	//! In converJ2u these boundaries will be replaced
	//! by "correct" inflow boundaries
	bool cp_plusJiBounds[6] = {1,0,1,0,1,0};
	cp_Boundaries(total_Ji_type, cp_plusJiBounds);

	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;
	
}


//!-----------------------------------------------------------
//! collect_rho_Ji_LAM_GAM
//!-----------------------------------------------------------
void CBlock::collect_rho_Ji_LAM_GAM(INT32 Ji_type,
				    INT32 Gam_type,
				    INT32 species_Ji_type)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	D_REAL *rho_species, *rho_species_array;
	D_REAL *JiX, *GamX, *species_JiX;
	
	JiX	          = Field_Type[Ji_type];
	GamX              = Field_Type[Gam_type];
	species_JiX       = Field_Type[species_Ji_type];
	rho_species_array = Field_Type[id_allRhoSpecies];
	
	
	for(short species=0; species<num_Particle_Species; species++)
	{
	
		//! set pointer to respective species
		rho_species = rho_species_array +species*num_nodes_in_block;
		memset(species_JiX    ,0,3*num_nodes_in_block*sizeof(D_REAL));
	
	
		//! Box
		if(is_gatherBlk)
		gather_from_parent_Block(species, getRHO, getJI, noVTH2, rho_species, species_JiX);
		else
		{
			for(INT32 oct=0; oct<8; oct++)
			if(!child_array[oct])
			gather_Oct(species, getRHO, getJI, noVTH2, EVERY_PARTICLE, oct, rho_species, species_JiX);
		}
	
		//! total rho and Lam are built in "sum_up_rhoSpecies" function
		//! then also charge is applied to rho
		for(i = 0; i < 3*num_nodes_in_block; i++)
		{
	
			JiX[i] +=   q_of[species] * species_JiX[i];
			GamX[i] += q2m_of[species] * species_JiX[i];
		}
	}


	//! Apply rho boundaries here AND (!!!) after ghost nodes are copied !!!
	//! 1)They have to be copied here to get interpolation from lower to higher
	//!   levels right
	//! 2)They have to be copied after GN update to get complete boundaries, else
	//!   edges will be incomplete 
 	bool cp_allBounds[6] = {0,0,0,0,0,0};
 	cp_Boundaries(id_allRhoSpecies, cp_allBounds);

	//! Also it is important to copy currents at PLUS boundaries, oth. GN at
	//! level boundaries will fail, since zero values in lower levels are used
	//! IT CAN BE SEEN AT BOX BOUNDARIES THAT ARE LEVEL BOUNDARIES AT THE SAME TIME !!!
 	bool cp_plusJiBounds[6] = {1,0,1,0,1,0};
 	cp_Boundaries(Ji_type, cp_plusJiBounds);


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}


//!-----------------------------------------------------------
//! collect_rho_Ji
//!-----------------------------------------------------------
void CBlock::collect_rho_Ji(INT32 species_to_collect,
			    INT32 species_rho_type,
			    INT32 species_Ji_type,
			    bool skip_unmarked_particle)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	D_REAL *rho, *JiX;
	
	rho = Field_Type[species_rho_type];
	JiX = Field_Type[species_Ji_type];

	
	memset(rho    ,0,1*num_nodes_in_block*sizeof(D_REAL));
	memset(JiX    ,0,3*num_nodes_in_block*sizeof(D_REAL));


	//! Box
	if(is_gatherBlk)
	gather_from_parent_Block(species_to_collect, getRHO, getJI, noVTH2, rho, JiX);
	else
	for(INT32 oct=0; oct<8; oct++)
	if(!child_array[oct])
	gather_Oct(species_to_collect, getRHO, getJI, noVTH2, skip_unmarked_particle, oct, rho, JiX);
	



	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}

//!----------------------------------------------------
//! sum_up_rhoSpecies:
//! rho total and Lam are build from all rho species 
//! so all rho species have to be update no !!!
//!----------------------------------------------------
void CBlock::sum_up_rhoSpecies(INT32 rho_type,
                           	INT32 Lam_type)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();



	D_REAL* rho_species;
	D_REAL* rho = Field_Type[rho_type];
	D_REAL* Lam = Field_Type[Lam_type];
	D_REAL* rho_species_array = Field_Type[id_allRhoSpecies];
	
	memset(rho,0,num_nodes_in_block*sizeof(D_REAL));
	memset(Lam,0,num_nodes_in_block*sizeof(D_REAL));




	//! It is important to use the correct boundary value for each rho species,
	//! otw. the gradPE will be calculated not correct which leads to an incorrect
	//! EField. As a consequence errors will rise from the simulation box boundaries.
	//! The following call copies at each BoxBounfary the penultimate layer to the 
	//! ultimate in order to eliminate gradients in EField calculation.
	//! Copy boundaries before summing up !!
 	bool cp_allBounds[6] = {0,0,0,0,0,0};
 	cp_Boundaries(id_allRhoSpecies, cp_allBounds);

	
	//! 2)sum up
	for(short species=0; species<num_Charged_Species; species++)
	{
	
		//! set pointer to respective species
		rho_species = rho_species_array +species*num_nodes_in_block;
	
	
	
		for(i = 0; i < num_nodes_in_block; i++)
		{
	
			//! artificially setting an inner density results in a 
			//! slightly asymetric density shape at obstacle surface.
			//! DO NOT SET INNER DENSITY HERE BUT AFTER J2U CONVERSION,
			//! OTHERWISE VELOCITIES ON OBSTALCE NEIGHBOUR NODES ARE 
			//! CONVERTED INCORRECT. 
			//! THIS ESPECIALLY EFFECTS VTH2 PROFILE	
// 			if(Flag[i])
// 			rho_species[i] = obstacle_rho[species];

			//! Convert number density to charge density
			rho_species[i] *= q_of[species];

			//! add to total rho
			rho[i]  +=  rho_species[i];
			
#if defined(use_dust_species_as_field)
			
			//! make sure that electron density does not become smaller than min charge density (or even negative)
			if(rho[i]<0 && Ion_Charges[species]<0)
			{
				rho_species[i] = q_of[species]*(rho[i] - rho_species[i] - MCD_BField);						
						
				rho[i] = MCD_BField;
			}
#endif			
			

			//! build Lam = q*q/m
			Lam[i]  += q_m_ratio[species] * rho_species[i];
		
			//! old:
// 			rho[i]  +=   q_of[species] * rho_species[i];
// 			Lam[i]  += q2m_of[species] * rho_species[i];
		}
	
	}

	//! NOTE:
	//! 1)
	//! No boundary procedure is applied to rho total, as boundary procedures 
	//! are already applied to every individual ion species.
	//! 2)
	//! No boundaries values are provided for Lam & Gam,
	//! since CAM is not performed on box boundaries


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


}

//!----------------------------------------------------
//! convert_Ji2Ui
//!----------------------------------------------------
void CBlock::convert_vth2rho_to_vth2(INT32 species,
				     INT32 rho_type,
				     INT32 vth2rho_type)

{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	D_REAL *rho, *vth2rho_X, *vth2rho_Y, *vth2rho_Z, numDensRho;
	
	rho = Field_Type[rho_type] +species*num_nodes_in_block;

	vth2rho_X = Field_Type[vth2rho_type] +0*num_nodes_in_block;
	vth2rho_Y = Field_Type[vth2rho_type] +1*num_nodes_in_block;
	vth2rho_Z = Field_Type[vth2rho_type] +2*num_nodes_in_block;
	
	
	D_REAL rez_rho_value = 0.;
	
	for(INT32 node=0; node < num_nodes_in_block; node++)
	{
	
		if(Flag[node]) //! Obstacle Cell 
		{

			vth2rho_X[node] = 0.;
			vth2rho_Y[node] = 0.;
			vth2rho_Z[node] = 0.;

		}
		else
		{
		
			//! total charge density is stored in rho, but we need number density:
			numDensRho = rho[node]/q_of[species];
		
			//! get rez_rho_value
			if(fabs(numDensRho) < MCD_J2U)
			rez_rho_value = 1./MCD_J2U;
			else
			rez_rho_value = 1./numDensRho;
		
			//! transform vth2rho to vth2
			vth2rho_X[node] *= rez_rho_value;
			vth2rho_Y[node] *= rez_rho_value;
			vth2rho_Z[node] *= rez_rho_value;
		
		}
	
	}
	
	//! since vth2 can not be collected on ghost cells, just copy them
	bool cp_allBounds[6] = {0,0,0,0,0,0};
	cp_Boundaries(vth2rho_type, cp_allBounds);


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}

//!----------------------------------------------------
//! convert_Ji2Ui
//!----------------------------------------------------
void CBlock::convert_Ji2Ui(INT32 species,
			   INT32 rho_type,
                           INT32 JiX_type,
			   INT32 UiX_type)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	D_REAL *rho;
	D_REAL *JiX, *JiY, *JiZ;
	D_REAL *UiX, *UiY, *UiZ;
	D_REAL *GamX, *GamY, *GamZ;
	
	rho = Field_Type[rho_type] +species*num_nodes_in_block;
	
	JiX = Field_Type[JiX_type];
	JiY = JiX +  num_nodes_in_block;
	JiZ = JiY +  num_nodes_in_block;
	
	UiX = Field_Type[UiX_type];
	UiY = UiX + num_nodes_in_block;
	UiZ = UiY + num_nodes_in_block;
	
	//! to set inner and  boundaries of gam
	GamX = Field_Type[id_Gam];
	GamY = GamX + num_nodes_in_block;
	GamZ = GamY + num_nodes_in_block;
	
	//! Check for MCD:
	//! This should ONLY be done in rho average, CAM and EField 
	//! procedure, otw. velocity will get to slow in MCD regions !!!
	
	//! (-> why does TB Code not produce errors when deviding 
	//! (   zero rho_values eg. when visulizing heavy ion densities ??)
	
	//! For comparison it is done like in TB code
	//    if(set_inner_Dens)
	//    for(INT32 i = 0; i < num_nodes_in_block; i++)
	//    if(rho[i] < MCD_BField) rho[i] = MCD_BField;


	//! NOTE:
	//! VISUALIZAION MAY LOOK FALSE and different to TB Code 
	//! for obstacle ions in case MCD_J2U is high (eg. 0.2).
	//! This is due to the fact, that TB Codes does not use
	//! a lower limit for Visualization but uses MCD
	//! for COMPUTATION of Fields.
	//! So all quantities (B,E,rho and U) are the same in both Codes,
	//! but in this Code visulizes the intern U,
	//! in TB Code intern and visulized U are different.



	D_REAL rez_rho_value = 0.;
	
	for(i = 0; i < BlkNds_X;i++)
	 for(j = 0; j < BlkNds_Y;j++)
	  for(k = 0; k < BlkNds_Z;k++)
	  {
	
		i_j_k =       i*BlkNds_Y*BlkNds_Z
			      +j*BlkNds_Z
			      +k;
	
		//! - use Flag is obsolete version
		//! - MCD / Inner Dens is never set to field 
		//!   directly but catched at impportant locations
		if(!Flag[i_j_k]) //! Plasma Cell 
		{
		
			if(rho[i_j_k]<MCD_J2U)
			rez_rho_value = 1./MCD_J2U;
			else
			rez_rho_value = 1./rho[i_j_k];

			UiX[i_j_k] =  JiX[i_j_k]*rez_rho_value;
			UiY[i_j_k] =  JiY[i_j_k]*rez_rho_value;
			UiZ[i_j_k] =  JiZ[i_j_k]*rez_rho_value;




		
		}
		else
		{
		
			//! should be zero anyway (no particle)
			//! Code also does work with out setting 
			//! U /GAM to zero inside obstacle, but 
			//! less "artefacs" around the obstacle 
			//! occur in case it is set 0.
			GamX[i_j_k] = 0.;
			GamY[i_j_k] = 0.;
			GamZ[i_j_k] = 0.;

			UiX[i_j_k] =  0.;
			UiY[i_j_k] =  0.;
			UiZ[i_j_k] =  0.;


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
	
				UiX[i_j_k] = v_tang[0];
				UiY[i_j_k] = v_tang[1];
				UiZ[i_j_k] = v_tang[2];
			}

		
		}
	}

	//! BELOW IS EXCLUSIVELY TREATEMENT OF BOUNDARIES,
	//! so in case this block is no box boundary,
	//! record particle calculation time and return
	if(!is_box_boundary[_ANY_])
	{

		time_finish = clock();
		time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;
		return;
	}


	//! Treatement of boundaries:
	cp_Boundaries(UiX_type, use_U_inflow_bounds);

	/*
	//! Treatement of boundaries:
	if(is_inflow_species[species])
	{

		INT32 cell_indices[3];
		D_REAL r_vec[3] = {0,0,0};
		D_REAL null_vec[3] = {0,0,0};

		//! NOTE:
		//! no boundaries values are provided for Lam & Gam,
		//! as CAM is not performed on box boundaries
	

		//! 1) Boundary values for total velocity at inflow
		if(is_box_boundary[0] && use_U_inflow_bounds[0]) //! set first plane to background value
		 for(i = 0; i < BlkNds_Y * BlkNds_Z; i++)
		 {



			cell_indices[0] =  i  /      (BlkNds_Y*BlkNds_Z);
			cell_indices[1] = (i -cell_indices[0]*(BlkNds_Y*BlkNds_Z))  /     BlkNds_Z;
			cell_indices[2] =  i -cell_indices[0]*(BlkNds_Y*BlkNds_Z) -cell_indices[1]*BlkNds_Z;


			intern2normedCoords(r_vec, null_vec, cell_indices);
			PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};

			D_REAL ext_v;
			estimate_extern_value(ext_v, 0, 1, x);

			UiX[i] =  ext_v;
			UiY[i] =  0.;
			UiZ[i] =  0.;


		 }
	
		//! 2) Boundary values for total velocity at inflow
		if(is_box_boundary[1] && use_U_inflow_bounds[1]) //! set first plane to background value
		 for(i = (BlkNds_X-1) *BlkNds_Y * BlkNds_Z; i < BlkNds_X *BlkNds_Y *BlkNds_Z; i++)
		 {


			cell_indices[0] =  i  /      (BlkNds_Y*BlkNds_Z);
			cell_indices[1] = (i -cell_indices[0]*(BlkNds_Y*BlkNds_Z))  /     BlkNds_Z;
			cell_indices[2] =  i -cell_indices[0]*(BlkNds_Y*BlkNds_Z) -cell_indices[1]*BlkNds_Z;


			intern2normedCoords(r_vec, null_vec, cell_indices);
			PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};

			D_REAL ext_v;
			estimate_extern_value(ext_v, 0, 1, x);

			UiX[i] =  ext_v;
			UiY[i] =  0.;
			UiZ[i] =  0.;
	
	
		 }
	
	
		
		//! 3) Front Y Boundary
		if(is_box_boundary[2] && use_U_inflow_bounds[2])
		 for(i = 0; i < BlkNds_X; i++)
		  for(k = 0; k < BlkNds_Z; k++)
		  {
			
			i_j_k =  i*BlkNds_Y*BlkNds_Z
				+0
				+k;

			cell_indices[0] =  i;
			cell_indices[1] =  0;
			cell_indices[2] =  k;


			intern2normedCoords(r_vec, null_vec, cell_indices);
			PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};

			D_REAL ext_v;
			estimate_extern_value(ext_v, 0, 1, x);

			UiX[i_j_k] =  ext_v;
			UiY[i_j_k] =  0.;
			UiZ[i_j_k] =  0.;	
	
		  }
	
		//! 4) Back Y Boundary
		if(is_box_boundary[3] && use_U_inflow_bounds[3])
		 for(i = 0; i < BlkNds_X; i++)
		  for(k = 0; k < BlkNds_Z; k++)
		  {
			
			i_j_k =  i*BlkNds_Y*BlkNds_Z
				+(BlkNds_Y-1)*BlkNds_Z
				+k;
		

			cell_indices[0] =  i;
			cell_indices[1] =  (BlkNds_Y-1);
			cell_indices[2] =  k;


			intern2normedCoords(r_vec, null_vec, cell_indices);
			PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};

			D_REAL ext_v;
			estimate_extern_value(ext_v, 0, 1, x);

			UiX[i_j_k] =  ext_v;
			UiY[i_j_k] =  0.;
			UiZ[i_j_k] =  0.;
	
		  }
		
		//! 5) Botton Z Boundary
		if(is_box_boundary[4] && use_U_inflow_bounds[4])
		 for(i = 0; i < BlkNds_X; i++)
		  for(j = 0; j < BlkNds_Y; j++)
		  {
			
			i_j_k =  i*BlkNds_Y*BlkNds_Z
				+j*BlkNds_Z
				+0;

		
			cell_indices[0] =  i;
			cell_indices[1] =  j;
			cell_indices[2] =  0;


			intern2normedCoords(r_vec, null_vec, cell_indices);
			PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};

			D_REAL ext_v;
			estimate_extern_value(ext_v, 0, 1, x);

			UiX[i_j_k] =  ext_v;
			UiY[i_j_k] =  0.;
			UiZ[i_j_k] =  0.;
		
	
		  }
	
		//! 6) Top Z Boundary
		if(is_box_boundary[5] && use_U_inflow_bounds[5])
		 for(i = 0; i < BlkNds_X; i++)
		  for(j = 0; j < BlkNds_Y; j++)
		  {
			
			i_j_k =  i*BlkNds_Y*BlkNds_Z
				+j*BlkNds_Z
				+BlkNds_Z-1;
		
			cell_indices[0] =  i;
			cell_indices[1] =  j;
			cell_indices[2] =  BlkNds_Z-1;


			intern2normedCoords(r_vec, null_vec, cell_indices);
			PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};

			D_REAL ext_v;
			estimate_extern_value(ext_v, 0, 1, x);

			UiX[i_j_k] =  ext_v;
			UiY[i_j_k] =  0.;
			UiZ[i_j_k] =  0.;
	
		  }

	
  	}//! end if is_inflow_species
*/
        
	//! Treatement of boundaries:

	if(is_inflow_species[species])
	{

                INT32 cell_indices[3];
                D_REAL r_vec[3] = {0,0,0};
                D_REAL null_vec[3] = {0,0,0};
            
		//! NOTE:
		//! no boundaries values are provided for Lam & Gam,
		//! as CAM is not performed on box boundaries
	
	
		//! 1) Boundary values for total velocity at inflow
                //! homogeneous
		if(is_box_boundary[0] && use_U_inflow_bounds[0] && use_hom_B_bounds[0]) //! set first plane to background value
		 for(i = 0; i < BlkNds_Y * BlkNds_Z; i++)
		 {

			UiX[i] = V_sw[0];
			UiY[i] = V_sw[1];
			UiZ[i] = V_sw[2];

		 }
		 
		 //! 1) Boundary values for total velocity at inflow
		 //! inhomogeneous
                if(is_box_boundary[0] && use_U_inflow_bounds[0] && !use_hom_B_bounds[0]) //! set first plane to background value
                for(i = 0; i < BlkNds_Y * BlkNds_Z; i++)
                {
                    cell_indices[0] =  i  /      (BlkNds_Y*BlkNds_Z);
                    cell_indices[1] = (i -cell_indices[0]*(BlkNds_Y*BlkNds_Z))  /     BlkNds_Z;
                    cell_indices[2] =  i -cell_indices[0]*(BlkNds_Y*BlkNds_Z) -cell_indices[1]*BlkNds_Z;
                    
                    
                    intern2normedCoords(r_vec, null_vec, cell_indices);
                    PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};
                    
                    PARTICLE_REAL vel[3];
                    set_inflow_velocity(vel,x,species);
                    
                    UiX[i] =  vel[0];
                    UiY[i] =  vel[1];
                    UiZ[i] =  vel[2];
                    
                }
                
		//! 2) Boundary values for total velocity at inflow
		//! homogeneous
		if(is_box_boundary[1] && use_U_inflow_bounds[1] && use_hom_B_bounds[1]) //! set first plane to background value
		 for(i = (BlkNds_X-1) *BlkNds_Y * BlkNds_Z; i < BlkNds_X *BlkNds_Y *BlkNds_Z; i++)
		 {
			UiX[i] = V_sw[0];
			UiY[i] = V_sw[1];
			UiZ[i] = V_sw[2];
	

		 }
                    
                //! 2) Boundary values for total velocity at inflow
               if(is_box_boundary[1] && use_U_inflow_bounds[1] && !use_hom_B_bounds[1]) //! set first plane to background value
                for(i = (BlkNds_X-1) *BlkNds_Y * BlkNds_Z; i < BlkNds_X *BlkNds_Y *BlkNds_Z; i++)
                {
                    cell_indices[0] =  i  /      (BlkNds_Y*BlkNds_Z);
                    cell_indices[1] = (i -cell_indices[0]*(BlkNds_Y*BlkNds_Z))  /     BlkNds_Z;
                    cell_indices[2] =  i -cell_indices[0]*(BlkNds_Y*BlkNds_Z) -cell_indices[1]*BlkNds_Z;
                    
                    
                    intern2normedCoords(r_vec, null_vec, cell_indices);
                    PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};
                    
                    PARTICLE_REAL vel[3];
                    set_inflow_velocity(vel,x,species);
                    
                    UiX[i] =  vel[0];
                    UiY[i] =  vel[1];
                    UiZ[i] =  vel[2];
                    
                    
                }
	
		
		//! 3) Front Y Boundary
		//! homogeneous
		if(is_box_boundary[2] && use_U_inflow_bounds[2] && use_hom_B_bounds[2])
		 for(i = 0; i < BlkNds_X; i++)
		  for(k = 0; k < BlkNds_Z; k++)
		  {
			
			i_j_k =  i*BlkNds_Y*BlkNds_Z
				+0
				+k;
	
			UiX[i_j_k] = V_sw[0];
			UiY[i_j_k] = V_sw[1];
			UiZ[i_j_k] = V_sw[2];
	
		  }
		  
                //! 3) Front Y Boundary
                if(is_box_boundary[2] && use_U_inflow_bounds[2] && !use_hom_B_bounds[2])
                    for(i = 0; i < BlkNds_X; i++)
                        for(k = 0; k < BlkNds_Z; k++)
                        {
                            
                            i_j_k =  i*BlkNds_Y*BlkNds_Z
                            +0
                            +k;
                            
                            cell_indices[0] =  i;
                            cell_indices[1] =  0;
                            cell_indices[2] =  k;
                            
                            
                            intern2normedCoords(r_vec, null_vec, cell_indices);
                            PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};
                            
                            PARTICLE_REAL vel[3];
                            set_inflow_velocity(vel,x,species);
                            
                            UiX[i_j_k] =  vel[0];
                            UiY[i_j_k] =  vel[1];
                            UiZ[i_j_k] =  vel[2];
                            
                        }
                          
            //! 4) Back Y Boundary
            //! homogeneous
            if(is_box_boundary[3] && use_U_inflow_bounds[3] && use_hom_B_bounds[3])
                for(i = 0; i < BlkNds_X; i++)
                for(k = 0; k < BlkNds_Z; k++)
                {
                    
                    i_j_k =  i*BlkNds_Y*BlkNds_Z
                            +(BlkNds_Y-1)*BlkNds_Z
                            +k;
            
                    UiX[i_j_k] = V_sw[0];
                    UiY[i_j_k] = V_sw[1];
                    UiZ[i_j_k] = V_sw[2];
    
                }
            
            //! 4) Back Y Boundary
            if(is_box_boundary[3] && use_U_inflow_bounds[3] && !use_hom_B_bounds[3])
                for(i = 0; i < BlkNds_X; i++)
                    for(k = 0; k < BlkNds_Z; k++)
                    {
                        
                        i_j_k =  i*BlkNds_Y*BlkNds_Z
                        +(BlkNds_Y-1)*BlkNds_Z
                        +k;
                        
                        
                        cell_indices[0] =  i;
                        cell_indices[1] =  (BlkNds_Y-1);
                        cell_indices[2] =  k;
                        
                        
                        intern2normedCoords(r_vec, null_vec, cell_indices);
                        PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};
                        
                        PARTICLE_REAL vel[3];
                        set_inflow_velocity(vel,x,species);
                        
                        UiX[i_j_k] =  vel[0];
                        UiY[i_j_k] =  vel[1];
                        UiZ[i_j_k] =  vel[2];
                    }
                    
            //! 5) Botton Z Boundary
            //! homogeneous
            if(is_box_boundary[4] && use_U_inflow_bounds[4] && use_hom_B_bounds[4])
                for(i = 0; i < BlkNds_X; i++)
                for(j = 0; j < BlkNds_Y; j++)
                {
                    
                    i_j_k =  i*BlkNds_Y*BlkNds_Z
                            +j*BlkNds_Z
                            +0;
            
                    UiX[i_j_k] = V_sw[0];
                    UiY[i_j_k] = V_sw[1];
                    UiZ[i_j_k] = V_sw[2];
    
                }

                //! 5) Botton Z Boundary
                if(is_box_boundary[4] && use_U_inflow_bounds[4] && !use_hom_B_bounds[4])
                    for(i = 0; i < BlkNds_X; i++)
                        for(j = 0; j < BlkNds_Y; j++)
                        {
                            
                            i_j_k =  i*BlkNds_Y*BlkNds_Z
                            +j*BlkNds_Z
                            +0;
                            
                            
                            cell_indices[0] =  i;
                            cell_indices[1] =  j;
                            cell_indices[2] =  0;
                            
                            
                            intern2normedCoords(r_vec, null_vec, cell_indices);
                            PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};
                            
                            PARTICLE_REAL vel[3];
                            set_inflow_velocity(vel,x,species);
                            
                            UiX[i_j_k] =  vel[0];
                            UiY[i_j_k] =  vel[1];
                            UiZ[i_j_k] =  vel[2];
                            
                        }
                        
            //! 6) Top Z Boundary
            //! homogeneous
            if(is_box_boundary[5] && use_U_inflow_bounds[5] && use_hom_B_bounds[5])
                for(i = 0; i < BlkNds_X; i++)
                for(j = 0; j < BlkNds_Y; j++)
                {
                    
                    i_j_k =  i*BlkNds_Y*BlkNds_Z
                            +j*BlkNds_Z
                            +BlkNds_Z-1;
            
                    UiX[i_j_k] = V_sw[0];
                    UiY[i_j_k] = V_sw[1];
                    UiZ[i_j_k] = V_sw[2];
    
                }
            //! 6) Top Z Boundary
            if(is_box_boundary[5] && use_U_inflow_bounds[5] && !use_hom_B_bounds[5])
                for(i = 0; i < BlkNds_X; i++)
                    for(j = 0; j < BlkNds_Y; j++)
                    {
                        
                        i_j_k =  i*BlkNds_Y*BlkNds_Z
                        +j*BlkNds_Z
                        +BlkNds_Z-1;
                        
                        cell_indices[0] =  i;
                        cell_indices[1] =  j;
                        cell_indices[2] =  BlkNds_Z-1;
                        
                        
                        intern2normedCoords(r_vec, null_vec, cell_indices);
                        PARTICLE_REAL x[3] = {r_vec[0], r_vec[1], r_vec[2]};
                        
                        PARTICLE_REAL vel[3];
                        set_inflow_velocity(vel,x,species);
                        
                        UiX[i_j_k] =  vel[0];
                        UiY[i_j_k] =  vel[1];
                        UiZ[i_j_k] =  vel[2];
                        
                    }
	
            }//! end if is_inflow_species

	//! record particle calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

	
}

//!-----------------------------------------------------------
//! average_Ui_rho_setREZrho
//!-----------------------------------------------------------
void CBlock::average_Ui_rho_setREZrho(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

	D_REAL *rho_n, *rho_np1;

	if(mesh_type==STAGGERED)
	{
		rho_n   = Field_Type[  id_rho_n_HIMesh];
		rho_np1 = Field_Type[id_rho_np1_HIMesh];
	}
	else
	{

		rho_n   = Field_Type[  id_rho_n];
		rho_np1 = Field_Type[id_rho_np1];
	}


	D_REAL* rez_rho = Field_Type[id_rho_rez];
	
	D_REAL* Ui_plus    = Field_Type[id_UI_plus];
	D_REAL* Ui_minus   = Field_Type[id_UI_minus];
	
	//! use id_UI_minus for result of average,
	//! as it not used after this any more
	//! (Ui_plus in contrast is used in CAM)
	for(i = 0; i < 3 *num_nodes_in_block; i++)
	Ui_minus[i] = 0.5 *(Ui_plus[i] + Ui_minus[i]);
	
	
	D_REAL average_rho = 0.;
	
	for(i = 0; i < num_nodes_in_block; i++)
	{


	
		//! Set rez_rho in Obstacle to zero in order to eliminate 
		//! hall term
		if(Flag[i])
		rez_rho[i] = 0.;
		else
		{
			average_rho = 0.5*(rho_np1[i] +rho_n[i]);
	
			if(average_rho<MCD_BField)
			rez_rho[i] = 1./MCD_BField;
			else
			rez_rho[i] = 1./average_rho;
	
		
		}
	
	
	
	}

	//! copy un-shifted fields
	rho_n   = Field_Type[  id_rho_n];
	rho_np1 = Field_Type[id_rho_np1];
	
	memcpy(rho_n, rho_np1, num_nodes_in_block *sizeof(D_REAL));

	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}

//!-------------------------------------------------------
//! gather_Oct
//!-------------------------------------------------------
void CBlock::gather_Oct(short species,
			bool collect_rho,
			bool collect_Ji,
			bool collect_vth2,
			bool skip_unmarked_particle,
			INT32 id_oct,
			D_REAL* rho,
			D_REAL* JiX)
{



	//! double precision is used for two reasons:
	//! - it is more precise
	//! - it is faster to typecast the pos and velocities into
	//!   double once and do the remaining calculations in
	//!   double precision
	D_REAL factor;
	D_REAL r[3],v[3], v2;
	particle* active_particle;
	
	
	D_REAL *JiY = JiX +num_nodes_in_block;
	D_REAL *JiZ = JiY +num_nodes_in_block;
	
	//! using a temp array for collecting is ~10% faster
	D_REAL rho_JX_JY_JZ_temp[32];
	
	
	D_REAL shape_func[8];
	D_REAL rez_CellVol = (1./CellVol_of_L[RLevel]);


	//! get index of oct (0 or 1 for each dimension)
	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);
	
	
	//! In case Block is Box boundary, collecting shall start at 0, else at 1
	//! As is_box_boundary returns 1 in case Box boundary and 0 else, the 
	//! result has to negated (using ! )
	INT32 start[3] = {!a*!is_box_boundary[0] +a*BlkNds_X/2,
			  !b*!is_box_boundary[2] +b*BlkNds_Y/2,
			  !c*!is_box_boundary[4] +c*BlkNds_Z/2};


	INT32 end[3]={BlkNds_X-1 -!a*((BlkNds_X-2)/2),
		     BlkNds_Y-1 -!b*((BlkNds_Y-2)/2),
		     BlkNds_Z-1 -!c*((BlkNds_Z-2)/2)};

	

   //! If Block is no Box-Boundary the following is valid:
   //! 1) start from 1 as 0 is not accesable by particles (GC)
   //! 2) just gather to Max-1 as no Max+1 does exist 

   //! THE FOLLOWING IS ONLY TRUE IN CASE BLOCK IS BOX BOUNDARY:
   //! As loop ends at BlkNds_X-2, all new inserted particles at BlkNds_X-1
   //! will NOT be gathered. So collected number is always smaller than
   //! total number of particles !!! Also the sgeometry is asymetric with
   //! respect to GN at 0 and BlkNds_X-1:
   //! 1) The particles attached to the GN 0 contribute 
   //!    to the denity of GN 0 and Node 1.
   //! 2) The particles attached to the GN BlkNds_X-1 DO ONLY contribute
   //!    to the denity of GN BlkNds_X-1.
   //! 3) As GN 0 and GN BlkNds_X-1 are overwritten with boundary values,
   //!    gathering procedure may stop at BlkNds_X-2.
   //! 4) Particles at BlkNds_X-1 may enter Box by thermal Motion and are 
   //!    necessary to avoid "plasma to vacuum" effect.

   //! It was tested that all particles are gatherd when switching off the 
   //! filling procedure at BlkNds_X-1 boundaries.
   for(i = start[0]; i < end[0]; i++)
    for(j = start[1]; j < end[1]; j++)
     for(k = start[2]; k < end[2]; k++)
     {

	i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;

      //! -----------------------------------------
      ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
      i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);


      //! ------------------------------------------
      ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
      ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
      i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);

      ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);


      memset(rho_JX_JY_JZ_temp,0,32*sizeof(D_REAL));



      //! standard shema in case rho or Ji is gatherd
      //! (catch this outside the loop for improved performance)
      if(!collect_vth2)
      	for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
#ifdef TRACK_PARTICLE
	if( !skip_unmarked_particle || pArray[species][i_j_k][part_index].number>=0)
#endif
      	{

	    active_particle = pArray[species][i_j_k] +part_index;


	    //! NOTE:
	    //! has to be set without using memcpy
	    //! FLOAT to DOUBLE
	    v[0] =  active_particle->v[0];
	    v[1] =  active_particle->v[1];
	    v[2] =  active_particle->v[2];

	    r[0] =  active_particle->rel_r[0];
	    r[1] =  active_particle->rel_r[1];
	    r[2] =  active_particle->rel_r[2];

	    factor = active_particle->weight * rez_CellVol;

            shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2])*factor;
            shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2])*factor;
            shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2])*factor;
            shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2])*factor;

            shape_func[4] = (   r[0])*(   r[1])*(1.-r[2])*factor;
            shape_func[5] = (   r[0])*(1.-r[1])*(   r[2])*factor;
            shape_func[6] = (1.-r[0])*(   r[1])*(   r[2])*factor;
            shape_func[7] = (   r[0])*(   r[1])*(   r[2])*factor;

	    //! --------- rho --------------------------------------------------
	    //! if requiered, charge/mass will be applied later but not here !!!
	    if(collect_rho)
	     for(INT32 corner=0; corner<8; corner++)
	     rho_JX_JY_JZ_temp[corner] += shape_func[corner];
	    //! ----------------------------------------------------------------

	    //! --------- Ji ---------------------------------------------------
	    //! if requiered, charge/mass will be applied later but not here !!!
	    if(collect_Ji)
	     for(INT32 corner=0; corner<8; corner++)
	     {
			rho_JX_JY_JZ_temp[corner + 8] += shape_func[corner]* v[0];
			rho_JX_JY_JZ_temp[corner +16] += shape_func[corner]* v[1];
			rho_JX_JY_JZ_temp[corner +24] += shape_func[corner]* v[2];
	     }
	    //! ----------------------------------------------------------------
	

	    num_total_particles_collected++;

	}//! end for particle list



	if(collect_vth2)
	 for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
#ifdef TRACK_PARTICLE
	if( !skip_unmarked_particle || pArray[species][i_j_k][part_index].number>=0)
#endif
	 {

	    	 active_particle = pArray[species][i_j_k] +part_index;


		//! NOTE:
		//! has to be  copied without using memcpy
		//! FLOAT to DOUBLE
		v[0] =  active_particle->v[0];
		v[1] =  active_particle->v[1];
		v[2] =  active_particle->v[2];
	
		r[0] =  active_particle->rel_r[0];
		r[1] =  active_particle->rel_r[1];
		r[2] =  active_particle->rel_r[2];

		factor = active_particle->weight * rez_CellVol;

		shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2])*factor;
		shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2])*factor;
		shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2])*factor;
		shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2])*factor;
	
		shape_func[4] = (   r[0])*(   r[1])*(1.-r[2])*factor;
		shape_func[5] = (   r[0])*(1.-r[1])*(   r[2])*factor;
		shape_func[6] = (1.-r[0])*(   r[1])*(   r[2])*factor;
		shape_func[7] = (   r[0])*(   r[1])*(   r[2])*factor;


		//! --------- vth --------------------------------------------------
		//! if requiered, charge/mass will be applied later but not here !!!
		 for(INT32 corner=0; corner<8; corner++)
		 {
			rho_JX_JY_JZ_temp[corner + 8] += shape_func[corner]* v[0]*v[0];
			rho_JX_JY_JZ_temp[corner +16] += shape_func[corner]* v[1]*v[1];
			rho_JX_JY_JZ_temp[corner +24] += shape_func[corner]* v[2]*v[2];
		 }
		//! ----------------------------------------------------------------

		//! Collect vp² to rho:
// 		v2 = vec_len2(v);

		//! distribute v² to cell corners
// 		for(INT32 corner=0; corner<8; corner++)
// 		rho_JX_JY_JZ_temp[corner] += shape_func[corner] *v2;

		num_total_particles_collected++;

      	}	//! end while particle->next


	if(collect_rho)
	{
		rho[i_j_k]   += rho_JX_JY_JZ_temp[0];
		rho[ip1_j_k] += rho_JX_JY_JZ_temp[1];
		rho[i_jp1_k] += rho_JX_JY_JZ_temp[2];
		rho[i_j_kp1] += rho_JX_JY_JZ_temp[3];

		rho[ip1_jp1_k]   += rho_JX_JY_JZ_temp[4];
		rho[ip1_j_kp1]   += rho_JX_JY_JZ_temp[5];
		rho[i_jp1_kp1]   += rho_JX_JY_JZ_temp[6];
		rho[ip1_jp1_kp1] += rho_JX_JY_JZ_temp[7];
	}



	//! --------- Ji ---------------------------
	//! gather Ji or roh_v2mean on Ji, depended on 
	//! which quantity is stored in rho_JX_JY_JZ_temp (see above).
	if(collect_Ji || collect_vth2)
	{

		//! --------- JiX ---------------------------
		JiX[i_j_k  ] += rho_JX_JY_JZ_temp[0 +8];
		JiX[ip1_j_k] += rho_JX_JY_JZ_temp[1 +8];
		JiX[i_jp1_k] += rho_JX_JY_JZ_temp[2 +8];
		JiX[i_j_kp1] += rho_JX_JY_JZ_temp[3 +8];

		JiX[ip1_jp1_k  ] += rho_JX_JY_JZ_temp[4 +8];
		JiX[ip1_j_kp1  ] += rho_JX_JY_JZ_temp[5 +8];
		JiX[i_jp1_kp1  ] += rho_JX_JY_JZ_temp[6 +8];
		JiX[ip1_jp1_kp1] += rho_JX_JY_JZ_temp[7 +8];

		//! --------- JiY ---------------------------
		JiY[i_j_k  ] += rho_JX_JY_JZ_temp[0 +16];
		JiY[ip1_j_k] += rho_JX_JY_JZ_temp[1 +16];
		JiY[i_jp1_k] += rho_JX_JY_JZ_temp[2 +16];
		JiY[i_j_kp1] += rho_JX_JY_JZ_temp[3 +16];

		JiY[ip1_jp1_k  ] += rho_JX_JY_JZ_temp[4 +16];
		JiY[ip1_j_kp1  ] += rho_JX_JY_JZ_temp[5 +16];
		JiY[i_jp1_kp1  ] += rho_JX_JY_JZ_temp[6 +16];
		JiY[ip1_jp1_kp1] += rho_JX_JY_JZ_temp[7 +16];

		//! --------- JiZ ---------------------------
		JiZ[i_j_k  ] += rho_JX_JY_JZ_temp[0 +24];
		JiZ[ip1_j_k] += rho_JX_JY_JZ_temp[1 +24];
		JiZ[i_jp1_k] += rho_JX_JY_JZ_temp[2 +24];
		JiZ[i_j_kp1] += rho_JX_JY_JZ_temp[3 +24];

		JiZ[ip1_jp1_k  ] += rho_JX_JY_JZ_temp[4 +24];
		JiZ[ip1_j_kp1  ] += rho_JX_JY_JZ_temp[5 +24];
		JiZ[i_jp1_kp1  ] += rho_JX_JY_JZ_temp[6 +24];
		JiZ[ip1_jp1_kp1] += rho_JX_JY_JZ_temp[7 +24];
	}
	



    } //! end for Block
}

//!-------------------------------------------------------
//! gather_Oct
//!-------------------------------------------------------
void CBlock::gather_from_parent_Block(short species,
			 	      bool collect_rho,
			  	      bool collect_Ji,
			  	      bool collect_vth2,
			  	      D_REAL* rho,
			  	      D_REAL* JiX)
{


	D_REAL* JiY = JiX + 1 * num_nodes_in_block;
	D_REAL* JiZ = JiX + 2 * num_nodes_in_block;


	particle* active_particle;
	
	//! double precision is used for to reasons:
	//! - it is more precise
	//! - it is faster to typecast the pos and velocities into
	//!   double once and do the remaining calculations in
	//!   double precision
	D_REAL factor = 1.;
	D_REAL r[3], v_part[3], shape_func[8], v2;
	WEIGHT_REAL rez_CellVol = (1./CellVol_of_L[RLevel]);
	
	
	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = Blk_Index[2] * (BlkNds_Z/2-1);


	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block
	for (INT32 a = 0; a < BlkNds_X/2-1; a++)
	 for (INT32 b = 0; b < BlkNds_Y/2-1; b++)
	  for (INT32 c = 0; c < BlkNds_Z/2-1; c++)
	  {
	
	
		//! i,j,k are indices of parent cell
		i = parent_Xstart +a+1;
		j = parent_Ystart +b+1;
		k = parent_Zstart +c+1;
		
		i_j_k =  i*BlkNds_Y*BlkNds_Z
			+j*BlkNds_Z
			+k;
	
	
		//       active_particle = parent->pArray[species][i_j_k];
		for(INT32 part_index=0; part_index<parent->num_MPiC[species][i_j_k]; part_index++)
		{
		
			active_particle = parent->pArray[species][i_j_k] +part_index;
		
		
			//! NOTE:
			//! has to be  copied without using memcpy
			//! FLOAT to DOUBLE
			v_part[0] = active_particle->v[0];
			v_part[1] = active_particle->v[1];
			v_part[2] = active_particle->v[2];
		
		
			//! NOTE:
			//! r has to be  copied since
			//! - float to double
			//! - it is changed during this function !
			r[0] = active_particle->rel_r[0];
			r[1] = active_particle->rel_r[1];
			r[2] = active_particle->rel_r[2];
		
		
			//! transform parent r to this Block r
			r[0] *= 2.;
			r[1] *= 2.;
			r[2] *= 2.;
		
			//! u,v,w are indices of child cell
			u = 2*a+1 +int(r[0]);
			v = 2*b+1 +int(r[1]);
			w = 2*c+1 +int(r[2]);
		
			//! r now is relative position in child cell (u,v,w)
			r[0] -= int(r[0]);
			r[1] -= int(r[1]);
			r[2] -= int(r[2]);
		
			u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
			up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC
			u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC
			up1_vp1_w =  (u+1)*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!LC
		
			u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
			up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC
			u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC
			up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN
		
		
		
			//! calc shape function
			factor = active_particle->weight * rez_CellVol;
		
			shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2])*factor;
			shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2])*factor;
			shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2])*factor;
			shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2])*factor;
		
			shape_func[4] = (   r[0])*(   r[1])*(1.-r[2])*factor;
			shape_func[5] = (   r[0])*(1.-r[1])*(   r[2])*factor;
			shape_func[6] = (1.-r[0])*(   r[1])*(   r[2])*factor;
			shape_func[7] = (   r[0])*(   r[1])*(   r[2])*factor;
		
			//! --------- rho ------------------
			if(collect_rho)
			{
		
				rho[u_v_w] += shape_func[0];
				rho[up1_v_w] += shape_func[1];
				rho[u_vp1_w] += shape_func[2];
				rho[u_v_wp1] += shape_func[3];
		
				rho[up1_vp1_w  ]   += shape_func[4];
				rho[up1_v_wp1  ]   += shape_func[5];
				rho[u_vp1_wp1  ]   += shape_func[6];
				rho[up1_vp1_wp1]   += shape_func[7];
		
			}
		
			if(collect_vth2)
			{
		
				//! Collect vp² to rho:
				v2 = vec_len2(v_part);
		
				rho[u_v_w]   += shape_func[0] *v2;
				rho[up1_v_w] += shape_func[1] *v2;
				rho[u_vp1_w] += shape_func[2] *v2;
				rho[u_v_wp1] += shape_func[3] *v2;
		
				rho[up1_vp1_w  ]   += shape_func[4] *v2;
				rho[up1_v_wp1  ]   += shape_func[5] *v2;
				rho[u_vp1_wp1  ]   += shape_func[6] *v2;
				rho[up1_vp1_wp1]   += shape_func[7] *v2;
			
			}
		
			//! --------- Ji ---------------------------
			if(collect_Ji && !collect_vth2)
			{
		
				//! --------- JiX ---------------------------
				JiX[u_v_w]   += shape_func[0] * v_part[0];
				JiX[up1_v_w] += shape_func[1] * v_part[0];
				JiX[u_vp1_w] += shape_func[2] * v_part[0];
				JiX[u_v_wp1] += shape_func[3] * v_part[0];
		
				JiX[up1_vp1_w]   += shape_func[4] * v_part[0];
				JiX[up1_v_wp1]   += shape_func[5] * v_part[0];
				JiX[u_vp1_wp1]   += shape_func[6] * v_part[0];
				JiX[up1_vp1_wp1] += shape_func[7] * v_part[0];
		
				//! --------- JiY ---------------------------
				JiY[u_v_w]   += shape_func[0] * v_part[1];
				JiY[up1_v_w] += shape_func[1] * v_part[1];
				JiY[u_vp1_w] += shape_func[2] * v_part[1];
				JiY[u_v_wp1] += shape_func[3] * v_part[1];
		
				JiY[up1_vp1_w]   += shape_func[4] * v_part[1];
				JiY[up1_v_wp1]   += shape_func[5] * v_part[1];
				JiY[u_vp1_wp1]   += shape_func[6] * v_part[1];
				JiY[up1_vp1_wp1] += shape_func[7] * v_part[1];
		
				//! --------- JiZ ---------------------------
				JiZ[u_v_w]   += shape_func[0] * v_part[2];
				JiZ[up1_v_w] += shape_func[1] * v_part[2];
				JiZ[u_vp1_w] += shape_func[2] * v_part[2];
				JiZ[u_v_wp1] += shape_func[3] * v_part[2];
		
				JiZ[up1_vp1_w]   += shape_func[4] * v_part[2];
				JiZ[up1_v_wp1]   += shape_func[5] * v_part[2];
				JiZ[u_vp1_wp1]   += shape_func[6] * v_part[2];
				JiZ[up1_vp1_wp1] += shape_func[7] * v_part[2];
			}
		
			num_particles_collected_from_parent++;
		
		
		}	//! for particle in cell
	} //! end for Block



}


//!-------------------------------------------------------------//
//! count number of cells with too much negative charges
//! if cell has rho<0, delete heaviest negative particle
//!-------------------------------------------------------------//
void CBlock::negative_particles(INT64 *num_negCells_delPart, INT32 level)
{

   	 D_REAL *rho = Field_Type[id_rho_n];
	
	//! cell indices
	INT32 i_j_k;
	INT32 ind[3];
	
	for(INT32 species=0; species<num_Charged_Species; species++)
	{

	//! check whether species has negative charge
	if(Ion_Charges[species]<0)
	 for(ind[0] = 1; ind[0] < BlkNds_X-1; ind[0]++)
	  for(ind[1] = 1; ind[1] < BlkNds_Y-1; ind[1]++)
	   for(ind[2] = 1; ind[2] < BlkNds_Z-1; ind[2]++)
	   {       
		
		i_j_k   =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];

		//! do not count cells within obstacle
		if(rho[i_j_k]<0 && !Flag[i_j_k])
		{
			
			short rho_neg = int(rho[i_j_k]*-10);
			
			if(num_MPiC[species][i_j_k]>1)
			{
				for(short i=0; i<=rho_neg; i++)
				{	
					//! first should always be heaviest
					INT32 part_index = 0;
					delete_particle(species, i_j_k, part_index);
					num_negCells_delPart[1]++;
					num_negCells_delPart[2+level]++;

					//! update statistics
					num_MPiC[species][i_j_k]--;
				}

			}
			num_negCells_delPart[0]++;

		}
	   }
	}

}