

#include "CHybrid.h"

#include "utils.h"
#include "parameters.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "absolute_Globals.h"
#include "unistd.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <math.h>

//! time to wait in micro-seconds
#define TIME_TO_WAIT 200000


//!--------------------------------------------------------
//!- inject_obstacle_ions:
//!  This function always is called by the main loop.
//!  Call any user defined function from this function
//!  In case no profile shall be used, commend out
//!  everything below.
//!--------------------------------------------------------
void CHybrid::inject_obstacle_ions(void)
{

	if(INJECT_PARTICLE_TO_TRACK && TL==TL_MARK_PARTICLE_TRACKS)
	  for(INT32 species=0; species<num_Particle_Species; species++)
	    if(num_mark_particle_in_species[species]>0)	
	      inject_callisto_particle_to_track(species);
	
	if(TestParticle_Simulation)
	return; 

	//#ifdef use_neutral_species_as_field
	//	insert_ions_from_neutral_profile();
	//#endif
	
// 	insert_sphere_cylinder_ions(1);
// 	inject_Moon_ions(1);
// 	insert_enceladus_ions(0);
	insert_Titan_Ions();
// 	insert_Titan_Ions2(2);insert_Titan_Ions2(3);insert_Titan_Ions2(4);
}


void CHybrid::insert_Titan_Ions2(INT32 species)
{
	clock_t start,finish;
	double time;
	start = clock();
	
	INT64* num_injected_particles_species = new INT64[num_ion_prod_fields*num_Neutral_Species];
        memset(num_injected_particles_species, 0,num_ion_prod_fields*num_Neutral_Species*sizeof(INT64));
	
	
	if(TL==2 || TL==TL_at_last_restore_state+2)
	{	
		INT32 total_particles_to_insert=0;
		for(INT32 species=0; species<num_Particle_Species; species++)
		total_particles_to_insert += obstacle_MP_num_each_TL[species];
		
		insert_every_x_TL= int(insert_every_x_TL*num_newlyionised_particles/total_particles_to_insert + 1);	
	}
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
// 			if(mpi_myRank == temp_Block->responsible_mpi_process)
// // 			  if(!temp_Block->child_array)
// 			    temp_Block->inject_extern_photorate(); 
			    //! If the process is responsible for block, proceed
			    if(mpi_myRank==temp_Block->responsible_mpi_process)
			      for(INT32 oct=0; oct<8; oct++)
			      {
				//! maybe just a oct of block is refined 
				if(!temp_Block->child_array[oct])
				  temp_Block->insert_ions_from_ionprod_file(insert_every_x_TL,oct,num_injected_particles_species);
				
			      }

			
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}
	
	INT64 local_info_values[num_ion_prod_fields*num_Neutral_Species];
	stringstream info_names[num_ion_prod_fields*num_Neutral_Species];
	
	for(INT32 neutralSpec=0; neutralSpec<num_Neutral_Species; neutralSpec++)
	 for(INT32 species=0; species<num_ion_prod_fields; species++)
	 {
		 
		local_info_values[species*num_Neutral_Species + neutralSpec] = num_injected_particles_species[species*num_Neutral_Species + neutralSpec];
		info_names[species*num_Neutral_Species + neutralSpec] << "   ->injected particles of species "<<species<< " from neutral species "<<neutralSpec<<" : ";		
	 }	

	//! Reduce via MPI if required
	show_information(local_info_values,
		info_names,
		num_ion_prod_fields*num_Neutral_Species,
		BUILD_SUM);
	
	//! get number of inserted particles to adjust this number in next time step
	if(TL==1 || TL==TL_at_last_restore_state+1)
	for(INT32 neutralSpec=0; neutralSpec<num_Neutral_Species; neutralSpec++)
	 for(INT32 species=0; species<num_ion_prod_fields; species++)
	  num_newlyionised_particles += local_info_values[species*num_Neutral_Species + neutralSpec];
 
 
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "  inject heavy ion time: " << time << "s." << endl << endl; 
  
  
  
}



void CHybrid::insert_Titan_Ions()
{
	clock_t start,finish;
	double time;
	start = clock();
	
	log_file<<"INJECTING IONS ..."<<endl;
	
	INT64* num_injected_particles_species = new INT64[num_ion_prod_fields];
        memset(num_injected_particles_species, 0,num_ion_prod_fields*sizeof(INT64));
	
	//! adjust number of particles that are injected
	//! do this only at second time step (first is used to get number of inserted particles)
	if(TL==2 || TL==TL_at_last_restore_state+2)
	{	

		INT32 total_particles_to_insert=0;

		for(INT32 species=0; species<num_Particle_Species; species++)
		total_particles_to_insert += obstacle_MP_num_each_TL[species];

		insert_every_x_TL= int(insert_every_x_TL*num_newlyionised_particles/total_particles_to_insert + 1);	

	}


	//! go over every level
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
// 			if(mpi_myRank == temp_Block->responsible_mpi_process)
// // 			  if(!temp_Block->child_array)
// 			    temp_Block->inject_extern_photorate(); 
			    //! If the process is responsible for block, proceed
			    if(mpi_myRank==temp_Block->responsible_mpi_process)
			      for(INT32 oct=0; oct<8; oct++)
			      {
				//! maybe just a oct of block is refined 
				if(!temp_Block->child_array[oct])
				  temp_Block->insert_ions_from_ionprod_file(insert_every_x_TL,oct,num_injected_particles_species);
				
			      }

			
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}

	//! Provide Information about how many partilce were injected in total
	INT64 local_info_values[num_ion_prod_fields];
	stringstream info_names[num_ion_prod_fields];
	
// 	for(INT32 neutralSpec=0; neutralSpec<num_ion_prod_fields; neutralSpec++)
	 for(INT32 species=0; species<num_ion_prod_fields; species++)
	 {
		 
		local_info_values[species] = num_injected_particles_species[species];
		info_names[species] << "   ->injected particles of species "<<species+num_Inflow_Species<< " from neutral species "<<species<<" : ";		
	 }	

	//! Reduce via MPI if required
	show_information(local_info_values,
		info_names,
		num_ion_prod_fields,
		BUILD_SUM);
	
	//! get number of inserted particles to adjust this number in next time step
	if(TL==1 || TL==TL_at_last_restore_state+1)
// 	for(INT32 neutralSpec=0; neutralSpec<num_Neutral_Species; neutralSpec++)
	 for(INT32 species=0; species<num_ion_prod_fields; species++)
	  num_newlyionised_particles += local_info_values[species];
 
 
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "  inject heavy ion time: " << time << "s." << endl << endl;

}




//!--------------------------------------------------------
//!- inject ions for simple geometries:
//!--------------------------------------------------------
void CHybrid::insert_sphere_cylinder_ions(INT32 species)
{
	
	log_file << " INSERTING PARTICLES IN TEST VOLUME ..." << endl;

	if(num_Particle_Species<=species)
	{
	    log_file << "  WARNING: " << endl
		     << "  species:           " << species << endl
		     << "  num_Particle_Species:   " << num_Particle_Species <<endl
		     << "  has to be at least " <<species+1<< endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	}

	if(obstacle_MP_num_each_TL[species]==0 || obstacle_MP_weight_each_t0[species]==0.)
	{
	    log_file << "  WARNING: " << endl
		     << "  obstacle_MP_num_each_TL["<<species<<"]   : " << obstacle_MP_num_each_TL[species] << endl
		     << "  obstacle_MP_weight_each_t0["<<species<<"]: " << obstacle_MP_weight_each_t0[species] << endl
		     << "  No Weight will be injected." << endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	}
	
	clock_t start,finish;
	double time;
	start = clock();

	//! -----------------------------------------------------
	//! Alloc Particle memory
	PARTICLE_REAL *posX, *posY, *posZ, *v_ionX, *v_ionY, *v_ionZ;

	posX = new PARTICLE_REAL[3*obstacle_MP_num_each_TL[species]];
	posY = posX +obstacle_MP_num_each_TL[species];
	posZ = posY +obstacle_MP_num_each_TL[species];
	memset(posX,0,3*obstacle_MP_num_each_TL[species]*sizeof(PARTICLE_REAL));

	v_ionX = new PARTICLE_REAL[3*obstacle_MP_num_each_TL[species]];
	v_ionY = v_ionX +obstacle_MP_num_each_TL[species];
	v_ionZ = v_ionY +obstacle_MP_num_each_TL[species];
	memset(v_ionX,0,3*obstacle_MP_num_each_TL[species]*sizeof(PARTICLE_REAL));
	//! -----------------------------------------------------

	D_REAL reject_val;
	
	//! for sphere
// 	D_REAL r, phi, theta;
// 	D_REAL inner_radius = 0;
// 	D_REAL outer_radius = 2.*R_Moon;
	
	//! for cylinder
	D_REAL r, phi, z;
	D_REAL inner_radius = 0;
	D_REAL outer_radius = 2.*R_Moon;
	D_REAL height_cylinder = 4.*R_Moon;	

	//! -----------------------------------------------------------------
	//! - 1) Estimate positions of particle and store them in arrays ----
	//! -----------------------------------------------------------------
	for (INT32 ion=0; ion<obstacle_MP_num_each_TL[species]; ion++)
	{
		//! for non isotropic ion volume
		do
		{

			//! spherical polar coordinates
// 			r     =	(outer_radius-inner_radius)*gsl_rng_uniform(randGen_of_species[species]) +inner_radius;
// 			theta =	0.5*M_PI		*gsl_rng_uniform(randGen_of_species[species]) + 0.5*M_PI ;
// 			phi   =	2.*M_PI		        *gsl_rng_uniform(randGen_of_species[species]);

				
			//! cylindrical coordinates
			r 	= (outer_radius-inner_radius)*gsl_rng_uniform(randGen_of_species[species])+inner_radius;
			phi 	= 2.* M_PI 		*gsl_rng_uniform(randGen_of_species[species]);
			z	= 1.*gsl_rng_uniform(randGen_of_species[species]);

			//! reject_val is given by det J/ det J|_max * f(x)/f_max(x)
 
			//! isotrop cylinder
			reject_val = r;

			//! isotrop sphere
// 			reject_val = r*r/(outer_radius*outer_radius)*sin(theta);
			
			//! sphere: decay with 1/r^2
//  			reject_val = sin(theta);			

		//! continue randomize in case reject value is to small 
		//! -> the smaller reject value,
		//!    the more likely it will be chosen.
		}while( 1.*random()/RAND_MAX > reject_val);
		
		//! spherical polar coordinates
// 		posX[ion] = r*cos(phi)*sin(theta);
// 		posY[ion] = r*sin(phi)*sin(theta);
// 		posZ[ion] = r*(cos(theta);
		
		//! cylindrical coordinates
		posX[ion] = ((outer_radius-inner_radius)*r+inner_radius)*cos(phi);
		posY[ion] = ((outer_radius-inner_radius)*r+inner_radius)*sin(phi);
		posZ[ion] = (2.*z - 1.)*height_cylinder;
		
		v_ionX[ion] = 0;
		v_ionY[ion] = 0;
		v_ionZ[ion] = 0;
		

	}
	
	//! -----------------------------------------------------------------
	//! - 2) Estimate start weight of macro particle --------------------
	//! -----------------------------------------------------------------
	PARTICLE_REAL start_weight_HI = (obstacle_MP_weight_each_t0[species] * dt)/(1.*obstacle_MP_num_each_TL[species]);

	//! -----------------------------------------------------------------
	//! - 3) Distribute the positions to Blocks -------------------------
	//! -    (particles finally get initialzed at respective Block) -----
	//! -----------------------------------------------------------------

	PARTICLE_REAL* positions[3]  = {  posX,   posY,   posZ};
	PARTICLE_REAL* velocities[3] = {v_ionX, v_ionY, v_ionZ};
	send_injected_Ions_to_Blks(species, start_weight_HI, positions, velocities, obstacle_MP_num_each_TL);

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " inject heavy ion time: " << time << "s." << endl << endl;


	delete[]   posX;
	delete[] v_ionX;

}




//!--------------------------------------------------------
//!- inject_Moon_ions to avoid vacuum in wake:
//!  - using initial velocity
//!--------------------------------------------------------
void CHybrid::inject_Moon_ions(INT32 species)
{

	log_file << " INJECTING GHOST IONS ..." << endl;


	if(num_Particle_Species<=species)
	{
	    log_file << "  WARNING: " << endl
		     << "  species:           " << species << endl
		     << "  num_Particle_Species:   " << num_Particle_Species <<endl
		     << "  has to be at least " <<species+1<< endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	}

	if(obstacle_MP_num_each_TL[species]==0 || obstacle_MP_weight_each_t0[species]==0.)
	{
	    log_file << "  WARNING: " << endl
		     << "  obstacle_MP_num_each_TL["<<species<<"]   : " << obstacle_MP_num_each_TL[species] << endl
		     << "  obstacle_MP_weight_each_t0["<<species<<"]: " << obstacle_MP_weight_each_t0[species] << endl
		     << "  No Weight will be injected." << endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	}

	clock_t start,finish;
	double time;
	start = clock();

	//! -----------------------------------------------------

	D_REAL r, z, phi, reject_val;

	PARTICLE_REAL *posX, *posY, *posZ, *v_ionX, *v_ionY, *v_ionZ;

	posX = new PARTICLE_REAL[3*obstacle_MP_num_each_TL[species]];
	posY = posX +obstacle_MP_num_each_TL[species];
	posZ = posY +obstacle_MP_num_each_TL[species];
	memset(posX,0,3*obstacle_MP_num_each_TL[species]*sizeof(PARTICLE_REAL));

	v_ionX = new PARTICLE_REAL[3*obstacle_MP_num_each_TL[species]];
	v_ionY = v_ionX +obstacle_MP_num_each_TL[species];
	v_ionZ = v_ionY +obstacle_MP_num_each_TL[species];
	memset(v_ionX,0,3*obstacle_MP_num_each_TL[species]*sizeof(PARTICLE_REAL));


	//! estimate maximal possible value for reject value in
	D_REAL Max_reject_val = R_Obstacle;

	//! -----------------------------------------------------------------
	//! - 1) Estimate positions of particle and store them in arrays ----
	//! -----------------------------------------------------------------
	for (INT32 ion=0; ion<obstacle_MP_num_each_TL[species]; ion++)
	{

		do
		{

			phi   	= 2.*M_PI 		   *gsl_rng_uniform(randGen_of_species[species]);
			z	= 2.*R_Obstacle		   *gsl_rng_uniform(randGen_of_species[species]);
			r	= R_Obstacle		   *gsl_rng_uniform(randGen_of_species[species]);
	
			//! uniform cylinder
			reject_val = r/Max_reject_val;
	
		//! continue randomize in case reject value is to small 
		//! -> the smaller reject value,
		//!    the more likely it will be chosen.
		}while( gsl_rng_uniform(randGen_of_species[species]) > reject_val);


		posX[ion] = z + 0.5 *R_Obstacle;
		posY[ion] = r*cos(phi);
		posZ[ion] = r*sin(phi);

		v_ionX[ion] = V_sw[0];
		v_ionY[ion] = V_sw[1];
		v_ionZ[ion] = V_sw[2];


	}

	//! -----------------------------------------------------------------
	//! - 2) Estimate start weight of macro particle --------------------
	//! -----------------------------------------------------------------
	PARTICLE_REAL start_weight_HI = (obstacle_MP_weight_each_t0[species] * dt)/(1.*obstacle_MP_num_each_TL[species]);


	//! -----------------------------------------------------------------
	//! - 3) Distribute the positions to Blocks -------------------------
	//! -    (particles finally get initialzed at respective Block) -----
	//! -----------------------------------------------------------------

	PARTICLE_REAL* positions[3]  = {  posX,   posY,   posZ};
	PARTICLE_REAL* velocities[3] = {v_ionX, v_ionY, v_ionZ};
	send_injected_Ions_to_Blks(species, start_weight_HI, positions, velocities, obstacle_MP_num_each_TL);

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " inject heavy ion time: " << time << "s." << endl << endl;


	delete[]   posX;
	delete[] v_ionX;

}


//!--------------------------------------------------------
//!- inject_Enceladus_plume_ions:
//!--------------------------------------------------------
void CHybrid::insert_enceladus_ions(INT32 species)
{

	
	log_file << " INSERTING ENCELADUS PLUME IONS/DUST ..." << endl;

	if(num_Particle_Species<=species)
	{
		log_file << "  WARNING: " << endl
		     << "  species:           " << species << endl
		     << "  num_Particle_Species:   " << num_Particle_Species <<endl
		     << "  has to be at least " <<species+1<< endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	}

	if(obstacle_MP_num_each_TL[species]==0 || obstacle_MP_weight_each_t0[species]==0.)
	{
		log_file << "  WARNING: " << endl
		     << "  obstacle_MP_num_each_TL["<<species<<"]   : " << obstacle_MP_num_each_TL[species] << endl
		     << "  obstacle_MP_weight_each_t0["<<species<<"]: " << obstacle_MP_weight_each_t0[species] << endl
		     << "  No Weight will be injected." << endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	}

	clock_t start,finish;
	double time;
	start = clock();

	//! -----------------------------------------------------

	const D_REAL H_theta_gas = 1./((opening_angle_gas*M_PI/180.)*(opening_angle_gas*M_PI/180.));
	const D_REAL H_theta_dust = 1./((opening_angle_dust*M_PI/180.)*(opening_angle_dust*M_PI/180.));
	const D_REAL cosphi0 = cos(phi0);
	const D_REAL costheta0 = cos(theta0);
	const D_REAL costheta0_ChEx = cos(M_PI - theta0);
	const D_REAL sinphi0 = sin(phi0);
	const D_REAL sintheta0 = sin(theta0);
	
	D_REAL r, theta, phi, reject_val;

	D_REAL x[3], rho[1];
	

	PARTICLE_REAL *posX, *posY, *posZ, *v_ionX, *v_ionY, *v_ionZ;

	posX = new PARTICLE_REAL[3*obstacle_MP_num_each_TL[species]];
	posY = posX +obstacle_MP_num_each_TL[species];
	posZ = posY +obstacle_MP_num_each_TL[species];
	memset(posX,0,3*obstacle_MP_num_each_TL[species]*sizeof(PARTICLE_REAL));

	v_ionX = new PARTICLE_REAL[3*obstacle_MP_num_each_TL[species]];
	v_ionY = v_ionX +obstacle_MP_num_each_TL[species];
	v_ionZ = v_ionY +obstacle_MP_num_each_TL[species];
	memset(v_ionX,0,3*obstacle_MP_num_each_TL[species]*sizeof(PARTICLE_REAL));

	//! -----------------------------------------------------------------
	//! - 1) Estimate positions of particle and store them in arrays ----
	//! -----------------------------------------------------------------
	for (INT32 ion=0; ion<obstacle_MP_num_each_TL[species]; ion++)
	{
		do
		{

			//! TODO:
			//! USE OWN GLOBAL RANDOM GENERATOR DIFFERENT FROM SW
			//! RANDOM GENERATOR IN ORDER TO PRODUCE HEAVY IONS
			//! -> ANY PROCESS PRODUCES EQUAL IONS
			//! -> NO NEED FOR COMMUNICATION

			r     =	10*R_Moon	*gsl_rng_uniform(randGen_of_species[species]) + R_Moon;
			theta =	0.5*M_PI		*gsl_rng_uniform(randGen_of_species[species]) + 0.5*M_PI ;
			phi   =	2.*M_PI		        *gsl_rng_uniform(randGen_of_species[species]);
	

			//! plume profile by J. Saur
			if(Ion_Charges[species]>0)
			reject_val=sin(theta)*exp(-(theta-M_PI)*(theta-M_PI)*H_theta_gas)*exp(-(r-R_Moon)/H_d);
			
			if(Ion_Charges[species]<0)
			reject_val=sin(theta)*r/R_Moon*exp(-(theta-M_PI)*(theta-M_PI)*H_theta_dust)*exp(-(r-R_Moon)/H_d_dust);
 

	
		//! continue randomize in case reject value is to small 
		//! -> the smaller reject value,
		//!    the more likely it will be chosen.
		}while( gsl_rng_uniform(randGen_of_species[species]) > reject_val);


		
		posX[ion]=r*(cos(phi)*cosphi0*costheta0*sin(theta)-sin(phi)*sinphi0*sin(theta)-cosphi0*cos(theta)*sintheta0);
		posY[ion] =2.*r*(cosphi0*sin(phi)*sin(theta)+sinphi0*(cos(phi)*costheta0*sin(theta)-cos(theta)*sintheta0)) ;
		posZ[ion] =r*(cos(theta)*costheta0+cos(phi)*sin(theta)*sintheta0);
		
// 		posX[ion] = 2.*RE*gsl_rng_uniform(randGen_of_species[species]);
// 		posY[ion] = 2.*RE*gsl_rng_uniform(randGen_of_species[species]);
// 		posZ[ion] = 2.*RE*gsl_rng_uniform(randGen_of_species[species]) - 5.*RE;
// 		
// 		v_ionX[ion] = 0.3;
	}
	

	//! -----------------------------------------------------------------
	//! - 2) Estimate start weight of macro particle --------------------
	//! -----------------------------------------------------------------
	PARTICLE_REAL start_weight_HI = (1.*obstacle_MP_weight_each_t0[species] * dt)/(1.*obstacle_MP_num_each_TL[species]);

	//! -----------------------------------------------------------------
	//! - 3) Distribute the positions to Blocks -------------------------
	//! -    (particles finally get initialzed at respective Block) -----
	//! -----------------------------------------------------------------

	PARTICLE_REAL* positions[3]  = {  posX,   posY,   posZ};
	PARTICLE_REAL* velocities[3] = {v_ionX, v_ionY, v_ionZ};
	send_injected_Ions_to_Blks(species, start_weight_HI, positions, velocities, obstacle_MP_num_each_TL);

	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " inject heavy ion time: " << time << "s." << endl << endl;


	delete[]   posX;
	delete[] v_ionX;

}




//!--------------------------------------------------------
//!- inject_particle_to_track:
//!--------------------------------------------------------
void CHybrid::inject_particle_to_track(INT32 species)
{

// 	if(TL!=TL_MARK_PARTICLE_TRACKS) return;
	
	log_file << " INJECTING PARTICLE TO TRACK ..." << endl;
	//! define which species shall be inserted


	if(num_Particle_Species<=species)
	{
		log_file << "  WARNING: " << endl
		     << "  species:           " << species << endl
		     << "  num_Particle_Species:   " << num_Particle_Species <<endl
		     << "  has to be at least " <<species+1<< endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	} 


	clock_t start,finish;
	double time;
	start = clock();

	//! -----------------------------------------------------
	D_REAL r, norm_r, theta, phi, reject_val;
	PARTICLE_REAL *posX, *posY, *posZ, *v_ionX, *v_ionY, *v_ionZ;

	posX = new PARTICLE_REAL[3*num_mark_particle_in_species[species]];
	posY = posX +num_mark_particle_in_species[species];
	posZ = posY +num_mark_particle_in_species[species];
	memset(posX,0,3*num_mark_particle_in_species[species]*sizeof(PARTICLE_REAL));

	v_ionX = new PARTICLE_REAL[3*num_mark_particle_in_species[species]];
	v_ionY = v_ionX +num_mark_particle_in_species[species];
	v_ionZ = v_ionY +num_mark_particle_in_species[species];
	memset(v_ionX,0,3*num_mark_particle_in_species[species]*sizeof(PARTICLE_REAL));

	//! -----------------------------------------------------------------
	//! - 1) Estimate positions of particle and store them in arrays ----
	//! -----------------------------------------------------------------
	for (INT32 ion=0; ion<num_mark_particle_in_species[species]; ion++)
	{
		
		//! NOTE be careful that each process generates same coordinates
		//! NOTE random generator for INFLOW species generates different coordinates for each process

		//! position particle into zCS starting at x = - 0.9 *Box_Origin[0];
		//posX[ion] = 0 + int(ion/50)*1./5.*R_Moon;
		//
		//posY[ion] = -1.*R_Moon + (ion%5)*1./2.*R_Moon;
		//
		//posZ[ion] = -1.5*R_Moon;
		//
		//v_ionX[ion] = 0 + int(0.1*ion)*200./SI_v0;
		//
		//log_file << "ion "<<ion<<"  ion%10 "<<ion%10<<" vx "<< v_ionX[ion]<<endl;
		//
		//v_ionY[ion] = 0;
		//v_ionZ[ion] = 0; 


	  posX[ion] = R_Obstacle;
	  posY[ion] = R_Obstacle;
	  posZ[ion] = R_Obstacle;
	  v_ionX[ion] = 0.000147283950676182;
	  v_ionY[ion] = 0;
	  v_ionZ[ion] = 0;
		
	}


	//! ----------------------------------------------------------------
	//! - 3) Distribute the positions to Blocks -------------------------
	//! -    (particles finally get initialzed at respective Block) -----
	//! -----------------------------------------------------------------
	PARTICLE_REAL start_weight = 1e-6;//(obstacle_MP_weight_each_t0[species] * dt)/(1.*num_mark_particle_in_species[species]);

	PARTICLE_REAL* positions[3]  = {  posX,   posY,   posZ};
	PARTICLE_REAL* velocities[3] = {v_ionX, v_ionY, v_ionZ};

	send_injected_Ions_to_Blks(species, start_weight, positions, velocities,num_mark_particle_in_species);
 
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " inject heavy ion time: " << time << "s." << endl << endl;


	delete[]   posX;
	delete[] v_ionX;

}


//!--------------------------------------------------------
//!- inject_callisto_particle_to_track:
//!--------------------------------------------------------
void CHybrid::inject_callisto_particle_to_track(INT32 species)
{

// 	if(TL!=TL_MARK_PARTICLE_TRACKS) return;
	
	log_file << "INJECTING CALLISTO PARTICLE TO TRACK ..." << endl;
	//! define which species shall be inserted


	if(num_Particle_Species<=species)
	{
		log_file << "  WARNING: " << endl
		     << "  species:           " << species << endl
		     << "  num_Particle_Species:   " << num_Particle_Species <<endl
		     << "  has to be at least " <<species+1<< endl
		     << " INJECTION CANCELED." << endl <<endl;
		return;

	} 


	clock_t start,finish;
	double time;
	start = clock();

	//! -----------------------------------------------------
	//! Define positions on sphere to inject
	//D_REAL theta, phi, v_theta, v_phi;
	D_REAL r = 1.0*R_Moon;
	D_REAL v_mag = 1.;

	//! -----------------------------------------------------
	//! -----------------------------------------------------
	//! TODO: Redefine positions and velocities to D_REAL.
	//! This whole allocation is fucked. I don't need arrays
	//! I just need a single double since the particle is injected
	//! at every value for a given v_phi. So no positions/velocities
	//! really need to be stored....
	//! -----------------------------------------------------
	//! --------- Allocate particle memory ------------------
	//! -----------------------------------------------------
	PARTICLE_REAL *posX, *posY, *posZ, *v_ionX, *v_ionY, *v_ionZ;
	//D_REAL posX, posY, posZ, v_ionX, v_ionY, v_ionZ;

	posX = new PARTICLE_REAL[3*num_mark_particle_in_species[species]];
	posY = posX +num_mark_particle_in_species[species];
	posZ = posY +num_mark_particle_in_species[species];
	memset(posX,0,3*num_mark_particle_in_species[species]*sizeof(PARTICLE_REAL));

	v_ionX = new PARTICLE_REAL[3*num_mark_particle_in_species[species]];
	v_ionY = v_ionX +num_mark_particle_in_species[species];
	v_ionZ = v_ionY +num_mark_particle_in_species[species];
	memset(v_ionX,0,3*num_mark_particle_in_species[species]*sizeof(PARTICLE_REAL));

	//! -----------------------------------------------------------------
	//! - 1) Estimate positions of particle and store them in arrays ----
	//! -----------------------------------------------------------------

	//! TODO: need to add a for loop here loping over every potision.
	//! want ALL ions in num_mark_particle at EVERY position on callisto

	//	D_REAL v_phi=44.*M_PI/180.;
	//	D_REAL v_theta=44.*M_PI/180.;

	INT32 ions = 0.;
	D_REAL theta_index = 10.;//2;
	D_REAL phi_index   = 20.;//4;

	for(D_REAL theta = 0.; theta < M_PI; theta=theta+theta_index*(M_PI/180.)){   //do i want theta + 1?
	  //log_file << "theta loop number: " << ion << endl;
	//  log_file << "theta: " << theta*180./M_PI << endl;
	  //ion = 0.;
	  for(D_REAL phi = 0.; phi < 2.0 * M_PI; phi=phi+phi_index*(M_PI/180.)){        //do i want phi + 1?
	    //log_file << "phi loop number: " << ion << endl;
	//    log_file << "phi: " << phi*180./M_PI << endl;
	   
	//all good above

	    for(D_REAL v_theta = 0.; v_theta < M_PI; v_theta=v_theta+theta_index*(M_PI/180.)){   //do i want theta + 1?
	      for(D_REAL v_phi = 0.; v_phi < 2.0 * M_PI; v_phi=v_phi+phi_index*(M_PI/180.)){        //do i want phi + 1?
//////
//////	//for(v_theta = 0; v_theta < 0.5 * M_PI; v_theta=v_theta+(M_PI/180.)){
//////	      //log_file << "v_theta loop number: " << ion << endl;
//////	//      log_file << "v_theta: " << v_theta*180./M_PI << endl;
//////	//      ion = 0;
//////	      //	theta = 0;
//////	      //	phi = 0;
//////	      //	v_theta = 0;
////////	      for(v_phi = 0; v_phi < 2.0 * M_PI; v_phi=v_phi+(M_PI/180.)){
////////		//log_file << "v_phi loop number: " << ion << endl;
////////		log_file << "v_phi: " << v_phi*180./M_PI << endl;
		//for(v_theta = 0; v_theta < M_PI; v_theta = v_theta + (2.*M_PI/180.)){
		//log_file << "theta: " << theta*180./M_PI << endl;
	  //ion = 0;
	  //for(v_phi = 0; v_phi < 2. * M_PI; v_phi = v_phi + (4.*M_PI/180.)){
		//log_file << "phi: " << phi*180./M_PI << endl;
		//log_file << "v_theta: " << v_theta*180./M_PI << endl;
		//log_file << "v_phi: " << v_phi*180./M_PI << endl;


	    //do these really need to be indexed at [ion]?
	    //i don't think so...
		posX[0] = r*cos(phi)*sin(theta);
		posY[0] = r*sin(phi)*sin(theta);
		posZ[0] = r*cos(theta);
//		//log_file << "here1" << endl;
		v_ionX[0] = v_mag*cos(v_phi)*sin(v_theta);
		v_ionY[0] = v_mag*sin(v_phi)*sin(v_theta);
		v_ionZ[0] = v_mag*cos(v_theta);

		/*
		if(theta + (2.*M_PI/180.) >= M_PI)
		  log_file << "theta max: " << theta*180./M_PI << endl;

		if(phi + (4.*M_PI/180.) >= 2.*M_PI)
		  log_file << "phi max: " << phi*180./M_PI << endl;

		if(v_theta + (2.*M_PI/180.) >= M_PI)
		  log_file << "v_theta max: " << v_theta*180./M_PI << endl;

		if(v_phi + (4.*M_PI/180.) >= 2.*M_PI)
		  log_file << "v_phi max: " << v_phi*180./M_PI << endl;
		*/

		//log_file << "here2" <<endl;
		//! Have to send positions and velocities
		//! to the blocks here, within the loop.
		//! THEN i can incriment ion, inject a new ion with
		//! new velocities (and resend pos/vel to blocks).
		//! After looping over all velocities at a specific phi position
		//! then incriment phi and start over.....
		PARTICLE_REAL start_weight = 1e-6; //(obstacle_MP_weight_each_t0[species] * dt)/(1.*num_mark_particle_in_species[species]);

		PARTICLE_REAL* positions[3]  = {  posX,   posY,   posZ}; 
		PARTICLE_REAL* velocities[3] = {v_ionX, v_ionY, v_ionZ};
		//log_file << "here3" << endl;
		send_injected_Ions_to_Blks(species, start_weight, positions, velocities, num_mark_particle_in_species);
		//log_file << "here4" << endl;

		//delete[] posX;
		//delete[] posY;
		//delete[] posZ;
		//delete[] v_ionX;		
		//delete[] v_ionY;
		//delete[] v_ionZ;
		//log_file << "hereA" << endl;
		ions++;
		//log_file << "hereB" << endl;
	  }
	}
		//}
		//}




	//i don't think i need to loop over num_mark.....
	//if it's done right, it'll be the same??
	//for (INT32 ion=0; ion<num_mark_particle_in_species[species]; ion++)
	  //{
		
		//! NOTE be careful that each process generates same coordinates
		//! NOTE random generator for INFLOW species generates different coordinates for each process

		//! position particle into zCS starting at x = - 0.9 *Box_Origin[0];
		//posX[ion] = 0 + int(ion/50)*1./5.*R_Moon;
		//
		//posY[ion] = -1.*R_Moon + (ion%5)*1./2.*R_Moon;
		//
		//posZ[ion] = -1.5*R_Moon;
		//
		//v_ionX[ion] = 0 + int(0.1*ion)*200./SI_v0;
		//
		//log_file << "ion "<<ion<<"  ion%10 "<<ion%10<<" vx "<< v_ionX[ion]<<endl;
		//
		//v_ionY[ion] = 0;
		//v_ionZ[ion] = 0; 


	  //posX[ion] = R_Obstacle;
	  //posY[ion] = R_Obstacle;
	  //posZ[ion] = R_Obstacle;
	  //v_ionX[ion] = 0.000147283950676182;
	  //v_ionY[ion] = 0;
	  //v_ionZ[ion] = 0;
		
	  }
	}


	//! THIS NEEDS TO BE DONE FOR EACH ION
	//! THEREFORE, SEND TO BLOCKS WITHIN THE FOR LOOP 

	//! ----------------------------------------------------------------
	//! - 3) Distribute the positions to Blocks -------------------------
	//! -    (particles finally get initialzed at respective Block) -----
	//! -----------------------------------------------------------------
	//PARTICLE_REAL start_weight = 1e-6;//(obstacle_MP_weight_each_t0[species] * dt)/(1.*num_mark_particle_in_species[species]);
	//
	//PARTICLE_REAL* positions[3]  = {  posX,   posY,   posZ};
	//PARTICLE_REAL* velocities[3] = {v_ionX, v_ionY, v_ionZ};
	//
	//send_injected_Ions_to_Blks(species, start_weight, positions, velocities,num_mark_particle_in_species);
 
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " inject CALLISTO TRACKED ion time: " << time << "s." << endl;
	log_file << "--> Number of ions of species " << species << " injected: " << ions << "." << endl << endl;

	//why aren't posY, posZ, v_ionY, v_ionZ deleted here??
	//delete[]   posX;
	//delete[] v_ionX;
	
}

//!--------------------------------------------------------
//!- inject_externRhoVelocity_ions:
//!--------------------------------------------------------
void CHybrid::insert_ions_from_neutral_profile(void)
{
  
	log_file << endl;
	log_file << " INJECTING OBSTACLE IONS FROM NEUTRAL DENSITY PROFILE(S):" << endl;	
	
	clock_t start,finish;
	double time;
	start = clock();
	
	INT64* num_injected_particles_species = new INT64[num_Particle_Species*num_Neutral_Species];
        memset(num_injected_particles_species, 0,num_Particle_Species*num_Neutral_Species*sizeof(INT64));
	
	//! adjust number of particles that are injected
	//! do this only at second time step (first is used to get number of inserted particles)
	if(TL==2 || TL==TL_at_last_restore_state+2)
	{	
		INT32 total_particles_to_insert=0;
		for(INT32 species=0; species<num_Particle_Species; species++)
		total_particles_to_insert += obstacle_MP_num_each_TL[species];
		
		insert_every_x_TL= int(insert_every_x_TL*num_newlyionised_particles/total_particles_to_insert + 1);	
	}
	
	
	if(ElectronionisationRate_Te_dependent)
	calc_electron_temperature(id_ElectronTemperature);
	
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
			if(mpi_myRank==temp_Block->responsible_mpi_process)
			for(INT32 oct=0; oct<8; oct++)
			{			
				if(!temp_Block->child_array[oct])
				temp_Block->insert_ions_from_neutral_profile(insert_every_x_TL,oct,num_injected_particles_species);
			}
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	

	// Provide Information about how many partilce were injected in total
	INT64 local_info_values[num_Particle_Species*num_Neutral_Species];
	stringstream info_names[num_Particle_Species*num_Neutral_Species];
	
	for(INT32 neutralSpec=0; neutralSpec<num_Neutral_Species; neutralSpec++)
	 for(INT32 species=0; species<num_Particle_Species; species++)
	 {
		 
		local_info_values[species*num_Neutral_Species + neutralSpec] = num_injected_particles_species[species*num_Neutral_Species + neutralSpec];
		info_names[species*num_Neutral_Species + neutralSpec] << "   ->injected particles of species "<<species<< " from neutral species "<<neutralSpec<<" : ";		
	 }	

	//! Reduce via MPI if required
	show_information(local_info_values,
		info_names,
		num_Particle_Species*num_Neutral_Species,
		BUILD_SUM);
	
	//! get number of inserted particles to adjust this number in next time step
	if(TL==1 || TL==TL_at_last_restore_state+1)
	for(INT32 neutralSpec=0; neutralSpec<num_Neutral_Species; neutralSpec++)
	 for(INT32 species=0; species<num_Particle_Species; species++)
	  num_newlyionised_particles += local_info_values[species*num_Neutral_Species + neutralSpec];
 
 
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << "  inject heavy ion time: " << time << "s." << endl << endl;

}



//!--------------------------------------------------------
//!- chemicalReactions:
//!--------------------------------------------------------
void CHybrid::chemical_Reactions(void)
{
            
        log_file << " PERFORMING CHEMICAL REACTIONS..." << endl;

	
        //! Preparation for recombination
        calc_RecombinationRate();

        INT64* num_exchanged_particle = new INT64[(num_Particle_Species)*(num_Particle_Species+1)];
        memset(num_exchanged_particle, 0,(num_Particle_Species)*(num_Particle_Species+1)*sizeof(INT64));
        
        clock_t start,finish;
        double time;
        start = clock();

	
        //! initial temp variable for charge exchange max beta
        max_ChemicalReaction_Probability = 0;
        memset(max_ChemicalReaction_Rates,0,(NUM_PARTICLE_SPECIES)*(NUM_PARTICLE_SPECIES+1)*sizeof(D_REAL));

	log_file << " start reactions ..." << endl;
	
        for(INT32 level=0; level<=MAX_LEVEL; level++)
        {
            CBlock* temp_Block = BlockList_of_Lev[level];
            while(temp_Block)
            {
                
                if(mpi_myRank == temp_Block->responsible_mpi_process)
                    temp_Block->chemical_Reactions(num_exchanged_particle);
                
                temp_Block = temp_Block->next_Blk_of_BlockList;
            }
        }
        
        log_file << "\tdone " << endl;
                  
	
        //! Provide Information about how many partilce were exchanged in each species
        INT64 local_info_values[(num_Particle_Species)*(num_Particle_Species+1)];
        stringstream info_names[(num_Particle_Species)*(num_Particle_Species+1)];
        
        INT32 counter=0;
        for(INT32 species=0; species<num_Particle_Species; species++)
	    if(species_does_chemical_reactions[species])
            for(INT32 destSpec=0; destSpec<(num_Particle_Species+1); destSpec++)
            {
		//! only display info about reactions that do really occur 
		//! NOTE TODO this already requires communication, otherwise not all processes have identical number of reactions occurring
// 		if(num_exchanged_particle[species*(num_Particle_Species+1)+destSpec]>0)
// 		{
			local_info_values[counter] = num_exchanged_particle[species*(num_Particle_Species+1)+destSpec];
			
			if(destSpec<num_Particle_Species)
			{	
				info_names[counter] << "   -> num particles "
				<< "S"    << species
				<< "->S"    << destSpec
				<< ": ";
			}
			else
			{
				info_names[counter] << "   -> num particles "
				<< "S"    << species
				<< " + e: ";
			}	
			counter++;
		
// 		}
            }
                         
        delete[] num_exchanged_particle;
        
	synchronize_allProcesses_MPI();
		

        //! Reduce via MPI if required
        show_information(local_info_values,
                            info_names,
                            counter,
                            BUILD_SUM);
       
	
        //! check if time step is ok for collision probability
        //! collision probability should be << 1, otherwise timestep inapproriate for describing collisions.
        //! But
        //! new ions have low velocity -> momentum loss due to several collisions within one time step is small
        //! inner region -> all ions have the same velocity
        if(check_max_reaction_probability || check_max_reaction_rates)
	{	
		D_REAL max_ChemicalReaction_Probability_global[(num_Particle_Species)*(num_Particle_Species+1)+1];
		memset(max_ChemicalReaction_Probability_global,0,((num_Particle_Species)*(num_Particle_Species+1)+1)*sizeof(D_REAL));
		D_REAL max_ChemicalReaction_Probability_local[(num_Particle_Species)*(num_Particle_Species+1)+1];
		memset(max_ChemicalReaction_Probability_local,0,((num_Particle_Species)*(num_Particle_Species+1)+1)*sizeof(D_REAL));

		//! first entry is total reaction probability
		max_ChemicalReaction_Probability_local[0] = max_ChemicalReaction_Probability;
		
		//! shift array entries by 1
		for(int i=0; i<(num_Particle_Species)*(num_Particle_Species+1);i++)
		max_ChemicalReaction_Probability_local[i+1]=max_ChemicalReaction_Rates[i];
	
		mpi_myComm.Allreduce(max_ChemicalReaction_Probability_local,
				max_ChemicalReaction_Probability_global,
				(num_Particle_Species)*(num_Particle_Species+1)+1,
				MPI_D_REAL,
				MPI_MAX);
		
		
		if(check_max_reaction_probability)
		log_file << " Max. Process Probability: " << max_ChemicalReaction_Probability_global[0] << endl;
		
		if(check_max_reaction_rates)
		{	
			log_file << " Max. Process Rates [1/s]: " << endl;
			log_file << " Species ";
			for(int species=0; species<num_Particle_Species; species++)
			{
				if(species==0)
				for(int destSpec=0; destSpec<num_Particle_Species+1; destSpec++)
				{					
						if(destSpec<num_Particle_Species)
						log_file<<" "<<destSpec<<" ";
						
						if(destSpec==num_Particle_Species)
						{
							log_file<<" e- ";
							log_file<<endl;
						}	
				}
				
				
				for(int destSpec=0; destSpec<num_Particle_Species+1; destSpec++)
				{
					if(destSpec==0)
					log_file<<"  "<<species<<"       ";
					
					log_file<<float(max_ChemicalReaction_Probability_global[species*(num_Particle_Species+1)+destSpec +1]/SI_t0)<<"   ";
				}
				log_file<<endl;
			}
		}
	}
	
        //! measure time
        finish = clock();
        time = (double(finish)-double(start))/CLOCKS_PER_SEC;
        log_file << " chemical reaction time: " << time << "s." << endl << endl;
}
