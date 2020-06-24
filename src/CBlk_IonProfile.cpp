


#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>



#include "CBlk.h"
#include "utils.h"
#include "parameters.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "absolute_Globals.h"

using namespace std;


//! extern Global Variables
//! (also used from functions outside this file)

extern D_REAL *CellVol_of_L;
extern PARTICLE_REAL  *thermal_velocity_of;
extern D_REAL *dt_particle_of_L;
extern D_REAL **delta_of_L;

extern WEIGHT_REAL startWeight_of_smallestCell;


//! Global variables of this File
extern particle* active_particle;

extern PARTICLE_REAL vth[3];

//! MPI related variables
extern ofstream log_file;

//! Random Generators
extern gsl_rng **randGen_of_species;

//! extern profiles
extern D_REAL **extern_1D_Field;

//!-------------------------------------------------------------//
//! modifications of ion production for numerical reasons
//!-------------------------------------------------------------//
bool CBlock::insert_at_rnormed(D_REAL *r_normed)
{
	//! comet
	if(vec_len2(r_normed) > R_Moon*R_Moon)
	return true;
	else
	return false;
}
//!-------------------------------------------------------------//


//!-------------------------------------------------------------//
//! ion production rate
//!-------------------------------------------------------------//
D_REAL CBlock::getIonProdRateAtPoint(PARTICLE_REAL *r_intern,INT32 * cell_indices, INT32 ion_species, INT32 neutral_species) {
        
	//! here you may include modifications of the neutral density profile if these are not written in the ion production field
	
	D_REAL IonProdRate=0;	
	
#if defined(use_ion_production_as_field)
	IonProdRate = pickScalarFieldValue(id_density_ionProdSpecies1 + neutral_species, r_intern ,  cell_indices);
#else	
	IonProdRate = PhotoionisationRate[ion_species][neutral_species]*pickScalarFieldValue(id_numberdensity_neutralSpecies1 + neutral_species, r_intern ,  cell_indices);
#endif	
        return IonProdRate;
}
//!-------------------------------------------------------------//

//!-------------------------------------------------------------//
void CBlock::getIonProdVelocityAtPoint(PARTICLE_REAL *r_intern,INT32 * cell_indices, INT32 ion_species, INT32 neutral_species,PARTICLE_REAL *vec_v) {

	pickVectorFieldValue(id_velocity_neutralSpecies1 + neutral_species, r_intern ,  cell_indices, vec_v);
}
//!-------------------------------------------------------------//


//!-------------------------------------------------------------//
//! Inject ions according to neutral profile and ionization frequency         //
//!     (theoretically, no modifications should be necessary)                 //
//!-------------------------------------------------------------/


void CBlock::insert_ions_from_ionprod_file(INT32 insert_every_x_TL,INT32 id_oct, INT64* num_injected_particles_species)
{
  
  D_REAL IonProdRate_at_point=0;
  //  D_REAL neutral_density_at_point=0;
  //  D_REAL ElecIoRate = 0;

  INT32 cell_indices[3];
  INT32 i_j_k;

  particle new_ion;
  PARTICLE_REAL x_part[3], v_part[3], weight_part, r_normed[3];
  
//   PARTICLE_REAL vth[3];  
//   memset(vth, 0, 3*sizeof(PARTICLE_REAL));
  
  //! get index of oct (0 or 1 for each dimension)
  INT32 a = id_oct/4;
  INT32 b = (id_oct -a*4)/2;
  INT32 c = (id_oct -a*4 -b*2);
	
  //! In case Block is Box boundary, filling shall start at 0, else at 1
  INT32 start[3] = {!a +a*(BlkNds_X/2),
			!b +b*(BlkNds_Y/2),
			!c +c*(BlkNds_Z/2)};



  INT32 end[3]={BlkNds_X-1 -!a*((BlkNds_X-2)/2),
		BlkNds_Y-1 -!b*((BlkNds_Y-2)/2),
		BlkNds_Z-1 -!c*((BlkNds_Z-2)/2)};


       for(INT32 i = start[0]; i < end[0]; i++)
        for(INT32 j = start[1]; j < end[1]; j++)
           for(INT32 k = start[2]; k < end[2]; k++)
            {        
                i_j_k  = i*BlkNds_Y*BlkNds_Z 
			+j*BlkNds_Z 
			+k;
			
                 //do not insert in obstacle cells
		//! NOTE this if condition may be removed if there are problems: for bad resolution with respect to obstacle
		//! cells which are half inside obstacle and half outside obstacle may be important...
		if(!Flag[i_j_k])
		for(INT32 ion_field=0; ion_field<num_ion_prod_fields; ion_field++)
		 {
			//! Füge nur zufällig alle insert_every_x_TL Zeitschritte ein MP ein
			D_REAL zufall =  1.*random()/RAND_MAX; 
			if(zufall<1./(1.*insert_every_x_TL)) 
			{
				//! Zufallsposition in der Zelle
				x_part[0] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
				x_part[1] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
				x_part[2] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
				
				cell_indices[0] = i;
				cell_indices[1] = j;
				cell_indices[2] = k;				
					
				//! since many people will introduce some more or less artificial modifications to the ion insertion,
				//! this will be considered here:
				//! Berechnung der Position in r_normed[3]
				intern2normedCoords(r_normed, x_part, cell_indices);
// 				if(insert_at_rnormed(r_normed))
// 				{	
				
					IonProdRate_at_point = getIonProdRateAtPoint(x_part,cell_indices,ion_field,ion_field);	
					//uncomment this if it doesn't work: we need to get neutral densities for elec impact:
					//	neutral_density_at_point = getNeutralDensityAtPoint(x_part,cell_indices,ion_field);
					/*
					if(ElectronionisationRate_Te_dependent)
					  ElecIoRate = pickScalarFieldValue(id_ElectronTemperature,x_part,cell_indices)*ElectronionisationRate[1][ion_field];
					else
					  ElecIoRate = ElectronionisationRate[1][ion_field];

					*/
					weight_part = insert_every_x_TL*(IonProdRate_at_point*dt*SI_t0/SI_n0);//
									 //+neutral_density_at_point*dt*SI_t0*ElecIoRate);
							
					//! only continue if weight larger than zero
					if(weight_part > 0.)
					{
					
						memcpy(new_ion.rel_r, x_part, 3*sizeof(PARTICLE_REAL));		
										
						// Get velocity				
						memset(v_part,0,3*sizeof(PARTICLE_REAL));
// 						getIonProdVelocityAtPoint(x_part,cell_indices,neutralSpec,neutralSpec,v_part);					
					
						PARTICLE_REAL vth[3];
						//! it is important to initialize vth, oth. error in case 
						//! temperature of species is zero.
						memset(vth, 0 ,3*sizeof(PARTICLE_REAL));
						
						//! Generate initial thermal velocity if required
						//! -> DO THIS AT EVERY PROCESS REGARDLESS OF WHETHER
						//!    THIS PROCESS WILL INSTERT THE ION.
						//! -> IF THIS IS DONE ONLY AT THE ION-INSERTING PROCESS
						//!    RANDOM GENERATORS WILL RUN OUT OF SYNCHRONISATION
						vth[0]= gsl_ran_gaussian_ziggurat(randGen_of_species[ion_field+1], Neutral_vth[ion_field]);
						vth[1]= gsl_ran_gaussian_ziggurat(randGen_of_species[ion_field+1], Neutral_vth[ion_field]);
						vth[2]= gsl_ran_gaussian_ziggurat(randGen_of_species[ion_field+1], Neutral_vth[ion_field]);

						
						new_ion.v[0] = v_part[0] +vth[0];
						new_ion.v[1] = v_part[1] +vth[1];
						new_ion.v[2] = v_part[2] +vth[2];
										
						
						new_ion.weight = weight_part *CellVol_of_L[RLevel];
					
		
						add_particle_to_pArray(ion_field+num_Inflow_Species, i_j_k, &new_ion);

						num_injected_particles_species[ion_field]++;
						num_total_particles             ++;
						num_total_particles_in_L[RLevel]++;
						//! num_MPiC is increased in add_particle_to_pArray function
					}						
// 				}
			}   
		}//! end species 
	}//! end octant              
 
}




void CBlock::insert_ions_from_neutral_profile(INT32 insert_every_x_TL,INT32 id_oct, INT64* num_injected_particles_species)
{
  
  D_REAL IonProdRate_at_point=0;
  D_REAL neutral_density_at_point=0;
  D_REAL ElecIoRate = 0;

  INT32 cell_indices[3];
  INT32 i_j_k;

  particle new_ion;
  PARTICLE_REAL x_part[3], v_part[3], weight_part, r_normed[3];
  
  PARTICLE_REAL vth[3];  
  memset(vth, 0, 3*sizeof(PARTICLE_REAL));
  
  //! get index of oct (0 or 1 for each dimension)
  INT32 a = id_oct/4;
  INT32 b = (id_oct -a*4)/2;
  INT32 c = (id_oct -a*4 -b*2);
	
  //! In case Block is Box boundary, filling shall start at 0, else at 1
  INT32 start[3] = {!a +a*(BlkNds_X/2),
			!b +b*(BlkNds_Y/2),
			!c +c*(BlkNds_Z/2)};



  INT32 end[3]={BlkNds_X-1 -!a*((BlkNds_X-2)/2),
		BlkNds_Y-1 -!b*((BlkNds_Y-2)/2),
		BlkNds_Z-1 -!c*((BlkNds_Z-2)/2)};


       for(INT32 i = start[0]; i < end[0]; i++)
        for(INT32 j = start[1]; j < end[1]; j++)
           for(INT32 k = start[2]; k < end[2]; k++)
            {        
                i_j_k  = i*BlkNds_Y*BlkNds_Z 
			+j*BlkNds_Z 
			+k;
			
                 //do not insert in obstacle cells
		//! NOTE this if condition may be removed if there are problems: for bad resolution with respect to obstacle
		//! cells which are half inside obstacle and half outside obstacle may be important...
		if(!Flag[i_j_k])
		for(INT32 neutralSpec=0; neutralSpec<num_Neutral_Species; neutralSpec++)
		 for(INT32 species=0; species<num_Particle_Species; species++)	 
		 {
			//! Füge nur zufällig alle insert_every_x_TL Zeitschritte ein MP ein
			D_REAL zufall =  1.*random()/RAND_MAX; 
			if(zufall<1./(1.*insert_every_x_TL)) 
			{
				//! Zufallsposition in der Zelle
				x_part[0] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
				x_part[1] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
				x_part[2] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
				
				cell_indices[0] = i;
				cell_indices[1] = j;
				cell_indices[2] = k;				
					
				//! since many people will introduce some more or less artificial modifications to the ion insertion,
				//! this will be considered here:
				//! Berechnung der Position in r_normed[3]
				intern2normedCoords(r_normed, x_part, cell_indices);
				if(insert_at_rnormed(r_normed))
				{	
				
					IonProdRate_at_point = getIonProdRateAtPoint(x_part,cell_indices,species,neutralSpec);
					neutral_density_at_point = getNeutralDensityAtPoint(x_part,cell_indices,neutralSpec);
					
					if(ElectronionisationRate_Te_dependent)
					ElecIoRate = pickScalarFieldValue(id_ElectronTemperature,x_part,cell_indices)*ElectronionisationRate[species][neutralSpec];
					else
					ElecIoRate = ElectronionisationRate[species][neutralSpec];	
					// comment out this to include only electron impact ionization
										weight_part = insert_every_x_TL*(IonProdRate_at_point*dt*SI_t0/SI_n0
														+neutral_density_at_point*dt*SI_t0*ElecIoRate);

					//					weight_part = insert_every_x_TL*(neutral_density_at_point*dt*SI_t0*ElecIoRate);
							
					//! only continue if weight larger than zero
					if(weight_part > 0.)
					{
					
						memcpy(new_ion.rel_r, x_part, 3*sizeof(PARTICLE_REAL));		
										
						// Get velocity				
						memset(v_part,0,3*sizeof(PARTICLE_REAL));
						getIonProdVelocityAtPoint(x_part,cell_indices,species,neutralSpec,v_part);					
					
						PARTICLE_REAL vth[3];
						//! it is important to initialize vth, oth. error in case 
						//! temperature of species is zero.
						memset(vth, 0 ,3*sizeof(PARTICLE_REAL));
						
						//! Generate initial thermal velocity if required
						//! -> DO THIS AT EVERY PROCESS REGARDLESS OF WHETHER
						//!    THIS PROCESS WILL INSTERT THE ION.
						//! -> IF THIS IS DONE ONLY AT THE ION-INSERTING PROCESS
						//!    RANDOM GENERATORS WILL RUN OUT OF SYNCHRONISATION
						
						// TO DO!:
						//The vth lines need to be commented out, not implemented properly, need to check initialization of Neutral_vth:

						//vth[0]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], Neutral_vth[neutralSpec]);
						//vth[1]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], Neutral_vth[neutralSpec]);
						//vth[2]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], Neutral_vth[neutralSpec]);

						
						new_ion.v[0] = v_part[0] +vth[0];
						new_ion.v[1] = v_part[1] +vth[1];
						new_ion.v[2] = v_part[2] +vth[2];
										
						
						new_ion.weight = weight_part *CellVol_of_L[RLevel];
					
		
						add_particle_to_pArray(species, i_j_k, &new_ion);

						num_injected_particles_species[species*num_Neutral_Species+neutralSpec]++;
						num_total_particles             ++;
						num_total_particles_in_L[RLevel]++;
						//! num_MPiC is increased in add_particle_to_pArray function
					}						
				}
			}   
		}//! end species 
	}//! end octant              
 
}


//!-------------------------------------------------------------//
//! getNeutralDensityAtPoint(D_REAL *r_normed,INT32 i_j_k)      //
//!-------------------------------------------------------------//
D_REAL CBlock::getNeutralDensityAtPoint(PARTICLE_REAL *r_intern,INT32 * cell_indices, INT32 neutral_species) {
        D_REAL NeutralDensity = pickScalarFieldValue(id_numberdensity_neutralSpecies1 + neutral_species, r_intern ,  cell_indices);
        return NeutralDensity;
}

//!-------------------------------------------------------------//
//! getNeutralVelocityAtPoint(D_REAL *r_normed,INT32 i_j_k,PARTICLE_REAL *vec_v)      //
//!-------------------------------------------------------------//
void CBlock::getNeutralVelocityAtPoint(PARTICLE_REAL *r_intern,INT32 * cell_indices, INT32 neutral_species,PARTICLE_REAL *vec_v) {
//         INT32 neutral_species = ion_species-1;
        pickVectorFieldValue(id_velocity_neutralSpecies1 + neutral_species, r_intern ,  cell_indices, vec_v);
}

