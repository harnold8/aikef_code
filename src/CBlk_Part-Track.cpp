#include <math.h>


#include <iostream>
#include <fstream>

#include "CBlk.h"
#include "defines.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"





//! ---------------------------------
//! Tracking the particles
//! ---------------------------------

//!-----------------------------------------------------------
//! mark_particle
//!-----------------------------------------------------------
void CBlock::mark_particle(INT64* particle_counter_of_species, INT64 num_partcle_in_list, INT64* Partilce_number_List)
{


//! only compile function in case particle tracing is use,
//! oterwise variable "number" does not exist in particle struct
#ifdef TRACK_PARTICLE

	PARTICLE_REAL x[3];
	PARTICLE_REAL v[3];

	particle *active_particle;

	for(INT32 i = 0; i < BlkNds_X; i++)
	 for(INT32 j = 0; j < BlkNds_Y; j++)
	  for(INT32 k = 0; k < BlkNds_Z; k++)
	  {
	
		INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
			     +j*BlkNds_Z 
			     +k;
	
	
		for(short species=0; species<num_Charged_Species; species++)
		  if(num_mark_particle_in_species[species] > 0)
		    {
		      for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
			{
			  
			  active_particle = pArray[species][i_j_k] +part_index;
			  INT32 cell_indices[3] = {i,j,k};
			  
			  //! Mark Particles that should be traced
			  intern2normedCoords(x, active_particle->rel_r, cell_indices);
			  // 			active_particle->number = -1;
			  
			  // 			only show in certain volume
			  
			  D_REAL propability = 0.1;
			  D_REAL dist = vec_len(x);
			  D_REAL vel = vec_len(v);
			  
			  log_file << "herehere " << species << endl;
			  
			  if(//   x[0] >  -LX && x[0] < -3.5*R_Obstacle
			     // && x[1] > -LY && x[1]  < +LY
			     // && x[2] > -LZ && x[2]  < +LZ
			     //dist > 1.3 * R_Obstacle && dist < 1.5 * R_Obstacle
			     (dist < 1.02*R_Obstacle*sqrt(3) || dist > 0.98*R_Obstacle*sqrt(3)) && (vel > 0.000147283 || vel < 0.000147285)
			     /*&& 1.*gsl_rng_uniform(randGen_general_asynchronous) < propability*/)
			    {
			      
			      // 				bool is_in_List = false;
			      // 
			      // 				for(INT64 counter=0; counter<num_partcle_in_list; counter++)
			      // 				 if(active_particle->number == Partilce_number_List[counter])
			      // 				 is_in_List = true;
			      // 
			      // 				if(!is_in_List)
			      // 				active_particle->number = -1;
			      // 				else
			      // 				{
			      // 					log_file << active_particle->number << "	";
			      // 					particle_counter_of_species[species]++;
			      // 
			      // 				}
			      
			      //! mark particle
			      active_particle->number = particle_counter_of_species[species];
			      particle_counter_of_species[species]++;
			      
			      
			    }
			  else
			    active_particle->number = -1;
			  
			}
		    }	
	  }//! for block knodes

#endif

}


//!-----------------------------------------------------------
//! mark_particle
//!-----------------------------------------------------------
void CBlock::write_particle_number(void)
{


//! only compile function in case particle tracing is use,
//! oterwise variable "number" does not exist in particle struct
#ifdef TRACK_PARTICLE

	PARTICLE_REAL x[3];


	char filename[200];
	ofstream ParticleNummber_FILE;

	sprintf(filename,"particle_numbers.txt");
	ParticleNummber_FILE.open(filename, ios_base::app);
	
	particle *active_particle;

	for(INT32 i = 0; i < BlkNds_X; i++)
	 for(INT32 j = 0; j < BlkNds_Y; j++)
	  for(INT32 k = 0; k < BlkNds_Z; k++)
	  {
	
		INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
			     +j*BlkNds_Z 
			     +k;
	
	
		for(short species=0; species<num_Charged_Species; species++)
		 for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
		 {
	
			active_particle = pArray[species][i_j_k] +part_index;
			INT32 cell_indices[3] = {i,j,k};
			

			intern2normedCoords(x, active_particle->rel_r, cell_indices);

			D_REAL propability = 0.1;
				
			if(   x[0] > -0.5*R_Obstacle && x[0] < 1.5*R_Obstacle
			    && x[1] > 1.5*R_Obstacle && x[1]  < 2.5*R_Obstacle
			    && x[2] > -1.*R_Obstacle && x[2]  < 1.*R_Obstacle
			    && 1.*gsl_rng_uniform(randGen_general_asynchronous) < propability)
			{
				//! mark particle
				ParticleNummber_FILE << active_particle->number << endl;

			}
	
		 }
		
	
	}//! for block knodes


	ParticleNummber_FILE.close();

#endif

}




//!-----------------------------------------------------------
//! trace_particle
//!-----------------------------------------------------------
void CBlock::trace_particle(INT64 *particle_counter_of_species)
{


//! only compile function in case particle tracing is use,
//! oterwise variable "number" does not exist in particle struct
#ifdef TRACK_PARTICLE




	char filename[200];
	PARTICLE_REAL x[3];
		
	particle *active_particle;


// 	if(binary_particle_tracks)
// 	{
// 		sprintf(filename,"%s/particle_tracks/Part_Tracks.bin",data_output_path);
// 		TraceParticle_FILE.open(filename, ios::binary | ios::app);
// 	}
	
	for(INT32 i = 0; i < BlkNds_X; i++)
	 for(INT32 j = 0; j < BlkNds_Y; j++)
	  for(INT32 k = 0; k < BlkNds_Z; k++)
	  {
	
		INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
			     +j*BlkNds_Z 
			     +k;

		for(short species=0; species<num_Charged_Species; species++)
		  if(num_mark_particle_in_species[species] > 0)
		    {
		      for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
			{
			  
			  active_particle = pArray[species][i_j_k] +part_index;
			  INT32 cell_indices[3] = {i,j,k};
			  
			  
			  intern2normedCoords(x, active_particle->rel_r, cell_indices);
			  
			  if(active_particle->number >= 0)
			    {
			      
			      
			      if(binary_particle_tracks)
				{
				  
				  sprintf(filename,"%s/particle_tracks/spec%d_part%07d.txt",data_output_path, species, active_particle->number);
				  TraceParticle_FILE.open(filename, ios_base::app);
				  TraceParticle_FILE.write(reinterpret_cast<char*> (x),3*sizeof(PARTICLE_REAL));
				  TraceParticle_FILE.close();
				  
				}
			      else
				{
				  sprintf(filename,"%s/particle_tracks/spec%d_part%07d.txt",data_output_path, species, active_particle->number);
				  TraceParticle_FILE.open(filename, ios_base::app);
				  
				  
				  //! write Time level 
				  TraceParticle_FILE << TL << "	";
				  
				  //! write particle number
				  //				  TraceParticle_FILE << active_particle->number << " ";

				  //! write position in normed units
				  TraceParticle_FILE << x[0] << " ";
				  TraceParticle_FILE << x[1] << " ";
				  TraceParticle_FILE << x[2] << " ";
				  
				  //! write velocity in normed units
				  TraceParticle_FILE << active_particle->v[0] << " ";
				  TraceParticle_FILE << active_particle->v[1] << " ";
				  TraceParticle_FILE << active_particle->v[2];
				  //				  TraceParticle_FILE << vec_len(active_particle->v);
				  
				  TraceParticle_FILE << endl;
				  
				  TraceParticle_FILE.close();
				}
			      
			      
			      particle_counter_of_species[species]++;
			      
			    }
			  
			}//! end for particle of cell
		    }
	  }//! end for knodes in block


// 	if(binary_particle_tracks)
// 	TraceParticle_FILE.close();

#endif
}
