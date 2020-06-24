


#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

#include "CBlk.h"
#include "parameters.h"
#include "absolute_Globals.h"
#include "utils.h"




//!-----------------------------------------------------------
//! check_weight_sorting
//!-----------------------------------------------------------
void CBlock::check_weight_sorting(void)
{
	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


  INT32 i,j,k, i_j_k;

  for(i = 1; i < BlkNds_X-1; i++)
   for(j = 1; j < BlkNds_Y-1; j++)
    for(k = 1; k < BlkNds_Z-1; k++)
    {

	i_j_k  = i*BlkNds_Y*BlkNds_Z 
		+j*BlkNds_Z 
		+k;


        for(INT32 species=0; species<num_Charged_Species; species++)
	//! finish at penultimate (no neighbur to compare with)
      	 for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]-1; part_index++)
	   if(  pArray[species][i_j_k][part_index].weight
	       > pArray[species][i_j_k][part_index+1].weight)
	      {

		log_file << "sorting error: " << endl;
		log_file << pArray[species][i_j_k][part_index].weight <<">" 
		     << pArray[species][i_j_k][part_index+1].weight<< endl;
		log_file << "i: " << i << endl;
		log_file << "j: " << j << endl;
		log_file << "k: " << k << endl;
		log_file << "part_index: " << part_index << endl;
		log_file << "num_MPiC[species][i_j_k]: " << num_MPiC[species][i_j_k] << endl;

		for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
		log_file << pArray[species][i_j_k][part_index].weight << endl;


		exit(1);
	      }





     }


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}




//!-----------------------------------------------------------
//! estimate_Cords_of_split_Particle_CoM_OFF
//!-----------------------------------------------------------
void CBlock::estimate_Cords_of_split_Particle_CoM_OFF(PARTICLE_REAL* x_new1,
						      PARTICLE_REAL* x_new2,
						      PARTICLE_REAL* x_old,
						      PARTICLE_REAL* v_part)
{


	PARTICLE_REAL ran1, ran2, ran3, ran4, rez_x_diff_length, rez_x_diff_length2;
	PARTICLE_REAL temp[3], x_diff[3], x_diff2[3];


// 	memcpy(temp,v_part,3*sizeof(PARTICLE_REAL));

	//!-------------------------------------------------------------
	//! generate random vector "x_diff" which is perpendicular to v
	//!-------------------------------------------------------------
	temp[0] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	temp[1] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	temp[2] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	vec_cross(x_diff,v_part,temp);
	rez_x_diff_length = 1. / vec_len(x_diff);




	//! normalize vector
	x_diff[0] *= rez_x_diff_length;
	x_diff[1] *= rez_x_diff_length;
	x_diff[2] *= rez_x_diff_length;

	//! build vector that is perpendicular to v AND x_diff
	vec_cross(x_diff2,v_part,x_diff);
	rez_x_diff_length2 = 1. / vec_len(x_diff2);



	//! normalize vector
	x_diff2[0] *= rez_x_diff_length2;
	x_diff2[1] *= rez_x_diff_length2;
	x_diff2[2] *= rez_x_diff_length2;

	//! x_diff & x_diff2 span a plane perpendicular to v 

	//!-------------------------------------------------------------
	//! add with random scale to old position
	//!-------------------------------------------------------------
	memcpy(x_new1,x_old,3*sizeof(PARTICLE_REAL));
	memcpy(x_new2,x_old,3*sizeof(PARTICLE_REAL));

// 	if(x_diff[0] > 1. || x_diff[1] > 1. || x_diff[2] > 1.)
// 	return;



	//! NOTE:
	//! Center of mass is only conserved in case
	//! one random number is used !!!
	//! use some constant offset such that ran1 > 0 !!!
	//! -> else particle can be at nearly same positions
	//!    which will make the code crash for some reason
	//!    (NaN velocities -> yet unclear why)4

	//! NOTE:
	//! When conserving centre of mass or macroscopic density, there will ALWAYS be
	//! one problem that cannot be overcome:
	//! - assume splitted particles are shifted by a certain dX that is > 0.1
	//! -> in general ions in cell centres are preferred rather than those at the edges
	//! -> cell centres tend to be emptied out while edge near regions are filled
	//! -> no problem at this level of refinement
	//! -> moving to the next level, "density stripes" will occur due to artificially
	//!    created high dens and low dens regions (sketch should make things clear)
	
// 	ran1 = 0.1 +0.9*gsl_rng_uniform(randGen_general_asynchronous);
// 
// 	x_new1[0] -= x_diff[0]*ran1;
// 	x_new1[1] -= x_diff[1]*ran1;
// 	x_new1[2] -= x_diff[2]*ran1;
// 	
// 	x_new2[0] += x_diff[0]*ran1;
// 	x_new2[1] += x_diff[1]*ran1;
// 	x_new2[2] += x_diff[2]*ran1;


	//! insert new particle at arbitrary position within the plane
	//! perpendicular to v
	ran1 = 2.*gsl_rng_uniform(randGen_general_asynchronous) -1.;
	ran2 = 2.*gsl_rng_uniform(randGen_general_asynchronous) -1.;
	ran3 = 2.*gsl_rng_uniform(randGen_general_asynchronous) -1.;
	ran4 = 2.*gsl_rng_uniform(randGen_general_asynchronous) -1.;

	for(int comp=0; comp<3; comp++)
	{

		x_new1[comp] = x_new1[comp] +x_diff[comp]*ran1 +x_diff2[comp]*ran2;
		x_new2[comp] = x_new2[comp] +x_diff[comp]*ran3 +x_diff2[comp]*ran4;
	}




	//! set x_new1 out of cell to intercept in are_in_cell
// 	if( vec_len(x_diff) > 1.+PART_REAL_PRECISION || vec_len(x_diff) < 1.-PART_REAL_PRECISION)
// 	x_new1[0] == 200.;


}



//!-----------------------------------------------------------
//! estimate_Cords_of_split_Particle_CoM_OFF
//!-----------------------------------------------------------
void CBlock::estimate_Cords_of_split_Particle_CoM_ON(PARTICLE_REAL* x_new1,
						     PARTICLE_REAL* x_new2,
						     PARTICLE_REAL* x_old,
						     PARTICLE_REAL* v_part)
{


	PARTICLE_REAL ran1, rez_x_diff_length;
	PARTICLE_REAL temp[3], x_diff[3];


// 	memcpy(temp,v_part,3*sizeof(PARTICLE_REAL));

	//!-------------------------------------------------------------
	//! generate random vector "x_diff" which is perpendicular to v
	//!-------------------------------------------------------------
	temp[0] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	temp[1] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	temp[2] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	vec_cross(x_diff,v_part,temp);
	rez_x_diff_length = 1. / vec_len(x_diff);



	//! normalize vector
	x_diff[0] *= rez_x_diff_length;
	x_diff[1] *= rez_x_diff_length;
	x_diff[2] *= rez_x_diff_length;



	//! x_diff & x_diff2 span a plane perpendicular to v 

	//!-------------------------------------------------------------
	//! add with random scale to old position
	//!-------------------------------------------------------------
	memcpy(x_new1,x_old,3*sizeof(PARTICLE_REAL));
	memcpy(x_new2,x_old,3*sizeof(PARTICLE_REAL));



	//! NOTE:
	//! Center of mass is only conserved in case
	//! one random number is used !!!
	//! use some constant offset such that ran1 > 0 !!!
	//! -> else particle can be at nearly same positions
	//!    which will make the code crash for some reason
	//!    (NaN velocities -> yet unclear why)4

	//! NOTE:
	//! When conserving centre of mass or macroscopic density, there will ALWAYS be
	//! one problem that cannot be overcome:
	//! - assume splitted particles are shifted by a certain dX that is > 0.1
	//! -> in general ions in cell centres are preferred rather than those at the edges
	//! -> cell centres tend to be emptied out while edge near regions are filled
	//! -> no problem at this level of refinement
	//! -> moving to the next level, "density stripes" will occur due to artificially
	//!    created high dens and low dens regions (sketch should make things clear)
	
	ran1 = 0.1 +0.9*gsl_rng_uniform(randGen_general_asynchronous);

	x_new1[0] -= x_diff[0]*ran1;
	x_new1[1] -= x_diff[1]*ran1;
	x_new1[2] -= x_diff[2]*ran1;
	
	x_new2[0] += x_diff[0]*ran1;
	x_new2[1] += x_diff[1]*ran1;
	x_new2[2] += x_diff[2]*ran1;





}



//!-----------------------------------------------------------
//! estimate_Cords_of_split_Particle_CoM_OFF
//!-----------------------------------------------------------
void CBlock::set_Cords_perp_to_velocity_CoM_OFF(PARTICLE_REAL* x_part,
						PARTICLE_REAL* v_part)
{


	PARTICLE_REAL ran1, rez_x_diff_length;
	PARTICLE_REAL temp[3], x_diff[3];




	//!-------------------------------------------------------------
	//! generate random vector "x_diff" which is perpendicular to v
	//!-------------------------------------------------------------
	temp[0] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	temp[1] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	temp[2] = 1.*gsl_rng_uniform(randGen_general_asynchronous);
	vec_cross(x_diff,v_part,temp);
	rez_x_diff_length = 1. / vec_len(x_diff);

	

	//! normalize vector
	x_diff[0] *= rez_x_diff_length;
	x_diff[1] *= rez_x_diff_length;
	x_diff[2] *= rez_x_diff_length;


	//! insert new particle at arbitrary position within the plane
	//! perpendicular to v
	ran1 = 2.*gsl_rng_uniform(randGen_general_asynchronous) -1.;

	for(int comp=0; comp<3; comp++)
	x_part[comp] += x_diff[comp]*ran1;

	
}


//!-----------------------------------------------------------
//! are_in_cell
//!-----------------------------------------------------------
inline bool is_in_cell(PARTICLE_REAL *x1)
{
	
      if(x1[0]<0.+PART_REAL_PRECISION || x1[0]>=1.-PART_REAL_PRECISION ||
	 x1[1]<0.+PART_REAL_PRECISION || x1[1]>=1.-PART_REAL_PRECISION ||
	 x1[2]<0.+PART_REAL_PRECISION || x1[2]>=1.-PART_REAL_PRECISION  )
      return false;



      //! for some reason the above expression does NOT intercept
      //! NaN numbers. So the interception below is required also.
      if(    vec_len2(x1) *0. !=0. || vec_len2(x1) > 10.)
      {

#ifdef DEBUG_MODUS
// 		log_file << endl;
// 		log_file << "----------------------------------------------------------" << endl;
// 		log_file << "- WARNING: NaN coordinates detected and rejected.        -" << endl;
// 		log_file << "-        This might result in errors on some compilers.  -" << endl;
// 		log_file << "----------------------------------------------------------" << endl;
// 		log_file << endl;
#endif
		return false;
      }




      return true;

}


//!-----------------------------------------------------------
//! are_in_cell
//!-----------------------------------------------------------
inline bool are_in_cell(PARTICLE_REAL *x1, PARTICLE_REAL *x2)
{
	
      if(x1[0]<0.+PART_REAL_PRECISION || x1[0]>=1.-PART_REAL_PRECISION ||
	 x1[1]<0.+PART_REAL_PRECISION || x1[1]>=1.-PART_REAL_PRECISION ||
	 x1[2]<0.+PART_REAL_PRECISION || x1[2]>=1.-PART_REAL_PRECISION ||

	 x2[0]<0.+PART_REAL_PRECISION || x2[0]>=1.-PART_REAL_PRECISION ||
	 x2[1]<0.+PART_REAL_PRECISION || x2[1]>=1.-PART_REAL_PRECISION ||
	 x2[2]<0.+PART_REAL_PRECISION || x2[2]>=1.-PART_REAL_PRECISION
	 )
      return false;



      //! for some reason the above expression does NOT intercept
      //! NaN numbers. So the interception below is required also.
      if(    vec_len2(x1) *0. !=0. || vec_len2(x1) > 10.
          || vec_len2(x2) *0. !=0. || vec_len2(x2) > 10.)
      {

#ifdef DEBUG_MODUS
// 		log_file << endl;
// 		log_file << "----------------------------------------------------------" << endl;
// 		log_file << "- WARNING: NaN coordinates detected and rejected.        -" << endl;
// 		log_file << "-        This might result in errors on some compilers.  -" << endl;
// 		log_file << "----------------------------------------------------------" << endl;
// 		log_file << endl;
#endif
		return false;
      }



      if(    x1[0] == x2[0]
          && x1[1] == x2[1]
          && x1[2] == x2[2])
      {
#ifdef DEBUG_MODUS
// 	log_file << "equal particle positions !!!! " << endl;
#endif
	return false;
      }



      return true;

}


//!-----------------------------------------------------------
//! split_lowest_weight_particle
//!-----------------------------------------------------------
void CBlock::split_heaviest_particle(void)
{


	//! NOTE: if a coarse mesh is used, for some reason this method produces "stripes" of density enhancement
	//! 	    at least in case a coarse mesh is used (delta cell > rg)


//    log_file << startWeight_of_smallestCell << "                     ";

  particle *part2Split;

  INT32 i,j,k, i_j_k;




  PARTICLE_REAL x_new1[3], x_new2[3];

  //! Instead of spliting 1->8 particle, it
  //! turned out to better splitting 1->2 and repeat 
  //! this 8 times.
  //! Doing so, every time the heaviest remaining particle 
  //! is newly estimated and split, which leads to better 
  //! weight balancing.
  for(i = 1; i < BlkNds_X-1; i++)
   for(j = 1; j < BlkNds_Y-1; j++)
    for(k = 1; k < BlkNds_Z-1; k++)
     {

	i_j_k  = i*BlkNds_Y*BlkNds_Z 
		+j*BlkNds_Z 
		+k;


      for(short species=0; species<num_Charged_Species; species++)
       if(do_split_in_species[species])
       {

        for(INT32 split_cycle=0; split_cycle<num_split_each_cell; split_cycle++)
	 if(       num_MPiC[species][i_j_k] > 2
 	    &&  1.*num_MPiC[species][i_j_k] < 1.*Blk_optimal_MPiC[species] *(0.75 - 0.5*gsl_rng_uniform(randGen_general_asynchronous)))
	      {



		
		//!-----------------------------------------------------------
		//!1) Search step:
		//!   Find Particle with maximal weight
		//!-----------------------------------------------------------
#ifdef SORT_BY_WEIGHT
		//! heaviest particle is last particle of pArray
		part2Split = pArray[species][i_j_k] +num_MPiC[species][i_j_k]-1;
#else


		WEIGHT_REAL max_Weight = -1.;
	
		//! estimate heaviest particle in List
		for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
		{
	
			part2Split = pArray[species][i_j_k] +part_index;
	
			if( part2Split->weight > max_Weight)
			{
				max_Weight = part2Split->weight;
				part2Split_index = part_index;
			}
	
	
		}
	
		part2Split = pArray[species][i_j_k] +part2Split_index;

#endif

	
		
		//!-----------------------------------------------------------
		//!2) Split step:
		//!   Split Particle into 2
		//!-----------------------------------------------------------
		//! check if particle weight to small to split
#ifdef TRACK_PARTICLE
		if(part2Split->weight < 1.e-7*startWeight_of_smallestCell)
#else
		if(part2Split->weight < 1.e-7*startWeight_of_smallestCell)
#endif
		{
			split_cycle = num_split_each_cell;
			num_split_canceled_low_weight++;
		}
		else
		{
	

		  //! The code below conserves the centre of mass,
		  //! but results in a terrble ppc balancing for
		  //! finer cells.
		  //! TODO:
		  //! Maybe x_diff should be scaled such that it always 
		  //! fits in the cell and for each particle seperately 
		  //! (of course this would not conserve the CoM any more,
		  //! but pe perp to velocity !)
	
		  //! generate x_diff perp to v
		
		

		  //! STIL REMAINING PROBLEMS WITH THE CODE BELOW:
		  //! Most particles are split when they enter
		  //! a Lp1 Blk, so they are at the entree 
		  //! of a cell.
		  //! When they are split fully randomly,
		  //! many of them will be put at the end 
		  //! of the cell which results in a lack
		  //! of density at the entree of the cell.
// 		  x_new1[0] = 1. *gsl_rng_uniform(randGen_general_asynchronous);
// 		  x_new1[1] = 1. *gsl_rng_uniform(randGen_general_asynchronous);
// 		  x_new1[2] = 1. *gsl_rng_uniform(randGen_general_asynchronous);
// 
// 		  x_new2[0] = 1. *gsl_rng_uniform(randGen_general_asynchronous);
// 		  x_new2[1] = 1. *gsl_rng_uniform(randGen_general_asynchronous);
// 		  x_new2[2] = 1. *gsl_rng_uniform(randGen_general_asynchronous);



		//! Estimate new coords end ckeck if they are in cell
		if(conserve_CoM_at_splitting)
		 for(INT32 num_try = 0; num_try<num_randomize_position_split; num_try++)
		 {


			estimate_Cords_of_split_Particle_CoM_ON(x_new1,
								x_new2,
								part2Split->rel_r,
								part2Split->v);

			if(are_in_cell(x_new1, x_new2))
			break;
		 }

		//! allow to set particle somewhere on plane without conserving CoM
		//! -> in order to avoid "stripes" in moments
		if(!conserve_CoM_at_splitting)
		{


			//! copy position of part2split to new positions
			memcpy(x_new1, part2Split->rel_r, 3*sizeof(PARTICLE_REAL));
			memcpy(x_new2, part2Split->rel_r, 3*sizeof(PARTICLE_REAL));

			//! randomize position for first particle
			for(INT32 num_try = 0; num_try<num_randomize_position_split; num_try++)
			{
				set_Cords_perp_to_velocity_CoM_OFF(x_new1, part2Split->v);

				if(is_in_cell(x_new1))
				num_try = num_randomize_position_split;
			}

			//! randomize position for second particle, if first is in cell
			if(is_in_cell(x_new1))
			for(INT32 num_try = 0; num_try<num_randomize_position_split; num_try++)
			{
	
				set_Cords_perp_to_velocity_CoM_OFF(x_new2, part2Split->v);

				if(is_in_cell(x_new2))
				num_try = num_randomize_position_split;

			}

		}




		//! check whether new cords are in old cell
		if(are_in_cell(x_new1, x_new2))
		{

		


			//! alloc further memory if required
			if(size_of_pArray[species][i_j_k] < num_MPiC[species][i_j_k]+1)
			{

				//! alloc new pArray of larger size
				resize_pArray(species, i_j_k, int(min_pArray_Size_factor *(num_MPiC[species][i_j_k]+1)));

				//! adress has changed during pArray re-allocation -> reset pointer !!!
				part2Split = pArray[species][i_j_k] +num_MPiC[species][i_j_k]-1;

			}


			//! half the weight for each particle
			part2Split->weight *= 0.5;


			//! SPLIT PARTICLE
			//! Store particle on temporary buffer
			particle part_buffer[2];

			memcpy(part_buffer +0, part2Split, sizeof(particle));
			memcpy(part_buffer +1, part2Split, sizeof(particle));

			memcpy(part_buffer[0].rel_r, x_new1, 3*sizeof(PARTICLE_REAL));
			memcpy(part_buffer[1].rel_r, x_new2, 3*sizeof(PARTICLE_REAL));



#ifdef SORT_BY_WEIGHT

			//! DELETE part2Split:
			//! part2Split is last of pArray.
			//! delete particle by decrementing num_MPiC
			num_MPiC[species][i_j_k]--;


			//! Sort particle in List
			//! Estimate first particle in list which weight is larger
			//! than part2Split->weight
			//! Loop from List end, start at position num_MPiC-2 to
			//! ignore part2Split which is at position num_MPiC-1
			//! Allow for temporary negative index (case list is empty)
			INT32 new_part_index;
			for(new_part_index=num_MPiC[species][i_j_k]-1; new_part_index>=0; new_part_index--)
			if(part2Split->weight >= pArray[species][i_j_k][new_part_index].weight)
			break;
	
			//! insert new particle behind estimated particle
			new_part_index++;


			//! shift particle_list to get space for 2 new particle.
			//! since memory overlaps, use memmove rather than memcpy.
			//! (memcpy may be undefined and IS NOT measurable faster)
			memmove(pArray[species][i_j_k]  +(new_part_index+2),
				pArray[species][i_j_k]  +(new_part_index),
				(num_MPiC[species][i_j_k] -new_part_index)*sizeof(particle));

			//! insert particle at position new_part_index
			memcpy(pArray[species][i_j_k]  +new_part_index,
			       part_buffer,
			       2*sizeof(particle));

#else


			//! shift particle_list to get space for 1 new particle right behind p2split
			memmove(pArray[species][i_j_k]  +(part2Split_index+2),
				pArray[species][i_j_k]  +(part2Split_index+1),
				(num_MPiC[species][i_j_k] -(part2Split_index+1))*sizeof(particle));

			//! insert both particles at position part2Split_index
			memcpy(pArray[species][i_j_k]  +part2Split_index,
			       part_buffer,
			       2*sizeof(particle));

			//! subtract one, add 2 as below
			num_MPiC[species][i_j_k]--;


#endif

		
			//! update statistics
			num_split_in_L[RLevel]++;
			num_MPiC[species][i_j_k]+=2;
			num_total_particles_in_L[RLevel]++;
			num_total_particles++;



		 }//! end if in cell
		 else
		 num_split_canceled_out_of_cell++;

       		}//! end if canceled by weight

       }//! end if criteria fulfilled
      }//! end for species
     }//! end foor Box

}





//!-----------------------------------------------------------
//! split_most_centred_particle
//!-----------------------------------------------------------
void CBlock::split_most_centred_particle(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

  particle *part2Split;


  INT32 i,j,k, i_j_k, part2Split_index=-1;



  PARTICLE_REAL x_new1[3], x_new2[3];

  //! Instead of spliting 1->8 particle, it
  //! turned out to better splitting 1->2 and repeat 
  //! this 8 times.
  //! Doing so, every time the heaviest remaining particle 
  //! is newly estimated and split, which leads to better 
  //! weight balancing.


  for(i = 1; i < BlkNds_X-1; i++)
   for(j = 1; j < BlkNds_Y-1; j++)
    for(k = 1; k < BlkNds_Z-1; k++)
     {

	i_j_k  = i*BlkNds_Y*BlkNds_Z 
		+j*BlkNds_Z 
		+k;



      for(short species=0; species<num_Charged_Species; species++)
       if(do_split_in_species[species])
       {

        for(INT32 split_cycle=0; split_cycle<num_split_each_cell; split_cycle++)
	 if(       num_MPiC[species][i_j_k] > 2
 	    &&  1.*num_MPiC[species][i_j_k] < Blk_optimal_MPiC[species] *(0.7 - 0.5*gsl_rng_uniform(randGen_general_asynchronous)))
	      {



		
		//!-----------------------------------------------------------
		//!1) Search step:
		//!   Find Particle that is best centred within the cell
		//!   (check for each component)
		//!-----------------------------------------------------------


		//! choose most centered particle

		int comp=-1;

		WEIGHT_REAL centered = 1.e9;

		//! NOTE: the selection method below for some reason produces "stripes" of density enhancement
		//! and therefore should not be used

// 		for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
// 		{
// 	
// 			part2Split = pArray[species][i_j_k] +part_index;
// 	
// 			//! check for each component
// 			for(INT32 test_comp=0; test_comp<3; test_comp++)
// 			if( fabs(part2Split->rel_r[test_comp] -0.5) < centered)
// 			{
// 				centered = fabs(part2Split->rel_r[test_comp] -0.5);
// 
// 				comp = test_comp;
// 				part2Split_index = part_index;
// 			}
// 		}

		PARTICLE_REAL center[3];
		for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
		{
	
			part2Split = pArray[species][i_j_k] +part_index;
	
			//! check for each component
			for(INT32 test_comp=0; test_comp<3; test_comp++)
			center[test_comp] = part2Split->rel_r[test_comp] -0.5;


			if( vec_len2(center) < centered)
			{
				centered = vec_len2(center);
				part2Split_index = part_index;
			}
		}


	
		part2Split = pArray[species][i_j_k] +part2Split_index;



	
		
		//!-----------------------------------------------------------
		//!2) Split step:
		//!   Split Particle into 2
		//!-----------------------------------------------------------
		//! check if particle weight to small to split
#ifdef TRACK_PARTICLE
		if(part2Split->number || part2Split->weight < 1.e-7*startWeight_of_smallestCell)
#else
		if(part2Split->weight < 1.e-7*startWeight_of_smallestCell)
#endif
		{
			split_cycle = num_split_each_cell;
			num_split_canceled_low_weight++;
		}
		else
		{


			bool coords_in_cell = false;
	
			memcpy(x_new1,part2Split->rel_r,3*sizeof(PARTICLE_REAL));
			memcpy(x_new2,part2Split->rel_r,3*sizeof(PARTICLE_REAL));
		
			//! Note:
			//! In case particle are split and put close to boundaries (either 0 or 1)
			//! stripe at "min" cells appear
			//! In case particle are split and put close to middle of cell (0.5)
			//! stripe at "plus" cells appear
			//! Best is to split particle that are centred in the cell and put them to
			//! +- 0.25 delta cell
			//! Use this to estimate optimal particle to Split
		

			//! TODO:
			//! SOLVED PROBLEM BELOW:
			//! 
			//! In case shift = 0 density with split method is identical to density without
			//! (obsolet: in case shift = 0. results should look like splitting turned off
			//! However, it does not slightly different)

			comp = int(3.*gsl_rng_uniform(randGen_general_asynchronous));

			const PARTICLE_REAL shift = 0.25;
		
			x_new1[comp] -= shift;
			x_new2[comp] += shift;


			if(  x_new1[comp]>=0.0 +PART_REAL_PRECISION && x_new2[comp]<1.-PART_REAL_PRECISION)
			coords_in_cell = true;

			//! the larger MPiC particle, the smaller interval can be chosen
			//! DON'T CHOOSE SMALLER THAN 0.04
			//! TODO:
			//! COUNT SPLIT IN RESPECTIVE DIRECTION X,Y,Z SEPERATELY !!!
			//! NOTE:
			//! IS IS IMPORTENT THAT PARTICLE HAVE ENOUGH TIME TO ADJUST
			//! SO dt MUST NOT EXEED COURANT CRETERIA BY MORE THAN FACTOR OF 5
			if(fabs(part2Split->rel_r[comp]-0.5) >= 0.04)
			coords_in_cell = false;


			//! check whether new cords are in old cell
			if(coords_in_cell)
			{
	
			
				//! alloc further memory if required
				if(size_of_pArray[species][i_j_k] < num_MPiC[species][i_j_k]+1)
				{
	
					//! alloc new pArray of larger size
					resize_pArray
					(
						  species,
						  i_j_k,
						  int(min_pArray_Size_factor *(num_MPiC[species][i_j_k]+1))
					);
	
					//! adress has changed -> reset pointer !!!
					part2Split = pArray[species][i_j_k] +part2Split_index;
	
				}
	
	
				//! halfen the weight of each particle
				part2Split->weight *= 0.5;
	
	
				//! SPLIT PARTICLE
				//! Store particle on temporary buffer
				particle part_buffer[2];
	
				memcpy(part_buffer +0, part2Split, sizeof(particle));
				memcpy(part_buffer +1, part2Split, sizeof(particle));
	
	
				memcpy(part_buffer[0].rel_r, x_new1, 3*sizeof(PARTICLE_REAL));
				memcpy(part_buffer[1].rel_r, x_new2, 3*sizeof(PARTICLE_REAL));
	
	
#ifdef SORT_BY_WEIGHT

				//! ----- DELETE part2Split: ----------------

				//! 1) shift particle_list over part2Split
				memmove(pArray[species][i_j_k]    +(part2Split_index),
					pArray[species][i_j_k]    +(part2Split_index+1),
					(num_MPiC[species][i_j_k] -(part2Split_index+1))*sizeof(particle));

				//! 2) decrement num_MPiC
				num_MPiC[species][i_j_k]--;
	
	
				//! Sort particle in List
				//! Estimate first particle in list which weight is larger
				//! than part2Split->weight

				//! Loop from List end to get small memmove of a lot weights equal

				//! Allow for negative index (case even zeroth particle weight is smaller)
				INT32 new_part_index;
				for(new_part_index=num_MPiC[species][i_j_k]-1; new_part_index>=0; new_part_index--)
				if(part_buffer[0].weight >= pArray[species][i_j_k][new_part_index].weight)
				break;
			
				//! new_part_index now is index of particle which weight is smaller,
				//! (or -1 in case new particle's weight is smallest in cell)
				//! so insert new particle behind estimated particle
				new_part_index++;
	
	
				//! shift particle_list to get space for 2 new particle.
				//! since memory overlaps, use memmove rather than memcpy.
				//! (memcpy may be undefined and IS NOT measurable faster)
				memmove(pArray[species][i_j_k]  +(new_part_index+2),
					pArray[species][i_j_k]  +(new_part_index),
					(num_MPiC[species][i_j_k] -new_part_index)*sizeof(particle));
	
				//! insert particle at position new_part_index
				memcpy(pArray[species][i_j_k]  +new_part_index,
				part_buffer,
				2*sizeof(particle));
	
#else
				
	
				//! shift particle_list to get space for 1 new particle right behind p2split
				memmove(pArray[species][i_j_k]  +(part2Split_index+2),
					pArray[species][i_j_k]  +(part2Split_index+1),
					(num_MPiC[species][i_j_k] -(part2Split_index+1))*sizeof(particle));
	
				//! insert both particles at position part2Split_index
				memcpy(pArray[species][i_j_k]  +part2Split_index,
				part_buffer,
				2*sizeof(particle));
	
				//! subtract one, add 2 see below
				num_MPiC[species][i_j_k]--;
	
	
#endif
	
			
				//! update statistics
				num_split_in_L[RLevel]++;
				num_MPiC[species][i_j_k]+=2;
				num_total_particles_in_L[RLevel]++;
				num_total_particles++;
	
	
	
			}//! end if in cell
			else
			num_split_canceled_out_of_cell++;

       		}//! end if canceled by weight

       }//! end if criteria fulfilled
      }//! end for species
     }//! end foor Box


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}


//!-----------------------------------------------------------
//! split_1to6particle
//!-----------------------------------------------------------
void CBlock::split_1to6particle(void)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();



	//    log_file << startWeight_of_smallestCell << "                     ";
	
	#ifndef SORT_BY_WEIGHT
	INT32 part2Split_index=0;
	#endif
	
	particle *part2Split;
	
	INT32 i,j,k, i_j_k;
	
	
	
	PARTICLE_REAL x_new[6][3];

	//! Instead of spliting 1->8 particle, it
	//! turned out to better splitting 1->2 and repeat 
	//! this 8 times.
	//! Doing so, every time the heaviest remaining particle 
	//! is newly estimated and split, which leads to better 
	//! weight balancing.



	for(i = 1; i < BlkNds_X-1; i++)
	 for(j = 1; j < BlkNds_Y-1; j++)
	  for(k = 1; k < BlkNds_Z-1; k++)
	  {
	
		i_j_k  = i*BlkNds_Y*BlkNds_Z 
			+j*BlkNds_Z 
			+k;
	
	
	
	
	      for(short species=0; species<num_Charged_Species; species++)
	       if(do_split_in_species[species])
	       {
	
		 for(INT32 split_cycle=0; split_cycle<num_split_each_cell; split_cycle++)
		  if(       num_MPiC[species][i_j_k] > 2
			&&  1.*num_MPiC[species][i_j_k] < Blk_optimal_MPiC[species] *(0.7 - 0.5*gsl_rng_uniform(randGen_general_asynchronous)))
		    {
	
	
	
			
			//!-----------------------------------------------------------
			//!1) Search step:
			//!   Find Particle with maximal weight
			//!-----------------------------------------------------------
#ifdef SORT_BY_WEIGHT
			//! heaviest particle is last particle of pArray
			part2Split = pArray[species][i_j_k] +num_MPiC[species][i_j_k]-1;
#else
	
	
			WEIGHT_REAL max_Weight = -1.;
		
			//! estimate heaviest particle in List
			for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
			{
		
				part2Split = pArray[species][i_j_k] +part_index;
		
				if( part2Split->weight > max_Weight)
				{
					max_Weight = part2Split->weight;
					part2Split_index = part_index;
				}
		
		
			}
		
			part2Split = pArray[species][i_j_k] +part2Split_index;
	
#endif
	
		
			
			//!-----------------------------------------------------------
			//!2) Split step:
			//!   Split Particle into 2
			//!-----------------------------------------------------------
			//! check if particle weight to small to split
			if(part2Split->weight < 1.e-7*startWeight_of_smallestCell)
			{
				split_cycle = num_split_each_cell;
				num_split_canceled_low_weight++;
			}
			else
			{
		
	
			int max_try = 1;
	
			bool coords_in_cell = false;
	
			//! Estimate new coords end ckeck if they are in cell
			for(INT32 num_try = 0; num_try<max_try; num_try++)
			{
	
	// 			INT32 fac = 1.*(max_try-num_try)/max_try;
	// 			PARTICLE_REAL shift = 0.1 +fac*gsl_rng_uniform(randGen_general_asynchronous);
	// 			PARTICLE_REAL shift = 0.2 +random;
	
				PARTICLE_REAL shift = 0.;//1./sqrt(num_MPiC[species][i_j_k]);
	
				for(INT32 part_ind=0; part_ind<6; part_ind++)
				memcpy(x_new[part_ind],
				part2Split->rel_r,
				3*sizeof(PARTICLE_REAL));
				
				shift=x_new[0][0];
	
				if(x_new[0][1] < shift) shift = x_new[0][1];
				if(x_new[0][2] < shift) shift = x_new[0][2];
	
				if( (1.-x_new[0][0]) < shift) shift = 1.-x_new[0][0];
				if( (1.-x_new[0][1]) < shift) shift = 1.-x_new[0][1];
				if( (1.-x_new[0][2]) < shift) shift = 1.-x_new[0][2];
	
				shift -=1.e-5;
				shift = 1.*gsl_rng_uniform(randGen_general_asynchronous);
				if(shift<0.02)
				break;
	
				x_new[0][0] -= shift;
				x_new[1][0] += shift;
	
				x_new[2][1] -= shift;
				x_new[3][1] += shift;
	
				x_new[4][2] -= shift;
				x_new[5][2] += shift;
	
				if(    x_new[0][0]>=0. && x_new[1][0]<1.
				&& x_new[2][1]>=0. && x_new[3][1]<1.
				&& x_new[4][2]>=0. && x_new[5][2]<1.)
				coords_in_cell = true;
	
				if(coords_in_cell)
				break;
			}
	
	
			//! check whether new cords are in old cell
			if(coords_in_cell)
			{
	
			
				//! alloc further memory if required
				if(size_of_pArray[species][i_j_k] < num_MPiC[species][i_j_k]+5)
				{
	
					//! alloc new pArray of larger size
					resize_pArray(species, i_j_k, int(min_pArray_Size_factor
									*(num_MPiC[species][i_j_k]+5)));
	
					//! adress has changed -> reset pointer !!!
					part2Split = pArray[species][i_j_k] +num_MPiC[species][i_j_k]-1;
	
				}
	
	
				//! half the weight for each particle
				part2Split->weight *= 1./6.;
	
	
				//! SPLIT PARTICLE
				//! Store particle on temporary buffer
				particle part_buffer[6];
	
				for(INT32 part_ind=0; part_ind<6; part_ind++)
				{
	
					memcpy(part_buffer +part_ind,   part2Split, sizeof(particle));
					memcpy(part_buffer[part_ind].rel_r, x_new[part_ind], 3*sizeof(PARTICLE_REAL));
	
				}
	
	
	
#ifdef SORT_BY_WEIGHT
	
				//! DELETE part2Split:
				//! part2Split is last of pArray.
				//! delete particle by decrementing num_MPiC
				num_MPiC[species][i_j_k]--;
	
	
				//! Sort particle in List
				//! Estimate first particle in list which weight is larger
				//! than part2Split->weight
				//! Loop from List end, start at position num_MPiC-2 to
				//! ignore part2Split which is at position num_MPiC-1
				//! Allow for temporary negative index (case list is empty)
				INT32 new_part_index;
				for(new_part_index=num_MPiC[species][i_j_k]-1; new_part_index>=0; new_part_index--)
				if(part2Split->weight >= pArray[species][i_j_k][new_part_index].weight)
				break;
		
				//! insert new particle behind estimated particle
				new_part_index++;
	
	
				//! shift particle_list to get space for 2 new particle.
				//! since memory overlaps, use memmove rather than memcpy.
				//! (memcpy may be undefined and IS NOT measurable faster)
				memmove(pArray[species][i_j_k]  +(new_part_index+6),
					pArray[species][i_j_k]  +(new_part_index),
					(num_MPiC[species][i_j_k] -new_part_index)*sizeof(particle));
	
				//! insert particle at position new_part_index
				memcpy(pArray[species][i_j_k]  +new_part_index,
				part_buffer,
				6*sizeof(particle));
	
#else
	
				log_file << " ACTIVATE SORT BY WEIGHT FOR 1->6 SPLIT" << endl << "Exiting ... " << endl;
				exit(0);
				//! shift particle_list to get space for 1 new particle right behind p2split
				memmove(pArray[species][i_j_k]  +(part2Split_index+2),
					pArray[species][i_j_k]  +(part2Split_index+1),
					(num_MPiC[species][i_j_k] -(part2Split_index+1))*sizeof(particle));
	
				//! insert both particles at position part2Split_index
				memcpy(pArray[species][i_j_k]  +part2Split_index,
				part_buffer,
				2*sizeof(particle));
	
				//! subtract one, add 2 as below
				num_MPiC[species][i_j_k]--;
	
	
#endif
	
			
				//! update statistics
				num_split_in_L[RLevel]++;
				num_MPiC[species][i_j_k]+=6;
				num_total_particles_in_L[RLevel]+=5;
				num_total_particles+=5;
	
	
	
			}//! end if in cell
			else
			num_split_canceled_out_of_cell++;
	
		}//! end if canceled by weight
	
	  }//! end if criteria fulfilled
	 }//! end for species
	}//! end foor Box


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}



//!-----------------------------------------------------------
//! merge_particle
//! NOTE:
//! This function gets awfully slow when to many particles are
//! within one cell -> N^2 performance 
//! Merging not often enough will result in many particles in 
//! one cell near level borders. This in turn leads to extremly long
//! merging computaion (finding 3 particles)
//!-----------------------------------------------------------
void CBlock::merge_particle(void)
{



	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


//! nur bei angeschalteter Magnetfeld routine "Streifen Muster", sonst nicht
//! -> Warum ???

//! "lacks" an gitter grenzen bei zu starken merging
//! -> auch ohne splitting !
//! -> gather from parent nicht ausreichend Particle


  INT32 num_merge_this_cell = 0;

  INT32 num_Compare;


  INT32 Small_Ekin = 0;
  INT32 Small_abs_L = 0;
  INT32 Small_abs_Plane = 0;
  INT32 Small_Plane_L = 0;
  D_REAL rho_of_cell = 0;
  D_REAL* rho_species;

  INT32 i,j,k, ip1_jp1_kp1;
  INT32 i_j_k, ip1_j_k, i_jp1_k ,i_j_kp1;
  INT32 ip1_j_kp1, ip1_jp1_k, i_jp1_kp1;

  particle *MergeList;

  //! indices for loops
  int a,b;
  int Part_Ind[3] = {-1,-1,-1};

  //! NOTE: weight related variables must be of MERGE_REAL precision
  WEIGHT_REAL deviation = 0.;
  WEIGHT_REAL min_deviation = 1.e18;
  WEIGHT_REAL weight[3],wA,wB,weight_total;
  WEIGHT_REAL lower_limit = 1.e-14;


  //! indices 0,1,2 are indices of initial particles
  //! indices A,B   are indices of merged  particles
  PARTICLE_REAL v0[3],v1[3],v2[3];
  PARTICLE_REAL x0[3],x1[3],x2[3];


  PARTICLE_REAL v_diff[3],x_diff[3];
  PARTICLE_REAL xA[3],xB[3], vA[3],vB[3];
  PARTICLE_REAL temp0[3],temp1[3],temp2[3],n_Plane[3],abs_Plane,abs_L,temp;

  PARTICLE_REAL E_kin_CM,v_magn,v_CM[3],r_CM[3],P_CM[3],L_CM[3],normed_L_CM[3];



   for(i = 1; i < BlkNds_X-1; i++)
    for(j = 1; j < BlkNds_Y-1; j++)
     for(k = 1; k < BlkNds_Z-1; k++)
     {

	i_j_k  = i*BlkNds_Y*BlkNds_Z 
		+j*BlkNds_Z 
		+k;



	//! -----------------------------------------
	ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
	i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
	i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
	
	//! ------------------------------------------
	ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
	ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
	i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
	ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
	

       for(short species=0; species<num_Charged_Species; species++)
        if(do_merge_in_species[species])
        {


	 //! use individual species density as reference density
	 rho_species = Field_Type[id_allRhoSpecies] +species *num_nodes_in_block;

	 //! density of cell i_j_k is average of 8 aligned nodes:
	 rho_of_cell =  0.125 *(  rho_species[i_j_k  ]
				 +rho_species[ip1_j_k]
				 +rho_species[i_jp1_k]
				 +rho_species[i_j_kp1]
		
				 +rho_species[ip1_jp1_k  ]
				 +rho_species[ip1_j_kp1  ]
				 +rho_species[i_jp1_kp1  ]
				 +rho_species[ip1_jp1_kp1]);


	//! In case of huge particle numbers repeat merge more often 
	if(num_MPiC[species][i_j_k] > 4*Blk_optimal_MPiC[species])
	num_merge_this_cell = (num_MPiC[species][i_j_k]/Blk_optimal_MPiC[species])*num_merge_each_cell;
	else
        num_merge_this_cell = num_merge_each_cell;
        
        //! Calc position of cell centre
        PARTICLE_REAL x_part[3],r_normed[3];
        INT32 cell_indices[3];
        
        x_part[0]=0.;
        x_part[1]=0.;
        x_part[2]=0.;
        
        cell_indices[0] = i;
        cell_indices[1] = j;
        cell_indices[2] = k;
        intern2normedCoords(r_normed, x_part, cell_indices);


	//! not to much merging eg. at shocks
        for(INT32 merge_cycle=0; merge_cycle<num_merge_this_cell; merge_cycle++)
         if(    //! 1) only merge in case num_MPiC>2 (3 particles are required for merging)
         	   num_MPiC[species][i_j_k] > 2 
		//! 2) never merge beyond oMPiC (else conflict with Split)
	     && 1.*num_MPiC[species][i_j_k] > Blk_optimal_MPiC[species]* (1.5 + 1.*gsl_rng_uniform(randGen_general_asynchronous))
		//! 3) - multiply optimal_MPiC with rho in order to avoid merging e.g. at shocks 
		//!    - fluctuations in Part number by factor of 2 are likely to occur, so
		//!      never merge unless PiC exeeds oPiC by at least a factor of 2 !!!
	     && (1.*num_MPiC[species][i_j_k] > Blk_optimal_MPiC[species]*rho_of_cell* (1.5 + 1.*gsl_rng_uniform(randGen_general_asynchronous))
             || r_normed[0] > merge_tail_distance + (1.*gsl_rng_uniform(randGen_general_asynchronous)))
        )
         {



		//! never try to compare more particle than in List

		//! Use const value specified in parameter.cpp 
		//! value should be similar to optimal_MPiC
		//! use sqrt of num_MPiC[species][i_j_k]
		//! doesnt perform significantly faster
		if(num_MPiC[species][i_j_k] < num_particle_in_MergeList)
		num_Compare = num_MPiC[species][i_j_k];
		else
		num_Compare = num_particle_in_MergeList;

                
		//! --------------------------------------------------------------
		//! Estimate triple of particle:
		//! First find pair of optimal fitting particles in velocity space
		//! --------------------------------------------------------------

		//! use MergeList as abbreviation for pArray
		MergeList = pArray[species][i_j_k];

		bool ind_set = false;

		min_deviation = 1.e18;
		for(a=0; a<num_Compare-1; a++)
		 for(b=a+1; b<num_Compare; b++)
		 {


		   v0[0] = MergeList[a].v[0] - MergeList[b].v[0];
		   v0[1] = MergeList[a].v[1] - MergeList[b].v[1];
		   v0[2] = MergeList[a].v[2] - MergeList[b].v[2];


		   deviation = (MergeList[a].weight + MergeList[b].weight) *vec_len2(v0);

		   if(deviation < min_deviation)
		   {
			min_deviation = deviation;
			Part_Ind[0] = a;
			Part_Ind[1] = b;
			ind_set = true;

		   }

		 }

#ifdef DEBUG_MODUS
		if(!ind_set)
		{
	
			log_file <<"ERROR in Merge:" << endl;
			log_file << "First comparison: indices have not been set. " << endl;

			log_file << "Exiting ... " << endl;
			exit(1);
		}
#endif

		//! now the first optimal pair is estimated, calculate v_CM
		weight[0]    = MergeList[Part_Ind[0]].weight;
		weight[1]    = MergeList[Part_Ind[1]].weight;
		weight_total = weight[0] +weight[1];

		for(a=0; a<3;a++)
		v_CM[a] = ( weight[0] *MergeList[Part_Ind[0]].v[a]
			   +weight[1] *MergeList[Part_Ind[1]].v[a]
			  )/weight_total;

		//! --------------------------------------------------------------
		//! Estimate a third particle matching optimal to 
		//! already estimateted pair
		//! --------------------------------------------------------------

		//! estimate a third particle matching optimal to already estimateted pair
		min_deviation = 1.e18;


		ind_set = false;
		for(a=0; a<num_Compare;a++)
		{

		   if(a!=Part_Ind[0] && a!=Part_Ind[1])
		   {
			v0[0] = MergeList[a].v[0]-v_CM[0];
			v0[1] = MergeList[a].v[1]-v_CM[1];
			v0[2] = MergeList[a].v[2]-v_CM[2];
			
			weight_total = MergeList[a].weight + weight[0] +weight[1];

			deviation = weight_total*vec_len2(v0);


			if(deviation < min_deviation)
		   	{
				min_deviation = deviation;
				Part_Ind[2] = a;
				ind_set = true;

		   	}

		   }

		}

#ifdef DEBUG_MODUS
		if(!ind_set)
		{
	
			log_file <<"ERROR in Merge:" << endl;
			log_file << "Second comparison: indices have not been set. " << endl;

			log_file << "Exiting ... " << endl;
			exit(1);
		}
#endif


		weight[2] = MergeList[Part_Ind[2]].weight;

// 		log_file << " (" << Part_Ind[0] <<","<<Part_Ind[1]<<","<<Part_Ind[2]<<")   ";
		//! ********************************************************
		//! 2) MERGE PARTICLES
		//! Now three Particle fitting optimal together have
		//! been estimated and indices are stored in Part_Ind[0,1,2]
		//! ********************************************************



		//! store positions to x[comp] to make code more readable
		memcpy(x0, MergeList[Part_Ind[0]].rel_r, 3*sizeof(PARTICLE_REAL));
		memcpy(x1, MergeList[Part_Ind[1]].rel_r, 3*sizeof(PARTICLE_REAL));
		memcpy(x2, MergeList[Part_Ind[2]].rel_r, 3*sizeof(PARTICLE_REAL));

		//! store velocities to v[comp] to make code more readable
		memcpy(v0, MergeList[Part_Ind[0]].v, 3*sizeof(PARTICLE_REAL));
		memcpy(v1, MergeList[Part_Ind[1]].v, 3*sizeof(PARTICLE_REAL));
		memcpy(v2, MergeList[Part_Ind[2]].v, 3*sizeof(PARTICLE_REAL));




		//! calculate center of mass, its momentum and velocity
		weight_total = weight[0] +weight[1] +weight[2];
		//! use half weight total for each new particle
		wA = wB = 0.5 * weight_total;

		for(a=0; a<3; a++)
		{
			r_CM[a] = (weight[0]*x0[a] +weight[1]*x1[a] +weight[2]*x2[a])/weight_total;
			P_CM[a] =  weight[0]*v0[a] +weight[1]*v1[a] +weight[2]*v2[a];

			//! velocity of centre of mass
			v_CM[a] = P_CM[a]/weight_total;
		}


#ifdef DEBUG_MODUS
		//! Check for to Large centre of velocity
		//! NOTE: Using the fprint c syntax results in errors 
		//!	    when using long doubles !!!
		if(vec_scalar( v_CM,v_CM) > 1.e4)
		log_file << "Large v_CM in cell ("<<i<<","<<j<<","<<k<<") "
		     << "with v_CM ("<<v_CM[0]<<","<<v_CM[1]<<","<<v_CM[2]<<") "
		     << "with P_CM ("<<P_CM[0]<<","<<P_CM[1]<<","<<P_CM[2]<<") "
		     << "with weight_total " << weight_total << endl;


		//! Check for to Large centre of mass
		if(vec_scalar( r_CM,r_CM) > 3.)
		log_file << "Large r_CM in cell ("<<i<<","<<j<<","<<k<<") "
		     << "with r_CM ("<<r_CM[0]<<","<<r_CM[1]<<","<<r_CM[2]<<") "
		     << "with P_CM ("<<P_CM[0]<<","<<P_CM[1]<<","<<P_CM[2]<<") "
		     << "with weight_total " << weight_total << endl;
#endif




		//! transform in center of mass system
		for(a=0; a<3; a++)
		{
			x0[a] -= r_CM[a];
			x1[a] -= r_CM[a];
			x2[a] -= r_CM[a];

			v0[a] -= v_CM[a];
			v1[a] -= v_CM[a];
			v2[a] -= v_CM[a];
		}



		//! Now velocities and positions are transformed in CM system
		//! Calculate new velocities of particle A and B
		E_kin_CM = 0.5*( weight[0]*vec_scalar(v0,v0)
			        +weight[1]*vec_scalar(v1,v1)
			        +weight[2]*vec_scalar(v2,v2));

		//! Angular Momentum
		vec_cross(temp0,x0,v0);
		vec_cross(temp1,x1,v1);
		vec_cross(temp2,x2,v2);

		for(a=0; a<3; a++)
		L_CM[a] = weight[0]*temp0[a] +weight[1]*temp1[a] +weight[2]*temp2[a];

		for(a=0; a<3; a++)
		{
			//! calcultion of plane which is defined by the
			//! positionn of all three particles
			temp0[a] = (x1[a] - x0[a]);
			temp1[a] = (x2[a] - x0[a]);

		}

		vec_cross(n_Plane,temp0,temp1);


		abs_Plane = vec_len(n_Plane);
		abs_L = vec_len(L_CM);




		if(E_kin_CM < lower_limit || abs_L < lower_limit)
		{
			//! TODO:
			//! HIER STIMMT WAS NICHT !!!
			//! IM KALTEN PLASMA ENTSTEHEH STREIFEN PARALLEL 
			//! ZU STROEMUNGSRICHTUNG 
			//! Above means kin Energy with respect to centre of mass
			//! is close to zero which means that velocities of particles
			//! are close to zero in CM-System, so just use v_CM for new particle.
			//! In order to get differing particle trajectories, the new 
			//! positions have to be different. So define arbitrary r_diff
			//! perp to v.

			//! 1. Case:
			//! Old Angular momentum L = r x mv is close to zero as |E| = 0.
			//! New L will be of course also zero in CM-System (as v_diff = 0).
			if(E_kin_CM  < lower_limit)
			{
				Small_Ekin++;
	
				for(a=0; a<3; a++)
				{
					vA[a] = v_CM[a];
					vB[a] = v_CM[a];
				}


				temp1[0] = 1. *gsl_rng_uniform(randGen_general_asynchronous);
				temp1[1] = 1. *gsl_rng_uniform(randGen_general_asynchronous);
				temp1[2] = 1. *gsl_rng_uniform(randGen_general_asynchronous);

				//! in case velocity is zero (eg. newborn heavy ions)
				//! add some offset, oth. |x_diff| equals zero
                                temp2[0] = v_CM[0]+1.e-5;
				temp2[1] = v_CM[1];
                                temp2[2] = v_CM[2];

				vec_cross(x_diff, temp1, temp2);

				temp = 0.25/ vec_len(x_diff);


#ifdef DEBUG_MODUS
				//! Check for to small x_diff
				if(temp > 1.e10)
				{
				
				log_file << "E_kin_CM -> Small x_diff in cell ("<<i<<","<<j<<","<<k<<") "
				     << "with v_CM ("<<v_CM[0]<<","<<v_CM[1]<<","<<v_CM[2]<<") "
				     << "with x_diff ("<<x_diff[0]<<","<<x_diff[1]<<","<<x_diff[2]<<") "
				     << "with weight_total " << weight_total << endl;

				log_file << " Part A v: ("
				     << MergeList[Part_Ind[0]].v[0] <<","
				     << MergeList[Part_Ind[0]].v[1] <<","
				     << MergeList[Part_Ind[0]].v[2] <<")" << endl;

				log_file << " Part B v: ("
				     << MergeList[Part_Ind[1]].v[0] <<","
				     << MergeList[Part_Ind[1]].v[1] <<","
				     << MergeList[Part_Ind[1]].v[2] <<")" << endl;

				log_file << " Part C v: ("
				     << MergeList[Part_Ind[2]].v[0] <<","
				     << MergeList[Part_Ind[2]].v[1] <<","
				     << MergeList[Part_Ind[2]].v[2] <<")" << endl;
				 

				}
#endif



				for(a=0; a<3; a++)
				{
					xA[a] =  r_CM[a] - x_diff[a]*temp;
					xB[a] =  r_CM[a] + x_diff[a]*temp;
				}

#ifdef DEBUG_MODUS
                                //! Check for to large x_diff
                                if( vec_len(x_diff) > 400.)
				    log_file << "E_kin_CM -> Large x_diff in cell ("<<i<<","<<j<<","<<k<<") "
					 << "with v_CM ("<<v_CM[0]<<","<<v_CM[1]<<","<<v_CM[2]<<") "
					 << "with x_diff ("<<x_diff[0]<<","<<x_diff[1]<<","<<x_diff[2]<<") "
					 << "with weight_total " << weight_total << endl;
#endif


			}
			//! in case |L| is close to zero either 
			//! 1) velocities are zero 
			//! 2) positions are at centre of mass 
			//! Condition 1) is already catched so condidition 
			//! 2) will be the case.

			//! In order to conserve angular momentum, positions are 
			//! set to CM-System (-> L = r x v  -r x v = 0).
			//! Plane normal |n| is used to define a v_diff.
			else if(abs_L  < lower_limit && abs_Plane > lower_limit) 
			{

				Small_abs_L++;
				v_magn = sqrt(E_kin_CM/wA);

				for(a=0; a<3; a++)
				{

					v_diff[a] = v_magn* (n_Plane[a]/abs_Plane);

					vA[a] = v_CM[a]-v_diff[a];
					vB[a] = v_CM[a]+v_diff[a];

					xA[a] =  r_CM[a];
					xB[a] =  r_CM[a];
				}


			}
			//! This will be a really rare case.
			//! For sake of simplicity just use 
			//! positions and velocities of old particle 0 & 1
			else if(abs_Plane < lower_limit && abs_L  < lower_limit) 
			{

				Small_Plane_L++;

				for(a=0; a<3; a++)
				{
					vA[a] = v_CM[a]+v0[a];
					vB[a] = v_CM[a]+v1[a];

					xA[a] =  r_CM[a]+x0[a];
					xB[a] =  r_CM[a]+x1[a];
				}


			}
			



		}
		else
		{

			//! Below means particle do not define a Plane
			//! e.g. they are at one point, one line.
			//! So just set n parallel to L and proceed
			//! as usual
			if(abs_Plane < lower_limit) 
			{
				Small_abs_Plane++;

				abs_Plane = abs_L;
				for(a=0; a<3; a++)
				n_Plane[a] = L_CM[a];
			}


			for(a=0; a<3; a++)
			{
				n_Plane[a]/=abs_Plane;
				normed_L_CM[a] = L_CM[a]/abs_L;
			}

			v_magn = sqrt(E_kin_CM/wA);
			temp = vec_scalar(normed_L_CM,n_Plane);

			for(a=0; a<3; a++)
			v_diff[a] = v_magn*(n_Plane[a] - temp*normed_L_CM[a]);

			if(vec_scalar(n_Plane,v_diff) < 0.)
			for(a=0; a<3; a++) v_diff[a] *= -1.;

			for(a=0; a<3; a++)
			{
				vA[a] = v_CM[a] - v_diff[a];
				vB[a] = v_CM[a] + v_diff[a];
			}


			vec_cross(x_diff,v_diff,L_CM);


			for(a=0; a<3; a++)
			x_diff[a] /= (2.*E_kin_CM);


	 		temp = vec_scalar(n_Plane,x_diff) / vec_scalar(n_Plane,v_diff);

			for(a=0; a<3; a++)
			{
  				x_diff[a] -= v_diff[a] * temp;

				xA[a] = r_CM[a] - x_diff[a];
				xB[a] = r_CM[a] + x_diff[a];
			}

#ifdef DEBUG_MODUS
			//! Check for to large x_diff
			if(vec_scalar( v_CM,v_CM) > 1.e4 || vec_scalar( x_diff,x_diff) > 200.)
			log_file << "STD -> Large v_diff/x_diff in cell ("<<i<<","<<j<<","<<k<<") "
			     << "with v_diff ("<<v_diff[0]<<","<<v_diff[1]<<","<<v_diff[2]<<") "
			     << "with x_diff ("<<x_diff[0]<<","<<x_diff[1]<<","<<x_diff[2]<<") "
			     << "with weight_total " << weight_total << endl;
#endif

		}
			

		//! allow to set particle somewhere on plane without conserving CoM
		//! -> in order to avoid "stripes" in moments
		if(!conserve_CoM_at_merging)
		{


			//! copy centre of mass position to new positions
			memcpy(xA, r_CM, 3*sizeof(PARTICLE_REAL));
			memcpy(xB, r_CM, 3*sizeof(PARTICLE_REAL));

			//! randomize position for first particle
			for(INT32 num_try = 0; num_try<num_randomize_position_merge; num_try++)
			{
				set_Cords_perp_to_velocity_CoM_OFF(xA, vA);

				if(is_in_cell(xA))
				break;
			}

			//! randomize position for second particle, if first is in cell
			if(is_in_cell(xA))
			 for(INT32 num_try = 0; num_try<num_randomize_position_merge; num_try++)
			 {
	
				set_Cords_perp_to_velocity_CoM_OFF(xB, vB);

				if(is_in_cell(xB))
				break;

			 }

		}



		//! --------------------------------------------------
		//! Check whether new positions are out of active cell
		//! Case yes, insert particle
		//! --------------------------------------------------
		if(are_in_cell(xA, xB))
		{



			//! ------------------------------------------
			//! ------ Delete Particle 1 -----------------
			//! ------------------------------------------
			//! shift mergelist to cover particle 1.
			//! since memory overlaps, use memmove rather than memcpy.
			//! (memcpy may be undefined and IS NOT measurable faster)
			memmove(MergeList +Part_Ind[0],
				MergeList +Part_Ind[0]+1,
				(num_MPiC[species][i_j_k] -(Part_Ind[0]+1))*sizeof(particle));


			//! adjust residual indices
			if(Part_Ind[1] > Part_Ind[0]) Part_Ind[1]--;
			if(Part_Ind[2] > Part_Ind[0]) Part_Ind[2]--;
			num_MPiC[species][i_j_k]--;


			//! ------------------------------------------
			//! ------ Delete Particle 2 -----------------
			//! ------------------------------------------
	
			//! shift mergelist to cover particle 2
			memmove(MergeList +Part_Ind[1],
				MergeList +Part_Ind[1]+1,
				(num_MPiC[species][i_j_k] -(Part_Ind[1]+1))*sizeof(particle));


			//! adjust residual indices
			if(Part_Ind[2] > Part_Ind[1]) Part_Ind[2]--;
			num_MPiC[species][i_j_k]--;


			//! ------------------------------------------
			//! ------ Delete Particle 3 -----------------
			//! ------------------------------------------
			//! shift mergelist to cover particle 3
			memmove(MergeList +Part_Ind[2],
				MergeList +Part_Ind[2]+1,
				(num_MPiC[species][i_j_k] -(Part_Ind[2]+1))*sizeof(particle));

			//! adjust residual indices
			num_MPiC[species][i_j_k]--;




			//! Sort particle in List
			//! Estimate first particle in list which weight is larger
			//! than wA
			//! Loop from List end
			//! Allow for temporary negative index (case list is empty)
			INT32 new_part_index;
			for(new_part_index=num_MPiC[species][i_j_k]-1; new_part_index>=0; new_part_index--)
			if(wA >= MergeList[new_part_index].weight)
			break;
	
			//! insert new particle behind estimated particle
			new_part_index++;


			//! shift particle_list to create space for 2 new particle.
			//! since memory overlaps, use memmove rather than memcpy.
			//! (memcpy may be undefined and IS NOT measurable faster)
			memmove(MergeList +(new_part_index+2),
				MergeList +(new_part_index),
				(num_MPiC[species][i_j_k] -new_part_index)*sizeof(particle));

			//! insert paticle at position new_part_index and new_part_index+1
			MergeList[new_part_index +0].weight = wA;
			MergeList[new_part_index +1].weight = wB;

			memcpy(MergeList[new_part_index +0].v,     vA, 3*sizeof(PARTICLE_REAL));
			memcpy(MergeList[new_part_index +1].v,     vB, 3*sizeof(PARTICLE_REAL));

			memcpy(MergeList[new_part_index +0].rel_r, xA, 3*sizeof(PARTICLE_REAL));
			memcpy(MergeList[new_part_index +1].rel_r, xB, 3*sizeof(PARTICLE_REAL));


			//! update statistics
			num_MPiC[species][i_j_k]+=2;
			num_total_particles--;
			num_total_particles_in_L[RLevel]--;


			//! --------------------------------------------------
			num_merge_in_L[RLevel]++;

		}
		else
		{
			num_merge_canceled_out_of_cell++;

			//! further retry when CoM conservation is set off
			if(conserve_CoM_at_merging)
			merge_cycle = num_merge_this_cell;
		}



	}//! end if merge particles ?
      }//! for species
     }//! end for all cells



	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


}

/*
//!-----------------------------------------------------------
//! save_alloc_particle
//!-----------------------------------------------------------
inline particle* save_alloc_particle(void)
{


    	particle *p;

	try
	{
		p = new particle;
	}
	catch (exception& e)
	{
// 		log_file << "Standard exception: " << e.what() << endl;
		bad_allocs++;
		p = 0;
	}
	return p;

}*/
