



#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

#include "CBlk.h"
#include "defines.h"
#include "parameters.h"
#include "absolute_Globals.h"
#include "utils.h"



using namespace std;



extern PARTICLE_REAL fastest_particle_v2[];


extern D_REAL *CellVol_of_L;
extern D_REAL *dt_particle_of_L, **delta_of_L;


const D_REAL squared_R_Obstacle = R_Obstacle*R_Obstacle;


//!--------------------------------------------------------------
//! delete_particle
//!--------------------------------------------------------------
void CBlock::delete_particle(INT32 species, INT32 i_j_k, INT32 &part_index)
{


	//! DELETE ACTIVE PARTICLE
	//! shift next up to last particle of list on top active_particle
	memmove(pArray[species][i_j_k] +(part_index),
		pArray[species][i_j_k] +(part_index+1),
		(num_MPiC[species][i_j_k] -(part_index+1))*sizeof(particle));
	
	//! part_index=-1 may temporary occcur here
	part_index--;
	
	//! update statistics
	num_total_particles--;
	num_total_particles_in_L[RLevel]--;
	num_MPiC[species][i_j_k]--;
	


}

//!--------------------------------------------------------------
//! block_intern_move_particle
//!--------------------------------------------------------------
bool CBlock::block_intern_move_particle(void)
{


	//! 1.)
	//! decide in advance, whether it is required to
	//! check for boundaries in move
	if(is_box_boundary[_ANY_] || is_box_boundary[_OBS_])
	{

		//! TODO:
		//! 1) distinguish inside move method ordinary, boundary & obstacle cell
		//!    using bitmask in Flag
		//! 2) combine flag and eta_flag by using bitmasks
		for(INT32 oct=0; oct<8; oct++)
		//!NOTE: below if statement will ignore moving of injected particle
		//!	   when boundary blocks are refined !!!
		//!	   update -> happens anyway !?!
		if(!child_array[oct])
		move_particle(oct);
	}
	else
	{
		for(INT32 oct=0; oct<8; oct++)
		//!NOTE: below if statement will ignore moving of injected particle
		//!	   when boundary blocks are refined !!!
		//!	   update -> happens anyway !?!
		if(!child_array[oct])
		move_particle_noObsBoxBoundCheck(oct);

	}



	//! 2.)
	if(!cell_reassign_particle())
	return false;


	//! 3.)
	//! TODO:
	//! THIS CAN CAUSE PROBLEMS IN CASE OF REFINEMENT:
	//! - L0 BLK BOUNDARY IS FILLED (-> PARTICLE IN 0 CELL)
	//! - BOX IS REFINED (-> ONLY PARTICLE INSIDE BLK ARE CONSIDERED)
	//! - SINCE L0 BLK IS COVERED BY L1 BLK, PARTICLE WILL NOT BE MOVED
	//! -> ERASE PARTICLE AT L0 BEFORE REFINEMENT
	//! -> FILL L1 BOUNDARY AFTER CREATION
	if(is_box_boundary[_ANY_])
	{


		for(INT32 oct=0; oct<8; oct++)
		if(!child_array[oct])
		inject_boundary_particle(oct);



	}

	return true;


}


//!--------------------------------------------------------------
//! move_particle:
//! NOTES:
//! - In case of higher temperature the move procedure gets 
//!   significant slower. This is certainly due to frequently 
//!   calling of the random function.
//! - In case unhook is called a particle can be inserted 
//!   in the new cells List sorted or unsorted. As this 
//!   sorting does not take much time, it can be always activated.
//!   (in case PPC = 20: t_sort < 1% of t_move)
//!--------------------------------------------------------------
void CBlock::move_particle(INT32 id_oct)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();



	INT32 new_i_j_k;
	INT32 index[3]; // allow for negativ numbers (which are deleted at once)
	
	PARTICLE_REAL r, x[3], temp_fastest_particle_v2;
	
	bool in_Obs = false;
	bool out_of_Box = false;
	
	
	
	//! temprary pointer to particle in list
	particle *active_particle;
	
	//! decide in init whether dt is different in indivudal levels
	PARTICLE_REAL dt_part = dt_particle_of_L[RLevel];
	
	
	//! con_base is needed as cell size for a particle (rel_r)
	//! ranges from 0 to 1 (= numerical units), wheras the physical
	//! cell size ranges from 0 to delta_cell (= physical units).
	//! As v is defined in physical units, con_base scales the dx 
	//! in moving step to numerical units.
	PARTICLE_REAL con_base[3] = {1./delta_of_L[RLevel][0],
				     1./delta_of_L[RLevel][1],
				     1./delta_of_L[RLevel][2]};




	//!	Some Particle remain unmoved for reassignement procedure
	//!	they could be moved by setting start to 0, but they are
	//!	not gathered anyway. Since they are only few (~1/1000000),
	//!	just leave it as before
	//!	(look xfig sketches for details)
	
	

	//! get index of oct (0 or 1 for each dimension)
	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);
	

	//! START:
	//! eg: BlkNds=6
	//! 1)
	//! a=0
	//! noBB
	//! -> start = 1
	//! 2)
	//! a=0
	//! isBB
	//! -> start = 0
	//! 3)
	//! a=1
	//! noBB
	//! -> start = 3
	//! 4)
	//! a=1
	//! isBB
	//! -> start = 3
	INT32 start[3] = {!a *!is_box_boundary[0] +a*BlkNds_X/2,
			 !b *!is_box_boundary[2] +b*BlkNds_Y/2,
			 !c *!is_box_boundary[4] +c*BlkNds_Z/2};


	//! END:
	//! eg: BlkNds=6
	//! 1)
	//! a=0
	//! noBB
	//! -> end = 3
	//! 2)
	//! a=0
	//! isBB
	//! -> end = 3
	//! 3)
	//! a=1
	//! noBB
	//! -> end = 5
	//! 4)
	//! a=1
	//! isBB
	//! -> end = 6
	INT32 end[3]={BlkNds_X -a*!is_box_boundary[1] -!a*(BlkNds_X/2),
		     BlkNds_Y -b*!is_box_boundary[3] -!b*(BlkNds_Y/2),
		     BlkNds_Z -c*!is_box_boundary[5] -!c*(BlkNds_Z/2)};

	



  //!-----------------------------------------------------------
  //!1) Move Step:
  //!   All particle weights are assumed to be > 0 now !!!
  //!-----------------------------------------------------------
  for(short species=0; species<num_Particle_Species; species++)
   for(INT32 i = start[0]; i < end[0]; i++)
    for(INT32 j = start[1]; j < end[1]; j++)
     for(INT32 k = start[2]; k < end[2]; k++)
     {

        INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
	        +j*BlkNds_Z 
	        +k;


	if(max_PiC<num_MPiC[species][i_j_k])
	max_PiC = num_MPiC[species][i_j_k];

      	for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
      	{


	 	active_particle = pArray[species][i_j_k] +part_index;


		//! ---- ADVANCE CARTESIEN POSITIONS ---------------------
		active_particle->rel_r[0] +=  dt_part * con_base[0] *active_particle->v[0];
		active_particle->rel_r[1] +=  dt_part * con_base[1] *active_particle->v[1];
		active_particle->rel_r[2] +=  dt_part * con_base[2] *active_particle->v[2];


		//! count particle as moved
		num_moved++;

{//! comments on interception below
		//! in order to avoid round off errors the interrogation below
		//! has to be performed in case PARTICLE_REAL is float.
		//! In case it is skipt, rel_r may get negative which may result
		//! in serious errors.'
		//! PART_REAL_PRECISION has to be chosen due to creteria below
		//! and is defined in defines.h

		//! To understand this, consider the following example:
		//! FLOAT Example:
		//! 1) Assume rel_r gets slighly negative which means particle left
		//!    cell at left hand side. New cell indices have to be calculated:
		//! 1a) 
		//! float rel_r  =    -1.e-7;
		//! float index  = 1. +rel_r;
		//! log_file << int(rel_r);
		//! log_file << int(index);
		//! OUTPUT:
		//! -1.e-7 ->CORRECT RESULT !!!!
		//! 0      ->CORRECT RESULT !!!!
		//! 	   -> paticle will leave cell and rel_r
		//! 	      will be reajuested


		//! 1b)
		//! float rel_r  =    -1.e-8;
		//! float index  = 1. +rel_r;
		//! log_file << int(rel_r);
		//! log_file << int(index);
		//! OUTPUT:
		//! -1.e-8 ->CORRECT RESULT !!!!
		//! 1      ->FALSE   RESULT !!!!
		//! 	   -> paticle will NOT leave cell and rel_r
		//! 	      will be remain negative

		//! In Example 1b) no cell crossing would be intercepted, hence
		//! rel_r remains negative which results in serious errors in
		//! gather and acceleration Method
		//! So, remember that float is only accurate up to -1.e-7
		//! The same happens for DOUBLE in case value below -1.e-16
		
		//! 
		//! 2a)
		//! float rel_r  = 1. -1.e-6;
		//! float index  = 1. +rel_r;
		//! log_file << int(rel_r);
		//! log_file << int(index);
		//! OUTPUT:
		//! 0.999999 ->CORRECT RESULT !!!!
		//! 1        ->CORRECT RESULT !!!!

		//! 
		//! 2b)
		//! float rel_r  = 1. -1.e-7;
		//! float index  = 1. +rel_r;
		//! log_file << int(rel_r);
		//! log_file << int(index);
		//! OUTPUT:
		//! 1  -> approx. CORRECT RESULT !!!!
		//! 1  -> FALSE RESULT !!!!

		//! 
		//! 2c)
		//! float rel_r  = 1. -1.e-8;
		//! float index  = 1. +rel_r;
		//! log_file << int(rel_r);
		//! log_file << int(index);
		//! OUTPUT:
		//! 1  -> approx. CORRECT RESULT !!!!
		//! 2  -> approx. CORRECT RESULT !!!!

		//! In Example 2b) no cell crossing would be intercepted, hence
		//! particle pos remains 1 and at the end of the cell.
		//! So it will be transported in next TL end result is
		//! in 1.e-6 accuriy correct.
		//! In Example 2c) particle will be transprted to next cell
		//! and position will be set to zero, which is correct in
		//! equal precision.
		//! -> no need to apply interrogation
}
//! end comments

		for(int comp=0; comp<3; comp++)
		{
			if(    active_particle->rel_r[comp] > 0. -PART_REAL_PRECISION
			    && active_particle->rel_r[comp] < 0. +PART_REAL_PRECISION)
			    active_particle->rel_r[comp] = 0.;
			else
			if(    active_particle->rel_r[comp] > 1. -PART_REAL_PRECISION
			    && active_particle->rel_r[comp] < 1. +PART_REAL_PRECISION)
			    active_particle->rel_r[comp] = 1.;

		}

#ifdef ESTIMATE_FASTEST_PARTICLE
		//! JUST FOR STATISTICS:
		//! FOR SOME REASON CAN SLOW SIGNIFICANTLY DOWN THE CODE !!!
		//! use squared velocity for comparison.
		//! only build sqrt for statistics (slightly faster)

		temp_fastest_particle_v2 = vec_len2(active_particle->v);

		if(temp_fastest_particle_v2 > fastest_particle_v2[RLevel])
		fastest_particle_v2[RLevel] = temp_fastest_particle_v2;
#endif



{//! comments to interception 
		//! NOTE: UPDATE:
		//! Interception below should not be necessay anymore when
		//! round off errors are avoided as above.

		//! NOTE:
		//! For some reason this "needless" if expression below has
		//! to be used, otw g++ compiled code fails. 
		//! Coompiling identical code with icc compiler this problem 
		//! does occur less frequently (or maybe it is not intercepted)

		//! !!! If rel_r = 1 than index must be i +1 (see above) !!!!
		//! For some reason this DOES NOT WORK RELYABLE (maybe round off) 
		//! So it may happen that rel_r = 1 and index = i
		//! -> particle is not inserted in next cell 
		//! -> rel_r[x] remains 1 after move method.
		//! -> This leads to int(2*rel_r[x]) = 2 in transport to Lp1
		//! -> This in turn leads to a Segmantation Fault.
		//! THIS IS NOT THE CASE IN DEBUG MODUS OR WHEN USING 
		//! DOUBLE PRECISION for PARTICLE_REAL !!!


#ifdef DEBUG_MODUS
		//! So the case below must not occur, however in rare cases
		//! it does and it HAS TO BE INTERCEPTED.
// 		 if(  (active_particle->rel_r[0]==1 && index[0] == i)
// 		   || (active_particle->rel_r[1]==1 && index[1] == j)
// 		   || (active_particle->rel_r[2]==1 && index[2] == k))
// 		{
// 
// 				log_file <<  "Particle rel_r: " 
// 				     << "(" <<active_particle->rel_r[0]<<" , "
// 				            <<active_particle->rel_r[1]<<" , "
// 				            <<active_particle->rel_r[2]<<") " << endl;
// 
// 
// 				log_file <<  "(i,j,k)  (" <<i<<","<<j<<","<<k<<") " << " <-> " ;
// 
// 				index[0] = i + int(active_particle->rel_r[0]);
// 				index[1] = j + int(active_particle->rel_r[1]);
// 				index[2] = k + int(active_particle->rel_r[2]);
// 
// 				log_file <<  "index[] (" <<index[0]<<","<<index[1]<<","
// 				     		      <<index[2]<<") " << endl;
// 				log_file << " corrected." << endl;
// 				exit(1);
// 
// 		}
#endif

//!---------------------- store entire phase function -------------------
#ifdef TL_PROTOCOL_ANY_PARTICLE
// 				if(TL_PROTOCOL_ANY_PARTICLE && TL%TL_PROTOCOL_ANY_PARTICLE==0)
// 				{
// 
// 
// 					INT32 cell_indices[3] = {i,j,k};
// 					intern2normedCoords(x, active_particle->rel_r, cell_indices);



// 					AnyParticle_FILE << part_counter<< "	";

					//! write weight
// 					AnyParticle_FILE << -active_particle->weight << "	";

					//! write position in normed units
// 					AnyParticle_FILE << x[0] << " ";
// 					AnyParticle_FILE << x[1] << " ";
// 					AnyParticle_FILE << x[2] << "	";


					//! write velocity in normed units
// 					AnyParticle_FILE << active_particle->v[0] << " ";
// 					AnyParticle_FILE << active_particle->v[1] << " ";
// 					AnyParticle_FILE << active_particle->v[2] << " ";

// 					AnyParticle_FILE << vec_len(active_particle->v);
// 					part_counter++;
// 	
// 					AnyParticle_FILE << endl;
// 				}
#endif
//!-----------------------------------------------------------------------
}
//!end comments



		//! reset out_of_Box
		in_Obs = false;
		out_of_Box = false;
		

		//! --------- OUT OF BOX ??? -------------------------------
		//! IF OUT OF BOX EREASE PARTICLE HERE TO AVOID REDUNDANT CALCULATIONS
		//! IN PARTICLE CELL REASSIGNEMENT
		//! check whether block is any box boundary
		if(is_box_boundary[_ANY_])
		{
			index[0] = int(1.*i + active_particle->rel_r[0]);
			index[1] = int(1.*j + active_particle->rel_r[1]);
			index[2] = int(1.*k + active_particle->rel_r[2]);
			
	
			//! 0: xMin, 1:xMax
			//! 2: yMin, 3:yMax
			//! 4: zMin, 5:zMax
			if(    (is_box_boundary[_Im1_] && index[0]<=0) 
			    || (is_box_boundary[_Ip1_] && index[0]>=BlkNds_X-1)
			    || (is_box_boundary[_Jm1_] && index[1]<=0) 
			    || (is_box_boundary[_Jp1_] && index[1]>=BlkNds_Y-1)
			    || (is_box_boundary[_Km1_] && index[2]<=0) 
			    || (is_box_boundary[_Kp1_] && index[2]>=BlkNds_Z-1))
			{

				delete_particle(species, i_j_k, part_index);
				num_out_of_Box[species]++;
				out_of_Box = true;
			}

		}





		//! --------- IN OBSTACLE ??? ------------------------------
		//! check whether block overlaps with obstacle
		//! IF IN OBS EREASE PARTICLE HERE TO AVOID REDUNDANT CALCULATIONS
		//! IN PARTICLE CELL REASSIGNEMENT
		if(is_box_boundary[_OBS_] && !out_of_Box)
		{
		
			INT32 cell_indices[3] = {i,j,k};
			intern2normedCoords(x, active_particle->rel_r, cell_indices);

			//! calculating squared distance is slightly faster
			r = vec_len2(x);
			
			D_REAL dist_to_orig_obs2 = -1;
			if(use_second_obstacle)
			dist_to_orig_obs2 = sqrt( (x[0]-Position_SecondObstacle[0])*(x[0]-Position_SecondObstacle[0])
						+(x[1]-Position_SecondObstacle[1])*(x[1]-Position_SecondObstacle[1]) 
						+ (x[2]-Position_SecondObstacle[2])*(x[2]-Position_SecondObstacle[2])
						);
			
		
			if(r < squared_R_Obstacle || (dist_to_orig_obs2 > 0 && dist_to_orig_obs2 < R_SecondObstacle) )
			{

//!---------------------- store particle impacting on obstacle ------------------------
#ifdef startTL_PROTOCOL_OBSTACLE_PARTICLE
				if(TL>=startTL_PROTOCOL_OBSTACLE_PARTICLE)
				{
					//! write TL
					ObstacleParticle_FILE << TL << "	";

					//! write weight
					ObstacleParticle_FILE << active_particle->weight << "	";

					//! write position in normed units
					ObstacleParticle_FILE << x[0] << " ";
					ObstacleParticle_FILE << x[1] << " ";
					ObstacleParticle_FILE << x[2] << "	";

					//! write velocity in normed units
					ObstacleParticle_FILE << active_particle->v[0] << " ";
					ObstacleParticle_FILE << active_particle->v[1] << " ";
					ObstacleParticle_FILE << active_particle->v[2] << " ";
	
					ObstacleParticle_FILE << endl;
				}
#endif
//!--------------------------------------------------------------------------------------

				//! DELETE ACTIVE PARTICLE
				delete_particle(species, i_j_k, part_index);
				num_in_obst[species]++;
				in_Obs = true;



			}//! end in obstacle
		}//! end is_box_boundary


#ifdef DELETE_PARTICLE_IF_V_EXCEEDS
		//! Delete particle at very high velocities which
		//! do not fulfill the courant criteria and thereby
		//! will make the move method crash
		//! -> this should rather be used for debugging
		//!    since ANY particle should fulfill the CC
		if(!out_of_Box && !in_Obs &&  vec_len2(active_particle->v) > DELETE_PARTICLE_IF_V_EXCEEDS*DELETE_PARTICLE_IF_V_EXCEEDS)
		{
			delete_particle(species, i_j_k, part_index);
			num_deleted_too_fast_particle++;
		}
#endif


	    }//! for particle of cell
	}//! for ijk


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


}

//!---------------------------------------------------------------
//! move_particle_noObsBoxBoundCheck:
//! same as move assuming block has no BB and does not overlap OBS
//! -> for more detailed comments see move
//!---------------------------------------------------------------
void CBlock::move_particle_noObsBoxBoundCheck(INT32 id_oct)
{



	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	INT32 new_i_j_k;
	INT32 index[3]; // allow for negativ numbers (which are deleted at once)
	
	PARTICLE_REAL r, x[3], temp_fastest_particle_v2;
	
	bool out_of_Box = false;
	
	
	
	//! temprary pointer to particle in list
	particle *active_particle;
	
	//! decide in init whether dt is different in indivudal levels
	PARTICLE_REAL dt_part = dt_particle_of_L[RLevel];
	
	
	//! con_base is needed as cell size for a particle (rel_r)
	//! ranges from 0 to 1 (= numerical units), wheras the physical
	//! cell size ranges from 0 to delta_cell (= physical units).
	//! As v is defined in physical units, con_base scales the dx 
	//! in moving step to numerical units.
	PARTICLE_REAL con_base[3] = {1./delta_of_L[RLevel][0],
				     1./delta_of_L[RLevel][1],
				     1./delta_of_L[RLevel][2]};




	//!	Some Particle remain unmoved for reassignement procedure
	//!	they could be moved by setting start to 0, but they are
	//!	not gathered anyway. Since they are only few (~1/1000000),
	//!	just leave it as before
	//!	(look xfig sketches for details)
	
	

	//! get index of oct (0 or 1 for each dimension)
	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);
	

	//! see move method for detailed comments
	INT32 start[3] = {!a *!is_box_boundary[0] +a*BlkNds_X/2,
			 !b *!is_box_boundary[2] +b*BlkNds_Y/2,
			 !c *!is_box_boundary[4] +c*BlkNds_Z/2};



	INT32 end[3]={BlkNds_X -a*!is_box_boundary[1] -!a*(BlkNds_X/2),
		     BlkNds_Y -b*!is_box_boundary[3] -!b*(BlkNds_Y/2),
		     BlkNds_Z -c*!is_box_boundary[5] -!c*(BlkNds_Z/2)};

	



  //!-----------------------------------------------------------
  //!1) Move Step:
  //!   All particle weights are assumed to be > 0 now !!!
  //!-----------------------------------------------------------
  for(short species=0; species<num_Charged_Species; species++)
   for(INT32 i = start[0]; i < end[0]; i++)
    for(INT32 j = start[1]; j < end[1]; j++)
     for(INT32 k = start[2]; k < end[2]; k++)
     {

        INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
	        +j*BlkNds_Z 
	        +k;


	if(max_PiC<num_MPiC[species][i_j_k])
	max_PiC = num_MPiC[species][i_j_k];
	

      	for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
      	{


		//! set particle pointer to active_particle
	 	active_particle = pArray[species][i_j_k] +part_index;


		//! ---- ADVANCE CARTESIEN POSITIONS ---------------------
		active_particle->rel_r[0] +=  dt_part * con_base[0] *active_particle->v[0];
		active_particle->rel_r[1] +=  dt_part * con_base[1] *active_particle->v[1];
		active_particle->rel_r[2] +=  dt_part * con_base[2] *active_particle->v[2];


		//! count particle as moved
		num_moved++;

		//!THIS HAS TO BE INTERCEPTED
		//! see move method for detailed comments
		for(int comp=0; comp<3; comp++)
		{
			if(    active_particle->rel_r[comp] > 0. -PART_REAL_PRECISION
			    && active_particle->rel_r[comp] < 0. +PART_REAL_PRECISION)
			    active_particle->rel_r[comp] = 0.;
			else
			if(    active_particle->rel_r[comp] > 1. -PART_REAL_PRECISION
			    && active_particle->rel_r[comp] < 1. +PART_REAL_PRECISION)
			    active_particle->rel_r[comp] = 1.;

		}

#ifdef ESTIMATE_FASTEST_PARTICLE
		//! JUST FOR STATISTICS:
		//! FOR SOME REASON CAN SLOW SIGNIFICANTLY DOWN THE CODE !!!
		//! use squared velocity for comparison.
		//! only build sqrt for statistics (slightly faster)

		temp_fastest_particle_v2 = vec_len2(active_particle->v);

                if(temp_fastest_particle_v2 > fastest_particle_v2[RLevel])
                    fastest_particle_v2[RLevel] = temp_fastest_particle_v2;
#endif

#ifdef DELETE_PARTICLE_IF_V_EXCEEDS
		//! Delete particle at very high velocities which
		//! do not fulfill the courant criteria and thereby
		//! will make the move method crash
		//! -> this should rather be used for debugging
		//!    since ANY particle should fulfill the CC
		if(vec_len2(active_particle->v) > DELETE_PARTICLE_IF_V_EXCEEDS*DELETE_PARTICLE_IF_V_EXCEEDS)
		{
			delete_particle(species, i_j_k, part_index);
			num_deleted_too_fast_particle++;
		}
#endif



	    }//! for particle of cell
	}//! for ijk


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


}

//!--------------------------------------------------------------
//! move_particle:
//! NOTES:
//!--------------------------------------------------------------
bool CBlock::cell_reassign_particle(void)
{


	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

	INT32 new_i_j_k;
	INT32 index[3]; //! allow for negativ numbers (which are deleted at once)

	//! temprary pointer to particle in list
	particle *active_particle;
	
	PARTICLE_REAL  x[3], temp_fastest_particle_v2;
	
	bool out_of_Box = false;
	

	//! In case Block is Box boundary, moving shall start at 0, else at 1
	INT32 start[3] = {!is_box_boundary[0],
			!is_box_boundary[2],
			!is_box_boundary[4]};
	
	
	INT32 end[3]   = {BlkNds_X -!is_box_boundary[1],
			 BlkNds_Y -!is_box_boundary[3],
			 BlkNds_Z -!is_box_boundary[5]};


	//!-----------------------------------------------------------
	//!2) reassign Step:
	//!-----------------------------------------------------------
	for(INT32 species=0; species<num_Charged_Species; species++)
	 for(INT32 i = start[0]; i < end[0]; i++)
	  for(INT32 j = start[1]; j < end[1]; j++)
	   for(INT32 k = start[2]; k < end[2]; k++)
	   {
	
		INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
			+j*BlkNds_Z 
			+k;


		for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
		{
	
			active_particle = pArray[species][i_j_k] +part_index;
	
	
			index[0] = int(1.*i + active_particle->rel_r[0]);
			index[1] = int(1.*j + active_particle->rel_r[1]);
			index[2] = int(1.*k + active_particle->rel_r[2]);
	
			//! --------- IN NEW CELL ??? ------------------------------
			if(index[0] != i || index[1] != j || index[2] != k)
			{
	
	
				new_i_j_k  = index[0]*BlkNds_Y*BlkNds_Z
					+index[1]*BlkNds_Z 
					+index[2];

{//! comments to interception 
// #ifdef DEBUG_MODUS
			
			if(   index[0] > BlkNds_X-1
			    || index[1] > BlkNds_Y-1
			    || index[2] > BlkNds_Z-1
			    || new_i_j_k<0)
			{


				log_file << "MALFUNCTION in cell_reassign_particle:" << endl;
				log_file << "New particle position out of Block:" << endl;
				log_file << "Block (" <<Blk_Index[0]<<","<<Blk_Index[1]<<","
						  <<Blk_Index[2] <<")" << endl;

				log_file <<  "Old cell (" <<i<<","<<j<<","<<k<<") " << endl;

				log_file <<  "New cell (" <<index[0]<<","<<index[1]<<","
				     		          <<index[2]<<") " << endl;


				INT32 Blk_Ind[3];
				if(use_SFC)
				{
					SFC_BlkNr_to_Indices(Block_Nr, Blk_Ind, RLevel);
	
				log_file << "In SimuBox (" <<Blk_Ind[0]<<
							","<<Blk_Ind[1]<<","
						  	   <<Blk_Ind[2] <<")" << endl;
				}



				log_file << "Level: " << RLevel << endl;

				log_file << "part_index: " << part_index << endl;

				log_file << "num_MPiC[species][i_j_k]: " << num_MPiC[species][i_j_k] << endl;

				log_file <<  "Particle velocity: " 
				     << "(" <<active_particle->v[0]<<","
				            <<active_particle->v[1]<<","
				            <<active_particle->v[2]<<") " << endl;

				log_file <<  "Particle position: " 
					 << "(" <<active_particle->rel_r[0]<<","
				         <<	  active_particle->rel_r[1]<<","
				         <<	  active_particle->rel_r[2]<<") " << endl;




				log_file <<  "exiting..." << endl;
				
				log_file <<  "THIS MAY BE FIXABLE BY: " << endl;
				log_file <<  " - using higher smoothing values." << endl;
				log_file <<  " - reducing the time step." << endl;
				log_file <<  " - increasing MCD." << endl;
// 				exit(1);
				return false;

			}

// #endif
}
//! end comments
				//! rel_r is either <0 or >1 -> in interval[-1:2].
				//! In order to get a positive number of [0:1] add 1
				//! and subtract the integer value
				for(int comp=0; comp<3; comp++)
				{
					active_particle->rel_r[comp] += 1.;
					active_particle->rel_r[comp] -= int(active_particle->rel_r[comp]);
				}

{//! comments to interception 
#ifdef DEBUG_MODUS
			if(   active_particle->rel_r[0] < 0. || active_particle->rel_r[0] >= 1.
			    || active_particle->rel_r[1] < 0. || active_particle->rel_r[1] >= 1.
			    || active_particle->rel_r[2] < 0. || active_particle->rel_r[2] >= 1.)
			{
				log_file << "ERROR in cell_reassign_particle - UNHOOK PROCEDURE:" << endl;


				log_file << "Species: " << species << endl;
	
				log_file <<  "----BEFORE MOVE:------" << endl;
				log_file <<  "(i,j,k)  (" <<i<<","<<j<<","<<k<<") " << endl;
	
	
				log_file <<  "----AFTER MOVE:-------" << endl;
	
				log_file <<  "index[] (" <<index[0]<<","
						     <<index[1]<<","
						     <<index[2]<<") " << endl;
	
				log_file << "(active_particle->rel_r[0] = "
				<<   active_particle->rel_r[0] << endl;
				log_file << "(active_particle->rel_r[1] = "
				<<   active_particle->rel_r[1] << endl;
				log_file << "(active_particle->rel_r[2] = "
				<<   active_particle->rel_r[2] << endl;
	
				log_file <<  "Particle velocity: " 
				<< "(" << active_particle->v[0]<<","
					<< active_particle->v[1]<<","
					<< active_particle->v[2]<<") " << endl;
	
	
				log_file << "Value is not allowed value for rel_r." << endl << endl;
				log_file << "USUALLY THIS MEANS PARTICLE TRAVELD FURTHER THAN ONE CELL PER TL." << endl;
				log_file << "SEE INFORMATION ABOVE FOR DETAILS (e.g. to high velocity?)" << endl;
				log_file << "exiting ..." << endl;
// 				exit(1);
				return false;
			}

#endif
}
//! end comments
				//! sort particle in new cells list
				swap_particles_pArray(species,
						      part_index,
						      i_j_k,
						      new_i_j_k,
						      active_particle);
			
				
				//! part_index=-1 may temporary occcur here
				part_index--;
				num_unhooked++;
		
			}//! end if new cell

{//! comments to interception 
#ifdef DEBUG_MODUS
		//! Even though this Error should not occur, it does 
		//! in case interception above is not used 
		//! -> below is just for double check.
		//! (see above for further details)
		if(     active_particle->rel_r[0] < 0. || active_particle->rel_r[0] >= 1.
		     || active_particle->rel_r[1] < 0. || active_particle->rel_r[1] >= 1.
		     || active_particle->rel_r[2] < 0. || active_particle->rel_r[2] >= 1.)
		{

			log_file << "ERROR in Move WITHOUT unhook:" << endl;


			log_file << "Species: " << species << endl;

			log_file <<  "----BEFORE MOVE:------" << endl;
			log_file <<  "(i,j,k)  (" <<i<<","<<j<<","<<k<<") " << endl;


			log_file <<  "----AFTER MOVE:-------" << endl;

			log_file <<  "index[] (" <<index[0]<<","<<index[1]<<","
				     	     <<index[2]<<") " << endl;

			log_file << "(active_particle->rel_r[0] = "
			<<   active_particle->rel_r[0] << endl;
			log_file << "(active_particle->rel_r[1] = "
			<<   active_particle->rel_r[1] << endl;
			log_file << "(active_particle->rel_r[2] = "
			<<   active_particle->rel_r[2] << endl;

			log_file <<  "Particle velocity: " 
			     << "(" << active_particle->v[0]<<","
				    << active_particle->v[1]<<","
				    << active_particle->v[2]<<") " << endl;



			log_file << "value is not allowed for rel_r." << endl;
			log_file << "exiting ..." << endl;
// 			exit(1);
			return false;
		}
#endif
}
//! end comments



	

      }//! end for particle list
    }//! end for Block 


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

	return true;

}

//!--------------------------------------------------------------
//! apply_inhom_particle_bounds:
//!--------------------------------------------------------------
void CBlock::apply_inhom_particle_bounds(const INT32 comp, INT32 species, const INT32 *start, const INT32 *end, const INT32 *src)
{

	INT32 ind[3];
	
	PARTICLE_REAL x[3];
	PARTICLE_REAL v_init[3] = {0., 0., 0.};

	PARTICLE_REAL rho_init;
	rho_init = rho_sw[species];

	for(ind[0] = start[0]; ind[0] < end[0]; ind[0]++)
	 for(ind[1] = start[1]; ind[1] < end[1]; ind[1]++)
	  for(ind[2] = start[2]; ind[2] < end[2]; ind[2]++)
	  {
		
		INT32 dest_i_j_k  =   ind[0]*BlkNds_Y*BlkNds_Z 
				    +ind[1]*BlkNds_Z 
				    +ind[2];

		INT32 src_i_j_k  =   (ind[0]+src[0])*BlkNds_Y*BlkNds_Z 
				   +(ind[1]+src[1])*BlkNds_Z 
				   +(ind[2]+src[2]);


		//! estimate coordinate of cell's centre
		cell_centre_normedCoords(x, ind);

		//! choose velocity depending on coordinate
		set_inflow_velocity(v_init, x, species);

		//! choose density depending on coordinate
		set_inflow_density(rho_init, x, species);
	
		//! copy cell when velocity is directed towards BB
		if(v_init[comp] *src[comp]<0.)
		copy_empty_Cell(species, dest_i_j_k, src_i_j_k);
		else
		fill_empty_Cell(species, dest_i_j_k, v_init, rho_init);

	
	   }

}

//!--------------------------------------------------------------
//! move_particle:
//! NOTES:
//!--------------------------------------------------------------
void CBlock::inject_boundary_particle(INT32 id_oct)
{



	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();

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





	//! Now all particle for all species are moved,
	//! Inject particles at boundaries.
	//! Do this for all inflow species
	for(INT32 inflow_species=0; inflow_species<num_Inflow_Species; inflow_species++)
	{

		INT32 species = index_Inflow_Species[inflow_species];


		//! check whether to inject respective species
		bool do_inject = false;

		//! always inject in case start is set to zero
		if(start_inject_species_each_xx_t0[species] == 0.)
		do_inject = true;
		//! else check wheter time is in interval start +duration
		else
		{

			D_REAL time_in_t0 = TL* dt;
			INT32 fac_start_t0 = time_in_t0 / start_inject_species_each_xx_t0[species];
		
			D_REAL start_inject = fac_start_t0 * start_inject_species_each_xx_t0[species];
			D_REAL  end_inject = start_inject +duration_inject_species_each_xx_t0[species];

			if(time_in_t0>=start_inject && time_in_t0<end_inject)
			do_inject = true;

		}


		if(do_inject)
		{

			// left X Boundary
			// "corrected" filling procedure:
			// start at j=0 & k=0 to fill entire i=0 layer
			//! simply leave out diagonal cell
			//! does not introduce any problems
			//! and is faster. also starting to insert at j,k=0
			//! somehow results in to high rho-boundary-corner values
	
			//!------------------------------------------------
			//!----------- I MIN Boundary ---------------------
			//!------------------------------------------------
	
			if(!a && is_box_boundary[Xmin_BB])
			{
	
				//! requiered for inhom boundaries
				const INT32 comp = 0;
	
				//! specify from which volume to copy
				const INT32 offset[3] = {1, 0, 0};
	
				//! specify which volume to fill
				const INT32 start[3]  = {0, oct_start[1], oct_start[2]};
				const INT32   end[3]  = {1,   oct_end[1],   oct_end[2]};
	
				if(use_hom_particle_bounds[Xmin_BB])
				{
					//! A) FILL CELL
					if(use_particle_inflow_bounds[Xmin_BB])
					fill_empty_Volume(species, start, end);
					else
					//! B) COPY PENULTIMATE CELL PARTICLE TO ULTIMATE
					if(use_particle_copy_bounds[Xmin_BB])
					copy_empty_Volume(species, start, end, offset);
					//! C) DO NOTHING, LEAVE CELL EMPTY
		
				}
				else
				apply_inhom_particle_bounds(comp, species, start, end, offset);
			}
			
	
	
	
			//!------------------------------------------------
			//!----------- I MAX Boundary ---------------------
			//!------------------------------------------------
			if( a && is_box_boundary[Xmax_BB])
			{
	
				//! requiered for inhom boundaries
				const INT32 comp = 0;
	
				//! specify from which volume to copy
				const INT32 offset[3] = {-1, 0, 0};
	
				//! specify which volume to fill
				const INT32 start[3] = {BlkNds_X-1,oct_start[1],oct_start[2]};
				const INT32 end[3]   = {BlkNds_X-0,  oct_end[1],  oct_end[2]};
	
	
				if(use_hom_particle_bounds[Xmax_BB])
				{
	
					//! A) FILL CELL
					if(use_particle_inflow_bounds[Xmax_BB])
					fill_empty_Volume(species, start, end);
					else
					//! B) COPY PENULTIMATE CELL PARTICLE TO ULTIMATE
					if(use_particle_copy_bounds[Xmax_BB])
					copy_empty_Volume(species, start, end, offset);
					//! C) DO NOTHING, LEAVE CELL EMPTY
	
				}
				else
				apply_inhom_particle_bounds(comp, species, start, end, offset);
	
			}
	
	
	
			//!------------------------------------------------
			//!----------- J MIN Boundary ---------------------
			//!------------------------------------------------
			if(!b && is_box_boundary[Ymin_BB])
			{
	
				//! requiered for inhom boundaries
				const INT32 comp = 1;
	
				//! specify from which volume to copy
				const INT32 offset[3] = {0, 1, 0};
	
				//! specify which volume to fill
				const INT32 start[3] = {oct_start[0], 0, oct_start[2]};
				const INT32 end[3]   = {  oct_end[0], 1,   oct_end[2]};
	
	
				if(use_hom_particle_bounds[Ymin_BB])
				{
	
					//! A) FILL CELL
					if(use_particle_inflow_bounds[Ymin_BB])
					fill_empty_Volume(species, start, end);
					else
					//! B) COPY PENULTIMATE CELL PARTICLE TO ULTIMATE
					if(use_particle_copy_bounds[Ymin_BB])
					copy_empty_Volume(species, start, end, offset);
					//! C) DO NOTHING, LEAVE CELL EMPTY
				}
				else
				apply_inhom_particle_bounds(comp, species, start, end, offset);
		
			}
	
	
	
			//!------------------------------------------------
			//!----------- J MAX Boundary --------------------
			//!------------------------------------------------
			if(b && is_box_boundary[Ymax_BB])
			{
	
				//! requiered for inhom boundaries
				const INT32 comp = 1;
	
				//! specify from which volume to copy
				const INT32 offset[3] = {0, -1, 0};
	
				//! specify which volume to fill
				const INT32 start[3] = {oct_start[0], BlkNds_Y-1, oct_start[2]};
				const INT32 end[3]   = {  oct_end[0], BlkNds_Y-0,   oct_end[2]};
	
	
				if(use_hom_particle_bounds[Ymax_BB])
				{
	
					//! A) FILL CELL
					if(use_particle_inflow_bounds[Ymax_BB])
					fill_empty_Volume(species, start, end);
					else
					//! B) COPY PENULTIMATE CELL PARTICLE TO ULTIMATE
					if(use_particle_copy_bounds[Ymax_BB])
					copy_empty_Volume(species, start, end, offset);
					//! C) DO NOTHING, LEAVE CELL EMPTY
				}
				else
				apply_inhom_particle_bounds(comp, species, start, end, offset);
	
	
			}
	
	
			//!------------------------------------------------
			//!----------- K MIN Boundary ---------------------
			//!------------------------------------------------
			if(!c && is_box_boundary[Zmin_BB])
			{
	
	
				//! requiered for inhom boundaries
				const INT32 comp = 2;
	
				//! specify from which volume to copy
				const INT32 offset[3] = {0, 0, 1};
	
				//! specify which volume to fill
				const INT32 start[3] = {oct_start[0], oct_start[1], 0};
				const INT32 end[3]   = {  oct_end[0],   oct_end[1], 1};
	
	
				if(use_hom_particle_bounds[Zmin_BB])
				{
	
					//! A) FILL CELL
					if(use_particle_inflow_bounds[Zmin_BB])
					fill_empty_Volume(species, start, end);
					else
					//! B) COPY PENULTIMATE CELL PARTICLE TO ULTIMATE
					if(use_particle_copy_bounds[Zmin_BB])
					copy_empty_Volume(species, start, end, offset);
					//! C) DO NOTHING, LEAVE CELL EMPTY
	
				}
				else
				apply_inhom_particle_bounds(comp, species, start, end, offset);
		
			}
	
	
	
			//!------------------------------------------------
			//!----------- K MAX Boundary --------------------
			//!------------------------------------------------
			if(c && is_box_boundary[Zmax_BB])
			{
	
	
				//! requiered for inhom boundaries
				const INT32 comp = 2;
	
				//! specify from which volume to copy
				const INT32 offset[3] = {0, 0, -1};
	
				//! specify which volume to fill
				const INT32 start[3] = {oct_start[0], oct_start[1], BlkNds_Z-1};
				const INT32 end[3]   = {  oct_end[0],   oct_end[1], BlkNds_Z-0};
	
	
				if(use_hom_particle_bounds[Zmax_BB])
				{
	
	
					//! A) FILL CELL
					if(use_particle_inflow_bounds[Zmax_BB])
					fill_empty_Volume(species, start, end);
					else
					//! B) COPY PENULTIMATE CELL PARTICLE TO ULTIMATE
					if(use_particle_copy_bounds[Zmax_BB])
					copy_empty_Volume(species, start, end, offset);
					//! C) DO NOTHING, LEAVE CELL EMPTY
	
				}
				else
				apply_inhom_particle_bounds(comp, species, start, end, offset);
	
			}
	
			//!------------------------------------------------
			//!----------- FILL EDGES AT MIN BB ---------------
			//!----- If this ist not done density lacks -------
			//!----- at 75% backround density will result -----
			//!----- at min boundries. ------------------------
			//!------------------------------------------------
	
			/*
			//!------------------------------------------------
			//!----------- X-Edge  ----------------------------
			//!------------------------------------------------
			if(!b && is_box_boundary[Ymin_BB] && !c && is_box_boundary[Zmin_BB])
			{
	
				//! specify which volume to fill
				const INT32 start[3]  = {oct_start[0], 0, 0};
				const INT32   end[3]  = {  oct_end[0], 1, 1};
				
	
				fill_empty_Volume(species, start, end);
	
	
			}
	
	
			//!------------------------------------------------
			//!----------- Y-Edge  ----------------------------
			//!------------------------------------------------
			if(!a && is_box_boundary[Xmin_BB] && !c && is_box_boundary[Zmin_BB])
			{
	
				//! specify which volume to fill
				const INT32 start[3]  = {0, oct_start[1], 0};
				const INT32   end[3]  = {1,   oct_end[1], 1};
	
	
				fill_empty_Volume(species, start, end);
	
	
			}
	
			//!------------------------------------------------
			//!----------- Z-Edge  ----------------------------
			//!------------------------------------------------
			if(!a && is_box_boundary[Xmin_BB] && !b && is_box_boundary[Ymin_BB])
			{
	
				//! specify which volume to fill
				const INT32 start[3]  = {0, 0, oct_start[2]};
				const INT32   end[3]  = {1, 1,   oct_end[2]};
	
	
				fill_empty_Volume(species, start, end);
	
	
			}
			*/
		} //! end do inject

 	}//! for inflow species


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


}



//!--------------------------------------------------------------
//! reset_particle_time_of_L:
//!--------------------------------------------------------------
void CBlock::reset_particle_time_of_L(void)
{


  for(short species=0; species<num_Charged_Species; species++)
   for(INT32 i = 0; i < BlkNds_X; i++)
    for(INT32 j = 0; j < BlkNds_Y; j++)
     for(INT32 k = 0; k < BlkNds_Z; k++)
     {

       INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
	       +j*BlkNds_Z 
	       +k;


//       	for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
// 	pArray[species][i_j_k][part_index].time = 0;


    }


}

//!-----------------------------------------------------------
//! calc_Energy_Momentum_Mass
//!-----------------------------------------------------------
void CBlock::calc_Energy_Momentum_Mass(void)
{



	particle *active_particle;
	

	D_REAL Cell_Volume = CellVol_of_L[RLevel]/(LX*LY*LZ);


	//! - only add up physical nodes that are not covered by refined octs
	//! - multiply with celvolume (integral)
	for(INT32 id_oct=0; id_oct<8; id_oct++)
	 if(child_array[id_oct] == NULL)
	 {
		//! get index of oct (0 or 1 for each dimension)
		INT32 a = id_oct/4;
		INT32 b = (id_oct -a*4)/2;
		INT32 c = (id_oct -a*4 -b*2);
		
	
		//! eg BlkNds=6
		//! -> 1-2
		//! -> 3-4

		//! see move method for detailed comments
		INT32 start[3] = {1+a*(BlkNds_X/2-1),
				  1+b*(BlkNds_Y/2-1),
				  1+c*(BlkNds_Z/2-1)};
	
	
	
		INT32 end[3]={BlkNds_X-1 -!a*(BlkNds_X/2-1),
			      BlkNds_Y-1 -!b*(BlkNds_Y/2-1),
			      BlkNds_Z-1 -!c*(BlkNds_Z/2-1)};

	
	
		for(INT32 i = start[0]; i < end[0]; i++)
		 for(INT32 j = start[1]; j < end[1]; j++)
		  for(INT32 k = start[2]; k < end[2]; k++)
		  {
		
			INT32 i_j_k  =   i*BlkNds_Y*BlkNds_Z 
					+j*BlkNds_Z 
					+k;
		
	
			//! --- BField ------------------------------
			D_REAL *B1 = Field_Type[id_BEven] +0*num_nodes_in_block;
			D_REAL *B2 = Field_Type[id_BEven] +1*num_nodes_in_block;
			D_REAL *B3 = Field_Type[id_BEven] +2*num_nodes_in_block;
	
			total_magnetic_field[0] += B1[i_j_k] *Cell_Volume;
			total_magnetic_field[1] += B2[i_j_k] *Cell_Volume;
			total_magnetic_field[2] += B3[i_j_k] *Cell_Volume;
	
			//! --- UField ------------------------------
			D_REAL *U1 = Field_Type[id_UI_plus] +0*num_nodes_in_block;
			D_REAL *U2 = Field_Type[id_UI_plus] +1*num_nodes_in_block;
			D_REAL *U3 = Field_Type[id_UI_plus] +2*num_nodes_in_block;
		
			total_collected_momentum[0] += U1[i_j_k] *Cell_Volume;
			total_collected_momentum[1] += U2[i_j_k] *Cell_Volume;
			total_collected_momentum[2] += U3[i_j_k] *Cell_Volume;
	
			//! --- TField ------------------------------
			D_REAL *T1 = Field_Type[id_PISpecies1] +0*num_nodes_in_block;
			D_REAL *T2 = Field_Type[id_PISpecies1] +1*num_nodes_in_block;
			D_REAL *T3 = Field_Type[id_PISpecies1] +2*num_nodes_in_block;
		
			total_collected_temperature[0] += T1[i_j_k] *Cell_Volume;
			total_collected_temperature[1] += T2[i_j_k] *Cell_Volume;
			total_collected_temperature[2] += T3[i_j_k] *Cell_Volume;
	
	
			//! --- Rho Field ---------------------------
			D_REAL *Rho = Field_Type[id_rho_np1];
			total_collected_density += Rho[i_j_k] *Cell_Volume;
	
	
		
			for(short species=0; species<num_Charged_Species; species++)
			for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
			{
			
					active_particle = pArray[species][i_j_k] +part_index;
				
					total_particle_mass   +=      Ion_Masses[species] *active_particle->weight;
					total_particle_energy += 0.5 *Ion_Masses[species] *active_particle->weight *vec_len2(active_particle->v);
				
					total_particle_momentum[0] += Ion_Masses[species] *active_particle->weight *active_particle->v[0];
					total_particle_momentum[1] += Ion_Masses[species] *active_particle->weight *active_particle->v[1];
					total_particle_momentum[2] += Ion_Masses[species] *active_particle->weight *active_particle->v[2];
			
					num_total_particles_collected++;
			
		  }//! end pArray
		}//! end for oct


	}//! end all octs

}

