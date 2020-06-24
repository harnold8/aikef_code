


#include "CBlk.h"
#include "utils.h"
#include "defines.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>





//!-------------------------------------------------------------//
//! refine_Blocks: - call by reference function			//
//!-------------------------------------------------------------//
bool CBlock::is_refinement_permitted(INT32 id_oct)
{


	if(child_array[id_oct]) return false;


	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);


	//!---------------------------------------------
	//!---- Check for Faces --- (6 cases) ----------
	//!---------------------------------------------

	if(    (!a && !Neighbour[_Im1_] && !is_box_boundary[0])
	    || ( a && !Neighbour[_Ip1_] && !is_box_boundary[1])
	    || (!b && !Neighbour[_Jm1_] && !is_box_boundary[2])
	    || ( b && !Neighbour[_Jp1_] && !is_box_boundary[3])
	    || (!c && !Neighbour[_Km1_] && !is_box_boundary[4])
	    || ( c && !Neighbour[_Kp1_] && !is_box_boundary[5]))
	return false;


	//!---------------------------------------------
	//!---- Check for Edges --- (12 cases) ---------
	//!---------------------------------------------

	//! ------ i-Direction -------------------------
	if(Neighbour[_Im1_])
	 if(   (!a && !b && !Neighbour[_Im1_]->Neighbour[_Jm1_] && !Neighbour[_Im1_]->is_box_boundary[2])
	     || (!a &&  b && !Neighbour[_Im1_]->Neighbour[_Jp1_] && !Neighbour[_Im1_]->is_box_boundary[3])
	     || (!a && !c && !Neighbour[_Im1_]->Neighbour[_Km1_] && !Neighbour[_Im1_]->is_box_boundary[4])
	     || (!a &&  c && !Neighbour[_Im1_]->Neighbour[_Kp1_] && !Neighbour[_Im1_]->is_box_boundary[5]))
	  return false;

	if(Neighbour[_Ip1_])
	 if(   (a && !b && !Neighbour[_Ip1_]->Neighbour[_Jm1_] && !Neighbour[_Ip1_]->is_box_boundary[2])
	     || (a &&  b && !Neighbour[_Ip1_]->Neighbour[_Jp1_] && !Neighbour[_Ip1_]->is_box_boundary[3])
	     || (a && !c && !Neighbour[_Ip1_]->Neighbour[_Km1_] && !Neighbour[_Ip1_]->is_box_boundary[4])
	     || (a &&  c && !Neighbour[_Ip1_]->Neighbour[_Kp1_] && !Neighbour[_Ip1_]->is_box_boundary[5]))
	  return false;

	//! ------ j-Direction -------------------------
	//! Remember:
	//! "Neighbour[_Jm1_]->Neighbour[_Im1_]"
	//! 	equals 
	//! "Neighbour[_Im1_]->Neighbour[_Jm1_]"

	if(Neighbour[_Jm1_])
	 if(    (!b && !c && !Neighbour[_Jm1_]->Neighbour[_Km1_] && !Neighbour[_Jm1_]->is_box_boundary[4])
	     || (!b &&  c && !Neighbour[_Jm1_]->Neighbour[_Kp1_] && !Neighbour[_Jm1_]->is_box_boundary[5]))
	  return false;

	if(Neighbour[_Jp1_])
	 if(    (b && !c && !Neighbour[_Jp1_]->Neighbour[_Km1_] && !Neighbour[_Jp1_]->is_box_boundary[4])
	     || (b &&  c && !Neighbour[_Jp1_]->Neighbour[_Kp1_] && !Neighbour[_Jp1_]->is_box_boundary[5]))
	  return false;

	//! ------ k-Direction -------------------------
	//! (all cases checkedf in i and j Direction)



	//!---------------------------------------------
	//!---- Check for Diagonals --- (8 cases) ------
	//!---------------------------------------------

	//! ------------i- min -Direction ---------------------------
	if(Neighbour[_Im1_] && Neighbour[_Jm1_])
	if(    (!a && !b && !c && !Neighbour[_Im1_]->Neighbour[_Jm1_]->Neighbour[_Km1_] && !Neighbour[_Im1_]->Neighbour[_Jm1_]->is_box_boundary[4])
	    || (!a && !b &&  c && !Neighbour[_Im1_]->Neighbour[_Jm1_]->Neighbour[_Kp1_] && !Neighbour[_Im1_]->Neighbour[_Jm1_]->is_box_boundary[5]))
	  return false;

	if(Neighbour[_Im1_] && Neighbour[_Jp1_])
	if(   (!a && b && !c && !Neighbour[_Im1_]->Neighbour[_Jp1_]->Neighbour[_Km1_] && !Neighbour[_Im1_]->Neighbour[_Jp1_]->is_box_boundary[4])
	    || (!a && b &&  c && !Neighbour[_Im1_]->Neighbour[_Jp1_]->Neighbour[_Kp1_] && !Neighbour[_Im1_]->Neighbour[_Jp1_]->is_box_boundary[5]))
	  return false;


	//! ------------i- plus -Direction ---------------------------
	if(Neighbour[_Ip1_] && Neighbour[_Jm1_])
	if(    (a && !b && !c && !Neighbour[_Ip1_]->Neighbour[_Jm1_]->Neighbour[_Km1_] && !Neighbour[_Ip1_]->Neighbour[_Jm1_]->is_box_boundary[4])
	    || (a && !b &&  c && !Neighbour[_Ip1_]->Neighbour[_Jm1_]->Neighbour[_Kp1_] && !Neighbour[_Ip1_]->Neighbour[_Jm1_]->is_box_boundary[5]))
	  return false;

	if(Neighbour[_Ip1_] && Neighbour[_Jp1_])
	if(   (a && b && !c && !Neighbour[_Ip1_]->Neighbour[_Jp1_]->Neighbour[_Km1_] && !Neighbour[_Ip1_]->Neighbour[_Jp1_]->is_box_boundary[4])
	    || ( a &&  b &&  c && !Neighbour[_Ip1_]->Neighbour[_Jp1_]->Neighbour[_Kp1_] && !Neighbour[_Ip1_]->Neighbour[_Jp1_]->is_box_boundary[5]))
	  return false;


	return true;


}

//!-------------------------------------------------------------//
//! is_removing_permitted							//
//!-------------------------------------------------------------//
bool CBlock::is_removing_permitted(INT32 id_oct)
{

	//! if any oct of child is refined return false
	for(INT32 child_of_child=0; child_of_child<8; child_of_child++)
	 if(child_array[id_oct]->child_array[child_of_child])
	 {
		num_rem_reject_child_is_refined++;
		return false;
	 }


	//! return in case block is flagged as initial_refinement
	if(child_array[id_oct]->initial_refined && never_remove_static_refinement)
	{
		num_rem_reject_initial_refined++;
		return false;
	}


	//! always allow to remove MAX_LEVEL Block
	if(RLevel==MAX_LEVEL-1 && !child_array[id_oct]->initial_refined)
	return true;

	//! method below not correct yet, for the moment
	//! do not allow for removing
	//! SHOULD WORK NOW
// 	return false;

	//! check whether refined neighbour exists
	if(child_array[id_oct]->has_refined_neighbour())
	{
		num_rem_reject_refined_neighbour++;
		return false;
	}

	//! no neighbour is refined -> allow for removing
	return true;
}


//!-------------------------------------------------------------//
//! flag_full_environment: - 							//
//! in order to avoid discontonuities on level boundaries,		//
//! also flag environment in THIS level for refinement		//
//!-------------------------------------------------------------//
void CBlock::flag_full_environment(INT32 oct, INT32 *num_flagged_neighbours)
{


	//! never refine gather Blocks
	if(is_gatherBlk)
	{
		log_file << " WARNING:" << endl
		         << " tried to flag ENVIRONMENT OF GATHER BLOCK."
		         << " Ignoring command and continuing ... " << endl;
		return;
	}


	//! never refine gather Blocks
	if(RLevel==MAX_LEVEL)
	{
		log_file << " WARNING:" << endl
		         << " tried to flag ENVIRONMENT IN MAX_LEVEL."
		         << " Ignoring command and continuing ... " << endl;
		return;
	}



	//! get index of oct (0 or 1 for each dimension)
	INT32 a =  oct/4;
	INT32 b = (oct -a*4)/2;
	INT32 c = (oct -a*4 -b*2);

	//! set position to respective CHILD of neighbour
	D_REAL centre_of_flagged_child[3] = {origin[0] +(0.5 +a) *0.5 *Blk_Length_of[RLevel][0],
					     origin[1] +(0.5 +b) *0.5 *Blk_Length_of[RLevel][1],
					     origin[2] +(0.5 +c) *0.5 *Blk_Length_of[RLevel][2]};


	//! will be set by seek_Block_at function
	INT32 child_id_at_position = 0;

	//! Try to refine the 26 (=3^3-1) Block around this Block
	//! in this level
	INT32 start[3] = {!is_box_boundary[0],
			 !is_box_boundary[2],
			 !is_box_boundary[4]};


	INT32   end[3] = {!is_box_boundary[1],
			 !is_box_boundary[3],
			 !is_box_boundary[5]};


	//! flag Neighbours around respective oct
	for(INT32 u=-start[0]; u<=end[0]; u++)
	 for(INT32 v=-start[1]; v<=end[1]; v++)
	  for(INT32 w=-start[2]; w<=end[2]; w++)
	  if(u||v||w)
	  {


		//! set position to respective CHILD of neighbour
		//! only use 0.5 times the Blk_Length of this level
		//! in order to get correct child_id
		D_REAL position[3] = {centre_of_flagged_child[0] +0.5 *u *Blk_Length_of[RLevel][0],
				      centre_of_flagged_child[1] +0.5 *v *Blk_Length_of[RLevel][1],
				      centre_of_flagged_child[2] +0.5 *w *Blk_Length_of[RLevel][2]};


		//! seek neighbour in this level or below
		CBlock *block_to_flag = seek_Block_at(RLevel, position, child_id_at_position);


		//! only flag if block in case it is not yet
		//! (if this is skipped, regular flagged blocks will be marked as neighbour
		//!   flagged which in turn will trigger less block refinement)
		if(!block_to_flag->child_flag_refinement[child_id_at_position])
		{


			//! remember that this block was not flagged by refinement criteria
			//! but by its neighbour
			//! If this is not done, this block will trigger aditional blocks
			//! for refinement
			block_to_flag->flagged_by_neighbour[child_id_at_position] = true;

			//! flag respective oct
			block_to_flag->child_flag_refinement[child_id_at_position] = true;
	
			//! if block_to_flag is in this Level,
			//! protocol it in first half of array
			if(block_to_flag->RLevel == RLevel)
			num_flagged_neighbours[RLevel]++;
	
			//! if block_to_flag is in smaller Level,
			//! protocol it in second half of array
			if(block_to_flag->RLevel < RLevel)
			num_flagged_neighbours[block_to_flag->RLevel +(MAX_LEVEL+1)]++;
		}


	  }

}

//!-------------------------------------------------------------//
//! has_refined_neighbour							//
//! e.g. to check whether removing is permitted			//
//!-------------------------------------------------------------//
CBlock* CBlock::has_refined_neighbour(void)
{


	//! if this block is within MAX_LEVEL, no neighbour can be refined,
	//! so always return false in this case
	if(RLevel==MAX_LEVEL)
	return false;


	//! (quite likely, this Block has no children when this function is called
	//! (yet it could have)

	//! check for all 8 potential children, whether a neighbour exist,
	//! 4^3 - 2^3 = 56 checks (= 7 checks each child =2^3-1)
	//! if any neighbour exist, return true

	//! parameter that is expected by seek_Block_at function
	//! but not used here
	INT32 dummy=0;

	for(INT32 a=0; a<2; a++)
	 for(INT32 b=0; b<2; b++)
	  for(INT32 c=0; c<2; c++)
	  {


		D_REAL child_origin[3] = {origin[0] +0.5 *a *Blk_Length_of[RLevel][0],
					  origin[1] +0.5 *b *Blk_Length_of[RLevel][1],
					  origin[2] +0.5 *c *Blk_Length_of[RLevel][2]};


		//! Start:
		//! a=0 && is  m1BB -> start at  0
		//! a=0 && not m1BB -> start at -1
		//! a=1 && any m1BB -> start at  0
		//! c,b respectively
		INT32 start[3] = {!a * !is_box_boundary[0],
				 !b * !is_box_boundary[2],
				 !c * !is_box_boundary[4]};
	
		//! End:
		//! a=1 && is  p1BB -> end at  0
		//! a=1 && not p1BB -> end at  1
		//! a=0 && any p1BB -> end at  0
		//! c,b respectively
		INT32   end[3] = {a * !is_box_boundary[1],
				 b * !is_box_boundary[3],
				 c * !is_box_boundary[5]};

		for(INT32 u=-start[0]; u<=end[0]; u++)
		 for(INT32 v=-start[1]; v<=end[1]; v++)
		  for(INT32 w=-start[2]; w<=end[2]; w++)
		  if(u||v||w)
		  {
	
			//! NOTE:
			//! shift position 50% of childs's length into child
			//! IF THIS IS NOT DONE ROUND OFF ERRORS WILL OCCUR !!!
			D_REAL position[3] = {child_origin[0] +(0.5 +u)*Blk_Length_of[RLevel+1][0],
					      child_origin[1] +(0.5 +v)*Blk_Length_of[RLevel+1][1],
					      child_origin[2] +(0.5 +w)*Blk_Length_of[RLevel+1][2]};
	
			
			CBlock *block_to_check = seek_Block_at(RLevel+1, position, dummy);

			if(block_to_check->RLevel == RLevel+1)
			return block_to_check;
	
		  }
	}


	//! else no neighbour is refined so return false
	return NULL;


}



//!-------------------------------------------------------------//
//! set_Blk_optimal_MPiC:					  		//
//!			  							//
//!-------------------------------------------------------------//
void CBlock::set_Blk_optimal_MPiC(void)
{

	//! set Blk_optimal_MPiC to optimal_MPiC by default
	for(INT32 species=0; species<num_Charged_Species; species++)
	Blk_optimal_MPiC[species] = optimal_MPiC[species];

	//! check whether this block is refinement boundary
// 	if(has_refined_neighbour())
// 	for(INT32 species=0; species<num_Charged_Species; species++)
// 	Blk_optimal_MPiC[species] =  int(fac_oMPiC_at_LevBoundBlks *optimal_MPiC[species]);


	//! below only effects level lower than MAX_LEVEL !!!
	if(RLevel==MAX_LEVEL)
	return;

	
	//! TODO:
	//! The maximum of all child should be taken here instead
	//! of the value of the first child that is found
// 	CBlock* child=NULL;
// 	if(child = has_refined_neighbour())
// 	for(INT32 species=0; species<num_Charged_Species; species++)
// 	Blk_optimal_MPiC[species] =  int(fac_oMPiC_at_LevBoundBlks*child->Blk_optimal_MPiC[species]);


        //! check whether any oct is refined:
//         for(INT32 oct=0; oct<8; oct++)
//         if(child_array[oct])
// 	for(INT32 species=0; species<num_Charged_Species; species++)
// 	if(Blk_optimal_MPiC[species] <  int(fac_oMPiC_at_LevBoundBlks*child_array[oct]->Blk_optimal_MPiC[species]))
//         Blk_optimal_MPiC[species] =  int(fac_oMPiC_at_LevBoundBlks*child_array[oct]->Blk_optimal_MPiC[species]);


	//! store pointer all children and sorrounding neighbours of RLevel+1
// 	CBlock* children_and_their_neighbours[64];
// 	memset(children_and_their_neighbours, 0, 64*sizeof(CBlock*));


	CBlock* temp_Block;

	//! "shifted_origin" is centered in child (0,0,0)
	D_REAL shifted_origin[3] = {origin[0] +0.25 *Blk_Length_of[RLevel][0],
				    origin[1] +0.25 *Blk_Length_of[RLevel][1],
				    origin[2] +0.25 *Blk_Length_of[RLevel][2]};

	//! loop across all potential children and their neighbours
	//! which is a 4x4x4 Block Volume
	INT32 dummy=0;
	for(INT32 a=is_box_boundary[_Im1_]; a<4 -is_box_boundary[_Ip1_]; a++)
	 for(INT32 b=is_box_boundary[_Jm1_]; b<4 -is_box_boundary[_Jp1_]; b++)
	  for(INT32 c=is_box_boundary[_Km1_]; c<4 -is_box_boundary[_Kp1_]; c++)
	  {


		D_REAL position[3] = {shifted_origin[0] +(a-1) *Blk_Length_of[RLevel+1][0],
				      shifted_origin[1] +(b-1) *Blk_Length_of[RLevel+1][1],
				      shifted_origin[2] +(c-1) *Blk_Length_of[RLevel+1][2]};


		temp_Block = seek_Block_at(RLevel+1, position, dummy);


		 //! Do not use method at the end of this function since to many
		 //! paricle are produced
		 if(temp_Block->RLevel == RLevel+1)
		 {
		 	for(INT32 species=0; species<num_Charged_Species; species++)
		 	Blk_optimal_MPiC[species] = int(fac_oMPiC_at_LevBoundBlks*optimal_MPiC[species]);
// 		 	if(Blk_optimal_MPiC[species] < int(fac_oMPiC_at_LevBoundBlks*temp_Block->Blk_optimal_MPiC[species]))
// 		 	    Blk_optimal_MPiC[species] = int(fac_oMPiC_at_LevBoundBlks*temp_Block->Blk_optimal_MPiC[species]);

		 }
		 
		if(eta_Alfven_Wing_boundary)
			if(position[2]>(LZ - Box_Origin[2] - resistive_bound_dist[4]) ||position[2]<( -Box_Origin[2] + resistive_bound_dist[5]) )
			{
				
					for(INT32 species=0; species<num_Charged_Species; species++)
					Blk_optimal_MPiC[species] = int(0.2*optimal_MPiC[species]) + 1;
			}

	}

	//! NOTE:
	//! Cildrens oMPiC is set at this point since method is called from MAX_LEVEL down to level 0
	//! find maximum of Blk_optimal_MPiC in children and their neighbours an set it to
	//! this blocks Blk_optimal_MPiC
// 	for(INT32 blk=0; blk<64; blk++)
// 	 if(children_and_their_neighbours[blk])
// 	  if(children_and_their_neighbours[blk]->RLevel == RLevel+1)
// 	   for(INT32 species=0; species<num_Charged_Species; species++)
// 	    if(Blk_optimal_MPiC[species] < int(fac_oMPiC_at_LevBoundBlks*children_and_their_neighbours[blk]->Blk_optimal_MPiC[species]))
// 	     Blk_optimal_MPiC[species] = int(fac_oMPiC_at_LevBoundBlks*children_and_their_neighbours[blk]->Blk_optimal_MPiC[species]);
		


		


}


//!-------------------------------------------------------------//
//! refine_oct_environment: - 						//
//! in order to allow for refinement
//!-------------------------------------------------------------//
void CBlock::refine_oct_environment(INT32 id_oct)
{



	//! never refine gather Blocks
	if(is_gatherBlk)
	{
		log_file << " WARNING:" << endl
		         << " tried to refine block in gather blocks."
		         << " Ignoring command and continuing ... " << endl;
		return;
	}

	//! get index of oct (0 or 1 for each dimension)
	INT32 a =  id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);

	D_REAL oct_centre[3] = { origin[0] +(0.25 +a*0.5)*Blk_Length_of[RLevel][0],
				 origin[1] +(0.25 +b*0.5)*Blk_Length_of[RLevel][1],
				 origin[2] +(0.25 +c*0.5)*Blk_Length_of[RLevel][2]};



	//! Try to refine the 7 (=2^3-1) Octs around this Block
	//! IN THIS LEVEL (one below oct to refine !!!)


	//! Start:
	//! a=0 && is  m1BB -> start at  0
	//! a=0 && not m1BB -> start at -1
	//! a=1 && any m1BB -> start at  0
	//! c,b respectively
	INT32 start[3] = {!a * !is_box_boundary[0],
			 !b * !is_box_boundary[2],
			 !c * !is_box_boundary[4]};

	//! End:
	//! a=1 && is  p1BB -> end at  0
	//! a=1 && not p1BB -> end at  1
	//! a=0 && any p1BB -> end at  0
	//! c,b respectively
	INT32   end[3] = {a * !is_box_boundary[1],
			 b * !is_box_boundary[3],
			 c * !is_box_boundary[5]};


	//! refine Neighbours around respective oct
	for(INT32 u=-start[0]; u<=end[0]; u++)
	 for(INT32 v=-start[1]; v<=end[1]; v++)
	  for(INT32 w=-start[2]; w<=end[2]; w++)
	  if(u||v||w)
	  {


		D_REAL position[3] = {oct_centre[0] +u*Blk_Length_of[RLevel][0],
				      oct_centre[1] +v*Blk_Length_of[RLevel][1],
				      oct_centre[2] +w*Blk_Length_of[RLevel][2]};



		refine_Oct_at(RLevel, position, false);

	  }

}



//!-------------------------------------------------------------//
//! refine_Blocks: - call by reference function			//
//!-------------------------------------------------------------//
void CBlock::refine_gatherEnvironment(void)
{


	//! never refine gather Blocks
	if(is_gatherBlk)
	{
		log_file << " WARNING:" << endl
		         << " tried to refine gather block in further gather blocks."
		         << " Ignoring command and continuing ... " << endl;
		return;
	}



	//! refine "min" Neighbours around this block
	//! (8-1=7 in total, leave out 0,0,0 block)
// 	for(INT32 u=0; u<=!is_box_boundary[0]; u++)
// 	 for(INT32 v=0; v<=!is_box_boundary[2]; v++)
// 	  for(INT32 w=0; w<=!is_box_boundary[4]; w++)
	for(INT32 u=-!is_box_boundary[1]; u<=!is_box_boundary[0]; u++)
	 for(INT32 v=-!is_box_boundary[3]; v<=!is_box_boundary[2]; v++)
	  for(INT32 w=-!is_box_boundary[5]; w<=!is_box_boundary[4]; w++)
	  if(u||v||w)
	  {


		//! NOTE:
		//! 0.1 is required to shift point 50% of block's length
		//! into block
		//! IF THIS IS NOT DONE ROUND OFF ERRORS WILL OCCUR !!!
		D_REAL position[3] = {origin[0] +(0.5 -u)*Blk_Length_of[RLevel][0],
				      origin[1] +(0.5 -v)*Blk_Length_of[RLevel][1],
				      origin[2] +(0.5 -w)*Blk_Length_of[RLevel][2]};


		refine_Oct_at(RLevel, position, true);

	  }



}



//!-------------------------------------------------------------//
//! refine_Oct_at:								//
//!-------------------------------------------------------------//
void CBlock::refine_Oct_at(INT32 req_level, D_REAL *position, bool ref_in_gatherOct)
{

	CBlock *temp_Block, *block_to_refine;
	INT32   existing_level, blk_index, *index_blk[3];


	//! alloc memory to store indices for each component/level
	for(INT32 comp=0; comp<3; comp++)
	{
 		index_blk[comp] = new INT32[MAX_LEVEL+1];
 		memset(index_blk[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
	}


	//!--------------------------------------------------
	//! 1) calculate block indices and check whether
	//!    position is in Simu Box.
	//!--------------------------------------------------
	//! 
	if(!calc_BlockIndex_of_Pos(index_blk, position))
	{
		log_file << " ERROR: " << endl
		 	 << " Block to refine out of Box: " << endl;


		for(INT32 level=0; level<=MAX_LEVEL; level++)
		log_file << " L"<<level<< " index_blk[0]: " << index_blk[0][level]
			 	       << " index_blk[1]: " << index_blk[1][level]
			 	       << " index_blk[2]: " << index_blk[2][level]
			 	       << endl;

		log_file << " Exiting ... " << endl;
		exit(0);

	}

	
	//!--------------------------------------------------
	//! 2) climb to highest possible level at respective
	//!    position.
	//!--------------------------------------------------


	INT32 block_indices[3] = {index_blk[0][0], index_blk[1][0], index_blk[2][0]};

	//! first calculate root block
	//! -> always level 0
	if(use_SFC)
	blk_index = SFC_Indices_to_BlkNr(block_indices, 0);
	else
	blk_index = LINEAR_Indices_to_BlkNr(block_indices);


	//! initialize to respective root block
	temp_Block      = Root_Block_Array +blk_index;
	block_to_refine = Root_Block_Array +blk_index;;

	//! descent to highest existing level
	for(INT32 level=1; level<=MAX_LEVEL; level++)
	{

		//! set block index to precalculated indices
		blk_index = index_blk[0][level] *2 *2
			   +index_blk[1][level] *2
			   +index_blk[2][level];
		

		//! Decent by gather child array, else
		//! gather octs are not visible
		//! They are removed before ordanary blocks are refined
		//! BUT THEY MIGHT BE CREATED SHORT BEFORE WITH IN "create_gather_block"
		//! function. SO GATHER CHILDS MUST BE USED HERE !!!
		temp_Block = temp_Block->gather_child_array[blk_index];

		//! cancel in case block does not exist
		if(!temp_Block)
		break;

		//! Set block_to_refine to temp_Block
		block_to_refine = temp_Block;

		
	}

	existing_level = block_to_refine->RLevel;


	//!--------------------------------------------------
	//! 3a) check whether estimated level is to small
	//!--------------------------------------------------
	if(existing_level < req_level-1)
	{
		
// 		log_file << " ERROR: " << endl
// 			 << " requested level is to small: " << endl
// 			 << " requested: "  << req_level << endl
// 			 << " estimated: "  << existing_level << endl
// 			 << " Exiting ... " << endl;

		return;
// 		exit(0);
	}

	//!-----------------------------------------------------
	//! 3a) check whether estimated level is equal or larger
	//!     -> block already does exist
	//!     -> do nothing
	//!-----------------------------------------------------
// 	log_file << "req_level: " << req_level << endl;
// 	log_file << "existing_level: " << existing_level << endl;
	if(existing_level >= req_level)
	{

		//! delete for each component/level
		for(INT32 comp=0; comp<3; comp++)
		delete[] index_blk[comp];

		return;
	}


	//!--------------------------------------------------
	//! 4) refine Block
	//!    due to interceptions above:
	//!    block_to_refine->RLevel == req_level
	//!--------------------------------------------------
	blk_index = index_blk[0][req_level] *2 *2
		   +index_blk[1][req_level] *2
		   +index_blk[2][req_level];


	if(ref_in_gatherOct)
	block_to_refine->refine_gatherOct(blk_index);
	else
	if(block_to_refine->is_refinement_permitted(blk_index))
	block_to_refine->refine_Oct(blk_index);



	//! delete for each component/level
	for(INT32 comp=0; comp<3; comp++)
 	delete[] index_blk[comp];


}

//!-------------------------------------------------------------//
//! seek_Block_at:								//
//! find Block at given Level and Position 				//
//!-------------------------------------------------------------//
CBlock* CBlock::seek_Block_at(INT32 req_level, D_REAL *position, INT32 &child_id_at_position)
{

	CBlock *temp_Block, *block_to_find;
	INT32   blk_index, *index_blk[3];


	//! alloc memory to store indices for each component/level
	for(INT32 comp=0; comp<3; comp++)
	{
 		index_blk[comp] = new INT32[MAX_LEVEL+1];
 		memset(index_blk[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
	}


	//!--------------------------------------------------
	//! 1) calculate block indices and check whether
	//!    position is in Simu Box.
	//!  - in case block is out of bocks return false
	//!    since no block available
	//!--------------------------------------------------
	if(!calc_BlockIndex_of_Pos(index_blk, position))
	{

		//! delete for each component/level
		for(INT32 comp=0; comp<3; comp++)
		delete[] index_blk[comp];


		log_file << " ERROR: " << endl
		 	 << " Block to seek out of Box: " << endl;


		for(INT32 level=0; level<=MAX_LEVEL; level++)
		log_file << " index_blk[0]["<<level<<"]: " << index_blk[0][level]
			 << " index_blk[1]["<<level<<"]: " << index_blk[1][level]
			 << " index_blk[2]["<<level<<"]: " << index_blk[2][level]
			 << endl;

		log_file << " Exiting ... " << endl;
		exit(0);
	}


	//!--------------------------------------------------
	//! 2) climb to highest possible level at respective
	//!    position.
	//!--------------------------------------------------
	INT32 block_indices[3] = {index_blk[0][0], index_blk[1][0], index_blk[2][0]};

	//! first calculate root block
	//! -> always level 0
	if(use_SFC)
	blk_index = SFC_Indices_to_BlkNr(block_indices, 0);
	else
	blk_index = LINEAR_Indices_to_BlkNr(block_indices);


	//! initialize to respective root block
	temp_Block    = Root_Block_Array +blk_index;
	block_to_find = Root_Block_Array +blk_index;

	//! descent to req_level if available
	for(INT32 level=1; level<=req_level; level++)
	{

		//! set block index to precalculated indices
		blk_index = index_blk[0][level] *2 *2
			   +index_blk[1][level] *2
			   +index_blk[2][level];
		

		//! Climp up by gather child array, else
		//! gather octs are not visible
		//! They are removed before ordanary blocks are refined
		//! BUT THEY MIGHT BE CREATED SHORT BEFORE WITHIN "create_gather_block"
		//! function. SO GATHER CHILDS MUST BE USED HERE !!!
		temp_Block = temp_Block->gather_child_array[blk_index];

		//! cancel in case block does not exist
		if(!temp_Block)
		break;

		//! Set block_to_find to temp_Block
		block_to_find = temp_Block;

		
	}


	//! RETURN BY REFERENCE child_id_at_position
	//! at level req_level+1
	if(block_to_find->RLevel<MAX_LEVEL)
	child_id_at_position = index_blk[0][block_to_find->RLevel+1] *2 *2
			      +index_blk[1][block_to_find->RLevel+1] *2
			      +index_blk[2][block_to_find->RLevel+1];
	else
	child_id_at_position = -1;


	//! delete for each component/level
	for(INT32 comp=0; comp<3; comp++)
 	delete[] index_blk[comp];
		

	//! return Block which is closed to requested level
	return block_to_find;


}


//!-------------------------------------------------------------//
//! refine_gatherly							//
//!-------------------------------------------------------------//
void CBlock::refine_gatherOct(INT32 id_oct)
{

	if(gather_child_array[id_oct])
	{
		log_file << " ERROR in refine_gatherOct: "<< endl
			<< " tried to refine child that already existed" << endl;
		exit(0);
// 		return;
	}

	//! get index of oct (0 or 1 for each dimension)
	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);

	INT32 indices[3] = {a,b,c};

	//! only gather_child_array is set.
	//! child_array remains unset.
	//! -> No gather Block is part of ordanary BlockList
	//! -> gather-Block-Neighbours are INVISIBLE for ordanary Blocks
	//! -> ordanary-Block-Neighbours are VISIBLE for gather Blocks
	gather_child_array[id_oct] = new CBlock;

	gather_child_array[id_oct]->init_Block(0, indices, RLevel+1, this);
	gather_child_array[id_oct]->alloc_GatherMemory();
	gather_child_array[id_oct]->is_gatherBlk = true;



	//! first all Blks have to be marked as gather before
	//! neighbour search

	//! REMEMBER:
	//! - gather_Neighbours of ordanary Blks are also set, hence
	//!   all ordanry Blks take part in the neighbour search !

	//! search includes ALL Blocks (ordanary and gather)
  	gather_child_array[id_oct]->find_all_gather_Neighbours();

	//! only gather - link is set
  	gather_child_array[id_oct]->link_all_gather_Neighbours();


	//! REMEMBER:
	//! - No gather Block is part of ordanary BlockList
	//! - GATHER_BlockList contains all ordanary Blks
	add_Block_to_GthrBlkList(id_oct);

	num_refined_Octs++;


}


//!-------------------------------------------------------------//
//! refine_Blk: - call by reference function 			//
//! 			- allocs memory and integrates the array in	//
//!			  the refninement level's linear list		//
//! 			- allocs fields and initializes coordinates	//
//!-------------------------------------------------------------//
void CBlock::refine_Oct(INT32 id_oct)
{

	//! use parents mpi process by default
	refine_Oct_set_mpi_process(id_oct,/* this->*/responsible_mpi_process);



}


//!-------------------------------------------------------------//
//! refine_Oct_skip_MempryAlloc:						//
//!-------------------------------------------------------------//
void CBlock::refine_Oct_set_mpi_process(INT32 id_oct, INT32 mpi_process)
{



	//! get index of oct (0 or 1 for each dimension)
	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);

	INT32 indices[3] = {a,b,c};



	if(child_array[id_oct])
	{
		log_file << " ERROR in refine_Oct: "<< endl
			<< " tried to refine child that already existed" << endl;
		exit(0);
	}


	//! increase global Block counter
	total_Blocks_L[RLevel+1]++;
	total_active_Blocks++;

	//! gather_child_array is ALWAYS set
	//! -> Every ordanary Block also is part of GATHER_BlockList.
	//! -> In turn no gather Block is part of ordanary BlockList.
	 
	//! -> gather-Block-Neighbours are INVISIBLE for ordanary Blocks
	//! -> ordanary-Block-Neighbours are VISIBLE for gather Blocks
	child_array[id_oct] = new CBlock;

	//! NOTE: DO NOT CALL MEMSET BELOW. IT WILL CRASH WHEN MEMORY IS DEALLOCATED
	//! 	delete child_array[id_child];
// 	memset(child_array[id_oct], 0, sizeof(CBlock));

	//! gather_child_array must be set, important eg. for neighbour search
	//! and add to GatherBlkList
	gather_child_array[id_oct] = child_array[id_oct];

	//! initialize oct
	child_array[id_oct]->init_Block(0, indices, RLevel+1,this);


	//! search includes only ordanary Blocks
  	child_array[id_oct]->find_all_ordanary_Neighbours();

	//! only ordanry links are set
  	child_array[id_oct]->link_all_ordanary_Neighbours();

	//! GATHER Neighbours have to be set , otw. density is 
	//! not exchanged between ANY Blk.
	//! search includes ALL Blocks (ordanary and gather)
	//! TODO: the above search is osolete, in case 
	//! 	    the search distinguishes between ord and gather 
	//! 	    and sets the respective Neighbour
  	child_array[id_oct]->find_all_gather_Neighbours();

	//! only gather - link is set
  	child_array[id_oct]->link_all_gather_Neighbours();

	//! FOR CLARIFICATION:
	//! - Every ordanary Block also is part of GATHER_BlockList.
	//! - No gather Block is part of ordanary BlockList
	add_Block_to_BlkList(id_oct);
	add_Block_to_GthrBlkList(id_oct);

	//! set to specified mpi process before memory allocation:
	child_array[id_oct]->responsible_mpi_process = mpi_process;
// 	child_array[id_oct]->responsible_mpi_process = RLevel+1;




	//! as usual exclusively alloc memory in case this procss is responsible
	if(mpi_myRank==child_array[id_oct]->responsible_mpi_process)
	alloc_set_OctFields(id_oct);


	//! do not send particle in first TL, since not available
	//! (what in case initial refinement is added?)
//	if(TL && TL-TL_at_last_restore_state>0)
	{
		//! create package with particle for child array and send via MPI if required
		//! in case my_rank processes THIS block
		if(mpi_myRank==responsible_mpi_process)
		send_particle_toChild(id_oct);
	
		//! insert package with particle from parent, recv via MPI if required
		//! in case my_rank processes CHILD block
		if(mpi_myRank==child_array[id_oct]->responsible_mpi_process)
		child_array[id_oct]->recv_particle_fromParent();
	}


	num_refined_Octs++;
	num_children++;

}


//!-------------------------------------------------------------//
//! alloc_set_OctFields:
//!-------------------------------------------------------------//
void CBlock::alloc_set_OctFields(INT32 id_oct)
{


	//! Alloc Memory and set fields
	child_array[id_oct]->alloc_process_specific_Memory();
  	child_array[id_oct]->set_Fields();


	//! used for field_from_parent
	INT32 start_ijk[3]={          0,          0,          0};
	INT32   end_ijk[3]={ BlkNds_X/2, BlkNds_Y/2, BlkNds_Z/2};

	//! - GN are included and need no update
	//! - Moments will be collected after refinement, no update required
	//! - do not interpolate in TL=0, use set fields instead
	if(TL /*&& TL-TL_at_last_restore_state>0*/)
	 //! EXCLUSIVELY INTERPOLATE NEWLY FIELDS FROM PARENT IF:
	 //!  --- PARENT-Process == CHILD-Process == mpi_myRank ----
	 //! -> THIS IS TRUE FOR ALL NEWLY GENERATED BLOCKS
	 //! -> CRASH OTHERWISE IN RESTORE
	 if(child_array[id_oct]->responsible_mpi_process == responsible_mpi_process)
	 {
		//! do not interpolate BTotal but rebuild BTotal instead at the end of refinement
		//! -> can cause serious issues in case of large dipolse else !!!
// 		child_array[id_oct]->field_from_parent(id_BTotal, id_BTotal, start_ijk, end_ijk);
// 		child_array[id_oct]->field_from_parent(id_Bcfbg,  id_Bcfbg,  start_ijk, end_ijk);
		child_array[id_oct]->field_from_parent(id_BEven,  id_BEven,  start_ijk, end_ijk);
	 } 

	//! set extern fields if available
// 	for(INT32 extern_rho=0; extern_rho<num_externRhoVelocityFields; extern_rho++)
// 	{
// 
// 		child_array[id_oct]->field_from_parent(id_externRho1 +extern_rho,
// 						       id_externRho1 +extern_rho,
// 						       start_ijk, end_ijk);
// 
// 		child_array[id_oct]->field_from_parent(id_extern_Ui1 +extern_rho,
// 						       id_extern_Ui1 +extern_rho,
// 						       start_ijk, end_ijk);
// 
// 	}



	
}




//!-------------------------------------------------------------//
//! remove_gather_Blk							//
//!-------------------------------------------------------------//
void CBlock::remove_gather_Blk(INT32 id_gChild)
{

	//! THIS HAS TO BE INTERCEPTED SINCE ORDANARY BLOCKS
	//! ARE LISTED IN gather_child_array ALSO !!!
	if(!gather_child_array[id_gChild]->is_gatherBlk)
	return;


	//! unlink all neighbours:
	gather_child_array[id_gChild]->unlink_all_Neighbours();

	//! GatherMemory is the only dynamically 
	//! allocated Memory of a Gather Block
	if(responsible_mpi_process == mpi_myRank)
	gather_child_array[id_gChild]->delete_GatherMemory();


	remove_GthrBlock_from_GthrBlkList(id_gChild);

	num_refined_Octs--;
	delete gather_child_array[id_gChild];

	gather_child_array[id_gChild] = 0;
	child_array[id_gChild] = 0;

	num_removed_gatherBlks+=1;

}


//!-------------------------------------------------------------//
//! remove_Child							//
//!-------------------------------------------------------------//
void CBlock::remove_Child(INT32 id_child)
{


	//! Do this before memory is deleted
	if(mpi_myRank==child_array[id_child]->responsible_mpi_process)
	child_array[id_child]->prepare_pPackage_toParent();

	//! insert particle of child
	if(mpi_myRank==responsible_mpi_process)
	insert_pPackage_fromChild(id_child);

	//! Unlink all neighbours:
	child_array[id_child]->unlink_all_Neighbours();

	//! Delete dynamically allocated Memory
	child_array[id_child]->delete_Memory();



	remove_Block_from_BlkList(id_child);
	remove_GthrBlock_from_GthrBlkList(id_child);

	total_Blocks_L[RLevel+1]--;
	total_active_Blocks--;

	num_children--;

	delete child_array[id_child];

	gather_child_array[id_child] = NULL;
	child_array[id_child] = NULL;



}


//!-------------------------------------------------------------//
//! remove_Child							//
//!-------------------------------------------------------------//
void CBlock::remove_GthrBlock_from_GthrBlkList(INT32 id_gChild)
{



	//! Unhook GATHERBlkArray of GATHERList
	if(gather_child_array[id_gChild]->prev_Blk_of_GatherBlockList)
	   gather_child_array[id_gChild]->prev_Blk_of_GatherBlockList->next_Blk_of_GatherBlockList
	=  gather_child_array[id_gChild]->next_Blk_of_GatherBlockList;
	else
	//! in case previous does not exist, this Array was first of List 
	//! -> set next to List start
	GATHER_BlockList_of_Lev[RLevel+1] =gather_child_array[id_gChild]->next_Blk_of_GatherBlockList;

	if(gather_child_array[id_gChild]->next_Blk_of_GatherBlockList)
	  gather_child_array[id_gChild]->next_Blk_of_GatherBlockList->prev_Blk_of_GatherBlockList
	= gather_child_array[id_gChild]->prev_Blk_of_GatherBlockList;


}


//!-------------------------------------------------------------//
//! remove_Child							//
//!-------------------------------------------------------------//
void CBlock::remove_Block_from_BlkList(INT32 id_child)
{


	//! Unhook BlkArray of List
	if(child_array[id_child]->prev_Blk_of_BlockList)
	  child_array[id_child]->prev_Blk_of_BlockList->next_Blk_of_BlockList
	= child_array[id_child]->next_Blk_of_BlockList;
	else
	//! in case previous does not exist, this Array was first of List 
	//! -> set next to List start
	BlockList_of_Lev[RLevel+1] = child_array[id_child]->next_Blk_of_BlockList;

	if(child_array[id_child]->next_Blk_of_BlockList)
	   child_array[id_child]->next_Blk_of_BlockList->prev_Blk_of_BlockList
	=  child_array[id_child]->prev_Blk_of_BlockList;

}






//!-------------------------------------------------------------//
//! add_Block_to_BlkList:							//
//!-------------------------------------------------------------//
void CBlock::add_Block_to_BlkList(INT32 id_oct)
{


	//! find correct rb indices:
	CBlock* temp = 0;
	CBlock* prev_of_temp = 0;
	//! NOTE: Just for clarification:
	//!	    - child_array is the same as &child_array[0]
	//!	    - child_array->next_Block_Array is the same as 
	//!	      child_array[0].next_Block_Array
	//!	 -> changes of linear level list are done in the first block
	//!	 -> prev_Block_Array & next_Block_Array of children 2-8 remain
	//!	    unchanged. (-> few bytes of memory are wasted).


	temp = BlockList_of_Lev[RLevel+1];
	//! - Is there at all a Block on this level ?
	//! - If yes is its number higher than the one of new BlkArray ?
	if(!temp || temp->Block_Nr > child_array[id_oct]->Block_Nr)
	{
		BlockList_of_Lev[RLevel+1] = child_array[id_oct];
		child_array[id_oct]->next_Blk_of_BlockList = temp;
		child_array[id_oct]->prev_Blk_of_BlockList = 0;

		if(temp)
		temp->prev_Blk_of_BlockList = child_array[id_oct];
	}
	else
	{


		//! Up to now:
		//! 1) BlockList_of_Lev NOT EMPTY
		//! 2) Block_Nr of child_array "<" present array
		//!	 -> first loop step will ALWAYS be performed
		//!	 -> temp & prev_of_temp get initialized !!!
		//! -> loop until next BA Nr. is higher
		while(    child_array[id_oct]->Block_Nr > temp->Block_Nr
			&& temp->next_Blk_of_BlockList)
		{
			prev_of_temp = temp;
			temp = temp->next_Blk_of_BlockList;

		}

		//! loop ended cause child_array Nr. "<" temp Nr. 
		if( child_array[id_oct]->Block_Nr < temp->Block_Nr)
		{
			//! inhook:
			//! prev <-> child_array <-> temp
			prev_of_temp->next_Blk_of_BlockList = child_array[id_oct];
			temp->prev_Blk_of_BlockList = child_array[id_oct];

			child_array[id_oct]->prev_Blk_of_BlockList = prev_of_temp;
			child_array[id_oct]->next_Blk_of_BlockList = temp;
		}
		//! loop ended cause temp->next_Blk_of_BlockList = 0
		//! -> child_array Nr. > temp Nr.
		else
		{
			child_array[id_oct]->next_Blk_of_BlockList = 0;
			child_array[id_oct]->prev_Blk_of_BlockList = temp;
			temp->next_Blk_of_BlockList = child_array[id_oct];
		}
	}

}

//!-------------------------------------------------------------//
//! add_Block_to_BlkList:						//
//!-------------------------------------------------------------//
void CBlock::add_Block_to_GthrBlkList(INT32 id_oct)
{

	//! find correct rb indices:
	CBlock* temp = 0;
	CBlock* prev_of_temp = 0;
	//! NOTE: Just for clarification:
	//!	    - child_array is the same as &child_array[0]
	//!	    - child_array->next_Block_Array is the same as 
	//!	      child_array[0].next_Block_Array
	//!	 -> changes of linear level list are done in the first block
	//!	 -> prev_Block_Array & next_Block_Array of children 2-8 remain
	//!	    unchanged. (-> few bytes of memory are wasted).

	temp = GATHER_BlockList_of_Lev[RLevel+1];
	//! - Is there at all a Block on this level ?
	//! - If yes is its number higher than the one of new BlkArray ?
	if(  !temp || temp->Block_Nr > gather_child_array[id_oct]->Block_Nr)
	{
		GATHER_BlockList_of_Lev[RLevel+1] = gather_child_array[id_oct];
		gather_child_array[id_oct]->next_Blk_of_GatherBlockList = temp;
		gather_child_array[id_oct]->prev_Blk_of_GatherBlockList = 0;

		if(temp)
		temp->prev_Blk_of_GatherBlockList = child_array[id_oct];
	}
	else
	{


		//! Up to now:
		//! 1) BlockList_of_Lev NOT EMPTY
		//! 2) Block_Nr of child_array "<" present array
		//!	 -> first loop step will ALWAYS be performed
		//!	 -> temp & prev_of_temp get initialized !!!
		//! -> loop until next BA Nr. is higher
		while(    gather_child_array[id_oct]->Block_Nr > temp->Block_Nr
			&& temp->next_Blk_of_GatherBlockList)
		{
			prev_of_temp = temp;
			temp = temp->next_Blk_of_GatherBlockList;

		}

		//! loop ended cause child_array Nr. "<" temp Nr.
		if( gather_child_array[id_oct]->Block_Nr < temp->Block_Nr)
		{

			//! inhook:
			//! prev <-> gather_child_array <-> temp
			prev_of_temp->next_Blk_of_GatherBlockList = gather_child_array[id_oct];
			temp->prev_Blk_of_GatherBlockList = gather_child_array[id_oct];

			gather_child_array[id_oct]->prev_Blk_of_GatherBlockList = prev_of_temp;
			gather_child_array[id_oct]->next_Blk_of_GatherBlockList = temp;
		}
		//! loop ended cause temp->next_Blk_of_BlockList = 0
		//! -> child_array Nr. > temp Nr.
		else
		{
			gather_child_array[id_oct]->next_Blk_of_GatherBlockList = 0;
			gather_child_array[id_oct]->prev_Blk_of_GatherBlockList = temp;
			temp->next_Blk_of_GatherBlockList = gather_child_array[id_oct];
		}

	}

}

/*
//!-------------------------------------------------------------//
//! estimate_refinement_efficiency:						//
//!-------------------------------------------------------------//
void CBlock::count_flagged_nodes(INT32 src_type, F_REAL *max_Value, INT32* num_flagged_nodes)
{

	INT32 i_j_k;

	D_REAL* AbsField = Field_Type[src_type];

	D_REAL  low_barrier = removeBlk_if_rating_below_L[RLevel] *max_Value[RLevel];
	D_REAL high_barrier = refineBlk_if_rating_above_L[RLevel] *max_Value[RLevel];


	//! check whether max_Value is set in this RLevel
	if(max_Value[RLevel]<0.)
	return;


	const INT32 i_offset = (BlkNds_X-2)/2;
	const INT32 j_offset = (BlkNds_Y-2)/2;
	const INT32 k_offset = (BlkNds_Z-2)/2;


	for(INT32 a = 0; a < 2; a++)
	 for(INT32 b = 0; b < 2; b++)
	  for(INT32 c = 0; c < 2; c++)
	  {

		//! define entire volume of oct (NO GC, BUT INSIDE PHYSICAL SPACE)
		INT32 start_ijk[3] = {(a+0)*i_offset+1, (b+0)*j_offset+1, (c+0)*k_offset+1};
		INT32 end_ijk[3] =   {(a+1)*i_offset+1, (b+1)*j_offset+1, (c+1)*k_offset+1};

		INT32 id_oct = a*2*2 +b*2 +c;

		//! if child_array exist, node counts as flagged unless its Abs is
		//! lager that low_barrier

		for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
		 for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
		  for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
		  {

			i_j_k =  i*BlkNds_Y*BlkNds_Z
				+j*BlkNds_Z
				+k;

			//! if child_array exist, node counts as flagged unless its Abs is
			//! lager that low_barrier
			if(child_array[id_oct] && AbsField[i_j_k] >= low_barrier)
			num_flagged_nodes[RLevel]++;

			//! if child_array does NOT exist, node counts as flagged unless its Abs is
			//! lager that high_barrier
			if(!child_array[id_oct] && AbsField[i_j_k] >= high_barrier)
			num_flagged_nodes[RLevel]++;
			

		}

	  }
	
}*/


//!-------------------------------------------------------------//
//! set_average_ref_value:							//
//! FOR 3 COMPONENTS ONLY							//
//!-------------------------------------------------------------//
void CBlock::set_average_ref_value(INT32 id_criteria)
{




	//! pointer to each field component
	D_REAL *Field[3];

	//!
	INT32 i_j_k;
	INT32 id_field = refcrit_field_IDs[id_criteria];
	INT32 num_field_comps = COMPs_FType[id_field];
	INT32 criteria_comp = refcrit_field_comp[id_criteria];


	//! distinguish scalar vs. vector field
	if(num_field_comps == 1)
	Field[0] = Field_Type[id_field];

	if(num_field_comps == 3)
	{
		Field[0] = Field_Type[id_field]+0*num_nodes_in_block;
		Field[1] = Field_Type[id_field]+1*num_nodes_in_block;
		Field[2] = Field_Type[id_field]+2*num_nodes_in_block;
	}



	const INT32 i_offset = (BlkNds_X-2)/2;
	const INT32 j_offset = (BlkNds_Y-2)/2;
	const INT32 k_offset = (BlkNds_Z-2)/2;


	//! loop accros 8 octs
	for(INT32 a = 0; a < 2; a++)
	 for(INT32 b = 0; b < 2; b++)
	  for(INT32 c = 0; c < 2; c++)
	  {

		//! define entire volume of oct (NO GC, BUT INSIDE PHYSICAL SPACE)
		INT32 start_ijk[3] = {(a+0)*i_offset+1, (b+0)*j_offset+1, (c+0)*k_offset+1};
		INT32 end_ijk[3] =   {(a+1)*i_offset+1, (b+1)*j_offset+1, (c+1)*k_offset+1};

		INT32 id_oct = a*2*2 +b*2 +c;

		//! reset average
                D_REAL average_value = 0.;


		//! loop accros nodes in oct
		for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
		 for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
		  for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
		  {
		
			i_j_k =  i*BlkNds_Y*BlkNds_Z
				+j*BlkNds_Z
				+k;

			if(criteria_comp < 3 || num_field_comps==1)
			average_value += Field[criteria_comp][i_j_k];
			else
			average_value += sqrt( Field[0][i_j_k]*Field[0][i_j_k]
					      +Field[1][i_j_k]*Field[1][i_j_k]
					      +Field[2][i_j_k]*Field[2][i_j_k]);
				
		  }


		  INT32 num_nodes_in_oct = (end_ijk[0] -start_ijk[0])
				         +(end_ijk[1] -start_ijk[1])
				         +(end_ijk[2] -start_ijk[2]);

		  average_value/=num_nodes_in_oct;
                  average_ref_value[id_oct] = average_value;


	  }
	

}
