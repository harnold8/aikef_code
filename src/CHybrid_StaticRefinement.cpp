

#include <math.h>
#include "CHybrid.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>



//!------------------------------------------------------------
//!- static_refinement_sphere: Generate initial Childs 
//!				- eg. near the obstacle
//!------------------------------------------------------------
void CHybrid::static_refinement_cuboid(void)
{




	D_REAL x_[3];
	INT64 num_blocks_createdL[MAX_LEVEL+2];
	memset(num_blocks_createdL,0,(MAX_LEVEL+2)*sizeof(INT64));



	log_file << " Creating children cuboid ...   " << endl;


	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];

		if(TL_REFINE_STATIC_CUBOID[level] >= 0 && TL_REFINE_STATIC_CUBOID[level] <= TL)
		while(temp_Block)
		{

			for(INT32 a=0; a<2; a++)
			 for(INT32 b=0; b<2; b++)
			  for(INT32 c=0; c<2; c++)
			  {


				x_[0] = temp_Block->origin[0]+(0.25 +0.5*a)*Blk_Length_of[level][0];
				x_[1] = temp_Block->origin[1]+(0.25 +0.5*b)*Blk_Length_of[level][1];
				x_[2] = temp_Block->origin[2]+(0.25 +0.5*c)*Blk_Length_of[level][2];
	
				INT32 oct = 2*2*a +2*b +c;

				//! ---- REFINE LEVEL [level] -----------------
				if(    x_[0] > minX_refBox_of_L[level] && x_[0] < maxX_refBox_of_L[level]
				    && x_[1] > minY_refBox_of_L[level] && x_[1] < maxY_refBox_of_L[level]
				    && x_[2] > minZ_refBox_of_L[level] && x_[2] < maxZ_refBox_of_L[level]
				)
				{
					//! in case state was restored existing child has to be marked
					//! as initial refined block
					if(temp_Block->child_array[oct])
					temp_Block->child_array[oct]->initial_refined = true;
					else if(temp_Block->is_refinement_permitted(oct))
					 {

						temp_Block->refine_Oct(oct);
						temp_Block->child_array[oct]->initial_refined = true;
			
						num_blocks_createdL[level+1]++;
						total_Blocks_L[level+1]++;
					 }
				}
			  }

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}


	log_file << " done." << endl;
	log_file << " Num blocks in ROOT Level: " << num_root_blocks << endl;
	INT32 j = 0;
	for(j=1; j<= MAX_LEVEL; j++)num_blocks_createdL[MAX_LEVEL+1] += num_blocks_createdL[j];
	for(j=1; j<= MAX_LEVEL; j++)
	log_file << " Num blocks created in Level "<<j<<": " << num_blocks_createdL[j] << endl;
	log_file << " Total child blocks created: " << num_blocks_createdL[MAX_LEVEL+1] << endl;
	log_file << " Total active Blocks: " << total_active_Blocks  << endl;


// 	log_file << endl;
// 	log_file << "Block File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*sizeof(CBlock*)+2 << " byte"  << endl;
// 	log_file << "BField File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*3*num_nodes_in_block+2 << " byte"  << endl;


}



//!------------------------------------------------------------
//!- static_refinement_sphere: Generate initial Childs 
//!				- eg. near the obstacle
//!------------------------------------------------------------
void CHybrid::static_refinement_sphere(void)
{




	D_REAL x_[3],radius;
	log_file << " Creating children sphere  ...   " << endl;


	INT64 num_blocks_createdL[MAX_LEVEL+2];
	memset(num_blocks_createdL,0,(MAX_LEVEL+2)*sizeof(INT64));


	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];

		if(TL_REFINE_STATIC_SPHERE[level] >=0 && TL_REFINE_STATIC_SPHERE[level] <= TL)
		while(temp_Block)
		{




			//! ---- REFINE LEVEL [level] -----------------
			//! in case state was restored existing child has to be markes
			//! as initial refined block

			for(INT32 a=0; a<2; a++)
			 for(INT32 b=0; b<2; b++)
			  for(INT32 c=0; c<2; c++)
			  {

				x_[0] = temp_Block->origin[0]+(0.25 +0.5*a)*Blk_Length_of[level][0];
				x_[1] = temp_Block->origin[1]+(0.25 +0.5*b)*Blk_Length_of[level][1];
				x_[2] = temp_Block->origin[2]+(0.25 +0.5*c)*Blk_Length_of[level][2];


				radius = vec_len(x_);

				INT32 oct = 2*2*a +2*b +c;

				if(radius < radius_refSphere_of_L[level]
				   && radius > SPHERE_INNER_BOUNDARY[level] * R_Obstacle
					/*&& (level <MAX_LEVEL-1 ||  (radius > 0.8*obstacle_core_fraction*R_Obstacle 
										&& (&& fabs(x_[2])<R_Obstacle) )*/
				    )
				{

					if(temp_Block->child_array[oct])
                                	temp_Block->child_array[oct]->initial_refined = true;

				  	else if(temp_Block->is_refinement_permitted(oct))
				  	{
						temp_Block->refine_Oct(oct);
						temp_Block->child_array[oct]->initial_refined = true;
		
						num_blocks_createdL[level+1]++;
				  	}

				}
			  }



			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}




	log_file << " done." << endl;
	log_file << " Num blocks in ROOT Level: " << num_root_blocks << endl;
	INT32 j = 0;
	for(j=1; j<= MAX_LEVEL; j++)num_blocks_createdL[MAX_LEVEL+1] += num_blocks_createdL[j];
	for(j=1; j<= MAX_LEVEL; j++)
	log_file << " Num blocks created in Level "<<j<<": " << num_blocks_createdL[j] << endl;
	log_file << " Total child blocks created: " << num_blocks_createdL[MAX_LEVEL+1] << endl;
	log_file << " Total blocks (including. ROOT): " << total_active_Blocks << endl;





	log_file << endl;
// 	log_file << "Block File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*sizeof(CBlock*)+2 << " byte"  << endl;
// 	log_file << "BField File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*3*num_nodes_in_block+2 << " byte"  << endl;


}


//!------------------------------------------------------------
//!- static_refinement_sphere: Generate initial Childs 
//!				- eg. near the obstacle
//!------------------------------------------------------------
void CHybrid::static_refinement_ZylinderX(void)
{




	D_REAL x_[3],radius;
	log_file << " Creating children Zylinder X-Axis  ...   " << endl;


	INT64 num_blocks_createdL[MAX_LEVEL+2];
	memset(num_blocks_createdL,0,(MAX_LEVEL+2)*sizeof(INT64));


	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];

		if(TL_REFINE_STATIC_ZYLINDER_X[level] >=0 && TL_REFINE_STATIC_ZYLINDER_X[level] <= TL)
		while(temp_Block)
		{




			//! ---- REFINE LEVEL [level] -----------------
			//! in case state was restored existing child has to be markes
			//! as initial refined block

			for(INT32 a=0; a<2; a++)
			 for(INT32 b=0; b<2; b++)
			  for(INT32 c=0; c<2; c++)
			  {

				x_[0] = temp_Block->origin[0]+(0.25 +0.5*a)*Blk_Length_of[level][0];
				x_[1] = temp_Block->origin[1]+(0.25 +0.5*b)*Blk_Length_of[level][1];
				x_[2] = temp_Block->origin[2]+(0.25 +0.5*c)*Blk_Length_of[level][2];


				radius = sqrt( x_[1]*x_[1] + x_[2]*x_[2]);

				INT32 oct = 2*2*a +2*b +c;

				if(      x_[0]>=minX_refZylinder_of_L[level] 
				    &&   x_[0]<maxX_refZylinder_of_L[level] 
				    &&    radius < radiusX_refZylinder_of_L[level])
				{

					if(temp_Block->child_array[oct])
                                	temp_Block->child_array[oct]->initial_refined = true;

				  	else if(temp_Block->is_refinement_permitted(oct))
				  	{
						temp_Block->refine_Oct(oct);
						temp_Block->child_array[oct]->initial_refined = true;
		
						num_blocks_createdL[level+1]++;
				  	}

				}
			  }



			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}




	log_file << " done." << endl;
	log_file << " Num blocks in ROOT Level: " << num_root_blocks << endl;
	INT32 j = 0;
	for(j=1; j<= MAX_LEVEL; j++)num_blocks_createdL[MAX_LEVEL+1] += num_blocks_createdL[j];
	for(j=1; j<= MAX_LEVEL; j++)
	log_file << " Num blocks created in Level "<<j<<": " << num_blocks_createdL[j] << endl;
	log_file << " Total child blocks created: " << num_blocks_createdL[MAX_LEVEL+1] << endl;
	log_file << " Total blocks (including. ROOT): " << total_active_Blocks << endl;





	log_file << endl;
// 	log_file << "Block File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*sizeof(CBlock*)+2 << " byte"  << endl;
// 	log_file << "BField File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*3*num_nodes_in_block+2 << " byte"  << endl;


}


//!------------------------------------------------------------
//!- static_refinement_sphere: Generate initial Childs 
//!				- eg. near the obstacle
//!------------------------------------------------------------
void CHybrid::static_refinement_ZylinderY(void)
{




	D_REAL x_[3],radius;
	log_file << " Creating children Zylinder Y-Axis  ...   " << endl;


	INT64 num_blocks_createdL[MAX_LEVEL+2];
	memset(num_blocks_createdL,0,(MAX_LEVEL+2)*sizeof(INT64));


	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];

		if(TL_REFINE_STATIC_ZYLINDER_Y[level] >=0 && TL_REFINE_STATIC_ZYLINDER_Y[level] <= TL)
		while(temp_Block)
		{




			//! ---- REFINE LEVEL [level] -----------------
			//! in case state was restored existing child has to be markes
			//! as initial refined block

			for(INT32 a=0; a<2; a++)
			 for(INT32 b=0; b<2; b++)
			  for(INT32 c=0; c<2; c++)
			  {

				x_[0] = temp_Block->origin[0]+(0.25 +0.5*a)*Blk_Length_of[level][0];
				x_[1] = temp_Block->origin[1]+(0.25 +0.5*b)*Blk_Length_of[level][1];
				x_[2] = temp_Block->origin[2]+(0.25 +0.5*c)*Blk_Length_of[level][2];


				radius = sqrt( x_[0]*x_[0] + x_[2]*x_[2]);

				INT32 oct = 2*2*a +2*b +c;

				if(      x_[1]>=minY_refZylinder_of_L[level] 
				    &&   x_[1]<maxY_refZylinder_of_L[level] 
				    &&    radius < radiusY_refZylinder_of_L[level])
				{

					if(temp_Block->child_array[oct])
                                	temp_Block->child_array[oct]->initial_refined = true;

				  	else if(temp_Block->is_refinement_permitted(oct))
				  	{
						temp_Block->refine_Oct(oct);
						temp_Block->child_array[oct]->initial_refined = true;
		
						num_blocks_createdL[level+1]++;
				  	}

				}
			  }



			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}




	log_file << " done." << endl;
	log_file << " Num blocks in ROOT Level: " << num_root_blocks << endl;
	INT32 j = 0;
	for(j=1; j<= MAX_LEVEL; j++)num_blocks_createdL[MAX_LEVEL+1] += num_blocks_createdL[j];
	for(j=1; j<= MAX_LEVEL; j++)
	log_file << " Num blocks created in Level "<<j<<": " << num_blocks_createdL[j] << endl;
	log_file << " Total child blocks created: " << num_blocks_createdL[MAX_LEVEL+1] << endl;
	log_file << " Total blocks (including. ROOT): " << total_active_Blocks << endl;





	log_file << endl;
// 	log_file << "Block File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*sizeof(CBlock*)+2 << " byte"  << endl;
// 	log_file << "BField File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*3*num_nodes_in_block+2 << " byte"  << endl;


}



//!------------------------------------------------------------
//!- static_refinement_sphere: Generate initial Childs 
//!				- eg. near the obstacle
//!------------------------------------------------------------
void CHybrid::static_refinement_ZylinderZ(void)
{




	D_REAL x_[3],radius;
	log_file << " Creating children Zylinder Z-Axis  ...   " << endl;


	INT64 num_blocks_createdL[MAX_LEVEL+2];
	memset(num_blocks_createdL,0,(MAX_LEVEL+2)*sizeof(INT64));


	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];

		if(TL_REFINE_STATIC_ZYLINDER_Z[level] >=0 && TL_REFINE_STATIC_ZYLINDER_Z[level] <= TL)
		while(temp_Block)
		{




			//! ---- REFINE LEVEL [level] -----------------
			//! in case state was restored existing child has to be markes
			//! as initial refined block

			for(INT32 a=0; a<2; a++)
			 for(INT32 b=0; b<2; b++)
			  for(INT32 c=0; c<2; c++)
			  {

				x_[0] = temp_Block->origin[0]+(0.25 +0.5*a)*Blk_Length_of[level][0];
				x_[1] = temp_Block->origin[1]+(0.25 +0.5*b)*Blk_Length_of[level][1];
				x_[2] = temp_Block->origin[2]+(0.25 +0.5*c)*Blk_Length_of[level][2];


				radius = sqrt( x_[0]*x_[0] + x_[1]*x_[1]);

				INT32 oct = 2*2*a +2*b +c;

				if(      x_[2]>=minZ_refZylinder_of_L[level]
				    &&   x_[2]<maxZ_refZylinder_of_L[level]
				    &&    radius < radiusZ_refZylinder_of_L[level])
				{

					if(temp_Block->child_array[oct])
                                	temp_Block->child_array[oct]->initial_refined = true;

				  	else if(temp_Block->is_refinement_permitted(oct))
				  	{
						temp_Block->refine_Oct(oct);
						temp_Block->child_array[oct]->initial_refined = true;
		
						num_blocks_createdL[level+1]++;
				  	}

				}
			  }



			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}




	log_file << " done." << endl;
	log_file << " Num blocks in ROOT Level: " << num_root_blocks << endl;
	INT32 j = 0;
	for(j=1; j<= MAX_LEVEL; j++)num_blocks_createdL[MAX_LEVEL+1] += num_blocks_createdL[j];
	for(j=1; j<= MAX_LEVEL; j++)
	log_file << " Num blocks created in Level "<<j<<": " << num_blocks_createdL[j] << endl;
	log_file << " Total child blocks created: " << num_blocks_createdL[MAX_LEVEL+1] << endl;
	log_file << " Total blocks (including. ROOT): " << total_active_Blocks << endl;





	log_file << endl;
// 	log_file << "Block File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*sizeof(CBlock*)+2 << " byte"  << endl;
// 	log_file << "BField File Size results in:  " <<
// 	 (num_blocks_createdL[MAX_LEVEL+1]+num_root_blocks)*3*num_nodes_in_block+2 << " byte"  << endl;


}
