


#include "CBlk.h"
#include "defines.h"
#include "parameters.h"
// #include "utils.h"


#include <iostream>
#include <fstream>
#include <math.h>
#include "absolute_Globals.h"


#define id_REQ_INFO_m1	0
#define id_REQ_INFO_p1	1
#define id_REQ_PART_m1	2
#define id_REQ_PART_p1	3

#define pID_TO_m1_NEIGHBOUR	0
#define pID_TO_p1_NEIGHBOUR	1

#define pID_FROM_m1_NEIGHBOUR	pID_TO_p1_NEIGHBOUR
#define pID_FROM_p1_NEIGHBOUR	pID_TO_m1_NEIGHBOUR





//!-----------------------------------------------------------
//! delete_package_memory
//!-----------------------------------------------------------
void CBlock::delete_package_memory(void)
{

	for(INT32 entry=0; entry<8; entry++)
	{

		//! delete memory that has been passed to "send_partPackage_MPI"  before.
		//! (the other memory buffers are released after "insert_add_pPackage")
		//! (so in principal this function only required for mpiruns).
// 		log_file << "particle_package " << entry <<":  "<< particle_package[entry] << endl;

		if(particle_package[entry])
		delete[] particle_package[entry];

		if(MPiC_Info_package[entry])
		delete[] MPiC_Info_package[entry];
	}

	memset(MPiC_Info_package, 0, 8*sizeof(INT32*));
	memset(particle_package , 0, 8*sizeof(particle*));

}

//!---------------------------------------------------------------------------
//! prepare_pPackage
//! Assign particle to i neighbour in equal level
//! in order to send all particle at once, they are
//! stored in a single array. In order to store the
//! information, which cell they belong to, MPiC_Info
//! is sent as well
//! 
//! start_ijk:
//! end_ijk:
//! pID: package ID
//!---------------------------------------------------------------------------
void CBlock::prepare_pPackage(INT32* start_ijk,
			      INT32* end_ijk,
			      INT32  pID)
{


	//! temporary variables
	INT32 index_of_particle_in_package, cell_index_of_block;
	
	
	INT32 num_nodes_in_package = (end_ijk[0] -start_ijk[0])
				   *(end_ijk[1] -start_ijk[1])
				   *(end_ijk[2] -start_ijk[2]);


	//! number of paticles that leave
	//! the Block in i-Direction
	INT32 num_partilce_in_package=0;

	//! estimate how many particle shall be sent:
	for(INT32 species=0; species < num_Charged_Species; species++)
	 for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
	  for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
	   for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
	   num_partilce_in_package += num_MPiC[species][i*BlkNds_Y*BlkNds_Z +j*BlkNds_Z +k];



	//! TODO:
	//! In principal "num_partilce_to_send" could be zero.
	//! (e.g. v=0 -> no reassignement required)
	//! Maybe this should be intercepted and no memory allocated ?!?

	//! alloc package array
	particle_package[pID] = new particle[num_partilce_in_package];

	//! alloc package info array
	MPiC_Info_package[pID] = new INT32[num_nodes_in_package *num_Charged_Species +1];

	//! store "num_partilce_in_package" at pen-ultimate position of MPiC_Info_package
	MPiC_Info_package[pID][num_nodes_in_package *num_Charged_Species] = num_partilce_in_package;


	//! set zero index_of_particle_in_package
	index_of_particle_in_package = 0;

	//! index of cell in layer that is sent
	//! (precisely it is:				   )
	//! (species* num_nodes_in_block +cell_index)
	INT32 cell_index_of_package = 0;

	//! fill package
	for(INT32 species=0; species < num_Charged_Species; species++)
	 for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
	  for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
	   for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
	   {
		//! Index of cell in layer of Block 
		cell_index_of_block = i*BlkNds_Y*BlkNds_Z + j*BlkNds_Z +k;


		//! copy all particle of cell to package
		memcpy(particle_package[pID] +index_of_particle_in_package,
			  pArray[species][cell_index_of_block],
			num_MPiC[species][cell_index_of_block] *sizeof(particle));

		//! Store size of cell's particle Array:
		//! An alternative method would be to add cell indices,
		//! than send, sort into new block and substract indices
		//! again. However, this might introduces aditional 
		//! computaion and maybe even rond off errors
		MPiC_Info_package[pID][cell_index_of_package] = num_MPiC[species][cell_index_of_block];


		//! increase package position
		index_of_particle_in_package += num_MPiC[species][cell_index_of_block];

		//! Finally set MPiC of ghost cell to zero
		//! (this means deleting particle information, even
		//!   though no physical memory is deleted,
		//!   but pArray[species][cell_index_of_block] of Block
		//!   is now assumed to be empty)
		num_MPiC[species][cell_index_of_block]=0;

		//! increase cell index in package-layer
		cell_index_of_package++;

	   }


}



//!---------------------------------------------------------------------------
//! insert_add_pPackage
//! Assign particle to i neighbour in equal level
//! in order to send all particle at once, they are
//! stored in a single array. In order to store the
//! information, which cell they belong to, MPiC_Info
//! is sent as well
//! 
//! i_dest:		      i_dest of layer
//! particle_package:  array for paticles that arrivefrom i-Direction
//! MPiC_Info_package: array to store how many MP are in each cell
//!---------------------------------------------------------------------------
void CBlock::insert_add_pPackage(INT32* start_ijk,
				 INT32* end_ijk,
				 particle* &particle_package,
				 INT32* &MPiC_Info_package)
{

	//! temporary variables
	INT32 index_of_particle_in_package, cell_index_of_block;

	//! set zero index_of_particle_in_package
	index_of_particle_in_package = 0;

	//! index of cell in layer which is received
	//! (precisely it is:				   )
	//! (species* num_nodes_in_block +cell_index)
	INT32 cell_index_of_package=0;


	//! insert package
	for(INT32 species=0; species < num_Charged_Species; species++)
	 for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
	  for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
	   for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
	   {
		//! Index of cell in layer of Block 
		cell_index_of_block = i*BlkNds_Y*BlkNds_Z + j*BlkNds_Z +k;
		
		//! sort particle of package to new cell
		for(INT32 part_index=0; part_index<MPiC_Info_package[cell_index_of_package]; part_index++)
		add_particle_to_pArray(species, cell_index_of_block, particle_package
				       +index_of_particle_in_package +part_index);

		//! increase package position
		index_of_particle_in_package += MPiC_Info_package[cell_index_of_package];
		
		//! increase cell index in package-layer
		cell_index_of_package++;
			
	   }


	delete[] particle_package;
	delete[] MPiC_Info_package;

	particle_package = 0;
	MPiC_Info_package = 0;

}

//!---------------------------------------------------------------------------
//! insert_pPackage_emptyBlock
//! use this method in case block is empty and
//! NO PARTICLE MEMORY IS ALLOCATED YET !!!!!
//!---------------------------------------------------------------------------
void CBlock::insert_pPackage_emptyBlock(INT32* start_ijk,
					INT32* end_ijk,
					particle* &particle_package,
					INT32* &MPiC_Info_package)
{

	//! temporary variables
	INT32 index_of_particle_in_package, cell_index_of_block;

	//! set zero index_of_particle_in_package
	index_of_particle_in_package = 0;

	//! index of cell in layer which is received
	//! (precisely it is:				   )
	//! (species* num_nodes_in_block +cell_index)
	INT32 cell_index_of_package=0;


	//! insert package
	for(INT32 species=0; species < num_Charged_Species; species++)
	 for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
	  for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
	   for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
	   {
		//! Index of cell in layer of Block 
		cell_index_of_block = i*BlkNds_Y*BlkNds_Z + j*BlkNds_Z +k;


		//! memory is not allocated yet, so use "resize_pArray" function
		resize_pArray(species, cell_index_of_block, int(min_pArray_Size_factor *MPiC_Info_package[cell_index_of_package]));

		//! copy particle from package into pArray
		memcpy(pArray[species][cell_index_of_block],
		       particle_package +index_of_particle_in_package,
		       MPiC_Info_package[cell_index_of_package] *sizeof(particle));


		//! protocol MPiC in array of block
		num_MPiC[species][cell_index_of_block] = MPiC_Info_package[cell_index_of_package];

		//! increase package position
		index_of_particle_in_package += MPiC_Info_package[cell_index_of_package];
		
		//! increase cell index in package-layer
		cell_index_of_package++;
			
	   }


	delete[] particle_package;
	delete[] MPiC_Info_package;

	particle_package = 0;
	MPiC_Info_package = 0;

}

//!---------------------------------------------------------------------------
//! insert_child_pPackage
//!---------------------------------------------------------------------------
void CBlock::insert_child_pPackage(INT32* start_ijk,
				   INT32* end_ijk,
				   INT32 child_index,
				   INT32 child_package_ID)
{


	//! temp variables
	INT32 dest_cell, dest_index[3], child_Blk_Index[3];
	PARTICLE_REAL r_offset[3];
	particle* active_particle;

	child_Blk_Index[0] =  child_index/4;
	child_Blk_Index[1] = (child_index -child_Blk_Index[0]*4)/2;
	child_Blk_Index[2] = (child_index -child_Blk_Index[0]*4 -child_Blk_Index[1]*2);



	//! set pointer to respective child memory
	particle* child_particle_package =  child_array[child_index]->particle_package[child_package_ID];
	INT32*     child_package_MPiC_Info = child_array[child_index]->MPiC_Info_package[child_package_ID];


	//! set zero index_of_particle_in_package
	INT32 index_of_particle_in_package = 0;
	
	//! index of cell in layer which is received
	//! (precisely it is:				   )
	//! (species* num_nodes_in_block +cell_index)
	INT32 cell_index_of_package = 0;
	

	//! loop over package cells
	for(INT32 species=0; species < num_Charged_Species; species++)
	 for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
	  for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
	   for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
	   {


		//!-----------------------------------------------------
		//! Calculate in which cell particle shall be inserted:
		//! (compare with sketch !!!)
		//! e.g. Blk_Index=0
		//! i=0 -> dest=0
		//! i=1 -> dest=1
		//! i=2 -> dest=1
		//! i=3 -> dest=2
		//! i=4 -> dest=2
		//!  ...	...
		//! e.g. Blk_Index=1, BlkNds=6
		//! i=0 -> dest=2
		//! i=1 -> dest=3
		//! i=2 -> dest=3
		//! i=3 -> dest=4
		//! i=4 -> dest=4
		//!  ...	...
		dest_index[0] = int(0.5* ((i-1)+2 +child_Blk_Index[0]*(BlkNds_X-2)));
		dest_index[1] = int(0.5* ((j-1)+2 +child_Blk_Index[1]*(BlkNds_Y-2)));
		dest_index[2] = int(0.5* ((k-1)+2 +child_Blk_Index[2]*(BlkNds_Z-2)));
		
		dest_cell =  dest_index[0]*BlkNds_Y*BlkNds_Z
			    +dest_index[1]*BlkNds_Z
			    +dest_index[2];
		//!-----------------------------------------------------
	
		//! particle from e.g. cell i=1 & i=2 are inserted both to dest 1, however
		//! an offset of 0.5 has to be addes to particle from i=2 cell
		//! -> add 0.5 to EVEN numbers
		r_offset[0] = 0.5 *((i+1)%2);
		r_offset[1] = 0.5 *((j+1)%2);
		r_offset[2] = 0.5 *((k+1)%2);

		//! process particle of received layer
		for(INT32 part_index=0; part_index < child_package_MPiC_Info[cell_index_of_package]; part_index++)
		{
				//! get pointer to active particle
			active_particle =  child_particle_package +index_of_particle_in_package +part_index;
			
			//! ajust positions of particle to level-1 cell
			active_particle->rel_r[0] = r_offset[0]+ 0.5 * active_particle->rel_r[0];
			active_particle->rel_r[1] = r_offset[1]+ 0.5 * active_particle->rel_r[1];
			active_particle->rel_r[2] = r_offset[2]+ 0.5 * active_particle->rel_r[2];

#ifdef DEBUG_MODUS
			if(    active_particle->rel_r[0] >=1.
			    || active_particle->rel_r[1] >=1.
			    || active_particle->rel_r[2] >=1.
/*			    ||(//! inside Block ?
   			     ( dest_index[0] >0 && dest_index[0] <BlkNds_X-1)
   			    && ( dest_index[1] >0 && dest_index[1] <BlkNds_Y-1)
   			    && ( dest_index[2] >0 && dest_index[2] <BlkNds_Z-1))*/)
			{



				log_file << "dest_index[0] " << dest_index[0] << endl;
				log_file << "dest_index[1] " << dest_index[1] << endl;
				log_file << "dest_index[2] " << dest_index[2] << endl;

				log_file << "start_ijk[0] " << start_ijk[0] << endl;
				log_file << "start_ijk[1] " << start_ijk[1] << endl;
				log_file << "start_ijk[2] " << start_ijk[2] << endl;


				log_file << "end_ijk[0] " << end_ijk[0] << endl;
				log_file << "end_ijk[1] " << end_ijk[1] << endl;
				log_file << "end_ijk[2] " << end_ijk[2] << endl;


				log_file << "ERROR in insert_child_pPackage" << endl;
// 				exit(1);

			}
			
#endif
			//! finally insert particle in new cell
			add_particle_to_pArray(species, dest_cell, active_particle);
		}

	
		//! increase package position
		index_of_particle_in_package += child_package_MPiC_Info[cell_index_of_package];
		
		//! increase cell index in i-layer
		cell_index_of_package++;

	}

	num_total_particles_in_L[RLevel  ]+= index_of_particle_in_package;
	num_total_particles_in_L[RLevel+1]-= index_of_particle_in_package;

}


//!----------------------------------------------------------------------------
//! insert_parent_pPackage_into_childArray
//! use this function to insert layers of parent or the entire pArray of parent
//! into the childArry
//!----------------------------------------------------------------------------
void CBlock::insert_parent_pPackage(INT32* start_ijk,
				    INT32* end_ijk,
				    INT32 id_package)
{


	//! temp variables
	particle* active_particle;
	INT32 u,v,w, cell_index_of_child;
	INT32 child_cell_index[3];
	
	//! index of cell in layer which is received
	//! (precisely it is:				   )
	//! (species* num_nodes_in_block +cell_index)
	INT32 cell_index_of_package = 0;

	//! set zero index_of_particle_in_package
	INT32 index_of_particle_in_package = 0;

	//! set pointer to respective child memory
	particle* parent_particle_package =  parent->particle_package[id_package];
	INT32*     parent_package_MPiC_Info = parent->MPiC_Info_package[id_package];


	//! Loop across parent package cells
	//! start_ijk has to be >= 1
	//! -> must not include GC
	for(INT32 species=0; species < num_Charged_Species; species++)
	 for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
	  for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
	   for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
	   {

	

		//! calculate corresponding cell index of child
		//! e.g.: BlkNds_X=8
		//! i=1 -> child_cell_index=1
		//! i=2 -> child_cell_index=3
		//! i=3 -> child_cell_index=5
		child_cell_index[0] = 2*i -1;
		child_cell_index[1] = 2*j -1;
		child_cell_index[2] = 2*k -1;



		//! process particle of received layer
		for(INT32 part_index=0; part_index < parent_package_MPiC_Info[cell_index_of_package]; part_index++)
		{
			//! get pointer to active particle
			active_particle =  parent_particle_package +index_of_particle_in_package +part_index;

			//! adjust positions of particle to level+1 cell size
			active_particle->rel_r[0] *= 2.;
			active_particle->rel_r[1] *= 2.;
			active_particle->rel_r[2] *= 2.;

			//! get integer for cell index calculation
		        u = int(active_particle->rel_r[0]);
			v = int(active_particle->rel_r[1]);
			w = int(active_particle->rel_r[2]);

			//! adjust positions of particle to interval [0;1[
			active_particle->rel_r[0] -= u;
			active_particle->rel_r[1] -= v;
			active_particle->rel_r[2] -= w;

#ifdef DEBUG_MODUS
			if(    active_particle->rel_r[0] <0.
			    || active_particle->rel_r[1] <0.
			    || active_particle->rel_r[2] <0.
			    || !(child_cell_index[0] +u)
			    || !(child_cell_index[1] +v)
			    || !(child_cell_index[2] +w))
			{
				log_file << "ERROR in insert_parent_pPackage_into_childArray" << endl;
				exit(1);

			}

		if(  !child_cell_index[0] || !child_cell_index[1] || !child_cell_index[2]
		   || child_cell_index[0] +u > BlkNds_X-2
		   || child_cell_index[1] +v > BlkNds_Y-2
		   || child_cell_index[2] +w > BlkNds_Z-2)
		{

			log_file << "must not" << endl; exit(1);
		}

			
#endif


			//! Index of cell in child block
			cell_index_of_child =  (child_cell_index[0] +u)*BlkNds_Y*BlkNds_Z
					      +(child_cell_index[1] +v)*BlkNds_Z
					      +(child_cell_index[2] +w);
			
			//! finally insert particle in new cell
			add_particle_to_pArray(species, cell_index_of_child, active_particle);
		}

		//! increase package position
		index_of_particle_in_package += parent_package_MPiC_Info[cell_index_of_package];
		
		//! increase cell index in i-layer
		cell_index_of_package++;


	}

	num_total_particles_in_L[RLevel  ]+= index_of_particle_in_package;
	num_total_particles_in_L[RLevel-1]-= index_of_particle_in_package;

}

//!----------------------------------------------------------------------------
//! insert_parent_pPackage_into_childArray
//! use this function to insert layers of parent or the entire pArray of parent
//! into the childArry
//!----------------------------------------------------------------------------
void CBlock::insert_parent_pPackage_into_childArray(INT32* start_ijk,
					       INT32* end_ijk,
					       INT32 id_package)
{


	//! temp variables
	CBlock *active_child;
	particle* active_particle;
	INT32 u,v,w,child, cell_index_of_child;
	INT32 child_cell_index[3], child_index[3];
	
	//! index of cell in layer which is received
	//! (precisely it is:				   )
	//! (species* num_nodes_in_block +cell_index)
	INT32 cell_index_of_package = 0;

	//! set zero index_of_particle_in_package
	INT32 index_of_particle_in_package = 0;

	//! set pointer to respective child memory
	particle* parent_particle_package =  parent->particle_package[id_package];
	INT32*     parent_package_MPiC_Info = parent->MPiC_Info_package[id_package];


	//! Loop across parent package cells
	//! even though it might be better cache optimized to loop over child cell,
	//! this is not possible (or at least very circuitous) since the parent's
	//! particle package is a lineary memory block and accesing this memory
	//! in a different order than it has been stored would involve to
	//! compute the respective parent cell-starting-address by summing
	//! over all num_MPiC.
	//! So keep the order of parent Block and sort particle into respective
	//! child
	//! NOTE:
	//! start_ijk hasto be >= 1
	//! -> must not include GC
	for(INT32 species=0; species < num_Charged_Species; species++)
	 for(INT32 i = start_ijk[0]; i < end_ijk[0]; i++)
	  for(INT32 j = start_ijk[1]; j < end_ijk[1]; j++)
	   for(INT32 k = start_ijk[2]; k < end_ijk[2]; k++)
	   {

		//! e.g.: BlkNds_X=8
		//! i=1 -> child_index=0
		//! i=2 -> child_index=0
		//! i=3 -> child_index=0
		//! i=4 -> child_index=1
		//! i=5 -> child_index=1
		//! i=6 -> child_index=1
		child_index[0] = (i-1)/( (BlkNds_X-2)/2 );
		child_index[1] = (j-1)/( (BlkNds_Y-2)/2 );
		child_index[2] = (k-1)/( (BlkNds_Z-2)/2 );


		//! calculate corresponding cell index of child
		//! e.g.: BlkNds_X=8
		//! i=1 -> child_cell_index=1
		//! i=2 -> child_cell_index=3
		//! i=3 -> child_cell_index=5
		//! i=4 -> child_cell_index=1
		//! i=5 -> child_cell_index=3
		//! i=6 -> child_cell_index=5
		child_cell_index[0] = 2*i -1 - child_index[0]* (BlkNds_X-2);
		child_cell_index[1] = 2*j -1 - child_index[1]* (BlkNds_Y-2);
		child_cell_index[2] = 2*k -1 - child_index[2]* (BlkNds_Z-2);


		//! calc which child of array
		child =  child_index[0]*2*2
			+child_index[1]*2
			+child_index[2];

		//! set pointer to arrays
		active_child = this +child;


		//! process particle of received layer
		for(INT32 part_index=0; part_index < parent_package_MPiC_Info[cell_index_of_package]; part_index++)
		{
			//! get pointer to active particle
			active_particle =  parent_particle_package +index_of_particle_in_package +part_index;

			//! adjust positions of particle to level+1 cell size
			active_particle->rel_r[0] *= 2.;
			active_particle->rel_r[1] *= 2.;
			active_particle->rel_r[2] *= 2.;

			//! get integer for cell index calculation
		        u = int(active_particle->rel_r[0]);
			v = int(active_particle->rel_r[1]);
			w = int(active_particle->rel_r[2]);

			//! adjust positions of particle to interval [0;1[
			active_particle->rel_r[0] -= u;
			active_particle->rel_r[1] -= v;
			active_particle->rel_r[2] -= w;

#ifdef DEBUG_MODUS
			if(    active_particle->rel_r[0] <0.
			    || active_particle->rel_r[1] <0.
			    || active_particle->rel_r[2] <0.
			    || !(child_cell_index[0] +u)
			    || !(child_cell_index[1] +v)
			    || !(child_cell_index[2] +w))
			{
				log_file << "ERROR in insert_parent_pPackage_into_childArray" << endl;
				exit(1);

			}
			
#endif


			//! Index of cell in child block
			cell_index_of_child =  (child_cell_index[0] +u)*BlkNds_Y*BlkNds_Z
					      +(child_cell_index[1] +v)*BlkNds_Z
					      +(child_cell_index[2] +w);
			
			//! finally insert particle in new cell
			active_child->add_particle_to_pArray(species, cell_index_of_child, active_particle);
		}

		//! increase package position
		index_of_particle_in_package += parent_package_MPiC_Info[cell_index_of_package];
		
		//! increase cell index in i-layer
		cell_index_of_package++;


	}

	num_total_particles_in_L[RLevel  ]+= index_of_particle_in_package;
	num_total_particles_in_L[RLevel-1]-= index_of_particle_in_package;

}


//!-----------------------------------------------------------
//! send_pPackage_toParent_MPI
//! six different directions are possible:
//! mX, pX, mY, pY, mZ, pZ
//!-----------------------------------------------------------
void CBlock::send_pPackage_toParent_MPI(INT32 id_direc, INT32 entries_in_MPiC_Info_package)
{


	INT32 INFO_request = id_direc;
	INT32 PART_request = id_direc +6;

	INT32 INFO_tag = mpi_tag +(id_direc +0)*total_num_mpi_tags;
	INT32 PART_tag = mpi_tag +(id_direc +6)*total_num_mpi_tags;


	send_MPiC_Package_MPI(INFO_request,
			      INFO_tag,
			      parent->responsible_mpi_process,
			      entries_in_MPiC_Info_package,
			      MPiC_Info_package[id_direc]);

	send_partPackage_MPI(PART_request,
			     PART_tag,
			     parent->responsible_mpi_process,
			     entries_in_MPiC_Info_package,
			     particle_package[id_direc],
			     MPiC_Info_package[id_direc]);
}


//!-----------------------------------------------------------
//! recv_pPackage_fromChild_MPI
//! six different directions are possible:
//! mX, pX, mY, pY, mZ, pZ
//!-----------------------------------------------------------
void CBlock::recv_pPackage_fromChild_MPI(INT32 id_child,
					 INT32 id_direc,
					 INT32 entries_in_MPiC_Info_package)
{



	INT32 INFO_tag = child_array[id_child]->mpi_tag +(id_direc +0)*total_num_mpi_tags;
	INT32 PART_tag = child_array[id_child]->mpi_tag +(id_direc +6)*total_num_mpi_tags;


	receive_infoPackage_MPI(INFO_tag,
				child_array[id_child]->responsible_mpi_process,
				entries_in_MPiC_Info_package,
				child_array[id_child]->MPiC_Info_package[id_direc]);

	receive_partPackage_MPI(PART_tag,
				child_array[id_child]->responsible_mpi_process,
				entries_in_MPiC_Info_package,
				child_array[id_child]->particle_package[id_direc],
				child_array[id_child]->MPiC_Info_package[id_direc]);
}


//!-----------------------------------------------------------
//! send_particle_toParent_LB
//!-----------------------------------------------------------
void CBlock::send_particle_toParent_LB(void)
{


	//! -------------------------------------------------------------------
	//! -- ALWAYS DELETE PACKAGE MEMORY BEFORE SEND CALL ------------------
	//! - HERE: delete memory from parent->children sent of last TL sent --
	//! -------------------------------------------------------------------
	delete_package_memory();



	//!---------------------------------------------------------------
	//! ---------- i - Direction -------------------------------------
	//!---------------------------------------------------------------

	//! Block has no neighbour and is not a im1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Im1_] && !is_box_boundary[_Im1_])
	{

		//! define X=0 layer
		INT32 start_ijk[3] = {0,	      0,       0};
		INT32 end_ijk[3] =   {1,BlkNds_Y,BlkNds_Z};

		//! in equal level
		prepare_pPackage(start_ijk, end_ijk, _Im1_);


		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		send_pPackage_toParent_MPI(_Im1_, BlkNds_Y *BlkNds_Z *num_Charged_Species +1);

	}

	//! Block has no neighbour and is not a ip1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Ip1_] && !is_box_boundary[_Ip1_])
	{

		//! define X=BlkNds_X-1 layer
		INT32 start_ijk[3] = {BlkNds_X-1,       0,        0};
		INT32 end_ijk[3]   = {BlkNds_X,  BlkNds_Y, BlkNds_Z};

		//! in equal level
		prepare_pPackage(start_ijk, end_ijk, _Ip1_);


		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		send_pPackage_toParent_MPI(_Ip1_, BlkNds_Y *BlkNds_Z *num_Charged_Species +1);

	}

	//!---------------------------------------------------------------
	//! ---------- j - Direction -------------------------------------
	//!---------------------------------------------------------------
	//! Block has no neighbour and is not a jm1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Jm1_] && !is_box_boundary[_Jm1_])
	{

		//! define Y=0 layer
		INT32 start_ijk[3] = {0,	      0,       0};
		INT32 end_ijk[3] =   {BlkNds_X,1,BlkNds_Z};

		//! in equal level
		prepare_pPackage(start_ijk, end_ijk, _Jm1_);

		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		send_pPackage_toParent_MPI(_Jm1_, BlkNds_X *BlkNds_Z *num_Charged_Species +1);

	}

	//! Block has no neighbour and is not a jp1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Jp1_] && !is_box_boundary[_Jp1_])
	{

		//! define Y=BlkNds_Y-1 layer
		INT32 start_ijk[3] = {0       , BlkNds_Y-1,       0};
		INT32 end_ijk[3] =   {BlkNds_X, BlkNds_Y  ,BlkNds_Z};

		//! in equal level
		prepare_pPackage(start_ijk, end_ijk, _Jp1_);

		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		send_pPackage_toParent_MPI(_Jp1_, BlkNds_X *BlkNds_Z *num_Charged_Species +1);


	}

	//!---------------------------------------------------------------
	//! ---------- k - Direction -------------------------------------
	//!---------------------------------------------------------------
	//! Block has no neighbour and is not a km1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Km1_] && !is_box_boundary[_Km1_])
	{

		//! define X=0 layer
		INT32 start_ijk[3] = {	    0,	     0,0};
		INT32 end_ijk[3] =   {BlkNds_X,BlkNds_Y,1};


		//! in equal level
		prepare_pPackage(start_ijk, end_ijk, _Km1_);

		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		send_pPackage_toParent_MPI(_Km1_, BlkNds_X *BlkNds_Y *num_Charged_Species +1);

	}

	//! Block has no neighbour and is not a jp1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Kp1_] && !is_box_boundary[_Kp1_])
	{

		//! define Y=BlkNds_Y-1 layer
		INT32 start_ijk[3] = {	    0,	     0,BlkNds_Z-1};
		INT32 end_ijk[3] =   {BlkNds_X,BlkNds_Y,  BlkNds_Z};


		//! in equal level
		prepare_pPackage(start_ijk, end_ijk, _Kp1_);

		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		send_pPackage_toParent_MPI(_Kp1_, BlkNds_X *BlkNds_Y *num_Charged_Species +1);


	}

}


//!-----------------------------------------------------------
//! wait_send_particle_toParent_LB
//!-----------------------------------------------------------
void CBlock::wait_send_particle_toParent_LB(void)
{

	//!---------------------------------------------------------------
	//! ---------- i - Direction -------------------------------------
	//!---------------------------------------------------------------

	//! Block has no neighbour and is not a im1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Im1_] && !is_box_boundary[_Im1_])
	{

		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		{
		  INT32 INFO_request = _Im1_;
		  INT32 PART_request = _Im1_ +6;
		  req_is_MPI_send_completed[INFO_request].Wait();
		  req_is_MPI_send_completed[PART_request].Wait();
		}
	}

	//! Block has no neighbour and is not a ip1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Ip1_] && !is_box_boundary[_Ip1_])
	{
		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		{
		  INT32 INFO_request = _Ip1_;
		  INT32 PART_request = _Ip1_ +6;
		  req_is_MPI_send_completed[INFO_request].Wait();
		  req_is_MPI_send_completed[PART_request].Wait();
		}
	}

	//!---------------------------------------------------------------
	//! ---------- j - Direction -------------------------------------
	//!---------------------------------------------------------------
	//! Block has no neighbour and is not a jm1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Jm1_] && !is_box_boundary[_Jm1_])
	{

		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		{
		  INT32 INFO_request = _Jm1_;
		  INT32 PART_request = _Jm1_ +6;
		  req_is_MPI_send_completed[INFO_request].Wait();
		  req_is_MPI_send_completed[PART_request].Wait();
		}

	}

	//! Block has no neighbour and is not a jp1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Jp1_] && !is_box_boundary[_Jp1_])
	{
		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		{
		  INT32 INFO_request = _Jp1_;
		  INT32 PART_request = _Jp1_ +6;
		  req_is_MPI_send_completed[INFO_request].Wait();
		  req_is_MPI_send_completed[PART_request].Wait();
		}

	}

	//!---------------------------------------------------------------
	//! ---------- k - Direction -------------------------------------
	//!---------------------------------------------------------------
	//! Block has no neighbour and is not a km1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Km1_] && !is_box_boundary[_Km1_])
	{
		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		{
		  INT32 INFO_request = _Km1_;
		  INT32 PART_request = _Km1_ +6;
		  req_is_MPI_send_completed[INFO_request].Wait();
		  req_is_MPI_send_completed[PART_request].Wait();
		}

	}

	//! Block has no neighbour and is not a jp1-Box Boundary, the only
	//! remaining possibility is that his parent has a neighbour
	if(!Neighbour[_Kp1_] && !is_box_boundary[_Kp1_])
	{
		//! check whether sent to parent via MPI is required
		if(parent->responsible_mpi_process != mpi_myRank)
		{
		  INT32 INFO_request = _Kp1_;
		  INT32 PART_request = _Kp1_ +6;
		  req_is_MPI_send_completed[INFO_request].Wait();
		  req_is_MPI_send_completed[PART_request].Wait();
		}
	}

}


//!-----------------------------------------------------------
//! recv_particle_fromChildren_LB
//!-----------------------------------------------------------
void CBlock::recv_particle_fromChildren_LB(void)
{


	//! 1)
	//! loop over all children
	//! check whether child exists

	//! 2) For each of 6 faces:
	//! in case child array has no neighbour and is no is_box_boundary,
	//! child_array particle have to travel to parent
	for(INT32 child=0; child<8; child++)
	if(child_array[child])
	{



		//! ---------------- Im1 Direction ---------------------------------
		if(!child_array[child]->Neighbour[_Im1_] && !child_array[child]->is_box_boundary[_Im1_])
		{

			//! adjust indices to child's i-min GC-Layer
			INT32 start_ijk[3] = {0,	      0,       0};
			INT32 end_ijk[3]   = {1,BlkNds_Y,BlkNds_Z};


			//! check whether recv from child via MPI is required
			if(child_array[child]->responsible_mpi_process != mpi_myRank)
			recv_pPackage_fromChild_MPI(child, _Im1_, BlkNds_Y *BlkNds_Z *num_Charged_Species +1);


			//! insert child's layer into Block
			insert_child_pPackage(start_ijk, end_ijk, child,_Im1_);

		}


		//! ---------------- Ip1 Direction ---------------------------------
		if(!child_array[child]->Neighbour[_Ip1_] && !child_array[child]->is_box_boundary[_Ip1_])
		{

			//! adjust indices to i-max layer
			INT32 start_ijk[3] = {BlkNds_X-1,       0,       0};
			INT32 end_ijk[3] =   {BlkNds_X  ,BlkNds_Y,BlkNds_Z};


			//! check whether recv from child via MPI is required
			if(child_array[child]->responsible_mpi_process != mpi_myRank)
			recv_pPackage_fromChild_MPI(child, _Ip1_, BlkNds_Y *BlkNds_Z *num_Charged_Species +1);


			//! insert child's layer into Block
			insert_child_pPackage(start_ijk,end_ijk,child,_Ip1_);
		}


		//! ---------------- Jm1 Direction ---------------------------------
		if(!child_array[child]->Neighbour[_Jm1_] && !child_array[child]->is_box_boundary[_Jm1_])
		{

			//! adjust indices to j-min layer
			INT32 start_ijk[3] = {       0,0,       0};
			INT32 end_ijk[3] =   {BlkNds_X,1,BlkNds_Z};


			//! check whether recv from child via MPI is required
			if(child_array[child]->responsible_mpi_process != mpi_myRank)
			recv_pPackage_fromChild_MPI(child, _Jm1_, BlkNds_X *BlkNds_Z *num_Charged_Species +1);


			//! insert child's layer into Block
			insert_child_pPackage(start_ijk,end_ijk,child,_Jm1_);
		}

		//! ---------------- Jp1 Direction ---------------------------------
		if(!child_array[child]->Neighbour[_Jp1_] && !child_array[child]->is_box_boundary[_Jp1_])
		{

			//! adjust indices to j-max layer
			INT32 start_ijk[3] = {       0,BlkNds_Y-1,       0};
			INT32 end_ijk[3] =   {BlkNds_X,BlkNds_Y,  BlkNds_Z};


			//! check whether recv from child via MPI is required
			if(child_array[child]->responsible_mpi_process != mpi_myRank)
			recv_pPackage_fromChild_MPI(child, _Jp1_, BlkNds_X *BlkNds_Z *num_Charged_Species +1);



			//! insert child's layer into Block
			insert_child_pPackage(start_ijk,end_ijk,child,_Jp1_);
		}

		//! ---------------- Km1 Direction ---------------------------------
		if(!child_array[child]->Neighbour[_Km1_] && !child_array[child]->is_box_boundary[_Km1_])
		{

			//! adjust indices to k-min layer
			INT32 start_ijk[3] = {       0,       0,0};
			INT32 end_ijk[3] =   {BlkNds_X,BlkNds_Y,1};


			//! check whether recv from child via MPI is required
			if(child_array[child]->responsible_mpi_process != mpi_myRank)
			recv_pPackage_fromChild_MPI(child, _Km1_, BlkNds_X *BlkNds_Y *num_Charged_Species +1);


			//! insert child's layer into Block
			insert_child_pPackage(start_ijk,end_ijk,child,_Km1_);
		}

		//! ---------------- Kp1 Direction ---------------------------------
		if(!child_array[child]->Neighbour[_Kp1_] && !child_array[child]->is_box_boundary[_Kp1_])
		{
	
			//! adjust indices to k-max layer
			INT32 start_ijk[3] = {       0,       0,BlkNds_Z-1};
			INT32 end_ijk[3] =   {BlkNds_X,BlkNds_Y,BlkNds_Z  };


			//! check whether recv from child via MPI is required
			if(child_array[child]->responsible_mpi_process != mpi_myRank)
			recv_pPackage_fromChild_MPI(child, _Kp1_, BlkNds_X *BlkNds_Y *num_Charged_Species +1);

	
			//! insert child's layer into Block
			insert_child_pPackage(start_ijk,end_ijk,child,_Kp1_);
		}


	}
	

}




//!-----------------------------------------------------------
//! send_particle_iDirection
//!-----------------------------------------------------------
void CBlock::prepare_pPackage_toNeighbour(INT32* m1_start_ijk,
					  INT32* m1_end_ijk,
					  INT32* p1_start_ijk,
					  INT32* p1_end_ijk,
					  INT32 _m1_,
					  INT32 _p1_)
{

	//! -----------------------------------------------------
	//! -- ALWAYS DELETE PACKAGE MEMORY BEFORE SEND CALL ----
	//! - HERE: delete memory from children->parent sent ----
	//! -----------------------------------------------------
	delete_package_memory();


	//!-----------------------------------------------------------
	//! send to Neighbour[_m1_]
	//!-----------------------------------------------------------
	//! copy particle/MPiC_Info
	if(Neighbour[_m1_])
	prepare_pPackage(m1_start_ijk, m1_end_ijk, pID_TO_m1_NEIGHBOUR);


	//!-----------------------------------------------------------
	//! send to Neighbour[_p1_]
	//!-----------------------------------------------------------
	//! copy particle/MPiC_Info
	if(Neighbour[_p1_])
	prepare_pPackage(p1_start_ijk, p1_end_ijk, pID_TO_p1_NEIGHBOUR);





}

//!-----------------------------------------------------------
//! send_MPiC_Info_toNeighbour
//!-----------------------------------------------------------
void CBlock::send_MPiC_Info_toNeighbour(INT32 entries_in_MPiC_Info_package,
					INT32 _m1_,
					INT32 _p1_)
{



	INT32 tag_info_to_m1 = mpi_tag +0*total_num_mpi_tags;
	INT32 tag_info_to_p1 = mpi_tag +1*total_num_mpi_tags;

	//!-----------------------------------------------------------
	//! send to Neighbour[_m1_]
	//!-----------------------------------------------------------
	if(Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	send_MPiC_Package_MPI(id_REQ_INFO_m1,
			      tag_info_to_m1,
		 	      Neighbour[_m1_]->responsible_mpi_process,
			      entries_in_MPiC_Info_package,
		 	      MPiC_Info_package[pID_TO_m1_NEIGHBOUR]);


	//!-----------------------------------------------------------
	//! send to Neighbour[_p1_]
	//!-----------------------------------------------------------
	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	send_MPiC_Package_MPI(id_REQ_INFO_p1,
		 	      tag_info_to_p1,
		 	      Neighbour[_p1_]->responsible_mpi_process,
			      entries_in_MPiC_Info_package,
		 	      MPiC_Info_package[pID_TO_p1_NEIGHBOUR]);


}


//!-----------------------------------------------------------
//! recv_MPiC_Info_fromNeighbour
//!-----------------------------------------------------------
void CBlock::recv_MPiC_Info_fromNeighbour(INT32 entries_in_MPiC_Info_package,
					  INT32 _m1_,
					  INT32 _p1_)
{





	//!-----------------------------------------------------------
	//! receive from Neighbour[_m1_]
	//!-----------------------------------------------------------
	if(Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	{

		//! initialize tag INSIDE brackets, else Neighbour might not
		//! exist !!! -> seghmentatiton fault
		INT32 tag_info_from_m1 = Neighbour[_m1_]->mpi_tag +1*total_num_mpi_tags;

		receive_infoPackage_MPI(tag_info_from_m1,
					Neighbour[_m1_]->responsible_mpi_process,
					entries_in_MPiC_Info_package,
					Neighbour[_m1_]->MPiC_Info_package[pID_FROM_m1_NEIGHBOUR]);
	}

	//!-----------------------------------------------------------
	//! receive from Neighbour[_p1_]
	//!-----------------------------------------------------------
	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	{

		//! initialize tag INSIDE brackets, else Neighbour might not
		//! exist !!! -> seghmentatiton fault
		INT32 tag_info_from_p1 = Neighbour[_p1_]->mpi_tag +0*total_num_mpi_tags;

		receive_infoPackage_MPI(tag_info_from_p1,
					Neighbour[_p1_]->responsible_mpi_process,
					entries_in_MPiC_Info_package,
					Neighbour[_p1_]->MPiC_Info_package[pID_FROM_p1_NEIGHBOUR]);
	}








}


//!-----------------------------------------------------------
//! send_particle_iDirection
//!-----------------------------------------------------------
void CBlock::send_particle_toNeighbour(INT32 entries_in_MPiC_Info_package,
			   INT32 _m1_,
			   INT32 _p1_)
{



	INT32 tag_part_to_m1 = mpi_tag +2*total_num_mpi_tags;
	INT32 tag_part_to_p1 = mpi_tag +3*total_num_mpi_tags;

	//!-----------------------------------------------------------
	//! send to Neighbour[_m1_]
	//!-----------------------------------------------------------
	if(Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	send_partPackage_MPI(id_REQ_PART_m1,
			     tag_part_to_m1,
			     Neighbour[_m1_]->responsible_mpi_process,
			     entries_in_MPiC_Info_package,
			     particle_package[pID_TO_m1_NEIGHBOUR],
			     MPiC_Info_package[pID_TO_m1_NEIGHBOUR]);


	//!-----------------------------------------------------------
	//! send to Neighbour[_p1_]
	//!-----------------------------------------------------------
	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	send_partPackage_MPI(id_REQ_PART_p1,
		 	     tag_part_to_p1,
		 	     Neighbour[_p1_]->responsible_mpi_process,
			     entries_in_MPiC_Info_package,
		 	     particle_package[pID_TO_p1_NEIGHBOUR],
		 	     MPiC_Info_package[pID_TO_p1_NEIGHBOUR]);

}





//!-----------------------------------------------------------
//! recv_particle_fromNeighbour
//!-----------------------------------------------------------
void CBlock::recv_particle_fromNeighbour(INT32 entries_in_MPiC_Info_package,
				   INT32 _m1_,
				   INT32 _p1_)
{




	//!-----------------------------------------------------------
	//! receive from Neighbour[_m1_]
	//!-----------------------------------------------------------
	if(Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	{

		//! initialize tag INSIDE brackets, else Neighbour might not
		//! exist !!! -> seghmentatiton fault
		INT32 tag_part_from_m1 = Neighbour[_m1_]->mpi_tag +3*total_num_mpi_tags;

		receive_partPackage_MPI(tag_part_from_m1,
					Neighbour[_m1_]->responsible_mpi_process,
					entries_in_MPiC_Info_package,
					Neighbour[_m1_]->particle_package[pID_FROM_m1_NEIGHBOUR],
					Neighbour[_m1_]->MPiC_Info_package[pID_FROM_m1_NEIGHBOUR]);
	}

	//!-----------------------------------------------------------
	//! receive from Neighbour[_p1_]
	//!-----------------------------------------------------------
	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	{

		//! initialize tag INSIDE brackets, else Neighbour might not
		//! exist !!! -> seghmentatiton fault
		INT32 tag_part_from_p1 = Neighbour[_p1_]->mpi_tag +2*total_num_mpi_tags;

		receive_partPackage_MPI(tag_part_from_p1,
					Neighbour[_p1_]->responsible_mpi_process,
					entries_in_MPiC_Info_package,
					Neighbour[_p1_]->particle_package[pID_FROM_p1_NEIGHBOUR],
					Neighbour[_p1_]->MPiC_Info_package[pID_FROM_p1_NEIGHBOUR]);
	}







	

}

//!---------------------------------------------------------------//
//! wait_for_completedSendParticle_MPI:						  //
//!---------------------------------------------------------------//
void CBlock::wait_for_completedSendParticle_MPI(INT32 _m1_, INT32 _p1_)
{

	//! force MPI to deallocate memory of Last Send
	//! Code will run without this command but result
	//! in a memory leak

	//! Minus will not be used for gather, intercept for this by
	//! specifyin negative value for _m1_

	if(_m1_>=0 && Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	{
	  req_is_MPI_send_completed[id_REQ_INFO_m1].Wait();
	  req_is_MPI_send_completed[id_REQ_PART_m1].Wait();
	}

	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	{
	  req_is_MPI_send_completed[id_REQ_INFO_p1].Wait();
	  req_is_MPI_send_completed[id_REQ_PART_p1].Wait();
	}

// 	MPI::Request::Waitall(NUM_REQUESTS, req_is_MPI_send_completed);

}


//!-----------------------------------------------------------
//! receive_particle_iDirection
//!-----------------------------------------------------------
void CBlock::insert_pPackage_fromNeighbour(INT32* m1_start_ijk,
					   INT32* m1_end_ijk,
					   INT32* p1_start_ijk,
					   INT32* p1_end_ijk,
					   INT32 _m1_,
					   INT32 _p1_)
{



	//!-----------------------------------------------------------
	//! receive from Neighbour[_m1_]
	//!-----------------------------------------------------------
	if(Neighbour[_m1_])
	insert_add_pPackage(m1_start_ijk,
			    m1_end_ijk,
			    Neighbour[_m1_]->particle_package[pID_FROM_m1_NEIGHBOUR],
			    Neighbour[_m1_]->MPiC_Info_package[pID_FROM_m1_NEIGHBOUR]);

	//!-----------------------------------------------------------
	//! receive from Neighbour[_p1_]
	//!-----------------------------------------------------------
	if(Neighbour[_p1_])
	insert_add_pPackage(p1_start_ijk,
			    p1_end_ijk,
			    Neighbour[_p1_]->particle_package[pID_FROM_p1_NEIGHBOUR],
			    Neighbour[_p1_]->MPiC_Info_package[pID_FROM_p1_NEIGHBOUR]);
}

/*
//!-----------------------------------------------------------
//! inserting_received_MPiC_Part
//!-----------------------------------------------------------
void CBlock::inserting_received_MPiC_Part(INT32* m1_start_ijk,
					  INT32* m1_end_ijk,
					  INT32* p1_start_ijk,
					  INT32* p1_end_ijk,
					  INT32 _m1_,
					  INT32 _p1_)
{

	//!-----------------------------------------------------------
	//! receive from Neighbour[_m1_]
	//!-----------------------------------------------------------
	if(Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	insert_add_pPackage(m1_start_ijk,
			    m1_end_ijk,
			    Neighbour[_m1_]->particle_package[pID_FROM_m1_NEIGHBOUR],
			    Neighbour[_m1_]->MPiC_Info_package[pID_FROM_m1_NEIGHBOUR]);
	
	//!-----------------------------------------------------------
	//! receive from Neighbour[_p1_]
	//!-----------------------------------------------------------
	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	insert_add_pPackage(p1_start_ijk,
			    p1_end_ijk,
			    Neighbour[_p1_]->particle_package[pID_FROM_p1_NEIGHBOUR],
			    Neighbour[_p1_]->MPiC_Info_package[pID_FROM_p1_NEIGHBOUR]);


}

*/






//!-----------------------------------------------------------
//! send_particle_toChild
//!-----------------------------------------------------------
void CBlock::send_particle_toChild(INT32 child)
{


	prepare_pPackage_toChild(child);

	//! for clarification:
	//! phys volume of oct = (BlkNds/2-1)³
	if(child_array[child]->responsible_mpi_process != mpi_myRank)
	send_pPackage_toChild_MPI(child, (BlkNds_X/2-1)*(BlkNds_Y/2-1)*(BlkNds_Z/2-1)*num_Charged_Species +1);


}

//!-----------------------------------------------------------
//! recv_particle_fromParent
//!-----------------------------------------------------------
void CBlock::recv_particle_fromParent(void)
{

	if(parent->responsible_mpi_process != mpi_myRank)
	recv_pPackage_fromParent_MPI((BlkNds_X/2-1)*(BlkNds_Y/2-1)*(BlkNds_Z/2-1)*num_Charged_Species +1);

	insert_pPackage_fromParent();

}

//!-----------------------------------------------------------
//! send_particle_toChildren_LB
//!-----------------------------------------------------------
void CBlock::send_particle_toChildren_LB(void)
{

	for(INT32 child=0; child<8; child++)
	if(do_send_particle_to_childArray[child])
	 send_particle_toChild(child);
	

}

//!-----------------------------------------------------------
//! wait_send_particle_toChildren_LB
//!-----------------------------------------------------------
void CBlock::wait_send_particle_toChildren_LB(void)
{
	for(INT32 child=0; child<8; child++)
	{
	  if(do_send_particle_to_childArray[child])
	  {
	    if(child_array[child]->responsible_mpi_process != mpi_myRank)
	      {
		INT32 INFO_request = child;
		INT32 PART_request = child +8;
		req_is_MPI_send_completed[INFO_request].Wait();
		req_is_MPI_send_completed[PART_request].Wait();
	      }
	  }
	}

}

//!-----------------------------------------------------------
//! send_pPackage_toParent_MPI
//! six different directions are possible:
//! mX, pX, mY, pY, mZ, pZ
//!-----------------------------------------------------------
void CBlock::send_pPackage_toChild_MPI(INT32 id_child, INT32 entries_in_MPiC_Info_package)
{


	INT32 INFO_request = id_child;
	INT32 PART_request = id_child +8;

	INT32 INFO_tag = mpi_tag +(id_child +0)*total_num_mpi_tags;
	INT32 PART_tag = mpi_tag +(id_child +8)*total_num_mpi_tags;


	send_MPiC_Package_MPI(INFO_request,
			      INFO_tag,
			      child_array[id_child]->responsible_mpi_process,
			      entries_in_MPiC_Info_package,
			      MPiC_Info_package[id_child]);

	send_partPackage_MPI(PART_request,
			     PART_tag,
			     child_array[id_child]->responsible_mpi_process,
			     entries_in_MPiC_Info_package,
			     particle_package[id_child],
			     MPiC_Info_package[id_child]);
}


//!-----------------------------------------------------------
//! recv_pPackage_fromChild_MPI
//! six different directions are possible:
//! mX, pX, mY, pY, mZ, pZ
//!-----------------------------------------------------------
void CBlock::recv_pPackage_fromParent_MPI(INT32 entries_in_MPiC_Info_package)
{


	INT32 id_child =  Blk_Index[0]*2*2
			+Blk_Index[1]*2
			+Blk_Index[2];

	INT32 INFO_tag = parent->mpi_tag +(id_child +0)*total_num_mpi_tags;
	INT32 PART_tag = parent->mpi_tag +(id_child +8)*total_num_mpi_tags;


	receive_infoPackage_MPI(INFO_tag,
				parent->responsible_mpi_process,
				entries_in_MPiC_Info_package,
				parent->MPiC_Info_package[id_child]);

	receive_partPackage_MPI(PART_tag,
				parent->responsible_mpi_process,
				entries_in_MPiC_Info_package,
				parent->particle_package[id_child],
				parent->MPiC_Info_package[id_child]);
}






//!-----------------------------------------------------------
//! recv_particle_fromParent_LB
//!-----------------------------------------------------------
void CBlock::recv_particle_fromParent_LB(void)
{


	if(do_receive_particle_from_parent)
	recv_particle_fromParent();


}


//!-----------------------------------------------------------
//! send_particle_to_child
//!-----------------------------------------------------------
void CBlock::prepare_pPackage_toChild(INT32 id_child)
{

	//! -----------------------------------------------------
	//! -- ALWAYS DELETE PACKAGE MEMORY BEFORE SEND CALL ----
	//! --------- BUT ONLY FOR RESPECTIVE CHILD -------------
	//! -----------------------------------------------------
	if(particle_package[id_child])
	delete particle_package[id_child];
	particle_package[id_child] = 0;

	if(MPiC_Info_package[id_child])
	delete MPiC_Info_package[id_child];
	MPiC_Info_package[id_child] = 0;


	INT32 a = id_child/4;
	INT32 b = (id_child -a*4)/2;
	INT32 c = (id_child -a*4 -b*2);


	INT32 i_offset = BlkNds_X/2-1;
	INT32 j_offset = BlkNds_Y/2-1;
	INT32 k_offset = BlkNds_Z/2-1;


	//! define entire volume of oct (NO GC, BUT INSIDE PHYSICAL SPACE)
	INT32 start_ijk[3] = {(a+0)*i_offset +1, (b+0)*j_offset+1, (c+0)*k_offset+1};
	INT32 end_ijk[3] =   {(a+1)*i_offset +1, (b+1)*j_offset+1, (c+1)*k_offset+1};

	//! prepare package for respective child
	prepare_pPackage(start_ijk, end_ijk, id_child);

}



//!-----------------------------------------------------------
//! receive_particle_from_parent
//!-----------------------------------------------------------
void CBlock::insert_pPackage_fromParent(void)
{




	//! define corresponding volume in parent idices always using
	//! the node (1,1,1) as origin
	//! (prepare package of parent ensures that correct volume is sent)
	INT32 start_ijk[3] = {		     1,		        1,		   1};
	INT32 end_ijk[3] =   {(BlkNds_X-2)/2 +1, (BlkNds_Y-2)/2 +1, (BlkNds_Z-2)/2 +1};

	INT32 pID =  Blk_Index[0]*2*2
		   +Blk_Index[1]*2
		   +Blk_Index[2];

	insert_parent_pPackage(start_ijk, end_ijk, pID);



}


//!-----------------------------------------------------------
//! send_particle_to_parent
//!-----------------------------------------------------------
void CBlock::prepare_pPackage_toParent(void)
{

	//! -------------------------------------------------------------------
	//! -- ALWAYS DELETE PACKAGE MEMORY BEFORE SEND CALL ------------------
	//! -------------------------------------------------------------------
	delete_package_memory();


	//! define entire volume (NO GC, BUT INSIDE PHYSICAL SPACE)
	INT32 start_ijk[3] = {         1,         1,         1};
	INT32 end_ijk[3] =   {BlkNds_X-1,BlkNds_Y-1,BlkNds_Z-1};

	prepare_pPackage(start_ijk, end_ijk, 0);

	if(parent->responsible_mpi_process!=mpi_myRank)
	send_pPackage_toParent_MPI(0, (BlkNds_X-2)*(BlkNds_Y-2)*(BlkNds_Z-2)*num_Charged_Species +1);

}

//!-----------------------------------------------------------
//! receive_particle_from_child
//!-----------------------------------------------------------
void CBlock::insert_pPackage_fromChild(INT32 id_child)
{


	//! define entire volume of child
	//! (NO GC, BUT INSIDE PHYSICAL SPACE)
	INT32 start_ijk[3] = {         1,         1,         1};
	INT32 end_ijk[3] =   {BlkNds_X-1,BlkNds_Y-1,BlkNds_Z-1};


	if(child_array[id_child]->responsible_mpi_process != mpi_myRank)
	recv_pPackage_fromChild_MPI(id_child, 0, (BlkNds_X-2)*(BlkNds_Y-2)*(BlkNds_Z-2)*num_Charged_Species +1);

	//! insert child's package into Block
	//! 0 is used as package id
	insert_child_pPackage(start_ijk, end_ijk, id_child, 0);


}


//!-----------------------------------------------------------
//! root_of_array is the first of 8 blocks.
//!-----------------------------------------------------------
void CBlock::check_if_receive_particle_from_parent_is_required(void)
{

	if(!parent) return;

	//! 1)
	//! in case all child array do exist, block is surrounded by
	//! full set of neighbours, so no flag will be set
	//! 2)
	//! in case e.g. right child refined but not left, left face
	//! might be level borde so do a full check
	INT32 id_child =  Blk_Index[0]*2*2
			+Blk_Index[1]*2
			+Blk_Index[2];

	//! New Version:
	//! check wether any of 26 neighbours has no child array,
	//! in case yes, send one package of entire Block,
	//! else to many exeptions have to be interceped each time since particle
	//! may arrive from 
	//! (6 faces, 8 corners, 12 edges)


	//! TODO:
	//! "is remove permitted" genauer spezifizieren!!!


	//!	 x    x    x
	//!	 |	|jm1 |
	//!	 |	|    |
	//!	 | im1| ip1|
	//!	 x----X----x
	//!	 |	|    |
	//!	 |	|    |
	//!	 |	|jp1 |
	//!	 x    x    x


	//! return in case any of 6 "face" neighbours does not exist
	if(   !Neighbour[_Im1_] || !Neighbour[_Ip1_]
	    || !Neighbour[_Jm1_] || !Neighbour[_Jp1_]
	    || !Neighbour[_Km1_] || !Neighbour[_Kp1_])
	{

		//! root_of_array is the first of 8 blocks.
		//! since only the first block is checked
		//! whether to send particle, ALWAYS set 
		//! first blocks flag to true !!!
		do_receive_particle_from_parent = true;
		parent->do_send_particle_to_childArray[id_child] = true;
		return;
	}


	//! return in case any of 12 "edge" neighbours does not exist
	    //! i->j
	if(   !Neighbour[_Im1_]->Neighbour[_Jm1_] ||  !Neighbour[_Im1_]->Neighbour[_Jp1_]
	    || !Neighbour[_Ip1_]->Neighbour[_Jm1_] ||  !Neighbour[_Ip1_]->Neighbour[_Jp1_]
	    //! i->k
	    || !Neighbour[_Im1_]->Neighbour[_Km1_] ||  !Neighbour[_Im1_]->Neighbour[_Kp1_]
	    || !Neighbour[_Ip1_]->Neighbour[_Km1_] ||  !Neighbour[_Ip1_]->Neighbour[_Kp1_]
	    //! j->k
	    || !Neighbour[_Jm1_]->Neighbour[_Km1_] ||  !Neighbour[_Jm1_]->Neighbour[_Kp1_]
	    || !Neighbour[_Jp1_]->Neighbour[_Km1_] ||  !Neighbour[_Jp1_]->Neighbour[_Kp1_])
	{

		//! root_of_array is the first of 8 blocks.
		//! since only the first block is checked
		//! whether to send particle, ALWAYS set 
		//! first blocks flag to true !!!
		do_receive_particle_from_parent=true;
		parent->do_send_particle_to_childArray[id_child] = true;
		return;
		
	}

	//! return in case any of 8 "corner" neighbours does not exist

	if(    !Neighbour[_Im1_]->Neighbour[_Jm1_]->Neighbour[_Km1_]
	     || !Neighbour[_Im1_]->Neighbour[_Jm1_]->Neighbour[_Kp1_]
	     || !Neighbour[_Im1_]->Neighbour[_Jp1_]->Neighbour[_Km1_]
	     || !Neighbour[_Im1_]->Neighbour[_Jp1_]->Neighbour[_Kp1_]

	     || !Neighbour[_Ip1_]->Neighbour[_Jm1_]->Neighbour[_Km1_]
	     || !Neighbour[_Ip1_]->Neighbour[_Jm1_]->Neighbour[_Kp1_]
	     || !Neighbour[_Ip1_]->Neighbour[_Jp1_]->Neighbour[_Km1_]
	     || !Neighbour[_Ip1_]->Neighbour[_Jp1_]->Neighbour[_Kp1_])
	{

		//! root_of_array is the first of 8 blocks.
		//! since only the first block is checked
		//! whether to send particle, ALWAYS set 
		//! first blocks flag to true !!!
		do_receive_particle_from_parent=true;
		parent->do_send_particle_to_childArray[id_child] = true;
		return;
		
	}




}



