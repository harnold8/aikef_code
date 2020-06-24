



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
particle* active_particle;

PARTICLE_REAL vth[3];

//! MPI related variables
extern ofstream log_file;


//! Random Generators
extern gsl_rng **randGen_of_species;

//! extern profiles
extern D_REAL **extern_1D_Field;

/*
//!--------------------------------------------------------------
//! get_a_b_c_from_abc:
//!--------------------------------------------------------------
void CBlock::get_a_b_c_from_abc(INT32 abc, INT32* a_b_c)
{
	INT32 temp = abc;

	a_b_c[0] = temp/(BlkNds_Y*BlkNds_Y);

	temp = abc%(BlkNds_Y*BlkNds_Z);

	a_b_c[1] = temp/(BlkNds_Z);

	a_b_c[2] = temp%(BlkNds_Z);


	if((a_b_c[0]*BlkNds_Y*BlkNds_Z +a_b_c[1]*BlkNds_Z +a_b_c[2]) != abc)
	{
			log_file << "false a_b_c" << endl;

			log_file << "abc:    " << abc << endl;
			log_file << "new_abc:    " << a_b_c[0]*BlkNds_Y*BlkNds_Z +a_b_c[1]*BlkNds_Z +a_b_c[2] << endl;

			log_file << "new_a:    " << a_b_c[0] << endl;
			log_file << "new_b:    " << a_b_c[1] << endl;
			log_file << "new_c:    " << a_b_c[2] << endl;
		exit(1);
	
	}
}*/





//!--------------------------------------------------------------
//! count_particle:
//!--------------------------------------------------------------
void CBlock::count_particle(void)
{

	//! reset particle nummer
	num_particle = 0;

	for(INT32 species=0; species<num_Charged_Species; species++)
	for(INT32 node=0; node<num_nodes_in_block; node++)
	num_particle += num_MPiC[species][node];

	num_particle_incChilds = num_particle;



}

//!--------------------------------------------------------------
//! resize_pArray:
//! increase or decrease
//!--------------------------------------------------------------
void CBlock::resize_pArray(INT32 species,
			   INT32 node,
			   INT32 new_size)
{
		//! - do not try allocate zero byte memory 
		//! - do nothing in case new_size==old_size
		if(!new_size ||  new_size == size_of_pArray[species][node])
		return;

		//! alloc new memory
		particle* buffer = new particle[new_size];
		num_pArrayReallocs++;

		if(num_MPiC[species][node] > new_size)
		{
			log_file << "EROOR in resize_pArray an node " << node << endl;
			log_file << "new_size: " << new_size << endl;
			log_file << "num_MPiC["<<species<<"]["<<node<<"]: " << num_MPiC[species][node]  << endl;
			exit(1);
		}


		//! copy and delete old memory in case it has been in use
		if(pArray[species][node])
		{
			memcpy(buffer,
			       pArray[species][node],
			       num_MPiC[species][node] *sizeof(particle));
	

			delete[] pArray[species][node];
		}

		num_particle_total_storage -= size_of_pArray[species][node];
		
		//! assign new to root of List
		pArray[species][node] = buffer;

		//! increase size_of_pArray
		size_of_pArray[species][node] = new_size;
		num_particle_total_storage += new_size;

}

//!--------------------------------------------------------------
//! resize_all_pArrays_of_Block:
//! increase or decres
//!--------------------------------------------------------------
void CBlock::resize_all_pArrays_of_Block(void)
{

	INT32 i,j,k,i_j_k;

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	INT32 min_pArray_Size = 2;

	for(INT32 species = 0; species<num_Charged_Species; species++)
// 	 for(i = is_box_boundary[0]; i < BlkNds_X-is_box_boundary[1]; i++)
// 	  for(j = is_box_boundary[2]; j < BlkNds_Y-is_box_boundary[3]; j++)
// 	   for(k = is_box_boundary[4]; k < BlkNds_Z-is_box_boundary[5]; k++)
	 for(i = 0; i < BlkNds_X; i++)
	  for(j = 0; j < BlkNds_Y; j++)
	   for(k = 0; k < BlkNds_Z; k++)
	   {
			
		i_j_k  = i*BlkNds_Y*BlkNds_Z 
			+j*BlkNds_Z 
			+k;


		//! decrease pArray Size in case MPiC is to low
		//! (use min_pArray_Size as lower limit)
		if( int(max_pArray_Size_factor*num_MPiC[species][i_j_k]) < size_of_pArray[species][i_j_k])
		{
			if(num_MPiC[species][i_j_k] > min_pArray_Size)
			resize_pArray(species, i_j_k, int(min_pArray_Size_factor*num_MPiC[species][i_j_k]));
// 			else
// 			if(num_MPiC[species][i_j_k]!=min_pArray_Size)
// 			resize_pArray(species, i_j_k, min_pArray_Size);
		}


		


	  }


	//! record particle calculation time
	time_finish = clock();
	time_process_particle += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}




//!--------------------------------------------------------------------------------
//! add_particle_to_pArray:
//! copy temp particle to pArray
//! species:	particle species
//! dest_cell:	destinatin cell of particle
//! particle_to_insert:		pointer or particle which will be copied to pArray
//!--------------------------------------------------------------------------------
void CBlock::add_particle_to_pArray(INT32 species,
				    INT32 dest_cell,
				    particle* particle_to_insert)
{

	//! check whether new list is large enough for additional particle
	if(size_of_pArray[species][dest_cell] < num_MPiC[species][dest_cell]+1)
	resize_pArray(species, dest_cell, int(min_pArray_Size_factor *(size_of_pArray[species][dest_cell]+1)));


#ifdef SORT_BY_WEIGHT
	//! estimate first particle in new list which weight is larger
	//! than particle_to_insert->weight
	//! START FROM END OF LIST
	//! -> Sort heavier or euqal weight Ions to end of list
	//!    allow for temporary negative index (case list is empty)
	INT32 new_part_index;
	

	for(new_part_index=num_MPiC[species][dest_cell]-1; new_part_index>=0; new_part_index--)
	 if( particle_to_insert->weight >= pArray[species][dest_cell][new_part_index].weight)
	 break;
	

	//! insert new particle behind estimated particle
	new_part_index++;

	//! --- INSERT IN NEW LIST -----------------------------------
	//! shift particle_list
	//! since memory overlaps, use memmove rather than memcpy.
	//! (memcpy may be undefined and IS NOT measurable faster)
	memmove(pArray[species][dest_cell] +(new_part_index+1),
		pArray[species][dest_cell] +(new_part_index),
		   (num_MPiC[species][dest_cell] -new_part_index)*sizeof(particle));

	//! copy temporary particle into pArray
	memcpy(pArray[species][dest_cell] +new_part_index,
	       particle_to_insert,
	       sizeof(particle));

#else

	//! copy temporary particle into pArray at ultimate position
	memcpy(pArray[species][dest_cell] +num_MPiC[species][dest_cell],
	       particle_to_insert,
	       sizeof(particle));

#endif



	//! update particle log for new cell
	num_MPiC[species][dest_cell]++;
		


}

//!-----------------------------------------------------------------------------
//! swap_particles_pArray:
//! species:	particle species
//! a_b_c:		old cell index
//! new_a_b_c:	new cell index
//! part_index:	particle position in old cells pArray
//! particle_to_insert:		pointer or particle which will be copied to pArray
//! 					has to be called by reference for correct update
//!-----------------------------------------------------------------------------
void CBlock::swap_particles_pArray(INT32 species,
				   INT32 part_index,
				   INT32 a_b_c,
				   INT32 new_a_b_c,
				   particle* &particle_to_insert)
{



	//! check whether new list is large enough for additional particle
	if(size_of_pArray[species][new_a_b_c] < num_MPiC[species][new_a_b_c]+1)
	resize_pArray(species,new_a_b_c,
			  int(min_pArray_Size_factor *(size_of_pArray[species][new_a_b_c]+1)));



#ifdef SORT_BY_WEIGHT
	//! estimate first particle in new list which weight is larger
	//! than particle_to_insert->weight
	//! START FROM END OF LIST
	//! -> Sort heavier or equal weight Ions to end of list
	//! allow for temporary negative index (case list is empty)
	INT32 new_part_index;
	
	for(new_part_index=num_MPiC[species][new_a_b_c]-1; new_part_index>=0; new_part_index--)
	 if( fabs(particle_to_insert->weight) >= fabs(pArray[species][new_a_b_c][new_part_index].weight))
	 break;
	

	//! insert new particle behind estimated particle
	new_part_index++;

	//! --- INSERT IN NEW LIST -----------------------------------
	//! since memory overlaps, use memmove rather than memcpy.
	//! (memcpy may be undefined and IS NOT measurable faster)
	memmove(pArray[species][new_a_b_c] +(new_part_index+1),
		pArray[species][new_a_b_c] +(new_part_index),
		   (num_MPiC[species][new_a_b_c] -(new_part_index))*sizeof(particle));

	memcpy(pArray[species][new_a_b_c] +new_part_index,
	       particle_to_insert,
	       sizeof(particle));

	//! update particle pointer
	particle_to_insert = pArray[species][new_a_b_c] +new_part_index;
		
	//! --- REMOVE FROM OLD LIST ---------------------------------
	//! --- !DO THIS AFTER PARTICLE HAS BEEN COPIED! -------------
	//! now the unhooked particle is copied, shift next up to last
	//! particle of list on top active_particle.
	//! since memory overlaps, use memmove rather than memcpy.
	//! (memcpy may be undefined and IS NOT measurable faster)
	memmove(pArray[species][a_b_c] +(part_index),
		pArray[species][a_b_c] +(part_index+1),
		   (num_MPiC[species][a_b_c] -(part_index+1))*sizeof(particle));

#else

	//! copy temporary particle into pArray at ultimate position
	memcpy(pArray[species][new_a_b_c] +num_MPiC[species][new_a_b_c],
	       particle_to_insert,
	       sizeof(particle));

	//! --- REMOVE FROM OLD LIST ---------------------------------
	//! --- !DO THIS AFTER PARTICLE HAS BEEN COPIED! -------------
	//! now the unhooked particle is copied, shift next up to last
	//! particle of list on top active_particle
	//! since memory overlaps, use memmove rather than memcpy.
	//! (memcpy may be undefined and IS NOT measurable faster)
	memmove(pArray[species][a_b_c] +(part_index),
		pArray[species][a_b_c] +(part_index+1),
		   (num_MPiC[species][a_b_c] -(part_index+1))*sizeof(particle));

#endif

	//! update particle logs for old and new cell
	num_MPiC[species][a_b_c]--;
	num_MPiC[species][new_a_b_c]++;


}


 
//!-------------------------------------------------------------//
//! estimate_extern_value: -								//
//!-------------------------------------------------------------//
void CBlock::estimate_extern_value(D_REAL &v, INT32 id_cord, INT32 id_field, PARTICLE_REAL *x)
{

	//! preset for extern profile
	const D_REAL extern_x0 = extern_1D_Field[id_cord][0];
	const D_REAL extern_dx = extern_1D_Field[id_cord][1] -extern_1D_Field[id_cord][0];


	//! set offset
	D_REAL temp_x = x[1] -extern_x0;

	//! estimate normed extern coordinate 
	D_REAL normed_ext_cord = temp_x/extern_dx;

	int ext_index = int(normed_ext_cord);
	D_REAL ext_cell_cord = normed_ext_cord - ext_index;

// 	v = extern_1D_Field[id_field][ext_index];
	//! interpolate extern velocity
	v = (1.-ext_cell_cord) *extern_1D_Field[id_field][ext_index]
	       +ext_cell_cord *extern_1D_Field[id_field][ext_index+1];

}

//!-------------------------------------------------------------//
//! fill_empty_Cell: -								//
//! faster than sorting particle into partially filled cell	//
//!-------------------------------------------------------------//
void CBlock::fill_empty_Cell(INT32 species, INT32 i_j_k, PARTICLE_REAL *v_init, PARTICLE_REAL rho_init)
{

	//! temp pointer to pArray position
	particle* new_particle;

	//! use the same weight for each particle
//         WEIGHT_REAL pWeight =  rho_sw[species]* CellVol_of_L[RLevel]/(1.*Blk_optimal_MPiC[species]);
        WEIGHT_REAL pWeight =  rho_init* CellVol_of_L[RLevel]/(1.*Blk_optimal_MPiC[species]);



	//! check whether cell is empty
	if(num_MPiC[species][i_j_k])
	{
		log_file << " Non empty cell must not be filled" << endl;
		log_file << " RLevel ... " << RLevel << endl;

		log_file << " i_j_k ... " << i_j_k << endl;
		log_file << " Exiting ... " << endl;
		exit(1);
// 		return;

	}	

	//! update statistics
	num_injected_particles		 += Blk_optimal_MPiC[species];
	num_MPiC[species][i_j_k]	 += Blk_optimal_MPiC[species];
	num_total_particles		 += Blk_optimal_MPiC[species];
	num_total_particles_in_L[RLevel] += Blk_optimal_MPiC[species];


	//! increase pArray size if required
	if(size_of_pArray[species][i_j_k] < num_MPiC[species][i_j_k])
	resize_pArray(species, i_j_k, int(min_pArray_Size_factor *num_MPiC[species][i_j_k]));

	if(size_of_pArray[species][i_j_k] < num_MPiC[species][i_j_k])
	{
		log_file << "num_MPiC[species][i_j_k]: " << num_MPiC[species][i_j_k] << endl;
		log_file << "size_of_pArray[species][i_j_k]: " << size_of_pArray[species][i_j_k] << endl;
		log_file << "int(min_pArray_Size_factor *num_MPiC[species][i_j_k]): " 
			 << int(min_pArray_Size_factor *num_MPiC[species][i_j_k]) << endl;
		exit(1);
	}





	//!------------------------------------------------------------
	//!------------------------------------------------------------
	//!------------------------------------------------------------
	INT32 cell_indices[3];
	PARTICLE_REAL abs_x1, x[3];
	PARTICLE_REAL sw_velocity[3];

	D_REAL width = 3.;

	cell_indices[0] =  i_j_k  /      (BlkNds_Y*BlkNds_Z);
	cell_indices[1] = (i_j_k -cell_indices[0]*(BlkNds_Y*BlkNds_Z))  /     BlkNds_Z;
	cell_indices[2] =  i_j_k -cell_indices[0]*(BlkNds_Y*BlkNds_Z) -cell_indices[1]*BlkNds_Z;



	memset(vth,0,3*sizeof(PARTICLE_REAL));

	//! fill set new particle
	//! NOTE: PARTICLE_LIST IS ASSUMED TO BE EMPTY
	//!      IF NOT OLD PARTICLES ARE DELETED
	for(INT32 z = 0; z  < Blk_optimal_MPiC[species];z++)
	{
		new_particle = pArray[species][i_j_k] +z;
		



		//! ---- old version ------------------------
// 		new_particle->rel_r[0] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
// 		new_particle->rel_r[1] =  1.*gsl_rng_uniform(randGen_general_asynchronous);
// 		new_particle->rel_r[2] =  1.*gsl_rng_uniform(randGen_general_asynchronous);

// 		get_maxwellian_v(species, vth);
		//! -----------------------------------------


		//! reset particle values (required in case of pTracing)
		memset(new_particle, 0, sizeof(particle));

		new_particle->rel_r[0] = gsl_rng_uniform(randGen_of_species[species]);
		new_particle->rel_r[1] = gsl_rng_uniform(randGen_of_species[species]);
		new_particle->rel_r[2] = gsl_rng_uniform(randGen_of_species[species]);


		get_maxwellian_GSL(species, vth);

		if(!special_Velocity_Distribution[species]) 
		{
		  //! Normal Maxwellian Distribution
		  new_particle->v[0] =  v_init[0] +vth[0];
		  new_particle->v[1] =  v_init[1] +vth[1];
		  new_particle->v[2] =  v_init[2] +vth[2];
		}
		else
		{
		  //! Special Distribution
		  //! This case Ring Distribution
		  D_REAL t = 2*M_PI*gsl_rng_uniform(randGen_of_species[species]);
		  D_REAL sinTheta = sin(B_angle*M_PI/180);
		  D_REAL cosTheta = cos(B_angle*M_PI/180);
                  new_particle->v[0] =  v_init[0]*sinTheta*(1-cos(t)) +vth[0];
                  new_particle->v[1] =  -1.*v_init[0]*sinTheta*cosTheta*(1-cos(t)) +vth[1];
                  new_particle->v[2] =  -1.*v_init[0]*sinTheta*(sin(t)) +vth[2];
		}

		new_particle->weight = pWeight;
		//! --------------------------------------

#ifdef TRACK_PARTICLE
		new_particle->number = -1;//local_mark_particle_counter;

		//local_mark_particle_counter ++;
#endif

// 		if(x[1]>5.)
// 		new_particle->v[0] = V_sw[0] +vth[0];
// 
// 		if(x[1]<-7.5)
// 		new_particle->v[0] = -V_sw[0] +vth[0];


// 		abs_x1 = fabs(x[1]);
// 
// 		if(abs_x1 < width)
// 		new_particle->v[0] =  sin(x[1]/width *0.5 *M_PI) *V_sw[0] +vth[0];
// 		else
// 		new_particle->v[0] =  x[1]/abs_x1 *V_sw[0] +vth[0];


// 		if(1.*gsl_rng_uniform(randGen_general_asynchronous) < 0.5)
// 		new_particle->v[0] +=5.;
// 		else
// 		new_particle->v[0] -=5.;

// 		new_particle->time = move_time[RLevel];

	}

}


//!-------------------------------------------------------------//
//! copy_empty_Cell: -								//
//! faster than sorting particle into partially filled cell	//
//!-------------------------------------------------------------//
void CBlock::copy_empty_Cell(INT32 species, INT32 dest_i_j_k, INT32 src_i_j_k)
{




	//! check whether cell is empty
	if(num_MPiC[species][dest_i_j_k])
	{
		log_file << " Non empty cell must not be filled" << endl;
		log_file << " RLevel ... " << RLevel << endl;

		log_file << " dest_i_j_k ... " << dest_i_j_k << endl;
		log_file << " Exiting ... " << endl;
		exit(1);

	}	

	//! update statistics
	num_injected_particles		 += num_MPiC[species][src_i_j_k];
	num_MPiC[species][dest_i_j_k]	 += num_MPiC[species][src_i_j_k];
	num_total_particles		 += num_MPiC[species][src_i_j_k];
	num_total_particles_in_L[RLevel] += num_MPiC[species][src_i_j_k];


	//! increase pArray size if required
	if(size_of_pArray[species][dest_i_j_k] < num_MPiC[species][dest_i_j_k])
	resize_pArray(species, dest_i_j_k, int(min_pArray_Size_factor *num_MPiC[species][dest_i_j_k]));

	if(size_of_pArray[species][dest_i_j_k] < num_MPiC[species][dest_i_j_k])
	{
		log_file << "num_MPiC[species][dest_i_j_k]: " << num_MPiC[species][dest_i_j_k] << endl;
		log_file << "size_of_pArray[species][dest_i_j_k]: " << size_of_pArray[species][dest_i_j_k] << endl;
		log_file << "int(min_pArray_Size_factor *num_MPiC[species][dest_i_j_k]): " 
			 << int(min_pArray_Size_factor *num_MPiC[species][dest_i_j_k]) << endl;
		exit(1);
	}

	INT32 temp1 = num_MPiC[species][src_i_j_k];
	INT32 temp2 = num_MPiC[species][dest_i_j_k];

	memcpy( pArray[species][dest_i_j_k],  pArray[species][src_i_j_k], num_MPiC[species][dest_i_j_k]*sizeof(particle));


}




//!-------------------------------------------------------------//
//! fill_empty_Volume: -							//
//!		//
//!-------------------------------------------------------------//
void CBlock::fill_empty_Volume(INT32 species, const INT32* start, const INT32* end)
{

		PARTICLE_REAL v_init[3];
		PARTICLE_REAL rho_init;

		v_init[0] = V_sw[0];
		v_init[1] = V_sw[1];
		v_init[2] = V_sw[2];

		rho_init = rho_sw[species];

		//! only use inhom boundaries when specified
		//! and at initial filling
		//! (inhom boundaries are filled by means of)
		//! (function: "apply_inhom_particle_bounds")
		if(use_hom_particle_bounds[_INIT_] || TL>0)
		{
			for(INT32 i = start[0]; i < end[0]; i++)
			 for(INT32 j = start[1]; j < end[1]; j++)
			  for(INT32 k = start[2]; k < end[2]; k++)
			  {
				
				INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
					+j*BlkNds_Z 
					+k;
			
				//! in case oct is plasma cell, fill cell with inflow species
				if(!Flag[i_j_k])
				fill_empty_Cell(species, i_j_k, v_init, rho_init);
			
			  }
		}
		else
		{

			INT32 ind[3];

			PARTICLE_REAL x[3];
			PARTICLE_REAL v_init[3] = {0., 0., 0.};

			for(ind[0] = start[0]; ind[0] < end[0]; ind[0]++)
			 for(ind[1] = start[1]; ind[1] < end[1]; ind[1]++)
			  for(ind[2] = start[2]; ind[2] < end[2]; ind[2]++)
			  {
				
				INT32 i_j_k  = ind[0]*BlkNds_Y*BlkNds_Z 
					     +ind[1]*BlkNds_Z 
					     +ind[2];

				cell_centre_normedCoords(x, ind);

				set_inflow_velocity(v_init,  x, species);
				set_inflow_density(rho_init, x, species);
			
				//! in case oct is plasma cell, fill cell with inflow species
				if(!Flag[i_j_k])
				fill_empty_Cell(species, i_j_k, v_init, rho_init);
			
			  }
		}
		
}

//!-------------------------------------------------------------//
//! copy_empty_Volume: -							//
//!		//
//!-------------------------------------------------------------//
void CBlock::copy_empty_Volume(INT32 species, const INT32* start, const INT32* end, const INT32* offset)
{

		for(INT32 i = start[0]; i < end[0]; i++)
		 for(INT32 j = start[1]; j < end[1]; j++)
		  for(INT32 k = start[2]; k < end[2]; k++)
		  {
			
			INT32 dest_i_j_k  =  i*BlkNds_Y*BlkNds_Z 
				           +j*BlkNds_Z 
				           +k;

			INT32 src_i_j_k  =   (i+offset[0])*BlkNds_Y*BlkNds_Z 
				           +(j+offset[1])*BlkNds_Z 
				           +(k+offset[2]);
		
			
			copy_empty_Cell(species, dest_i_j_k, src_i_j_k);
		
		  }
}


//!-------------------------------------------------------------//
//! fill_Block: -							//
//!		//
//!-------------------------------------------------------------//
void CBlock::fill_Oct(INT32 id_oct)
{


	
	//! get index of oct (0 or 1 for each dimension)
	INT32 a = id_oct/4;
	INT32 b = (id_oct -a*4)/2;
	INT32 c = (id_oct -a*4 -b*2);
	
	
	//! In case Block is Box boundary, filling shall start at 0, else at 1
	INT32 start[3] = {!a*!is_box_boundary[_Im1_] +a*(BlkNds_X/2),
			 !b*!is_box_boundary[_Jm1_] +b*(BlkNds_Y/2),
			 !c*!is_box_boundary[_Km1_] +c*(BlkNds_Z/2)};


	INT32 end[3]={BlkNds_X-1 -!a*((BlkNds_X-2)/2),
		     BlkNds_Y-1 -!b*((BlkNds_Y-2)/2),
		     BlkNds_Z-1 -!c*((BlkNds_Z-2)/2)};
	
	//! If Block is no Box-Boundary the following is valid:
	//! 1) start from 1 as 0 is not accesable by particles (GC)
	//! 2) Insert up to Max-1, as no cell is associated with Max
	for(INT32 inflow_species = 0; inflow_species<num_Inflow_Species; inflow_species++)
	fill_empty_Volume(index_Inflow_Species[inflow_species], start, end);

}


 
//!-------------------------------------------------------------//
//! intern2normedUnits: -							//
//!-------------------------------------------------------------//
void CBlock::cell_centre_normedCoords(double* x, INT32* cell_indices)
{


	//!---------------------------------------------------
	//! Grid Generation of Adaptive Mesh Code
	//!---------------------------------------------------
	x[0] = origin[0] + (0.5 +cell_indices[0] -1)*delta_of_L[RLevel][0];
	x[1] = origin[1] + (0.5 +cell_indices[1] -1)*delta_of_L[RLevel][1];
	x[2] = origin[2] + (0.5 +cell_indices[2] -1)*delta_of_L[RLevel][2];


}



//!-------------------------------------------------------------//
//! intern2normedUnits: -							//
//!-------------------------------------------------------------//
void CBlock::intern2normedCoords(double* x, double* cell_intern_r, INT32* cell_indices)
{


	//!---------------------------------------------------
	//! Grid Generation of Adaptive Mesh Code
	//!---------------------------------------------------
	x[0] = origin[0] + (cell_intern_r[0] +cell_indices[0] -1)*delta_of_L[RLevel][0];
	x[1] = origin[1] + (cell_intern_r[1] +cell_indices[1] -1)*delta_of_L[RLevel][1];
	x[2] = origin[2] + (cell_intern_r[2] +cell_indices[2] -1)*delta_of_L[RLevel][2];

	//!---------------------------------------------------
	//! Grid Generation as in original TB Code
	//!---------------------------------------------------
// 	x[0] = LX*((cell_intern_r[0] +1.*cell_indices[0])/(BlkNds_X-2.)-0.5);
// 	x[1] = LY*((cell_intern_r[1] +1.*cell_indices[1])/(BlkNds_Y-2.)-0.5);
// 	x[2] = LZ*((cell_intern_r[2] +1.*cell_indices[2])/(BlkNds_Z-2.)-0.5);

}

//!-------------------------------------------------------------//
//! intern2normedUnits: -							//
//!-------------------------------------------------------------//
void CBlock::intern2normedCoords_HIMesh(double* x, double* cell_intern_r, INT32* cell_indices)
{


	//!---------------------------------------------------
	//! Grid Generation of Adaptive Mesh Code
	//!---------------------------------------------------
	x[0] = origin[0] + (1.*cell_intern_r[0] +1.*cell_indices[0] -0.5)*delta_of_L[RLevel][0];
	x[1] = origin[1] + (1.*cell_intern_r[1] +1.*cell_indices[1] -0.5)*delta_of_L[RLevel][1];
	x[2] = origin[2] + (1.*cell_intern_r[2] +1.*cell_indices[2] -0.5)*delta_of_L[RLevel][2];


}

//!-------------------------------------------------------------//
//! normed2internUnits: -							//
//!-------------------------------------------------------------//
void CBlock::normed2internCoords(PARTICLE_REAL* intern_r, PARTICLE_REAL* x)
{



	//!---------------------------------------------------
	//! Grid Generation of Adaptive Mesh Code
	//!---------------------------------------------------

	//! Get Position in Block
	intern_r[0] = x[0] - origin[0];
	intern_r[1] = x[1] - origin[1];
	intern_r[2] = x[2] - origin[2];

 	//! Estimate intern index
	intern_r[0] = intern_r[0]/delta_of_L[RLevel][0] +1.;
	intern_r[1] = intern_r[1]/delta_of_L[RLevel][1] +1.;
	intern_r[2] = intern_r[2]/delta_of_L[RLevel][2] +1.;

			
	//!---------------------------------------------------
	//! Grid Generation as in original TB Code
	//!---------------------------------------------------
// 	intern_r[0] = (x[0]/LX  +0.5) * (BlkNds_X-2.);
// 	intern_r[1] = (x[1]/LY  +0.5) * (BlkNds_Y-2.);
// 	intern_r[2] = (x[2]/LX  +0.5) * (BlkNds_Z-2.);





}




