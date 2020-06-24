

#include <math.h>

#include "utils.h"
#include "CHybrid.h"
#include "parameters.h"
#include "hilbert.h"
#include "absolute_Globals.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"




//!------------------------------------------------------------------------
//!- LINEAR_Indices_to_BlkNr:
//!------------------------------------------------------------------------
INT32 LINEAR_Indices_to_BlkNr(INT32 ind[3])
{ 

	return   ind[0] *RB_Y *RB_Z
		+ind[1] *RB_Z
		+ind[2];
}


//!------------------------------------------------------------------------
//!- Linear_blk_nr_to_indices:
//!------------------------------------------------------------------------
void LINEAR_BlkNr_to_Indices(INT32 blk_nr, INT32 ind[3])
{ 
	ind[0] =  blk_nr  /      (RB_Y*RB_Z);
	ind[1] = (blk_nr -ind[0]*(RB_Y*RB_Z))  /     RB_Z;
	ind[2] =  blk_nr -ind[0]*(RB_Y*RB_Z) -ind[1]*RB_Z;
}


//!--------------------------------------------------------
//!- SFC_BlkNr_to_Indices:
//!--------------------------------------------------------
INT32 SFC_Indices_to_BlkNr(INT32 *indices, INT32 level)
{

	//! specify number of dimensions
	const int num_dims=3;
        bitmask_t Coord[3] = {indices[0], indices[1], indices[2]};


	//! NOTE:
	//! using "SFC_RB_power +level" instead of "RB_X" HAS TO BE TESTED
	return hilbert_c2i(num_dims, SFC_RB_power +level, Coord);

}

//!--------------------------------------------------------
//!- SFC_BlkNr_to_Indices:
//!--------------------------------------------------------
void SFC_BlkNr_to_Indices(INT32 BlkNr, INT32 *indices, INT32 level)
{


	//! specify number of dimensions
	const int num_dims=3;

        bitmask_t Coord[3];

	//! only one RB exist -> set indices to zero
	//! (RB_X=RB_Y=RB_Z has been checked before)
	//! (no SFC calculation required)
	if(RB_X==1)
	{
		memset(indices, 0, 3*sizeof(INT32));
		return;
	}


	//! NOTE:
	//! using "SFC_RB_power +level" instead of "RB_X" HAS TO BE TESTED
	//! transform index to coordinates
        hilbert_i2c(num_dims,
        	    SFC_RB_power +level,
        	    BlkNr,
        	    Coord);

	//! cp local copy to funtion parameter
	indices[0] = Coord[0];
	indices[1] = Coord[1];
	indices[2] = Coord[2];

}




//!------------------------------------------------------------------------
//!- vec_len:
//!------------------------------------------------------------------------
double vec_len(double *v)
{ 
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double vec_len(const double *v)
{ 
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

float vec_len(float *v)
{ 
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

//!------------------------------------------------------------------------
//!- vec_len2:
//!------------------------------------------------------------------------
float vec_len2(float *v)
{ 
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

//!------------------------------------------------------------------------
//!- vec_len2:
//!------------------------------------------------------------------------
double vec_len2(double *v)
{ 
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

//!------------------------------------------------------------------------
//!- vec_scalar:
//!------------------------------------------------------------------------

// template <typename T>
// T vec_scalar(T *v1, T *v2)
// { 
//   return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
// }

double vec_scalar(double *v1, double *v2)
{ 
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

// double vec_scalar(const double *v1, double *v2)
// { 
//   return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
// }

float vec_scalar(float *v1, float *v2)
{ 
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/*-------------------------------------------------------------------------*/
/* vec_cross                                                               */
/*-------------------------------------------------------------------------*/

void vec_cross(float *res, float *v1, float *v2)
{

  res[0] = v1[1]*v2[2]-v1[2]*v2[1];
  res[1] = v1[2]*v2[0]-v1[0]*v2[2];
  res[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

void vec_cross(double *res, double *v1, double *v2)
{

  res[0] = v1[1]*v2[2]-v1[2]*v2[1];
  res[1] = v1[2]*v2[0]-v1[0]*v2[2];
  res[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

void vec_cross(double *res, float *v1, double *v2)
{

  res[0] = v1[1]*v2[2]-v1[2]*v2[1];
  res[1] = v1[2]*v2[0]-v1[0]*v2[2];
  res[2] = v1[0]*v2[1]-v1[1]*v2[0];
}


void vec_cross(long double *res, long double *v1, long double *v2)
{

  res[0] = v1[1]*v2[2]-v1[2]*v2[1];
  res[1] = v1[2]*v2[0]-v1[0]*v2[2];
  res[2] = v1[0]*v2[1]-v1[1]*v2[0];
}


//!----------------------------------------------------
//!- calc_cells_to_cut: precalculate the XCell which is
//!    is within the cross section. Doing this once for
//!    for every level once in advance saves some time

//!   x_box = [0;1[   normalized position in Simu-Box
//! x_block = [0;1[   normalized position in block
//!  x_cell = [0;1[   normalized position in cell
//!----------------------------------------------------
bool calc_cells_to_cut(INT32 CS, INT32 *indedx_cell_L)
{


     D_REAL Length[3] = {LX,LY,LZ};
     INT32  num_RB[3]  = {RB_X,RB_Y,RB_Z};
     INT32  num_ND[3] = {BlkNds_X,BlkNds_Y,BlkNds_Z};

     INT32 index_RB;
     D_REAL x_box, x_blk, x_cell;

     //! get normalized position in Simu-Box
     x_box = (CROSS_SECTION[CS]+Box_Origin[CS])/Length[CS];
     if(x_box<0. || x_box >= 1.) return false;


     //!-------- Level 0 ---------------------------
     //! estimate block index to cut in root array and
     //! normalized position in root block
     index_RB = int(x_box*num_RB[CS]);
        x_blk =     x_box*num_RB[CS] - index_RB;



     //! estimate index of cell to cut in root block and
     //! normalized position in cell
     indedx_cell_L[0] = int(x_blk*(num_ND[CS]-2) +1);
               x_cell =    (x_blk*(num_ND[CS]-2) +1) -indedx_cell_L[0];

     //! since x_blk=[0;1[, only values [1,..,num_ND-2] can exist
     //! -> no GhostNodes can be hit in this round off style!!!
     //! in nearest neighbour style [1,..,num_ND-1] can exist 
     //! -> "plus" GhostNodes can be hit.
     //! Using nearest neighbour style is more convenient:
     if(x_cell>0.5)
     indedx_cell_L[0]++;


//     log_file << "index_RB: " << index_RB
// 	 << "  x_blk: " << x_blk
// 	 << "  indedx_cell_L[0]: " << indedx_cell_L[0]
// 	 << "  x_cell: " << x_cell << endl;

     //! avoid "plus" GhostNodes in nearest neighbour fashion
     if(avoid_GhostCXS && indedx_cell_L[0]==num_ND[CS]-1)
     indedx_cell_L[0]--;


     //! estimate cut-CELL for every level
     for(INT32 level=1; level<=MAX_LEVEL; level++)
     {

	//! x_blk=[0.0;0.5[ -> scale x_blk back to [0;1[
	if(x_blk<0.5)
	x_blk*=2.;
	//! x_blk=[0.5;1.0[ -> scale x_blk back to [0;1[
     	else 
	x_blk= 2.*x_blk -1.;

	//! estimate index of cell to cut in block of respective level
	//! and normalized position in cell
	indedx_cell_L[level] = int(x_blk*(num_ND[CS]-2) +1);
                      x_cell =    (x_blk*(num_ND[CS]-2) +1) -indedx_cell_L[level];

	//! Using nearest neighbour style is more convenient:
	if(x_cell>0.5)
	indedx_cell_L[level]++;

	//! avoid "plus" GhostNodes in nearest neighbour fashion
	if(avoid_GhostCXS && indedx_cell_L[level]==num_ND[CS]-1)
	indedx_cell_L[level]--;

     }

     return true;

}

//!-----------------------------------------------------
//!- calc_BlockNodeXCell_of_Pos:
//!    precalculate:
//!     index_blk[CS][level]
//!
//!   x_box = [  0;1[   normalized position in Simu-Box
//! x_block = [  0;1[   normalized position in block
//!   r_vec = [O-L:O[   normalized position in Simu-Box
//!-----------------------------------------------------
bool calc_BlockIndex_of_Pos(INT32**   index_blk,
		            D_REAL*  r_vec)
{


     D_REAL x_box, x_blk;

     D_REAL Length[3] = {LX,LY,LZ};
     INT32  num_RB[3] = {RB_X,RB_Y,RB_Z};


     for(INT32 comp=0; comp<3; comp++)
     {

		//! scale r_vec from normalized units to absolute positon in box [0;1[.
		x_box = (r_vec[comp]+Box_Origin[comp])/Length[comp];


		//! May occur and be correct in case of periodic boundaries
		if(x_box<0. || x_box >= 1.)
		{

			if(x_box < 0. && use_periodic_bounds[comp])
			x_box += 1.;

			if(x_box >= 1. && use_periodic_bounds[comp])
			x_box -= 1.;
			
			if(!use_periodic_bounds[comp])
			{

// 				log_file << " WARNING in 'calc_BlockIndex_of_Pos':" << endl;
// 				log_file << " requested position out of box (ignoring x=("<<r_vec[0]<<","
// 										<<r_vec[1]<<","
// 										<<r_vec[2]<<"))." << endl;
				return false;
			}
		}
		
		
		//!-------- Level 0 ---------------------------
		//! estimate block index of root array and
		//! normalized position in root block
		index_blk[comp][0] = int(x_box*num_RB[comp]);
			   x_blk =     x_box*num_RB[comp] - index_blk[comp][0];
		
		
		
		//! estimate block index for every level
		for(INT32 level=1; level<=MAX_LEVEL; level++)
		{
		
			//! x_blk=[0.0;0.5[ first half of parent block
			if(x_blk<0.5)
			{
				//! at level>0 first blk index = 0
				index_blk[comp][level] = 0;
				
				//! scale x_blk back to [0;1[
				x_blk*=2.;
			}
			//! x_blk=[0.5;1.0[ second half of parent block
			else 
			{
				//! at level>0 second blk index = 1
				index_blk[comp][level] = 1;
		
				//! scale x_blk back to [0;1[
				x_blk =2.*x_blk -1.;
			}
		}

     }

     return true;

}



//!----------------------------------------------------
//!- calc_BlockNodeXCell_of_Pos:
//!    precalculate:
//!     index_blk[CS][level]
//!    index_cell[CS][level]
//!        x_cell[CS][level]
//!
//!  requires:
//!   r_vec = [0;1[   normalized position in Simu-Box
//!
//!   x_box = [0;1[   normalized position in Simu-Box
//! x_block = [0;1[   normalized position in block
//!  x_cell = [0;1[   normalized position in cell
//!----------------------------------------------------
bool calc_BlockNodeXCell_of_Pos(INT32**  index_blk,
				INT32**  index_cell,
				D_REAL** x_cell,
				D_REAL*  r_vec)
{


	//! set values to arrays
	INT32  num_RB[3] = {RB_X,RB_Y,RB_Z};
	INT32  num_ND[3] = {BlkNds_X,BlkNds_Y,BlkNds_Z};
	
	
	D_REAL x_box, x_blk;
	
	for(INT32 CS=0; CS<3; CS++)
	{

		//! r_vec is normalized position in Simu-Box
//      	D_REAL Length[3] = {LX,LY,LZ};
// 		x_box = (r_vec[CS]+Box_Origin[CS])/Length[CS];
		x_box = r_vec[CS];
		if(x_box<0. || x_box >= 1.)
		{
// 			log_file << " WARNING in 'calc_BlockNodeXCell_of_Pos':" << endl;
// 			log_file << " requested position out of box (ignoring x=("<<r_vec[0]<<","
// 									<<r_vec[1]<<","
// 									<<r_vec[2]<<"))." << endl;
			return false;
		}
		
		
		//!-------- Level 0 ---------------------------
		//! estimate block index to cut in root array and
		//! normalized position in root block
		index_blk[CS][0] = int(x_box*num_RB[CS]);
			   x_blk =     x_box*num_RB[CS] - index_blk[CS][0];
		
		
		
		//! estimate index of cell to cut in root block and
		//! normalized position in cell
		index_cell[CS][0] = int(x_blk*(num_ND[CS]-2) +1);
		    x_cell[CS][0] =    (x_blk*(num_ND[CS]-2) +1) -index_cell[CS][0];
		
// 		    log_file << "   index_blk[0][0]: " << index_blk[CS][0]
// 			 << "   index_cell[0][0]: " <<  index_cell[CS][0]
// 			 << "   x_cell[0][0]: " << x_cell[CS][0] << endl;
		
		//! estimate cut-CELL for every level
		for(INT32 level=1; level<=MAX_LEVEL; level++)
		{
		
			//! x_blk=[0.0;0.5[ first half of parent block
			if(x_blk<0.5)
			{
				//! at level>0 first blk index = 0
				index_blk[CS][level] = 0;
				
				//! scale x_blk back to [0;1[
				x_blk*=2.;
			}
			//! x_blk=[0.5;1.0[ second half of parent block
			else 
			{
				//! at level>0 second blk index = 1
				index_blk[CS][level] = 1;
		
				//! scale x_blk back to [0;1[
				x_blk =2.*x_blk -1.;
		
			}
		
			//! estimate index of cell to cut in block of respective level
			//! and normalized position in cell
			index_cell[CS][level] = int(x_blk*(num_ND[CS]-2) +1);
			    x_cell[CS][level] =    (x_blk*(num_ND[CS]-2) +1) -index_cell[CS][level];

// 			log_file << "   index_blk["<<CS<<"]["<<level<<"]: " << index_blk[CS][level]
// 			     << "   index_cell["<<CS<<"]["<<level<<"]: " <<  index_cell[CS][level]
// 			     << "   x_cell["<<CS<<"]["<<level<<"]: " << x_cell[CS][level] << endl;
		}

     }

     return true;

}


//!-------------------------------------------------------------//
//! get_maxwellian_GSL:
//! faster maxwellian generation
//!-------------------------------------------------------------//
void get_maxwellian_GSL(INT32 species, PARTICLE_REAL *vth)
{

	if(!thermal_velocity_para_of[species] && !thermal_velocity_perp_of[species])
	return;

	if(thermal_velocity_perp_of[species])
	{	
		vth[0]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], thermal_velocity_perp_of[species]);
		vth[1]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], thermal_velocity_perp_of[species]);
	}
	else 
	{
		vth[0]=0;
		vth[1]=0;
	}
	
	if(thermal_velocity_para_of[species])
	vth[2]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], thermal_velocity_para_of[species]);
	else
	vth[2];

}

//!-------------------------------------------------------------//
//! get_maxwellian_GSL:
//! faster maxwellian generation
//!-------------------------------------------------------------//
void get_maxwellian_GSL(INT32 species, D_REAL ion_beta, PARTICLE_REAL *vth)
{

	const D_REAL norm_value = 1.;

	//! see S.Simon's PhD Thesis pp. 261
// 	norm_value = 1. /(Ion_Masses[species]*rho_sw[species]);
		

	D_REAL v2mean = 1.5 * ion_beta *norm_value;
	D_REAL thermal_velocity = sqrt(v2mean);

	D_REAL sigma = thermal_velocity * (1./sqrt(3.));

	vth[0]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], sigma);
	vth[1]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], sigma);
	vth[2]= gsl_ran_gaussian_ziggurat(randGen_of_species[species], sigma);

// 	const D_REAL max_vth = 10.;
// 
// 	if(vth[0] > max_vth || vth[0] < -max_vth) vth[0] = max_vth;
// 	if(vth[1] > max_vth || vth[1] < -max_vth) vth[1] = max_vth;
// 	if(vth[2] > max_vth || vth[2] < -max_vth) vth[2] = max_vth;

}


//!----------------------------------------------------------------
//! get_maxwellian_v:
//! maxwellian generation as in original TB code
//!----------------------------------------------------------------
// void get_maxwellian_v(INT32 species, PARTICLE_REAL *vth)
// {
// 
// 	if(!thermal_velocity_of[species])
// 	return;
// 
// 	int  comp = 0;
// 	PARTICLE_REAL f_maxwell, v2mean, vmean_x_MAX_VTH;
// 
// 	v2mean = thermal_velocity_of[species]*thermal_velocity_of[species];
// 	vmean_x_MAX_VTH = thermal_velocity_of[species]*max_vth;
// 
// 
//  	for (comp=0; comp<3; comp++)
// 	 do
// 	 {
// 
// 		//! [-MAX_VTH:+MAX_VTH]
// 		vth[comp]=(2.*random()/RAND_MAX-1.)*vmean_x_MAX_VTH;
// 
// 		f_maxwell=exp(-1.5*vth[comp]*vth[comp]/(v2mean));
// 
// 	 }
// 	 while (1.*random()/RAND_MAX > f_maxwell);
// 
// }





//!--------------------------------------------------------------
//!- copy_Field:
//!--------------------------------------------------------------
void CHybrid::copy_Field_L(INT32 level, INT32 dest, INT32 src)
{


	CBlock *temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->copy_Field(dest, src);
	
		temp_Block = temp_Block->next_Blk_of_BlockList;
	}

}

//!--------------------------------------------------------------
//!- copy_Field:
//!--------------------------------------------------------------
void CHybrid::copy_Field(INT32 dest, INT32 src)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
  	copy_Field_L(level, dest, src);


}

//!--------------------------------------------------------------
//!- copy_Field:
//!--------------------------------------------------------------
void CHybrid::add_Vector2Field(INT32 dest, INT32 src, D_REAL* vec)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			temp_Block->add_Vector2Field(dest, src, vec);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------------
//!- copy_Field:
//!--------------------------------------------------------------
void CHybrid::add_multipliedField(INT32 dest, INT32 src, D_REAL factor)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->add_multipliedField(dest, src, factor);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------------
//!- sub_squaredField:
//!- dest: scalr field 
//!- src: vector field
//!--------------------------------------------------------------
void CHybrid::sub_squaredField_Multiply(INT32 dest, INT32 src, D_REAL factor)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->sub_squaredField_Multiply(dest, src, factor);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------------
//!- Multiply_fields:
//!--------------------------------------------------------------
void CHybrid::Multiply_fields(INT32 dest, INT32 src1, INT32 src2)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->Multiply_fields(dest, src1, src2);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------------
//!- Multiply_field:
//!--------------------------------------------------------------
void CHybrid::Multiply_field(INT32 dest, INT32 src1, D_REAL src2)
{

    for(INT32 level=0; level<= MAX_LEVEL; level++)
    {

        CBlock *temp_Block = BlockList_of_Lev[level];
        while(temp_Block)
        {

            if(mpi_myRank == temp_Block->responsible_mpi_process)
                temp_Block->Multiply_field(dest, src1, src2);
            
            temp_Block = temp_Block->next_Blk_of_BlockList;
        }
    }

}

//!--------------------------------------------------------------
//!- Add_fields:
//!--------------------------------------------------------------
void CHybrid::Add_fields(INT32 dest, INT32 src1, INT32 src2)
{
    
    for(INT32 level=0; level<= MAX_LEVEL; level++)
    {
        
        CBlock *temp_Block = BlockList_of_Lev[level];
        while(temp_Block)
        {
            
            if(mpi_myRank == temp_Block->responsible_mpi_process)
                temp_Block->Add_fields(dest, src1, src2);
            
            temp_Block = temp_Block->next_Blk_of_BlockList;
        }
    }
    
}

//!--------------------------------------------------------------
//!- Devide_fields:
//!--------------------------------------------------------------
void CHybrid::Devide_fields(INT32 dest, INT32 src1, INT32 src2, D_REAL min)
{
    
    for(INT32 level=0; level<= MAX_LEVEL; level++)
    {
        
        CBlock *temp_Block = BlockList_of_Lev[level];
        while(temp_Block)
        {
            
            if(mpi_myRank == temp_Block->responsible_mpi_process)
                temp_Block->Devide_fields(dest, src1, src2, min);
            
            temp_Block = temp_Block->next_Blk_of_BlockList;
        }
    }
    
}



//!--------------------------------------------------------------
//!- set_zero_field:
//!--------------------------------------------------------------
void CHybrid::set_zero_field(INT32 dest)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			memset(temp_Block->Field_Type[dest], 0, COMPs_FType[dest]*num_nodes_in_block *sizeof(D_REAL));
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------------
//!- set_zero_field_inCore:
//!--------------------------------------------------------------
void CHybrid::set_zero_field_inCore(INT32 dest)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->set_zero_Field_inCore(dest);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------------
//!- Multiply_fields:
//!--------------------------------------------------------------
void CHybrid::set_zero_field_incGather(INT32 dest)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{

// 			if(mpi_myRank == temp_Block->responsible_mpi_process)
			//! some blocks that are not processed by myrank have field buffers
			//! that must be set to zero as well.

			if(temp_Block->Field_Type[dest])
			memset(temp_Block->Field_Type[dest], 0, COMPs_FType[dest]*num_nodes_in_block *sizeof(D_REAL));

			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}
	}

}

//!--------------------------------------------------------------
//!- square_Field:
//!--------------------------------------------------------------
void CHybrid::square_Field(INT32 dest, INT32 src)
{

  	for(INT32 level=0; level<= MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->square_Field(dest, src);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}
