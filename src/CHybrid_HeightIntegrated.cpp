#include "parameters.h"
#include "CHybrid.h"
#include "absolute_Globals.h"
#include "utils.h"
#include "CBlk.h"  
#include <cmath>

extern CHybrid* Hybrid;

//! -----------------------------------------------------------------------------------------------
//! -----------------------------------------------------------------------------------------------
//! The three functions
// void getField_on_uniform_grid(INT32 Field_ID)
// void CHybrid::get_height_integrated_Field(INT32 src_field_id, INT32 dest_field_id)
// void CBlock::set_Hint_to_field(INT32 src_field_id,
//			       INT32 dest_field_id,
// 			      short  *num_uniform_Nodes,
// 			      D_REAL *extern_Origin,
// 			      D_REAL *extern_Length,
// 			      D_REAL *Hint_array)
//! write the height integrated value of src_field_id into all cross sections of
//! dest_field_id
//! NOTE: prototypes for last two have to be defined in CHybrid.h and CBlk.h, respectively
//! (last in CBlk_HeightIntegrated.cpp )
//! Usage example:
// 	set_zero_field(id_rotB);
// 	get_height_integrated_Field(id_BTotal,id_rotB);
//  	silo_3DWrite_Field(id_rotB);
//! BField is integrated and written into rotB
//! functions work as following: (1) assuming a uniform cartesian mesh, field is
//! 				     written into one array
//!				 (2) summation along axis
//!				 (3) summed field is written to simu mesh
//! -----------------------------------------------------------------------------------------------
//! -----------------------------------------------------------------------------------------------

//! -----------------------------------------------------------
//! parameters for height integrated field
//! -----------------------------------------------------------
//! choose axis for integration
//! NOTE: only one entry must be true
bool integration_axis[3] = {false,false,true};  

//! if integration shall only be performed for e.g. z>0, set origin[2] = 0;
//! or if integration shall only be performed for e.g. z<0, set origin[2] = Box_Origin[2];
//! if integration shall only be performed for whole Box length, set origin[2] = Box_Origin[2];
D_REAL origin_for_hInt[3] = {Box_Origin[0],Box_Origin[1],Box_Origin[2]};

//! length of integration axis
//! if integration shall only be performed for e.g. z>0, set LBox_for_hInt[2] = LZ-Box_Origin[2];
//! or if integration shall only be performed for e.g. z<0, set LBox_for_hInt[2] = Box_Origin[2];
//! if integration shall only be performed for whole Box length, set LBox_for_hInt[2] = LZ;
D_REAL LBox_for_hInt[3] = {LX,LY,Box_Origin[2]}; 

//! -----------------------------------------------------------
			    
			    
short num_nds_uniform[3] = {RB_X*(BlkNds_X-2)*pow(2,MAX_LEVEL),
				RB_Y*(BlkNds_Y-2)*pow(2,MAX_LEVEL),
				RB_Z*(BlkNds_Z-2)*pow(2,MAX_LEVEL)};
//! values larger than the number of cell in the respective direction don't make sense...
//! -----------------------------------------------------------


const INT32 num_grid_nds_uniform =  num_nds_uniform[0] * num_nds_uniform[1] * num_nds_uniform[2];

D_REAL* field_array_height_int;

//!------------------------------------------------------------
//!- get field on uniform mesh (in array):
//!------------------------------------------------------------
void getField_on_uniform_grid(INT32 Field_ID)
{

	log_file.flush();
	log_file << "     Get Field "<< Field_Name[Field_ID] << " ...";
	log_file.flush();


	INT32 u_v_w;

        short num_comps = COMPs_FType[Field_ID];

	//! indices for interpolation
	INT32  i,j,k;
	INT32  i_j_k, ip1_jp1_kp1;
	INT32  ip1_j_k, i_jp1_k,  i_j_kp1;
	INT32  ip1_jp1_k, ip1_j_kp1,  i_jp1_kp1;

	INT32  *index_blk[3], *index_cell[3], blk_index;
	D_REAL *x_cell[3], r_vec[3], A[3], shape_func[8], *BlkField[3];

	CBlock *temp_Block, *top_level_Block;

	//!-----------------------------------
	//! 1) alloc variable memory
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		 index_blk[comp] = new INT32[MAX_LEVEL+1];
		index_cell[comp] = new INT32[MAX_LEVEL+1];
		    x_cell[comp] = new D_REAL[MAX_LEVEL+1];
	}



	//!-----------------------------------
	//! 2) trace
	//!-----------------------------------

	for(INT32 u=0; u < num_nds_uniform[0]  ; u++)
	  for(INT32 v=0; v < num_nds_uniform[1] ; v++)
	    for(INT32 w=0; w < num_nds_uniform[2] ; w++)
	    {
		
		u_v_w  = u*num_nds_uniform[1]*num_nds_uniform[2]
			+v*num_nds_uniform[2]
			+w;

			
		//! same as for grid nodes
		//! TODO 0.99 weil sonst außerhalb der Box
		r_vec[0] = - origin_for_hInt[0]  + u * LBox_for_hInt[0]/(num_nds_uniform[0]-0.99) ;
		r_vec[1] = - origin_for_hInt[1]  + v * LBox_for_hInt[1]/(num_nds_uniform[1]-0.99) ;
		r_vec[2] = - origin_for_hInt[2]  + w * LBox_for_hInt[2]/(num_nds_uniform[2]-0.99) ;


		//! normalize coordinate to box length
		//! get normalized position in Simu-Box
		r_vec[0] = (r_vec[0]+Box_Origin[0])/LX;
		r_vec[1] = (r_vec[1]+Box_Origin[1])/LY;
		r_vec[2] = (r_vec[2]+Box_Origin[2])/LZ;
	
		
		//!--------------------------------------------------
		//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
		//!     (even though target blk might not be refined)
		//!--------------------------------------------------
	
		//! calculate indices and check whether position is in Simu Box
		if(calc_BlockNodeXCell_of_Pos(index_blk, index_cell, x_cell, r_vec))
		{


			//!--------------------------------------------------
			//! 2b) climb to highest possible level at respective
			//!    position.
			//!--------------------------------------------------
		
	
			//! copy indices to temp array
			INT32 blk_indices[3] = {index_blk[0][0],
					       index_blk[1][0],
					       index_blk[2][0]};

			
			//! first calculate root block
			//! -> always level 0
			if(use_SFC)
			blk_index = SFC_Indices_to_BlkNr(blk_indices, 0);
			else
			blk_index = LINEAR_Indices_to_BlkNr(blk_indices);
		
		
			//! initialize to respective root block
			temp_Block      = Root_Block_Array +blk_index;
			top_level_Block = Root_Block_Array +blk_index;;
		
			//! descent to highest existing level
			for(INT32 level=1; level<=MAX_LEVEL; level++)
			{
		
				//! set block index to precalculated indices
				blk_index = index_blk[0][level] *2 *2
					   +index_blk[1][level] *2
					   +index_blk[2][level];
				
		
				//! Decent by ordanary child array
				temp_Block = temp_Block->child_array[blk_index];
		
				//! cancel in case block does not exist
				if(!temp_Block)
				break;
		
				//! Set top_level_Block to temp_Block
				top_level_Block = temp_Block;
		
				
			}
	

			//!-----------------------------------------------------
			//! 2c) Now target Block and Cell are known, interpolate
			//!     Field at respective position.
			//!-----------------------------------------------------
	
			if(mpi_myRank == top_level_Block->responsible_mpi_process)
			{
				//! use pointers to block's field 
				for(INT32 comp=0; comp<COMPs_FType[Field_ID]; comp++)
				BlkField[comp] = top_level_Block->Field_Type[Field_ID] +comp *num_nodes_in_block;
				
		
				INT32 level=top_level_Block->RLevel;
		
				//! use common i,j,k notation:
				i = index_cell[0][level];
				j = index_cell[1][level];
				k = index_cell[2][level];
		
				//! ------------------------------------------
				i_j_k   =      i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
				
				//! ------------------------------------------
				ip1_j_k = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z     +k;
				i_jp1_k =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
				i_j_kp1 =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
				
				//! -------------------------------------------
				ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z     +k;
				ip1_j_kp1 = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +(k+1);
				i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
				
				ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +(k+1);
		
		
				//! use now r_vec for psition in respective cell
				r_vec[0] = x_cell[0][level];
				r_vec[1] = x_cell[1][level];
				r_vec[2] = x_cell[2][level];
		
				shape_func[0] = (1.-r_vec[0])*(1.-r_vec[1])*(1.-r_vec[2]);
				shape_func[1] = (   r_vec[0])*(1.-r_vec[1])*(1.-r_vec[2]);
				shape_func[2] = (1.-r_vec[0])*(   r_vec[1])*(1.-r_vec[2]);
				shape_func[3] = (1.-r_vec[0])*(1.-r_vec[1])*(   r_vec[2]);
				
				shape_func[4] = (   r_vec[0])*(   r_vec[1])*(1.-r_vec[2]);
				shape_func[5] = (   r_vec[0])*(1.-r_vec[1])*(   r_vec[2]);
				shape_func[6] = (1.-r_vec[0])*(   r_vec[1])*(   r_vec[2]);
				shape_func[7] = (   r_vec[0])*(   r_vec[1])*(   r_vec[2]);
	

				//! interpolate and store field in variable A
				for(INT32 comp=0; comp<COMPs_FType[Field_ID]; comp++)
				{
		
				      A[comp] = BlkField[comp][  i_j_k] * shape_func[0]
						+BlkField[comp][ip1_j_k] * shape_func[1]
						+BlkField[comp][i_jp1_k] * shape_func[2]
						+BlkField[comp][i_j_kp1] * shape_func[3]
						+BlkField[comp][  ip1_jp1_k] * shape_func[4]
						+BlkField[comp][  ip1_j_kp1] * shape_func[5]
						+BlkField[comp][  i_jp1_kp1] * shape_func[6]
						+BlkField[comp][ip1_jp1_kp1] * shape_func[7];
		

					field_array_height_int[comp*num_grid_nds_uniform+u_v_w] = A[comp];
				}
			
			}//! end if responsible_mpi_process

		}

	     }//! end for positions

	log_file << " done." <<endl;

}


//!------------------------------------------------------------
//!- get_height_integrated_Field
//!------------------------------------------------------------
void CHybrid::get_height_integrated_Field(INT32 src_field_id, INT32 dest_field_id)
{


	short axis;
	if(integration_axis[0])
		axis = 0;
	else if(integration_axis[1])
		axis = 1;
	else if(integration_axis[2])
		axis = 2;
	
	if(integration_axis[0]+integration_axis[1]+integration_axis[2] != 1)
	{
		log_file << " Error: can only integrate along ONE axis ... "<<endl;
		return;
	}	
	
	log_file << endl << " Integrating field  "<< Field_Name[src_field_id]<<" along axis " << axis << " ..." << endl;

	//! alloc array
	INT32 num_values = COMPs_FType[src_field_id]*num_grid_nds_uniform;
	field_array_height_int = new D_REAL[num_values];
	memset(field_array_height_int,0,num_values*sizeof(D_REAL));

	//! calc field on uniform grid
	getField_on_uniform_grid(src_field_id);
	
	
	Hybrid->synchronize_allProcesses_MPI();


	//! below, the height integration is performed
	
	INT32 u_v_w;
	D_REAL uniform_delta_g[3];
	
	for(INT32 comp=0; comp<3; comp++)
	uniform_delta_g[comp] = LBox_for_hInt[comp]/(num_nds_uniform[comp]-1);
		
	for(short comps=0; comps<COMPs_FType[src_field_id]; comps++)
	{	
		if(integration_axis[0])
		{	
			//! calculate coordinates of mesh nodes
			for(INT32 v=0; v < num_nds_uniform[1] ; v++)
			 for(INT32 w=0; w < num_nds_uniform[2] ; w++)
			 {	 
			
			 D_REAL sum_along_u = 0;	
			
			 //! get height_integrated quantity
			 for(INT32 u=0;u < num_nds_uniform[0] ; u++)
			 {
		
				u_v_w  = u*num_nds_uniform[1]*num_nds_uniform[2]
					+v*num_nds_uniform[2]
					+w;
							
				sum_along_u += field_array_height_int[u_v_w + comps*num_grid_nds_uniform ];	
			 }
			 //! set height integrated quantity to all cross sections
			 for(INT32 u=0; u < num_nds_uniform[0] ; u++)
			 {
		
				u_v_w  = u*num_nds_uniform[1]*num_nds_uniform[2]
					+v*num_nds_uniform[2]
					+w;
				
				//! if in arbitrary units	
// 				field_array_height_int[u_v_w + comps*num_grid_nds_uniform] = sum_along_u;
				//! if in physical units:
				field_array_height_int[u_v_w + comps*num_grid_nds_uniform] = sum_along_u*uniform_delta_g[0];
			 }
			}
		}	
		
		if(integration_axis[1])
		{	
			//! calculate coordinates of mesh nodes
			for(INT32 u=0; u < num_nds_uniform[0] ; u++)
			 for(INT32 w=0; w < num_nds_uniform[2] ; w++)
			 {	 
			
			 D_REAL sum_along_v = 0;	
			
			 //! get height_integrated quantity
			 for(INT32 v=0;v < num_nds_uniform[1] ; v++)
			 {
		
				u_v_w  = u*num_nds_uniform[1]*num_nds_uniform[2]
					+v*num_nds_uniform[2]
					+w;
							
				sum_along_v += field_array_height_int[u_v_w + comps*num_grid_nds_uniform ];	
			 }
			 //! set height integrated quantity to all cross sections
			 for(INT32 v=0; v < num_nds_uniform[1] ; v++)
			 {
		
				u_v_w  = u*num_nds_uniform[1]*num_nds_uniform[2]
					+v*num_nds_uniform[2]
					+w;
				
				//! if in arbitrary units	
// 				field_array_height_int[u_v_w + comps*num_grid_nds_uniform] = sum_along_v;
				//! if in physical units:
				field_array_height_int[u_v_w + comps*num_grid_nds_uniform] = sum_along_v*uniform_delta_g[1];
			 }
			}
		}
		
		if(integration_axis[2])
		{	
			//! calculate coordinates of mesh nodes
			for(INT32 u=0; u < num_nds_uniform[0] ; u++)
			 for(INT32 v=0; v < num_nds_uniform[1] ; v++)
			 {	 
			
			 D_REAL sum_along_w = 0;	
			
			 //! get height_integrated quantity
			 for(INT32 w=0; w < num_nds_uniform[2] ; w++)
			 {
		
				u_v_w  = u*num_nds_uniform[1]*num_nds_uniform[2]
					+v*num_nds_uniform[2]
					+w;
							
				sum_along_w += field_array_height_int[u_v_w + comps*num_grid_nds_uniform ];	
			 }
			 //! set height integrated quantity to all cross sections
			 for(INT32 w=0; w < num_nds_uniform[2] ; w++)
			 {
		
				u_v_w  = u*num_nds_uniform[1]*num_nds_uniform[2]
					+v*num_nds_uniform[2]
					+w;
				
				//! if in arbitrary units	
// 				field_array_height_int[u_v_w + comps*num_grid_nds_uniform] = sum_along_w;
				//! if in physical units:
				field_array_height_int[u_v_w + comps*num_grid_nds_uniform] = sum_along_w*uniform_delta_g[2];
			 }
			}
		}
	}

	
	mpi_build_sum(field_array_height_int,field_array_height_int,num_values);				     
		     
	
	//!---------------------------------------------------------//
	//! set height integrated field to simulation mesh for plotting
	//!---------------------------------------------------------//
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->set_Hint_to_field(src_field_id,
						      dest_field_id, 
						       num_nds_uniform,
						       origin_for_hInt,
						       LBox_for_hInt,
						       field_array_height_int,
						      integration_axis);

		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	

	delete[] field_array_height_int;

	log_file  << " done." << endl << endl;

}




