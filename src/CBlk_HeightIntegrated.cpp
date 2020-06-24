#include "CBlk.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"
#include "defines.h"
#include <math.h>

//!-------------------------------------------------------------//
//! set_Hint_to_field: slightly modified version of set_RhoUi_extern								//
//!-------------------------------------------------------------//
void CBlock::set_Hint_to_field(INT32 src_field_id,
			       INT32 dest_field_id,
			      short  *num_uniform_Nodes,
			      D_REAL *origin_for_hInt,
			      D_REAL *LBox_for_hInt,
			      D_REAL *Hint_array,
			      bool *integration_axis)
{

	INT32  ind[3];
	INT32  num_uniform_box_nodes;
	INT32 a,b,c, a_b_c, i_j_k;
	INT32 ap1_b_c, a_bp1_c, a_b_cp1;
	INT32 ap1_bp1_c, ap1_b_cp1, a_bp1_cp1, ap1_bp1_cp1;

	D_REAL *Hint_dest_field, *Hint_dest_fieldX, *Hint_dest_fieldY, *Hint_dest_fieldZ;
// 	D_REAL *extern_UiX, *extern_UiY, *extern_UiZ;
	D_REAL r[3], extern_delta[3], shape_func[8];

	PARTICLE_REAL x_BlockNode[3], x_extern[3], cell_intern_r[3];
	memset(cell_intern_r,0,3*sizeof(PARTICLE_REAL));


	for(INT32 comp=0; comp<3; comp++)
	extern_delta[comp] = LBox_for_hInt[comp]/(num_uniform_Nodes[comp]-1);


	//! Set pointer to extern field
	num_uniform_box_nodes =   num_uniform_Nodes[0]
				*num_uniform_Nodes[1]
				*num_uniform_Nodes[2];

	
	
	//! Set pointer to intern (Block) field
	if(COMPs_FType[src_field_id]==1)
	Hint_dest_field = Field_Type[dest_field_id];

	if(COMPs_FType[src_field_id]==3)
	{
		Hint_dest_fieldX = Field_Type[dest_field_id];
		Hint_dest_fieldY = Hint_dest_fieldX + num_nodes_in_block;
		Hint_dest_fieldZ = Hint_dest_fieldY + num_nodes_in_block;
	}	


	//! use a,b,c for extern mesh 
	//! use i,j,k for intern mesh
	for (ind[0]=0; ind[0] < BlkNds_X; ind[0]++)
	 for (ind[1]=0; ind[1] < BlkNds_Y; ind[1]++)
	  for (ind[2]=0; ind[2] < BlkNds_Z; ind[2]++)
	  {
	
	
		i_j_k =  ind[0]*BlkNds_Y*BlkNds_Z
			+ind[1]*BlkNds_Z
			+ind[2];
		
		//! get Coordinate of intern Block Node 
		intern2normedCoords(x_BlockNode, cell_intern_r, ind);
	
		//! find lower next neighbour node in read Mesh
		x_extern[0] = (x_BlockNode[0] +origin_for_hInt[0])/extern_delta[0];
		x_extern[1] = (x_BlockNode[1] +origin_for_hInt[1])/extern_delta[1];
		x_extern[2] = (x_BlockNode[2] +origin_for_hInt[2])/extern_delta[2];
		

		 a  = int(x_extern[0]);
		 b  = int(x_extern[1]);
		 c  = int(x_extern[2]);

		
		r[0] = x_extern[0]-a;
		r[1] = x_extern[1]-b;
		r[2] = x_extern[2]-c;
		
		shape_func[0] = (1.-r[0])*(1.-r[1])*(1.-r[2]);
		shape_func[1] = (   r[0])*(1.-r[1])*(1.-r[2]);
		shape_func[2] = (1.-r[0])*(   r[1])*(1.-r[2]);
		shape_func[3] = (1.-r[0])*(1.-r[1])*(   r[2]);
		
		shape_func[4] = (   r[0])*(   r[1])*(1.-r[2]);
		shape_func[5] = (   r[0])*(1.-r[1])*(   r[2]);
		shape_func[6] = (1.-r[0])*(   r[1])*(   r[2]);
		shape_func[7] = (   r[0])*(   r[1])*(   r[2]);
		
		//! -----------------------------------------------
		a_b_c   =      a*num_uniform_Nodes[1]*num_uniform_Nodes[2]    +b*num_uniform_Nodes[2]      +c;
		
		//! -----------------------------------------------
		ap1_b_c = (a+1)*num_uniform_Nodes[1]*num_uniform_Nodes[2]     +b*num_uniform_Nodes[2]     +c;
		a_bp1_c =     a*num_uniform_Nodes[1]*num_uniform_Nodes[2] +(b+1)*num_uniform_Nodes[2]     +c;
		a_b_cp1 =     a*num_uniform_Nodes[1]*num_uniform_Nodes[2]     +b*num_uniform_Nodes[2] +(c+1);
		
		//! ------------------------------------------------
		ap1_bp1_c = (a+1)*num_uniform_Nodes[1]*num_uniform_Nodes[2] +(b+1)*num_uniform_Nodes[2]     +c;
		ap1_b_cp1 = (a+1)*num_uniform_Nodes[1]*num_uniform_Nodes[2]     +b*num_uniform_Nodes[2] +(c+1);
		a_bp1_cp1 =     a*num_uniform_Nodes[1]*num_uniform_Nodes[2] +(b+1)*num_uniform_Nodes[2] +(c+1);
		
		ap1_bp1_cp1 = (a+1)*num_uniform_Nodes[1]*num_uniform_Nodes[2] +(b+1)*num_uniform_Nodes[2] +(c+1);

		
		//! Set pointer to intern (Block) field
		if(COMPs_FType[src_field_id]==1)
		{
			Hint_dest_field[i_j_k]  =  Hint_array[  a_b_c] * shape_func[0]
					+Hint_array[ap1_b_c] * shape_func[1]
					+Hint_array[a_bp1_c] * shape_func[2]
					+Hint_array[a_b_cp1] * shape_func[3]
			
					+Hint_array[  ap1_bp1_c] * shape_func[4]
					+Hint_array[  ap1_b_cp1] * shape_func[5]
					+Hint_array[  a_bp1_cp1] * shape_func[6]
					+Hint_array[ap1_bp1_cp1] * shape_func[7];	
					
	
			//! NOTE set zero for regions where no integration is performed	
			//! I am not perfectly sure if this works also if integration starts at z=const 
			//! instead of z=0
			for(INT32 axis=0; axis<3; axis++)
			if(integration_axis[axis])
			{
				if( (x_BlockNode[axis]<origin_for_hInt[axis] && origin_for_hInt[axis]<Box_Origin[axis]) || 
					(x_BlockNode[axis]>origin_for_hInt[axis]-Box_Origin[axis] && origin_for_hInt[axis]==LBox_for_hInt[axis]))
				{
					Hint_dest_field[i_j_k] = 0;
				}
			}

		}	

		if(COMPs_FType[src_field_id]==3)
		{
			
			Hint_dest_fieldX[i_j_k]  =  Hint_array[  a_b_c] * shape_func[0]
				+Hint_array[ap1_b_c] * shape_func[1]
				+Hint_array[a_bp1_c] * shape_func[2]
				+Hint_array[a_b_cp1] * shape_func[3]
				+Hint_array[  ap1_bp1_c] * shape_func[4]
				+Hint_array[  ap1_b_cp1] * shape_func[5]
				+Hint_array[  a_bp1_cp1] * shape_func[6]
				+Hint_array[ap1_bp1_cp1] * shape_func[7];	
				
			Hint_dest_fieldY[i_j_k]  =  Hint_array[  a_b_c +num_uniform_box_nodes] * shape_func[0]
				+Hint_array[ap1_b_c + num_uniform_box_nodes] * shape_func[1]
				+Hint_array[a_bp1_c + num_uniform_box_nodes] * shape_func[2]
				+Hint_array[a_b_cp1 + num_uniform_box_nodes] * shape_func[3]
				+Hint_array[  ap1_bp1_c + num_uniform_box_nodes] * shape_func[4]
				+Hint_array[  ap1_b_cp1 + num_uniform_box_nodes] * shape_func[5]
				+Hint_array[  a_bp1_cp1 + num_uniform_box_nodes] * shape_func[6]
				+Hint_array[ap1_bp1_cp1 + num_uniform_box_nodes] * shape_func[7];	
				
			Hint_dest_fieldZ[i_j_k]  =  Hint_array[  a_b_c + 2*num_uniform_box_nodes] * shape_func[0]
				+Hint_array[ap1_b_c + 2*num_uniform_box_nodes] * shape_func[1]
				+Hint_array[a_bp1_c + 2*num_uniform_box_nodes] * shape_func[2]
				+Hint_array[a_b_cp1 + 2*num_uniform_box_nodes] * shape_func[3]
				+Hint_array[  ap1_bp1_c + 2*num_uniform_box_nodes] * shape_func[4]
				+Hint_array[  ap1_b_cp1 + 2*num_uniform_box_nodes] * shape_func[5]
				+Hint_array[  a_bp1_cp1 + 2*num_uniform_box_nodes] * shape_func[6]
				+Hint_array[ap1_bp1_cp1 + 2*num_uniform_box_nodes] * shape_func[7];
				
			//! NOTE set zero for regions where no integration is performed	
			//! I am not perfectly sure if this works also if integration starts at z=const 
			//! instead of z=0
			for(INT32 axis=0; axis<3; axis++)
			if(integration_axis[axis])
			{
				if( (x_BlockNode[axis]<origin_for_hInt[axis] && origin_for_hInt[axis]<Box_Origin[axis]) || 
					(x_BlockNode[axis]>origin_for_hInt[axis]-Box_Origin[axis] && origin_for_hInt[axis]==LBox_for_hInt[axis]))
				{
					Hint_dest_fieldX[i_j_k] = 0;	
					Hint_dest_fieldY[i_j_k] = 0;	
					Hint_dest_fieldZ[i_j_k] = 0;	
				}
			}	
		}

	  }

}