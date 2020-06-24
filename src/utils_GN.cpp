



#include "CHybrid.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"

#include <math.h>





//! define array of function pointers for
//! cp GN procedure
//! (they are initialzed in CHybrid.h)
void (*get_GN_equal_process[3])(CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32    field_type);


//!************************************************************************************************//
//! ******************************** COPY Ghost Nodes *********************************************//
//!************************************************************************************************//

//!---------------------------------------------------------------//
//! get_Im1_GN_equal_process:  					  //
//!---------------------------------------------------------------//
void get_I_GN_equal_process(bool direc,
			    CBlock* dest_Blk,
			    CBlock* src_Blk,
			    INT32    field_type)
{

	//! src refers to neighboring blk, dest to calling blk.
	INT32 XCS_src, XCS_dest;

	//! abbreviations
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];
	D_REAL *src_Field  = src_Blk->Field_Type[field_type];

	//! decide whether to send first or last physical cells
	if(!direc)
	{	
		XCS_src = (BlkNds_X-2) * BlkNds_Y * BlkNds_Z;
		XCS_dest =            0 * BlkNds_Y * BlkNds_Z;

	}
	else
	{
		XCS_src =              1* BlkNds_Y *BlkNds_Z;
		XCS_dest =  (BlkNds_X-1) *BlkNds_Y *BlkNds_Z;
	}


	//! cp Full Cross Section
	//! do this for each component
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	memcpy(dest_Field +comp*num_nodes_in_block + XCS_dest,
	       src_Field  +comp*num_nodes_in_block + XCS_src,
	       BlkNds_Y *BlkNds_Z *sizeof(D_REAL));
}





//!---------------------------------------------------------------//
//! get_Jm1_GN_equal_process:  						  //
//!---------------------------------------------------------------//
void get_J_GN_equal_process(bool direc,
			    CBlock* dest_Blk,
			    CBlock* src_Blk,
			    INT32    field_type)
{


	INT32 k_src_row, k_dest_row;

	//! abbreviations
	INT32 XCS_pos;
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];
	D_REAL *src_Field  =  src_Blk->Field_Type[field_type];



	//! decide whether to send first or last physical cells
	if(!direc)
	{	
		//! Start Pos of penultimate k-row of a given XCS
		k_src_row  =   (BlkNds_Y-2) *BlkNds_Z;
		k_dest_row =             0  *BlkNds_Z;


	}
	else
	{
		//! Start Pos of penultimate k-row of a given XCS
		k_src_row  =             1  *BlkNds_Z;
		k_dest_row =   (BlkNds_Y-1) *BlkNds_Z;

	}



	//! Do not cp Full Cross Section
	//! leaving out k_row's at BlkNds_X=0 & BlkNds_X-2 doesnt matter.
	//! BUT: Full k_row has to be copied else.
	//! this loop copies the k-rows
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X-1; a++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		XCS_pos = a *BlkNds_Y  *BlkNds_Z;
	
		memcpy(dest_Field +comp*num_nodes_in_block +XCS_pos +k_dest_row,
		        src_Field +comp*num_nodes_in_block +XCS_pos +k_src_row,
		        BlkNds_Z *sizeof(D_REAL));
	 }


}



//!---------------------------------------------------------------//
//! cpGC_k_to_kp1_equal_process:  //
//!---------------------------------------------------------------//
void get_K_GN_equal_process(bool direc,
			    CBlock* dest_Blk,
			    CBlock* src_Blk,
			    INT32    field_type)
{

	//! Do not cp Full Cross Section
	//! Leaving out border cells doesnt matter
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];
	D_REAL *src_Field  = src_Blk->Field_Type[field_type];


	INT32 Y_row, XCS_pos, k_src_cell, k_dest_cell;

	//! decide whether to send first or last physical cells
	if(!direc)
	{	
		//! Start Pos of penultimate k-cell
		k_src_cell  =   (BlkNds_Z-2);
		k_dest_cell =             0;

	}
	else
	{
		//! Start Pos of penultimate k-cell
		k_src_cell  =             1;
		k_dest_cell =   (BlkNds_Z-1);

	}


	//! this loop copies the cells
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X-1; a++)
	 {
		XCS_pos = a *BlkNds_Y *BlkNds_Z;

		for(INT32 b=1; b<BlkNds_Y-1; b++)
		{
			Y_row = b *BlkNds_Z;

			dest_Field[comp*num_nodes_in_block +XCS_pos +Y_row + k_dest_cell] =
			 src_Field[comp*num_nodes_in_block +XCS_pos +Y_row + k_src_cell];

		}
	 }


}



//!************************************************************************************************//
//! ******************************** ADD Ghost Nodes **********************************************//
//!************************************************************************************************//

//! array of function pointers for
//! add GN procedure


void (*get_add_GN_equal_process[3])(CBlock* dest_Blk,
				    CBlock* src_Blk,
				    INT32    field_type);


//!---------------------------------------------------------------//
//! add_I_GN_same_process: necessary for moment gathering         //
//!---------------------------------------------------------------//
void add_I_GN_same_process(CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type)
{

	//! abreviations
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];
	D_REAL *src_Field  = src_Blk->Field_Type[field_type];


	INT32 Y_row = 0;
	INT32 XCS_planeSRC = (BlkNds_X-1) * BlkNds_Y * BlkNds_Z;
	INT32 XCS_planeDEST = 1 * BlkNds_Y * BlkNds_Z;

	//!NOTE: be carefull with global indices a,b,c

	//! this loop copies the cells
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	for(INT32 b=1; b<BlkNds_Y; b++)
	{
		Y_row = b *BlkNds_Z;

		for(INT32 c=1; c<BlkNds_Z; c++)
		dest_Field[comp*num_nodes_in_block +XCS_planeDEST +Y_row + c] +=
		 src_Field[comp*num_nodes_in_block +XCS_planeSRC +Y_row + c];

	}

}



//!---------------------------------------------------------------//
//! add_Plane_j_to_jp1: necessary for moment gathering            //
//!---------------------------------------------------------------//
void add_J_GN_same_process(CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type)
{


	//! abreviations
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];
	D_REAL *src_Field  = src_Blk->Field_Type[field_type];


	INT32 X_plane = 0;
	INT32 Y_SRC =  (BlkNds_Y-1) * BlkNds_Z;
	INT32 Y_DEST = 1 * BlkNds_Z;

	//!NOTE: be carefull with global indices a,b,c

	//! this loop copies the cells
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	for(INT32 a=1; a<BlkNds_X; a++)
	{
		X_plane = a *BlkNds_Y *BlkNds_Z;

		for(INT32 c=1; c<BlkNds_Z; c++)
		{

		    dest_Field[comp*num_nodes_in_block +X_plane +Y_DEST + c] 
		   += src_Field[comp*num_nodes_in_block +X_plane +Y_SRC  + c];

		}
	}

}


//!---------------------------------------------------------------//
//! add_K_GN_same_process: necessary for moment gathering            //
//!---------------------------------------------------------------//
void add_K_GN_same_process(CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type)
{

	//! abreviations
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];
	D_REAL *src_Field  = src_Blk->Field_Type[field_type];

	//! src is GC layer
	INT32 Z_SRC =  BlkNds_Z-1;

	//! dest is first active layer
	INT32 Z_DEST = 1;

	INT32 Y_row, XCS_pos;

	//!NOTE: be carefull with global indices a,b,c
	
	//! this loop adds the cells
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X; a++)
	 {

		XCS_pos = a *BlkNds_Y *BlkNds_Z;
		for(INT32 b=1; b<BlkNds_Y; b++)
		{
		   Y_row = b *BlkNds_Z;

		   dest_Field[comp*num_nodes_in_block +XCS_pos + Y_row + Z_DEST]
		  += src_Field[comp*num_nodes_in_block +XCS_pos + Y_row + Z_SRC];
		}
	 }

}












//!************************************************************************************************//
//! ******************************** GET Ghost Nodes from Parent **********************************//
//!************************************************************************************************//
//! define array of function pointers for
//! cp GN procedure
//! (they are initialzed in CHybrid.h)
void (*BV_from_parent[6])(CBlock* dest_Blk,
			  INT32    field_type);



//!-------------------------------------------------------------//
//! im1_BV_from_parent:						//
//! TODO: Evtl nur äußere von Parent nehmen, 			//
//!       innere von Nachbarblock des selben Levels			//
//!-------------------------------------------------------------//
void im1_BV_from_parent(CBlock* child_block, INT32 field_type)
{

	D_REAL* child_field  = child_block->Field_Type[field_type];
	D_REAL* parent_field = child_block->parent->Field_Type[field_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	//! TODO: Check this for odd number of BlkNds_
	INT32 parent_Xstart = child_block->Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = child_block->Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = child_block->Blk_Index[2] * (BlkNds_Z/2-1);

	INT32 a,b,c;
	//! a,b,c have no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block

       //!----------------------------------------------------------------------
       //!------ X=0 Crossections ----------------------------------------------
       //!----------------------------------------------------------------------
       a=0;
//        if(Blk_Index[0]==0)
       for (b = 0; b < BlkNds_Y/2; b++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		INT32 u = 2*a;
		INT32 v = 2*b;
		INT32 w = 2*c;

		INT32 u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
		INT32 u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC

		INT32 u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
		INT32 u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC

		INT32 i = parent_Xstart +a;
		INT32 j = parent_Ystart +b;
		INT32 k = parent_Zstart +c;

		INT32 i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		INT32 ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		INT32 i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		INT32 ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		INT32 i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		INT32 ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		INT32 i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
		{
		 INT32 comp_start = comp*num_nodes_in_block;

		 //!1) cell centered: u_v_w Node
		 child_field[comp_start+u_v_w] = 0.125 * //! k-Plane
					( parent_field[comp_start+    i_j_k]
					 +parent_field[comp_start+  ip1_j_k]
					 +parent_field[comp_start+  i_jp1_k]
					 +parent_field[comp_start+ip1_jp1_k]
					      	//! kp1-Plane
					 +parent_field[comp_start+    i_j_kp1]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);


		 //!2) face centered: u_vp1_w Node jp1-Plane
		 child_field[comp_start+u_vp1_w] = 0.25 *
					( parent_field[comp_start+    i_jp1_k]
					 +parent_field[comp_start+  ip1_jp1_k]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);
	
		 //!3) face centered: u_v_wp1 Node kp1-Plane
		 child_field[comp_start+u_v_wp1] = 0.25 *
					 (parent_field[comp_start+    i_j_kp1]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);


		 //!4) line centered: u_vp1_wp1 Node 
		 child_field[comp_start+u_vp1_wp1] = 0.5 *
					 (parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);
		 

		}

	}




}






//!-------------------------------------------------------------//
//! ip1_BV_from_parent:								//
//!-------------------------------------------------------------//
void ip1_BV_from_parent(CBlock* child_block, INT32 field_type)
{

	D_REAL* child_field  = child_block->Field_Type[field_type];
	D_REAL* parent_field = child_block->parent->Field_Type[field_type];

	//! start from nodes in parent field depends
	//! on the Block indices:
	//! TODO: Check this for odd number of BlkNds_
	INT32 parent_Xstart = child_block->Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = child_block->Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = child_block->Blk_Index[2] * (BlkNds_Z/2-1);


	INT32 a,b,c;
	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block


       //!----------------------------------------------------------------------
       //!------ X=X_MAX Crossections ------------------------------------------
       //!----------------------------------------------------------------------
       a = BlkNds_X/2-1;
//        if(Blk_Index[0]==1)
       for (b = 0; b < BlkNds_Y/2; b++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		INT32 u = 2*a;
		INT32 v = 2*b;
		INT32 w = 2*c;

		INT32 up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC
		INT32 up1_vp1_w = (u+1)*BlkNds_Y*BlkNds_Z   +(v+1)*BlkNds_Z +w; //!LC

		INT32 up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC
		INT32 up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		INT32 i = parent_Xstart +a;
		INT32 j = parent_Ystart +b;
		INT32 k = parent_Zstart +c;


		INT32 ip1_j_k     = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		INT32 ip1_jp1_k   = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
		INT32 ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
		{
		 INT32 comp_start = comp*num_nodes_in_block;


		 //!1) face centered: up1_v_w Node in ip1-Plane
		 child_field[comp_start+up1_v_w] = 0.25 *
					( parent_field[comp_start+    ip1_j_k]
					 +parent_field[comp_start+  ip1_jp1_k]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		  //!2) line centered: up1_vp1_w Node 
		 child_field[comp_start+up1_vp1_w] = 0.5 *
					 (parent_field[comp_start+ip1_jp1_k]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!3) line centered: up1_v_wp1 Node 
		 child_field[comp_start+up1_v_wp1] = 0.5 *
					 (parent_field[comp_start+ip1_j_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!4) identical node: up1_vp1_wp1 Node 
		 child_field[comp_start+up1_vp1_wp1] = parent_field[comp_start+ip1_jp1_kp1];
		 

		}

	}


}

//!-------------------------------------------------------------//
//! jm1_BV_from_parent:						//
//!-------------------------------------------------------------//
void jm1_BV_from_parent(CBlock* child_block, INT32 field_type)
{

	D_REAL* child_field  = child_block->Field_Type[field_type];
	D_REAL* parent_field = child_block->parent->Field_Type[field_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = child_block->Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = child_block->Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = child_block->Blk_Index[2] * (BlkNds_Z/2-1);


	INT32 a,b,c;
	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block


       //!----------------------------------------------------------------------
       //!------ Y=0 Crossections ----------------------------------------------
       //!----------------------------------------------------------------------
       b = 0;
//        if(Blk_Index[1]==0)
       for (a = 0; a < BlkNds_X/2; a++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		INT32 u = 2*a;
		INT32 v = 2*b;
		INT32 w = 2*c;

		INT32 u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
		INT32 up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC


		INT32 u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
		INT32 up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC

		INT32 i = parent_Xstart +a;
		INT32 j = parent_Ystart +b;
		INT32 k = parent_Zstart +c;

		INT32 i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		INT32 ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		INT32 i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		INT32 ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		INT32 i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		INT32 ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		INT32 i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
		{
		 INT32 comp_start = comp*num_nodes_in_block;


		 //!1) cell centered: u_v_w Node
		 child_field[comp_start+u_v_w] = 0.125 * //! k-Plane
					( parent_field[comp_start+    i_j_k]
					 +parent_field[comp_start+  ip1_j_k]
					 +parent_field[comp_start+  i_jp1_k]
					 +parent_field[comp_start+ip1_jp1_k]
					      	//! kp1-Plane
					 +parent_field[comp_start+    i_j_kp1]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!2) face centered: up1_v_w Node in ip1-Plane
		 child_field[comp_start+up1_v_w] =  0.25 *
					( parent_field[comp_start+    ip1_j_k]
					 +parent_field[comp_start+  ip1_jp1_k]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);


		 //!3) face centered: u_v_wp1 Node kp1-Plane
		 child_field[comp_start+u_v_wp1] = 0.25 *
					 (parent_field[comp_start+    i_j_kp1]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);



		 //!4) line centered: up1_v_wp1 Node 
		 child_field[comp_start+up1_v_wp1] = 0.5 *
					 (parent_field[comp_start+ip1_j_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		}

	}



}

//!-------------------------------------------------------------//
//! jp1_BV_from_parent:						//
//!-------------------------------------------------------------//
void jp1_BV_from_parent(CBlock* child_block, INT32 field_type)
{

	D_REAL* child_field  = child_block->Field_Type[field_type];
	D_REAL* parent_field = child_block->parent->Field_Type[field_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	//! TODO: Check this for odd number of BlkNds_
	INT32 parent_Xstart = child_block->Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = child_block->Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = child_block->Blk_Index[2] * (BlkNds_Z/2-1);


	INT32 a,b,c;
	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block


       //!----------------------------------------------------------------------
       //!------ Y=Y_MAX Crossections ------------------------------------------
       //!----------------------------------------------------------------------
       b=BlkNds_Y/2-1;
//        if(Blk_Index[1]==1)
       for (a = 0; a < BlkNds_X/2; a++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		INT32 u = 2*a;
		INT32 v = 2*b;
		INT32 w = 2*c;

		INT32 u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC
		INT32 up1_vp1_w = (u+1)*BlkNds_Y*BlkNds_Z   +(v+1)*BlkNds_Z +w; //!LC

		INT32 u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC
		INT32 up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		INT32 i = parent_Xstart +a;
		INT32 j = parent_Ystart +b;
		INT32 k = parent_Zstart +c;


      		INT32   i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		INT32 ip1_jp1_k   = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		INT32   i_jp1_kp1 =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
		{
		 INT32 comp_start = comp*num_nodes_in_block;




		 //!1) face centered: u_vp1_w Node jp1-Plane
		 child_field[comp_start+u_vp1_w] = 0.25 *
					( parent_field[comp_start+    i_jp1_k]
					 +parent_field[comp_start+  ip1_jp1_k]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);
	
		  //!2) line centered: up1_vp1_w Node 
		 child_field[comp_start+up1_vp1_w] = 0.5 *
					 (parent_field[comp_start+ip1_jp1_k]
					 +parent_field[comp_start+ip1_jp1_kp1]);


		 //!3) line centered: u_vp1_wp1 Node 
		 child_field[comp_start+u_vp1_wp1] = 0.5 *
					 (parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!4) identical node: up1_vp1_wp1 Node 
		 child_field[comp_start+up1_vp1_wp1] = parent_field[comp_start+ip1_jp1_kp1];
		 

		}

	}



}

//!-------------------------------------------------------------//
//! km1_BV_from_parent:						//
//!-------------------------------------------------------------//
void km1_BV_from_parent(CBlock* child_block, INT32 field_type)
{


	D_REAL* child_field  = child_block->Field_Type[field_type];
	D_REAL* parent_field = child_block->parent->Field_Type[field_type];

	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = child_block->Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = child_block->Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = child_block->Blk_Index[2] * (BlkNds_Z/2-1);


	INT32 a,b,c;
	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block


       //!----------------------------------------------------------------------
       //!------ Z=0 Crossections ----------------------------------------------
       //!----------------------------------------------------------------------
       c=0;
//        if(Blk_Index[2]==0)
       for (a = 0; a < BlkNds_X/2; a++)
        for (b = 0; b < BlkNds_Y/2; b++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		INT32 u = 2*a;
		INT32 v = 2*b;
		INT32 w = 2*c;

		INT32 u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
		INT32 up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC
		INT32 u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC
		INT32 up1_vp1_w =  (u+1)*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!LC

		INT32 i = parent_Xstart +a;
		INT32 j = parent_Ystart +b;
		INT32 k = parent_Zstart +c;

		INT32 i_j_k       =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		INT32 ip1_j_k     = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		INT32 i_jp1_k     =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		INT32 ip1_jp1_k   = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		INT32 i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		INT32 ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		INT32 i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
		{
			INT32 comp_start = comp*num_nodes_in_block;
	
			//!1) cell centered: u_v_w Node
			child_field[comp_start+u_v_w] = 0.125 * //! k-Plane
						( parent_field[comp_start+    i_j_k]
						+parent_field[comp_start+  ip1_j_k]
						+parent_field[comp_start+  i_jp1_k]
						+parent_field[comp_start+ip1_jp1_k]
							//! kp1-Plane
						+parent_field[comp_start+    i_j_kp1]
						+parent_field[comp_start+  ip1_j_kp1]
						+parent_field[comp_start+  i_jp1_kp1]
						+parent_field[comp_start+ip1_jp1_kp1]);
	
			//!2) face centered: up1_v_w Node in ip1-Plane
			child_field[comp_start+up1_v_w] = 0.25 *
						( parent_field[comp_start+    ip1_j_k]
						+parent_field[comp_start+  ip1_jp1_k]
						+parent_field[comp_start+  ip1_j_kp1]
						+parent_field[comp_start+ip1_jp1_kp1]);
	
			//!3) face centered: u_vp1_w Node jp1-Plane
			child_field[comp_start+u_vp1_w] = 0.25 *
						( parent_field[comp_start+    i_jp1_k]
						+parent_field[comp_start+  ip1_jp1_k]
						+parent_field[comp_start+  i_jp1_kp1]
						+parent_field[comp_start+ip1_jp1_kp1]);
		
			//!4) line centered: up1_vp1_w Node 
			child_field[comp_start+up1_vp1_w] = 0.5 *
						(parent_field[comp_start+ip1_jp1_k]
						+parent_field[comp_start+ip1_jp1_kp1]);

		 

		}

	}

}


//!-------------------------------------------------------------//
//! kp1_BV_from_parent:						//
//!-------------------------------------------------------------//
void kp1_BV_from_parent(CBlock* child_block, INT32 field_type)
{



	D_REAL* child_field  = child_block->Field_Type[field_type];
	D_REAL* parent_field = child_block->parent->Field_Type[field_type];

	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = child_block->Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = child_block->Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = child_block->Blk_Index[2] * (BlkNds_Z/2-1);


	INT32 a,b,c;
	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block


       //!----------------------------------------------------------------------
       //!------ Z=MAX_Z Crossections ------------------------------------------
       //!----------------------------------------------------------------------
       c=BlkNds_Z/2-1;

       for (a = 0; a < BlkNds_X/2; a++)
        for (b = 0; b < BlkNds_Y/2; b++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		INT32 u = 2*a;
		INT32 v = 2*b;
		INT32 w = 2*c;

		INT32 u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
		INT32 up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC
		INT32 u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC
		INT32 up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		INT32 i = parent_Xstart +a;
		INT32 j = parent_Ystart +b;
		INT32 k = parent_Zstart +c;


		INT32 i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		INT32 ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		INT32 i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		INT32 ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
		{
		 INT32 comp_start = comp*num_nodes_in_block;



	
		 //!1) face centered: u_v_wp1 Node kp1-Plane
		 child_field[comp_start+u_v_wp1] = 0.25 *
					 (parent_field[comp_start+    i_j_kp1]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);


		 //!2) line centered: up1_v_wp1 Node 
		 child_field[comp_start+up1_v_wp1] = 0.5 *
					 (parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!3) line centered: u_vp1_wp1 Node 
		 child_field[comp_start+u_vp1_wp1] = 0.5 *
					 (parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!4) identical node: up1_vp1_wp1 Node 
		 child_field[comp_start+up1_vp1_wp1] =  parent_field[comp_start+ip1_jp1_kp1];
		 

		}

	}

}


