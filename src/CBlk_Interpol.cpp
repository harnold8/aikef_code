



#include "CBlk.h"
#include "parameters.h"


#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include "absolute_Globals.h"

extern D_REAL ***smooth_coeff;
extern bool* is_inflow_species;


using namespace std;


//! --- Some general Comments: ---------------------------------------
//! 1) While ghost cells of eg. jm1 lay in between to parent "rows"s,
//!	 the jp1 ghost cells lay directly on a parent "row". So the 
//!	 copying is not symmetrically


//!------------------------------------------------------------
//!---------- CBlock: Global Interpolation Solver Variables ---
//!------------------------------------------------------------
extern INT32 *COMPs_FType;


INT32 u,v,w, u_v_w, up1_v_w, u_vp1_w, u_v_wp1;
INT32 up1_vp1_w, up1_v_wp1, u_vp1_wp1, up1_vp1_wp1;

INT32 comp_start, temp_RLevel, target_RLevel;
// INT32 parent_Xstart, parent_Ystart, parent_Zstart;



CBlock* temp_parent;

extern INT32 *temp_i, *temp_j, *temp_k;

extern INT32 a, b, c;
extern INT32 i,j,k, ip1_jp1_kp1;
extern INT32 i_j_k, ip1_j_k, im1_j_k, i_jp1_k , i_jm1_k ,i_j_kp1 ,i_j_km1;
extern INT32 ip1_j_kp1, ip1_jp1_k, i_jp1_kp1;




//!-------------------------------------------------------------//
//! field_from_parent:						//
//!-------------------------------------------------------------//
void CBlock::field_from_parent(INT32 parent_type, INT32 child_type, INT32* start_ijk, INT32* end_ijk)
{
	INT32 type = child_type;
	
	//! of course num_comps of parent_type and child_type must equal
	D_REAL* child_field  = Field_Type[child_type];
	D_REAL* parent_field = parent->Field_Type[parent_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = Blk_Index[2] * (BlkNds_Z/2-1);



	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block
	for (INT32 a = start_ijk[0]; a < end_ijk[0]; a++)
         for (INT32 b = start_ijk[1]; b < end_ijk[1]; b++)
          for (INT32 c = start_ijk[2]; c < end_ijk[2]; c++)
	  {

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		u = 2*a;
		v = 2*b;
		w = 2*c;

		u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
		up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC
		u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC
		up1_vp1_w = (u+1)*BlkNds_Y*BlkNds_Z   +(v+1)*BlkNds_Z +w; //!LC

		u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
		up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC
		u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC
		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;



		for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;

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
	
		 //!4) face centered: u_v_wp1 Node kp1-Plane
		 child_field[comp_start+u_v_wp1] = 0.25 *	
					 (parent_field[comp_start+    i_j_kp1]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		  //!5) line centered: up1_vp1_w Node 
		 child_field[comp_start+up1_vp1_w] = 0.5 *	
					 (parent_field[comp_start+ip1_jp1_k]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!6) line centered: up1_v_wp1 Node 
		 child_field[comp_start+up1_v_wp1] = 0.5 *	
					 (parent_field[comp_start+ip1_j_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!7) line centered: u_vp1_wp1 Node 
		 child_field[comp_start+u_vp1_wp1] = 0.5 *	
					 (parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);
		 //!8) identical node: up1_vp1_wp1 Node 
		 child_field[comp_start+up1_vp1_wp1] =
					 parent_field[comp_start+ip1_jp1_kp1];
		 
		}
	   }


}



/*
//!-------------------------------------------------------------//
//! BV_from_parent_EvenNP:					//
//! TODO: Evtl nur äußere von Parent nehmen, 			//
//!       innere von Nachbarblock des selben Levels		//
//!-------------------------------------------------------------//
void CBlock::BV_from_parent(INT32 parent_type, INT32 child_type)
{
	INT32 type = child_type;
	//! of course num_comps of parent_type and child_type must equal
	D_REAL* child_field  = Field_Type[child_type];
	D_REAL* parent_field = parent->Field_Type[parent_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	//! TODO: Check this for odd number of BlkNds_
	parent_Xstart = Blk_Index[0] * (BlkNds_X/2-1);
	parent_Ystart = Blk_Index[1] * (BlkNds_Y/2-1);
	parent_Zstart = Blk_Index[2] * (BlkNds_Z/2-1);


	INT32 a,b,c;
	//! a,b,c no meaning (arbitrary counter index)
	//! i,j,k meaning Indices of Parent
	//! u,v,w meaning Indices of Block

       //!----------------------------------------------------------------------
       //!------ X=0 Crossections ----------------------------------------------
       //!----------------------------------------------------------------------
       a=0;
       if(Blk_Index[0]==0)
       for (b = 0; b < BlkNds_Y/2; b++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		u = 2*a;
		v = 2*b;
		w = 2*c;

		u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
		u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC

		u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
		u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;

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


       //!----------------------------------------------------------------------
       //!------ X=X_MAX Crossections ------------------------------------------
       //!----------------------------------------------------------------------
       a = BlkNds_X/2-1;
       if(Blk_Index[0]==1)
       for (b = 0; b < BlkNds_Y/2; b++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		u = 2*a;
		v = 2*b;
		w = 2*c;

		up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC
		up1_vp1_w = (u+1)*BlkNds_Y*BlkNds_Z   +(v+1)*BlkNds_Z +w; //!LC

		up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC
		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;


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
		 child_field[comp_start+up1_vp1_wp1] =
					 parent_field[comp_start+ip1_jp1_kp1];
		 

		}

	}

       //!----------------------------------------------------------------------
       //!------ Y=0 Crossections ----------------------------------------------
       //!----------------------------------------------------------------------
       b = 0;
       if(Blk_Index[1]==0)
       for (a = 0; a < BlkNds_X/2; a++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		u = 2*a;
		v = 2*b;
		w = 2*c;

		u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
		up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC


		u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
		up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;

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


       //!----------------------------------------------------------------------
       //!------ Y=Y_MAX Crossections ------------------------------------------
       //!----------------------------------------------------------------------
       b=BlkNds_Y/2-1;
       if(Blk_Index[1]==1)
       for (a = 0; a < BlkNds_X/2; a++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		u = 2*a;
		v = 2*b;
		w = 2*c;

		u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC
		up1_vp1_w = (u+1)*BlkNds_Y*BlkNds_Z   +(v+1)*BlkNds_Z +w; //!LC

		u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC
		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;




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
		 child_field[comp_start+up1_vp1_wp1] =
					 parent_field[comp_start+ip1_jp1_kp1];
		 

		}

	}

       //!----------------------------------------------------------------------
       //!------ Z=0 Crossections ----------------------------------------------
       //!----------------------------------------------------------------------
       c=0;
       if(Blk_Index[2]==0)
       for (a = 0; a < BlkNds_X/2; a++)
        for (b = 0; b < BlkNds_Y/2; b++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		u = 2*a;
		v = 2*b;
		w = 2*c;

		u_v_w     =      u*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!CC
		up1_v_w   =  (u+1)*BlkNds_Y*BlkNds_Z      +v*BlkNds_Z +w; //!FC
		u_vp1_w   =      u*BlkNds_Y*BlkNds_Z  +(v+1)*BlkNds_Z +w; //!FC
		up1_vp1_w = (u+1)*BlkNds_Y*BlkNds_Z   +(v+1)*BlkNds_Z +w; //!LC

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;

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

       //!----------------------------------------------------------------------
       //!------ Z=MAX_Z Crossections ------------------------------------------
       //!----------------------------------------------------------------------
       c=BlkNds_Y/2-1;
       if(Blk_Index[2]==1)
       for (a = 0; a < BlkNds_X/2; a++)
        for (b = 0; b < BlkNds_Y/2; b++)
	{

		//! Note:
		//! use u,v,w as loop parameters and write
		//! i = parent_Xstart +u/2 instead results in errors
		u = 2*a;
		v = 2*b;
		w = 2*c;

		u_v_wp1 =         u*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!FC
		up1_v_wp1   = (u+1)*BlkNds_Y*BlkNds_Z     +v*BlkNds_Z +w+1;//!LC
		u_vp1_wp1   =     u*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!LC
		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		i_j_k     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
		ip1_j_k   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k;
      		i_jp1_k   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;
      		ip1_jp1_k = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k;

		i_j_kp1     =     i*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
		ip1_j_kp1   = (i+1)*BlkNds_Y*BlkNds_Z     +j*BlkNds_Z +k+1;
      		i_jp1_kp1   =     i*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
      		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		for(INT32 comp=0; comp<COMPs_FType[type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;



	
		 //!1) face centered: u_v_wp1 Node kp1-Plane
		 child_field[comp_start+u_v_wp1] = 0.25 *	
					 (parent_field[comp_start+    i_j_kp1]
					 +parent_field[comp_start+  ip1_j_kp1]
					 +parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);


		 //!2) line centered: up1_v_wp1 Node 
		 child_field[comp_start+up1_v_wp1] = 0.5 *	
					 (parent_field[comp_start+ip1_j_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);

		 //!3) line centered: u_vp1_wp1 Node 
		 child_field[comp_start+u_vp1_wp1] = 0.5 *	
					 (parent_field[comp_start+  i_jp1_kp1]
					 +parent_field[comp_start+ip1_jp1_kp1]);
		 //!4) identical node: up1_vp1_wp1 Node 
		 child_field[comp_start+up1_vp1_wp1] =
					 parent_field[comp_start+ip1_jp1_kp1];
		 

		}

	}


}

*/

//!-------------------------------------------------------------//
//! cp_child_to_parent:							//
//! inverse function of field_FROM_parent node #8			//
//!-------------------------------------------------------------//
void CBlock::add_field_to_parent(INT32 parent_type, INT32 child_type)
{

	D_REAL* child_field  = Field_Type[child_type];
	D_REAL* parent_field = parent->Field_Type[parent_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = Blk_Index[2] * (BlkNds_Z/2-1);

	//! NOTE: !!!
	//! NOTE: !!!
	//! NOTE: !!!
	//! Parent indices ranging from:
	//! [       1 <= BlkNds/2] including BlkNds/2 !!!
	//! [BlkNds/2 < BlkNds     ]
	//! -> region of two distinct octs overlaps !!!
	//! NOTE: !!!
	//! NOTE: !!!
	//! NOTE: !!!

	//! e.g.: BlkNds_X = 6:

	//! Blk_Index=0
	//! -> parent_Xstart=0
	//! -> a=0,1,2
	//! -> v=0,2,4
	//! -> i=0,1,2
	//! -> ip1=1,2,3
	//! -> vp1=1,3,5


	//! Blk_Index=1
	//! -> parent_Xstart=2
	//! -> a=0,1,2
	//! -> v=0,2,4
	//! -> i=2,3,4
	//! -> ip1=3,4,5
	//! -> vp1=1,3,5


       INT32 a,b,c;
      //! former version: Ghost knodes are not considered
      //! -> loop up to BlkNds_X/2-1
      for (a = 0; a < BlkNds_X/2; a++)
       for (b = 0; b < BlkNds_Y/2; b++)
        for (c = 0; c < BlkNds_Z/2; c++)
	{

		u = 2*a;
		v = 2*b;
		w = 2*c;

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		


		for(INT32 comp=0; comp<COMPs_FType[child_type]; comp++)
		{
			comp_start = comp*num_nodes_in_block;
			parent_field[comp_start+ip1_jp1_kp1] += child_field[comp_start+up1_vp1_wp1];
		}

	}


}


//!-------------------------------------------------------------//
//! cp_child_to_parent:								//
//! inverse function of field_FROM_parent				//
//!-------------------------------------------------------------//
void CBlock::field_to_parent_smooth(INT32 parent_type, INT32 child_type)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	INT32 a,b,c;
	INT32 child_elmt;

	D_REAL* child_field  = Field_Type[child_type];
	D_REAL* parent_field = parent->Field_Type[parent_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = Blk_Index[2] * (BlkNds_Z/2-1);

	//! former version: Ghost knodes are not considered
	//! -> loop up to BlkNds_X/2-1
	for (a = 0; a < BlkNds_X/2-1; a++)
	 for (b = 0; b < BlkNds_Y/2-1; b++)
	  for (c = 0; c < BlkNds_Z/2-1; c++)
	  {

		//! e.g.: BlkNds_X = 6:

		//! Blk_Index=0
		//! -> parent_Xstart=0
		//! -> a=0,1
		//! -> v=0,2
		//! -> i=0,1
		//! -> ip1=1,2
		//! -> vp1=1,3


		//! Blk_Index=1
		//! -> parent_Xstart=2
		//! -> a=0,1
		//! -> v=0,2
		//! -> i=2,3
		//! -> ip1=3,4
		//! -> vp1=1,3

		u = 2*a;
		v = 2*b;
		w = 2*c;

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z
			     +(j+1)*BlkNds_Z
			     + k+1;


	      //! using a averaging procedure to copy child values to parents
	      //! Even though less numerical artefacts occur, the solution 
	      //! sometimes looks a bit "unphysical, mesh influenced".
	      //! Especially upstream disturbence is generated.
	      //! Maybe at the end the best way is to use the copy procedure 
	      //! and use a std smoothing procedure.
	      for(INT32 comp=0; comp<COMPs_FType[child_type]; comp++)
	      {
	        comp_start = comp*num_nodes_in_block;

		//! Set parent field to zero for smooth child to parent
		parent_field[comp_start+ip1_jp1_kp1] = 0;


		for (short e=-1; e<=1; e++)
	 	 for (short f=-1; f<=1; f++)
	  	  for (short g=-1; g<=1; g++)
		  {
		   child_elmt = comp*num_nodes_in_block +(u+e+1)*BlkNds_Y*BlkNds_Z 
							+(v+f+1)*BlkNds_Z
							+(w+g+1);
          	   parent_field[comp_start+ip1_jp1_kp1] 
				+= smooth_coeff[e+1][f+1][g+1] * child_field[child_elmt];
		  }
	      }
	}


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}

//!-------------------------------------------------------------//
//! cp_child_to_parent:							//
//! inverse function of field_FROM_parent node #8			//
//!-------------------------------------------------------------//
void CBlock::moments_to_parent_smooth(INT32 field_type)
{

	//! record processing time of this block
	clock_t time_start, time_finish;
	time_start = clock();


	INT32 a,b,c;
	INT32 child_elmt;

	D_REAL* child_field  = Field_Type[field_type];
	D_REAL* parent_field = parent->Field_Type[field_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	INT32 parent_Xstart = Blk_Index[0] * (BlkNds_X/2-1);
	INT32 parent_Ystart = Blk_Index[1] * (BlkNds_Y/2-1);
	INT32 parent_Zstart = Blk_Index[2] * (BlkNds_Z/2-1);

	INT32 gstart[3];
	if(gather_Neighbour[_Im1_])
	gstart[0] = 0;
	else
	gstart[0] = 1;

	if(gather_Neighbour[_Jm1_])
	gstart[1] = 0;
	else
	gstart[1] = 1;

	if(gather_Neighbour[_Km1_])
	gstart[2] = 0;
	else
	gstart[2] = 1;
	

      //! former version: Ghost knodes are not considered
      //! -> loop up to BlkNds_X/2-1
      for (a = gstart[0]; a < BlkNds_X/2-1; a++)
       for (b = gstart[1]; b < BlkNds_Y/2-1; b++)
        for (c = gstart[2]; c < BlkNds_Z/2-1; c++)
	{

		//! e.g.: BlkNds_X = 6, Blk_Index=1
		//! -> parent_Xstart=2
		//! -> a=0,1
		//! -> v=0,2
		//! -> i=2,3
		//! -> ip1=3,4
		//! -> vp1=1,3

		u = 2*a;
		v = 2*b;
		w = 2*c;

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z
			     +(j+1)*BlkNds_Z
			     + k+1;


		
		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z
		             +(v+1)*BlkNds_Z
		             + w+1;
// 
// 
// 
// 		for(INT32 comp=0; comp<COMPs_FType[child_type]; comp++)
// 		{
// 		 comp_start = comp*num_nodes_in_block;
// 		 parent_field[comp_start+ip1_jp1_kp1] = child_field[comp_start+up1_vp1_wp1];
// 		}


	      //! using a averaging procedure to copy child values to parents
	      //! Even though less numerical artefacts occur, the solution 
	      //! sometimes looks a bit "unphysical, mesh influenced".
	      //! Especially upstream disturbence is generated.
	      //! Maybe at the end the best way is to use the copy procedure 
	      //! and use a std smoothing procedure.
	      for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	      {
	        comp_start = comp*num_nodes_in_block;

		//! Set parent field to zero for smooth child to parent
		parent_field[comp_start+ip1_jp1_kp1] = 0;


		for (short e=-1; e<=1; e++)
	 	 for (short f=-1; f<=1; f++)
	  	  for (short g=-1; g<=1; g++)
		  {
		   child_elmt = comp*num_nodes_in_block +(u+e+1)*BlkNds_Y*BlkNds_Z 
							+(v+f+1)*BlkNds_Z
							+(w+g+1);
          	   parent_field[comp_start+ip1_jp1_kp1] 
				+= smooth_coeff[e+1][f+1][g+1] * child_field[child_elmt];
		  }
	      }
	}


	//! record field calculation time
	time_finish = clock();
	time_process_fields += (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;

}
/*
//!-------------------------------------------------------------//
//! cp_child_to_parent:							//
//! inverse function of field_FROM_parent node #8			//
//!-------------------------------------------------------------//
void CBlock::add_inject_GC_to_parent(INT32 parent_type, INT32 child_type)
{
	INT32 a,b,c;
	short e,f,g;
	INT32 child_elmt;

	child_field  = Field_Type[child_type];
	parent_field = parent->Field_Type[parent_type];


	//! start from nodes in parent field depends
	//! on the Block indices:
	parent_Xstart = Blk_Index[0] * (BlkNds_X/2-1);
	parent_Ystart = Blk_Index[1] * (BlkNds_Y/2-1);
	parent_Zstart = Blk_Index[2] * (BlkNds_Z/2-1);




       //!----------------------------------------------------------------------
       //!------ X=BlkNds_X/2-1; Crossections -----------------
       //!----------------------------------------------------------------------


	//! -0 and -1 in loop ends symbolies the "loop end offset":
       //! Offset in ijk-loops is necessary to avoid some 
       //! nodes (diagonal & Block edge) are added to their
       //! parent GC several times.

       a = BlkNds_X/2-1;
       for (b = 0; b < BlkNds_Y/2-0; b++)
	for (c = 0; c < BlkNds_Z/2-0; c++)
	{


		u = 2*a;
		v = 2*b;
		w = 2*c;

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		//! e.g.: BlkNds_X = 6, Blk_Index=1
		//! -> parent_Xstart=2
		//! -> a=2
		//! -> u=4
		//! -> i=4
		//! -> ip1=5
		//! -> up1=5
		

		//! INJECTION VERSION
		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;
		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;

		for(INT32 comp=0; comp<COMPs_FType[child_type]; comp++)
		{
		comp_start = comp*num_nodes_in_block;
		parent_field[comp_start+ip1_jp1_kp1] += (  child_field[comp_start+up1_vp1_wp1]);
		}


		//! SMOOTHED VERSION DOES NOT WORK AT LEVEL BOUNDARIES,
		//! REALLY DIFFICULT, ONLY IMLEMENT IN CASE SERIOUS
		//! PROBLEMS OCCUR DUE TO INJECTION METHOD
		
// 	      for(INT32 comp=0; comp<COMPs_FType[child_type]; comp++)
// 	      {
// 	        comp_start = comp*num_nodes_in_block;
// 		//! Set parent field to zero for smooth child to parent
// // 		parent_field[comp_start+ip1_jp1_kp1] = 0;
// 		for (short e=-1; e<=0; e++)
// 	 	 for (short f=-1; f<=1; f++)
// 	  	  for (short g=-1; g<=1; g++)
// 		  {
// 		   child_elmt = comp*num_nodes_in_block +(u+e+1)*BlkNds_Y*BlkNds_Z 
// 							+(v+f+1)*BlkNds_Z +(w+g+1);
// 
// // 		   if( ((u+e+1)*BlkNds_Y*BlkNds_Z +(v+f+1)*BlkNds_Z +(w+g+1)) < num_nodes_in_block)
//           	   parent_field[comp_start+ip1_jp1_kp1] 
// 				+= smooth_coeff[e+1][f+1][g+1] * child_field[child_elmt];
// 		  }
// 	      }


	}

       //!----------------------------------------------------------------------
       //!------ Y=BlkNds_Y/2-1; Crossections -----------------
       //!----------------------------------------------------------------------
       b = BlkNds_Y/2-1;
       for (a = 0; a < BlkNds_X/2-1; a++)
        for (c = 0; c < BlkNds_Z/2-0; c++)
	{

		u = 2*a;
		v = 2*b;
		w = 2*c;

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		for(INT32 comp=0; comp<COMPs_FType[child_type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;
		 parent_field[comp_start+ip1_jp1_kp1] += child_field[comp_start+up1_vp1_wp1];
		}

	}

	//!----------------------------------------------------------------------
       //!------ Z=BlkNds_Z/2-1; Crossections -----------------
       //!----------------------------------------------------------------------
       c = BlkNds_Z/2-1;
       for (a = 0; a < BlkNds_X/2-1; a++)
        for (b = 0; b < BlkNds_Y/2-1; b++)
	{

		u = 2*a;
		v = 2*b;
		w = 2*c;

		i = parent_Xstart +a;
		j = parent_Ystart +b;
		k = parent_Zstart +c;

		ip1_jp1_kp1 = (i+1)*BlkNds_Y*BlkNds_Z +(j+1)*BlkNds_Z +k+1;


		up1_vp1_wp1 = (u+1)*BlkNds_Y*BlkNds_Z +(v+1)*BlkNds_Z +w+1;//!IN

		for(INT32 comp=0; comp<COMPs_FType[child_type]; comp++)
		{
		 comp_start = comp*num_nodes_in_block;
		 parent_field[comp_start+ip1_jp1_kp1] += child_field[comp_start+up1_vp1_wp1];
		}

	}


}
*/
//!-------------------------------------------------------------//
//! find_im1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_im1_Neighbour(void)
{


	INT32 temp_oct_id=0;
	
	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array
	//! (neighbour may exit or may not)
	if(Blk_Index[0]==1)
	Neighbour[_Im1_] = parent->child_array[ 0*2*2 +Blk_Index[1]*2+Blk_Index[2] ];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[0]==0)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[0]==0 && temp_parent->RLevel>0)
		{

			//! store index of y,z dimension 
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;

		 }

		temp_parent = temp_parent->Neighbour[_Im1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{


			temp_oct_id = +1*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +temp_k[temp_parent->RLevel+1];

			if(!temp_parent->child_array[temp_oct_id])
			{
				//! No im1 neighbour with the same Level
				//! of refinement does exist
				Neighbour[_Im1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step will be stored)
			Neighbour[_Im1_] = temp_parent;


		 }
		 else Neighbour[_Im1_] = 0;

	}

}

//!-------------------------------------------------------------//
//! find_ip1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_ip1_Neighbour(void)
{
	

	INT32 temp_oct_id=0;


	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array
	if(Blk_Index[0]==0)
	Neighbour[_Ip1_] = parent->child_array[ 1*2*2 +Blk_Index[1]*2+Blk_Index[2] ];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[0]==1)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[0]==1 && temp_parent->RLevel>0)
		{

			//! store index of y,z dimension 
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->Neighbour[_Ip1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{

			temp_oct_id = +0*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +temp_k[temp_parent->RLevel+1];


			if(!temp_parent->child_array[temp_oct_id])
			{
				//! No ip1 neighbour with the same Level
				//! of refinement does exist
				Neighbour[_Ip1_] = 0;
				break;
			}



			//! climb to next Level 
			temp_parent=temp_parent->child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			Neighbour[_Ip1_] = temp_parent;


		 }
		 else Neighbour[_Ip1_] = 0;

	}



}



//!-------------------------------------------------------------//
//! find_jm1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_jm1_Neighbour(void)
{
	
	

	INT32 temp_oct_id=0;

	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array
	if(Blk_Index[1]==1)
	Neighbour[_Jm1_] = parent->child_array[ Blk_Index[0]*2*2 +0*2 +Blk_Index[2] ];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[1]==0)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[1]==0 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->Neighbour[_Jm1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{


			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +1*2
				      +temp_k[temp_parent->RLevel+1];


			if(!temp_parent->child_array[temp_oct_id])
			{
				//! No jm1 neighbour with the same Level
				//! of refinement does exist
				Neighbour[_Jm1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			Neighbour[_Jm1_] = temp_parent;


		 }
		 else Neighbour[_Jm1_] = 0;

	}



}

//!-------------------------------------------------------------//
//! find_jp1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_jp1_Neighbour(void)
{
	

	INT32 temp_oct_id=0;

	
	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array
	if(Blk_Index[1]==0)
	Neighbour[_Jp1_] = parent->child_array[ Blk_Index[0]*2*2 +1*2 +Blk_Index[2] ];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[1]==1)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[1]==1 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->Neighbour[_Jp1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{

			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +0*2
				      +temp_k[temp_parent->RLevel+1];


			if(!temp_parent->child_array[temp_oct_id])
			{
				//! No jp1 neighbour with the same Level
				//! of refinement does exist
				Neighbour[_Jp1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			Neighbour[_Jp1_] = temp_parent;


		 }
		 else Neighbour[_Jp1_] = 0;

	}



}



//!-------------------------------------------------------------//
//! find_km1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_km1_Neighbour(void)
{
	
	INT32 temp_oct_id=0;


	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[2]==1)
	Neighbour[_Km1_] = parent->child_array[ Blk_Index[0]*2*2 +Blk_Index[1]*2 +0];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[2]==0)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[2]==0 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->Neighbour[_Km1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{

			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +1;


			if(!temp_parent->child_array[temp_oct_id])
			{
				//! No km1 neighbour with the same Level
				//! of refinement does exist
				Neighbour[_Km1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			Neighbour[_Km1_] = temp_parent;


		 }
		 else Neighbour[_Km1_] = 0;

	}

}

//!-------------------------------------------------------------//
//! find_kp1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_kp1_Neighbour(void)
{
	

	INT32 temp_oct_id=0;


	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[2]==0)
	Neighbour[_Kp1_] = parent->child_array[ Blk_Index[0]*2*2 +Blk_Index[1]*2+1];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[2]==1)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[2]==1 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->Neighbour[_Kp1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{

			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +0;


			if(!temp_parent->child_array[temp_oct_id])
			{
				//! No kp1 neighbour with the same Level
				//! of refinement does exist
				Neighbour[_Kp1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			Neighbour[_Kp1_] = temp_parent;


		 }
		 else Neighbour[_Kp1_] = 0;

	}



}

//!-------------------------------------------------------------//
//! find_im1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_gather_im1_Neighbour(void)
{
	
	INT32 temp_oct_id=0;

	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[0]==1)
	gather_Neighbour[_Im1_] = parent->gather_child_array[ 0*2*2 +Blk_Index[1]*2+Blk_Index[2]];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[0]==0)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[0]==0 && temp_parent->RLevel>0)
		{

			//! store index of y,z dimension 
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;

		 }

		temp_parent = temp_parent->gather_Neighbour[_Im1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{

			temp_oct_id = +1*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +temp_k[temp_parent->RLevel+1];


			if(!temp_parent->gather_child_array[temp_oct_id])
			{
				//! No im1 neighbour with the same Level
				//! of refinement does exist
				gather_Neighbour[_Im1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->gather_child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			gather_Neighbour[_Im1_] = temp_parent;


		 }
		 else gather_Neighbour[_Im1_] = 0;

	}


}

//!-------------------------------------------------------------//
//! find_ip1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_gather_ip1_Neighbour(void)
{
	

	INT32 temp_oct_id=0;

	
	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[0]==0)
	gather_Neighbour[_Ip1_] = parent->gather_child_array[ 1*2*2 +Blk_Index[1]*2+Blk_Index[2]];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[0]==1)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[0]==1 && temp_parent->RLevel>0)
		{

			//! store index of y,z dimension 
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->gather_Neighbour[_Ip1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{


			temp_oct_id = +0*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +temp_k[temp_parent->RLevel+1];

			if(!temp_parent->gather_child_array[temp_oct_id])
			{
				//! No ip1 neighbour with the same Level
				//! of refinement does exist
				gather_Neighbour[_Ip1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->gather_child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			gather_Neighbour[_Ip1_] = temp_parent;


		 }
		 else gather_Neighbour[_Ip1_] = 0;

	}



}



//!-------------------------------------------------------------//
//! find_jm1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_gather_jm1_Neighbour(void)
{
	
	
	INT32 temp_oct_id=0;


	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[1]==1)
	gather_Neighbour[_Jm1_] = parent->gather_child_array[ Blk_Index[0]*2*2 +0*2 +Blk_Index[2]];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[1]==0)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[1]==0 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->gather_Neighbour[_Jm1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{

			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +1*2
				      +temp_k[temp_parent->RLevel+1];


			if(!temp_parent->gather_child_array[temp_oct_id])
			{
				//! No jm1 neighbour with the same Level
				//! of refinement does exist
				gather_Neighbour[_Jm1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->gather_child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			gather_Neighbour[_Jm1_] = temp_parent;


		 }
		 else gather_Neighbour[_Jm1_] = 0;

	}



}

//!-------------------------------------------------------------//
//! find_jp1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_gather_jp1_Neighbour(void)
{

	
	INT32 temp_oct_id=0;

	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[1]==0)
	gather_Neighbour[_Jp1_] = parent->gather_child_array[ Blk_Index[0]*2*2 +1*2 +Blk_Index[2]];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[1]==1)
	{




		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[1]==1 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_k[temp_parent->RLevel] = temp_parent->Blk_Index[2];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->gather_Neighbour[_Jp1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{


			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +0*2
				      +temp_k[temp_parent->RLevel+1];
	
			if(!temp_parent->gather_child_array[temp_oct_id])
			{
				//! No jp1 neighbour with the same Level
				//! of refinement does exist
				gather_Neighbour[_Jp1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->gather_child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			gather_Neighbour[_Jp1_] = temp_parent;


		 }
		 else gather_Neighbour[_Jp1_] = 0;

	}



}



//!-------------------------------------------------------------//
//! find_km1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_gather_km1_Neighbour(void)
{
	

	INT32 temp_oct_id=0;


	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[2]==1)
	gather_Neighbour[_Km1_] = parent->gather_child_array[ Blk_Index[0]*2*2 +Blk_Index[1]*2 +0];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[2]==0)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[2]==0 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->gather_Neighbour[_Km1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{

			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +1;

			if(!temp_parent->gather_child_array[temp_oct_id])
			{
				//! No km1 neighbour with the same Level
				//! of refinement does exist
				gather_Neighbour[_Km1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->gather_child_array[temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			gather_Neighbour[_Km1_] = temp_parent;


		 }
		 else gather_Neighbour[_Km1_] = 0;

	}

}

//!-------------------------------------------------------------//
//! find_kp1_Neighbour:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_gather_kp1_Neighbour(void)
{


	INT32 temp_oct_id=0;
	
	
	temp_parent = this;
	target_RLevel = RLevel;


	//! neighbour is within same Block Array and does exist
	if(Blk_Index[2]==0)
	gather_Neighbour[_Kp1_] = parent->gather_child_array[ Blk_Index[0]*2*2 +Blk_Index[1]*2+1];

	//! neighbour is NOT within the same parent-array
	//! and may not exist
	if(Blk_Index[2]==1)
	{

		//! descent to corporate parent array,
		//! this may also be the root array
		while(temp_parent->Blk_Index[2]==1 && temp_parent->RLevel>0)
		{

			//! store index of x,z dimension 
			temp_i[temp_parent->RLevel] = temp_parent->Blk_Index[0];
			temp_j[temp_parent->RLevel] = temp_parent->Blk_Index[1];

			temp_parent=temp_parent->parent;


		 }

		temp_parent = temp_parent->gather_Neighbour[_Kp1_];

		//! now climb back to target_RLevel if possible
		if(temp_parent)
		while(temp_parent->RLevel<target_RLevel)
		{


			temp_oct_id = +temp_i[temp_parent->RLevel+1]*2*2 
				      +temp_j[temp_parent->RLevel+1]*2
				      +0;


			if(!temp_parent->gather_child_array[temp_oct_id])
			{
				//! No kp1 neighbour with the same Level
				//! of refinement does exist
				gather_Neighbour[_Kp1_] = 0;
				break;
			}

			//! climb to next Level 
			temp_parent=temp_parent->gather_child_array [temp_oct_id];

			//! always set to active Block 
			//! (only temp_parent of last loop step wiil be stored)
			gather_Neighbour[_Kp1_] = temp_parent;


		 }
		 else gather_Neighbour[_Kp1_] = 0;

	}



}


//!-------------------------------------------------------------//
//! find_all_ordanary_Neighbours:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_all_ordanary_Neighbours(void)
{

	find_im1_Neighbour();
	find_ip1_Neighbour();

	find_jm1_Neighbour();
	find_jp1_Neighbour();

	find_km1_Neighbour();
	find_kp1_Neighbour();

}

//!-------------------------------------------------------------//
//! find_all_ordanary_Neighbours:							//
//! Only to use when Block is no root block			//
//!-------------------------------------------------------------//
void CBlock::find_all_gather_Neighbours(void)
{

	find_gather_im1_Neighbour();
	find_gather_ip1_Neighbour();

	find_gather_jm1_Neighbour();
	find_gather_jp1_Neighbour();

	find_gather_km1_Neighbour();
	find_gather_kp1_Neighbour();

}

//!-------------------------------------------------------------//
//! link_Neighbours_pts_to_this:					  //
//! Only to use when Block is no root block				  //
//!-------------------------------------------------------------//
void CBlock::link_all_ordanary_Neighbours(void)
{

	if(Neighbour[_Im1_]) Neighbour[_Im1_]->Neighbour[_Ip1_] = this;
	if(Neighbour[_Ip1_]) Neighbour[_Ip1_]->Neighbour[_Im1_] = this;

	if(Neighbour[_Jm1_]) Neighbour[_Jm1_]->Neighbour[_Jp1_] = this;
	if(Neighbour[_Jp1_]) Neighbour[_Jp1_]->Neighbour[_Jm1_] = this;

	if(Neighbour[_Km1_]) Neighbour[_Km1_]->Neighbour[_Kp1_] = this;
	if(Neighbour[_Kp1_]) Neighbour[_Kp1_]->Neighbour[_Km1_] = this;


}




//!-------------------------------------------------------------//
//! link_Neighbours_pts_to_this:					  //
//! Only to use when Block is no root block				  //
//!-------------------------------------------------------------//
void CBlock::link_all_gather_Neighbours(void)
{

	if(gather_Neighbour[_Im1_]) gather_Neighbour[_Im1_]->gather_Neighbour[_Ip1_] = this;
	if(gather_Neighbour[_Ip1_]) gather_Neighbour[_Ip1_]->gather_Neighbour[_Im1_] = this;

	if(gather_Neighbour[_Jm1_]) gather_Neighbour[_Jm1_]->gather_Neighbour[_Jp1_] = this;
	if(gather_Neighbour[_Jp1_]) gather_Neighbour[_Jp1_]->gather_Neighbour[_Jm1_] = this;

	if(gather_Neighbour[_Km1_]) gather_Neighbour[_Km1_]->gather_Neighbour[_Kp1_] = this;
	if(gather_Neighbour[_Kp1_]) gather_Neighbour[_Kp1_]->gather_Neighbour[_Km1_] = this;


}

//!-------------------------------------------------------------//
//! unlink_all_Neighbours:						//
//! Only to use when Block is no root block				//
//!-------------------------------------------------------------//
void CBlock::unlink_all_Neighbours(void)
{

	if(Neighbour[_Im1_]) Neighbour[_Im1_]->Neighbour[_Ip1_] = 0;
	if(Neighbour[_Ip1_]) Neighbour[_Ip1_]->Neighbour[_Im1_] = 0;

	if(Neighbour[_Jm1_]) Neighbour[_Jm1_]->Neighbour[_Jp1_] = 0;
	if(Neighbour[_Jp1_]) Neighbour[_Jp1_]->Neighbour[_Jm1_] = 0;

	if(Neighbour[_Km1_]) Neighbour[_Km1_]->Neighbour[_Kp1_] = 0;
	if(Neighbour[_Kp1_]) Neighbour[_Kp1_]->Neighbour[_Km1_] = 0;

	if(gather_Neighbour[_Im1_]) gather_Neighbour[_Im1_]->gather_Neighbour[_Ip1_] = 0;
	if(gather_Neighbour[_Ip1_]) gather_Neighbour[_Ip1_]->gather_Neighbour[_Im1_] = 0;

	if(gather_Neighbour[_Jm1_]) gather_Neighbour[_Jm1_]->gather_Neighbour[_Jp1_] = 0;
	if(gather_Neighbour[_Jp1_]) gather_Neighbour[_Jp1_]->gather_Neighbour[_Jm1_] = 0;

	if(gather_Neighbour[_Km1_]) gather_Neighbour[_Km1_]->gather_Neighbour[_Kp1_] = 0;
	if(gather_Neighbour[_Kp1_]) gather_Neighbour[_Kp1_]->gather_Neighbour[_Km1_] = 0;

}


