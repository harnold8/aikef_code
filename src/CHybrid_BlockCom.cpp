

//! iostream for memcpy
#include "CHybrid.h"
#include "defines.h"
#include "parameters.h"
#include "utils.h"
#include "absolute_Globals.h"

#include <sstream>

#define ADD_GN		true
#define DO_NOT_ADD_GN	false
#define NO_minDIREC_SEND	-1

#define send	true
#define recv	false

#define compose_field	true
#define    copy_field	false




extern void (*BV_from_parent[6])(CBlock* dest_Blk,
				  INT32    field_type);

//!---------------------------------------------------------------//
//! cp_GN_equal_process: 					  //
//!---------------------------------------------------------------//
void CHybrid::cp_GN_equal_process(INT32 level, INT32 _m1_, INT32 _p1_, INT32 field_type)
{


	//! decide wether to loop through ordinary or gather blocklist
	if(_p1_ < gathNEIB)
	{

		//! ORDINARY LIST
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->cp_GN_equal_process(_m1_, _p1_, field_type);
		

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}

	}
	else
	{

		//! GATHER LIST
		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->cp_GN_equal_process(_m1_, _p1_, field_type);
		

			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}

	}


}



//!---------------------------------------------------------------//
//! cp_GN_send_MPI: 						  //
//!---------------------------------------------------------------//
void CHybrid::cp_GN_send_MPI(INT32 level, INT32 _m1_, INT32 _p1_, INT32 field_type)
{



	//! decide wether to loop through ordinary or gather blocklist
	if(_p1_ < gathNEIB)
	{

		//! ORDINARY LIST
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->cp_GN_send_MPI(_m1_, _p1_, field_type);
		

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}

	}
	else
	{

		//! GATHER LIST
		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->cp_GN_send_MPI(_m1_, _p1_, field_type);
		

			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}

	}

}

//!---------------------------------------------------------------//
//! cp_GN_recv_MPI: 						  //
//!---------------------------------------------------------------//
void CHybrid::cp_GN_recv_MPI(INT32 level, INT32 _m1_, INT32 _p1_, INT32 field_type)
{



	//! decide wether to loop through ordinary or gather blocklist
	if(_p1_ < gathNEIB)
	{

		//! ORDINARY LIST
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->cp_GN_recv_MPI(_m1_, _p1_, field_type);
		

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}

	}
	else
	{

		//! GATHER LIST
		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->cp_GN_recv_MPI(_m1_, _p1_, field_type);
		

			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}

	}

}


//!---------------------------------------------------------------//
//! wait_all_MPI: 				  //
//!---------------------------------------------------------------//
void CHybrid::wait_all_MPI(INT32 level)
{

	if(level <0)
	return;

	//! GATHER LIST
	CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
	while(temp_Block)
	{
	

// 		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->wait_all_MPI();
	

		temp_Block = temp_Block->next_Blk_of_GatherBlockList;

	}


}


//!---------------------------------------------------------------//
//! wait_SendFields_MPI: 				  //
//!---------------------------------------------------------------//
void CHybrid::wait_SendFields_MPI(INT32 level, INT32 _m1_, INT32 _p1_)
{


	//! decide wether to loop through ordinary or gather blocklist
	if(_p1_ < gathNEIB)
	{

		//! ORDINARY LIST
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->wait_SendFields_MPI(_m1_, _p1_);
		

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}

	}
	else
	{

		//! GATHER LIST
		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		while(temp_Block)
		{
		
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->wait_SendFields_MPI(_m1_, _p1_);
		

			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}

	}

}



//!------------------------------------------------------------------------------------//
//! get_GN_from_parent: 											//
//!------------------------------------------------------------------------------------//
void CHybrid::get_GN_from_parent(INT32 level,
				 INT32 _m1_,
				 INT32 _p1_,
				 INT32 field_type,
				 INT32* end_ijk,
				 bool incl_first_phys_face)
{


	INT32 null_vec[3] = {0,0,0};

	CBlock* temp_Block;

	//! get GN at level borders
	//! (not parallelized yet)
	temp_Block = BlockList_of_Lev[level];
	if(level)
	while(temp_Block)
	{
	
		//! Boundary Blocks don't receive values from neighbour/parent
		//! minus direction (might include first physical row when gathering)
		if(mpi_myRank == temp_Block->responsible_mpi_process)
		{
			if(!temp_Block->is_box_boundary[_m1_] && !temp_Block->Neighbour[_m1_])
			{
				if(incl_first_phys_face)
				temp_Block->field_from_parent(field_type, field_type, null_vec, end_ijk);
				else
				BV_from_parent[_m1_](temp_Block, field_type);
			}

			//! plus direction
			if(!temp_Block->is_box_boundary[_p1_] && !temp_Block->Neighbour[_p1_])
			BV_from_parent[_p1_](temp_Block, field_type);
		}
	
	
		temp_Block = temp_Block->next_Blk_of_BlockList;
	
	}



}




//!-------------------------------------------------------------//
//! zero_parent_add_children_field_MPI: 					//
//! send injected field to parent via MPI				//
//!-------------------------------------------------------------//
void CHybrid::zero_parent_add_children_field_MPI(INT32 level, INT32 id_field)
{


	child_send_recv_field_MPI( level,    id_field, send);
	parent_send_recv_field_MPI(level-1, id_field, recv, compose_field);
	if(sync_mpi_send_rec) {
	  wait_child_send_recv_field_MPI( level,    id_field);
	}
	else {
	  wait_all_MPI(level);
	}

}

//!-------------------------------------------------------------//
//! zero_parent_add_children_field_MPI: 					//
//! send injected field to parent via MPI				//
//!-------------------------------------------------------------//
void CHybrid::copy_children_field_to_parent_MPI(INT32 level, INT32 id_field)
{

	child_send_recv_field_MPI( level,    id_field, send);
	parent_send_recv_field_MPI(level-1, id_field, recv, copy_field);
	if(sync_mpi_send_rec) {
	  wait_child_send_recv_field_MPI( level,    id_field);
	}
	else {
	  wait_all_MPI(level);
	}

}


//!-------------------------------------------------------------//
//! field_to_children_MPI: 						//
//!-------------------------------------------------------------//
void CHybrid::field_to_children_MPI(INT32 level, INT32 id_field)
{

	//! get field from parents to allow interpolation
	//! at block boundaries
	//! NOTE: WAIT HAS TO BE CALLED IN PARENT LEVEL
	parent_send_recv_field_MPI(level-1, id_field, send, compose_field);
	child_send_recv_field_MPI( level  , id_field, recv);
	if(sync_mpi_send_rec) {
	  wait_parent_send_recv_field_MPI(level-1, id_field);
	}
	else {
	  wait_all_MPI(level-1);
	}


}



//!-------------------------------------------------------------//
//! parent_send_recv_field_MPI: 						//
//! - This method handles blocks whose					//
//!   children ARE NOT on this process					//
//!   -> receive up to eight children fields each block		//
//!-------------------------------------------------------------//
void CHybrid::parent_send_recv_field_MPI(INT32 level, INT32 id_field, bool do_send, bool compose)
{


	if(level <0) return;



	//! make sure not to send/recv field of parent more than once.
	bool* process_done = new bool[mpi_num_processes];


	//! there can be 8 distinct processors each oct
	//! (at maximum)
	D_REAL *buffer_of_proc[8];


	INT32* process_buffer_id =  new INT32[mpi_num_processes];
	
	//! only allocate when recv is called
	if(!do_send)
	 for(INT32 proc=0; proc<8; proc++)
	 buffer_of_proc[proc] = new D_REAL[COMPs_FType[id_field] *num_nodes_in_block];


	//! loop through parent blocks
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{


		//! reset all processes
		memset(process_done, 0, mpi_num_processes*sizeof(bool));


		//! reset process counter each block
		INT32 buffer_id = 0;

		//! check for every child:
		//! -> has different responsible process?
		//! -> if yes, sent parent field
		//! -> make sure that field is only sent once to
		//!    every process, since all 8 octs share
		//!    the same parent field
		if(mpi_myRank==temp_Block->responsible_mpi_process)
		 for(INT32 oct=0; oct<8; oct++)
		  if(temp_Block->child_array[oct])
		  {

			//! get child process
			INT32 child_process = temp_Block->child_array[oct]->responsible_mpi_process;


			//! Make sure to send / recv field only once each process
			if(child_process!=mpi_myRank && !process_done[child_process])
			{
	


				//! SEND TO MAXIMAL 8 CHILDREN PROCESSES
				if(do_send)
				{
					//! use oct nr. as request id
					//! maximal 8 request may be used
					//! in case each children is
					//! on another process
					INT32 id_request = oct;

					//! use parent's mpi_tag
					send_MPI_Field(id_field,
						       id_request,
						       child_process,
						       temp_Block->mpi_tag,
						       temp_Block);

					process_done[child_process] = true;

				}
				//! RECEIVE FROM MAXIMAL 8 CHILDREN PROCESSES
				else
				{

// 					log_file << " receiving: " << temp_Block->mpi_tag   +9*total_num_mpi_tags << endl;

					//! - one field may arrive from each child
					//!   they all have to be buffered and results
					//!   composed together in the final parent field
					//! - use parent's mpi_tag
					//! - recv to buffer
					//! - recv field of FULL BLOCK SIZE, even though
					//!   oct is smaller, BUT:
					//!   IN CASE ALL CHILDREN ARE AT THE SAME PROCESS,
					//!   THEY ARE SEND AT ONCE WHICH REQUIRES THE FULL
					//!   BLOCK SIZE.
					recv_MPI_Field(id_field,
						       child_process,
						       temp_Block->mpi_tag,
						       buffer_of_proc[buffer_id]);



					process_buffer_id[child_process] = buffer_id;


					//! count number of mpi processes
					//! from all eight octs
					buffer_id++;


					//! mark process as completed
					process_done[child_process] = true;


				}//! end decide send/recv
			}//! end if oct !my_rank
		 }//! if has child


// 		cout << "num_distinct_mpi_procs: " << num_distinct_mpi_procs << endl;

		//! in case field has been received from at least one other process, copy
		//! received fields of respective process to parent field memory
		if(buffer_id)
		{
			if(compose)
			temp_Block->compose_recvd_fields(id_field,
							 buffer_id,
							 buffer_of_proc);
			else
			temp_Block->copy_recvd_fields(id_field,
						      process_buffer_id,
						      buffer_of_proc);

		}
		

		temp_Block = temp_Block->next_Blk_of_BlockList;

	}



	delete[] process_done;
	delete[] process_buffer_id;

	if(!do_send)
	for(INT32 proc=0; proc<8; proc++)
	 delete[] buffer_of_proc[proc];


}

//!-------------------------------------------------------------//
//! wait_parent_send_recv_field_MPI: 				//
//! - This method handles blocks whose				//
//!   children ARE NOT on this process				//
//!   -> receive up to eight children fields each block		//
//! - Complete Send Recive Process				//
//!-------------------------------------------------------------//
void CHybrid::wait_parent_send_recv_field_MPI(INT32 level, INT32 id_field)
{


	if(level <0) return;



	//! make sure not to send/recv field of parent more than once.
	bool* process_done = new bool[mpi_num_processes];


	//! there can be 8 distinct processors each oct
	//! (at maximum)
	D_REAL *buffer_of_proc[8];

	INT32* process_buffer_id =  new INT32[mpi_num_processes];

	//! loop through parent blocks
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{


		//! reset all processes
		memset(process_done, 0, mpi_num_processes*sizeof(bool));


		//! reset process counter each block
		INT32 buffer_id = 0;

		//! check for every child:
		//! -> has different responsible process?
		//! -> if yes, sent parent field
		//! -> make sure that field is only sent once to
		//!    every process, since all 8 octs share
		//!    the same parent field
		if(mpi_myRank==temp_Block->responsible_mpi_process)
		 for(INT32 oct=0; oct<8; oct++)
		  if(temp_Block->child_array[oct])
		  {

			//! get child process
			INT32 child_process = temp_Block->child_array[oct]->responsible_mpi_process;


			//! Make sure to send / recv field only once each process
			if(child_process!=mpi_myRank && !process_done[child_process])
			{
	


				//! use oct nr. as request id
				//! maximal 8 request may be used
				//! in case each children is
				//! on another process
				INT32 id_request = oct;

				temp_Block->req_is_MPI_send_completed[id_request].Wait();

				process_done[child_process] = true;

			}//! end if oct !my_rank
		 }//! if has child



		temp_Block = temp_Block->next_Blk_of_BlockList;

	}

	delete[] process_done;
	delete[] process_buffer_id;

}


//!-------------------------------------------------------//
//! child_send_recv_field_MPI: 					//
//! - This method handles blocks whose				//
//!   parent IS NOT on this process				//
//!-------------------------------------------------------//
void CHybrid::child_send_recv_field_MPI(INT32 level, INT32 id_field, bool do_send)
{


	//! no parent available in level 0
	if(!level) return;


	//! make sure not to send field of one block
	//! to another process twice.
	bool already_processed = false;


	//! LOOP THROUGH PARENT LEVEL RATHER THAN
	//! CHILDREN LIST
	CBlock* temp_Block = BlockList_of_Lev[level-1];
	while(temp_Block)
	{
	
		already_processed = false;

		//! check for every parent:
		//! -> has parent different responsible process?
		//! -> if yes, recv from parent field
		//! -> make sure that field is only recv once from
		//!    every process, since all 8 octs share
		//!    the same parent field
		if(mpi_myRank!=temp_Block->responsible_mpi_process)
		 for(INT32 oct=0; oct<8; oct++)
		  if(temp_Block->child_array[oct] && !already_processed)
		  {


			//! set abbreviations to respective processes
			INT32 parent_process = temp_Block->responsible_mpi_process;
			INT32  child_process = temp_Block->child_array[oct]->responsible_mpi_process;
			
			if(child_process==mpi_myRank)
			{

				//! mark this block as finished
				already_processed = true;

				//! SEND TO EXACTLY ONCE
				if(do_send)
				{

					//! use oct nr. as request id
					//! maximal 8 request may be used
					//! in case each children is
					//! on another process
					INT32 id_request = oct;

// 					log_file << " sending to p" << 

					//! - use parent's mpi_tag
					//! - do not send childrens's field but
					//!   field that has been smoothed to parent already
					//!   -> use "temp_Block" rather than  temp_Block->child_array[oct]
					send_MPI_Field(id_field,
						       id_request,
						       parent_process,
						       temp_Block->mpi_tag,
						       temp_Block);
					

				}
				//! RECEIVE FROM EXACTLY ONE PARENT PROCESS
				else
				{
					//! - use parent's mpi_tag
					//! - recv parents field into parent memory
					//!   -> use "temp_Block->Field_Type[id_field]"
					//!      rather than  "temp_Block->child_array[oct]->Field_Type[id_field]"
					recv_MPI_Field(id_field,
						       parent_process,
						       temp_Block->mpi_tag,
						       temp_Block->Field_Type[id_field]);

// 					cout << " recv " << id_field << endl;
	

				}
			}

		  }

		temp_Block = temp_Block->next_Blk_of_BlockList;

	}

// 	if(mpi_myRank)
// 	cout  << endl<< endl<< endl;

}

//!-------------------------------------------------------//
//! wait_child_send_recv_field_MPI: 			//
//! - This method handles blocks whose			//
//!   parent IS NOT on this process			//
//! - Wait to finish send receive process		//
//!-------------------------------------------------------//
void CHybrid::wait_child_send_recv_field_MPI(INT32 level, INT32 id_field)
{


	//! no parent available in level 0
	if(!level) return;


	//! make sure not to send field of one block
	//! to another process twice.
	bool already_processed = false;


	//! LOOP THROUGH PARENT LEVEL RATHER THAN
	//! CHILDREN LIST
	CBlock* temp_Block = BlockList_of_Lev[level-1];
	while(temp_Block)
	{
	
		already_processed = false;

		//! check for every parent:
		//! -> has parent different responsible process?
		//! -> if yes, recv from parent field
		//! -> make sure that field is only recv once from
		//!    every process, since all 8 octs share
		//!    the same parent field
		if(mpi_myRank!=temp_Block->responsible_mpi_process)
		 for(INT32 oct=0; oct<8; oct++)
		  if(temp_Block->child_array[oct] && !already_processed)
		  {


			//! set abbreviations to respective processes
			INT32 parent_process = temp_Block->responsible_mpi_process;
			INT32  child_process = temp_Block->child_array[oct]->responsible_mpi_process;
			
			if(child_process==mpi_myRank)
			{

				//! mark this block as finished
				already_processed = true;

				//! use oct nr. as request id
				//! maximal 8 request may be used
				//! in case each children is
				//! on another process
				INT32 id_request = oct;

				temp_Block->req_is_MPI_send_completed[id_request].Wait();

			}

		  }

		temp_Block = temp_Block->next_Blk_of_BlockList;

	}

}

//!---------------------------------------------------------------//
//! cp_rootArray_GC: 						  //
//!---------------------------------------------------------------//
void CHybrid::FULL_GN_UPDATE(INT32 field)
{
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	update_GN_of_level(level, field, SKIP_PHYS_FACE, DO_PARENT_BUFFER_UPDATE);
}


//!---------------------------------------------------------------//
//! update_GN_of_level: 					  //
//!---------------------------------------------------------------//
void CHybrid::update_GN_of_level(INT32 level,
				 INT32 id_field,
				 bool incl_first_phys_face,
				 bool update_parent_buffer)
{

	
	//! for loop
	CBlock* temp_Block;
	
	
	INT32 null_vec[3]={0,0,0};
	
	//! NOTE: inversed order of copying boundary CS 
	//!		1)k 
	//!		2)j 
	//!		3)i
	
	//! indices for field_from_parent function
	INT32 Imin_indices[3]={         1, BlkNds_Y/2, BlkNds_Z/2};
	INT32 Jmin_indices[3]={BlkNds_X/2,          1, BlkNds_Z/2};
	INT32 Kmin_indices[3]={BlkNds_X/2, BlkNds_Y/2,          1};



	//! get field from parents to allow interpolation
	//! at block boundaries
	if(update_parent_buffer)
	field_to_children_MPI(level, id_field);

	

	//!----------------------------------------------------------------
	//! 1) get cells in k direction
	//!----------------------------------------------------------------
	 cp_GN_send_MPI     (level, _Km1_, _Kp1_, id_field);
	 cp_GN_equal_process(level, _Km1_, _Kp1_, id_field);
	get_GN_from_parent  (level, _Km1_, _Kp1_, id_field, Kmin_indices, incl_first_phys_face);
	 cp_GN_recv_MPI     (level, _Km1_, _Kp1_, id_field);
	if(sync_mpi_send_rec) {
	  wait_SendFields_MPI(level, _Km1_, _Kp1_);
	}
	else {
	  wait_all_MPI(level);
	}
	
	


	//!----------------------------------------------------------------
	//! 2) get cells in j direction
	//!   (all i-direction GN must be received by now)
	//!----------------------------------------------------------------
	 cp_GN_send_MPI     (level, _Jm1_, _Jp1_, id_field);
	 cp_GN_equal_process(level, _Jm1_, _Jp1_, id_field);
	get_GN_from_parent  (level, _Jm1_, _Jp1_, id_field, Jmin_indices, incl_first_phys_face);
	 cp_GN_recv_MPI     (level, _Jm1_, _Jp1_, id_field);
	if(sync_mpi_send_rec) {
	  wait_SendFields_MPI(level, _Jm1_, _Jp1_);
	}
	else {
	  wait_all_MPI(level);
	}
	
	
	


	//!----------------------------------------------------------------
	//! 3) get GN in i direction
	//!   (all j-direction GN must be received by now)
	//!----------------------------------------------------------------
	 cp_GN_send_MPI     (level, _Im1_, _Ip1_, id_field);
	 cp_GN_equal_process(level, _Im1_, _Ip1_, id_field);
	get_GN_from_parent  (level, _Im1_, _Ip1_, id_field, Imin_indices, incl_first_phys_face);
	 cp_GN_recv_MPI     (level, _Im1_, _Ip1_, id_field);
	if(sync_mpi_send_rec) {
	  wait_SendFields_MPI (level, _Im1_, _Ip1_);
	}
	else {
	  wait_all_MPI(level);
	}
	
	

}


//!---------------------------------------------------------------//
//! update_GN_of_level: 					  //
//!---------------------------------------------------------------//
void CHybrid::update_GATHER_GC_of_level(INT32 level, INT32 field_type)
{


	//! NOTE: inversed order of copying boundary CS 
	//!		1)k 
	//!		2)j 
	//!		3)i
	
	//!----------------------------------------------------------------
	//! 1) get cells in k direction
	//!----------------------------------------------------------------
	cp_GN_send_MPI		        (level, _Km1_ +gathNEIB, _Kp1_ +gathNEIB, field_type);
	cp_GN_equal_process	        (level, _Km1_ +gathNEIB, _Kp1_ +gathNEIB, field_type);
	cp_GN_recv_MPI		        (level, _Km1_ +gathNEIB, _Kp1_ +gathNEIB, field_type);
	if(sync_mpi_send_rec) {
	  wait_SendFields_MPI(level, _Km1_ +gathNEIB, _Kp1_ +gathNEIB);
	}
	else {
	  wait_all_MPI(level);
	}
// 	

	
	
	//!----------------------------------------------------------------
	//! 2) get cells in j direction
	//!----------------------------------------------------------------
	cp_GN_send_MPI		        (level, _Jm1_ +gathNEIB, _Jp1_ +gathNEIB, field_type);
	cp_GN_equal_process	        (level, _Jm1_ +gathNEIB, _Jp1_ +gathNEIB, field_type);
	cp_GN_recv_MPI		        (level, _Jm1_ +gathNEIB, _Jp1_ +gathNEIB, field_type);
	if(sync_mpi_send_rec) {
	  wait_SendFields_MPI(level, _Jm1_ +gathNEIB, _Jp1_ +gathNEIB);
	}
	else {
	  wait_all_MPI(level);
	}

	
	//!----------------------------------------------------------------
	//! 3) get cells in i direction
	//!----------------------------------------------------------------
	cp_GN_send_MPI		        (level, _Im1_ +gathNEIB, _Ip1_ +gathNEIB, field_type);
	cp_GN_equal_process	        (level, _Im1_ +gathNEIB, _Ip1_ +gathNEIB, field_type);
	cp_GN_recv_MPI		        (level, _Im1_ +gathNEIB, _Ip1_ +gathNEIB, field_type);
	if(sync_mpi_send_rec) {
	  wait_SendFields_MPI(level, _Im1_ +gathNEIB, _Ip1_ +gathNEIB);
	}
	else {
	  wait_all_MPI(level);
	}



}

//!-------------------------------------------------------------//
//! update_GN_of_level:								//
//!-------------------------------------------------------------//
void CHybrid::moments_to_parent_smooth(INT32 field_type)
{

	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{

		//! TODO:
		//! Optimize this
		//! introduce bool "smooth_gathBlk_to_parent" and set after refinement
		//! (diagonale? 6->8->10?)

		update_GATHER_GC_of_level(level, field_type);

		CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
		if(level)
		while(temp_Block)
		{
			//! smooth field to parent
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			if(!temp_Block->is_gatherBlk)
			 temp_Block->moments_to_parent_smooth(field_type);
			 else if(   (temp_Block->Neighbour[6] && !temp_Block->Neighbour[6]->is_gatherBlk)
			          || (temp_Block->Neighbour[8] && !temp_Block->Neighbour[8]->is_gatherBlk)
			          || (temp_Block->Neighbour[10] && !temp_Block->Neighbour[10]->is_gatherBlk))
			 temp_Block->moments_to_parent_smooth(field_type);
			 else
			  if(temp_Block->Neighbour[6] && temp_Block->Neighbour[6]->Neighbour[8] && !temp_Block->Neighbour[6]->Neighbour[8]->is_gatherBlk)
			  temp_Block->moments_to_parent_smooth(field_type);
			   else if(temp_Block->Neighbour[6] && temp_Block->Neighbour[6]->Neighbour[10] && !temp_Block->Neighbour[6]->Neighbour[10]->is_gatherBlk)
			   temp_Block->moments_to_parent_smooth(field_type);
			    else if(temp_Block->Neighbour[8] && temp_Block->Neighbour[8]->Neighbour[10] && !temp_Block->Neighbour[8]->Neighbour[10]->is_gatherBlk)
			    temp_Block->moments_to_parent_smooth(field_type);
			    else if(    temp_Block->Neighbour[6] && temp_Block->Neighbour[6]->Neighbour[8]
			              &&  temp_Block->Neighbour[6]->Neighbour[8]->Neighbour[10]
			              && !temp_Block->Neighbour[6]->Neighbour[8]->Neighbour[10]->is_gatherBlk)
			    temp_Block->moments_to_parent_smooth(field_type);


			temp_Block = temp_Block->next_Blk_of_GatherBlockList;

		}
	}
	
}





//!---------------------------------------------------------------//
//! cp_rootArray_GC: 							  //
//!---------------------------------------------------------------//
void CHybrid::FULL_MOMENT_UPDATE(INT32 field)
{


	//! this must not be gather (otw. GC in level > 0 incorrect)
	//! TODO:
	//! ALSO UPDATE ULTIMATE PHYSICAL NODE LINE AT LEVEL BORDERS
	//! -> NO GATHER BLOCKS SHOULD BE REQUIRED


	//! adding the first physical line is only required, in
	//! case gather blocks are NOT used
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	update_GN_of_level(level, field, !use_gather_blocks, DO_PARENT_BUFFER_UPDATE);

	//! optinonally apply field_to_parent_smooth to Moments
	if(MAX_LEVEL>0 && smooth_moments_to_parent)
	 for(INT32 level=MAX_LEVEL; level>0; level--)
	 {

		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			//! smooth field to parent
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->field_to_parent_smooth(field, field);


			temp_Block = temp_Block->next_Blk_of_BlockList;

		}

		//! update GN of parent
		update_GN_of_level(level-1, field, false, DO_PARENT_BUFFER_UPDATE);
	 }

	

}


//!---------------------------------------------------------------//
//! add_GN_equal_process: 					  //
//!---------------------------------------------------------------//
void CHybrid::add_GN_equal_process(INT32 level, INT32 _m1_, INT32 field_type)
{

	CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->add_GN_equal_process(_m1_, field_type);

		temp_Block = temp_Block->next_Blk_of_GatherBlockList;

	}

}


//!---------------------------------------------------------------//
//! add_GN_send_MPI: 						  //
//!---------------------------------------------------------------//
void CHybrid::add_GN_send_MPI(INT32 level, INT32 _p1_, INT32 field_type)
{

	CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
	while(temp_Block)
	{

		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->add_GN_send_MPI( _p1_, field_type);

		temp_Block = temp_Block->next_Blk_of_GatherBlockList;

	}

}


//!---------------------------------------------------------------//
//! add_GN_recv_MPI: 								  //
//!---------------------------------------------------------------//
void CHybrid::add_GN_recv_MPI(INT32 level, INT32 _m1_, INT32 field_type)
{


	CBlock* temp_Block = GATHER_BlockList_of_Lev[level];
	while(temp_Block)
	{


		if(mpi_myRank == temp_Block->responsible_mpi_process)
		temp_Block->add_GN_recv_MPI( _m1_, field_type);

		temp_Block = temp_Block->next_Blk_of_GatherBlockList;

	}

}



//!---------------------------------------------------------------//
//! Add_Boundary_Moments: 						  //
//!---------------------------------------------------------------//
void CHybrid::Add_Boundary_Moments_of_L(INT32 level, INT32 field_type)
{

   //!----------------------------------------------------------------
   //! 1) get cells in k direction
   //!----------------------------------------------------------------
   add_GN_send_MPI     (level, _Kp1_, field_type);
   add_GN_equal_process(level, _Km1_, field_type);
   add_GN_recv_MPI     (level, _Km1_, field_type);
   if(sync_mpi_send_rec) {
    wait_SendFields_MPI(level, NO_minDIREC_SEND, _Kp1_);
    synchronize_allProcesses_MPI();
   }
   else {
     wait_all_MPI(level);
   }


   //!----------------------------------------------------------------
   //! 2) get cells in J direction
   //!----------------------------------------------------------------
   add_GN_send_MPI     (level, _Jp1_, field_type);
   add_GN_equal_process(level, _Jm1_, field_type);
   add_GN_recv_MPI     (level, _Jm1_, field_type);
   if(sync_mpi_send_rec) {
     wait_SendFields_MPI(level, NO_minDIREC_SEND, _Jp1_);
     synchronize_allProcesses_MPI();
   }
   else {
     wait_all_MPI(level);
   }

   //!----------------------------------------------------------------
   //! 3) get cells in I direction
   //!----------------------------------------------------------------
   add_GN_send_MPI     (level, _Ip1_, field_type);
   add_GN_equal_process(level, _Im1_, field_type);
   add_GN_recv_MPI     (level, _Im1_, field_type);
   if(sync_mpi_send_rec) {
    wait_SendFields_MPI(level, NO_minDIREC_SEND, _Ip1_);
    synchronize_allProcesses_MPI();
   }
   else {
     wait_all_MPI(level);
   } 

}


/*old version by ck
//!--------------------------------------------------------
//!- send_injected_Ions_to_Blks:
//!--------------------------------------------------------
void CHybrid::send_injected_Ions_to_Blks(INT32            species,
					 PARTICLE_REAL   start_weight_HI,
					 PARTICLE_REAL** positions,
					 PARTICLE_REAL** velocities)
{


	if(num_Charged_Species <= species)
	{
	    log_file << endl
		     << " ERROR:" <<endl
		     << " Cannot inject Obstalce Ions-Species " << species <<endl
		     << " num_Charged_Species = " << num_Charged_Species << endl 
		     << " Ions-Species must be smaller than "<<num_Charged_Species <<"." << endl << endl;
		exit(1);
	}


	particle new_ion;


	PARTICLE_REAL vth[3];
	D_REAL r_ion[3], *x_cell[3];
	D_REAL Length[3] = {LX,LY,LZ};
	CBlock *temp_Block, *top_level_Block;
	INT32  i_j_k, blk_index, *index_blk[3], *index_cell[3];

	num_injected_particles = 0;



	INT64* num_inserted_in_L = new INT64[MAX_LEVEL+1];
  	memset(num_inserted_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));




	//! it is important to initialize vth, oth. error in case 
	//! temperature of species is zero.
	memset(vth, 0 ,3*sizeof(PARTICLE_REAL));



	//!-----------------------------------
	//! 1) alloc variable memory
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		 index_blk[comp] = new INT32[MAX_LEVEL+1];
		index_cell[comp] = new INT32[MAX_LEVEL+1];
		    x_cell[comp] = new D_REAL[MAX_LEVEL+1];

 		memset( index_blk[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
 		memset(index_cell[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
 		memset(    x_cell[comp], 0, (MAX_LEVEL+1)* sizeof(D_REAL));
	}




	//! loop across all ions that have to be inserted
	for(INT32 ion=0; ion<obstacle_MP_num_each_TL[species]; ion++)
	{


		//!---------------------------------------------------
		//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
		//!     (even though target blk might not be refined)
		//!---------------------------------------------------

		//! transform r_ion to normalized position in Simu-Box [0;1[
		r_ion[0] = (positions[0][ion]+Box_Origin[0])/Length[0];
		r_ion[1] = (positions[1][ion]+Box_Origin[1])/Length[1];
		r_ion[2] = (positions[2][ion]+Box_Origin[2])/Length[2];


		//! calculate indices and check whether position is in Simu Box
		if(calc_BlockNodeXCell_of_Pos(index_blk, index_cell, x_cell, r_ion))
		{


			//!--------------------------------------------------
			//! 2b) climb to highest possible level at respective
			//!     position.
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


			//! Generate initial thermal velocity if required
			//! -> DO THIS AT EVERY PROCESS REGARDLESS OF WHETHER
			//!    THIS PROCESS WILL INSTERT THE ION.
			//! -> IF THIS IS DONE ONLY AT THE ION-INSERTING PROCESS
			//!    RANDOM GENERATORS WILL RUN OUT OF SYNCHRONISATION
			get_maxwellian_GSL(species, vth);

			//!-----------------------------------------------------
			//! 2c) Now target Block and Cell are known, insert
			//!     ion in respective Block.
			//!-----------------------------------------------------
			if(mpi_myRank == top_level_Block->responsible_mpi_process)
			{
	

				//! record processing time of this block
				clock_t time_start, time_finish;
				time_start = clock();


				//! set top_level_Block's level 
				INT32 level = top_level_Block->RLevel;

				//! set cell index
				i_j_k =  index_cell[0][level] *BlkNds_Z *BlkNds_Y
				        +index_cell[1][level] *BlkNds_Z
				        +index_cell[2][level];
					
				//! for negative charges:
				//! get charge density of cell
				D_REAL rho_of_cell = 0;
				if(Ion_Charges[species]<0)
				{
				  D_REAL *rho = top_level_Block->Field_Type[id_rho_n];
				  rho_of_cell = rho[i_j_k];
				}

				//! only insert new ion, if cell has
				//! charge density larger than or equal to zero
				if(rho_of_cell>=0)
				{
					
					//! reset all particle member variables
					memset(&new_ion, 0, sizeof(particle));


					//! apply velocities / coordinates
					for(INT32 comp=0; comp<3; comp++)
					{
						new_ion.rel_r[comp] = x_cell[comp][level];
						new_ion.v[comp]     = velocities[comp][ion] +vth[comp];

					}

					//! set initital weight
					new_ion.weight = start_weight_HI;


	#ifdef	TRACK_PARTICLE
					//! enumerate particle for marking in case
					//! particle tracking is activated
					new_ion.number = ion;
	#endif

					//! sort particle in new cells list
					top_level_Block->add_particle_to_pArray(species, i_j_k, &new_ion);

					//! update statistics
					num_injected_particles		 ++;
					num_total_particles		 ++;
					num_total_particles_in_L[level]  ++;
					num_inserted_in_L[level]++;

				}

				//! record particle calculation time
				time_finish = clock();
				top_level_Block->time_process_particle
					+= (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


			}
		}
	}


	INT64 local_info_values[1 +(MAX_LEVEL+1)];
	stringstream info_names[INFO_ARRAY_SIZE];

	local_info_values[0]=0;

	for(int z=0; z<= MAX_LEVEL; z++)
	local_info_values[0] += num_inserted_in_L[z];


	info_names[0] << "   ->heavy ions species " << species << " injected: ";
	

	//! particle in respective level:
	for(int z=0; z<= MAX_LEVEL; z++)
	{

		info_names[1+z] << "     -> L[" <<z<< "]:  ";
		local_info_values[1+z] = num_inserted_in_L[z];
	}



	show_information(local_info_values,
			 info_names,
			 1 +(MAX_LEVEL+1), BUILD_SUM);


	delete[] num_inserted_in_L;





}*/

//!--------------------------------------------------------
//!- send_injected_Ions_to_Blks:
//!--------------------------------------------------------
void CHybrid::send_injected_Ions_to_Blks(INT32            species,
					 PARTICLE_REAL   start_weight_HI,
					 PARTICLE_REAL** positions,
					 PARTICLE_REAL** velocities,
					 const INT32* num_particles_to_inject)
{
	if(num_Particle_Species <= species)
	{
	    log_file << endl
		     << " ERROR:" <<endl
		     << " Cannot inject Obstalce Ions-Species " << species <<endl
		     << " num_Charged_Species = " << num_Particle_Species << endl 
		     << " Ions-Species must be smaller than "<<num_Particle_Species <<"." << endl << endl;
		exit(1);
	}

	particle new_ion;


	PARTICLE_REAL vth[3];
	D_REAL r_ion[3], *x_cell[3];
	D_REAL Length[3] = {LX,LY,LZ};
	CBlock *temp_Block, *top_level_Block;
	INT32  i_j_k, blk_index, *index_blk[3], *index_cell[3];

	num_injected_particles = 0;

	INT64* num_inserted_in_L = new INT64[MAX_LEVEL+1];
  	memset(num_inserted_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));

	//! it is important to initialize vth, oth. error in case 
	//! temperature of species is zero.
	memset(vth, 0 ,3*sizeof(PARTICLE_REAL));
	//!-----------------------------------
	//! 1) alloc variable memory
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		 index_blk[comp] = new INT32[MAX_LEVEL+1];
		index_cell[comp] = new INT32[MAX_LEVEL+1];
		    x_cell[comp] = new D_REAL[MAX_LEVEL+1];

 		memset( index_blk[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
 		memset(index_cell[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
 		memset(    x_cell[comp], 0, (MAX_LEVEL+1)* sizeof(D_REAL));
	}

	//! loop across all ions that have to be inserted
	for(INT32 ion=0; ion<num_particles_to_inject[species]; ion++)
	{
 
		//!---------------------------------------------------
		//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
		//!     (even though target blk might not be refined)
		//!---------------------------------------------------

		//! transform r_ion to normalized position in Simu-Box [0;1[
		r_ion[0] = (positions[0][ion]+Box_Origin[0])/Length[0];
		r_ion[1] = (positions[1][ion]+Box_Origin[1])/Length[1];
		r_ion[2] = (positions[2][ion]+Box_Origin[2])/Length[2];


		//! calculate indices and check whether position is in Simu Box
		if(calc_BlockNodeXCell_of_Pos(index_blk, index_cell, x_cell, r_ion))
		{


			//!--------------------------------------------------
			//! 2b) climb to highest possible level at respective
			//!     position.
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


			//! Generate initial thermal velocity if required
			//! -> DO THIS AT EVERY PROCESS REGARDLESS OF WHETHER
			//!    THIS PROCESS WILL INSTERT THE ION.
			//! -> IF THIS IS DONE ONLY AT THE ION-INSERTING PROCESS
			//!    RANDOM GENERATORS WILL RUN OUT OF SYNCHRONISATION
			get_maxwellian_GSL(species, vth);

			//!-----------------------------------------------------
			//! 2c) Now target Block and Cell are known, insert
			//!     ion in respective Block.
			//!-----------------------------------------------------
			if(mpi_myRank == top_level_Block->responsible_mpi_process)
			{
	

				//! record processing time of this block
				clock_t time_start, time_finish;
				time_start = clock();


				//! set top_level_Block's level 
				INT32 level = top_level_Block->RLevel;

				//! set cell index
				i_j_k =  index_cell[0][level] *BlkNds_Z *BlkNds_Y
				        +index_cell[1][level] *BlkNds_Z
				        +index_cell[2][level];
					
				//! for negative charges:
				//! get charge density of cell
				D_REAL rho_of_cell = 0;
				if(Ion_Charges[species]<0)
				{
				  D_REAL *rho = top_level_Block->Field_Type[id_rho_n];
				  rho_of_cell = rho[i_j_k];
				}

				//! only insert new ion, if cell has
				//! charge density larger than or equal to zero
				if(rho_of_cell>=0)
				{
					
					
					
					//! reset all particle member variables
					memset(&new_ion, 0, sizeof(particle));


					//! apply velocities / coordinates
					for(INT32 comp=0; comp<3; comp++)
					{
						new_ion.rel_r[comp] = x_cell[comp][level];
	#ifdef TRACK_PARTICLE
						new_ion.v[comp]     = velocities[comp][ion];
	#else
						new_ion.v[comp]     = velocities[comp][ion] +vth[comp];
	#endif
					}

					//! set initital weight
					new_ion.weight = start_weight_HI;


	#ifdef	TRACK_PARTICLE
					//! enumerate particle for marking in case
					//! particle tracking is activated
					//! +1 to start with one to have number=0 as unmarked particles 
					//! and num_particles_to_inject[species] marked particles
					new_ion.number = ion+1;
	#endif

					//! sort particle in new cells list
					top_level_Block->add_particle_to_pArray(species, i_j_k, &new_ion);

					//! update statistics
					num_injected_particles		 ++;
					num_total_particles		 ++;
					num_total_particles_in_L[level]  ++;
					num_inserted_in_L[level]++;

				}

				//! record particle calculation time
				time_finish = clock();
				top_level_Block->time_process_particle
					+= (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


			}
		}
	}

	// Commented out the following....
	// Inject Callisto Particle to Track took
	// too much time to write all this info to proc0 logfile:
       	
	//	INT64 local_info_values[1 +(MAX_LEVEL+1)];
	//	stringstream info_names[INFO_ARRAY_SIZE];

	//	local_info_values[0]=0;
	//	for(int z=0; z<= MAX_LEVEL; z++)
	  //	local_info_values[0] += num_inserted_in_L[z];


	//	info_names[0] << "   ->heavy ions species " << species << " injected: ";
	
	////! particle in respective level:
	//	for(int z=0; z<= MAX_LEVEL; z++)
	//	{
	//
	//		info_names[1+z] << "     -> L[" <<z<< "]:  ";
	//		local_info_values[1+z] = num_inserted_in_L[z];
	//	}

	//	show_information(local_info_values,
	//			 info_names,
	//			 1 +(MAX_LEVEL+1), BUILD_SUM);

	delete[] num_inserted_in_L;



}

//!--------------------------------------------------------
//!- send_injected_Ions_to_Blks_calc_weight:
//!--------------------------------------------------------
void CHybrid::send_injected_Ions_to_Blks_calc_weight(INT32            species,
						    PARTICLE_REAL** positions,
						    PARTICLE_REAL** velocities)
{


	if(num_Charged_Species <= species)
	{
	    log_file << endl
		     << " ERROR:" <<endl
		     << " Cannot inject Obstalce Ions-Species " << species <<endl
		     << " num_Charged_Species = " << num_Charged_Species << endl 
		     << " Ions-Species must be smaller than "<<num_Charged_Species <<"." << endl << endl;
		exit(1);
	}


	particle new_ion;


	PARTICLE_REAL vth[3];
	D_REAL r_ion[3], *x_cell[3];
	D_REAL Length[3] = {LX,LY,LZ};
	CBlock *temp_Block, *top_level_Block;
	INT32  i_j_k, blk_index, *index_blk[3], *index_cell[3];

	//! Distance of particle to Box Origin
	D_REAL distance_sqr;
	//! total weight of all non-normalised particle weights
	WEIGHT_REAL total_non_normalised_weight = 0;
	  
	num_injected_particles = 0;



	INT64* num_inserted_in_L = new INT64[MAX_LEVEL+1];
  	memset(num_inserted_in_L,0,(MAX_LEVEL+1)*sizeof(INT64));




	//! it is important to initialize vth, oth. error in case 
	//! temperature of species is zero.
	memset(vth, 0 ,3*sizeof(PARTICLE_REAL));



	//!-----------------------------------
	//! 1) alloc variable memory
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		 index_blk[comp] = new INT32[MAX_LEVEL+1];
		index_cell[comp] = new INT32[MAX_LEVEL+1];
		    x_cell[comp] = new D_REAL[MAX_LEVEL+1];

 		memset( index_blk[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
 		memset(index_cell[comp], 0, (MAX_LEVEL+1)* sizeof(INT32));
 		memset(    x_cell[comp], 0, (MAX_LEVEL+1)* sizeof(D_REAL));
	}
	//! Alloc memory for weight of MP
	WEIGHT_REAL *weight_ions;
	weight_ions = new WEIGHT_REAL[obstacle_MP_num_each_TL[species]];
	memset(weight_ions,0,obstacle_MP_num_each_TL[species]*sizeof(WEIGHT_REAL));

	//!-----------------------------------
	//! 1b) calc weight of ions
	//!-----------------------------------
	
	//! Calc non normalised weight for each ion
	//! Sum over all non-normalised weight
	for(INT32 ion=0; ion<obstacle_MP_num_each_TL[species]; ion++)
	{
	  distance_sqr = positions[0][ion]*positions[0][ion]+positions[1][ion]*positions[1][ion]+positions[2][ion]*positions[2][ion];
	  weight_ions[ion] = 1./distance_sqr;
	  total_non_normalised_weight+=weight_ions[ion];
	}
	
	//! loop across all ions that have to be inserted
	for(INT32 ion=0; ion<obstacle_MP_num_each_TL[species]; ion++)
	{


		//!---------------------------------------------------
		//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
		//!     (even though target blk might not be refined)
		//!---------------------------------------------------

		//! transform r_ion to normalized position in Simu-Box [0;1[
		r_ion[0] = (positions[0][ion]+Box_Origin[0])/Length[0];
		r_ion[1] = (positions[1][ion]+Box_Origin[1])/Length[1];
		r_ion[2] = (positions[2][ion]+Box_Origin[2])/Length[2];


		//! calculate indices and check whether position is in Simu Box
		if(calc_BlockNodeXCell_of_Pos(index_blk, index_cell, x_cell, r_ion))
		{


			//!--------------------------------------------------
			//! 2b) climb to highest possible level at respective
			//!     position.
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


			//! Generate initial thermal velocity if required
			//! -> DO THIS AT EVERY PROCESS REGARDLESS OF WHETHER
			//!    THIS PROCESS WILL INSTERT THE ION.
			//! -> IF THIS IS DONE ONLY AT THE ION-INSERTING PROCESS
			//!    RANDOM GENERATORS WILL RUN OUT OF SYNCHRONISATION
			get_maxwellian_GSL(species, vth);

			//!-----------------------------------------------------
			//! 2c) Now target Block and Cell are known, insert
			//!     ion in respective Block.
			//!-----------------------------------------------------
			if(mpi_myRank == top_level_Block->responsible_mpi_process)
			{
	

				//! record processing time of this block
				clock_t time_start, time_finish;
				time_start = clock();


				//! set top_level_Block's level 
				INT32 level = top_level_Block->RLevel;

				//! set cell index
				i_j_k =  index_cell[0][level] *BlkNds_Z *BlkNds_Y
				        +index_cell[1][level] *BlkNds_Z
				        +index_cell[2][level];

				//! for negative charges:
				//! get charge density of cell
				D_REAL rho_of_cell = 0;
				if(Ion_Charges[species]<0)
				{
				  D_REAL *rho = top_level_Block->Field_Type[id_rho_n];
				  rho_of_cell = rho[i_j_k];
				}

				//! only insert new ion, if cell has
				//! charge density larger than or equal to zero
				if(rho_of_cell>=0)
				{

					//! reset all particle member variables
					memset(&new_ion, 0, sizeof(particle));


					//! apply velocities / coordinates
					for(INT32 comp=0; comp<3; comp++)
					{
						new_ion.rel_r[comp] = x_cell[comp][level];
						new_ion.v[comp]     = velocities[comp][ion] +vth[comp];

					}

					//! set initital weight
					//! Set nomalised weight for each ion
					new_ion.weight = weight_ions[ion] * obstacle_MP_weight_each_t0[species] * dt / total_non_normalised_weight;


	#ifdef	TRACK_PARTICLE
					//! enumerate particle for marking in case
					//! particle tracking is activated
					new_ion.number = ion;
	#endif

					//! sort particle in new cells list
					top_level_Block->add_particle_to_pArray(species, i_j_k, &new_ion);

					//! update statistics
					num_injected_particles		 ++;
					num_total_particles		 ++;
					num_total_particles_in_L[level]  ++;
					num_inserted_in_L[level]++;

				}

				//! record particle calculation time
				time_finish = clock();
				top_level_Block->time_process_particle
					+= (double(time_finish)-double(time_start))/CLOCKS_PER_SEC;


			}
		}
	}


	INT64 local_info_values[1 +(MAX_LEVEL+1)];
	stringstream info_names[INFO_ARRAY_SIZE];

	local_info_values[0]=0;

	for(int z=0; z<= MAX_LEVEL; z++)
	local_info_values[0] += num_inserted_in_L[z];


	info_names[0] << "   ->heavy ions species " << species << " injected: ";
	

	//! particle in respective level:
	for(int z=0; z<= MAX_LEVEL; z++)
	{

		info_names[1+z] << "     -> L[" <<z<< "]:  ";
		local_info_values[1+z] = num_inserted_in_L[z];
	}



	show_information(local_info_values,
			 info_names,
			 1 +(MAX_LEVEL+1), BUILD_SUM);


	delete[] num_inserted_in_L;
	delete[] weight_ions;

}
