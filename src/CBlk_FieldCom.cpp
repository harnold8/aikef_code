/***************************************************************************
 *   Copyright (C) 2009 by Joachim Müller   *
 *   joa.mueller@tu-bs.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/



#include "CBlk.h"
#include "defines.h"
#include "parameters.h"
#include "absolute_Globals.h"
#include "utils.h"

#include <iostream>
#include <fstream>

#include <math.h>

#define minDirec	0
#define plsDirec	1

using namespace std;

extern INT32 *COMPs_FType;

//! MPI related variables
extern INT32 mpi_myRank;
extern ofstream log_file;





//!************************************************************************************************//
//! ******************************** Copy Ghost Nodes *********************************************//
//!************************************************************************************************//

//! extern array of function pointers for
//! cp GN procedure
extern void (*send_GN_MPI[3])(bool direc,
			      CBlock* dest_Blk,
			      CBlock* src_Blk,
			      INT32 id_package,
			      INT32 field_type);

extern void (*recv_GN_MPI[3])(bool direc,
			      CBlock* dest_Blk,
			      CBlock* src_Blk,
			      INT32 id_package,
			      INT32 field_type);


extern void (*get_GN_equal_process[3])(bool direc,
				       CBlock* dest_Blk,
				       CBlock* src_Blk,
				       INT32    field_type);

//!---------------------------------------------------------------//
//! cp_GN_equal_process:					  //
//!---------------------------------------------------------------//
void CBlock::cp_GN_equal_process(INT32 _m1_, INT32 _p1_, INT32 field_type)
{


	//! Explenation how to choose indices of function pointer:
	//! _m1_/2 and _m1_/2 set i,j or k component
	//! e.g. Im1=0 -> =0
	//! e.g. Ip1=1 -> =0
	//! e.g. Jm1=2 -> =1
	//! e.g. Kp1=5 -> =2

	//! - use function pointer to decide which function to call
	//!   definition in utils.h, initialization in CHybrid.h
	//! - use %gathNEIB to distinguish ordinary vs. gather Neighbours
	if(Neighbour[_m1_] && mpi_myRank==Neighbour[_m1_]->responsible_mpi_process)
	get_GN_equal_process[(_m1_%gathNEIB)/2](minDirec, this, Neighbour[_m1_], field_type);


	if(Neighbour[_p1_] && mpi_myRank==Neighbour[_p1_]->responsible_mpi_process)
	get_GN_equal_process[(_p1_%gathNEIB)/2](plsDirec, this, Neighbour[_p1_], field_type);

}



//!---------------------------------------------------------------//
//! cp_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void CBlock::cp_GN_send_MPI(INT32 _m1_, INT32 _p1_, INT32 field_type)
{


	//! set outgoing package ids
	//! it should be fine to use same id's in i,j,k direction since messages
	//! cannot overtake




	//! - Use function pointer to decide which function to call
	//!   definition in utils.h, initialization in CHybrid.h
	//! - Explenation how to choose indices of function pointer:
	//!   _m1_/2 and _m1_/2 set i,j or k component
	//!   e.g. Im1=0 -> =0
	//!   e.g. Kp1=5 -> =2
	if(Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	{

		//! - use %gathNEIB to distinguish ordinary vs. gather Neighbours
		INT32 tag_to_m1 = mpi_tag +0*total_num_mpi_tags;
		send_GN_MPI[(_m1_%gathNEIB)/2](minDirec, Neighbour[_m1_], this, tag_to_m1, field_type);
	}

	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	{

		//! initialize tag
		INT32 tag_to_p1 = mpi_tag +1*total_num_mpi_tags;
		send_GN_MPI[(_p1_%gathNEIB)/2](plsDirec, Neighbour[_p1_], this, tag_to_p1, field_type);
	}

}


//!---------------------------------------------------------------//
//! cp_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void CBlock::cp_GN_recv_MPI(INT32 _m1_, INT32 _p1_, INT32 field_type)
{

	//! Explenation how to choose indices of function pointer:
	//! _m1_/2 and _m1_/2 set i,j or k component
	//! e.g. Im1=0 -> =0
	//! e.g. Kp1=5 -> =2


	//! use function pointer to decide which function to call
	//! definition in utils.h, initialization in CHybrid.h
	if(Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	{

		//! - initialize tag INSIDE brackets, else Neighbour might not
		//!    exist !!! -> seghmentatiton fault
		//! - use %gathNEIB to distinguish ordinary vs. gather Neighbours
		INT32 tag_from_m1 = Neighbour[_m1_]->mpi_tag +1*total_num_mpi_tags;
		recv_GN_MPI[(_m1_%gathNEIB)/2](minDirec, this, Neighbour[_m1_], tag_from_m1, field_type);

	}


	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	{

		//! - initialize tag INSIDE brackets, else Neighbour might not
		//!   exist !!! -> seghmentatiton fault
		//! - use %gathNEIB to distinguish ordinary vs. gather Neighbours
		INT32 tag_from_p1 = Neighbour[_p1_]->mpi_tag +0*total_num_mpi_tags;
		recv_GN_MPI[(_p1_%gathNEIB)/2](plsDirec, this, Neighbour[_p1_], tag_from_p1, field_type);
	}

}


//!---------------------------------------------------------------//
//! wait_SendFields_MPI:						  //
//!---------------------------------------------------------------//
void CBlock::wait_SendFields_MPI(INT32 _m1_, INT32 _p1_)
{

	//! force MPI to deallocate memory of Last Send
	//! Code will run without this command but result
	//! in a memory leak

	//! Minus will not be used for gather, intercept for this by
	//! specifyin negative value for _m1_
	if(_m1_>=0 && Neighbour[_m1_] && mpi_myRank!=Neighbour[_m1_]->responsible_mpi_process)
	req_is_MPI_send_completed[minDirec].Wait();

	if(Neighbour[_p1_] && mpi_myRank!=Neighbour[_p1_]->responsible_mpi_process)
	req_is_MPI_send_completed[plsDirec].Wait();

// 	MPI::Request::Waitall(NUM_REQUESTS, req_is_MPI_send_completed);


}



//!---------------------------------------------------------------//
//! wait_all_MPI:									  //
//!---------------------------------------------------------------//
void CBlock::wait_all_MPI(void)
{




	MPI::Request::Waitall(NUM_REQUESTS, req_is_MPI_send_completed);
	//!	Description
	//! 
	//!	This subroutine blocks until all operations associated with active handles
	//!	in the list complete, and returns the status of each operation.
	//!	array_of_requests and array_of statuses contain count entries.
	//! 
	//!	The ith entry in array_of_statuses is set to the return status of the ith
	//!	operation. Requests created by nonblocking operations are deallocated and
	//!	the corresponding handles in the array are set to MPI_REQUEST_NULL. If
	//!	array_of_requests contains null or inactive handles, MPI_WAITALL sets the
	//!	status of each one to empty.
	//! 
	//!	MPI_WAITALL(count, array_of_requests, array_of_statuses) has the same
	//!	effect as the execution of MPI_WAIT(array_of_requests[i],
	//!	array_of_statuses[i]) for i = 0, 1, ..., (count-1), in some arbitrary
	//!	order. MPI_WAITALL with an array of length one is equivalent to MPI_WAIT.
	//! 
	//!	The error fields are never modified unless the function gives a return code
	//!	of MPI_ERR_IN_STATUS. In which case, the error field of every MPI_Status is
	//!	modified to reflect the result of the corresponding request.
	//! 
	//!	Passing MPI_STATUSES_IGNORE for the array_of statuses argument causes PE
	//!	MPI to skip filling in the status fields. By passing this value for
	//!	array_of statuses, you can avoid having to allocate a status object array
	//!	in programs that do not need to examine the status fields.


}



//!************************************************************************************************//
//!********************************* ADD Ghost Nodes **********************************************//
//!************************************************************************************************//

//! extern array of function pointers for
//! cp GN procedure
extern void (*send_add_GN_MPI[3])(CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32 id_package,
				INT32 field_type);

extern void (*recv_add_GN_MPI[3])(CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32 id_package,
				INT32 field_type);


extern void (*get_add_GN_equal_process[3])(CBlock* dest_Blk,
					 CBlock* src_Blk,
					 INT32    field_type);

//!---------------------------------------------------------------//
//! add_GN_equal_process:					  //
//!---------------------------------------------------------------//
void CBlock::add_GN_equal_process(INT32 _m1_, INT32 field_type)
{

	//! Explenation how to choose indices of function pointer:
	//! _m1_/2 and _m1_/2 set i,j or k component
	//! e.g. Im1=0 -> =0
	//! e.g. Kp1=5 -> =2

	//! use function pointer to decide which function to call
	//! definition in utils.h, initialization in CHybrid.h
	if(gather_Neighbour[_m1_] && mpi_myRank==gather_Neighbour[_m1_]->responsible_mpi_process)
	get_add_GN_equal_process[_m1_/2](this, gather_Neighbour[_m1_], field_type);

}


//!---------------------------------------------------------------//
//! add_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void CBlock::add_GN_send_MPI(INT32 _p1_, INT32 field_type)
{


	//! Explenation how to choose indices of function pointer:
	//! _m1_/2 and _m1_/2 set i,j or k component
	//! e.g. Ip1=0 -> =0
	//! e.g. Kp1=5 -> =2


	//! use function pointer to decide which function to call
	//! definition in utils.h, initialization in CHybrid.h
	if(gather_Neighbour[_p1_] && mpi_myRank!=gather_Neighbour[_p1_]->responsible_mpi_process)
	{

		//! set outgoing package ids
		INT32 tag_to_p1 = mpi_tag;
		send_add_GN_MPI[_p1_/2](gather_Neighbour[_p1_], this, tag_to_p1, field_type);
	}
}


//!---------------------------------------------------------------//
//! add_GN_recv_MPI:						  //
//!---------------------------------------------------------------//
void CBlock::add_GN_recv_MPI(INT32 _m1_, INT32 field_type)
{


	//! Explenation how to choose indices of function pointer:
	//! _m1_/2 and _m1_/2 set i,j or k component
	//! e.g. Im1=0 -> =0
	//! e.g. Km1=4 -> =2

	//! use function pointer to decide which function to call
	//! definition in utils.h, initialization in CHybrid.h
	if(gather_Neighbour[_m1_] && mpi_myRank!=gather_Neighbour[_m1_]->responsible_mpi_process)
	{

		//! initialize tag INSIDE brackets, else Neighbour might not
		//! exist !!! -> seghmentatiton fault
		INT32 tag_from_p1 = gather_Neighbour[_m1_]->mpi_tag;

		recv_add_GN_MPI[_m1_/2](this, gather_Neighbour[_m1_], tag_from_p1, field_type);
	}

}
