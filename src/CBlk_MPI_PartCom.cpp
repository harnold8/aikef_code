



#include "CBlk.h"
#include "defines.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>




//!-----------------------------------------------------------
//! send_MPiC_Package_MPI
//!-----------------------------------------------------------
void CBlock::send_MPiC_Package_MPI(INT32  id_request,
				   INT32  id_package,
				   INT32  dest_process,
				   INT32  num_entries_in_InfoPackage,
				   INT32* &MPiC_Info_package)
{


	//! force MPI to deallocate memory of Last Send
	//! Code will run without this command but result
	//! in a memory leak
	//! MEANWHILE THIS SHOULD BE INTERCEPTED BY "wait_SendFields_MPI".
// 	if(req_is_MPI_send_completed[direc])
// 	req_is_MPI_send_completed[direc].Wait();


	//! maybe it's better just to use MPI's BYTE datatype and sizeof() lile memcpy
	//! -> easy to try
	req_is_MPI_send_completed[id_request] = mpi_myComm.IB_OR_IS_SEND(MPiC_Info_package,
									 num_entries_in_InfoPackage *sizeof(INT32),
									 MPI::BYTE,// INT,
									 dest_process,
									 id_package);

}



//!-----------------------------------------------------------
//! send_partPackage_MPI
//!-----------------------------------------------------------
void CBlock::send_partPackage_MPI(INT32 id_request,
				  INT32 id_package,
				  INT32  dest_process,
				  INT32 num_entries_in_InfoPackage,
				  particle* &particle_package,
				  INT32* &MPiC_Info_package)
{




	//! "num_entries_in_PartPackage" is stored at ultimate position of MPiC_Info_package
	INT32 num_particle_in_PartPackage = MPiC_Info_package[num_entries_in_InfoPackage-1];


	req_is_MPI_send_completed[id_request]
		= mpi_myComm.IB_OR_IS_SEND(particle_package,
					   num_particle_in_PartPackage*sizeof(particle),
					   MPI::BYTE, //PARTICLE_MPI,
					   dest_process,
					   id_package);




}


//!-----------------------------------------------------------
//! receive_infoPackage_MPI
//!-----------------------------------------------------------
void CBlock::receive_infoPackage_MPI(INT32  id_package,
				     INT32  src_process,
				     INT32  num_entries_in_InfoPackage,
				     INT32* &MPiC_Info_package)
{



	//! - since neighbour is not handled by this process, memory has not been allocated yet
	//! - "1" is stored at ultimate position of MPiC_Info_package, indicating for MPI that
	//!   array has been received and may be sorted to block (make sure it is zero before)
	MPiC_Info_package = new INT32[num_entries_in_InfoPackage];


	//!  receive MPiC Info package
	mpi_myComm.Recv(MPiC_Info_package,
			num_entries_in_InfoPackage *sizeof(INT32),
			MPI::BYTE,
			src_process,
			id_package);


}


//!-----------------------------------------------------------
//! receive_partPackage_MPI
//!-----------------------------------------------------------
void CBlock::receive_partPackage_MPI(INT32  id_package,
				     INT32  src_process,
				     INT32  num_entries_in_InfoPackage,
				     particle* &particle_package,
				     INT32* &MPiC_Info_package)
{



	//! "num_particle_in_PartPackage" is stored at pen ultimate position of MPiC_Info_package
	INT32 num_particle_in_PartPackage = MPiC_Info_package[num_entries_in_InfoPackage-1];

	//! since neighbour is not handled by this process, memory has not been allocated yet
	particle_package = new particle[num_particle_in_PartPackage];


	//!  receive particle package
	mpi_myComm.Recv(particle_package,
			num_particle_in_PartPackage *sizeof(particle),
			MPI::BYTE,
			src_process,
			id_package);


}

//!-----------------------------------------------------------
//! compose_recvd_fields
//!-----------------------------------------------------------
void CBlock::compose_recvd_fields(INT32 id_field,
				  INT32 num_distinct_mpi_procs,
				  D_REAL** buffer_proc)
{


	//! set pointer do destination field
	D_REAL* parent_field = Field_Type[id_field];


	//! 1)
	//! in case moments are smoothed to parent, entire field can be set to zero

	//! 2)
	//! in case EM Fields are smoothed to parent, only

// 	memset(parent_field, 0, COMPs_FType[id_field]*num_nodes_in_block*sizeof(D_REAL));


	//! use a,b,c as oct indices
	for(INT32 proc=0; proc<num_distinct_mpi_procs; proc++)
	{


		D_REAL* received_field = buffer_proc[proc];


		for(INT32 u=0; u<BlkNds_X; u++)
		 for(INT32 v=0; v<BlkNds_Y; v++)
		  for(INT32 w=0; w<BlkNds_Z; w++)
		  {

			INT32 u_v_w =  u*BlkNds_Z*BlkNds_Y
					+v*BlkNds_Z
					+w;

			//! NOTE:
			//! since parent field region of two different octs may
			//! overlab in gather method,fields have to be added !!
			//! -> see "add_field_to_parent" method !!!
			for(INT32 comp=0; comp<COMPs_FType[id_field]; comp++)
				parent_field[u_v_w +comp*num_nodes_in_block]
			   += received_field[u_v_w +comp*num_nodes_in_block];

		   }
	  }
}


//!-----------------------------------------------------------
//! compose_recvd_fields
//!-----------------------------------------------------------
void CBlock::copy_recvd_fields(INT32 id_field,
			       INT32* process_buffer_id,
			       D_REAL** buffer_proc)
{


	//! in case EM Fields are smoothed to parent, nothing has to be set zero
	//! since values are copied

	//! set pointer do destination field
	D_REAL* parent_field = Field_Type[id_field];


	//! use a,b,c as oct indices
	for(INT32 a = 0; a < 1; a++)
	 for(INT32 b = 0; b < 1; b++)
	  for(INT32 c = 0; c < 1; c++)
	  {
	

		INT32 id_oct = a*2*2 +b*2 +c;


		//! only copy when
		//! - child exists
		//! - is not my rank
		if(child_array[id_oct] && child_array[id_oct]->responsible_mpi_process!=mpi_myRank)
		{

			INT32 id_buffer = process_buffer_id[child_array[id_oct]->responsible_mpi_process];

			D_REAL* received_field = buffer_proc[id_buffer];
	
	
	
			//! START:
			//! eg: BlkNds=6
			//! 1)
			//! a=0
			//! -> start = 1
			//! 2)
			//! a=1
			//! -> start = 3
			const INT32 oct_start[3] = { 1 +a*(BlkNds_X/2-1),
						    1 +b*(BlkNds_Y/2-1),
						    1 +c*(BlkNds_Z/2-1)};
		
		
			//! END:
			//! eg: BlkNds=6
			//! a=0
			//! -> end = 3
			//! a=1
			//! -> end = 5
			const INT32 oct_end[3]={ BlkNds_X/2 +a*(BlkNds_X/2-1),
						BlkNds_Y/2 +b*(BlkNds_Y/2-1),
						BlkNds_Z/2 +c*(BlkNds_Z/2-1)};
		
			for(INT32 i = oct_start[0]; i < oct_end[0]; i++)
			 for(INT32 j = oct_start[1]; j < oct_end[1]; j++)
			  for(INT32 k = oct_start[2]; k < oct_end[2]; k++)
			  {
		
	
				INT32 i_j_k =   i*BlkNds_Y*BlkNds_Z
					      +j*BlkNds_Z
					      +k;
		
				for(INT32 comp=0; comp<COMPs_FType[id_field]; comp++)
				parent_field[i_j_k +comp*num_nodes_in_block] = received_field[i_j_k +comp*num_nodes_in_block];
		
			 }

		}
	}
}






