

#include "defines.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <ctime>

// #include <stdio.h>
#include <mpi.h>





//!************************************************************************************************//
//! ******************************** FULL BLOCK ***************************************************//
//!************************************************************************************************//


//!---------------------------------------------------------------//
//! send_MPI_Field:						  //
//!---------------------------------------------------------------//
void send_MPI_Field(INT32 id_field,
		    INT32 id_request,
		    INT32 id_dest_process,
		    INT32 tag,
		    CBlock* src_Blk)
{



	//! calculate size of outgoing buffer
	INT32 buffer_size = COMPs_FType[id_field] *num_nodes_in_block;


	//! - send buffer to responsible process
	//! (may try Bsend, IB_OR_IS_SEND or Ibsend also)
	src_Blk->req_is_MPI_send_completed[id_request]

		= mpi_myComm.IB_OR_IS_SEND(src_Blk->Field_Type[id_field],
					   buffer_size,
					   MPI_D_REAL,
					   id_dest_process,
					   tag);
}




//!---------------------------------------------------------------//
//! recv_MPI_Field:								  //
//!---------------------------------------------------------------//
void recv_MPI_Field(INT32 id_field,
		    INT32 id_src_process,
		    INT32 tag,
		    D_REAL* Field)
{

	//! calculate size of incoming buffer
	INT32 buffer_size = COMPs_FType[id_field] *num_nodes_in_block;

	//! - receive buffer from responsible process
	mpi_myComm.Recv(Field,
			buffer_size,
			MPI_D_REAL,
			id_src_process,
			tag);


}

//!************************************************************************************************//
//! ******************************** COPY SEND  ***************************************************//
//!************************************************************************************************//


//!---------------------------------------------------------------//
//! send_GN_MPI:						  //
//!---------------------------------------------------------------//
//! define array of function pointers for
//! cp GN procedure
//! (they are initialzed in CHybrid_Init.cpp)
void (*send_GN_MPI[3])(bool direc,
		       CBlock* dest_Blk,
			CBlock* src_Blk,
			INT32 id_package,
			INT32 field_type);



//!---------------------------------------------------------------//
//! cp_Ip1_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void cp_I_GN_send_MPI(bool direc,
		      CBlock* dest_Blk,
		      CBlock* src_Blk,
		      INT32 id_package,
		      INT32 field_type)
{
	//! Do not cp Full Cross Section
	//! leaving out k_row's at BlkNds_X=0 & BlkNds_X-2 doesnt matter.
	//! BUT: Full k_row has to be copied else.


	//! abbreviations
	INT32 XCS_src;
	const INT32 scalar_buffer_size  =  BlkNds_Y *BlkNds_Z;

	//! decide whether to send first or last physical cells
	if(!direc) XCS_src =             1* BlkNds_Y *BlkNds_Z;
	else       XCS_src = (BlkNds_X-2) * BlkNds_Y *BlkNds_Z;


	//! delete memory from previous send
	delete[] src_Blk->send_buffer[direc];


	//! alloc buffer to store send data
	src_Blk->send_buffer[direc] = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *src_Field = src_Blk->Field_Type[field_type];

	//! cp Full Cross Section
	//! do this for each component
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	memcpy(src_Blk->send_buffer[direc] +comp*scalar_buffer_size,
	       src_Field   +comp*num_nodes_in_block + XCS_src,
	       BlkNds_Y *BlkNds_Z *sizeof(D_REAL));


	//! - send buffer to responsible process
	//! (may try Bsend, IB_OR_IS_SEND or Ibsend also)
	src_Blk->req_is_MPI_send_completed[direc] =
	 mpi_myComm.IB_OR_IS_SEND(src_Blk->send_buffer[direc],
			 scalar_buffer_size*COMPs_FType[field_type],
			 MPI_D_REAL,
			 dest_Blk->responsible_mpi_process,
			 id_package);



}


//!---------------------------------------------------------------//
//! cp_Jm1_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void cp_J_GN_send_MPI(bool direc,
		      CBlock* dest_Blk,
		      CBlock* src_Blk,
		      INT32 id_package,
		      INT32 field_type)
{

	//! Do not cp Full Cross Section
	//! leaving out k_row's at BlkNds_X=0 & BlkNds_X-2 doesnt matter.
	//! BUT: Full k_row has to be copied else.

	//! abbreviations
	INT32 XCS_pos, XCSb_pos, k_src_row;
	const INT32 scalar_buffer_size = (BlkNds_X-2)*BlkNds_Z;

	//! decide whether to send first or last physical cells
	if(!direc) k_src_row = 1 *BlkNds_Z;
	else       k_src_row = (BlkNds_Y-2) *BlkNds_Z;
	

	//! delete memory from previous send
	delete[] src_Blk->send_buffer[direc];

	//! TODO:
	//! check whether vector masking performs faster

	//! alloc buffer to store send data
	src_Blk->send_buffer[direc] = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *src_Field = src_Blk->Field_Type[field_type];



	//! this loop copies the k-rows into send-buffer
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X-1; a++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		XCSb_pos = (a-1) *BlkNds_Z;
		XCS_pos  =  a    *BlkNds_Y  *BlkNds_Z;
	
		memcpy(src_Blk->send_buffer[direc] +comp*scalar_buffer_size +XCSb_pos,
		        	    src_Field  +comp*num_nodes_in_block +XCS_pos  +k_src_row,
		        BlkNds_Z *sizeof(D_REAL));
	 }	


	//! - send buffer to responsible process
	//! (may try Bsend, IB_OR_IS_SEND or Ibsend also)
	src_Blk->req_is_MPI_send_completed[direc] = mpi_myComm.IB_OR_IS_SEND(src_Blk->send_buffer[direc],
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			dest_Blk->responsible_mpi_process,
			id_package );


}



//!---------------------------------------------------------------//
//! cp_Km1_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void cp_K_GN_send_MPI(bool direc,
		      CBlock* dest_Blk,
		      CBlock* src_Blk,
		      INT32 id_package,
		      INT32 field_type)
{

	//! abbreviations
	INT32 XCS_pos, XCSb_pos, k_src_cell;
	const INT32 scalar_buffer_size = (BlkNds_X-2)*(BlkNds_Y-2);

	//! decide whether to send first or last physical cells
	if(!direc) k_src_cell = 1;
	else       k_src_cell = (BlkNds_Z-2);


	//! delete memory from previous send
	delete[] src_Blk->send_buffer[direc];

	//! alloc buffer to store received data
	src_Blk->send_buffer[direc] = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *src_Field = src_Blk->Field_Type[field_type];


	//! Do not cp Full Cross Section
	//! Leaving out border cells doesnt matter
	//! this loop sorts the nodes into send_buffer
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X-1; a++)
	 {
		XCS_pos  =  a    *BlkNds_Y *BlkNds_Z;
		XCSb_pos = (a-1) *(BlkNds_Y-2);

		for(INT32 b=1; b<BlkNds_Y-1; b++)
		src_Blk->send_buffer[direc][comp*scalar_buffer_size + XCSb_pos + (b-1)] =
			          src_Field[comp*num_nodes_in_block +  XCS_pos +  b  *BlkNds_Z  +k_src_cell];


	 }

	//! - send buffer to responsible process
	//! (may try Bsend, IB_OR_IS_SEND or Ibsend also)
	src_Blk->req_is_MPI_send_completed[direc] = mpi_myComm.IB_OR_IS_SEND(src_Blk->send_buffer[direc],
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			dest_Blk->responsible_mpi_process,
			id_package);


}




//!************************************************************************************************//
//! ****************************** RECEIVE ********************************************************//
//!************************************************************************************************//


//!---------------------------------------------------------------//
//! send_GN_MPI:						  //
//!---------------------------------------------------------------//
//! define array of function pointers for
//! cp GN procedure
//! (they are initialzed in CHybrid.h)
void (*recv_GN_MPI[3])(bool direc,
			CBlock* dest_Blk,
			CBlock* src_Blk,
			INT32 id_package,
			INT32 field_type);

//!---------------------------------------------------------------//
//! cp_I_GN_receive_MPI:					  //
//!---------------------------------------------------------------//
void cp_I_GN_recv_MPI(bool direc,
			CBlock* dest_Blk,
			CBlock* src_Blk,
			INT32 id_package,
			INT32 field_type)
{


	//! abbreviations
	INT32 XCS_dest;
	const INT32 scalar_buffer_size = BlkNds_Y*BlkNds_Z;

	//! decide whether to receive first or last cells
	if(!direc) XCS_dest =            0 * BlkNds_Y * BlkNds_Z;
	else       XCS_dest =  (BlkNds_X-1) *BlkNds_Y *BlkNds_Z;

	//! alloc buffer to store received data
	D_REAL *receive_buffer = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];

	//! - receive buffer from responsible process
	mpi_myComm.Recv(receive_buffer,
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			src_Blk->responsible_mpi_process,
			id_package);

	//! cp Full Cross Section
	//! do this for each component
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	memcpy(dest_Field     +comp*num_nodes_in_block +XCS_dest,
	       receive_buffer +comp*scalar_buffer_size,
	       BlkNds_Y *BlkNds_Z *sizeof(D_REAL));


	delete[] receive_buffer;

}


//!---------------------------------------------------------------//
//! cp_Jm1_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void cp_J_GN_recv_MPI(bool direc,
			CBlock* dest_Blk,
			CBlock* src_Blk,
			INT32 id_package,
			INT32 field_type)
{


	//! abbreviations
	INT32 XCS_pos, XCSb_pos, k_dest_row;
	const INT32 scalar_buffer_size = (BlkNds_X-2)*BlkNds_Z;

	//! decide whether to receive first or last cells
	if(!direc) k_dest_row =             0 *BlkNds_Z;
	else       k_dest_row =  (BlkNds_Y-1) *BlkNds_Z;

	//! alloc buffer to store received data
	D_REAL *receive_buffer = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];

	//! - receive buffer from responsible process
	mpi_myComm.Recv(receive_buffer,
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			src_Blk->responsible_mpi_process,
			id_package);


	//! this loop copies the k-rows into destination Field
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X-1; a++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		XCSb_pos = (a-1) *BlkNds_Z;
		XCS_pos  =  a    *BlkNds_Y  *BlkNds_Z;
	
		memcpy(dest_Field      +comp*num_nodes_in_block +XCS_pos  +k_dest_row,
		        receive_buffer +comp*scalar_buffer_size +XCSb_pos,
		        BlkNds_Z *sizeof(D_REAL));
	 }



	delete[] receive_buffer;

}


//!---------------------------------------------------------------//
//! cp_Kp1_GN_recv_MPI:						  //
//!---------------------------------------------------------------//
void cp_K_GN_recv_MPI(bool direc,
		      CBlock* dest_Blk,
		      CBlock* src_Blk,
		      INT32 id_package,
		      INT32 field_type)
{


	//! abbreviations
	INT32 XCS_pos, XCSb_pos, k_dest_cell;
	const INT32 scalar_buffer_size = (BlkNds_X-2)*(BlkNds_Y-2);

	//! decide whether to receive first or last cells
	if(!direc) k_dest_cell =  0;
	else       k_dest_cell =  BlkNds_Z-1;

	//! alloc buffer to store received data
	D_REAL *receive_buffer = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];

	//! - receive buffer from responsible process
	mpi_myComm.Recv(receive_buffer,
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			src_Blk->responsible_mpi_process,
			id_package);


	//! Do not cp Full Cross Section
	//! Leaving out border cells doesnt matter
	//! this loop sorts the nodes into send_buffer
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X-1; a++)
	 {
		XCS_pos  =  a    *BlkNds_Y *BlkNds_Z;
		XCSb_pos = (a-1) *(BlkNds_Y-2);

		for(INT32 b=1; b<BlkNds_Y-1; b++)
		dest_Field[	 comp*num_nodes_in_block +XCS_pos  +  b * BlkNds_Z  +k_dest_cell] =
		  receive_buffer[comp*scalar_buffer_size +XCSb_pos + (b-1)];

	 }

	delete[] receive_buffer;

}



//!************************************************************************************************//
//! ******************************** ADD SEND *****************************************************//
//!************************************************************************************************//


//! array of function pointers for
//! add GN procedure
void (*send_add_GN_MPI[3])(CBlock* dest_Blk,
				CBlock* src_Blk,
				INT32 id_package,
				INT32 field_type);


//!---------------------------------------------------------------//
//! add_I_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void add_I_GN_send_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type)
{

	//! for sake of simplicity send (BlkNds_Y *BlkNds_Z)
	//! rather than (BlkNds_Y-1) *(BlkNds_Z-1)


	//! abbreviations
	const INT32 XCS_src  = (BlkNds_X-1) * BlkNds_Y *BlkNds_Z;
	const INT32 scalar_buffer_size  =     BlkNds_Y *BlkNds_Z;


	//! delete memory from previous send
	delete[] src_Blk->send_buffer[1];


	//! alloc buffer to store send data
	src_Blk->send_buffer[1] = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *src_Field = src_Blk->Field_Type[field_type];

	//! cp Full Cross Section
	//! do this for each component
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	memcpy(src_Blk->send_buffer[1] +comp*scalar_buffer_size,
	     	    	 src_Field  +comp*num_nodes_in_block + XCS_src,
	       BlkNds_Y *BlkNds_Z *sizeof(D_REAL));


	//! - send buffer to responsible process
	src_Blk->req_is_MPI_send_completed[1]=mpi_myComm.IB_OR_IS_SEND(src_Blk->send_buffer[1],
			 scalar_buffer_size*COMPs_FType[field_type],
			 MPI_D_REAL,
			 dest_Blk->responsible_mpi_process,
			 id_package);


}


//!---------------------------------------------------------------//
//! add_J_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void add_J_GN_send_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type)
{


	//! abbreviations
	INT32 XCS_pos, XCSb_pos;
	const INT32 src_Plane = (BlkNds_Y-1);

	//!+1 below is to to skip the first value of k row
	const INT32 k_src_row  =  (src_Plane *BlkNds_Z) +1;
	const INT32 scalar_buffer_size = (BlkNds_X-1)*(BlkNds_Z-1);

	//! delete memory from previous send
	delete[] src_Blk->send_buffer[1];

	//! alloc buffer to store send data
	src_Blk->send_buffer[1] = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *src_Field = src_Blk->Field_Type[field_type];



	//! this loop copies the k-rows into send-buffer
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X; a++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		XCSb_pos = (a-1) *(BlkNds_Z-1);
		XCS_pos  =  a    *BlkNds_Y  *BlkNds_Z;
	
		memcpy(src_Blk->send_buffer[1] +comp*scalar_buffer_size +XCSb_pos,
				    src_Field  +comp*num_nodes_in_block +XCS_pos  +k_src_row,
		        (BlkNds_Z-1)*sizeof(D_REAL));
	 }	


	//! - send buffer to responsible process
	//! (may try Bsend, IB_OR_IS_SEND or Ibsend also)
	src_Blk->req_is_MPI_send_completed[1] = mpi_myComm.IB_OR_IS_SEND(src_Blk->send_buffer[1],
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			dest_Blk->responsible_mpi_process,
			id_package);


}



//!---------------------------------------------------------------//
//! cp_K_GN_send_MPI:						  //
//!---------------------------------------------------------------//
void add_K_GN_send_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type)
{


	//! abbreviations
	INT32 XCS_pos, XCSb_pos;

	const INT32 src_Plane = BlkNds_Z-1;
	const INT32 scalar_buffer_size = (BlkNds_X-1)*(BlkNds_Y-1);


	//! delete memory from previous send
	delete[] src_Blk->send_buffer[1];

	//! alloc buffer to store received data
	src_Blk->send_buffer[1] = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *src_Field = src_Blk->Field_Type[field_type];


	//! Do not cp Full Cross Section
	//! Leaving out border cells doesnt matter
	//! this loop sorts the nodes into send_buffer
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X; a++)
	 {
		XCS_pos  =  a    *BlkNds_Y *BlkNds_Z;
		XCSb_pos = (a-1) *(BlkNds_Y-1);

		for(INT32 b=1; b<BlkNds_Y; b++)
		src_Blk->send_buffer[1][comp*scalar_buffer_size + XCSb_pos + (b-1)] =
			      src_Field[comp*num_nodes_in_block +  XCS_pos +  b  *BlkNds_Z    +src_Plane];


	 }

	//! - send buffer to responsible process
	//! (may try Bsend, IB_OR_IS_SEND or Ibsend also)
	src_Blk->req_is_MPI_send_completed[1] = mpi_myComm.IB_OR_IS_SEND(src_Blk->send_buffer[1],
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			dest_Blk->responsible_mpi_process,
			id_package);



}

//!************************************************************************************************//
//! ******************************** ADD RECEIVE **************************************************//
//!************************************************************************************************//

//! array of function pointers for
//! add GN procedure
void (*recv_add_GN_MPI[3])(CBlock* dest_Blk,
			   CBlock* src_Blk,
			   INT32 id_package,
			   INT32 field_type);


//!---------------------------------------------------------------//
//! cp_I_GN_receive_MPI:					  //
//!---------------------------------------------------------------//
void add_I_GN_receive_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type)
{


	//! abbreviations
	INT32 Y_row;
	const INT32 XCS_dest = 1 * BlkNds_Y * BlkNds_Z;
	const INT32 scalar_buffer_size = BlkNds_Y*BlkNds_Z;

	//! TODO:
	//! check whether vector masking performs faster

	//! alloc buffer to store received data
	D_REAL *receive_buffer = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];

	//! - receive buffer from responsible process
	mpi_myComm.Recv(receive_buffer,
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			src_Blk->responsible_mpi_process,
			id_package);



	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 b=1; b<BlkNds_Y; b++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		Y_row  =  b *BlkNds_Z;

		for(INT32 c=1; c<BlkNds_Z; c++)
		    dest_Field[comp*num_nodes_in_block +XCS_dest +Y_row +c] +=
		receive_buffer[comp*scalar_buffer_size           +Y_row +c];

	 }



	delete[] receive_buffer;

}


//!---------------------------------------------------------------//
//! add_J_GN_receive_MPI:					  //
//!---------------------------------------------------------------//
void add_J_GN_receive_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type)
{


	//! abbreviations
	INT32 XCS_pos, XCSb_pos;

	const INT32 dest_Plane=1;
	const INT32 k_dest_row =  dest_Plane *BlkNds_Z;
	const INT32 scalar_buffer_size = (BlkNds_X-1)*(BlkNds_Z-1);

	//! TODO:
	//! check whether vector masking performs faster

	//! alloc buffer to store received data
	D_REAL *receive_buffer = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];

	//! - receive buffer from responsible process
	mpi_myComm.Recv(receive_buffer,
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			src_Blk->responsible_mpi_process,
			id_package);

	//! add buffer to field
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X; a++)
	 {
		//! start Pos of X=a Cross Section (XCS)
		XCSb_pos = (a-1) *(BlkNds_Z-1);
		XCS_pos  =  a    *BlkNds_Y  *BlkNds_Z;

		for(INT32 c=1; c<BlkNds_Z; c++)
		    dest_Field[comp*num_nodes_in_block  +XCS_pos +k_dest_row     + c] +=
		receive_buffer[comp*scalar_buffer_size +XCSb_pos             + (c-1)];

	 }


	delete[] receive_buffer;

}



//!---------------------------------------------------------------//
//! add_K_GN_receive_MPI:					  //
//!---------------------------------------------------------------//
void add_K_GN_receive_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type)
{


	//! abbreviations
	INT32 XCS_pos, XCSb_pos;

	const INT32 dest_Plane = 1;
	const INT32 scalar_buffer_size = (BlkNds_X-1)*(BlkNds_Y-1);

	//! alloc buffer to store received data
	D_REAL *receive_buffer = new D_REAL[scalar_buffer_size*COMPs_FType[field_type]];
	D_REAL *dest_Field = dest_Blk->Field_Type[field_type];

	//! - receive buffer from responsible process
	mpi_myComm.Recv(receive_buffer,
			scalar_buffer_size*COMPs_FType[field_type],
			MPI_D_REAL,
			src_Blk->responsible_mpi_process,
			id_package);


	//! add buffer to field
	for(INT32 comp=0; comp<COMPs_FType[field_type]; comp++)
	 for(INT32 a=1; a<BlkNds_X; a++)
	 {
		XCS_pos  =  a    *BlkNds_Y *BlkNds_Z;
		XCSb_pos = (a-1) *(BlkNds_Y-1);

		for(INT32 b=1; b<BlkNds_Y; b++)
		dest_Field[	 comp*num_nodes_in_block +XCS_pos  +  b * BlkNds_Z   +dest_Plane] +=
		  receive_buffer[comp*scalar_buffer_size +XCSb_pos + (b-1)];

	 }

	delete[] receive_buffer;

}

