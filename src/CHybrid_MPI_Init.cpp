

#include "utils.h"
#include "CHybrid.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <math.h>

// #include <stdio.h>
#include <mpi.h>


//!--------------------------------------------------------
//!- init_MPI
//!  get MPI Globals
//!--------------------------------------------------------
void CHybrid::init_MPI(int argc, char *argv[])
{

	//! set com to global comm
	mpi_myComm = MPI::COMM_WORLD;

	//! be carefull not to confuse init of Hybrid and Init of MPI
	MPI::Init(argc,argv);
	mpi_myRank        = mpi_myComm.Get_rank();
	mpi_num_processes = mpi_myComm.Get_size();



	//! to store number of blocks at each process,
	//! will be initialized in init_Block or restore state
// 	num_rootBlocks_at_process = new INT32[mpi_num_processes];
// 	memset(num_rootBlocks_at_process, 0, mpi_num_processes *sizeof(INT32));

	
	char filename[200];
	sprintf(filename,"process_%d.log", mpi_myRank);

	//! use std out here !!!
	cout << "See "<<filename<<" for further information." << endl;

	if(!show_proc0_logfile_only || !mpi_myRank)
	log_file.open(filename, ios_base::app);

	log_file << endl << endl;
	log_file << "/******************************************/" << endl;
	log_file << "/*           A.I.K.E.F.mpi V4.2           */" << endl;
	log_file << "/******************************************/" << endl;
// 	log_file << endl;

	time_t rawtime;
	time ( &rawtime );
	log_file << " Start of Simulation: " <<  ctime (&rawtime);

	log_file << " System Clock uses " <<  CLOCKS_PER_SEC/1000 << "k clockticks/second." << endl <<endl;

	log_file << endl;
	log_file << "STARTING PROCESS "<< mpi_myRank << " LOG_FILE ..." << endl;
	log_file << "mpi_myRank:        " <<        mpi_myRank << endl;
	log_file << "mpi_num_processes: " << mpi_num_processes << endl;
	log_file << endl;


	log_file << "sizeof(MPI::Request): " <<  sizeof(MPI::Request) << endl;

	//! only attach mpi_buffer in case parallel run is started
	if(mpi_num_processes > 1)
	{
		//! attaching a buffer of resonable size doesn't seem to do any difference at HPCx
		//! estimate 6 face send for particle for send & receive
		INT32 buffer_size = (6*(BlkNds_X *BlkNds_Y+2)
					*RB_X*RB_Y*RB_Z)*optimal_MPiC[0] / mpi_num_processes;
	
// 		INT32 buffer_size =  (NUM_SCALAR_FIELDS_EACH_BLK * num_nodes_in_block * sizeof(D_REAL)
// 				   +optimal_MPiC[0] * num_nodes_in_block *sizeof(particle)) *(RB_X*RB_Y*RB_Z)/ mpi_num_processes;

// 		INT32 buffer_size = 0;


// 		log_file << "Attaching buffer to MPI of " 
// 			<< 1.*buffer_size *sizeof(particle)/1.e6 <<"MB." << endl;
// 		log_file << endl;
// 	
// 		particle* buffer = new particle[buffer_size];
// 		MPI::Attach_buffer(buffer, buffer_size);
	}

	log_file << "done." << endl;




}

	//-----------------------------------------
	//--- PARTICLE_MPI (CBlk_Globaals.h) ------
	//-----------------------------------------
// 	log_file << "Creating particle_mpi structure:" << endl;
// 
// 	int lengths[3]={3,3,1};
// 	Datatype types[3]={FLOAT,FLOAT,MPI_D_REAL};
// 	Aint displace[3];
// 	Aint float_size_mpi, double_size_mpi;
// 	
// 
// 	float_size_mpi = (Aint) FLOAT.Get_size();
// 	double_size_mpi= (Aint) MPI_D_REAL.Get_size();
// 
// 	displace[0] = (Aint) 0;
// 	displace[1] = 3*float_size_mpi;
// 	displace[2] = 6*float_size_mpi;
// 
// 
// 	PARTICLE_MPI=Datatype::Create_struct(3,lengths,displace,types);
// 	PARTICLE_MPI.Commit();

// 	log_file << "float_size_mpi: " << float_size_mpi << endl;
// 	log_file << "double_size_mpi: " << double_size_mpi << endl;
// 
// 	log_file << "displace[0]: " << displace[0] << endl;
// 	log_file << "displace[1]: " << displace[1] << endl;
// 	log_file << "displace[2]: " << displace[2] << endl;

//!--------------------------------------------------------
//!- init_MPI
//!  get MPI Globals
//!--------------------------------------------------------
void CHybrid::synchronize_allProcesses_MPI(void)
{

// 	log_file << " synchronizing ... " << endl;
	mpi_myComm.Barrier();

}

//!--------------------------------------------------------
//!- init_MPI
//!  get MPI Globals
//!--------------------------------------------------------
void CHybrid::show_information(INT64* local_values,
			       stringstream* info_names,
			       const INT32 num_values,
			       INT32 option)
{


	//! plot informations
	if(TL_LOGFILE_globalINFO  && TL%TL_LOGFILE_globalINFO==0)
	{

		INT64 global_values[num_values];

		if(option==0)
		mpi_myComm.Allreduce(local_values,
				global_values,
				num_values,
				MPI_INT64,
				MPI_SUM);
		else
		mpi_myComm.Allreduce(local_values,
				global_values,
				num_values,
				MPI_INT64,
				MPI_MAX);


		log_file << "  Info (global)" << endl;
		for(int i=0; i<num_values; i++)
		log_file <<  info_names[i].str() << global_values[i] << endl;


		memcpy(local_values, global_values, num_values *sizeof(INT64));


	}
	else
	{

		log_file << "  Info (local)" << endl;
		for(int i=0; i<num_values; i++)
		log_file <<  info_names[i].str() << local_values[i] << endl;

	}



}

//!--------------------------------------------------------
//!- init_MPI
//!  get MPI Globals
//!--------------------------------------------------------
void CHybrid::mpi_build_sum(D_REAL *local_values,
			    D_REAL *global_values,
			    INT32 num_values)
{

			mpi_myComm.Allreduce(local_values,
					     global_values,
					     num_values,
					     MPI_DOUBLE,
					     MPI_SUM);

}



//!--------------------------------------------------------
//!- check_for_NaN_MPI
//!--------------------------------------------------------
void CHybrid::check_for_NaN_MPI(bool is_NaN_this_process)
{

	//! code will crash in case in and out arrays are identical
	int *NaN_in_process_in = new int[mpi_num_processes];
	memset(NaN_in_process_in, 0, mpi_num_processes*sizeof(int));

	int *NaN_in_process_out = new int[mpi_num_processes];
	memset(NaN_in_process_out, 0, mpi_num_processes*sizeof(int));


	//! set 0/1 of this processes in array
	//! (all others are zero)
	if(is_NaN_this_process)
	NaN_in_process_in[mpi_myRank] = 1;

	//! sum up all
	//! -> if 1. at any element, respective process at NaN
	mpi_myComm.Allreduce(NaN_in_process_in,
			     NaN_in_process_out,
			     mpi_num_processes,
			     MPI_INT,
			     MPI_SUM);

	bool NaN_at_any_process = false;
	for(INT32 process=0; process<mpi_num_processes; process++)
	{
		if(NaN_in_process_out[process])
		{
			log_file << "  NaN at process " << process << " (see respective process.log for further information)"<< endl;
			NaN_at_any_process = true;
		}
	
	}

	if(NaN_at_any_process)
	{

		log_file << "  NaN at at least one process"  << endl;
		log_file << "  Exiting ..."  << endl;

		finalize_MPI();
		exit(1);
		
	}

	delete[] NaN_in_process_in;
	delete[] NaN_in_process_out;

}



//!--------------------------------------------------------
//!- init_MPI
//!  get MPI Globals
//!--------------------------------------------------------
void CHybrid::finalize_MPI(void)
{

// 	log_file << " synchronizing ... " << endl;
	MPI::Finalize();

}



