



#include "CHybrid.h"

#include "utils.h"
#include "parameters.h"
#include "output_gnuplot.h"
#include "output_silo.h"
#include "output_uniform_grid.h"
#include "absolute_Globals.h"
// #include "CBlk_Globals.h"
#include "unistd.h"

#define TIME_FIELD	0
#define TIME_PART	1
#define TIME_TOTAL	2

#define ERROR_SMALL_STATEFILE  "Missmatch in blocksize of state file and active run \n"\
			       "Statefile to SMALL\n"\
			       "Exiting ...."

#define ERROR_LARGE_STATEFILE  "Missmatch in blocksize of state file and active run \n"\
			       "Statefile to LARGE\n"\
			       "Exiting ...."


#define TIME_TO_WAIT 200000

#include <iostream>
#include <fstream>
#include <sstream>


INT64 num_Blks_in_CrossSection[3];

ofstream Particle_File;
ifstream Particle_inFile;




//!------------------------------------------------------------
//!- output_all:
//!------------------------------------------------------------
void CHybrid::output_all(INT32 TL)
{



	//! average fields if required
	handle_average_fields(TL);
        

// 	sum_up_block_timing();

// 	if(TL_OUTPUT_2D_SILO && TL%TL_OUTPUT_2D_SILO==0)
// 	collect_RHOnp1_UIplus_LAM_GAM();


// 	if(TL_OUTPUT_2D_NATIVE && TL%TL_OUTPUT_2D_NATIVE==0)
// 	collect_RHOnp1_UIplus_LAM_GAM();



	//! ---- marked particles ------------------------------
	//! mark before output so they are available in first TL
	if(TL_MARK_PARTICLE_TRACKS==TL)
	mark_particle();

	//! -------------------------------------
	//! -- OUTPUT ---------------------------
	//! -------------------------------------

	//! ---- Silo ---------------------------
	if(TL_OUTPUT_2D_SILO && TL%TL_OUTPUT_2D_SILO==0)
	silo_output_2D();

	if(TL_OUTPUT_3D_SILO && TL%TL_OUTPUT_3D_SILO==0)
	silo_output_3D();


	//! ---- Native -------------------------
	if(TL_OUTPUT_2D_NATIVE && TL%TL_OUTPUT_2D_NATIVE==0)
	native_output_2D(TL);

	if(TL_OUTPUT_3D_NATIVE && TL%TL_OUTPUT_3D_NATIVE==0)
	native_output_3D(TL);


	//! ---- Gnuplot ------------------------
	if(TL_OUTPUT_2D_GNUPLOT && TL%TL_OUTPUT_2D_GNUPLOT==0)
	gnuplot_output(TL);


	//! ---- ErgMomMass ---------------------
	if(TL_OUTPUT_ErgMomMass && TL%TL_OUTPUT_ErgMomMass==0)
	collect_Energy_Momentum_Mass();


	//! ---- Trajectory ---------------------
	if(TL_OUTPUT_TRAJECTORY && TL%TL_OUTPUT_TRAJECTORY==0)
	{
		trajectory_output();
	}
	
	//! ---- LineOut ------------------------
	if(TL_OUTPUT_LINEOUT && TL%TL_OUTPUT_LINEOUT==0)
	{

		for(INT32 id_line=0; id_line<num_lines; id_line++)
		lineout_output(id_line);
	}
	
	//! ---- GetMesh --------------------
	if(TL_OUTPUT_GETMESH && TL%TL_OUTPUT_GETMESH==0)
        {
                get_Mesh();
        }
	
	//! ---- uniform_grid ------------------------
	if(TL_OUTPUT_3D_uniform_grid && TL%TL_OUTPUT_3D_uniform_grid==0)
	uniform_grid_output(TL);

	//! ---- ASCII ------------------------
	if(TL_OUTPUT_3D_ASCII && TL%TL_OUTPUT_3D_ASCII==0)
	ascii_output_3D(TL);


	
	if(TL_OUTPUT_PARTICLE_TRACKS && TL%TL_OUTPUT_PARTICLE_TRACKS==0)
	trace_particle_output();

	//! ---- particleDetector ------------------------
	if(TL_OUTPUT_PARTICLE_DETECTOR && TL%TL_OUTPUT_PARTICLE_DETECTOR==0)
	{
	  for(INT32 id_detector=0; id_detector<num_particle_detector; id_detector++)
	    particle_detector_output(id_detector);
	}

// 	if(TL_OUTPUT_PARTICLE_TRACKS && TL==TL_MAX)
// 	parallel_particle_track_to_single_file();

	if(TL_CONVERT_PARTICLE_TRACKS_TO_SILO_MESH && TL%TL_CONVERT_PARTICLE_TRACKS_TO_SILO_MESH==0 )
	{
		synchronize_allProcesses_MPI();
		convert_particle_tracks_to_silo_mesh();
	}

}

//!--------------------------------------------------------
//! handle_average_fields
//!--------------------------------------------------------
void CHybrid::handle_average_fields(INT32 TL)
{

	if(!TL_FINISH_AVERAGE_FIELDS || !num_average_fields)
	return;

	clock_t start,finish;
	double time;
	

	//! set average fields to zero if required
	if((TL +num_TL_average_Fields) %TL_FINISH_AVERAGE_FIELDS == 0)
	{

		log_file << " RESETTING AVERAGE FIELDS ... " << endl;
		reset_average_fields();
		log_file << " done. " << endl;

	}
	

	//! add fields to average field array if required
	if((TL +num_TL_average_Fields) %TL_FINISH_AVERAGE_FIELDS < num_TL_average_Fields)
	{
		start = clock();
		log_file << " ADDING FIELDS TO AVERAGE FIELD ARRAY ... " << endl;
		
		add_fields_to_average_fields();
		
		finish = clock();
		time = double(finish-start)/CLOCKS_PER_SEC;
		log_file << " done.\n";
		log_file << " adding fields to avrg fields time : " << time << "s." << endl;
	
	}


	//! devide field by number of average TL
	if((TL+1) %TL_FINISH_AVERAGE_FIELDS == 0)
	{
		start = clock();
		log_file << " NORMALIZING AVERAGE FIELDS ... " << endl;
		
		normalize_average_fields();
		
		finish = clock();
		time = double(finish-start)/CLOCKS_PER_SEC;
		log_file << " done.\n";
		log_file << " normalizing avrg fields time : " << time << "s." << endl;
	}


}

//!--------------------------------------------------------
//! reset_average_fields
//!--------------------------------------------------------
void CHybrid::reset_average_fields(void)
{

	//!--------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->reset_average_fields();
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}


//!--------------------------------------------------------
//! average_fields
//!--------------------------------------------------------
void CHybrid::add_fields_to_average_fields(void)
{


	bool rotB_in_use = false;

	//! calc vTherm
	//! calculate vt2
	for(INT32 average_field=0; average_field<num_average_fields; average_field++)
	 for(int species=0; species<num_Charged_Species; species++)
	 if(IDs_of_Fields_to_average[average_field]==id_PISpecies1)
	 {

		rotB_in_use = true;

		//! provide rho
		//! gather Ui 
		collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);

		//! provide rho
		//! provide Ui
		//! gather vth2
		collect_Ui_vth2_of_species(species, id_UI_Species1, id_PISpecies1 ,getVTH2);

	 }


	for(INT32 average_field=0; average_field<num_average_fields; average_field++)
	if(IDs_of_Fields_to_average[average_field]==id_rotB)
	{


		if(mesh_type == UNIFORM || rotB_in_use==false)
		calc_rot(id_BEven, id_rotB);
		else
		{
			log_file << " cannot average curlB" << endl;
			log_file << " no unused field available." << endl;

		}

	}




	//!--------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->add_fields_to_average_fields();
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}

//!--------------------------------------------------------
//! average_fields
//!--------------------------------------------------------
void CHybrid::normalize_average_fields(void)
{

	//!--------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->normalize_average_fields();
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

}


//!--------------------------------------------------------
//! Trace Particles
//!--------------------------------------------------------
void CHybrid::mark_particle(void)
{
	clock_t start,finish;
	double time;
	start = clock();
	
	log_file << " MARKING PARTICLES ... " << endl;


	bool write_particle_number_into_file = true;


	//! ----------------------------------------
	//! In case particles alredy are enumerated
	//! and only a given set of particles shall
	//! be traced, read list of numbers which
	//! to trace
	//! ----------------------------------------

	char filename[200];
	INT64 num_partcle_in_list=0;
	INT64* Partilce_number_List;
	sprintf(filename,"particle_numbers.txt");
	ifstream ParticleNummber_FILE;
	ParticleNummber_FILE.open(filename);

	if(ParticleNummber_FILE)
	{

		INT64 dummy;
		INT64 num_partcle_in_list = 0;
		for(;;)
		{
	
			ParticleNummber_FILE >> dummy;
	
			if(ParticleNummber_FILE.eof())
			break;
	
	
			num_partcle_in_list++;
	
	
		}
		ParticleNummber_FILE.close();
	
		if(num_partcle_in_list)
		Partilce_number_List = new INT64[num_partcle_in_list];
	
	
		num_partcle_in_list = 0;
		ParticleNummber_FILE.open(filename, ios_base::app);
		for(;;)
		{
	
			ParticleNummber_FILE >> Partilce_number_List[num_partcle_in_list];
	
			if(ParticleNummber_FILE.eof())
			break;
	
			num_partcle_in_list++;
	
	
		}
		ParticleNummber_FILE.close();
	}
	else
	log_file << " No Particle List File Found " << endl;
	//! ----------------------------------------



	
	//! now begin marking the particles
	INT64* particle_counter_of_species = new INT64[num_Charged_Species];
	memset(particle_counter_of_species, 0, num_Charged_Species*sizeof(INT64));


	//! In order to avoid that all processes write at the same time,
	//! let p0 write first, than send message to p1 to trigger writing,
	//! then p2 .... up to pN
	INT32 tag = 0;
	INT32 is_finished = 0;
	INT32 num_entries = num_Charged_Species;
	bool serialize_marking_of_Particles = true;

	num_marked_particle = 0;

	//! - every process must receive except for first process
	//! - use Blocking receive


	//! HAS TO BE SERIALIZED AND NUMBER MUST BE SENT, OTHERWISE NEXT 
	//! PROCESS DOES NOT KNOW WHERE TO CONTINUE COUNTING !!!
	if(serialize_marking_of_Particles && mpi_myRank)
	mpi_myComm.Recv(particle_counter_of_species, num_entries, MPI_INT64, mpi_myRank-1, tag );

	//! mark particles
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->mark_particle(particle_counter_of_species, num_partcle_in_list, Partilce_number_List);


// 			if(write_particle_number_into_file)
// 			if(mpi_myRank == temp_Block->responsible_mpi_process)
// 			temp_Block->write_particle_number();
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	//! - every process must send except for last process
	//! - use Blocking Send
	if(serialize_marking_of_Particles && mpi_myRank < mpi_num_processes-1)
	{	
		//! wait time in micro seconds
		usleep(TIME_TO_WAIT);
		mpi_myComm.Send(particle_counter_of_species, num_entries, MPI_INT64, mpi_myRank+1, tag);

	}


	//! Show information:
	INT64* local_info_values = new INT64[num_Charged_Species];
	stringstream info_names[INFO_ARRAY_SIZE];

	log_file << "  Marked particle this process: " << endl;
	for(INT32 species=0; species<num_Charged_Species; species++)
	{

		info_names[species] << "   ->species " << species << ": ";
		local_info_values[species] = particle_counter_of_species[species];

		log_file << info_names[species].str() <<  local_info_values[species] << endl;
	}


	delete[] particle_counter_of_species;


	if(num_partcle_in_list)
	delete[] Partilce_number_List;


	show_information(local_info_values,
			 info_names,
			 num_Charged_Species,
			 BUILD_MAX);




	//! store number of marked paricles
	for(INT32 species=0; species<num_Charged_Species; species++)
	num_marked_particle += local_info_values[species];





	//! measure time
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " marking particles time: " << time << "s." << endl << endl;
}






//!--------------------------------------------------------
//! Trace Particles
//!--------------------------------------------------------
void CHybrid::trace_particle_output(void)
{
  //  if(num_mark_particle_in_species(species) > 0
  
	clock_t start,finish;
	double time;
	start = clock();
	
	log_file << " TRACING PARTICLES ... "<<endl;


	INT64* particle_counter_of_species = new INT64[num_Charged_Species];
	memset(particle_counter_of_species, 0, num_Charged_Species*sizeof(INT64));


	//! In order to avoid that all processes write at the same time,
	//! let p0 write first, than send message to p1 to trigger writing,
	//! then p2 .... up to pN
	INT32 tag = 0;
	INT32 is_finished = 0;
	INT32 length_in_byte = 1;

// 	bool serialize_writing_of_ParticleTracks = true;

	//! - every process must receive except for first process
	//! - use Blocking receive
// 	if(serialize_writing_of_ParticleTracks && mpi_myRank)
// 	mpi_myComm.Recv(&is_finished, length_in_byte, MPI_INT32, mpi_myRank-1, tag );

	//!--------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->trace_particle(particle_counter_of_species);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}

	is_finished = 1;

	//! - every process must send except for last process
	//! - use Blocking Send
// 	if(serialize_writing_of_ParticleTracks && mpi_myRank < mpi_num_processes-1)
// 	{	
// 		//! wait time in micro seconds
// 		usleep(TIME_TO_WAIT);
// 		mpi_myComm.Send(&is_finished, length_in_byte, MPI_INT32, mpi_myRank+1, tag);
// 
// 	}





	//! Show information:
	INT64* local_info_values = new INT64[num_Charged_Species];
	stringstream info_names[INFO_ARRAY_SIZE];

	log_file << "  Traced particle this process: " << endl;
	for(INT32 species=0; species<num_Charged_Species; species++)
	{



		info_names[species] << "   ->species " << species << ": ";
		local_info_values[species] = particle_counter_of_species[species];

		log_file << info_names[species].str() <<  local_info_values[species] << endl;



	}


	show_information(local_info_values,
			 info_names,
			 num_Charged_Species,
			 BUILD_SUM);



	 num_marked_particle = 0;

	//! store number of marked paricles
	for(INT32 species=0; species<num_Charged_Species; species++)
	num_marked_particle += local_info_values[species];



	delete[] local_info_values;
	delete[] particle_counter_of_species;
	//! measure time
	finish = clock();
	time = (double(finish)-double(start))/CLOCKS_PER_SEC;
	log_file << " tracing particle time: " << time << "s." << endl << endl;
}


//!--------------------------------------------------------
//! parallel_particle_track_to_single_file
//!--------------------------------------------------------
void CHybrid::parallel_particle_track_to_single_file(void)
{

	log_file << " Packing particle data to single files ... " << endl;
	log_file << " ---------------------------------------------------------------------------------------------------- " << endl;
	log_file << " ";



	//! assume particles to be in s0
	INT32 species = 0;

	char buffer[300];
	char Particle_outFileName[300];
	char Particle_inFileName[300];

	INT64 time_level = 0;
	PARTICLE_REAL x[3], v[4];

	//! only master shall pack files
	if(!mpi_myRank)
	//! each ion has its on file
	for (INT32 particle=0; particle<obstacle_MP_num_each_TL[species]; particle++)
	{


		sprintf(Particle_outFileName,"%s/particle_tracks/spec%d_part%07d.txt",
		data_output_path, species, particle);

		//! open file in which particle data shall
		//! be written (of all processes)
		ofstream Particle_outFile;
		Particle_outFile.open(Particle_outFileName);


		for (INT32 proc=0; proc<mpi_num_processes; proc++)
		{

			sprintf(Particle_inFileName,"%s/particle_tracks/spec%d_part%07d_p%05d.txt",
			data_output_path, species, particle, proc);

			//! open file from which particle data shall
			//! be read (one specific process)
			ifstream Particle_inFile;
			Particle_inFile.open(Particle_inFileName);

			//! process only writes Particle_inFile in case particle has been
			//! at process, so in general not every proces will rite file
			if(Particle_inFile)
			{

				for(;;)
				{
			
					//! get time level
					Particle_inFile >> time_level;

			
					if(Particle_inFile.eof())
					break;

					//! get particle position
					Particle_inFile >> x[0];
					Particle_inFile >> x[1];
					Particle_inFile >> x[2];

					//! get particle velocitiy
// 					Particle_inFile >> v[0];
// 					Particle_inFile >> v[1];
// 					Particle_inFile >> v[2];


					Particle_outFile << time_level << "	";

					Particle_outFile << x[0] << "	";
					Particle_outFile << x[1] << "	";
					Particle_outFile << x[2] << endl;

					
				}//! end of for all positions in file

				Particle_outFile.flush();

				Particle_inFile.close();
		
				sprintf(buffer,"rm %s", Particle_inFileName);
				system(buffer);

			}//! end if file


		}//! end for all processes


		//! close file for this particle
		Particle_outFile.close();


		if(particle % 100 ==0 )
		log_file << "-";
		

	}//! end for all particle



	log_file << " done. "  << endl;

	synchronize_allProcesses_MPI();

}





//!------------------------------------------------------------
//!- trajectory_count_values:
//!
//!------------------------------------------------------------
void CHybrid::trajectory_count_values(INT32 id_trajectory, INT32 &num_values, INT32 &num_positions)
{


	char buffer1[200];
	ifstream Trajectory_File;
	Trajectory_File.open(Trajectory_FileName[id_trajectory]);

	num_values=0;
	num_positions=0;


	if(!Trajectory_File)
	{
		log_file << " no "<<Trajectory_FileName[id_trajectory]<<" File found." << endl;
		return;
	}
	
	
	log_file << endl << " Reading trajectory of file "<<Trajectory_FileName[id_trajectory]<<" ..." << endl;


	while(!Trajectory_File.eof())
	{
		Trajectory_File >> buffer1;
		num_values++;

// 		if(num_values%10000==0)
// 		log_file << num_values << endl;
// 		log_file << buffer << endl;
	}

	Trajectory_File.close();

	if(do_read_timestring && !do_read_SpaceCraftField)
	{
		num_positions = num_values/4;
		log_file <<"  ->Values: "<< num_values<<"."<< endl;
		log_file <<"  ->Positions: "<< num_positions<<"."<< endl;

	}
	if(do_read_timestring && do_read_SpaceCraftField)
	{
		num_positions = num_values/7;
		log_file <<"  ->Values: "<< num_values<<"."<< endl;
		log_file <<"  ->Positions: "<< num_positions<<"."<< endl;

	}
	if(!do_read_timestring && !do_read_SpaceCraftField)
	{
		num_positions = num_values/3;
		log_file <<"  ->Values: "<< num_values<<"."<< endl;
		log_file <<"  ->Positions: "<< num_positions<<"."<< endl;

	}

}


//!------------------------------------------------------------
//!- trajectory_count_values:
//!
//!------------------------------------------------------------
void CHybrid::read_trajectory(INT32 id_trajectory, INT32 &num_positions, D_REAL** &positions, char** &time_string_array, D_REAL** &SCField)
{



	ifstream Trajectory_File;
	INT32 num_pos_out_of_box=0;


	//! alloc time string memory
	if(do_read_timestring)
	{
		time_string_array = new char*[num_positions];
		for(INT32 pos=0; pos<num_positions; pos++)
		time_string_array[pos] = new char[30];
	}


	positions = new D_REAL*[3];

	//! alloc coordinates memory
	positions[0] = new D_REAL[num_positions];
	positions[1] = new D_REAL[num_positions];
	positions[2] = new D_REAL[num_positions];


	//! alloc SpaceCraft field memory
	if(do_read_SpaceCraftField)
	{
		SCField = new D_REAL*[3];

		SCField[0] = new D_REAL[num_positions];
		SCField[1] = new D_REAL[num_positions];
		SCField[2] = new D_REAL[num_positions];

	}
	
	Trajectory_File.open(Trajectory_FileName[id_trajectory]);
	for(INT32 counter=0; counter<num_positions; counter++)
	{


		if(do_read_timestring)
		Trajectory_File >> time_string_array[counter];
		

		//! read coordinates
		Trajectory_File >> positions[0][counter];
		Trajectory_File >> positions[1][counter];
		Trajectory_File >> positions[2][counter];
		
		
		if(do_read_SpaceCraftField)
		{
			Trajectory_File >> SCField[0][counter];
			Trajectory_File >> SCField[1][counter];
			Trajectory_File >> SCField[2][counter];
		}


		//! conversion km -> normalized units
		positions[0][counter] = positions[0][counter]/trajectory_pos_conversion_fact;
		positions[1][counter] = positions[1][counter]/trajectory_pos_conversion_fact;
		positions[2][counter] = positions[2][counter]/trajectory_pos_conversion_fact;

		//! normalize coordinate to box length
		//! get normalized position in Simu-Box
		positions[0][counter] = (positions[0][counter]+Box_Origin[0])/LX;
		positions[1][counter] = (positions[1][counter]+Box_Origin[1])/LY;
		positions[2][counter] = (positions[2][counter]+Box_Origin[2])/LZ;

		if(    positions[0][counter]<0. || positions[0][counter] >= 1.
		    || positions[1][counter]<0. || positions[1][counter] >= 1.
		    || positions[2][counter]<0. || positions[2][counter] >= 1.)
		{


// 			log_file << " requested position "<<counter<<" out of box:"<< endl
// 			         << " ignoring x=("<<posX[counter]<<","
// 					       <<posY[counter]<<","
// 					       <<posZ[counter]<<") in units of Box Length." << endl;

			num_positions--;
			counter--;
			num_pos_out_of_box++;
		}

	}


	log_file << "  ->valid positions:       " << num_positions      << endl;
	log_file << "  ->positions outside box: " << num_pos_out_of_box << endl;
	Trajectory_File.close();


}



//!------------------------------------------------------------
//!- trajectory_output:
//!
//!------------------------------------------------------------
void CHybrid::trajectory_output(void)
{

	//! alloc array for trace density or only Bfield...
	bool** trace_field_flyby;
	trace_field_flyby = new bool*[num_trajectories];
	for(INT32 id_trajectory=0; id_trajectory<num_trajectories; id_trajectory++)
	trace_field_flyby[id_trajectory] = new bool[num_fields_to_trace];

	for(INT32 flyby=0; flyby<num_trajectories; flyby++)
	  for(INT32 field=0; field<num_fields_to_trace; field++)
	  {	
		trace_field_flyby[flyby][field]=true;
		
		//! only write neutral density at TL 0
// 		if(TL)
// 		{
// 			//! BField for all
// 			if(ID_of_Fields_to_trace[field]==id_BTotal 
// 				|| ID_of_Fields_to_trace[field]==id_BEven
// 				|| ID_of_Fields_to_trace[field]==id_average_Field1 +3)
// 			trace_field_flyby[flyby][field]=true;
// 			
// 			//! Density for E2, E3, E5, E7, E17, E15
// 			if(flyby==2 || flyby==3 || flyby==5 || flyby==7 || flyby==15 || flyby==17)
// 			trace_field_flyby[flyby][field]=true;
// 		}
// 		else
// 		{
// 			if(ID_of_Fields_to_trace[field]==id_externRho1 && (flyby==2 || flyby==3 || flyby==5 || flyby==7) )
// 			trace_field_flyby[flyby][field]=true;	
// 		}	
	  }


	for(INT32 id_trajectory=0; id_trajectory<num_trajectories; id_trajectory++)
	{	
		if(trajectory_2D_out_only[id_trajectory])
		return;
	

		D_REAL** positions;
		char** time_string_array;
		D_REAL** SCField;
		INT32 num_values, num_positions;

		
		//!-----------------------------------
		//! 1) count number of values in file
		//!-----------------------------------
		trajectory_count_values(id_trajectory, num_values, num_positions);


		//!-----------------------------------
		//! 2) read values in file
		//!-----------------------------------
		read_trajectory(id_trajectory, num_positions, positions, time_string_array, SCField);


		//!-----------------------------------
		//! 3) estimate field data
		//!-----------------------------------
		char buffer2[200];
		sprintf(buffer2,"mkdir %s/trajectories/temp", data_output_path);

		if(!mpi_myRank)
		system(buffer2);
		//! wait until directory is created
		synchronize_allProcesses_MPI();


		if(!num_positions)
		log_file << " no valid positions." << endl;
		else
		{	
			
			bool  calc_field = false;
		
			for(INT32 field=0; field<num_fields_to_trace; field++)
			if(trace_field_flyby[id_trajectory][field])
 			{
			
				
				
	
				if(ID_of_Fields_to_trace[field]!=id_PMagnetic &&  ID_of_Fields_to_trace[field]!=id_ExB)
				{
					if(collect_only_ion_species)
					{
						collect_only_ion_species=false;
						collect_RHOnp1_UIplus_LAM_GAM();
					}	
				}
					
				//! calculate current
				if(ID_of_Fields_to_trace[field]==id_rotB)
				calc_rot(id_BEven, id_rotB);
		
				//! calculate divergence
				if(ID_of_Fields_to_trace[field]==id_divB)
				calc_div(id_BEven, id_divB);
				
				//! provide rho
				//! gather Ui 
				for(INT32 species=0; species<num_Charged_Species; species++)
				if(ID_of_Fields_to_trace[field]==id_UI_Species1 +species)
				collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);
		
		
				//! calculate vt2
				for(int species=0; species<num_Charged_Species; species++)
				if(ID_of_Fields_to_trace[field]==id_PISpecies1 +species)
				{
		
					//! provide rho
					//! gather Ui 
					collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);
		
					//! provide rho
					//! provide Ui
					//! gather vth2
					collect_Ui_vth2_of_species(species, id_UI_Species1, id_PISpecies1 ,getVTH2);
		
				}
						
				
				//! total ion density
				if(ID_of_Fields_to_trace[field]==id_PMagnetic && !calc_field)
				{
					collect_only_ion_species=true;
					collect_RHOnp1_UIplus_LAM_GAM();
					copy_Field(id_PMagnetic,id_rho_np1);
					calc_field = true;
				}
				
				//! total ion velocity ????
				if(ID_of_Fields_to_trace[field]==id_ExB)
				{
					if(!calc_field)
					{
						collect_only_ion_species=true;
						collect_RHOnp1_UIplus_LAM_GAM();
						copy_Field(id_ExB,id_UI_plus);
						calc_field = true;
					}
					if(calc_field)
					copy_Field(id_ExB,id_UI_plus);	
				}
				
				
				if(ID_of_Fields_to_trace[field]==ENERGY_SPECTRUM)
				trace_trajectory_energy_spectrum(id_trajectory, ID_of_Fields_to_trace[field], num_positions, positions, time_string_array, SCField);
				else
				trace_trajectory(id_trajectory, ID_of_Fields_to_trace[field], num_positions, positions, time_string_array, SCField);
		
				
				if(pack_parallel_output_to_single_file)
				{
		
					//! All processors must have written files before master can
					//! begin packing
					synchronize_allProcesses_MPI();
		
					if(!mpi_myRank)
					{
						if(ID_of_Fields_to_trace[field]==ENERGY_SPECTRUM)
						parallel_output_to_single_file_energy_spectrum(id_trajectory, ID_of_Fields_to_trace[field]);
						else
						parallel_output_to_single_file(id_trajectory, ID_of_Fields_to_trace[field]);
					}
				}
		
			}




		}
		
		if(collect_only_ion_species)
		{
			collect_only_ion_species=false;
			collect_RHOnp1_UIplus_LAM_GAM();
		}	

		//!-----------------------------------
		//! 4) clean up
		//!-----------------------------------
		//! delete time_string_array memory	
		if(do_read_timestring)
		{
			for(INT32 pos=0; pos<num_positions; pos++)
			delete[] time_string_array[pos];

			delete[] time_string_array;
		}


		//! delete coordinates memory
		delete[] positions[0];
		delete[] positions[1];
		delete[] positions[2];

		delete[] positions;
		
		if(do_read_SpaceCraftField)
		{
			delete[] SCField[0];
			delete[] SCField[1];
			delete[] SCField[2];

			delete[] SCField;
		}

		sprintf(buffer2,"rm -r %s/trajectories/temp", data_output_path);

		if(!mpi_myRank)
		system(buffer2);
	}

	for(INT32 id_trajectory=0; id_trajectory<num_trajectories; id_trajectory++)
	delete[] trace_field_flyby[id_trajectory];

	delete[] trace_field_flyby;

}


//!------------------------------------------------------------
//!trace_trajectory:
//!
//! - Field_id: id of Field that shall be traced along trajectory
//! - positions: x,y,z = [0;1[ normalized position in Simu-Box
//!------------------------------------------------------------
void CHybrid::parallel_output_to_single_file(INT32 id_trajectory, INT32 id_field)
{




	ifstream Trajectory_File;
	Trajectory_File.open(Trajectory_FileName[id_trajectory]);


	if(!Trajectory_File)
	{
		log_file << " no " << Trajectory_FileName[id_trajectory] << " File found." << endl;
		return;
	}

	Trajectory_File.close();


	log_file << "     Packing files of field '" << Field_Name[id_field] << "' to single file... ";


	char buffer[200];
	double Field[3], r[3], SC_Field[3];;
	char outfile_name[200];
	char time_string[200];

	ofstream outfile;


	//! decide whether to attach data to file or open new one
	if(create_new_file_each_TL)
	{
		sprintf(outfile_name,"%s/trajectories/%s_%s_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[id_field], TL);
		outfile.open(outfile_name);


	}
	else
	{

		sprintf(outfile_name,"%s/trajectories/%s_%s.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[id_field]);

		if(!TL)
		outfile.open(outfile_name);
		else
		outfile.open(outfile_name, ios_base::app);


	}






	int counter_total = 0;
	for(int process=0; process<mpi_num_processes; process++)
	{


		char infile_process[200];

		if(create_new_file_each_TL)
		sprintf(infile_process,"%s/trajectories/temp/%s_%s_p%05d_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[id_field], process, TL);
		else
		sprintf(infile_process,"%s/trajectories/temp/%s_%s_p%05d.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[id_field], process);


		ifstream infile;
		infile.open(infile_process);


		if(!infile)
		{
			log_file << "Error in 'parallel_output_to_single_file':" << endl;
			log_file << "no file: " << infile_process << endl;
			log_file << "Returning ... ." << endl;
			return;
			

		}

	

		int counter_file = 0;
		for(;;)
		{


			if(do_read_timestring)
			infile >> time_string;

			//! get position
			infile >> r[0];
			infile >> r[1];
			infile >> r[2];
			
			if(do_read_SpaceCraftField)
			{
				infile >> SC_Field[0];
				infile >> SC_Field[1];
				infile >> SC_Field[2];
			}

			if(infile.eof())
			break;
	

			if(do_read_timestring)
			outfile << time_string <<"	";

			r[0]*=factor_scale_mesh;
			r[1]*=factor_scale_mesh;
			r[2]*=factor_scale_mesh;

			//! write outfile
			outfile 	<< r[0] <<"		";
			outfile 	<< r[1] <<"		";
			outfile 	<< r[2] <<"		";
			outfile <<  vec_len(r)  <<"		";

			if(do_read_SpaceCraftField)
			{
				outfile 	<< SC_Field[0] <<"		";
				outfile 	<< SC_Field[1] <<"		";
				outfile 	<< SC_Field[2] <<"		";
// 				outfile << sqrt(SC_Field[0]*SC_Field[0] +SC_Field[1]*SC_Field[1] +SC_Field[2]*SC_Field[2]) <<"   	";
			}
			
			
			if(COMPs_FType[id_field]==1)
			{
				//! get Field
				infile >> Field[0];
		
				//! write BField
				outfile << Field[0] <<"   	";
				outfile << endl;
			}

			if(COMPs_FType[id_field]==3)
			{
				//! get Field
				infile >> Field[0];
				infile >> Field[1];
				infile >> Field[2];
		
				//! write BField
				outfile << Field[0] <<"   	";
				outfile << Field[1] <<"   	";
				outfile << Field[2] <<"   	";
				outfile << sqrt(Field[0]*Field[0] +Field[1]*Field[1] +Field[2]*Field[2]) <<"   	";
				outfile << endl;
			}
	
			counter_file++;

		}
		counter_total +=counter_file;

		infile.close();


// 		sprintf(buffer,"rm %s", infile_process);
// 		system(buffer);

// 		log_file << "  " << counter_file << " values" << endl;


	}


	log_file << " done. "  << endl;
	log_file << "     ->Filename: '" << outfile_name <<"'"<< endl;
	log_file << "     ->Values read/written: " << counter_total << endl;

	log_file << endl;


	outfile.close();




}


//!------------------------------------------------------------
//!trace_trajectory:
//!
//! - Field_id: id of Field that shall be traced along trajectory
//! - positions: x,y,z = [0;1[ normalized position in Simu-Box
//!------------------------------------------------------------
void CHybrid::trace_trajectory(INT32 id_trajectory, INT32 Field_id, INT32 num_positions, D_REAL** positions, char** time_string_array, D_REAL** SCField)
{


	log_file <<"     Tracing "<<Field_Name[Field_id]<<" ...  ";

	//! indices for interpolation
	INT32  i,j,k;
	INT32  i_j_k, ip1_jp1_kp1;
	INT32  ip1_j_k, i_jp1_k,  i_j_kp1;
	INT32  ip1_jp1_k, ip1_j_kp1,  i_jp1_kp1;

	ofstream trajectory_outfile;
	INT32  *index_blk[3], *index_cell[3], blk_index;
	D_REAL *x_cell[3], r_vec[3], A[3], shape_func[8], *BlkField[3];

	CBlock *temp_Block, *top_level_Block;



	//! decide whether to attach data to file or open new one
	if(create_new_file_each_TL)
	{

		char filename[200];

		if(pack_parallel_output_to_single_file)
		sprintf(filename,"%s/trajectories/temp/%s_%s_p%05d_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[Field_id], mpi_myRank, TL);
		else
		sprintf(filename,"%s/trajectories/%s_%s_p%05d_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[Field_id], mpi_myRank, TL);


		trajectory_outfile.open(filename);
	}
	else
	{

		char filename[200];

		if(pack_parallel_output_to_single_file)
		sprintf(filename,"%s/trajectories/temp/%s_%s_p%05d.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[Field_id], mpi_myRank);
		else
		sprintf(filename,"%s/trajectories/%s_%s_p%05d.txt", data_output_path, Trajectory_FileName[id_trajectory], Field_Name[Field_id], mpi_myRank);

		if(!TL)
		trajectory_outfile.open(filename);
		else
		trajectory_outfile.open(filename, ios_base::app);

// 		trajectory_outfile << TL*dt << "	";

	}

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

	for(INT32 counter=0; counter<num_positions; counter++)
	{

	
		//!--------------------------------------------------
		//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
		//!     (even though target blk might not be refined)
		//!--------------------------------------------------
		r_vec[0] = positions[0][counter];
		r_vec[1] = positions[1][counter];
		r_vec[2] = positions[2][counter];


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
				for(INT32 comp=0; comp<COMPs_FType[Field_id]; comp++)
				BlkField[comp] = top_level_Block->Field_Type[Field_id] +comp *num_nodes_in_block;
				
		
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
	


				if(do_read_timestring)
				trajectory_outfile << time_string_array[counter] <<"	";
	
				trajectory_outfile << positions[0][counter] *LX -Box_Origin[0] <<"		"
						   << positions[1][counter] *LY -Box_Origin[1] <<"		"
						   << positions[2][counter] *LZ -Box_Origin[2] <<"		";
						   
				if(do_read_SpaceCraftField)
				{
					trajectory_outfile 	<< SCField[0][counter] <<"		"
								<< SCField[1][counter] <<"		"
								<< SCField[2][counter] <<"		";
				}

		
				//! interpolate and store field in variable A
				for(INT32 comp=0; comp<COMPs_FType[Field_id]; comp++)
// 				for(INT32 comp=0; comp<1; comp++)
				{
		
				      A[comp] = BlkField[comp][  i_j_k] * shape_func[0]
						+BlkField[comp][ip1_j_k] * shape_func[1]
						+BlkField[comp][i_jp1_k] * shape_func[2]
						+BlkField[comp][i_j_kp1] * shape_func[3]
				
						+BlkField[comp][  ip1_jp1_k] * shape_func[4]
						+BlkField[comp][  ip1_j_kp1] * shape_func[5]
						+BlkField[comp][  i_jp1_kp1] * shape_func[6]
						+BlkField[comp][ip1_jp1_kp1] * shape_func[7];
		
					trajectory_outfile << A[comp] << "	";
				}
				trajectory_outfile << endl;
			}//! end if responsible_mpi_process
		}

	}//! end for positions
	//! NOTE: remove after DFT
// 	trajectory_outfile << endl;


	//!-----------------------------------
	//! ) clean up
	//!-----------------------------------
	for(INT32 comp=0; comp<3; comp++)
	{
		delete index_blk[comp];
		delete index_cell[comp];
		delete x_cell[comp] ;
	}


	trajectory_outfile.flush();
	trajectory_outfile.close();

	log_file <<"done." <<endl;


}

//!------------------------------------------------------------
//!- native_output_2D: Save all important infomrmation for one TL
//!		(which is basically all 3D Fields)
//!	NOTE: Fucntion "write_num_of_Blks_in_XS" has to be 
//!		called before any field is plotted !!!
//!------------------------------------------------------------
void CHybrid::native_output_2D(INT32 TL)
{



	log_file << endl << " Writing Native 2D Data ... " << endl;

	char filename[200];

	sprintf(filename,"%s/native/%s_2d_TL%d.dat",data_output_path, Run_Name, TL);
	Outfile2D.open( filename ,ofstream::binary);


	//! This has to be called before any field is plotted !
	write_num_of_Blks_in_XS();



	show_time(TIME_FIELD, id_FieldTime);
 	write_all_native_CS(id_FieldTime);

	show_time(TIME_PART, id_ParticleTime);
 	write_all_native_CS(id_ParticleTime);

// 	show_time(TIME_TOTAL, id_scratch_vector);
//  	write_all_native_CS(id_scratch_vector);

	show_responsible_Proc(id_scratch_scalar);
 	write_all_native_CS(id_scratch_scalar);

	//! ---- BField Cuts -----------------
 	write_all_native_CS(id_BTotal);

	//! ---- EField Cuts -----------------
 	write_all_native_CS(id_EField);

	//! ---- UField Cuts -----------------
 	write_all_native_CS(id_UI_plus);

	//! ---- rho Cuts -----------------
 	write_all_native_CS(id_rho_np1);


	//! ---- Eta Cuts --------------------
 	write_all_native_CS(id_Eta);


	//! ---- rotB Cuts -------------------
	calc_rot(id_BEven, id_rotB);
   	FULL_GN_UPDATE(id_rotB);
 	write_all_native_CS(id_rotB);


	//! ---- divB Cuts -------------------
	calc_div(id_BEven, id_divB);
   	FULL_GN_UPDATE(id_divB);
	write_all_native_CS(id_divB);


	//! ---- divE Cuts -------------------
	calc_div(id_EField, id_divE);
   	FULL_GN_UPDATE(id_divE);
	write_all_native_CS(id_divE);


	//! ---- id_PMagnetic -------------------
	square_Field(id_PMagnetic, id_BTotal);
// 	calc_macro_Force(id_PMagnetic, id_EField, id_UI_plus, id_BTotal);
	write_all_native_CS(id_PMagnetic);


	//! set zero PTotal here !
	set_zero_field(id_PTotal);
	add_multipliedField(id_PTotal, id_PMagnetic, 1.);

	//! ---- rhoSpecies/UI_Species Cuts -------------------
	for(int species=0; species<num_Charged_Species; species++)
	{


		//! write rho of species
		write_all_native_CS(id_rhoSpecies1 +species);


		//! provide rho
		//! gather Ui 
		collect_Ui_vth2_of_species(species, id_UI_Species1, 0,noVTH2);
		write_all_native_CS(id_UI_Species1 +species);

		calc_macro_Force(id_ForceSpecies1, id_EField, id_UI_Species1 +species, id_BTotal);
		write_all_native_CS(id_ForceSpecies1 +species);

		//! provide rho
		//! provide Ui
		//! gather vth2
		collect_Ui_vth2_of_species(species, id_UI_Species1, id_PISpecies1 ,getVTH2);
		write_all_native_CS(id_PISpecies1 +species);
		
		show_PE(species);
		write_all_native_CS(id_PESpecies1 +species);

		//! add pressures to total pressure
		add_multipliedField(id_PTotal, id_PISpecies1+species, 1.);
		add_multipliedField(id_PTotal, id_PESpecies1+species, 1.);


	}

	//! finally write total presure
	write_all_native_CS(id_PTotal);



	//! ---- extern_DensityFields Cuts -------------------
	for(INT32 rho_extern=0; rho_extern<num_externRhoVelocityFields; rho_extern++)
	{
		write_XCrossSection(id_externRho1 +rho_extern);
		write_YCrossSection(id_externRho1 +rho_extern);
		write_ZCrossSection(id_externRho1 +rho_extern);

		write_XCrossSection(id_extern_Ui1 +rho_extern);
		write_YCrossSection(id_extern_Ui1 +rho_extern);
		write_ZCrossSection(id_extern_Ui1 +rho_extern);
	}

	//! close and flush output
	Outfile2D.flush();
	Outfile2D.close();

	log_file << " writing finished." << endl << endl;


}


//!-------------------------------------------------------------//
//! silo_output_2D									//
//!-------------------------------------------------------------//
void CHybrid::silo_2D_Trajectory(INT32 id_trajectory)
{


	D_REAL** positions;
	char** time_string_array;
	D_REAL** SCField;
	INT32 num_values, num_positions;



	//!-----------------------------------
	//! 1) count number of values in file
	//!-----------------------------------
	trajectory_count_values(id_trajectory, num_values, num_positions);


	//!-----------------------------------
	//! 2) read values in file
	//!-----------------------------------
	read_trajectory(id_trajectory, num_positions, positions, time_string_array, SCField);


	//! write Orbit for each CrossSection
	silo_2DWrite_Trajectory(id_trajectory, num_positions, positions, 0);
	silo_2DWrite_Trajectory(id_trajectory, num_positions, positions, 1);
	silo_2DWrite_Trajectory(id_trajectory, num_positions, positions, 2);


	//!-----------------------------------
	//! 4) clean up
	//!-----------------------------------
	//! delete time_string_array memory	
	if(do_read_timestring)
	{
		for(INT32 pos=0; pos<num_positions; pos++)
		delete[] time_string_array[pos];

		delete[] time_string_array;
	}

	//! delete coordinates memory
	delete[] positions[0];
	delete[] positions[1];
	delete[] positions[2];

	delete[] positions;

	if(do_read_SpaceCraftField)
	{
		delete[] SCField[0];
		delete[] SCField[1];
		delete[] SCField[2];

		delete[] SCField;
	}

}


//!-------------------------------------------------------------//
//! silo_output_2D									//
//!-------------------------------------------------------------//
void CHybrid::silo_3D_Trajectory(INT32 id_trajectory)
{


	D_REAL** positions;
	char** time_string_array;
	D_REAL** SCField;
	INT32 num_values, num_positions;



	//!-----------------------------------
	//! 1) count number of values in file
	//!-----------------------------------
	trajectory_count_values(id_trajectory, num_values, num_positions);


	//!-----------------------------------
	//! 2) read values in file
	//!-----------------------------------
	read_trajectory(id_trajectory, num_positions, positions, time_string_array, SCField);


	//! write Orbit
	silo_3DWrite_Trajectory(id_trajectory, num_positions, positions);


	//!-----------------------------------
	//! 4) clean up
	//!-----------------------------------
	//! delete time_string_array memory	
	if(do_read_timestring)
	{
		for(INT32 pos=0; pos<num_positions; pos++)
		delete[] time_string_array[pos];

		delete[] time_string_array;
	}

	//! delete coordinates memory
	delete[] positions[0];
	delete[] positions[1];
	delete[] positions[2];

	delete[] positions;

	if(do_read_SpaceCraftField)
	{
		delete[] SCField[0];
		delete[] SCField[1];
		delete[] SCField[2];

		delete[] SCField;
	}
}



//!-------------------------------------------------------------//
//! silo_output_2D									//
//!-------------------------------------------------------------//
void CHybrid::silo_output_2D(void)
{


	if(TestParticle_Simulation)
	collect_RHOnp1_UIplus_LAM_GAM();

	//! open file, alloc memory etc.
	silo_2DWrite_prepare();

	log_file<<"prepare sucess"<<endl;
	
	//! write trajectories as point meshes
	for(INT32 id_trajectory=0; id_trajectory<num_trajectories; id_trajectory++)
	silo_2D_Trajectory(id_trajectory);
	//! write lineout as point meshes
	for(INT32 id_line=0; id_line<num_lines; id_line++)
	silo_2D_LineOut(id_line);



	//! write mesh always has to be called
	//! before fields can be written.
	silo_2DWrite_MeshCS(0);
	silo_2DWrite_MeshCS(1);
	silo_2DWrite_MeshCS(2);


	

	//! ---- BField Cuts -----------------
	silo_2DWrite_allCS(id_BTotal);

	//! ---- BField Cuts -----------------
	if(!TestParticle_Simulation)
	silo_2DWrite_allCS(id_BEven);

	//! ---- EField Cuts -----------------
	silo_2DWrite_allCS(id_EField);
	
	
#ifdef TRACK_PARTICLE

	collect_rho_Ji(0, id_scratch_scalar, id_scratch_vector, true);
	silo_2DWrite_allCS(id_scratch_scalar);
	silo_2DWrite_allCS(id_scratch_vector);

#endif

	if(!TestParticle_Simulation)
	{
        
		//! ---- UField Cuts -----------------
                silo_2DWrite_allCS(id_UI_plus);
                
                //! ---- divU Cuts ----------------
                calc_div(id_UI_plus, id_divU);
                FULL_GN_UPDATE(id_divU);
                silo_2DWrite_allCS(id_divU);
                
                //! ---- rho Cuts -----------------
                silo_2DWrite_allCS(id_rho_np1);
        
                //! ---- recombined rho Cuts -----------------
                silo_2DWrite_allCS(id_rho_np1_recombined);
                
                //! ---- Eta Cuts --------------------
                silo_2DWrite_allCS(id_Eta);
		
// 		silo_2DWrite_allCS(id_rho_rez);

        
                //! ---- rotB Cuts -------------------
                calc_rot(id_BEven, id_rotB);
                FULL_GN_UPDATE(id_rotB);
                silo_2DWrite_allCS(id_rotB);
        
                //! ---- divrotB Cuts ----------------
                calc_div(id_rotB, id_divrotB);
                FULL_GN_UPDATE(id_divrotB);
                silo_2DWrite_allCS(id_divrotB);
        
                //! ---- divB Cuts -------------------
                calc_div(id_BEven, id_divB);
                FULL_GN_UPDATE(id_divB);
                silo_2DWrite_allCS(id_divB);
        
                //! ---- id_gradB Cuts -------------------
                calc_grad(id_BTotal, id_gradB);
                FULL_GN_UPDATE(id_gradB);
                silo_2DWrite_allCS(id_gradB);
                
// 		log_file << " ja? "<<endl;		
		
                //! ---- show Num_PiC --------------
                show_MPiC(id_Refine_Rating);
                silo_2DWrite_allCS(id_Refine_Rating);		

// 		log_file << " ja? "<<endl;		
// 
		
//                 //! ---- Magnetic Pressure 
// 		square_Field(id_PMagnetic, id_BTotal);
//                 silo_2DWrite_allCS(id_PMagnetic);
// 		
//                 //! ---- Dynamic Pressure
//                 //! Vector quantity
//                 //! p_dyn^* = SUM_species ((n*)_species * (m*)_species * (u*)_species)
//                 set_zero_field(id_scratch_vector); //! Save Total Value
// 		set_zero_field(id_gradPI0);
// 
// 		  
// 		collect_Ui_vth2_of_species(0, id_UI_Species1, 0, noVTH2);
// 		Multiply_fields(id_gradPI0,id_UI_Species1,id_UI_Species1);
// 		Multiply_fields(id_scratch_vector,id_gradPI0,id_rhoSpecies1);
// 		  
                
               
//                 for(int species=0; species<1; species++) {
//                     //! Dynamic Pressue
//                     //! p_dyn^* = SUM_species ((n*)_species * (m*)_species * (u*)_species)
//                     //! density of species: id_rhoSpecies1 +species
//                     //! mass of species: Ion_Masses[species]
//                     //! velocitiy of species: id_UI_Species1 +species
//                     
//                     //! Collect Velocity of Species
//                    set_zero_field(id_gradPI0); //! Calculate u^2
//                    set_zero_field(id_gradB); //! Calulate n*u^2
//                     
//                    collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);
// 
// // 		   Multiply_fields(id_scratch_vector,id_UI_Species1 +species,id_UI_Species1 +species);
//                    Multiply_fields(id_gradPI0,id_UI_Species1 +species,id_UI_Species1 +species);
//                    Multiply_fields(id_gradB,id_gradPI0,id_rhoSpecies1 +species);
//                    add_multipliedField(id_scratch_vector,id_gradB,Ion_Masses[species]);
// // 		   add_multipliedField(id_scratch_vector,id_rhoSpecies1 +species,Ion_Charges[species]);
// // 		   FULL_GN_UPDATE(id_scratch_vector);
//                 }
//                 silo_2DWrite_allCS(id_scratch_vector);

                //! Mean Electron Temperature
                //! Electrons Adiabatic -> Density 
                //! Density of Species in field: id_rhoSpecies1 +species
                calc_electron_temperature(id_ElectronTemperature);
                silo_2DWrite_allCS(id_ElectronTemperature);
                
                
		//!NOTE: make sure not to use the same field twice
		
	// 	show_time(TIME_TOTAL, id_scratch_vector);
	//  	silo_2DWrite_allCS(id_scratch_vector);
	
// 		show_responsible_Proc(id_scratch_scalar);
// 		silo_2DWrite_allCS(id_scratch_scalar);

		silo_2DWrite_allCS(id_PhiDC);

	}
	
	for(int species=0; species<num_Particle_Species; species++)
	{

		//! write rho
		silo_2DWrite_allCS(id_rhoSpecies1 +species);

		//! provide rho
		//! gather Ui 
                collect_Ui_vth2_of_species(species, id_UI_Species1 +species, 0, noVTH2);
		silo_2DWrite_allCS(id_UI_Species1 +species);

                //! calc gyroradius ions
                calc_local_gyro_radius(id_gyro_Species1 + species,species,id_BTotal,id_UI_Species1 +species);
                silo_2DWrite_allCS(id_gyro_Species1 +species);
                //! calc gyroradius electrons
                calc_local_electron_gyro_radius(id_gyro_el_Species1 + species,species,id_BTotal,id_UI_Species1 +species);
                silo_2DWrite_allCS(id_gyro_el_Species1 +species);
                
 		calc_macro_Force(id_ForceSpecies1, id_EField, id_UI_Species1 +species, id_BTotal);
		silo_2DWrite_allCS(id_ForceSpecies1 +species);


		//! provide rho
		//! provide Ui
		//! gather vth2
		collect_Ui_vth2_of_species(species, id_UI_Species1+species, id_PISpecies1+species ,getVTH2);
		//! for temperature in eV
		add_multipliedField(id_PISpecies1 +species, id_PISpecies1 +species, Ion_Masses[species]*SI_m0/e*SI_v0*SI_v0);
		silo_2DWrite_allCS(id_PISpecies1 +species);
		
		
		show_PE(species);
		silo_2DWrite_allCS(id_PESpecies1 + species);

// 		calc_grad(id_PISpecies1 +species, id_gradPI0);
//  		FULL_GN_UPDATE(id_gradPI0);
// 		silo_2DWrite_allCS(id_gradPI0);

// 		calc_rot(id_UI_Species1 +species, id_ForceSpecies1 +species);
// 		FULL_GN_UPDATE(id_ForceSpecies1 +species);
// 		silo_2DWrite_allCS(id_ForceSpecies1 +species);


// 		calc_div(id_UI_plus, id_PISpecies1 +0);
// 		FULL_GN_UPDATE(id_PISpecies1 +0);
// 		silo_2DWrite_allCS(id_PISpecies1 +0);



	}
	
	for(int species=num_Particle_Species; species<num_Charged_Species; species++)
	{
		silo_2DWrite_allCS(id_rhoSpecies1 +species);
		//! provide rho

	}	
	
#if defined nonadiabatic_gradPE_TERM
	//! write total electron pressure
	silo_2DWrite_allCS(id_PEtotal);
#else
	log_file << "  Compute total electron pressure... " << endl;
	set_zero_field(id_PEtotal);
	for(int species=0; species<num_Charged_Species; species++)
	{
	  add_multipliedField(id_PEtotal, id_PESpecies1+species, 1.);
	}
	//! write total electron pressure
	silo_2DWrite_allCS(id_PEtotal);
#endif
	
	
	

#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
	for(INT32 species=0; species<num_Neutral_Species; ++species)
	{

		//! write number density
		silo_2DWrite_allCS(id_numberdensity_neutralSpecies1     + species);
		
		//! write bulk velocity
		silo_2DWrite_allCS(id_velocity_neutralSpecies1          + species);
		
		//! write pressure
		silo_2DWrite_allCS(id_pressure_neutralSpecies1          + species);
		
		//! write new electron beta
		silo_2DWrite_allCS(id_new_electron_beta_neutralSpecies1 + species);
		
	}
#endif
	
#if defined(use_ion_production_as_field)
	for(INT32 species=0; species<num_ion_prod_fields; ++species)
	{

		//! write number density
		silo_2DWrite_allCS(id_density_ionProdSpecies1     + species);
		
		//! write bulk velocity
		silo_2DWrite_allCS(id_velocity_ionProdSpecies1          + species);		
	
	}
#endif
	
#if defined(use_dust_species_as_field)
	for(INT32 species=0; species<num_Dust_Species; ++species)
	{

		//! write number density
		silo_2DWrite_allCS(id_density_dustSpecies1     + species);
		
		//! write bulk velocity
		silo_2DWrite_allCS(id_velocity_dustSpecies1    + species);		
	
	}
#endif
	
	//! ---- precalced_DensityFields Cuts -------------------
	for(INT32 rho_extern=0; rho_extern<num_externRhoVelocityFields; rho_extern++)
	{
		silo_2DWrite_allCS(id_externRho1 +rho_extern);
		silo_2DWrite_allCS(id_extern_Ui1 +rho_extern);
	}

	//! ---- averaged fields -------------------
	if(!TestParticle_Simulation)
	for(INT32 average_field=0; average_field<num_average_fields; average_field++)
	silo_2DWrite_allCS(id_average_Field1 +average_field);


	//! close file, delete memory etc
	silo_2DFlush_CleanUp();


}

//!-------------------------------------------------------------//
//! silo_output_3D								//
//!-------------------------------------------------------------//
void CHybrid::silo_output_3D(void)
{


	//! open file, alloc memory etc.
	silo_3DWrite_prepare();


	//! write trajectories as point meshes
	for(INT32 id_trajectory=0; id_trajectory<num_trajectories; id_trajectory++)
	silo_3D_Trajectory(id_trajectory);


	//! write lineout as point meshes
	for(INT32 id_line=0; id_line<num_lines; id_line++)
	silo_3D_LineOut(id_line);

	//! write mesh always has to be called
	//! before fields can be written.
	silo_3DWrite_Mesh();

	//! ---- BField Cuts -----------------
 	silo_3DWrite_Field(id_BTotal);

	//! ---- EField Cuts -----------------
	silo_3DWrite_Field(id_EField);


#ifdef TRACK_PARTICLE
	collect_rho_Ji(0, id_scratch_scalar, id_scratch_vector, true);
	silo_3DWrite_Field(id_scratch_scalar);
	silo_3DWrite_Field(id_scratch_vector);

#endif

	if(!TestParticle_Simulation)
	{
		//! ---- UField Cuts -----------------
		silo_3DWrite_Field(id_UI_plus);
	
		//! ---- rho Cuts -----------------
		silo_3DWrite_Field(id_rho_np1);
                
                //! ---- recombined rho Cuts -----------------
                silo_3DWrite_Field(id_rho_np1_recombined);

		//! ---- Eta Cuts --------------------
		silo_3DWrite_Field(id_Eta);
	
		//! ---- rotB Cuts -------------------
		calc_rot(id_BEven, id_rotB);
		FULL_GN_UPDATE(id_rotB);
		silo_3DWrite_Field(id_rotB);
	
		//! ---- divB Cuts -------------------
		calc_div(id_BEven, id_divB);
		FULL_GN_UPDATE(id_divB);
		silo_3DWrite_Field(id_divB);
	
		//! ---- divE Cuts -------------------
		calc_div(id_EField, id_divE);
		FULL_GN_UPDATE(id_divE);
		silo_3DWrite_Field(id_divE);
	

                //! ---- Magnetic Pressure 
		square_Field(id_PMagnetic, id_BTotal);
                silo_3DWrite_Field(id_PMagnetic);
                
                //! ---- show Num_PiC --------------
                show_MPiC(id_Refine_Rating);
                silo_3DWrite_Field(id_Refine_Rating);           
                
                //! ---- Dynamic Pressure
                //! Vector quantity
                //! p_dyn^* = SUM_species ((n*)_species * (m*)_species * (u*)_species)
                set_zero_field(id_scratch_vector); //! Save Total Value
		set_zero_field(id_gradPI0);

		  
		collect_Ui_vth2_of_species(0, id_UI_Species1, 0, noVTH2);
		Multiply_fields(id_gradPI0,id_UI_Species1,id_UI_Species1);
		Multiply_fields(id_scratch_vector,id_gradPI0,id_rhoSpecies1);
                
                silo_3DWrite_Field(id_scratch_vector);
                
                //! Mean Electron Temperature
                //! Electrons Adiabatic -> Density 
                //! Density of Species in field: id_rhoSpecies1 +species
                calc_electron_temperature(id_ElectronTemperature);
                //                 if(TL==10)
                //                 calc_RecombinationAlpha(id_ElectronTemperature);
                silo_3DWrite_Field(id_ElectronTemperature);
                
                /*
                for(int species=0; species<1; species++) {
                    //! Dynamic Pressue
                    //! p_dyn^* = SUM_species ((n*)_species * (m*)_species * (u*)_species)
                    //! density of species: id_rhoSpecies1 +species
                    //! mass of species: Ion_Masses[species]
                    //! velocitiy of species: id_UI_Species1 +species
                    
                    //! Collect Velocity of Species
                   set_zero_field(id_gradPI0); //! Calculate u^2
                   set_zero_field(id_gradB); //! Calulate n*u^2
                    
                   collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);

// 		   Multiply_fields(id_scratch_vector,id_UI_Species1 +species,id_UI_Species1 +species);
                   Multiply_fields(id_gradPI0,id_UI_Species1 +species,id_UI_Species1 +species);
                   Multiply_fields(id_gradB,id_gradPI0,id_rhoSpecies1 +species);
                   add_multipliedField(id_scratch_vector,id_gradB,Ion_Masses[species]);
// 		   add_multipliedField(id_scratch_vector,id_rhoSpecies1 +species,Ion_Charges[species]);
// 		   FULL_GN_UPDATE(id_scratch_vector);
                }*/
                
	}


	//! ---- Rho Cuts -------------------
	for(int species=0; species<num_Charged_Species; species++)
	{

		//! write rho
		silo_3DWrite_Field(id_rhoSpecies1 +species);

		//! provide rho
		//! gather Ui 
		collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);
		silo_3DWrite_Field(id_UI_Species1 +species);

                calc_macro_Force(id_ForceSpecies1, id_EField, id_UI_Species1 +species, id_BTotal);
                silo_3DWrite_Field(id_ForceSpecies1 +species);
                
		//! provide rho
		//! provide Ui
		//! gather vth2
		collect_Ui_vth2_of_species(species, id_UI_Species1, id_PISpecies1 ,getVTH2);
		silo_3DWrite_Field(id_PISpecies1 +species);
		
		show_PE(species);
		silo_3DWrite_Field(id_PESpecies1 + species);

	}

#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
	if(num_Neutral_Species>0)
	for(INT32 neutralSpec=0; neutralSpec<num_Neutral_Species; neutralSpec++)
	{
		silo_3DWrite_Field(id_numberdensity_neutralSpecies1+neutralSpec);
		silo_3DWrite_Field(id_velocity_neutralSpecies1+neutralSpec);
	}	
#endif
	
	
#ifdef use_dust_species_as_field	
	if(num_Dust_Species>0)
	for(INT32 dustSpec=0; dustSpec<num_Dust_Species; dustSpec++)
	{
		silo_3DWrite_Field(id_density_dustSpecies1+dustSpec);
		silo_3DWrite_Field(id_velocity_dustSpecies1+dustSpec);
	}
#endif	

	//! ---- precalced_DensityFields Cuts -------------------
	for(INT32 rho_extern=0; rho_extern<num_externRhoVelocityFields; rho_extern++)
	{
		silo_3DWrite_Field(id_externRho1 +rho_extern);
		silo_3DWrite_Field(id_extern_Ui1 +rho_extern);
	}
	
	//! close file, delete memory etc
	silo_3DFlush_CleanUp();



}
//!------------------------------------------------------------
//!- write_num_of_Blks_in_XS
//!------------------------------------------------------------
void CHybrid::write_num_of_Blks_in_XS(void)
{


   memset(num_Blks_in_CrossSection,0,3*sizeof(INT64));

   //!-----------------------------------------//
   //! 	Loop over all Blks			 //
   //!----------------------------------------//
   for(INT32 level=0; level<=MAX_LEVEL; level++)
   {
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{



		 for(INT32 comp=0; comp<3; comp++)
		  if(   CROSS_SECTION[comp] >= temp_Block->origin[comp]
		     && CROSS_SECTION[comp] <  temp_Block->origin[comp]
						+Blk_Length_of[level][comp])
		  num_Blks_in_CrossSection[comp]++;

		  temp_Block = temp_Block->next_Blk_of_BlockList;
	}
   }

   //! use XSEC to identify correct CrossSection
   Outfile2D.write(reinterpret_cast<char*> (num_Blks_in_CrossSection),3*sizeof(INT64));



	
}

//!------------------------------------------------------------
//! write_all_native_CS
//!------------------------------------------------------------
void CHybrid::write_all_native_CS(INT32 type)
{
	write_XCrossSection(type);
	write_YCrossSection(type);
	write_ZCrossSection(type);
}




//!------------------------------------------------------------
//!- write_XCrossSection: Writes out the Field values for
//!	a given Cross Section Coordinate specified in defines.h.
//!   - Whether a Block is within the CS is decided by its 
//!	  Interval [orig,orig + length]
//!   - No Interpolation is performed. Data is written out in 
//!     a nearest neighbour style
//!------------------------------------------------------------
void CHybrid::write_XCrossSection(INT32 type)
{


   INT32 XSEC  = 0;
   INT32 COMP   = COMPs_FType[type];
   INT32 VisuNr = VisuNr_FType[type];
   //! ----- Initialize File-Stream ------------------------------


   //! use XSEC to identify correct CrossSection
   Outfile2D.write(reinterpret_cast<char*> (&XSEC),sizeof(INT32));
   //! use type to apply correct dimension in Visu
   Outfile2D.write(reinterpret_cast<char*> (&type),sizeof(INT32));
   //! use COMP to define number of components
   Outfile2D.write(reinterpret_cast<char*> (&COMP),sizeof(INT32));
   //! use VisuNr to order field in Visus DropDown Menu
   Outfile2D.write(reinterpret_cast<char*> (&VisuNr),sizeof(INT32));
   //! use Field_Name to define Name which appears in Visu
   Outfile2D.write(reinterpret_cast<char*> (Field_Name[type]),Field_Name_Size*sizeof(char));

   //! ----- Cut Preparation -----------------------------------
   //! in order to draw obstacle norm_pos (in Box) has
   //! to be stored:

   calc_cells_to_cut(0, cell_indedx_L);

   FILE_REAL norm_pos = (CROSS_SECTION[0]+Box_Origin[0])/LX;
   Outfile2D.write(reinterpret_cast<char*> (&norm_pos),sizeof(FILE_REAL));

   //!-----------------------------------------//
   //! 	Loop over Blks				//
   //!-----------------------------------------//
   for(INT32 level=0; level<=MAX_LEVEL; level++)
   {
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{


		  if(      CROSS_SECTION[0] >= temp_Block->origin[0]
			&& CROSS_SECTION[0] <  temp_Block->origin[0]+Blk_Length_of[level][0])
		  write_XBlock_Cut(level,type,temp_Block);

		  temp_Block = temp_Block->next_Blk_of_BlockList;
	}
   }



	
}

//!------------------------------------------------------------
//!- write_YCrossSection: Writes out the Field values for
//!	a given Cross Section Coordinate specified in defines.h.
//!   - Whether a Block is within the CS is decided by its 
//!	  Interval [orig,orig + length]
//!   - No Interpolation is performed. Data is written out in 
//!     a nearest neighbour style
//!------------------------------------------------------------
void CHybrid::write_YCrossSection(INT32 type)
{


   INT32 XSEC  = 1;
   INT32 COMP   = COMPs_FType[type];
   INT32 VisuNr = VisuNr_FType[type];
   //! ----- Initialize File-Stream ------------------------------

   //! use XSEC to identify correct CrossSection
   Outfile2D.write(reinterpret_cast<char*> (&XSEC),sizeof(INT32));
   //! use type to apply correct dimension in Visu
   Outfile2D.write(reinterpret_cast<char*> (&type),sizeof(INT32));
   //! use COMP to define number of components
   Outfile2D.write(reinterpret_cast<char*> (&COMP),sizeof(INT32));
   //! use VisuNr to order field in Visus DropDown Menu
   Outfile2D.write(reinterpret_cast<char*> (&VisuNr),sizeof(INT32));
   //! use Field_Name to define Name which appears in Visu
   Outfile2D.write(reinterpret_cast<char*> (Field_Name[type]),Field_Name_Size*sizeof(char));



   //! ----- Cut Preparation -----------------------------------
   //! in order to draw obstacle norm_pos (in Box) has
   //! to be stored:

   calc_cells_to_cut(1, cell_indedx_L);

   FILE_REAL norm_pos = (CROSS_SECTION[1]+Box_Origin[1])/LY;
   Outfile2D.write(reinterpret_cast<char*> (&norm_pos),sizeof(FILE_REAL));


   //!-----------------------------------------//
   //! 	Loop over all Blks			 //
   //!----------------------------------------//
   for(INT32 level=0; level<=MAX_LEVEL; level++)
   {
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{

		  if(      CROSS_SECTION[1] >= temp_Block->origin[1]
			&& CROSS_SECTION[1] <  temp_Block->origin[1]+ Blk_Length_of[level][1])
		  //! Write down child Field
		  write_YBlock_Cut(level,type,temp_Block);

		  temp_Block = temp_Block->next_Blk_of_BlockList;
	}

    }



	
}


//!------------------------------------------------------------
//!- write_ZCrossSection: Writes out the Field values for
//!	a given Cross Section Coordinate specified in defines.h.
//!   - Whether a Block is within the CS is decided by its 
//!	  Interval [orig,orig + length]
//!   - No Interpolation is performed. Data is written out in 
//!     a nearest neighbour style
//!------------------------------------------------------------
void CHybrid::write_ZCrossSection(INT32 type)
{


   INT32 XSEC  = 2;
   INT32 COMP   = COMPs_FType[type];
   INT32 VisuNr = VisuNr_FType[type];
   //! ----- Initialize File-Stream ------------------------------

   //! use XSEC to identify correct CrossSection
   Outfile2D.write(reinterpret_cast<char*> (&XSEC),sizeof(INT32));
   //! use type to apply correct dimension in Visu
   Outfile2D.write(reinterpret_cast<char*> (&type),sizeof(INT32));
   //! use COMP to define number of components
   Outfile2D.write(reinterpret_cast<char*> (&COMP),sizeof(INT32));
   //! use VisuNr to order field in Visus DropDown Menu
   Outfile2D.write(reinterpret_cast<char*> (&VisuNr),sizeof(INT32));
   //! use Field_Name to define Name which appears in Visu
   Outfile2D.write(reinterpret_cast<char*> (Field_Name[type]),Field_Name_Size*sizeof(char));



   //! ----- Cut Preparation -----------------------------------
   //! in order to draw obstacle norm_pos (in Box) has
   //! to be stored:

   calc_cells_to_cut(2, cell_indedx_L);

   FILE_REAL norm_pos = (CROSS_SECTION[2]+Box_Origin[2])/LZ;
   Outfile2D.write(reinterpret_cast<char*> (&norm_pos),sizeof(FILE_REAL));


   //!-----------------------------------------//
   //! 	Loop over all Blks			 //
   //!----------------------------------------//
   for(INT32 level=0; level<=MAX_LEVEL; level++)
   {
	CBlock* temp_Block = BlockList_of_Lev[level];
	while(temp_Block)
	{


		  if(      CROSS_SECTION[2] >= temp_Block->origin[2]
			&& CROSS_SECTION[2] <  temp_Block->origin[2] +Blk_Length_of[level][2])
		  write_ZBlock_Cut(level,type,temp_Block);

		  temp_Block = temp_Block->next_Blk_of_BlockList;
	}
    }


	
}


//!------------------------------------------------------
//!- write_XBlock_Cut: The function "calc_cells_to_cut" has to 
//!    be called in advance to set "cell_indedx_L[Level]".
//!    Otherwise wrong cell indeces will be used..
//!------------------------------------------------------
void CHybrid::write_XBlock_Cut(INT32 Level, INT32 F_Type, CBlock* Block)
{

	 INT32 COMP = COMPs_FType[F_Type];

	 FILE_REAL Block_Cords[3] = {Block->origin[1]+Box_Origin[1],
				     Block->origin[2]+Box_Origin[2],
				     Block->RLevel};
	//! Write down root Block Cords
	 Outfile2D.write(reinterpret_cast<char*> (Block_Cords), 3*sizeof(FILE_REAL));


	INT32 temp4 = cell_indedx_L[Level] * BlkNds_Y * BlkNds_Z;
	for(int X=0; X <COMP; X++)
	{
		INT32 comp = X * BlkNds_Y *BlkNds_Z;
		for(int b=0; b<BlkNds_Y; b++)
		{
	   	  INT32 temp5 = b*BlkNds_Z;

	   	  for(int c=0; c<BlkNds_Z; c++)
	   	  XBlock_CrossSection[comp  +temp5 +c] = Block->Field_Type[F_Type][X*num_nodes_in_block +temp4 +temp5 +c];

		}
	}

	Outfile2D.write(reinterpret_cast<char*> (XBlock_CrossSection),COMP*BlkNds_Y *BlkNds_Z*sizeof(FILE_REAL));
	num_Blocks_transfered_L[Level]++;

}

//!------------------------------------------------------------
//!- write_YBlock_Cut: The function "calc_cells_to_cut" has to 
//!    be called in advance to set "cell_indedx_L[Level]".
//!    Otherwise wrong cell indeces will be used..
//!------------------------------------------------------------
void CHybrid::write_YBlock_Cut(INT32 Level, INT32 F_Type,CBlock* Block)
{

	 INT32 COMP = COMPs_FType[F_Type];

	 FILE_REAL Block_Cords[3] = {Block->origin[0]+Box_Origin[0],
				     Block->origin[2]+Box_Origin[2],
				     Block->RLevel};
	//! Write down root Block Cords
	Outfile2D.write(reinterpret_cast<char*> (Block_Cords), 3*sizeof(FILE_REAL));

	INT32 temp4,temp6;
	INT32 temp5 = cell_indedx_L[Level] * BlkNds_Z;
	for(int X=0; X <COMP; X++)
	{
		INT32 comp = X * BlkNds_X *BlkNds_Z;
		for(int a=0; a<BlkNds_X; a++)
		{

		  temp4 = a*BlkNds_Y*BlkNds_Z;
		  temp6 = a*BlkNds_Z;

	   	  for(int c=0; c<BlkNds_Z; c++)
	   	  YBlock_CrossSection[comp +temp6 +c] = Block->Field_Type[F_Type][X*num_nodes_in_block +temp4 +temp5 +c];

		}
	}

	Outfile2D.write(reinterpret_cast<char*> (YBlock_CrossSection),COMP*BlkNds_X *BlkNds_Z*sizeof(FILE_REAL));
	num_Blocks_transfered_L[Level]++;

}


//!------------------------------------------------------------
//!- write_ZBlock_Cut: The function "calc_cells_to_cut" has to 
//!    be called in advance to set "cell_indedx_L[Level]".
//!    Otherwise wrong cell indeces will be used..
//!------------------------------------------------------------
void CHybrid::write_ZBlock_Cut(INT32 Level, INT32 F_Type,CBlock* Block)
{

	INT32 COMP = COMPs_FType[F_Type];

	 FILE_REAL Block_Cords[3] = {Block->origin[0]+Box_Origin[0],
				     Block->origin[1]+Box_Origin[1],
				     Block->RLevel};
	//! Write down root Block Cords
	Outfile2D.write(reinterpret_cast<char*> (Block_Cords), 3*sizeof(FILE_REAL));

	INT32 temp4, temp5, temp6;
	int c = cell_indedx_L[Level];
	for(int X=0; X <COMP; X++)
	{
		INT32 comp = X * BlkNds_X *BlkNds_Y;
		for(int a=0; a<BlkNds_X; a++)
		{

		  temp4 = a*BlkNds_Y*BlkNds_Z;
		  temp6 = a*BlkNds_Y;

	   	  for(int b=0; b<BlkNds_Y; b++)
		  {
		   temp5 = b * BlkNds_Z;
	   	   ZBlock_CrossSection[comp +temp6 +b] = Block->Field_Type[F_Type][X*num_nodes_in_block +temp4 +temp5 +c];
		  }
		
		}
	}

	Outfile2D.write(reinterpret_cast<char*> (ZBlock_CrossSection),COMP*BlkNds_X *BlkNds_Y*sizeof(FILE_REAL));
	num_Blocks_transfered_L[Level]++;

}



//!------------------------------------------------------------
//!- native_output_3D: Save all important infomrmation for one TL
//!		(which is basically all 3D Fields)
//!------------------------------------------------------------
void CHybrid::native_output_3D(INT32 TL)
{


	write_Field3D(id_BEven,TL, false);


}

//!------------------------------------------------------------
//!- gnuplot_output:
//!------------------------------------------------------------
void CHybrid::gnuplot_output(INT32 TL)
{

	if(RB_X>1 || RB_Y>1 || RB_Z>1)
	return;

	log_file << endl << " Writing Gnuplot 2D Data ..." << endl;

	INT32 lay[3];

	for(INT32 cs=0; cs<3; cs++)
	{
		calc_cells_to_cut(cs, cell_indedx_L);
		lay[cs] = cell_indedx_L[0];
	}


	//! For gnuplot output it is assumed that there is only 
	//! a single Block in level 0 and no further refinements
	CBlock* gnu_Block = BlockList_of_Lev[0];

	//! Write constant fields
	write_field(gnu_Block, lay, id_Eta, 0);
	gnuplot_write_inner(gnu_Block->Flag, lay[0], lay[1], lay[2]);
	gnuplot_write_grid(gnu_Block, 	     lay[0], lay[1], lay[2]);




	//! Write simulation results

	//! EM-Fields
	write_field(gnu_Block, lay, id_BTotal, TL);
	write_field(gnu_Block, lay, id_EField, TL);


	//! rotB & divB have to be calculated first
	gnu_Block->calc_rot(id_BEven, id_rotB);
	write_field(gnu_Block, lay, id_rotB, TL);

	gnu_Block->calc_div(id_BEven, id_divB);
	write_field(gnu_Block, lay, id_divB, TL);


	//! write total moments 
	write_field(gnu_Block, lay, id_UI_plus, TL);
	write_field(gnu_Block, lay, id_rho_np1, TL);


	//! collect and write moments of each species
	for(INT32 species=0; species<num_Charged_Species; species++)
	write_field(gnu_Block, lay, id_rhoSpecies1+species, TL);


	//! NOTE:
	//! VISUALIZAION MAY LOOK FALSE and different to TB Code 
	//! for obstacle ions in case MCD_J2U is high (eg. 0.2).
	//! This is due to the fact, that TB Codes does not use
	//! a lower limit for Visualization but uses MCD
	//! for COMPUTATION of Fields.
	//! So all quantities (B,E,rho and U) are the same in both Codes,
	//! but in this Code visulizes the intern U,
	//! in TB Code intern and visulized U are different.

	for(INT32 species=0; species<num_Charged_Species; species++)
	{
		collect_Ui_vth2_of_species(species, id_UI_Species1+species, 0, noVTH2);
		write_field(gnu_Block, lay, id_UI_Species1+species, TL);

	}

	//! ---- precalced_DensityFields Cuts -------------------
	for(INT32 rho_extern=0; rho_extern<num_externRhoVelocityFields; rho_extern++)
	write_field(gnu_Block, lay, id_externRho1 +rho_extern, TL);

	log_file << " writing finished." << endl;
}


//!--------------------------------------------------------
//!- secure_state_files:
//!- move Statefile after successful restart
//!--------------------------------------------------------
void CHybrid::secure_state_files(void)
{

   	char buffer[200];

	log_file << "  secure state ... ";

	//! move all file from State to State_TL"TimeLevel"
	//! only master moves new to old state
	if(!mpi_myRank)
	{	 
		sprintf(buffer,"mv -f State State_TL%d",TL);
		system(buffer);
	}


	log_file << "done." << endl;
	
	//! use barrier here otherwise
	//! it may happen that som processes
	//! start to save state  while another is moving
	//! their state files 
	synchronize_allProcesses_MPI();

}


//!--------------------------------------------------------
//!- save state
//!- USE TIME LEVEL AS TRIGGER FOR SAVE STATE RATHER
//!  THAN RUN TIME SINCE RUN TIME ALWAYS RESULTED IN
//!  SYNCHRONIZATION ERRORS !!!
//!--------------------------------------------------------
void CHybrid::save_state(void)
{


	//! NOTE:
	//! In very rare cases the states of some processes (2 out of 128) have
	//! not been present in the state directory even though it was protocolled
	//! within their log_files that they have been written
	//! -> it seem to be some filesytem synchronization problem, so aditional
	//!    barriers in write state were set


// 	if(!TL_SAVE_STATE)
// 	return;


	
	
// 	if( !(TL % TL_SAVE_STATE ==0) )
// 	{
// 		log_file << " Next save state in "<< TL_SAVE_STATE -TL%TL_SAVE_STATE <<" TL." << endl;
// 		return;
// 	}
	
	
	log_file << " Saving state..." << endl;
	synchronize_allProcesses_MPI();	
	
	//!- move any statefile BEFORE writing state starts, else
	//!  it may happen when crash in writing that some
	//!  state file are up to date and others are not !!!
	backup_state_files();
	synchronize_allProcesses_MPI();	
	


	//! Write block file including respective mpi_process
	//! It is sufficient to write Blockfile once since it
	//! is identical on every process
	//! -> only Master will write block_file
	if(!mpi_myRank)
	write_BlockFile();

	//! write density and velocity state files
	write_Field3D(id_rho_np1,TL,true);
	synchronize_allProcesses_MPI();	
	
	for(int species=0; species<num_Charged_Species; species++)
	{
	  write_Field3D(id_rhoSpecies1 +species,TL,true);
	  synchronize_allProcesses_MPI();	
	  collect_Ui_vth2_of_species(species, id_UI_Species1, 0, noVTH2);
	  write_Field3D(id_UI_Species1 +species,TL,true);
	  synchronize_allProcesses_MPI();	
	}
	
	write_Field3D(id_UI_plus,TL,true);
	synchronize_allProcesses_MPI();	
	
	synchronize_allProcesses_MPI();	
	//! - carefull, id_BTotal and id_BEven may have the same names
	//! - id_BTotal should not be required
	//    write_Field3D(id_BTotal,TL,true);
	write_Field3D(id_BEven,TL,true);
	synchronize_allProcesses_MPI();	

	write_Field3D(id_EField,TL,true);
	synchronize_allProcesses_MPI();
        
        write_Field3D(id_PEtotal,TL,true);
        synchronize_allProcesses_MPI();

	if(!TestParticle_Simulation)
	write_Field3D(id_PhiDC,TL,true);
	synchronize_allProcesses_MPI();

	if(!TestParticle_Simulation)
	for(INT32 average_field=0; average_field<num_average_fields; average_field++)
	write_Field3D(id_average_Field1 +average_field, TL, true);

	
	write_Particle();
	synchronize_allProcesses_MPI();	
	
	//! restart timer for save state
	
	log_file << endl << " writing finished." << endl;
	
	//    synchronize_allProcesses_MPI();
	//    time(&last_save_state_time);
	
	//    synchronize_allProcesses_MPI();
	//    global_MPI_lastState_time = MPI_Wtime();
	

}


//!--------------------------------------------------------
//!- backup_state_files:
//!- move any statefile BEFORE writing state starts,
//!  else it may happens when crash in writing that
//!  some state file are up to date and others are
//!  not.
//!--------------------------------------------------------
void CHybrid::backup_state_files(void)
{

   	char buffer[200];

	log_file << "  backing up last state ... ";

	//! move all file from State to Last_State
	//! only master moves new to old state
	if(!mpi_myRank)
	{
		sprintf(buffer,"mkdir Last_State");
		system(buffer);

		sprintf(buffer,"mv -f State/* Last_State/");
		system(buffer);
	}


	log_file << "done." << endl;
	
	//! use barrier here otherwise
	//! it may happen that som processes
	//! start to save state  while another is moving
	//! their state files 
	synchronize_allProcesses_MPI();


}


//!--------------------------------------------------------
//!- count_BlockTree_mpi_processes
//!--------------------------------------------------------
INT32 CHybrid::count_BlockTree_mpi_processes(void)
{

	INT32 number_of_BlockTree_mpi_processes = 0;

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{


		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(temp_Block->responsible_mpi_process > number_of_BlockTree_mpi_processes)
			number_of_BlockTree_mpi_processes = temp_Block->responsible_mpi_process;


			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}

	return number_of_BlockTree_mpi_processes +1;


}


//!--------------------------------------------------------
//!- reorganize_memory_for_assigned_blocks
//!--------------------------------------------------------
void CHybrid::reorganize_memory_for_assigned_blocks(INT32 restore_start_id, INT32 restore_end_id)
{


	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			//! in general mpi_myRank will not be within the
			//! interval [restore_start_id:restore_end_id[
			//! so first erease memory that was allocated
			//! during restore
			if(temp_Block->responsible_mpi_process == mpi_myRank)
			temp_Block->delete_process_specific_Memory();

			for(INT32 restore_process_id = restore_start_id;  restore_process_id<restore_end_id; restore_process_id++)
			{

				//! reallocate memory in case mpi_myRank takes responsibility 
				if(temp_Block->responsible_mpi_process == restore_process_id)
				{

					temp_Block->alloc_process_specific_Memory();
					temp_Block->set_Fields();
				}
			}

			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}
}


//!--------------------------------------------------------
//!- reorganize_memory_for_assigned_blocks
//!--------------------------------------------------------
void CHybrid::reset_responsible_mpi_process(INT32 number_of_processes_in_BlockFile)
{


	INT32 processes_to_restore_each_process = number_of_processes_in_BlockFile / mpi_num_processes;
	INT32 additional_processes = number_of_processes_in_BlockFile%mpi_num_processes;


	for(INT32 rank=0; rank<mpi_num_processes; rank++)
	{

		INT32 restore_start_id = (rank +0) * processes_to_restore_each_process +additional_processes;
		INT32 restore_end_id   = (rank +1) * processes_to_restore_each_process +additional_processes;

		if(rank < additional_processes)
		{

			restore_start_id -= (additional_processes-rank-0);
			restore_end_id   -= (additional_processes-rank-1);

		}


		//! loop across all files to restore
		for(INT32 restore_process_id = restore_start_id;  restore_process_id<restore_end_id; restore_process_id++)
		{
	
			for(INT32 level=0; level<=MAX_LEVEL; level++)
			{
			
				CBlock *temp_Block = BlockList_of_Lev[level];
				while(temp_Block)
				{
					if(temp_Block->responsible_mpi_process == restore_process_id)
					temp_Block->responsible_mpi_process = rank;
	
					temp_Block = temp_Block->next_Blk_of_BlockList;
				}
			}
		}
	}

}

//!--------------------------------------------------------
//!- restore_state
//!--------------------------------------------------------
bool CHybrid::restore_state(void)
{




	//! start timer for save state
	time(&last_save_state_time);
	
	particle_restored = false;
	
	log_file << endl << " Searching for state files ...     ";
	
	
	if(!read_BlockFile()) return false;


	//!-----------------------------------------------------------


	//! set default values for restore assuming number of processes
	//! this run equals number of of processes in statefile
	INT32 processes_to_restore_each_process = 1;

	INT32 restore_start_id = mpi_myRank +0;
	INT32 restore_end_id   = mpi_myRank +1;


	//! count number of involved mpi_processes
	INT32 number_of_processes_in_BlockFile = count_BlockTree_mpi_processes();



	log_file << " number of processes specified in Block File: " << number_of_processes_in_BlockFile << endl;
	log_file << " number of processes specified for this job: "  << mpi_num_processes << endl;
	log_file << endl;


	//! THREE CASES MAY OCCUR:

	//! 1) number_of_processes_in_BlockFile  < 'mpi_num_processes'
	//!	 -> DO RESTORE AS USUAL
	//!	 -> REASSIGNE BLOCKS IN REDISTRIBUTION STEP

	//! 2) number_of_processes_in_BlockFile == 'mpi_num_processes'
	//!	 -> DO RESTORE AS USUAL

	//! 3) number_of_processes_in_BlockFile  > 'mpi_num_processes'
	//!	 -> ASSIGN EACH PROCESS AN INTERVAL OF PROCESSES TO RESTORE
	//!	 -> REORGANIZE ALLOCATED FIELD MEMORY
	//!	 -> RESTORE FILES
	//!	 -> FINALLY SET 'mpi_myRank' TO EACH BLOCK OF INTERVAL

	//! check whether to few processes are available
	if(mpi_num_processes < number_of_processes_in_BlockFile)
	{

		log_file << " number of processes specified in Block File exeeds" << endl;
		log_file << " number of processes specified for this simulation" << endl;

		processes_to_restore_each_process = number_of_processes_in_BlockFile / mpi_num_processes;

		INT32 additional_processes = number_of_processes_in_BlockFile%mpi_num_processes;

		log_file << " Processes_to_restore_each_process: " <<   processes_to_restore_each_process << endl;
		log_file << " The first " << additional_processes  << " processes will restore one additional process." << endl;



		restore_start_id = (mpi_myRank +0) * processes_to_restore_each_process +additional_processes;
		restore_end_id   = (mpi_myRank +1) * processes_to_restore_each_process +additional_processes;

		if(mpi_myRank < additional_processes)
		{

			restore_start_id -= (additional_processes-mpi_myRank-0);
			restore_end_id   -= (additional_processes-mpi_myRank-1);

		}
		log_file << " process " << mpi_myRank << " will handle ["<<restore_start_id<<":"<<restore_end_id<< "[" << endl;



		reorganize_memory_for_assigned_blocks(restore_start_id, restore_end_id);


	}


	//!-----------------------------------------------------------

	//! loop across all files to restore
	for(INT32 restore_process_id = restore_start_id;  restore_process_id<restore_end_id; restore_process_id++)
	{
	
		log_file << " **********************************" << endl;
		log_file << " * Restoring data of process " << restore_process_id <<  endl;
		log_file << " * This Process: " << mpi_myRank <<  endl;
		log_file << " **********************************" << endl;

		//! - carefull, id_BTotal and id_BEven may have the same name
		//! - restore of id_BTotal should not be required
		read_Field3D(  id_BEven, restore_process_id);
		read_Field3D( id_EField, restore_process_id);
                read_Field3D( id_PEtotal, restore_process_id);

		read_Field3D(  id_PhiDC, restore_process_id);
		
		for(INT32 average_field=0; average_field<num_average_fields; average_field++)
		read_Field3D(id_average_Field1 +average_field, restore_process_id);
			
		read_Particle(restore_process_id);
	
	}

	if(mpi_num_processes < number_of_processes_in_BlockFile)
	reset_responsible_mpi_process(number_of_processes_in_BlockFile);

	//! build BTotal rather than restoring it
	//! (done in post mesh refinement)
	build_B_total();
	
	log_file << " done." <<endl;
	synchronize_allProcesses_MPI();
	
	return true;

}


//!------------------------------------------------------------
//!- write_BlockFile:
//!------------------------------------------------------------
void CHybrid::write_BlockFile(void)
{


	system("mkdir State");
	system("mkdir Last_State");


	INT32 level;
	CBlock *temp_Block;

	bool has_child[8];
	INT32 mpi_process[8];
	INT32 *root_mpi_process = new INT32[num_root_blocks];

	memset(has_child,        0,               8*sizeof(bool));
	memset(mpi_process,      0,               8*sizeof(INT32));
	memset(root_mpi_process, 0, num_root_blocks*sizeof(INT32));
	

	char filename[200];
	sprintf(filename,"State/State_Blocks_p%05d",mpi_myRank);
	
	ofstream Block_File;
	Block_File.open(filename, ofstream::binary);
	

	
	
	memset(num_Blocks_transfered_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
	log_file <<  endl;
	log_file << " --------------------------------------------------" << endl;
        log_file << " <<<<<<<<<<<<<<<<Writing Block Tree>>>>>>>>>>>>>>>>" << endl;
	

	//! write update TL into Blockfile
	Block_File.write(reinterpret_cast<char*> (&TL), sizeof(INT32));


	//! at store mpi_processes of root array:
	for(INT32 blk=0; blk<num_root_blocks; blk++)
	root_mpi_process[blk] = Root_Block_Array[blk].responsible_mpi_process;

	//! write processes
	Block_File.write(reinterpret_cast<char*>   (root_mpi_process), num_root_blocks*sizeof(INT32));
	num_Blocks_transfered_L[0] = num_root_blocks;

	for(level=0; level<MAX_LEVEL; level++)
	{
	
		//! - set bool true when child active
		//! - doing climbing to MAX_LEVEL-1 is sufficient
		//! - count number of active children
		temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			memset(has_child,   0, 8*sizeof(bool));
			memset(mpi_process, 0, 8*sizeof(INT32));

			for(INT32 oct=0; oct<8; oct++)
			 if(temp_Block->child_array[oct])
			 {
		
				has_child[oct] = true;
				mpi_process[oct] = temp_Block->child_array[oct]->responsible_mpi_process;

				num_Blocks_transfered_L[level+1]++;
			 }
	
			Block_File.write(reinterpret_cast<char*>   (has_child), 8*sizeof(bool));
			Block_File.write(reinterpret_cast<char*> (mpi_process), 8*sizeof(INT32));
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	
	Block_File.flush();
	Block_File.close();

	//! show information
	INT32 total_blocks_written = 0;
   	for(INT32 level=0; level<=MAX_LEVEL; level++)
	total_blocks_written +=num_Blocks_transfered_L[level];

	
	delete[] root_mpi_process;


	log_file << " Block Tree written ..." << endl;
	log_file << "  Info (global==local, Blocks exits on each proc)" << endl;
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	log_file << " -> L" << level << ":  " << num_Blocks_transfered_L[level] << endl;
	log_file << " total : " << total_blocks_written  << endl;
	log_file << " --------------------------------------------------" << endl << endl;


	
}



//!------------------------------------------------------------
//!- read_BlockFile:
//!------------------------------------------------------------
bool CHybrid::read_BlockFile(void)
{

	INT32 level;
	bool has_child[8];
	INT32 mpi_process[8];
	INT32 *root_mpi_process = new INT32[num_root_blocks];

	memset(has_child,        0,               8*sizeof(bool));
	memset(mpi_process,      0,               8*sizeof(INT32));
	memset(root_mpi_process, 0, num_root_blocks*sizeof(INT32));
	
	
	char filename[200];
// 	sprintf(filename,"State/State_Blocks_p%05d",mpi_myRank);

	
	sprintf(filename,"State/State_Blocks_p00000");
	
	
	ifstream Block_File;
	Block_File.open(filename, ifstream::binary);
	
	
	
	
	if(!Block_File)
	{

		log_file << "no Block File found." << endl;
		return false;

	}
	
	log_file << " found." << endl;
	//! first sizeof(INT32) BYTES are TL

	//! read time level
	Block_File.read(reinterpret_cast<char*> (&TL), sizeof(INT32));
	
	//! in order to build an average time of time level, time since prog start
	//! has to be devided by (TL-TL_at_last_restore_state)
	TL_at_last_restore_state = TL;
	
	log_file << " Continuing at TL " << TL << "." << endl << endl;
	
	log_file << " --------------------------------------------------" << endl;
	log_file << " >>>>>>>>>>>>>>>>Restoring Block Tree<<<<<<<<<<<<<<" << endl;
	

	log_file << "  -> restoring root_mpi_process" << endl << endl;

	//! restore mpi_processes of root array:
	Block_File.read(reinterpret_cast<char*> (root_mpi_process), num_root_blocks*sizeof(INT32));

	//! alloc root array using restored mpi processes
	init_RootBlocks_set_mpi_process(root_mpi_process);


	log_file << "  ->  restoring children " << endl;

	for(level=0; level<MAX_LEVEL; level++)
	{

	
		log_file << "  ->  L" << level << endl;

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			memset(has_child,   0, 8*sizeof(bool));
			memset(mpi_process, 0, 8*sizeof(INT32));

			Block_File.read(reinterpret_cast<char*>   (has_child), 8*sizeof(bool));
			Block_File.read(reinterpret_cast<char*> (mpi_process), 8*sizeof(INT32));
	
			for(INT32 oct=0; oct<8; oct++)
			 if(has_child[oct])
			 temp_Block->refine_Oct_set_mpi_process(oct, mpi_process[oct]);


			temp_Block = temp_Block->next_Blk_of_BlockList;
		}

	}
	
	
	Block_File.close();
	

	delete[] root_mpi_process;



	synchronize_allProcesses_MPI();





	//! Note:
	//! no MPI reduce required here since all blocks exist on
	//! every process
	
	//! Note:
	//! total_active_Blocks is incresed in refine_Oct(oct), so no additional
	//! summation is required here
	log_file << " Block Tree reconstructed ..." << endl;
	log_file << "  Info (global==local, Blocks exits on each proc)" << endl;
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	log_file << " -> L" << level << ":  " << total_Blocks_L[level] << endl;
	log_file << " total active Blocks: " << total_active_Blocks  << endl;
	log_file << " --------------------------------------------------" << endl << endl;
	log_file << endl;


	

	return true;


}

//!------------------------------------------------------------
//!- write_Field3D:
//!------------------------------------------------------------
void CHybrid::write_Field3D(INT32 FType, INT32 TL, bool write_state)
{



	//! In order to avoid that all processes write at the same time,
	//! let p0 write first, than send message to p1 to trigger writing,
	//! then p2 .... up to pN
	INT32 tag = 0;
	INT32 is_finished = 0;
	INT32 length_in_byte = 1;


	//! - every process must receive except for first process
	//! - use Blocking receive
	if(serialize_writing_of_StateFile && mpi_myRank)
	mpi_myComm.Recv(&is_finished, length_in_byte, MPI_INT32, mpi_myRank-1, tag );



	INT32 level;
	CBlock *temp_Block;
	
	char buffer[200];
	
	ofstream Field3D_File;
	
	if(write_state)
	sprintf(buffer,"State/State_%s_p%05d",Field_Name[FType], mpi_myRank);
	else
	sprintf(buffer,"%s_of_%s_TL%d.dat",Field_Name[FType],Run_Name,TL);
	
	Field3D_File.open(buffer, ofstream::binary);
	
	memset(num_Blocks_transfered_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
	
	sprintf(buffer,"%s",Field_Name[FType]);
	
	log_file << " -------------------------------" << endl;
	log_file << " <<<<<<<<Writing "<< buffer << ">>>>>>>>>>" << endl;

	for(level=0; level<=MAX_LEVEL; level++)
	{
		temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
	
	
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{
	
				Field3D_File.write(reinterpret_cast<char*>
					(temp_Block->Field_Type[FType]),
					COMPs_FType[FType]*num_nodes_in_block*sizeof(D_REAL));
				num_Blocks_transfered_L[level]++;
			}
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	
	Field3D_File.flush();
	Field3D_File.close();


	is_finished = 1;

	//! - every process must send except for last process
	//! - use Blocking Send
	if(serialize_writing_of_StateFile && mpi_myRank < mpi_num_processes-1)
	{	
		//! wait time in micro seconds
		usleep(TIME_TO_WAIT);
		mpi_myComm.Send(&is_finished, length_in_byte, MPI_INT32, mpi_myRank+1, tag);

	}


	log_file << " Block's " << buffer << " written ..." << endl;
	//! show information
	INT64 local_info_values[(MAX_LEVEL+1)];
	stringstream info_names[INFO_ARRAY_SIZE];

	//! blocks in respective level:
   	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		info_names[level] << " -> L" << level <<":  ";
		local_info_values[level] = num_Blocks_transfered_L[level];
	}

	show_information(local_info_values,
			 info_names,
			 (MAX_LEVEL+1), BUILD_SUM);

	INT32 total_Blocks_transfered = 0;

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	total_Blocks_transfered +=  local_info_values[level];
	log_file << " total: " << total_Blocks_transfered << endl;
	log_file << " -----------------------------" << endl << endl;

	
}


//!------------------------------------------------------------
//!- write_Particle:
//!------------------------------------------------------------
void CHybrid::write_Particle()
{

	INT32 level;
	CBlock *temp_Block;
	
	log_file << endl << " writing Particle ..." << endl << endl;
	
	log_file << " *******************************" << endl;
	log_file << " Particle to write (all species):" <<  endl;

	//! show information
	INT64 local_info_values[(MAX_LEVEL+1)];
	stringstream info_names[INFO_ARRAY_SIZE];

	//! particles in respective level:
   	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		info_names[level] << " -> Particle in L" << level << ":  ";
		local_info_values[level] = num_total_particles_in_L[level];
	}

	show_information(local_info_values,
			 info_names,
			 (MAX_LEVEL+1), BUILD_SUM);

	INT64 total_particle_to_write = 0;

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	total_particle_to_write +=  local_info_values[level];
	
	log_file << " -> total Particles: " << total_particle_to_write << endl;
	log_file << " *******************************" << endl << endl;




   for(INT32 species=0; species<num_Charged_Species; species++)
   {



		//! ---------------- Serialize Code -------------------------------
		//! In order to avoid that all processes write at the same time,
		//! let p0 write first, than send message to p1 to trigger writing,
		//! then p2 .... up to pN
		INT32 tag = 0;
		INT32 is_finished = 0;
		INT32 length_in_byte = 1;

		//! - every process must receive except for first process
		//! - use Blocking receive
		if(serialize_writing_of_StateFile && mpi_myRank)
		mpi_myComm.Recv(&is_finished, length_in_byte, MPI_INT32, mpi_myRank-1, tag );
		//! --------------------------------------------------------------

		char buffer[200];
	
		//! Particle_File is global ofstream
		//! (defined at top of this file)
		sprintf(buffer,"State/State_PartSpecies%d_p%05d", species, mpi_myRank);
		Particle_File.open(buffer, ofstream::binary);
		
		memset(num_Blocks_transfered_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
		memset(num_particle_processed_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
		
		log_file << " -----------------------------------" << endl;
		log_file << " <<<<<<<<<<<<<Species "<< species << ">>>>>>>>>>>>>>" << endl;
		for(level=0; level<=MAX_LEVEL; level++)
		{
			temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{
		
	
				if(mpi_myRank == temp_Block->responsible_mpi_process)
				{
		
					//! TODO:
					//! In general it should be sufficient 
					//! only to use top level Blks
					//! Each oct has to be checked so for sake of simplicity
					//! write all Blocks
	// 				if(!temp_Block->child_array)
					{
						write_Particle_of_Block(species, temp_Block);
						num_Blocks_transfered_L[level]++;
					}
				}
		
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}
		}
		
		Particle_File.flush();
		Particle_File.close();


		//! ---------------- Serialize Code ----------------------------
		//! - every process must send except for last process
		//! - use Blocking Send
		is_finished = 1;
		if(serialize_writing_of_StateFile && mpi_myRank < mpi_num_processes-1)
		{
			//! wait time in micro seconds
			usleep(TIME_TO_WAIT);
			mpi_myComm.Send(&is_finished, length_in_byte, MPI_INT32, mpi_myRank+1, tag);
		}
		//! ------------------------------------------------------------
	


		//! show information
		INT64 local_info_values[(MAX_LEVEL+1)];
		stringstream info_names[INFO_ARRAY_SIZE];
	
		log_file << " -> Particle written: " << endl;
		//! blocks in respective level:
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{
	
			info_names[level] << " -> L"<< level <<":  ";
			local_info_values[level] = num_particle_processed_L[level];
		}
	
		show_information(local_info_values,
				info_names,
				(MAX_LEVEL+1), BUILD_SUM);
		
	
		INT64 total_particle_written = 0;
	
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		total_particle_written +=  local_info_values[level];
	
		log_file << " -> Particle in species "<< species <<": " << total_particle_written << endl;
		log_file << " -----------------------------------" << endl << endl;
	
	
	}

	
}

//!------------------------------------------------------------
//!- write_Particle_of_Block:
//!------------------------------------------------------------
void CHybrid::write_Particle_of_Block(INT32 species, CBlock *active_block)
{


	for(INT32 cell=0; cell < num_nodes_in_block; cell++)
	{

		Particle_File.write(reinterpret_cast<char*>
				    (active_block->num_MPiC[species]+cell),
				    sizeof(INT32)
				   );

		Particle_File.write(reinterpret_cast<char*>
				    (active_block->pArray[species][cell]),
				     active_block->num_MPiC[species][cell] *sizeof(particle)
				   );


		num_particle_processed_L[active_block->RLevel]+=active_block->num_MPiC[species][cell];

	}


}



//!------------------------------------------------------------
//!- read_Field3D:
//!------------------------------------------------------------
void CHybrid::read_Field3D(INT32 FType, INT32 restorer_mpi_process)
{
	
	INT32 level;
	CBlock *temp_Block;
	

	bool missmatch_in_file_size = false;
	
	char buffer[200];
	sprintf(buffer,"State/State_%s_p%05d", Field_Name[FType], restorer_mpi_process);
	
	ifstream Field3D_File;
	Field3D_File.open(buffer, ifstream::binary);
	
	sprintf(buffer,"%s",Field_Name[FType]);
	if(!Field3D_File)
	{
		log_file << "no "<< buffer <<" File found." << endl;
// 		return;
		//! Do not return, otherwise some processes will
		//! reach collective mpi messages below but
		//! others won't -> deadlock
		//! (may only happen when number of procs was increased after
		//! (restore -> some procs dont have Fields/Particles till then)
	}
	else
	{
	
		memset(num_Blocks_transfered_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
		
		
		log_file << " -----------------------------" << endl;
		log_file << " >>>>>>Restoring "<< buffer << "<<<<<<<<<" << endl;
		for(level=0; level<=MAX_LEVEL; level++)
		{
			temp_Block = BlockList_of_Lev[level];
			while(temp_Block)
			{
		
				if(restorer_mpi_process == temp_Block->responsible_mpi_process)
				{
		
					//! check whether end of file is already reached
					if(Field3D_File.eof())
					{
						log_file     << ERROR_SMALL_STATEFILE << endl;
						log_file << ERROR_SMALL_STATEFILE << endl;
						missmatch_in_file_size = true;

					}
		
		
					Field3D_File.read(reinterpret_cast<char*> 
						(temp_Block->Field_Type[FType]),
						COMPs_FType[FType]*num_nodes_in_block*sizeof(D_REAL));
					num_Blocks_transfered_L[level]++;
				}
		
				temp_Block = temp_Block->next_Blk_of_BlockList;
			}
		}
		
	
		//! check whether end of file is reached
		double final;
		Field3D_File.read(reinterpret_cast<char*> (&final), sizeof(D_REAL));
	
		if(!Field3D_File.eof())
		{
	
			log_file << ERROR_LARGE_STATEFILE << endl;
			missmatch_in_file_size = true;

		}
	
	
		
		Field3D_File.close();

	}

	log_file << " checking for missmatch in file size  ..." << endl;
	log_file << " -> (will be marked as 'NaN')" << endl;
	check_for_NaN_MPI(missmatch_in_file_size);
	log_file << " done." << endl;




	log_file << " Block's " << buffer << " restored ..." << endl;
	//! show information
	INT64 local_info_values[(MAX_LEVEL+1)];
	stringstream info_names[INFO_ARRAY_SIZE];

	//! blocks in respective level:
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		info_names[level] << " -> L" << level << ":  ";
		local_info_values[level] = num_Blocks_transfered_L[level];
	}

	show_information(local_info_values,
			info_names,
			(MAX_LEVEL+1), BUILD_SUM);
	

	INT32 total_Blocks_transfered = 0;

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	total_Blocks_transfered +=  local_info_values[level];
	log_file << " total: " << total_Blocks_transfered << endl;
	log_file << " -----------------------------" << endl << endl;



	
}

//!------------------------------------------------------------
//!- read_Particle:
//!------------------------------------------------------------
void CHybrid::read_Particle(INT32 restorer_mpi_process)
{

	INT32 level;
	CBlock *temp_Block;
	
	char buffer[200];
	num_pArrayReallocs = 0;

	bool missmatch_in_file_size = false;

	log_file << endl;
	log_file << " -----------------------------"<< endl;
	log_file << " >>>>>>Restoring Particle<<<<<" << endl;
	log_file << " Modus: " << endl;
	switch(convert_TRACK_PARTICLE_State)
	{

		case 0: log_file << " ->LEAVE_UNCHANGED" << endl;
			 break;

		case 1: log_file << " ->CONVERT_TO_TRACK" << endl;
			 break;

		case 2: log_file << " ->CONVERT_TO_ORDINARY" << endl;
			 break;


	}


	for(INT32 species=0; species<num_Charged_Species; species++)
	{
		//! Particle_inFile is global ifstream
		//! (defined at top of this file)
		sprintf(buffer,"State/State_PartSpecies%d_p%05d", species, restorer_mpi_process);
		Particle_inFile.open(buffer, ofstream::binary);
		if(!Particle_inFile)
		{
			log_file << "no Particle Species "<< species << " File found." << endl;
// 			return;
			//! Do not return, otherwise some processes will
			//! reach collective mpi messages below but
			//! others won't -> deadlock
			//! (may only happen when number of procs was increased after
			//! (restore -> some procs dont have Fields/Particles till then)
		}
		else
		{
		
		
			memset(num_Blocks_transfered_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
			memset(num_particle_processed_L, 0,(MAX_LEVEL+2)*sizeof(INT64));
			
			
		// 	log_file << " restoring Particle of species " << species << " ..." << endl;
			for(level=0; level<=MAX_LEVEL; level++)
			{
				temp_Block = BlockList_of_Lev[level];
				while(temp_Block)
				{
			
					if(restorer_mpi_process == temp_Block->responsible_mpi_process)
					{
			
						//! TODO:
						//! In general it should be sufficient 
						//! only to use top level Blks, however 
						//! not all particles are written down 
						//! when doing so !
						//! -> some refined Blks contain particle ?!?
		// 				if(!temp_Block->child_array)
						{

							switch(convert_TRACK_PARTICLE_State)
							{

								case 0: if(!read_Particle_of_Block(species, temp_Block))
									 missmatch_in_file_size = true;
									 break;

								case 1: if(!read_Particle_of_Block_convert_to_TRACK_particle(species, temp_Block))
									 missmatch_in_file_size = true;
									 break;

								case 2: if(!read_Particle_of_Block_convert_to_ORDINARY_particle(species, temp_Block))
									 missmatch_in_file_size = true;
									 break;
	
	
							}

							num_Blocks_transfered_L[level]++;
						}
					}
			
					temp_Block = temp_Block->next_Blk_of_BlockList;
				}
			}
			
			//! check whether end of file is reached
			double final;
			Particle_inFile.read(reinterpret_cast<char*> (&final), sizeof(D_REAL));
		
			if(!Particle_inFile.eof())
			{
		
				log_file     << ERROR_LARGE_STATEFILE << endl;
				log_file << ERROR_LARGE_STATEFILE << endl;
				missmatch_in_file_size = true;
// 				exit(1);
			}
		
			
			Particle_inFile.close();

		}


		log_file << " checking for missmatch in file size  ..." << endl;
		log_file << " -> (will be marked as 'NaN')" << endl;
		check_for_NaN_MPI(missmatch_in_file_size);
		log_file << " done." << endl;



		//! set restored to true, otherwise box will be filled
		particle_restored = true;


	
	// 	log_file << "num_pArrayReallocs: "<< num_pArrayReallocs << endl;
		
	
	// 	for(INT32 level=0; level<=MAX_LEVEL; level++)
	// 	log_file << " Top Level Blks in L" << level <<":  "
	// 	         << num_Blocks_transfered_L[level] << endl;
	// 	log_file << endl;
		
	
		log_file << " -----------------------------" << endl;
		log_file << " >>>>>>>>>>Species "<< species << "<<<<<<<<<<<" << endl;
		log_file << "  particle read: " << endl;
	
		//! show information
		INT64 local_info_values[(MAX_LEVEL+1)];
		stringstream info_names[INFO_ARRAY_SIZE];
	
		//! blocks in respective level:
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		{
	
			info_names[level] << " -> L" << level << ":  ";
			local_info_values[level] = num_particle_processed_L[level];
		}
	
		show_information(local_info_values,
				info_names,
				(MAX_LEVEL+1), BUILD_SUM);
		
	
		INT64 total_particle_read = 0;
	
		for(INT32 level=0; level<=MAX_LEVEL; level++)
		total_particle_read +=  local_info_values[level];
		log_file << " Particle in species "<<species<< ": " << total_particle_read  << endl;
		log_file << " -----------------------------" << endl << endl;




    	}


	log_file << " *****************************" << endl;
	log_file << " Particle read (all species): " << endl;
	//! show information
	INT64 local_info_values[(MAX_LEVEL+1)];
	stringstream info_names[INFO_ARRAY_SIZE];

	//! particles in respective level:
   	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		info_names[level] << " -> Particle in L" << level << ":  ";
		local_info_values[level] = num_total_particles_in_L[level];
	}

	show_information(local_info_values,
			 info_names,
			 (MAX_LEVEL+1), BUILD_SUM);

	INT64 total_particle_read = 0;

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	total_particle_read +=  local_info_values[level];

	log_file << " -> total Particles: " << total_particle_read << endl;
	log_file << " *****************************" << endl << endl;


	
}


//!------------------------------------------------------------
//!- read_Particle_of_Block:
//!------------------------------------------------------------
bool CHybrid::read_Particle_of_Block(INT32 species, CBlock *active_block)
{


	for(INT32 cell=0; cell < num_nodes_in_block; cell++)
	{

		//! check whether end of file is already reached
		if(Particle_inFile.eof())
		{
			log_file     << ERROR_SMALL_STATEFILE << endl;
			log_file << ERROR_SMALL_STATEFILE << endl;
			return false;
// 			exit(1);
		}


		Particle_inFile.read( reinterpret_cast<char*>
				      (active_block->num_MPiC[species]+cell),
				      sizeof(INT32)
				     );

		//! increase size/set new size if required
		if(active_block->size_of_pArray[species][cell] < active_block->num_MPiC[species][cell])
		active_block->resize_pArray(species, cell, int( min_pArray_Size_factor *active_block->num_MPiC[species][cell]));


		Particle_inFile.read(reinterpret_cast<char*>
				    (active_block->pArray[species][cell]),
				     active_block->num_MPiC[species][cell] *(sizeof(particle))
				   );

		num_total_particles += 				  active_block->num_MPiC[species][cell];
		num_total_particles_in_L[active_block->RLevel] += active_block->num_MPiC[species][cell];
		num_particle_processed_L[active_block->RLevel] += active_block->num_MPiC[species][cell];
	

	}


	return true;
}


//!----------------------------------------------------------------------------------
//!- Resubmit job if the simulation is not finished at the end of the walltime
//!----------------------------------------------------------------------------------
void CHybrid::resubjob(void)
{
	 INT32 stopping_tl=TL_MAX+1;
	 stringstream info_resubmit[2];
	 INT64 local_info_values[2];
 
	log_file << endl;
	info_resubmit[0] << " Time until walltime ends [min]: ";
	info_resubmit[1] << " TL until walltime ends: ";
 
	INT32 TL_processed = (TL -TL_at_last_restore_state);
	//! measure mpi time
	double active_time_mpi = MPI_Wtime();
	run_time_mpi = active_time_mpi -run_start_mpi;
	double average_time_mpi = run_time_mpi / TL_processed;
	
	if (resubmit_enable)
	{
		synchronize_allProcesses_MPI();
	
		if (run_time_mpi >= resubmit_walltime-resubmit_security_factor)
		{
			stopping_tl=TL;
			info_resubmit[0] << " Saving and stopping simulation because of walltime ending ";
		} 
		local_info_values[0]=1.*(resubmit_walltime-run_time_mpi)/60.;
		local_info_values[1]=1.*(resubmit_walltime-run_time_mpi)/(1.*average_time_mpi);

		//!TODO make nicer
		show_information(local_info_values,
		    info_resubmit,
		    2, BUILD_MAX);
		
		synchronize_allProcesses_MPI();
	
		if (TL==stopping_tl)
		{
			if(TL_SAVE_STATE!=-1)
			save_state();
			
			if(!mpi_myRank)
			system(resubmit_jobscript); 
				
			exit(1); 
		} 
	} 
}

//!------------------------------------------------------------
//!- read_Particle_of_Block:
//!------------------------------------------------------------
bool CHybrid::read_Particle_of_Block_convert_to_TRACK_particle(INT32 species, CBlock *active_block)
{

#ifdef TRACK_PARTICLE


	ORDINARY_particle* ORDINARY_particle_Array = NULL;

	for(INT32 cell=0; cell < num_nodes_in_block; cell++)
	{

		//! check whether end of file is already reached
		if(Particle_inFile.eof())
		{
			log_file     << ERROR_SMALL_STATEFILE << endl;
			log_file << ERROR_SMALL_STATEFILE << endl;
			return false;
// 			exit(1);
		}


		Particle_inFile.read( reinterpret_cast<char*>
				      (active_block->num_MPiC[species]+cell),
				      sizeof(INT32)
				     );




		//! increase size/set new size if required
		if(active_block->size_of_pArray[species][cell] < active_block->num_MPiC[species][cell])
		active_block->resize_pArray(species, cell, int( min_pArray_Size_factor *active_block->num_MPiC[species][cell]));



		if(active_block->size_of_pArray[species][cell]>0)
		{
			ORDINARY_particle_Array = new ORDINARY_particle[active_block->size_of_pArray[species][cell]];
	
			Particle_inFile.read(reinterpret_cast<char*>
					(ORDINARY_particle_Array),
					active_block->num_MPiC[species][cell] *(sizeof(ORDINARY_particle))
					);
	
			for(INT32 part=0; part<active_block->size_of_pArray[species][cell]; part++)
			{
	
	
				active_block->pArray[species][cell][part].rel_r[0] = ORDINARY_particle_Array[part].rel_r[0];
				active_block->pArray[species][cell][part].rel_r[1] = ORDINARY_particle_Array[part].rel_r[1];
				active_block->pArray[species][cell][part].rel_r[2] = ORDINARY_particle_Array[part].rel_r[2];
	
				active_block->pArray[species][cell][part].v[0] = ORDINARY_particle_Array[part].v[0];
				active_block->pArray[species][cell][part].v[1] = ORDINARY_particle_Array[part].v[1];
				active_block->pArray[species][cell][part].v[2] = ORDINARY_particle_Array[part].v[2];
	
				active_block->pArray[species][cell][part].weight = ORDINARY_particle_Array[part].weight;
	
	
				active_block->pArray[species][cell][part].number   =  -1;
				active_block->pArray[species][cell][part].exchange =  -1;
	
	
			}
	
			delete[] ORDINARY_particle_Array;
	
	
			num_total_particles += 				  active_block->num_MPiC[species][cell];
			num_total_particles_in_L[active_block->RLevel] += active_block->num_MPiC[species][cell];
			num_particle_processed_L[active_block->RLevel] += active_block->num_MPiC[species][cell];
		}
	

	}

	return true;

#endif

}


//!------------------------------------------------------------
//!- read_Particle_of_Block:
//!------------------------------------------------------------
bool CHybrid::read_Particle_of_Block_convert_to_ORDINARY_particle(INT32 species, CBlock *active_block)
{


	TRACK_particle* TRACK_particle_Array = NULL;

	for(INT32 cell=0; cell < num_nodes_in_block; cell++)
	{

		//! check whether end of file is already reached
		if(Particle_inFile.eof())
		{
			log_file     << ERROR_SMALL_STATEFILE << endl;
			log_file << ERROR_SMALL_STATEFILE << endl;
			return false;
// 			exit(1);
		}


		Particle_inFile.read( reinterpret_cast<char*>
				      (active_block->num_MPiC[species]+cell),
				      sizeof(INT32)
				     );

			//! increase size/set new size if required
			if(active_block->size_of_pArray[species][cell] < active_block->num_MPiC[species][cell])
			active_block->resize_pArray(species, cell, int( min_pArray_Size_factor *active_block->num_MPiC[species][cell]));
	
			if(active_block->size_of_pArray[species][cell]>0)
			{
	
			TRACK_particle_Array = new TRACK_particle[active_block->size_of_pArray[species][cell]];
	
			Particle_inFile.read(reinterpret_cast<char*>
					(TRACK_particle_Array),
					active_block->num_MPiC[species][cell] *(sizeof(TRACK_particle))
					);
	
			for(INT32 part=0; part<active_block->size_of_pArray[species][cell]; part++)
			{
	
	
				active_block->pArray[species][cell][part].rel_r[0] = TRACK_particle_Array[part].rel_r[0];
				active_block->pArray[species][cell][part].rel_r[1] = TRACK_particle_Array[part].rel_r[1];
				active_block->pArray[species][cell][part].rel_r[2] = TRACK_particle_Array[part].rel_r[2];
	
				active_block->pArray[species][cell][part].v[0] = TRACK_particle_Array[part].v[0];
				active_block->pArray[species][cell][part].v[1] = TRACK_particle_Array[part].v[1];
				active_block->pArray[species][cell][part].v[2] = TRACK_particle_Array[part].v[2];
	
				active_block->pArray[species][cell][part].weight = TRACK_particle_Array[part].weight;

	
			}
	
			delete[] TRACK_particle_Array;
	
	
			Particle_inFile.read(reinterpret_cast<char*>
					(active_block->pArray[species][cell]),
					active_block->num_MPiC[species][cell] *(sizeof(particle))
					);
	
			num_total_particles += 				  active_block->num_MPiC[species][cell];
			num_total_particles_in_L[active_block->RLevel] += active_block->num_MPiC[species][cell];
			num_particle_processed_L[active_block->RLevel] += active_block->num_MPiC[species][cell];

		}
	

	}


	return true;
}

//! Get_Mesh
//! Returns the mesh in CESQ Coordinate system

void CHybrid::get_Mesh(void)
{
    
    log_file << " Save Mesh" << endl;
    //!-----------------------------------------//
    //!     Prepare                             //
    //!----------------------------------------//
    char buffer2[200];
    sprintf(buffer2,"mkdir %s/mesh/temp", data_output_path);
    
    if(!mpi_myRank)
        system(buffer2);
    
    synchronize_allProcesses_MPI();
    
    char filename[200];
    ofstream mesh_outfile;
    sprintf(filename,"%s/mesh/temp/Mesh_p%05d.txt", data_output_path,mpi_myRank);
    
    mesh_outfile.open(filename);
    
    //open File
    //!-----------------------------------------//
    //!     Loop over all Blks                       //
    //!----------------------------------------//
    for(INT32 level=0; level<=MAX_LEVEL; level++)
    {
        
        CBlock* temp_Block = BlockList_of_Lev[level];
        while(temp_Block)
        {
            
            if(mpi_myRank==temp_Block->responsible_mpi_process)
                for(INT32 oct=0; oct<8; oct++)
                {
                    
                    if(!temp_Block->child_array[oct])
                                               
//                         temp_Block->inject_cometary_profile_ions(oct,species,insert_every_x_TL);
                        // Auslesen der Position und geschreiben 
                            temp_Block->getMesh(oct,mesh_outfile);
//                         mesh_outfile <<  temp_Block->origin[0] << "\t" << temp_Block->origin[1] << "\t" <<temp_Block->origin[2] << endl;
                }
                
                temp_Block = temp_Block->next_Blk_of_BlockList;
            
        }
    }
    
    synchronize_allProcesses_MPI();
 
    log_file << "     Packing files of mesh to single file... ";
    
    if(!mpi_myRank) 
    {
        double r[3];
        //! Create Single File
        ofstream outfile;
        char outfile_name[200];
        sprintf(outfile_name,"%s/mesh/mesh.txt",data_output_path);
        outfile.open(outfile_name);
        
        int counter_total = 0;
        
        for(int process=0; process<mpi_num_processes; process++)
        {
            char infile_process[200];
            sprintf(infile_process,"%s/mesh/temp/Mesh_p%05d.txt", data_output_path,process);
            
            ifstream infile;
            infile.open(infile_process);
            
            if(!infile)
            {
                log_file << "Error in 'get_Mesh':" << endl;
                log_file << "no file: " << infile_process << endl;
                log_file << "Returning ... ." << endl;
                return;
                
            }
            
            int counter_file = 0;
            for(;;)
            {
                //! get position
                infile >> r[0];
                infile >> r[1];
                infile >> r[2];
                
                if(infile.eof())
                    break;
                
                outfile << -1*r[0]*SI_x0 << "\t";
                outfile << r[1]*SI_x0 << "\t";
                outfile << r[2]*SI_x0 << endl;
                
                counter_file++;
            }
            counter_total +=counter_file;
            infile.close();
        }
        
        log_file << " done. "  << endl;
        log_file << "     ->Filename: '" << outfile_name <<"'"<< endl;
        log_file << "     ->Values read/written: " << counter_total << endl;
        
        log_file << endl;
        
        
        outfile.close();
	
    }
    
    sprintf(buffer2,"rm -r %s/mesh/temp", data_output_path);

    if(!mpi_myRank)
      system(buffer2);
	
}
