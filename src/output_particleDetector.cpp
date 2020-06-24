
#include "parameters.h"
#include "CHybrid.h"
#include "absolute_Globals.h"
#include "utils.h"

#include <iostream>
#include <fstream>
#include <cmath>


void CHybrid::particle_detector_output(INT32 id_detector)
{
  //!-----------------------------------
  //! 1) prepare output
  //! estimate field data
  //!-----------------------------------
  log_file 	<< " Prepare particle detector ..." << endl;
  char buffer[200];
  sprintf(buffer,"mkdir %s/particle_detector/temp", data_output_path);
  if(!mpi_myRank)
    system(buffer);
  
  //! wait until directory is created
  synchronize_allProcesses_MPI();
  
  //!-----------------------------------------------------------
  //!1) Perpare Write
  //!-----------------------------------------------------------
//   log_file 	<< " Open file ..." << endl;
  ofstream particle_detector_file;
  char particle_detector_file_name[200];
  sprintf(particle_detector_file_name,"%s/particle_detector/temp/%s_p%05d_TL%06d.txt", data_output_path, Detector_FileName[id_detector], mpi_myRank, TL);
  particle_detector_file.open(particle_detector_file_name);

  //!-----------------------------------
  //! 3) Loop over all blocks
  //!-----------------------------------
  log_file 	<< " Loop over all blocks ..." << endl;

  for(INT32 level=0; level<=MAX_LEVEL; level++)
    {
      CBlock* temp_Block = BlockList_of_Lev[level];
      while(temp_Block)
      {
	if(mpi_myRank == temp_Block->responsible_mpi_process)
	  temp_Block->particle_detector_output(id_detector,particle_detector_file);

	temp_Block = temp_Block->next_Blk_of_BlockList;
      }
    }
  
  //!-----------------------------------------------------------
  //!1) Close File
  //!-----------------------------------------------------------
//   log_file 	<< " Close File ..." << endl;
  particle_detector_file.close();

  //!-----------------------------------
  //! 4) Move parallel output to single sile
  //!-----------------------------------
  //! All processors must have written files before master can
  //! begin packing
  synchronize_allProcesses_MPI();
  log_file 	<< " Pack to single file ..." << endl;
  if(!mpi_myRank)
    parallel_particle_detector_to_single_file(id_detector);
  //!-----------------------------------
  //! 5) Clean up 
  //!-----------------------------------
  log_file 	<< " Clean up ..." << endl;
  sprintf(buffer,"rm -r %s/particle_detector/temp", data_output_path);

      if(!mpi_myRank)
      system(buffer);
}



void CHybrid::parallel_particle_detector_to_single_file(INT32 id_detector)
{
  log_file << " Packing particle data to single files ... " << endl;
  char outfile_name[200];
  short det_species;
  D_REAL position[3];
  D_REAL velocity[3];
  D_REAL weight;
  D_REAL speed;
  
  ofstream outfile;
  sprintf(outfile_name,"%s/particle_detector/%s_TL%06d.txt", data_output_path,  Detector_FileName[id_detector], TL);
  outfile.open(outfile_name);
  
  int counter_total = 0;
  
  for(int process=0; process<mpi_num_processes; process++)
  {


    char infile_process[200];

    sprintf(infile_process,"%s/particle_detector/temp/%s_p%05d_TL%06d.txt", data_output_path, Detector_FileName[id_detector], process, TL);

    ifstream infile;
    infile.open(infile_process);


    if(!infile)
    {
	    log_file << "Error in 'parallel_output_to_single_file_lineout':" << endl;
	    log_file << "no file: " << infile_process << endl;
	    log_file << "Returning ... ." << endl;
	    return;
	    

    }
    
    for(;;)
    {

	    //! get position
	    infile >> det_species >> weight >> position[0] >> position[1] >> position[2] >> velocity[0] >> velocity[1] >> velocity[2] >> speed;

	    if(infile.eof())
	    break;
	    
	    outfile << det_species<< "\t" << weight << "\t" << position[0]*factor_scale_mesh << "\t" << position[1]*factor_scale_mesh << "\t" << position[2]*factor_scale_mesh << "\t" 
		    << velocity[0] << "\t" << velocity[1] << "\t" << velocity[2] << "\t" << speed;
	    outfile << endl;
	    
	    counter_total++;
    }

    infile.close();

  }


  log_file << " done. "  << endl;
  log_file << "     ->Filename: '" << outfile_name <<"'"<< endl;
  log_file << "     ->Values read/written: " << counter_total << endl;

  log_file << endl;

  outfile.close();
  
}


void CBlock::particle_detector_output(INT32 id_detector, ostream &particle_detector_file)
{
  //! temprary pointer to particle in list
  particle *active_particle;
  PARTICLE_REAL x[3];
  
  for(INT32 oct=0; oct<8; oct++)
  {
    if(!child_array[oct])
    {
      //!-----------------------------------------------------------
      //!1) Move Step:
      //!   All particle weights are assumed to be > 0 now !!!
      //!-----------------------------------------------------------
      //! get index of oct (0 or 1 for each dimension)
      INT32 a = oct/4;
      INT32 b = (oct -a*4)/2;
      INT32 c = (oct -a*4 -b*2);
      
      INT32 start[3] = {!a *!is_box_boundary[0] +a*BlkNds_X/2,
		      !b *!is_box_boundary[2] +b*BlkNds_Y/2,
		      !c *!is_box_boundary[4] +c*BlkNds_Z/2};

      INT32 end[3]={BlkNds_X -a*!is_box_boundary[1] -!a*(BlkNds_X/2),
		  BlkNds_Y -b*!is_box_boundary[3] -!b*(BlkNds_Y/2),
		  BlkNds_Z -c*!is_box_boundary[5] -!c*(BlkNds_Z/2)};

      //!-----------------------------------------------------------
      //!1) Loop over all Particle
      //!-----------------------------------------------------------
      for(short species=0; species<num_Charged_Species; species++)
      {

	
	for(INT32 i = start[0]; i < end[0]; i++)
	  for(INT32 j = start[1]; j < end[1]; j++)
	    for(INT32 k = start[2]; k < end[2]; k++)
	    {

	    INT32 i_j_k  = i*BlkNds_Y*BlkNds_Z 
		    +j*BlkNds_Z 
		    +k;

	    for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k]; part_index++)
	    {

		active_particle = pArray[species][i_j_k] +part_index;

		INT32 cell_indices[3] = {i,j,k};
		intern2normedCoords(x, active_particle->rel_r, cell_indices);

		//!-----------------------------------------------------------
		//!1) Check if conditions are full filled
		//!-----------------------------------------------------------
		if(detector_box_xmin[id_detector] < x[0] && x[0] < detector_box_xmax[id_detector]) 
		  if(detector_box_ymin[id_detector] < x[1] && x[1] < detector_box_ymax[id_detector]) 
		    if(detector_box_zmin[id_detector] < x[2] && x[2] < detector_box_zmax[id_detector]) 
		    {
		      //! write species
		      particle_detector_file << species << "\t";
		      //! write weight
		      particle_detector_file << active_particle->weight << "\t";

		      //! write position in normed units
		      particle_detector_file << x[0] << "\t";
		      particle_detector_file << x[1] << "\t";
		      particle_detector_file << x[2] << "\t";


		      //! write velocity in normed units
		      particle_detector_file << active_particle->v[0] << "\t";
		      particle_detector_file << active_particle->v[1] << "\t";
		      particle_detector_file << active_particle->v[2] << "\t";

		      particle_detector_file << vec_len(active_particle->v);

		      particle_detector_file << endl;
		    }

		}//! for particle of cell
	    }//! for ijk
      }
    }
  }//! end oct
      
}