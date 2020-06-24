#include "CHybrid.h"

#include "utils.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <sstream>



//!------------------------------------------------------------
//!trace_trajectory:
//!
//! - Field_id: id of Field that shall be traced along trajectory
//! - positions: x,y,z = [0;1[ normalized position in Simu-Box
//!------------------------------------------------------------
void CHybrid::trace_trajectory_energy_spectrum(INT32 id_trajectory, INT32 Field_id, INT32 num_positions, D_REAL** positions, char** time_string_array, D_REAL** SCField)
{


	log_file <<"     Tracing energy spectrum data ...  ";

	//! indices for interpolation
	INT32  i,j,k;
	INT32  i_j_k, ip1_jp1_kp1;
	INT32  ip1_j_k, i_jp1_k,  i_j_kp1;
	INT32  ip1_jp1_k, ip1_j_kp1,  i_jp1_kp1;

	ofstream trajectory_outfile;
	INT32  *index_blk[3], *index_cell[3], blk_index;
	D_REAL *x_cell[3], r_vec[3], A[3], shape_func[8], *BlkField[3];

	CBlock *temp_Block, *top_level_Block;


	particle* active_particle;

	//! decide whether to attach data to file or open new one
	if(create_new_file_each_TL)
	{

		char filename[200];

		if(pack_parallel_output_to_single_file)
		sprintf(filename,"%s/trajectories/temp/%s_Energy_p%05d_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory], mpi_myRank, TL);
		else
		sprintf(filename,"%s/trajectories/%s_Energy_p%05d_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory],mpi_myRank, TL);


		trajectory_outfile.open(filename);
	}
	else
	{

		char filename[200];

		if(pack_parallel_output_to_single_file)
		sprintf(filename,"%s/trajectories/temp/%s_Energy_p%05d.txt", data_output_path, Trajectory_FileName[id_trajectory], mpi_myRank);
		else
		sprintf(filename,"%s/trajectories/%s_Energy_p%05d.txt", data_output_path, Trajectory_FileName[id_trajectory],  mpi_myRank);

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
// 				for(INT32 comp=0; comp<COMPs_FType[Field_id]; comp++)
// 				BlkField[comp] = top_level_Block->Field_Type[Field_id] +comp *num_nodes_in_block;
				
		
				INT32 level=top_level_Block->RLevel;
		
				//! use common i,j,k notation:
				i = index_cell[0][level];
				j = index_cell[1][level];
				k = index_cell[2][level];
								
				i_j_k   =  i*BlkNds_Y*BlkNds_Z	+	j*BlkNds_Z	+k;
				

				INT32* number_MPiC = new INT32[num_Particle_Species];
				
				D_REAL** velocity = new D_REAL*[3];
				
				
				for(INT32 species=0; species<num_Particle_Species; species++)
				{
// 					if(Ion_Charges[species]>0)
// 					{
					 
						number_MPiC[species]=top_level_Block->num_MPiC[species][i_j_k];	
					
						velocity[0] = new D_REAL[number_MPiC[species]];
						velocity[1] = new D_REAL[number_MPiC[species]];
						velocity[2] = new D_REAL[number_MPiC[species]];
						
						D_REAL* gewicht = new D_REAL[number_MPiC[species]];
						
						top_level_Block->get_velocity(velocity, species, i_j_k, gewicht);
						
						for(INT32 part_in_cell=0;part_in_cell<number_MPiC[species]; part_in_cell++)
						{
											
							if(do_read_timestring)
							trajectory_outfile << time_string_array[counter] <<"	";
				
							trajectory_outfile << positions[0][counter] *LX -Box_Origin[0] <<"		"
									<< positions[1][counter] *LY -Box_Origin[1] <<"		"
									<< positions[2][counter] *LZ -Box_Origin[2] <<"		";
								
							trajectory_outfile << velocity[0][part_in_cell] << "	  "
									<< velocity[1][part_in_cell] << "	  "
									<< velocity[2][part_in_cell] << "	  ";
									
							trajectory_outfile << gewicht[part_in_cell] << "	  ";	
							
							trajectory_outfile << species << "	  ";		
									
// 							trajectory_outfile << level << "	  ";		
									
									
							trajectory_outfile << endl;

						}
						
						//! delete coordinates memory
						delete[] velocity[0];
						delete[] velocity[1];
						delete[] velocity[2];
					
// 					}
				
					
				}

				delete[] velocity;
				delete[] number_MPiC;
				
				trajectory_outfile << endl;
// 				trajectory_outfile << endl;
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
//!trace_line:
//!
//! - Field_id: id of Field that shall be traced along line
//! - positions: x,y,z = [0;1[ normalized position in Simu-Box
//!------------------------------------------------------------
void CHybrid::trace_line_energy_spectrum(INT32 id_line, INT32 Field_id, INT32 num_positions, D_REAL** positions)
{


	log_file <<"     Tracing energy spectrum data ...  ";

	//! indices for interpolation
	INT32  i,j,k;
	INT32  i_j_k, ip1_jp1_kp1;
	INT32  ip1_j_k, i_jp1_k,  i_j_kp1;
	INT32  ip1_jp1_k, ip1_j_kp1,  i_jp1_kp1;

	ofstream lineout_outfile;
	INT32  *index_blk[3], *index_cell[3], blk_index;
	D_REAL *x_cell[3], r_vec[3], A[3], shape_func[8], *BlkField[3];

	CBlock *temp_Block, *top_level_Block;



	//! decide whether to attach data to file or open new one
	if(create_new_file_each_TL_lineout)
	{

		char filename[200];

		if(pack_parallel_output_to_single_file)
		{
		  sprintf(filename,"%s/lineout/temp/%s_Energy_p%05d_TL%06d.txt", data_output_path, LineOut_FileName[id_line], mpi_myRank, TL);
		}
		else
		{
		  sprintf(filename,"%s/lineout/%s_Energy_p%05d_TL%06d.txt", data_output_path, LineOut_FileName[id_line],  mpi_myRank, TL);
		}

		lineout_outfile.open(filename);
	}
	else
	{

		char filename[200];

		if(pack_parallel_output_to_single_file_lineout)
		{
		  sprintf(filename,"%s/lineout/temp/%s_Energy_p%05d.txt", data_output_path, LineOut_FileName[id_line],  mpi_myRank);
		}
		else
		{
		  sprintf(filename,"%s/lineout/%s_Energy_p%05d.txt", data_output_path, LineOut_FileName[id_line], mpi_myRank);
		}

		if(!TL)
		lineout_outfile.open(filename);
		else
		lineout_outfile.open(filename, ios_base::app);

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

				INT32* number_MPiC = new INT32[num_Particle_Species];
				
				D_REAL** velocity = new D_REAL*[3];
				
				
				for(INT32 species=0; species<num_Particle_Species; species++)
				{
// 					if(Ion_Charges[species]>0)
// 					{
					 
						number_MPiC[species]=top_level_Block->num_MPiC[species][i_j_k];	
					
						velocity[0] = new D_REAL[number_MPiC[species]];
						velocity[1] = new D_REAL[number_MPiC[species]];
						velocity[2] = new D_REAL[number_MPiC[species]];
						
						D_REAL* gewicht = new D_REAL[number_MPiC[species]];
						
						top_level_Block->get_velocity(velocity, species, i_j_k, gewicht);
						
						for(INT32 part_in_cell=0;part_in_cell<number_MPiC[species]; part_in_cell++)
						{
											
				
							lineout_outfile << positions[0][counter] *LX -Box_Origin[0] <<"		"
									<< positions[1][counter] *LY -Box_Origin[1] <<"		"
									<< positions[2][counter] *LZ -Box_Origin[2] <<"		";
								
							lineout_outfile << velocity[0][part_in_cell] << "	  "
									<< velocity[1][part_in_cell] << "	  "
									<< velocity[2][part_in_cell] << "	  ";
									
							lineout_outfile << gewicht[part_in_cell] << "	  ";	
							
							lineout_outfile << species << "	  ";		
									
// 							trajectory_outfile << level << "	  ";		
									
									
							lineout_outfile << endl;

						}
						
						//! delete coordinates memory
						delete[] velocity[0];
						delete[] velocity[1];
						delete[] velocity[2];
					
// 					}
				
					
				}

				delete[] velocity;
				delete[] number_MPiC;
				
				lineout_outfile << endl;
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


	lineout_outfile.flush();
	lineout_outfile.close();

	log_file <<"done." <<endl;


}


void CHybrid::parallel_output_to_single_file_energy_spectrum(INT32 id_trajectory, INT32 id_field)
{




// 	ifstream Trajectory_File;
// 	Trajectory_File.open(Trajectory_FileName[id_trajectory]);
// 
// 
// 	if(!Trajectory_File)
// 	{
// 		log_file << " no " << Trajectory_FileName[id_trajectory] << " File found." << endl;
// 		return;
// 	}
// 
// 	Trajectory_File.close();


	log_file << "     Packing energy spectrum files to single file... ";


	char buffer[200];
	double Field[3], r[3], SC_Field[3], weight;
	INT32 level, species;
	char outfile_name[200];
// 	char time_string[200];
	string time_string, time_string_temp;

	ofstream outfile;


	//! decide whether to attach data to file or open new one
	if(create_new_file_each_TL)
	{
		sprintf(outfile_name,"%s/trajectories/%s_Energy_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory], TL);
		outfile.open(outfile_name);


	}
	else
	{

		sprintf(outfile_name,"%s/trajectories/%s_Energy.txt", data_output_path, Trajectory_FileName[id_trajectory]);

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
		sprintf(infile_process,"%s/trajectories/temp/%s_Energy_p%05d_TL%06d.txt", data_output_path, Trajectory_FileName[id_trajectory],  process, TL);
		else
		sprintf(infile_process,"%s/trajectories/temp/%s_Energy_p%05d.txt", data_output_path, Trajectory_FileName[id_trajectory],  process);


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

			time_string_temp = time_string;
			
			if(do_read_timestring)
			infile >> time_string;
			


			//! get position
			infile >> r[0];
			infile >> r[1];
			infile >> r[2];
			
// 			if(do_read_SpaceCraftField)
// 			{
// 				infile >> SC_Field[0];
// 				infile >> SC_Field[1];
// 				infile >> SC_Field[2];
// 			}

			if(infile.eof())
			break;
	
			if(time_string_temp != time_string)
			outfile << endl;
			
			if(do_read_timestring)
			outfile << time_string <<"	";

			r[0]*=factor_scale_mesh;
			r[1]*=factor_scale_mesh;
			r[2]*=factor_scale_mesh;

			//! write outfile
			outfile 	<< r[0] <<"		";
			outfile 	<< r[1] <<"		";
			outfile 	<< r[2] <<"		";
// 			outfile <<  vec_len(r)  <<"		";

// 			if(do_read_SpaceCraftField)
// 			{
// 				outfile 	<< SC_Field[0] <<"		";
// 				outfile 	<< SC_Field[1] <<"		";
// 				outfile 	<< SC_Field[2] <<"		";
// // 				outfile << sqrt(SC_Field[0]*SC_Field[0] +SC_Field[1]*SC_Field[1] +SC_Field[2]*SC_Field[2]) <<"   	";
// 			}
			
			
			infile >> Field[0];
			infile >> Field[1];
			infile >> Field[2];
	
			//! write BField
// 			outfile << Field[0] <<"   	";
// 			outfile << Field[1] <<"   	";
// 			outfile << Field[2] <<"   	";
			
			outfile << Field[0] * SI_v0<<"   	";
			outfile << Field[1] * SI_v0<<"   	";
			outfile << Field[2] * SI_v0<<"   	";
			
			outfile <<  Field[0]*SI_v0*Field[0]*SI_v0 +   Field[1]*SI_v0*Field[1]*SI_v0 + Field[2]*SI_v0*Field[2]*SI_v0 <<" 	  ";
			
// 				outfile << sqrt(Field[0]*Field[0] +Field[1]*Field[1] +Field[2]*Field[2]) <<"   	";
			
			infile >> weight;
			infile >> species;
// 			infile >> level;
			
			outfile << weight<<" 	  ";
			outfile << species<<" 	  ";
// 			outfile << level;
			
			outfile << endl;

		
// 			if(COMPs_FType[id_field]==1)
// 			{
// 				//! get Field
// 				infile >> Field[0];
// 		
// 				//! write BField
// 				outfile << Field[0] <<"   	";
// 				outfile << endl;
// 			}
// 
// 			if(COMPs_FType[id_field]==3)
// 			{
// 				//! get Field
// 				infile >> Field[0];
// 				infile >> Field[1];
// 				infile >> Field[2];
// 		
// 				//! write BField
// 				outfile << Field[0] <<"   	";
// 				outfile << Field[1] <<"   	";
// 				outfile << Field[2] <<"   	";
// // 				outfile << sqrt(Field[0]*Field[0] +Field[1]*Field[1] +Field[2]*Field[2]) <<"   	";
// 				outfile << endl;
// 			}
	
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


void CHybrid::parallel_output_to_single_lineout_file_energy_spectrum(INT32 id_line, INT32 id_field)
{




// 	ifstream Trajectory_File;
// 	Trajectory_File.open(Trajectory_FileName[id_trajectory]);
// 
// 
// 	if(!Trajectory_File)
// 	{
// 		log_file << " no " << Trajectory_FileName[id_trajectory] << " File found." << endl;
// 		return;
// 	}
// 
// 	Trajectory_File.close();


	log_file << "     Packing energy spectrum files to single file... ";


	char buffer[200];
	double Field[3], r[3], SC_Field[3], weight;
	INT32 level, species;
	char outfile_name[200];
// 	char time_string[200];
	string time_string, time_string_temp;

	ofstream outfile;


	//! decide whether to attach data to file or open new one
	if(create_new_file_each_TL_lineout)
	{
		sprintf(outfile_name,"%s/lineout/%s_Energy_TL%06d.txt", data_output_path, LineOut_FileName[id_line], TL);
		outfile.open(outfile_name);


	}
	else
	{

		sprintf(outfile_name,"%s/lineout/%s_Energy.txt", data_output_path, LineOut_FileName[id_line]);

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
		sprintf(infile_process,"%s/lineout/temp/%s_Energy_p%05d_TL%06d.txt", data_output_path, LineOut_FileName[id_line],  process, TL);
		else
		sprintf(infile_process,"%s/lineout/temp/%s_Energy_p%05d.txt", data_output_path, LineOut_FileName[id_line],  process);


		ifstream infile;
		infile.open(infile_process);


		if(!infile)
		{
			log_file << "Error in 'parallel_output_to_single_linefile':" << endl;
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
	
	

			r[0]*=factor_scale_mesh;
			r[1]*=factor_scale_mesh;
			r[2]*=factor_scale_mesh;

			//! write outfile
			outfile 	<< r[0] <<"		";
			outfile 	<< r[1] <<"		";
			outfile 	<< r[2] <<"		";
// 			outfile <<  vec_len(r)  <<"		";


			
			
			infile >> Field[0];
			infile >> Field[1];
			infile >> Field[2];
	
			//! write BField
// 			outfile << Field[0] <<"   	";
// 			outfile << Field[1] <<"   	";
// 			outfile << Field[2] <<"   	";
			
			outfile << Field[0] * SI_v0<<"   	";
			outfile << Field[1] * SI_v0<<"   	";
			outfile << Field[2] * SI_v0<<"   	";
			
			outfile <<  Field[0]*SI_v0*Field[0]*SI_v0 +   Field[1]*SI_v0*Field[1]*SI_v0 + Field[2]*SI_v0*Field[2]*SI_v0 <<" 	  ";
			
// 				outfile << sqrt(Field[0]*Field[0] +Field[1]*Field[1] +Field[2]*Field[2]) <<"   	";
			
			infile >> weight;
			infile >> species;
// 			infile >> level;
			
			outfile << weight<<" 	  ";
			outfile << species<<" 	  ";
// 			outfile << level;
			
			outfile << endl;

		
// 			if(COMPs_FType[id_field]==1)
// 			{
// 				//! get Field
// 				infile >> Field[0];
// 		
// 				//! write BField
// 				outfile << Field[0] <<"   	";
// 				outfile << endl;
// 			}
// 
// 			if(COMPs_FType[id_field]==3)
// 			{
// 				//! get Field
// 				infile >> Field[0];
// 				infile >> Field[1];
// 				infile >> Field[2];
// 		
// 				//! write BField
// 				outfile << Field[0] <<"   	";
// 				outfile << Field[1] <<"   	";
// 				outfile << Field[2] <<"   	";
// // 				outfile << sqrt(Field[0]*Field[0] +Field[1]*Field[1] +Field[2]*Field[2]) <<"   	";
// 				outfile << endl;
// 			}
	
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

//!------------------------------------------------------------------------
//!- add fields
//!------------------------------------------------------------------------
void CHybrid::calc_scalar_fields_sum(INT32 id_Result, INT32 id_Field1, INT32 id_Field2)
{


	log_file << "  build sum of scalar fields: "<< Field_Name[id_Result]
		      <<" = "<< Field_Name[id_Field1]
		      <<" + "<< Field_Name[id_Field2];
	

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{

				D_REAL *R = temp_Block->Field_Type[id_Result] ;
				D_REAL *A = temp_Block->Field_Type[id_Field1] ;
				D_REAL *B = temp_Block->Field_Type[id_Field2] ;

				for(INT32 elm=0; elm<num_nodes_in_block; elm++)
				R[elm] = A[elm] + B[elm];

			}

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	
	log_file << " done." << endl;

}

//!------------------------------------------------------------------------
//!- add fields
//!------------------------------------------------------------------------
void CHybrid::calc_vector_fields_sum(INT32 id_Result, INT32 id_Field1, INT32 id_Field2)
{


	log_file << "  build sum of vector fields: "<< Field_Name[id_Result]
		      <<" = "<< Field_Name[id_Field1]
		      <<" + "<< Field_Name[id_Field2];
	
	//! in order to advance odd BField derivatives of
	//! even BField are necessary. So first cp even GC.

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{

				D_REAL *R1 = temp_Block->Field_Type[id_Result] +0*num_nodes_in_block;
				D_REAL *R2 = temp_Block->Field_Type[id_Result] +1*num_nodes_in_block;
				D_REAL *R3 = temp_Block->Field_Type[id_Result] +2*num_nodes_in_block;

				D_REAL *A1 = temp_Block->Field_Type[id_Field1] +0*num_nodes_in_block;
				D_REAL *A2 = temp_Block->Field_Type[id_Field1] +1*num_nodes_in_block;
				D_REAL *A3 = temp_Block->Field_Type[id_Field1] +2*num_nodes_in_block;

				D_REAL *B1 = temp_Block->Field_Type[id_Field2] +0*num_nodes_in_block;
				D_REAL *B2 = temp_Block->Field_Type[id_Field2] +1*num_nodes_in_block;
				D_REAL *B3 = temp_Block->Field_Type[id_Field2] +2*num_nodes_in_block;


				for(INT32 elm=0; elm<num_nodes_in_block; elm++)
				{

					R1[elm] = A1[elm] + B1[elm];
					R2[elm] = A2[elm] + B2[elm];
					R3[elm] = A3[elm] + B3[elm];
	
				}

			}

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	
	log_file << " done." << endl;

}

