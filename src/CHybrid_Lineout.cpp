
#include "CHybrid.h"
#include "output_silo.h"
#include "utils.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <cmath>

//!-----------------------------------------------------------
//! get_positions_from_analytical_formula
//! The methode generates the positions of the line
//!-----------------------------------------------------------

void CHybrid::get_positions_from_analytical_formula(INT32 id_line, INT32 &num_positions, D_REAL** &positions)
{

	//!-----------------------------
	//!----Line Out Config----------
	//!-----------------------------
	D_REAL start_pos[3];
	D_REAL line_vec[3];
	D_REAL stop_pos[3];
	D_REAL delta_r[3];
	
	//! --- X-Achse
	if(id_line==0) {
	  //! Start Position
	  start_pos[0] = -Box_Origin[0];
	  start_pos[1] = 0.;
	  start_pos[2] = 0.;
	  //! Direction Vector
	  //! normalised
	  line_vec[0] = 1.;
	  line_vec[1] = 0.;
	  line_vec[2] = 0.;
	  //! Stop Position
	  stop_pos[0] = (LX-Box_Origin[0])-1;
	  stop_pos[1] = 0.;
	  stop_pos[2] = 0.;
	}
	
		//! --- 300kmX-Achse
	if(id_line==1) {
	  //! Start Position
	  start_pos[0] = -Box_Origin[0];
	  start_pos[1] = 1.38;
	  start_pos[2] = 0.;
	  //! Direction Vector
	  //! normalised
	  line_vec[0] = 0.;
	  line_vec[1] = 0.;
	  line_vec[2] = 1.;
	  //! Stop Position
	  stop_pos[0] = (LX-Box_Origin[0])-1;
	  stop_pos[1] = 1.38;
	  stop_pos[2] = 0.;
	}
	
	//! --- Y-Achse
	if(id_line==3) {
	  //! Start Position
	  start_pos[0] = 0.;
	  start_pos[1] = -Box_Origin[1];
	  start_pos[2] = 0.;
	  //! Direction Vector
	  //! normalised
	  line_vec[0] = 0.;
	  line_vec[1] = 1.;
	  line_vec[2] = 0.;
	  //! Stop Position
	  stop_pos[0] = 0.;
	  stop_pos[1] = (LX-Box_Origin[1])-1;
	  stop_pos[2] = 0.;
	}
	
	//! --- Z-Achse
	if(id_line==2) {
	  //! Start Position
	  start_pos[0] = 0.;
	  start_pos[1] = 0.;
	  start_pos[2] = -Box_Origin[2];
	  //! Direction Vector
	  //! normalised
	  line_vec[0] = 0.;
	  line_vec[1] = 0.;
	  line_vec[2] = 1.;
	  //! Stop Position
	  stop_pos[0] = 0.;
	  stop_pos[1] = 0.;
	  stop_pos[2] = (LX-Box_Origin[2])-1;
	}
	

	
	//! Stepsize
	//! 1/2 * Smallest Box Size 
	delta_r[0] = LX/(2.*pow(2.,(double)MAX_LEVEL)*(BlkNds_X-2)*RB_X);
 	delta_r[1] = LY/(2.*pow(2.,(double)MAX_LEVEL)*(BlkNds_Y-2)*RB_Y);
 	delta_r[2] = LZ/(2.*pow(2.,(double)MAX_LEVEL)*(BlkNds_Z-2)*RB_Z);
// 	delta_r[0] = LX/2.;
// 	delta_r[1] = LY/2.;
// 	delta_r[2] = LZ/2.;
	//! --- End LineOut Config
	
	//! Calculate number of steps
	
	D_REAL temp_pos[3];
	
	num_positions = 1800;
	
	//! alloc coordinates memory
	positions = new D_REAL*[3];
	positions[0] = new D_REAL[num_positions];
	positions[1] = new D_REAL[num_positions];
	positions[2] = new D_REAL[num_positions];
	
	INT32 counter = 0;
	
	while(counter<num_positions)
	{
		//! ----------------------------------------------------------
		//!-----------------------------------------------------------
		//! -> HERE the equations for the positions are defined
		//!	and the positions are written into the array
		//!
		positions[0][counter] = start_pos[0] + line_vec[0]*delta_r[0]*counter;
		positions[1][counter] = start_pos[1] + line_vec[1]*delta_r[1]*counter;
		positions[2][counter] = start_pos[2] + line_vec[2]*delta_r[2]*counter;
		//! ----------------------------------------------------------
		//! ----------------------------------------------------------
	

		//! if generated in km and not in simu units
		//! conversion km -> normalized units
 		positions[0][counter] = positions[0][counter]/lineout_pos_conversion_fact;
 		positions[1][counter] = positions[1][counter]/lineout_pos_conversion_fact;
 		positions[2][counter] = positions[2][counter]/lineout_pos_conversion_fact;

		//! normalize coordinate to box length
		//! get normalized position in Simu-Box
		positions[0][counter] = (positions[0][counter]+Box_Origin[0])/LX;
		positions[1][counter] = (positions[1][counter]+Box_Origin[1])/LY;
		positions[2][counter] = (positions[2][counter]+Box_Origin[2])/LZ;
				
// 		log_file << "Position: " << positions[0][counter] << " | " <<  positions[1][counter] << " | " << positions[2][counter] << " Counter " << counter << endl;
		counter++;
	}
	log_file << "  ->valid positions:       " << num_positions      << endl;

}

//!------------------------------------------------------------
//!- lineout_output:
//!
//!------------------------------------------------------------
void CHybrid::lineout_output(INT32 id_line)
{


	D_REAL** positions;
	INT32 num_values, num_positions;


	//!-----------------------------------
	//! 1) get positions on line
	//!-----------------------------------
	get_positions_from_analytical_formula(id_line, num_positions, positions);
	//!-----------------------------------
	//! 2) estimate field data
	//!-----------------------------------
	char buffer2[200];
	sprintf(buffer2,"mkdir %s/lineout/temp", data_output_path);

	if(!mpi_myRank)
	system(buffer2);
	//! wait until directory is created
	synchronize_allProcesses_MPI();


	if(!num_positions)
	log_file << " no valid positions." << endl;
	else
	for(INT32 field=0; field<num_fields_to_trace_lineout; field++)
	{
		//! calculate current
		if(ID_of_Fields_to_trace_lineout[field]==id_rotB)
		calc_rot(id_BEven, id_rotB);

		//! calculate divergence
		if(ID_of_Fields_to_trace_lineout[field]==id_divB)
		calc_div(id_BEven, id_divB);
		
		//! calculate pressure magnetic
		if(ID_of_Fields_to_trace_lineout[field]==id_PMagnetic)
		square_Field(id_PMagnetic, id_BTotal);
		
		//! calculate ram magnetic
		if(ID_of_Fields_to_trace_lineout[field]==id_scratch_vector) {
		  //! ---- Dynamic Pressure
		  //! Vector quantity
		  //! p_dyn^* = SUM_species ((n*)_species * (m*)_species * (u*)_species)
		  set_zero_field(id_scratch_vector); //! Save Total Value
		  set_zero_field(id_gradPI0);

		  
		  collect_Ui_vth2_of_species(0, id_UI_Species1, 0, noVTH2);
		  Multiply_fields(id_gradPI0,id_UI_Species1,id_UI_Species1);
		  Multiply_fields(id_scratch_vector,id_gradPI0,id_rhoSpecies1);
		    
		}
		
		
		//! calculate divergence
		if(ID_of_Fields_to_trace_lineout[field]==id_divU)
		calc_div(id_UI_plus, id_divU);

		for(INT32 species=0; species<num_Charged_Species; species++)
                    if(ID_of_Fields_to_trace[field]==id_UI_Species1 +species)
                        collect_Ui_vth2_of_species(species, id_UI_Species1+species, 0, noVTH2);
                    
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
                        
		if(ID_of_Fields_to_trace_lineout[field]==ENERGY_SPECTRUM)
		trace_line_energy_spectrum(id_line, ID_of_Fields_to_trace_lineout[field], num_positions, positions);
		else
		trace_line(id_line, ID_of_Fields_to_trace_lineout[field], num_positions, positions);
	      
		if(pack_parallel_output_to_single_file_lineout)
		{

			//! All processors must have written files before master can
			//! begin packing
			synchronize_allProcesses_MPI();
			
			if(!mpi_myRank)
			{
				if(ID_of_Fields_to_trace_lineout[field]==ENERGY_SPECTRUM)
				parallel_output_to_single_lineout_file_energy_spectrum(id_line, ID_of_Fields_to_trace_lineout[field]);
				else
				parallel_output_to_single_file_lineout(id_line, ID_of_Fields_to_trace_lineout[field]);
			}
			
		}

		
		
	}


	//!-----------------------------------
	//! 3) clean up
	//!-----------------------------------

	//! delete coordinates memory
	delete[] positions[0];
	delete[] positions[1];
	delete[] positions[2];

	delete[] positions;

	sprintf(buffer2,"rm -r %s/lineout/temp", data_output_path);

	if(!mpi_myRank)
	system(buffer2);
	
}

//!------------------------------------------------------------
//!parallel_output_to_single_file_lineout(INT32 id_line, INT32 id_field)
//!
//! - Field_id: id of Field that shall be traced along lineout
//! - positions: x,y,z = [0;1[ normalized position in Simu-Box
//!------------------------------------------------------------
void CHybrid::parallel_output_to_single_file_lineout(INT32 id_line, INT32 id_field)
{

	log_file << "     Packing files of field '" << Field_Name[id_field] << "' to single file... ";


	char buffer[200];
	double Field[3], r[3];
	char outfile_name[200];

	ofstream outfile;


	//! decide whether to attach data to file or open new one
	if(create_new_file_each_TL_lineout)
	{
		sprintf(outfile_name,"%s/lineout/%s_%s_TL%06d.txt", data_output_path, LineOut_FileName[id_line], Field_Name[id_field], TL);
		outfile.open(outfile_name);


	}
	else
	{

		sprintf(outfile_name,"%s/lineout/%s_%s.txt", data_output_path, LineOut_FileName[id_line], Field_Name[id_field]);

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
		sprintf(infile_process,"%s/lineout/temp/%s_%s_p%05d_TL%06d.txt", data_output_path, LineOut_FileName[id_line], Field_Name[id_field], process, TL);
		else
		sprintf(infile_process,"%s/lineout/temp/%s_%s_p%05d.txt", data_output_path, LineOut_FileName[id_line], Field_Name[id_field], process);


		ifstream infile;
		infile.open(infile_process);


		if(!infile)
		{
			log_file << "Error in 'parallel_output_to_single_file_lineout':" << endl;
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
			outfile <<  vec_len(r)  <<"		";


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
//!trace_line:
//!
//! - Field_id: id of Field that shall be traced along line
//! - positions: x,y,z = [0;1[ normalized position in Simu-Box
//!------------------------------------------------------------
void CHybrid::trace_line(INT32 id_line, INT32 Field_id, INT32 num_positions, D_REAL** positions)
{


	log_file <<"     Tracing "<<Field_Name[Field_id]<<" ...  " << endl;

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
		  sprintf(filename,"%s/lineout/temp/%s_%s_p%05d_TL%06d.txt", data_output_path, LineOut_FileName[id_line], Field_Name[Field_id], mpi_myRank, TL);
		}
		else
		{
		  sprintf(filename,"%s/lineout/%s_%s_p%05d_TL%06d.txt", data_output_path, LineOut_FileName[id_line], Field_Name[Field_id], mpi_myRank, TL);
		}

		lineout_outfile.open(filename);
	}
	else
	{

		char filename[200];

		if(pack_parallel_output_to_single_file_lineout)
		{
		  sprintf(filename,"%s/lineout/temp/%s_%s_p%05d.txt", data_output_path, LineOut_FileName[id_line], Field_Name[Field_id], mpi_myRank);
		}
		else
		{
		  sprintf(filename,"%s/lineout/%s_%s_p%05d.txt", data_output_path, LineOut_FileName[id_line], Field_Name[Field_id], mpi_myRank);
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
	



				lineout_outfile << positions[0][counter] *LX -Box_Origin[0] <<"		"
						<< positions[1][counter] *LY -Box_Origin[1] <<"		"
						<< positions[2][counter] *LZ -Box_Origin[2] <<"		";

		
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
		
					lineout_outfile << A[comp] << "	";
				}
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

//!-------------------------------------------------------------//
//! silo_2D_LineOut						//
//!-------------------------------------------------------------//
void CHybrid::silo_2D_LineOut(INT32 id_line)
{

	log_file << " Writing line " << LineOut_FileName[id_line] <<endl;

	D_REAL** positions;
	INT32 num_values, num_positions;

	//!-----------------------------------
	//! 1) get positions on line
	//!-----------------------------------
	get_positions_from_analytical_formula(id_line, num_positions, positions);
	
	//! write Orbit for each CrossSection
	silo_2DWrite_Line(id_line, num_positions, positions, 0);
	silo_2DWrite_Line(id_line, num_positions, positions, 1);
	silo_2DWrite_Line(id_line, num_positions, positions, 2);


	//!-----------------------------------
	//! 4) clean up
	//!-----------------------------------
	
	//! delete coordinates memory
	delete[] positions[0];
	delete[] positions[1];
	delete[] positions[2];

	delete[] positions;


}

//!-------------------------------------------------------------//
//! silo_3D_LineOut(INT32 id_line)				//
//!-------------------------------------------------------------//
void CHybrid::silo_3D_LineOut(INT32 id_line)
{
	
	log_file << " Writing line " << LineOut_FileName[id_line] <<endl;

	D_REAL** positions;
	INT32 num_values, num_positions;
	
	
	//!-----------------------------------
	//! 1) get positions on line
	//!-----------------------------------
	get_positions_from_analytical_formula(id_line, num_positions, positions);


	//! write Orbit
	silo_3DWrite_Line(id_line, num_positions, positions);


	//!-----------------------------------
	//! 4) clean up
	//!-----------------------------------
	//! delete time_string_array memory	

	//! delete coordinates memory
	delete[] positions[0];
	delete[] positions[1];
	delete[] positions[2];

	delete[] positions;


}
