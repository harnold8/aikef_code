//! As stated in parameters.cpp, uniform grid interpolates the values of the
//! adaptive mesh to an equidistant cartesian mesh. Of course, you either loose
//! resolution or produce data without information. However, this could be 
//! neccessary if you want to distribute your results to someone, who is not
//! using our code or to use the results of one simulation as an input for 
//! another, while both simulations might have a different mesh (i.e. number
//! of root blocks, level of refinement, whatever).

//! resulting file is binary (otherwise writing and reading of the most
//! propably occuring large files would take very long)

//! uniform_grid_output() starts writing the file with the general information
//! given in parameters.cpp and then writes the positions of the cartesian
//! mesh nodes.

//! uniform_grid_write() calls uniform_grid_writeField() for every process
//! and uniform_grid_parallel_output_to_single_file() for process0

//! in uniform_grid_writeField() each block writes the values of the field in that
//! block in a separate file
//! The way in which the arrays are used, might be confusing, but since the number
//! of grid nodes might be different for each process, at first all values are
//! written in a large array (box_sized_array) and then transfered to a smaller array
//! (blk_sized_array)

//! uniform_grid_parallel_output_to_single_file() reads these files and attaches
//! the field for the whole box to the main output file

// I am pretty sure that some functions could be optimized, but since this works
// and writing output is not done very often, i will not spend more time on this.


#include "parameters.h"
#include "output_uniform_grid.h"
#include "CHybrid.h"
#include "absolute_Globals.h"
#include "utils.h"

#include <iostream>
#include <fstream>



// D_REAL shrink_domain = 1.;

D_REAL LX_uniform_grid = LX/uniform_grid_shrink_domain[0];
D_REAL LY_uniform_grid = LY/uniform_grid_shrink_domain[1];
D_REAL LZ_uniform_grid = LZ/uniform_grid_shrink_domain[2];

D_REAL origin_uniform_grid[3] = {Box_Origin[0]/uniform_grid_shrink_domain[0],
			    Box_Origin[1]/uniform_grid_shrink_domain[1],
			    Box_Origin[2]/uniform_grid_shrink_domain[2]};

extern CHybrid* Hybrid;
ofstream uniform_grid_Outfile;
ofstream uniform_grid_Outfile_field;

const INT32 num_grid_nds =  numNds_uniform_gridMesh[0] * numNds_uniform_gridMesh[1] * numNds_uniform_gridMesh[2];
FILE_REAL* field_array_blk_boxsize;
INT32* index_node_box;
//!-------------------------------------------------------------------------
//! init_BoxInfo_Structure:
//!
//!-------------------------------------------------------------------------
void init_GridInfo_Structure(SGrid_Info* Grid_Info, int TL)
{
 

	//! set SI Values 

	//! B in T
        Grid_Info->B0 = uniform_grid_B0SI;

	//! rho in N/m³
        Grid_Info->n0 = uniform_grid_n0SI;

	//! x in m
        Grid_Info->x0 = uniform_grid_x0SI;

	//! v in m/s
        Grid_Info->v0 = uniform_grid_v0SI;

        sprintf(Grid_Info->Run_Name,Run_Name);

        Grid_Info->TL = TL;

        Grid_Info->NP[0] = numNds_uniform_gridMesh[0];
        Grid_Info->NP[1] = numNds_uniform_gridMesh[1];
        Grid_Info->NP[2] = numNds_uniform_gridMesh[2];

        Grid_Info->Length[0] = LX_uniform_grid;
        Grid_Info->Length[1] = LY_uniform_grid;
        Grid_Info->Length[2] = LZ_uniform_grid;

	//! for visualization, origin should be set to zero
	//! since this is the case for EVERY simulation,
	//! this could be done in uniform_grid source code
	Grid_Info->Origin[0] = 0.;
	Grid_Info->Origin[1] = 0.;
	Grid_Info->Origin[2] = 0.;

	//! for writing the plume, origin must be box_origin
	//! (for correct reading of uniform_grid output as input)
// 	Grid_Info->Origin[0] = origin_uniform_grid[0] -0.5*LX;
// 	Grid_Info->Origin[1] = origin_uniform_grid[1] -0.5*LY;
// 	Grid_Info->Origin[2] = origin_uniform_grid[2] -0.5*LZ;


        Grid_Info->R_Obstacle = R_Obstacle;


}

//!------------------------------------------------------------
//!- uniform_grid_writeField:
//!------------------------------------------------------------
void uniform_grid_writeField(INT32 Field_ID, INT32 uniform_grid_TYPE)
{

	log_file.flush();
	log_file << "     Write field "<< Field_Name[Field_ID] << " ...";
	log_file.flush();


	INT32 u_v_w;

        short num_comps = COMPs_FType[Field_ID];

	//! indices for interpolation
	INT32  i,j,k;
	INT32  i_j_k, ip1_jp1_kp1;
	INT32  ip1_j_k, i_jp1_k,  i_j_kp1;
	INT32  ip1_jp1_k, ip1_j_kp1,  i_jp1_kp1;

	INT32  *index_blk[3], *index_cell[3], blk_index, num_uniform_gridNds_in_Blk=0;
	D_REAL *x_cell[3], r_vec[3], A[3], shape_func[8], *BlkField[3];

	CBlock *temp_Block, *top_level_Block;

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

	for(INT32 u=0; u < numNds_uniform_gridMesh[0]  ; u++)
	  for(INT32 v=0; v < numNds_uniform_gridMesh[1] ; v++)
	    for(INT32 w=0; w < numNds_uniform_gridMesh[2] ; w++)
	    {
		
		u_v_w  = u*numNds_uniform_gridMesh[1]*numNds_uniform_gridMesh[2]
			+v*numNds_uniform_gridMesh[2]
			+w;

			
		//! same as for grid nodes
		//! TODO 0.99 weil sonst außerhalb der Box
		r_vec[0] = - origin_uniform_grid[0]  + u * LX_uniform_grid/(numNds_uniform_gridMesh[0]-0.99) ;
		r_vec[1] = - origin_uniform_grid[1]  + v * LY_uniform_grid/(numNds_uniform_gridMesh[1]-0.99) ;
		r_vec[2] = - origin_uniform_grid[2]  + w * LZ_uniform_grid/(numNds_uniform_gridMesh[2]-0.99) ;


		//! normalize coordinate to box length
		//! get normalized position in Simu-Box
		r_vec[0] = (r_vec[0]+Box_Origin[0])/LX;
		r_vec[1] = (r_vec[1]+Box_Origin[1])/LY;
		r_vec[2] = (r_vec[2]+Box_Origin[2])/LZ;
	
		//!--------------------------------------------------
		//! 2a) estimate blk Cell/Blk-Indices up to MAX_LEVEL
		//!     (even though target blk might not be refined)
		//!--------------------------------------------------
	
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
				for(INT32 comp=0; comp<COMPs_FType[Field_ID]; comp++)
				BlkField[comp] = top_level_Block->Field_Type[Field_ID] +comp *num_nodes_in_block;
				
		
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
	


				//! this could only be done in loop for respective proc
				index_node_box[u_v_w] = u_v_w;
				//! since only nodes for respective proc should be written and all other
				//! nodes in this array will be dropped, the "real" zero has to be marked
				if(u_v_w==0)
				index_node_box[u_v_w] = -1;

				//! interpolate and store field in variable A
				for(INT32 comp=0; comp<COMPs_FType[Field_ID]; comp++)
				{
		
				      A[comp] = BlkField[comp][  i_j_k] * shape_func[0]
						+BlkField[comp][ip1_j_k] * shape_func[1]
						+BlkField[comp][i_jp1_k] * shape_func[2]
						+BlkField[comp][i_j_kp1] * shape_func[3]
				
						+BlkField[comp][  ip1_jp1_k] * shape_func[4]
						+BlkField[comp][  ip1_j_kp1] * shape_func[5]
						+BlkField[comp][  i_jp1_kp1] * shape_func[6]
						+BlkField[comp][ip1_jp1_kp1] * shape_func[7];
		

					field_array_blk_boxsize[comp*num_grid_nds+num_uniform_gridNds_in_Blk] = A[comp];
				}
				num_uniform_gridNds_in_Blk++;
			
			}//! end if responsible_mpi_process

		}

	     }//! end for positions

	INT32* index_node_blk = new INT32[num_uniform_gridNds_in_Blk];
	
	//! drop all entries which are not in this blk/proc
	INT32 z=0;
	for(INT32 j=0; j<num_grid_nds; j++)
	{
		if(index_node_box[j]!=0)
		{
			index_node_blk[z]=index_node_box[j];
			if(index_node_box[j]==-1)
			index_node_blk[z]=0;

			z++;
		}
	}

	FILE_REAL* field_array_blk = new FILE_REAL[num_comps*num_uniform_gridNds_in_Blk];
	for(INT32 k=0;k<num_uniform_gridNds_in_Blk;k++)
	{
		for(INT32 comp=0; comp<COMPs_FType[Field_ID]; comp++)
		{
			field_array_blk[comp*num_uniform_gridNds_in_Blk + k] = field_array_blk_boxsize[comp*num_grid_nds+k];
		}
	}	


        char filename[100];

        sprintf(filename,"%s/uniform_grid/temp/%s_uniform_grid_%s_TL%d_mpi%d.dat",data_output_path,Run_Name,Field_Name[Field_ID],TL,mpi_myRank);
        uniform_grid_Outfile_field.open(filename, ofstream::binary);
	
        uniform_grid_Outfile_field.write( reinterpret_cast<char*> (&num_comps),sizeof(short));
	//! write field_array_blk
        uniform_grid_Outfile_field.write( reinterpret_cast<char*> (&num_uniform_gridNds_in_Blk), sizeof(INT32));
        uniform_grid_Outfile_field.write( reinterpret_cast<char*> (index_node_blk), num_uniform_gridNds_in_Blk*sizeof(INT32));
        uniform_grid_Outfile_field.write( reinterpret_cast<char*> (field_array_blk), num_comps *num_uniform_gridNds_in_Blk*sizeof(FILE_REAL));


	uniform_grid_Outfile_field.flush();
        uniform_grid_Outfile_field.close();

	//! clean up memory

	delete[] index_node_blk;
	delete[] field_array_blk;

	log_file << " done." <<endl;

}

//!------------------------------------------------------------
//!- uniform_grid_output:
//!------------------------------------------------------------
void uniform_grid_output(INT32 TL)
{

	log_file << endl << " Writing uniform_grid 3D Data of TL " << TL <<" ..." << endl;


        SGrid_Info* Grid_Info = new SGrid_Info;
        memset(Grid_Info,0,sizeof(SGrid_Info));

        //! Write Box info & Si constants to Info Structure
        init_GridInfo_Structure(Grid_Info, TL);

        char filename[100];

	Hybrid->synchronize_allProcesses_MPI();

	if(!mpi_myRank)
	{
		sprintf(filename,"%s/uniform_grid/%s_uniform_grid_TL%d.dat",data_output_path,Run_Name,TL);
	
		uniform_grid_Outfile.open(filename, ofstream::binary);
	
		//!  write down Info structure
		uniform_grid_Outfile.write( reinterpret_cast<char*> (Grid_Info),sizeof(SGrid_Info));
		delete(Grid_Info);
		
		//!--------- Write Nodes
		FILE_REAL* mesh_nds_pos = new FILE_REAL[3*num_grid_nds];
		INT32 u_v_w;

		//! calculate coordinates of mesh nodes
		for(INT32 u=0; u < numNds_uniform_gridMesh[0] ; u++)
		 for(INT32 v=0; v < numNds_uniform_gridMesh[1] ; v++)
		  for(INT32 w=0; w < numNds_uniform_gridMesh[2] ; w++)
		  {
	
			u_v_w  = u*numNds_uniform_gridMesh[1]*numNds_uniform_gridMesh[2]
				+v*numNds_uniform_gridMesh[2]
				+w;	
	
			mesh_nds_pos[0*num_grid_nds +u_v_w] = - origin_uniform_grid[0] + 0.5* LX_uniform_grid/numNds_uniform_gridMesh[0] + u * LX_uniform_grid/(numNds_uniform_gridMesh[0]-1);
			mesh_nds_pos[1*num_grid_nds +u_v_w] = - origin_uniform_grid[1] + 0.5* LY_uniform_grid/numNds_uniform_gridMesh[1] + v * LY_uniform_grid/(numNds_uniform_gridMesh[1]-1);
			mesh_nds_pos[2*num_grid_nds +u_v_w] = - origin_uniform_grid[2] + 0.5* LZ_uniform_grid/numNds_uniform_gridMesh[2] + w * LZ_uniform_grid/(numNds_uniform_gridMesh[2]-1);
	
		  }

		//! First Field has to be Grid Node
		uniform_grid_Outfile.write( reinterpret_cast<char*> (mesh_nds_pos),3*num_grid_nds*sizeof(FILE_REAL));
	
		delete[] mesh_nds_pos;
	
		char buffer2[200];
		sprintf(buffer2,"mkdir %s/uniform_grid/temp", data_output_path);
		system(buffer2);

	}//! end loop for process 0

	//! wait until directory is created
	Hybrid->synchronize_allProcesses_MPI();

	

	index_node_box = new INT32[num_grid_nds];
	memset(index_node_box,0,num_grid_nds*sizeof(INT32));

	//! write Fields
	for(INT32 i=0; i<uniform_grid_num_fields_to_trace; i++)
	{
		//! allocate field_array_blk to store field of type FILE_REAL
		field_array_blk_boxsize = new FILE_REAL[COMPs_FType[uniform_grid_ID_of_Fields_to_trace[i]]*num_grid_nds];

		uniform_grid_write(uniform_grid_ID_of_Fields_to_trace[i], uniform_grid_Field_type[i]);

		delete[] field_array_blk_boxsize;

	}

	delete[] index_node_box;

	//! ---- rhoSpecies/UI_Species Cuts -------------------
// 	for(int species=0; species<num_Charged_Species; species++)
// 	{

		//! write rho of species
// 		uniform_grid_write(id_rhoSpecies1  +species, uniform_grid_TYPE_DENSITY_FIELD);

		//! provide rho
		//! gather Ui 
// 		Hybrid->collect_Ui_vth2_of_species(species, id_UI_Species1, 0,noVTH2);
// 		uniform_grid_write(id_UI_Species1 +species, uniform_grid_TYPE_VELOCITY_FIELD);

// 	}

	if(!mpi_myRank)
	{
		uniform_grid_Outfile.flush();
		uniform_grid_Outfile.close();
		sprintf(filename,"rm -r %s/uniform_grid/temp", data_output_path);
		system(filename);
	}

	log_file  << " done." << endl << endl;

}


void uniform_grid_write(INT32 Field_ID, INT32 uniform_grid_TYPE)
{
	//! following has to be done for each field

        //! write Field 
	uniform_grid_writeField(Field_ID, uniform_grid_TYPE);
	//! All processors must have written files before master can
	//! begin packing
	Hybrid->synchronize_allProcesses_MPI();
	if(!mpi_myRank)
	uniform_grid_parallel_output_to_single_file(Field_ID, uniform_grid_TYPE, num_grid_nds);

// 	Hybrid->synchronize_allProcesses_MPI();
}



void uniform_grid_parallel_output_to_single_file(INT32 Field_ID, INT32 uniform_grid_TYPE, INT32 num_grid_nds)
{
	log_file << "     Packing files of field " << Field_Name[Field_ID] << " to single file ...";

	char  buffer[200];
	INT32 uniform_grid_Field_Name_Size = 50;
	int ftype = uniform_grid_TYPE;
        short num_comps = COMPs_FType[Field_ID];
	short buff_num_comps;
	INT32 num_uniform_gridNds_in_Blk;
	FILE_REAL* field_array_box = new FILE_REAL[num_comps*num_grid_nds];

	for(int process=0; process<mpi_num_processes; process++)
	{

		char infile_process[200];
       		sprintf(infile_process,"%s/uniform_grid/temp/%s_uniform_grid_%s_TL%d_mpi%d.dat",data_output_path,Run_Name,Field_Name[Field_ID],TL,process);

		ifstream infile;
		infile.open(infile_process, ifstream::binary);

		if(!infile)
		{
			log_file << "Error in 'parallel_output_to_single_file':" << endl;
			log_file << "no file: " << infile_process << endl;
			log_file << "Returning ... ." << endl;
			return;			
		}
	
		infile.read( reinterpret_cast<char*> (&buff_num_comps),sizeof(short));

		//! check if reading file works properly
		if(buff_num_comps != num_comps)
		{
			log_file << "Error in 'parallel_output_to_single_file':" << endl;
			log_file << "Mismatch in number of field components: " << infile_process << endl;
			log_file << "Returning ... ." << endl;
			return;
		}

		//! first read number of nodes written by each process
		infile.read(reinterpret_cast<char*> (&num_uniform_gridNds_in_Blk), sizeof(INT32));

		//! alloc arrays for reading
		INT32* ijk_array_proc = new INT32[num_uniform_gridNds_in_Blk];
		FILE_REAL* field_array_proc = new FILE_REAL[buff_num_comps*num_uniform_gridNds_in_Blk];

		//! read indices of nodes written by this process
		infile.read(reinterpret_cast<char*> (ijk_array_proc),  num_uniform_gridNds_in_Blk*sizeof(INT32));

		//! read field
		infile.read(reinterpret_cast<char*> (field_array_proc),  buff_num_comps*num_uniform_gridNds_in_Blk*sizeof(FILE_REAL));

		//! write field in one array for all processes
		for(INT32 i=0; i<num_uniform_gridNds_in_Blk; i++)
		  for(INT32 comp=0; comp<COMPs_FType[Field_ID]; comp++)
		{		
			field_array_box[comp*num_grid_nds +ijk_array_proc[i]] = field_array_proc[comp*num_uniform_gridNds_in_Blk + i];
		}

		infile.close();

		delete[] ijk_array_proc;
		delete[] field_array_proc;
		
		//! remove temporary files
		sprintf(buffer,"rm %s", infile_process);
		system(buffer);

	}

	uniform_grid_Outfile.write( reinterpret_cast<char*> (&num_comps),sizeof(short));
	uniform_grid_Outfile.write( reinterpret_cast<char*> (&ftype),sizeof(int));
        uniform_grid_Outfile.write( reinterpret_cast<char*> (Field_Name[Field_ID]), uniform_grid_Field_Name_Size*sizeof(char));
        uniform_grid_Outfile.write( reinterpret_cast<char*> (field_array_box), num_comps*num_grid_nds*sizeof(FILE_REAL));

	delete[] field_array_box;

	log_file << " done. "  << endl;


}
