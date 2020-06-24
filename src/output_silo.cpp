

#include "CBlk.h"
#include "utils.h"
#include "parameters.h"
#include "output_silo.h"
#include "absolute_Globals.h"


#include "silo.h"
#include <iostream>
#include <fstream>


extern D_REAL **Blk_Length_of;



//! file pointer to silo file
DBfile *silo_file2D;
DBfile *silo_file3D;

//! number of blocks to write
INT32 num_TopLevelBlks;

INT32 num_TopLevelBlks_CS[3];
INT32 num_MaxTopLevelBlks2D;



//! Number of dimensions and number
//! of nodes each dimension
INT32 ind[3];
INT32 num_dims;



INT32 max_name_length = 100;


INT32 num_FieldValues;
INT32 num_FieldValues_in_dim[3];

FILE_REAL *FieldX, *FieldY, *FieldZ;
FILE_REAL *FieldData[3];

INT32 num_FieldValuesCS[3];
INT32 num_FieldValuesCS_in_dim[3][2];


//! pGN should be visible, mGN should be hidden
//! set 0 to show
//! set 1 to hide

//! (higher resolution in plus direction)
//! (see sketches in document for details)
//! (default: hide minus; show plus)
INT32 hide_pGN = 0;
INT32 hide_mGN = 1;
INT32 CENTER_STYLE;


//! MultiMesh/MultiVar arrays
int   *MultiObject_FieldTypes;

char **MultiObject_FieldNames;
char **MultiObject_FieldNames_XComp;
char **MultiObject_FieldNames_YComp;
char **MultiObject_FieldNames_ZComp;

char silo_filename[200];

//!-------------------------------------------------------------//
//! silo_2DWrite_allCS: -								//
//!-------------------------------------------------------------//
void silo_2DWrite_allCS(INT32 Field_ID)
{

	log_file << "  Write " << Field_Name[Field_ID] <<"... ";


	silo_2DWrite_FieldCS(Field_ID,0);
	silo_2DWrite_FieldCS(Field_ID,1);
	silo_2DWrite_FieldCS(Field_ID,2);

	log_file << "done." << endl;


}
//!-------------------------------------------------------------//
//! silo_2DWrite_prepare: -								//
//!-------------------------------------------------------------//
void silo_2DWrite_prepare(void)
{

	


	//!--------------------------------------------------
	//! 1) - Open Silo file
	//! 	 - Create Blocks directory
	//!--------------------------------------------------

	log_file << " Writing Silo 2DFile ..." << endl;
	silo_file2D = NULL;

	sprintf(silo_filename,"%s/silo/%s_2d_p%05d_TL%05d.silo",data_output_path, Run_Name, mpi_myRank, TL);
	//! Open the Silo file
	//! In case visit is compiled with HDF5 support, DB_PDB can be 
	//! replaced by DB_HDF5
	
	if (COMPRESS_SILO)
	{
	  DBSetCompression("METHOD=GZIP LEVEL=8");
	}
	silo_file2D = DBCreate(silo_filename, DB_CLOBBER, DB_LOCAL, "2d", DB_HDF5);
	if(silo_file2D == NULL)
	{
		log_file <<  "Could not create Silo file!" << endl;
		return;
	}
	

	
	//! create new directory in silo file
	DBMkDir(silo_file2D,"/XCS_Blocks");
	DBMkDir(silo_file2D,"/YCS_Blocks");
	DBMkDir(silo_file2D,"/ZCS_Blocks");

	DIR_of_trajectories_exist = false;
	DIR_of_lines_exist = false;
	DIR_of_particleTracks2D_exist = false;
	
	memset(DIR_of_FieldID_exists,0,NUM_FIELDS*sizeof(bool));

	//!--------------------------------------------------
	//! 2) Specify dimensions to write depending on
	//!    whether or not GN in plus and minus
	//!    direction (pGN, mGN) are to written out
	//!--------------------------------------------------


	//! set style for DB FileOut
	if(SILO_STYLE_ZONAL)
	CENTER_STYLE = DB_ZONECENT;
	else
	CENTER_STYLE= DB_NODECENT;

	//! Set number of nodes each Block
	num_dims = 2;

	num_FieldValuesCS_in_dim[0][0] = (BlkNds_Y -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);
	num_FieldValuesCS_in_dim[0][1] = (BlkNds_Z -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);


	num_FieldValuesCS_in_dim[1][0] = (BlkNds_X -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);
	num_FieldValuesCS_in_dim[1][1] = (BlkNds_Z -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);


	num_FieldValuesCS_in_dim[2][0] = (BlkNds_X -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);
	num_FieldValuesCS_in_dim[2][1] = (BlkNds_Y -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);



	//! Total number of Nodes to write each Block
	num_FieldValuesCS[0] = num_FieldValuesCS_in_dim[0][0] *num_FieldValuesCS_in_dim[0][1];
	num_FieldValuesCS[1] = num_FieldValuesCS_in_dim[1][0] *num_FieldValuesCS_in_dim[1][1];
	num_FieldValuesCS[2] = num_FieldValuesCS_in_dim[2][0] *num_FieldValuesCS_in_dim[2][1];

	//! estimate CS of maximum mesh nodes
	INT32 num_MaxFieldValues2D = num_FieldValuesCS[0];

	if(num_MaxFieldValues2D<num_FieldValuesCS[1])
	num_MaxFieldValues2D = num_FieldValuesCS[1];

	if(num_MaxFieldValues2D<num_FieldValuesCS[2])
	num_MaxFieldValues2D = num_FieldValuesCS[2];

	//!-------------------------------------------------------
	//! - Allocate one array for each component
	//! - Use CS of maximal mesh nodes
	//!-------------------------------------------------------

	FieldX = new FILE_REAL[num_MaxFieldValues2D];
	FieldY = new FILE_REAL[num_MaxFieldValues2D];
	FieldZ = new FILE_REAL[num_MaxFieldValues2D];

	FieldData[0] = FieldX;
	FieldData[1] = FieldY;
	FieldData[2] = FieldZ;





	//!----------------------------------------//
	//!3) Count top level Blks in each CS	  //
	//!----------------------------------------//
	num_TopLevelBlks_CS[0] = 0;
	num_TopLevelBlks_CS[1] = 0;
	num_TopLevelBlks_CS[2] = 0;
	

	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			for(INT32 comp=0; comp<3; comp++)
			if(    CROSS_SECTION[comp] >= temp_Block->origin[comp]
			    && CROSS_SECTION[comp] <  temp_Block->origin[comp] +Blk_Length_of[level][comp])
			num_TopLevelBlks_CS[comp]++;
	
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}
	
	//! estimate CS of maximum Blocks
	num_MaxTopLevelBlks2D = num_TopLevelBlks_CS[0];

	if(num_MaxTopLevelBlks2D<num_TopLevelBlks_CS[1])
	num_MaxTopLevelBlks2D = num_TopLevelBlks_CS[1];

	if(num_MaxTopLevelBlks2D<num_TopLevelBlks_CS[2])
	num_MaxTopLevelBlks2D = num_TopLevelBlks_CS[2];

	//!--------------------------------------------------
	//! 4) free store memory below is used in write Mesh
	//!    and write Field method and will be deallocated
	//!    in function "silo_3DFlush_CleanUp".
	//!--------------------------------------------------
	MultiObject_FieldTypes = new int[num_MaxTopLevelBlks2D];
        
	MultiObject_FieldNames = new char*[num_MaxTopLevelBlks2D];
	MultiObject_FieldNames_XComp = new char*[num_MaxTopLevelBlks2D];
	MultiObject_FieldNames_YComp = new char*[num_MaxTopLevelBlks2D];
	MultiObject_FieldNames_ZComp = new char*[num_MaxTopLevelBlks2D];

	for(INT32 blk=0; blk<num_MaxTopLevelBlks2D; blk++)
	{
		MultiObject_FieldNames[blk] = new char[max_name_length];
		MultiObject_FieldNames_XComp[blk] = new char[max_name_length];
		MultiObject_FieldNames_YComp[blk] = new char[max_name_length];
		MultiObject_FieldNames_ZComp[blk] = new char[max_name_length];
	}
	
}


//!-------------------------------------------------------------//
//! silo_2D_writeTrajectory: -						//
//!-------------------------------------------------------------//
void silo_2DWrite_Trajectory(INT32 Trajectory_ID, INT32 num_positions, D_REAL** positions, INT32 CS)
{


	INT32 lp_ind0, lp_ind1;


	if(!DIR_of_trajectories_exist)
	{
		DBMkDir(silo_file2D,"/Trajectories");
		DBMkDir(silo_file2D,"/Trajectories/X_CrossSection");
		DBMkDir(silo_file2D,"/Trajectories/Y_CrossSection");
		DBMkDir(silo_file2D,"/Trajectories/Z_CrossSection");
		DIR_of_trajectories_exist = true;
	}

	char trajectory_name[200];



	//! Depending of CS, decide which loop indices have to be used
	switch(CS)
	{
		//! X-Cross Section: loop across Y and Z dimension
		case 0: lp_ind0 = 1;
			 lp_ind1 = 2;
			 sprintf(trajectory_name,"/Trajectories/X_CrossSection/%s", Trajectory_FileName[Trajectory_ID]);
			 break;

		//! Y-Cross Section: loop across X and Z dimension
		case 1: lp_ind0 = 0;
			 lp_ind1 = 2;
			 sprintf(trajectory_name,"/Trajectories/Y_CrossSection/%s", Trajectory_FileName[Trajectory_ID]);
			 break;

		//! Z-Cross Section: loop across X and Y dimension
		case 2: lp_ind0 = 0;
			 lp_ind1 = 1;
			 sprintf(trajectory_name,"/Trajectories/Z_CrossSection/%s", Trajectory_FileName[Trajectory_ID]);
			 break;

	}




	double *x = new double[num_positions];
	double *y = new double[num_positions];


	double *coords[] = {x,y};

	 double L[3] = {LX, LY, LZ};

	for(INT32 pos=0; pos<num_positions; pos++)
	{
		x[pos] = (positions[lp_ind0][pos] *L[lp_ind0] -Box_Origin[lp_ind0]) *factor_scale_mesh;
		y[pos] = (positions[lp_ind1][pos] *L[lp_ind1] -Box_Origin[lp_ind1]) *factor_scale_mesh;

	}


	INT32 ndims = 2;
	DBPutPointmesh(silo_file2D, trajectory_name, ndims, coords, num_positions, DB_DOUBLE, NULL);


	delete[] x;
	delete[] y;

}


//!-------------------------------------------------------------//
//!silo_2DWrite_ParticleTrack: -						//
//!-------------------------------------------------------------//
void silo_2DWrite_ParticleTrack(INT32 GROUP_ID, INT32 num_positions, double** positions, INT32 CS)
{


	INT32 lp_ind0, lp_ind1;



	char particle_track_name[200];

	DBoptlist *optlist = DBMakeOptlist(4);
	DBAddOption(optlist, DBOPT_XUNITS, (void *) visit_length_unit);
	DBAddOption(optlist, DBOPT_YUNITS, (void *) visit_length_unit);



	//! Depending of CS, decide which loop indices have to be used
	switch(CS)
	{
		//! X-Cross Section: loop across Y and Z dimension
		case 0: lp_ind0 = 1;
			 lp_ind1 = 2;
			 DBAddOption(optlist, DBOPT_XLABEL, (void *) visit_ylabel);
			 DBAddOption(optlist, DBOPT_YLABEL, (void *) visit_zlabel);
			 sprintf(particle_track_name, "/ParticleTracks/X_CrossSection/Group_%04d", GROUP_ID);
			 break;

		//! Y-Cross Section: loop across X and Z dimension
		case 1: lp_ind0 = 0;
			 lp_ind1 = 2;
			 DBAddOption(optlist, DBOPT_XLABEL, (void *) visit_xlabel);
			 DBAddOption(optlist, DBOPT_YLABEL, (void *) visit_zlabel);
			 sprintf(particle_track_name, "/ParticleTracks/Y_CrossSection/Group_%04d", GROUP_ID);
			 break;

		//! Z-Cross Section: loop across X and Y dimension
		case 2: lp_ind0 = 0;
			 lp_ind1 = 1;
			 DBAddOption(optlist, DBOPT_XLABEL, (void *) visit_xlabel);
			 DBAddOption(optlist, DBOPT_YLABEL, (void *) visit_ylabel);
			 sprintf(particle_track_name, "/ParticleTracks/Z_CrossSection/Group_%04d", GROUP_ID);
			 break;

	}


	double* positions_2D[2] = {positions[lp_ind0], positions[lp_ind1]};


	INT32 ndims = 2;
	DBPutPointmesh(silo_file2D, particle_track_name, ndims, positions_2D, num_positions, DB_DOUBLE, optlist);



}

//!-------------------------------------------------------------//
//! silo_2D_writeLine: -						//
//!-------------------------------------------------------------//
void silo_2DWrite_Line(INT32 Line_ID, INT32 num_positions, D_REAL** positions, INT32 CS)
{


	INT32 lp_ind0, lp_ind1;
	
	if(!DIR_of_lines_exist)
	{	
		DBMkDir(silo_file2D,"/LineOut");
		DBMkDir(silo_file2D,"/LineOut/X_CrossSection");
		DBMkDir(silo_file2D,"/LineOut/Y_CrossSection");
		DBMkDir(silo_file2D,"/LineOut/Z_CrossSection");
		DIR_of_lines_exist = true;
	}
	

	char line_name[200];



	//! Depending of CS, decide which loop indices have to be used
	switch(CS)
	{
		//! X-Cross Section: loop across Y and Z dimension
		case 0: lp_ind0 = 1;
			 lp_ind1 = 2;
			 sprintf(line_name,"/LineOut/X_CrossSection/%s", LineOut_FileName[Line_ID]);
			 break;

		//! Y-Cross Section: loop across X and Z dimension
		case 1: lp_ind0 = 0;
			 lp_ind1 = 2;
			 sprintf(line_name,"/LineOut/Y_CrossSection/%s", LineOut_FileName[Line_ID]);
			 break;

		//! Z-Cross Section: loop across X and Y dimension
		case 2: lp_ind0 = 0;
			 lp_ind1 = 1;
			 sprintf(line_name,"/LineOut/Z_CrossSection/%s", LineOut_FileName[Line_ID]);
			 break;

	}




	double *x = new double[num_positions];
	double *y = new double[num_positions];


	double *coords[] = {x,y};

	double L[3] = {LX, LY, LZ};

	for(INT32 pos=0; pos<num_positions; pos++)
	{
		x[pos] = (positions[lp_ind0][pos] *L[lp_ind0] -Box_Origin[lp_ind0]) *factor_scale_mesh;
		y[pos] = (positions[lp_ind1][pos] *L[lp_ind1] -Box_Origin[lp_ind1]) *factor_scale_mesh;

	}


	INT32 ndims = 2;
	DBPutPointmesh(silo_file2D, line_name, ndims, coords, num_positions, DB_DOUBLE, NULL);


	delete[] x;
	delete[] y;

}

//!-------------------------------------------------------------//
//! silo_2DFlush_CleanUp: -					//
//!-------------------------------------------------------------//
void silo_2DFlush_CleanUp(void)
{

	//! Close the Silo file.
	DBClose(silo_file2D);

	//! free Field Data
	delete[] FieldX;
	delete[] FieldY;
	delete[] FieldZ;

	for(INT32 blk=0; blk<num_MaxTopLevelBlks2D; blk++)
	{
		delete[] MultiObject_FieldNames[blk];
		delete[] MultiObject_FieldNames_XComp[blk];
		delete[] MultiObject_FieldNames_YComp[blk];
		delete[] MultiObject_FieldNames_ZComp[blk];
	}

	delete[] MultiObject_FieldTypes;

	delete[] MultiObject_FieldNames;
	delete[] MultiObject_FieldNames_XComp;
	delete[] MultiObject_FieldNames_YComp;
	delete[] MultiObject_FieldNames_ZComp;

	log_file << " finished." << endl << endl;

}

//!------------------------------------------------------------
//!- write_XCrossSection: Writes out the Field values for
//!	a given Cross Section Coordinate specified in defines.h.
//!   - Whether a Block is within the CS is decided by its 
//!	  Interval [orig,orig + length]
//!   - No Interpolation is performed. Data is written out in 
//!     a nearest neighbour style
//!------------------------------------------------------------
void silo_2DWrite_FieldCS(INT32 Field_ID, INT32 CS)
{

if (COMPRESS_SILO)
	{
	  DBSetCompression("METHOD=GZIP LEVEL=8");
	}

	//!------------------------------------
	//! 1) Variable Declaration
	//!------------------------------------
	//! loop variables
	INT32 i_j_k, i_j_k_DATA, lp_ind0=-1, lp_ind1=-1, lp_ind2=-1, active_TopLevelBlk;
	char active_Meshname[50], temp_filename[50];
	char active_VariableName[50];
	D_REAL *VFieldX, *VFieldY, *VFieldZ;

	//! specify names for variable components
	char comp0[50], comp1[50], comp2[50];
	char* comp_names[3] = {comp0, comp1, comp2};
        
	sprintf(comp_names[0],"%sX",Field_Name[Field_ID]);
	sprintf(comp_names[1],"%sY",Field_Name[Field_ID]);
	sprintf(comp_names[2],"%sZ",Field_Name[Field_ID]);


	//! path depending on which CS o write
	char cs_path[50];

	//! Depending of CS, decide which loop indices have to be used
	switch(CS)
	{
		//! X-Cross Section: loop across Y and Z dimension
		case 0: lp_ind0 = 1;
			 lp_ind1 = 2;
			 lp_ind2 = 0;
			 sprintf(cs_path,"/XCS_Blocks");
			 break;

		//! Y-Cross Section: loop across X and Z dimension
		case 1: lp_ind0 = 0;
			 lp_ind1 = 2;
			 lp_ind2 = 1;
			 sprintf(cs_path,"/YCS_Blocks");
			 break;

		//! Z-Cross Section: loop across X and Y dimension
		case 2: lp_ind0 = 0;
			 lp_ind1 = 1;
			 lp_ind2 = 2;
			 sprintf(cs_path,"/ZCS_Blocks");
			 break;

	}

	//! calc blk_indedx_L0 and cell_indedx_L
	calc_cells_to_cut(CS, cell_indedx_L);

	//!------------------------------------------------------------------//
	//! 	3a) Master create MultiMesh names				     //
	//!	    This loop must only be called by master process !!!	     //
	//!	    since every process is writing to its own file,		     //
	//!	    the name of the file has to be specified within		     //
	//!	    the MULTIMESH OBJECT prior to the respective meshname !!!  //
	//!	    TODO:								     //
	//!	    Use limited numbre of shared files using MPIO !!!	     //
	//!------------------------------------------------------------------//
	active_TopLevelBlk=0;
	if(!mpi_myRank)
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

	
			if(   CROSS_SECTION[CS] >= temp_Block->origin[CS]
				&& CROSS_SECTION[CS] <  temp_Block->origin[CS]+Blk_Length_of[level][CS])
			{


				//! Don't specify filename in case master has to process Block
				//! since multivar and var are in the same file
				//! FOR SOME STRANGE REASON THIS HIDES ALL THE OTHER BLOCK DIRECTORIES
				//! THAT ARE VISIBLE IN CASE FILENAME IS SPECIFIED FOR MASTER ALSO !!!
				//! NOTE:
				//! IN CONTRAST TO 3D OUTPUT, THE SLASH HAS TO BE SUPPRESSED ALSO !!!!
				if(!temp_Block->responsible_mpi_process)
				sprintf(temp_filename,"");
				else
		   		sprintf(temp_filename,"%s_2d_p%05d_TL%05d.silo:/",
		   				      Run_Name,
		   				      temp_Block->responsible_mpi_process,
		   				      TL);


				MultiObject_FieldTypes[active_TopLevelBlk] = DB_QUADVAR;

				//!-----------------------------------------//
				//!  write scalar field			     //
				//!----------------------------------------//
				if(COMPs_FType[Field_ID]==1)
				{
	
					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames[active_TopLevelBlk],"%s%s/Block%d/%s%d",
											temp_filename,
											cs_path,
											active_TopLevelBlk,
											Field_Name[Field_ID],
											active_TopLevelBlk);
		
	
				}



				//!-----------------------------------------//
				//!  write vector field			     //
				//!----------------------------------------//
				if(COMPs_FType[Field_ID]==3)
				{
	

					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames_XComp[active_TopLevelBlk],"%s%s/Block%d/%s%d_X",
									     temp_filename,
									     cs_path,
									     active_TopLevelBlk,
									     Field_Name[Field_ID],
									     active_TopLevelBlk);
	
					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames_YComp[active_TopLevelBlk],"%s%s/Block%d/%s%d_Y",
									     temp_filename,
									     cs_path,
									     active_TopLevelBlk,
									     Field_Name[Field_ID],
									     active_TopLevelBlk);
	
					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames_ZComp[active_TopLevelBlk],"%s%s/Block%d/%s%d_Z",
									     temp_filename,
									     cs_path,
									     active_TopLevelBlk,
									     Field_Name[Field_ID],
									     active_TopLevelBlk);


	
				}

				//! increase active_TopLevelBlk independent of COMPs_FType
				active_TopLevelBlk++;

			}
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}//! end while
	}//! end for level


	active_TopLevelBlk = 0;
	//!-----------------------------------------//
	//! 	Loop across Blks		    //
	//!-----------------------------------------//
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


		  if(   CROSS_SECTION[CS] >= temp_Block->origin[CS]
			&& CROSS_SECTION[CS] <  temp_Block->origin[CS]+Blk_Length_of[level][CS])
		  {

		      if(mpi_myRank==temp_Block->responsible_mpi_process)
		      {


			//! Names/Types of Block Variable
			sprintf(active_Meshname, "%s/Block%d/Mesh%d",
						 cs_path,
						 active_TopLevelBlk,
						 active_TopLevelBlk);

			sprintf(active_VariableName, "%s/Block%d/%s%d",
						     cs_path,
						     active_TopLevelBlk,
						     Field_Name[Field_ID],
						     active_TopLevelBlk);



			//! set constant cell index to cell CS of respective level
			ind[lp_ind2] = cell_indedx_L[level];

			//!-----------------------------------------//
			//!  write scalar field				   //
			//!----------------------------------------//
			if(COMPs_FType[Field_ID]==1)
			{

				
				//! Loop across CS:
				//! Finish at num_FieldValuesCS_in_dim[CS][0] +hide_mGN, since
				//! hide_mGN is substracted prior to this loop.
				for(ind[lp_ind0]=hide_mGN; ind[lp_ind0] <num_FieldValuesCS_in_dim[CS][0] +hide_mGN; ind[lp_ind0]++)
				 for(ind[lp_ind1]=hide_mGN; ind[lp_ind1] <num_FieldValuesCS_in_dim[CS][1] +hide_mGN; ind[lp_ind1]++)
				 {
				
	
					i_j_k  = ind[0]*BlkNds_Y*BlkNds_Z 
						+ind[1]*BlkNds_Z
						+ind[2];
	
					i_j_k_DATA  = (ind[lp_ind1]-hide_mGN)*num_FieldValuesCS_in_dim[CS][0]
						     +(ind[lp_ind0]-hide_mGN);
	
					FieldX[i_j_k_DATA] = temp_Block->Field_Type[Field_ID][i_j_k];
					

				 }


				DBPutQuadvar1(silo_file2D,
					      active_VariableName,
					      active_Meshname,
					      FieldX,num_FieldValuesCS_in_dim[CS],
					      num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);




		   	 }//! end if scalar field


			//!-----------------------------------------//
			//!  write vector field			     //
			//!----------------------------------------//
			if(COMPs_FType[Field_ID]==3)
			{

				//! Write Vector
				VFieldX = temp_Block->Field_Type[Field_ID];
				VFieldY = VFieldX +num_nodes_in_block;
				VFieldZ = VFieldY +num_nodes_in_block;

				//! Loop across CS:
				//! Finish at num_FieldValuesCS_in_dim[CS][0] +hide_mGN, since
				//! hide_mGN is substracted prior to this loop.
				for(ind[lp_ind0]=hide_mGN; ind[lp_ind0] <num_FieldValuesCS_in_dim[CS][0] +hide_mGN; ind[lp_ind0]++)
				 for(ind[lp_ind1]=hide_mGN; ind[lp_ind1] <num_FieldValuesCS_in_dim[CS][1] +hide_mGN; ind[lp_ind1]++)
				 {

					i_j_k  = ind[0]*BlkNds_Y*BlkNds_Z 
						+ind[1]*BlkNds_Z
						+ind[2];
	
					i_j_k_DATA  = (ind[lp_ind1]-hide_mGN)*num_FieldValuesCS_in_dim[CS][0]
						     +(ind[lp_ind0]-hide_mGN);

				
					FieldX[i_j_k_DATA] = VFieldX[i_j_k];
					FieldY[i_j_k_DATA] = VFieldY[i_j_k];
					FieldZ[i_j_k_DATA] = VFieldZ[i_j_k];
					
				 }

			
				//! ---------- write X Component of Block ------------------------------------
				sprintf(active_VariableName,"%s/Block%d/%s%d_X",
							    cs_path,
							    active_TopLevelBlk,
							    Field_Name[Field_ID],
							    active_TopLevelBlk);

				DBPutQuadvar1(silo_file2D,
					active_VariableName,
					active_Meshname,
					FieldX,num_FieldValuesCS_in_dim[CS],
					num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);


				//! ---------- write Y Component of Block ------------------------------------
				sprintf(active_VariableName,"%s/Block%d/%s%d_Y",
							    cs_path,
							    active_TopLevelBlk,
							    Field_Name[Field_ID],
							    active_TopLevelBlk);

				DBPutQuadvar1(silo_file2D,
					active_VariableName,
					active_Meshname,
					FieldY,num_FieldValuesCS_in_dim[CS],
					num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);


				//! ---------- write Z Component of Block ------------------------------------
				sprintf(active_VariableName,"%s/Block%d/%s%d_Z",
							    cs_path,
							    active_TopLevelBlk,
							    Field_Name[Field_ID],
							    active_TopLevelBlk);

				DBPutQuadvar1(silo_file2D,
					active_VariableName,
					active_Meshname,
					FieldZ,num_FieldValuesCS_in_dim[CS],
					num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);



			}//! end if vec field

		     }//! end this process

		    active_TopLevelBlk++;
		  }


		  temp_Block = temp_Block->next_Blk_of_BlockList;
	    }
	}


	//!------------------------------------
	//! 4) Build Multi-Objects
	//!------------------------------------
	DBSetDir(silo_file2D, "/");

	//! only create dir in case it does not exist 
	//! (could have been created in write CS-call different to this)
	if(!DIR_of_FieldID_exists[Field_ID])
	{
		DIR_of_FieldID_exists[Field_ID] = true;
		DBMkDir(silo_file2D,Field_Name[Field_ID]);
	}

	DBSetDir(silo_file2D, Field_Name[Field_ID]);

	char CS_name[200];
	if(CS==0) sprintf(CS_name,"X_CrossSection");
	if(CS==1) sprintf(CS_name,"Y_CrossSection");
	if(CS==2) sprintf(CS_name,"Z_CrossSection");

	DBMkDir(silo_file2D,CS_name);
	DBSetDir(silo_file2D, CS_name);

	

	//! scalar field
	if(!mpi_myRank && COMPs_FType[Field_ID]==1)
	{

	
		DBPutMultivar(silo_file2D,
			      Field_Name[Field_ID],
			      num_TopLevelBlks_CS[CS],
			      MultiObject_FieldNames,
			      MultiObject_FieldTypes,
			      NULL);
	}


	//! vector field
	if(!mpi_myRank && COMPs_FType[Field_ID]==3)
	{

		char FieldCompNameX[200];
		char FieldCompNameY[200];
		char FieldCompNameZ[200];


		sprintf(FieldCompNameX,"%s_X",Field_Name[Field_ID]);
		sprintf(FieldCompNameY,"%s_Y",Field_Name[Field_ID]);
		sprintf(FieldCompNameZ,"%s_Z",Field_Name[Field_ID]);




		//! write x comp
		DBPutMultivar(silo_file2D,
			      FieldCompNameX,
			      num_TopLevelBlks_CS[CS],
			      MultiObject_FieldNames_XComp,
			      MultiObject_FieldTypes,
			      NULL);

		//! write y comp
		DBPutMultivar(silo_file2D,
			      FieldCompNameY,
			      num_TopLevelBlks_CS[CS],
			      MultiObject_FieldNames_YComp,
			      MultiObject_FieldTypes,
			      NULL);

		//! write z comp
		DBPutMultivar(silo_file2D,
			      FieldCompNameZ,
			      num_TopLevelBlks_CS[CS],
			      MultiObject_FieldNames_ZComp,
			      MultiObject_FieldTypes,
			      NULL);

		int    types[2];
		char   vnames[2][200];
		char   defns[2][600];

		char *pvnames[2] = {vnames[0], vnames[1]};
		char *pdefns[2] = {defns[0], defns[1]};


		//! ---------------- Def Vars -------------------------------------------------
		if(CS==0)
		{
			//! for some reason absolute path WITHOUT "/" at the beginning has to be given
			sprintf(vnames[0],"%s/%s/%s_YZvec",Field_Name[Field_ID],CS_name,Field_Name[Field_ID]);

			//! define 2D Vector
			sprintf(defns[0],"{<%s/%s/%s>,<%s/%s/%s>}",Field_Name[Field_ID],CS_name,FieldCompNameY,
								   Field_Name[Field_ID],CS_name,FieldCompNameZ);
		}

		if(CS==1)
		{
			//! for some reason absolute path WITHOUT "/" at the beginning has to be given
			sprintf(vnames[0],"%s/%s/%s_XZvec",Field_Name[Field_ID],CS_name,Field_Name[Field_ID]);

			//! define 2D Vector
			sprintf(defns[0],"{<%s/%s/%s>,<%s/%s/%s>}",Field_Name[Field_ID],CS_name,FieldCompNameX,
								   Field_Name[Field_ID],CS_name,FieldCompNameZ);
		};

		if(CS==2)
		{
			//! for some reason absolute path WITHOUT "/" at the beginning has to be given
			sprintf(vnames[0],"%s/%s/%s_XYvec",Field_Name[Field_ID],CS_name,Field_Name[Field_ID]);

			//! define 2D Vector
			sprintf(defns[0],"{<%s/%s/%s>,<%s/%s/%s>}",Field_Name[Field_ID],CS_name,FieldCompNameX,
								   Field_Name[Field_ID],CS_name,FieldCompNameY);
		}


		types[0] = DB_VARTYPE_VECTOR;



		//! ---------------- def1 -------------------------------------------------
		//! define Magnitude of 3D Vector
		sprintf(defns[1],"sqrt(<%s/%s/%s>^2+<%s/%s/%s>^2+<%s/%s/%s>^2)",
					Field_Name[Field_ID],CS_name,FieldCompNameX,
					Field_Name[Field_ID],CS_name,FieldCompNameY,
					Field_Name[Field_ID],CS_name,FieldCompNameZ);

		//! for some reason absolute path WITHOUT "/" at the beginning has to be given
		sprintf(vnames[1],"%s/%s/%s_Magnitude",Field_Name[Field_ID],CS_name,Field_Name[Field_ID]);
		types[1] = DB_VARTYPE_SCALAR;


		//! write defs
		char def_name[200];
		sprintf(def_name,"%s_defvars",Field_Name[Field_ID]);
		DBPutDefvars(silo_file2D, def_name, 2, pvnames, types,pdefns, NULL);



	}




}


//!-----------------------------------------------------------------------------//
//! silo_2DWrite_MeshCS: -
//!-----------------------------------------------------------------------------//
void silo_2DWrite_MeshCS(INT32 CS)
{



	log_file << "  Writing Mesh... ";


	//!------------------------------------
	//! 1) Variable Declaration
	//!------------------------------------
	//! loop variables
	INT32 lp_ind0, lp_ind1, active_TopLevelBlk;
	char dirName[200];
	char active_Meshname[50];


	//! variable for NodeIndex to Coord Transformation
	double r_vec[3] = {0.,0.,0.};
	double cell_intern_r[3] = {0.,0.,0.};


	INT32 num_MeshValuesCS_in_dim[3][2];


	num_MeshValuesCS_in_dim[0][0] = (BlkNds_Y -hide_mGN -hide_pGN);
	num_MeshValuesCS_in_dim[0][1] = (BlkNds_Z -hide_mGN -hide_pGN);


	num_MeshValuesCS_in_dim[1][0] = (BlkNds_X -hide_mGN -hide_pGN);
	num_MeshValuesCS_in_dim[1][1] = (BlkNds_Z -hide_mGN -hide_pGN);


	num_MeshValuesCS_in_dim[2][0] = (BlkNds_X -hide_mGN -hide_pGN);
	num_MeshValuesCS_in_dim[2][1] = (BlkNds_Y -hide_mGN -hide_pGN);


	double *Mesh0, *Mesh1;
	//! Allocate one array for each component
	Mesh0 = new double[num_MeshValuesCS_in_dim[CS][0]];
	Mesh1 = new double[num_MeshValuesCS_in_dim[CS][1]];


	memset(Mesh0, 0, num_MeshValuesCS_in_dim[CS][0] *sizeof(double));
	memset(Mesh1, 0, num_MeshValuesCS_in_dim[CS][1] *sizeof(double));



	//! group component arrays to VecData array
	double *MeshData[2] = {Mesh0, Mesh1};


	char cs_path[50];
	char temp_filename[200];

	DBoptlist *optlist = DBMakeOptlist(4);
	
	DBAddOption(optlist, DBOPT_XUNITS, (void *) visit_length_unit);
	DBAddOption(optlist, DBOPT_YUNITS, (void *) visit_length_unit);

	//! Depending of CS, decide which loop indices have to be used
	switch(CS)
	{
		//! X-Cross Section: loop across Y and Z dimension
		case 0: lp_ind0 = 1;
			 lp_ind1 = 2;
			 sprintf(cs_path,"/XCS_Blocks");
			 DBAddOption(optlist, DBOPT_XLABEL, (void *) visit_ylabel);
			 DBAddOption(optlist, DBOPT_YLABEL, (void *) visit_zlabel);
			 break;

		//! Y-Cross Section: loop across X and Z dimension
		case 1: lp_ind0 = 0;
			 lp_ind1 = 2;
			 sprintf(cs_path,"/YCS_Blocks");
			 DBAddOption(optlist, DBOPT_XLABEL, (void *) visit_xlabel);
			 DBAddOption(optlist, DBOPT_YLABEL, (void *) visit_zlabel);
			 break;

		//! Z-Cross Section: loop across X and Y dimension
		case 2: lp_ind0 = 0;
			 lp_ind1 = 1;
			 sprintf(cs_path,"/ZCS_Blocks");
			 DBAddOption(optlist, DBOPT_XLABEL, (void *) visit_xlabel);
			 DBAddOption(optlist, DBOPT_YLABEL, (void *) visit_ylabel);
			 break;

	}




	//!------------------------------------------------------------------//
	//! 	3a) Master create MultiMesh names				     //
	//!	    This loop must only be called by master process !!!	     //
	//!	    since every process is writing to its own file,		     //
	//!	    the name of the file has to be specified within		     //
	//!	    the MULTIMESH OBJECT prior to the respective meshname !!!  //
	//!	    TODO:								     //
	//!	    Use limited numbre of shared files using MPIO !!!	     //
	//!------------------------------------------------------------------//
	active_TopLevelBlk=0;
	if(!mpi_myRank)
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

	
			if(   CROSS_SECTION[CS] >= temp_Block->origin[CS]
				&& CROSS_SECTION[CS] <  temp_Block->origin[CS]+Blk_Length_of[level][CS])
			{



				//! Don't specify filename in case master has to process Block
				//! since multivar and var are in the same file
				//! FOR SOME STRANGE REASON THIS HIDES ALL THE OTHER BLOCK DIRECTORIES
				//! THAT ARE VISIBLE IN CASE FILENAME IS SPECIFIED FOR MASTER ALSO !!!
				//! NOTE:
				//! IN CONTRAST TO 3D OUTPUT, THE SLASH HAS TO BE SUPPRESSED ALSO !!!!
				//! ELSE E.G. Y-CS DATA IS ATTACHED TO X-CS MESH, SFC FAILS ETC. ETC ....


				if(!temp_Block->responsible_mpi_process)
				sprintf(temp_filename,"");
				else
		   		sprintf(temp_filename,"%s_2d_p%05d_TL%05d.silo:",
		   					     Run_Name,
		   					     temp_Block->responsible_mpi_process,
		   					     TL);


				sprintf(MultiObject_FieldNames[active_TopLevelBlk],"%s%s/Block%d/Mesh%d",
												temp_filename,
												cs_path,
												active_TopLevelBlk,
												active_TopLevelBlk);
				MultiObject_FieldTypes[active_TopLevelBlk] = DB_QUAD_RECT;

				active_TopLevelBlk++;

			}
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}//! end while
	}//! end for level


	//!----------------------------------------//
	//! 	3) Write down Mesh of each block   //
	//! 	   and count number of top level   //
	//! 	   blocks	  		   //
	//!----------------------------------------//
	active_TopLevelBlk=0;
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

	
			if(        CROSS_SECTION[CS] >= temp_Block->origin[CS]
				&& CROSS_SECTION[CS] <  temp_Block->origin[CS]+Blk_Length_of[level][CS])
			{
	
			  if(mpi_myRank==temp_Block->responsible_mpi_process)
			  {

				//! use below formulation for DB_COLLINEAR and resortphysical fields
				for(ind[lp_ind0]=hide_mGN; ind[lp_ind0]<num_MeshValuesCS_in_dim[CS][0] +hide_mGN; ind[lp_ind0]++)
				{
					temp_Block->intern2normedCoords(r_vec, cell_intern_r,ind);
					Mesh0[ind[lp_ind0]-hide_mGN] = factor_scale_mesh *r_vec[lp_ind0] /*+ofst*/;
				}
				
				for(ind[lp_ind1]=hide_mGN; ind[lp_ind1]<num_MeshValuesCS_in_dim[CS][1] +hide_mGN; ind[lp_ind1]++)
				{
					temp_Block->intern2normedCoords(r_vec, cell_intern_r,ind);
					Mesh1[ind[lp_ind1]-hide_mGN] = factor_scale_mesh *r_vec[lp_ind1];
				}
			

				//! create new directory in silo file
				sprintf(dirName,"%s/Block%d", cs_path, active_TopLevelBlk);
				DBMkDir(silo_file2D,dirName);
	
				//! set output to active blocks directory
				DBSetDir(silo_file2D, dirName);
				sprintf(active_Meshname,"Mesh%d",active_TopLevelBlk);

				
				//! write mesh
				DBPutQuadmesh(silo_file2D,
						active_Meshname,
						NULL,
						MeshData,
						num_MeshValuesCS_in_dim[CS],
						num_dims,
						DB_DOUBLE,
						DB_COLLINEAR,
						optlist);

			  }
			  active_TopLevelBlk++;


		        }//! end inside CS

		  temp_Block = temp_Block->next_Blk_of_BlockList;
		}//! end while
	}//! end for level



	//!--------------------------------------------//
	//! 	4) Create Multimesh Object:		//
	//!	   - create list of mesh names/types	//
	//!        for multimesh/multivar objects	//
	//! 	   - define OptionList			//
	//! 	   - write Multimesh			//
	//!--------------------------------------------//

	//! choose field name depending on which CS to write
	char mesh_name[50];
	if(CS==0) sprintf(mesh_name,"/Mesh2D_XCS");
	if(CS==1) sprintf(mesh_name,"/Mesh2D_YCS");
	if(CS==2) sprintf(mesh_name,"/Mesh2D_ZCS");




	//! create a Multi-Block Object
	if(!mpi_myRank)
	DBPutMultimesh (silo_file2D, mesh_name, num_TopLevelBlks_CS[CS],
			MultiObject_FieldNames,
			MultiObject_FieldTypes,
			NULL);



	DBFreeOptlist(optlist);

	delete[] Mesh0;
	delete[] Mesh1;




	log_file << "done." << endl;

}




//!-------------------------------------------------------------//
//! silo_3DWrite_prepare: -								//
//!-------------------------------------------------------------//
void silo_3DWrite_prepare(void)
{




	//!--------------------------------------------------
	//! 1) - Open Silo file
	//! 	 - Create Blocks directory
	//!--------------------------------------------------

	log_file << " Writing Silo File 3D ..." << endl;
	silo_file3D = NULL;

	sprintf(silo_filename,"%s/silo_3D/%s_3d_p%05d_TL%05d.silo",data_output_path, Run_Name, mpi_myRank, TL);
	//! Open the Silo file
	//! In case visit is compiled with HDF5 support, DB_PDB can be 
	//! replaced by DB_HDF5
	
	if (COMPRESS_SILO)
	{
	  DBSetCompression("METHOD=GZIP LEVEL=8");
	}
	silo_file3D = DBCreate(silo_filename, DB_CLOBBER, DB_LOCAL,"3d" , DB_HDF5);
	if(silo_file3D == NULL)
	{
		log_file <<  "Could not create Silo file!" << endl;
		exit(0);
// 		return;
	}
	
	//! create new directory in silo file
	DBMkDir(silo_file3D,"/Blocks");

	DIR_of_trajectories_exist = false;
	DIR_of_lines_exist = false;
	DIR_of_particleTracks3D_exist = false;

	//!--------------------------------------------------
	//! 2) Specify dimensions to write depending on
	//!    whether or not GN in plus and minus
	//!    direction (pGN, mGN) are to written out
	//!--------------------------------------------------

	//! set style for DB FileOut
	if(SILO_STYLE_ZONAL)
	CENTER_STYLE = DB_ZONECENT;
	else
	CENTER_STYLE= DB_NODECENT;

	//! Set number of nodes each Block
	num_dims = 3;

	num_FieldValues_in_dim[0] = (BlkNds_X -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);
	num_FieldValues_in_dim[1] = (BlkNds_Y -hide_mGN -hide_pGN -SILO_STYLE_ZONAL),
	num_FieldValues_in_dim[2] = (BlkNds_Z -hide_mGN -hide_pGN -SILO_STYLE_ZONAL);

	//! Total number of Nodes to write each Block

	num_FieldValues = num_FieldValues_in_dim[0]*num_FieldValues_in_dim[1]*num_FieldValues_in_dim[2];


	//! Allocate one array for each component
	FieldX = new FILE_REAL[num_FieldValues];
	FieldY = new FILE_REAL[num_FieldValues];
	FieldZ = new FILE_REAL[num_FieldValues];

	FieldData[0] = FieldX;
	FieldData[1] = FieldY;
	FieldData[2] = FieldZ;



	//!-----------------------------------------//
	//!3) Count top level Blks		     //
	//!----------------------------------------//
	num_TopLevelBlks = 0;
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


			if(temp_Block->num_children<8)
			num_TopLevelBlks++;
	
			temp_Block = temp_Block->next_Blk_of_BlockList;

		}
	}

	//!--------------------------------------------------
	//! 4) free store memory below is used in write Mesh
	//!    and write Field method and will be deallocated
	//!    in function "silo_3DFlush_CleanUp".
	//!--------------------------------------------------
	MultiObject_FieldTypes = new int[num_TopLevelBlks];


	MultiObject_FieldNames = new char*[num_TopLevelBlks];
	MultiObject_FieldNames_XComp = new char*[num_TopLevelBlks];
	MultiObject_FieldNames_YComp = new char*[num_TopLevelBlks];
	MultiObject_FieldNames_ZComp = new char*[num_TopLevelBlks];

	for(INT32 blk=0; blk<num_TopLevelBlks; blk++)
	{
		MultiObject_FieldNames[blk] = new char[max_name_length];
		MultiObject_FieldNames_XComp[blk] = new char[max_name_length];
		MultiObject_FieldNames_YComp[blk] = new char[max_name_length];
		MultiObject_FieldNames_ZComp[blk] = new char[max_name_length];
	}



}



//!-------------------------------------------------------------//
//! silo_3DFlush_CleanUp: -								//
//!-------------------------------------------------------------//
void silo_3DFlush_CleanUp(void)
{





	//! Close the Silo file
	log_file << "  Flushing data ... " <<endl;
	DBClose(silo_file3D);
	log_file << "  flushed." << endl;

	//! free Field Data
	delete[] FieldX;
	delete[] FieldY;
	delete[] FieldZ;

	for(INT32 blk=0; blk<num_TopLevelBlks; blk++)
	{
		delete[] MultiObject_FieldNames[blk];
		delete[] MultiObject_FieldNames_XComp[blk];
		delete[] MultiObject_FieldNames_YComp[blk];
		delete[] MultiObject_FieldNames_ZComp[blk];
	}

	delete[] MultiObject_FieldTypes;

	delete[] MultiObject_FieldNames;
	delete[] MultiObject_FieldNames_XComp;
	delete[] MultiObject_FieldNames_YComp;
	delete[] MultiObject_FieldNames_ZComp;



	log_file << " finished." << endl <<endl;

}


//!-------------------------------------------------------------//
//! silo_3D_writeTrajectory: -						//
//!-------------------------------------------------------------//
void silo_3DWrite_Trajectory(INT32 Trajectory_ID, INT32 num_positions, D_REAL** positions)
{


	if(!DIR_of_trajectories_exist)
	{
	  DBMkDir(silo_file3D,"/Trajectories");
	  DIR_of_trajectories_exist = true;
	}


	INT32 ndims = 3;


	double *x = new double[num_positions];
	double *y = new double[num_positions];
	double *z = new double[num_positions];

	double *coords[] = {x,y,z};


	for(INT32 pos=0; pos<num_positions; pos++)
	{
		x[pos] = (positions[0][pos] *LX -Box_Origin[0]) *factor_scale_mesh;
		y[pos] = (positions[1][pos] *LY -Box_Origin[1]) *factor_scale_mesh;
		z[pos] = (positions[2][pos] *LZ -Box_Origin[2]) *factor_scale_mesh;

	}


	char trajectory_name[200];
	sprintf(trajectory_name,"/Trajectories/%s",Trajectory_FileName[Trajectory_ID]);


	DBPutPointmesh(silo_file3D, trajectory_name, ndims, coords, num_positions, DB_DOUBLE, NULL);


	delete[] x;
	delete[] y;
	delete[] z;
}

//!-------------------------------------------------------------//
//! silo_3D_writeLine: -						//
//!-------------------------------------------------------------//
void silo_3DWrite_Line(INT32 Line_ID, INT32 num_positions, D_REAL** positions)
{


	if(!DIR_of_lines_exist) 
	{
	  DBMkDir(silo_file3D,"/LineOut");
	  DIR_of_lines_exist = true;
	}


	INT32 ndims = 3;


	double *x = new double[num_positions];
	double *y = new double[num_positions];
	double *z = new double[num_positions];

	double *coords[] = {x,y,z};


	for(INT32 pos=0; pos<num_positions; pos++)
	{
		x[pos] = (positions[0][pos] *LX -Box_Origin[0]) *factor_scale_mesh;
		y[pos] = (positions[1][pos] *LY -Box_Origin[1]) *factor_scale_mesh;
		z[pos] = (positions[2][pos] *LZ -Box_Origin[2]) *factor_scale_mesh;

	}


	char line_name[200];
	sprintf(line_name,"/LineOut/%s",LineOut_FileName[Line_ID]);


	DBPutPointmesh(silo_file3D, line_name, ndims, coords, num_positions, DB_DOUBLE, NULL);


	delete[] x;
	delete[] y;
	delete[] z;
}

//!-------------------------------------------------------------//
//! silo_3DWrite_ParticleTrack: -						//
//!-------------------------------------------------------------//
void silo_3DWrite_ParticleTrack(INT32 GROUP_ID, INT32 num_positions, double** positions)
{


	INT32 lp_ind0, lp_ind1;


	char particle_track_name[200];
 	sprintf(particle_track_name, "/ParticleTracks/Group_%04d", GROUP_ID);

// 	DBoptlist *optlist = DBMakeOptlist(4);
// 	DBAddOption(optlist, DBOPT_XUNITS, (void *) visit_length_unit);
// 	DBAddOption(optlist, DBOPT_YUNITS, (void *) visit_length_unit);


	INT32 ndims = 3;
	DBPutPointmesh(silo_file3D, particle_track_name, ndims, positions, num_positions, DB_DOUBLE, NULL);



}







//!-------------------------------------------------------------//
//! silo_3DWrite_Field: -							//
//!										//
//! If no subdirecotires are used for each Block, attention has // 
//! to be drawn to the following format convention:			//
//! 										//
//! For Scalar Fields (SF) you have to provide:			//
//! 1) variable name							//
//! 2)  associated  mesh name						//
//! E.g. SF1 to Mesh1,							//
//!      SF2 to Mesh2, etc						//
//! 										//
//! For Vector Fields (VF) you have to provide			//
//! 1) variable name							//
//! 2) associated  mesh name						//
//! 3) Besides several other parameter you have to specify the	//
//! names for each  component, say vx, vy and vz.			//
//! E.g. VF1(vx,vy,vz) to Mesh1, 					//
//!      VF2(vx,vy,vz) to Mesh2, etc					//
//!										//
//! Even though this should be enough information for VISIT, it	//
//! does not work. In order to get a correct result, each		//
//! component has to be indexed with the domain number:		//
//! E.g. VF1(vx1,vy1,vz1) to Mesh1					//
//!       VF2(vx2,vy2,vz2) to Mesh2					//
//!-------------------------------------------------------------//
void silo_3DWrite_Field(INT32 Field_ID)
{


	log_file << "  Write " << Field_Name[Field_ID] <<"... ";


	//!------------------------------------
	//! 1) Variable Declaration
	//!------------------------------------
	//! loop variables
	INT32 i_j_k, i_j_k_DATA, active_TopLevelBlk;
	char active_Meshname[50];
	char active_VariableName[50];
	char temp_filename[50];
	D_REAL *VFieldX, *VFieldY, *VFieldZ;



	//! specify names for variable components
	char comp0[50], comp1[50], comp2[50];
	char* comp_names[3] = {comp0, comp1, comp2};
	sprintf(comp_names[0],"%sX",Field_Name[Field_ID]);
	sprintf(comp_names[1],"%sY",Field_Name[Field_ID]);
	sprintf(comp_names[2],"%sZ",Field_Name[Field_ID]);



	//!-------------------------------------------------------------
	//!2) Meta Data for Visit optimization
	//!-------------------------------------------------------------
	//! Min/Max Values to accelerate contur plots
	double **extents = new double*[num_TopLevelBlks];
	for(INT32 blk=0; blk<num_TopLevelBlks; blk++)
	extents[blk] = new double[2];



	//!------------------------------------------------------------------//
	//! 	3a) Master create MultiMesh names				     //
	//!	    This loop must only be called by master process !!!	     //
	//!	    since every process is writing to its own file,		     //
	//!	    the name of the file has to be specified within		     //
	//!	    the MULTIMESH OBJECT prior to the respective meshname !!!  //
	//!	    TODO:								     //
	//!	    Use limited numbre of shared files using MPIO !!!	     //
	//!------------------------------------------------------------------//
	active_TopLevelBlk=0;
	if(!mpi_myRank)
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

	   CBlock *temp_Block = BlockList_of_Lev[level];
	   while(temp_Block)
	   {

		  if(temp_Block->num_children<8)
		  {



			//! Don't specify filename in case master has to process Block
			//! since multivar and var are in the same file
			//! FOR SOME STRANGE REASON THIS HIDES ALL THE OTHER BLOCK DIRECTORIES
			//! THAT ARE VISIBLE IN CASE FILENAME IS SPECIFIED FOR MASTER ALSO !!!
			if(!temp_Block->responsible_mpi_process)
			sprintf(temp_filename,"");
			else
		   	sprintf(temp_filename,"%s_3d_p%05d_TL%05d.silo:",
		   					     Run_Name,
		   					     temp_Block->responsible_mpi_process,
		   					     TL);


			//! Types of Block MultiVar
			MultiObject_FieldTypes[active_TopLevelBlk] = DB_QUADVAR;
	
			//!-----------------------------------------//
			//!  write scalar field			   //
			//!-----------------------------------------//
			if(COMPs_FType[Field_ID]==1)
			{
	
	
					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames[active_TopLevelBlk],"%s/Blocks/Block%d/%s%d",
											temp_filename,
											active_TopLevelBlk,
											Field_Name[Field_ID],
											active_TopLevelBlk);
		
		
	
		
					active_TopLevelBlk++;
	
			}

	
			//!-----------------------------------------//
			//!  write vector field				   //
			//!----------------------------------------//
			if(COMPs_FType[Field_ID]==3)
			{
	
	
					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames_XComp[active_TopLevelBlk],"%s/Blocks/Block%d/%s%d_X",
												temp_filename,
												active_TopLevelBlk,
												Field_Name[Field_ID],
												active_TopLevelBlk);
	
					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames_YComp[active_TopLevelBlk],"%s/Blocks/Block%d/%s%d_Y",
												temp_filename,
												active_TopLevelBlk,
												Field_Name[Field_ID],
												active_TopLevelBlk);
	
					//! MultiVarList of Names/Types of block variable
					sprintf(MultiObject_FieldNames_ZComp[active_TopLevelBlk],"%s/Blocks/Block%d/%s%d_Z",
												temp_filename,
												active_TopLevelBlk,
												Field_Name[Field_ID],
												active_TopLevelBlk);
	
					active_TopLevelBlk++;
	
			}


		  }

		  temp_Block = temp_Block->next_Blk_of_BlockList;
	   }

	}



	//!-----------------------------------------//
	//!3) Loop over all Blks			   //
	//!----------------------------------------//
	active_TopLevelBlk=0;
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{



		  if(temp_Block->num_children<8)
		  {

		    if(mpi_myRank==temp_Block->responsible_mpi_process)
		    {



			//! Name of Block Variables Mesh
			sprintf(active_Meshname,"/Blocks/Block%d/Mesh%d",active_TopLevelBlk, active_TopLevelBlk);


			//!-----------------------------------------//
			//!  write scalar field				   //
			//!----------------------------------------//
			if(COMPs_FType[Field_ID]==1)
			{
				D_REAL min_Value =  1.e9;
				D_REAL max_Value = -1.e9;
				
				for(ind[0]=hide_mGN; ind[0]<BlkNds_X -hide_pGN -SILO_STYLE_ZONAL; ind[0]++)
				 for(ind[1]=hide_mGN; ind[1]<BlkNds_Y -hide_pGN -SILO_STYLE_ZONAL; ind[1]++)
				  for(ind[2]=hide_mGN; ind[2]<BlkNds_Z -hide_pGN -SILO_STYLE_ZONAL; ind[2]++)
				  {
				
	
					i_j_k  = ind[0]*BlkNds_Y*BlkNds_Z 
						+ind[1]*BlkNds_Z 
						+ind[2];
	
					i_j_k_DATA  = (ind[2]-hide_mGN)*num_FieldValues_in_dim[1]*num_FieldValues_in_dim[0]
						     +(ind[1]-hide_mGN)*num_FieldValues_in_dim[0]
						     +(ind[0]-hide_mGN);
	
					FieldX[i_j_k_DATA] = temp_Block->Field_Type[Field_ID][i_j_k];
					
					if(FieldX[i_j_k_DATA]>max_Value) max_Value = FieldX[i_j_k_DATA];
					if(FieldX[i_j_k_DATA]<min_Value) min_Value = FieldX[i_j_k_DATA];
				  }
	

				//! extents are applied to multivar rather than Quadvar,
				//! so no optlist below
				extents[active_TopLevelBlk][0] = min_Value;
				extents[active_TopLevelBlk][1] = max_Value;


				sprintf(active_VariableName,"/Blocks/Block%d/%s%d",
							   active_TopLevelBlk,
							   Field_Name[Field_ID],
							   active_TopLevelBlk);
	


				DBPutQuadvar1(silo_file3D,
					      active_VariableName,
					      active_Meshname,
					      FieldX,num_FieldValues_in_dim,
					      num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);



		   	 }//! end if scalar field


			//!-----------------------------------------//
			//!  write vector field				   //
			//!----------------------------------------//
			if(COMPs_FType[Field_ID]==3)
			{

				D_REAL min_Value =  1.e9;
				D_REAL max_Value = -1.e9;
				
				D_REAL vec[3], magnitude;

				//! Write Vector
				VFieldX = temp_Block->Field_Type[Field_ID];
				VFieldY = VFieldX +num_nodes_in_block;
				VFieldZ = VFieldY +num_nodes_in_block;

				for(ind[0]=hide_mGN; ind[0]<BlkNds_X -hide_pGN -SILO_STYLE_ZONAL; ind[0]++)
				 for(ind[1]=hide_mGN; ind[1]<BlkNds_Y -hide_pGN -SILO_STYLE_ZONAL; ind[1]++)
				  for(ind[2]=hide_mGN; ind[2]<BlkNds_Z -hide_pGN -SILO_STYLE_ZONAL; ind[2]++)
				  {

					i_j_k  = ind[0]*BlkNds_Y*BlkNds_Z 
						+ind[1]*BlkNds_Z 
						+ind[2];
	
					i_j_k_DATA  = (ind[2]-hide_mGN)*num_FieldValues_in_dim[1]*num_FieldValues_in_dim[0]
						     +(ind[1]-hide_mGN)*num_FieldValues_in_dim[0]
						     +(ind[0]-hide_mGN);

				
					FieldX[i_j_k_DATA] = vec[0] = VFieldX[i_j_k];
					FieldY[i_j_k_DATA] = vec[1] = VFieldY[i_j_k];
					FieldZ[i_j_k_DATA] = vec[2] = VFieldZ[i_j_k];
					
					magnitude = vec_len(vec);
					if(magnitude>max_Value) max_Value = magnitude;
					if(magnitude<min_Value) min_Value = magnitude;
					
				  }

				//! extents are applied to multivar rather than Quadvar,
				//! so no optlist below
				extents[active_TopLevelBlk][0] = min_Value;
				extents[active_TopLevelBlk][1] = max_Value;


				//! ---------- write X Component of Block ------------------------------------
				sprintf(active_VariableName,"/Blocks/Block%d/%s%d_X",
							    active_TopLevelBlk,
							    Field_Name[Field_ID],
							    active_TopLevelBlk);

				DBPutQuadvar1(silo_file3D,
					active_VariableName,
					active_Meshname,
					FieldX, num_FieldValues_in_dim,
					num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);



				//! ---------- write Y Component of Block ------------------------------------
				sprintf(active_VariableName,"/Blocks/Block%d/%s%d_Y",
							    active_TopLevelBlk,
							    Field_Name[Field_ID],
							    active_TopLevelBlk);

				DBPutQuadvar1(silo_file3D,
					active_VariableName,
					active_Meshname,
					FieldY, num_FieldValues_in_dim,
					num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);

				//! ---------- write Z Component of Block ------------------------------------
				sprintf(active_VariableName,"/Blocks/Block%d/%s%d_Z",
							    active_TopLevelBlk,
							    Field_Name[Field_ID],
							    active_TopLevelBlk);


				DBPutQuadvar1(silo_file3D,
					active_VariableName,
					active_Meshname,
					FieldZ, num_FieldValues_in_dim,
					num_dims, NULL, 0, DB_FLOAT, CENTER_STYLE, NULL);


			}//! end if vec field
			

		     }//! mpi_myRank==temp_Block->responsible_mpi_process
		     active_TopLevelBlk++;


		  }//! end for block array
		  temp_Block = temp_Block->next_Blk_of_BlockList;

		}//! end while
	}//! end for level



	//!------------------------------------
	//! 4) Build Multi-Objects
	//!------------------------------------

	const int two = 2;
	DBoptlist *MultiVarOptlist = NULL;
	MultiVarOptlist = DBMakeOptlist(2);
	DBAddOption(MultiVarOptlist, DBOPT_EXTENTS_SIZE, (void *)&two);
	DBAddOption(MultiVarOptlist, DBOPT_EXTENTS, (void *)extents);

	DBSetDir(silo_file3D, "/");

	//! MASTER && scalar field
	if(!mpi_myRank && COMPs_FType[Field_ID]==1)
	{

		DBPutMultivar(silo_file3D,
			Field_Name[Field_ID],
			num_TopLevelBlks,
			MultiObject_FieldNames,
			MultiObject_FieldTypes,
			MultiVarOptlist);

	}

	//! MASTER && vector field
	if(!mpi_myRank && COMPs_FType[Field_ID]==3)
	{

		//! create new directory in silo file
		DBMkDir(silo_file3D, Field_Name[Field_ID]);
		DBSetDir(silo_file3D, Field_Name[Field_ID]);

		char FieldCompNameX[200];
		char FieldCompNameY[200];
		char FieldCompNameZ[200];


		sprintf(FieldCompNameX,"%s_X",Field_Name[Field_ID]);
		sprintf(FieldCompNameY,"%s_Y",Field_Name[Field_ID]);
		sprintf(FieldCompNameZ,"%s_Z",Field_Name[Field_ID]);


		//! use vector min/max extents for components also (upper limit for components)

		//! write x comp
		DBPutMultivar(silo_file3D,
			      FieldCompNameX,
			      num_TopLevelBlks,
			      MultiObject_FieldNames_XComp,
			      MultiObject_FieldTypes,
			      MultiVarOptlist);

		//! write y comp
		DBPutMultivar(silo_file3D,
			      FieldCompNameY,
			      num_TopLevelBlks,
			      MultiObject_FieldNames_YComp,
			      MultiObject_FieldTypes,
			      MultiVarOptlist);

		//! write z comp
		DBPutMultivar(silo_file3D,
			      FieldCompNameZ,
			      num_TopLevelBlks,
			      MultiObject_FieldNames_ZComp,
			      MultiObject_FieldTypes,
			      MultiVarOptlist);

		int    types[2];
		char   vnames[2][200];
		char   defns[2][600];

		char *pvnames[2] = {vnames[0], vnames[1]};
		char *pdefns[2] = {defns[0], defns[1]};

		types[0] = DB_VARTYPE_VECTOR;

		//! for some reason absolute path WITHOUT "/" at the beginning has to be given
		sprintf(vnames[0],"%s/%s_vec",Field_Name[Field_ID],Field_Name[Field_ID]);

		//! define 2D Vector
		sprintf(defns[0],"{<%s/%s>,<%s/%s>,<%s/%s>}",Field_Name[Field_ID],FieldCompNameX,
							     Field_Name[Field_ID],FieldCompNameY,
							     Field_Name[Field_ID],FieldCompNameZ);

		//! write defs
		char def_name[200];
		sprintf(def_name,"%s_defvars",Field_Name[Field_ID]);
		DBPutDefvars(silo_file3D, def_name, 1, pvnames, types,pdefns, NULL);


	}

	//! Free dynamical memory
	DBFreeOptlist(MultiVarOptlist);

	for(INT32 blk=0; blk<num_TopLevelBlks; blk++)
	delete extents[blk];
	delete extents;


	log_file << "done." << endl;

}


//!-----------------------------------------------------------------------------//
//! silo_3DWrite_Mesh: -
//!-----------------------------------------------------------------------------//
void silo_3DWrite_Mesh(void)
{



	log_file << "  Writing Mesh... ";


	//!------------------------------------
	//! 1) Variable Declaration
	//!------------------------------------
	//! loop variables
	INT32 active_TopLevelBlk;
	char dirName[200];
	char temp_filename[200];
	
	char active_Meshname[50];


	//! variable for NodeIndex to Coord Transformation
	double r_vec[3] = {0.,0.,0.};;
	double cell_intern_r[3] = {0.,0.,0.};



	INT32 num_MeshValues_in_dim[3];
	num_MeshValues_in_dim[0] = (BlkNds_X -hide_mGN -hide_pGN);
	num_MeshValues_in_dim[1] = (BlkNds_Y -hide_mGN -hide_pGN),
	num_MeshValues_in_dim[2] = (BlkNds_Z -hide_mGN -hide_pGN);

	//! Total number of Nodes to write each Block
	INT32 num_MeshValues;
	num_MeshValues  = num_MeshValues_in_dim[0] *num_MeshValues_in_dim[1] *num_MeshValues_in_dim[2];


	double *MeshX, *MeshY, *MeshZ;
	//! Allocate one array for each component
	MeshX = new double[num_MeshValues];
	MeshY = new double[num_MeshValues];
	MeshZ = new double[num_MeshValues];

	//! group component arrays to VecData array
	double *MeshData[3] = {MeshX, MeshY, MeshZ};



	//!------------------------------------------------------------------//
	//! 	3a) Master create MultiMesh names				     //
	//!	    This loop must only be called by master process !!!	     //
	//!	    since every process is writing to its own file,		     //
	//!	    the name of the file has to be specified within		     //
	//!	    the MULTIMESH OBJECT prior to the respective meshname !!!  //
	//!	    TODO:								     //
	//!	    Use limited numbre of shared files using MPIO !!!	     //
	//!------------------------------------------------------------------//
	active_TopLevelBlk=0;
	if(!mpi_myRank)
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

	   CBlock *temp_Block = BlockList_of_Lev[level];
	   while(temp_Block)
	   {


		  if(temp_Block->num_children<8)
		  {


			//! don't specify filename in case master has to process Block
			//! since multivar and var are in the same file
			//! FOR SOME STRANGE REASON THIS HIDES ALL THE OTHER BLOCK DIRECTORIES
			//! THAT ARE VISIBLE IN CASE FILENAME IS SPECIFIED FOR MASTER ALSO !!!
			if(!temp_Block->responsible_mpi_process)
			sprintf(temp_filename,"");
			else
		   	sprintf(temp_filename,"%s_3d_p%05d_TL%05d.silo:",
		   					     Run_Name,
		   					     temp_Block->responsible_mpi_process,
		   					     TL);


			//! Master set write multimesh names
			sprintf(MultiObject_FieldNames[active_TopLevelBlk],"%s/Blocks/Block%d/Mesh%d",
											temp_filename,
											active_TopLevelBlk,
											active_TopLevelBlk);


		  	MultiObject_FieldTypes[active_TopLevelBlk] = DB_QUAD_RECT;

		  	active_TopLevelBlk++;

		  }

		  temp_Block = temp_Block->next_Blk_of_BlockList;
	   }

	}

	//!----------------------------------------//
	//! 	3) Write down Mesh of each block    //
	//! 	   and count number of top level    //
	//! 	   blocks	  			    //
	//!---------------------------------------//
	active_TopLevelBlk=0;
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{

		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{


		  if(temp_Block->num_children<8)
		  {
	
		    if(mpi_myRank==temp_Block->responsible_mpi_process)
		    {

			//! use below formulation for DB_COLLINEAR and resortphysical fields
			for(ind[0]=hide_mGN; ind[0]<BlkNds_X-hide_pGN; ind[0]++)
			{
				temp_Block->intern2normedCoords(r_vec, cell_intern_r,ind);
				MeshX[ind[0]-hide_mGN] = factor_scale_mesh * r_vec[0];
			}
			
			for(ind[1]=hide_mGN; ind[1]<BlkNds_Y-hide_pGN; ind[1]++)
			{
				temp_Block->intern2normedCoords(r_vec, cell_intern_r,ind);
				MeshY[ind[1]-hide_mGN] = factor_scale_mesh * r_vec[1];
			}
			
			for(ind[2]=hide_mGN; ind[2]<BlkNds_Z-hide_pGN; ind[2]++)
			{
				temp_Block->intern2normedCoords(r_vec, cell_intern_r,ind);
				MeshZ[ind[2]-hide_mGN] = factor_scale_mesh * r_vec[2];
			}



			//! create new directory in silo file
			sprintf(dirName,"/Blocks/Block%d",active_TopLevelBlk);
			DBMkDir(silo_file3D,dirName);

			//! set output to active blocks directory
			DBSetDir(silo_file3D,dirName);
			sprintf(active_Meshname,"Mesh%d",active_TopLevelBlk);


			//! write mesh
			DBPutQuadmesh(silo_file3D, active_Meshname, NULL, MeshData, num_MeshValues_in_dim, num_dims,DB_DOUBLE, DB_COLLINEAR, NULL);



		   }//! end if responsible_mpi_process
		   active_TopLevelBlk++;

		  }//! end for block array

		  temp_Block = temp_Block->next_Blk_of_BlockList;
		}//! end while
	}//! end for level



	//!--------------------------------------------//
	//! 	2) Create Multimesh Object:			//
	//!	   - create list of mesh names/types	//
	//!        for multimesh/multivar objects	//
	//! 	   - define OptionList				//
	//! 	   - write Multimesh				//
	//!--------------------------------------------//


	//! Create an option list to contain labels and units.
// 	DBoptlist *optlist = DBMakeOptlist(6);
// 	DBAddOption(optlist, DBOPT_XLABEL, (void *)"test");
// 	DBAddOption(optlist, DBOPT_XUNITS, (void *)"test");
// 	DBAddOption(optlist, DBOPT_YLABEL, (void *)"testtest");
// 	DBAddOption(optlist, DBOPT_YUNITS, (void *)"test");
// 	DBAddOption(optlist, DBOPT_ZLABEL, (void *)"test");
// 	DBAddOption(optlist, DBOPT_ZUNITS, (void *)"test");


// 	const int six = 6;
// 	DBoptlist *optlist = DBMakeOptlist(2);
// 	DBAddOption(optlist, DBOPT_EXTENTS_SIZE, (void *)&six);
// 	DBAddOption(optlist, DBOPT_EXTENTS, (void *)spatial_extents);

	//! create a Multi-Block Object
	if(!mpi_myRank)
	DBPutMultimesh (silo_file3D, "/Mesh3D", num_TopLevelBlks,
			MultiObject_FieldNames,
			MultiObject_FieldTypes,
			NULL);


	delete[] MeshX;
	delete[] MeshY;
	delete[] MeshZ;


	//! Free the option list.
// 	DBFreeOptlist(optlist);

// 	for(INT32 blk=0; blk<num_TopLevelBlks; blk++)
// 	delete spatial_extents[blk];
// 	delete spatial_extents;

	log_file << "done." << endl;

}


//!--------------------------------------------------------
//! convert_particle_tracks_to_silo_mesh
//!--------------------------------------------------------
void convert_particle_tracks_to_silo_mesh(void)
{

if (COMPRESS_SILO)
	{
	  DBSetCompression("METHOD=GZIP LEVEL=8");
	}

	
	log_file << endl;
	log_file << " CONVERTING PARTICLE-DATA TO SILO-MESH ... " << endl;

	if(!binary_particle_tracks)
	{

		log_file << " binary_particle_tracks required ... " << endl;
		log_file << " canceling." << endl;
		return;

	}





	char silo_filename[200];

	//! open file, alloc memory etc.
	//! Open the Silo file
	//! In case visit is compiled with HDF5 support, DB_PDB can be 
	//! replaced by DB_HDF5
	

	//! ------- 2D ---------------------------------------------------------
	if(TL_OUTPUT_2D_SILO && TL%TL_OUTPUT_2D_SILO==0)
	{
		sprintf(silo_filename,"%s/silo/%s_2d_p%05d_TL%05d.silo",data_output_path, Run_Name, mpi_myRank, TL);
		silo_file2D = DBOpen(silo_filename, DB_HDF5, DB_APPEND);
		if(silo_file2D == NULL)
		{
			silo_file2D = DBCreate(silo_filename, DB_CLOBBER, DB_LOCAL,"2d" , DB_HDF5);


			if(silo_file2D == NULL)
			{
				log_file <<  "Could not open nor create Silo file!" << endl;
				return;
			}
		
		}

		
	}

	//! ------- 3D ---------------------------------------------------------
	if(TL_OUTPUT_3D_SILO && TL%TL_OUTPUT_3D_SILO==0)
	{
		sprintf(silo_filename,"%s/silo_3D/%s_3d_p%05d_TL%05d.silo",data_output_path, Run_Name, mpi_myRank, TL);
		silo_file3D = DBOpen(silo_filename, DB_HDF5, DB_APPEND);
		if(silo_file3D == NULL)
		{

			silo_file3D = DBCreate(silo_filename, DB_CLOBBER, DB_LOCAL,"3d" , DB_HDF5);


			if(silo_file3D == NULL)
			{
	
				log_file <<  "Could not open nor create Silo file!" << endl;
				return;
			}
		
		}
	}

	if(!DIR_of_particleTracks2D_exist && TL_OUTPUT_2D_SILO && TL%TL_OUTPUT_2D_SILO==0)
	{
		DBMkDir(silo_file2D,"/ParticleTracks");
		DBMkDir(silo_file2D,"/ParticleTracks/X_CrossSection");
		DBMkDir(silo_file2D,"/ParticleTracks/Y_CrossSection");
		DBMkDir(silo_file2D,"/ParticleTracks/Z_CrossSection");
		
		DIR_of_particleTracks2D_exist = true;
		
	}
	
	if(!DIR_of_particleTracks3D_exist && TL_OUTPUT_3D_SILO && TL%TL_OUTPUT_3D_SILO==0)
	{
		DBMkDir(silo_file3D,"/ParticleTracks");

		DIR_of_particleTracks3D_exist = true;
	
	}
	
	//! assume particles to be in s0
	INT32 species = 0;

	char buffer[300];
	char Particle_inFileName[300];
	INT64 num_bytes_this_track, num_pos_this_track, num_pos_this_group;
	ifstream Particle_inFile;

// 	INT64 time_level = 0;
	PARTICLE_REAL x[3];


	double *x_silo;
	double *y_silo;
	double *z_silo;




	INT32 num_particle_tracks = num_marked_particle;

	log_file << "num_particle_tracks: " << num_particle_tracks << endl;

	INT32 num_groups = num_particle_tracks/num_tracks_each_group;

	log_file << "num_groups: " << num_groups << endl;

	//! only master shall convert files
	if(!mpi_myRank)
	//! each ion has its on file
	for (INT32 group=0; group < num_groups-1; group++)
	{

		num_pos_this_group = 0;

		for(INT32 track=0; track<num_tracks_each_group; track++)
		{

			sprintf(Particle_inFileName,"%s/particle_tracks/spec%d_part%07d.txt",data_output_path, species,  track +group*num_tracks_each_group);
			Particle_inFile.open(Particle_inFileName);


			//! open file in which particle data shall
			//! be written (of all processes)

// 			if(!Particle_inFile);
// 			{
// 
// 				log_file << " File not found: " << endl;
// 				log_file << " Particle_inFileName: " << Particle_inFileName << endl;
// 
// 			}


			//! count number of bytes in respective Particle File
			Particle_inFile.seekg (0, ios::end);
			num_bytes_this_track = Particle_inFile.tellg();
			Particle_inFile.close();



			num_pos_this_track = num_bytes_this_track / (3*sizeof(PARTICLE_REAL));

			if(num_bytes_this_track % (3*sizeof(PARTICLE_REAL)) != 0)
			{
				log_file << " ERROR in convert particle track to silo mesh " << endl;
				log_file << " file size no multiple of 3*sizeof(PARTICLE_REAL)" << endl;

				log_file << endl;
				log_file << " Particle_inFileName: " << Particle_inFileName << endl;
				log_file << " num_bytes_this_track: " << num_bytes_this_track << endl;
				log_file << " 3*sizeof(PARTICLE_REAL): " << 3*sizeof(PARTICLE_REAL) << endl;


				log_file << " EXITING ...." << endl;
				exit(1);
			}

			num_pos_this_group += num_pos_this_track;

		}
	

// 		log_file << " group:              " << group << endl;
// 		log_file << " num_pos_this_group: " << num_pos_this_group << endl;
// 		log_file << endl;


		x_silo = new double[num_pos_this_group];
		y_silo = new double[num_pos_this_group];
		z_silo = new double[num_pos_this_group];
		double* positions[3] = {x_silo, y_silo, z_silo};

		INT64 pos_in_group = 0;

		for(INT32 track=0; track<num_tracks_each_group; track++)
		{

// 			log_file << " track: "  << track << endl;


			sprintf(Particle_inFileName,"%s/particle_tracks/spec%d_part%07d.txt",data_output_path, species, track +group*num_tracks_each_group);
			Particle_inFile.open(Particle_inFileName);

			Particle_inFile.seekg (0, ios::end);
			num_bytes_this_track = Particle_inFile.tellg();

			if(num_bytes_this_track % (3*sizeof(PARTICLE_REAL)) != 0 || num_bytes_this_track==0)
			{
				log_file << " ERROR in convert particle track to silo mesh " << endl;
				log_file << " file size no multiple of 3*sizeof(PARTICLE_REAL)" << endl;

				log_file << endl;
				log_file << " Particle_inFileName: " << Particle_inFileName << endl;
				log_file << " num_bytes_this_track: " << num_bytes_this_track << endl;
				log_file << " 3*sizeof(PARTICLE_REAL): " << 3*sizeof(PARTICLE_REAL) << endl;


				log_file << " EXITING ...." << endl;
				exit(1);
			}

			num_pos_this_track = num_bytes_this_track / (3*sizeof(PARTICLE_REAL));

			Particle_inFile.seekg (0, ios::beg);
			for(INT32 pos=0; pos<num_pos_this_track; pos++)
			{



				Particle_inFile.read(reinterpret_cast<char*> (x),3*sizeof(PARTICLE_REAL));

				x_silo[pos_in_group] = x[0]*factor_scale_mesh;
				y_silo[pos_in_group] = x[1]*factor_scale_mesh;
				z_silo[pos_in_group] = x[2]*factor_scale_mesh;

				pos_in_group++;
			}

			Particle_inFile.close();


		}




		if(pos_in_group  != num_pos_this_group)
		{
			log_file << " ERROR in convert particle track to silo mesh " << endl;
			log_file << " pos_in_group  != num_pos_this_group" << endl;
			log_file << " EXITING ...." << endl;
			exit(1);
		}



		if(TL_OUTPUT_2D_SILO && TL%TL_OUTPUT_2D_SILO==0)
		{
			silo_2DWrite_ParticleTrack(group, num_pos_this_group, positions, 0);
			silo_2DWrite_ParticleTrack(group, num_pos_this_group, positions, 1);
			silo_2DWrite_ParticleTrack(group, num_pos_this_group, positions, 2);
		}

		if(TL_OUTPUT_3D_SILO && TL%TL_OUTPUT_3D_SILO==0)
		silo_3DWrite_ParticleTrack(group, num_pos_this_group, positions);

		delete[] x_silo;
		delete[] y_silo;
		delete[] z_silo;
		

	}//! end for group




	//! close file, delete memory etc



	if(TL_OUTPUT_2D_SILO && TL%TL_OUTPUT_2D_SILO==0)
	DBClose(silo_file2D);
	
	if(TL_OUTPUT_3D_SILO && TL%TL_OUTPUT_3D_SILO==0)
	DBClose(silo_file3D);

	log_file << " Done." << endl;
	log_file << endl;

}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
		/*
			//! Save mesh nodes DB_NONCOLLINEAR:
			//! - Usually a DB_COLLINEAR mesh would be sufficient. However, using it
			//!   needs reordering of all Var-Values (it seems as silo sorts different
			//!   than x-Plane,y-Colums, z-Values).
			//! - The NONCOLLINEAR way choosen below little overhed data is produced which
			//!   should be negligible compared to all the Physical Fields data.
			//! - Especially it seems Visit can process COLLINEAR and NONCOLLINEAR Mesh-Data
			//!   equally fast.
			//! - See end of this file for collinear version.
			for(ind[0]=hide_mGN; ind[0]<BlkNds_X-hide_pGN; ind[0]++)
			 for(ind[1]=hide_mGN; ind[1]<BlkNds_Y-hide_pGN; ind[1]++)
			  for(ind[2]=hide_mGN; ind[2]<BlkNds_Z-hide_pGN; ind[2]++)
			  {
			

				//! logical order in silo format differs from order 
				//! in adaptive_hybrid code:

				//! adaptive_hybrid:
				//! z-values, y-clolumns, x-planes

				//! silo:
				//! x-values, y-clolumns, z-planes



				i_j_k  = ind[0]*BlkNds_Y*BlkNds_Z
					+ind[1]*BlkNds_Z
					+ind[2];


				i_j_k_DATA  = (ind[2]-hide_mGN)*num_MeshValues_in_dim[1]*num_MeshValues_in_dim[0]
					     +(ind[1]-hide_mGN)*num_MeshValues_in_dim[0]
					     +(ind[0]-hide_mGN);

				//! calc coords of respective node 
				temp_Block->intern2normedCoords(r_vec, cell_intern_r,ind);

				MeshX[i_j_k_DATA] = r_vec[0];
				MeshY[i_j_k_DATA] = r_vec[1];
				MeshZ[i_j_k_DATA] = r_vec[2];

			  }*/


	//! using expression below does not allow to acces individual componets in visu
				/*
				DBPutQuadvar(silo_file,			//! Database ﬁle pointer
					      active_VariableName,	//! Name of the variable
					      active_Meshname,		//! Name of the mesh associated with this variable
					      				//! (written with DBPutQuadmesh or DBPutUcdmesh). If
					      				//! no association is to be made, this value should
					      				//! be NULL.
					      3,			//! nvars: number of sub-variables which comprise
					      				//! this variable. For a scalar array, this is 
					      				//! one. If writing a vector quantity, however, this
					      				//! would be two for a 2-D vector and three for a
					      				//! 3-D vector.
					      comp_names,		//! Array of length "nvars" containing pointers to
					       				//! character strings deﬁning the names associated
					       				//! with each sub-variable.
					      FieldData,		//! Array of length "nvars" containing pointers to
									//! arrays deﬁning the values associated with each
									//! subvariable
					      num_FieldValuesCS_in_dim[CS], //! Array of length ndims which describes the
					      				//! dimensionality of the variable.
					      num_dims,			//! Number of dimensions.
					      NULL,			//! Array of length nvars containing pointers to arrays
					      				//! deﬁning the mixed-data values associated with each
					      				//! subvariable. If no mixed values are present, this
					      				//! should be NULL.
					      0,			//! Length of mixed data arrays, if provided.
					      DB_FLOAT,			//! Datatype of the variable. One of the predeﬁned
									//! Silo data types.
					      CENTER_STYLE,		//! Centering of the sub-variables on the associated
					      				//! mesh.
					      NULL);		//! Pointer to an option list structure
				*/
