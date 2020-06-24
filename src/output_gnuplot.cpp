
#include "CBlk.h"
#include "parameters.h"
#include "output_gnuplot.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>
#include <math.h>


extern D_REAL **delta_of_L;




//!-------------------------------------------------------------//
//! gnuplot_xlayer: -							//
//!-------------------------------------------------------------//
void ascii_output_3D(INT32 TL)
{


	log_file << endl << " Writing ASCII 3D Data ... " << endl;

	//! For ascii_output_3D output it is assumed that there is only
	//! a single Block in level 0 and no further refinements
	CBlock* ascii3d_Block = BlockList_of_Lev[0];

	write_ascii3D_field(ascii3d_Block, id_BTotal, TL);
	write_ascii3D_field(ascii3d_Block, id_EField, TL);

	log_file << " writing finished." << endl << endl;

}


//!-------------------------------------------------------------//
//! gnuplot_xlayer: -							//
//!-------------------------------------------------------------//
void write_ascii3D_field(CBlock* block, INT32 Field_Type, INT32 TL)
{

	INT32 i_j_k;
	char filename[200];

	INT32 start_elemet;
	D_REAL* Field = block->Field_Type[Field_Type];

	ofstream ascii3D_outfile;

	for(INT32 comp=0; comp<COMPs_FType[Field_Type]; comp++)
	{

		start_elemet = comp *num_nodes_in_block;


		char filename[200];
		sprintf(filename,"%s/gnuplot/%s_%s_comp%d_3d_TL%d.ascii3d",data_output_path, Run_Name, Field_Name[Field_Type], comp, TL);
		ascii3D_outfile.open( filename);

		for (INT32 k=0; k<BlkNds_Z; k++)
		{
			for (INT32 j=0; j<BlkNds_Y; j++)
			{
				for (INT32 i=0; i<BlkNds_X; i++)
				{
			
			
					i_j_k  = i*BlkNds_Y*BlkNds_Z 
						+j*BlkNds_Z 
						+k;
			
					ascii3D_outfile << Field[start_elemet +i_j_k] <<"	";
			
				}
				ascii3D_outfile << endl;
				
			}
			ascii3D_outfile << endl;
		}

		ascii3D_outfile.close();
	}
}

//!-------------------------------------------------------------//
//! gnuplot_xlayer: -							//
//!-------------------------------------------------------------//
void write_field(CBlock* block, INT32* lay, INT32 Field_Type, INT32 TL)
{


	for(INT32 comp=0; comp<COMPs_FType[Field_Type]; comp++)
	{

		gnuplot_xlayer(block->Field_Type[Field_Type],Field_Name[Field_Type],comp,lay[0],TL);
		gnuplot_ylayer(block->Field_Type[Field_Type],Field_Name[Field_Type],comp,lay[1],TL);
		gnuplot_zlayer(block->Field_Type[Field_Type],Field_Name[Field_Type],comp,lay[2],TL);

	}

}


//!-------------------------------------------------------------//
//! gnuplot_xlayer: -								//
//!-------------------------------------------------------------//
void gnuplot_write_inner(bool* Flag, INT32 xlay, INT32 ylay, INT32 zlay)
{

	D_REAL* inner_field = new D_REAL[num_nodes_in_block];


	for(INT32 i=0; i<num_nodes_in_block; i++)
	inner_field[i] = D_REAL(Flag[i]);


	char field_name[100] = "inner";
	gnuplot_xlayer(inner_field,field_name,0,xlay,0);
	gnuplot_ylayer(inner_field,field_name,0,ylay,0);
	gnuplot_zlayer(inner_field,field_name,0,zlay,0);


	delete inner_field;



}

//!-------------------------------------------------------------//
//! gnuplot_write_grid: -							//
//!-------------------------------------------------------------//
void gnuplot_write_grid(CBlock *temp_blk, INT32 xlay, INT32 ylay, INT32 zlay)
{


	INT32 i_j_k;
	ofstream Outfile;
	char filename[200];


	D_REAL* GridX = new D_REAL[num_nodes_in_block];
	D_REAL* GridY = new D_REAL[num_nodes_in_block];
	D_REAL* GridZ = new D_REAL[num_nodes_in_block];


	for (INT32 i=0; i<BlkNds_X; i++)
	 for (INT32 j=0; j<BlkNds_Y; j++)
	  for (INT32 k=0; k<BlkNds_Z; k++)
	  {


		i_j_k  = i*BlkNds_Y*BlkNds_Z 
			+j*BlkNds_Z 
			+k;

		GridX[i_j_k] = temp_blk->origin[0] + (i-1) *delta_of_L[temp_blk->RLevel][0];
		GridY[i_j_k] = temp_blk->origin[1] + (j-1) *delta_of_L[temp_blk->RLevel][1];
		GridZ[i_j_k] = temp_blk->origin[2] + (k-1) *delta_of_L[temp_blk->RLevel][2];

	  }


	//!--------------------------------------------------
	//!--------- X-LAYER GRID NODES ---------------------
	//!--------------------------------------------------


	sprintf (filename,"%s/gnuplot/%s_xygrid",data_output_path,Run_Name);
	Outfile.open(filename);
	INT32 i = xlay;
	for (INT32 k=0; k<BlkNds_Z; k++)
	{
		for (INT32 j=0; j<BlkNds_Y; j++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << GridY[i_j_k] << "   ";
		}
		Outfile << endl;
	}
	Outfile.close();

	sprintf (filename,"%s/gnuplot/%s_xzgrid",data_output_path,Run_Name);
	Outfile.open(filename);
	for (INT32 k=0; k<BlkNds_Z; k++)
	{
		for (INT32 j=0; j<BlkNds_Y; j++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << GridZ[i_j_k] << "   ";
		}
		Outfile << endl;
	}
	Outfile.close();

	//!--------------------------------------------------
	//!--------- Y-LAYER GRID NODES ---------------------
	//!--------------------------------------------------
	sprintf (filename,"%s/gnuplot/%s_yxgrid",data_output_path, Run_Name);
	Outfile.open(filename);
	INT32 j = ylay;
	for (INT32 k=0; k<BlkNds_Z; k++)
	{
		for (INT32 i=0; i<BlkNds_X; i++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << GridX[i_j_k] << "   ";
		}
		Outfile << endl;
	}
	Outfile.close();

	sprintf (filename,"%s/gnuplot/%s_yzgrid",data_output_path, Run_Name);
	Outfile.open(filename);
	for (INT32 k=0; k<BlkNds_Z; k++)
	{
		for (INT32 i=0; i<BlkNds_X; i++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << GridZ[i_j_k] << "   ";
		}
		Outfile << endl;
	}
	Outfile.close();


	//!--------------------------------------------------
	//!--------- Z-LAYER GRID NODES ---------------------
	//!--------------------------------------------------
	sprintf (filename,"%s/gnuplot/%s_zxgrid",data_output_path, Run_Name);
	Outfile.open(filename);
	INT32 k = zlay;
	for (INT32 j=0; j<BlkNds_Y; j++)
	{
		for (INT32 i=0; i<BlkNds_X; i++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << GridX[i_j_k] << "   ";
		}
		Outfile << endl;
	}
	Outfile.close();

	sprintf (filename,"%s/gnuplot/%s_zygrid",data_output_path, Run_Name);
	Outfile.open(filename);
	for (INT32 j=0; j<BlkNds_Y; j++)
	{
		for (INT32 i=0; i<BlkNds_X; i++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << GridY[i_j_k] << "   ";
		}
		Outfile << endl;
	}
	Outfile.close();




	delete GridX;
	delete GridY;
	delete GridZ;



}



//!-------------------------------------------------------------//
//! gnuplot_xlayer: -								//
//!-------------------------------------------------------------//
void gnuplot_xlayer(D_REAL* Field, char *field_name, INT32 c, INT32 xlay, INT32 step)
{

	INT32 i_j_k;
	ofstream Outfile;

	char filename[200];
	
	sprintf (filename,"%s/gnuplot/%s_%s%d_x%d_%d",data_output_path, Run_Name,field_name,c,xlay,step);
// 	log_file << " Plotting GNUPLOT xlayer" << xlay
// 	     << " of component "<< c
// 	     <<	" of vfield "<< field_name
// 	     <<	" in file " << filename << endl;



	Outfile.open(filename);
	
	INT32 i = xlay;
	for (INT32 k=0; k<BlkNds_Z; k++)
	{
		for (INT32 j=0; j<BlkNds_Y; j++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << Field[i_j_k +c*num_nodes_in_block] << "   ";
		}
		Outfile << endl;
	}
	
	Outfile.close();
// 	log_file << "done. " << endl;

}


//!-------------------------------------------------------------//
//! gnuplot_ylayer: -								//
//!-------------------------------------------------------------//
void gnuplot_ylayer(D_REAL* Field, char *field_name, INT32 c, INT32 ylay, INT32 step)
{

	INT32 i_j_k;
	ofstream Outfile;

	char filename[100];
	
	sprintf (filename,"%s/gnuplot/%s_%s%d_y%d_%d",data_output_path, Run_Name,field_name,c,ylay,step);
// 	log_file << " Plotting GNUPLOT ylayer" <<ylay
// 	     << " of component "<< c
// 	     <<	" of vfield "<< field_name
// 	     <<	" in file " << filename << endl;



	Outfile.open(filename);
	
	INT32 j = ylay;
	for (INT32 k=0; k<BlkNds_Z; k++)
	{
		for (INT32 i=0; i<BlkNds_X; i++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << Field[i_j_k +c*num_nodes_in_block] << "   ";
		}
		Outfile << endl;
	}
	
	Outfile.close();
// 	log_file << "done. " << endl;

}

//!-------------------------------------------------------------//
//! gnuplot_xlayer: -								//
//!-------------------------------------------------------------//
void gnuplot_zlayer(D_REAL* Field, char *field_name, INT32 c, INT32 zlay, INT32 step)
{

	INT32 i_j_k;
	ofstream Outfile;

	char filename[100];
	
	sprintf (filename,"%s/gnuplot/%s_%s%d_z%d_%d",data_output_path, Run_Name,field_name,c,zlay,step);
// 	log_file << " Plotting GNUPLOT zlayer" << zlay
// 	     << " of component "<< c
// 	     <<	" of vfield "<< field_name
// 	     <<	" in file " << filename << endl;



	Outfile.open(filename);
	
	INT32 k = zlay;
	for (INT32 j=0; j<BlkNds_Y; j++)
	{
		for (INT32 i=0; i<BlkNds_X; i++)
		{
		
			i_j_k =       i*BlkNds_Y*BlkNds_Z    +j*BlkNds_Z      +k;
		
			Outfile << Field[i_j_k +c*num_nodes_in_block] << "   ";
		}
		Outfile << endl;
	}
	
	Outfile.close();
// 	log_file << "done. " << endl;

}


