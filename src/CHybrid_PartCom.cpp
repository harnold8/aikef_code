/***************************************************************************
 *   Copyright (C) 2009 by Joachim Mueller   *
 *   joachim@champ   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#include "CHybrid.h"


#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>





//!--------------------------------------------------------
//!- set_Blk_optimal_MPiC:
//!--------------------------------------------------------
void CHybrid::set_Blk_optimal_MPiC(void)
{

	log_file << "  Setting oMPiC ...  " << endl;
	for(INT32 level=MAX_LEVEL; level>=0; level--)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			temp_Block->set_Blk_optimal_MPiC();
	
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
	log_file << "  done." << endl;

}



//!-----------------------------------------------------------
//! reassign_particle_iDirection
//!-----------------------------------------------------------
void CHybrid::reassign_particle_iDirection(void)
{



	//! define X=x_Min layer
	INT32 SENTm1_start_ijk[3] = {         0,        0,        0};
	INT32 SENTm1_end_ijk[3] =   {         1, BlkNds_Y, BlkNds_Z};

	//! define X=x_Max layer
	INT32 SENTp1_start_ijk[3] = {BlkNds_X-1,        0,        0};
	INT32 SENTp1_end_ijk[3] =   {BlkNds_X-0, BlkNds_Y, BlkNds_Z};

	//! define X=x_Min+1 layer
	INT32 RECVm1_start_ijk[3] = {         1,        0,        0};
	INT32 RECVm1_end_ijk[3] =   {         2, BlkNds_Y, BlkNds_Z};

	//! define X=x_Max-1 layer
	INT32 RECVp1_start_ijk[3] = {BlkNds_X-2,        0,        0};
	INT32 RECVp1_end_ijk[3] =   {BlkNds_X-1, BlkNds_Y, BlkNds_Z};

	INT32 entries_in_MPiC_Info_package = BlkNds_Y *BlkNds_Z *num_Charged_Species +1;

	reassign_particle_equal_level( _Im1_, _Ip1_,
				      SENTm1_start_ijk,
				      SENTm1_end_ijk,
				      SENTp1_start_ijk,
				      SENTp1_end_ijk,
				      RECVm1_start_ijk,
				      RECVm1_end_ijk,
				      RECVp1_start_ijk,
				      RECVp1_end_ijk,
				      entries_in_MPiC_Info_package);



}

//!-----------------------------------------------------------
//! reassign_particle_jDirection
//!-----------------------------------------------------------
void CHybrid::reassign_particle_jDirection(void)
{



	//! define Y=y_Min layer
	INT32 SENTm1_start_ijk[3] = {       0,          0,        0};
	INT32 SENTm1_end_ijk[3] =   {BlkNds_X,          1, BlkNds_Z};

	//! define Y=y_Max layer
	INT32 SENTp1_start_ijk[3] = {       0, BlkNds_Y-1,        0};
	INT32 SENTp1_end_ijk[3] =   {BlkNds_X, BlkNds_Y-0, BlkNds_Z};

	//! define Y=y_Min+1 layer
	INT32 RECVm1_start_ijk[3] = {       0,          1,        0};
	INT32 RECVm1_end_ijk[3] =   {BlkNds_X,          2, BlkNds_Z};

	//! define Y=y_Max-1 layer
	INT32 RECVp1_start_ijk[3] = {       0, BlkNds_Y-2,        0};
	INT32 RECVp1_end_ijk[3] =   {BlkNds_X, BlkNds_Y-1, BlkNds_Z};

	INT32 entries_in_MPiC_Info_package = BlkNds_X *BlkNds_Z *num_Charged_Species +1;

	reassign_particle_equal_level( _Jm1_, _Jp1_,
				      SENTm1_start_ijk,
				      SENTm1_end_ijk,
				      SENTp1_start_ijk,
				      SENTp1_end_ijk,
				      RECVm1_start_ijk,
				      RECVm1_end_ijk,
				      RECVp1_start_ijk,
				      RECVp1_end_ijk,
				      entries_in_MPiC_Info_package);



}

//!-----------------------------------------------------------
//! reassign_particle_kDirection
//!-----------------------------------------------------------
void CHybrid::reassign_particle_kDirection(void)
{


	//! define Z=z_Min layer
	INT32 SENTm1_start_ijk[3] = {       0,        0,          0};
	INT32 SENTm1_end_ijk[3] =   {BlkNds_X, BlkNds_Y,          1};

	//! define Z=z_Max layer
	INT32 SENTp1_start_ijk[3] = {       0,        0, BlkNds_Z-1};
	INT32 SENTp1_end_ijk[3] =   {BlkNds_X, BlkNds_Y, BlkNds_Z-0};

	//! define Z=z_Min+1 layer
	INT32 RECVm1_start_ijk[3] = {       0,        0,          1};
	INT32 RECVm1_end_ijk[3] =   {BlkNds_X, BlkNds_Y,          2};

	//! define Z=z_Max-1 layer
	INT32 RECVp1_start_ijk[3] = {       0,        0, BlkNds_Z-2};
	INT32 RECVp1_end_ijk[3] =   {BlkNds_X, BlkNds_Y, BlkNds_Z-1};

	INT32 entries_in_MPiC_Info_package = BlkNds_X *BlkNds_Y *num_Charged_Species +1;

	reassign_particle_equal_level( _Km1_, _Kp1_,
				      SENTm1_start_ijk,
				      SENTm1_end_ijk,
				      SENTp1_start_ijk,
				      SENTp1_end_ijk,
				      RECVm1_start_ijk,
				      RECVm1_end_ijk,
				      RECVp1_start_ijk,
				      RECVp1_end_ijk,
				      entries_in_MPiC_Info_package);

}



//!-----------------------------------------------------------
//! reassign_particle_iDirection
//!-----------------------------------------------------------
void CHybrid::reassign_particle_equal_level(INT32 _m1NB_, INT32 _p1NB_,
					    INT32* SENTm1_start_ijk,
					    INT32* SENTm1_end_ijk,
					    INT32* SENTp1_start_ijk,
					    INT32* SENTp1_end_ijk,
					    INT32* RECVm1_start_ijk,
					    INT32* RECVm1_end_ijk,
					    INT32* RECVp1_start_ijk,
					    INT32* RECVp1_end_ijk,
					    INT32  entries_in_MPiC_Info_package)
{




	//!-------------------------- (1) ----------------------------------------
	//!--------------- SEND INFO & PARTICLE PACKAGE --------------------------
	//!-----------------------------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock*  temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			//! send in case my_rank is responsible
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			{

				temp_Block->prepare_pPackage_toNeighbour(SENTm1_start_ijk,
									 SENTm1_end_ijk,
									 SENTp1_start_ijk,
									 SENTp1_end_ijk,
									 _m1NB_,_p1NB_);

				
				temp_Block->send_MPiC_Info_toNeighbour(entries_in_MPiC_Info_package,_m1NB_, _p1NB_);
				temp_Block->send_particle_toNeighbour(entries_in_MPiC_Info_package, _m1NB_, _p1NB_);
			}


			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}







	//!-------------------------- (2) ----------------------------------------
	//!--------------- RECEIVE INFO & PARTICLE PACKAGE -----------------------
	//!-----------------------------------------------------------------------
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank==temp_Block->responsible_mpi_process)
			{
				temp_Block->recv_MPiC_Info_fromNeighbour(entries_in_MPiC_Info_package,_m1NB_, _p1NB_);
				temp_Block->recv_particle_fromNeighbour(entries_in_MPiC_Info_package, _m1NB_, _p1NB_);


				temp_Block->insert_pPackage_fromNeighbour(RECVm1_start_ijk,
									  RECVm1_end_ijk,
									  RECVp1_start_ijk,
									  RECVp1_end_ijk,
									  _m1NB_, _p1NB_);

			}

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}




	//!-------------------------- (3) ----------------------------------------
	//!--------------- RELEASE MPI MEMORY -----------------------------------
	//!----------------------------------------------------------------------
	if(sync_mpi_send_rec) {
	  for(INT32 level=0; level<=MAX_LEVEL; level++)
	  {
	    CBlock* temp_Block = BlockList_of_Lev[level];
	    while(temp_Block)
	    {

	      if(mpi_myRank==temp_Block->responsible_mpi_process)
	      {
		temp_Block->wait_for_completedSendParticle_MPI(_m1NB_,_p1NB_);

	      }

	      temp_Block = temp_Block->next_Blk_of_BlockList;
	    }
	  }
	}
	else {
	  for(INT32 level=0; level<=MAX_LEVEL; level++)
	    wait_all_MPI(level);
	}


}





//!----------------------------------------------------------------------------------------
//!- reassign_particle_child_to_parent_LB:
//!----------------------------------------------------------------------------------------
void CHybrid::reassign_particle_child_to_parent_LB(void)
{


	//!-------------------------- (1) ----------------------------------------
	//!--------------- SEND INFO & PARTICLE PACKAGE TO PARENT ----------------
	//!-----------------------------------------------------------------------
	for(INT32 level=1; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->send_particle_toParent_LB();

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}


	//!-------------------------- (2) ----------------------------------------
	//!--------------- RECEIVE INFO & PARTICLE PACKAGE FROM CHILDREN ---------
	//!-----------------------------------------------------------------------
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			
			if(mpi_myRank==temp_Block->responsible_mpi_process)
			temp_Block->recv_particle_fromChildren_LB();



			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}




	//!-------------------------- (3) ----------------------------------------
	//!--------------- RELEASE MPI MEMORY -----------------------------------
	//!----------------------------------------------------------------------
	if(sync_mpi_send_rec) {
	  for(INT32 level=1; level<=MAX_LEVEL; level++)
	  {
	    CBlock* temp_Block = BlockList_of_Lev[level];
	    while(temp_Block)
	    {

	      if(mpi_myRank == temp_Block->responsible_mpi_process)
	      temp_Block->wait_send_particle_toParent_LB();

	      temp_Block = temp_Block->next_Blk_of_BlockList;
	    }
	}
	}
	else {
	  for(INT32 level=0; level<=MAX_LEVEL; level++)
	    wait_all_MPI(level);
	}


}





//!--------------------------------------------------------
//!- reassign_particle_parent_to_child_LB:
//!--------------------------------------------------------
void CHybrid::reassign_particle_parent_to_child_LB(void)
{


	//!-------------------------- (1) ----------------------------------------
	//!--------------- SEND INFO & PARTICLE PACKAGE TO CHILDREN --------------
	//!-----------------------------------------------------------------------
	for(INT32 level=0; level<MAX_LEVEL; level++)
	{
		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->send_particle_toChildren_LB();

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}


	//!-------------------------- (2) ----------------------------------------
	//!--------------- RECEIVE INFO & PARTICLE PACKAGE FROM PARENT -----------
	//!-----------------------------------------------------------------------
	for(INT32 level=1; level<=MAX_LEVEL; level++)
	{
		CBlock *temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{

			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->recv_particle_fromParent_LB();

			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}




	//!-------------------------- (3) ----------------------------------------
	//!--------------- RELEASE MPI MEMORY -----------------------------------
	//!----------------------------------------------------------------------
	if(sync_mpi_send_rec) {
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	  CBlock* temp_Block = BlockList_of_Lev[level];
	  while(temp_Block)
	  {
	    if(mpi_myRank == temp_Block->responsible_mpi_process)
	    temp_Block->wait_send_particle_toChildren_LB();

	    temp_Block = temp_Block->next_Blk_of_BlockList;
	  }
	}
	}
	else {
	  for(INT32 level=0; level<=MAX_LEVEL; level++)
	    wait_all_MPI(level);
	}



}





//!--------------------------------------------------------
//!- reassign_particle:
//!--------------------------------------------------------
void CHybrid::reassign_particle(D_REAL &reassign_time)
{



	clock_t start,finish;
	//! No start assigning particle of Ghost Cells to new Blocks
	start = clock();


	//! children->parent
	reassign_particle_child_to_parent_LB();


	//! neighbour->neighbour
	reassign_particle_iDirection();
	reassign_particle_jDirection();
	reassign_particle_kDirection();


	//! parent->children
	reassign_particle_parent_to_child_LB();




	finish = clock();
	reassign_time = (double(finish)-double(start))/CLOCKS_PER_SEC;


}
