
#include "defines.h"



//! forward declaration
class vfield;

//!-----------------------------------------------------------
//!SBox_Info: Structure to store general Box information 
//!		 - If variables are added, the "SBox_Info Stucure
//!		   of the visalization has to be syncronized
//!		 - All quantities are stored in normalized units
//!-----------------------------------------------------------


struct SGrid_Info
{
   
        //! Simulation Name
    char Run_Name[42];   // BUG FIX: Changed from Run_Name[40].


    //! number of mesh points
    short NP[3];

    //! time level of 3D Data
    int TL;


    //! Box Origin
    FILE_REAL Origin[3];

    //! Length of Box
    FILE_REAL Length[3];

    //! Raidius of Obstacle 
    FILE_REAL R_Obstacle;


    //! SI Units for re-Normalization:
    //! inertia lenth
    FILE_REAL x0;
    //! Alfven velocity
    FILE_REAL v0;
    //! IMF
    FILE_REAL B0;
    //! Backround plasma density
    FILE_REAL n0;


};


void uniform_grid_output(int TL);
void uniform_grid_writeField(INT32 Field_ID, INT32 uniform_grid_TYPE);
void uniform_grid_write(INT32 Field_ID, INT32 uniform_grid_TYPE);
void init_GridInfo_Structure(SGrid_Info* Grid_Info, int TL);
void uniform_grid_parallel_output_to_single_file(INT32 Field_ID, INT32 uniform_grid_TYPE, INT32 num_grid_nds);

