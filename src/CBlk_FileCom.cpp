#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctime>



#include "CBlk.h"
#include "utils.h"
#include "parameters.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "absolute_Globals.h"



void CBlock::getMesh(INT32 id_oct,ofstream& outfile)
{
    
    //! get index of oct (0 or 1 for each dimension)
    INT32 a = id_oct/4;
    INT32 b = (id_oct -a*4)/2;
    INT32 c = (id_oct -a*4 -b*2);
    
    
    //! In case Block is Box boundary, filling shall start at 0, else at 1
    INT32 start[3] = {!a*!is_box_boundary[_Im1_] +a*(BlkNds_X/2),
    !b*!is_box_boundary[_Jm1_] +b*(BlkNds_Y/2),
    !c*!is_box_boundary[_Km1_] +c*(BlkNds_Z/2)};
    
    
    INT32 end[3]={BlkNds_X-1 -!a*((BlkNds_X-2)/2),
    BlkNds_Y-1 -!b*((BlkNds_Y-2)/2),
    BlkNds_Z-1 -!c*((BlkNds_Z-2)/2)};
    
    
    INT32 i,j,k,i_j_k;
    INT32 cell_indices[3];
    
  
    PARTICLE_REAL x_part[3], r_normed[3];

    for(i = start[0]; i < end[0]; i++)
        for(j = start[1]; j < end[1]; j++)
            for(k = start[2]; k < end[2]; k++)
            {        
                i_j_k  = i*BlkNds_Y*BlkNds_Z 
                +j*BlkNds_Z 
                +k;
                
                //! Position 0 in der Zelle
                x_part[0] =  0.;
                x_part[1] =  0.;
                x_part[2] =  0.;
                    
                cell_indices[0] = i;
                cell_indices[1] = j;
                cell_indices[2] = k;
                    
                //! Berechnung der internen Position
                intern2normedCoords(r_normed, x_part, cell_indices);
                
                outfile << r_normed[0] << "\t" << r_normed[1] << "\t" <<  r_normed[2] << endl;
            }                
}
