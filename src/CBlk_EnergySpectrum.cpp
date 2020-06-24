#include "CBlk.h"
#include "utils.h"
#include <math.h>
#include "parameters.h"
#include "absolute_Globals.h"

using namespace std;

extern D_REAL *CellVol_of_L;


void CBlock::get_velocity(D_REAL** &velocity, INT32 species,  INT32 i_j_k, D_REAL* &gewicht)
{
	

	particle* active_particle;

	//! position of particle
	PARTICLE_REAL x[3];

	for(INT32 part_index=0; part_index<num_MPiC[species][i_j_k];part_index++)
	{

		//! set pointer to particle at position part_index
		active_particle = pArray[species][i_j_k] +part_index;
		
		for(INT32 comp=0; comp<3; comp++)
		velocity[comp][part_index] =  active_particle->v[comp];
		
		gewicht[part_index] = active_particle->weight/CellVol_of_L[RLevel];
		
#ifdef DUST_PARTICLE_RADIUS		
		if(Ion_Charges[species]<0)
		gewicht[part_index] *= 1./(1.*pow(active_particle->radius,mu_slope));
#endif		
		
	}

}