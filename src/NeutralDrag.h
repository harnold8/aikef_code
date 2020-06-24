


#include "parameters.h"
// #include "utils.h"

#include <math.h>
#include <iostream>
#include <fstream>

// #define  COM_G		8.0e+25  // ions per sec. (SI units)
// #define  SI_nuni	1.7e-15  // Coll. freq. neutrals/ions in m^3/s
// #define  SI_t0		5.22045
// #define  SI_x0		192472.
// #define  COM_v0         1.e+3    // outgasing velocity

//! apply neutral drag like in TB CODE (Nikos Comet Version)
// const D_REAL ND_k =1.*COM_G*SI_nuni*SI_t0/(4.*M_PI*COM_v0*SI_x0*SI_x0); //Bagdonat;

using namespace std;
//!--------------------------------------------------------------------//
//! apply_NeutralDrag: -					       //
//!	r: position of particle in normed coords	(no pointer    //
//! 	v: pointer to velocity of particle			       //
//!	use Ion_Masses[species]  as mass of paricle		       //
//!	use Ion_Charges[species] as charge of paricle		       //
//!--------------------------------------------------------------------//
inline void apply_NeutralDrag(/*INT32 species,*/ D_REAL *r, D_REAL *v, D_REAL dt_particle)
{

	//! get squared distance to origin
	D_REAL rsqr = vec_len2(r);

	//! in case rsqr gets 0, force gets infinitly high
	if(rsqr<1.e-4)
	rsqr = 1.e-4;
	//! New Neutral by Christoph Koenders with Euler
	//! decelerate particle velocity

	//! Neutral drag is switched off see parameter.h paramers.cpp
	const D_REAL ND_k=0.;
	v[0] = v[0] -(v[0] - r[0]* COM_v0 /(sqrt(rsqr)*SI_v0) )*ND_k/rsqr*dt_particle;
	v[1] = v[1] -(v[1] - r[1]* COM_v0 /(sqrt(rsqr)*SI_v0) )*ND_k/rsqr*dt_particle;
	v[2] = v[2] -(v[2] - r[2]* COM_v0 /(sqrt(rsqr)*SI_v0) )*ND_k/rsqr*dt_particle;
    //! Old Neutral by Nikolaos Gortsas with Euler
	//! decelerate particle velocity
// 	v[0] = v[0] -(v[0])*ND_k/rsqr*dt_particle;
// 	v[1] = v[1] -(v[1])*ND_k/rsqr*dt_particle;
// 	v[2] = v[2] -(v[2])*ND_k/rsqr*dt_particle;



}
