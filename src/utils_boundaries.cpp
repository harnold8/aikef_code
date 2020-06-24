

#include <math.h>
#include <iostream>

#include "utils.h"
#include "CHybrid.h"
#include "parameters.h"
#include "hilbert.h"
#include "absolute_Globals.h"

extern ofstream log_file;

const INT32   TL_start = 12000;
const D_REAL barrier_start = -0.3*LY;
const D_REAL barrier_speed = V_sw[1];

//     const D_REAL a = ;
//     const D_REAL b = ;
//     const D_REAL c = ;
//     const D_REAL d = ;
//     const D_REAL e = ;

D_REAL velocity(D_REAL r)
{   
    const D_REAL a = 0.697276;
    const D_REAL b = 3.0442e-9;
    const D_REAL c = -4.50639e-17;
    const D_REAL d = 5.3092e-25;
    const D_REAL e = -3.00721e-33;
    return a + b*r + c*pow(r,2) + d*pow(r,3) + e*pow(r,4);
    //     
}


D_REAL densityCometaryIons(D_REAL r)
{   
    const D_REAL a = 0.021562;
    const D_REAL b = -3.4494e-10;
    const D_REAL c = 6.43525e-18;
    const D_REAL d = -8.49322e-26;
    const D_REAL e = 5.07107e-34;
    return a + b*r + c*pow(r,2) + d*pow(r,3) + e*pow(r,4);
    //      
}

D_REAL densitySolarWindIons(D_REAL r)
{   
    const D_REAL a = 1.43415;
    const D_REAL b = -6.21666e-9;
    const D_REAL c = 1.11861e-16;
    const D_REAL d = -1.45827e-24;
    const D_REAL e = 8.67022e-33;
    return a + b*r + c*pow(r,2) + d*pow(r,3) + e*pow(r,4);
    //    
}



//!-------------------------------------------------------------//
//! set_inflow_BField:
//!-------------------------------------------------------------//
void set_inflow_BField(D_REAL *B, PARTICLE_REAL *x)
{

// 	//! use field 0
// 	INT32 field=0;
// 
// 
// 	D_REAL* coord = values_EBFs[field] +0 *num_rows_EBFs[field];
// 	D_REAL*    BX = values_EBFs[field] +1 *num_rows_EBFs[field];
// 	D_REAL*    BY = values_EBFs[field] +2 *num_rows_EBFs[field];
// 	D_REAL*    BZ = values_EBFs[field] +3 *num_rows_EBFs[field];
// 
// 	B[0] = B_sw[0];
// 	B[1] = B_sw[1];
// 	B[2] = B_sw[2] * (2./(1.+exp(0.4 *(x[1]-TL*dt*barrier_speed -barrier_start) ))-1.);
// 
// 	D_REAL barrier_shift = (TL-TL_start)*dt*barrier_speed +barrier_start;
// 
// 	const D_REAL a = 0.368416;
// 	const D_REAL b = 27.687-30.;
// 	B[2] = 2./(1.+exp(a*(x[1]-b -barrier_shift))) -1.;
        //! Quadratischer Abstand zur x Achse 
        D_REAL rad = (x[1]*x[1]+x[2]*x[2])*SI_x0*SI_x0/(1000*1000);
        //! Änderung der Geschwindigkeit
        D_REAL vel = velocity(rad);
//         log_file << x[0] << " " << x[1] << " "<< x[2] << " - " << rad << " vel " << vel << "\n";
        
        B[0] = B_sw[0]/vel;
        B[1] = B_sw[1]/vel;
        B[2] = B_sw[2];
        // B_sw[2] bleibt unverändert.
        
}


//!-------------------------------------------------------------//
//! set_inflow_velocity:
//!-------------------------------------------------------------//
void set_inflow_velocity(PARTICLE_REAL *v, PARTICLE_REAL *x, INT32 species)
{

/*
	D_REAL barrier_shift = (TL-TL_start)*dt*barrier_speed +barrier_start;

	if(x[1]-barrier_shift<0)
	{
		const D_REAL a1=-478515.;
		const D_REAL b1=395175.;
		const D_REAL c1=-0.0693621;
		const D_REAL d1=0.54794;
		v[0] = a1/(b1+exp(d1*(x[1]+29. -barrier_shift))) -c1;

	}
	else
	{
		const D_REAL a2=-1.12428;
		const D_REAL b2=2.03927;
		const D_REAL c2=33.0949 -29.;
		const D_REAL d2=1.11713 -0.02;
		v[0] = a2/(1.+b2*exp(x[1]-c2 -barrier_shift)) +d2;

	}



// 	v[0] = -V_sw[0] * (2./(1.+exp(0.3 *(x[1]-TL*dt*barrier_speed -barrier_start) ))-1.);
//	v[1] =  V_sw[1];
//	v[2] =  V_sw[2];*/

        //! Quadratischer Abstand zur x Achse 
        D_REAL rad = (x[1]*x[1]+x[2]*x[2])*SI_x0*SI_x0/(1000*1000);
        //! Änderung der Geschwindigkeit
        D_REAL vel = velocity(rad);
        
        v[0] = V_sw[0]*vel;
        v[1] = V_sw[1]*vel;
        v[2] = V_sw[2]*vel;

}

//!-------------------------------------------------------------//
//! set_inflow_density:
//!-------------------------------------------------------------//
void set_inflow_density(PARTICLE_REAL &rho, PARTICLE_REAL *x, INT32 species)
{       
        //! Quadratischer Abstand zur x Achse 
        D_REAL rad = (x[1]*x[1]+x[2]*x[2])*SI_x0*SI_x0/(1000*1000);
        if(species==0) 
            rho = densitySolarWindIons(rad);
        else 
        {
            rho = densityCometaryIons(rad)*COM_Q[species-1]/(COM_Q_Total);
        }

}
