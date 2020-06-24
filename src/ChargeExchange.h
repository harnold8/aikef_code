


#include "utils.h"
#include "parameters.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include <stdlib.h>
#include <math.h>
#include <ctime>

//! Random Generators
using namespace std;
extern gsl_rng **randGen_of_species;


//!--------------------------------------------------------------------//
//! calculate_charge_exchange calculates computes a probability for
//! charge exchange to occurs and returns "true" or "false".
//!--------------------------------------------------------------------//


//!-------------------------------------------------------------------------------//
//! calculate_charge_exchange:							 		//
//! calcute_charge_exchange returns true, if charge exchange should take place	//
//! 	x: position of particle in normed coords	(no pointer)		 		//
//!------------------------------------------------------------------------------//
inline bool calculate_charge_exchange(PARTICLE_REAL *x, PARTICLE_REAL *v, D_REAL dt_particle)
{

			
	D_REAL dichte, theta, absx, absv;
	//! Obstacle Radius neccessary for some density profiles, even if no obstacle exists
	const D_REAL R_Obstacle2 = 2.216;	
	
	//! constant parameters for neutral profile
	const D_REAL theta0 = 0.;
	const D_REAL H_theta = 10.;
	const D_REAL n0 = 1.418e+10;
		
	//! sigma = cross section, see explanation file for details
	const D_REAL sigma = 1.8772e-06;

	//! beta = survival probability against charge exchange 
	D_REAL beta=1.;

	//! below some checks, in order to save computational time (exp is expensive!),
	//! depending on the used neutral density profile
	if(x[2]<0)
	{ 
		absx = vec_len(x);
		//! in obstacle?
		if(absx>R_Obstacle2)	
		{
			theta = atan2( sqrt( x[0]*x[0]+x[1]*x[1]),x[2] );
			//! further selection 
			if(theta>0.9*M_PI && theta<M_PI)
			{
				//! "dichte" is the neutral density depending on x[]
				dichte = n0*(R_Obstacle2*R_Obstacle2/(absx*absx))*exp(-(theta-M_PI)*(theta-M_PI)*H_theta)*exp(-(absx-R_Obstacle2)/14.0);
				//! one more selections, still to determine parameter
				if(dichte>100.)
				{

					absv = vec_len(v);
					beta = exp(- dichte*absv*sigma*dt_particle);
						
					//! here is determined, if charge exchange occurs
					if(beta < 1.*random()/RAND_MAX)
					{
						//! set velocity to zero
						memset(v, 0, 3*sizeof(PARTICLE_REAL));
						return true;
					}

				}
			}
		}
	}

	//! in case any if interception id false, return false
	return false;
}

