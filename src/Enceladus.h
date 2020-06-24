



//! calculated values
const D_REAL mass_constant = 4./3.*M_PI*1.e-24;
const D_REAL charge_constant = 4.*M_PI*eps_0*1.e-09*dust_potential/e;
const D_REAL q_m_constant = -3.*eps_0*dust_potential/1000.*average_ion_mass*m_p/e*1.e+18;
 
//! calculations from values
const D_REAL H_theta_gas = 1./((opening_angle_gas*M_PI/180.)*(opening_angle_gas*M_PI/180.));
const D_REAL H_theta_dust = 1./((opening_angle_dust*M_PI/180.)*(opening_angle_dust*M_PI/180.));
const D_REAL cosphi0 = cos(phi0);
const D_REAL costheta0 = cos(theta0);
const D_REAL costheta0_ChEx = cos(M_PI - theta0);
const D_REAL sinphi0 = sin(phi0);
const D_REAL sintheta0 = sin(theta0);
const D_REAL dichte_normiert_ChEX = base_density * 1.e-06 * (average_ion_mass*m_p/(e*B_SI));



 
