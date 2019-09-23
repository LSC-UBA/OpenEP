#ifndef ELEC_H_
#define ELEC_H_

#include "par.h"

/* --------------------- Electrical variable calculation ------------------- */

/* Calculate initial distribution of Phi */
void init_phi(ScalarField & Phi);

/* Calculate Phi: div( sigma grad(Phi) ) = 0 */
void calc_nonlinear_phi( ScalarField & Phi, ScalarField & Phi_aux,
                         ScalarField & Sigma, long double curr_time,
                         int it_number, int max_sub_it = par::max_sub_it_Phi);

/* Calculate electric field: E = - grad Phi */
void calc_electric_field(VectorField & ElectricField, ScalarField & Phi);

/* Calculate electric conductivity: sigma = sigma0 * (1 + alpha (T − temp) */
void calc_sigma( ScalarField & Sigma, ScalarField & Temperature,
                 long double curr_time, int it_number );

/* Calculate electric current: integral_surface J * ds */
long double calc_electric_current(VectorField & CurrentDensity);

/* Set boundary condition at the electrodes */
void set_in_electrodes(ScalarField & scalarField, const long double val,
                       int ** electrodes, int electrodes_number);


#endif
