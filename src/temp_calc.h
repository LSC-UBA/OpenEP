#ifndef TEMP_H_
#define TEMP_H_

#include "par.h"

/* ------------------------- Temperature calculation ----------------------- */

void calc_temperature(  ScalarField & Temperature, ScalarField & Temperature_old,
                        VectorField & ElectricField, ScalarField & K, 
                        ScalarField & Sigma, long double curr_time, int it_number,
                        long double dt);

#endif
