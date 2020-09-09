#ifndef SAVE_H_
#define SAVE_H_

#include <iomanip>

#include "par.h"
#include "scalar_field.h"
#include "vector_field.h"
#include "electrics_calc.h"

/* ------------------------ Save and log functions ------------------------- */

void save(  long double curr_time, int& save_counter, ScalarField & Phi,
            VectorField & ElectricField, ScalarField & Sigma,
            VectorField & CurrentDensity, ScalarField & Temperature);

void log(   long double curr_time, int pulse, ScalarField & Temperature,
            ScalarField & Phi, ScalarField & Sigma, VectorField & ElectricField, VectorField & CurrentDensity,
            long double q_accum, bool & log_created );

void save_parameters();


#endif
