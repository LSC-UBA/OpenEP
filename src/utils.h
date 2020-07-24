#ifndef UTILS_H_
#define UTILS_H_

#include "par.h"

/* Validate parameters in par.h */
void validate_par(){
	int no_volt = sizeof(par::volt_to_dist) / sizeof(par::volt_to_dist[0]);
	int no_on = sizeof(par::on_pulse_times) / sizeof(par::on_pulse_times[0]);
	int no_off = sizeof(par::off_pulse_times) / sizeof(par::off_pulse_times[0]);
	int no_pulse = sizeof(par::pulse_repetitions) / sizeof(par::pulse_repetitions[0]);

	bool treatment_par_ok;
	treatment_par_ok = (no_volt == no_on);
	treatment_par_ok &= (no_on == no_off);
	treatment_par_ok &= (no_off == no_pulse);
	treatment_par_ok &= (no_pulse == par::no_elems_per_cycle);

	if (!treatment_par_ok)
		throw "Error. Treatment parameters are inconsistent.";
}

#endif