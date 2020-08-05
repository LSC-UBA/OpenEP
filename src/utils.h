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

	if (!treatment_par_ok)
		throw "Error. Treatment parameters are inconsistent.";
}

// int calc_no_pulses()
// {
// 	int n_repetitions = sizeof(par::pulse_repetitions) / sizeof(par::pulse_repetitions[0]);
// 	int sum_repetitions = 0;
// 	for (int i = 0; i < n_repetitions; i++) 
// 	{
// 		sum_repetitions += par::pulse_repetitions[i]; 
// 	}
// 	return sum_repetitions * par::no_cycles;
// }

// long double calc_max_time()
// {
// 	double parcial_time = 0;
// 	for (int i = 0; i < par::no_elems_per_cycle; i++)
// 	{
// 		parcial_time += (par::on_pulse_times[i] + par::off_pulse_times[i]) * par::pulse_repetitions[i];
// 	}

// 	return  parcial_time * par::no_cycles;
// }

// long double * calc_max_voltages()
// {
// 	long double * max_voltages = new long double[par::no_elems_per_cycle];
// 	for (int i = 0; i < par::no_elems_per_cycle; i++)
// 	{
// 		max_voltages[i] = par::volt_to_dist[i] * gap_anode_cathode;
// 	}	
// 	return (long double *) max_voltages;
// }

// long double * calc_dts_on(int dt_fractions)
// {
// 	long double * dts_on = new long double[par::no_elems_per_cycle];
// 	for (int i = 0; i < par::no_elems_per_cycle; i++)
// 	{
// 		dts_on[i] = par::on_pulse_times[i] / dt_fractions;
// 	}	
// 	return (long double *) dts_on;

// }

// long double * calc_dts_off(int dt_fractions)
// {
// 	long double * dts_off = new long double[par::no_elems_per_cycle];
// 	for (int i = 0; i < par::no_elems_per_cycle; i++)
// 	{
// 		dts_off[i] = par::off_pulse_times[i] / dt_fractions;
// 	}	
// 	return (long double *) dts_off;

// }

// int * calc_save_steps()
// {
// 	int * save_steps = new int[par::no_elems_per_cycle];	
// 	for (int i = 0; i < par::no_elems_per_cycle; i++)
// 	{
// 		save_steps[i] = 1. / par::dts_on_pulse[i];
// 	}	
// 	return (int *) save_steps;
// }

#endif