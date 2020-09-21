#ifndef PAR_H_
#define PAR_H_

#include "mesh.h"
#include "scalar_field.h"
#include "vector_field.h"

/* ---------- Here you will find the parameters of the simulation ---------- */

namespace par
{

/* ----------------------- Declaration of functions defined at the end ----------------------------- */ 

int calc_no_pulses();
long double calc_max_time();
long double * calc_max_voltages();
long double * calc_dts_on(int);
long double * calc_dts_off(int);
int * calc_save_steps();

/* ----------------------- Universal constants ----------------------------- */

/* Tissue temperature [K] */
const long double temp = 310.15;
/* Air temperature. [K] */
const long double temp_air = 298.15;                                    

/* ---------------------------- Electrodes ---------------------------------- */

/* Electrode length [m] */
const long double electrode_length = 7e-3;
/* Electrode width [m] */
const long double electrode_width = 0.7e-3;
/* Electrode thickness [m] */
const long double electrode_thickness = 0.7e-3;
/* Space between anode-cathode, horizontally [m] */
const long double gap_anode_cathode = 8e-3;
/* Type of electrode: "plates" or "needles" */
const std::string electrode_type = "needles";
/* Space between anode-anode and cathode-cathode. [m] */
/* Use only with needle electrodes */
const long double gap_elect_elect = 5e-3;
/* Use only with needle electrodes. No. of electrodes. */
const int no_electrodes = 1;

/* ------------------------- Domain parameters ----------------------------- */

/* Maximum length of the x axis [m] */
const long double x_max = 15e-3;
/* Maximum length of the y axis [m] */
const long double y_max = 15e-3;
/* Maximum length of the z axis [m] */
const long double z_max = 10e-3;

/* ----------------------- Treatment parameters ---------------------------- */

/* Voltages to distance ratio [V m^-1] */
const double volt_to_dist[] = {25000};
/* On times [s] */
const double on_pulse_times[] = {0.05};
/* Off times [s] */
const double off_pulse_times[] = {0.950};
/* How many repetitions of each pulse */
const int pulse_repetitions[] = {8};
/* No. of cycle repetition. A cycle is a sequence of pulses. */
const int no_cycles = 1;



/* ----------- Other parameters, calculated from the above ----------------- */

/* No. of pulses */
const int no_pulses = calc_no_pulses();
/* No. of elements of each cycle. Each element could correspond to several pulses */
const int no_elems_per_cycle = sizeof(par::volt_to_dist) / sizeof(par::volt_to_dist[0]);
/* Simulation total time [s] */
const long double max_time = calc_max_time();
/* Maximum voltages [V] */
inline const long double * max_voltages = calc_max_voltages();
/* Potentials in the anode [V] */
inline const long double * phi_anodes = max_voltages;
/* Potential in the cathode [V] */
const long double phi_cathode = 0.;
/* Initial electric potential [V] */
const long double phi_init = 0.;
/* Initial electric field [V m^-1] */
const long double ef_init = 0.;

/* -------------------------- Tissue paramters ------------------------------ */

/* Initial temperature [K] */
const long double temperature_init = temp;
/* Thermal conductivity [ W m^-1 K^-1 ] */
const long double k_domain_init = 0.512;
/* Electrical conductivity [ S m^-1 ] */
const long double sigma_domain_init = 0.504;
/* Tissue density [ kg m^3 ] */
const long double rho = 1050.;
/* Heat capacity [J kg^-1 K^-1] */
const long double cp = 3600.;
/* Metabolic heat generation [ W m^-3 ] */
const long double qm = 420.;                                            

/* --------------------------- Blood parameters ------------------------------ */

/* Temperature of the arterial blood [K] */
const long double temperature_b = 310.15;
/* Density of the blood [ kg m^3 ] */
const long double rho_b = 1060;
/* Heat capacity of the blood [J kg^-1 K^-1] */
const long double c_b = 3600;
/* Blood profusion [ s^-1 ] */
const long double w_b = 0.0044;

/* -------------------------- Electrode parameters --------------------------- */

/* Thermal conductivity [ W m^-1 K^-1 ] */
const long double k_electrode_init = 16.3;
/* Electrical conductivity [ S m^-1 ] */
const long double sigma_electrode_init = 1.398e6;
/* Electrode density [ kg m^-3 ] */
const long double rho_electrode = 1050.;
/* Heat capacity [J kg^-1 K^-1] */                             
const long double cp_electrode = 490.;
                              
/* --------------------------- Boundary parameters ---------------------------- */

/* Convective heat-transfer coefficient [ W m^-2 K^-1 ]  */
const long double h = 10.;
/* Electrical conductivity increment with regard to temperature */
const long double alpha0 = 0.015;

/* ----------------------- Geometry and Mesh Parameters ----------------------- */

/* No. of nodes between anodes and cathodes */
const int resolution = 34;                                              
/* Mesh */
const Mesh mesh(    x_max, y_max, z_max, electrode_length, electrode_width,
                    electrode_thickness, gap_anode_cathode, gap_elect_elect,
                    no_electrodes, resolution, electrode_type );
/* Anode index positions */
int ** const anodes = mesh.get_anodes();
/* Cathode index positions */
int ** const cathodes = mesh.get_cathodes();
/* No. of x axis divisions */
const int ii = mesh.get_ii();
/* No. of y axis divisions */
const int jj = mesh.get_jj();
/* No. of z axis divisions */
const int kk = mesh.get_kk();
/* Space between nodes in x axes */
const long double dx = mesh.get_dx();
/* Space between nodes in y axes */
const long double dy = mesh.get_dy();
/* Space between nodes in z axes */
const long double dz = mesh.get_dz();

/* ----------------------- Other numerical parameters ------------------------- */

/* Convergence tolerance */
const double tol = 1.e-9;
/* Max. sub iteration number for Phi */
const int max_sub_it_Phi = 50;
/* Max. sub iteration number for intial condition of Phi */
const int max_sub_it_init_Phi = 1000;
/* Relaxation */
const long double omega =  1.;

/* ----------------------------- Pulse parameters ----------------------------- */

/* Time step during on pulse [s] */
const long double dt_fractions_on = 5000;
inline const long double * dts_on_pulse =  calc_dts_on(dt_fractions_on);
/* Time step during off pulse [s] */
const long double dt_fractions_off = 5000;
inline const long double * dts_off_pulse = calc_dts_off(dt_fractions_off);

/* -------------------------- Save and log parameters  ------------------------ */

/* Save each ... iterations */
inline const int * save_steps = calc_save_steps();
/* Log each ... iterations */
const int log_step = 50;
/* "vtk" or "csv" */
const std::string save_format = "vtk";

/* -------------------------- Calculation of some parameters  ------------------------ */

inline int calc_no_pulses()
{
    int n_repetitions = sizeof(pulse_repetitions) / sizeof(pulse_repetitions[0]);
    int sum_repetitions = 0;
    for (int i = 0; i < n_repetitions; i++) 
    {
        sum_repetitions += pulse_repetitions[i]; 
    }
    return sum_repetitions * no_cycles;
}

inline long double calc_max_time()
{
    double parcial_time = 0;
    for (int i = 0; i < no_elems_per_cycle; i++)
    {
        parcial_time += (on_pulse_times[i] + off_pulse_times[i]) * pulse_repetitions[i];
    }

    return  parcial_time * no_cycles;
}

inline long double * calc_max_voltages()
{
    long double * max_voltages = new long double[no_elems_per_cycle];
    for (int i = 0; i < no_elems_per_cycle; i++)
    {
        max_voltages[i] = volt_to_dist[i] * gap_anode_cathode;
    }   
    return (long double *) max_voltages;
}

inline long double * calc_dts_on(int dt_fractions)
{
    long double * dts_on = new long double[no_elems_per_cycle];
    for (int i = 0; i < no_elems_per_cycle; i++)
    {
        dts_on[i] = on_pulse_times[i] / dt_fractions;
    }   
    return (long double *) dts_on;

}

inline long double * calc_dts_off(int dt_fractions)
{
    long double * dts_off = new long double[no_elems_per_cycle];
    for (int i = 0; i < no_elems_per_cycle; i++)
    {
        dts_off[i] = off_pulse_times[i] / dt_fractions;
    }   
    return (long double *) dts_off;

}

inline int * calc_save_steps()
{
    int * save_steps = new int[no_elems_per_cycle]; 
    for (int i = 0; i < no_elems_per_cycle; i++)
    {
        save_steps[i] = 1. / dts_on_pulse[i];
    }   
    return (int *) save_steps;
}



}

#endif
