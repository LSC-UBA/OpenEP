#ifndef PAR_H_
#define PAR_H_

#include "mesh.h"
#include "scalar_field.h"
#include "vector_field.h"

/* ---------- Here you will find the parameters of the simulation ---------- */

namespace par
{

/* ----------------------- Universal constants ----------------------------- */

/* Tissue temperature [K] */
const long double temp = 310.15;
/* Air temperature. [K] */
const long double temp_air = 298.15;                                    

/* ---------------------------- Electrodes ---------------------------------- */

/* Electrode length [m] */
const long double electrode_length = 19e-3;
/* Electrode width [m] */
const long double electrode_width = 7e-3;
/* Electrode thickness [m] */
const long double electrode_thickness = 1e-3;
/* Space between anode-cathode, horizontally [m] */
const long double gap_anode_cathode = 8e-3;
/* Type of electrode: "plates" or "needles" */
const std::string electrode_type = "plates";
/* Space between anode-anode and cathode-cathode. [m] */
/* Use only with needle electrodes */
const long double gap_elect_elect = 5e-3;
/* Use only with needle electrodes. No. of electrodes. */
const int no_electrodes = 1;

/* ------------------------- Domain parameters ----------------------------- */

/* Maximum length of the x axis [m] */
const long double x_max = 2e-2;
/* Maximum length of the y axis [m] */
const long double y_max = 1.5e-2;
/* Maximum length of the z axis [m] */
const long double z_max = 2.5e-2;

/* ----------------------- Treatment parameters ---------------------------- */

/* Voltage to distance ratio [V m^-1] */
const long double volt_to_dist = 25000.;
/* Frequency [Hz] */
const long double freq = 1.;
/* On time [s] */
const long double on_pulse_time = 0.05;
/* No. of pulses */
const long double nbr_pulses = 8;

/* ----------- Other parameters, calculated from the above ----------------- */

/* Total pulse time [s] */
const long double total_pulse_time = 1. / freq;
/* Off pulse time [s] */
const long double off_pulse_time = total_pulse_time - on_pulse_time;
/* Simulation total time [s] */
const long double max_time = nbr_pulses * total_pulse_time;
/* Maximum voltage [V] */
const long double max_voltage = volt_to_dist * gap_anode_cathode;
/* Potential in the anode [V] */
const long double phi_anode = max_voltage;
/* Potential in the cathode [V] */
const long double phi_cathode =  0.;
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
const int resolution = 30;                                              
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
const int max_sub_it_init_Phi = 10000;
/* Relaxation */
const long double omega =  1.;

/* ----------------------------- Pulse parameters ----------------------------- */

/* Time step during on pulse [s] */
const long double dt_on_pulse =  on_pulse_time / 500.;
/* Time step during off pulse [s] */
const long double dt_off_pulse = on_pulse_time / 500.;

/* -------------------------- Save and log parameters  ------------------------ */

/* Save each ... iterations */
const int save_step = 1. / dt_on_pulse;
/* Log each ... iterations */
const int log_step = 500;
/* "vtk" or "csv" */
const std::string save_format = "vtk";


}

#endif
