#include "save.h"

/* ------------------------ Save and log functions ------------------------- */

void save(  long double curr_time, int& save_counter, ScalarField & Phi,
            VectorField & ElectricField, ScalarField & Sigma,
            VectorField & CurrentDensity, ScalarField & Temperature)
{   

    /* Save scalar fields */
    std::ostringstream stream;
    stream  << std::setfill('0') << std::setw(10) << std::fixed 
            << std::setprecision(0) << save_counter;
    Temperature.save("Temp_" + stream.str(), par::save_format);
    Phi.save("Phi_" + stream.str(), par::save_format);
    ElectricField.save("EField_" + stream.str(), par::save_format);
    ElectricField.get_norm().save("EFieldNorm_" + stream.str(), par::save_format);
    Sigma.save("Sigma_" + stream.str(), par::save_format);
    CurrentDensity.save("CurrentDensity_" + stream.str(), par::save_format);
    CurrentDensity.get_norm().save("CurrentDensityNorm_" + stream.str(), par::save_format);

    /* Save save_counter vs time */
    std::ofstream file;
    file.open("save_counter_vs_time.csv", std::ios_base::app);
    file << save_counter << "," << curr_time << std::endl;
    file.close();

    /* Update save counter */
    save_counter++;

}

void log(   long double curr_time, int pulse, ScalarField & Temperature,
            ScalarField & Phi, ScalarField & Sigma, VectorField & ElectricField, VectorField & CurrentDensity,
            long double q_accum, bool & log_created )

{

    /* Open file */
    std::ofstream file("log.csv", std::ios::app);

    /* Write header */
    if (!log_created)
    {
        file    <<  "Time [s],"
                    "Pulse number,"
                    "Potential difference [V],"
                    "Electric conductivity (middle) [S/m],"
                    "Electric current [A],"
                    "Electric charge [C],"
                    "Temperature (middle) [K],"
                    "Electric field (middle) [V/m]"
                << std::endl;
        log_created = true;
    }

    /* Write content */

    long double time = curr_time;
    long double potential_difference = 
                (   Phi( par::anodes[0][0], par::anodes[0][1], par::anodes[0][2]) - 
                    Phi( par::cathodes[0][0], par::cathodes[0][1], par::cathodes[0][2]) );
    long double electricConductivity = 
                Sigma( par::ii / 2, par::jj / 2, par::kk / 2);
    long double electricCurrent =
                calc_electric_current( CurrentDensity );
    long double temperature_middle = Temperature( par::ii / 2, par::jj / 2, par::kk / 2);
    long double ef_middle = ElectricField.get_norm().get(par::ii / 2, par::jj / 2, par::kk / 2);
    
    file    << time << ", "
            << pulse << ", "
            << potential_difference << ", "
            << electricConductivity << ", "
            << electricCurrent << ", "
            << q_accum  << ", "
            << temperature_middle << ", "
            << ef_middle
            << std::endl;

    /* Close file */
    file.close();
}


void save_parameters() {

    std::ofstream file ;

    /* Open file */
    file .open("parameters.dat");

    /* Write content */

    file  << "/* -- Universal constants --------------------- */" << std::endl;

    /* Tissue temperature [K] */
    file  << "  temp: " << par::temp << std::endl;  
    /* Air temperature. [K] */
    file  << "  temp_air: " << par::temp_air << std::endl;  

    file  << "/* -- Electrodes ------------------------------ */" << std::endl;

    /* Electrode length [m] */
    file  << "  electrode_length: " << par::electrode_length << std::endl;  
    /* Electrode width [m] */
    file  << "  electrode_width: " << par::electrode_width << std::endl;  
    /* Electrode thickness [m] */
    file  << "  electrode_thickness: " << par::electrode_thickness << std::endl;  
    /* Space between anode-cathode, horizontally [m] */
    file  << "  gap_anode_cathode: " << par::gap_anode_cathode << std::endl;  
    /* Type of electrode: "plates" or "needles" */
    file  << "  electrode_type: " << par::electrode_type << std::endl;  
    /* Space between anode-anode and cathode-cathode. [m] */
    /* Use only with needle electrodes */
    file  << "  gap_elect_elect: " << par::gap_elect_elect << std::endl;  
    /* Use only with needle electrodes. No. of electrodes. */
    file  << "  no_electrodes: " << par::no_electrodes << std::endl;  

    file  << "/* -- Domain parameters ----------------------- */" << std::endl;

    /* Maximum length of the x axis [m] */
    file  << "  x_max: " << par::x_max << std::endl;  
    /* Maximum length of the y axis [m] */
    file  << "  y_max: " << par::y_max << std::endl;  
    /* Maximum length of the z axis [m] */
    file  << "  z_max: " << par::z_max << std::endl;  

    file  << "/* -- Treatment parameters -------------------- */" << std::endl;

    /* Voltage to distance ratio [V m^-1] */
    file  << "  volt_to_dist: {";
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::volt_to_dist[i] << "}";
        else
            file << par::volt_to_dist[i] << ", ";
    }
    file << std::endl;  
    /* Pulse repetitions*/
    // file << par::pulse_repetitions << std::endl;
    /* On time [s] */
    file  << "  on_pulse_times: {";
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::on_pulse_times[i] << "}";
        else
            file << par::on_pulse_times[i] << ", ";
    }
    file << std::endl;  
    /* No. of pulses */
    file  << "  no_pulses: " << par::no_pulses << std::endl;  
    
    file  << "  no_elems_per_cycle: " << par::no_elems_per_cycle << std::endl;  

    file  << "/* - Other parameters, calculated from the above */" << std::endl;

    /* Total pulse time [s] */
    // file  << "  total_pulse_time: " << par::total_pulse_time << std::endl;  
    /* Off pulse times [s] */
    file  << "  off_pulse_times: {"; 
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::off_pulse_times[i] << "}";
        else
            file << par::off_pulse_times[i] << ", ";
    }
    file << std::endl;   
    /* Simulation total time [s] */
    file  << "  max_time: " << par::max_time << std::endl;  
    /* Maximum voltages [V] */
    file  << "  max_voltages: {"; 
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::max_voltages[i] << "}";
        else
            file << par::max_voltages[i] << ", ";
    }
    file << std::endl; 
    /* Potential in the anode [V] */
    file  << "  phi_anodes: {"; 
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::phi_anodes[i] << "}";
        else
            file << par::phi_anodes[i] << ", ";
    }
    file << std::endl; 
    /* Potential in the cathode [V] */
    file  << "  phi_cathode: " << par::phi_cathode << std::endl;  
    /* Initial electric potential [V] */
    file  << "  phi_init: " << par::phi_init << std::endl;  
    /* Initial electric filed [V m^-1] */
    file  << "  ef_init: " << par::ef_init << std::endl;  

    file  << "/* -- Tissue paramters ------------------------ */" << std::endl;

    /* Initial temperature [K] */
    file  << "  temperature_init: " << par::temperature_init << std::endl;  
    /* Thermal conductivity [ W m^-1 K^-1 ] */
    file  << "  k_domain_init: " << par::k_domain_init << std::endl;  
    /* Electrical conductivity [ S m^-1 ] */
    file  << "  sigma_domain_init: " << par::sigma_domain_init << std::endl;  
    /* Tissue density [ kg m^3 ] */
    file  << "  rho: " << par::rho << std::endl;  
    /* Heat capacity [J kg^-1 K^-1] */
    file  << "  cp: " << par::cp << std::endl;  
    /* Metabolic heat generation [ W m^-3 ] */
    file  << "  qm: " << par::qm << std::endl;  

    file  << "/* -- Blood parameters ------------------------ */" << std::endl;

    /* Temperature of the arterial blood [K] */
    file  << "  temperature_b: " << par::temperature_b << std::endl;  
    /* Density of the blood [ kg m^3 ] */
    file  << "  rho_b: " << par::rho_b << std::endl;  
    /* Heat capacity of the blood [J kg^-1 K^-1] */
    file  << "  c_b: " << par::c_b << std::endl;  
    /* Blood profusion [ s^-1 ] */
    file  << "  w_b: " << par::w_b << std::endl;  

    file  << "/* -- Electrode parameters -------------------- */" << std::endl;

    /* Thermal conductivity [ W m^-1 K^-1 ] */
    file  << "  k_electrode_init: " << par::k_electrode_init << std::endl;  
    /* Electrical conductivity [ S m^-1 ] */
    file  << "  sigma_electrode_init: " << par::sigma_electrode_init << std::endl;  
    /* Electrode density [ kg m^-3 ] */
    file  << "  rho_electrode: " << par::rho_electrode << std::endl;  
    /* Heat capacity [J kg^-1 K^-1] */                             
    file  << "  cp_electrode: " << par::cp_electrode << std::endl;  
                                  
    file  << "/* -- Boundary parameters --------------------- */" << std::endl;

    /* Convective heat-transfer coefficient [ W m^-2 K^-1 ]  */
    file  << "  h: " << par::h << std::endl;  
    /* Electrical conductivity increment with regard to temperature */
    file  << "  alpha0: " << par::alpha0 << std::endl;  

    file  << "/* -- Geometry and Mesh Parameters ------------ */" << std::endl;

    /* No. of nodes between anodes and cathodes */
    file  << "  resolution: " << par::resolution << std::endl;  
    /* Mesh */

    /* No. of x axis divisions */
    file  << "  ii: " << par::ii << std::endl;  
    /* No. of y axis divisions */
    file  << "  jj: " << par::jj << std::endl;  
    /* No. of z axis divisions */
    file  << "  kk: " << par::kk << std::endl;  
    /* Space between nodes in x axes */
    file  << "  dx: " << par::dx << std::endl;  
    /* Space between nodes in y axes */
    file  << "  dy: " << par::dy << std::endl;  
    /* Space between nodes in z axes */
    file  << "  dz: " << par::dz << std::endl;  

    file  << "/* -- Other numerical parameters -------------- */" << std::endl;

    /* Convergence tolerance */
    file  << "  tol: " << par::tol << std::endl;  
    /* Max. sub iteration number for Phi */
    file  << "  max_sub_it_Phi: " << par::max_sub_it_Phi << std::endl;  
    /* Max. sub iteration number for intial condition of Phi */
    file  << "  max_sub_it_init_Phi: " << par::max_sub_it_init_Phi << std::endl;  
    /* Relaxation */
    file  << "  omega: " << par::omega << std::endl;  

    file  << "/* -- Pulse parameters ------------------------ */" << std::endl;

    /* Time step during on pulse [s] */
    file  << "  dts_on_pulse: {";
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::dts_on_pulse[i] << "}";
        else
            file << par::dts_on_pulse[i] << ", ";
    }
    file << std::endl;  
    /* Time step during off pulse [s] */
    file  << "  dts_off_pulse: {";
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::dts_off_pulse[i] << "}";
        else
            file << par::dts_off_pulse[i] << ", ";
    }
    file << std::endl;  

    file  << "/* -- Save and log parameters  ---------------- */" << std::endl;

    /* Save each ... iterations */
    file  << "  save_steps: {";
    for (int i = 0; i < par::no_elems_per_cycle; i++)
    {
        if(i == par::no_elems_per_cycle - 1)
            file << par::save_steps[i] << "}";
        else
            file << par::save_steps[i] << ", ";
    }
    file << std::endl;  
    /* Log each ... iterations */
    file  << "  log_step: " << par::log_step << std::endl;  
    /* "vtk" or "csv" */
    file  << "  save_format: " << par::save_format << std::endl;  


    file  << "/* -- Anode and cathode positions  ------------ */" << std::endl;

    file  << "  anodes" << std::endl;
    for( int i = 0; i < par::mesh.get_anode_number(); i++ )
        file  << par::anodes[i][0] << ",: " << par::anodes[i][1] << ",: " << par::anodes[i][2] << std::endl;

    file  << "  cathodes" << std::endl;
    for( int i = 0; i < par::mesh.get_cathode_number(); i++ )
        file  << par::cathodes[i][0] << ",: " <<  par::cathodes[i][1] << ",: " <<  par::cathodes[i][2] << std::endl;

    file  << std::endl;


    /* Close file */
    file .close();

}
