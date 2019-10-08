#include "par.h"
#include "mesh.h"
#include "scalar_field.h"
#include "vector_field.h"
#include "temp_calc.h"
#include "electrics_calc.h"
#include "save.h"

/* Main simulation */
int main()
{

    /* Starting time calc. */
    long double begin = time(NULL);

    /* Declare dependent variables */
    /* Phi in time N+1, new iteration */
    ScalarField Phi(par::mesh, par::phi_init); 
    /* Phi in time N+1, last iteration. Needed in semi-implicit method */
    ScalarField Phi_aux(par::mesh, par::phi_init); 
    /* Phi in time N */
    ScalarField Phi_old(par::mesh, par::phi_init); 
    /* Zero Electric Potential in time N, it is used for an optimization in the off pulse */
    ScalarField PhiZero(par::mesh, 0.); 
    /* Electric Potential in time N, it is used for an optimization in the off pulse */
    ScalarField * pPhi = NULL;  
    /* Electric Field in time N */
    VectorField ElectricField(par::mesh, par::ef_init); 
    /* Electric Field in time N, it is used for an optimization in the off pulse */
    VectorField ElectricFieldZero(par::mesh, 0.); 
    /* Electric Field pointer, it is used for an optimization in the off pulse */
    VectorField * pElectricField = NULL; 
    /* Sigma in time N+1 */
    ScalarField Sigma(par::mesh, par::sigma_domain_init); 
    /* Current density, J in time N+1 */
    VectorField CurrentDensity(par::mesh, 0.); 
    /* Temperature in time N+1, new iteration */
    ScalarField Temperature(par::mesh, par::temperature_init); 
    /* Temperature in time N */
    ScalarField Temperature_old(par::mesh, par::temperature_init); 
    /* Thermal conductivity in time N */
    ScalarField K(par::mesh, par::k_domain_init); 

    (void)par::anodes;
    (void)par::cathodes;
    (void)par::ii;
    (void)par::jj;
    (void)par::kk;
    (void)par::dx;
    (void)par::dy;
    (void)par::dz;

    /* Declare state variables */
    int it_number = 0;
    int save_counter = 0;
    long double curr_time = 0.;
    long double pulse_time_acum = 0.;
    bool pulse_on = true;
    int pulse = 0;
    bool log_created = false;
    long double q_accum = 0;
    long double dt = par::dt_on_pulse;
    
    /* Parameters */
    save_parameters();

    /* ----------------- Initial and boundary conditions ------------------- */

    std::cout << "Initiating variables..." << std::endl;

    /* Set Sigma in anode and cathode */
    set_in_electrodes(Sigma, par::sigma_electrode_init, par::mesh.get_anodes(),
                      par::mesh.get_anode_number() );
    set_in_electrodes(Sigma, par::sigma_electrode_init, par::mesh.get_cathodes(),
                      par::mesh.get_cathode_number() );

    /* Set K in anode and cathode */
    set_in_electrodes(K, par::k_electrode_init, par::mesh.get_anodes(),
                     par::mesh.get_anode_number() );
    set_in_electrodes(K, par::k_electrode_init, par::mesh.get_cathodes(),
                    par::mesh.get_cathode_number() );
    
    /* Calc. initial Phi. Phi is constant along the simulation */
    init_phi(Phi);
    Phi_aux = Phi;
    calc_nonlinear_phi( Phi, Phi_aux, Sigma, curr_time, it_number,
                        par::max_sub_it_init_Phi);
    pPhi = &Phi;

    /* Calc. initial electric field */
    calc_electric_field(ElectricField, Phi);
    pElectricField = &ElectricField;

    /* Calc. Current Density = ElectricField * Sigma */
    CurrentDensity.set_coordinate(0, pElectricField->get_coordinate(0) * Sigma );
    CurrentDensity.set_coordinate(1, pElectricField->get_coordinate(1) * Sigma );
    CurrentDensity.set_coordinate(2, pElectricField->get_coordinate(2) * Sigma );

    /* --------------------------- Simulation ------------------------------ */

    std::cout << "Simulation begins...";

    /* Temporal iterations */
    while( curr_time <= par::max_time )
    {
        
        /* Save data */
            
        if( it_number % par::save_step == 0 )
        {
            save(curr_time, save_counter, *pPhi, *pElectricField, Sigma, 
                 CurrentDensity, Temperature);
        }

        /* Log iteration data */

        if( it_number % par::log_step == 0 )
        {   
            log( curr_time, pulse, Temperature, *pPhi, Sigma, CurrentDensity,
                 q_accum, log_created);
        }

        /* Calc. in domain boundaries */

        /* Pulse change: On --> Off */
        if( pulse_on && ( pulse_time_acum > par::on_pulse_time) )
        {
            /* Save and log */
            save(curr_time, save_counter, *pPhi,*pElectricField, Sigma, 
                 CurrentDensity, Temperature);
            log( curr_time, pulse,  Temperature, *pPhi, Sigma, CurrentDensity,
                 q_accum, log_created);

            /* Some control variables */
            pulse_on = false;
            pulse_time_acum = 0.;
            
            /* Increment dt due to off pulse is very long */
            dt = par::dt_off_pulse;
            
            /* Set Electric potential and Electric Field to zero. When pulse */
            /* is off, electric potential and field are zero at any point. */
            pPhi = &PhiZero;
            pElectricField = &ElectricFieldZero;

            /* Calc. Current Density */
            /* CurrentDensity = ElectricField * Sigma; */
            CurrentDensity.set_coordinate(0, 
            	pElectricField->get_coordinate(0) * Sigma );
            CurrentDensity.set_coordinate(1,
            	pElectricField->get_coordinate(1) * Sigma );
            CurrentDensity.set_coordinate(2,
            	pElectricField->get_coordinate(2) * Sigma );
        }
        /* Pulse change: Off --> On */
        else if( !pulse_on && ( pulse_time_acum > par::off_pulse_time) ) 
        {
            /* Save and log */
            save( curr_time, save_counter, *pPhi,*pElectricField, Sigma,
                  CurrentDensity, Temperature);
            log( curr_time, pulse,  Temperature, *pPhi, Sigma, CurrentDensity,
                 q_accum, log_created);

            /* Some control variables */
            pulse_on = true;
            pulse_time_acum = 0.;
            
            /* Decrease dt due to on pulse is very short */
            dt = par::dt_on_pulse;
            
            /* Incremente number of pulses */
            pulse++;
            
            /* Set Electric potential and Field to the values before pulse off. */
            pPhi = &Phi;
            pElectricField = &ElectricField;

            /* Calc. Current Density */
            /* CurrentDensity = ElectricField * Sigma; */
            CurrentDensity.set_coordinate(0, pElectricField->get_coordinate(0) * Sigma );
            CurrentDensity.set_coordinate(1, pElectricField->get_coordinate(1) * Sigma );
            CurrentDensity.set_coordinate(2, pElectricField->get_coordinate(2) * Sigma );
        }
        
        /* Calc. Electric Potential and Electric Field based on Sigma. 
        *  Only if pulse is on. */ 
        if( pulse_on )
        {
            /* Calc. Electric Potential */
            calc_nonlinear_phi( Phi, Phi_aux, Sigma, curr_time, it_number);
                        
            /* Calc. Electric Field */
            calc_electric_field(ElectricField, Phi);
        }
        
        /* Calc. Sigma (Electrical Conductivity) */
        calc_sigma( Sigma, Temperature, curr_time, it_number);


        /* Calc. and update Temperature */
        calc_temperature(Temperature, Temperature_old, *pElectricField, K,
        				 Sigma, curr_time, it_number, dt);
        Temperature_old = Temperature;

        /* Calc. accum. elec. charge */
        q_accum += calc_electric_current( CurrentDensity ) * dt;

        /* Increment time and counters */
        it_number++;
        curr_time += dt;
        pulse_time_acum += dt;

    }

    /* Ending time calc. */
    std::cout << "Simulation ends" << std::endl;
    long double time_elapsed = time(NULL) - begin;
    std::cout << "Time elapsed: " << time_elapsed / 60. << " minutes." << std::endl;

    return 0;
}


