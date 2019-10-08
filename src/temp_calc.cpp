#include "temp_calc.h"

/* ------------------------- Temperature calculation ----------------------- */

void calc_temperature(  ScalarField & Temperature, ScalarField & Temperature_old,
                        VectorField & ElectricField, ScalarField & K, 
                        ScalarField & Sigma, long double curr_time, int it_number,
                        long double dt)
{

    long double mid1, mid2, dk_by_gradtemp_dx, dk_by_gradtemp_dy,
                dk_by_gradtemp_dz;
    long double aux, aux_1, aux_2;
    (void)par::anodes;
    (void)par::cathodes;
  
    #pragma omp parallel for collapse(3) private(mid1, mid2, dk_by_gradtemp_dx, dk_by_gradtemp_dy,dk_by_gradtemp_dz, aux, aux_1, aux_2 ) schedule(static)
    for( int k = 1; k < par::kk - 1; k++ )
    {
        for( int j = 1; j < par::jj - 1; j++ )
        {
            for( int i = 1; i < par::ii - 1; i++ )
            {
                if( !par::mesh.is_electrode(i,j,k) )
                {

                    // Precalculating elements
                    mid1 = ( K(i+1,j,k) + K(i,j,k) ) / 2.;
                    aux_1 = mid1
                            * ( Temperature_old(i+1,j,k) - Temperature_old(i,j,k) ) / par::dx;
                    mid2 = ( K(i,j,k) + K(i-1,j,k) ) / 2.;
                    aux_2 = mid2
                            * ( Temperature_old(i,j,k) - Temperature_old(i-1,j,k) ) / par::dx;
                    dk_by_gradtemp_dx = ( aux_1 - aux_2 ) / par::dx ;
                    
                    mid1 = ( K(i,j+1,k) + K(i,j,k) ) / 2.;
                    aux_1 = mid1
                            * ( Temperature_old(i,j+1,k) - Temperature_old(i,j,k) ) / par::dy;
                    mid2 = ( K(i,j,k) + K(i,j-1,k) ) / 2.;
                    aux_2 = mid2
                            * ( Temperature_old(i,j,k) - Temperature_old(i,j-1,k) ) / par::dy;
                    dk_by_gradtemp_dy = ( aux_1 - aux_2 ) / par::dy;
                    
                    mid1 = ( K(i,j,k+1) + K(i,j,k) ) / 2.;
                    aux_1 = mid1
                            * ( Temperature_old(i,j,k+1) - Temperature_old(i,j,k) ) / par::dz;
                    mid2 = ( K(i,j,k) + K(i,j,k-1) ) / 2.;
                    aux_2 = mid2
                            * ( Temperature_old(i,j,k) - Temperature_old(i,j,k-1) ) / par::dz;
                    dk_by_gradtemp_dz = ( aux_1 - aux_2 ) / par::dz;

                    aux = 0.;

                    // Diffusive term: P9 * div( K * grad(T) )
                    aux +=  1. / (par::rho * par::cp )
                            * ( dk_by_gradtemp_dx + dk_by_gradtemp_dy + dk_by_gradtemp_dz );

                    // Blood impact in temperature change
                    aux -=  1. / (par::rho * par::cp )
                            * par::rho_b * par::w_b * par::c_b * (  Temperature_old(i,j,k) - par::temperature_b );

                    // Metabolic heat generation: qm
                    aux +=  1. / (par::rho * par::cp ) * par::qm;

                    // Joule effect term: sigma * abs( grad(Phi) )^2
                    aux +=  1. / (par::rho * par::cp ) 
                            * Sigma(i,j,k) 
                            * abs(  pow(ElectricField(0,i,j,k),2)
                                    + pow(ElectricField(1,i,j,k),2) 
                                    + pow(ElectricField(2,i,j,k),2) );

                    // Temporal term
                    aux = aux * dt + Temperature_old(i,j,k);

                    if( aux < 0. ){
                        std::cout   << "domain temperature less than zero:"
                                    << i << ", " << j << ", " << k << std::endl;
                        aux = ( Temperature_old(i+1,j,k) + Temperature_old(i-1,j,k) ) / 2.;
                    }
                                        
                    // Set new value
                    Temperature(i,j,k) = aux;

                    // If is not a valid number (NaN) stop the simulation
                    if( isnan(aux) )
                    {
                        std::cout << "  Time: " << curr_time << ", iteration: " << it_number << std::endl;
                        std::cout << "  Exit from 'calc_temperature' due to NaNs" << std::endl;
                        std::cout << "  i, j, k: " << i << "," << j << "," << k << std::endl;

                        exit(EXIT_FAILURE);
                    }
                }
                else
                {
                    // Precalculating elements
                    mid1 = ( K(i+1,j,k) + K(i,j,k) ) / 2.;
                    aux_1 = mid1
                            * ( Temperature_old(i+1,j,k) - Temperature_old(i,j,k) ) / par::dx;
                    mid2 = ( K(i,j,k) + K(i-1,j,k) ) / 2.;
                    aux_2 = mid2
                            * ( Temperature_old(i,j,k) - Temperature_old(i-1,j,k) ) / par::dx;
                    dk_by_gradtemp_dx = ( aux_1 - aux_2 ) / par::dx;
                    
                    mid1 = ( K(i,j+1,k) + K(i,j,k) ) / 2.;
                    aux_1 = mid1
                            * ( Temperature_old(i,j+1,k) - Temperature_old(i,j,k) ) / par::dy;
                    mid2 = ( K(i,j,k) + K(i,j-1,k) ) / 2.;
                    aux_2 = mid2
                            * ( Temperature_old(i,j,k) - Temperature_old(i,j-1,k) ) / par::dy;
                    dk_by_gradtemp_dy = ( aux_1 - aux_2 ) / par::dy ;
                    
                    mid1 = ( K(i,j,k+1) + K(i,j,k) ) / 2.;
                    aux_1 = mid1
                            * ( Temperature_old(i,j,k+1) - Temperature_old(i,j,k) ) / par::dz;
                    mid2 = ( K(i,j,k) + K(i,j,k-1) ) / 2.;
                    aux_2 = mid2
                            * ( Temperature_old(i,j,k) - Temperature_old(i,j,k-1) ) / par::dz;
                    dk_by_gradtemp_dz = ( aux_1 - aux_2 ) / par::dz;

                    aux = 0.;

                    // Diffusive term: P9 * div( K * grad(T) )
                    aux +=  1. / ( par::rho_electrode * par::cp_electrode )
                            * ( dk_by_gradtemp_dx +
                                dk_by_gradtemp_dy +
                                dk_by_gradtemp_dz );

                    // Joule effect is neglectable, E = 0 in electrode. 

                    // Temporal term
                    aux = aux * dt + Temperature_old(i,j,k);

                    if( aux < 0. ){
                        std::cout   << "domain temperature less than zero:"
                                    << i << ", " << j << ", " << k << std::endl;
                        aux = ( Temperature_old(i+1,j,k) + Temperature_old(i-1,j,k) ) / 2.;
                    }

                    // Set new value
                    Temperature(i,j,k) = aux;
                
                    // If is not a valid number (NaN) stop the simulation
                    if( isnan(aux) )
                    {
                        std::cout << "  Time: " << curr_time << ", iteration: " << it_number << std::endl;
                        std::cout << "  Exit from 'calc_temperature' due to NaNs" << std::endl;
                        std::cout << "  i, j, k: " << i << "," << j << "," << k << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }

              }
          }
      }

    // Main computation: neumann boundary conditions
    int i = par::ii - 1;
    #pragma omp parallel for collapse(2) private(aux) schedule(static)
    for( int j = 0; j < par::jj; j++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {
              aux = Temperature(i-1,j,k);
              Temperature(i,j,k) = aux;
        }
    }

    i = 0;
    #pragma omp parallel for collapse(2) private(aux) schedule(static)
    for( int j = 0; j < par::jj; j++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {
            aux = Temperature(i+1,j,k);
            Temperature(i,j,k) = aux;
        }
    }

    int j = par::jj - 1;
    #pragma omp parallel for collapse(2) private(aux) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {
            aux = Temperature(i,j-1,k);
            Temperature(i,j,k) = aux;
        }
    }

    j = 0;
    #pragma omp parallel for collapse(2) private(aux) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {
              aux = Temperature(i,j+1,k);
              Temperature(i,j,k) = aux;
        }
    }

    int k = par::kk - 1;
    #pragma omp parallel for collapse(2) private(aux) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int j = 0; j < par::jj; j++ )
        {
            aux = Temperature(i,j,k-1);
            Temperature(i,j,k) = aux;
        }
    }

    k = 0;
    #pragma omp parallel for collapse(2) private(aux) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int j = 0; j < par::jj; j++ )
        {
            if(par::mesh.is_electrode(i,j,k)) 
            {
                long double r = - par::h * par::dz / K(i,j,k); 
                aux = ( r * par::temp_air - Temperature(i,j,k+1) ) / ( r - 1 );
                Temperature(i,j,k) = aux;
            } 
            else 
            {
                aux = Temperature(i,j,k+1);
                Temperature(i,j,k) = aux;
            }
        }
    }

}

