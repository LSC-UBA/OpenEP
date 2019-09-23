#include "electrics_calc.h"

/* --------------------- Electrical variable calculation ------------------- */

/* Calculate initial distribution of Phi */
void init_phi(ScalarField & Phi)
{
    #pragma omp parallel for collapse(3) schedule(static)
    for( int k = 1; k < par::kk - 1; k++ )
    {
        for( int j = 1; j < par::jj - 1; j++ )
        {
            for( int i = 1; i < par::ii - 1; i++ )
            {
                if( i < par::ii / 2 )
                    Phi.set(i,j,k,par::max_voltage);
                else
                    Phi.set(i,j,k,0.);
            }
        }
    }
}

/* Calculate Phi: div( sigma grad(Phi) ) = 0 */
void calc_nonlinear_phi(    ScalarField & Phi, ScalarField & Phi_aux,
                            ScalarField & Sigma, long double curr_time,
                            int it_number, int max_sub_it)
{

    long double sigmai0, sigmai1, sigmaj0, sigmaj1, sigmak0, sigmak1, 
                num, den, res;
    long int subiteration_number = 0;
        
    /* Subiteration for solving eq. system in domain */
    while( subiteration_number < max_sub_it )
    {
        #pragma omp parallel for collapse(3) private(sigmai0, sigmai1, sigmaj0, sigmaj1, sigmak0, sigmak1, num, den, res) schedule(static)
        for( int k = 1; k < par::kk - 1; k++ )
        {
            for( int j = 1; j < par::jj - 1; j++ )
            {
                for( int i = 1; i < par::ii - 1; i++ )
                {

                    if( par::mesh.is_elec_isolating_material(i,j,k) )
                    {
                        Phi.set(i,j,k,0.);
                    }
                    else if( !par::mesh.is_electrode(i,j,k) )
                    {
                        
                        num = 0.;
                        
                        sigmai1 = ( Sigma(i+1,j,k) + Sigma(i,j,k) ) / 2.;
                        sigmai0 = ( Sigma(i,j,k) + Sigma(i-1,j,k) ) / 2.;
                        num += par::dx * ( sigmai1 * Phi_aux(i+1,j,k) \
                        	+ sigmai0 * Phi_aux(i-1,j,k) ) ;
                        
                        sigmaj1 = ( Sigma(i,j+1,k) + Sigma(i,j,k) ) / 2.;
                        sigmaj0 = ( Sigma(i,j,k) + Sigma(i,j-1,k) ) / 2.;
                        num += par::dy * ( sigmaj1 * Phi_aux(i,j+1,k) \
                        	+ sigmaj0 * Phi_aux(i,j-1,k) );
                        
                        sigmak1 = ( Sigma(i,j,k+1) + Sigma(i,j,k) ) / 2.;
                        sigmak0 = ( Sigma(i,j,k) + Sigma(i,j,k-1) ) / 2.;
                        num += par::dz * ( sigmak1 * Phi_aux(i,j,k+1) \
                        	+ sigmak0 * Phi_aux(i,j,k-1) ) ;
                        
                        den = 0.;
                        
                        den += par::dx * (sigmai1 + sigmai0 );
                        den += par::dy * (sigmaj1 + sigmaj0 );
                        den += par::dz * (sigmak1 + sigmak0 );
                        
                        res = num / den;

                        res = res * par::omega 
                            + Phi_aux(i,j,k) * ( 1. - par::omega );
                        
                        Phi.set(i,j,k,res);

                        /* If is not a valid number (NaN) stop the simulation */
                        if( isnan(res) )        
                        {
                                std::cout << "    Time: " << curr_time
                                << ", iteration: " << it_number << std::endl;
                                std::cout 
                                << "    Exit from 'calc_nonlinear_phi' due to NaNs"
                                << std::endl;
                                exit(EXIT_FAILURE);
                        }

                    }

                }
            }
        }

        /* Neumann boundary conditions */

        int i = par::ii - 1;
        #pragma omp parallel for collapse(2) private(res) schedule(static)
        for( int j = 0; j < par::jj; j++ )
        {
            for( int k = 0; k < par::kk; k++ )
            {
                res = Phi(i-1,j,k);
                Phi.set(i,j,k,res);
            }
        }

        i = 0;
        #pragma omp parallel for collapse(2) private(res) schedule(static)
        for( int j = 0; j < par::jj; j++ )
        {
            for( int k = 0; k < par::kk; k++ )
            {
                res = Phi(i+1,j,k);
                Phi.set(i,j,k,res);
            }
        }

        int j = par::jj - 1;
        #pragma omp parallel for collapse(2) private(res) schedule(static)
        for( int i = 0; i < par::ii; i++ )
        {
            for( int k = 0; k < par::kk; k++ )
            {
                res = Phi(i,j-1,k);
                Phi.set(i,j,k,res);
            }
        }

        j = 0;
        #pragma omp parallel for collapse(2) private(res) schedule(static)
        for( int i = 0; i < par::ii; i++ )
        {
            for( int k = 0; k < par::kk; k++ )
            {
                res = Phi(i,j+1,k);
                Phi.set(i,j,k,res);
            }
        }

        int k = par::kk - 1;
        #pragma omp parallel for collapse(2) private(res) schedule(static)
        for( int i = 0; i < par::ii; i++ )
        {
            for( int j = 0; j < par::jj; j++ )
            {
                    res = Phi(i,j,k-1);
                    Phi.set(i,j,k,res);
            }
        }

        k = 0;
        #pragma omp parallel for collapse(2) private(res) schedule(static)
        for( int i = 0; i < par::ii; i++ )
        {
            for( int j = 0; j < par::jj; j++ )
            {
                    res = Phi(i,j,k+1);
                    Phi.set(i,j,k,res);
            }
        }

        /* Update subiteration concentratio and number */

        Phi_aux = Phi;
        subiteration_number++;

    }

}

/* Calculate electric field: E = - grad Phi */
void calc_electric_field(VectorField & ElectricField, ScalarField & Phi)
{
        long double aux;
        int i, j, k;
    
        #pragma omp parallel for collapse(3) private(aux) schedule(static)
        for( int k = 1; k < par::kk - 1; k++ )
        {
            for( int j = 1; j < par::jj - 1; j++ )
            {
                for( int i = 1; i < par::ii - 1; i++ )
                {
                    if( !par::mesh.is_electrode(i,j,k) 
                        && !par::mesh.is_elec_isolating_material(i,j,k) )
                    {
                    
                        aux = - ( Phi(i+1,j,k) - Phi(i-1,j,k) ) / ( 2 * par::dx );
                        ElectricField.set(0, i, j, k, aux );
                        
                        aux = - (Phi(i,j+1,k) -    Phi(i,j-1,k)) / ( 2 * par::dy );
                        ElectricField.set(1, i, j, k, aux );
                        
                        aux = - (Phi(i,j,k+1) -    Phi(i,j,k-1)) / ( 2 * par::dz );
                        ElectricField.set(2, i, j, k, aux );

                    }
                }
            }
        }
        
        
        i = 0;
        #pragma omp parallel for collapse(2) private(aux) schedule(static)
        for( int k = 0; k < par::kk; k++ )
        {
            for( int j = 0; j < par::jj; j++ )
            {

                aux = ElectricField(0,i+1,j,k);
                ElectricField.set(0, i, j, k, aux);
                
                aux = ElectricField(1,i+1,j,k);
                ElectricField.set(1, i, j, k, aux);
                
                aux = ElectricField(2,i+1,j,k);
                ElectricField.set(2, i, j, k, aux);

            }
        }
        
        i = par::ii - 1;
        #pragma omp parallel for collapse(2) private(aux) schedule(static)
        for( int k = 0; k < par::kk; k++ )
        {
            for( int j = 0; j < par::jj; j++ )
            {

                aux = ElectricField(0,i-1,j,k);
                ElectricField.set(0, i, j, k, aux);
                
                aux = ElectricField(1,i-1,j,k);
                ElectricField.set(1, i, j, k, aux);
                
                aux = ElectricField(2,i-1,j,k);
                ElectricField.set(2, i, j, k, aux);

            }
        }
        
        
        j = 0;
        #pragma omp parallel for collapse(2) private(aux) schedule(static)
        for( int k = 0; k < par::kk; k++ )
        {
            for( int i = 0; i < par::ii; i++ )
            {

                aux = ElectricField(0,i,j+1,k);
                ElectricField.set(0, i, j, k, aux);
                
                aux = ElectricField(1,i,j+1,k);
                ElectricField.set(1, i, j, k, aux);
                
                aux = ElectricField(2,i,j+1,k);
                ElectricField.set(2, i, j, k, aux);

            }
        }
        
        j = par::jj - 1;
        #pragma omp parallel for collapse(2) private(aux) schedule(static)
        for( int k = 0; k < par::kk; k++ )
        {
            for( int i = 0; i < par::ii; i++ )
            {

                aux = ElectricField(0,i,j-1,k);
                ElectricField.set(0, i, j, k, aux);
                
                aux = ElectricField(1,i,j-1,k);
                ElectricField.set(1, i, j, k, aux);
                
                aux = ElectricField(2,i,j-1,k);
                ElectricField.set(2, i, j, k, aux);

            }
        }
        
        k = 0;
        #pragma omp parallel for collapse(2) private(aux) schedule(static)
        for( int j = 0; j < par::jj; j++ )
        {
            for( int i = 0; i < par::ii; i++ )
            {

                aux = ElectricField(0,i,j,k+1);
                ElectricField.set(0, i, j, k, aux);
                
                aux = ElectricField(1,i,j,k+1);
                ElectricField.set(1, i, j, k, aux);
                
                aux = ElectricField(2,i,j,k+1);
                ElectricField.set(2, i, j, k, aux);

            }
        }
        
        k = par::kk - 1;
        #pragma omp parallel for collapse(2) private(aux) schedule(static)
        for( int j = 0; j < par::jj; j++ )
        {
            for( int i = 0; i < par::ii; i++ )
            {

                aux = ElectricField(0,i,j,k-1);
                ElectricField.set(0, i, j, k, aux);
                
                aux = ElectricField(1,i,j,k-1);
                ElectricField.set(1, i, j, k, aux);
                
                aux = ElectricField(2,i,j,k-1);
                ElectricField.set(2, i, j, k, aux);

            }
        }
        
}

/* Calculate electric conductivity: sigma = sigma0 * (1 + alpha (T − temp) */
void calc_sigma( ScalarField & Sigma, ScalarField & Temperature,
                 long double curr_time, int it_number )
{
    long double res = 0;

    #pragma omp parallel for private(res) schedule(static)
    for( int k = 1; k < par::kk - 1; k++ )
    {
        for( int j = 1; j < par::jj - 1; j++ )
        {
            for( int i = 1; i < par::ii - 1; i++ )
            {
                if( !par::mesh.is_electrode(i,j,k)
                    && !par::mesh.is_elec_isolating_material(i,j,k) )
                {

                    res = par::sigma_domain_init 
                        * (1 + par::alpha0 *(Temperature(i,j,k) - par::temp));
                    Sigma.set(i,j,k,res);

                    /* If is not a valid number (NaN) stop the simulation */
                    if( isnan(res) )
                    {
                        std::cout << "    Time: " << curr_time 
                        << ", iteration: " << it_number << std::endl;
                        std::cout <<
                        "    Exit from 'calc sigma' due to NaNs" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                }                 

            }
        }
    }

    /* Main computation: neumann boundary conditions */
    int i = par::ii - 1;
    #pragma omp parallel for collapse(2) private(res) schedule(static)
    for( int j = 0; j < par::jj; j++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {

            res = Sigma(i-1,j,k);
            Sigma.set(i,j,k,res);

        }
    }

    i = 0;
    #pragma omp parallel for collapse(2) private(res) schedule(static)
    for( int j = 0; j < par::jj; j++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {

            res = Sigma(i+1,j,k);
            Sigma.set(i,j,k,res);

        }
    }

    int j = par::jj - 1;
    #pragma omp parallel for collapse(2) private(res) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {

            res = Sigma(i,j-1,k);
            Sigma.set(i,j,k,res);

        }
    }

    j = 0;
    #pragma omp parallel for collapse(2) private(res) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int k = 0; k < par::kk; k++ )
        {

            res = Sigma(i,j+1,k);
            Sigma.set(i,j,k,res);

        }
    }

    int k = par::kk - 1;
    #pragma omp parallel for collapse(2) private(res) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int j = 0; j < par::jj; j++ )
        {

            res = Sigma(i,j,k-1);
            Sigma.set(i,j,k,res);

        }
    }

    k = 0;
    #pragma omp parallel for collapse(2) private(res) schedule(static)
    for( int i = 0; i < par::ii; i++ )
    {
        for( int j = 0; j < par::jj; j++ )
        {

            res = Sigma(i,j,k+1);
            Sigma.set(i,j,k,res);

        }
    }
}


/* Calculate electric current: integral_surface J * ds */
long double calc_electric_current(VectorField & CurrentDensity)
{
    /* The middle between anode and cathode x-distance */
    int off_i = par::anodes[0][0] 
                + floor(abs(par::anodes[0][0] - par::cathodes[0][0]) / 2);

    long double electricCurrent = 0.;

    /* Right anf left sides */
    for(int j = 0; j < par::jj; j++)
    {
        for(int k = 0; k < par::kk; k++)
        {

            /* n = (1,0,0) */
            electricCurrent += CurrentDensity(0, off_i, j, k ) 
                            * par::dy * par::dz;
            /* n = (1,0,0) */
            electricCurrent += - CurrentDensity(0, 0, j, k )
                            * par::dy * par::dz;

        }
    }

    /* Back and forward sides */
    for(int i = 0; i < par::ii; i++)
    {
        for(int k = 0; k < par::kk; k++)
        {

            /* n = (0,1,0) */
            electricCurrent += CurrentDensity(1, i, par::jj-1, k )
                            * par::dx * par::dz;
            /* n = (0,-1,0) */
            electricCurrent += - CurrentDensity(1, i, 0, k )
                            * par::dx * par::dz; 

        }
    }

    /* Down */
    for(int j = 0; j < par::jj; j++)
    {
        for(int i = 0; i < off_i; i++)
        {
            if( !par::mesh.is_electrode(i, j, par::kk-1 ) )
            {
                /* n = (0,0,-1) */
                electricCurrent += - CurrentDensity(2, i, j, par::kk-1 )
                                * par::dx * par::dy;

            }
        }
    }

    /* Up */
    for(int j = 0; j < par::jj; j++)
    {
        for(int i = 0; i < off_i; i++)            
        {
            if( !par::mesh.is_electrode(i, j, 0) )
            {

                /* n = (0,0,1)  */
                electricCurrent += CurrentDensity(2, i, j, 0 )
                                * par::dx * par::dy;

            }
        }
    }

    return electricCurrent;
}

/* Set boundary condition at the electrodes */
void set_in_electrodes(ScalarField & scalarField, const long double val,
                       int ** electrodes, int electrodes_number)
{
    int i, j, k;
    for( int l = 0; l < electrodes_number; l++ )
    {
        i = electrodes[l][0];
        j = electrodes[l][1];
        k = electrodes[l][2];
        scalarField.set(i,j,k,val);
    }
}
