// ----------------------------------------------------------------------------
//! \file Tools.cpp
//
//! \brief This class contains a set of tools for numerical aspects
//
// ----------------------------------------------------------------------------

#include "Tools.h"

// ----------------------------------------------------------------------------
//! \brief Distribute equally 1D array into chunks
//! This function returns tables of indexes and length for each chunk
//
//! \param nb_chunks Total number of chunks
//! \param nb_elems Total number of element to be distributed
//! \param imin_table Index of the first element for each chunk
//! \param length_table Number of element for each chunk
// ----------------------------------------------------------------------------
void Tools::distributeArray(
    int nb_chunks,
    int nb_elems,
    int *imin_table,
    int *length_table )
{

    // If more chunks than elements,
    // only a part of the processes will work
    if( nb_chunks >= nb_elems ) {
        #pragma omp simd
        for( int chunk = 0 ; chunk < nb_elems ; chunk ++ ) {
            imin_table[chunk] = chunk;
            length_table[chunk] = 1;
        }
        #pragma omp simd
        for( int chunk = nb_elems ; chunk < nb_chunks ; chunk ++ ) {
            imin_table[chunk] = nb_elems;
            length_table[chunk] = 0;
        }
    } else {
    
        int quotient;
        int remainder;
        
        // Part of the load equally distributed
        quotient = nb_elems/nb_chunks;
        
        // Remaining load to be distributed after balanced repartition
        remainder = nb_elems%nb_chunks;
        
        #pragma omp simd
        for( int chunk = 0 ; chunk < remainder ; chunk ++ ) {
            imin_table[chunk] =  chunk*quotient+chunk;
            length_table[chunk] = quotient + 1;
        }
        #pragma omp simd
        for( int chunk = remainder ; chunk < nb_chunks ; chunk ++ ) {
            imin_table[chunk] = remainder + chunk*quotient;
            length_table[chunk] = quotient;
        }
    }
}


// ----------------------------------------------------------------------------
//! \brief Computation of the roots and weights for the
//! Gauss-Legendre Integration between x_min and x_max.
//
//! \param x_min minimum integration boundary
//! \param x_max maximum integration boundary
//! \param roots array of abscissa
//! \param weights array of weight for integration
//! \param number_of_roots number of iteration for integration (array size)
//! \param eps accuracy threshold for coef computation
//! For more information:
//! - https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
//! - Numerical recipes
// ----------------------------------------------------------------------------
void Tools::GaussLegendreQuadrature( double x_min, double x_max, double *roots,
        double *weights, int number_of_roots, double eps )
{
    
    // Checks
    if( number_of_roots <= 0 ) {
        ERROR("Gauss-Legendre quadrature: Number of iteration <= 1");
    }
    if( x_max < x_min ) {
        ERROR("Gauss-Legendre quadrature: x_max < x_min");
    }
    if( eps <= 0 ) {
        ERROR("Gauss-Legendre quadrature: accuracy threshold epsilon <= 0");
    }
    
    // Legendre Polynomials
    double P_prev;
    double P_curr;
    double P_next;
    
    // Derivative of P
    double P_derivative;
    
    // Newton-Raphson integration
    double root;
    double root_prev;

    // The roots are symmetric, so we only find half of them.
    int half_number_of_roots=( number_of_roots+1 )/2;
    double x_average=0.5*( x_min+x_max );
    double x_half_length=0.5*( x_max-x_min );
    
    // Loop over the roots
    for( int i_root=0; i_root<=half_number_of_roots-1; i_root++ ) { /* Loop over the desired roots. */
        
        // First guess for the i-th root of a n-order polynomial Pn for the Newton-Raphson iteration method
        root=cos( M_PI*( i_root+1.0-0.25 )/( number_of_roots+0.5 ) );
        

        do {
            
            // Initialization of the first polynomials
            P_next=root;
            P_curr=1.0;
            
            // We use the recurrence of Bonnet to compute the Legendre polynomial at order n+1
            for( int order=2; order<=number_of_roots; order++ ) {
                
                P_prev=P_curr;
                P_curr=P_next;
                
                // Recurrence
                P_next=( ( 2.0*order-1.0 )*root*P_curr-( order-1.0 )*P_prev )/order;
            }

            // Derivatibve computation
            P_derivative=number_of_roots*( root*P_next-P_curr )/( root*root-1.0 );
            
            root_prev=root;
            
            // Newton-Raphson ierative method to get closer to the i-th root
            root=root_prev-P_next/P_derivative;
            
        } while( fabs( root-root_prev ) > eps );
        
        // Update of the root
        // The algorithm is between -1 and 1 so we rescale between x_min and x_max
        roots[i_root]=x_average-x_half_length*root;
        roots[number_of_roots-1-i_root]=x_average+x_half_length*root;
        
        // Update of the Weight
        weights[i_root]=2.0*x_half_length/( ( 1.0-root*root )*P_derivative*P_derivative );
        weights[number_of_roots-1-i_root]=weights[i_root];

    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function returns true/flase whether the file exists or not
//! \param file file name to test
// ---------------------------------------------------------------------------------------------------------------------
bool Tools::fileCreated( const std::string &filename )
{
    std::ifstream file( filename.c_str() );
    return !file.fail();
}

// ----------------------------------------------------------------------------
//! \brief This function computes the second kind modified Bessel function K
// ----------------------------------------------------------------------------
double Tools::BesselK(double nu, double z)
{
    double K;
    if (z > 3000.0) {
        K = Tools::asymptoticBesselK(nu, z);
    } else {
        K = boost::math::cyl_bessel_k( nu, z);
        //K = std::cyl_bessel_k( nu, z);
        //K = gsl_sf_bessel_Kn( nu, z);
    }
    return K;
}

// ----------------------------------------------------------------------------
//! \brief Asymptotic Behavior of Bessel K
//! See Abramowitz and Stegun
// ----------------------------------------------------------------------------
double Tools::asymptoticBesselK(double nu, double z)
{
    double mu = 4 * nu*nu;
    double iz = 1/z;
    double D = iz*0.125;
    double C = std::sqrt(M_PI*0.5*iz)*std::exp(-z);
    double K;
    K = 1 + (mu-1)*0.125*iz
          + (mu-1)*(mu-9)*0.5*std::pow(D,2)
          + (mu-1)*(mu-9)*(mu-25)*std::pow(D,3)/6.0;
    K *= C;
    
    return K;
}
