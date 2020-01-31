// ----------------------------------------------------------------------------
//! \file Tools.h
//
//! \brief This class contains a set of tools for numerical aspects
//
// ----------------------------------------------------------------------------

#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>
#include <stdio.h>
#include <boost/math/special_functions/bessel.hpp> 
// #include <gsl/gsl_sf_bessel.h>

#define __header(__msg,__txt) std::cout << "\t[" << __msg << "] " << __FILE__ << ":" << __LINE__ << " (" \
<< __FUNCTION__ << ") " << __txt << std::endl

#define WARNING(__txt) {int __rk; MPI_Comm_rank( MPI_COMM_WORLD, &__rk ); if (__rk==0) {std::cout << "\t[WARNING] " << __txt << std::endl;}}

#define ERROR(__txt) {__header("ERROR", __txt); MPI_Finalize(); exit(EXIT_FAILURE);}

class Tools
{
    public:
        
        //! Load repartition in 1d between MPI processes.
        //! This function returns tables of indexes and length for all rank
        static void distributeArray(
            int nb_chunks,
            int nb_elems,
            int *imin_table,
            int *length_table );
            
        //! Computation of the Gauss-Legendre abscissa and weight
        // static void GaussLegendreQuadrature( double x_min, double x_max, double * roots,
        //                                  double *weights, int number_of_roots, double eps );
        static void GaussLegendreQuadrature( double xmin, double xmax, double *x,
                                                 double *w, int nb_iterations, double eps );
        
        //! This function returns true/flase whether the file exists or not
        //! \param file file name to test
        static bool fileCreated( const std::string &filename ) ;
        
        //! This function computes the second kind modified Bessel function K
        static double BesselK(double nu, double z);
        
        //! This method provides the asymptotic behaviors for large
        //! z of the second kind modified Bessel function K
        static double asymptoticBesselK(double nu, double z);
                                         
};

#endif
