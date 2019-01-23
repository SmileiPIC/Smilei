#include "tabulatedFunctions.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>


using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Inverse error function
// ---------------------------------------------------------------------------------------------------------------------

// method used to load the tabulated function
// ------------------------------------------
void erfinv::prepare()
{
    if( erfinv_tab_.size()==0 ) {
    
        erfinv_tabSize_  = 1000;
        erfinv_xmin_     = 0.0001;
        erfinv_xmax_     = 0.9999;
        erfinv_alpha_    = log( erfinv_xmax_/erfinv_xmin_ )/( double )( erfinv_tabSize_-1 );
        
        erfinv_tab_.resize( erfinv_tabSize_ );
        erfinv_x_.resize( erfinv_tabSize_ );
        
        // -----------------------------------
        // TABULATE THE INVERSE ERROR FUNCTION
        // (using a logarithmic scale)
        // -----------------------------------
        double tiny=1.e-10; // numerical parameter (~precision)
        
        for( unsigned int n=0; n<erfinv_tabSize_; n++ ) {
        
            erfinv_x_[n] = 1.0-erfinv_xmax_ * exp( -( double )( n )*erfinv_alpha_ );
            
            double vl = 0.0;
            double vr = 20.0;
            double vm = 10.0;
            
            while( abs( erfinv_x_[n]-erf( vm ) )>tiny ) {
                vm=0.5*( vl+vr );
                if( erfinv_x_[n]>erf( vm ) ) {
                    vl=vm;
                } else {
                    vr=vm;
                }
            }
            erfinv_tab_[n] = 0.5*( vl+vr );
            
        }//n
        
    } else {
        DEBUG( "trying to call this again!" );
    }//needLoad
    
}

// method used to compute the value for a given x (use linear interpolation on log. scale axis)
// --------------------------------------------------------------------------------------------
double erfinv::call( double x )
{

    double pi  = M_PI;
    double val = 0.0;
    
    if( x <= erfinv_x_[0] ) {
        val = 0.5*sqrt( pi )*x + pi/24.0 *pow( x, 3 );
    } else if( x >= erfinv_x_.back() ) {
        double eta = -log( sqrt( pi )*( 1.0-x ) );
        double log_eta = log( eta );
        val = sqrt( eta - 0.5*log_eta + ( 0.25*log_eta-0.5 )/eta );
    } else {
        unsigned int n = floor( log( erfinv_xmax_/( 1.0-x ) )/erfinv_alpha_ );
        double wl = ( erfinv_x_[n+1]-x )/( erfinv_x_[n+1]-erfinv_x_[n] );
        double wr = ( x-erfinv_x_[n] )  /( erfinv_x_[n+1]-erfinv_x_[n] );
        val = wl*erfinv_tab_[n] + wr*erfinv_tab_[n+1];
    }
    
    return val;
}






