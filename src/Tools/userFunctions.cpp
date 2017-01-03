#include <limits>
#include <math.h>
#include "userFunctions.h"


//! inverse error function is taken from NIST
double userFunctions::erfinv (double x)
{
    if (x < -1. || x > 1.)
        return std::numeric_limits<double>::quiet_NaN();
    
    if (x == 0)
        return 0;
    
    int sign_x=(x > 0? 1 : -1);
    
    double r;
    if (x <= 0.686) {
        double x2 = x * x;
        r  = x * (((-0.140543331 * x2 + 0.914624893) * x2 + -1.645349621) * x2 + 0.886226899);
        r /= (((0.012229801 * x2 + -0.329097515) * x2 + 1.442710462) * x2 + -2.118377725) * x2 + 1.;
    } else {
        double y = sqrt (-log((1. - x) / 2.));
        r  = (((1.641345311 * y + 3.429567803) * y + -1.62490649) * y + -1.970840454);
        r /= (1.637067800 * y + 3.543889200) * y + 1.;
    }
    
    r *= (double)sign_x;
    x *= (double)sign_x;
    
    r -= (erf(r) - x) / (2. / sqrt (M_PI) * exp(-r*r));
    
    return r;
}

//! inverse error function is taken from M.B. Giles. 'Approximating the erfinv function'. In GPU Computing Gems, volume 2, Morgan Kaufmann, 2011.
double userFunctions::erfinv2 (double x)
{
    double w, p; 
    w = -log((1.0-x)*(1.0+x));

    if ( w < 5.000000 ) {
        w = w - 2.500000; 
        p = +2.81022636000e-08      ;
        p = +3.43273939000e-07 + p*w;
        p = -3.52338770000e-06 + p*w;
        p = -4.39150654000e-06 + p*w;
        p = +0.00021858087e+00 + p*w;
        p = -0.00125372503e+00 + p*w;
        p = -0.00417768164e+00 + p*w;
        p = +0.24664072700e+00 + p*w;
        p = +1.50140941000e+00 + p*w;
    } else {
        w = sqrt(w) - 3.000000; 
        p = -0.000200214257      ;
        p = +0.000100950558 + p*w;
        p = +0.001349343220 + p*w;
        p = -0.003673428440 + p*w;
        p = +0.005739507730 + p*w;
        p = -0.007622461300 + p*w;
        p = +0.009438870470 + p*w;
        p = +1.001674060000 + p*w;
        p = +2.832976820000 + p*w;
    }  
    return p*x;
}



