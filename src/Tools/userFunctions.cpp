#include <cmath>
#include "userFunctions.h"


//! inverse error function is taken from NIST
double userFunctions::erfinv (double x)
{
    if (x < -1 || x > 1)
        return NAN;
    
    if (x == 0)
        return 0;
    
    int  sign_x;
    if (x > 0) {
        sign_x = 1;
    } else {
        sign_x = -1;
        x = -x;
    }
    
    double r;
    if (x <= 0.686) {
        double x2 = x * x;
        r = x * (((-0.140543331 * x2 + 0.914624893) * x2 + -1.645349621) * x2 + 0.886226899);
        r /= (((0.012229801 * x2 + -0.329097515) * x2 + 1.442710462) * x2 + -2.118377725) * x2 + 1;
    } else {
        double y = sqrt (-log((1 - x) / 2));
        r = (((1.641345311 * y + 3.429567803) * y + -1.62490649) * y + -1.970840454);
        r /= ((1.637067800 * y + 3.543889200) * y + 1);
    }
    
    r *= sign_x;
    x *= sign_x;
    
    r -= (erf(r) - x) / (2 / sqrt (M_PI) * exp(-r*r));
    
    return r;
}

