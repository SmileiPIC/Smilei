#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

class userFunctions {
    
public:

    static double erfinv(double x);
    static double erfinv2(double x);

    //! Modified Bessel function of first and second kind
    static void modified_bessel_IK(double n, double x, double & I, double & dI, 
                                   double & K, double & dK, unsigned int maxit, double eps);

    //! Chebychev evaluation
    static double chebychev_eval(const double * c, const int m, const double x);
 
private:
    
    
};
#endif

