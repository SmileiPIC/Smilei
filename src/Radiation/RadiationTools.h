// ----------------------------------------------------------------------------
//! \file RadiationTools.h
//
//! \brief This class contains tools
//! for the Nonlinear Inverse Compton Scattering
//
// ----------------------------------------------------------------------------

#ifndef RADIATIONTOOLS_H
#define RADIATIONTOOLS_H

#include <cmath>

class RadiationTools {

    public:

        // -----------------------------------------------------------------------------
        //! Return the value of the function h(particle_chi) of Niel et al.
        //! from a polynomial numerical fit at order 10
        //! Valid between particle_chi in 1E-3 and 1E1
        //! \param particle_chi particle quantum parameter
        // -----------------------------------------------------------------------------
        static inline double getHNielFitOrder10(double particle_chi)
        {
            // Max relative error ~2E-4
            double logchi1 = log(particle_chi);
            double logchi2 = logchi1 * logchi1;
            double logchi3 = logchi2 * logchi1;
            double logchi4 = logchi3 * logchi1;
            double logchi5 = logchi4 * logchi1;
            return exp(-3.231764974833856e-08 * logchi5*logchi5
                       -7.574417415366786e-07 * logchi5*logchi4
                       -5.437005218419013e-06 * logchi5*logchi3
                       -4.359062260446135e-06 * logchi5*logchi2
                       +5.417842511821415e-05 * logchi5*logchi1
                       -1.263905701127627e-04 * logchi5
                       +9.899812622393002e-04 * logchi4
                       +1.076648497464146e-02 * logchi3
                       -1.624860613422593e-01 * logchi2
                       +1.496340836237785e+00 * logchi1
                       -2.756744141581370e+00);
        }

        // -----------------------------------------------------------------------------
        //! Return the value of the function h(particle_chi) of Niel et al.
        //! from a polynomial numerical fit at order 5
        //! Valid between particle_chi in 1E-3 and 1E1
        //! \param particle_chi particle quantum parameter
        // -----------------------------------------------------------------------------
        static double inline getHNielFitOrder5(double particle_chi)
        {

            double logchi1 = log(particle_chi);
            double logchi2 = logchi1 * logchi1;
            double logchi3 = logchi2 * logchi1;
            // Max relative error ~0.02
            return exp(+1.399937206900322e-04 * logchi3*logchi2
                       +3.123718241260330e-03 * logchi3*logchi1
                       +1.096559086628964e-02 * logchi3
                       -1.733977278199592e-01 * logchi2
                       +1.492675770100125e+00 * logchi1
                       -2.748991631516466e+00 );
        }

        // -----------------------------------------------------------------------------
        //! Return the value of the function h(particle_chi) of Niel et al.
        //! using the numerical fit of Ridgers in
        //! Ridgers et al., ArXiv 1708.04511 (2017)
        //! \param particle_chi particle quantum parameter
        // -----------------------------------------------------------------------------
        static double inline getHNielFitRidgers(double particle_chi)
        {
            double chi2 = particle_chi * particle_chi;
            double chi3 = chi2 * particle_chi;
            return chi3*1.9846415503393384
                *std::pow(
                    1.0 + (1. + 4.528*particle_chi)*log(1.+12.29*particle_chi) + 4.632*chi2
                    ,-7./6.
                );
        }

        //! Computation of the function g of Erber using the Ridgers
        //! approximation formulae
        //! \param particle_chi particle quantum parameter
        //#pragma omp declare simd
        static double inline computeGRidgers(double particle_chi)
        {
            return std::pow(1. + 4.8*(1.0+particle_chi)*log(1. + 1.7*particle_chi)
                       + 2.44*particle_chi*particle_chi,-2./3.);
        };

        // -----------------------------------------------------------------------------
        //! Return f1(nu) = Int_nu^\infty K_{5/3}(y) dy
        //! used in computed synchrotron power spectrum
        // -----------------------------------------------------------------------------
        static double inline computeF1Nu(double nu)
        {
            if (nu<0.1)      return 2.149528241483088*std::pow(nu,-0.6666666666666667) - 1.813799364234217;
            else if (nu>10)  return 1.253314137315500*std::pow(nu,-0.5)*exp(-nu);
            else {
                double lognu = log(nu);
                double lognu_power_n = lognu;
                double f = -4.341018460806052e-01 - 1.687909081004528e+00 * lognu_power_n;
                lognu_power_n *= lognu;
                f -= 4.575331390887448e-01 * lognu_power_n; // n=2

                lognu_power_n *= lognu;
                f -= 1.570476212230771e-01 * lognu_power_n; // n=3

                lognu_power_n *= lognu;
                f -= 5.349995695960174e-02 * lognu_power_n; // n=4

                lognu_power_n *= lognu;
                f -= 1.042081355552157e-02 * lognu_power_n; // n=5

                return exp(f);

                /*return exp(-1.042081355552157e-02 * pow(lognu,5)
                           -5.349995695960174e-02 * pow(lognu,4)
                           -1.570476212230771e-01 * pow(lognu,3)
                           -4.575331390887448e-01 * pow(lognu,2)
                           -1.687909081004528e+00 * lognu
                           -4.341018460806052e-01) ;*/
            }
        }

        // -----------------------------------------------------------------------------
        //! Return f2(nu) = BesselK_{2/3}(nu)
        //! used in computed synchrotron power spectrum
        // -----------------------------------------------------------------------------
        static double inline computeF2Nu(double nu)
        {
            if (nu<0.05)     return 1.074764120720013*std::pow(nu,-0.6666666666666667);
            else if (nu>10)  return 1.253314137315500*std::pow(nu,-0.5)*exp(-nu);
            else {
                double lognu = log(nu);
                double lognu_power_n = lognu;
                double f = -7.121012104149862e-01 - 1.539212709860801e+00 * lognu_power_n;
                lognu_power_n *= lognu;
                f -= 4.589601096726573e-01 * lognu_power_n; //n=2

                lognu_power_n *= lognu;
                f -= 1.782660550734939e-01 * lognu_power_n; //n=3

                lognu_power_n *= lognu;
                f -= 5.412029310872778e-02 * lognu_power_n; //n=4

                lognu_power_n *= lognu;
                f -= 7.694562217592761e-03 * lognu_power_n; //n=5

                return exp(f);

                /*return exp(-7.694562217592761e-03 * pow(lognu,5)
                           -5.412029310872778e-02 * pow(lognu,4)
                           -1.782660550734939e-01 * pow(lognu,3)
                           -4.589601096726573e-01 * pow(lognu,2)
                           -1.539212709860801e+00 * pow(lognu,1)
                           -7.121012104149862e-01) ;*/
            }
        }

        // -----------------------------------------------------------------------------
        //! Return f1(nu) + cst * f2(nu)
        //! = Int_nu^\infty K_{5/3}(y) dy + cst * BesselK_{2/3}(nu)
        //! used in computed synchrotron power spectrum
        // -----------------------------------------------------------------------------
        static double inline computeBesselPartsRadiatedPower(double nu, double cst)
        {
            double f1, f2;
            if (nu<0.1)
            {
                f2 = 1.074764120720013 / cbrt(nu*nu);
                f1 = 2*f2 - 1.813799364234217;
                return f1 + cst*f2;
            }
            else if (nu>10)
            {
                return (1.+cst)*1.253314137315500*exp(-nu)/sqrt(nu);
            }
            else
            {
                double lognu = log(nu);
                double lognu_power_n = lognu;

                f1 = - 4.364684279797524e-01;
                f2 = - 7.121012104149862e-01;
                f1 -= 1.670543589881836e+00 * lognu_power_n; //n=1
                f2 -= 1.539212709860801e+00 * lognu_power_n; //n=1

                lognu_power_n *= lognu;
                f1 -= 4.533108925728350e-01 * lognu_power_n; //n=2
                f2 -= 4.589601096726573e-01 * lognu_power_n; //n=2

                lognu_power_n *= lognu;
                f1 -= 1.723519212869859e-01 * lognu_power_n; //n=3
                f2 -= 1.782660550734939e-01 * lognu_power_n; //n=3

                lognu_power_n *= lognu;
                f1 -= 5.431864123685266e-02 * lognu_power_n; //n=4
                f2 -= 5.412029310872778e-02 * lognu_power_n; //n=4

                lognu_power_n *= lognu;
                f1 -= 7.892740572869308e-03 * lognu_power_n; //n=5
                f2 -= 7.694562217592761e-03 * lognu_power_n; //n=5


                return exp(f1)+cst*exp(f2);
            }
        }

    private:


};
#endif
