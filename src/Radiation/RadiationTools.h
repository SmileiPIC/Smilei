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
        //! Return the value of the function h(chipa) of Niel et al.
        //! from a polynomial numerical fit at order 10
        //! Valid between chipa in 1E-3 and 1E1
        //! \param chipa particle quantum parameter
        // -----------------------------------------------------------------------------
        static double inline get_h_Niel_from_fit_order10(double chipa)
        {
            // Max relative error ~2E-4
            return exp(-3.231764974833856e-08 * pow(log(chipa),10)
                       -7.574417415366786e-07 * pow(log(chipa),9)
                       -5.437005218419013e-06 * pow(log(chipa),8)
                       -4.359062260446135e-06 * pow(log(chipa),7)
                       + 5.417842511821415e-05 * pow(log(chipa),6)
                       -1.263905701127627e-04 * pow(log(chipa),5)
                       + 9.899812622393002e-04 * pow(log(chipa),4)
                       + 1.076648497464146e-02 * pow(log(chipa),3)
                       -1.624860613422593e-01 * pow(log(chipa),2)
                       + 1.496340836237785e+00 * log(chipa)
                       -2.756744141581370e+00);
        }

        // -----------------------------------------------------------------------------
        //! Return the value of the function h(chipa) of Niel et al.
        //! from a polynomial numerical fit at order 5
        //! Valid between chipa in 1E-3 and 1E1
        //! \param chipa particle quantum parameter
        // -----------------------------------------------------------------------------
        static double inline get_h_Niel_from_fit_order5(double chipa)
        {
            // Max relative error ~0.02
            return exp(1.399937206900322e-04 * pow(log(chipa),5)
                       + 3.123718241260330e-03 * pow(log(chipa),4)
                       + 1.096559086628964e-02 * pow(log(chipa),3)
                       -1.733977278199592e-01 * pow(log(chipa),2)
                       + 1.492675770100125e+00 * log(chipa)
                       -2.748991631516466e+00 );
        }

        // -----------------------------------------------------------------------------
        //! Return the value of the function h(chipa) of Niel et al.
        //! using the numerical fit of Ridgers in
        //! Ridgers et al., ArXiv 1708.04511 (2017)
        //! \param chipa particle quantum parameter
        // -----------------------------------------------------------------------------
        static double inline get_h_Niel_from_fit_Ridgers(double chipa)
        {
            return pow(chipa,3)*1.9846415503393384*pow(1. +
                    (1. + 4.528*chipa)*log(1.+12.29*chipa) + 4.632*pow(chipa,2),-7./6.);
        }

        //! Computation of the function g of Erber using the Ridgers
        //! approximation formulae
        //! \param chipa particle quantum parameter
        //#pragma omp declare simd
        static double inline compute_g_Ridgers(double chipa)
        {
            return pow(1. + 4.8*(1.+chipa)*log(1. + 1.7*chipa)
                       + 2.44*chipa*chipa,-2./3.);
        };

        // -----------------------------------------------------------------------------
        //! Return f1(nu) = Int_nu^\infty K_{5/3}(y) dy
        //! used in computed synchrotron power spectrum
        // -----------------------------------------------------------------------------
        static double inline compute_f1_nu(double nu)
        {
            if (nu<0.1)      return 2.149528241483088*pow(nu,-0.6666666666666667) - 1.813799364234217;
            else if (nu>10)  return 1.253314137315500*pow(nu,-0.5)*exp(-nu);
            else {
                return exp(-1.042081355552157e-02 * pow(log(nu),5)
                           -5.349995695960174e-02 * pow(log(nu),4)
                           -1.570476212230771e-01 * pow(log(nu),3)
                           -4.575331390887448e-01 * pow(log(nu),2)
                           -1.687909081004528e+00 * pow(log(nu),1)
                           -4.341018460806052e-01) ;
            }
        }

        // -----------------------------------------------------------------------------
        //! Return f2(nu) = BesselK_{2/3}(nu)
        //! used in computed synchrotron power spectrum
        // -----------------------------------------------------------------------------
        static double inline compute_f2_nu(double nu)
        {
            if (nu<0.05)     return 1.074764120720013*pow(nu,-0.6666666666666667);
            else if (nu>10)  return 1.253314137315500*pow(nu,-0.5)*exp(-nu);
            else {
                return exp(-7.694562217592761e-03 * pow(log(nu),5)
                           -5.412029310872778e-02 * pow(log(nu),4)
                           -1.782660550734939e-01 * pow(log(nu),3)
                           -4.589601096726573e-01 * pow(log(nu),2)
                           -1.539212709860801e+00 * pow(log(nu),1)
                           -7.121012104149862e-01) ;
            }
        }

    private:


};
#endif
