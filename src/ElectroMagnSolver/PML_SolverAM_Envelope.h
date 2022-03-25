#ifndef PML_SOLVERAM_ENVELOPE_H
#define PML_SOLVERAM_ENVELOPE_H

#include "SolverAM.h"
class ElectroMagn;
class LaserEnvelope;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PML_SolverAM_Envelope : public SolverAM
{

public:
    PML_SolverAM_Envelope( Params &params );
    virtual ~PML_SolverAM_Envelope();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) ;

    void setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch );

    void compute_A_from_G( LaserEnvelope *envelope, int iDim, int min_or_max, int solver_min, int solver_max );

protected:
    double alpha_l_max ;
    double alpha_cl ;
    double sigma_l_max ;
    double kappa_l_max ;
    double power_pml_kappa_l ;
    double power_pml_sigma_l ;
    double power_pml_alpha_l ;
    double alpha_r_max ;
    double alpha_cr ;
    double sigma_r_max ;
    double kappa_r_max ;
    double power_pml_kappa_r ;
    double power_pml_sigma_r ;
    double power_pml_alpha_r ;

    std::vector<double> kappa_l_p;
    std::vector<double> sigma_l_p;
    std::vector<double> alpha_l_p;
    std::vector<double> kappa_prime_l_p;
    std::vector<double> sigma_prime_l_p;
    std::vector<double> alpha_prime_l_p;
    std::vector<double> kappa_r_p;
    std::vector<double> sigma_r_p;
    std::vector<double> alpha_r_p;
    std::vector<double> kappa_prime_r_p;
    std::vector<double> sigma_prime_r_p;
    std::vector<double> alpha_prime_r_p;
    std::vector<double> integrate_kappa_r_p ;
    std::vector<double> integrate_sigma_r_p ;
    std::vector<double> integrate_alpha_r_p ;

    //double xmin;
    //double ymin;

    double length_l_pml;
    double length_r_pml;
    double length_l_pml_lmin;
    double length_l_pml_lmax;

    int j_glob_pml;
    double rmax;
    double r0;

};//END class

#endif
