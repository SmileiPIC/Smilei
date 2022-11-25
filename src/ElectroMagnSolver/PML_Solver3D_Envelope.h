#ifndef PML_SOLVER3D_ENVELOPE_H
#define PML_SOLVER3D_ENVELOPE_H

#include "Solver3D.h"
class ElectroMagn;
class LaserEnvelope;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PML_Solver3D_Envelope : public Solver3D
{

public:
    PML_Solver3D_Envelope( Params &params );
    virtual ~PML_Solver3D_Envelope();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields ) ;

    void setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch );

    void compute_A_from_G( LaserEnvelope *envelope, int iDim, int min_or_max, unsigned int solver_min, unsigned int solver_max );

protected:
    double alpha_x_max ;
    double sigma_x_max ;
    double kappa_x_max ;
    double power_pml_kappa_x ;
    double power_pml_sigma_x ;
    double power_pml_alpha_x ;
    double alpha_cx ;
    double alpha_y_max ;
    double sigma_y_max ;
    double kappa_y_max ;
    double power_pml_kappa_y ;
    double power_pml_sigma_y ;
    double power_pml_alpha_y ;
    double alpha_cy ;
    double alpha_z_max ;
    double sigma_z_max ;
    double kappa_z_max ;
    double power_pml_kappa_z ;
    double power_pml_sigma_z ;
    double power_pml_alpha_z ;
    double alpha_cz ;

    std::vector<double> kappa_x_p;
    std::vector<double> sigma_x_p;
    std::vector<double> alpha_x_p;
    std::vector<double> kappa_prime_x_p;
    std::vector<double> sigma_prime_x_p;
    std::vector<double> alpha_prime_x_p;
    std::vector<double> kappa_y_p;
    std::vector<double> sigma_y_p;
    std::vector<double> alpha_y_p;
    std::vector<double> kappa_prime_y_p;
    std::vector<double> sigma_prime_y_p;
    std::vector<double> alpha_prime_y_p;
    std::vector<double> kappa_z_p;
    std::vector<double> sigma_z_p;
    std::vector<double> alpha_z_p;
    std::vector<double> kappa_prime_z_p;
    std::vector<double> sigma_prime_z_p;
    std::vector<double> alpha_prime_z_p;

    //double xmin;
    //double ymin;

    double length_x_pml;
    double length_y_pml;
    double length_z_pml;
    double length_x_pml_xmin;
    double length_x_pml_xmax;
    double length_y_pml_ymin;
    double length_y_pml_ymax;

};//END class

#endif
