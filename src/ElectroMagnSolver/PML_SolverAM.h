#ifndef PML_SOLVERAM_H
#define PML_SOLVERAM_H

#include "SolverAM.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PML_SolverAM : public SolverAM
{

public:
    PML_SolverAM( Params &params );
    virtual ~PML_SolverAM();

    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );

    void setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml_domain, int startpml, int* ncells_pml_lmin, int* ncells_pml_lmax, Patch* patch );

    void compute_E_from_D( ElectroMagn *fields, int iDim, int min_or_max, unsigned int solvermin, unsigned int solvermax );
    void compute_H_from_B( ElectroMagn *fields, int iDim, int min_or_max, unsigned int solvermin, unsigned int solvermax );

protected:
    double sigma_r_max;
    double kappa_r_max;
    double power_pml_r;
    double sigma_l_max;
    double kappa_l_max;
    double power_pml_l;

    std::vector<double> kappa_r_p;
    std::vector<double> sigma_r_p;
    std::vector<double> kappa_r_d;
    std::vector<double> sigma_r_d;
    std::vector<double> integrate_kappa_r_p;
    std::vector<double> integrate_sigma_r_p;
    std::vector<double> integrate_kappa_r_d;
    std::vector<double> integrate_sigma_r_d;
    std::vector<double> kappa_l_p;
    std::vector<double> sigma_l_p;
    std::vector<double> kappa_l_d;
    std::vector<double> sigma_l_d;

    std::vector<double> c1_p_lfield;
    std::vector<double> c2_p_lfield;
    std::vector<double> c3_p_lfield;
    std::vector<double> c4_p_lfield;
    std::vector<double> c5_p_lfield;
    std::vector<double> c6_p_lfield;
    std::vector<double> c1_d_lfield;
    std::vector<double> c2_d_lfield;
    std::vector<double> c3_d_lfield;
    std::vector<double> c4_d_lfield;
    std::vector<double> c5_d_lfield;
    std::vector<double> c6_d_lfield;
    std::vector<double> c1_p_rfield;
    std::vector<double> c2_p_rfield;
    std::vector<double> c3_p_rfield;
    std::vector<double> c4_p_rfield;
    std::vector<double> c5_p_rfield;
    std::vector<double> c6_p_rfield;
    std::vector<double> c1_d_rfield;
    std::vector<double> c2_d_rfield;
    std::vector<double> c3_d_rfield;
    std::vector<double> c4_d_rfield;
    std::vector<double> c5_d_rfield;
    std::vector<double> c6_d_rfield;
    std::vector<double> c1_p_tfield;
    std::vector<double> c2_p_tfield;
    std::vector<double> c3_p_tfield;
    std::vector<double> c4_p_tfield;
    std::vector<double> c5_p_tfield;
    std::vector<double> c6_p_tfield;
    std::vector<double> c1_d_tfield;
    std::vector<double> c2_d_tfield;
    std::vector<double> c3_d_tfield;
    std::vector<double> c4_d_tfield;
    std::vector<double> c5_d_tfield;
    std::vector<double> c6_d_tfield;

    // double lmax;
    // double l0;
    int j_glob_pml;
    double rmax;
    double r0;
    double length_r_pml;
    double length_l_pml;
    double length_l_pml_lmin;
    double length_l_pml_lmax;

    bool isMin;
    bool isMax;

    std::complex<double> Dl_pml_old;
    std::complex<double> Dr_pml_old;
    std::complex<double> Dt_pml_old;
    std::complex<double> Bl_pml_old;
    std::complex<double> Br_pml_old;
    std::complex<double> Bt_pml_old;

};//END class

#endif
