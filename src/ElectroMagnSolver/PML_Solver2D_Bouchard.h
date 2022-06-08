#ifndef PML_SOLVER2D_BOUCHARD_H
#define PML_SOLVER2D_BOUCHARD_H

#include "Solver2D.h"
class ElectroMagn;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class PML_Solver2D_Bouchard : public Solver2D
{

public:
    PML_Solver2D_Bouchard( Params &params );
    virtual ~PML_Solver2D_Bouchard();
    
    //! Overloading of () operator
    virtual void operator()( ElectroMagn *fields );

    void setDomainSizeAndCoefficients( int iDim, int min_or_max, int ncells_pml_domain, int startpml, int* ncells_pml_min, int* ncells_pml_max, Patch* patch );

    void compute_E_from_D( ElectroMagn *fields, int iDim, int min_or_max, int solver_min, int solver_max );
    void compute_H_from_B( ElectroMagn *fields, int iDim, int min_or_max, int solver_min, int solver_max );

protected:

    // Parameters for the Maxwell-Faraday solver
    
    double delta_x;
    double delta_y;
    double beta_xy;
    double beta_yx;
    double alpha_x;
    double alpha_y;
    double Ax;
    double Ay;
    double Bx;
    double By;
    double Dx;
    double Dy;

    double sigma_x_max;
    double kappa_x_max;
    double sigma_power_pml_x;
    double kappa_power_pml_x;
    double sigma_y_max;
    double kappa_y_max;
    double sigma_power_pml_y;
    double kappa_power_pml_y;   
 
    std::vector<double> kappa_x_p;
    std::vector<double> sigma_x_p;
    std::vector<double> kappa_x_d;
    std::vector<double> sigma_x_d;
    std::vector<double> kappa_y_p;
    std::vector<double> sigma_y_p;
    std::vector<double> kappa_y_d;
    std::vector<double> sigma_y_d;
    std::vector<double> kappa_z_p;
    std::vector<double> sigma_z_p;
    std::vector<double> kappa_z_d;
    std::vector<double> sigma_z_d;

    std::vector<double> c1_p_xfield;
    std::vector<double> c2_p_xfield;
    std::vector<double> c3_p_xfield;
    std::vector<double> c4_p_xfield;
    std::vector<double> c5_p_xfield;
    std::vector<double> c6_p_xfield;
    std::vector<double> c1_d_xfield;
    std::vector<double> c2_d_xfield;
    std::vector<double> c3_d_xfield;
    std::vector<double> c4_d_xfield;
    std::vector<double> c5_d_xfield;
    std::vector<double> c6_d_xfield;
    std::vector<double> c1_p_yfield;
    std::vector<double> c2_p_yfield;
    std::vector<double> c3_p_yfield;
    std::vector<double> c4_p_yfield;
    std::vector<double> c5_p_yfield;
    std::vector<double> c6_p_yfield;
    std::vector<double> c1_d_yfield;
    std::vector<double> c2_d_yfield;
    std::vector<double> c3_d_yfield;
    std::vector<double> c4_d_yfield;
    std::vector<double> c5_d_yfield;
    std::vector<double> c6_d_yfield;
    std::vector<double> c1_p_zfield;
    std::vector<double> c2_p_zfield;
    std::vector<double> c3_p_zfield;
    std::vector<double> c4_p_zfield;
    std::vector<double> c5_p_zfield;
    std::vector<double> c6_p_zfield;
    std::vector<double> c1_d_zfield;
    std::vector<double> c2_d_zfield;
    std::vector<double> c3_d_zfield;
    std::vector<double> c4_d_zfield;
    std::vector<double> c5_d_zfield;
    std::vector<double> c6_d_zfield;    

    //double xmin;
    //double ymin;
    
    double length_x_pml;
    double length_y_pml;
    double length_x_pml_xmin;
    double length_x_pml_xmax;

    bool isMin ;
    bool isMax ;
    double Dx_pml_old ;
    double Dy_pml_old ;
    double Dz_pml_old ;
    double Bx_pml_old ;
    double By_pml_old ;
    double Bz_pml_old ;
    
};//END class

#endif

