#ifndef ELECTROMAGN1D_H
#define ELECTROMAGN1D_H

#include "ElectroMagn.h"

class PicParams;

//! class ElectroMagn1D containing all information on the electromagnetic fields & currents for 1d3v simulations
class ElectroMagn1D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn1D
    ElectroMagn1D(PicParams &params,  LaserParams &laser_params, Patch* patch);

    //! Destructor for ElectroMagn1D
    ~ElectroMagn1D();

    //! Oversize
    unsigned int oversize_;

    // --------------------------------------
    //  --------- PATCH IN PROGRESS ---------
    // --------------------------------------
    void initPoisson(Patch *patch);
    double compute_r();
    void compute_Ap(Patch *patch);
    //Access to Ap
    double compute_pAp();
    void update_pand_r(double r_dot_r, double p_dot_Ap);
    void update_p(double rnew_dot_rnew, double r_dot_r);
    void initE(Patch *patch);
    void centeringE( std::vector<double> E_Add );

    double getEx_West() { return (*Ex_)(index_bc_min[0]);}//(*Ex_)     (0); }
    double getEx_East() { return (*Ex_)(index_bc_max[0]);}//(*Ex_)(nx_d-1); }

    double getEx_WestNorth() { return 0.; }
    double getEy_WestNorth() { return 0.; }
    double getEx_EastSouth() { return 0.; }
    double getEy_EastSouth() { return 0.; }

    // --------------------------------------
    //  --------- PATCH IN PROGRESS ---------
    // --------------------------------------
    
    //! Method used to solve Maxwell-Ampere equation
    void solveMaxwellAmpere();

    //! Method used to solve Maxwell-Faraday equation
    void solveMaxwellFaraday();

    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields();

    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields();
    
    //! Method used to reset/increment the averaged fields
    void incrementAvgFields(unsigned int time_step, unsigned int ntime_step_avg);

    //! Method used to restart the total charge densities and currents
    void restartRhoJ();
    //! Method used to restart the total charge densities and currents
    void restartRhoJs();

    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ();

    //! \todo Create properties the laser time-profile (MG & TV)
 
    //! Number of nodes on the primal grid
    unsigned int nx_p;

     //! Number of nodes on the dual grid
    unsigned int nx_d;

    //! Spatial step dx for 1d3v cartesian simulations
    double dx;

    //! Ratio of the time-step by the spatial-step dt/dx for 1d3v cartesian simulations
    double dt_ov_dx;

    //! Ratio of the spatial-step by the time-step dx/dt for 1d3v cartesian simulations
    double dx_ov_dt;

    //! compute Poynting on borders
    void computePoynting();

private:
    //! from patch is west
    const bool isWestern;
    
    //! from patch is east
    const bool isEastern;
};

#endif
