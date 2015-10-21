#ifndef ELECTROMAGN2D_H
#define ELECTROMAGN2D_H

#include "ElectroMagn.h"
#include "Field.h"
#include "Field2D.h"
class PicParams;

//! class ElectroMagn2D containing all information on the electromagnetic fields & currents for 2d3v simulations
class ElectroMagn2D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn2D
    ElectroMagn2D(PicParams &params, LaserParams &laser_params, Patch* patch);

    //! Destructor for ElectroMagn2D
    ~ElectroMagn2D();
    
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

    double getEx_West() { return 0.; }
    double getEx_East() { return 0.; }

    double getEx_WestNorth() { return (*Ex_)(0,ny_p-1); }
    double getEy_WestNorth() { return (*Ey_)(0,ny_d-1); }
    double getEx_EastSouth() { return (*Ex_)(nx_d-1,0); }
    double getEy_EastSouth() { return (*Ey_)(nx_p-1,0); }

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
    
    //! Method used to initialize the total charge densities and currents
    void restartRhoJ();
    //! Method used to initialize the total charge densities and currents of species
    void restartRhoJs();

    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ();
    void addToGlobalRho(int ispec, unsigned int clrw);
    void computeTotalRhoJs(unsigned int clrw);
    //! Method used to gather species densities and currents on a single array
    void synchronizePatch(unsigned int clrw);
    void finalizePatch(unsigned int clrw);

    //! \todo Create properties the laser time-profile (MG & TV)

    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;

    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;

    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;

    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;

    //! Spatial step dx for 2D3V cartesian simulations
    double dx;

    //! Spatial step dy for 2D3V cartesian simulations
    double dy;

    //! Ratio of the time-step by the spatial-step dt/dx for 2D3V cartesian simulations
    double dt_ov_dx;

    //! Ratio of the time-step by the spatial-step dt/dy for 2D3V cartesian simulations
    double dt_ov_dy;

    //! Ratio of the spatial-step by the time-step dx/dt for 2D3V cartesian simulations
    double dx_ov_dt;

    //! Ratio of the spatial-step by the time-step dy/dt for 2D3V cartesian simulations
    double dy_ov_dt;

    //! compute Poynting on borders
    void computePoynting();

private:
    
    //! from smpi is west
    const bool isWestern;
    
    //! from smpi is east
    const bool isEastern;
    
    //! from smpi is north
    const bool isSouthern;
    
    //! from smpi is south
    const bool isNorthern;
};

#endif
