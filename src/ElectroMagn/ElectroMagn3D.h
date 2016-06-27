#ifndef ELECTROMAGN3D_H
#define ELECTROMAGN3D_H

#include "ElectroMagn.h"
#include "Field.h"
#include "Field3D.h"

class Params;

//! class ElectroMagn3D containing all information on the electromagnetic fields & currents for 2d3v simulations
class ElectroMagn3D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn3D
    ElectroMagn3D(Params &params, std::vector<Species*>& vecSpecies, Patch* patch);

    //! Destructor for ElectroMagn3D
    ~ElectroMagn3D();
    
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

    double getEx_WestNorth() { return 0.; }
    double getEy_WestNorth() { return 0.; }
    double getEx_EastSouth() { return 0.; }
    double getEy_EastSouth() { return 0.; }

#ifdef _PATCH3D_TODO
#endif

    // --------------------------------------
    //  --------- PATCH IN PROGRESS ---------
    // --------------------------------------
    
    //! Method used to solve Maxwell-Ampere equation
    void solveMaxwellAmpere();

    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields();

    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields();

    //! Method used to reset/increment the averaged fields
    void incrementAvgFields(unsigned int time_step);
    
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

    //! Number of nodes on the primal grid in the z-direction
    unsigned int nz_p;

    //! Number of nodes on the dual grid in the z-direction
    unsigned int nz_d;

    //! Spatial step dx for 3D3V cartesian simulations
    double dx;

    //! Spatial step dy for 3D3V cartesian simulations
    double dy;

    //! Spatial step dz for 3D3V cartesian simulations
    double dz;

    //! Ratio of the time-step by the spatial-step dt/dx for 3D3V cartesian simulations
    double dt_ov_dx;

    //! Ratio of the time-step by the spatial-step dt/dy for 3D3V cartesian simulations
    double dt_ov_dy;

    //! Ratio of the time-step by the spatial-step dt/dz for 3D3V cartesian simulations
    double dt_ov_dz;

    //! Ratio of the spatial-step by the time-step dx/dt for 3D3V cartesian simulations
    double dx_ov_dt;

    //! Ratio of the spatial-step by the time-step dy/dt for 3D3V cartesian simulations
    double dy_ov_dt;

    //! Ratio of the spatial-step by the time-step dz/dt for 3D3V cartesian simulations
    double dz_ov_dt;

    //! compute Poynting on borders
    void computePoynting();

    //! Method used to impose external fields
    void applyExternalField(Field*, Profile*, Patch*);
    
    void initAntennas(Patch* patch);
        
private:
    
    //! from smpi is west
    const bool isWestern;
    
    //! from smpi is east
    const bool isEastern;
    
    //! from smpi is north
    const bool isSouthern;
    
    //! from smpi is south
    const bool isNorthern;

    //! from smpi is top
    const bool isTop;
    
    //! from smpi is bottom
    const bool isBottom;
};

#endif
