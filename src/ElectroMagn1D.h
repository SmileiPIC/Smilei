#ifndef ELECTROMAGN1D_H
#define ELECTROMAGN1D_H

#include "ElectroMagn.h"

class PicParams;

//! class ElectroMagn1D containing all information on the electromagnetic fields & currents for 1d3v simulations
class ElectroMagn1D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn1D
    ElectroMagn1D(PicParams* params, SmileiMPI* smpi);

    //! Destructor for ElectroMagn1D
    ~ElectroMagn1D();

    //! Oversize
    unsigned int oversize_;

    //! Method used for initializing Maxwell solver
    void solvePoisson(SmileiMPI* smpi);

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
    //! from smpi is west
    const bool isWestern;
    
    //! from smpi is east
    const bool isEastern;
};

#endif
