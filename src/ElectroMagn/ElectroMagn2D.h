#ifndef ELECTROMAGN2D_H
#define ELECTROMAGN2D_H

#include "ElectroMagn.h"

class PicParams;

//! class ElectroMagn2D containing all information on the electromagnetic fields & currents for 2d3v simulations
class ElectroMagn2D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn2D
    ElectroMagn2D(PicParams &params, LaserParams &laser_params, SmileiMPI* smpi);

    //! Destructor for ElectroMagn2D
    ~ElectroMagn2D();

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
    
    //! Method used to initialize the total charge densities and currents
    void restartRhoJ();
    //! Method used to initialize the total charge densities and currents of species
    void restartRhoJs(int ispec, bool currents);

    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ();
    //! Method used to gather species densities and currents on a single array
    void sumtwins();

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
