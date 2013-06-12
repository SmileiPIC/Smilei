#ifndef ELECTROMAGN2D_H
#define ELECTROMAGN2D_H

#include "ElectroMagn.h"

class PicParams;

//! class ElectroMagn2D containing all information on the electromagnetic fields & currents for 2d3v simulations
class ElectroMagn2D : public ElectroMagn
{
public:
	//! Constructor for ElectroMagn2D
	ElectroMagn2D(PicParams* params, SmileiMPI* smpi);
	//! Destructor for ElectroMagn2D
	~ElectroMagn2D();

	//! Constant used for the Silver-Mueller boundary conditions
	double A_;
	//! Constant used for the Silver-Mueller boundary conditions
	double B_;
	//! Constant used for the Silver-Mueller boundary conditions
	double C_;
	    
	//! Method used for initializing Maxwell solver
	void solvePoisson(SmileiMPI* smpi);
    
	//! Method used to initialize the total charge densities and currents
	void initRhoJ();
    
	//! Method used to solve Maxwell's equations
	void solveMaxwell(double time_dual, SmileiMPI* smpi);
	void solveMaxwellAmpere();
	void solveMaxwellFaraday();
	void boundaryConditions(double time_dual, SmileiMPI* smpi);

	//! \todo Create properties the laser time-profile (MG & TV)

	//! Number of nodes on the primal grid
	unsigned int nx_p;
    
	//! Number of nodes on the dual grid
	unsigned int nx_d;
    
	//! Spatial step dx for 2D3v cartesian simulations
	double dx;
	double dy;
    
	//! Ratio of the time-step by the spatial-step dt/dx for 2d3v cartesian simulations
	double dt_ov_dx;
	double dt_ov_dy;
    
	//! Ratio of the spatial-step by the time-step dx/dt for 2d3v cartesian simulations
	double dx_ov_dt;
	double dy_ov_dt;
    
private:
};

#endif
