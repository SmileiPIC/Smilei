#include "ElectroMagn2D.h"
#include "PicParams.h"
#include "Field2D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart2D.h"

#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::ElectroMagn2D(PicParams* params, SmileiMPI* smpi)
	: ElectroMagn(params, smpi)
{
	MESSAGE( "to be implemented" );
		
} // END constructor Electromagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::~ElectroMagn2D()
{
}


// ---------------------------------------------------------------------------------------------------------------------
// Solve Poisson
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solvePoisson(SmileiMPI* smpi)
{
	MESSAGE( "to be implemented" );
	
}//END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwell(double time_dual, SmileiMPI* smpi)
{
	MESSAGE( "to be implemented" );
	
}

void ElectroMagn2D::solveMaxwellAmpere()
{
	MESSAGE( "to be implemented" );

}


void ElectroMagn2D::solveMaxwellFaraday()
{
	MESSAGE( "to be implemented" );

}


void ElectroMagn2D::boundaryConditions(double time_dual, SmileiMPI* smpi)
{
	MESSAGE( "to be implemented" );

}


// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge density and transverse currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::initRhoJ()
{
	MESSAGE( "to be implemented" );
    
}
