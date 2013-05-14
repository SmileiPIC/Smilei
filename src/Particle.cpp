#include "Particle.h"
#include "PicParams.h"

#include <iostream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Particle (takes the dimension of the particle as input parameter)
// ---------------------------------------------------------------------------------------------------------------------
Particle::Particle(int nDim)
{
	buf  = new double[nDim+3+1+nDim];

	Position = new (buf)        double[nDim];
	Momentum = new (buf+nDim)   double[3];
	Weight	 = new (buf+nDim+3) double[1];

	Position_old = new (buf+nDim+3+1) double[nDim];

	Position[0]     = 0.;
	Position_old[0] = 0.;
	for (unsigned int i=0 ; i<3 ; i++ ) Momentum[i] = 0.;
	Weight[0] = 0.;

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Particle
// ---------------------------------------------------------------------------------------------------------------------
Particle::~Particle()
{
	delete [] buf;
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used to print the Particle properties
// ---------------------------------------------------------------------------------------------------------------------
void Particle::Print(PicParams* params)
{
	for (unsigned int i=0 ; i<params->nDim_field ; i++ ) cout << Position[i] << " ";
	for (unsigned int i=0 ; i<3 ; i++ )                  cout << Momentum[i] << " ";
	cout <<  Weight << " " << endl;
}
