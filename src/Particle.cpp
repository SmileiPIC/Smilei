#include "Particle.h"
#include "PicParams.h"

#include <iostream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Particle (takes the dimension of the particle as input parameter)
// ---------------------------------------------------------------------------------------------------------------------
Particle::Particle(int nDim)
{
	Position     = new double[nDim];
    Position_old = new double[nDim];
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Particle
// ---------------------------------------------------------------------------------------------------------------------
Particle::~Particle()
{
	delete [] Position;
    delete [] Position_old;
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
