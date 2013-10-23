#ifndef IONIZATION_H
#define IONIZATION_H

#include "Tools.h"
#include "PicParams.h"
#include "Field.h"
#include <map>

class Particle;

//! Class Ionization: generic class allowing to define Ionization physics
class Ionization
{
    
public:
	//! Constructor for Ionization
	Ionization(PicParams *params, int ispec);
	virtual ~Ionization();
    
    //! Overloading of () operator
	virtual void operator() (Particle* part, LocalFields Epart) = 0;
	
	
protected:
	std::vector<double> Potential;
	std::vector<int> azimuthal_quantum_number;
	
private:
	double dt, dts2;
	// mass_ and charge_ relative to Species but used in the particle pusher
	double mass_;
	double charge_;
	double charge_over_mass_;
	
	int atomic_number_;
	int nDim_;

};

#endif
