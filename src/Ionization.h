#ifndef IONIZATION_H
#define IONIZATION_H

#include "Tools.h"
#include "PicParams.h"
#include "Field.h"


class Particle;

class IonizationData {
	std::vector<double> Potential;
	std::vector<int> azimuthal_quantum_number;
};

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
//	std::map<int,IonizationData> IonizationDataValues;
	
private:
	
};

#endif
