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
    
    //! Overloading of () operator
	virtual void operator() (Particle* part, LocalFields Epart, LocalFields Jion) = 0;
	
    std::vector<Particle*> new_electrons;
	
protected:
	std::vector<double> Potential;
	std::vector<double> Azimuthal_quantum_number;
    
    double eV_to_au;
    double EC_to_au;
    double au_to_w0;
	
    double wavelength_SI;
	double dt;
	unsigned int nDim_field;
	unsigned int nDim_particle;
	unsigned int atomic_number_;
    unsigned int ionized_species_mass;
    
private:


};

#endif
