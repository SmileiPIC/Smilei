#include "ElectroMagn.h"
#include "ElectroMagn1D.h"
#include "PicParams.h"
#include "Species.h"
#include "Projector.h"
#include "Laser.h"
#include "Field.h"
#include <limits>

#include <iostream>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn(PicParams* params)
{

    // check for laser conditions
	laser_.resize(params->n_laser);
    
	for (unsigned int i=0; i<laser_.size(); i++)
    {
		DEBUG(5,"Initializing Laser "<<i);
		laser_[i] = new Laser(params->laser_param[i]);
	}
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::~ElectroMagn()
{
	delete Ex_;
	delete Ey_;
	delete Ez_;
	delete Bx_;
	delete By_;
	delete Bz_;
	delete Bx_m;
	delete By_m;
	delete Bz_m;
	delete Jx_;
	delete Jy_;
	delete Jz_;
	delete rho_;
	delete rho_o;
	for (unsigned int i=0; i< laser_.size(); i++){
        delete laser_[i];
    }
	
}//END Destructer



// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump of the data contained in ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::dump(PicParams* params)
{
    //!\todo Check for none-cartesian grid & for generic grid (neither all dual or all primal) (MG & JD)
    
    std::vector<unsigned int> dimPrim; dimPrim.resize(1); dimPrim[0] = params->n_space[0]+1;
    std::vector<unsigned int> dimDual; dimDual.resize(1); dimDual[0] = params->n_space[0]+2;
    
    // dump of the electromagnetic fields
	Ex_->dump(dimDual);
	Ey_->dump(dimPrim);
	Ez_->dump(dimPrim);
	Bx_->dump(dimPrim);
	By_->dump(dimDual);
	Bz_->dump(dimDual);
    // dump of the total charge density & currents
	rho_->dump(dimPrim);
	Jx_->dump(dimDual);
	Jy_->dump(dimPrim);
	Jz_->dump(dimPrim); 
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used initialize the total charge density
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::initRho(vector<Species*> vecSpecies, Projector* Proj)
{
    //! \todo Check that one uses only none-test particles
    // number of (none-test) used in the simulation
	unsigned int n_species = vecSpecies.size();

    //loop on all (none-test) Species
	for (unsigned int iSpec=0 ; iSpec<n_species; iSpec++ )
    {
		std::vector<Particle*> cuParticles = vecSpecies[iSpec]->getParticlesList();
		unsigned int n_particles = vecSpecies[iSpec]->getNbrOfParticles();
		for (unsigned int iPart=0 ; iPart<n_particles; iPart++ )
        {
			(*Proj)( rho_ , cuParticles[iPart] );
		}
	}//iSpec

}
