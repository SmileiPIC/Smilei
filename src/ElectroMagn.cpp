
#include "ElectroMagn.h"
#include "ElectroMagn1D.h"

#include "PicParams.h"
#include "Species.h"
#include "Projector.h"

#include "Laser.h"
#include "Field.h"

#include <iostream>
using namespace std;

#include <limits>

//! Creator for the virtual class ElectroMagn
ElectroMagn::ElectroMagn(PicParams* params, SmileiMPI* smpi)
{

	laser_.resize(params->n_laser);
	for (unsigned int i=0; i<laser_.size(); i++) {
		DEBUG(5,"Initializing Laser "<<i);
		laser_[i] = new Laser(params->laser_param[i]);
	}

}

//! Destructor for the virtual class ElectroMagn
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


void ElectroMagn::dump()
{
	Ex_->dump(dimDual);
	Ey_->dump(dimPrim);
	Ez_->dump(dimPrim);
	Bx_->dump(dimPrim);
	By_->dump(dimDual);
	Bz_->dump(dimDual);

	rho_->dump(dimPrim);
	Jx_->dump(dimDual);
	Jy_->dump(dimPrim);
	Jz_->dump(dimPrim);
}


void ElectroMagn::initRho(vector<Species*> vecSpecies, Projector* Proj)
{
	unsigned int n_species = vecSpecies.size();
	unsigned int n_particles;
	
	for (unsigned int iSpec=0 ; iSpec<n_species; iSpec++ ) {
		std::vector<Particle*> cuParticles = vecSpecies[iSpec]->getParticlesList();
		n_particles = vecSpecies[iSpec]->getNbrOfParticles();
		for (unsigned int iPart=0 ; iPart<n_particles; iPart++ ) {
			(*Proj)( rho_ , cuParticles[iPart] );
		}
	}

}


