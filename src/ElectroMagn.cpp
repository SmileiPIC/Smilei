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
ElectroMagn::ElectroMagn(PicParams* params, SmileiMPI* smpi)
{

	// take useful things from params
	cell_volume=params->cell_volume;
	n_space=params->n_space;
	
	oversize.resize(3,0);
	for (unsigned int i=0; i<params->oversize.size(); i++) {
		oversize[i]=params->oversize[i];
	}

	for (unsigned int i=0; i<3; i++) {
		DEBUG("____________________ OVERSIZE: " <<i << " " << oversize[i]);
	}
	
	if (n_space.size() != 3) ERROR("this should not happend");
		
	// check for laser conditions
	laser_.resize(params->n_laser);
    
	for (unsigned int i=0; i<laser_.size(); i++) {
		DEBUG(5,"Initializing Laser "<<i);
		laser_[i] = new Laser(params->sim_time, params->laser_param[i]);
	}

	Ex_=NULL;
	Ey_=NULL;
	Ez_=NULL;
	Bx_=NULL;
	By_=NULL;
	Bz_=NULL;
	Bx_m=NULL;
	By_m=NULL;
	Bz_m=NULL;
	Jx_=NULL;
	Jy_=NULL;
	Jz_=NULL;
	rho_=NULL;
	rho_o=NULL;		
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
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
/*void ElectroMagn::solveMaxwell(double time_dual, SmileiMPI* smpi)
{
	//solve Maxwell's equations
	solveMaxwellAmpere();
	//smpi->exchangeE( EMfields );
	solveMaxwellFaraday();
	smpi->exchangeB( this );
	boundaryConditions(time_dual, smpi);

}*/


// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump of the data contained in ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::dump(PicParams* params)
{
    //!\todo Check for none-cartesian grid & for generic grid (neither all dual or all primal) (MG & JD)
    
    std::vector<unsigned int> dimPrim; dimPrim.resize(1); dimPrim[0] = params->n_space[0]+2*params->oversize[0]+1;
    std::vector<unsigned int> dimDual; dimDual.resize(1); dimDual[0] = params->n_space[0]+2*params->oversize[0]+2;
    
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
// Method used to initialize the total charge density
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
        DEBUG(n_particles<<" species "<<iSpec);
		for (unsigned int iPart=0 ; iPart<n_particles; iPart++ )
        {
			(*Proj)( rho_ , cuParticles[iPart] );
		}
        DEBUG("projection done for initRho");
	}//iSpec

}

void ElectroMagn::computeScalars()
{
	
	std::vector<Field*> fields;

	fields.push_back(Ex_);
	fields.push_back(Ey_);
	fields.push_back(Ez_);
	fields.push_back(Bx_m);
	fields.push_back(By_m);
	fields.push_back(Bz_m);

	for (vector<Field*>::iterator field=fields.begin(); field!=fields.end(); field++) {

		map<string,vector<double> > scalars_map;
		
		vector<double> Etot(1);
		
		for (unsigned int k=oversize[2]; k<n_space[2]-oversize[2]; k++) {
			for (unsigned int j=oversize[1]; j<n_space[1]-oversize[1]; j++) {
				for (unsigned int i=oversize[0]; i<n_space[0]-oversize[0]; i++) {
					unsigned int ii=i+j*n_space[0]+k*n_space[0]*n_space[1];
					Etot[0]+=pow((**field)(ii),2);
				}
			}
		}
		Etot[0]*=0.5*cell_volume;
		scalars_map["Etot"]=Etot;
		
		scalars[(*field)->name+"_U"]=scalars_map;
	}
	
	fields.push_back(Jx_);
	fields.push_back(Jy_);
	fields.push_back(Jz_);
	fields.push_back(rho_);	
	
	

	for (vector<Field*>::iterator field=fields.begin(); field!=fields.end(); field++) {
		
		map<string,vector<double> > scalars_map;
	
		
// this does not work!		
		vector<double> minVec(4);
		vector<double> maxVec(4);
		
		minVec[0]=(**field)(0);
		maxVec[0]=(**field)(0);
		minVec[1]=maxVec[1]=0;
		minVec[2]=maxVec[2]=0;
		minVec[3]=maxVec[3]=0;
		
		for (unsigned int k=oversize[2]; k<n_space[2]-oversize[2]; k++) {
			for (unsigned int j=oversize[1]; j<n_space[1]-oversize[1]; j++) {
				for (unsigned int i=oversize[0]; i<n_space[0]-oversize[0]; i++) {
					unsigned int ii=i+j*n_space[0]+k*n_space[0]*n_space[1];
					if (minVec[0]>(**field)(ii)) {
						minVec[0]=(**field)(ii);
						minVec[1]=i;
						minVec[2]=j;
						minVec[3]=k;
					}
					if (maxVec[0]<(**field)(ii)) {
						maxVec[0]=(**field)(ii);
						maxVec[1]=i;
						maxVec[2]=j;
						maxVec[3]=k;
					}
				}
			}
		}
		minVec.resize(1+(*field)->dims_.size());
		maxVec.resize(1+(*field)->dims_.size());

		
//		unsigned int tot_size_field=1;
//		for (unsigned int i =0; i<(*field)->dims_.size(); i++) {
//			tot_size_field*=(*field)->dims_[i];
//		}
//		vector<double> minVec(1+(*field)->dims_.size());
//		vector<double> maxVec(1+(*field)->dims_.size());
//		
//		minVec[0]=(**field)(0);
//		maxVec[0]=(**field)(0);
//		for (unsigned int ii =0; ii<tot_size_field; ii++){
//			
//			if (minVec[0]>(**field)(ii)) {
//				minVec[0]=(**field)(ii);
//				for (unsigned int i =0; i<(*field)->dims_.size(); i++) {
//					minVec[1+i]=ii;
//				}
//			}
//			if (maxVec[0]<(**field)(ii)) {
//				maxVec[0]=(**field)(ii);
//				for (unsigned int i =0; i<(*field)->dims_.size(); i++) {
//					maxVec[1+i]=ii;
//				}
//			}
//			
//			
//		}
					
		// we just store the values that change
		scalars_map["min"]=minVec;
		scalars_map["max"]=maxVec;
		scalars[(*field)->name]=scalars_map;
	}
}


