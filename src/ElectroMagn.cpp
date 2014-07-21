#include "ElectroMagn.h"

#include <limits>
#include <iostream>

#include "ElectroMagn1D.h"
#include "PicParams.h"
#include "Species.h"
#include "Projector.h"
#include "Laser.h"
#include "Field.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn(PicParams* params, SmileiMPI* smpi)
{

    poynting[0].resize(params->nDim_field,0.0);
    poynting[1].resize(params->nDim_field,0.0);

    // take useful things from params
    cell_volume=params->cell_volume;
    n_space=params->n_space;

    oversize=params->oversize;

    for (unsigned int i=0; i<3; i++) {
        DEBUG("____________________ OVERSIZE: " <<i << " " << oversize[i]);
    }

    if (n_space.size() != 3) ERROR("this should not happend");

    // check for laser conditions
    laser_.resize(params->n_laser);

    for (unsigned int i=0; i<laser_.size(); i++) {
        DEBUG(5,"Initializing Laser "<<i);
        laser_[i] = new Laser(params->sim_time, params->sim_length[1], params->laser_param[i]);
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
    
    Ex_avg=NULL;
    Ey_avg=NULL;
    Ez_avg=NULL;
    Bx_avg=NULL;
    By_avg=NULL;
    Bz_avg=NULL;

    // Species charge currents and density
    n_species = params->n_species;
    Jx_s.resize(n_species);
    Jy_s.resize(n_species);
    Jz_s.resize(n_species);
    rho_s.resize(n_species);
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = NULL;
        Jy_s[ispec]  = NULL;
        Jz_s[ispec]  = NULL;
        rho_s[ispec] = NULL;
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
    for (unsigned int i=0; i< laser_.size(); i++) {
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

    vector<unsigned int> dimPrim;
    dimPrim.resize(1);
    dimPrim[0] = params->n_space[0]+2*params->oversize[0]+1;
    vector<unsigned int> dimDual;
    dimDual.resize(1);
    dimDual[0] = params->n_space[0]+2*params->oversize[0]+2;

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
void ElectroMagn::initRhoJ(vector<Species*> vecSpecies, Projector* Proj)
{
    //! \todo Check that one uses only none-test particles
    // number of (none-test) used in the simulation
    unsigned int n_species = vecSpecies.size();

    //loop on all (none-test) Species
    for (unsigned int iSpec=0 ; iSpec<n_species; iSpec++ ) {
        Particles cuParticles = vecSpecies[iSpec]->getParticlesList();
        unsigned int n_particles = vecSpecies[iSpec]->getNbrOfParticles();

        DEBUG(n_particles<<" species "<<iSpec);
        for (unsigned int iPart=0 ; iPart<n_particles; iPart++ ) {
            // project charge & current densities
            (*Proj)(Jx_s[iSpec], Jy_s[iSpec], Jz_s[iSpec], rho_s[iSpec], cuParticles, iPart,
                    cuParticles.lor_fac(iPart));
        }

    }//iSpec

    computeTotalRhoJ();
    DEBUG("projection done for initRhoJ");

}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to compute EM fields related scalar quantities used in diagnostics
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::computeScalars()
{

    vector<Field*> fields;

    fields.push_back(Ex_);
    fields.push_back(Ey_);
    fields.push_back(Ez_);
    fields.push_back(Bx_m);
    fields.push_back(By_m);
    fields.push_back(Bz_m);

    for (vector<Field*>::iterator field=fields.begin(); field!=fields.end(); field++) {

        map<string,vector<double> > scalars_map;

        vector<double> Etot(1);

	unsigned int iFieldStart[3], iFieldSize[3];
	for ( int i=0 ; i<(*field)->isDual_.size() ; i++ ) {
	    iFieldStart[i] = istart [i][(*field)->isDual(i)];
	    iFieldSize [i] = bufsize[i][(*field)->isDual(i)];
	}

	for (unsigned int k=iFieldStart[2]; k<iFieldSize[2]; k++) {
	    for (unsigned int j=iFieldStart[1]; j<iFieldSize[1]; j++) {
		for (unsigned int i=iFieldStart[0]; i<iFieldSize[0]; i++) {
		    unsigned int ii=i+j*n_space[0]+k*n_space[0]*n_space[1];
		    Etot[0]+=pow((**field)(ii),2);
		}
	    }
        }
        Etot[0]*=0.5*cell_volume;
        scalars_map["sum"]=Etot;

        scalars[(*field)->name+"_U"]=scalars_map;
    }

    fields.push_back(Jx_);
    fields.push_back(Jy_);
    fields.push_back(Jz_);
    fields.push_back(rho_);



    for (vector<Field*>::iterator field=fields.begin(); field!=fields.end(); field++) {

        map<string,vector<double> > scalars_map;

        vector<double> minVec(4);
        vector<double> maxVec(4);

        minVec[0]=(**field)(0);
        maxVec[0]=(**field)(0);
        minVec[1]=maxVec[1]=0;
        minVec[2]=maxVec[2]=0;
        minVec[3]=maxVec[3]=0;

	unsigned int iFieldStart[3], iFieldSize[3];
	for ( int i=0 ; i<(*field)->isDual_.size() ; i++ ) {
	    iFieldStart[i] = istart [i][(*field)->isDual(i)];
	    iFieldSize [i] = bufsize[i][(*field)->isDual(i)];
	}

	for (unsigned int k=iFieldStart[2]; k<iFieldSize[2]; k++) {
	    for (unsigned int j=iFieldStart[1]; j<iFieldSize[1]; j++) {
		for (unsigned int i=iFieldStart[0]; i<iFieldSize[0]; i++) {
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


        // we just store the values that change
        scalars_map["min"]=minVec;
        scalars_map["max"]=maxVec;
        scalars[(*field)->name]=scalars_map;
    }
    
    // poynting stuff
    map<string,vector<double> > poynting_map_inf;
    poynting_map_inf["sum"]=poynting[0];
    scalars["Poy_inf"]=poynting_map_inf;

    map<string,vector<double> > poynting_map_sup;
    poynting_map_sup["sum"]=poynting[1];
    scalars["Poy_sup"]=poynting_map_sup;
    
    
}
