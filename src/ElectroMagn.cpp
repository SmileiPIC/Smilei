#include "ElectroMagn.h"

#include <limits>
#include <iostream>

#include "ElectroMagn1D.h"
#include "PicParams.h"
#include "Species.h"
#include "Projector.h"
#include "Laser.h"
#include "Field.h"
#include "FieldsBC.h"
#include "FieldsBC_Factory.h"
#include "SimWindow.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn(PicParams* params, SmileiMPI* smpi)
{
    params_ = params;

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

    fieldsBoundCond = FieldsBC_Factory::create(*params);

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

    int nBC = fieldsBoundCond.size();
    for ( int i=0 ; i<nBC ;i++ )
      delete fieldsBoundCond[i];

}//END Destructer

// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::solveMaxwell(int itime, double time_dual, SmileiMPI* smpi, PicParams &params, SimWindow* simWindow)
{
    saveMagneticFields();

    // Compute Ex_, Ey_, Ez_
    solveMaxwellAmpere();
    // Exchange Ex_, Ey_, Ez_
    smpi->exchangeE( this );

    // Compute Bx_, By_, Bz_
    solveMaxwellFaraday();

    // Update Bx_, By_, Bz_
    if ((!simWindow) || (!simWindow->isMoving(itime)) )
	fieldsBoundCond[0]->apply(this, time_dual, smpi);
    if ( (!params.use_transverse_periodic) && (fieldsBoundCond.size()>1) )
	fieldsBoundCond[1]->apply(this, time_dual, smpi);
 
    // Exchange Bx_, By_, Bz_
    smpi->exchangeB( this );

    // Compute Bx_m, By_m, Bz_m
    centerMagneticFields();

}


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


void ElectroMagn::movingWindow_x(unsigned int shift, SmileiMPI *smpi)
{
    Ex_->shift_x(shift);
    Ey_->shift_x(shift);
    Ez_->shift_x(shift);
    //! \ Comms to optimize, only in x, east to west 
    //! \ Implement SmileiMPI::exchangeE( EMFields*, int nDim, int nbNeighbors );
    smpi->exchangeE( this );

    Bx_->shift_x(shift);
    By_->shift_x(shift);
    Bz_->shift_x(shift);
    smpi->exchangeB( this );

    //Here you might want to apply some new boundary conditions on the +x boundary. For the moment, all fields are set to 0. 
}

