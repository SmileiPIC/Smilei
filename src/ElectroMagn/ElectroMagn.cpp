#include "ElectroMagn.h"

#include <limits>
#include <iostream>

#include "PicParams.h"
#include "Species.h"
#include "Projector.h"
#include "Field.h"
#include "ElectroMagnBC.h"
#include "ElectroMagnBC_Factory.h"
#include "SimWindow.h"
#include "ExtFieldProfile1D.h"
#include "ExtFieldProfile2D.h"
#include "SolverFactory.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn(PicParams &params, InputData &input_data, SmileiMPI* smpi) :
laser_params(params, input_data),
extfield_params(params, input_data),
timestep(params.timestep),
cell_length(params.cell_length),
n_species(params.n_species),
nDim_field(params.nDim_field),
cell_volume(params.cell_volume),
n_space(params.n_space),
oversize(params.oversize)
{
    // initialize poynting vector
    poynting[0].resize(nDim_field,0.0);
    poynting[1].resize(nDim_field,0.0);
    poynting_inst[0].resize(nDim_field,0.0);
    poynting_inst[1].resize(nDim_field,0.0);
    
    // take useful things from params
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
    
    Ex_avg=NULL;
    Ey_avg=NULL;
    Ez_avg=NULL;
    Bx_avg=NULL;
    By_avg=NULL;
    Bz_avg=NULL;
    
    // Species charge currents and density
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

    
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<2; j++) {
            istart[i][j]=0;
            bufsize[i][j]=0;
        }
    }    

    emBoundCond = ElectroMagnBC_Factory::create(params, laser_params);
    
    MaxwellFaradaySolver_ = SolverFactory::create(params);

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

    if (Ex_avg!=NULL) {
        delete Ex_avg;
        delete Ey_avg;
        delete Ez_avg;
        delete Bx_avg;
        delete By_avg;
        delete Bz_avg;
    }

    for (unsigned int ispec=0; ispec<n_species; ispec++) {
      delete Jx_s[ispec];
      delete Jy_s[ispec];
      delete Jz_s[ispec];
      delete rho_s[ispec];
    }
  
    int nBC = emBoundCond.size();
    for ( int i=0 ; i<nBC ;i++ )
      if (emBoundCond[i]!=NULL) delete emBoundCond[i];

    delete MaxwellFaradaySolver_;

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
void ElectroMagn::solveMaxwell(int itime, double time_dual, SmileiMPI* smpi, PicParams &params, SimWindow* simWindow)
{
#pragma omp parallel
{
    // saving magnetic fields (to compute centered fields used in the particle pusher)
    saveMagneticFields();

    // Compute Ex_, Ey_, Ez_
    solveMaxwellAmpere();

#pragma omp single
{
    // Exchange Ex_, Ey_, Ez_
    smpi->exchangeE( this );
}// end single

    // Compute Bx_, By_, Bz_
    (*MaxwellFaradaySolver_)(this);

#pragma omp single
{
    // Update Bx_, By_, Bz_
    if ((!simWindow) || (!simWindow->isMoving(time_dual)) )
        if (emBoundCond[0]!=NULL) // <=> if !periodic
	    emBoundCond[0]->apply(this, time_dual, smpi);
    if ( (emBoundCond.size()>1) )
        if (emBoundCond[1]!=NULL) // <=> if !periodic
	    emBoundCond[1]->apply(this, time_dual, smpi);
 
    // Exchange Bx_, By_, Bz_
    smpi->exchangeB( this );
}// end single

    // Compute Bx_m, By_m, Bz_m
    centerMagneticFields();
} // end parallel
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump of the data contained in ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::dump()
{
    //!\todo Check for none-cartesian grid & for generic grid (neither all dual or all primal) (MG & JD)
    
    vector<unsigned int> dimPrim;
    dimPrim.resize(1);
    dimPrim[0] = n_space[0]+2*oversize[0]+1;
    vector<unsigned int> dimDual;
    dimDual.resize(1);
    dimDual[0] = n_space[0]+2*oversize[0]+2;
    
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
void ElectroMagn::initRhoJ(vector<Species*>& vecSpecies, Projector* Proj)
{
    //! \todo Check that one uses only none-test particles
    // number of (none-test) used in the simulation
    //! \todo fix this: n_species is already a member of electromagn, is it this confusing? what happens if n_species grows (i.e. with ionization)?
    unsigned int n_species = vecSpecies.size();
    
    //loop on all (none-test) Species
    for (unsigned int iSpec=0 ; iSpec<n_species; iSpec++ ) {
        Particles &cuParticles = vecSpecies[iSpec]->getParticlesList();
        unsigned int n_particles = vecSpecies[iSpec]->getNbrOfParticles();
        
        DEBUG(n_particles<<" species "<<iSpec);
        for (unsigned int iPart=0 ; iPart<n_particles; iPart++ ) {
            // project charge & current densities
            (*Proj)(Jx_s[iSpec], Jy_s[iSpec], Jz_s[iSpec], rho_s[iSpec], cuParticles, iPart,
                    cuParticles.lor_fac(iPart));
        }
        
    }//iSpec
    DEBUG("before computeTotalRhoJ");    
    computeTotalRhoJ();
    DEBUG("projection done for initRhoJ");
    
}




void ElectroMagn::movingWindow_x(unsigned int shift, SmileiMPI *smpi)
{
    //! \todo{ Why the if test ? Remove it ? (AB for JD)}
    if (emBoundCond[0]!=NULL)
        emBoundCond[0]->laserDisabled();

    // For nrj balance
    nrj_mw_lost += computeNRJ(shift, smpi);

    smpi->exchangeE( this, shift );

    smpi->exchangeB( this, shift );
    
    smpi->exchangeBm( this, shift );

    if (Ex_avg!=NULL) {
        Ex_avg->shift_x(shift);
        Ey_avg->shift_x(shift);
        Ez_avg->shift_x(shift);
        Bx_avg->shift_x(shift);
        By_avg->shift_x(shift);
        Bz_avg->shift_x(shift);
        smpi->exchangeAvg( this );
    }

    // For now, fields introduced with moving window set to 0 
    nrj_new_fields =+ 0.;

    
    //Here you might want to apply some new boundary conditions on the +x boundary. For the moment, all fields are set to 0.
}

double ElectroMagn::computeNRJ(unsigned int shift, SmileiMPI *smpi) {
    double nrj(0.);

    if ( smpi->isWestern() ) {
	nrj += Ex_->computeNRJ(shift, istart, bufsize);
	nrj += Ey_->computeNRJ(shift, istart, bufsize);
	nrj += Ez_->computeNRJ(shift, istart, bufsize);

	nrj += Bx_m->computeNRJ(shift, istart, bufsize);
	nrj += By_m->computeNRJ(shift, istart, bufsize);
	nrj += Bz_m->computeNRJ(shift, istart, bufsize);

    }

    return nrj;
}

bool ElectroMagn::isRhoNull(SmileiMPI* smpi)
{
    double norm2(0.);
    double locnorm2(0.);

    // rho_->isDual(i) = 0 for all i
    // istart[i][0] & bufsize[i][0]

    vector<unsigned int> iFieldStart(3,0), iFieldEnd(3,1), iFieldGlobalSize(3,1);
    for (unsigned int i=0 ; i<rho_->isDual_.size() ; i++ ) {
	iFieldStart[i] = istart[i][0];
	iFieldEnd [i] = iFieldStart[i] + bufsize[i][0];
	iFieldGlobalSize [i] = rho_->dims_[i];
    }

    for (unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++) {
	for (unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++) {
	    for (unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++) {
		unsigned int ii=k+ j*iFieldGlobalSize[2] +i*iFieldGlobalSize[1]*iFieldGlobalSize[2];
		locnorm2 += (*rho_)(ii)*(*rho_)(ii);
	    }
	}
    }

    MPI_Allreduce(&locnorm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return (norm2<=0.);

}

string LowerCase(string in){
    string out=in;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

void ElectroMagn::applyExternalFields(SmileiMPI* smpi) {    
    vector<Field*> my_fields;
    my_fields.push_back(Ex_);
    my_fields.push_back(Ey_);
    my_fields.push_back(Ez_);
    my_fields.push_back(Bx_);
    my_fields.push_back(By_);
    my_fields.push_back(Bz_);
    
    for (vector<Field*>::iterator field=my_fields.begin(); field!=my_fields.end(); field++) {
        if (*field) {
            for (vector<ExtFieldStructure>::iterator extfield=extfield_params.structs.begin(); extfield!=extfield_params.structs.end(); extfield++ ) {
                ExtFieldProfile *my_ExtFieldProfile=NULL;
                if (extfield_params.geometry == "1d3v") {
                    my_ExtFieldProfile = (ExtFieldProfile*) (new ExtFieldProfile1D(*extfield));
                } else if (extfield_params.geometry == "2d3v") {
                    my_ExtFieldProfile = (ExtFieldProfile*) (new ExtFieldProfile2D(*extfield));
                }
                if (my_ExtFieldProfile) {
                    for (vector<string>::iterator fieldName=(*extfield).fields.begin();fieldName!=(*extfield).fields.end();fieldName++) {
                        if (LowerCase((*field)->name)==LowerCase(*fieldName)) {
                            applyExternalField(*field,my_ExtFieldProfile, smpi);
                        }
                    }
                    delete my_ExtFieldProfile;
                    my_ExtFieldProfile=NULL;
                } else{
		    ERROR("Could not initialize external field Profile");
		}
            }
        }
    }
    Bx_m->copyFrom(Bx_);
    By_m->copyFrom(By_);
    Bz_m->copyFrom(Bz_);
}    
