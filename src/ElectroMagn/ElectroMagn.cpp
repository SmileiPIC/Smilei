#include "ElectroMagn.h"

#include <limits>
#include <iostream>

#include "Params.h"
#include "Species.h"
#include "Projector.h"
#include "Field.h"
#include "ElectroMagnBC.h"
#include "ElectroMagnBC_Factory.h"
#include "SimWindow.h"
#include "Patch.h"
#include "Profile.h"
#include "SolverFactory.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn(Params &params, vector<Species*>& vecSpecies, Patch* patch) :
timestep(params.timestep),
cell_length(params.cell_length),
n_species(vecSpecies.size()),
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
    
    if (n_space.size() != 3) ERROR("this should not happen");
    
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
    
    
    emBoundCond = ElectroMagnBC_Factory::create(params, patch);
    
    MaxwellFaradaySolver_ = SolverFactory::create(params);
    
}


void ElectroMagn::finishInitialization(int nspecies, Patch* patch)
{

    initAntennas(patch);
    
    // Fill allfields
    allFields.push_back(Ex_ );
    allFields.push_back(Ey_ );
    allFields.push_back(Ez_ );
    allFields.push_back(Bx_ );
    allFields.push_back(By_ );
    allFields.push_back(Bz_ );
    allFields.push_back(Bx_m);
    allFields.push_back(By_m);
    allFields.push_back(Bz_m);
    allFields.push_back(Jx_ );
    allFields.push_back(Jy_ );
    allFields.push_back(Jz_ );
    allFields.push_back(rho_);

    for (int ispec=0; ispec<nspecies; ispec++) {
        allFields.push_back(Jx_s[ispec] );
        allFields.push_back(Jy_s[ispec] );
        allFields.push_back(Jz_s[ispec] );
        allFields.push_back(rho_s[ispec]);
    }
                
    allFields_avg.push_back(Ex_avg);
    allFields_avg.push_back(Ey_avg);
    allFields_avg.push_back(Ez_avg);
    allFields_avg.push_back(Bx_avg);
    allFields_avg.push_back(By_avg);
    allFields_avg.push_back(Bz_avg);
    
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
    
    //antenna cleanup
    for (vector<Antenna>::iterator antenna=antennas.begin(); antenna!=antennas.end(); antenna++ ) {
        delete antenna->field;
        antenna->field=NULL;
    }
    
}//END Destructer


void ElectroMagn::clean()
{
    for ( int i=0 ; i<emBoundCond.size() ;i++ )
        if (emBoundCond[i]!=NULL) emBoundCond[i]->clean();
}

// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
// In the main program 
//     - saveMagneticFields
//     - solveMaxwellAmpere
//     - solveMaxwellFaraday
//     - boundaryConditions
//     - vecPatches::exchangeB (patch & MPI sync)
//     - centerMagneticFields



void ElectroMagn::boundaryConditions(int itime, double time_dual, Patch* patch, Params &params, SimWindow* simWindow)
{
    // Compute EM Bcs
    if ( (!simWindow) || (!simWindow->isMoving(time_dual)) ) {
        if (emBoundCond[0]!=NULL) { // <=> if !periodic
            emBoundCond[0]->apply_xmin(this, time_dual, patch);
            emBoundCond[1]->apply_xmax(this, time_dual, patch);
        }
    }
    if (emBoundCond.size()>2) {
        if (emBoundCond[2]!=NULL) {// <=> if !periodic
            emBoundCond[2]->apply_ymin(this, time_dual, patch);
            emBoundCond[3]->apply_ymax(this, time_dual, patch);
        }
    }
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
    // number of (none-test) used in the simulation
    //! \todo fix this: n_species is already a member of electromagn, is it this confusing? what happens if n_species grows (i.e. with ionization)?
    unsigned int n_species = vecSpecies.size();
    
    //loop on all (none-test) Species
    for (unsigned int iSpec=0 ; iSpec<n_species; iSpec++ ) {
        Particles &cuParticles = vecSpecies[iSpec]->getParticlesList();
        unsigned int n_particles = vecSpecies[iSpec]->getNbrOfParticles();
        
        DEBUG(n_particles<<" species "<<iSpec);
        if (!cuParticles.isTest) {
            for (unsigned int iPart=0 ; iPart<n_particles; iPart++ ) {
                // project charge & current densities
                (*Proj)(Jx_s[iSpec], Jy_s[iSpec], Jz_s[iSpec], rho_s[iSpec], cuParticles, iPart, cuParticles.lor_fac(iPart));
            }
        }
        
    }//iSpec
    DEBUG("before computeTotalRhoJ");    
    computeTotalRhoJ();
    DEBUG("projection done for initRhoJ");
    
}


void ElectroMagn::movingWindow_x(unsigned int shift)
{
    //! \todo{ Why the if test ? Remove it ? (AB for JD)}
    if (emBoundCond[0]!=NULL)
        emBoundCond[0]->laserDisabled();

    // For nrj balance
    //nrj_mw_lost += computeNRJ(); // Integreated in SimWindow::operate

    // For now, fields introduced with moving window set to 0 
    nrj_new_fields =+ 0.;

    
    //Here you might want to apply some new boundary conditions on the +x boundary. For the moment, all fields are set to 0.
}

void ElectroMagn::laserDisabled()
{
    if ( emBoundCond.size() )
	emBoundCond[0]->laserDisabled();
}

double ElectroMagn::computeNRJ() {
    double nrj(0.);

    nrj += Ex_->norm2(istart, bufsize);
    nrj += Ey_->norm2(istart, bufsize);
    nrj += Ez_->norm2(istart, bufsize);

    nrj += Bx_m->norm2(istart, bufsize);
    nrj += By_m->norm2(istart, bufsize);
    nrj += Bz_m->norm2(istart, bufsize);

    return nrj;
}

string LowerCase(string in){
    string out=in;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

void ElectroMagn::applyExternalFields(Patch* patch) {    
    Field * field;
    bool found=false;
    for (vector<ExtField>::iterator extfield=extFields.begin(); extfield!=extFields.end(); extfield++ ) {
        for (vector<string>::iterator fieldName=extfield->fields.begin();fieldName!=extfield->fields.end();fieldName++) {
            string name = LowerCase(*fieldName);
            if      ( Ex_ && name==LowerCase(Ex_->name) ) field = Ex_;
            else if ( Ey_ && name==LowerCase(Ey_->name) ) field = Ey_;
            else if ( Ez_ && name==LowerCase(Ez_->name) ) field = Ez_;
            else if ( Bx_ && name==LowerCase(Bx_->name) ) field = Bx_;
            else if ( By_ && name==LowerCase(By_->name) ) field = By_;
            else if ( Bz_ && name==LowerCase(Bz_->name) ) field = Bz_;
            else field = NULL;
            
            if( field ) {
                applyExternalField( field, extfield->profile, patch );
                found=true;
            }
        }
    }
    if (patch->isMaster()) {
        if (found) {
            MESSAGE(1,"Finish");
        } else {
            MESSAGE(1,"Nothing to do");
        }
    }
    Bx_m->copyFrom(Bx_);
    By_m->copyFrom(By_);
    Bz_m->copyFrom(Bz_);
}

void ElectroMagn::applyAntennas(SmileiMPI* smpi, double time) {
    Field * field;
    
    for (vector<Antenna>::iterator antenna=antennas.begin(); antenna!=antennas.end(); antenna++ ) {
        if (antenna->field) {
            double intensity = antenna->time_profile->valueAt(time);
            
            if     ( antenna->field->name == "Jx" ) field = Jx_;
            else if( antenna->field->name == "Jy" ) field = Jy_;
            else if( antenna->field->name == "Jz" ) field = Jz_;
            
            if (antenna->field->globalDims_ == field->globalDims_) // to do (TV): is this check really necessary ?
                for (unsigned int i=0; i< field->globalDims_ ; i++)
                    (*field)(i) += intensity * (*antenna->field)(i);
        }
    }
}

