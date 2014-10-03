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
ElectroMagn::ElectroMagn(PicParams &params, LaserParams &laser_params, SmileiMPI* smpi) :
timestep(params.timestep),
cell_length(params.cell_length),
nDim_field(params.nDim_field),
cell_volume(params.cell_volume),
n_species(params.n_species),
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
    rho_o=NULL;
    
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

    fieldsBoundCond = FieldsBC_Factory::create(params, laser_params);

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
      if (fieldsBoundCond[i]!=NULL) delete fieldsBoundCond[i];

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
    saveMagneticFields();

    // Compute Ex_, Ey_, Ez_
    solveMaxwellAmpere();
    // Exchange Ex_, Ey_, Ez_
    smpi->exchangeE( this );

    // Compute Bx_, By_, Bz_
    solveMaxwellFaraday();

    // Update Bx_, By_, Bz_
    if ((!simWindow) || (!simWindow->isMoving(time_dual)) )
        if (fieldsBoundCond[0]!=NULL) // <=> if !periodic
	    fieldsBoundCond[0]->apply(this, time_dual, smpi);
    if ( (fieldsBoundCond.size()>1) )
        if (fieldsBoundCond[1]!=NULL) // <=> if !periodic
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
void ElectroMagn::initRhoJ(vector<Species*> vecSpecies, Projector* Proj)
{
    //! \todo Check that one uses only none-test particles
    // number of (none-test) used in the simulation
    //! \todo fix this: n_species is already a member of electromagn, is it this confusing? what happens if n_species grows (i.e. with ionization)?
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




void ElectroMagn::movingWindow_x(unsigned int shift, SmileiMPI *smpi)
{
    if (fieldsBoundCond[0]!=NULL)
        fieldsBoundCond[0]->laserDisabled();

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

