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
#include "DomainDecompositionFactory.h"
#include "LaserEnvelope.h"

#include "PatchAM.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn( Params &params, vector<Species *> &vecSpecies, Patch *patch ) :
    timestep( params.timestep ),
    cell_length( params.cell_length ),
    n_species( vecSpecies.size() ),
    nDim_field( params.nDim_field ),
    cell_volume( params.cell_volume ),
    isXmin( patch->isXmin() ),
    isXmax( patch->isXmax() ),
    is_pxr( params.is_pxr ),
    nrj_mw_out( 0. ),
    nrj_mw_inj( 0. ),
    filter_( NULL ),
    use_BTIS3(params.use_BTIS3)
{
    size_ = patch->size_;
    oversize = patch->oversize;
    // Dimension of the primal and dual grids
    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    for( unsigned int i=0 ; i < nDim_field ; i++ ) {
        dimPrim[i] = size_[i] + 2*oversize[i] + 1;
        dimDual[i] = size_[i] + 2*oversize[i] + 2 - (params.is_pxr);
    }
    
    // take useful things from params
    initElectroMagnQuantities();
    emBoundCond = ElectroMagnBC_Factory::create( params, patch );
    MaxwellAmpereSolver_  = SolverFactory::createMA( params );
    MaxwellFaradaySolver_ = SolverFactory::createMF( params );
    
    envelope = NULL;
    
}

// ---------------------------------------------------------------------------------------------------------------------
// ElectroMagn constructor for patches
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::ElectroMagn( ElectroMagn *emFields, Params &params, Patch *patch ) :
    timestep( emFields->timestep ),
    cell_length( emFields->cell_length ),
    n_species( emFields->n_species ),
    nDim_field( emFields->nDim_field ),
    cell_volume( emFields->cell_volume ),
    size_( emFields->size_ ),
    oversize( emFields->oversize ),
    isXmin( patch->isXmin() ),
    isXmax( patch->isXmax() ),
    is_pxr( emFields->is_pxr ),
    nrj_mw_out( 0. ),
    nrj_mw_inj( 0. ),
    filter_( NULL ),
    use_BTIS3(emFields->use_BTIS3)
{
    dimPrim = emFields->dimPrim;
    dimDual = emFields->dimDual;
    
    initElectroMagnQuantities();
    
    emBoundCond = ElectroMagnBC_Factory::create( params, patch );
    
    MaxwellAmpereSolver_  = SolverFactory::createMA( params );
    MaxwellFaradaySolver_ = SolverFactory::createMF( params );
    
    envelope = NULL;
}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize quantities used in ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::initElectroMagnQuantities()
{
    // initialize poynting vector
    poynting[0].resize( nDim_field, 0.0 );
    poynting[1].resize( nDim_field, 0.0 );
    poynting_inst[0].resize( nDim_field, 0.0 );
    poynting_inst[1].resize( nDim_field, 0.0 );
    
    Ex_=NULL;
    Ey_=NULL;
    Ez_=NULL;
    Bx_=NULL;
    By_=NULL;
    Bz_=NULL;
    Bx_m=NULL;
    By_m=NULL;
    Bz_m=NULL;
    By_mBTIS3=NULL;
    Bz_mBTIS3=NULL;
    Jx_=NULL;
    Jy_=NULL;
    Jz_=NULL;
    rho_=NULL;
    Env_A_abs_=NULL;
    Env_Chi_  =NULL;
    Env_E_abs_=NULL;
    Env_Ex_abs_=NULL;
    
    
    // Species charge currents and density
    Jx_s.resize( n_species );
    Jy_s.resize( n_species );
    Jz_s.resize( n_species );
    rho_s.resize( n_species );
    
    Env_Chi_s.resize( n_species );
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        Jx_s[ispec]  = NULL;
        Jy_s[ispec]  = NULL;
        Jz_s[ispec]  = NULL;
        rho_s[ispec] = NULL;
        Env_Chi_s[ispec] = NULL;
    }
    
    for( unsigned int i=0; i<3; i++ ) {
        for( unsigned int j=0; j<2; j++ ) {
            istart[i][j]=0;
            bufsize[i][j]=0;
        }
    }
}


void ElectroMagn::finishInitialization( int nspecies, Patch * )
{

    // Fill allfields
    allFields.push_back( Ex_ );
    allFields.push_back( Ey_ );
    allFields.push_back( Ez_ );
    allFields.push_back( Bx_ );
    allFields.push_back( By_ );
    allFields.push_back( Bz_ );
    allFields.push_back( Bx_m );
    allFields.push_back( By_m );
    allFields.push_back( Bz_m );
    allFields.push_back( Jx_ );
    allFields.push_back( Jy_ );
    allFields.push_back( Jz_ );
    allFields.push_back( rho_ );
    if( Env_A_abs_ != NULL ) {
        allFields.push_back( Env_A_abs_ );
        allFields.push_back( Env_Chi_ );
        allFields.push_back( Env_E_abs_ );
        allFields.push_back( Env_Ex_abs_ );
    }
    if ( By_mBTIS3 != NULL ){
        allFields.push_back( By_mBTIS3 );
        allFields.push_back( Bz_mBTIS3 );
    }
    
    // For species-related fields
    // The order is necessary in DiagnosticProbes - DO NOT CHANGE -
    species_starts.resize( 0 );
    for( int ispec=0; ispec<nspecies; ispec++ ) {
        species_starts.push_back( allFields.size() );
        allFields.push_back( Jx_s[ispec] );
        allFields.push_back( Jy_s[ispec] );
        allFields.push_back( Jz_s[ispec] );
        allFields.push_back( rho_s[ispec] );
        if( Env_A_abs_ != NULL ) {
            allFields.push_back( Env_Chi_s[ispec] );
        }
    }
    species_starts.push_back( allFields.size() );
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for the virtual class ElectroMagn
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn::~ElectroMagn()
{

    if( Ex_ != NULL ) {
        delete Ex_;
    }
    if( Ey_ != NULL ) {
        delete Ey_;
    }
    if( Ez_ != NULL ) {
        delete Ez_;
    }
    if( Bx_ != NULL ) {
        delete Bx_;
    }
    if( By_ != NULL ) {
        delete By_;
    }
    if( Bz_ != NULL ) {
        delete Bz_;
    }
    if( !is_pxr ) {
        if( Bx_m != NULL ) {
            delete Bx_m;
        }
        if( By_m != NULL ) {
            delete By_m;
        }
        if( Bz_m != NULL ) {
            delete Bz_m;
        }
    }
    if( Jx_ != NULL ) {
        delete Jx_;
    }
    if( Jy_ != NULL ) {
        delete Jy_;
    }
    if( Jz_ != NULL ) {
        delete Jz_;
    }
    if( rho_ != NULL ) {
        delete rho_;
    }
    
    if( Env_A_abs_ != NULL ) {
        delete Env_A_abs_;
    }
    if( Env_Chi_   != NULL ) {
        delete Env_Chi_;
    }
    if( Env_E_abs_ != NULL ) {
        delete Env_E_abs_;
    }
    if( Env_Ex_abs_ != NULL ) {
        delete Env_Ex_abs_;
    }
    
    for( unsigned int idiag=0; idiag<allFields_avg.size(); idiag++ )
        for( unsigned int ifield=0; ifield<allFields_avg[idiag].size(); ifield++ ) {
            delete allFields_avg[idiag][ifield];
        }
        
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        if( Jx_s [ispec] ) {
            delete Jx_s [ispec];
        }
        if( Jy_s [ispec] ) {
            delete Jy_s [ispec];
        }
        if( Jz_s [ispec] ) {
            delete Jz_s [ispec];
        }
        if( rho_s[ispec] ) {
            delete rho_s[ispec];
        }
        if( Env_Chi_s[ispec] ) {
            delete Env_Chi_s[ispec];
        }
    }
    
    if( filter_ ) {
        delete filter_;
    }
    
    int nBC = emBoundCond.size();
    for( int i=0 ; i<nBC ; i++ )
        if( emBoundCond[i]!=NULL ) {
            delete emBoundCond[i];
        }
        
    delete MaxwellAmpereSolver_;
    delete MaxwellFaradaySolver_;
    
    if( envelope != NULL ) {
        delete envelope;
    }
    
    //antenna cleanup
    for( vector<Antenna>::iterator antenna=antennas.begin(); antenna!=antennas.end(); antenna++ ) {
        delete antenna->field;
        antenna->field=NULL;
    }
    
    if (use_BTIS3){
        if(By_mBTIS3 != NULL){
            delete By_mBTIS3;
        }
    
        if(Bz_mBTIS3 != NULL){
            delete Bz_mBTIS3;
        }
    }

//     for ( unsigned int iExt = 0 ; iExt < prescribedFields.size() ; iExt++ ) {
//     	delete prescribedFields[iExt].savedField;
//     	#pragma omp single
//         if (prescribedFields[iExt].profile!=NULL) {
// 	    	delete prescribedFields[iExt].profile;
// 	    }
//     }
//     prescribedFields.clear();
    
    /*for ( unsigned int iExt = 0 ; iExt < extFields.size() ; iExt++ ) {
        if (extFields[iExt].profile!=NULL) {
            delete extFields[iExt].profile;
            extFields[iExt].profile = NULL;
        } // Pb wih clones
    }*/
    
}//END Destructer


void ElectroMagn::updateGridSize( Params &params, Patch *patch )
{
    isXmin = patch->isXmin();
    isXmax = patch->isXmax();
    
    unsigned int i=0;
    {
        for( int isDual=0 ; isDual<2 ; isDual++ ) {
            bufsize[i][isDual] = size_[i] + 1;
        }
        
        for( int isDual=0 ; isDual<2 ; isDual++ ) {
            bufsize[i][isDual] += isDual;
            if( params.number_of_patches[i]!=1 ) {
            
                if( ( !isDual ) ) {
                    bufsize[i][isDual]--;
                } else if( isDual ) {
                    bufsize[i][isDual]--;
                    bufsize[i][isDual]--;
                }
                
            } // if ( params.number_of_patches[i]!=1 )
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0)
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



void ElectroMagn::boundaryConditions( double time_dual, Patch *patch, SimWindow *simWindow )
{
    // Compute EM Bcs
    if( !( simWindow && simWindow->isMoving( time_dual ) ) ) { //Boundary conditions are applied after moving the window.
        if( emBoundCond[0]!=NULL ) { // <=> if !periodic
            emBoundCond[0]->apply( this, time_dual, patch );
            emBoundCond[1]->apply( this, time_dual, patch );
        }
    }
    if( emBoundCond.size()>2 ) {
        if( emBoundCond[2]!=NULL ) { // <=> if !periodic
            emBoundCond[2]->apply( this, time_dual, patch );
        }
        if( emBoundCond[3]!=NULL ) { // <=> if !periodic
            emBoundCond[3]->apply( this, time_dual, patch );
        }
    }
    if( emBoundCond.size()>4 ) {
        if( emBoundCond[4]!=NULL ) { // <=> if !periodic
            emBoundCond[4]->apply( this, time_dual, patch );
            emBoundCond[5]->apply( this, time_dual, patch );
        }
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge densities and currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::restartRhoJ()
{
    Jx_ ->put_to( 0. );
    Jy_ ->put_to( 0. );
    Jz_ ->put_to( 0. );
    rho_->put_to( 0. );
}

void ElectroMagn::restartEnvChi()
{
    Env_Chi_->put_to( 0. );
}

void ElectroMagn::restartRhoJs()
{
    for( unsigned int ispec=0 ; ispec < n_species ; ispec++ ) {
        if( Jx_s [ispec] ) {
            Jx_s [ispec]->put_to( 0. );
        }
        if( Jy_s [ispec] ) {
            Jy_s [ispec]->put_to( 0. );
        }
        if( Jz_s [ispec] ) {
            Jz_s [ispec]->put_to( 0. );
        }
        if( rho_s[ispec] ) {
            rho_s[ispec]->put_to( 0. );
        }
    }
    restartRhoJ();
}

void ElectroMagn::restartEnvChis()
{
    for( unsigned int ispec=0 ; ispec < n_species ; ispec++ ) {
        if( Env_Chi_s [ispec] ) {
            Env_Chi_s [ispec]->put_to( 0. );
        }
    }
    Env_Chi_->put_to( 0. );
}

// ---------------------------------------------------------------------------------------------------------------------
// Increment an averaged field
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn::incrementAvgField( Field *field, Field *field_avg )
{
    for( unsigned int i=0; i<field->number_of_points_; i++ ) {
        ( *field_avg )( i ) += ( *field )( i );
    }
}//END incrementAvgField



void ElectroMagn::laserDisabled()
{
    for( unsigned int i=0; i<emBoundCond.size(); i++ ) {
        if( emBoundCond[i] ) {
            emBoundCond[i]->laserDisabled();
        }
    }
}

double ElectroMagn::computeEnergy()
{
    double nrj( 0. );
    
    nrj += Ex_->norm2( istart, bufsize );
    nrj += Ey_->norm2( istart, bufsize );
    nrj += Ez_->norm2( istart, bufsize );
    
    nrj += Bx_m->norm2( istart, bufsize );
    nrj += By_m->norm2( istart, bufsize );
    nrj += Bz_m->norm2( istart, bufsize );
    
    return nrj;
}

void ElectroMagn::applyExternalFields( Patch *patch )
{
    std::vector<bool> input, copy;
    input = {false, false, false};
    copy = {true, true, true};
    for( vector<ExtField>::iterator extfield=extFields.begin(); extfield!=extFields.end(); extfield++ ) {
        if( extfield->index < allFields.size() ) {
            applyExternalField( allFields[extfield->index], extfield->profile, patch );
        }
        string name = LowerCase( extfield->field );
        if( name==LowerCase( Bx_->name ) ) input[0] = true; 
        if( name==LowerCase( By_->name ) ) input[1] = true; 
        if( name==LowerCase( Bz_->name ) ) input[2] = true; 
        if( name==LowerCase( Bx_m->name ) ) copy[0] = false; 
        if( name==LowerCase( By_m->name ) ) copy[1] = false; 
        if( name==LowerCase( Bz_m->name ) ) copy[2] = false; 
    }
    // These should be additions and not copies. In case B is given but not B_m, B is assumed constant over time and B_m=B.
    if (input[0] && copy[0]) Bx_m->copyFrom( Bx_ );
    if (input[1] && copy[1]) By_m->copyFrom( By_ );
    if (input[2] && copy[2]) Bz_m->copyFrom( Bz_ );
}

void ElectroMagn::saveExternalFields( Patch *patch )
{
    for( vector<ExtField>::iterator extfield=extFields.begin(); extfield!=extFields.end(); extfield++ ) {
        if( extfield->index < allFields.size() ) {
            for( auto &embc: emBoundCond ) {
                if( embc ) {
                    embc->save_fields( allFields[extfield->index], patch );
                }
            }
        }
    }
}

void ElectroMagn::applyPrescribedFields( Patch *patch, double time )
{
    for( vector<PrescribedField>::iterator pf=prescribedFields.begin(); pf!=prescribedFields.end(); pf++ ) {
        if( pf->index < allFields.size() ) {
            pf->savedField->copyFrom(allFields[pf->index]);
            applyPrescribedField( allFields[pf->index], pf->profile, patch, time );
        }
    }
}

void ElectroMagn::resetPrescribedFields()
{
    for( vector<PrescribedField>::iterator pf=prescribedFields.begin(); pf!=prescribedFields.end(); pf++ ) {
        if( pf->index < allFields.size() ) {
            allFields[pf->index]->copyFrom(pf->savedField);
        }
    }
}




void ElectroMagn::applyAntenna( unsigned int iAntenna, double intensity )
{
    Field *field=nullptr;
    Field *antennaField = antennas[iAntenna].field;
    
    if( antennaField && antennas[iAntenna].index<allFields.size() ) {
    
        field = allFields[antennas[iAntenna].index];
        
        for( unsigned int i=0; i< field->number_of_points_ ; i++ ) {
            ( *field )( i ) += intensity * ( *antennaField )( i );
        }
        
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Compute the total density and currents from species density and currents on Device
//! This function is valid wathever the geometry
// ---------------------------------------------------------------------------------------------------------------------
#if defined( SMILEI_ACCELERATOR_GPU )
void ElectroMagn::computeTotalRhoJOnDevice()
{

    double *const __restrict__ Jxp = Jx_->data();
    double *const __restrict__ Jyp = Jy_->data();
    double *const __restrict__ Jzp = Jz_->data();
    double *const __restrict__ rhop = rho_->data();

    unsigned int Jx_size = Jx_->number_of_points_;
    unsigned int Jy_size = Jy_->number_of_points_;
    unsigned int Jz_size = Jz_->number_of_points_;
    unsigned int rho_size = rho_->number_of_points_;

    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {

        double *const __restrict__ Jxsp = Jx_s[ispec] ? Jx_s[ispec]->data() : nullptr;
        double *const __restrict__ Jysp = Jy_s[ispec] ? Jy_s[ispec]->data() : nullptr;
        double *const __restrict__ Jzsp = Jz_s[ispec] ? Jz_s[ispec]->data() : nullptr;
        double *const __restrict__ rhosp = rho_s[ispec] ? rho_s[ispec]->data() : nullptr;


#if defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc parallel present( \
                                          Jxp[0:Jx_size],     \
                                          Jyp[0:Jy_size],     \
                                          Jzp[0:Jz_size],     \
                                          rhop[0:rho_size],   \
                                          Jxsp[0:Jx_size],    \
                                          Jysp[0:Jy_size],    \
                                          Jzsp[0:Jz_size],    \
                                          rhosp[0:rho_size]   \
                                          )  
            {
#endif
        if (Jxsp) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<Jx_size; i++ ) {
                Jxp[i] += Jxsp[i];
            }
        }
        if (Jysp) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<Jy_size; i++ ) {
                Jyp[i] += Jysp[i];
            }
        }
        if (Jzsp) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<Jz_size; i++ ) {
                Jzp[i] += Jzsp[i];
            }
        }
        if (rhosp) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
            #pragma omp target
            #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
            #pragma acc loop gang worker vector
#endif
            for( unsigned int i=0 ; i<rho_size; i++ ) {
                rhop[i] += rhosp[i];
            }
        }
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
        } // end parallel region
#endif

        //smilei::tools::gpu::HostDeviceMemoryManagement::CopyDeviceToHost( rhosp, rho_size );

        //double sum = 0;
        //for (int i = 0 ; i < rho_size ; i++) {
        //    sum += rhosp[i];
        //}
        //std::cerr << "sum rhos"<<ispec<<" in total: " << sum << std::endl;

        //cerr << Jxsp << " " << Jysp << " " << " " << Jzsp << " " << rhosp << std::endl;
    } // end loop species

    //smilei::tools::gpu::HostDeviceMemoryManagement::CopyDeviceToHost( Jyp, Jy_size );

    //double sum = 0;
    //for (int i = 0 ; i < rho_size ; i++) {
    //   sum += rhop[i];
    //}
    //std::cerr << "sum rho in total: " << sum << std::endl;

    //cerr << "end:computeTotalRho" << std::endl;
} //END computeTotalRhoJOnDevice
#endif
