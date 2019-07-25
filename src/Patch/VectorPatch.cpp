#include "VectorPatch.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <math.h>
//#include <string>

#include "Collisions.h"
#include "DomainDecompositionFactory.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"
#include "PeekAtSpecies.h"
#include "SimWindow.h"
#include "SolverFactory.h"
#include "DiagnosticFactory.h"
#include "LaserEnvelope.h"

#include "SyncVectorPatch.h"
#include "interface.h"
#include "Timers.h"

using namespace std;


VectorPatch::VectorPatch()
{
    domain_decomposition_ = NULL ;
}


VectorPatch::VectorPatch( Params &params )
{
    domain_decomposition_ = DomainDecompositionFactory::create( params );
}


VectorPatch::~VectorPatch()
{
    if( domain_decomposition_ != NULL ) {
        delete domain_decomposition_;
    }
}


void VectorPatch::close( SmileiMPI *smpiData )
{
    closeAllDiags( smpiData );


    if( diag_timers.size() ) {
        MESSAGE( "\n\tDiagnostics profile :" );
    }
    for( unsigned int idiag = 0 ;  idiag < diag_timers.size() ; idiag++ ) {
        double sum( 0 );
        MPI_Reduce( &diag_timers[idiag]->time_acc_, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        MESSAGE( "\t\t" << setw( 20 ) << diag_timers[idiag]->name_ << "\t" << sum/( double )smpiData->getSize() );
    }

    for( unsigned int idiag = 0 ;  idiag < diag_timers.size() ; idiag++ ) {
        delete diag_timers[idiag];
    }
    diag_timers.clear();


    for( unsigned int idiag=0 ; idiag<localDiags.size(); idiag++ ) {
        delete localDiags[idiag];
    }
    localDiags.clear();

    for( unsigned int idiag=0 ; idiag<globalDiags.size(); idiag++ ) {
        delete globalDiags[idiag];
    }
    globalDiags.clear();

    for( unsigned int ipatch=0 ; ipatch<size(); ipatch++ ) {
        delete patches_[ipatch];
    }

    patches_.clear();
}

void VectorPatch::createDiags( Params &params, SmileiMPI *smpi, OpenPMDparams &openPMD )
{
    globalDiags = DiagnosticFactory::createGlobalDiagnostics( params, smpi, *this );
    localDiags  = DiagnosticFactory::createLocalDiagnostics( params, smpi, *this, openPMD );

    // Delete all unused fields
    for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
        if( params.geometry!="AMcylindrical" ) {
            for( unsigned int ifield=0 ; ifield<( *this )( ipatch )->EMfields->Jx_s.size(); ifield++ ) {
                if( ( *this )( ipatch )->EMfields->Jx_s[ifield]->data_ == NULL ) {
                    delete( *this )( ipatch )->EMfields->Jx_s[ifield];
                    ( *this )( ipatch )->EMfields->Jx_s[ifield]=NULL;
                }

            }
            for( unsigned int ifield=0 ; ifield<( *this )( ipatch )->EMfields->Jy_s.size(); ifield++ ) {
                if( ( *this )( ipatch )->EMfields->Jy_s[ifield]->data_ == NULL ) {
                    delete( *this )( ipatch )->EMfields->Jy_s[ifield];
                    ( *this )( ipatch )->EMfields->Jy_s[ifield]=NULL;
                }
            }
            for( unsigned int ifield=0 ; ifield<( *this )( ipatch )->EMfields->Jz_s.size(); ifield++ ) {
                if( ( *this )( ipatch )->EMfields->Jz_s[ifield]->data_ == NULL ) {
                    delete( *this )( ipatch )->EMfields->Jz_s[ifield];
                    ( *this )( ipatch )->EMfields->Jz_s[ifield]=NULL;
                }
            }
            for( unsigned int ifield=0 ; ifield<( *this )( ipatch )->EMfields->rho_s.size(); ifield++ ) {
                if( ( *this )( ipatch )->EMfields->rho_s[ifield]->data_ == NULL ) {
                    delete( *this )( ipatch )->EMfields->rho_s[ifield];
                    ( *this )( ipatch )->EMfields->rho_s[ifield]=NULL;
                }
            }

        } else {
            ElectroMagnAM *EMfields = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
            for( unsigned int ifield=0 ; ifield<EMfields->Jl_s.size(); ifield++ ) {
                if( EMfields->Jl_s[ifield]->cdata_ == NULL ) {
                    delete EMfields->Jl_s[ifield];
                    EMfields->Jl_s[ifield]=NULL;
                }
            }
            for( unsigned int ifield=0 ; ifield<EMfields->Jr_s.size(); ifield++ ) {
                if( EMfields->Jr_s[ifield]->cdata_ == NULL ) {
                    delete EMfields->Jr_s[ifield];
                    EMfields->Jr_s[ifield]=NULL;
                }
            }
            for( unsigned int ifield=0 ; ifield<EMfields->Jt_s.size(); ifield++ ) {
                if( EMfields->Jt_s[ifield]->cdata_ == NULL ) {
                    delete EMfields->Jt_s[ifield];
                    EMfields->Jt_s[ifield]=NULL;
                }
            }

            for( unsigned int ifield=0 ; ifield<EMfields->rho_AM_s.size(); ifield++ ) {
                if( EMfields->rho_AM_s[ifield]->cdata_ == NULL ) {
                    delete EMfields->rho_AM_s[ifield];
                    EMfields->rho_AM_s[ifield]=NULL;
                }
            }

        }


        if( params.Laser_Envelope_model ) {
            for( unsigned int ifield=0 ; ifield<( *this )( ipatch )->EMfields->Env_Chi_s.size(); ifield++ ) {
                if( ( *this )( ipatch )->EMfields->Env_Chi_s[ifield]->data_ == NULL ) {
                    delete( *this )( ipatch )->EMfields->Env_Chi_s[ifield];
                    ( *this )( ipatch )->EMfields->Env_Chi_s[ifield]=NULL;
                }
            }
        }



    }

    for( unsigned int idiag = 0 ;  idiag < globalDiags.size() ; idiag++ ) {
        diag_timers.push_back( new Timer( globalDiags[idiag]->filename ) );
    }
    for( unsigned int idiag = 0 ;  idiag < localDiags.size() ; idiag++ ) {
        diag_timers.push_back( new Timer( localDiags[idiag]->filename ) );
    }

    for( unsigned int idiag = 0 ;  idiag < diag_timers.size() ; idiag++ ) {
        diag_timers[idiag]->init( smpi );
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------       INTERFACES        ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------------------------------
// Configure all patches for the new time step
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::configuration( Params &params, Timers &timers, int itime )
{

    //if (params.has_adaptive_vectorization)
    //{

    timers.reconfiguration.restart();

    unsigned int npatches = this->size();

    // Clean buffers
    #pragma omp master
    {
        for( unsigned int ipatch=0 ; ipatch< npatches; ipatch++ ) {
            // For all species
            for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
                ( *this )( ipatch )->cleanMPIBuffers( ispec, params );
            }
        }
    }
    #pragma omp barrier

    // Species configuration according to the default mode

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<npatches ; ipatch++ ) {
        // Particle importation for all species
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            species( ipatch, ispec )->initial_configuration( params, ( *this )( ipatch ) );
        }
    }

    timers.reconfiguration.update( params.printNow( itime ) );
    //}

}

// ---------------------------------------------------------------------------------------------------------------------
// Reconfigure all patches for the new time step
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::reconfiguration( Params &params, Timers &timers, int itime )
{
    //if (params.has_adaptive_vectorization)
    //{

    timers.reconfiguration.restart();

    unsigned int npatches = this->size();

    // Clean buffers
    #pragma omp master
    {
        for( unsigned int ipatch=0 ; ipatch < npatches ; ipatch++ ) {
            // For all species
            for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
                ( *this )( ipatch )->cleanMPIBuffers( ispec, params );
            }
        }
    }
    #pragma omp barrier

    // Species reconfiguration for best performance
    // Change the status to use vectorized or not-vectorized operators
    // as a function of the metrics
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch < npatches ; ipatch++ ) {
        // Particle importation for all species
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            species( ipatch, ispec )->reconfiguration( params, ( *this )( ipatch ) );
        }
    }

    timers.reconfiguration.update( params.printNow( itime ) );
    //}
}


// ---------------------------------------------------------------------------------------------------------------------
// Reconfigure all patches for the new time step
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::sort_all_particles( Params &params )
{
#ifdef _VECTO
    if((  params.vectorization_mode != "off" ) || (params.cell_sorting) ) {
        //Need to sort because particles are not well sorted at creation
        for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
            for( unsigned int ispec=0 ; ispec<patches_[ipatch]->vecSpecies.size(); ispec++ ) {
                patches_[ipatch]->vecSpecies[ispec]->compute_part_cell_keys( params );
                patches_[ipatch]->vecSpecies[ispec]->sort_part( params );
            }
        }
    }
#endif
}

// ---------------------------------------------------------------------------------------------------------------------
// For all patches, move particles (restartRhoJ(s), dynamics and exchangeParticles)
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::dynamics( Params &params,
                            SmileiMPI *smpi,
                            SimWindow *simWindow,
                            RadiationTables &RadiationTables,
                            MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                            double time_dual, Timers &timers, int itime )
{

    #pragma omp single
    diag_flag = needsRhoJsNow( itime );

    timers.particles.restart();
    ostringstream t;
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->restartRhoJ();
        //MESSAGE("restart rhoj");
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            Species *spec = species( ipatch, ispec );
            if( spec->ponderomotive_dynamics ) {
                continue;
            }
            if( spec->isProj( time_dual, simWindow ) || diag_flag ) {
                // Dynamics with vectorized operators
                if( spec->vectorized_operators || params.cell_sorting ) {
                    spec->dynamics( time_dual, ispec,
                                    emfields( ipatch ),
                                    params, diag_flag, partwalls( ipatch ),
                                    ( *this )( ipatch ), smpi,
                                    RadiationTables,
                                    MultiphotonBreitWheelerTables,
                                    localDiags );
                }
                // Dynamics with scalar operators
                else {
                    if( params.vectorization_mode == "adaptive" ) {
                        spec->scalar_dynamics( time_dual, ispec,
                                               emfields( ipatch ),
                                               params, diag_flag, partwalls( ipatch ),
                                               ( *this )( ipatch ), smpi,
                                               RadiationTables,
                                               MultiphotonBreitWheelerTables,
                                               localDiags );
                    } else {
                        spec->Species::dynamics( time_dual, ispec,
                                                 emfields( ipatch ),
                                                 params, diag_flag, partwalls( ipatch ),
                                                 ( *this )( ipatch ), smpi,
                                                 RadiationTables,
                                                 MultiphotonBreitWheelerTables,
                                                 localDiags );
                    }
                } // end if condition on envelope dynamics
            } // end if condition on species
        } // end loop on species
        //MESSAGE("species dynamics");
    } // end loop on patches


    timers.particles.update( params.printNow( itime ) );
#ifdef __DETAILED_TIMERS
    timers.interpolator.update( *this, params.printNow( itime ) );
    timers.pusher.update( *this, params.printNow( itime ) );
    timers.projector.update( *this, params.printNow( itime ) );
    timers.cell_keys.update( *this, params.printNow( itime ) );
    timers.ionization.update( *this, params.printNow( itime ) );
    timers.radiation.update( *this, params.printNow( itime ) );
    timers.multiphoton_Breit_Wheeler_timer.update( *this, params.printNow( itime ) );
#endif

    timers.syncPart.restart();
    for( unsigned int ispec=0 ; ispec<( *this )( 0 )->vecSpecies.size(); ispec++ ) {
        Species *spec = species( 0, ispec );
        if( !spec->ponderomotive_dynamics && spec->isProj( time_dual, simWindow ) ) {
            SyncVectorPatch::exchangeParticles( ( *this ), ispec, params, smpi, timers, itime ); // Included sort_part
        } // end condition on species
    } // end loop on species
    //MESSAGE("exchange particles");
    timers.syncPart.update( params.printNow( itime ) );
#ifdef __DETAILED_TIMERS
    timers.sorting.update( *this, params.printNow( itime ) );
    timers.merging.update( *this, params.printNow( itime ) );
#endif
} // END dynamics

// ---------------------------------------------------------------------------------------------------------------------
// For all patches, project charge and current densities with standard scheme for diag purposes at t=0
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::projection_for_diags( Params &params,
                                        SmileiMPI *smpi,
                                        SimWindow *simWindow,
                                        double time_dual, Timers &timers, int itime )
{

    #pragma omp single
    diag_flag = needsRhoJsNow( itime );

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->restartRhoJ();
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            if( ( *this )( ipatch )->vecSpecies[ispec]->isProj( time_dual, simWindow ) || diag_flag ) {
                species( ipatch, ispec )->projection_for_diags( time_dual, ispec,
                        emfields( ipatch ),
                        params, diag_flag,
                        ( *this )( ipatch ), smpi );
            }
        }

    }

    // if Envelope is used, project the susceptibility of the particles interacting with the envelope
    if( params.Laser_Envelope_model ) {
        #pragma omp for schedule(runtime)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->restartEnvChi();
            for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
                if( ( *this )( ipatch )->vecSpecies[ispec]->isProj( time_dual, simWindow ) || diag_flag ) {
                    if( species( ipatch, ispec )->ponderomotive_dynamics ) {
                        species( ipatch, ispec )->ponderomotive_project_susceptibility( time_dual, ispec,
                                emfields( ipatch ),
                                params, diag_flag,
                                ( *this )( ipatch ), smpi,
                                localDiags );
                    } // end condition on ponderomotive dynamics
                } // end diagnostic or projection if condition on species
            } // end loop on species
        } // end loop on patches
    }

} // END projection for diags

void VectorPatch::finalize_and_sort_parts( Params &params, SmileiMPI *smpi, SimWindow *simWindow,
        double time_dual, Timers &timers, int itime )
{
    timers.syncPart.restart();


    // Particle synchronization and sorting
    // ----------------------------------------

    for( unsigned int ispec=0 ; ispec<( *this )( 0 )->vecSpecies.size(); ispec++ ) {
        if( ( *this )( 0 )->vecSpecies[ispec]->isProj( time_dual, simWindow ) ) {
            SyncVectorPatch::finalize_and_sort_parts( ( *this ), ispec, params, smpi, timers, itime ); // Included sort_part
        }

    }

    // Particle importation
    // ----------------------------------------

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        // Particle importation for all species
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            if( ( *this )( ipatch )->vecSpecies[ispec]->isProj( time_dual, simWindow ) || diag_flag ) {
                species( ipatch, ispec )->dynamics_import_particles( time_dual, ispec,
                        params,
                        ( *this )( ipatch ), smpi,
                        localDiags );
            }
        }
    }

    // Particle merging
    // ---------------------------------------

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        // Particle importation for all species
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            // Check if the particle merging is activated for this species
            if (species( ipatch, ispec )->has_merging_) {

                // Check the time selection
                if( species( ipatch, ispec )->merging_time_selection_->theTimeIsNow( itime ) ) {
                    species( ipatch, ispec )->mergeParticles( time_dual, ispec,
                            params,
                            ( *this )( ipatch ), smpi,
                            localDiags );
                }
            }
        }
    }


    // Species reconfiguration for best performance
    // Change the status to use vectorized or not-vectorized operators
    // as a function of the metrics
    /*#pragma omp for schedule(runtime)
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
        // Particle importation for all species
        for (unsigned int ispec=0 ; ispec<(*this)(ipatch)->vecSpecies.size() ; ispec++) {
            species(ipatch, ispec)->reconfiguration(params, (*this)(ipatch));
        }
    }*/

    if( itime%params.every_clean_particles_overhead==0 ) {
        #pragma omp master
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->cleanParticlesOverhead( params );
        }
        #pragma omp barrier
    }

    timers.syncPart.update( params.printNow( itime ) );

} // END finalize_and_sort_parts


void VectorPatch::computeCharge()
{
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->restartRhoJ();
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            if( ( *this )( ipatch )->vecSpecies[ispec]->vectorized_operators ) {
                species( ipatch, ispec )->computeCharge( ispec, emfields( ipatch ) );
            } else {
                species( ipatch, ispec )->Species::computeCharge( ispec, emfields( ipatch ) );
            }
        }
    }

} // END computeRho

void VectorPatch::computeChargeRelativisticSpecies( double time_primal )
{
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->restartRhoJ();
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            // project only if species needs relativistic initialization and it is the right time to initialize its fields
            if( ( species( ipatch, ispec )->relativistic_field_initialization ) &&
                    ( time_primal == species( ipatch, ispec )->time_relativistic_initialization ) ) {
                if( ( *this )( ipatch )->vecSpecies[ispec]->vectorized_operators ) {
                    species( ipatch, ispec )->computeCharge( ispec, emfields( ipatch ) );
                } else {
                    species( ipatch, ispec )->Species::computeCharge( ispec, emfields( ipatch ) );
                }
            }
        }
    }
} // END computeRho

void VectorPatch::resetRhoJ()
{
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->restartRhoJ();
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// For all patch, sum densities on ghost cells (sum per species if needed, sync per patch and MPI sync)
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::sumDensities( Params &params, double time_dual, Timers &timers, int itime, SimWindow *simWindow, SmileiMPI *smpi )
{
    bool some_particles_are_moving = false;
    unsigned int n_species( ( *this )( 0 )->vecSpecies.size() );
    for( unsigned int ispec=0 ; ispec < n_species ; ispec++ ) {
        if( ( *this )( 0 )->vecSpecies[ispec]->isProj( time_dual, simWindow ) ) {
            some_particles_are_moving = true;
        }
    }
    if( !some_particles_are_moving  && !diag_flag ) {
        return;
    }

    timers.densities.restart();
    if( diag_flag ) {
        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // Per species in global, Attention if output -> Sync / per species fields
            ( *this )( ipatch )->EMfields->computeTotalRhoJ();
        }
    } //MESSAGE ("bug before");
    timers.densities.update();

    timers.syncDens.restart();
    if( params.geometry != "AMcylindrical" ) {
        SyncVectorPatch::sumRhoJ( params, ( *this ), smpi, timers, itime ); // MPI
    } else {
        for( unsigned int imode = 0 ; imode < static_cast<ElectroMagnAM *>( patches_[0]->EMfields )->Jl_.size() ; imode++ ) {
            SyncVectorPatch::sumRhoJ( params, ( *this ), imode, smpi, timers, itime );
        }
    }
    //MESSAGE ("bug after");
    if( diag_flag ) {
        for( unsigned int ispec=0 ; ispec<( *this )( 0 )->vecSpecies.size(); ispec++ ) {
            if( !( *this )( 0 )->vecSpecies[ispec]->particles->is_test ) {
                update_field_list( ispec, smpi );
                if( params.geometry != "AMcylindrical" ) {
                    SyncVectorPatch::sumRhoJs( params, ( *this ), ispec, smpi, timers, itime ); // MPI
                } else {
                    for( unsigned int imode = 0 ; imode < static_cast<ElectroMagnAM *>( patches_[0]->EMfields )->Jl_.size() ; imode++ ) {
                        SyncVectorPatch::sumRhoJs( params, ( *this ), imode, ispec, smpi, timers, itime );
                    }
                }
            }
        }
    }
    timers.syncDens.update( params.printNow( itime ) );
} // End sumDensities


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::sumSusceptibility( Params &params, double time_dual, Timers &timers, int itime, SimWindow *simWindow, SmileiMPI *smpi )
{
    bool some_particles_are_moving = false;
    unsigned int n_species( ( *this )( 0 )->vecSpecies.size() );
    for( unsigned int ispec=0 ; ispec < n_species ; ispec++ ) {
        if( ( *this )( 0 )->vecSpecies[ispec]->isProj( time_dual, simWindow ) ) {
            some_particles_are_moving = true;
        }
    }
    if( !some_particles_are_moving  && !diag_flag ) {
        return;
    }

    timers.susceptibility.restart();
    if( diag_flag ) {
        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // Per species in global, Attention if output -> Sync / per species fields
            ( *this )( ipatch )->EMfields->computeTotalEnvChi();
        }
    }

    if( diag_flag ) {
        for( unsigned int ispec=0 ; ispec<( *this )( 0 )->vecSpecies.size(); ispec++ ) {
            if( !( *this )( 0 )->vecSpecies[ispec]->particles->is_test ) {
                if( species( 0, ispec )->ponderomotive_dynamics ) {
                    update_field_list( ispec, smpi );
                    SyncVectorPatch::sumEnvChis( params, ( *this ), ispec, smpi, timers, itime );
                } // MPI
            }
        }
    }

    timers.susceptibility.update();

    timers.susceptibility.restart();
    if( ( params.geometry == "1Dcartesian" ) or ( params.geometry == "2Dcartesian" ) or ( params.geometry == "3Dcartesian" ) ) {
        SyncVectorPatch::sumEnvChi( params, ( *this ), smpi, timers, itime ); // MPI
    } else {
        ERROR( "Envelope model not yet implemented in this geometry" );
        // for (unsigned int imode = 0 ; imode < static_cast<ElectroMagnAM*>(patches_[0]->EMfields)->Jl_.size() ; imode++  ) {
        //     SyncVectorPatch::sumRhoJ( params, (*this), imode, timers, itime );
        // }
    }

    timers.susceptibility.update();

} // End sumSusceptibility



// ---------------------------------------------------------------------------------------------------------------------
// For all patch, update E and B (Ampere, Faraday, boundary conditions, exchange B and center B)
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::solveMaxwell( Params &params, SimWindow *simWindow, int itime, double time_dual, Timers &timers, SmileiMPI *smpi )
{
    timers.maxwell.restart();

    for( unsigned int ipassfilter=0 ; ipassfilter<params.currentFilter_passes ; ipassfilter++ ) {
        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // Current spatial filtering
            ( *this )( ipatch )->EMfields->binomialCurrentFilter();
        }
        if (params.geometry != "AMcylindrical"){
            SyncVectorPatch::exchange_along_all_directions<double,Field>( listJx_, *this, smpi );
            SyncVectorPatch::finalize_exchange_along_all_directions( listJx_, *this );
            SyncVectorPatch::exchange_along_all_directions<double,Field>( listJy_, *this, smpi );
            SyncVectorPatch::finalize_exchange_along_all_directions( listJy_, *this );
            SyncVectorPatch::exchange_along_all_directions<double,Field>( listJz_, *this, smpi );
            SyncVectorPatch::finalize_exchange_along_all_directions( listJz_, *this );
        } else {
            for (unsigned int imode=0 ; imode < params.nmodes; imode++) {
                SyncVectorPatch::exchange_along_all_directions<complex<double>,cField>( listJl_[imode], *this, smpi );
                SyncVectorPatch::finalize_exchange_along_all_directions( listJl_[imode], *this );
                SyncVectorPatch::exchange_along_all_directions<complex<double>,cField>( listJr_[imode], *this, smpi );
                SyncVectorPatch::finalize_exchange_along_all_directions( listJr_[imode], *this );
                SyncVectorPatch::exchange_along_all_directions<complex<double>,cField>( listJt_[imode], *this, smpi );
                SyncVectorPatch::finalize_exchange_along_all_directions( listJt_[imode], *this );
            }
        }
    }
    #pragma omp for schedule(static)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        if( !params.is_spectral ) {
            // Saving magnetic fields (to compute centered fields used in the particle pusher)
            // Stores B at time n in B_m.
            ( *this )( ipatch )->EMfields->saveMagneticFields( params.is_spectral );
        }
        // Computes Ex_, Ey_, Ez_ on all points.
        // E is already synchronized because J has been synchronized before.
        ( *( *this )( ipatch )->EMfields->MaxwellAmpereSolver_ )( ( *this )( ipatch )->EMfields );
        //MESSAGE("SOLVE MAXWELL AMPERE");
        // Computes Bx_, By_, Bz_ at time n+1 on interior points.
        //for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
        ( *( *this )( ipatch )->EMfields->MaxwellFaradaySolver_ )( ( *this )( ipatch )->EMfields );
        //MESSAGE("SOLVE MAXWELL FARADAY");
    }
    //Synchronize B fields between patches.
    timers.maxwell.update( params.printNow( itime ) );


    timers.syncField.restart();
    if( params.geometry != "AMcylindrical" ) {
        if( params.is_spectral ) {
            SyncVectorPatch::exchangeE( params, ( *this ), smpi );
        }
        SyncVectorPatch::exchangeB( params, ( *this ), smpi );
    } else {
        for( unsigned int imode = 0 ; imode < static_cast<ElectroMagnAM *>( patches_[0]->EMfields )->El_.size() ; imode++ ) {
            SyncVectorPatch::exchangeB( params, ( *this ), imode, smpi );
            SyncVectorPatch::finalizeexchangeB( params, ( *this ), imode ); // disable async, because of tags which is the same for all modes
        }
    }
    timers.syncField.update( params.printNow( itime ) );


#ifdef _PICSAR
    //if ( (params.is_spectral) && (itime!=0) && ( time_dual > params.time_fields_frozen ) ) {
    if( ( itime!=0 ) && ( time_dual > params.time_fields_frozen ) ) {
        timers.syncField.restart();
        if( params.is_spectral ) {
            SyncVectorPatch::finalizeexchangeE( params, ( *this ) );
        }

        SyncVectorPatch::finalizeexchangeB( params, ( *this ) );
        timers.syncField.update( params.printNow( itime ) );

        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // Applies boundary conditions on B
            ( *this )( ipatch )->EMfields->boundaryConditions( itime, time_dual, ( *this )( ipatch ), params, simWindow );
            // Computes B at time n using B and B_m.
            if( !params.is_spectral ) {
                ( *this )( ipatch )->EMfields->centerMagneticFields();
            } else {
                ( *this )( ipatch )->EMfields->saveMagneticFields( params.is_spectral );
            }
        }
        if( params.is_spectral ) {
            save_old_rho( params );
        }
    }
#endif


} // END solveMaxwell

void VectorPatch::solveEnvelope( Params &params, SimWindow *simWindow, int itime, double time_dual, Timers &timers, SmileiMPI *smpi )
{

    if( ( *this )( 0 )->EMfields->envelope!=NULL ) {

        timers.envelope.restart();
        // Exchange susceptibility
        SyncVectorPatch::exchangeEnvChi( params, ( *this ), smpi );

        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {

            // Saving Phi and GradPhi fields
            // (to compute centered quantities used in the particle position ponderomotive pusher)
            // Stores Phi at time n in Phi_m, GradPhi at time n in GradPhi_m
            ( *this )( ipatch )->EMfields->envelope->savePhi_and_GradPhi();

            // Computes A in all points
            ( *this )( ipatch )->EMfields->envelope->compute( ( *this )( ipatch )->EMfields );
            ( *this )( ipatch )->EMfields->envelope->boundaryConditions( itime, time_dual, ( *this )( ipatch ), params, simWindow );

            // Compute ponderomotive potential Phi=|A|^2/2
            ( *this )( ipatch )->EMfields->envelope->compute_Phi( ( *this )( ipatch )->EMfields );

        }

        // Exchange envelope A
        SyncVectorPatch::exchangeA( params, ( *this ), smpi );
        SyncVectorPatch::finalizeexchangeA( params, ( *this ) );

        // Exchange Phi
        SyncVectorPatch::exchangePhi( params, ( *this ), smpi );
        SyncVectorPatch::finalizeexchangePhi( params, ( *this ) );


        // Compute gradients of Phi
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->envelope->compute_gradient_Phi( ( *this )( ipatch )->EMfields );
            // Computes Phi and GradPhi at time n+1/2 using their values at timestep n+1 and n (the latter already in Phi_m and GradPhi_m)
            ( *this )( ipatch )->EMfields->envelope->centerPhi_and_GradPhi();
        }

        // Exchange GradPhi
        SyncVectorPatch::exchangeGradPhi( params, ( *this ), smpi );
        SyncVectorPatch::finalizeexchangeGradPhi( params, ( *this ) );
        timers.envelope.update();
    }

} // END solveEnvelope

void VectorPatch::finalize_sync_and_bc_fields( Params &params, SmileiMPI *smpi, SimWindow *simWindow,
        double time_dual, Timers &timers, int itime )
{
#ifndef _PICSAR
    if( ( !params.is_spectral ) && ( itime!=0 ) && ( time_dual > params.time_fields_frozen ) ) {
        if( params.geometry != "AMcylindrical" ) {
            timers.syncField.restart();
            SyncVectorPatch::finalizeexchangeB( params, ( *this ) );
            timers.syncField.update( params.printNow( itime ) );
        }

        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // Applies boundary conditions on B
            ( *this )( ipatch )->EMfields->boundaryConditions( itime, time_dual, ( *this )( ipatch ), params, simWindow );
            // Computes B at time n using B and B_m.
            ( *this )( ipatch )->EMfields->centerMagneticFields();
        }
    }
#endif

} // END finalize_sync_and_bc_fields


void VectorPatch::initExternals( Params &params )
{
    // Init all lasers
    for( unsigned int ipatch=0; ipatch<size(); ipatch++ ) {
        if( ( *this )( ipatch )->isXmin() && ( *this )( ipatch )->EMfields->emBoundCond[0] != NULL ) {
            unsigned int nlaser = ( *this )( ipatch )->EMfields->emBoundCond[0]->vecLaser.size();
            for( unsigned int ilaser = 0; ilaser < nlaser; ilaser++ ) {
                ( *this )( ipatch )->EMfields->emBoundCond[0]->vecLaser[ilaser]->initFields( params, ( *this )( ipatch ) );
            }
        }

        if( ( *this )( ipatch )->isXmax() && ( *this )( ipatch )->EMfields->emBoundCond[1] != NULL ) {
            unsigned int nlaser = ( *this )( ipatch )->EMfields->emBoundCond[1]->vecLaser.size();
            for( unsigned int ilaser = 0; ilaser < nlaser; ilaser++ ) {
                ( *this )( ipatch )->EMfields->emBoundCond[1]->vecLaser[ilaser]->initFields( params, ( *this )( ipatch ) );
            }
        }
    }

    // Init all antennas
    for( unsigned int ipatch=0; ipatch<size(); ipatch++ ) {
        ( *this )( ipatch )->EMfields->initAntennas( ( *this )( ipatch ) );
    }
}


void VectorPatch::initAllDiags( Params &params, SmileiMPI *smpi )
{
    // Global diags: scalars + particles
    for( unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++ ) {
        globalDiags[idiag]->init( params, smpi, *this );
        // MPI master creates the file
        if( smpi->isMaster() ) {
            globalDiags[idiag]->openFile( params, smpi, true );
        }
    }

    // Local diags : fields, probes, tracks
    for( unsigned int idiag = 0 ; idiag < localDiags.size() ; idiag++ ) {
        localDiags[idiag]->init( params, smpi, *this );
    }

} // END initAllDiags


void VectorPatch::closeAllDiags( SmileiMPI *smpi )
{
    // MPI master closes all global diags
    if( smpi->isMaster() )
        for( unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++ ) {
            globalDiags[idiag]->closeFile();
        }

    // All MPI close local diags
    for( unsigned int idiag = 0 ; idiag < localDiags.size() ; idiag++ ) {
        localDiags[idiag]->closeFile();
    }
}


void VectorPatch::openAllDiags( Params &params, SmileiMPI *smpi )
{
    // MPI master opens all global diags
    if( smpi->isMaster() )
        for( unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++ ) {
            globalDiags[idiag]->openFile( params, smpi, false );
        }

    // All MPI open local diags
    for( unsigned int idiag = 0 ; idiag < localDiags.size() ; idiag++ ) {
        localDiags[idiag]->openFile( params, smpi, false );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// For all patch, Compute and Write all diags
//   - Scalars, Probes, Phases, TrackParticles, Fields, Average fields
//   - set diag_flag to 0 after write
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::runAllDiags( Params &params, SmileiMPI *smpi, unsigned int itime, Timers &timers, SimWindow *simWindow )
{
    // Global diags: scalars + particles
    timers.diags.restart();
    for( unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++ ) {
        diag_timers[idiag]->restart();

        #pragma omp single
        globalDiags[idiag]->theTimeIsNow = globalDiags[idiag]->prepare( itime );
        #pragma omp barrier
        if( globalDiags[idiag]->theTimeIsNow ) {
            // All patches run
            #pragma omp for schedule(runtime)
            for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
                globalDiags[idiag]->run( ( *this )( ipatch ), itime, simWindow );
            }
            // MPI procs gather the data and compute
            #pragma omp single
            smpi->computeGlobalDiags( globalDiags[idiag], itime );
            // MPI master writes
            #pragma omp single
            globalDiags[idiag]->write( itime, smpi );
        }

        diag_timers[idiag]->update();
    }

    // Local diags : fields, probes, tracks
    for( unsigned int idiag = 0 ; idiag < localDiags.size() ; idiag++ ) {
        diag_timers[globalDiags.size()+idiag]->restart();

        #pragma omp single
        localDiags[idiag]->theTimeIsNow = localDiags[idiag]->prepare( itime );
        #pragma omp barrier
        // All MPI run their stuff and write out
        if( localDiags[idiag]->theTimeIsNow ) {
            localDiags[idiag]->run( smpi, *this, itime, simWindow, timers );
        }

        diag_timers[globalDiags.size()+idiag]->update();
    }

    // Manage the "diag_flag" parameter, which indicates whether Rho and Js were used
    if( diag_flag ) {
        #pragma omp barrier
        #pragma omp single
        diag_flag = false;
        #pragma omp for
        for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->restartRhoJs();
        }
    }
    timers.diags.update();

} // END runAllDiags


// ---------------------------------------------------------------------------------------------------------------------
// Check if rho is null (MPI & patch sync)
// ---------------------------------------------------------------------------------------------------------------------
bool VectorPatch::isRhoNull( SmileiMPI *smpi )
{
    double norm2( 0. );
    double locnorm2( 0. );
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        locnorm2 += ( *this )( ipatch )->EMfields->computeRhoNorm2();
    }

    MPI_Allreduce( &locnorm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    return ( norm2<=0. );
} // END isRhoNull


// ---------------------------------------------------------------------------------------------------------------------
// Solve Poisson to initialize E
//   - all steps are done locally, sync per patch, sync per MPI process
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::solvePoisson( Params &params, SmileiMPI *smpi )
{
    Timer ptimer( "global" );
    ptimer.init( smpi );
    ptimer.restart();


    unsigned int iteration_max = params.poisson_max_iteration;
    double           error_max = params.poisson_max_error;
    unsigned int iteration=0;

    // Init & Store internal data (phi, r, p, Ap) per patch
    double rnew_dot_rnew_local( 0. );
    double rnew_dot_rnew( 0. );
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->initPoisson( ( *this )( ipatch ) );
        rnew_dot_rnew_local += ( *this )( ipatch )->EMfields->compute_r();
    }
    MPI_Allreduce( &rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    std::vector<Field *> Ex_;
    std::vector<Field *> Ap_;

    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        Ex_.push_back( ( *this )( ipatch )->EMfields->Ex_ );
        Ap_.push_back( ( *this )( ipatch )->EMfields->Ap_ );
    }

    unsigned int nx_p2_global = ( params.n_space_global[0]+1 );
    if( Ex_[0]->dims_.size()>1 ) {
        nx_p2_global *= ( params.n_space_global[1]+1 );
        if( Ex_[0]->dims_.size()>2 ) {
            nx_p2_global *= ( params.n_space_global[2]+1 );
        }
    }

    // compute control parameter
    double ctrl = rnew_dot_rnew / ( double )( nx_p2_global );

    // ---------------------------------------------------------
    // Starting iterative loop for the conjugate gradient method
    // ---------------------------------------------------------
    if( smpi->isMaster() ) {
        DEBUG( "Starting iterative loop for CG method" );
    }
    while( ( ctrl > error_max ) && ( iteration<iteration_max ) ) {
        iteration++;
        if( smpi->isMaster() ) {
            DEBUG( "iteration " << iteration << " started with control parameter ctrl = " << ctrl*1.e14 << " x 1e-14" );
        }

        // scalar product of the residual
        double r_dot_r = rnew_dot_rnew;

        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->compute_Ap( ( *this )( ipatch ) );
        }

        // Exchange Ap_ (intra & extra MPI)
        SyncVectorPatch::exchange_along_all_directions<double,Field>( Ap_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions( Ap_, *this );

        // scalar product p.Ap
        double p_dot_Ap       = 0.0;
        double p_dot_Ap_local = 0.0;
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            p_dot_Ap_local += ( *this )( ipatch )->EMfields->compute_pAp();
        }
        MPI_Allreduce( &p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );


        // compute new potential and residual
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->update_pand_r( r_dot_r, p_dot_Ap );
        }

        // compute new residual norm
        rnew_dot_rnew       = 0.0;
        rnew_dot_rnew_local = 0.0;
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            rnew_dot_rnew_local += ( *this )( ipatch )->EMfields->compute_r();
        }
        MPI_Allreduce( &rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if( smpi->isMaster() ) {
            DEBUG( "new residual norm: rnew_dot_rnew = " << rnew_dot_rnew );
        }

        // compute new directio
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->update_p( rnew_dot_rnew, r_dot_r );
        }

        // compute control parameter
        ctrl = rnew_dot_rnew / ( double )( nx_p2_global );
        if( smpi->isMaster() ) {
            DEBUG( "iteration " << iteration << " done, exiting with control parameter ctrl = " << ctrl );
        }

    }//End of the iterative loop


    // --------------------------------
    // Status of the solver convergence
    // --------------------------------
    if( iteration_max>0 && iteration == iteration_max ) {
        if( smpi->isMaster() )
            WARNING( "Poisson solver did not converge: reached maximum iteration number: " << iteration
                     << ", relative err is ctrl = " << 1.0e14*ctrl << " x 1e-14" );
    } else {
        if( smpi->isMaster() )
            MESSAGE( 1, "Poisson solver converged at iteration: " << iteration
                     << ", relative err is ctrl = " << 1.0e14*ctrl << " x 1e-14" );
    }

    // ------------------------------------------
    // Compute the electrostatic fields Ex and Ey
    // ------------------------------------------
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->initE( ( *this )( ipatch ) );
    }

    SyncVectorPatch::exchangeE( params, *this, smpi );
    SyncVectorPatch::finalizeexchangeE( params, *this );

    // Centering of the electrostatic fields
    // -------------------------------------
    vector<double> E_Add( Ex_[0]->dims_.size(), 0. );
    if( Ex_[0]->dims_.size()==3 ) {
        double Ex_avg_local( 0. ), Ex_avg( 0. ), Ey_avg_local( 0. ), Ey_avg( 0. ), Ez_avg_local( 0. ), Ez_avg( 0. );
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            Ex_avg_local += ( *this )( ipatch )->EMfields->computeExSum();
            Ey_avg_local += ( *this )( ipatch )->EMfields->computeEySum();
            Ez_avg_local += ( *this )( ipatch )->EMfields->computeEzSum();
        }

        MPI_Allreduce( &Ex_avg_local, &Ex_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &Ey_avg_local, &Ey_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &Ez_avg_local, &Ez_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        E_Add[0] = -Ex_avg/( ( params.n_space[0]+2 )*( params.n_space[1]+1 )*( params.n_space[2]+1 ) );
        E_Add[1] = -Ey_avg/( ( params.n_space[0]+1 )*( params.n_space[1]+2 )*( params.n_space[2]+1 ) );;
        E_Add[2] = -Ez_avg/( ( params.n_space[0]+1 )*( params.n_space[1]+1 )*( params.n_space[2]+2 ) );;
    } else if( Ex_[0]->dims_.size()==2 ) {
        double Ex_XminYmax = 0.0;
        double Ey_XminYmax = 0.0;
        double Ex_XmaxYmin = 0.0;
        double Ey_XmaxYmin = 0.0;

        //The YmaxXmin patch has Patch coordinates X=0, Y=2^m1-1= number_of_patches[1]-1.
        std::vector<int> xcall( 2, 0 );
        xcall[0] = 0;
        xcall[1] = params.number_of_patches[1]-1;
        int patch_YmaxXmin = domain_decomposition_->getDomainId( xcall );
        //The MPI rank owning it is
        int rank_XminYmax = smpi->hrank( patch_YmaxXmin );
        //The YminXmax patch has Patch coordinates X=2^m0-1= number_of_patches[0]-1, Y=0.
        //Its hindex is
        xcall[0] = params.number_of_patches[0]-1;
        xcall[1] = 0;
        int patch_YminXmax = domain_decomposition_->getDomainId( xcall );
        //The MPI rank owning it is
        int rank_XmaxYmin = smpi->hrank( patch_YminXmax );


        //cout << patch_YmaxXmin << " " << rank_XminYmax << " " << patch_YminXmax << " " << rank_XmaxYmin << endl;

        if( smpi->getRank() == rank_XminYmax ) {
            Ex_XminYmax = ( *this )( patch_YmaxXmin-( this->refHindex_ ) )->EMfields->getEx_XminYmax();
            Ey_XminYmax = ( *this )( patch_YmaxXmin-( this->refHindex_ ) )->EMfields->getEy_XminYmax();
        }

        // Xmax-Ymin corner
        if( smpi->getRank() == rank_XmaxYmin ) {
            Ex_XmaxYmin = ( *this )( patch_YminXmax-( this->refHindex_ ) )->EMfields->getEx_XmaxYmin();
            Ey_XmaxYmin = ( *this )( patch_YminXmax-( this->refHindex_ ) )->EMfields->getEy_XmaxYmin();
        }

        MPI_Bcast( &Ex_XminYmax, 1, MPI_DOUBLE, rank_XminYmax, MPI_COMM_WORLD );
        MPI_Bcast( &Ey_XminYmax, 1, MPI_DOUBLE, rank_XminYmax, MPI_COMM_WORLD );

        MPI_Bcast( &Ex_XmaxYmin, 1, MPI_DOUBLE, rank_XmaxYmin, MPI_COMM_WORLD );
        MPI_Bcast( &Ey_XmaxYmin, 1, MPI_DOUBLE, rank_XmaxYmin, MPI_COMM_WORLD );

        //This correction is always done, independantly of the periodicity. Is this correct ?
        E_Add[0] = -0.5*( Ex_XminYmax+Ex_XmaxYmin );
        E_Add[1] = -0.5*( Ey_XminYmax+Ey_XmaxYmin );

#ifdef _3D_LIKE_CENTERING
        double Ex_avg_local( 0. ), Ex_avg( 0. ), Ey_avg_local( 0. ), Ey_avg( 0. );
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            Ex_avg_local += ( *this )( ipatch )->EMfields->computeExSum();
            Ey_avg_local += ( *this )( ipatch )->EMfields->computeEySum();
        }

        MPI_Allreduce( &Ex_avg_local, &Ex_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &Ey_avg_local, &Ey_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        E_Add[0] = -Ex_avg/( ( params.n_space[0]+2 )*( params.n_space[1]+1 ) );
        E_Add[1] = -Ey_avg/( ( params.n_space[0]+1 )*( params.n_space[1]+2 ) );;
#endif

    } else if( Ex_[0]->dims_.size()==1 ) {
        double Ex_Xmin = 0.0;
        double Ex_Xmax = 0.0;

        unsigned int rankXmin = 0;
        if( smpi->getRank() == 0 ) {
            //Ex_Xmin = (*Ex1D)(index_bc_min[0]);
            Ex_Xmin = ( *this )( ( 0 )-( this->refHindex_ ) )->EMfields->getEx_Xmin();
        }
        MPI_Bcast( &Ex_Xmin, 1, MPI_DOUBLE, rankXmin, MPI_COMM_WORLD );

        unsigned int rankXmax = smpi->getSize()-1;
        if( smpi->getRank() == smpi->getSize()-1 ) {
            //Ex_Xmax = (*Ex1D)(index_bc_max[0]);
            Ex_Xmax = ( *this )( ( params.number_of_patches[0]-1 )-( this->refHindex_ ) )->EMfields->getEx_Xmax();
        }
        MPI_Bcast( &Ex_Xmax, 1, MPI_DOUBLE, rankXmax, MPI_COMM_WORLD );
        E_Add[0] = -0.5*( Ex_Xmin+Ex_Xmax );

#ifdef _3D_LIKE_CENTERING
        double Ex_avg_local( 0. ), Ex_avg( 0. );
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            Ex_avg_local += ( *this )( ipatch )->EMfields->computeExSum();
        }

        MPI_Allreduce( &Ex_avg_local, &Ex_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        E_Add[0] = -Ex_avg/( ( params.n_space[0]+2 ) );
#endif

    }

    // Centering electrostatic fields
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->centeringE( E_Add );
    }


    // Compute error on the Poisson equation
    double deltaPoisson_max = 0.0;
    int i_deltaPoisson_max  = -1;

#ifdef _A_FINALISER
    for( unsigned int i=0; i<nx_p; i++ ) {
        double deltaPoisson = abs( ( ( *Ex1D )( i+1 )-( *Ex1D )( i ) )/dx - ( *rho1D )( i ) );
        if( deltaPoisson > deltaPoisson_max ) {
            deltaPoisson_max   = deltaPoisson;
            i_deltaPoisson_max = i;
        }
    }
#endif

    //!\todo Reduce to find global max
    if( smpi->isMaster() ) {
        MESSAGE( 1, "Poisson equation solved. Maximum err = " << deltaPoisson_max << " at i= " << i_deltaPoisson_max );
    }

    ptimer.update();
    MESSAGE( "Time in Poisson : " << ptimer.getTime() );

} // END solvePoisson


void VectorPatch::runRelativisticModule( double time_prim, Params &params, SmileiMPI* smpi,  Timers &timers )
{
    // Compute rho only for species needing relativistic field Initialization
    computeChargeRelativisticSpecies( time_prim );

    if (params.geometry != "AMcylindrical"){
        SyncVectorPatch::sum<double,Field>( listrho_, (*this), smpi, timers, 0 );
    } else {
        for( unsigned int imode=0 ; imode<params.nmodes ; imode++ ) {
            SyncVectorPatch::sumRhoJ( params, (*this), imode, smpi, timers, 0 );
        }
    }

    #pragma omp master
    {
        // Initialize the fields for these species
        if( !isRhoNull( smpi ) ) {
            TITLE( "Initializing relativistic species fields" );
            if (params.geometry != "AMcylindrical"){
                solveRelativisticPoisson( params, smpi, time_prim );
            } else {
                solveRelativisticPoissonAM( params, smpi, time_prim );
            }
        }
    }
    #pragma omp barrier

    // Reset rho and J and return to PIC loop
    resetRhoJ();

}


void VectorPatch::solveRelativisticPoisson( Params &params, SmileiMPI *smpi, double time_primal )
{


    //Timer ptimer("global");
    //ptimer.init(smpi);
    //ptimer.restart();

    // Assumption: one or more species move in vacuum with mean lorentz gamma factor gamma_mean in the x direction,
    // with low energy spread.
    // The electromagnetic fields of this species can be initialized solving a Poisson-like problem (here informally
    // referred to as "relativistic Poisson problem") and then performing a Lorentz back-transformation to find the
    // electromagnetic fields of the species in the lab frame.
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    // In case of non-monoenergetic relativistic distribution (NOT IMPLEMENTED AT THE MOMENT), the linearity of Maxwell's equations can be exploited:
    // divide the species in quasi-monoenergetic bins with gamma_i and repeat the same procedure for described above
    // for all bins. Finally, in the laboratory frame sum all the fields of the various energy-bin ensembles of particles.

    // All the parameters for the Poisson problem (e.g. maximum iteration) are the same used in the namelist
    // for the traditional Poisson problem

    // compute gamma_mean for the species for which the field is initialized
    double s_gamma( 0. );
    uint64_t nparticles( 0 );
    for( unsigned int ispec=0 ; ispec<( *this )( 0 )->vecSpecies.size() ; ispec++ ) {
        if( species( 0, ispec )->relativistic_field_initialization ) {
            for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
                if( time_primal==species( ipatch, ispec )->time_relativistic_initialization ) {
                    s_gamma += species( ipatch, ispec )->sum_gamma();
                    nparticles += species( ipatch, ispec )->getNbrOfParticles();
                }
            }
        }
    }
    double gamma_global( 0. );
    MPI_Allreduce( &s_gamma, &gamma_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    uint64_t nparticles_global( 0 );
    MPI_Allreduce( &nparticles, &nparticles_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
    MESSAGE( "GAMMA = " << gamma_global/( double )nparticles_global );

    //Timer ptimer("global");
    //ptimer.init(smpi);
    //ptimer.restart();

    double gamma_mean = gamma_global/( double )nparticles_global;

    unsigned int iteration_max = params.relativistic_poisson_max_iteration;
    double           error_max = params.relativistic_poisson_max_error;
    unsigned int iteration=0;

    // Init & Store internal data (phi, r, p, Ap) per patch
    double rnew_dot_rnew_local( 0. );
    double rnew_dot_rnew( 0. );
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->initPoisson( ( *this )( ipatch ) );
        rnew_dot_rnew_local += ( *this )( ipatch )->EMfields->compute_r();
        //cout << std::scientific << "rnew_dot_rnew_local = " << rnew_dot_rnew_local << endl;
        ( *this )( ipatch )->EMfields->initRelativisticPoissonFields( ( *this )( ipatch ) );
    }
    //cout << std::scientific << "rnew_dot_rnew_local = " << rnew_dot_rnew_local << endl;
    MPI_Allreduce( &rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    std::vector<Field *> Ex_;
    std::vector<Field *> Ey_;
    std::vector<Field *> Ez_;
    std::vector<Field *> Bx_;
    std::vector<Field *> By_;
    std::vector<Field *> Bz_;
    std::vector<Field *> Bx_m;
    std::vector<Field *> By_m;
    std::vector<Field *> Bz_m;

    std::vector<Field *> Ex_rel_;
    std::vector<Field *> Ey_rel_;
    std::vector<Field *> Ez_rel_;
    std::vector<Field *> Bx_rel_;
    std::vector<Field *> By_rel_;
    std::vector<Field *> Bz_rel_;

    std::vector<Field *> Bx_rel_t_plus_halfdt_;
    std::vector<Field *> By_rel_t_plus_halfdt_;
    std::vector<Field *> Bz_rel_t_plus_halfdt_;
    std::vector<Field *> Bx_rel_t_minus_halfdt_;
    std::vector<Field *> By_rel_t_minus_halfdt_;
    std::vector<Field *> Bz_rel_t_minus_halfdt_;


    std::vector<Field *> Ap_;

    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        Ex_.push_back( ( *this )( ipatch )->EMfields->Ex_ );
        Ey_.push_back( ( *this )( ipatch )->EMfields->Ey_ );
        Ez_.push_back( ( *this )( ipatch )->EMfields->Ez_ );
        Bx_.push_back( ( *this )( ipatch )->EMfields->Bx_ );
        By_.push_back( ( *this )( ipatch )->EMfields->By_ );
        Bz_.push_back( ( *this )( ipatch )->EMfields->Bz_ );
        Bx_m.push_back( ( *this )( ipatch )->EMfields->Bx_m );
        By_m.push_back( ( *this )( ipatch )->EMfields->By_m );
        Bz_m.push_back( ( *this )( ipatch )->EMfields->Bz_m );
        Ex_rel_.push_back( ( *this )( ipatch )->EMfields->Ex_rel_ );
        Ey_rel_.push_back( ( *this )( ipatch )->EMfields->Ey_rel_ );
        Ez_rel_.push_back( ( *this )( ipatch )->EMfields->Ez_rel_ );
        Bx_rel_.push_back( ( *this )( ipatch )->EMfields->Bx_rel_ );
        By_rel_.push_back( ( *this )( ipatch )->EMfields->By_rel_ );
        Bz_rel_.push_back( ( *this )( ipatch )->EMfields->Bz_rel_ );
        Bx_rel_t_plus_halfdt_.push_back( ( *this )( ipatch )->EMfields->Bx_rel_t_plus_halfdt_ );
        By_rel_t_plus_halfdt_.push_back( ( *this )( ipatch )->EMfields->By_rel_t_plus_halfdt_ );
        Bz_rel_t_plus_halfdt_.push_back( ( *this )( ipatch )->EMfields->Bz_rel_t_plus_halfdt_ );
        Bx_rel_t_minus_halfdt_.push_back( ( *this )( ipatch )->EMfields->Bx_rel_t_minus_halfdt_ );
        By_rel_t_minus_halfdt_.push_back( ( *this )( ipatch )->EMfields->By_rel_t_minus_halfdt_ );
        Bz_rel_t_minus_halfdt_.push_back( ( *this )( ipatch )->EMfields->Bz_rel_t_minus_halfdt_ );

        Ap_.push_back( ( *this )( ipatch )->EMfields->Ap_ );
    }

    unsigned int nx_p2_global = ( params.n_space_global[0]+1 );
    //if ( Ex_[0]->dims_.size()>1 ) {
    if( Ex_rel_[0]->dims_.size()>1 ) {
        nx_p2_global *= ( params.n_space_global[1]+1 );
        if( Ex_rel_[0]->dims_.size()>2 ) {
            nx_p2_global *= ( params.n_space_global[2]+1 );
        }
    }


    // compute control parameter
    double norm2_source_term = sqrt( rnew_dot_rnew );
    //double ctrl = rnew_dot_rnew / (double)(nx_p2_global);
    double ctrl = sqrt( rnew_dot_rnew ) / norm2_source_term; // initially is equal to one

    // ---------------------------------------------------------
    // Starting iterative loop for the conjugate gradient method
    // ---------------------------------------------------------
    if( smpi->isMaster() ) {
        DEBUG( "Starting iterative loop for CG method" );
    }
    //cout << std::scientific << ctrl << "\t" << error_max << "\t" << iteration << "\t" << iteration_max << endl;
    while( ( ctrl > error_max ) && ( iteration<iteration_max ) ) {
        iteration++;

        if( ( smpi->isMaster() ) && ( iteration%1000==0 ) ) {
            MESSAGE( "iteration " << iteration << " started with control parameter ctrl = " << 1.0e22*ctrl << " x 1.e-22" );
        }

        // scalar product of the residual
        double r_dot_r = rnew_dot_rnew;

        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->compute_Ap_relativistic_Poisson( ( *this )( ipatch ), gamma_mean );
        }

        // Exchange Ap_ (intra & extra MPI)
        SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Ap_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ap_, *this );


        // scalar product p.Ap
        double p_dot_Ap       = 0.0;
        double p_dot_Ap_local = 0.0;
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            p_dot_Ap_local += ( *this )( ipatch )->EMfields->compute_pAp();
        }
        MPI_Allreduce( &p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );


        // compute new potential and residual
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->update_pand_r( r_dot_r, p_dot_Ap );
        }

        // compute new residual norm
        rnew_dot_rnew       = 0.0;
        rnew_dot_rnew_local = 0.0;
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            rnew_dot_rnew_local += ( *this )( ipatch )->EMfields->compute_r();
        }
        MPI_Allreduce( &rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if( smpi->isMaster() ) {
            DEBUG( "new residual norm: rnew_dot_rnew = " << rnew_dot_rnew );
        }

        // compute new directio
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->update_p( rnew_dot_rnew, r_dot_r );
        }

        // compute control parameter
        
        ctrl = sqrt( rnew_dot_rnew )/norm2_source_term;
        if( smpi->isMaster() ) {
            DEBUG( "iteration " << iteration << " done, exiting with control parameter ctrl = " << 1.0e22*ctrl << " x 1.e-22" );
        }

    }//End of the iterative loop


    // --------------------------------
    // Status of the solver convergence
    // --------------------------------
    if( iteration_max>0 && iteration == iteration_max ) {
        if( smpi->isMaster() )
            WARNING( "Relativistic Poisson solver did not converge: reached maximum iteration number: " << iteration
                     << ", relative err is ctrl = " << 1.0e22*ctrl << "x 1.e-22" );
    } else {
        if( smpi->isMaster() )
            MESSAGE( 1, "Relativistic Poisson solver converged at iteration: " << iteration
                     << ", relative err is ctrl = " << 1.0e22*ctrl << " x 1.e-22" );
    }

    // ------------------------------------------
    // Compute the electromagnetic fields E and B
    // ------------------------------------------

    // sync the potential
    //SyncVectorPatch::exchange( (*this)(ipatch)->EMfields->phi_, *this, smpi );
    //SyncVectorPatch::finalizeexchange( (*this)(ipatch)->EMfields->phi_, *this );

    // compute E and sync
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        // begin loop on patches
        ( *this )( ipatch )->EMfields->initE_relativistic_Poisson( ( *this )( ipatch ), gamma_mean );
    } // end loop on patches
    
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Ex_rel_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ex_rel_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Ey_rel_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ey_rel_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Ez_rel_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ez_rel_, *this );
    //SyncVectorPatch::exchangeE( params, *this, smpi );
    //SyncVectorPatch::finalizeexchangeE( params, *this );

    // Force to zero the average value of electric field, as in traditional Poisson solver
    //// -------------------------------------
    vector<double> E_Add( Ex_rel_[0]->dims_.size(), 0. );
    if( Ex_rel_[0]->dims_.size()==3 ) {
        double Ex_avg_local( 0. ), Ex_avg( 0. ), Ey_avg_local( 0. ), Ey_avg( 0. ), Ez_avg_local( 0. ), Ez_avg( 0. );
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            Ex_avg_local += ( *this )( ipatch )->EMfields->computeExrelSum();
            Ey_avg_local += ( *this )( ipatch )->EMfields->computeEyrelSum();
            Ez_avg_local += ( *this )( ipatch )->EMfields->computeEzrelSum();
        }

        MPI_Allreduce( &Ex_avg_local, &Ex_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &Ey_avg_local, &Ey_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &Ez_avg_local, &Ez_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        E_Add[0] = -Ex_avg/( ( params.n_space[0]+2 )*( params.n_space[1]+1 )*( params.n_space[2]+1 ) );
        E_Add[1] = -Ey_avg/( ( params.n_space[0]+1 )*( params.n_space[1]+2 )*( params.n_space[2]+1 ) );;
        E_Add[2] = -Ez_avg/( ( params.n_space[0]+1 )*( params.n_space[1]+1 )*( params.n_space[2]+2 ) );;
    } else if( Ex_rel_[0]->dims_.size()==2 ) {
        double Ex_XminYmax = 0.0;
        double Ey_XminYmax = 0.0;
        double Ex_XmaxYmin = 0.0;
        double Ey_XmaxYmin = 0.0;

        //The YmaxXmin patch has Patch coordinates X=0, Y=2^m1-1= number_of_patches[1]-1.
        std::vector<int> xcall( 2, 0 );
        xcall[0] = 0;
        xcall[1] = params.number_of_patches[1]-1;
        int patch_YmaxXmin = domain_decomposition_->getDomainId( xcall );
        //The MPI rank owning it is
        int rank_XminYmax = smpi->hrank( patch_YmaxXmin );
        //The YminXmax patch has Patch coordinates X=2^m0-1= number_of_patches[0]-1, Y=0.
        //Its hindex is
        xcall[0] = params.number_of_patches[0]-1;
        xcall[1] = 0;
        int patch_YminXmax = domain_decomposition_->getDomainId( xcall );
        //The MPI rank owning it is
        int rank_XmaxYmin = smpi->hrank( patch_YminXmax );


        //cout << patch_YmaxXmin << " " << rank_XminYmax << " " << patch_YminXmax << " " << rank_XmaxYmin << endl;

        if( smpi->getRank() == rank_XminYmax ) {
            Ex_XminYmax = ( *this )( patch_YmaxXmin-( this->refHindex_ ) )->EMfields->getExrel_XminYmax();
            Ey_XminYmax = ( *this )( patch_YmaxXmin-( this->refHindex_ ) )->EMfields->getEyrel_XminYmax();
        }

        // Xmax-Ymin corner
        if( smpi->getRank() == rank_XmaxYmin ) {
            Ex_XmaxYmin = ( *this )( patch_YminXmax-( this->refHindex_ ) )->EMfields->getExrel_XmaxYmin();
            Ey_XmaxYmin = ( *this )( patch_YminXmax-( this->refHindex_ ) )->EMfields->getEyrel_XmaxYmin();
        }

        MPI_Bcast( &Ex_XminYmax, 1, MPI_DOUBLE, rank_XminYmax, MPI_COMM_WORLD );
        MPI_Bcast( &Ey_XminYmax, 1, MPI_DOUBLE, rank_XminYmax, MPI_COMM_WORLD );

        MPI_Bcast( &Ex_XmaxYmin, 1, MPI_DOUBLE, rank_XmaxYmin, MPI_COMM_WORLD );
        MPI_Bcast( &Ey_XmaxYmin, 1, MPI_DOUBLE, rank_XmaxYmin, MPI_COMM_WORLD );

        //This correction is always done, independantly of the periodicity. Is this correct ?
        E_Add[0] = -0.5*( Ex_XminYmax+Ex_XmaxYmin );
        E_Add[1] = -0.5*( Ey_XminYmax+Ey_XmaxYmin );

#ifdef _3D_LIKE_CENTERING
        double Ex_avg_local( 0. ), Ex_avg( 0. ), Ey_avg_local( 0. ), Ey_avg( 0. );
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            Ex_avg_local += ( *this )( ipatch )->EMfields->computeExrelSum();
            Ey_avg_local += ( *this )( ipatch )->EMfields->computeEyrelSum();
        }

        MPI_Allreduce( &Ex_avg_local, &Ex_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( &Ey_avg_local, &Ey_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        E_Add[0] = -Ex_avg/( ( params.n_space[0]+2 )*( params.n_space[1]+1 ) );
        E_Add[1] = -Ey_avg/( ( params.n_space[0]+1 )*( params.n_space[1]+2 ) );;
#endif

    }

    else if( Ex_rel_[0]->dims_.size()==1 ) {
        double Ex_Xmin = 0.0;
        double Ex_Xmax = 0.0;

        unsigned int rankXmin = 0;
        if( smpi->getRank() == 0 ) {
            //Ex_Xmin = (*Ex1D)(index_bc_min[0]);
            Ex_Xmin = ( *this )( ( 0 )-( this->refHindex_ ) )->EMfields->getExrel_Xmin();
        }
        MPI_Bcast( &Ex_Xmin, 1, MPI_DOUBLE, rankXmin, MPI_COMM_WORLD );

        unsigned int rankXmax = smpi->getSize()-1;
        if( smpi->getRank() == smpi->getSize()-1 ) {
            //Ex_Xmax = (*Ex1D)(index_bc_max[0]);
            Ex_Xmax = ( *this )( ( params.number_of_patches[0]-1 )-( this->refHindex_ ) )->EMfields->getExrel_Xmax();
        }
        MPI_Bcast( &Ex_Xmax, 1, MPI_DOUBLE, rankXmax, MPI_COMM_WORLD );
        E_Add[0] = -0.5*( Ex_Xmin+Ex_Xmax );

#ifdef _3D_LIKE_CENTERING
        double Ex_avg_local( 0. ), Ex_avg( 0. );
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            Ex_avg_local += ( *this )( ipatch )->EMfields->computeExrelSum();
        }

        MPI_Allreduce( &Ex_avg_local, &Ex_avg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        E_Add[0] = -Ex_avg/( ( params.n_space[0]+2 ) );
#endif

    }

    // Centering electrostatic fields
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->centeringErel( E_Add );
    }

    // compute B and sync
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        // begin loop on patches
        ( *this )( ipatch )->EMfields->initB_relativistic_Poisson( ( *this )( ipatch ), gamma_mean );
    } // end loop on patches
    
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bx_rel_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bx_rel_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( By_rel_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( By_rel_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bz_rel_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bz_rel_, *this );


    // Proper spatial centering of the B fields in the Yee Cell through interpolation
    // (from B_rel to B_rel_t_plus_halfdt and B_rel_t_minus_halfdt)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        // begin loop on patches
        ( *this )( ipatch )->EMfields->center_fields_from_relativistic_Poisson( ( *this )( ipatch ) );
    } // end loop on patches

    // Re-exchange the properly spatially centered B field
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bx_rel_t_plus_halfdt_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bx_rel_t_plus_halfdt_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( By_rel_t_plus_halfdt_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( By_rel_t_plus_halfdt_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bz_rel_t_plus_halfdt_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bz_rel_t_plus_halfdt_, *this );
    
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bx_rel_t_minus_halfdt_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bx_rel_t_minus_halfdt_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( By_rel_t_minus_halfdt_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( By_rel_t_minus_halfdt_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bz_rel_t_minus_halfdt_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bz_rel_t_minus_halfdt_, *this );



    MESSAGE( 0, "Summing fields of relativistic species to the grid fields" );
    // sum the fields found  by relativistic Poisson solver to the existing em fields
    // E  = E  + E_rel
    // B  = B  + B_rel_t_plus_halfdt
    // Bm = Bm + B_rel_t_minus_halfdt

    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        // begin loop on patches
        ( *this )( ipatch )->EMfields->sum_rel_fields_to_em_fields( ( *this )( ipatch ) );
    } // end loop on patches

    // Exchange the fields after the addition of the relativistic species fields
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Ex_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ex_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Ey_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ey_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Ez_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ez_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bx_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bx_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( By_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( By_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bz_, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bz_, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bx_m, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bx_m, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( By_m, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( By_m, *this );
    SyncVectorPatch::exchange_along_all_directions_noomp<double,Field>( Bz_m, *this, smpi );
    SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bz_m, *this );

    MESSAGE( 0, "Fields of relativistic species initialized" );
    //!\todo Reduce to find global max
    //if (smpi->isMaster())
    //  MESSAGE(1,"Relativistic Poisson equation solved. Maximum err = ");

    //ptimer.update();
    //MESSAGE("Time in Relativistic Poisson : " << ptimer.getTime() );


    //ptimer.update();
    //MESSAGE("Time in Relativistic Poisson : " << ptimer.getTime() );
    MESSAGE( "Relativistic Poisson finished" );

} // END solveRelativisticPoisson



void VectorPatch::solveRelativisticPoissonAM( Params &params, SmileiMPI *smpi, double time_primal )
{

    //Timer ptimer("global");
    //ptimer.init(smpi);
    //ptimer.restart();
    
    // Assumption: one or more species move in vacuum with mean lorentz gamma factor gamma_mean in the x direction,
    // with low energy spread.
    // The electromagnetic fields of this species can be initialized solving a Poisson-like problem (here informally
    // referred to as "relativistic Poisson problem") and then performing a Lorentz back-transformation to find the
    // electromagnetic fields of the species in the lab frame.
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    // In case of non-monoenergetic relativistic distribution (NOT IMPLEMENTED AT THE MOMENT), the linearity of Maxwell's equations can be exploited:
    // divide the species in quasi-monoenergetic bins with gamma_i and repeat the same procedure for described above
    // for all bins. Finally, in the laboratory frame sum all the fields of the various energy-bin ensembles of particles.
    
    // All the parameters for the Poisson problem (e.g. maximum iteration) are the same used in the namelist
    // for the traditional Poisson problem
    
    // compute gamma_mean for the species for which the field is initialized
    double s_gamma( 0. );
    uint64_t nparticles( 0 );
    for( unsigned int ispec=0 ; ispec<( *this )( 0 )->vecSpecies.size() ; ispec++ ) {
        if( species( 0, ispec )->relativistic_field_initialization ) {
            for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
                if( time_primal==species( ipatch, ispec )->time_relativistic_initialization ) {
                    s_gamma += species( ipatch, ispec )->sum_gamma();
                    nparticles += species( ipatch, ispec )->getNbrOfParticles();
                }
            }
        }
    }
    double gamma_global( 0. );
    MPI_Allreduce( &s_gamma, &gamma_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    uint64_t nparticles_global( 0 );
    MPI_Allreduce( &nparticles, &nparticles_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
    MESSAGE( "GAMMA = " << gamma_global/( double )nparticles_global );
    
    //Timer ptimer("global");
    //ptimer.init(smpi);
    //ptimer.restart();
    
    double gamma_mean = gamma_global/( double )nparticles_global;
    
    unsigned int iteration_max = params.relativistic_poisson_max_iteration;
    double           error_max = params.relativistic_poisson_max_error;
    unsigned int iteration=0;
    
    // Init & Store internal data (phi, r, p, Ap) per patch
    double rnew_dot_rnew_localAM_( 0. );
    double rnew_dot_rnewAM_( 0. );


    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->initPoisson( ( *this )( ipatch ) );
        //cout << std::scientific << "rnew_dot_rnew_local = " << rnew_dot_rnew_local << endl;
        ( *this )( ipatch )->EMfields->initRelativisticPoissonFields( ( *this )( ipatch ) );
    }
    // //cout << std::scientific << "rnew_dot_rnew_local = " << rnew_dot_rnew_local << endl;

    std::vector<Field *> El_;
    std::vector<Field *> Er_;
    std::vector<Field *> Et_;
    std::vector<Field *> Bl_;
    std::vector<Field *> Br_;
    std::vector<Field *> Bt_;
    std::vector<Field *> Bl_m;
    std::vector<Field *> Br_m;
    std::vector<Field *> Bt_m;
    
    std::vector<Field *> El_rel_;
    std::vector<Field *> Er_rel_;
    std::vector<Field *> Et_rel_;
    std::vector<Field *> Bl_rel_;
    std::vector<Field *> Br_rel_;
    std::vector<Field *> Bt_rel_;
    
    std::vector<Field *> Bl_rel_t_plus_halfdt_;
    std::vector<Field *> Br_rel_t_plus_halfdt_;
    std::vector<Field *> Bt_rel_t_plus_halfdt_;
    std::vector<Field *> Bl_rel_t_minus_halfdt_;
    std::vector<Field *> Br_rel_t_minus_halfdt_;
    std::vector<Field *> Bt_rel_t_minus_halfdt_;
    
    std::vector<Field *> Ap_AM_;
    
    // For each mode, repeat the initialization procedure
    // (the relativistic Poisson equation is linear, so it can be decomposed in azimuthal modes)
    for( unsigned int imode=0 ; imode<params.nmodes_rel_field_init ; imode++ ) {
        
        // init Phi, r, p values
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
            emAM->initPoisson_init_phi_r_p_Ap( ( *this )( ipatch ), imode );
            rnew_dot_rnew_localAM_ += emAM->compute_r();
        }
        
        MPI_Allreduce( &rnew_dot_rnew_localAM_, &rnew_dot_rnewAM_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
            El_.push_back( emAM->El_[imode] );
            Er_.push_back( emAM->Er_[imode] );
            Et_.push_back( emAM->Et_[imode] );
            Bl_.push_back( emAM->Bl_[imode] );
            Br_.push_back( emAM->Br_[imode] );
            Bt_.push_back( emAM->Bt_[imode] );
            Bl_m.push_back( emAM->Bl_m[imode] );
            Br_m.push_back( emAM->Br_m[imode] );
            Bt_m.push_back( emAM->Bt_m[imode] );
            El_rel_.push_back( emAM->El_rel_ );
            Er_rel_.push_back( emAM->Er_rel_ );
            Et_rel_.push_back( emAM->Et_rel_ );
            Bl_rel_.push_back( emAM->Bl_rel_ );
            Br_rel_.push_back( emAM->Br_rel_ );
            Bt_rel_.push_back( emAM->Bt_rel_ );
            Bl_rel_t_plus_halfdt_.push_back( emAM->Bl_rel_t_plus_halfdt_ );
            Br_rel_t_plus_halfdt_.push_back( emAM->Br_rel_t_plus_halfdt_ );
            Bt_rel_t_plus_halfdt_.push_back( emAM->Bt_rel_t_plus_halfdt_);
            Bl_rel_t_minus_halfdt_.push_back( emAM->Bl_rel_t_minus_halfdt_ );
            Br_rel_t_minus_halfdt_.push_back( emAM->Br_rel_t_minus_halfdt_ );
            Bt_rel_t_minus_halfdt_.push_back( emAM->Bt_rel_t_minus_halfdt_ );
            
            Ap_AM_.push_back( emAM->Ap_AM_ );
        }

        unsigned int nx_p2_global = ( params.n_space_global[0]+1 );
        //if ( Ex_[0]->dims_.size()>1 ) {
        if( El_rel_[0]->dims_.size()>1 ) {
            nx_p2_global *= ( params.n_space_global[1]+1 );
            if( El_rel_[0]->dims_.size()>2 ) {
                nx_p2_global *= ( params.n_space_global[2]+1 );
            }
        }

        // compute control parameter
        double norm2_source_term = sqrt( std::abs(rnew_dot_rnewAM_) );
        //double ctrl = rnew_dot_rnew / (double)(nx_p2_global);
        double ctrl = sqrt( std::abs(rnew_dot_rnewAM_) ) / norm2_source_term; // initially is equal to one
        
        // ---------------------------------------------------------
        // Starting iterative loop for the conjugate gradient method
        // ---------------------------------------------------------
        if( smpi->isMaster() ) {
            DEBUG( "Starting iterative loop for CG method for the mode "<<imode );
        }
        
        iteration = 0;//MESSAGE("Initial error parameter (must be 1) : "<<ctrl);
        //cout << std::scientific << ctrl << "\t" << error_max << "\t" << iteration << "\t" << iteration_max << endl;
        while( ( ctrl > error_max ) && ( iteration<iteration_max ) ) {
            iteration++;
        
            if( ( smpi->isMaster() ) && ( iteration%1000==0 ) ) {
                MESSAGE( "iteration " << iteration << " started with control parameter ctrl = " << 1.0e22*ctrl << " x 1.e-22" );
            }
        
            // scalar product of the residual
            double r_dot_rAM_ = rnew_dot_rnewAM_;
        
            for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
                emAM->compute_Ap_relativistic_Poisson_AM( ( *this )( ipatch ), gamma_mean, imode );
            }
        
            // Exchange Ap_ (intra & extra MPI)
            SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Ap_AM_, *this, smpi );
            SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Ap_AM_, *this );
        
        
            // scalar product p.Ap
            std::complex<double> p_dot_ApAM_       = 0.0;
            std::complex<double> p_dot_Ap_localAM_ = 0.0;
            for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
                p_dot_Ap_localAM_ += emAM->compute_pAp_AM();
            }
            MPI_Allreduce( &p_dot_Ap_localAM_, &p_dot_ApAM_, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        
        
            // compute new potential and residual
            for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
                emAM->update_pand_r_AM( r_dot_rAM_, p_dot_ApAM_ );
            }
        
            // compute new residual norm
            rnew_dot_rnewAM_       = 0.0;
            rnew_dot_rnew_localAM_ = 0.0;
            for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
                rnew_dot_rnew_localAM_ += emAM->compute_r();
            }
            MPI_Allreduce( &rnew_dot_rnew_localAM_, &rnew_dot_rnewAM_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            if( smpi->isMaster() ) {
                DEBUG( "new residual norm: rnew_dot_rnew = " << rnew_dot_rnewAM_ );
            }
        
            // compute new directio
            for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
                emAM->update_p( rnew_dot_rnewAM_, r_dot_rAM_ );
            }
        
            // compute control parameter
          
            ctrl = sqrt( std::abs(rnew_dot_rnewAM_) )/norm2_source_term;
            if( smpi->isMaster() ) {
                DEBUG( "iteration " << iteration << " done, exiting with control parameter ctrl = " << 1.0e22*ctrl << " x 1.e-22" );
            }
        
        }//End of the iterative loop


        // --------------------------------
        // Status of the solver convergence
        // --------------------------------
        if( iteration_max>0 && iteration == iteration_max ) {
            if( smpi->isMaster() )
                WARNING( "Relativistic Poisson solver did not converge: reached maximum iteration number: " << iteration
                         << ", relative err is ctrl = " << 1.0e22*ctrl << "x 1.e-22" );
        } else {
            if( smpi->isMaster() )
                MESSAGE( 1, "Relativistic Poisson solver converged at iteration: " << iteration
                         << ", relative err is ctrl = " << 1.0e22*ctrl << " x 1.e-22" );
        }

        // ------------------------------------------
        // Compute the electromagnetic fields E and B
        // ------------------------------------------
        
        // sync the potential
        //SyncVectorPatch::exchange( (*this)(ipatch)->EMfields->phi_, *this, smpi );
        //SyncVectorPatch::finalizeexchange( (*this)(ipatch)->EMfields->phi_, *this );
        
        // compute E and sync
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
            // begin loop on patches
            emAM->initE_relativistic_Poisson_AM( ( *this )( ipatch ), gamma_mean, imode );
        } // end loop on patches
        
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( El_rel_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( El_rel_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Er_rel_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Er_rel_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Et_rel_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Et_rel_, *this );

        // compute B and sync
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // begin loop on patches
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
            emAM->initB_relativistic_Poisson_AM( ( *this )( ipatch ), gamma_mean );
        } // end loop on patches
      
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Bl_rel_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bl_rel_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Br_rel_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Br_rel_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Bt_rel_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bt_rel_, *this );
        
        
        // Proper spatial centering of the B fields in the Yee Cell through interpolation
        // (from B_rel to B_rel_t_plus_halfdt and B_rel_t_minus_halfdt)
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // begin loop on patches
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
            emAM->center_fields_from_relativistic_Poisson_AM( ( *this )( ipatch ) );
        } // end loop on patches
        
        // Re-exchange the properly spatially centered B field
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Bl_rel_t_plus_halfdt_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bl_rel_t_plus_halfdt_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Br_rel_t_plus_halfdt_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Br_rel_t_plus_halfdt_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Bt_rel_t_plus_halfdt_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bt_rel_t_plus_halfdt_, *this );
        
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Bl_rel_t_minus_halfdt_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bl_rel_t_minus_halfdt_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Br_rel_t_minus_halfdt_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Br_rel_t_minus_halfdt_, *this );
        SyncVectorPatch::exchange_along_all_directions_noomp<complex<double>,cField>( Bt_rel_t_minus_halfdt_, *this, smpi );
        SyncVectorPatch::finalize_exchange_along_all_directions_noomp( Bt_rel_t_minus_halfdt_, *this );
        
        
        MESSAGE( 0, "Summing fields of relativistic species to the grid fields" );
        // sum the fields found  by relativistic Poisson solver to the existing em fields
        // E  = E  + E_rel
        // B  = B  + B_rel_t_plus_halfdt
        // Bm = Bm + B_rel_t_minus_halfdt
        
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            // begin loop on patches
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
            emAM->sum_rel_fields_to_em_fields_AM( ( *this )( ipatch ), params, imode );
        } // end loop on patches
        
        
        // clean the auxiliary vectors for the present azimuthal mode
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            
            El_.pop_back();
            Er_.pop_back();
            Et_.pop_back();
            Bl_.pop_back();
            Br_.pop_back();
            Bt_.pop_back();
            Bl_m.pop_back();
            Br_m.pop_back();
            Bt_m.pop_back();
            El_rel_.pop_back();
            Er_rel_.pop_back();
            Et_rel_.pop_back();
            Bl_rel_.pop_back();
            Br_rel_.pop_back();
            Bt_rel_.pop_back();
            Bl_rel_t_plus_halfdt_.pop_back();
            Br_rel_t_plus_halfdt_.pop_back();
            Bt_rel_t_plus_halfdt_.pop_back();
            Bl_rel_t_minus_halfdt_.pop_back();
            Br_rel_t_minus_halfdt_.pop_back();
            Bt_rel_t_minus_halfdt_.pop_back();
            
            Ap_AM_.pop_back();

        }
        
    }  // end loop on the modes

    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( ( *this )( ipatch )->EMfields );
        emAM->delete_phi_r_p_Ap( ( *this )( ipatch ) );
        emAM->delete_relativistic_fields( ( *this )( ipatch ) );
    }

    // // Exchange the fields after the addition of the relativistic species fields
    for( unsigned int imode = 0 ; imode < params.nmodes_rel_field_init ; imode++ ) {
        SyncVectorPatch::exchangeE( params, ( *this ), imode, smpi );
        SyncVectorPatch::finalizeexchangeE( params, ( *this ), imode ); // disable async, because of tags which is the same for all modes
    }
    for( unsigned int imode = 0 ; imode < params.nmodes_rel_field_init ; imode++ ) {
        SyncVectorPatch::exchangeB( params, ( *this ), imode, smpi );
        SyncVectorPatch::finalizeexchangeB( params, ( *this ), imode ); // disable async, because of tags which is the same for all modes
    }
    
    MESSAGE( 0, "Fields of relativistic species initialized" );
    //!\todo Reduce to find global max
    //if (smpi->isMaster())
    //  MESSAGE(1,"Relativistic Poisson equation solved. Maximum err = ");
    
    //ptimer.update();
    //MESSAGE("Time in Relativistic Poisson : " << ptimer.getTime() );
    
    
    //ptimer.update();
    //MESSAGE("Time in Relativistic Poisson : " << ptimer.getTime() );
    MESSAGE( "Relativistic Poisson finished" );



}

// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------    BALANCING METHODS    ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------


void VectorPatch::load_balance( Params &params, double time_dual, SmileiMPI *smpi, SimWindow *simWindow, unsigned int itime )
{

    // Compute new patch distribution
    smpi->recompute_patch_count( params, *this, time_dual );

    // Create empty patches according to this new distribution
    this->createPatches( params, smpi, simWindow );

    // Proceed to patch exchange, and delete patch which moved
    this->exchangePatches( smpi, params );

    // Tell that the patches moved this iteration (needed for probes)
    lastIterationPatchesMoved = itime;

}


// ---------------------------------------------------------------------------------------------------------------------
// Explicits patch movement regarding new patch distribution stored in smpi->patch_count
//   - compute send_patch_id_
//   - compute recv_patch_id_
//   - create empty (not really, created like at t0) new patch in recv_patches_
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::createPatches( Params &params, SmileiMPI *smpi, SimWindow *simWindow )
{
    unsigned int n_moved( 0 );
    recv_patches_.resize( 0 );

    // Set Index of the 1st patch of the vector yet on current MPI rank
    // Is this really necessary ? It should be done already ...
    refHindex_ = ( *this )( 0 )->Hindex();

    // Current number of patch
    int nPatches_now = this->size() ;

    // When going to openMP, these two vectors must be stored by patch and not by vectorPatch.
    recv_patch_id_.clear();
    send_patch_id_.clear();

    // istart = Index of the futur 1st patch
    int istart( 0 );
    for( int irk=0 ; irk<smpi->getRank() ; irk++ ) {
        istart += smpi->patch_count[irk];
    }

    // recv_patch_id_ = vector of the hindex this process must own at the end of the exchange.
    for( int ipatch=0 ; ipatch<smpi->patch_count[smpi->getRank()] ; ipatch++ ) {
        recv_patch_id_.push_back( istart+ipatch );
    }


    // Loop on current patches to define patch to send
    for( int ipatch=0 ; ipatch < nPatches_now ; ipatch++ ) {
        //if  current hindex     <  future refHindex   OR      current hindex > future last hindex...
        if( ( refHindex_+ipatch < recv_patch_id_[0] ) || ( refHindex_+ipatch > recv_patch_id_.back() ) ) {
            // Put this patch in the send list.
            send_patch_id_.push_back( ipatch );
        }
    }


    // Backward loop on future patches to define suppress patch in receive list
    // before this loop, recv_patch_id_ stores all patches index define in SmileiMPI::patch_count
    int existing_patch_id = -1;
    for( int ipatch=recv_patch_id_.size()-1 ; ipatch>=0 ; ipatch-- ) {
        //if    future patch hindex  >= current refHindex AND  future patch hindex <= current last hindex
        if( ( recv_patch_id_[ipatch]>=refHindex_ ) && ( recv_patch_id_[ipatch] <= refHindex_ + nPatches_now - 1 ) ) {
            //Store an existing patch id for cloning.
            existing_patch_id = recv_patch_id_[ipatch];
            //Remove this patch from the receive list because I already own it.
            recv_patch_id_.erase( recv_patch_id_.begin()+ipatch );
        }
    }


    // Get an existing patch that will be used for cloning
    if( existing_patch_id<0 ) {
        ERROR( "No patch to clone. This should never happen!" );
    }
    Patch *existing_patch = ( *this )( existing_patch_id-refHindex_ );


    // Create new Patches
    n_moved = simWindow->getNmoved();
    // Store in local vector future patches
    // Loop on the patches I have to receive and do not already own.
    for( unsigned int ipatch=0 ; ipatch < recv_patch_id_.size() ; ipatch++ ) {
        // density profile is initializes as if t = 0 !
        // Species will be cleared when, nbr of particles will be known
        // Creation of a new patch, ready to receive its content from MPI neighbours.
        Patch *newPatch = PatchesFactory::clone( existing_patch, params, smpi, domain_decomposition_, recv_patch_id_[ipatch], n_moved, false );
        newPatch->finalizeMPIenvironment( params );
        //Store pointers to newly created patch in recv_patches_.
        recv_patches_.push_back( newPatch );
    }


} // END createPatches


// ---------------------------------------------------------------------------------------------------------------------
// Exchange patches, based on createPatches initialization
//   take care of reinitialize patch master and diag file managment
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::exchangePatches( SmileiMPI *smpi, Params &params )
{

    //int newMPIrankbis, oldMPIrankbis, tmp;
    int newMPIrank = smpi->getRank() -1;
    int oldMPIrank = smpi->getRank() -1;
    int istart = 0;
    int nmessage = nrequests;

    for( int irk=0 ; irk<smpi->getRank() ; irk++ ) {
        istart += smpi->patch_count[irk];
    }


    // Send particles
    for( unsigned int ipatch=0 ; ipatch < send_patch_id_.size() ; ipatch++ ) {
        // locate rank which will own send_patch_id_[ipatch]
        // We assume patches are only exchanged with neighbours.
        // Once all patches supposed to be sent to the left are done, we send the rest to the right.
        // if hindex of patch to be sent      >  future hindex of the first patch owned by this process
        if( send_patch_id_[ipatch]+refHindex_ > istart ) {
            newMPIrank = smpi->getRank() + 1;
        }

        smpi->isend_species( ( *this )( send_patch_id_[ipatch] ), newMPIrank, ( refHindex_+send_patch_id_[ipatch] )*nmessage, params );
    }

    for( unsigned int ipatch=0 ; ipatch < recv_patch_id_.size() ; ipatch++ ) {
        //if  hindex of patch to be received > first hindex actually owned, that means it comes from the next MPI process and not from the previous anymore.
        if( recv_patch_id_[ipatch] > refHindex_ ) {
            oldMPIrank = smpi->getRank() + 1;
        }

        smpi->recv_species( recv_patches_[ipatch], oldMPIrank, recv_patch_id_[ipatch]*nmessage, params );
    }


    for( unsigned int ipatch=0 ; ipatch < send_patch_id_.size() ; ipatch++ ) {
        smpi->waitall( ( *this )( send_patch_id_[ipatch] ) );
    }

    smpi->barrier();


    // Split the exchangePatches process to avoid deadlock with OpenMPI (observed with OpenMPI on Irene and Poicnare, not with IntelMPI)
    newMPIrank = smpi->getRank() -1;
    oldMPIrank = smpi->getRank() -1;


    // Send fields
    for( unsigned int ipatch=0 ; ipatch < send_patch_id_.size() ; ipatch++ ) {
        // locate rank which will own send_patch_id_[ipatch]
        // We assume patches are only exchanged with neighbours.
        // Once all patches supposed to be sent to the left are done, we send the rest to the right.
        // if hindex of patch to be sent      >  future hindex of the first patch owned by this process
        if( send_patch_id_[ipatch]+refHindex_ > istart ) {
            newMPIrank = smpi->getRank() + 1;
        }

        smpi->isend_fields( ( *this )( send_patch_id_[ipatch] ), newMPIrank, ( refHindex_+send_patch_id_[ipatch] )*nmessage, params );
    }

    for( unsigned int ipatch=0 ; ipatch < recv_patch_id_.size() ; ipatch++ ) {
        //if  hindex of patch to be received > first hindex actually owned, that means it comes from the next MPI process and not from the previous anymore.
        if( recv_patch_id_[ipatch] > refHindex_ ) {
            oldMPIrank = smpi->getRank() + 1;
        }

        smpi->recv_fields( recv_patches_[ipatch], oldMPIrank, recv_patch_id_[ipatch]*nmessage, params );
    }


    for( unsigned int ipatch=0 ; ipatch < send_patch_id_.size() ; ipatch++ ) {
        smpi->waitall( ( *this )( send_patch_id_[ipatch] ) );
    }

    smpi->barrier();


    //Delete sent patches
    int nPatchSend( send_patch_id_.size() );
    for( int ipatch=nPatchSend-1 ; ipatch>=0 ; ipatch-- ) {
        //Ok while at least 1 old patch stay inon current CPU
        delete( *this )( send_patch_id_[ipatch] );
        patches_[ send_patch_id_[ipatch] ] = NULL;
        patches_.erase( patches_.begin() + send_patch_id_[ipatch] );

    }

#ifdef _VECTO
    if( params.vectorization_mode == "adaptive_mixed_sort" ) {
        // adaptive vectorization -- mixed sort
        // Recompute the cell keys before the next step and configure operators
        for( unsigned int ipatch=0 ; ipatch<recv_patch_id_.size() ; ipatch++ ) {
            for( unsigned int ispec=0 ; ispec< recv_patches_[ipatch]->vecSpecies.size() ; ispec++ ) {
                if( dynamic_cast<SpeciesVAdaptiveMixedSort *>( recv_patches_[ipatch]->vecSpecies[ispec] ) ) {
                    dynamic_cast<SpeciesVAdaptiveMixedSort *>( recv_patches_[ipatch]->vecSpecies[ispec] )->compute_part_cell_keys( params );
                    dynamic_cast<SpeciesVAdaptiveMixedSort *>( recv_patches_[ipatch]->vecSpecies[ispec] )->reconfigure_operators( params, recv_patches_[ipatch] );
                }
            }
        }
    } else if( params.vectorization_mode == "adaptive" ) {
        // adaptive vectorization --  always sort
        // Recompute the cell keys before the next step and configure operators
        for( unsigned int ipatch=0 ; ipatch<recv_patch_id_.size() ; ipatch++ ) {
            for( unsigned int ispec=0 ; ispec< recv_patches_[ipatch]->vecSpecies.size() ; ispec++ ) {
                if( dynamic_cast<SpeciesVAdaptive *>( recv_patches_[ipatch]->vecSpecies[ispec] ) ) {
                    dynamic_cast<SpeciesVAdaptive *>( recv_patches_[ipatch]->vecSpecies[ispec] )->compute_part_cell_keys( params );
                    dynamic_cast<SpeciesVAdaptive *>( recv_patches_[ipatch]->vecSpecies[ispec] )->reconfigure_operators( params, recv_patches_[ipatch] );
                }
            }
        }
    }
#endif

    //Put received patches in the global vecPatches
    for( unsigned int ipatch=0 ; ipatch<recv_patch_id_.size() ; ipatch++ ) {
        if( recv_patch_id_[ipatch] > refHindex_ ) {
            patches_.push_back( recv_patches_[ipatch] );
        } else {
            patches_.insert( patches_.begin()+ipatch, recv_patches_[ipatch] );
        }
    }
    recv_patches_.clear();


    for( unsigned int ipatch=0 ; ipatch<patches_.size() ; ipatch++ ) {
        ( *this )( ipatch )->updateMPIenv( smpi );
        if( ( *this )( ipatch )->has_an_MPI_neighbor() ) {
            ( *this )( ipatch )->createType( params );
        } else {
            ( *this )( ipatch )->cleanType();
        }
    }
    this->set_refHindex() ;
    update_field_list( smpi ) ;

} // END exchangePatches

// ---------------------------------------------------------------------------------------------------------------------
// Write in a file patches communications
//   - Send/Recv MPI rank
//   - Send/Recv patch Id
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::output_exchanges( SmileiMPI *smpi )
{
    ofstream output_file;
    ostringstream name( "" );
    name << "debug_output"<<smpi->getRank()<<".txt" ;
    output_file.open( name.str().c_str(), std::ofstream::out | std::ofstream::app );
    int newMPIrank, oldMPIrank;
    newMPIrank = smpi->getRank() -1;
    oldMPIrank = smpi->getRank() -1;
    int istart( 0 );
    for( int irk=0 ; irk<smpi->getRank() ; irk++ ) {
        istart += smpi->patch_count[irk];
    }
    for( unsigned int ipatch=0 ; ipatch < send_patch_id_.size() ; ipatch++ ) {
        if( send_patch_id_[ipatch]+refHindex_ > istart ) {
            newMPIrank = smpi->getRank() + 1;
        }
        output_file << "Rank " << smpi->getRank() << " sending patch " << send_patch_id_[ipatch]+refHindex_ << " to " << newMPIrank << endl;
    }
    for( unsigned int ipatch=0 ; ipatch < recv_patch_id_.size() ; ipatch++ ) {
        if( recv_patch_id_[ipatch] > refHindex_ ) {
            oldMPIrank = smpi->getRank() + 1;
        }
        output_file << "Rank " << smpi->getRank() << " receiving patch " << recv_patch_id_[ipatch] << " from " << oldMPIrank << endl;
    }
    output_file << "NEXT" << endl;
    output_file.close();
} // END output_exchanges

//! Resize vector of field*
void VectorPatch::update_field_list( SmileiMPI *smpi )
{
    int nDim( 0 );
    if( !dynamic_cast<ElectroMagnAM *>( patches_[0]->EMfields ) ) {
        nDim = patches_[0]->EMfields->Ex_->dims_.size();
    } else {
        nDim = static_cast<ElectroMagnAM *>( patches_[0]->EMfields )->El_[0]->dims_.size();
    }
    densities.resize( 3*size() ) ; // Jx + Jy + Jz

    //                          1D  2D  3D
    Bs0.resize( 2*size() ) ; //  2   2   2
    Bs1.resize( 2*size() ) ; //  0   2   2
    Bs2.resize( 2*size() ) ; //  0   0   2

    densitiesLocalx.clear();
    densitiesLocaly.clear();
    densitiesLocalz.clear();
    densitiesMPIx.clear();
    densitiesMPIy.clear();
    densitiesMPIz.clear();
    LocalxIdx.clear();
    LocalyIdx.clear();
    LocalzIdx.clear();
    MPIxIdx.clear();
    MPIyIdx.clear();
    MPIzIdx.clear();

    if( !dynamic_cast<ElectroMagnAM *>( patches_[0]->EMfields ) ) {

        listJx_.resize( size() ) ;
        listJy_.resize( size() ) ;
        listJz_.resize( size() ) ;
        listrho_.resize( size() ) ;
        listEx_.resize( size() ) ;
        listEy_.resize( size() ) ;
        listEz_.resize( size() ) ;
        listBx_.resize( size() ) ;
        listBy_.resize( size() ) ;
        listBz_.resize( size() ) ;

        if( patches_[0]->EMfields->envelope != NULL ) {
            listA_.resize( size() ) ;
            listA0_.resize( size() ) ;
            listPhi_.resize( size() ) ;
            listPhi0_.resize( size() ) ;
            listGradPhix_.resize( size() ) ;
            listGradPhiy_.resize( size() ) ;
            listGradPhiz_.resize( size() ) ;
            listGradPhix0_.resize( size() ) ;
            listGradPhiy0_.resize( size() ) ;
            listGradPhiz0_.resize( size() ) ;
            listEnv_Chi_.resize( size() ) ;
        }

        for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
            listJx_[ipatch] = patches_[ipatch]->EMfields->Jx_ ;
            listJy_[ipatch] = patches_[ipatch]->EMfields->Jy_ ;
            listJz_[ipatch] = patches_[ipatch]->EMfields->Jz_ ;
            listrho_[ipatch] =patches_[ipatch]->EMfields->rho_;
            listEx_[ipatch] = patches_[ipatch]->EMfields->Ex_ ;
            listEy_[ipatch] = patches_[ipatch]->EMfields->Ey_ ;
            listEz_[ipatch] = patches_[ipatch]->EMfields->Ez_ ;
            listBx_[ipatch] = patches_[ipatch]->EMfields->Bx_ ;
            listBy_[ipatch] = patches_[ipatch]->EMfields->By_ ;
            listBz_[ipatch] = patches_[ipatch]->EMfields->Bz_ ;
        }
        if( patches_[0]->EMfields->envelope != NULL ) {
            for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
                listA_[ipatch]         = patches_[ipatch]->EMfields->envelope->A_ ;
                listA0_[ipatch]        = patches_[ipatch]->EMfields->envelope->A0_ ;
                listPhi_[ipatch]       = patches_[ipatch]->EMfields->envelope->Phi_ ;
                listPhi0_[ipatch]      = patches_[ipatch]->EMfields->envelope->Phi_m ;
                listGradPhix_[ipatch]  = patches_[ipatch]->EMfields->envelope->GradPhix_ ;
                listGradPhiy_[ipatch]  = patches_[ipatch]->EMfields->envelope->GradPhiy_ ;
                listGradPhiz_[ipatch]  = patches_[ipatch]->EMfields->envelope->GradPhiz_ ;
                listGradPhix0_[ipatch] = patches_[ipatch]->EMfields->envelope->GradPhix_m ;
                listGradPhiy0_[ipatch] = patches_[ipatch]->EMfields->envelope->GradPhiy_m ;
                listGradPhiz0_[ipatch] = patches_[ipatch]->EMfields->envelope->GradPhiz_m ;
                listEnv_Chi_[ipatch]   = patches_[ipatch]->EMfields->Env_Chi_ ;
            }
        }

    } else {
        unsigned int nmodes = static_cast<ElectroMagnAM *>( patches_[0]->EMfields )->El_.size();
        listJl_.resize( nmodes ) ;
        listJr_.resize( nmodes ) ;
        listJt_.resize( nmodes ) ;
        listrho_AM_.resize( nmodes ) ;
        listJls_.resize( nmodes ) ;
        listJrs_.resize( nmodes ) ;
        listJts_.resize( nmodes ) ;
        listrhos_AM_.resize( nmodes ) ;
        listEl_.resize( nmodes ) ;
        listEr_.resize( nmodes ) ;
        listEt_.resize( nmodes ) ;
        listBl_.resize( nmodes ) ;
        listBr_.resize( nmodes ) ;
        listBt_.resize( nmodes ) ;

        for( unsigned int imode=0 ; imode < nmodes ; imode++ ) {
            listJl_[imode].resize( size() );
            listJr_[imode].resize( size() );
            listJt_[imode].resize( size() );
            listrho_AM_[imode].resize( size() );
            listEl_[imode].resize( size() );
            listEr_[imode].resize( size() );
            listEt_[imode].resize( size() );
            listBl_[imode].resize( size() );
            listBr_[imode].resize( size() );
            listBt_[imode].resize( size() );
            for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
                listJl_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Jl_[imode] ;
                listJr_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Jr_[imode] ;
                listJt_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Jt_[imode] ;
                listrho_AM_[imode][ipatch] =static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->rho_AM_[imode];
                listEl_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->El_[imode] ;
                listEr_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Er_[imode] ;
                listEt_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Et_[imode] ;
                listBl_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Bl_[imode] ;
                listBr_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Br_[imode] ;
                listBt_[imode][ipatch]     = static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields )->Bt_[imode] ;
            }
        }
    }

    B_localx.clear();
    B_MPIx.clear();

    B1_localy.clear();
    B1_MPIy.clear();

    B2_localz.clear();
    B2_MPIz.clear();

    for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
        densities[ipatch         ] = patches_[ipatch]->EMfields->Jx_ ;
        densities[ipatch+  size()] = patches_[ipatch]->EMfields->Jy_ ;
        densities[ipatch+2*size()] = patches_[ipatch]->EMfields->Jz_ ;

        Bs0[ipatch       ] = patches_[ipatch]->EMfields->By_ ;
        Bs0[ipatch+size()] = patches_[ipatch]->EMfields->Bz_ ;

        // TO DO , B size depend of nDim
        // Pas grave, au pire inutil
        Bs1[ipatch       ] = patches_[ipatch]->EMfields->Bx_ ;
        Bs1[ipatch+size()] = patches_[ipatch]->EMfields->Bz_ ;

        // TO DO , B size depend of nDim
        // Pas grave, au pire inutil
        Bs2[ipatch       ] = patches_[ipatch]->EMfields->Bx_ ;
        Bs2[ipatch+size()] = patches_[ipatch]->EMfields->By_ ;
    }

    for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
        if( ( *this )( ipatch )->has_an_MPI_neighbor( 0 ) ) {
            MPIxIdx.push_back( ipatch );
        }
        if( ( *this )( ipatch )->has_an_local_neighbor( 0 ) ) {
            LocalxIdx.push_back( ipatch );
        }
    }
    if( nDim>1 ) {
        for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
            if( ( *this )( ipatch )->has_an_MPI_neighbor( 1 ) ) {
                MPIyIdx.push_back( ipatch );
            }
            if( ( *this )( ipatch )->has_an_local_neighbor( 1 ) ) {
                LocalyIdx.push_back( ipatch );
            }
        }
        if( nDim>2 ) {
            for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {

                if( ( *this )( ipatch )->has_an_MPI_neighbor( 2 ) ) {
                    MPIzIdx.push_back( ipatch );
                }
                if( ( *this )( ipatch )->has_an_local_neighbor( 2 ) ) {
                    LocalzIdx.push_back( ipatch );
                }
            }
        }
    }

    B_MPIx.resize( 2*MPIxIdx.size() );
    B_localx.resize( 2*LocalxIdx.size() );
    B1_MPIy.resize( 2*MPIyIdx.size() );
    B1_localy.resize( 2*LocalyIdx.size() );
    B2_MPIz.resize( 2*MPIzIdx.size() );
    B2_localz.resize( 2*LocalzIdx.size() );

    densitiesMPIx.resize( 3*MPIxIdx.size() );
    densitiesLocalx.resize( 3*LocalxIdx.size() );
    densitiesMPIy.resize( 3*MPIyIdx.size() );
    densitiesLocaly.resize( 3*LocalyIdx.size() );
    densitiesMPIz.resize( 3*MPIzIdx.size() );
    densitiesLocalz.resize( 3*LocalzIdx.size() );

    int mpix( 0 ), locx( 0 ), mpiy( 0 ), locy( 0 ), mpiz( 0 ), locz( 0 );

    for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {

        if( ( *this )( ipatch )->has_an_MPI_neighbor( 0 ) ) {
            B_MPIx[mpix               ] = patches_[ipatch]->EMfields->By_;
            B_MPIx[mpix+MPIxIdx.size()] = patches_[ipatch]->EMfields->Bz_;

            densitiesMPIx[mpix                 ] = patches_[ipatch]->EMfields->Jx_;
            densitiesMPIx[mpix+  MPIxIdx.size()] = patches_[ipatch]->EMfields->Jy_;
            densitiesMPIx[mpix+2*MPIxIdx.size()] = patches_[ipatch]->EMfields->Jz_;
            mpix++;
        }
        if( ( *this )( ipatch )->has_an_local_neighbor( 0 ) ) {
            B_localx[locx                 ] = patches_[ipatch]->EMfields->By_;
            B_localx[locx+LocalxIdx.size()] = patches_[ipatch]->EMfields->Bz_;

            densitiesLocalx[locx                   ] = patches_[ipatch]->EMfields->Jx_;
            densitiesLocalx[locx+  LocalxIdx.size()] = patches_[ipatch]->EMfields->Jy_;
            densitiesLocalx[locx+2*LocalxIdx.size()] = patches_[ipatch]->EMfields->Jz_;
            locx++;
        }
    }
    if( nDim>1 ) {
        for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
            if( ( *this )( ipatch )->has_an_MPI_neighbor( 1 ) ) {
                B1_MPIy[mpiy               ] = patches_[ipatch]->EMfields->Bx_;
                B1_MPIy[mpiy+MPIyIdx.size()] = patches_[ipatch]->EMfields->Bz_;

                densitiesMPIy[mpiy                 ] = patches_[ipatch]->EMfields->Jx_;
                densitiesMPIy[mpiy+  MPIyIdx.size()] = patches_[ipatch]->EMfields->Jy_;
                densitiesMPIy[mpiy+2*MPIyIdx.size()] = patches_[ipatch]->EMfields->Jz_;
                mpiy++;
            }
            if( ( *this )( ipatch )->has_an_local_neighbor( 1 ) ) {
                B1_localy[locy                 ] = patches_[ipatch]->EMfields->Bx_;
                B1_localy[locy+LocalyIdx.size()] = patches_[ipatch]->EMfields->Bz_;

                densitiesLocaly[locy                   ] = patches_[ipatch]->EMfields->Jx_;
                densitiesLocaly[locy+  LocalyIdx.size()] = patches_[ipatch]->EMfields->Jy_;
                densitiesLocaly[locy+2*LocalyIdx.size()] = patches_[ipatch]->EMfields->Jz_;
                locy++;
            }
        }
        if( nDim>2 ) {
            for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
                if( ( *this )( ipatch )->has_an_MPI_neighbor( 2 ) ) {
                    B2_MPIz[mpiz               ] = patches_[ipatch]->EMfields->Bx_;
                    B2_MPIz[mpiz+MPIzIdx.size()] = patches_[ipatch]->EMfields->By_;

                    densitiesMPIz[mpiz                 ] = patches_[ipatch]->EMfields->Jx_;
                    densitiesMPIz[mpiz+  MPIzIdx.size()] = patches_[ipatch]->EMfields->Jy_;
                    densitiesMPIz[mpiz+2*MPIzIdx.size()] = patches_[ipatch]->EMfields->Jz_;
                    mpiz++;
                }
                if( ( *this )( ipatch )->has_an_local_neighbor( 2 ) ) {
                    B2_localz[locz                 ] = patches_[ipatch]->EMfields->Bx_;
                    B2_localz[locz+LocalzIdx.size()] = patches_[ipatch]->EMfields->By_;

                    densitiesLocalz[locz                   ] = patches_[ipatch]->EMfields->Jx_;
                    densitiesLocalz[locz+  LocalzIdx.size()] = patches_[ipatch]->EMfields->Jy_;
                    densitiesLocalz[locz+2*LocalzIdx.size()] = patches_[ipatch]->EMfields->Jz_;
                    locz++;
                }
            }
        }

    }

    if( !dynamic_cast<ElectroMagnAM *>( patches_[0]->EMfields ) ) {
        for( unsigned int ipatch = 0 ; ipatch < size() ; ipatch++ ) {
            listJx_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 1 );
            listJy_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 2 );
            listJz_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 3 );
            listBx_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 6 );
            listBy_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 7 );
            listBz_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 8 );
            listrho_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 4 );
        }
        if( patches_[0]->EMfields->envelope != NULL ) {
            for( unsigned int ipatch = 0 ; ipatch < size() ; ipatch++ ) {
                listA_ [ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listA0_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listPhi_ [ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listPhi0_ [ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listGradPhix_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listGradPhiy_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listGradPhiz_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listGradPhix0_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listGradPhiy0_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listGradPhiz0_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
                listEnv_Chi_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 ) ;
            }

        }
    } else {
        unsigned int nmodes = static_cast<ElectroMagnAM *>( patches_[0]->EMfields )->El_.size();
        for( unsigned int imode=0 ; imode < nmodes ; imode++ ) {
            for( unsigned int ipatch = 0 ; ipatch < size() ; ipatch++ ) {
                listJl_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                listJr_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                listJt_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                listBl_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                listBr_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                listBt_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                listrho_AM_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
            }
        }
    }
}



void VectorPatch::update_field_list( int ispec, SmileiMPI *smpi )
{
    #pragma omp barrier
    if( !dynamic_cast<ElectroMagnAM *>( patches_[0]->EMfields ) ) {
        #pragma omp single
        {
            if( patches_[0]->EMfields->Jx_s [ispec] )
            {
                listJxs_.resize( size() ) ;
            } else
            {
                listJxs_.clear();
            }
            if( patches_[0]->EMfields->Jy_s [ispec] )
            {
                listJys_.resize( size() ) ;
            } else
            {
                listJys_.clear();
            }
            if( patches_[0]->EMfields->Jz_s [ispec] )
            {
                listJzs_.resize( size() ) ;
            } else
            {
                listJzs_.clear();
            }
            if( patches_[0]->EMfields->rho_s[ispec] )
            {
                listrhos_.resize( size() ) ;
            } else
            {
                listrhos_.clear();
            }

            if( patches_[0]->EMfields->envelope != NULL )
            {
                if( patches_[0]->EMfields->Env_Chi_s[ispec] ) {
                    listEnv_Chis_.resize( size() ) ;
                } else {
                    listEnv_Chis_.clear();
                }
            }
        }

        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
            if( patches_[ipatch]->EMfields->Jx_s [ispec] ) {
                listJxs_ [ipatch] = patches_[ipatch]->EMfields->Jx_s [ispec];
                listJxs_ [ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
            }
            if( patches_[ipatch]->EMfields->Jy_s [ispec] ) {
                listJys_ [ipatch] = patches_[ipatch]->EMfields->Jy_s [ispec];
                listJys_ [ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
            }
            if( patches_[ipatch]->EMfields->Jz_s [ispec] ) {
                listJzs_ [ipatch] = patches_[ipatch]->EMfields->Jz_s [ispec];
                listJzs_ [ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
            }
            if( patches_[ipatch]->EMfields->rho_s[ispec] ) {
                listrhos_[ipatch] = patches_[ipatch]->EMfields->rho_s[ispec];
                listrhos_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
            }

            if( patches_[0]->EMfields->envelope != NULL ) {
                if( patches_[ipatch]->EMfields->Env_Chi_s[ispec] ) {
                    listEnv_Chis_[ipatch] = patches_[ipatch]->EMfields->Env_Chi_s[ispec];
                    listEnv_Chis_[ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                }
            }
        }
    } else { // if ( dynamic_cast<ElectroMagnAM*>(patches_[0]->EMfields) )
        ElectroMagnAM *emAM =  static_cast<ElectroMagnAM *>( patches_[0]->EMfields );
        unsigned int nmodes = emAM->El_.size();
        unsigned int n_species = emAM->n_species;
        #pragma omp single
        {
            for( unsigned int imode=0 ; imode < nmodes ; imode++ ) {
                unsigned int ifield = imode*n_species + ispec ;
                if( emAM->Jl_s [ifield] ) {
                    listJls_[imode].resize( size() ) ;
                } else {
                    listJls_[imode].clear();
                }
                if( emAM->Jr_s [ifield] ) {
                    listJrs_[imode].resize( size() ) ;
                } else {
                    listJrs_[imode].clear();
                }
                if( emAM->Jt_s [ifield] ) {
                    listJts_[imode].resize( size() ) ;
                } else {
                    listJts_[imode].clear();
                }
                if( emAM->rho_AM_s [ifield] ) {
                    listrhos_AM_[imode].resize( size() ) ;
                } else {
                    listrhos_AM_[imode].clear();
                }
            }
        }
        for( unsigned int imode=0 ; imode < nmodes ; imode++ ) {
            unsigned int ifield = imode*n_species + ispec ;
            #pragma omp for schedule(static)
            for( unsigned int ipatch=0 ; ipatch < size() ; ipatch++ ) {
                emAM =  static_cast<ElectroMagnAM *>( patches_[ipatch]->EMfields );
                if( emAM->Jl_s [ifield] ) {
                    listJls_[imode][ipatch] = emAM->Jl_s [ifield];
                    listJls_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                }
                if( emAM->Jr_s [ifield] ) {
                    listJrs_[imode][ipatch] = emAM->Jr_s [ifield];
                    listJrs_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                }
                if( emAM->Jt_s [ifield] ) {
                    listJts_[imode][ipatch] = emAM->Jt_s [ifield];
                    listJts_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                }
                if( emAM->rho_AM_s [ifield] ) {
                    listrhos_AM_[imode][ipatch] = emAM->rho_AM_s [ifield];
                    listrhos_AM_[imode][ipatch]->MPIbuff.defineTags( patches_[ipatch], smpi, 0 );
                }
            }
        }

    }



}


void VectorPatch::applyAntennas( double time )
{
#ifdef  __DEBUG
    if( nAntennas>0 ) {
        #pragma omp single
        TITLE( "Applying antennas at time t = " << time );
    }
#endif

    // Loop antennas
    for( unsigned int iAntenna=0; iAntenna<nAntennas; iAntenna++ ) {

        // Get intensity from antenna of the first patch
        #pragma omp single
        antenna_intensity = patches_[0]->EMfields->antennas[iAntenna].time_profile->valueAt( time );

        // Loop patches to apply
        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
            patches_[ipatch]->EMfields->applyAntenna( iAntenna, antenna_intensity );
        }

    }
}

// For each patch, apply the collisions
void VectorPatch::applyCollisions( Params &params, int itime, Timers &timers )
{
    timers.collisions.restart();

    if( Collisions::debye_length_required )
        #pragma omp for schedule(runtime)
        for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
            Collisions::calculate_debye_length( params, patches_[ipatch] );
        }

    unsigned int ncoll = patches_[0]->vecCollisions.size();

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ )
        for( unsigned int icoll=0 ; icoll<ncoll; icoll++ ) {
            patches_[ipatch]->vecCollisions[icoll]->collide( params, patches_[ipatch], itime, localDiags );
        }

    #pragma omp single
    for( unsigned int icoll=0 ; icoll<ncoll; icoll++ ) {
        Collisions::debug( params, itime, icoll, *this );
    }
    #pragma omp barrier

    timers.collisions.update();
}

// For all patches, allocate a field if not allocated
void VectorPatch::allocateField( unsigned int ifield, Params &params )
{
    for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
        if( params.geometry != "AMcylindrical" ) {
            Field *field = emfields( ipatch )->allFields[ifield];
            if( field->data_ != NULL ) {
                continue;
            }
            if( field->name.substr( 0, 2 )=="Jx" ) {
                field->allocateDims( 0, false );
            } else if( field->name.substr( 0, 2 )=="Jy" ) {
                field->allocateDims( 1, false );
            } else if( field->name.substr( 0, 2 )=="Jz" ) {
                field->allocateDims( 2, false );
            } else if( field->name.substr( 0, 2 )=="Rh" ) {
                field->allocateDims();
            }
            //MESSAGE("HNA4");
        } else {
            cField2D *field = static_cast<cField2D *>( emfields( ipatch )->allFields[ifield] );
            if( field->cdata_ != NULL ) {
                continue;
            }
            if( field->name.substr( 0, 2 )=="Jl" ) {
                field->allocateDims( 0, false );
            } else if( field->name.substr( 0, 2 )=="Jr" ) {
                field->allocateDims( 1, false );
            } else if( field->name.substr( 0, 2 )=="Jt" ) {
                field->allocateDims( 2, false );
            } else if( field->name.substr( 0, 2 )=="Rh" ) {
                field->allocateDims();
            }
        }
    }
}


// For each patch, apply external fields
void VectorPatch::applyExternalFields()
{
    for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
        patches_[ipatch]->EMfields->applyExternalFields( ( *this )( ipatch ) );    // Must be patch
    }
}


// For each patch, apply external fields
void VectorPatch::saveExternalFields( Params &params )
{
    if( params.save_magnectic_fields_for_SM ) {
        for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
            patches_[ipatch]->EMfields->saveExternalFields( ( *this )( ipatch ) );    // Must be patch
        }
    }
}


// Print information on the memory consumption
void VectorPatch::check_memory_consumption( SmileiMPI *smpi )
{
    long int particlesMem( 0 );
    for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ )
        for( unsigned int ispec=0 ; ispec<patches_[ipatch]->vecSpecies.size(); ispec++ ) {
            particlesMem += patches_[ipatch]->vecSpecies[ispec]->getMemFootPrint();
        }
    MESSAGE( 1, "(Master) Species part = " << ( int )( ( double )particlesMem / 1024./1024. ) << " MB" );

    long double dParticlesMem = ( double )particlesMem / 1024./1024./1024.;
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dParticlesMem, &dParticlesMem, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( 1, setprecision( 3 ) << "Global Species part = " << dParticlesMem << " GB" );

    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&particlesMem, &particlesMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( 1, "Max Species part = " << ( int )( ( double )particlesMem / 1024./1024. ) << " MB" );

    // fieldsMem contains field per species and average fields
    long int fieldsMem( 0 );
    for( unsigned int ipatch=0 ; ipatch<size() ; ipatch++ ) {
        fieldsMem += patches_[ipatch]->EMfields->getMemFootPrint();
    }
    MESSAGE( 1, "(Master) Fields part = " << ( int )( ( double )fieldsMem / 1024./1024. ) << " MB" );

    long double dFieldsMem = ( double )fieldsMem / 1024./1024./1024.;
    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dFieldsMem, &dFieldsMem, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MESSAGE( 1, setprecision( 3 ) << "Global Fields part = " << dFieldsMem << " GB" );

    MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&fieldsMem, &fieldsMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
    MESSAGE( 1, "Max Fields part = " << ( int )( ( double )fieldsMem / 1024./1024. ) << " MB" );


    for( unsigned int idiags=0 ; idiags<globalDiags.size() ; idiags++ ) {
        // fieldsMem contains field per species
        long int diagsMem( 0 );
        diagsMem += globalDiags[idiags]->getMemFootPrint();

        long double dDiagsMem = ( double )diagsMem / 1024./1024./1024.;
        MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dDiagsMem, &dDiagsMem, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        if( dDiagsMem>0. ) {
            MESSAGE( 1, "(Master) " <<  globalDiags[idiags]->filename << "  = " << ( int )( ( double )diagsMem / 1024./1024. ) << " MB" );
            MESSAGE( 1, setprecision( 3 ) << "Global " <<  globalDiags[idiags]->filename << " = " << dDiagsMem << " GB" );
        }

        MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&diagsMem, &diagsMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
        if( dDiagsMem>0. ) {
            MESSAGE( 1, "Max " <<  globalDiags[idiags]->filename << " = " << ( int )( ( double )diagsMem / 1024./1024. ) << " MB" );
        }
    }

    for( unsigned int idiags=0 ; idiags<localDiags.size() ; idiags++ ) {
        // fieldsMem contains field per species
        long int diagsMem( 0 );
        diagsMem += localDiags[idiags]->getMemFootPrint();

        long double dDiagsMem = ( double )diagsMem / 1024./1024./1024.;
        MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&dDiagsMem, &dDiagsMem, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        if( dDiagsMem>0. ) {
            MESSAGE( 1, "(Master) " <<  localDiags[idiags]->filename << "  = " << ( int )( ( double )diagsMem / 1024./1024. ) << " MB" );
            MESSAGE( 1, setprecision( 3 ) << "Global " <<  localDiags[idiags]->filename << " = " << dDiagsMem << " GB" );
        }

        MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&diagsMem, &diagsMem, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
        if( dDiagsMem>0. ) {
            MESSAGE( 1, "Max " <<  localDiags[idiags]->filename << " = " << ( int )( ( double )diagsMem / 1024./1024. ) << " MB" );
        }
    }

    // Read value in /proc/pid/status
    //Tools::printMemFootPrint( "End Initialization" );
}


void VectorPatch::save_old_rho( Params &params )
{
    int n=0;
    #pragma omp for schedule(static)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        n = ( *this )( ipatch )->EMfields->rhoold_->dims_[0]*( *this )( ipatch )->EMfields->rhoold_->dims_[1]; //*(*this)(ipatch)->EMfields->rhoold_->dims_[2];
        if( params.nDim_field ==3 ) {
            n*=( *this )( ipatch )->EMfields->rhoold_->dims_[2];
        }
        std::memcpy( ( *this )( ipatch )->EMfields->rhoold_->data_, ( *this )( ipatch )->EMfields->rho_->data_, sizeof( double )*n );
    }
}



// Print information on the memory consumption
void VectorPatch::check_expected_disk_usage( SmileiMPI *smpi, Params &params, Checkpoint &checkpoint )
{
    if( smpi->isMaster() ) {

        MESSAGE( 1, "WARNING: disk usage by non-uniform particles maybe strongly underestimated," );
        MESSAGE( 1, "   especially when particles are created at runtime (ionization, pair generation, etc.)" );
        MESSAGE( 1, "" );

        // Find the initial and final timesteps for this simulation
        int istart = 0, istop = params.n_time;
        // If restarting simulation define the starting point
        if( params.restart ) {
            istart = checkpoint.this_run_start_step+1;
        }
        // If leaving the simulation after dump, define the stopping point
        if( checkpoint.dump_step > 0 && checkpoint.exit_after_dump ) {
            int ncheckpoint = ( istart/( int )checkpoint.dump_step ) + 1;
            int nextdumptime = ncheckpoint * ( int )checkpoint.dump_step;
            if( nextdumptime < istop ) {
                istop = nextdumptime;
            }
        }

        MESSAGE( 1, "Expected disk usage for diagnostics:" );
        // Calculate the footprint from local then global diagnostics
        uint64_t diagnostics_footprint = 0;
        for( unsigned int idiags=0 ; idiags<localDiags.size() ; idiags++ ) {
            uint64_t footprint = localDiags[idiags]->getDiskFootPrint( istart, istop, patches_[0] );
            diagnostics_footprint += footprint;
            MESSAGE( 2, "File " << localDiags[idiags]->filename << ": " << Tools::printBytes( footprint ) );
        }
        for( unsigned int idiags=0 ; idiags<globalDiags.size() ; idiags++ ) {
            uint64_t footprint = globalDiags[idiags]->getDiskFootPrint( istart, istop, patches_[0] );
            diagnostics_footprint += footprint;
            MESSAGE( 2, "File " << globalDiags[idiags]->filename << ": " << Tools::printBytes( footprint ) );
        }
        MESSAGE( 1, "Total disk usage for diagnostics: " << Tools::printBytes( diagnostics_footprint ) );
        MESSAGE( 1, "" );

        // If checkpoints to be written, estimate their size
        if( checkpoint.dump_step > 0 || checkpoint.dump_minutes > 0 ) {
            MESSAGE( 1, "Expected disk usage for each checkpoint:" );

            // - Contribution from the grid
            ElectroMagn *EM = patches_[0]->EMfields;
            //     * Calculate first the number of grid points in total
            uint64_t n_grid_points = 1;
            for( unsigned int i=0; i<params.nDim_field; i++ ) {
                n_grid_points *= ( params.n_space[i] + 2*params.oversize[i]+1 );
            }
            n_grid_points *= params.tot_number_of_patches;
            //     * Now calculate the total number of fields
            unsigned int n_fields = 9
                                    + EM->Exfilter.size() + EM->Eyfilter.size() + EM->Ezfilter.size()
                                    + EM->Bxfilter.size() + EM->Byfilter.size() + EM->Bzfilter.size();
            for( unsigned int idiag=0; idiag<EM->allFields_avg.size(); idiag++ ) {
                n_fields += EM->allFields_avg[idiag].size();
            }
            //     * Conclude the total field disk footprint
            uint64_t checkpoint_fields_footprint = n_grid_points * ( uint64_t )( n_fields * sizeof( double ) );
            MESSAGE( 2, "For fields: " << Tools::printBytes( checkpoint_fields_footprint ) );

            // - Contribution from particles
            uint64_t checkpoint_particles_footprint = 0;
            for( unsigned int ispec=0 ; ispec<patches_[0]->vecSpecies.size() ; ispec++ ) {
                Species *s = patches_[0]->vecSpecies[ispec];
                Particles *p = s->particles;
                //     * Calculate the size of particles' individual parameters
                uint64_t one_particle_size = 0;
                one_particle_size += ( p->Position.size() + p->Momentum.size() + 1 ) * sizeof( double );
                one_particle_size += 1 * sizeof( short );
                if( p->tracked ) {
                    one_particle_size += 1 * sizeof( uint64_t );
                }
                //     * Calculate an approximate number of particles
                PeekAtSpecies peek( params, ispec );
                uint64_t number_of_particles = peek.totalNumberofParticles();
                //     * Calculate the size of the first_index and last_index arrays
                uint64_t b_size = ( s->first_index.size() + s->last_index.size() ) * params.tot_number_of_patches * sizeof( int );
                //     * Conclude the disk footprint of this species
                checkpoint_particles_footprint += one_particle_size*number_of_particles + b_size;
            }
            MESSAGE( 2, "For particles: " << Tools::printBytes( checkpoint_particles_footprint ) );

            // - Contribution from diagnostics
            uint64_t checkpoint_diags_footprint = 0;
            //     * Averaged field diagnostics
            n_fields = 0;
            for( unsigned int idiag=0; idiag<EM->allFields_avg.size(); idiag++ ) {
                n_fields += EM->allFields_avg[idiag].size();
            }
            checkpoint_diags_footprint += n_grid_points * ( uint64_t )( n_fields * sizeof( double ) );
            //     * Screen diagnostics
            for( unsigned int idiag=0; idiag<globalDiags.size(); idiag++ )
                if( DiagnosticScreen *screen = dynamic_cast<DiagnosticScreen *>( globalDiags[idiag] ) ) {
                    checkpoint_diags_footprint += screen->data_sum.size() * sizeof( double );
                }
            MESSAGE( 2, "For diagnostics: " << Tools::printBytes( checkpoint_diags_footprint ) );

            uint64_t checkpoint_footprint = checkpoint_fields_footprint + checkpoint_particles_footprint + checkpoint_diags_footprint;
            MESSAGE( 1, "Total disk usage for one checkpoint: " << Tools::printBytes( checkpoint_footprint ) );
        }

    }
}

void VectorPatch::runEnvelopeModule( Params &params,
        SmileiMPI *smpi,
        SimWindow *simWindow,
        double time_dual, Timers &timers, int itime )
{
    // interpolate envelope for susceptibility deposition, project susceptibility for envelope equation, momentum advance
    ponderomotive_update_susceptibility_and_momentum( params, smpi, simWindow, time_dual, timers, itime );

    // comm and sum susceptibility
    sumSusceptibility( params, time_dual, timers, itime, simWindow, smpi );

    // solve envelope equation and comm envelope
    solveEnvelope( params, simWindow, itime, time_dual, timers, smpi );

    // interp updated envelope for position advance, update positions and currents for Maxwell's equations
    ponderomotive_update_position_and_currents( params, smpi, simWindow, time_dual, timers, itime );

}

// ---------------------------------------------------------------------------------------------------------------------
// For all patch, update momentum for particles interacting with envelope
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::ponderomotive_update_susceptibility_and_momentum( Params &params,
        SmileiMPI *smpi,
        SimWindow *simWindow,
        double time_dual, Timers &timers, int itime )
{

    #pragma omp single
    diag_flag = needsRhoJsNow( itime );

    timers.particles.restart();

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        ( *this )( ipatch )->EMfields->restartEnvChi();
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            if( ( *this )( ipatch )->vecSpecies[ispec]->isProj( time_dual, simWindow ) || diag_flag ) {
                if( species( ipatch, ispec )->ponderomotive_dynamics ) {
                    if( ( *this )( ipatch )->vecSpecies[ispec]->vectorized_operators || params.cell_sorting )
                        species( ipatch, ispec )->ponderomotive_update_susceptibility_and_momentum( time_dual, ispec,
                                emfields( ipatch ),
                                params, diag_flag,
                                ( *this )( ipatch ), smpi,
                                localDiags );
                    else {
                        if( params.vectorization_mode == "adaptive" ) {
                            species( ipatch, ispec )->scalar_ponderomotive_update_susceptibility_and_momentum( time_dual, ispec,
                                    emfields( ipatch ),
                                    params, diag_flag,
                                    ( *this )( ipatch ), smpi,
                                    localDiags );
                        } else {
                            species( ipatch, ispec )->Species::ponderomotive_update_susceptibility_and_momentum( time_dual, ispec,
                                    emfields( ipatch ),
                                    params, diag_flag,
                                    ( *this )( ipatch ), smpi,
                                    localDiags );
                        }
                    }

                } // end condition on ponderomotive dynamics
            } // end diagnostic or projection if condition on species
        } // end loop on species
    } // end loop on patches

    timers.particles.update( );
#ifdef __DETAILED_TIMERS
    timers.interp_fields_env.update( *this, params.printNow( itime ) );
    timers.proj_susceptibility.update( *this, params.printNow( itime ) );
    timers.push_mom.update( *this, params.printNow( itime ) );
#endif

} // END ponderomotive_update_susceptibility_and_momentum

void VectorPatch::ponderomotive_update_position_and_currents( Params &params,
        SmileiMPI *smpi,
        SimWindow *simWindow,
        double time_dual, Timers &timers, int itime )
{

    #pragma omp single
    diag_flag = needsRhoJsNow( itime );

    timers.particles.restart();

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
        for( unsigned int ispec=0 ; ispec<( *this )( ipatch )->vecSpecies.size() ; ispec++ ) {
            if( ( *this )( ipatch )->vecSpecies[ispec]->isProj( time_dual, simWindow ) || diag_flag ) {
                if( species( ipatch, ispec )->ponderomotive_dynamics ) {
                    if( ( *this )( ipatch )->vecSpecies[ispec]->vectorized_operators || params.cell_sorting )
                        species( ipatch, ispec )->ponderomotive_update_position_and_currents( time_dual, ispec,
                                emfields( ipatch ),
                                params, diag_flag, partwalls( ipatch ),
                                ( *this )( ipatch ), smpi,
                                localDiags );
                    else {

                        if( params.vectorization_mode == "adaptive" ) {
                            species( ipatch, ispec )->scalar_ponderomotive_update_position_and_currents( time_dual, ispec,
                                    emfields( ipatch ),
                                    params, diag_flag, partwalls( ipatch ),
                                    ( *this )( ipatch ), smpi,
                                    localDiags );
                        } else {
                            species( ipatch, ispec )->Species::ponderomotive_update_position_and_currents( time_dual, ispec,
                                    emfields( ipatch ),
                                    params, diag_flag, partwalls( ipatch ),
                                    ( *this )( ipatch ), smpi,
                                    localDiags );
                        }
                    }

                } // end condition on ponderomotive dynamics
            } // end diagnostic or projection if condition on species
        } // end loop on species
    } // end loop on patches

    timers.particles.update( params.printNow( itime ) );
#ifdef __DETAILED_TIMERS
    timers.interp_env_old.update( *this, params.printNow( itime ) );
    timers.proj_currents.update( *this, params.printNow( itime ) );
    timers.push_pos.update( *this, params.printNow( itime ) );
    timers.cell_keys.update( *this, params.printNow( itime ) );
#endif

    timers.syncPart.restart();
    for( unsigned int ispec=0 ; ispec<( *this )( 0 )->vecSpecies.size(); ispec++ ) {
        if( ( *this )( 0 )->vecSpecies[ispec]->ponderomotive_dynamics ) {
            if( ( *this )( 0 )->vecSpecies[ispec]->isProj( time_dual, simWindow ) ) {
                SyncVectorPatch::exchangeParticles( ( *this ), ispec, params, smpi, timers, itime ); // Included sort_part
            } // end condition on species
        } // end condition on envelope dynamics
    } // end loop on species
    timers.syncPart.update( params.printNow( itime ) );



} // END ponderomotive_update_position_and_currents


void VectorPatch::init_new_envelope( Params &params )
{
    if( ( *this )( 0 )->EMfields->envelope!=NULL ) {
        // for all patches, init new envelope from input namelist parameters
        for( unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++ ) {
            ( *this )( ipatch )->EMfields->envelope->initEnvelope( ( *this )( ipatch ), ( *this )( ipatch )->EMfields );
        } // end loop on patches
    }
} // END init_new_envelope
