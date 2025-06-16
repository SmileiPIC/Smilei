/*
 * Checkpoint.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "Checkpoint.h"

#include <sstream>
#include <iomanip>
#include <string>

#include <mpi.h>

#include "Params.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "Region.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC2D_SM.h"
#include "ElectroMagnBC3D_SM.h"
#include "ElectroMagnBC2D_PML.h"
#include "ElectroMagnBC3D_PML.h"
#include "ElectroMagnBCAM_PML.h"
#include "EnvelopeBCAM_PML.h"
#include "EnvelopeBC2D_PML.h"
#include "EnvelopeBC3D_PML.h"
#include "Laser.h"
#include "Species.h"
#include "DiagnosticProbes.h"
#include "DiagnosticScreen.h"
#include "DiagnosticTrack.h"
#include "LaserEnvelope.h"
#include "BinaryProcesses.h"
#include "CollisionalNuclearReaction.h"

using namespace std;

// static varable must be defined and initialized here
int Checkpoint::signal_received=0;

Checkpoint::Checkpoint( Params &params, SmileiMPI *smpi ) :
    dump_number( 0 ),
    this_run_start_step( 0 ),
    exit_asap( false ),
    dump_step( 0 ),
    dump_minutes( 0.0 ),
    exit_after_dump( true ),
    time_reference( MPI_Wtime() ),
    keep_n_dumps( 2 ),
    keep_n_dumps_max( 10000 ),
    dump_deflate( 0 ),
    file_grouping( 0 )
{

    if( PyTools::nComponents( "Checkpoints" ) > 0 ) {

        PyTools::extract( "dump_step", dump_step, "Checkpoints"  );
        if( dump_step > 0 ) {
            MESSAGE( 1, "Code will dump after " << dump_step << " steps" );
        }

        PyTools::extract( "dump_minutes", dump_minutes, "Checkpoints"  );
        if( dump_minutes > 0 ) {
            MESSAGE( 1, "Code will dump after " << dump_minutes << " minutes" );
        }

        PyTools::extract( "keep_n_dumps", keep_n_dumps, "Checkpoints"  );
        if( keep_n_dumps<1 ) {
            keep_n_dumps=1;
        }

        if( keep_n_dumps > keep_n_dumps_max ) {
            WARNING( "Smilei supports a maximum of keep_n_dumps of "<< keep_n_dumps_max );
            keep_n_dumps = keep_n_dumps_max;
        }

        PyTools::extract( "exit_after_dump", exit_after_dump, "Checkpoints"  );

        PyTools::extract( "dump_deflate", dump_deflate, "Checkpoints"  );

        PyTools::extract( "file_grouping", file_grouping, "Checkpoints"  );
        if( file_grouping > 0 ) {
            if( file_grouping > ( unsigned int )( smpi->getSize() ) ) {
                file_grouping = smpi->getSize();
            }
            MESSAGE( 1, "Code will group checkpoint files by "<< file_grouping );
        }

        smpi->barrier();

        if( params.restart ) {
            std::vector<std::string> restart_files;
            if( ! PyTools::extractV( "restart_files", restart_files, "Checkpoints" ) ) {
                ERROR( "Internal parameter `restart_files` not understood. This should not happen" );
            }

            // This will open all dumps and pick the last one
            restart_file = "";
            for( unsigned int num_dump=0; num_dump<restart_files.size(); num_dump++ ) {
                string dump_name = restart_files[num_dump];
                H5Read f( dump_name, NULL, false );
                if( f.valid() ) {
                    unsigned int dump_step = 0;
                    f.attr( "dump_step", dump_step );
                    if( dump_step > this_run_start_step ) {
                        this_run_start_step = dump_step;
                        restart_file = dump_name;
                        dump_number = num_dump;
                        f.attr( "dump_number", dump_number );
                    }
                }
            }

            if( restart_file == "" ) {
                ERROR( "Cannot find a valid restart file for rank "<<smpi->getRank() );
            }

            // Make sure all ranks have the same dump number
            // Different numbers can be due to corrupted restart files
            if( ! smpi->test_mode ) {
                unsigned int prev_number;
                MPI_Status status;
                MPI_Sendrecv(
                    &dump_number, 1, MPI_UNSIGNED, (smpi->getRank()+1) % smpi->getSize(), smpi->getRank(),
                    &prev_number, 1, MPI_UNSIGNED, (smpi->getRank()+smpi->getSize()-1) % smpi->getSize(), (smpi->getRank()+smpi->getSize()-1)%smpi->getSize(),
                    smpi->world(), &status
                );
                int problem = (prev_number != dump_number);
                int any_problem;
                MPI_Allreduce( &problem, &any_problem, 1, MPI_INT, MPI_LOR, smpi->world() );
                if( any_problem ) {
                    if( problem ) {
                        ostringstream t("");
                        t << "\t[ERROR]: Issue with restart file on rank " << smpi->getRank() << endl;
                        cout << t.str();
                    }
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
            }

#ifdef  __DEBUG
            MESSAGEALL( 2, " : Restarting fields and particles, dump_number = " << dump_number << " step=" << this_run_start_step << "\n\t\t" << restart_file );
#else
            MESSAGE( 2, "Restarting fields and particles at step: " << this_run_start_step );
            MESSAGE( 2, "                            master file: " << restart_file );
#endif

        }
    }

    if( dump_step>0 || dump_minutes>0. ) {
        if( exit_after_dump ) {
            MESSAGE( 1, "Code will exit after dump" );
        } else {
            ostringstream message( "" );
            message << "Code will dump";
            if( dump_step>0 ) {
                message << " every "<< dump_step << " steps,";
            }
            if( dump_minutes>0. ) {
                message << " every "<<dump_minutes<< " min,";
            }
            message << " keeping "<< keep_n_dumps << " dumps at maximum";
            MESSAGE( 1, message.str() );
        }
    }

    // registering signal handler
    if( SIG_ERR == signal( SIGUSR1, Checkpoint::signal_callback_handler ) ) {
        WARNING( "Cannot catch signal SIGUSR1" );
    }
    if( SIG_ERR == signal( SIGUSR2, Checkpoint::signal_callback_handler ) ) {
        WARNING( "Cannot catch signal SIGUSR2" );
    }

    nDim_particle=params.nDim_particle;
}

Checkpoint::~Checkpoint() {}

void Checkpoint::dump( VectorPatch &vecPatches, Region &region, unsigned int itime, SmileiMPI *smpi, SimWindow *simWindow, Params &params )
{
    bool dump_now = false;
    
    // Find out whether we should make a checkpoint due to dump_minutes
    if( dump_minutes != 0. ) {
        // master checks whenever we passed the time limit
        if( smpi->isMaster() &&  MPI_Wtime() - time_reference > dump_minutes * 60. ) {
            dump_now = true;
        }
        // Broadcast the result
        MPI_Bcast( &dump_now, 1, MPI_CXX_BOOL, 0, smpi->world() );
    }
    
    // Dump if at requested timestep
    dump_now = dump_now || ( dump_step != 0 && ( ( itime-this_run_start_step ) % dump_step == 0 ) );
    
    if( signal_received != 0 || dump_now ) {
        dumpAll( vecPatches, region, itime,  smpi, simWindow, params );
        if( exit_after_dump || ( ( signal_received!=0 ) && ( signal_received != SIGUSR2 ) ) ) {
            exit_asap = true;
        }
        signal_received = 0;
        time_reference = MPI_Wtime();
    }
}

void Checkpoint::dumpAll( VectorPatch &vecPatches, Region &region, unsigned int itime,  SmileiMPI *smpi, SimWindow *simWin,  Params &params )
{
    unsigned int num_dump=dump_number % keep_n_dumps;

    ostringstream nameDumpTmp( "" );
    nameDumpTmp << "checkpoints" << PATH_SEPARATOR;
    if( file_grouping>0 ) {
        nameDumpTmp << setfill( '0' ) << setw( int( 1+log10( smpi->getSize()/file_grouping+1 ) ) ) << smpi->getRank()/file_grouping << PATH_SEPARATOR;
    }

    nameDumpTmp << "dump-" << setfill( '0' ) << setw( 5 ) << num_dump << "-" << setfill( '0' ) << setw( 10 ) << smpi->getRank() << ".h5" ;
    std::string dumpName=nameDumpTmp.str();


    H5Write f( dumpName );
    dump_number++;

#ifdef  __DEBUG
    //MESSAGEALL( "Step " << itime << " : DUMP fields and particles " << dumpName );
    MESSAGEALL( " Checkpoint #" << dumpName << " at iteration " << itime << " dumped" );
#else
    MESSAGE( " Checkpoint #" << num_dump << " at iteration " << itime << " dumped" );
#endif

#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_ACCELERATOR_GPU_OACC )
    MESSAGE( " Copying device data in main memory" );
    // TODO(Etienne M): This may very well be redundant if we did a diagnostic
    // during the last iteration. Indeed, we copy everything from the device to
    // the main memory every time we do even the smallest diagnostic.
    // if( !vecPatches.diag_flag ) could be used (?) to avoid copying everything
    // a second time.
    vecPatches.copyParticlesFromDeviceToHost();
    vecPatches.copyFieldsFromDeviceToHost();
    //vecPatches.copyDeviceStateToHost();
#endif

    // Write basic attributes
    f.attr( "Version", string( __VERSION ) );

    f.attr( "dump_step", itime );
    f.attr( "dump_number", dump_number );

    f.vect( "patch_count", smpi->patch_count );

    // Write diags scalar data
    DiagnosticScalar *scalars = static_cast<DiagnosticScalar *>( vecPatches.globalDiags[0] );
    f.attr( "latest_timestep",   scalars->latest_timestep );
    // Scalars only by master
    if( smpi->isMaster() ) {
        f.attr( "Energy_time_zero",  scalars->Energy_time_zero );
        f.attr( "EnergyUsedForNorm", scalars->EnergyUsedForNorm );
    }
    // Poynting scalars
    unsigned int k=0;
    for( unsigned int j=0; j<2; j++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
        for( unsigned int i=0; i<params.nDim_field; i++ ) { //axis 0=x, 1=y, 2=z
            if( scalars->necessary_poy[k] ) {
                string poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
                double poy_val = 0.;
                for( unsigned ipatch=0; ipatch<vecPatches.size(); ipatch++ ) {
                    poy_val += vecPatches( ipatch )->EMfields->poynting[j][i];
                }
                f.attr( poy_name, poy_val );
                k++;
            }
        }
    }

    // Write the diags screen data
    if( smpi->isMaster() ) {
        unsigned int iscreen = 0;
        for( unsigned int idiag=0; idiag<vecPatches.globalDiags.size(); idiag++ ) {
            if( DiagnosticScreen *screen = dynamic_cast<DiagnosticScreen *>( vecPatches.globalDiags[idiag] ) ) {
                ostringstream diagName( "" );
                diagName << "DiagScreen" << iscreen;
                f.vect( diagName.str(), *(screen->getData()) );
                iscreen++;
            }
        }
    }

    // Write all the patch data
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++ ) {

        // Open a group
        ostringstream patch_name( "" );
        patch_name << setfill( '0' ) << setw( 6 ) << vecPatches( ipatch )->Hindex();
        string patchName=Tools::merge( "patch-", patch_name.str() );
        H5Write g = f.group( patchName.c_str() );

        dumpPatch( vecPatches( ipatch ), params, g );

        // Random number generator state
        g.attr( "xorshift32_state", vecPatches( ipatch )->rand_->xorshift32_state );

    }

    if (params.multiple_decomposition) {
        // Open a group
        ostringstream patch_name( "" );
        patch_name << setfill( '0' ) << setw( 6 ) << region.patch_->Hindex();
        string patchName=Tools::merge( "region-", patch_name.str() );
        H5Write g = f.group( patchName.c_str() );
        dumpPatch( region.patch_, params, g );
    }

    // Write the latest Id that the MPI processes have given to each species
    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ ) {
        if( DiagnosticTrack *track = dynamic_cast<DiagnosticTrack *>( vecPatches.localDiags[idiag] ) ) {
            ostringstream n( "" );
            n<< "latest_ID_" << track->species_name_;
            f.attr( n.str(), track->latest_Id, H5T_NATIVE_UINT64 );
        }
    }

    // Write the moving window status
    if( simWin!=NULL ) {
        dumpMovingWindow( f, simWin );
    }

}


void Checkpoint::dumpPatch( Patch *patch, Params &params, H5Write &g )
{
    ElectroMagn * EMfields = patch->EMfields;
    if (  params.geometry != "AMcylindrical" ) {
        dumpFieldsPerProc( g, EMfields->Ex_ );
        dumpFieldsPerProc( g, EMfields->Ey_ );
        dumpFieldsPerProc( g, EMfields->Ez_ );
        dumpFieldsPerProc( g, EMfields->Bx_ );
        dumpFieldsPerProc( g, EMfields->By_ );
        dumpFieldsPerProc( g, EMfields->Bz_ );
        dumpFieldsPerProc( g, EMfields->Bx_m );
        dumpFieldsPerProc( g, EMfields->By_m );
        dumpFieldsPerProc( g, EMfields->Bz_m );
        if (params.use_BTIS3){
            dumpFieldsPerProc( g, EMfields->By_mBTIS3 );
            dumpFieldsPerProc( g, EMfields->Bz_mBTIS3 );  
        }
        for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
            if( dynamic_cast<ElectroMagnBC2D_PML *>( EMfields->emBoundCond[bcId] )){
                ElectroMagnBC2D_PML *embc = static_cast<ElectroMagnBC2D_PML *>( EMfields->emBoundCond[bcId] );
                if (embc->Hx_) dump_PML(embc, g);
            } else if( dynamic_cast<ElectroMagnBC3D_PML *>( EMfields->emBoundCond[bcId] )){
                ElectroMagnBC3D_PML *embc = static_cast<ElectroMagnBC3D_PML *>( EMfields->emBoundCond[bcId] );
                if (embc->Hx_) dump_PML(embc, g);
            }
        }
    }
    else {
        for ( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            dump_cFieldsPerProc( g, emAM->El_[imode] );
            dump_cFieldsPerProc( g, emAM->Er_[imode] );
            dump_cFieldsPerProc( g, emAM->Et_[imode] );
            dump_cFieldsPerProc( g, emAM->Bl_[imode] );
            dump_cFieldsPerProc( g, emAM->Br_[imode] );
            dump_cFieldsPerProc( g, emAM->Bt_[imode] );
            dump_cFieldsPerProc( g, emAM->Bl_m[imode] );
            dump_cFieldsPerProc( g, emAM->Br_m[imode] );
            dump_cFieldsPerProc( g, emAM->Bt_m[imode] );
            if (params.use_BTIS3){
                dump_cFieldsPerProc( g, emAM->Br_mBTIS3[imode] );
                dump_cFieldsPerProc( g, emAM->Bt_mBTIS3[imode] );
            }
            if( params.is_pxr ) {
                dump_cFieldsPerProc( g, emAM->rho_old_AM_[imode] );
            }
            for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
                if( dynamic_cast<ElectroMagnBCAM_PML *>( EMfields->emBoundCond[bcId] )){
                    ElectroMagnBCAM_PML *embc = static_cast<ElectroMagnBCAM_PML *>( EMfields->emBoundCond[bcId] );
                    if (embc->Hl_[imode]) dump_PML(embc, g, imode);
                }
            }
        }
    }

    if( EMfields->envelope!=NULL ) {
        dump_cFieldsPerProc( g, EMfields->envelope->A_ );
        dump_cFieldsPerProc( g, EMfields->envelope->A0_ );
        dumpFieldsPerProc( g, EMfields->Env_Chi_ );
        if (  params.geometry != "AMcylindrical" ) {
            for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
                if( dynamic_cast<EnvelopeBC2D_PML *>( EMfields->envelope->EnvBoundCond[bcId] )){
                    EnvelopeBC2D_PML *envbc = static_cast<EnvelopeBC2D_PML *>( EMfields->envelope->EnvBoundCond[bcId] );
                    if (envbc->A_n_) dump_PMLenvelope(envbc, g, bcId);
                } else if( dynamic_cast<EnvelopeBC3D_PML *>( EMfields->envelope->EnvBoundCond[bcId] )){
                    EnvelopeBC3D_PML *envbc = static_cast<EnvelopeBC3D_PML *>( EMfields->envelope->EnvBoundCond[bcId] );
                    if (envbc->A_n_) dump_PMLenvelope(envbc, g, bcId);
                }
            }

        } else {
            for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
                if( dynamic_cast<EnvelopeBCAM_PML *>( EMfields->envelope->EnvBoundCond[bcId] )){
                    EnvelopeBCAM_PML *envbc = static_cast<EnvelopeBCAM_PML *>( EMfields->envelope->EnvBoundCond[bcId] );
                    if (envbc->A_n_) dump_PMLenvelopeAM(envbc, g, bcId);
                }
            }
        }
    }

    // filtered Electric fields
    if( EMfields->filter_ ) {
        if (params.geometry!="AMcylindrical"){
            for( unsigned int i=0; i<EMfields->filter_->Ex_.size(); i++ ) {
                dumpFieldsPerProc( g, EMfields->filter_->Ex_[i] );
            }
            for( unsigned int i=0; i<EMfields->filter_->Ey_.size(); i++ ) {
                dumpFieldsPerProc( g, EMfields->filter_->Ey_[i] );
            }
            for( unsigned int i=0; i<EMfields->filter_->Ez_.size(); i++ ) {
                dumpFieldsPerProc( g, EMfields->filter_->Ez_[i] );
            }
        } else{
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            for ( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {
                for( unsigned int i=0; i<EMfields->filter_->El_[imode].size(); i++ ) {
                    dumpFieldsPerProc( g, emAM->filter_->El_[imode][i] );
                }
                for( unsigned int i=0; i<EMfields->filter_->Er_[imode].size(); i++ ) {
                    dumpFieldsPerProc( g, emAM->filter_->Er_[imode][i] );
                }
                for( unsigned int i=0; i<EMfields->filter_->Et_[imode].size(); i++ ) {
                    dumpFieldsPerProc( g, emAM->filter_->Et_[imode][i] );
                }
            } // end loop on modes        
        } // end if condition on geometry
      
    }
    
    // Fields required for DiagFields
    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
        ostringstream group_name( "" );
        group_name << "FieldsForDiag" << idiag;
        H5Write diag = g.group( group_name.str() );

        for( unsigned int ifield=0; ifield<EMfields->allFields_avg[idiag].size(); ifield++ ) {
            dumpFieldsPerProc( diag, EMfields->allFields_avg[idiag][ifield] );
        }
    }

    // Fields required for DiagProbes with time integral
    for( unsigned int iprobe=0; iprobe<patch->probes.size(); iprobe++ ) {
        unsigned int nFields = patch->probes[iprobe]->integrated_data.size();
        if( nFields > 0 ) {
            ostringstream group_name( "" );
            group_name << "DataForProbes" << iprobe;
            H5Write diag = g.group( group_name.str() );
            for( unsigned int ifield=0; ifield<nFields; ifield++ ) {
                ostringstream field( "" );
                field << "field" << ifield;
                diag.vect( field.str(), patch->probes[iprobe]->integrated_data[ifield] );
            }
        }
    }

    if( ( EMfields->extFields.size()>0 ) && ( params.save_magnectic_fields_for_SM ) ) {
        for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
            if( ! EMfields->emBoundCond[bcId] ) {
                continue;
            }
            if( dynamic_cast<ElectroMagnBC1D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC1D_SM *embc = static_cast<ElectroMagnBC1D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                H5Write b = g.group( groupName );
                b.attr( "By_val", embc->By_val_ );
                b.attr( "Bz_val", embc->Bz_val_ );
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                H5Write b = g.group( groupName );
                b.vect( "Bx_val", embc->B_val[0] );
                b.vect( "By_val", embc->B_val[1] );
                b.vect( "Bz_val", embc->B_val[2] );
            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                H5Write b = g.group( groupName );
                if( embc->B_val[0] ) {
                    dumpFieldsPerProc( b, embc->B_val[0] );
                }
                if( embc->B_val[1] ) {
                    dumpFieldsPerProc( b, embc->B_val[1] );
                }
                if( embc->B_val[2] ) {
                    dumpFieldsPerProc( b, embc->B_val[2] );
                }
            }
        }
    }

    g.flush();
    g.attr( "species", patch->vecSpecies.size() );

    for( unsigned int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++ ) {
        Species *spec = patch->vecSpecies[ispec];

        ostringstream name( "" );
        name << setfill( '0' ) << setw( 2 ) << ispec;
        string groupName=Tools::merge( "species-", name.str(), "-", spec->name_ );
        H5Write s = g.group( groupName );

        s.attr( "partCapacity", spec->getParticlesCapacity() );
        s.attr( "partSize", spec->getNbrOfParticles() );

        s.attr( "nrj_bc_lost", spec->nrj_bc_lost );
        s.attr( "nrj_mw_inj", spec->nrj_mw_inj );
        s.attr( "nrj_mw_out", spec->nrj_mw_out );
        s.attr( "nrj_new_part", spec->nrj_new_part_ );
        s.attr( "radiatedEnergy", spec->nrj_radiated_ );

        if( spec->getNbrOfParticles()>0 ) {
            dumpParticles( s, *spec->particles );
            s.vect( "first_index", spec->particles->first_index );
            s.vect( "last_index", spec->particles->last_index );
        }
        
        // Copy birth records that haven't been written yet
        if( spec->birth_records_ ) {
            H5Write b = s.group( "birth_records" );
            b.vect( "birth_time", spec->birth_records_->birth_time_ );
            dumpParticles( b, spec->birth_records_->p_ );
        }
    } // End for ispec

    // Save some scalars
    g.attr( "nrj_mw_inj", EMfields->nrj_mw_inj );
    g.attr( "nrj_mw_out", EMfields->nrj_mw_out );
    // Manage some collisions parameters
    std::vector<double> rate_multiplier( patch->vecBPs.size() );
    for( unsigned int icoll = 0; icoll < patch->vecBPs.size(); icoll++ ) {
        if( CollisionalNuclearReaction * NR = patch->vecBPs[icoll]->nuclear_reactions_ ) {
            rate_multiplier[icoll] =  NR->rate_multiplier_;
        }
    }
    g.vect( "nuclear_reaction_multiplier", rate_multiplier );

    // Save data for LaserProfileFile (i.e. LaserOffset)
    for( unsigned int ii = 0; ii < 2; ii++ ) {
        if( ! EMfields->emBoundCond[ii] ) continue;
        std::vector<Laser *> * veclaser = & EMfields->emBoundCond[ii]->vecLaser;
        for( unsigned int ilas = 0; ilas < veclaser->size(); ilas++ ) {
            Laser * las = (*veclaser)[ilas];
            for( unsigned int iprof = 0; iprof < las->profiles.size(); iprof++ ) {
                LaserProfile * prof = las->profiles[iprof];
                if( dynamic_cast<LaserProfileFile *>( prof ) ) {
                    LaserProfileFile *p = static_cast<LaserProfileFile *>( prof );
                    if( p->magnitude && p->phase ) {
                        ostringstream t1, t2, t3, t4;
                        t1 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_mag";
                        t2 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_phase";
                        t3 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_omega";
                        t4 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_dims";
                        p->magnitude->name = t1.str();
                        p->phase->name = t2.str();
                        dumpFieldsPerProc( g, p->magnitude );
                        dumpFieldsPerProc( g, p->phase );
                        g.vect( t3.str(), p->omega );
                        g.vect( t4.str(), p->magnitude->dims_ );
                    }
                }
            }
        }
    }
};


void Checkpoint::readPatchDistribution( SmileiMPI *smpi, SimWindow *simWin )
{
    H5Read f( restart_file );

    // Read basic attributes
    string dump_version;
    f.attr( "Version", dump_version );

    if( dump_version != string( __VERSION ) ) {
        WARNING( "The code version that dumped the file is " << dump_version );
        WARNING( "                while running version is " << string( __VERSION ) );
    }

    vector<int> patch_count( smpi->getSize() );
    f.vect( "patch_count", patch_count );
    smpi->patch_count = patch_count;

    smpi->patch_refHindexes.resize( smpi->patch_count.size(), 0 );
    smpi->patch_refHindexes[0] = 0;
    for( int rk=1 ; rk<smpi->smilei_sz ; rk++ ) {
        smpi->patch_refHindexes[rk] = smpi->patch_refHindexes[rk-1] + smpi->patch_count[rk-1];
    }

    // load window status : required to know the patch movement
    restartMovingWindow( f, simWin );
}


void Checkpoint::restartAll( VectorPatch &vecPatches, Region &region, SmileiMPI *smpi, Params &params )
{
    MESSAGE( 1, "READING fields and particles for restart" );

    H5Read f( restart_file );

    // Write diags scalar data
    DiagnosticScalar *scalars = static_cast<DiagnosticScalar *>( vecPatches.globalDiags[0] );
    f.attr( "latest_timestep", scalars->latest_timestep );
    // Scalars only by master
    if( smpi->isMaster() ) {
        f.attr( "Energy_time_zero",  scalars->Energy_time_zero );
        f.attr( "EnergyUsedForNorm", scalars->EnergyUsedForNorm );
    }
    // Poynting scalars
    unsigned int k=0;
    for( unsigned int j=0; j<2; j++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
        for( unsigned int i=0; i<params.nDim_field; i++ ) { //axis 0=x, 1=y, 2=z
            string poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
            if( f.hasAttr( poy_name ) ) {
                f.attr( poy_name, vecPatches( 0 )->EMfields->poynting[j][i] );
            }
            k++;
        }
    }

    // Read the diags screen data
    if( smpi->isMaster() ) {
        unsigned int iscreen = 0;
        for( unsigned int idiag=0; idiag<vecPatches.globalDiags.size(); idiag++ ) {
            if( DiagnosticScreen *screen = dynamic_cast<DiagnosticScreen *>( vecPatches.globalDiags[idiag] ) ) {
                ostringstream diagName( "" );
                diagName << "DiagScreen" << iscreen;
                int target_size = screen->getData()->size();
                int vect_size = f.vectSize( diagName.str() );
                if( vect_size == target_size ) {
                    f.vect( diagName.str(), *(screen->getData()) );
                } else {
                    WARNING( "Restart: DiagScreen[" << iscreen << "] size mismatch. Previous data discarded" );
                }
                iscreen++;
            }
        }
    }

    // Read all the patch data
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++ ) {

        ostringstream patch_name( "" );
        patch_name << setfill( '0' ) << setw( 6 ) << vecPatches( ipatch )->Hindex();
        string patchName = Tools::merge( "patch-", patch_name.str() );
        H5Read g = f.group( patchName );

        restartPatch( vecPatches( ipatch ), params, g );

        // Random number generator state
        g.attr( "xorshift32_state", vecPatches( ipatch )->rand_->xorshift32_state );

    }

    if (params.multiple_decomposition) {
        ostringstream patch_name( "" );
        patch_name << setfill( '0' ) << setw( 6 ) << region.patch_->Hindex();
        string patchName = Tools::merge( "region-", patch_name.str() );
        H5Read g = f.group( patchName );
        restartPatch( region.patch_, params, g );
    }

    // Read the latest Id that the MPI processes have given to each species
    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ ) {
        if( DiagnosticTrack *track = dynamic_cast<DiagnosticTrack *>( vecPatches.localDiags[idiag] ) ) {
            ostringstream n( "" );
            n<< "latest_ID_" << track->species_name_;
            if( f.hasAttr( n.str() ) ) {
                f.attr( n.str(), track->latest_Id, H5T_NATIVE_UINT64 );
            } else {
                track->IDs_done=false;
            }
        }
    }

}


void Checkpoint::readRegionDistribution( Region &region )
{
    int read_hindex( -1 );

    hid_t file = H5Fopen(restart_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t grp = H5Gopen(file,"/",H5P_DEFAULT);

    hsize_t nobj;
    H5Gget_num_objs(grp, &nobj);
    char memb_name[1024];
    for (int i = 0; i < (int)(nobj); i++) {
        H5Gget_objname_by_idx(grp, (hsize_t)i, memb_name, (size_t)1024 );
        string test( memb_name );
        if ( test.find("region") != std::string::npos ) {
            //patch_name << setfill( '0' ) << setw( 6 ) << region.patch_->Hindex(); -> 6
            //string patchName=Tools::merge( "region-", patch_name.str() );         -> 7
            read_hindex = std::stoi( test.substr(7,6) );
        }
    }

    H5Gclose(grp);
    H5Fclose(file);

    region.vecPatch_.refHindex_ = read_hindex;

}


void Checkpoint::restartPatch( Patch *patch, Params &params, H5Read &g )
{
    ElectroMagn * EMfields = patch->EMfields;

    if ( params.geometry != "AMcylindrical" ) {
        restartFieldsPerProc( g, EMfields->Ex_ );
        restartFieldsPerProc( g, EMfields->Ey_ );
        restartFieldsPerProc( g, EMfields->Ez_ );
        restartFieldsPerProc( g, EMfields->Bx_ );
        restartFieldsPerProc( g, EMfields->By_ );
        restartFieldsPerProc( g, EMfields->Bz_ );
        restartFieldsPerProc( g, EMfields->Bx_m );
        restartFieldsPerProc( g, EMfields->By_m );
        restartFieldsPerProc( g, EMfields->Bz_m );
        if (params.use_BTIS3){
            restartFieldsPerProc( g, EMfields->By_mBTIS3 );
            restartFieldsPerProc( g, EMfields->Bz_mBTIS3 );  
        }
        for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
            if( dynamic_cast<ElectroMagnBC2D_PML *>( EMfields->emBoundCond[bcId] )){
                ElectroMagnBC2D_PML *embc = static_cast<ElectroMagnBC2D_PML *>( EMfields->emBoundCond[bcId] );
                if (embc->Hx_) restart_PML(embc, g);
            } else if( dynamic_cast<ElectroMagnBC3D_PML *>( EMfields->emBoundCond[bcId] )){
                ElectroMagnBC3D_PML *embc = static_cast<ElectroMagnBC3D_PML *>( EMfields->emBoundCond[bcId] );
                if (embc->Hx_) restart_PML(embc, g);
            }
        }

    }
    else {
        for ( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            restart_cFieldsPerProc( g, emAM->El_[imode] );
            restart_cFieldsPerProc( g, emAM->Er_[imode] );
            restart_cFieldsPerProc( g, emAM->Et_[imode] );
            restart_cFieldsPerProc( g, emAM->Bl_[imode] );
            restart_cFieldsPerProc( g, emAM->Br_[imode] );
            restart_cFieldsPerProc( g, emAM->Bt_[imode] );
            restart_cFieldsPerProc( g, emAM->Bl_m[imode] );
            restart_cFieldsPerProc( g, emAM->Br_m[imode] );
            restart_cFieldsPerProc( g, emAM->Bt_m[imode] );

            if (params.use_BTIS3){
                restart_cFieldsPerProc( g, emAM->Br_mBTIS3[imode] );
                restart_cFieldsPerProc( g, emAM->Bt_mBTIS3[imode] );  
            }
            if( params.is_pxr ) {
                restart_cFieldsPerProc( g, emAM->rho_old_AM_[imode] );
            }
            for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
                if( dynamic_cast<ElectroMagnBCAM_PML *>( EMfields->emBoundCond[bcId] )){
                    ElectroMagnBCAM_PML *embc = static_cast<ElectroMagnBCAM_PML *>( EMfields->emBoundCond[bcId] );
                    if (embc->Hl_[0]) restart_PML(embc, g, imode);
                }
            }

        }
    }


    if( EMfields->envelope!=NULL ) {
        DEBUG( "restarting envelope" );
        restart_cFieldsPerProc( g, EMfields->envelope->A_ );
        restart_cFieldsPerProc( g, EMfields->envelope->A0_ );
        restartFieldsPerProc( g, EMfields->Env_Chi_ );
        if (  params.geometry != "AMcylindrical" ) {
            for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
                if( dynamic_cast<EnvelopeBC2D_PML *>( EMfields->envelope->EnvBoundCond[bcId] )){
                    EnvelopeBC2D_PML *envbc = static_cast<EnvelopeBC2D_PML *>( EMfields->envelope->EnvBoundCond[bcId] );
                    if (envbc->A_n_) restart_PMLenvelope(envbc, g, bcId);
                } else if( dynamic_cast<EnvelopeBC3D_PML *>( EMfields->envelope->EnvBoundCond[bcId] )){
                    EnvelopeBC3D_PML *envbc = static_cast<EnvelopeBC3D_PML *>( EMfields->envelope->EnvBoundCond[bcId] );
                    if (envbc->A_n_) restart_PMLenvelope(envbc, g, bcId);
                }
            }
        } else {
            for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
                if( dynamic_cast<EnvelopeBCAM_PML *>( EMfields->envelope->EnvBoundCond[bcId] )){
                    EnvelopeBCAM_PML *envbc = static_cast<EnvelopeBCAM_PML *>( EMfields->envelope->EnvBoundCond[bcId] );
                    if (envbc->A_n_) restart_PMLenvelopeAM(envbc, g, bcId);
                }
            }
        }

    } else {
        DEBUG( "envelope is null" );
    }

    if( EMfields->filter_ ) {
        if (params.geometry!="AMcylindrical"){
            // filtered Electric fields
            for( unsigned int i=0; i<EMfields->filter_->Ex_.size(); i++ ) {
                restartFieldsPerProc( g, EMfields->filter_->Ex_[i] );
            }
            for( unsigned int i=0; i<EMfields->filter_->Ey_.size(); i++ ) {
                restartFieldsPerProc( g, EMfields->filter_->Ey_[i] );
            }
            for( unsigned int i=0; i<EMfields->filter_->Ez_.size(); i++ ) {
                restartFieldsPerProc( g, EMfields->filter_->Ez_[i] );
            }
       } else {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            for ( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {    
                for( unsigned int i=0; i<EMfields->filter_->El_[imode].size(); i++ ) {
                    restart_cFieldsPerProc( g, emAM->filter_->El_[imode][i] );
                }
                for( unsigned int i=0; i<EMfields->filter_->Er_[imode].size(); i++ ) {
                    restart_cFieldsPerProc( g, emAM->filter_->Er_[imode][i] );
                }
                for( unsigned int i=0; i<EMfields->filter_->Et_[imode].size(); i++ ) {
                    restart_cFieldsPerProc( g, emAM->filter_->Et_[imode][i] );
                }
            } // end imode loop  
        } // end if condition on geometry
    }

    // Fields required for DiagFields
    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
        ostringstream group_name( "" );
        group_name << "FieldsForDiag" << idiag;
        if( g.has( group_name.str() ) ) {

            H5Read d = g.group( group_name.str() );
            for( unsigned int ifield=0; ifield<EMfields->allFields_avg[idiag].size(); ifield++ ) {
                restartFieldsPerProc( d, EMfields->allFields_avg[idiag][ifield] );
            }

        } else if( EMfields->allFields_avg[idiag].size() > 0 ) {
            // When the restart occurs in the middle of an average, and the field diag is new,
            // there is missing data that will cause wrong results on the first output after restart
            WARNING( "New average Field diag "<<idiag<<" may produce wrong first output after restart" );
        }
    }

    // Fields required for DiagProbes with time integral
    for( unsigned int iprobe=0; iprobe<patch->probes.size(); iprobe++ ) {
        ostringstream group_name( "" );
        group_name << "DataForProbes" << iprobe;
        if(  g.has( group_name.str() ) ) {
            H5Read diag = g.group( group_name.str() );
            for( unsigned int ifield = 0; true; ifield++ ) {
                ostringstream field( "" );
                field << "field" << ifield;
                if( ! diag.has( field.str() ) ) {
                    break;
                }
                patch->probes[iprobe]->integrated_data.resize(ifield+1);
                diag.vect( field.str(), patch->probes[iprobe]->integrated_data[ifield], true );
            }
        }
    }

    if( ( EMfields->extFields.size()>0 ) && ( params.save_magnectic_fields_for_SM ) ) {
        for( unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
            if( ! EMfields->emBoundCond[bcId] ) {
                continue;
            }
            if( dynamic_cast<ElectroMagnBC1D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC1D_SM *embc = static_cast<ElectroMagnBC1D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName = Tools::merge( "EM_boundary-species-", name.str() );
                H5Read b = g.group( groupName );
                b.attr( "By_val", embc->By_val_ );
                b.attr( "Bz_val", embc->Bz_val_ );
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName = Tools::merge( "EM_boundary-species-", name.str() );
                H5Read b = g.group( groupName );
                b.vect( "Bx_val", embc->B_val[0] );
                b.vect( "By_val", embc->B_val[1] );
                b.vect( "Bz_val", embc->B_val[2] );
            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName = Tools::merge( "EM_boundary-species-", name.str() );
                H5Read b = g.group( groupName );
                if( embc->B_val[0] ) {
                    restartFieldsPerProc( b, embc->B_val[0] );
                }
                if( embc->B_val[1] ) {
                    restartFieldsPerProc( b, embc->B_val[1] );
                }
                if( embc->B_val[2] ) {
                    restartFieldsPerProc( b, embc->B_val[2] );
                }
            }
        }
    }

    unsigned int vecSpeciesSize=0;
    g.attr( "species", vecSpeciesSize );

    if( vecSpeciesSize != patch->vecSpecies.size() ) {
        ERROR_NAMELIST( "Number of species differs between dump (" << vecSpeciesSize << ") and namelist ("<<patch->vecSpecies.size()<<")",
        "https://smileipic.github.io/Smilei/namelist.html#checkpoints");
    }


    for( unsigned int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++ ) {
        Species * spec = patch->vecSpecies[ispec];

        ostringstream name( "" );
        name << setfill( '0' ) << setw( 2 ) << ispec;
        string groupName = Tools::merge( "species-", name.str(), "-", spec->name_ );
        H5Read s = g.group( groupName );

        unsigned int partCapacity=0;
        s.attr( "partCapacity", partCapacity );
        //spec->particles->reserve( partCapacity, nDim_particle );

        unsigned int partSize=0;
        s.attr( "partSize", partSize );
        spec->particles->initialize( partSize, nDim_particle, params.keep_position_old );

        s.attr( "nrj_bc_lost", spec->nrj_bc_lost );
        s.attr( "nrj_mw_inj", spec->nrj_mw_inj );
        s.attr( "nrj_mw_out", spec->nrj_mw_out );
        s.attr( "nrj_new_part", spec->nrj_new_part_ );
        s.attr( "radiatedEnergy", spec->nrj_radiated_ );

        if( partSize>0 ) {
            restartParticles( s, *spec->particles );
            
            if( ! params.cell_sorting_ ) {
                s.vect( "first_index", spec->particles->first_index, true );
                s.vect( "last_index", spec->particles->last_index, true );
            }
            // When cell sorting is activated, indexes are recomputed directly after the restart.
        }
        
        // Read birth records that haven't been written yet
        if( spec->birth_records_ && s.has( "birth_records" ) ) {
            H5Read b = s.group( "birth_records" );
            b.vect( "birth_time", spec->birth_records_->birth_time_, true );
            spec->birth_records_->p_.initialize( spec->birth_records_->birth_time_.size(), *spec->particles );
            restartParticles( b, spec->birth_records_->p_ );
        }
        
    }

    // Load some scalars
    g.attr( "nrj_mw_inj", EMfields->nrj_mw_inj );
    g.attr( "nrj_mw_out", EMfields->nrj_mw_out );
    // Manage some collisions parameters
    if( g.vectSize( "nuclear_reaction_multiplier" ) > 0 ) {
        std::vector<double> rate_multiplier;
        g.vect( "nuclear_reaction_multiplier", rate_multiplier, true );
        for( unsigned int icoll = 0; icoll<rate_multiplier.size(); icoll++ ) {
            if( CollisionalNuclearReaction * NR = patch->vecBPs[icoll]->nuclear_reactions_ ) {
                NR->rate_multiplier_ = rate_multiplier[icoll];
            }
        }
    }

    // Load data for LaserProfileFile (i.e. LaserOffset)
    for( unsigned int ii = 0; ii < 2; ii++ ) {
        if( ! EMfields->emBoundCond[ii] ) continue;
        std::vector<Laser *> * veclaser = & EMfields->emBoundCond[ii]->vecLaser;
        for( unsigned int ilas = 0; ilas < veclaser->size(); ilas++ ) {
            Laser * las = (*veclaser)[ilas];
            for( unsigned int iprof = 0; iprof < las->profiles.size(); iprof++ ) {
                LaserProfile * prof = las->profiles[iprof];
                if( dynamic_cast<LaserProfileFile *>( prof ) ) {
                    LaserProfileFile *p = static_cast<LaserProfileFile *>( prof );
                    if( p->magnitude && p->phase ) {
                        ostringstream t1, t2, t3, t4;
                        t1 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_mag";
                        t2 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_phase";
                        t3 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_omega";
                        t4 << "LaserFile_" << ii << "_" << ilas << "_" << iprof << "_dims";
                        p->magnitude->name = t1.str();
                        p->phase->name = t2.str();
                        vector<unsigned int> dims;
                        g.vect( t4.str(), dims, true );
                        p->magnitude->allocateDims( dims );
                        p->phase->allocateDims( dims );
                        restartFieldsPerProc( g, p->magnitude );
                        restartFieldsPerProc( g, p->phase );
                        g.vect( t3.str(), p->omega, true );
                    }
                }
            }
        }
    }
}

void Checkpoint::dumpFieldsPerProc( H5Write &g, Field *field )
{
    g.vect( field->name, *field->data_, field->number_of_points_, H5T_NATIVE_DOUBLE );
}

void Checkpoint::dump_cFieldsPerProc( H5Write &g, Field *field )
{
    cField *cfield = static_cast<cField *>( field );
    g.vect( field->name, *cfield->cdata_, 2*field->number_of_points_, H5T_NATIVE_DOUBLE );
}

void Checkpoint::restartFieldsPerProc( H5Read &g, Field *field )
{
    g.vect( field->name, *field->data_, H5T_NATIVE_DOUBLE );
}

void Checkpoint::restart_cFieldsPerProc( H5Read &g, Field *field )
{
    cField *cfield = static_cast<cField *>( field );
    g.vect( field->name, *cfield->cdata_, H5T_NATIVE_DOUBLE );
}

void Checkpoint::dumpParticles( H5Write& s, Particles &p )
{
    for( unsigned int i=0; i<p.Position.size(); i++ ) {
        ostringstream my_name( "" );
        my_name << "Position-" << i;
        s.vect( my_name.str(), p.Position[i] );//, dump_deflate );
    }
    
    for( unsigned int i=0; i<p.Momentum.size(); i++ ) {
        ostringstream my_name( "" );
        my_name << "Momentum-" << i;
        s.vect( my_name.str(),p.Momentum[i] );//, dump_deflate );
    }
    
    s.vect( "Weight", p.Weight );//, dump_deflate );
    s.vect( "Charge", p.Charge );//, dump_deflate );
    
    if( p.tracked ) {
        s.vect( "Id", p.Id, H5T_NATIVE_UINT64 );//, dump_deflate );
    }
    
    // Monte-Carlo process
    if( p.has_Monte_Carlo_process ) {
        s.vect( "Tau", p.Tau );//, dump_deflate );
    }
    
    // Copy interpolated fields that must be accumulated over time
    if( p.interpolated_fields_ ) {
        if( p.interpolated_fields_->mode_[6] == 2 ) {
            s.vect( "Wx", p.interpolated_fields_->F_[6] );
        }
        if( p.interpolated_fields_->mode_[7] == 2 ) {
            s.vect( "Wy", p.interpolated_fields_->F_[7] );
        }
        if( p.interpolated_fields_->mode_[8] == 2 ) {
            s.vect( "Wz", p.interpolated_fields_->F_[8] );
        }
    }
}

void Checkpoint::restartParticles( H5Read& s, Particles &p )
{
    for( unsigned int i=0; i<p.Position.size(); i++ ) {
        ostringstream my_name( "" );
        my_name << "Position-" << i;
        s.vect( my_name.str(), p.Position[i] );
    }
    
    for( unsigned int i=0; i<p.Momentum.size(); i++ ) {
        ostringstream my_name( "" );
        my_name << "Momentum-" << i;
        s.vect( my_name.str(), p.Momentum[i] );
    }
    
    s.vect( "Weight", p.Weight );
    s.vect( "Charge", p.Charge );
    
    if( p.tracked && s.has( "Id" ) ) {
        s.vect( "Id", p.Id, H5T_NATIVE_UINT64 );
    }
    
    if( p.has_Monte_Carlo_process ) {
        s.vect( "Tau", p.Tau );
    }
    
    // Retrieve interpolated fields that must be accumulated over time
    if( p.interpolated_fields_ ) {
        if( p.interpolated_fields_->mode_[6] == 2 ) {
            s.vect( "Wx", p.interpolated_fields_->F_[6] );
        }
        if( p.interpolated_fields_->mode_[7] == 2 ) {
            s.vect( "Wy", p.interpolated_fields_->F_[7] );
        }
        if( p.interpolated_fields_->mode_[8] == 2 ) {
            s.vect( "Wz", p.interpolated_fields_->F_[8] );
        }
    }
}

void Checkpoint::dumpMovingWindow( H5Write &f, SimWindow *simWin )
{
    f.attr( "x_moved", simWin->getXmoved() );
    f.attr( "n_moved", simWin->getNmoved() );
}
void Checkpoint::restartMovingWindow( H5Read &f, SimWindow *simWin )
{

    double x_moved=0.;
    f.attr( "x_moved", x_moved );
    simWin->setXmoved( x_moved );

    unsigned int n_moved=0;
    f.attr( "n_moved", n_moved );
    simWin->setNmoved( n_moved );

}
template <typename Tpml> //ElectroMagnBC2D_PML or ElectroMagnBC3D_PML
void  Checkpoint::dump_PML(Tpml embc, H5Write &g ){
    dumpFieldsPerProc( g, embc->Hx_ );
    dumpFieldsPerProc( g, embc->Hy_ );
    dumpFieldsPerProc( g, embc->Hz_ );
    dumpFieldsPerProc( g, embc->Bx_ );
    dumpFieldsPerProc( g, embc->By_ );
    dumpFieldsPerProc( g, embc->Bz_ );
    dumpFieldsPerProc( g, embc->Ex_ );
    dumpFieldsPerProc( g, embc->Ey_ );
    dumpFieldsPerProc( g, embc->Ez_ );
    dumpFieldsPerProc( g, embc->Dx_ );
    dumpFieldsPerProc( g, embc->Dy_ );
    dumpFieldsPerProc( g, embc->Dz_ );
}
void  Checkpoint::dump_PML( ElectroMagnBCAM_PML *embc, H5Write &g, unsigned int imode ){
    dump_cFieldsPerProc( g, embc->Hl_[imode] );
    dump_cFieldsPerProc( g, embc->Hr_[imode] );
    dump_cFieldsPerProc( g, embc->Ht_[imode] );
    dump_cFieldsPerProc( g, embc->Bl_[imode] );
    dump_cFieldsPerProc( g, embc->Br_[imode] );
    dump_cFieldsPerProc( g, embc->Bt_[imode] );
    dump_cFieldsPerProc( g, embc->El_[imode] );
    dump_cFieldsPerProc( g, embc->Er_[imode] );
    dump_cFieldsPerProc( g, embc->Et_[imode] );
    dump_cFieldsPerProc( g, embc->Dl_[imode] );
    dump_cFieldsPerProc( g, embc->Dr_[imode] );
    dump_cFieldsPerProc( g, embc->Dt_[imode] );
}
template <typename Tpml> //EnvelopBC2D_PML or EnvelopeBC3D_PML
void  Checkpoint::dump_PMLenvelope(Tpml envbc, H5Write &g, unsigned int bcId ){
    dumpFieldsPerProc( g, envbc->Chi_ );
    dump_cFieldsPerProc( g, envbc->A_n_ );
    dump_cFieldsPerProc( g, envbc->A_nm1_ );
    dump_cFieldsPerProc( g, envbc->u1_nm1_x_ );
    dump_cFieldsPerProc( g, envbc->u2_nm1_x_ );
    dump_cFieldsPerProc( g, envbc->u3_nm1_x_ );
    if ( bcId > 1) {
        dump_cFieldsPerProc( g, envbc->u1_nm1_y_ );
        dump_cFieldsPerProc( g, envbc->u2_nm1_y_ );
        dump_cFieldsPerProc( g, envbc->u3_nm1_y_ );
    }
    if ( std::is_same<Tpml, EnvelopeBC3D_PML>::value and bcId > 3) {
        EnvelopeBC3D_PML *envbc3d = dynamic_cast<EnvelopeBC3D_PML *>( envbc );
        dump_cFieldsPerProc( g, envbc3d->u1_nm1_z_ );
        dump_cFieldsPerProc( g, envbc3d->u2_nm1_z_ );
        dump_cFieldsPerProc( g, envbc3d->u3_nm1_z_ );
    }
}
void  Checkpoint::dump_PMLenvelopeAM(EnvelopeBCAM_PML *envbc, H5Write &g, unsigned int bcId ){
    dumpFieldsPerProc( g, envbc->Chi_ );
    dump_cFieldsPerProc( g, envbc->A_n_ );
    dump_cFieldsPerProc( g, envbc->A_nm1_ );
    dump_cFieldsPerProc( g, envbc->G_n_ );
    dump_cFieldsPerProc( g, envbc->G_nm1_ );
    dump_cFieldsPerProc( g, envbc->u1_nm1_l_ );
    dump_cFieldsPerProc( g, envbc->u2_nm1_l_ );
    dump_cFieldsPerProc( g, envbc->u3_nm1_l_ );
    if ( bcId == 3) {
        dump_cFieldsPerProc( g, envbc->u1_nm1_r_ );
        dump_cFieldsPerProc( g, envbc->u2_nm1_r_ );
        dump_cFieldsPerProc( g, envbc->u3_nm1_r_ );
    }
}

template <typename Tpml> //ElectroMagnBC2D_PML or ElectroMagnBC3D_PML
void  Checkpoint::restart_PML(Tpml embc, H5Read &g ){
    restartFieldsPerProc( g, embc->Hx_ );
    restartFieldsPerProc( g, embc->Hy_ );
    restartFieldsPerProc( g, embc->Hz_ );
    restartFieldsPerProc( g, embc->Bx_ );
    restartFieldsPerProc( g, embc->By_ );
    restartFieldsPerProc( g, embc->Bz_ );
    restartFieldsPerProc( g, embc->Ex_ );
    restartFieldsPerProc( g, embc->Ey_ );
    restartFieldsPerProc( g, embc->Ez_ );
    restartFieldsPerProc( g, embc->Dx_ );
    restartFieldsPerProc( g, embc->Dy_ );
    restartFieldsPerProc( g, embc->Dz_ );
}
void  Checkpoint::restart_PML(ElectroMagnBCAM_PML *embc, H5Read &g, unsigned int imode ){
    restart_cFieldsPerProc( g, embc->Hl_[imode] );
    restart_cFieldsPerProc( g, embc->Hr_[imode] );
    restart_cFieldsPerProc( g, embc->Ht_[imode] );
    restart_cFieldsPerProc( g, embc->Bl_[imode] );
    restart_cFieldsPerProc( g, embc->Br_[imode] );
    restart_cFieldsPerProc( g, embc->Bt_[imode] );
    restart_cFieldsPerProc( g, embc->El_[imode] );
    restart_cFieldsPerProc( g, embc->Er_[imode] );
    restart_cFieldsPerProc( g, embc->Et_[imode] );
    restart_cFieldsPerProc( g, embc->Dl_[imode] );
    restart_cFieldsPerProc( g, embc->Dr_[imode] );
    restart_cFieldsPerProc( g, embc->Dt_[imode] );
}
template <typename Tpml> //EnvelopBC2D_PML or EnvelopeBC3D_PML
void  Checkpoint::restart_PMLenvelope(Tpml envbc, H5Read &g, unsigned int bcId ){
    restartFieldsPerProc( g, envbc->Chi_ );
    restart_cFieldsPerProc( g, envbc->A_n_ );
    restart_cFieldsPerProc( g, envbc->A_nm1_ );
    restart_cFieldsPerProc( g, envbc->u1_nm1_x_ );
    restart_cFieldsPerProc( g, envbc->u2_nm1_x_ );
    restart_cFieldsPerProc( g, envbc->u3_nm1_x_ );
    if ( bcId > 1) {
        restart_cFieldsPerProc( g, envbc->u1_nm1_y_ );
        restart_cFieldsPerProc( g, envbc->u2_nm1_y_ );
        restart_cFieldsPerProc( g, envbc->u3_nm1_y_ );
    }
    if ( std::is_same<Tpml, EnvelopeBC3D_PML>::value and bcId > 3) {
        EnvelopeBC3D_PML *envbc3d = dynamic_cast<EnvelopeBC3D_PML *>( envbc );
        restart_cFieldsPerProc( g, envbc3d->u1_nm1_z_ );
        restart_cFieldsPerProc( g, envbc3d->u2_nm1_z_ );
        restart_cFieldsPerProc( g, envbc3d->u3_nm1_z_ );
    }
}
void  Checkpoint::restart_PMLenvelopeAM(EnvelopeBCAM_PML *envbc, H5Read &g, unsigned int bcId){
    restartFieldsPerProc( g, envbc->Chi_ );
    restart_cFieldsPerProc( g, envbc->A_n_ );
    restart_cFieldsPerProc( g, envbc->A_nm1_ );
    restart_cFieldsPerProc( g, envbc->G_n_ );
    restart_cFieldsPerProc( g, envbc->G_nm1_ );
    restart_cFieldsPerProc( g, envbc->u1_nm1_l_ );
    restart_cFieldsPerProc( g, envbc->u2_nm1_l_ );
    restart_cFieldsPerProc( g, envbc->u3_nm1_l_ );
    if ( bcId == 3) {
        restart_cFieldsPerProc( g, envbc->u1_nm1_r_ );
        restart_cFieldsPerProc( g, envbc->u2_nm1_r_ );
        restart_cFieldsPerProc( g, envbc->u3_nm1_r_ );
    }
}
