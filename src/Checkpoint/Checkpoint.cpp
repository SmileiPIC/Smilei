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
#include "OpenPMDparams.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "Region.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC2D_SM.h"
#include "ElectroMagnBC3D_SM.h"
#include "Laser.h"
#include "Species.h"
#include "PatchesFactory.h"
#include "DiagnosticScreen.h"
#include "DiagnosticTrack.h"
#include "LaserEnvelope.h"
#include "Collisions.h"

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
    time_dump_step( 0 ),
    keep_n_dumps( 2 ),
    keep_n_dumps_max( 10000 ),
    dump_deflate( 0 ),
    dump_request( smpi->getSize() ),
    file_grouping( 0 )
{

    if( PyTools::nComponents( "Checkpoints" ) > 0 ) {
    
        PyTools::extract( "dump_step", dump_step, "Checkpoints"  );
        if( dump_step > 0 ) {
            MESSAGE( 1, "Code will dump after " << dump_step << " steps" );
        }
        
        PyTools::extract( "dump_minutes", dump_minutes, "Checkpoints"  );
        if( dump_minutes > 0 ) {
            MESSAGE( 1, "Code will stop after " << dump_minutes << " minutes" );
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
                H5Read f( dump_name, false, false );
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
                    smpi->SMILEI_COMM_WORLD, &status
                );
                int problem = (prev_number != dump_number);
                int any_problem;
                MPI_Allreduce( &problem, &any_problem, 1, MPI_INT, MPI_LOR, smpi->SMILEI_COMM_WORLD );
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

void Checkpoint::dump( VectorPatch &vecPatches, Region &region, unsigned int itime, SmileiMPI *smpi, SimWindow *simWindow, Params &params )
{

    // check for excedeed time
    if( dump_minutes != 0.0 ) {
        int tagUB = smpi->getTagUB();
        // master checks whenever we passed the time limit
        if( smpi->isMaster() && time_dump_step==0 ) {
            double elapsed_time = ( MPI_Wtime() - time_reference )/60.;
            if( elapsed_time > dump_minutes ) {
                time_dump_step = itime+1; // we will dump at next timestep (in case non-master already passed)
                MESSAGE( "Reached time limit : " << elapsed_time << " minutes. Dump timestep : " << time_dump_step );
                // master does a non-blocking send
                for( unsigned int dest=0; dest < ( unsigned int ) smpi->getSize(); dest++ ) {
                    MPI_Isend( &time_dump_step, 1, MPI_UNSIGNED, dest, tagUB, smpi->SMILEI_COMM_WORLD, &dump_request[dest] );
                }
            }
        } else { // non master nodes receive the time_dump_step (non-blocking)
            int todump=0;
            MPI_Iprobe( 0, tagUB, MPI_COMM_WORLD, &todump, &dump_status_prob );
            if( todump ) {
                MPI_Recv( &time_dump_step, 1, MPI_UNSIGNED, 0, tagUB, smpi->SMILEI_COMM_WORLD, &dump_status_recv );
            }
        }
        smpi->barrier();
    }
    
    if( signal_received!=0 ||
            ( dump_step != 0 && ( ( itime-this_run_start_step ) % dump_step == 0 ) ) ||
            ( time_dump_step!=0 && itime==time_dump_step ) ) {
        dumpAll( vecPatches, region, itime,  smpi, simWindow, params );
        if( exit_after_dump || ( ( signal_received!=0 ) && ( signal_received != SIGUSR2 ) ) ) {
            exit_asap=true;
        }
        signal_received=0;
        time_dump_step=0;
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
    MESSAGEALL( "Step " << itime << " : DUMP fields and particles " << dumpName );
#else
    MESSAGE( "Step " << itime << " : DUMP fields and particles " << num_dump );
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
        // Poynting scalars
        unsigned int k=0;
        for( unsigned int j=0; j<2; j++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
            for( unsigned int i=0; i<params.nDim_field; i++ ) { //axis 0=x, 1=y, 2=z
                if( scalars->necessary_poy[k] ) {
                    string poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
                    f.attr( poy_name, ( double )*( scalars->poy[k] ) );
                    k++;
                }
            }
        }
    }
    
    // Write the diags screen data
    ostringstream diagName( "" );
    if( smpi->isMaster() ) {
        unsigned int iscreen = 0;
        for( unsigned int idiag=0; idiag<vecPatches.globalDiags.size(); idiag++ ) {
            if( DiagnosticScreen *screen = dynamic_cast<DiagnosticScreen *>( vecPatches.globalDiags[idiag] ) ) {
                diagName.str( "" );
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
        
        dumpPatch( vecPatches( ipatch )->EMfields, vecPatches( ipatch )->vecSpecies, vecPatches( ipatch )->vecCollisions, params, g );
        
        // Random number generator state
        g.attr( "xorshift32_state", vecPatches( ipatch )->rand_->xorshift32_state );
        
    }

    if (params.uncoupled_grids) {
        // Open a group
        ostringstream patch_name( "" );
        patch_name << setfill( '0' ) << setw( 6 ) << region.patch_->Hindex();
        string patchName=Tools::merge( "region-", patch_name.str() );
        H5Write g = f.group( patchName.c_str() );
        dumpPatch( region.patch_->EMfields, region.patch_->vecSpecies, region.patch_->vecCollisions, params, g );
    }

    // Write the latest Id that the MPI processes have given to each species
    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ ) {
        if( DiagnosticTrack *track = dynamic_cast<DiagnosticTrack *>( vecPatches.localDiags[idiag] ) ) {
            ostringstream n( "" );
            n<< "latest_ID_" << vecPatches( 0 )->vecSpecies[track->speciesId_]->name_;
            f.attr( n.str(), track->latest_Id, H5T_NATIVE_UINT64 );
        }
    }
    
    // Write the moving window status
    if( simWin!=NULL ) {
        dumpMovingWindow( f, simWin );
    }
    
}


void Checkpoint::dumpPatch( ElectroMagn *EMfields, std::vector<Species *> vecSpecies, std::vector<Collisions *> &vecCollisions, Params &params, H5Write &g )
{
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
            
            if(params.is_pxr == true)
                dump_cFieldsPerProc( g, emAM->rho_old_AM_[imode] );
            
        }
    }
    
    if( EMfields->envelope!=NULL ) {
        dump_cFieldsPerProc( g, EMfields->envelope->A_ );
        dump_cFieldsPerProc( g, EMfields->envelope->A0_ );
        dumpFieldsPerProc( g, EMfields->Env_Chi_ );
    }
    
    // filtered Electric fields
    for( unsigned int i=0; i<EMfields->Exfilter.size(); i++ ) {
        dumpFieldsPerProc( g, EMfields->Exfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Eyfilter.size(); i++ ) {
        dumpFieldsPerProc( g, EMfields->Eyfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Ezfilter.size(); i++ ) {
        dumpFieldsPerProc( g, EMfields->Ezfilter[i] );
    }
    // filtered Magnetic fields
    for( unsigned int i=0; i<EMfields->Bxfilter.size(); i++ ) {
        dumpFieldsPerProc( g, EMfields->Bxfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Byfilter.size(); i++ ) {
        dumpFieldsPerProc( g, EMfields->Byfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Bzfilter.size(); i++ ) {
        dumpFieldsPerProc( g, EMfields->Bzfilter[i] );
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
                b.attr( "By_val", embc->By_val );
                b.attr( "Bz_val", embc->Bz_val );
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                H5Write b = g.group( groupName );
                g.vect( "Bx_val", embc->Bx_val );
                g.vect( "By_val", embc->By_val );
                g.vect( "Bz_val", embc->Bz_val );
            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                H5Write b = g.group( groupName );
                if( embc->Bx_val ) {
                    dumpFieldsPerProc( b, embc->Bx_val );
                }
                if( embc->By_val ) {
                    dumpFieldsPerProc( b, embc->By_val );
                }
                if( embc->Bz_val ) {
                    dumpFieldsPerProc( b, embc->Bz_val );
                }
            }
        }
    }
    
    g.flush();
    g.attr( "species", vecSpecies.size() );
    
    for( unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++ ) {
        ostringstream name( "" );
        name << setfill( '0' ) << setw( 2 ) << ispec;
        string groupName=Tools::merge( "species-", name.str(), "-", vecSpecies[ispec]->name_ );
        H5Write s = g.group( groupName );
        
        s.attr( "partCapacity", vecSpecies[ispec]->particles->capacity() );
        s.attr( "partSize", vecSpecies[ispec]->particles->size() );
        s.attr( "nrj_radiation", vecSpecies[ispec]->getNrjRadiation() );
        
        if( vecSpecies[ispec]->particles->size()>0 ) {
        
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++ ) {
                ostringstream my_name( "" );
                my_name << "Position-" << i;
                s.vect( my_name.str(), vecSpecies[ispec]->particles->Position[i] );//, dump_deflate );
            }
            
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++ ) {
                ostringstream my_name( "" );
                my_name << "Momentum-" << i;
                s.vect( my_name.str(), vecSpecies[ispec]->particles->Momentum[i] );//, dump_deflate );
            }
            
            s.vect( "Weight", vecSpecies[ispec]->particles->Weight );//, dump_deflate );
            s.vect( "Charge", vecSpecies[ispec]->particles->Charge );//, dump_deflate );
            
            if( vecSpecies[ispec]->particles->tracked ) {
                s.vect( "Id", vecSpecies[ispec]->particles->Id, H5T_NATIVE_UINT64 );//, dump_deflate );
            }
            
            s.vect( "first_index", vecSpecies[ispec]->particles->first_index );
            s.vect( "last_index", vecSpecies[ispec]->particles->last_index );
            
        } // End if partSize
        
    } // End for ispec
    
    // Manage some collisions parameters
    std::vector<double> rate_multiplier( vecCollisions.size() );
    for( unsigned int icoll = 0; icoll<vecCollisions.size(); icoll++ ) {
        rate_multiplier[icoll] = vecCollisions[icoll]->NuclearReaction->rate_multiplier_;
    }
    g.vect( "collisions_rate_multiplier", rate_multiplier );
    
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
    
    string dump_date;
    f.attr( "CommitDate", dump_date );
    
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


void Checkpoint::restartAll( VectorPatch &vecPatches, Region &region, SmileiMPI *smpi, SimWindow *simWin, Params &params, OpenPMDparams &openPMD )
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
    }
    
    // Read the diags screen data
    ostringstream diagName( "" );
    if( smpi->isMaster() ) {
        unsigned int iscreen = 0;
        for( unsigned int idiag=0; idiag<vecPatches.globalDiags.size(); idiag++ ) {
            if( DiagnosticScreen *screen = dynamic_cast<DiagnosticScreen *>( vecPatches.globalDiags[idiag] ) ) {
                diagName.str( "" );
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
        
        restartPatch( vecPatches( ipatch )->EMfields, vecPatches( ipatch )->vecSpecies, vecPatches( ipatch )->vecCollisions, params, g );
        
        // Random number generator state
        g.attr( "xorshift32_state", vecPatches( ipatch )->rand_->xorshift32_state );
        
    }

    if (params.uncoupled_grids) {
        ostringstream patch_name( "" );
        patch_name << setfill( '0' ) << setw( 6 ) << region.patch_->Hindex();
        string patchName = Tools::merge( "region-", patch_name.str() );
        H5Read g = f.group( patchName );
        restartPatch( region.patch_->EMfields, region.patch_->vecSpecies, region.patch_->vecCollisions, params, g );
    }
    
    // Read the latest Id that the MPI processes have given to each species
    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ ) {
        if( DiagnosticTrack *track = dynamic_cast<DiagnosticTrack *>( vecPatches.localDiags[idiag] ) ) {
            ostringstream n( "" );
            n<< "latest_ID_" << vecPatches( 0 )->vecSpecies[track->speciesId_]->name_;
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
    for (int i = 0; i < nobj; i++) {
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


void Checkpoint::restartPatch( ElectroMagn *EMfields, std::vector<Species *> &vecSpecies, std::vector<Collisions *> &vecCollisions, Params &params, H5Read &g )
{
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
            
            if(params.is_pxr == true)
                restart_cFieldsPerProc( g, emAM->rho_old_AM_[imode] );
            
        }
    }
    

    if( EMfields->envelope!=NULL ) {
        DEBUG( "restarting envelope" );
        restart_cFieldsPerProc( g, EMfields->envelope->A_ );
        restart_cFieldsPerProc( g, EMfields->envelope->A0_ );
        restartFieldsPerProc( g, EMfields->Env_Chi_ );
    } else {
        DEBUG( "envelope is null" );
    }
    
    
    // filtered Electric fields
    for( unsigned int i=0; i<EMfields->Exfilter.size(); i++ ) {
        restartFieldsPerProc( g, EMfields->Exfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Eyfilter.size(); i++ ) {
        restartFieldsPerProc( g, EMfields->Eyfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Ezfilter.size(); i++ ) {
        restartFieldsPerProc( g, EMfields->Ezfilter[i] );
    }
    // filtered Magnetic fields
    for( unsigned int i=0; i<EMfields->Bxfilter.size(); i++ ) {
        restartFieldsPerProc( g, EMfields->Bxfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Byfilter.size(); i++ ) {
        restartFieldsPerProc( g, EMfields->Byfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Bzfilter.size(); i++ ) {
        restartFieldsPerProc( g, EMfields->Bzfilter[i] );
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
                b.attr( "By_val", embc->By_val );
                b.attr( "Bz_val", embc->Bz_val );
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName = Tools::merge( "EM_boundary-species-", name.str() );
                H5Read b = g.group( groupName );
                b.vect( "Bx_val", embc->Bx_val );
                b.vect( "By_val", embc->By_val );
                b.vect( "Bz_val", embc->Bz_val );
            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName = Tools::merge( "EM_boundary-species-", name.str() );
                H5Read b = g.group( groupName );
                if( embc->Bx_val ) {
                    restartFieldsPerProc( b, embc->Bx_val );
                }
                if( embc->By_val ) {
                    restartFieldsPerProc( b, embc->By_val );
                }
                if( embc->Bz_val ) {
                    restartFieldsPerProc( b, embc->Bz_val );
                }
            }
        }
    }
    
    unsigned int vecSpeciesSize=0;
    g.attr( "species", vecSpeciesSize );
    
    if( vecSpeciesSize != vecSpecies.size() ) {
        ERROR( "Number of species differs between dump (" << vecSpeciesSize << ") and namelist ("<<vecSpecies.size()<<")" );
    }
    
    
    for( unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++ ) {
        ostringstream name( "" );
        name << setfill( '0' ) << setw( 2 ) << ispec;
        string groupName = Tools::merge( "species-", name.str(), "-", vecSpecies[ispec]->name_ );
        H5Read s = g.group( groupName );
        
        unsigned int partCapacity=0;
        s.attr( "partCapacity", partCapacity );
        vecSpecies[ispec]->particles->reserve( partCapacity, nDim_particle );
        
        unsigned int partSize=0;
        s.attr( "partSize", partSize );
        vecSpecies[ispec]->particles->initialize( partSize, nDim_particle, params.keep_position_old );
        
        double nrj_radiation;
        if( s.hasAttr( "nrj_radiation" ) ) {
            s.attr( "nrj_radiation", nrj_radiation );
            vecSpecies[ispec]->setNrjRadiation( nrj_radiation );
        }
        
        if( partSize>0 ) {
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++ ) {
                ostringstream namePos( "" );
                namePos << "Position-" << i;
                s.vect( namePos.str(), vecSpecies[ispec]->particles->Position[i] );
            }
            
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++ ) {
                ostringstream namePos( "" );
                namePos << "Momentum-" << i;
                s.vect( namePos.str(), vecSpecies[ispec]->particles->Momentum[i] );
            }
            
            s.vect( "Weight", vecSpecies[ispec]->particles->Weight );
            
            s.vect( "Charge", vecSpecies[ispec]->particles->Charge );
            
            if( vecSpecies[ispec]->particles->tracked ) {
                s.vect( "Id", vecSpecies[ispec]->particles->Id, H5T_NATIVE_UINT64 );
            }
            
            if( params.vectorization_mode == "off" || params.vectorization_mode == "on" || params.cell_sorting ) {
                s.vect( "first_index", vecSpecies[ispec]->particles->first_index, true );
                s.vect( "last_index", vecSpecies[ispec]->particles->last_index, true );
            }
            // In the adaptive vectorization case, the bins will be recomputed
            // latter in the patch reconfiguration
            
        }
    }
    
    // Manage some collisions parameters
    if( g.vectSize( "collisions_rate_multiplier" ) > 0 ) {
        std::vector<double> rate_multiplier;
        g.vect( "collisions_rate_multiplier", rate_multiplier, true );
        for( unsigned int icoll = 0; icoll<rate_multiplier.size(); icoll++ ) {
            vecCollisions[icoll]->NuclearReaction->rate_multiplier_ = rate_multiplier[icoll];
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
    g.vect( field->name, *field->data_, field->globalDims_, H5T_NATIVE_DOUBLE );
}

void Checkpoint::dump_cFieldsPerProc( H5Write &g, Field *field )
{
    cField *cfield = static_cast<cField *>( field );
    g.vect( field->name, *cfield->cdata_, 2*field->globalDims_, H5T_NATIVE_DOUBLE );
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
