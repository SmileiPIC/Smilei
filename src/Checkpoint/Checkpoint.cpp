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
                string dump_name=restart_files[num_dump];
                hid_t fid = H5::Fopen( dump_name );
                if( fid < 0 ) {
                    continue;
                }
                unsigned int stepStartTmp=0;
                H5::getAttr( fid, "dump_step", stepStartTmp );
                if( stepStartTmp>this_run_start_step ) {
                    this_run_start_step=stepStartTmp;
                    restart_file=dump_name;
                    dump_number=num_dump;
                    H5::getAttr( fid, "dump_number", dump_number );
                }
                H5Fclose( fid );
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

void Checkpoint::dump( VectorPatch &vecPatches, unsigned int itime, SmileiMPI *smpi, SimWindow *simWindow, Params &params )
{

    // check for excedeed time
    if( dump_minutes != 0.0 ) {
        // master checks whenever we passed the time limit
        if( smpi->isMaster() && time_dump_step==0 ) {
            double elapsed_time = ( MPI_Wtime() - time_reference )/60.;
            if( elapsed_time > dump_minutes ) {
                time_dump_step = itime+1; // we will dump at next timestep (in case non-master already passed)
                MESSAGE( "Reached time limit : " << elapsed_time << " minutes. Dump timestep : " << time_dump_step );
                // master does a non-blocking send
                for( unsigned int dest=0; dest < ( unsigned int ) smpi->getSize(); dest++ ) {
                    MPI_Isend( &time_dump_step, 1, MPI_UNSIGNED, dest, SMILEI_COMM_DUMP_TIME, smpi->SMILEI_COMM_WORLD, &dump_request[dest] );
                }
            }
        } else { // non master nodes receive the time_dump_step (non-blocking)
            int todump=0;
            MPI_Iprobe( 0, SMILEI_COMM_DUMP_TIME, MPI_COMM_WORLD, &todump, &dump_status_prob );
            if( todump ) {
                MPI_Recv( &time_dump_step, 1, MPI_UNSIGNED, 0, SMILEI_COMM_DUMP_TIME, smpi->SMILEI_COMM_WORLD, &dump_status_recv );
            }
        }
        smpi->barrier();
    }
    
    if( signal_received!=0 ||
            ( dump_step != 0 && ( ( itime-this_run_start_step ) % dump_step == 0 ) ) ||
            ( time_dump_step!=0 && itime==time_dump_step ) ) {
        dumpAll( vecPatches, itime,  smpi, simWindow, params );
        if( exit_after_dump || ( ( signal_received!=0 ) && ( signal_received != SIGUSR2 ) ) ) {
            exit_asap=true;
        }
        signal_received=0;
        time_dump_step=0;
        time_reference = MPI_Wtime();
    }
}

void Checkpoint::dumpAll( VectorPatch &vecPatches, unsigned int itime,  SmileiMPI *smpi, SimWindow *simWin,  Params &params )
{
    unsigned int num_dump=dump_number % keep_n_dumps;
    
    ostringstream nameDumpTmp( "" );
    nameDumpTmp << "checkpoints" << PATH_SEPARATOR;
    if( file_grouping>0 ) {
        nameDumpTmp << setfill( '0' ) << setw( int( 1+log10( smpi->getSize()/file_grouping+1 ) ) ) << smpi->getRank()/file_grouping << PATH_SEPARATOR;
    }
    
    nameDumpTmp << "dump-" << setfill( '0' ) << setw( 5 ) << num_dump << "-" << setfill( '0' ) << setw( 10 ) << smpi->getRank() << ".h5" ;
    std::string dumpName=nameDumpTmp.str();
    
    
    hid_t fid = H5Fcreate( dumpName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    if (fid<0) {
        ERROR("Can't open file for writing checkpoint " << dumpName.c_str())
    } else {
        dump_number++;
    }
    
#ifdef  __DEBUG
    MESSAGEALL( "Step " << itime << " : DUMP fields and particles " << dumpName );
#else
    MESSAGE( "Step " << itime << " : DUMP fields and particles " << num_dump );
#endif
    
    
    // Write basic attributes
    H5::attr( fid, "Version", string( __VERSION ) );
    
    H5::attr( fid, "dump_step", itime );
    H5::attr( fid, "dump_number", dump_number );
    
    H5::vect( fid, "patch_count", smpi->patch_count );
    
    // Write diags scalar data
    DiagnosticScalar *scalars = static_cast<DiagnosticScalar *>( vecPatches.globalDiags[0] );
    H5::attr( fid, "latest_timestep",   scalars->latest_timestep );
    // Scalars only by master
    if( smpi->isMaster() ) {
        H5::attr( fid, "Energy_time_zero",  scalars->Energy_time_zero );
        H5::attr( fid, "EnergyUsedForNorm", scalars->EnergyUsedForNorm );
        // Poynting scalars
        unsigned int k=0;
        for( unsigned int j=0; j<2; j++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
            for( unsigned int i=0; i<params.nDim_field; i++ ) { //axis 0=x, 1=y, 2=z
                if( scalars->necessary_poy[k] ) {
                    string poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
                    H5::attr( fid, poy_name, ( double )*( scalars->poy[k] ) );
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
                H5::vect( fid, diagName.str(), *(screen->getData()) );
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
        hid_t patch_gid = H5::group( fid, patchName.c_str() );
        
        dumpPatch( vecPatches( ipatch )->EMfields, vecPatches( ipatch )->vecSpecies, vecPatches( ipatch )->vecCollisions, params, patch_gid );
        
        // Random number generator state
        H5::attr( patch_gid, "xorshift32_state", vecPatches( ipatch )->rand_->xorshift32_state );
        
        // Close a group
        H5Gclose( patch_gid );
        
    }
    
    // Write the latest Id that the MPI processes have given to each species
    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ ) {
        if( DiagnosticTrack *track = dynamic_cast<DiagnosticTrack *>( vecPatches.localDiags[idiag] ) ) {
            ostringstream n( "" );
            n<< "latest_ID_" << vecPatches( 0 )->vecSpecies[track->speciesId_]->name_;
            H5::attr( fid, n.str().c_str(), track->latest_Id, H5T_NATIVE_UINT64 );
        }
    }
    
    // Write the moving window status
    if( simWin!=NULL ) {
        dumpMovingWindow( fid, simWin );
    }
    
    herr_t tclose = H5Fclose( fid );
    if (tclose < 0) {
        ERROR("Can't close file " << dumpName.c_str())
    }
    
}

void Checkpoint::dumpPatch( ElectroMagn *EMfields, std::vector<Species *> vecSpecies, std::vector<Collisions *> &vecCollisions, Params &params, hid_t patch_gid )
{
    if (  params.geometry != "AMcylindrical" ) {
        dumpFieldsPerProc( patch_gid, EMfields->Ex_ );
        dumpFieldsPerProc( patch_gid, EMfields->Ey_ );
        dumpFieldsPerProc( patch_gid, EMfields->Ez_ );
        dumpFieldsPerProc( patch_gid, EMfields->Bx_ );
        dumpFieldsPerProc( patch_gid, EMfields->By_ );
        dumpFieldsPerProc( patch_gid, EMfields->Bz_ );
        dumpFieldsPerProc( patch_gid, EMfields->Bx_m );
        dumpFieldsPerProc( patch_gid, EMfields->By_m );
        dumpFieldsPerProc( patch_gid, EMfields->Bz_m );
    }
    else {
        for ( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            dump_cFieldsPerProc( patch_gid, emAM->El_[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Er_[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Et_[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Bl_[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Br_[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Bt_[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Bl_m[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Br_m[imode] );
            dump_cFieldsPerProc( patch_gid, emAM->Bt_m[imode] );
            
            if(params.is_pxr == true)
                dump_cFieldsPerProc( patch_gid, emAM->rho_old_AM_[imode] );
            
        }
    }
    
    if( EMfields->envelope!=NULL ) {
        dump_cFieldsPerProc( patch_gid, EMfields->envelope->A_ );
        dump_cFieldsPerProc( patch_gid, EMfields->envelope->A0_ );
        dumpFieldsPerProc( patch_gid, EMfields->Env_Chi_ );
    }
    
    // filtered Electric fields
    for( unsigned int i=0; i<EMfields->Exfilter.size(); i++ ) {
        dumpFieldsPerProc( patch_gid, EMfields->Exfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Eyfilter.size(); i++ ) {
        dumpFieldsPerProc( patch_gid, EMfields->Eyfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Ezfilter.size(); i++ ) {
        dumpFieldsPerProc( patch_gid, EMfields->Ezfilter[i] );
    }
    // filtered Magnetic fields
    for( unsigned int i=0; i<EMfields->Bxfilter.size(); i++ ) {
        dumpFieldsPerProc( patch_gid, EMfields->Bxfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Byfilter.size(); i++ ) {
        dumpFieldsPerProc( patch_gid, EMfields->Byfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Bzfilter.size(); i++ ) {
        dumpFieldsPerProc( patch_gid, EMfields->Bzfilter[i] );
    }
    
    // Fields required for DiagFields
    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
        ostringstream group_name( "" );
        group_name << "FieldsForDiag" << idiag;
        hid_t diag_gid = H5::group( patch_gid, group_name.str() );
        
        for( unsigned int ifield=0; ifield<EMfields->allFields_avg[idiag].size(); ifield++ ) {
            dumpFieldsPerProc( diag_gid, EMfields->allFields_avg[idiag][ifield] );
        }
        
        H5Gclose( diag_gid );
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
                hid_t gid = H5::group( patch_gid, groupName );
                H5::attr( gid, "By_val", embc->By_val );
                H5::attr( gid, "Bz_val", embc->Bz_val );
                H5Gclose( gid );
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                hid_t gid = H5::group( patch_gid, groupName );
                H5::vect( gid, "Bx_val", embc->Bx_val );
                H5::vect( gid, "By_val", embc->By_val );
                H5::vect( gid, "Bz_val", embc->Bz_val );
                H5Gclose( gid );
            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                
                hid_t gid = H5::group( patch_gid, groupName );
                
                if( embc->Bx_val ) {
                    dumpFieldsPerProc( gid, embc->Bx_val );
                }
                if( embc->By_val ) {
                    dumpFieldsPerProc( gid, embc->By_val );
                }
                if( embc->Bz_val ) {
                    dumpFieldsPerProc( gid, embc->Bz_val );
                }
                
                H5Gclose( gid );
            }
        }
    }
    
    H5Fflush( patch_gid, H5F_SCOPE_GLOBAL );
    H5::attr( patch_gid, "species", vecSpecies.size() );
    
    for( unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++ ) {
        ostringstream name( "" );
        name << setfill( '0' ) << setw( 2 ) << ispec;
        string groupName=Tools::merge( "species-", name.str(), "-", vecSpecies[ispec]->name_ );
        hid_t gid = H5::group( patch_gid, groupName );
        
        H5::attr( gid, "partCapacity", vecSpecies[ispec]->particles->capacity() );
        H5::attr( gid, "partSize", vecSpecies[ispec]->particles->size() );
        H5::attr( gid, "nrj_radiation", vecSpecies[ispec]->getNrjRadiation() );
        
        
        if( vecSpecies[ispec]->particles->size()>0 ) {
        
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++ ) {
                ostringstream my_name( "" );
                my_name << "Position-" << i;
                H5::vect( gid, my_name.str(), vecSpecies[ispec]->particles->Position[i], dump_deflate );
            }
            
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++ ) {
                ostringstream my_name( "" );
                my_name << "Momentum-" << i;
                H5::vect( gid, my_name.str(), vecSpecies[ispec]->particles->Momentum[i], dump_deflate );
            }
            
            H5::vect( gid, "Weight", vecSpecies[ispec]->particles->Weight, dump_deflate );
            H5::vect( gid, "Charge", vecSpecies[ispec]->particles->Charge, dump_deflate );
            
            if( vecSpecies[ispec]->particles->tracked ) {
                H5::vect( gid, "Id", vecSpecies[ispec]->particles->Id, H5T_NATIVE_UINT64, dump_deflate );
            }
            
            
            H5::vect( gid, "first_index", vecSpecies[ispec]->particles->first_index );
            H5::vect( gid, "last_index", vecSpecies[ispec]->particles->last_index );
            
        } // End if partSize
        
        H5Gclose( gid );
        
    } // End for ispec
    
    // Manage some collisions parameters
    std::vector<double> rate_multiplier( vecCollisions.size() );
    for( unsigned int icoll = 0; icoll<vecCollisions.size(); icoll++ ) {
        rate_multiplier[icoll] = vecCollisions[icoll]->NuclearReaction->rate_multiplier_;
    }
    H5::vect( patch_gid, "collisions_rate_multiplier", rate_multiplier );
    
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
                        dumpFieldsPerProc( patch_gid, p->magnitude );
                        dumpFieldsPerProc( patch_gid, p->phase );
                        H5::vect( patch_gid, t3.str(), p->omega );
                        H5::vect( patch_gid, t4.str(), p->magnitude->dims_ );
                    }
                }
            }
        }
    }
};


void Checkpoint::readPatchDistribution( SmileiMPI *smpi, SimWindow *simWin )
{
    hid_t fid = H5Fopen( restart_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
    if( fid < 0 ) {
        ERROR( restart_file << " is not a valid HDF5 file" );
    }
    
    // Read basic attributes
    string dump_version;
    H5::getAttr( fid, "Version", dump_version );
    
    string dump_date;
    H5::getAttr( fid, "CommitDate", dump_date );
    
    if( dump_version != string( __VERSION ) ) {
        WARNING( "The code version that dumped the file is " << dump_version );
        WARNING( "                while running version is " << string( __VERSION ) );
    }
    
    vector<int> patch_count( smpi->getSize() );
    H5::getVect( fid, "patch_count", patch_count );
    smpi->patch_count = patch_count;
    
    smpi->patch_refHindexes.resize( smpi->patch_count.size(), 0 );
    smpi->patch_refHindexes[0] = 0;
    for( int rk=1 ; rk<smpi->smilei_sz ; rk++ ) {
        smpi->patch_refHindexes[rk] = smpi->patch_refHindexes[rk-1] + smpi->patch_count[rk-1];
    }
    
    // load window status : required to know the patch movement
    restartMovingWindow( fid, simWin );
    
    H5Fclose( fid );
}


void Checkpoint::restartAll( VectorPatch &vecPatches,  SmileiMPI *smpi, SimWindow *simWin, Params &params, OpenPMDparams &openPMD )
{
    MESSAGE( 1, "READING fields and particles for restart" );
    
    hid_t fid = H5Fopen( restart_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
    if( fid < 0 ) {
        ERROR( restart_file << " is not a valid HDF5 file" );
    }
    
    // Write diags scalar data
    DiagnosticScalar *scalars = static_cast<DiagnosticScalar *>( vecPatches.globalDiags[0] );
    H5::getAttr( fid, "latest_timestep", scalars->latest_timestep );
    // Scalars only by master
    if( smpi->isMaster() ) {
        H5::getAttr( fid, "Energy_time_zero",  scalars->Energy_time_zero );
        H5::getAttr( fid, "EnergyUsedForNorm", scalars->EnergyUsedForNorm );
        // Poynting scalars
        unsigned int k=0;
        for( unsigned int j=0; j<2; j++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
            for( unsigned int i=0; i<params.nDim_field; i++ ) { //axis 0=x, 1=y, 2=z
                string poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
                if( H5Aexists( fid, poy_name.c_str() )>0 ) {
                    H5::getAttr( fid, poy_name, vecPatches( 0 )->EMfields->poynting[j][i] );
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
                int vect_size = H5::getVectSize( fid, diagName.str() );
                if( vect_size == target_size ) {
                    H5::getVect( fid, diagName.str(), *(screen->getData()) );
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
        string patchName=Tools::merge( "patch-", patch_name.str() );
        hid_t patch_gid = H5Gopen( fid, patchName.c_str(), H5P_DEFAULT );
        
        restartPatch( vecPatches( ipatch )->EMfields, vecPatches( ipatch )->vecSpecies, vecPatches( ipatch )->vecCollisions, params, patch_gid );
        
        // Random number generator state
        H5::getAttr( patch_gid, "xorshift32_state", vecPatches( ipatch )->rand_->xorshift32_state );
        
        H5Gclose( patch_gid );
        
    }
    
    // Read the latest Id that the MPI processes have given to each species
    for( unsigned int idiag=0; idiag<vecPatches.localDiags.size(); idiag++ ) {
        if( DiagnosticTrack *track = dynamic_cast<DiagnosticTrack *>( vecPatches.localDiags[idiag] ) ) {
            ostringstream n( "" );
            n<< "latest_ID_" << vecPatches( 0 )->vecSpecies[track->speciesId_]->name_;
            if( H5::hasAttr( fid, n.str() ) ) {
                H5::getAttr( fid, n.str(), track->latest_Id, H5T_NATIVE_UINT64 );
            } else {
                track->IDs_done=false;
            }
        }
    }
    
    H5Fclose( fid );
    
}


void Checkpoint::restartPatch( ElectroMagn *EMfields, std::vector<Species *> &vecSpecies, std::vector<Collisions *> &vecCollisions, Params &params, hid_t patch_gid )
{
    if ( params.geometry != "AMcylindrical" ) {
        restartFieldsPerProc( patch_gid, EMfields->Ex_ );
        restartFieldsPerProc( patch_gid, EMfields->Ey_ );
        restartFieldsPerProc( patch_gid, EMfields->Ez_ );
        restartFieldsPerProc( patch_gid, EMfields->Bx_ );
        restartFieldsPerProc( patch_gid, EMfields->By_ );
        restartFieldsPerProc( patch_gid, EMfields->Bz_ );
        restartFieldsPerProc( patch_gid, EMfields->Bx_m );
        restartFieldsPerProc( patch_gid, EMfields->By_m );
        restartFieldsPerProc( patch_gid, EMfields->Bz_m );
    }
    else {
        for ( unsigned int imode = 0 ; imode < params.nmodes ; imode++ ) {
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            restart_cFieldsPerProc( patch_gid, emAM->El_[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Er_[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Et_[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Bl_[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Br_[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Bt_[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Bl_m[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Br_m[imode] );
            restart_cFieldsPerProc( patch_gid, emAM->Bt_m[imode] );
            
            if(params.is_pxr == true)
                restart_cFieldsPerProc( patch_gid, emAM->rho_old_AM_[imode] );
            
        }
    }
    

    if( EMfields->envelope!=NULL ) {
        DEBUG( "restarting envelope" );
        restart_cFieldsPerProc( patch_gid, EMfields->envelope->A_ );
        restart_cFieldsPerProc( patch_gid, EMfields->envelope->A0_ );
        restartFieldsPerProc( patch_gid, EMfields->Env_Chi_ );
    } else {
        DEBUG( "envelope is null" );
    }
    
    
    // filtered Electric fields
    for( unsigned int i=0; i<EMfields->Exfilter.size(); i++ ) {
        restartFieldsPerProc( patch_gid, EMfields->Exfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Eyfilter.size(); i++ ) {
        restartFieldsPerProc( patch_gid, EMfields->Eyfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Ezfilter.size(); i++ ) {
        restartFieldsPerProc( patch_gid, EMfields->Ezfilter[i] );
    }
    // filtered Magnetic fields
    for( unsigned int i=0; i<EMfields->Bxfilter.size(); i++ ) {
        restartFieldsPerProc( patch_gid, EMfields->Bxfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Byfilter.size(); i++ ) {
        restartFieldsPerProc( patch_gid, EMfields->Byfilter[i] );
    }
    for( unsigned int i=0; i<EMfields->Bzfilter.size(); i++ ) {
        restartFieldsPerProc( patch_gid, EMfields->Bzfilter[i] );
    }
    
    // Fields required for DiagFields
    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
        ostringstream group_name( "" );
        group_name << "FieldsForDiag" << idiag;
        htri_t status = H5Lexists( patch_gid, group_name.str().c_str(), H5P_DEFAULT );
        if( status > 0 ) {
            hid_t diag_gid = H5Gopen( patch_gid, group_name.str().c_str(), H5P_DEFAULT );
            
            for( unsigned int ifield=0; ifield<EMfields->allFields_avg[idiag].size(); ifield++ ) {
                restartFieldsPerProc( diag_gid, EMfields->allFields_avg[idiag][ifield] );
            }
            
            H5Gclose( diag_gid );
            
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
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                hid_t gid = H5Gopen( patch_gid, groupName.c_str(), H5P_DEFAULT );
                H5::getAttr( gid, "By_val", embc->By_val );
                H5::getAttr( gid, "Bz_val", embc->Bz_val );
                H5Gclose( gid );
                
            } else if( dynamic_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC2D_SM *embc = static_cast<ElectroMagnBC2D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                hid_t gid = H5Gopen( patch_gid, groupName.c_str(), H5P_DEFAULT );
                H5::getVect( gid, "Bx_val", embc->Bx_val );
                H5::getVect( gid, "By_val", embc->By_val );
                H5::getVect( gid, "Bz_val", embc->Bz_val );
                H5Gclose( gid );
            } else if( dynamic_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] ) ) {
                ElectroMagnBC3D_SM *embc = static_cast<ElectroMagnBC3D_SM *>( EMfields->emBoundCond[bcId] );
                ostringstream name( "" );
                name << setfill( '0' ) << setw( 2 ) << bcId;
                string groupName=Tools::merge( "EM_boundary-species-", name.str() );
                hid_t gid = H5Gopen( patch_gid, groupName.c_str(), H5P_DEFAULT );
                
                if( embc->Bx_val ) {
                    restartFieldsPerProc( gid, embc->Bx_val );
                }
                if( embc->By_val ) {
                    restartFieldsPerProc( gid, embc->By_val );
                }
                if( embc->Bz_val ) {
                    restartFieldsPerProc( gid, embc->Bz_val );
                }
                H5Gclose( gid );
            }
        }
    }
    
    unsigned int vecSpeciesSize=0;
    H5::getAttr( patch_gid, "species", vecSpeciesSize );
    
    if( vecSpeciesSize != vecSpecies.size() ) {
        ERROR( "Number of species differs between dump (" << vecSpeciesSize << ") and namelist ("<<vecSpecies.size()<<")" );
    }
    
    
    for( unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++ ) {
        ostringstream name( "" );
        name << setfill( '0' ) << setw( 2 ) << ispec;
        string groupName=Tools::merge( "species-", name.str(), "-", vecSpecies[ispec]->name_ );
        hid_t gid = H5Gopen( patch_gid, groupName.c_str(), H5P_DEFAULT );
        
        unsigned int partCapacity=0;
        H5::getAttr( gid, "partCapacity", partCapacity );
        vecSpecies[ispec]->particles->reserve( partCapacity, nDim_particle );
        
        unsigned int partSize=0;
        H5::getAttr( gid, "partSize", partSize );
        vecSpecies[ispec]->particles->initialize( partSize, nDim_particle );
        
        double nrj_radiation;
        if( H5::hasAttr( gid, "nrj_radiation" ) ) {
            H5::getAttr( gid, "nrj_radiation", nrj_radiation );
            vecSpecies[ispec]->setNrjRadiation( nrj_radiation );
        }
        
        if( partSize>0 ) {
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++ ) {
                ostringstream namePos( "" );
                namePos << "Position-" << i;
                H5::getVect( gid, namePos.str(), vecSpecies[ispec]->particles->Position[i] );
            }
            
            for( unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++ ) {
                ostringstream namePos( "" );
                namePos << "Momentum-" << i;
                H5::getVect( gid, namePos.str(), vecSpecies[ispec]->particles->Momentum[i] );
            }
            
            H5::getVect( gid, "Weight", vecSpecies[ispec]->particles->Weight );
            
            H5::getVect( gid, "Charge", vecSpecies[ispec]->particles->Charge );
            
            if( vecSpecies[ispec]->particles->tracked ) {
                H5::getVect( gid, "Id", vecSpecies[ispec]->particles->Id, H5T_NATIVE_UINT64 );
            }
            
            if( params.vectorization_mode == "off" || params.vectorization_mode == "on" || params.cell_sorting ) {
                H5::getVect( gid, "first_index", vecSpecies[ispec]->particles->first_index, true );
                H5::getVect( gid, "last_index", vecSpecies[ispec]->particles->last_index, true );
            }
            // In the adaptive vectorization case, the bins will be recomputed
            // latter in the patch reconfiguration
            
        }
        
        H5Gclose( gid );
    }
    
    // Manage some collisions parameters
    if( H5::getVectSize( patch_gid, "collisions_rate_multiplier" ) > 0 ) {
        std::vector<double> rate_multiplier;
        H5::getVect( patch_gid, "collisions_rate_multiplier", rate_multiplier, true );
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
                        H5::getVect( patch_gid, t4.str(), dims, true );
                        p->magnitude->allocateDims( dims );
                        p->phase->allocateDims( dims );
                        restartFieldsPerProc( patch_gid, p->magnitude );
                        restartFieldsPerProc( patch_gid, p->phase );
                        H5::getVect( patch_gid, t3.str(), p->omega, true );
                    }
                }
            }
        }
    }
}

void Checkpoint::dumpFieldsPerProc( hid_t fid, Field *field )
{
    hsize_t dims[1]= {field->globalDims_};
    hid_t sid = H5Screate_simple( 1, dims, NULL );
    hid_t did = H5Dcreate( fid, field->name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &field->data_[0] );
    H5Dclose( did );
    H5Sclose( sid );
}

void Checkpoint::dump_cFieldsPerProc( hid_t fid, Field *field )
{
    cField *cfield = static_cast<cField *>( field );
    hsize_t dims[1]= {2*field->globalDims_}; //*2 : to manage complex data
    hid_t sid = H5Screate_simple( 1, dims, NULL );
    hid_t did = H5Dcreate( fid, field->name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    H5Dwrite( did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cfield->cdata_[0] );
    H5Dclose( did );
    H5Sclose( sid );
}

void Checkpoint::restartFieldsPerProc( hid_t fid, Field *field )
{
    hid_t did = H5Dopen( fid, field->name.c_str(), H5P_DEFAULT );
    H5Dread( did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &field->data_[0] );
    H5Dclose( did );
}

void Checkpoint::restart_cFieldsPerProc( hid_t fid, Field *field )
{
    cField *cfield = static_cast<cField *>( field );
    hid_t did = H5Dopen( fid, field->name.c_str(), H5P_DEFAULT );
    H5Dread( did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cfield->cdata_[0] );
    H5Dclose( did );
}

void Checkpoint::dumpMovingWindow( hid_t fid, SimWindow *simWin )
{
    H5::attr( fid, "x_moved", simWin->getXmoved() );
    H5::attr( fid, "n_moved", simWin->getNmoved() );
    
}
void Checkpoint::restartMovingWindow( hid_t fid, SimWindow *simWin )
{

    double x_moved=0.;
    H5::getAttr( fid, "x_moved", x_moved );
    simWin->setXmoved( x_moved );
    
    unsigned int n_moved=0;
    H5::getAttr( fid, "n_moved", n_moved );
    simWin->setNmoved( n_moved );
    
}
