#include <mpi.h>
#include <iomanip>

#include "Timers.h"

#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

Timers::Timers( SmileiMPI *smpi ) :
    global( "Global" ),           // The entire time loop
    particles( "Particles" ),     // Call dynamics + restartRhoJ(s)
    maxwell( "Maxwell" ),         // Maxwell
    diags( "Diagnostics" ),       // Diags.runAllDiags + MPI & Patch sync
    densities( "Densities" ),     // Local sum of rho, Jxyz
    collisions( "Collisions" ),             // Call to Collisions methods
    movWindow( "Mov window" ),              // Moving Window
    loadBal( "Load balancing" ),            // Load balancing
    syncPart( "Sync Particles" ),           // Call exchangeParticles (MPI & Patch sync)
    syncField( "Sync Fields" ),             // Call sumRhoJ(s), exchangeB (MPI & Patch sync)
    syncDens( "Sync Densities" ),           // If necessary the following timers can be reintroduced
    particleMerging( "Part Merging" ),      // Particle merging
    particleInjection( "Part Injection" ),  // Particle injection
    diagsNEW( "DiagnosticsNEW" ),           // Diags.runAllDiags + MPI & Patch sync
    reconfiguration( "Reconfiguration" ),   // Patch reconfiguration
    envelope( "Envelope" ),
    susceptibility( "Sync_Susceptibility" ),
    grids("Grids")
#ifdef __DETAILED_TIMERS
    // Details of Dynamic
    , interpolator( "Interpolator" ),
    pusher( "Pusher" ),
    projector( "Projector" ),
    cell_keys( "Cell_keys" ),
    ionization( "Ionization" ),
    radiation( "Radiation" ),
    multiphoton_Breit_Wheeler_timer( "Multiphoton_Breit-Wheeler" ),
    // Details of Envelop
    interp_fields_env( "Interp_Fields_Env" ),
    proj_susceptibility( "Proj_Susceptibility" ),
    push_mom( "Push_Momentum" ),
    interp_env_old( "Interp_Env_Old" ),
    proj_currents( "Proj_Currents" ),
    push_pos( "Push_Pos" ),
    // Details of Sync Particles
    sorting( "Sorting" )
#endif
{
    timers.resize( 0 );
    timers.push_back( &global );
    timers.push_back( &particles );
    timers.push_back( &maxwell );
    timers.push_back( &diags );
    timers.push_back( &densities );
    timers.push_back( &collisions );
    timers.push_back( &movWindow );
    timers.push_back( &loadBal );
    timers.push_back( &syncPart );
    timers.push_back( &syncField );
    timers.push_back( &syncDens );
    timers.push_back( &particleMerging );
    timers.push_back( &particleInjection );
    timers.push_back( &diagsNEW );
    timers.push_back( &reconfiguration );
    timers.push_back( &envelope );
    timers.push_back( &susceptibility );
    timers.push_back( &grids );
    patch_timer_id_start = timers.size()-1;
#ifdef __DETAILED_TIMERS
    timers.push_back( &interpolator );
    timers.back()->patch_timer_id = 0;
    timers.push_back( &pusher );
    timers.back()->patch_timer_id = 1;
    timers.push_back( &projector );
    timers.back()->patch_timer_id = 2;
    timers.push_back( &cell_keys );
    timers.back()->patch_timer_id = 3;
    timers.push_back( &ionization );
    timers.back()->patch_timer_id = 4;
    timers.push_back( &radiation );
    timers.back()->patch_timer_id = 5;
    timers.push_back( &multiphoton_Breit_Wheeler_timer );
    timers.back()->patch_timer_id = 6;
    
    timers.push_back( &interp_fields_env );
    timers.back()->patch_timer_id = 7;
    timers.push_back( &proj_susceptibility );
    timers.back()->patch_timer_id = 8;
    timers.push_back( &push_mom );
    timers.back()->patch_timer_id = 9;
    
    timers.push_back( &interp_env_old );
    timers.back()->patch_timer_id = 10;
    timers.push_back( &push_pos );
    timers.back()->patch_timer_id = 11;
    timers.push_back( &proj_currents ) ;
    timers.back()->patch_timer_id = 12;
    
    // Details of Sync Particles
    timers.push_back( &sorting ) ;
    timers.back()->patch_timer_id = 13;
#endif
    
    for( unsigned int i=0; i<timers.size(); i++ ) {
        timers[i]->init( smpi );
    }
    
    if( smpi->getRank()==0 && ! smpi->test_mode ) {
        remove( "profil.txt" );
        ofstream fout;
        fout.open( "profil.txt" );
        fout.close();
    }
    
}

Timers::~Timers()
{
}

void Timers::reboot()
{
    for( unsigned int i=0; i<timers.size(); i++ ) {
        timers[i]->reboot();
    }
}

//! Output the timer profile
void Timers::profile( SmileiMPI *smpi )
{
    std::vector<Timer *> avg_timers = consolidate( smpi, true );
    
    if( smpi->isMaster() ) {
        double coverage( 0. );
        // Computation of the coverage: it only takes into account
        // the main timers (14)
        for( unsigned int i=1 ; i<patch_timer_id_start+1 ; i++ ) {
            coverage += timers[i]->getTime();
        }
        
        MESSAGE( "Time_in_time_loop\t" << global.getTime() << "\t"<<coverage/global.getTime()*100.<< "% coverage" );
        
#ifdef __DETAILED_TIMERS
        
        for( unsigned int i=0 ; i<patch_timer_id_start ; i++ ) {
            avg_timers[i]->print( global.getTime() );
        }
        
        MESSAGE( "\n Patch average timers:" );
        
        for( unsigned int i=patch_timer_id_start ; i<avg_timers.size() ; i++ ) {
            avg_timers[i]->print( global.getTime() );
        }
        
        MESSAGE( 0, "\n\t Printed times are averaged per MPI process" );
        MESSAGE( 0, "\t\t See advanced metrics in profil.txt" );
        
#else
        
        for( unsigned int i=0 ; i<avg_timers.size() ; i++ ) {
            avg_timers[i]->print( global.getTime() );
        }
        MESSAGE( 0, "\n\t Printed times are averaged per MPI process" );
        MESSAGE( 0, "\t\t See advanced metrics in profil.txt" );
        
#endif
        
    }
    for( unsigned int i=0 ; i<avg_timers.size() ; i++ ) {
        delete avg_timers[i];
    }
}

//! Perform the required processing on the timers for output
std::vector<Timer *> Timers::consolidate( SmileiMPI *smpi, bool final_profile )
{
    std::vector<Timer *> avg_timers;
    int sz = smpi->getSize(), rk = smpi->getRank();
    
    ofstream fout;
    if( rk==0 && ! smpi->test_mode ) {
        fout.open( "profil.txt", ofstream::out | ofstream::app );
        fout << endl << endl << "--- Timestep = " << ( timers[1]->register_timers.size()-1 ) << " x Main.print_every = " <<  " ---" << endl;
            fout << setw(14) << scientific << setprecision(3)
                 << "Time \t " << "Min   "
                 << "\t\t " << "Avg  "
                 << "\t\t " << "Max   "
                 << "\t\t " << "SD "
                 << endl;
    }
    
    // timers[0] is the global PIC loop timer, naturally synchronized
    for( unsigned int itimer = 1 ; itimer < timers.size() ; itimer++ ) {
        int nrecords( 0 );
        if( timers[itimer]->register_timers.size()>0 ) {
            nrecords = timers[itimer]->register_timers.size();
        } else {
            nrecords = 1;
        }
        double *tmp = new double[sz*nrecords];
        
        if( timers[itimer]->register_timers.size()>0 )
            MPI_Gather( &( timers[itimer]->register_timers[0] ), nrecords, MPI_DOUBLE,
                        &( tmp[0] ), nrecords, MPI_DOUBLE,
                        0, MPI_COMM_WORLD );
        else
            MPI_Gather( &( timers[itimer]->time_acc_ ), 1, MPI_DOUBLE,
                        &( tmp[0] ), 1, MPI_DOUBLE,
                        0, MPI_COMM_WORLD );
                        
                        
        // Mean on last records
        int idx_records=nrecords-1; // = last record
        
        double min( tmp[idx_records] ), max( tmp[idx_records] ), avg( tmp[idx_records] );
        double sum( tmp[idx_records] ), sig( tmp[idx_records]*tmp[idx_records]/( double )sz );
        if( rk==0 )
            for( int i=1 ; i<sz ; i++ ) {
                if( tmp[idx_records+i*( nrecords )] < min ) {
                    min = tmp[idx_records+i*( nrecords )];
                }
                if( tmp[idx_records+i*( nrecords )] > max ) {
                    max = tmp[idx_records+i*( nrecords )];
                }
                sum += tmp[idx_records+i*( nrecords )];
                sig += tmp[idx_records+i*( nrecords )]*tmp[idx_records+i*( nrecords )]/( double )sz;
            }
        avg = sum / sz;
        
        delete [] tmp;
        
        if( ( max>0. ) && ( rk==0 ) && ! smpi->test_mode ) {
            fout.setf( ios::fixed,  ios::floatfield );
            if (avg/timers[0]->time_acc_>0.001) {
                fout << setw(14) << scientific << setprecision(3)
                     << timers[itimer]->name_ 
                     << "\t " << min
                     << "\t " << avg
                     << "\t " << max
                     << "\t " << sqrt( sig-sum*sum/(double)(sz)/(double)(sz) )
                     << endl;
            }
        }
        if( ( rk==0 ) && ( final_profile ) ) {
            Timer *newTimer = new Timer( "" );
            newTimer->time_acc_ = avg;
            newTimer->name_ = timers[itimer]->name_;
            avg_timers.push_back( newTimer );
        }
        
    }
    if( rk==0 && ! smpi->test_mode ) {
        fout.close();
    }
    
    return avg_timers;
    
}
