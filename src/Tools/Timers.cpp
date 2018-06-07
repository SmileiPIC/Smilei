#include <mpi.h>
#include <iomanip>

#include "Timers.h"

#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

Timers::Timers( SmileiMPI * smpi ) :
    global    ("Global"        ), // The entire time loop
    particles ("Particles"     ), // Call dynamics + restartRhoJ(s)
    maxwell   ("Maxwell"       ), // Maxwell
    diags     ("Diagnostics"   ), // Diags.runAllDiags + MPI & Patch sync
    densities ("Densities"     ), // Local sum of rho, Jxyz
    collisions("Collisions"    ), // Call to Collisions methods
    movWindow ("Mov window"    ), // Moving Window
    loadBal   ("Load balancing" ), // Load balancing
    syncPart  ("Sync Particles"), // Call exchangeParticles (MPI & Patch sync)
    syncField ("Sync Fields"   ), // Call sumRhoJ(s), exchangeB (MPI & Patch sync)
    syncDens  ("Sync Densities"),  // If necessary the following timers can be reintroduced
    diagsNEW  ("DiagnosticsNEW" ), // Diags.runAllDiags + MPI & Patch sync
    reconfiguration("Reconfiguration")
#ifdef __DETAILED_TIMERS
    ,interpolator("Interpolator"),
    pusher("Pusher"             ),
    projector("Projector"       ),
    particles_boundaries("Particles boundaries"),
    ionization("Ionization"       ),
    radiation("Radiation"       ),
    multiphoton_Breit_Wheeler_timer("Multiphoton Breit-Wheeler"       )
#endif
{
    timers.resize(0);
    timers.push_back( &global     );
    timers.push_back( &particles  );
    timers.push_back( &maxwell    );
    timers.push_back( &diags      );
    timers.push_back( &densities  );
    timers.push_back( &collisions );
    timers.push_back( &movWindow  );
    timers.push_back( &loadBal    );
    timers.push_back( &syncPart   );
    timers.push_back( &syncField  );
    timers.push_back( &syncDens   );
    timers.push_back( &diagsNEW   );
    timers.push_back( &reconfiguration   );
#ifdef __DETAILED_TIMERS
    patch_timer_id_start = timers.size()-1;
    timers.push_back( &interpolator   );
    timers.back()->patch_timer_id = 0;
    timers.push_back( &pusher   );
    timers.back()->patch_timer_id = 1;
    timers.push_back( &projector   );
    timers.back()->patch_timer_id = 2;
    timers.push_back( &particles_boundaries   );
    timers.back()->patch_timer_id = 3;
    timers.push_back( &ionization   );
    timers.back()->patch_timer_id = 4;
    timers.push_back( &radiation   );
    timers.back()->patch_timer_id = 5;
    timers.push_back( &multiphoton_Breit_Wheeler_timer   );
    timers.back()->patch_timer_id = 6;
#endif

    for( unsigned int i=0; i<timers.size(); i++)
        timers[i]->init(smpi);

    if (smpi->getRank()==0) {
        remove ("profil.txt");
        ofstream fout;
        fout.open ("profil.txt");
        fout.close();
    }

}

Timers::~Timers()
{
}

void Timers::reboot()
{
    for( unsigned int i=0; i<timers.size(); i++)
        timers[i]->reboot();
}

//! Output the timer profile
void Timers::profile(SmileiMPI * smpi)
{
    std::vector<Timer*> avg_timers = consolidate(smpi);

    if ( smpi->isMaster() ) {
        double coverage(0.);
        for (unsigned int i=1 ; i<timers.size() ; i++)
            coverage += timers[i]->getTime();

        MESSAGE("Time in time loop :\t" << global.getTime() << "\t"<<coverage/global.getTime()*100.<< "% coverage" );

#ifdef __DETAILED_TIMERS

        for (unsigned int i=0 ; i<patch_timer_id_start ; i++)
            avg_timers[i]->print(global.getTime());

        MESSAGE(0, "\n\t Patch timers:" );

        for (unsigned int i=patch_timer_id_start ; i<avg_timers.size() ; i++)
            avg_timers[i]->print(global.getTime());

        MESSAGE(0, "\n\t Printed times are averaged per MPI process" );
        MESSAGE(0, "\t\t See advanced metrics in profil.txt");

#else

        for (unsigned int i=0 ; i<avg_timers.size() ; i++)
            avg_timers[i]->print(global.getTime());
        MESSAGE(0, "\n\t Printed times are averaged per MPI process" );
        MESSAGE(0, "\t\t See advanced metrics in profil.txt");

#endif

    }
    for (unsigned int i=0 ; i<avg_timers.size() ; i++)
        delete avg_timers[i];
}

//! Perform the required processing on the timers for output
std::vector<Timer*> Timers::consolidate(SmileiMPI * smpi)
{
    std::vector<Timer*> avg_timers;

    int sz = smpi->getSize(), rk = smpi->getRank();

    ofstream fout;
    if (rk==0) fout.open ("profil.txt", ofstream::out | ofstream::app );
    fout << endl << endl << "--- Timestep = " << (timers[1]->register_timers.size()-1) << " x Main.print_every = " <<  " ---" << endl;

    // timers[0] is the global PIC loop timer, naturally synchronized
    for ( unsigned int itimer = 1 ; itimer < timers.size() ; itimer++ ) {
        int nrecords(0);
        if (timers[itimer]->register_timers.size()>0)
            nrecords = timers[itimer]->register_timers.size();
        else
            nrecords = 1;
        double* tmp = new double[sz*nrecords];

        if (timers[itimer]->register_timers.size()>0)
            MPI_Gather( &(timers[itimer]->register_timers[0]), nrecords, MPI_DOUBLE,
                        &(tmp[0]), nrecords, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
        else
            MPI_Gather( &(timers[itimer]->time_acc_), 1, MPI_DOUBLE,
                        &(tmp[0]), 1, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);


        // Mean on last records
        int idx_records=nrecords-1; // = last record

        double min( tmp[idx_records] ), max( tmp[idx_records] ), avg( tmp[idx_records] );
        double sum( tmp[idx_records] ), sig( tmp[idx_records]*tmp[idx_records]/(double)sz );
        if (rk==0)
            for ( int i=1 ; i<sz ; i++) {
                if ( tmp[idx_records+i*(nrecords)] < min ) min = tmp[idx_records+i*(nrecords)];
                if ( tmp[idx_records+i*(nrecords)] > max ) max = tmp[idx_records+i*(nrecords)];
                sum += tmp[idx_records+i*(nrecords)];
                sig += tmp[idx_records+i*(nrecords)]*tmp[idx_records+i*(nrecords)]/(double)sz;
            }
        avg = sum / sz;

        delete [] tmp;

        if ((max>0.) && (rk==0)) {
            fout.setf( ios::fixed,  ios::floatfield );
            fout << setw(14) << scientific << setprecision(3)
                 << timers[itimer]->name_ << "\t : " << "Min time =  " << min
                 << "\t - \t" <<  "Avg time =  " << avg
                 << "\t - \t" <<  "Max time =  " << max
                 << "\t - \t" <<  "SD time =  " << sqrt( sig-sum*sum/(double)(sz)/(double)(sz) )
                 << endl;
        }
        if (rk==0) {
            Timer * newTimer = new Timer("");
            newTimer->time_acc_ = avg;
            newTimer->name_ = timers[itimer]->name_;
            avg_timers.push_back( newTimer );
        }

    }
    if (rk==0) fout.close();

    return avg_timers;

}
