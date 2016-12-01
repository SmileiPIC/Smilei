#include "Timer.h"

#include <iomanip>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

Timer::Timer() :
name_(""),
time_acc_(0.0),
smpi_(NULL)
{
    register_timers.resize(0,0.);    
}

Timer::~Timer()
{
}

void Timer::init(SmileiMPI *smpi, string name)
{
    smpi_ = smpi;
    smpi_->barrier();
    last_start_ = MPI_Wtime();
    name_ = name;
}


void Timer::update(bool store)
{
    #pragma omp barrier
    #pragma omp master
    {
        time_acc_ +=  MPI_Wtime()-last_start_;
        last_start_ = MPI_Wtime();
        if (store) register_timers.push_back( time_acc_ );
    }
}

void Timer::restart()
{
    #pragma omp barrier
    #pragma omp master
    {
        last_start_ = MPI_Wtime();
    }
}

void Timer::reboot()
{
    smpi_->barrier();
    last_start_ =  MPI_Wtime();
    time_acc_ = 0.;
}

void Timer::print(double tot)
{
    if ((time_acc_>0.) && (name_!="")) {
      //cout << "\t" << setw(12) << name_ << "\t" << time_acc_  << "\t(" << 100.0*time_acc_/tot << "%)" << endl ;
        //MESSAGE(0, "\t" << setw(12) << name_ << "\t" << time_acc_  << "\t(" << 100.0*time_acc_/tot << "%)");
        double perc=100.0*time_acc_/tot;
        if (perc<0.001) return;
        if (perc<1) {
            MESSAGE(0, "\t" << setw(14) << name_ << "\t" << time_acc_  << "\t" << "<1%");
        } else {
            MESSAGE(0, "\t" << setw(14) << name_ << "\t" << time_acc_  << "\t" << perc << "%");
        }
    }
}

std::vector<Timer> Timer::consolidate_timers( std::vector<Timer> timers )
{
    std::vector<Timer> avg_timers(timers.size());

    int sz( 1 ), rk( 0 );
    MPI_Comm_size( MPI_COMM_WORLD, &sz );
    MPI_Comm_rank( MPI_COMM_WORLD, &rk );

    ofstream fout;
    if (rk==0) fout.open ("profil.txt");

    // timers[0] is the global PIC loop timer, naturally synchronized
    for ( unsigned int itimer = 1 ; itimer < timers.size() ; itimer++ ) {
        int nrecords(0);
        if (timers[itimer].register_timers.size()>0)
            nrecords = timers[itimer].register_timers.size();
        else
            nrecords = 1;
        double* tmp = new double[sz*nrecords];

        if (timers[itimer].register_timers.size()>0)
            MPI_Gather( &(timers[itimer].register_timers[0]), nrecords, MPI_DOUBLE,
                        &(tmp[0]), nrecords, MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
        else
            MPI_Gather( &(timers[itimer].time_acc_), 1, MPI_DOUBLE,
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
                 << timers[itimer].name_ << "\t : " << "Min time =  " << min
                 << "\t - \t" <<  "Avg time =  " << avg
                 << "\t - \t" <<  "Max time =  " << max
                 << "\t - \t" <<  "SD time =  " << sqrt( sig-sum*sum/(double)(sz)/(double)(sz) )
                 << endl;
        }
        if (rk==0) {
            avg_timers[itimer].time_acc_ = avg;
            avg_timers[itimer].name_ = timers[itimer].name_;
        }

        
    }
    if (rk==0) fout.close();

    return avg_timers;

}

std::vector<Timer> Timer::initialize_timers( SmileiMPI *smpi ) {
    // GC IDRIS : "Timer timer[ntimer];" to "Timer timer[8];"
    int ntimer(13);
    vector<Timer> timer(ntimer);
    // The entire time loop
    timer[0].init(smpi, "Global");
    // Call dynamics + restartRhoJ(s)
    timer[1].init(smpi, "Particles");
    // Maxwell
    timer[2].init(smpi, "Maxwell");
    // Diags.runAllDiags + MPI & Patch sync
    timer[3].init(smpi, "Diagnostics");
    // Local sum of rho, Jxyz
    timer[4].init(smpi, "Densities");
    // Moving Window
    timer[5].init(smpi, "Mov window");
    // Dump fields (including average)
    timer[6].init(smpi, "Diag fields");
    // Load balancing
    timer[7].init(smpi, "Load balacing");
    // Call exchangeParticles (MPI & Patch sync)
    timer[8].init(smpi, "Sync Particles");
    // Call sumRhoJ(s), exchangeB (MPI & Patch sync)
    timer[9].init(smpi, "Sync Fields");
    // Call to Collisions methods
    timer[10].init(smpi, "Collisions");
    // If necessary the following timers can be reintroduced
    timer[11].init(smpi, "Sync Densities");
    //timer[12].init(smpi, "AvgFields");
    
    return timer;
}
