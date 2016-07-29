#include "Timer.h"

#include <iomanip>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

Timer::Timer() :
time_acc_(0.0),
name_(""),
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

void Timer::consolidate_timers( std::vector<Timer> timers )
{
    int sz( 1 ), rk( 0 );
    MPI_Comm_size( MPI_COMM_WORLD, &sz );
    MPI_Comm_rank( MPI_COMM_WORLD, &rk );

    ofstream fout;
    if (rk==0) fout.open ("profil.txt");

    for ( unsigned int itimer = 0 ; itimer < timers.size() ; itimer++ ) {
        int nrecords( timers[itimer].register_timers.size() );
        double* tmp = new double[sz*nrecords];

        MPI_Gather( &(timers[itimer].register_timers[0]), nrecords, MPI_DOUBLE,
                    &(tmp[0]), nrecords, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);


        // Mean on last records
        int idx_records=nrecords-1; // = last record

        double min( tmp[idx_records] ), max( tmp[idx_records] ), avg( tmp[idx_records] );
        if (rk==0)
            for ( int i=1 ; i<sz ; i++) {
                if ( tmp[idx_records+i*(nrecords)] < min ) min = tmp[idx_records+i*(nrecords)];
                if ( tmp[idx_records+i*(nrecords)] > max ) max = tmp[idx_records+i*(nrecords)];
                avg += tmp[idx_records+i*(nrecords)];
            }
        avg /= sz;

        delete [] tmp;

        if ((max>0.) && (rk==0)) {   
            fout.setf( ios::fixed,  ios::floatfield );
            fout << setw(14) << scientific << setprecision(3)
                 << timers[itimer].name_ << "\t : " << "Min time =  " << min
                 << "\t - \t" <<  "Avg time =  " << avg
                 << "\t - \t" <<  "Max time =  " << max
                 << endl;
        }
        
    }
    if (rk==0) fout.close();

}
