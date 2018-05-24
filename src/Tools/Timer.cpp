#include "Timer.h"

#include <iomanip>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

Timer::Timer( string name ) :
name_(name),
time_acc_(0.0),
smpi_(NULL)
{
    register_timers.resize(0,0.);    
}

Timer::~Timer()
{
}

void Timer::init(SmileiMPI *smpi)
{
    smpi_ = smpi;
    smpi_->barrier();
    last_start_ = MPI_Wtime();
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
    register_timers.clear();
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
