#include "Timer.h"

#include <iomanip>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"
#include "Tools.h"
#include "VectorPatch.h"

using namespace std;

Timer::Timer( string name ) :
    name_( name ),
    time_acc_( 0.0 ),
    smpi_( NULL )
{
    register_timers.resize( 0, 0. );
}

Timer::~Timer()
{
}

void Timer::init( SmileiMPI *smpi )
{
    smpi_ = smpi;
    smpi_->barrier();
    last_start_ = MPI_Wtime();
}

//! Accumulate time couting from last init/restart
void Timer::update( bool store )
{
    #pragma omp barrier
    #pragma omp master
    {
        time_acc_ +=  MPI_Wtime()-last_start_;
        last_start_ = MPI_Wtime();
        if( store )
        {
            register_timers.push_back( time_acc_ );
        }
    }
}

#ifdef __DETAILED_TIMERS
//!Accumulate time couting from last init/restart using patch detailed timers
void Timer::update( VectorPatch &vecPatches, bool store )
{
    #pragma omp barrier
    #pragma omp master
    {
        // Reduce the time spent in all patches in time_tmp
        double time_tmp = 0.;
        for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ )
        {
            time_tmp += vecPatches( ipatch )->patch_timers[this->patch_timer_id];
            vecPatches( ipatch )->patch_timers[this->patch_timer_id] = 0;
        }
        
        // Get the number of threads per MPI in order to evaluate the mean per patch
        int thread_number = 0.;
#ifdef _OPENMP
        thread_number = omp_get_num_threads();
#else
        thread_number = 1;
#endif
        
        // Average over all patches
        this->time_acc_  += time_tmp / ( double )( thread_number );
        
        if( store )
        {
            register_timers.push_back( time_acc_ );
        }
    }
}
#endif

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

void Timer::print( double tot )
{
    if( ( time_acc_>0. ) && ( name_!="" ) ) {
        //cout << "\t" << setw(12) << name_ << "\t" << time_acc_  << "\t(" << 100.0*time_acc_/tot << "%)" << endl ;
        //MESSAGE(0, "\t" << setw(12) << name_ << "\t" << time_acc_  << "\t(" << 100.0*time_acc_/tot << "%)");
        double perc=100.0*time_acc_/tot;
        if( perc<0.05 ) {
            return;
        }
        if( perc<1 ) {
            MESSAGE( 0, "\t" << setw( 20 ) << name_ << "\t" << setprecision(6) << time_acc_ << setprecision(1) << "\t" << "    <1%" );
        } else {
            MESSAGE( 0, "\t" << setw( 20 ) << name_ << "\t" << setprecision(6) << time_acc_ << setprecision(1) << "\t" << perc << "%" );
        }
    }
}
