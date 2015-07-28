#include "Timer.h"

#include <iomanip>
#include <string>

#include <mpi.h>

#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

Timer::Timer()
{
    smpi_ = NULL;
    name_ = "";
    time_acc_=0.0;
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


void Timer::update()
{
    smpi_->barrier();
    time_acc_ +=  MPI_Wtime()-last_start_;
    last_start_ = MPI_Wtime();
}

void Timer::restart()
{
    smpi_->barrier();
    last_start_ = MPI_Wtime();
}

void Timer::print(double tot)
{
    if ((time_acc_>0.) && (name_!=""))
        MESSAGE(0, "\t" << setw(12) << name_ << "\t" << time_acc_  << "\t(" << 100.0*time_acc_/tot << "%)");
}

