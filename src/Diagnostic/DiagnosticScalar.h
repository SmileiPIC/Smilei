#ifndef DIAGNOSTICSCALAR_H
#define DIAGNOSTICSCALAR_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Tools.h"
#include "Species.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;


//! double-int structure to communicate min/max and location trough MPI 
struct val_index
{
    //! min/max 
    double val;
    //! cell location index
    int index;
};


//! class that calculates scalars and writes them on a file

//! the user who wants to implement a scalar diagnostic, can fill the scalars map in species::computeScalar
class DiagnosticScalar {

public:
    //! creator (called from Diagnostic)
    DiagnosticScalar(SmileiMPI* smpi);
    //! destructor
    ~DiagnosticScalar(){};
    
    //! close the file
    void close();

    //! calls the compute_proc_gather, compute and write
    void run(int timestep, ElectroMagn* EMfields, std::vector<Species*>&, SmileiMPI *smpi);

    //! ask to each processor to compute the scalars and gather them in the map mpi_spec_scalars[cpu][species]
    void compute(ElectroMagn* EMfields, std::vector<Species*>&, SmileiMPI *smpi);

    //! write the out_list data onto a file
    void write(int timestep);

    //! get a particular scalar
    double getScalar(std::string name);

    //! every step to calculate scalars
    unsigned int every;
    
    //! this is copied from params
    double res_time;
    
    double tmin;
    double tmax;
    double dt;

    //! write precision
    unsigned int precision;
    
    //! this is a list to keep variable name and value
    std::vector<std::pair<std::string,double> > out_list;
    
    //! list of keys for scalars to be written
    std::vector<std::string> vars;
    
    //! copied from params
    double cell_volume;
    
private:
    //! check if proc is master (from smpi)
    const bool isMaster;
    
    //! tot number of cpus (from smpi)
    const unsigned int cpuSize;

    //! initial energy (kinetic + EM)
    double Energy_time_zero;
    
    //! energy used for the normalization of energy balance (former total energy)
    double EnergyUsedForNorm;
    
    
    //! output stream
    std::ofstream fout;
    
    //! append to outlist
    void append(std::string, double);

    //! prepend to outlist
    void prepend(std::string, double);

    //! check if key is allowed
    bool allowedKey(std::string);

};

#endif
