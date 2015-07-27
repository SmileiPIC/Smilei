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
class Patch;


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
    friend class VectorPatch;
    friend class SmileiMPI;
public:
    //! creator (called from Diagnostic)
    DiagnosticScalar(PicParams &params, DiagParams &diagParams, Patch* patch);
    //! destructor
    ~DiagnosticScalar(){};
    
    //! close the file
    void close();

    //! calls the compute_proc_gather, compute and write
    void run(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);

    //! ask to each processor to compute the scalars and gather them in the map mpi_spec_scalars[cpu][species]
    void compute(ElectroMagn* EMfields, std::vector<Species*>&);

    //! write the out_list data onto a file
    void write(int timestep);

    //! get a particular scalar
    double getScalar(std::string name);
    //! get a particular scalar
    void setScalar(std::string name, double value);

    std::vector<std::pair<std::string,double> >::iterator itDiagScalar;

private:
    //! check if patch is master (from patch)
    const bool isMaster;
    
    //! initial energy (kinetic + EM)
    double Energy_time_zero;
    
    //! energy used for the normalization of energy balance (former total energy)
    double EnergyUsedForNorm;
    
    //! this is copied from params
    const double res_time;
    
    //! every step to calculate scalars
    const unsigned int every;
    
    //! output stream
    std::ofstream fout;
    
    //! copied from params
    double cell_volume;
    
    //! write precision
    unsigned int precision;
    
    //! this is a list to keep variable name and value
    std::vector<std::pair<std::string,double> > out_list;
        
    //! append to outlist
    void append(std::string, double);

    //! prepend to outlist
    void prepend(std::string, double);

    //! list of keys for scalars to be written
    std::vector<std::string> vars;

    //! check if key is allowed
    bool allowedKey(std::string);

};

#endif
