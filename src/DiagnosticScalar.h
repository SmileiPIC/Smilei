#ifndef DIAGNOSTICSCALAR_H
#define DIAGNOSTICSCALAR_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Tools.h"
//#include "Diagnostic.h"
#include "Species.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;


struct val_index
{
    double val;
    int index;
};


//! class that calculates scalars and writes them on a file

//! the user who wants to implement a scalar diagnostic, can fill the scalars map in species::computeScalar
class DiagnosticScalar {

public:
    //! creator (called from Diagnostic)
    DiagnosticScalar(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi);
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

    //! this is a list to keep variable name and value
    std::vector<std::pair<std::string,double> > out_list;
    
    //! get a particular scalar
    double getScalar(std::string name);

private:
    //! check if proc is master (from smpi)
    const bool isMaster;
    
    //! tot number of cpus (from smpi)
    const unsigned int cpuSize;

    //! initial energy (kinetic + EM)
    double Energy_time_zero;
    
    //! this is copied from params
    const double res_time;
    
    //! every step to calculate scalars
    const unsigned int every;
    
    //! output stream
    std::ofstream fout;
    
    //! mpi_spec_scalars [iCpu][iSpec]
    std::vector<std::vector<std::map<std::string, double> > > mpi_spec_scalars;

    //! mpi_EM_scalars [iCpu]["Field name"]["key"]<values>
    std::vector<std::map<std::string,std::map<std::string,std::vector<double> > > > mpi_EM_scalars;
    
    //! copied from params
    double cell_volume;
    
    //! write precision
    unsigned int precision;
};

#endif
