#ifndef DIAGNOSTICSCALAR_H
#define DIAGNOSTICSCALAR_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Tools.h"
#include "Species.h"
#include "TimeSelection.h"

class Params;
class ElectroMagn;
class DiagParams;
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
class DiagnosticScalar {    friend class VectorPatch;
    friend class SmileiMPI;
    friend class SimWindow;
public:
    //! creator (called from Diagnostic)
    DiagnosticScalar(Params &params, Patch* patch);
    //! destructor
    ~DiagnosticScalar(){};
    
    //! close the file
    void closeFile();


    //! close the file
    void open();
    void open(bool append);

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

    //std::vector<std::pair<std::string,double> >::iterator itDiagScalar;
    std::vector<std::string>::iterator itDiagScalarName;
    std::vector<double>::iterator itDiagScalarValue;

    //! check if patch is master (from patch)
    bool isMaster;
    
//    //! every step to calculate scalars
//    unsigned int every;

    //! Time selection
    TimeSelection * timeSelection;
    
    //! this is copied from params
    double res_time;
    
//    double tmin;
//    double tmax;
    double dt;

    //! write precision
    unsigned int precision;

    //! this is a list to keep variable name and value
    std::vector<std::pair<std::string,double> > out_list;
    
    //! these are lists to keep variable names and values
    std::vector<std::string> out_key;
    std::vector<double>      out_value;
    //! width of each field
    std::vector<unsigned int> out_width;

    //! list of keys for scalars to be written
    std::vector<std::string> vars;
    
    //! copied from params
    double cell_volume;
    
    //! initial energy (kinetic + EM)
    double Energy_time_zero;
    
    //! energy used for the normalization of energy balance (former total energy)
    double EnergyUsedForNorm;
    
        
private:    
    //! append to outlist
    void append(std::string, double);

    //! prepend to outlist
    void prepend(std::string, double);

    //! check if key is allowed
    bool allowedKey(std::string);

protected :    
    //! output stream
    std::ofstream fout;

};

#endif
