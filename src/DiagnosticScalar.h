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

//! class that calculates scalars and writes them on a file

//! the user who wants to implement a scalar diagnostic, can fill the scalars map in species::computeScalar
class DiagnosticScalar {

public:
    // creator (called from Diagnostic)
    DiagnosticScalar(PicParams* params, SmileiMPI* smpi);
    //destructor
    ~DiagnosticScalar();

    //! calls the compute_proc_gather, compute and write
    void run(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);

    //! ask to each processor to compute the scalars and gather them in the map mpi_spec_scalars[cpu][species]
    void compute_proc_gather(ElectroMagn* EMfields, std::vector<Species*>&);

    //! fill the out_list based on mpi_spec_scalars[cpu][species] (user should just change this)
    void compute();

    //! write the out_list data onto a file
    void write(int timestep);

private:
    std::ofstream fout;
    SmileiMPI* smpi_;
    //! mpi_spec_scalars [iCpu][iSpec]
    std::vector<std::vector<std::map<std::string, double> > > mpi_spec_scalars;

    //! mpi_EM_scalars [iCpu]["Field name"]["key"]<values>
    std::vector<std::map<std::string,std::map<std::string,std::vector<double> > > > mpi_EM_scalars;

    //! this is a list to keep variable name and value
    std::vector<std::pair<std::string,double> > out_list;
};

#endif
