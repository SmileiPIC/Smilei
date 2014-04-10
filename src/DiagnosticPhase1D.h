#ifndef DiagnosticPhase1D_H
#define DiagnosticPhase1D_H

#include <cmath>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <hdf5.h>

#include "Tools.h"

#include "Species.h"
#include "Interpolator.h"
#include "Particles.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

class DiagnosticPhase1D {

public:

    DiagnosticPhase1D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi);
    ~DiagnosticPhase1D();

    void set_proc();
    void set_file_name();

    void run(int timestep, std::vector<Species*>&);

    void open_file();
    void close();
	
private:
    SmileiMPI* smpi_;
	
	double momentum_min;
	std::vector<double> momentumDistr;

};
#endif
