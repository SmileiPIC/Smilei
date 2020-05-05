#ifndef INTERFACE_H
#define INTERFACE_H
#include "Field3D.h"
#include "Params.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include "ElectroMagn3D.h"
#include <iostream>

namespace picsar
{

extern  "C"
{


    void init_params_picsar( int *, int *, int *, double *, double *, double *, double *,
                             int *, int *, int *, int *, int *, int *, bool *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *, double *, double *, int * );
    void init_params_picsar_AM( int*, int *, int *, int *, double *, double *, double *,
                                int *, int *, int *, int *, 
                                std::complex<double> *, std::complex<double> *, std::complex<double> *,
                                std::complex<double> *, std::complex<double> *, std::complex<double> *, 
                                std::complex<double> *, std::complex<double> *, std::complex<double> *, 
                                std::complex<double> *, std::complex<double> * );
    void free_params_picsar_AM();
    void push_psatd_ebfield_();
    void solve_maxwell_fdtd_pxr();
    void rotational_cleaning();
    void densities_correction();
};

}

void copy_field_3d( Field3D *out, Field3D *in );
void copy_field_2d( Field2D *out, Field2D *in );
void duplicate_field_into_pxr( ElectroMagn * );
void duplicate_field_into_smilei( ElectroMagn * );


#endif
