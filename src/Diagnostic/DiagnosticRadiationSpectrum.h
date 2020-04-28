#ifndef DIAGNOSTICRADIATIONSPECTRUM_H
#define DIAGNOSTICRADIATIONSPECTRUM_H

#include "DiagnosticParticleBinning.h"

class DiagnosticRadiationSpectrum : public DiagnosticParticleBinning
{
    friend class SmileiMPI;

public :

    //! Default constructor
    DiagnosticRadiationSpectrum( Params &params, SmileiMPI *smpi, Patch *patch, RadiationTables *radiation_tables_, int diagId );
    //! Default destructor
    ~DiagnosticRadiationSpectrum();
    
    void openFile( Params &params, SmileiMPI *smpi, bool newfile ) override;
    
    void run( Patch *patch, int timestep, SimWindow *simWindow ) override;

private :

    //! constant 2/3
    double two_third;
    
    //! Extra axis for photon energies
    HistogramAxis *photon_axis;
    
    //! Minimum photon energy for radiation spectrum
    double minimum_chi_continuous_;
    
    //! axis containing the values of the binned photon_energies
    std::vector<double> photon_energies;
    
    //! axis containing the values of the delta on the binned photon_energies
    std::vector<double> delta_energies;

};

#endif

