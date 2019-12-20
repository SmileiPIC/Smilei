#ifndef DIAGNOSTICRADIATIONSPECTRUM_H
#define DIAGNOSTICRADIATIONSPECTRUM_H

#include "Diagnostic.h"

#include "Histogram.h"

class DiagnosticRadiationSpectrum : public Diagnostic {
    friend class SmileiMPI;

public :

    //! Default constructor
    DiagnosticRadiationSpectrum( Params &params, SmileiMPI* smpi, Patch* patch, int diagId );
    //! Cloning constructor
    DiagnosticRadiationSpectrum( DiagnosticRadiationSpectrum* );
    //! Default destructor
    ~DiagnosticRadiationSpectrum() override;

    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;

    void closeFile() override;

    bool prepare( int timestep ) override;

    void run( Patch* patch, int timestep, SimWindow* simWindow ) override;

    void write(int timestep, SmileiMPI* smpi) override;

    //! Clear the array
    void clear();

     //! Get memory footprint of current diagnostic
    int getMemFootPrint() override {
        int size = output_size*sizeof(double);
        // + data_array + index_array +  axis_array
        // + nparts_max * (sizeof(double)+sizeof(int)+sizeof(double))
        return size;
    };

    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint(int istart, int istop, Patch* patch) override;

private :

    //! constant 2/3
    double two_third;

    //! normalization factor for the emitted power spectrum
    double factor;

    //! number of timesteps during which outputs are averaged
    int time_average;

    //! list of the species that will be accounted for
    std::vector<unsigned int> species;

    //! vector for saving the output array for time-averaging
    std::vector<double> data_sum;

    //! Histogram object
    Histogram * histogram;

    unsigned int output_size;

    //! Minimum and maximum spatial coordinates that are useful for this diag
    std::vector<double> spatial_min, spatial_max;

    //! Minimum photon energy for radiation spectrum
    double minimum_chi_continuous_;

    //! Minimum photon energy for radiation spectrum
    double photon_energy_min;

    //! Maximum photon energy for radiation spectrum
    double photon_energy_max;

    //! Spacing in the linear or log10 scale
    double spacing;

    //! Number of energy bins for radiation spectrum
    int photon_energy_nbins;

    //! is logscale used for the photon energy axis
    bool photon_energy_logscale;

    //! are edges inclusive for the photon energy axis
    bool photon_energy_edge_inclusive;

    //! axis containing the values of the binned photon_energies
    std::vector<double> photon_energies;

    //! axis containing the values of the delta on the binned photon_energies
    std::vector<double> delta_energies;

};

#endif
