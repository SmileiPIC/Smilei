#include "Histogram.h"

using namespace std;


Histogram::Histogram(vector<HistogramAxis*> axes) :
    axes(axes)
{}

Histogram::~Histogram()
{}

void Histogram::digitize(  Species *s,
    std::vector<double> &double_buffer,
    std::vector<int>    &int_buffer,
    std::vector<double> &output_array    )
{
    unsigned int ipart, npart=s->particles->size();
    int ind;
    
    fill(int_buffer.begin(), int_buffer.end(), 0);
    
    // 1 - loop on the different axes requested and compute the output index of each particle
    // --------------------------------------------------------------------------------------
    for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
        
        // first loop on particles to store the indexing (axis) quantity
        axes[iaxis]->digitize( s, double_buffer, npart );
        // Now, double_buffer has the location of each particle along the axis
        
        // if log scale, loop again and convert to log
        if (axes[iaxis]->logscale) {
            for (ipart = 0 ; ipart < npart ; ipart++)
                double_buffer[ipart] = log10(double_buffer[ipart]);
        }
        
        // The indexes are "reshaped" in one dimension.
        // For instance, in 3d, the index has the form  i = i3 + n3*( i2 + n2*i1 )
        // Here we do the multiplication by n3 or n2 (etc.)
        if (iaxis>0) {
            for (ipart = 0 ; ipart < npart ; ipart++)
                int_buffer[ipart] *= axes[iaxis]->nbins;
        }
        
        // loop again on the particles and calculate the index
        // This is separated in two cases: edge_inclusive and edge_exclusive
        if (!axes[iaxis]->edge_inclusive) { // if the particles out of the "box" must be excluded
            
            for (ipart = 0 ; ipart < npart ; ipart++) {
                // skip already discarded particles
                if (int_buffer[ipart] < 0) continue; 
                // calculate index
                ind = floor( (double_buffer[ipart]-axes[iaxis]->actual_min) * axes[iaxis]->coeff );
                // index valid only if in the "box"
                if (ind >= 0  &&  ind < axes[iaxis]->nbins) int_buffer[ipart] += ind;
                else int_buffer[ipart] = -1; // discard particle
            }
            
        } else { // if the particles out of the "box" must be included
            
            for (ipart = 0 ; ipart < npart ; ipart++) {
                // skip already discarded particles
                if (int_buffer[ipart] < 0) continue; 
                // calculate index
                ind = floor( (double_buffer[ipart]-axes[iaxis]->actual_min) * axes[iaxis]->coeff );
                // move out-of-range indexes back into range
                if (ind < 0) ind = 0;
                if (ind >= axes[iaxis]->nbins) ind = axes[iaxis]->nbins-1;
                int_buffer[ipart] += ind;
            }
            
        }
        
    } // loop axes
    
    // 2 - fill the buffer with the particle data
    // -----------------------------------------------
    data_filling( s, double_buffer, npart );
    
    // 3 - sum the data into the data_sum according to the indexes
    // ---------------------------------------------------------------
    for (ipart = 0 ; ipart < npart ; ipart++) {
        ind = int_buffer[ipart];
        if (ind<0) continue; // skip discarded particles
        #pragma omp atomic
        output_array[ind] += double_buffer[ipart];
    }
    
}



// Functions that are used to compute the data_array
void Histogram_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart];
}
void Histogram_charge_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart]);
}
void Histogram_jx_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart])
                     * s->particles->Momentum[0][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_jy_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart])
                     * s->particles->Momentum[1][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_jz_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart])
                     * s->particles->Momentum[2][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_ekin_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * ( sqrt(1. + pow(s->particles->Momentum[0][ipart],2)
                                 + pow(s->particles->Momentum[1][ipart],2)
                                 + pow(s->particles->Momentum[2][ipart],2)) - 1.);
}
void Histogram_p_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * sqrt(pow(s->particles->Momentum[0][ipart],2)
                          + pow(s->particles->Momentum[1][ipart],2)
                          + pow(s->particles->Momentum[2][ipart],2));
}
void Histogram_px_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart] * s->particles->Momentum[0][ipart];
}
void Histogram_py_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart] * s->particles->Momentum[1][ipart];
}
void Histogram_pz_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart] * s->particles->Momentum[2][ipart];
}
void Histogram_pressure_xx::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * pow(s->particles->Momentum[0][ipart],2)
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_pressure_yy::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * pow(s->particles->Momentum[1][ipart],2)
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_pressure_zz::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * pow(s->particles->Momentum[2][ipart],2)
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_pressure_xy::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[0][ipart]
                     * s->particles->Momentum[1][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_pressure_xz::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[0][ipart]
                     * s->particles->Momentum[2][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_pressure_yz::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[1][ipart]
                     * s->particles->Momentum[2][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void Histogram_ekin_vx_density::data_filling(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[0][ipart]
                     * (1. - 1./sqrt(1. + pow(s->particles->Momentum[0][ipart],2)
                                        + pow(s->particles->Momentum[1][ipart],2)
                                        + pow(s->particles->Momentum[2][ipart],2)));
}




HistogramAxis::HistogramAxis(
    std::string    type,
    double         min,
    double         max,
    int            nbins,
    bool           logscale,
    bool           edge_inclusive,
    vector<double> coefficients
) :
    type           ( type           ),
    min            ( min            ),
    max            ( max            ),
    nbins          ( nbins          ),
    logscale       ( logscale       ),
    edge_inclusive ( edge_inclusive ),
    coefficients   ( coefficients   )
{
    if( logscale ) {
        actual_min = log10(min);
        actual_max = log10(max);
    } else {
        actual_min = min;
        actual_max = max;
    }
    
    coeff = ((double)nbins)/(actual_max-actual_min);
}

HistogramAxis::~HistogramAxis()
{}

void HistogramAxis_x::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for ( unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Position[0][ipart];
}

void HistogramAxis_y::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Position[1][ipart];
}

void HistogramAxis_z::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Position[2][ipart];
}

void HistogramAxis_px::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Momentum[0][ipart];
}

void HistogramAxis_py::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Momentum[1][ipart];
}

void HistogramAxis_pz::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Momentum[2][ipart];
}

void HistogramAxis_p::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * sqrt(pow(s->particles->Momentum[0][ipart],2)
                                    + pow(s->particles->Momentum[1][ipart],2)
                                    + pow(s->particles->Momentum[2][ipart],2));
}

void HistogramAxis_gamma::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}

void HistogramAxis_ekin::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * (sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                           + pow(s->particles->Momentum[1][ipart],2)
                                           + pow(s->particles->Momentum[2][ipart],2) ) - 1.);
}

void HistogramAxis_vx::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Momentum[0][ipart]
                       / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                  + pow(s->particles->Momentum[1][ipart],2)
                                  + pow(s->particles->Momentum[2][ipart],2) );
}

void HistogramAxis_vy::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Momentum[1][ipart]
                       / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                  + pow(s->particles->Momentum[1][ipart],2)
                                  + pow(s->particles->Momentum[2][ipart],2) );
}

void HistogramAxis_vz::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Momentum[2][ipart]
                       / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                  + pow(s->particles->Momentum[1][ipart],2)
                                  + pow(s->particles->Momentum[2][ipart],2) );
}

void HistogramAxis_v::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = pow( 1. + 1./(pow(s->particles->Momentum[0][ipart],2)
                                   + pow(s->particles->Momentum[1][ipart],2)
                                   + pow(s->particles->Momentum[2][ipart],2)) , -0.5);
}

void HistogramAxis_vperp2::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = (  pow(s->particles->Momentum[1][ipart],2)
                        + pow(s->particles->Momentum[2][ipart],2)
                       ) / (1. + pow(s->particles->Momentum[0][ipart],2)
                               + pow(s->particles->Momentum[1][ipart],2)
                               + pow(s->particles->Momentum[2][ipart],2) );
}

void HistogramAxis_charge::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = (double) s->particles->Charge[ipart];
}

void HistogramAxis_chi::digitize(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Chi[ipart];
}

void HistogramAxis_composite::digitize(Species * s, vector<double> &array, unsigned int npart) {
    unsigned int idim, ndim = coefficients.size();
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++) {
        array[ipart] = 0.;
        for (idim = 0 ; idim < ndim ; idim++)
            array[ipart] += coefficients[idim] * s->particles->Position[idim][ipart];
    }
}
