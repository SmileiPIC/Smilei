#include "PyTools.h"
#include "Histogram.h"
#include "Patch.h"
#include "ParticleData.h"

#include <algorithm>

using namespace std;

// Loop on the different axes requested and compute the output index of each particle
void Histogram::digitize( vector<Species *> species,
                          vector<double> &double_buffer,
                          vector<int>    &int_buffer,
                          SimWindow *simWindow )
{
    unsigned int npart = double_buffer.size();
    
    for( unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++ ) {
        
        HistogramAxis * axis = axes[iaxis];
        
        // first loop on particles to store the indexing (axis) quantity
        unsigned int istart = 0;
        for( unsigned int ispec=0; ispec < species.size(); ispec++ ) {
            unsigned int npart = species[ispec]->getNbrOfParticles();
            axis->calculate_locations( species[ispec], &double_buffer[istart], &int_buffer[istart], npart, simWindow );
            istart += npart;
        }
        // Now, double_buffer has the location of each particle along the axis
        
        // if log scale, loop again and convert to log
        if( axis->logscale ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( int_buffer[ipart] < 0 ) {
                    continue;
                }
                double_buffer[ipart] = log10( abs( double_buffer[ipart] ) );
            }
        }
        
        double actual_min = axis->logscale ? log10( axis->global_min ) : axis->global_min;
        double actual_max = axis->logscale ? log10( axis->global_max ) : axis->global_max;
        double coeff = ( ( double ) axis->nbins )/( actual_max - actual_min );
        
        // The indexes are "reshaped" in one dimension.
        // For instance, in 3d, the index has the form  i = i3 + n3*( i2 + n2*i1 )
        // Here we do the multiplication by n3 or n2 (etc.)
        if( iaxis>0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                int_buffer[ipart] *= axis->nbins;
            }
        }
        
        // loop again on the particles and calculate the index
        // This is separated in two cases: edge_inclusive and edge_exclusive
        if( !axis->edge_inclusive ) { // if the particles out of the "box" must be excluded
        
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                // skip already discarded particles
                if( int_buffer[ipart] < 0 ) {
                    continue;
                }
                // calculate index
                int ind = floor( ( double_buffer[ipart]-actual_min ) * coeff );
                // index valid only if in the "box"
                if( ind >= 0  &&  ind < axis->nbins ) {
                    int_buffer[ipart] += ind;
                } else {
                    int_buffer[ipart] = -1;    // discard particle
                }
            }
            
        } else { // if the particles out of the "box" must be included
        
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                // skip already discarded particles
                if( int_buffer[ipart] < 0 ) {
                    continue;
                }
                // calculate index
                int ind = floor( ( double_buffer[ipart]-actual_min ) * coeff );
                // move out-of-range indexes back into range
                if( ind < 0 ) {
                    ind = 0;
                }
                if( ind >= axis->nbins ) {
                    ind = axis->nbins-1;
                }
                int_buffer[ipart] += ind;
            }
            
        }
        
    } // loop axes
}

void Histogram::distribute(
    std::vector<double> &double_buffer,
    std::vector<int>    &int_buffer,
    std::vector<double> &output_array )
{

    unsigned int ipart, npart=double_buffer.size();
    int ind;
    
    // Sum the data into the data_sum according to the indexes
    // ---------------------------------------------------------------
    for( ipart = 0 ; ipart < npart ; ipart++ ) {
        ind = int_buffer[ipart];
        if( ind<0 ) {
            continue;    // skip discarded particles
        }
        #pragma omp atomic
        output_array[ind] += double_buffer[ipart];
    }
    
}



void HistogramAxis::init( string type_, double min_, double max_, int nbins_, bool logscale_, bool edge_inclusive_, vector<double> coefficients_ )
{
    type           = type_          ;
    min            = min_           ;
    max            = max_           ;
    nbins          = nbins_         ;
    logscale       = logscale_      ;
    edge_inclusive = edge_inclusive_;
    coefficients   = coefficients_  ;
    global_min     = min;
    global_max     = max;
}
