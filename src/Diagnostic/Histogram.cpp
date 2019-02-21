#include "PyTools.h"
#include "Histogram.h"
#include "Patch.h"
#include "ParticleData.h"


using namespace std;

// Loop on the different axes requested and compute the output index of each particle
void Histogram::digitize( Species *s,
                          std::vector<double> &double_buffer,
                          std::vector<int>    &int_buffer,
                          SimWindow *simWindow )
{
    unsigned int ipart, npart=s->particles->size();
    int ind;
    
    for( unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++ ) {
    
        // first loop on particles to store the indexing (axis) quantity
        axes[iaxis]->digitize( s, double_buffer, int_buffer, npart, simWindow );
        // Now, double_buffer has the location of each particle along the axis
        
        // if log scale, loop again and convert to log
        if( axes[iaxis]->logscale ) {
            for( ipart = 0 ; ipart < npart ; ipart++ ) {
                if( int_buffer[ipart] < 0 ) {
                    continue;
                }
                double_buffer[ipart] = log10( abs( double_buffer[ipart] ) );
            }
        }
        
        // The indexes are "reshaped" in one dimension.
        // For instance, in 3d, the index has the form  i = i3 + n3*( i2 + n2*i1 )
        // Here we do the multiplication by n3 or n2 (etc.)
        if( iaxis>0 ) {
            for( ipart = 0 ; ipart < npart ; ipart++ ) {
                int_buffer[ipart] *= axes[iaxis]->nbins;
            }
        }
        
        // loop again on the particles and calculate the index
        // This is separated in two cases: edge_inclusive and edge_exclusive
        if( !axes[iaxis]->edge_inclusive ) { // if the particles out of the "box" must be excluded
        
            for( ipart = 0 ; ipart < npart ; ipart++ ) {
                // skip already discarded particles
                if( int_buffer[ipart] < 0 ) {
                    continue;
                }
                // calculate index
                ind = floor( ( double_buffer[ipart]-axes[iaxis]->actual_min ) * axes[iaxis]->coeff );
                // index valid only if in the "box"
                if( ind >= 0  &&  ind < axes[iaxis]->nbins ) {
                    int_buffer[ipart] += ind;
                } else {
                    int_buffer[ipart] = -1;    // discard particle
                }
            }
            
        } else { // if the particles out of the "box" must be included
        
            for( ipart = 0 ; ipart < npart ; ipart++ ) {
                // skip already discarded particles
                if( int_buffer[ipart] < 0 ) {
                    continue;
                }
                // calculate index
                ind = floor( ( double_buffer[ipart]-axes[iaxis]->actual_min ) * axes[iaxis]->coeff );
                // move out-of-range indexes back into range
                if( ind < 0 ) {
                    ind = 0;
                }
                if( ind >= axes[iaxis]->nbins ) {
                    ind = axes[iaxis]->nbins-1;
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
    
    if( logscale ) {
        actual_min = log10( min );
        actual_max = log10( max );
    } else {
        actual_min = min;
        actual_max = max;
    }
    
    coeff = ( ( double )nbins )/( actual_max-actual_min );
}
