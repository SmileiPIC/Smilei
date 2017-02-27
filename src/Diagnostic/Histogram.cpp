#include "Histogram.h"
#include "PyTools.h"
#include "Patch.h"

using namespace std;


void Histogram::init( Params &params, vector<PyObject*> pyAxes, vector<unsigned int> species,
    string errorPrefix, Patch * patch, vector<string> excluded_types )
{
    
    string type;
    double min, max;
    int nbins;
    bool logscale, edge_inclusive;
    
    if (pyAxes.size() == 0)
        ERROR(errorPrefix << ": axes must contain something");
    
    // Loop axes and extract their format
    for (unsigned int iaxis=0; iaxis<pyAxes.size(); iaxis++ ) {
        PyObject *pyAxis=pyAxes[iaxis];
        
        // Axis must be a list
        if (!PyTuple_Check(pyAxis) && !PyList_Check(pyAxis))
            ERROR(errorPrefix << ": axis #" << iaxis << " must be a list");
        PyObject* seq = PySequence_Fast(pyAxis, "expected a sequence");
        
        // Axis must have 4 elements or more
        unsigned int lenAxisArgs=PySequence_Size(seq);
        if (lenAxisArgs<4)
            ERROR(errorPrefix << ": axis #" << iaxis << " contain at least 4 arguments");
        
        // Try to extract first element: type
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 0), type))
            ERROR(errorPrefix << ", axis #" << iaxis << ": First item must be a string (axis type)");
        for( unsigned int i=0; i<excluded_types.size(); i++ )
            if( type == excluded_types[i] )
                ERROR(errorPrefix << ", axis #" << iaxis << ": type " << type << " unknown");
        
        // Try to extract second element: axis min
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 1), min)) {
            ERROR(errorPrefix << ", axis #" << iaxis << ": Second item must be a double (axis min)");
        }
        
        // Try to extract third element: axis max
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 2), max)) {
            ERROR(errorPrefix << ", axis #" << iaxis << ": Third item must be a double (axis max)");
        }
        
        // Try to extract fourth element: axis nbins
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 3), nbins)) {
            ERROR(errorPrefix << ", axis #" << iaxis << ": Fourth item must be an int (number of bins)");
        }
        
        // Check for  other keywords such as "logscale" and "edge_inclusive"
        logscale = false;
        edge_inclusive = false;
        for(unsigned int i=4; i<lenAxisArgs; i++) {
            string my_str("");
            PyTools::convert(PySequence_Fast_GET_ITEM(seq, i),my_str);
            if(my_str=="logscale" ||  my_str=="log_scale" || my_str=="log")
                logscale = true;
            else if(my_str=="edges" ||  my_str=="edge" ||  my_str=="edge_inclusive" ||  my_str=="edges_inclusive")
                edge_inclusive = true;
            else
                ERROR(errorPrefix << ": keyword `" << my_str << "` not understood");
        }
        
        HistogramAxis * axis;
        vector<double> coefficients(0);
        if        (type == "x" ) {
            axis = new HistogramAxis_x();
        } else if (type == "y" ) {
            if (params.nDim_particle <2)
                ERROR(errorPrefix << ": axis y cannot exist in <2D");
            axis = new HistogramAxis_y();
        } else if (type == "z" ) {
            if (params.nDim_particle <3)
                ERROR(errorPrefix << ": axis z cannot exist in <3D");
            axis = new HistogramAxis_z();
        } else if (type == "a" ) {
            if (params.nDim_particle <2)
                ERROR(errorPrefix << ": axis a cannot exist in <2D");
            axis = new HistogramAxis_vector();
        } else if (type == "b" ) {
            if (params.nDim_particle <3)
                ERROR(errorPrefix << ": axis b cannot exist in <3D");
            axis = new HistogramAxis_vector();
        } else if (type == "theta" ) {
            if (params.nDim_particle == 1) {
                ERROR(errorPrefix << ": axis theta cannot exist in 1D");
            } else if (params.nDim_particle == 2) {
                axis = new HistogramAxis_theta2D();
            } else if (params.nDim_particle == 3) {
                axis = new HistogramAxis_theta3D();
            } else{
                ERROR(errorPrefix << ": impossible");
            }
        } else if (type == "phi" ) {
            if (params.nDim_particle <3)
                ERROR(errorPrefix << ": axis phi cannot exist in <3D");
            axis = new HistogramAxis_phi();
        } else if (type == "px" ) {
            axis = new HistogramAxis_px();
        } else if (type == "py" ) {
            axis = new HistogramAxis_py();
        } else if (type == "pz" ) {
            axis = new HistogramAxis_pz();
        } else if (type == "p" ) {
            axis = new HistogramAxis_p();
        } else if (type == "gamma" ) {
            axis = new HistogramAxis_gamma();
        } else if (type == "ekin" ) {
            axis = new HistogramAxis_ekin();
        } else if (type == "vx" ) {
            axis = new HistogramAxis_vx();
        } else if (type == "vy" ) {
            axis = new HistogramAxis_vy();
        } else if (type == "vz" ) {
            axis = new HistogramAxis_vz();
        } else if (type == "v" ) {
            axis = new HistogramAxis_v();
        } else if (type == "vperp2" ) {
            axis = new HistogramAxis_vperp2();
        } else if (type == "charge" ) {
            axis = new HistogramAxis_charge();
        } else if (type == "chi" ) {
            // The requested species must be radiating
            for (unsigned int ispec=0 ; ispec < species.size() ; ispec++)
                if( ! patch->vecSpecies[species[ispec]]->particles->isRadReaction )
                    ERROR(errorPrefix << ": axis #" << iaxis << " 'chi' requires all species to be 'radiating'");
            axis = new HistogramAxis_chi();
        } else if (type == "composite") {
            ERROR(errorPrefix << ": axis type cannot be 'composite'");
        
        } else {
            // If not "usual" type, try to find composite type
            for( unsigned int i=1; i<=type.length(); i++ )
                if( type.substr(i,1) == " " )
                    ERROR(errorPrefix << ": axis #" << iaxis << " type cannot contain whitespace");
            if( type.length()<2 )
                ERROR(errorPrefix << ": axis #" << iaxis << " type not understood");
            
            // Analyse character by character
            coefficients.resize( params.nDim_particle , 0. );
            unsigned int previ=0;
            double sign=1.;
            type += "+";
            for( unsigned int i=1; i<=type.length(); i++ ) {
                // Split string at "+" or "-" location
                if( type.substr(i,1) == "+" || type.substr(i,1) == "-" ) {
                    // Get one segment of the split string
                    string segment = type.substr(previ,i-previ);
                    // Get the last character, which should be one of x, y, or z
                    unsigned int j = segment.length();
                    string direction = j>0 ? segment.substr(j-1,1) : "";
                    unsigned int direction_index;
                    if     ( direction == "x" ) direction_index = 0;
                    else if( direction == "y" ) direction_index = 1;
                    else if( direction == "z" ) direction_index = 2;
                    else { ERROR(errorPrefix << ": axis #" << iaxis << " type not understood"); }
                    if( direction_index >= params.nDim_particle )
                        ERROR(errorPrefix << ": axis #" << iaxis << " type " << direction << " cannot exist in " << params.nDim_particle << "D");
                    if( coefficients[direction_index] != 0. )
                        ERROR(errorPrefix << ": axis #" << iaxis << " type " << direction << " appears twice");
                    // Get the remaining characters, which should be a number
                    coefficients[direction_index] = j>1 ? ::atof(segment.substr(0,j-1).c_str()) : 1.;
                    coefficients[direction_index] *= sign;
                    // Save sign and position for next segment
                    sign = type.substr(i,1) == "+" ? 1. : -1;
                    previ = i+1;
                }
            }
            
            type = "composite:"+type.substr(0,type.length()-1);
            axis = new HistogramAxis_composite();
        }
        
        Py_DECREF(seq);
        
        axis->init( type, min, max, nbins, logscale, edge_inclusive, coefficients );
        axes.push_back( axis );
    }
    
}


// Loop on the different axes requested and compute the output index of each particle
void Histogram::digitize(  Species *s,
    std::vector<double> &double_buffer,
    std::vector<int>    &int_buffer    )
{
    unsigned int ipart, npart=s->particles->size();
    int ind;
    
    for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
        
        // first loop on particles to store the indexing (axis) quantity
        axes[iaxis]->digitize( s, double_buffer, int_buffer, npart );
        // Now, double_buffer has the location of each particle along the axis
        
        // if log scale, loop again and convert to log
        if (axes[iaxis]->logscale) {
            for (ipart = 0 ; ipart < npart ; ipart++) {
                if( int_buffer[ipart] < 0 ) continue;
                double_buffer[ipart] = log10(double_buffer[ipart]);
            }
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
}

void Histogram::distribute(
    std::vector<double> &double_buffer,
    std::vector<int>    &int_buffer,
    std::vector<double> &output_array    )
{
    
    unsigned int ipart, npart=double_buffer.size();
    int ind;
    
    // Sum the data into the data_sum according to the indexes
    // ---------------------------------------------------------------
    for (ipart = 0 ; ipart < npart ; ipart++) {
        ind = int_buffer[ipart];
        if (ind<0) continue; // skip discarded particles
        #pragma omp atomic
        output_array[ind] += double_buffer[ipart];
    }
    
}



void HistogramAxis::init( string type_, double min_, double max_, int nbins_, bool logscale_, bool edge_inclusive_, vector<double> coefficients_)
{
    type           = type_          ;
    min            = min_           ;
    max            = max_           ;
    nbins          = nbins_         ;
    logscale       = logscale_      ;
    edge_inclusive = edge_inclusive_;
    coefficients   = coefficients_  ;
    
    if( logscale ) {
        actual_min = log10(min);
        actual_max = log10(max);
    } else {
        actual_min = min;
        actual_max = max;
    }
    
    coeff = ((double)nbins)/(actual_max-actual_min);
}
