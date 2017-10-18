#ifndef HISTOGRAMFACTORY_H
#define HISTOGRAMFACTORY_H

#include "PyTools.h"
#include "Histogram.h"
#include "Params.h"
#include "Patch.h"
#include "ParticleData.h"
#include <sstream>

class HistogramFactory {
public:
    static Histogram* create(
        Params &params,
        PyObject* deposited_quantity_object,
        std::vector<PyObject*> &pyAxes,
        std::vector<unsigned int> &species,
        Patch* patch,
        std::vector<std::string> &excluded_axes,
        std::string errorPrefix
    )
    {
        
        Histogram * histogram;
        std::string deposited_quantity = "";
        std::ostringstream deposited_quantityName( "" );
        deposited_quantityName << errorPrefix << ": parameter `deposited_quantity`";
        std::string deposited_quantityPrefix = deposited_quantityName.str();
        
        // By default, deposited_quantity=None, but should not
        if( deposited_quantity_object == Py_None ) {
            ERROR(deposited_quantityPrefix << " required");
            
        // If string, then ok
        } else if( PyTools::convert(deposited_quantity_object, deposited_quantity) ) {
            
            if ( deposited_quantity == "user_function" ) {
                ERROR(deposited_quantityPrefix << " = " << deposited_quantity <<" not understood");
            } else if (deposited_quantity == "weight"        ) {
                histogram = new Histogram_density        ();
            } else if (deposited_quantity == "weight_charge" ) {
                histogram = new Histogram_charge_density ();
            } else if (deposited_quantity == "weight_charge_vx" ) {
                histogram = new Histogram_jx_density     ();
            } else if (deposited_quantity == "weight_charge_vy" ) {
                histogram = new Histogram_jy_density     ();
            } else if (deposited_quantity == "weight_charge_vz" ) {
                histogram = new Histogram_jz_density     ();
            } else if (deposited_quantity == "weight_ekin"   ) {
                histogram = new Histogram_ekin_density   ();
            } else if (deposited_quantity == "weight_p"      ) {
                histogram = new Histogram_p_density      ();
            } else if (deposited_quantity == "weight_px"     ) {
                histogram = new Histogram_px_density     ();
            } else if (deposited_quantity == "weight_py"     ) {
                histogram = new Histogram_py_density     ();
            } else if (deposited_quantity == "weight_pz"     ) {
                histogram = new Histogram_pz_density     ();
            } else if (deposited_quantity == "weight_vx_px"    ) {
                histogram = new Histogram_pressure_xx    ();
            } else if (deposited_quantity == "weight_vy_py"    ) {
                histogram = new Histogram_pressure_yy    ();
            } else if (deposited_quantity == "weight_vz_pz"    ) {
                histogram = new Histogram_pressure_zz    ();
            } else if (deposited_quantity == "weight_vx_py"    ) {
                histogram = new Histogram_pressure_xy    ();
            } else if (deposited_quantity == "weight_vx_pz"    ) {
                histogram = new Histogram_pressure_xz    ();
            } else if (deposited_quantity == "weight_vy_pz"    ) {
                histogram = new Histogram_pressure_yz    ();
            } else if (deposited_quantity == "weight_ekin_vx") {
                histogram = new Histogram_ekin_vx_density();
            } else if (deposited_quantity == "weight_chi" ) {
                histogram = new Histogram_chi_density    (patch, species, errorPrefix);
            } else {
                ERROR(deposited_quantityPrefix << " not understood");
            }
            histogram->deposited_quantity = deposited_quantity;
            
        // If numpy supported, also accept deposited_quantity = any function
        } else {
#ifdef SMILEI_USE_NUMPY
            // Test the function with temporary, "fake" particles
            double * dummy = NULL;
            ParticleData test( params.nDim_particle, deposited_quantity_object, deposited_quantityPrefix, dummy );
            histogram = new Histogram_user_function(deposited_quantity_object);
            histogram->deposited_quantity = "user_function";
#else
            ERROR(deposited_quantityPrefix << " should be a string");
#endif
        }
        
        // Now setup each axis
        std::string type;
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
            PyObject * type_object = PySequence_Fast_GET_ITEM(seq, 0);
            if ( PyTools::convert(type_object, type) ) {
                if ( type == "user_function" )
                    ERROR(errorPrefix << ", axis #" << iaxis << ": type " << type << " unknown");
                for( unsigned int i=0; i<excluded_axes.size(); i++ )
                    if( type == excluded_axes[i] )
                        ERROR(errorPrefix << ", axis #" << iaxis << ": type " << type << " unknown");
            // If numpy supported, also accept type = any function
            } else {
                std::ostringstream typePrefix("");
                typePrefix << errorPrefix << ", axis #" << iaxis << ": type";
#ifdef SMILEI_USE_NUMPY
                // Test the function with temporary, "fake" particles
                double * dummy = NULL;
                ParticleData test( params.nDim_particle, type_object, typePrefix.str(), dummy );
                type = "user_function";
#else
                ERROR(errorPrefix << ", axis #" << iaxis << ": First item must be a string (axis type)");
#endif
            }
            
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
                std::string my_str("");
                PyTools::convert(PySequence_Fast_GET_ITEM(seq, i),my_str);
                if(my_str=="logscale" ||  my_str=="log_scale" || my_str=="log")
                    logscale = true;
                else if(my_str=="edges" ||  my_str=="edge" ||  my_str=="edge_inclusive" ||  my_str=="edges_inclusive")
                    edge_inclusive = true;
                else
                    ERROR(errorPrefix << ": keyword `" << my_str << "` not understood");
            }
            
            HistogramAxis * axis;
            std::vector<double> coefficients(0);
            if        (type == "x" ) {
                axis = new HistogramAxis_x();
            } else if (type == "moving_x" ) {
                axis = new HistogramAxis_moving_x();
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
                    if( ! patch->vecSpecies[species[ispec]]->particles->isQuantumParameter )
                        ERROR(errorPrefix << ": axis #" << iaxis << " 'chi' requires all species to be 'radiating'");
                axis = new HistogramAxis_chi();
            } else if (type == "composite") {
                ERROR(errorPrefix << ": axis type cannot be 'composite'");
            } else if (type == "user_function") {
                axis = new HistogramAxis_user_function( type_object );
            
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
                        std::string segment = type.substr(previ,i-previ);
                        // Get the last character, which should be one of x, y, or z
                        unsigned int j = segment.length();
                        std::string direction = j>0 ? segment.substr(j-1,1) : "";
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
                        coefficients[direction_index] = 1.;
                        if( j>1 ) {
                            std::stringstream number("");
                            number << segment.substr(0,j-1);
                            number >> coefficients[direction_index];
                            if( ! coefficients[direction_index] )
                                ERROR(errorPrefix << ": axis #" << iaxis << " type not understood");
                        }
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
            histogram->axes.push_back( axis );
        }
        
        return histogram;
    }
};

#endif
