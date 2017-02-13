
#include "DiagnosticParticles.h"

#include <iomanip>

using namespace std;



// Functions that are used to compute the axis_array
void axis_x        (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for ( unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Position[0][ipart];
}
void axis_y        (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Position[1][ipart];
}
void axis_z        (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Position[2][ipart];
}
void axis_px       (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Momentum[0][ipart];
}
void axis_py       (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Momentum[1][ipart];
}
void axis_pz       (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Momentum[2][ipart];
}
void axis_p        (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * sqrt(pow(s->particles->Momentum[0][ipart],2)
                                    + pow(s->particles->Momentum[1][ipart],2)
                                    + pow(s->particles->Momentum[2][ipart],2));
}
void axis_gamma    (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void axis_ekin     (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * (sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                           + pow(s->particles->Momentum[1][ipart],2)
                                           + pow(s->particles->Momentum[2][ipart],2) ) - 1.);
}
void axis_vx       (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Momentum[0][ipart]
                       / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                  + pow(s->particles->Momentum[1][ipart],2)
                                  + pow(s->particles->Momentum[2][ipart],2) );
}
void axis_vy       (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Momentum[1][ipart]
                       / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                  + pow(s->particles->Momentum[1][ipart],2)
                                  + pow(s->particles->Momentum[2][ipart],2) );
}
void axis_vz       (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Momentum[2][ipart]
                       / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                  + pow(s->particles->Momentum[1][ipart],2)
                                  + pow(s->particles->Momentum[2][ipart],2) );
}
void axis_v        (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = pow( 1. + 1./(pow(s->particles->Momentum[0][ipart],2)
                                   + pow(s->particles->Momentum[1][ipart],2)
                                   + pow(s->particles->Momentum[2][ipart],2)) , -0.5);
}
void axis_vperp2   (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = (  pow(s->particles->Momentum[1][ipart],2)
                        + pow(s->particles->Momentum[2][ipart],2)
                       ) / (1. + pow(s->particles->Momentum[0][ipart],2)
                               + pow(s->particles->Momentum[1][ipart],2)
                               + pow(s->particles->Momentum[2][ipart],2) );
}
void axis_charge   (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = (double) s->particles->Charge[ipart];
}
void axis_chi      (Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Chi[ipart];
}
void axis_composite(Species * s, vector<double> &array, unsigned int npart, DiagnosticParticlesAxis &axis) {
    unsigned int idim, ndim = axis.coefficients.size();
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++) {
        array[ipart] = 0.;
        for (idim = 0 ; idim < ndim ; idim++)
            array[ipart] += axis.coefficients[idim] * s->particles->Position[idim][ipart];
    }
}


// Functions that are used to compute the data_array
void data_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart];
}
void data_charge_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart]);
}
void data_jx_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart])
                     * s->particles->Momentum[0][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_jy_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart])
                     * s->particles->Momentum[1][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_jz_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->particles->Weight[ipart] * (double)(s->particles->Charge[ipart])
                     * s->particles->Momentum[2][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_ekin_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * ( sqrt(1. + pow(s->particles->Momentum[0][ipart],2)
                                 + pow(s->particles->Momentum[1][ipart],2)
                                 + pow(s->particles->Momentum[2][ipart],2)) - 1.);
}
void data_p_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * sqrt(pow(s->particles->Momentum[0][ipart],2)
                          + pow(s->particles->Momentum[1][ipart],2)
                          + pow(s->particles->Momentum[2][ipart],2));
}
void data_px_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart] * s->particles->Momentum[0][ipart];
}
void data_py_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart] * s->particles->Momentum[1][ipart];
}
void data_pz_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart] * s->particles->Momentum[2][ipart];
}
void data_pressure_xx(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * pow(s->particles->Momentum[0][ipart],2)
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_pressure_yy(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * pow(s->particles->Momentum[1][ipart],2)
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_pressure_zz(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * pow(s->particles->Momentum[2][ipart],2)
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_pressure_xy(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[0][ipart]
                     * s->particles->Momentum[1][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_pressure_xz(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[0][ipart]
                     * s->particles->Momentum[2][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_pressure_yz(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[1][ipart]
                     * s->particles->Momentum[2][ipart]
                     / sqrt( 1. + pow(s->particles->Momentum[0][ipart],2)
                                + pow(s->particles->Momentum[1][ipart],2)
                                + pow(s->particles->Momentum[2][ipart],2) );
}
void data_ekin_vx_density(Species * s, vector<double> &array, unsigned int npart) {
    for (unsigned int ipart = 0 ; ipart < npart ; ipart++)
        array[ipart] = s->mass * s->particles->Weight[ipart]
                     * s->particles->Momentum[0][ipart]
                     * (1. - 1./sqrt(1. + pow(s->particles->Momentum[0][ipart],2)
                                        + pow(s->particles->Momentum[1][ipart],2)
                                        + pow(s->particles->Momentum[2][ipart],2)));
}


// Constructor
DiagnosticParticles::DiagnosticParticles( Params &params, SmileiMPI* smpi, Patch* patch, int diagId )
{
    fileId_ = 0;
    
    int n_diag_particles = diagId;
    
    ostringstream name("");
    name << "Diagnotic Particles #" << n_diag_particles;
    string errorPrefix = name.str();
    
    // get parameter "output" that determines the quantity to sum in the output array
    output = "";
    if (!PyTools::extract("output",output,"DiagParticles",n_diag_particles))
        ERROR(errorPrefix << ": parameter `output` required");
    
    if        (output == "density"        ) {
        data_filling = &data_density;
    } else if (output == "charge_density" ) {
        data_filling = &data_charge_density;
    } else if (output == "jx_density"     ) {
        data_filling = &data_jx_density;
    } else if (output == "jy_density"     ) {
        data_filling = &data_jy_density;
    } else if (output == "jz_density"     ) {
        data_filling = &data_jz_density;
    } else if (output == "ekin_density"   ) {
        data_filling = &data_ekin_density;
    } else if (output == "p_density"      ) {
        data_filling = &data_p_density;
    } else if (output == "px_density"     ) {
        data_filling = &data_px_density;
    } else if (output == "py_density"     ) {
        data_filling = &data_py_density;
    } else if (output == "pz_density"     ) {
        data_filling = &data_pz_density;
    } else if (output == "pressure_xx"    ) {
        data_filling = &data_pressure_xx;
    } else if (output == "pressure_yy"    ) {
        data_filling = &data_pressure_yy;
    } else if (output == "pressure_zz"    ) {
        data_filling = &data_pressure_zz;
    } else if (output == "pressure_xy"    ) {
        data_filling = &data_pressure_xy;
    } else if (output == "pressure_xz"    ) {
        data_filling = &data_pressure_xz;
    } else if (output == "pressure_yz"    ) {
        data_filling = &data_pressure_yz;
    } else if (output == "ekin_vx_density") {
        data_filling = &data_ekin_vx_density;
    } else {
        ERROR(errorPrefix << ": parameter `output = "<< output <<"` not understood");
    }

    // get parameter "every" which describes a timestep selection
    timeSelection = new TimeSelection(
        PyTools::extract_py("every", "DiagParticles", n_diag_particles),
        name.str()
    );
    
    // get parameter "flush_every" which describes a timestep selection for flushing the file
    flush_timeSelection = new TimeSelection(
        PyTools::extract_py("flush_every", "DiagParticles", n_diag_particles),
        name.str()
    );
    
    // get parameter "time_average" that determines the number of timestep to average the outputs
    time_average = 1;
    PyTools::extract("time_average",time_average,"DiagParticles",n_diag_particles);
    if ( time_average < 1 ) time_average=1;
    if ( time_average > timeSelection->smallestInterval() )
        ERROR(errorPrefix << ": `time_average` is incompatible with `every`");
    
    // get parameter "species" that determines the species to use (can be a list of species)
    vector<string> species_names;
    if (!PyTools::extract("species",species_names,"DiagParticles",n_diag_particles))
        ERROR(errorPrefix << ": parameter `species` required");
    // verify that the species exist, remove duplicates and sort by number
    species = params.FindSpecies(patch->vecSpecies, species_names);
    
    
    // get parameter "axes" that adds axes to the diagnostic
    // Each axis should contain several items:
    //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
    vector<PyObject*> allAxes=PyTools::extract_pyVec("axes","DiagParticles",n_diag_particles);
    
    if (allAxes.size() == 0)
        ERROR(errorPrefix << ": axes must contain something");
    
    // Loop axes and extract their format
    for (unsigned int iaxis=0; iaxis<allAxes.size(); iaxis++ ) {
        DiagnosticParticlesAxis tmpAxis;
        PyObject *oneAxis=allAxes[iaxis];
        
        // Axis must be a list
        if (!PyTuple_Check(oneAxis) && !PyList_Check(oneAxis))
            ERROR(errorPrefix << ": axis #" << iaxis << " must be a list");
        PyObject* seq = PySequence_Fast(oneAxis, "expected a sequence");
        
        // Axis must have 4 elements or more
        unsigned int lenAxisArgs=PySequence_Size(seq);
        if (lenAxisArgs<4)
            ERROR(errorPrefix << ": axis #" << iaxis << " contain at least 4 arguments");
        
        // Try to extract first element: type
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 0),tmpAxis.type))
            ERROR(errorPrefix << ", axis #" << iaxis << ": First item must be a string (axis type)");
        
        if        (tmpAxis.type == "x" ) {
            tmpAxis.axis_binning = &axis_x;
        } else if (tmpAxis.type == "y" ) {
            if (params.nDim_particle <2)
                ERROR(errorPrefix << ": axis y cannot exist in <2D");
            tmpAxis.axis_binning = &axis_y;
        } else if (tmpAxis.type == "z" ) {
            if (params.nDim_particle <3)
                ERROR(errorPrefix << ": axis z cannot exist in <3D");
            tmpAxis.axis_binning = &axis_z;
        } else if (tmpAxis.type == "px" ) {
            tmpAxis.axis_binning = &axis_px;
        } else if (tmpAxis.type == "py" ) {
            tmpAxis.axis_binning = &axis_py;
        } else if (tmpAxis.type == "pz" ) {
            tmpAxis.axis_binning = &axis_pz;
        } else if (tmpAxis.type == "p" ) {
            tmpAxis.axis_binning = &axis_p;
        } else if (tmpAxis.type == "gamma" ) {
            tmpAxis.axis_binning = &axis_gamma;
        } else if (tmpAxis.type == "ekin" ) {
            tmpAxis.axis_binning = &axis_ekin;
        } else if (tmpAxis.type == "vx" ) {
            tmpAxis.axis_binning = &axis_vx;
        } else if (tmpAxis.type == "vy" ) {
            tmpAxis.axis_binning = &axis_vy;
        } else if (tmpAxis.type == "vz" ) {
            tmpAxis.axis_binning = &axis_vz;
        } else if (tmpAxis.type == "v" ) {
            tmpAxis.axis_binning = &axis_v;
        } else if (tmpAxis.type == "vperp2" ) {
            tmpAxis.axis_binning = &axis_vperp2;
        } else if (tmpAxis.type == "charge" ) {
            tmpAxis.axis_binning = &axis_charge;
        } else if (tmpAxis.type == "chi" ) {
            // The requested species must be radiating
            for (unsigned int ispec=0 ; ispec < species.size() ; ispec++)
                if( ! patch->vecSpecies[species[ispec]]->particles->isRadReaction )
                    ERROR(errorPrefix << ": axis #" << iaxis << " 'chi' requires all species to be 'radiating'");
            tmpAxis.axis_binning = &axis_chi;
        } else if (tmpAxis.type == "composite") {
            ERROR(errorPrefix << ": axis type cannot be 'composite'");
        
        } else {
            // If not "usual" type, try to find composite type
            for( unsigned int i=1; i<=tmpAxis.type.length(); i++ )
                if( tmpAxis.type.substr(i,1) == " " )
                    ERROR(errorPrefix << ": axis #" << iaxis << " type cannot contain whitespace");
            if( tmpAxis.type.length()<2 )
                ERROR(errorPrefix << ": axis #" << iaxis << " type not understood");
            
            // Analyse character by character
            tmpAxis.coefficients.resize( params.nDim_particle , 0. );
            unsigned int previ=0;
            double sign=1.;
            tmpAxis.type += "+";
            for( unsigned int i=1; i<=tmpAxis.type.length(); i++ ) {
                // Split string at "+" or "-" location
                if( tmpAxis.type.substr(i,1) == "+" || tmpAxis.type.substr(i,1) == "-" ) {
                    // Get one segment of the split string
                    string segment = tmpAxis.type.substr(previ,i-previ);
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
                    if( tmpAxis.coefficients[direction_index] != 0. )
                        ERROR(errorPrefix << ": axis #" << iaxis << " type " << direction << " appears twice");
                    // Get the remaining characters, which should be a number
                    tmpAxis.coefficients[direction_index] = j>1 ? ::atof(segment.substr(0,j-1).c_str()) : 1.;
                    tmpAxis.coefficients[direction_index] *= sign;
                    // Save sign and position for next segment
                    sign = tmpAxis.type.substr(i,1) == "+" ? 1. : -1;
                    previ = i+1;
                }
            }
            
            tmpAxis.type = "composite:"+tmpAxis.type.substr(0,tmpAxis.type.length()-1);
            tmpAxis.axis_binning = &axis_composite;
        }
        
        // Try to extract second element: axis min
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 1),tmpAxis.min)) {
            ERROR(errorPrefix << ", axis #" << iaxis << ": Second item must be a double (axis min)");
        }
        
        // Try to extract third element: axis max
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 2),tmpAxis.max)) {
            ERROR(errorPrefix << ", axis #" << iaxis << ": Third item must be a double (axis max)");
        }
        
        // Try to extract fourth element: axis nbins
        if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 3),tmpAxis.nbins)) {
            ERROR(errorPrefix << ", axis #" << iaxis << ": Fourth item must be an int (number of bins)");
        }
        
        // 5 - Check for  other keywords such as "logscale" and "edge_inclusive"
        tmpAxis.logscale = false;
        tmpAxis.edge_inclusive = false;
        for(unsigned int i=4; i<lenAxisArgs; i++) {
            string my_str("");
            PyTools::convert(PySequence_Fast_GET_ITEM(seq, i),my_str);
            if(my_str=="logscale" ||  my_str=="log_scale" || my_str=="log")
                tmpAxis.logscale = true;
            else if(my_str=="edges" ||  my_str=="edge" ||  my_str=="edge_inclusive" ||  my_str=="edges_inclusive")
                tmpAxis.edge_inclusive = true;
            else
                ERROR(errorPrefix << ": keyword `" << my_str << "` not understood");
        }
        
        axes.push_back(tmpAxis);
        
        Py_DECREF(seq);
        
    }
    
    // calculate the total size of the output array
    output_size = 1;
    for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++)
        output_size *= axes[iaxis].nbins;
    
    // Output info on diagnostics
    if ( smpi->isMaster() ) {
        ostringstream mystream("");
        mystream.str("");
        mystream << species_names[0];
        for(unsigned int i=1; i<species_names.size(); i++)
            mystream << "," << species_names[i];
        MESSAGE(1,"Created particle diagnostic #" << n_diag_particles << ": species " << mystream.str());
        for(unsigned int i=0; i<axes.size(); i++) {
            mystream.str("");
            mystream << "Axis ";
            if( axes[i].type.substr(0,9) == "composite" ) {
                bool first = true;
                for( unsigned int idim=0; idim<axes[i].coefficients.size(); idim++ ) {
                    if( axes[i].coefficients[idim]==0. ) continue;
                    bool negative = axes[i].coefficients[idim]<0.;
                    double coeff = (negative?-1.:1.)*axes[i].coefficients[idim];
                    mystream << (negative?"-":(first?"":"+"));
                    if( coeff!=1. ) mystream << coeff;
                    mystream << (idim==0?"x":(idim==1?"y":"z"));
                    first = false;
                }
            } else {
                mystream << axes[i].type;
            }
            mystream << " from " << axes[i].min << " to " << axes[i].max << " in " << axes[i].nbins << " steps";
            if( axes[i].logscale       ) mystream << " [LOGSCALE] ";
            if( axes[i].edge_inclusive ) mystream << " [EDGE INCLUSIVE]";
            MESSAGE(2,mystream.str());
        }
        
        // init HDF files (by master, only if it doesn't yet exist)
        mystream.str(""); // clear
        mystream << "ParticleDiagnostic" << n_diag_particles << ".h5";
        filename = mystream.str();
    }

} // END DiagnosticParticles::DiagnosticParticles


DiagnosticParticles::~DiagnosticParticles()
{
    delete timeSelection;
    delete flush_timeSelection;
} // END DiagnosticParticles::~DiagnosticParticles


// Called only by patch master of process master
void DiagnosticParticles::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if (!smpi->isMaster()) return;
    
    if( fileId_>0 ) return;
    
    if ( newfile ) {
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        // write all parameters as HDF5 attributes
        H5::attr(fileId_, "Version", string(__VERSION));
        H5::attr(fileId_, "output" , output);
        H5::attr(fileId_, "time_average"  , time_average);
        // write all species
        ostringstream mystream("");
        mystream.str(""); // clear
        for (unsigned int i=0 ; i < species.size() ; i++)
            mystream << species[i] << " ";
        H5::attr(fileId_, "species", mystream.str());
        // write each axis
        for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
            mystream.str(""); // clear
            mystream << "axis" << iaxis;
            string str1 = mystream.str();
            mystream.str(""); // clear
            mystream << axes[iaxis].type << " " << axes[iaxis].min << " " << axes[iaxis].max << " "
                     << axes[iaxis].nbins << " " << axes[iaxis].logscale << " " << axes[iaxis].edge_inclusive << " [";
            for( unsigned int idim=0; idim<axes[iaxis].coefficients.size(); idim++) {
                mystream << axes[iaxis].coefficients[idim];
                if(idim<axes[iaxis].coefficients.size()-1) mystream << ",";
            }
            mystream << "]";
            string str2 = mystream.str();
            H5::attr(fileId_, str1, str2);
        }
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
    else {
        fileId_ = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
}


void DiagnosticParticles::closeFile()
{
    if (fileId_!=0) {
        H5Fclose(fileId_);
        fileId_ = 0;
    }

} // END closeFile


bool DiagnosticParticles::prepare( int timestep )
{
    // Get the previous timestep of the time selection
    int previousTime = timeSelection->previousTime(timestep);
    
    // Leave if the timestep is not the good one
    if (timestep - previousTime >= time_average) return false;
    
    // Allocate memory for the output array (already done if time-averaging)
    data_sum.resize(output_size);
    
    // if first time, erase output array
    if (timestep == previousTime)
        fill(data_sum.begin(), data_sum.end(), 0.);
    
    return true;
    
} // END prepare


// run one particle diagnostic
void DiagnosticParticles::run( Patch* patch, int timestep )
{
    
    vector<int> index_array;
    vector<double> axis_array, data_array;
    int axissize, ind;
    double axismin, axismax, coeff;
    unsigned int ipart, npart;
    
    // loop species
    for (unsigned int ispec=0 ; ispec < species.size() ; ispec++) {
        
        Species * s = patch->vecSpecies[species[ispec]];
        npart = s->particles->size();
        axis_array .resize(npart); // array to store particle axis data
        index_array.resize(npart); // array to store particle output index
        data_array .resize(npart); // array to store particle output data
        
        fill(index_array.begin(), index_array.end(), 0);
        
        // 1 - loop on the different axes requested and compute the output index of each particle
        // --------------------------------------------------------------------------------------
        for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
            
            axismin  = axes[iaxis].min;
            axismax  = axes[iaxis].max;
            axissize = axes[iaxis].nbins;
            
            // first loop on particles to store the indexing (axis) quantity
            (*(axes[iaxis].axis_binning))( s, axis_array, npart, axes[iaxis] );
            
            // Now, axis_array has the index of each particle
            // This index tells where, along the axis, its data should go
            
            // if log scale
            if (axes[iaxis].logscale) {
                // then loop again and convert to log
                for (ipart = 0 ; ipart < npart ; ipart++)
                    axis_array[ipart] = log10(axis_array[ipart]);
                // also convert other quantities
                axismin = log10(axismin);
                axismax = log10(axismax);
            }
            
            // The indexes are "reshaped" in one dimension.
            // For instance, in 3d, the index has the form  i = i3 + n3*( i2 + n2*i1 )
            // Here we do the multiplication by n3 or n2 (etc.)
            if (iaxis>0) {
                for (ipart = 0 ; ipart < npart ; ipart++)
                    index_array[ipart] *= axissize;
            }
            
            // loop again on the particles and calculate index
            // This is separated in two cases: edge_inclusive and edge_exclusive
            coeff = ((double)axissize)/(axismax-axismin);
            if (!axes[iaxis].edge_inclusive) { // if the particles out of the "box" must be excluded
                
                for (ipart = 0 ; ipart < npart ; ipart++) {
                    // skip already discarded particles
                    if (index_array[ipart] < 0) continue; 
                    
                    // calculate index
                    ind = floor( (axis_array[ipart]-axismin) * coeff );
                    
                    // index valid only if in the "box"
                    if (ind >= 0  &&  ind < axissize) index_array[ipart] += ind;
                    else index_array[ipart] = -1; // discard particle
                }
                
            } else { // if the particles out of the "box" must be included
                
                for (ipart = 0 ; ipart < npart ; ipart++) {
                    // skip already discarded particles
                    if (index_array[ipart] < 0) continue; 
                    
                    // calculate index
                    ind = floor( (axis_array[ipart]-axismin) * coeff );
                    
                    // move out-of-range indexes back into range
                    if (ind < 0) ind = 0;
                    if (ind >= axissize) ind = axissize-1;
                    index_array[ipart] += ind;
                }
                
            }
            
        } // loop axes
        
        // 2 - fill the data array with the particle data
        // -----------------------------------------------
        (*data_filling)( s, data_array, npart );
        
        // 3 - sum the data into the data_sum according to the indexes
        // ---------------------------------------------------------------
        for (ipart = 0 ; ipart < npart ; ipart++) {
            ind = index_array[ipart];
            if (ind<0) continue; // skip discarded particles
            #pragma omp atomic
            data_sum[ind] += data_array[ipart];
        }
        
    } // loop species
    
} // END run


// Now the data_sum has been filled
// if needed now, store result to hdf file
void DiagnosticParticles::write(int timestep, SmileiMPI* smpi)
{
    if ( !smpi->isMaster() ) return;
    
    if (timestep - timeSelection->previousTime() != time_average-1) return;
    
    double coeff;
    // if time_average, then we need to divide by the number of timesteps
    if (time_average > 1) {
        coeff = 1./((double)time_average);
        for (int i=0; i<output_size; i++)
            data_sum[i] *= coeff;
    }
    
    // make name of the array
    ostringstream mystream("");
    mystream.str("");
    mystream << "timestep" << setw(8) << setfill('0') << timestep;
    
    // write the array if it does not exist already
    if (! H5Lexists( fileId_, mystream.str().c_str(), H5P_DEFAULT ) ) {
        // Prepare array dimensions
        unsigned int naxes = axes.size();
        hsize_t dims[naxes];
        for( unsigned int iaxis=0; iaxis<naxes; iaxis++) dims[iaxis] = axes[iaxis].nbins;
        // Create file space
        hid_t sid = H5Screate_simple(naxes, &dims[0], NULL);
        hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
        // create dataset
        hid_t did = H5Dcreate(fileId_, mystream.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid, H5P_DEFAULT);
        // write vector in dataset
        H5Dwrite(did, H5T_NATIVE_DOUBLE, sid, sid, H5P_DEFAULT, &data_sum[0]);
        // close all
        H5Dclose(did);
        H5Pclose(pid);
        H5Sclose(sid);
    }
    
    if( flush_timeSelection->theTimeIsNow(timestep) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    
    // Clear the array
    clear();
    data_sum.resize(0);
} // END write


//! Clear the array
void DiagnosticParticles::clear() {
    data_sum.resize(0);
    vector<double>().swap( data_sum );
}
