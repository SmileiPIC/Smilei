
#include "DiagnosticScreen.h"

#include <iomanip>

using namespace std;


// Constructor
DiagnosticScreen::DiagnosticScreen( Params &params, SmileiMPI* smpi, Patch* patch, int diagId )
{
    fileId_ = 0;
    
    screen_id = diagId;
    dt = params.timestep;
    
    ostringstream name("");
    name << "Diagnotic Screen #" << screen_id;
    string errorPrefix = name.str();
    
    // get parameter "shape" that determines if the screen is a plane or a sphere
    if (!PyTools::extract("shape",screen_shape,"DiagScreen",screen_id))
        ERROR(errorPrefix << ": parameter `shape` required");
    if( screen_shape == "plane" )
        screen_type = 0;
    else if( screen_shape != "sphere" )
        screen_type = 1;
    else
        ERROR(errorPrefix << ": parameter `shape` must be 'plane' or 'sphere'");
    
    // get parameter "point" that determines the plane reference point, or the sphere center
    if (!PyTools::extract("point",screen_point,"DiagScreen",screen_id))
        ERROR(errorPrefix << ": parameter `point` required");
    if( screen_point.size() != params.nDim_particle )
        ERROR(errorPrefix << ": parameter `point` must have "<<params.nDim_particle<<" elements in "<<params.nDim_particle<<"D");
    
    // get parameter "vector" that determines the plane normal, or the sphere radius
    if (!PyTools::extract("vector",screen_vector,"DiagScreen",screen_id)) {
        if( params.nDim_particle == 1 ) {
            screen_vector.resize( 1, 1. );
        } else {
            ERROR(errorPrefix << ": parameter `vector` required");
        }
    }
    if( screen_vector.size() != params.nDim_particle )
        ERROR(errorPrefix << ": parameter `vector` must have "<<params.nDim_particle<<" elements in "<<params.nDim_particle<<"D");
    // Calculate the unit vector
    screen_unitvector.resize( params.nDim_particle, 0. );
    screen_vectornorm = 0.;
    for( unsigned int i=0; i<params.nDim_particle; i++ ) screen_vectornorm += screen_vector[i]*screen_vector[i];
    if( screen_vectornorm == 0. )
        ERROR(errorPrefix << ": parameter `vector` must not be a null vector");
    screen_vectornorm = sqrt(screen_vectornorm);
    for( unsigned int i=0; i<params.nDim_particle; i++ ) screen_unitvector[i] = screen_vector[i]/screen_vectornorm;
    // Calculate other unit vectors to form an orthogonal base
    screen_vector_a.resize( params.nDim_particle, 0. );
    screen_vector_b.resize( params.nDim_particle, 0. );
    if( params.nDim_particle > 1 ) {
        screen_vector_a[0] = -screen_unitvector[1];
        screen_vector_a[1] =  screen_unitvector[0];
        double norm = sqrt(pow(screen_vector_a[0],2) + pow(screen_vector_a[1],2));
        if( norm < 1.e-8 ) {
            screen_vector_a[0] = 0.;
            screen_vector_a[1] = 1.;
        } else {
            screen_vector_a[0] /= norm;
            screen_vector_a[1] /= norm;
        }
    }
    if( params.nDim_particle > 2 ) {
        screen_vector_b[0] = screen_unitvector[1]*screen_vector_a[2] - screen_unitvector[2]*screen_vector_a[1];
        screen_vector_b[1] = screen_unitvector[2]*screen_vector_a[0] - screen_unitvector[0]*screen_vector_a[2];
        screen_vector_b[2] = screen_unitvector[0]*screen_vector_a[1] - screen_unitvector[1]*screen_vector_a[0];
    }
    
    // get parameter "output" that determines the quantity to sum in the output array
    if (!PyTools::extract("output",output,"DiagScreen",screen_id))
        ERROR(errorPrefix << ": parameter `output` required");
    
    // get parameter "every" which describes a timestep selection
    timeSelection = new TimeSelection(
        PyTools::extract_py("every", "DiagScreen", screen_id),
        name.str()
    );
    
    // get parameter "flush_every" which describes a timestep selection for flushing the file
    flush_timeSelection = new TimeSelection(
        PyTools::extract_py("flush_every", "DiagScreen", screen_id),
        name.str()
    );
    
    // get parameter "species" that determines the species to use (can be a list of species)
    vector<string> species_names;
    if (!PyTools::extract("species",species_names,"DiagScreen",screen_id))
        ERROR(errorPrefix << ": parameter `species` required");
    // verify that the species exist, remove duplicates and sort by number
    species = params.FindSpecies(patch->vecSpecies, species_names);
    
    // get parameter "axes" that adds axes to the diagnostic
    // Each axis should contain several items:
    //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
    vector<PyObject*> pyAxes=PyTools::extract_pyVec("axes","DiagScreen",screen_id);
    
    vector<string> excluded_types(0);
    if( screen_shape == "plane" ) {
        excluded_types.push_back( "theta_yx" );
        excluded_types.push_back( "theta_zx" );
    } else {
        excluded_types.push_back( "a" );
        excluded_types.push_back( "b" );
    }
    
    // Create the Histogram object
    if        (output == "density"        ) {
        histogram = new Histogram_density        ();
    } else if (output == "charge_density" ) {
        histogram = new Histogram_charge_density ();
    } else if (output == "jx_density"     ) {
        histogram = new Histogram_jx_density     ();
    } else if (output == "jy_density"     ) {
        histogram = new Histogram_jy_density     ();
    } else if (output == "jz_density"     ) {
        histogram = new Histogram_jz_density     ();
    } else if (output == "ekin_density"   ) {
        histogram = new Histogram_ekin_density   ();
    } else if (output == "p_density"      ) {
        histogram = new Histogram_p_density      ();
    } else if (output == "px_density"     ) {
        histogram = new Histogram_px_density     ();
    } else if (output == "py_density"     ) {
        histogram = new Histogram_py_density     ();
    } else if (output == "pz_density"     ) {
        histogram = new Histogram_pz_density     ();
    } else if (output == "pressure_xx"    ) {
        histogram = new Histogram_pressure_xx    ();
    } else if (output == "pressure_yy"    ) {
        histogram = new Histogram_pressure_yy    ();
    } else if (output == "pressure_zz"    ) {
        histogram = new Histogram_pressure_zz    ();
    } else if (output == "pressure_xy"    ) {
        histogram = new Histogram_pressure_xy    ();
    } else if (output == "pressure_xz"    ) {
        histogram = new Histogram_pressure_xz    ();
    } else if (output == "pressure_yz"    ) {
        histogram = new Histogram_pressure_yz    ();
    } else if (output == "ekin_vx_density") {
        histogram = new Histogram_ekin_vx_density();
    } else {
        ERROR(errorPrefix << ": parameter `output = "<< output <<"` not understood");
    }
    
    histogram->init(params, pyAxes, species, errorPrefix, patch, excluded_types);
    
    // If axes are "a", "b", "theta" or "phi", they need some coefficients
    unsigned int idim;
    for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
        string type = histogram->axes[i]->type;
        if( type=="a" || type=="b" || type=="theta" || type=="phi" ) {
            vector<double> coefficients(params.nDim_particle * 2);
            for( idim=0; idim<params.nDim_particle; idim++ )
                coefficients[idim] = screen_point[idim];
            if( type == "a" )
                for( idim=0; idim<params.nDim_particle; idim++ )
                    coefficients[params.nDim_particle+idim] = screen_vector_a[idim];
            else if( type == "b" )
                for( idim=0; idim<params.nDim_particle; idim++ )
                    coefficients[params.nDim_particle+idim] = screen_vector_b[idim];
            else if( type == "theta" )
                for( idim=0; idim<params.nDim_particle; idim++ )
                    coefficients[params.nDim_particle+idim] = screen_vector[idim] / pow(screen_vectornorm,2);
            else if( type == "phi" ) {
                coefficients.resize(params.nDim_particle * 3);
                for( idim=0; idim<params.nDim_particle; idim++ )
                    coefficients[params.nDim_particle+idim] = screen_vector_a[idim];
                for( idim=0; idim<params.nDim_particle; idim++ )
                    coefficients[2*params.nDim_particle+idim] = screen_vector_b[idim];
            }
            histogram->axes[i]->coefficients = coefficients;
        }
    }
    
    // Calculate the size of the output array
    output_size = 1;
    for( unsigned int i=0; i<histogram->axes.size(); i++ )
        output_size *= histogram->axes[i]->nbins;
    data_sum.resize( output_size, 0. );
    
    // Output info on diagnostics
    if ( smpi->isMaster() ) {
        ostringstream mystream("");
        mystream.str("");
        mystream << species_names[0];
        for(unsigned int i=1; i<species_names.size(); i++)
            mystream << "," << species_names[i];
        MESSAGE(1,"Created screen diagnostic #" << screen_id << ": species " << mystream.str());
        for(unsigned int i=0; i<histogram->axes.size(); i++) {
            HistogramAxis * axis = histogram->axes[i];
            mystream.str("");
            mystream << "Axis ";
            if( axis->type.substr(0,9) == "composite" ) {
                bool first = true;
                for( unsigned int idim=0; idim<axis->coefficients.size(); idim++ ) {
                    if( axis->coefficients[idim]==0. ) continue;
                    bool negative = axis->coefficients[idim]<0.;
                    double coeff = (negative?-1.:1.)*axis->coefficients[idim];
                    mystream << (negative?"-":(first?"":"+"));
                    if( coeff!=1. ) mystream << coeff;
                    mystream << (idim==0?"x":(idim==1?"y":"z"));
                    first = false;
                }
            } else {
                mystream << axis->type;
            }
            mystream << " from " << axis->min << " to " << axis->max << " in " << axis->nbins << " steps";
            if( axis->logscale       ) mystream << " [LOGSCALE] ";
            if( axis->edge_inclusive ) mystream << " [EDGE INCLUSIVE]";
            MESSAGE(2,mystream.str());
        }
        
        // init HDF files (by master, only if it doesn't yet exist)
        mystream.str(""); // clear
        mystream << "Screen" << screen_id << ".h5";
        filename = mystream.str();
    }

} // END DiagnosticScreen::DiagnosticScreen


DiagnosticScreen::~DiagnosticScreen()
{
    delete timeSelection;
    delete flush_timeSelection;
} // END DiagnosticScreen::~DiagnosticScreen


// Called only by patch master of process master
void DiagnosticScreen::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if (!smpi->isMaster()) return;
    
    if( fileId_>0 ) return;
    
    if ( newfile ) {
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        // write all parameters as HDF5 attributes
        H5::attr(fileId_, "Version", string(__VERSION));
        H5::attr(fileId_, "output" , output);
        // write all species
        ostringstream mystream("");
        mystream.str(""); // clear
        for (unsigned int i=0 ; i < species.size() ; i++)
            mystream << species[i] << " ";
        H5::attr(fileId_, "species", mystream.str());
        // write each axis
        for (unsigned int iaxis=0 ; iaxis < histogram->axes.size() ; iaxis++) {
            mystream.str(""); // clear
            mystream << "axis" << iaxis;
            string str1 = mystream.str();
            mystream.str(""); // clear
            mystream << histogram->axes[iaxis]->type << " " << histogram->axes[iaxis]->min << " " << histogram->axes[iaxis]->max << " "
                     << histogram->axes[iaxis]->nbins << " " << histogram->axes[iaxis]->logscale << " " << histogram->axes[iaxis]->edge_inclusive << " [";
            for( unsigned int idim=0; idim<histogram->axes[iaxis]->coefficients.size(); idim++) {
                mystream << histogram->axes[iaxis]->coefficients[idim];
                if(idim<histogram->axes[iaxis]->coefficients.size()-1) mystream << ",";
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


void DiagnosticScreen::closeFile()
{
    if (fileId_!=0) {
        H5Fclose(fileId_);
        fileId_ = 0;
    }

} // END closeFile


bool DiagnosticScreen::prepare( int timestep )
{
    
    // This diag always runs, but the output is not done at every timestep
    return true;
    
} // END prepare


// run one screen diagnostic
void DiagnosticScreen::run( Patch* patch, int timestep )
{
    
    vector<int> int_buffer;
    vector<double> double_buffer;
    unsigned int npart, ndim = screen_point.size(), ipart, idim;
    double side, side_old, dtg;
    
    // Verify that this patch is in a useful region for this diag
    if( screen_type == 0 ) { // plane
        double distance_to_plane = 0.;
        for( idim=0; idim<ndim; idim++ )
            distance_to_plane += (patch->center[idim] - screen_point[idim]) * screen_unitvector[idim];
        if( abs(distance_to_plane) > patch->radius ) return;
    } else { // sphere
        double distance_to_center = 0.;
        for( idim=0; idim<ndim; idim++ )
            distance_to_center += pow(patch->center[idim] - screen_point[idim], 2);
        distance_to_center = sqrt(distance_to_center);
        if( abs(screen_vectornorm - distance_to_center) > patch->radius ) return;
    }
    
    // loop species
    for (unsigned int ispec=0 ; ispec < species.size() ; ispec++) {
        
        Species * s = patch->vecSpecies[species[ispec]];
        npart = s->particles->size();
        int_buffer   .resize(npart);
        double_buffer.resize(npart);
        
        // Fill the int_buffer with -1 (not crossing screen) and 0 (crossing screen)
        if( screen_type == 0 ) { // plane
            for( ipart=0; ipart<npart; ipart++ ) {
                side = 0.;
                side_old = 0.;
                dtg = dt / sqrt( 1. + pow(s->particles->Momentum[0][ipart], 2)
                                    + pow(s->particles->Momentum[1][ipart], 2)
                                    + pow(s->particles->Momentum[2][ipart], 2) );
                for( idim=0; idim<ndim; idim++ ) {
                    side += (s->particles->Position[idim][ipart] - screen_point[idim]) * screen_unitvector[idim];
                    side_old += (s->particles->Position[idim][ipart] - dtg*(s->particles->Momentum[idim][ipart]) - screen_point[idim]) * screen_unitvector[idim];
                }
                if( side*side_old < 0. )
                    int_buffer[ipart] = 0;
                else
                    int_buffer[ipart] = -1;
            }
        } else { // sphere
            for( ipart=0; ipart<npart; ipart++ ) {
                side = 0.;
                side_old = 0.;
                dtg = dt / sqrt( 1. + pow(s->particles->Momentum[0][ipart], 2)
                                    + pow(s->particles->Momentum[1][ipart], 2)
                                    + pow(s->particles->Momentum[2][ipart], 2) );
                for( idim=0; idim<ndim; idim++ ) {
                    side += pow(s->particles->Position[idim][ipart] - screen_point[idim], 2);
                    side_old += pow(s->particles->Position[idim][ipart] - dtg*(s->particles->Momentum[idim][ipart]) - screen_point[idim], 2);
                }
                if( (screen_vectornorm-sqrt(side))*(screen_vectornorm-sqrt(side_old)) < 0. )
                    int_buffer[ipart] = 0;
                else
                    int_buffer[ipart] = -1;
            }
        }
        
        // digitize
        histogram->digitize( s, double_buffer, int_buffer, data_sum );
        
    }
    
} // END run


// if needed now, store result to hdf file
void DiagnosticScreen::write(int timestep, SmileiMPI* smpi)
{
    if ( !smpi->isMaster() ) return;
    
    if (!timeSelection->theTimeIsNow(timestep) ) return;
    
    // make name of the array
    ostringstream mystream("");
    mystream.str("");
    mystream << "timestep" << setw(8) << setfill('0') << timestep;
    
    // write the array if it does not exist already
    if (! H5Lexists( fileId_, mystream.str().c_str(), H5P_DEFAULT ) ) {
        // Prepare array dimensions
        unsigned int naxes = histogram->axes.size();
        hsize_t dims[naxes];
        for( unsigned int iaxis=0; iaxis<naxes; iaxis++) dims[iaxis] = histogram->axes[iaxis]->nbins;
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
} // END write


//! Zero the array
void DiagnosticScreen::clear() {
    fill( data_sum.begin(), data_sum.end(), 0. );
}
