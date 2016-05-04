
#include "DiagnosticParticles.h"

#include <iomanip>

using namespace std;


DiagnosticParticles::DiagnosticParticles( Params &params, SmileiMPI* smpi, Patch* patch, int diagId )
{
    fileId_ = 0;
    
    int n_diag_particles = diagId;

    // n_diag_particles ...

    std::vector<Species*>& vecSpecies = patch->vecSpecies;

    ostringstream name("");
    name << "Diagnotic Particles #" << n_diag_particles;
    string errorPrefix = name.str();
    
    // get parameter "output" that determines the quantity to sum in the output array
    output = "";
    if (!PyTools::extract("output",output,"DiagParticles",n_diag_particles))
        ERROR(errorPrefix << ": parameter `output` required");
    
    // get parameter "every" which describes a timestep selection
    timeSelection = new TimeSelection(
        PyTools::extract_py("every", "DiagParticles", n_diag_particles),
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
    species = params.FindSpecies(vecSpecies, species_names);
    
    
    // get parameter "axes" that adds axes to the diagnostic
    // Each axis should contain several items:
    //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
    vector<PyObject*> allAxes=PyTools::extract_pyVec("axes","DiagParticles",n_diag_particles);
    
    if (allAxes.size() == 0)
        ERROR(errorPrefix << ": axes must contain something");
    
    for (unsigned int iaxis=0; iaxis<allAxes.size(); iaxis++ ) {
        DiagnosticParticlesAxis tmpAxis;
        PyObject *oneAxis=allAxes[iaxis];
        if (PyTuple_Check(oneAxis) || PyList_Check(oneAxis)) {
            PyObject* seq = PySequence_Fast(oneAxis, "expected a sequence");
            unsigned int lenAxisArgs=PySequence_Size(seq);
            if (lenAxisArgs<4)
                ERROR(errorPrefix << ": axis #" << iaxis << " contain at least 4 arguments");
            
            if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 0),tmpAxis.type)) {
                ERROR(errorPrefix << ", axis #" << iaxis << ": First item must be a string (axis type)");
            } else {
                if (   (tmpAxis.type == "z" && params.nDim_particle <3)
                    || (tmpAxis.type == "y" && params.nDim_particle <2) )
                    ERROR(errorPrefix << ": axis " << tmpAxis.type << " cannot exist in " << params.nDim_particle << "D");
            }
            
            if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 1),tmpAxis.min)) {
                ERROR(errorPrefix << ", axis #" << iaxis << ": Second item must be a double (axis min)");
            }
            
            if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 2),tmpAxis.max)) {
                ERROR(errorPrefix << ", axis #" << iaxis << ": Third item must be a double (axis max)");
            }
            
            
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
            mystream << "Axis " << axes[i].type << " from " << axes[i].min << " to " << axes[i].max << " in " << axes[i].nbins << " steps";
            if( axes[i].logscale       ) mystream << " [LOGSCALE] ";
            if( axes[i].edge_inclusive ) mystream << " [EDGE INCLUSIVE]";
            MESSAGE(2,mystream.str());
        }

        // init HDF files (by master, only if it doesn't yet exist)
        mystream.str(""); // clear
        mystream << "ParticleDiagnostic" << n_diag_particles << ".h5";
        filename = mystream.str();
    }

    type_ = "Particles";

} // END DiagnosticParticles::DiagnosticParticles


// Cloning constructor // NOT USED
DiagnosticParticles::DiagnosticParticles( DiagnosticParticles* diag)
{
    
    output        = diag->output;
    time_average  = diag->time_average;
    species       = diag->species ;
    timeSelection = new TimeSelection(diag->timeSelection);
    
    for (unsigned int iaxis=0; iaxis<diag->axes.size(); iaxis++ ) {
        DiagnosticParticlesAxis tmpAxis;
        tmpAxis.type           = diag->axes[iaxis].type;
        tmpAxis.min            = diag->axes[iaxis].min;
        tmpAxis.max            = diag->axes[iaxis].max;
        tmpAxis.nbins          = diag->axes[iaxis].nbins;
        tmpAxis.logscale       = diag->axes[iaxis].logscale;
        tmpAxis.edge_inclusive = diag->axes[iaxis].edge_inclusive;
        axes.push_back(tmpAxis);
    }
    
    output_size = diag->output_size;
    
    type_ = "Particles";
}



DiagnosticParticles::~DiagnosticParticles()
{
    delete timeSelection;
    
} // END DiagnosticParticles::~DiagnosticParticles


// Called only by patch master of process master
void DiagnosticParticles::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if (!smpi->isMaster()) return;

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
                     << axes[iaxis].nbins << " " << axes[iaxis].logscale << " " << axes[iaxis].edge_inclusive;
            string str2 = mystream.str();
            H5::attr(fileId_, str1, str2);
        }
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

    std::vector<Species*>& vecSpecies = patch->vecSpecies;

    Species *s;
    Particles *p;
    vector<int> index_array;
    vector<double> *x, *y, *z, *px, *py, *pz, *w, *chi=NULL, axis_array, data_array;
    vector<short> *q;
    int nbins = vecSpecies[0]->bmin.size(); // number of bins in the particles binning (openMP)
    int bmin, bmax, axissize, ind;
    double axismin, axismax, mass, coeff;
    string axistype;
    ostringstream mystream("");
    
    // loop species
    for (unsigned int ispec=0 ; ispec < species.size() ; ispec++) {
        
        // make shortcuts
        s    = vecSpecies[species[ispec]];  // current species
        p    = (s->particles);              // current particles array
        x    = &(p->Position[0]);           // -+
        y    = &(p->Position[1]);           //  |- position
        z    = &(p->Position[2]);           // -+
        px   = &(p->Momentum[0]);           // -+
        py   = &(p->Momentum[1]);           //  |- momentum
        pz   = &(p->Momentum[2]);           // -+
        q    = &(p->Charge);                // charge
        w    = &(p->Weight);                // weight
        if (s->dynamics_type == "rrll")
            chi  = &(p->Chi);               // chi (for rad reaction particles particles)
        mass = s->mass;       // mass
        
        axis_array .resize(p->size()); // array to store particle axis data
        index_array.resize(p->size()); // array to store particle output index
        data_array .resize(p->size()); // array to store particle output data
        
        fill(index_array.begin(), index_array.end(), 0);
        
        // loop each openMP bin
        //! \todo Make OpenMP parallelization
        for (int ibin=0 ; ibin<nbins ; ibin++) {
            
            bmin = s->bmin[ibin];
            bmax = s->bmax[ibin];
            
            // 1 - loop on the different axes requested and compute the output index of each particle
            // --------------------------------------------------------------------------------------
            for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
                
                axismin  = axes[iaxis].min;
                axismax  = axes[iaxis].max;
                axissize = axes[iaxis].nbins;
                axistype = axes[iaxis].type;
                
                // first loop on particles to store the indexing (axis) quantity
                
                if      (axistype == "x"     )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*x)[ipart];
                
                else if (axistype == "y"     )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*y)[ipart];
                
                else if (axistype == "z"     )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*z)[ipart];
                
                else if (axistype == "px"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (*px)[ipart];
                
                else if (axistype == "py"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (*py)[ipart];
                
                else if (axistype == "pz"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (*pz)[ipart];
                
                else if (axistype == "p"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * sqrt(pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2));
                
                else if (axistype == "gamma" )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "ekin" )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) ) - 1.);
                
                else if (axistype == "vx"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*px)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "vy"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*py)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "vz"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*pz)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                        
                else if (axistype == "v"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = pow( 1. + 1./(pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2)) , -0.5);
                
                else if (axistype == "vperp2"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (pow((*py)[ipart],2) + pow((*pz)[ipart],2)) / (1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "charge")
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (double) (*q)[ipart];
                
                else if (axistype == "chi")
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (double) (*chi)[ipart];
                
                else   ERROR("In particle diagnostics, axis `" << axistype << "` unknown");
                
                // Now, axis_array points to the array that contains the particles data (for indexing)
                
                // if log scale
                if (axes[iaxis].logscale) {
                    // then loop again and convert to log
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = log10(axis_array[ipart]);
                    // also convert other quantities
                    axismin = log10(axismin);
                    axismax = log10(axismax);
                }
                
                // The indexes are "reshaped" in one dimension.
                // For instance, in 3d, the index has the form  i = i3 + n3*( i2 + n2*i1 )
                // Here we do the multiplication by n3 or n2 (etc.)
                if (iaxis>0) {
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        index_array[ipart] *= axissize;
                }
                
                // loop again on the particles and calculate index
                // This is separated in two cases: edge_inclusive and edge_exclusive
                coeff = ((double)axissize)/(axismax-axismin);
                if (!axes[iaxis].edge_inclusive) { // if the particles out of the "box" must be excluded
                    
                    for (int ipart = bmin ; ipart < bmax ; ipart++) {
                        // skip already discarded particles
                        if (index_array[ipart] < 0) continue; 
                        
                        // calculate index
                        ind = floor( (axis_array[ipart]-axismin) * coeff );
                        
                        // index valid only if in the "box"
                        if (ind >= 0  &&  ind < axissize) index_array[ipart] += ind;
                        else index_array[ipart] = -1; // discard particle
                    }
                    
                } else { // if the particles out of the "box" must be included
                    
                    for (int ipart = bmin ; ipart < bmax ; ipart++) {
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
            
            // 2 - prepare the data to output
            // ------------------------------
            if      (output == "density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart];
            
            else if (output == "charge_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]);
            
            else if (output == "jx_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]) * (*px)[ipart] / sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "jy_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]) * (*py)[ipart] / sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "jz_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]) * (*pz)[ipart] / sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "ekin_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (sqrt(1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart]) - 1.);
            
            else if (output == "p_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * sqrt((*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart]);
            
            else if (output == "px_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*px)[ipart];
            
            else if (output == "py_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*py)[ipart];
            
            else if (output == "pz_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*pz)[ipart];
            
            else if (output == "pressure_xx")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*px)[ipart]*(*px)[ipart]/ sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "pressure_yy")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*py)[ipart]*(*py)[ipart]/ sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "pressure_zz")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*pz)[ipart]*(*pz)[ipart]/ sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "pressure_xy")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*px)[ipart]*(*py)[ipart]/ sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "pressure_xz")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*px)[ipart]*(*pz)[ipart]/ sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            else if (output == "pressure_yz")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*py)[ipart]*(*pz)[ipart]/ sqrt( 1. + (*px)[ipart]*(*px)[ipart] + (*py)[ipart]*(*py)[ipart] + (*pz)[ipart]*(*pz)[ipart] );
            
            // 3 - sum the data into the data_sum according to the indexes
            // ---------------------------------------------------------------
            for (int ipart = bmin ; ipart < bmax ; ipart++) {
                ind = index_array[ipart];
                if (ind<0) continue; // skip discarded particles
                data_sum[ind] += data_array[ipart];
            }
            
        } // loop bins
        
    } // loop species

} // END run


// Now the data_sum has been filled
// if needed now, store result to hdf file
// called by MPI master only, when time-average has finished
void DiagnosticParticles::write(int timestep)
{
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
    // write the array
    if (! H5Lexists( fileId_, mystream.str().c_str(), H5P_DEFAULT ) )
        H5::vect(fileId_, mystream.str(), data_sum);
    else
        WARNING("DIAG PARTICLES COULD NOT WRITE");
    
    // Clear the array
    data_sum.resize(0); 

} // END write
