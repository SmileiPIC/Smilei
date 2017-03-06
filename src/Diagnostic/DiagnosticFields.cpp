
#include <string>

#include "DiagnosticFields.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, int ndiag )
{
    fileId_ = 0;
    data_group_id = 0;
    tmp_dset_id = 0;
    diag_n = ndiag;
    timestep = params.timestep;
    
    filespace_firstwrite = 0;
    memspace_firstwrite = 0;
    filespace_reread = 0;
    memspace_reread = 0;
    
    // Extract the time_average parameter
    time_average = 1;
    PyTools::extract("time_average", time_average, "DiagFields", ndiag);
    if( time_average < 1 ) time_average = 1;
    time_average_inv = 1./((double)time_average);
    
    // Define the filename
    ostringstream fn("");
    fn << "Fields"<< ndiag <<".h5";
    filename = fn.str();
    
    // Extract the requested fields
    vector<string> fieldsToDump(0);
    PyTools::extract("fields", fieldsToDump, "DiagFields", ndiag);
    
    // List all fields that are requested
    std::vector<Field*> allFields (0);
    ostringstream ss("");
    fields_indexes.resize(0);
    fields_names  .resize(0);
    hasRhoJs = false;
    // Loop fields
    for( unsigned int i=0; i<vecPatches(0)->EMfields->allFields.size(); i++ ) {
        string field_name = vecPatches(0)->EMfields->allFields[i]->name;
        bool RhoJ = field_name.at(0)=='J' || field_name.at(0)=='R';
        bool species_field = (field_name.at(0)=='J' && field_name.length()>2) || (field_name.at(0)=='R' && field_name.length()>3);
        // If field in list of fields to dump, then add it
        if( hasField(field_name, fieldsToDump) ) {
            ss << field_name << " ";
            fields_indexes.push_back( i );
            fields_names  .push_back( field_name );
            if( RhoJ ) hasRhoJs = true;
            // If field specific to a species, then allocate it
            if( species_field ) {
                for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                    Field * field = vecPatches(ipatch)->EMfields->allFields[i];
                    if( field->data_ != NULL ) continue;
                    if     ( field_name.substr(0,2)=="Jx" ) field->allocateDims(0,false);
                    else if( field_name.substr(0,2)=="Jy" ) field->allocateDims(1,false);
                    else if( field_name.substr(0,2)=="Jz" ) field->allocateDims(2,false);
                    else if( field_name.substr(0,2)=="Rh" ) field->allocateDims();
                }
            }
        }
    }
    
    // Some output
    ostringstream p("");
    p << "(time average = " << time_average << ")";
    MESSAGE(1,"Diagnostic Fields #"<<ndiag<<" "<<(time_average>1?p.str():"")<<" :");
    MESSAGE(2, ss.str() );
    
    // Create new fields in each patch, for time-average storage
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->EMfields->allFields_avg.resize( diag_n+1 );
        if( time_average > 1 ) {
            for( unsigned int ifield=0; ifield<fields_names.size(); ifield++)
                vecPatches(ipatch)->EMfields->allFields_avg[diag_n].push_back(
                    vecPatches(ipatch)->EMfields->createField(fields_names[ifield])
                );
        }
    }
    
    // Extract the time selection
    timeSelection = new TimeSelection( PyTools::extract_py( "every", "DiagFields", ndiag ), "DiagFields" );
    
    // If the time selection contains intervals smaller than the time average, then error
    if( timeSelection->smallestInterval() < time_average )
        ERROR("Diagnostic Fields #"<<ndiag<<" has a time average too large compared to its time-selection interval ('every')");
    
    // Extract the flush time selection
    flush_timeSelection = new TimeSelection( PyTools::extract_py( "flush_every", "DiagFields", ndiag ), "DiagFields flush_every" );
    
    // Copy the total number of patches
    tot_number_of_patches = params.tot_number_of_patches;
    
    // Prepare the property list for HDF5 output
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);
    
    // Prepare openPMD attributes
    string xyz = "xyz";
    axisLabels      .resize( params.nDim_field );
    gridGlobalOffset.resize( params.nDim_field );
    gridOffset      .resize( params.nDim_field );
    gridSpacing = params.cell_length;
    for( unsigned int idim=0; idim<params.nDim_field; idim++ ) {
        axisLabels      [idim] = xyz.substr(idim, 1);
        gridGlobalOffset[idim] = 0.;
        gridOffset      [idim] = 0.;
    }
    unitDimension.resize( fields_names.size() );
    for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ ) {
        unitDimension[ifield].resize(7, 0.);
        string first_char = fields_names[ifield].substr(0,1);
        if        ( first_char == "E" ) {
            unitDimension[ifield][0] = 1.;
            unitDimension[ifield][1] = 1.;
            unitDimension[ifield][2] = -3.;
            unitDimension[ifield][3] = -1.;
        } else if ( first_char == "B" ) {
            unitDimension[ifield][1] = 1.;
            unitDimension[ifield][2] = -2.;
            unitDimension[ifield][3] = -1.;
        } else if ( first_char == "J" ) {
            unitDimension[ifield][0] = -2.;
            unitDimension[ifield][3] = 1.;
        } else if ( first_char == "R" ) {
            unitDimension[ifield][0] = -3.;
        } else {
            ERROR(" impossible field name ");
        }
    }
}


DiagnosticFields::~DiagnosticFields()
{
    H5Pclose( write_plist );
    
    delete timeSelection;
    delete flush_timeSelection;
}


bool DiagnosticFields::hasField(string field_name, vector<string> fieldsToDump)
{
    bool hasfield;
    if( fieldsToDump.size()==0 ) {
        hasfield = true;
    } else {
        hasfield = false;
        for( unsigned int j=0; j<fieldsToDump.size(); j++ ) {
            if( field_name == fieldsToDump[j] ) {
                hasfield = true;
                break;
            }
        }
    }
    return hasfield;
}

void DiagnosticFields::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if( fileId_>0 ) return;
    
    if ( newfile ) {
        // Create file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Attributes for openPMD
//        uint32_t extension = 1; // ED-PIC extension. Not supported yet
        uint32_t extension = 0;
        H5::attr( fileId_, "openPMDextension", extension, H5T_NATIVE_UINT32);
        H5::attr( fileId_, "openPMD", "1.0.0");
        H5::attr( fileId_, "basePath", "/data/%T/");
        H5::attr( fileId_, "meshesPath", "");
        H5::attr( fileId_, "particlesPath", "");
        H5::attr( fileId_, "software", "Smilei");
        H5::attr( fileId_, "softwareVersion", __VERSION);
        H5::attr( fileId_, "date", params.getLocalTime());
        H5::attr( fileId_, "iterationEncoding", "groupBased");
        H5::attr( fileId_, "iterationFormat", "/data/%T/");
        
        // Make main "data" group where everything will be stored (required by openPMD)
        data_group_id = H5::group( fileId_, "data" );
        
        // Prepare some attributes for openPMD compatibility
        fieldSolverParameters = "";
        if       ( params.maxwell_sol == "Yee" ) {
            fieldSolver = "Yee";
        } else if( params.maxwell_sol == "Lehe" ) {
            fieldSolver = "Lehe";
        } else {
            fieldSolver = "other";
            fieldSolverParameters = params.maxwell_sol;
        }
        fieldBoundary.resize(params.nDim_field * 2);
        fieldBoundaryParameters.resize(params.nDim_field * 2);
        em_bc(params.bc_em_type_x[0], fieldBoundary[0], fieldBoundaryParameters[0]);
        em_bc(params.bc_em_type_x[1], fieldBoundary[1], fieldBoundaryParameters[1]);
        if( params.nDim_field > 1 ) {
            em_bc(params.bc_em_type_y[0], fieldBoundary[2], fieldBoundaryParameters[2]);
            em_bc(params.bc_em_type_y[1], fieldBoundary[3], fieldBoundaryParameters[3]);
            if( params.nDim_field > 2 ) {
                em_bc(params.bc_em_type_z[0], fieldBoundary[4], fieldBoundaryParameters[4]);
                em_bc(params.bc_em_type_z[1], fieldBoundary[5], fieldBoundaryParameters[5]);
            }
        }
        particleBoundary.resize(params.nDim_field * 2, "");
        particleBoundaryParameters.resize(params.nDim_field * 2, "");
        currentSmoothing = "none";
        currentSmoothingParameters = "";
        if( params.currentFilter_int > 0 ) {
            currentSmoothing = "Binomial";
            ostringstream t("");
            t << "numPasses="<<params.currentFilter_int;
            currentSmoothingParameters = t.str();
        }
    }
    else {
        // Open the existing file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose(pid);
        data_group_id = H5Gopen( fileId_, "data", H5P_DEFAULT );
    }
}

void DiagnosticFields::closeFile()
{
    if ( filespace_firstwrite>0 ) H5Sclose( filespace_firstwrite );
    if ( memspace_firstwrite >0 ) H5Sclose( memspace_firstwrite );
    if ( filespace_reread    >0 ) H5Sclose( filespace_reread );
    if ( memspace_reread     >0 ) H5Sclose( memspace_reread );
    if ( filespace           >0 ) H5Sclose( filespace );
    if ( memspace            >0 ) H5Sclose( memspace );
    if ( tmp_dset_id         >0 ) H5Dclose( tmp_dset_id );
    
    if( data_group_id>0 ) H5Gclose( data_group_id );
    data_group_id = 0;
    if( fileId_>0 ) H5Fclose( fileId_ );
    fileId_ = 0;
}



void DiagnosticFields::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}

bool DiagnosticFields::prepare( int itime )
{
    
    // Leave if the iteration is not the good one
    if (itime - timeSelection->previousTime(itime) >= time_average) return false;
    
    return true;
}


void DiagnosticFields::run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime )
{
    // If time-averaging, increment the average
    if( time_average>1 ) {
        #pragma omp for schedule(static)
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            for( unsigned int ifield=0; ifield<fields_names.size(); ifield++) {
                vecPatches(ipatch)->EMfields->incrementAvgField(
                    vecPatches(ipatch)->EMfields->allFields[fields_indexes[ifield]], // instantaneous field
                    vecPatches(ipatch)->EMfields->allFields_avg[diag_n][ifield]      // averaged field
                );
            }
        }
    }
    
    // If is not writing iteration, leave
    if (itime - timeSelection->previousTime(itime) != time_average-1) return;
    
    #pragma omp master
    {
        // Calculate the structure of the file depending on 1D, 2D, ...
        refHindex = (unsigned int)(vecPatches.refHindex_);
        setFileSplitting( smpi, vecPatches );
        
        // Create group for this iteration
        ostringstream name_t;
        name_t.str("");
        name_t << setfill('0') << setw(10) << itime;
        status = H5Lexists(data_group_id, name_t.str().c_str(), H5P_DEFAULT);
        if( status==0 )
           iteration_group_id = H5::group(data_group_id, name_t.str().c_str());
        // Warning if file unreachable
        if( status < 0 ) WARNING("Fields diagnostics could not write");
        // Add openPMD attributes ( "basePath" )
        H5::attr( iteration_group_id, "time", (double)(itime*timestep));
        H5::attr( iteration_group_id, "dt", (double)timestep);
        H5::attr( iteration_group_id, "timeUnitSI", 0.); // not relevant
        // Add openPMD attributes ( "meshesPath" )
        H5::attr( iteration_group_id, "fieldSolver", fieldSolver);
        H5::attr( iteration_group_id, "fieldSolverParameters", fieldSolverParameters);
        H5::attr( iteration_group_id, "fieldBoundary", fieldBoundary);
        H5::attr( iteration_group_id, "fieldBoundaryParameters", fieldBoundaryParameters);
        H5::attr( iteration_group_id, "particleBoundary", particleBoundary);
        H5::attr( iteration_group_id, "particleBoundaryParameters", particleBoundaryParameters);
        H5::attr( iteration_group_id, "currentSmoothing", currentSmoothing);
        H5::attr( iteration_group_id, "currentSmoothingParameters", currentSmoothingParameters);
        H5::attr( iteration_group_id, "chargeCorrection", "none");
        H5::attr( iteration_group_id, "chargeCorrectionParameters", "");
        H5::attr( iteration_group_id, "fieldSmoothing", "none");
        H5::attr( iteration_group_id, "fieldSmoothingParameters", "");
    }
    #pragma omp barrier
    
    // Do not output diag if this iteration has already been written or if problem with file
    if( status != 0 ) return;
    
    unsigned int nPatches( vecPatches.size() );
    
    // For each field, combine all patches and write out
    for( unsigned int ifield=0; ifield < fields_indexes.size(); ifield++ ) {
        
        // Copy the patch field to the buffer
        #pragma omp barrier
        #pragma omp for schedule(static)
        for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++)
            getField( vecPatches(ipatch), ifield );
        
        #pragma omp master
        {
            // Create field dataset in HDF5
            hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
            hid_t dset_id  = H5Dcreate( iteration_group_id, fields_names[ifield].c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
            H5Pclose(plist_id);
            
            // Write
            writeField(dset_id, itime);
            
            // Attributes for openPMD
            H5::attr(dset_id, "geometry", "cartesian");
            H5::attr(dset_id, "dataOrder", "C");
            H5::attr(dset_id, "axisLabels", axisLabels);
            H5::attr(dset_id, "gridSpacing", gridSpacing);
            H5::attr(dset_id, "gridGlobalOffset", gridGlobalOffset);
            H5::attr(dset_id, "gridOffset", gridOffset);
            H5::attr(dset_id, "gridUnitSI", 0.);      
            H5::attr(dset_id, "unitSI", 0.);      
            H5::attr(dset_id, "unitDimension", unitDimension[ifield]);
            H5::attr(dset_id, "timeOffset", 0.);      
            
            // Close dataset
            H5Dclose( dset_id );
        }
    }
    
    #pragma omp master
    {
        H5Gclose(iteration_group_id);
        if( tmp_dset_id>0 ) H5Dclose( tmp_dset_id );
        tmp_dset_id=0;
        if( flush_timeSelection->theTimeIsNow(itime) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
}

void DiagnosticFields::em_bc(string bc, string& bc_0, string& bc_1)
{
    if( bc == "periodic" ) {
        bc_0 = "periodic";
        bc_1 = "periodic";
    } else if( bc == "reflective" ) {
        bc_0 = "reflecting";
        bc_1 = "reflecting";
    } else if( bc == "silver-muller" ) {
        bc_0 = "open";
        bc_1 = "silver-muller";
    } else {
        ERROR(" impossible boundary condition ");
    }
}

bool DiagnosticFields::needsRhoJs(int itime) {
    return hasRhoJs && timeSelection->theTimeIsNow(itime);
}
