
#include <string>

#include "DiagnosticFields.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI* smpi, Patch* patch, int ndiag )
{
    fileId_ = 0;
    tmp_dset_id = 0;
    
    // Extract the time_average parameter
    time_average = 1;
    PyTools::extract("time_average", time_average, "DiagFields", ndiag);
    if( time_average < 1 )
        time_average = 1;
    
    // Verify that only one diag of this type exists
    int tavg;
    for( int idiag=0; idiag<ndiag; idiag++ ) {
        tavg = 1;
        PyTools::extract("time_average", tavg, "DiagFields", idiag);
        if( tavg < 1 )
            tavg = 1;
        if( tavg*time_average == 1 || (tavg>1 && time_average>1) )
            ERROR("Cannot have two DiagFields with time_average "<<(tavg==1?"=":">")<<" 1");
    }
    
    // Define the filename and get the list of fields
    std::vector<Field*> * allFields;
    if ( time_average==1 ) {
        filename = "Fields.h5";
        allFields = &(patch->EMfields->allFields);
    }
    else {
        filename = "Fields_avg.h5";
        allFields = &(patch->EMfields->allFields_avg);
    }
    
    // Extract the requested fields
    std::vector<std::string> fieldsToDump(0);
    PyTools::extract("fields", fieldsToDump, "DiagFields", ndiag);
    
    // List all fields that are requested
    ostringstream ss("");
    fields_indexes.resize(0);
    fields_names  .resize(0);
    bool hasfield;
    for( unsigned int i=0; i<allFields->size(); i++ ) {
        string field_name = (*allFields)[i]->name;
        if( field_name.find("_avg") < string::npos ) field_name.erase(field_name.find("_avg"));
        
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
        
        if( hasfield ) {
            ss << field_name << " ";
            fields_indexes.push_back( i );
            fields_names  .push_back( field_name );
        }
    }
    MESSAGE(1,"EM fields dump "<<(time_average>1?"(avg)":"     ")<<" :");
    MESSAGE(2, ss.str() );
    
    // Extract the time selection
    timeSelection = new TimeSelection( PyTools::extract_py( "every", "DiagFields", ndiag ), "DiagFields" );
    
    // Copy the total number of patches
    tot_number_of_patches = params.tot_number_of_patches;
    
    // Prepare the property list for HDF5 output
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);
    
    type_ = "Fields";
}


DiagnosticFields::~DiagnosticFields()
{
    H5Pclose( write_plist );
    H5Sclose( filespace );
    H5Sclose( memspace );
    
    delete timeSelection;
}


void DiagnosticFields::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if( fileId_>0 ) return;
    
    if ( newfile ) {
        // ----------------------------
        // Management of global IO file
        // ----------------------------
        MPI_Info info  = MPI_INFO_NULL;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
        
        // Fields.h5
        // ---------
        fileId_  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        
        // Create property list for collective dataset write: for Fields.h5
        vector<double> my_cell_length=params.cell_length;
        my_cell_length.resize(params.nDim_field);
        H5::attr(fileId_, "res_time"    , params.res_time);
        H5::attr(fileId_, "res_space"   , params.res_space);
        H5::attr(fileId_, "cell_length" , my_cell_length);
        H5::attr(fileId_, "sim_length"  , params.sim_length);
        
        H5Pclose(plist_id);
        
    }
    else {
        // Open the existing file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose(pid);
    }
    
}

void DiagnosticFields::closeFile()
{
    if( tmp_dset_id>0 ) H5Dclose( tmp_dset_id );
    tmp_dset_id=0;
    if( fileId_>0 ) H5Fclose( fileId_ );
    fileId_ = 0;
}



void DiagnosticFields::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    // create the file
    openFile( params, smpi, true );
    closeFile();
}

bool DiagnosticFields::prepare( int timestep )
{
    
    // Leave if the timestep is not the good one
    if (timestep - timeSelection->previousTime(timestep) >= time_average) return false;
    
    return true;
}


void DiagnosticFields::run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep )
{
    
    refHindex = (unsigned int)(vecPatches.refHindex_);
    
    setFileSplitting( smpi, vecPatches );
    
    // If time-averaging, increment the average
    if( time_average>1 )
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            vecPatches(ipatch)->EMfields->incrementAvgFields(timestep);
    
    // If is not writing timestep, leave
    if (timestep - timeSelection->previousTime(timestep) != time_average-1) return;
    
    // Create group for this timestep
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << timestep;
    
    htri_t status = H5Lexists(fileId_, name_t.str().c_str(), H5P_DEFAULT);
    // Do not output diag if this timestep has already been written
    if( status > 0 ) return;
    // Warning if file unreachable
    if( status < 0 ) {
        WARNING("Fields diagnostics could not write");
        return;
    }
    
    timestep_group_id = H5Gcreate(fileId_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // For each field, combine all patches and write out
    for( unsigned int ifield=0; ifield < fields_indexes.size(); ifield++ ) {
        
        // Copy the patch field to the buffer
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            getField( vecPatches(ipatch), fields_indexes[ifield] );
        
        // Create or open field dataset in HDF5
        hid_t dset_id;
        htri_t status = H5Lexists( timestep_group_id, fields_names[ifield].c_str(), H5P_DEFAULT );
        if (!status) {
            hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
            dset_id  = H5Dcreate( timestep_group_id, fields_names[ifield].c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
            H5Pclose(plist_id);
        } else {
            dset_id = H5Dopen( timestep_group_id, fields_names[ifield].c_str(), H5P_DEFAULT);                
        }
        
        // Write
        writeField(dset_id, timestep);
        
        // Close dataset
        H5Dclose( dset_id );
    }
    
    H5Gclose(timestep_group_id);
    
    // Final loop on patches to zero RhoJs
    if (fields_indexes.size()>0)
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            vecPatches(ipatch)->EMfields->restartRhoJs();
    
}

 