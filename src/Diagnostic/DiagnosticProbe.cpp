#include "DiagnosticProbe.h"

#include <iomanip>
#include <string>
#include <iomanip>
#include <sstream>

#include "Params.h"
#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"
#include "SmileiMPI_Cart2D.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field.h"
#include "H5.h"

using namespace std;

DiagnosticProbe::DiagnosticProbe():
fileId(0)
{
}


DiagnosticProbe::~DiagnosticProbe()
{
    for ( int np=0 ; np < probesArray.size() ; np++ )
    delete probesArray[np];
}


void DiagnosticProbe::close() {
    if (fileId>0) {
        H5Fclose(fileId);
    }
}

string DiagnosticProbe::probeName(int p) {
    ostringstream prob_name("");
    prob_name << "p" << setfill('0') << setw(4) << p;
    return prob_name.str();
}

void DiagnosticProbe::run(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp) {

    double time = (double)timestep * dt;
    
    // Loop probes
    for (unsigned int np=0; np<every.size(); np++) {
        // skip if current timestep is not requested
        if ( (every[np]  && timestep % every[np] == 0) &&
             (time <= tmax[np]) && (time >= tmin[np]) ) {
            
            // Loop probe ("fake") particles
            unsigned int nPart_local = probeParticles[np].size();
            for (int iprob=0; iprob<nPart_local; iprob++) {
                
                // Interpolate fields at the location of the fake particles
                (*interp)(EMfields,probeParticles[np],iprob,&Eloc_fields,&Bloc_fields,&Jloc_fields,&Rloc_fields);
                
                //! here we fill the probe data!!!
                probesArray[np]->data_2D[fieldlocation[np][0]][iprob]=Eloc_fields.x;
                probesArray[np]->data_2D[fieldlocation[np][1]][iprob]=Eloc_fields.y;
                probesArray[np]->data_2D[fieldlocation[np][2]][iprob]=Eloc_fields.z;
                probesArray[np]->data_2D[fieldlocation[np][3]][iprob]=Bloc_fields.x;
                probesArray[np]->data_2D[fieldlocation[np][4]][iprob]=Bloc_fields.y;
                probesArray[np]->data_2D[fieldlocation[np][5]][iprob]=Bloc_fields.z;
                probesArray[np]->data_2D[fieldlocation[np][6]][iprob]=Jloc_fields.x;
                probesArray[np]->data_2D[fieldlocation[np][7]][iprob]=Jloc_fields.y;
                probesArray[np]->data_2D[fieldlocation[np][8]][iprob]=Jloc_fields.z;
                probesArray[np]->data_2D[fieldlocation[np][9]][iprob]=Rloc_fields;
                
            }
            
            // Make the name for the array
            ostringstream name_t;
            name_t.str("");
            name_t << "/" << probeName(np).c_str() << "/" << setfill('0') << setw(10) << timestep;
            
            // Open the existing HDF5 group for that probe
            hid_t did = H5Gopen2(fileId, probeName(np).c_str(), H5P_DEFAULT);
            // Write the positions array into the current HDF5 group
            H5::matrix_MPI(did, name_t.str(), probesArray[np]->data_2D[0][0], nPart_total[np], nFields[np], probesStart[np], nPart_local);
            // Close the group
            H5Gclose(did);
            
        }
    }
    if (fileId) H5Fflush(fileId, H5F_SCOPE_GLOBAL );
}
