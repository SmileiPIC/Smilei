#include "DiagnosticProbe.h"

#include <iomanip>
#include <string>
#include <iomanip>
#include <sstream>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"
#include "SmileiMPI_Cart2D.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field.h"
#include "DiagParams.h"

using namespace std;

DiagnosticProbe::DiagnosticProbe(PicParams &params, DiagParams &diagParams, SmileiMPI* smpi):
cpuRank((int)smpi->getRank()),
probeSize(10), 
fileId(0) {
    
    // diagParams.probeStruc is an array of several probeStructure objects.
    // They were collected in DiagParams.cpp by reading the input file.
    unsigned int nprobes = diagParams.probeStruc.size();
    if ( nprobes>0) {
        
        // Create the HDF5 file that will contain all the probes
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId = H5Fcreate( "Probes.h5", H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Write the version of the code as an attribute
        string ver(__VERSION);
        H5::attr(fileId, "Version", ver);
        
        // Initialize arrays with 1 element for each probe
        dt = params.timestep;
        every         .resize(nprobes);
        tmin          .resize(nprobes);
        tmax          .resize(nprobes);
        probeParticles.resize(nprobes);
        nPart_total   .resize(nprobes);
        probesArray   .resize(nprobes);
        probesStart   .resize(nprobes);
        
        // Loop each probe and collect the data from the probeStructure
        for (unsigned int np=0; np<diagParams.probeStruc.size(); np++) {
            every[np]=diagParams.probeStruc[np].every;
            tmin[np] =diagParams.probeStruc[np].tmin;
            tmax[np] =diagParams.probeStruc[np].tmax;
            
            // dimension of the probe grid
            unsigned int dimProbe=diagParams.probeStruc[np].dim;
            // dimension of the simulation
            unsigned int ndim=params.nDim_particle;
            // "shape" of the probe grid
            vector<unsigned int> vecNumber=diagParams.probeStruc[np].number;
            // list of the positions of the points "pos", "pos_first", etc.
            vector< vector<double> > * p = &(diagParams.probeStruc[np].pos);
            
            // Calculate the total number of points in the grid
            // Each point is actually a "fake" macro-particle
            nPart_total[np]=1;
            for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                nPart_total[np] *= vecNumber[iDimProbe];
            }
            
            // Initialize the list of "fake" particles the same way are actual macro-particles
            probeParticles[np].initialize(nPart_total[np], ndim);
            
            // For each grid point, calculate its position and assign that position to the particle
            // The particle position is a linear combination of the `pos` with `pos_first` or `pos_second`, etc.
            double partPos, dx;
            vector<unsigned int> ipartND (dimProbe);
            for(unsigned int ipart=0; ipart<nPart_total[np]; ++ipart) { // for each particle
                // first, convert the index `ipart` into N-D indexes
                unsigned int i = ipart;
                for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                    ipartND[iDimProbe] = i%vecNumber[iDimProbe];
                    i = i/vecNumber[iDimProbe]; // integer division
                }
                // Now assign the position of the particle
                for(unsigned int iDim=0; iDim!=ndim; ++iDim) { // for each dimension of the simulation
                    partPos = (*p)[0][iDim]; // position of `pos`
                    for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) { // for each of `pos`, `pos_first`, etc.
                        dx = ((*p)[iDimProbe+1][iDim]-(*p)[0][iDim])/(vecNumber[iDimProbe]-1); // distance between 2 gridpoints
                        partPos += ipartND[iDimProbe] * dx;
                    }
                    probeParticles[np].position(iDim,ipart) = partPos;
                }
            }
            
            // Remove particles out of the domain
            for ( int ipb=nPart_total[np]-1 ; ipb>=0 ; ipb--) {
                if (!probeParticles[np].is_part_in_domain(ipb, smpi))
                    probeParticles[np].erase_particle(ipb);
            }
            unsigned int nPart_local = probeParticles[np].size(); // number of fake particles for this proc
            
            // Make the array that will contain the data
            // probesArray : 10 x nPart_tot
            vector<unsigned int> probesArraySize(2);
            probesArraySize[0] = nPart_local; // number of particles
            probesArraySize[1] = probeSize; // number of fields (Ex, Ey, etc)
            probesArray[np] = new Field2D(probesArraySize);
            
            // Exchange data between MPI cpus so that they can figure out which part
            // of the grid they have to manage
            MPI_Status status;
            // Receive the location where to start from the previous node
            probesStart[np] = 0;
            if (cpuRank>0) MPI_Recv( &(probesStart[np]), 1, MPI_INTEGER, cpuRank-1, 0, MPI_COMM_WORLD, &status );
            // Send the location where to end to the next node
            int probeEnd = probesStart[np]+nPart_local;
            if (cpuRank!=smpi->getSize()-1) MPI_Send( &probeEnd, 1, MPI_INTEGER, cpuRank+1, 0, MPI_COMM_WORLD );
            
            // Create group for the current probe
            hid_t did = H5Gcreate(fileId, probeName(np).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
            // Create an array to hold the positions of local probe particles
            double posArray [nPart_local][ndim];
            for (int ipb=0 ; ipb<nPart_local ; ipb++)
                for (int idim=0 ; idim<ndim  ; idim++)
                    posArray[ipb][idim] = probeParticles[np].position(idim,ipb);
            
            // Add array "positions" into the current HDF5 group
            H5::matrix_MPI(did, "positions", posArray[0][0], nPart_total[np], ndim, probesStart[np], nPart_local);
            
            // Add arrays "p0", "p1", ... to the current group
            ostringstream pk;
            for (unsigned int iDimProbe=0; iDimProbe<=dimProbe; iDimProbe++) {
                pk.str("");
                pk << "p" << iDimProbe;
                H5::vector(did, pk.str(), diagParams.probeStruc[np].pos[iDimProbe][0], ndim);
            }
            
            // Add array "number" to the current group
            H5::vector(did, "number", diagParams.probeStruc[np].number[0], dimProbe);
            
            // Add attribute every to the current group
            H5::attr(did, "every", every[np]);
            // Add attribute "dimension" to the current group
            H5::attr(did, "dimension", diagParams.probeStruc[np].dim);
            
            // Close current group
            H5Gclose(did);
        }
    
    }
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
                probesArray[np]->data_2D[iprob][0]=Eloc_fields.x;
                probesArray[np]->data_2D[iprob][1]=Eloc_fields.y;
                probesArray[np]->data_2D[iprob][2]=Eloc_fields.z;
                probesArray[np]->data_2D[iprob][3]=Bloc_fields.x;
                probesArray[np]->data_2D[iprob][4]=Bloc_fields.y;
                probesArray[np]->data_2D[iprob][5]=Bloc_fields.z;
                probesArray[np]->data_2D[iprob][6]=Jloc_fields.x;
                probesArray[np]->data_2D[iprob][7]=Jloc_fields.y;
                probesArray[np]->data_2D[iprob][8]=Jloc_fields.z;
                probesArray[np]->data_2D[iprob][9]=Rloc_fields;
                
            }
            
            // Make the name for the array
            ostringstream name_t;
            name_t.str("");
            name_t << "/" << probeName(np).c_str() << "/" << setfill('0') << setw(10) << timestep;
            
            // Open the existing HDF5 group for that probe
            hid_t did = H5Gopen2(fileId, probeName(np).c_str(), H5P_DEFAULT);
            // Write the positions array into the current HDF5 group
            H5::matrix_MPI(did, name_t.str(), probesArray[np]->data_2D[0][0], nPart_total[np], probeSize, probesStart[np], nPart_local);
            // Close the group
            H5Gclose(did);
            
        }
    }
    if (fileId) H5Fflush(fileId, H5F_SCOPE_GLOBAL );
}
