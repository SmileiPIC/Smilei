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

DiagnosticProbe::DiagnosticProbe(Params& params, SmileiMPI *smpi) {
    bool ok;
    
    // loop all "diagnostic probe" groups in the input file
    unsigned  numProbes=PyTools::nComponents("DiagProbe");
    
    if (numProbes>0) {
        // Create the HDF5 file that will contain all the probes
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        hid_t fileId = H5Fcreate("Probes.h5", H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Write the version of the code as an attribute
        H5::attr(fileId, "Version", string(__VERSION));
        H5::attr(fileId, "CommitDate", string(__COMMITDATE));
        
        dt = params.timestep;
        every         .clear();
        tmin          .clear();
        tmax          .clear();
        probeParticles.clear();
        nPart_total   .clear();
        probesArray   .clear();
        probesStart   .clear();
        fieldname     .clear();
        fieldlocation .clear();
        nFields       .clear();
        
        for (unsigned int n_probe = 0; n_probe < numProbes; n_probe++) {
        
            // Extract "every" (number of timesteps between each output)
            unsigned int my_every=0;
            ok=PyTools::extract("every",my_every,"DiagProbe",n_probe);
            if (!ok) my_every=params.global_every;
            every.push_back(my_every);
            
            // Extract "time_range" (tmin and tmax of the outputs)
            vector<double> time_range(2,0.);
            double my_tmin,my_tmax;
            ok=PyTools::extract("time_range",time_range,"DiagProbe",n_probe);
            if (!ok) {
                my_tmin = 0.;
                my_tmax = params.sim_time;
            } else {
                my_tmin = time_range[0];
                my_tmax = time_range[1];
            }
            tmin.push_back(my_tmin);
            tmax.push_back(my_tmax);
            
            // Extract "number" (number of points you have in each dimension of the probe,
            // which must be smaller than the code dimensions)
            vector<unsigned int> vecNumber;
            PyTools::extract("number",vecNumber,"DiagProbe",n_probe);
            
            // Dimension of the probe grid
            unsigned int dimProbe=vecNumber.size();
            if (dimProbe > params.nDim_particle) {
                ERROR("Probe #"<<n_probe<<": probe dimension is greater than simulation dimension")
            }
            
            // If there is no "number" argument provided, then it corresponds to
            // a zero-dimensional probe (one point). In this case, we say the probe
            // has actually one dimension with only one point.
            unsigned int dim=vecNumber.size();
            if (vecNumber.size() == 0) {
                vecNumber.resize(1);
                vecNumber[0]=1;
            }
            
            // Dimension of the simulation
            unsigned int ndim=params.nDim_particle;
            
            // Extract "pos", "pos_first", "pos_second" and "pos_third"
            // (positions of the vertices of the grid)
            vector< vector<double> > allPos;
            vector<double> pos;
            
            if (PyTools::extract("pos",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            if (PyTools::extract("pos_first",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos_first size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            if (PyTools::extract("pos_second",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos_second size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            if (PyTools::extract("pos_third",pos,"DiagProbe",n_probe)) {
                if (pos.size()!=ndim) {
                    ERROR("Probe #"<<n_probe<<": pos_third size(" << pos.size() << ") != ndim(" << ndim<< ")");
                }
                allPos.push_back(pos);
            }
            
            // Extract the list of requested fields
            vector<string> fs;
            if(!PyTools::extract("fields",fs,"DiagProbe",n_probe)) {
                fs.resize(10);
                fs[0]="Ex"; fs[1]="Ey"; fs[2]="Ez";
                fs[3]="Bx"; fs[4]="By"; fs[5]="Bz";
                fs[6]="Jx"; fs[7]="Jy"; fs[8]="Jz"; fs[9]="Rho";
            }
            vector<unsigned int> locations;
            locations.resize(10);
            for( unsigned int i=0; i<10; i++) locations[i] = fs.size();
            for( unsigned int i=0; i<fs.size(); i++) {
                for( unsigned int j=0; j<i; j++) {
                    if( fs[i]==fs[j] ) {
                        ERROR("Probe #"<<n_probe<<": field "<<fs[i]<<" appears twice");
                    }
                }
                if     ( fs[i]=="Ex" ) locations[0] = i;
                else if( fs[i]=="Ey" ) locations[1] = i;
                else if( fs[i]=="Ez" ) locations[2] = i;
                else if( fs[i]=="Bx" ) locations[3] = i;
                else if( fs[i]=="By" ) locations[4] = i;
                else if( fs[i]=="Bz" ) locations[5] = i;
                else if( fs[i]=="Jx" ) locations[6] = i;
                else if( fs[i]=="Jy" ) locations[7] = i;
                else if( fs[i]=="Jz" ) locations[8] = i;
                else if( fs[i]=="Rho") locations[9] = i;
                else {
                    ERROR("Probe #"<<n_probe<<": unknown field "<<fs[i]);
                }
            }
            fieldlocation.push_back(locations);
            fieldname.push_back(fs);
            nFields.push_back(fs.size());
            
            // Calculate the total number of points in the grid
            // Each point is actually a "fake" macro-particle
            unsigned int my_nPart=1;
            for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                my_nPart *= vecNumber[iDimProbe];
            }
            nPart_total.push_back(my_nPart);
            
            
            // Initialize the list of "fake" particles just as actual macro-particles
            Particles my_parts;
            my_parts.initialize(my_nPart, params.nDim_particle);
            
            // For each grid point, calculate its position and assign that position to the particle
            // The particle position is a linear combination of the `pos` with `pos_first` or `pos_second`, etc.
            double partPos, dx;
            vector<unsigned int> ipartND (dimProbe);
            for(unsigned int ipart=0; ipart<my_nPart; ++ipart) { // for each particle
                // first, convert the index `ipart` into N-D indexes
                unsigned int i = ipart;
                for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
                    ipartND[iDimProbe] = i%vecNumber[iDimProbe];
                    i = i/vecNumber[iDimProbe]; // integer division
                }
                // Now assign the position of the particle
                for(unsigned int iDim=0; iDim!=ndim; ++iDim) { // for each dimension of the simulation
                    partPos = allPos[0][iDim]; // position of `pos`
                    for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) { // for each of `pos`, `pos_first`, etc.
                        dx = (allPos[iDimProbe+1][iDim]-allPos[0][iDim])/(vecNumber[iDimProbe]-1); // distance between 2 gridpoints
                        partPos += ipartND[iDimProbe] * dx;
                    }
                    my_parts.position(iDim,ipart) = partPos;
                }
            }
            
            
            // Remove particles out of the domain
            for ( int ipb=my_nPart-1 ; ipb>=0 ; ipb--) {
                if (!my_parts.is_part_in_domain(ipb, smpi))
                    my_parts.erase_particle(ipb);
            }
            probeParticles.push_back(my_parts);
            
            unsigned int nPart_local = my_parts.size(); // number of fake particles for this proc
            
            // Make the array that will contain the data
            // probesArray : 10 x nPart_tot
            vector<unsigned int> probesArraySize(2);
            probesArraySize[1] = nPart_local; // number of particles
            probesArraySize[0] = nFields[n_probe] + 1; // number of fields (Ex, Ey, etc) +1 for garbage
            Field2D *myfield = new Field2D(probesArraySize);
            probesArray.push_back(myfield);
            
            // Exchange data between MPI cpus so that they can figure out which part
            // of the grid they have to manage
            MPI_Status status;
            // Receive the location where to start from the previous node
            int my_Start = 0;
            if (smpi->getRank()>0) MPI_Recv( &(my_Start), 1, MPI_INTEGER, smpi->getRank()-1, 0, MPI_COMM_WORLD, &status );
            // Send the location where to end to the next node
            int probeEnd = my_Start+nPart_local;
            if (smpi->getRank()!=smpi->getSize()-1) MPI_Send( &probeEnd, 1, MPI_INTEGER, smpi->getRank()+1, 0, MPI_COMM_WORLD );
            
            // Create group for the current probe
            ostringstream prob_name("");
            prob_name << "p" << setfill('0') << setw(4) << n_probe;
            hid_t gid = H5::group(fileId, prob_name.str());
            
            // Create an array to hold the positions of local probe particles
            Field2D fieldPosProbe;
            fieldPosProbe.allocateDims(ndim,nPart_local);
            
            for (unsigned int ipb=0 ; ipb<nPart_local ; ipb++)
                for (unsigned int idim=0 ; idim<ndim  ; idim++)
                    fieldPosProbe(idim,ipb) = my_parts.position(idim,ipb);
            
            // Add array "positions" into the current HDF5 group
            H5::matrix_MPI(gid, "positions", fieldPosProbe.data_2D[0][0], my_nPart, ndim, my_Start, nPart_local);
            
            probesStart.push_back(my_Start);
            
            // Add arrays "p0", "p1", ... to the current group
            ostringstream pk;
            for (unsigned int iDimProbe=0; iDimProbe<=dimProbe; iDimProbe++) {
                pk.str("");
                pk << "p" << iDimProbe;
                H5::vect(gid, pk.str(), allPos[iDimProbe]);
            }
            
            // Add array "number" to the current group
            H5::vect(gid, "number", vecNumber);
            
            // Add attribute every to the current group
            H5::attr(gid, "every", my_every);
            // Add attribute "dimension" to the current group
            H5::attr(gid, "dimension", dim);
            
            // Add "fields" to the current group
            ostringstream fields("");
            fields << fs[0];
            for( unsigned int i=1; i<fs.size(); i++) fields << "," << fs[i];
            H5::attr(gid, "fields", fields.str());
            
            // Close current group
            H5Gclose(gid);
            
        }
    
        H5Fclose(fileId);
    }
    
}

DiagnosticProbe::~DiagnosticProbe()
{
    for (unsigned int np=0 ; np < probesArray.size() ; np++ )
    delete probesArray[np];
}


string DiagnosticProbe::probeName(int p) {
    ostringstream prob_name("");
    prob_name << "p" << setfill('0') << setw(4) << p;
    return prob_name.str();
}

void DiagnosticProbe::run(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp) {

    double time = (double)timestep * dt;
    
    // Loop probes
    
    hid_t fileId=0;
    for (unsigned int np=0; np<every.size(); np++) {
        // skip if current timestep is not requested
        if ( (every[np]  && timestep % every[np] == 0) &&
             (time <= tmax[np]) && (time >= tmin[np]) ) {
            
            if (fileId==0) {
                // open the HDF5 file that will contain all the probes
                hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
                H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
                fileId = H5Fopen( "Probes.h5", H5F_ACC_RDWR, pid);
                H5Pclose(pid);
            }

            // Loop probe ("fake") particles
            unsigned int nPart_local = probeParticles[np].size();
            for (unsigned int iprob=0; iprob<nPart_local; iprob++) {
                
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
    if (fileId) H5Fclose(fileId);
}
