#include "Diagnostic.h"

#include <string>
#include <iomanip>

#include <hdf5.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "DiagnosticPhaseMomMom.h"
#include "DiagnosticPhasePosLor.h"
#include "DiagnosticPhasePosMom.h"


using namespace std;

Diagnostic::Diagnostic(PicParams& params, InputData &ifile, SmileiMPI *smpi) :
dtimer(4)
{
    dtimer[0].init(smpi, "scalars");
    dtimer[1].init(smpi, "probes");
    dtimer[2].init(smpi, "phases");
    dtimer[3].init(smpi, "particles");
    
    // defining default values & reading diagnostic every-parameter
    // ------------------------------------------------------------
    print_every=params.n_time/10;
    ifile.extract("print_every", print_every);
    
    fieldDump_every=0;
    if (!ifile.extract("fieldDump_every", fieldDump_every)) fieldDump_every=params.global_every;
    
    if (!ifile.extract("fieldDump_every", fieldDump_every)) {
        DEBUG("activating all fields to dump");
    }    
    
    avgfieldDump_every=params.res_time*10;
    if (!ifile.extract("avgfieldDump_every", avgfieldDump_every)) avgfieldDump_every=params.global_every;
    
    //!\todo Define default behaviour : 0 or params.res_time
    //ntime_step_avg=params.res_time;
    ntime_step_avg=0;
    ifile.extract("ntime_step_avg", ntime_step_avg);
    
    particleDump_every=0;
    if (ifile.extract("particleDump_every", particleDump_every))
        WARNING("Option particleDump_every disabled");
    
    // scalars initialization
    initScalars(params,ifile, smpi);
    
    // probes initialization
    initProbes(params,ifile,smpi);
    
    // phasespaces initialization    
    initPhases(params,ifile,smpi);
    
    // particles initialization    
    initParticles(params,ifile);
    
}

void Diagnostic::closeAll (SmileiMPI* smpi) {
    
    scalars.closeFile(smpi);
    probes.close();
    phases.close();
    
    for (int i=0; i<vecDiagnosticParticles.size(); i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->close();
    
}

void Diagnostic::printTimers (SmileiMPI *smpi, double tottime) {
    
    double coverage(0.);
    if ( smpi->isMaster() ) {
        for (int i=0 ; i<dtimer.size() ; i++) {
            coverage += dtimer[i].getTime();
        }
    }
    MESSAGE(0, "\nTime in diagnostics : \t"<< tottime <<"\t(" << coverage/tottime*100. << "% coverage)" );    
    if ( smpi->isMaster() ) {
        for (int i=0 ; i<dtimer.size() ; i++) {
            dtimer[i].print(tottime) ;
        }
    }
}

double Diagnostic::getScalar(string name){
    return scalars.getScalar(name);
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator *interp, SmileiMPI *smpi) {
    dtimer[0].restart();
    scalars.run(timestep, EMfields, vecSpecies, smpi);
    dtimer[0].update();
    
    dtimer[1].restart();
    probes.run(timestep, EMfields, interp);
    dtimer[1].update();
    
    dtimer[2].restart();
    phases.run(timestep, vecSpecies);
    dtimer[2].update();
    
    // run all the particle diagnostics
    dtimer[3].restart();
    for (int i=0; i<vecDiagnosticParticles.size(); i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->run(timestep, vecSpecies, smpi);
    dtimer[3].update();

}

void Diagnostic::initScalars(PicParams& params, InputData &ifile, SmileiMPI *smpi) {
    //open file scalars.txt
    scalars.openFile(smpi);
    
    MESSAGE(1,"Initilizing scalars");
    scalars.every=0;
    bool ok=ifile.extract("every",scalars.every,"DiagScalar");
    if (!ok) scalars.every=params.global_every;
    
    vector<double> scalar_time_range(2,0.);
    
    ok=ifile.extract("time_range",scalar_time_range,"DiagScalar");        
    if (!ok) { 
        scalars.tmin = 0.;
        scalars.tmax = params.sim_time;
    }
    else {
        scalars.tmin = scalar_time_range[0]*params.conv_fac;
        scalars.tmax = scalar_time_range[1]*params.conv_fac;
    }
    
    scalars.precision=10;
    ifile.extract("precision",scalars.precision,"DiagScalar");
    ifile.extract("vars",scalar_vars,"DiagScalar");
    
    // copy from params remaining stuff
    scalars.res_time=params.res_time;
    scalars.dt=params.timestep;
    scalars.cell_volume=params.cell_volume;
}

void Diagnostic::initProbes(PicParams& params, InputData &ifile, SmileiMPI *smpi) {
    
    bool ok;
    
    // loop all "diagnostic probe" groups in the input file
    unsigned  numProbes=ifile.nComponents("DiagProbe");
    MESSAGE(1,"Initializing Probes (" << numProbes << ")");
    for (unsigned int n_probe = 0; n_probe < numProbes; n_probe++) {
        
        if (n_probe==0) {
            // Create the HDF5 file that will contain all the probes
            hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
            probes.fileId = H5Fcreate( "Probes.h5", H5F_ACC_TRUNC, H5P_DEFAULT, pid);
            H5Pclose(pid);
            
            // Write the version of the code as an attribute
            string ver(__VERSION);
            H5::attr(probes.fileId, "Version", ver);
            
            probes.dt = params.timestep;
            probes.every         .resize(0);
            probes.tmin          .resize(0);
            probes.tmax          .resize(0);
            probes.probeParticles.resize(0);
            probes.nPart_total   .resize(0);
            probes.probesArray   .resize(0);
            probes.probesStart   .resize(0);
        }
        
        
        // Extract "every" (number of timesteps between each output)
        unsigned int every=0;
        ok=ifile.extract("every",every,"DiagProbe",n_probe);        
        if (!ok) every=params.global_every;
        probes.every.push_back(every);
        
        // Extract "time_range" (tmin and tmax of the outputs)
        vector<double> time_range(2,0.);
        double tmin,tmax;
        ok=ifile.extract("time_range",time_range,"DiagProbe",n_probe);        
        if (!ok) {
            tmin = 0.;
            tmax = params.sim_time;
        } else {
            tmin = time_range[0]*params.conv_fac;
            tmax = time_range[1]*params.conv_fac;
        }
        probes.tmin.push_back(tmin);
        probes.tmax.push_back(tmax);
        
        // Extract "number" (number of points you have in each dimension of the probe,
        // which must be smaller than the code dimensions)
        vector<unsigned int> vecNumber; 
        ifile.extract("number",vecNumber,"DiagProbe",n_probe);
        
        // Dimension of the probe grid
        unsigned int dimProbe=vecNumber.size();
        if (dimProbe > params.nDim_particle) {
            ERROR("probe dimension is greater than simulation dimension")
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
        ifile.extract("pos",pos,"DiagProbe",n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= params.conv_fac;
        if (pos.size()>0) allPos.push_back(pos);
        
        ifile.extract("pos_first",pos,"DiagProbe",n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= params.conv_fac;
        if (pos.size()>0) allPos.push_back(pos);
        
        ifile.extract("pos_second",pos,"DiagProbe",n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= params.conv_fac;
        if (pos.size()>0) allPos.push_back(pos);
        
        ifile.extract("pos_third",pos,"DiagProbe",n_probe);
        for (unsigned int i=0; i<pos.size(); i++)
            pos[i] *= params.conv_fac;
        if (pos.size()>0) allPos.push_back(pos);
        
        
        // Calculate the total number of points in the grid
        // Each point is actually a "fake" macro-particle
        unsigned int nPart_total=1;
        for (unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++) {
            nPart_total *= vecNumber[iDimProbe];
        }
        probes.nPart_total.push_back(nPart_total);
        
        
        // Initialize the list of "fake" particles just as actual macro-particles
        Particles probeParticles;
        probeParticles.initialize(nPart_total, ndim);
        
        // For each grid point, calculate its position and assign that position to the particle
        // The particle position is a linear combination of the `pos` with `pos_first` or `pos_second`, etc.
        double partPos, dx;
        vector<unsigned int> ipartND (dimProbe);
        for(unsigned int ipart=0; ipart<nPart_total; ++ipart) { // for each particle
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
                probeParticles.position(iDim,ipart) = partPos;
            }
        }
        
        
        // Remove particles out of the domain
        for ( int ipb=nPart_total-1 ; ipb>=0 ; ipb--) {
            if (!probeParticles.is_part_in_domain(ipb, smpi))
                probeParticles.erase_particle(ipb);
        }
        probes.probeParticles.push_back(probeParticles);
        
        unsigned int nPart_local = probeParticles.size(); // number of fake particles for this proc
        
        // Make the array that will contain the data
        // probesArray : 10 x nPart_tot
        vector<unsigned int> probesArraySize(2);
        probesArraySize[0] = nPart_local; // number of particles
        probesArraySize[1] = probes.probeSize; // number of fields (Ex, Ey, etc)
        Field2D *myfield = new Field2D(probesArraySize);
        probes.probesArray.push_back(myfield);
        
        // Exchange data between MPI cpus so that they can figure out which part
        // of the grid they have to manage
        MPI_Status status;
        // Receive the location where to start from the previous node
        int probesStart = 0;
        if (smpi->getRank()>0) MPI_Recv( &(probesStart), 1, MPI_INTEGER, smpi->getRank()-1, 0, MPI_COMM_WORLD, &status );
        // Send the location where to end to the next node
        int probeEnd = probesStart+nPart_local;
        if (smpi->getRank()!=smpi->getSize()-1) MPI_Send( &probeEnd, 1, MPI_INTEGER, smpi->getRank()+1, 0, MPI_COMM_WORLD );
        
        // Create group for the current probe
        ostringstream prob_name("");
        prob_name << "p" << setfill('0') << setw(4) << n_probe;
        hid_t did = H5Gcreate(probes.fileId, prob_name.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        // Create an array to hold the positions of local probe particles
        double posArray [nPart_local][ndim];
        for (unsigned int ipb=0 ; ipb<nPart_local ; ipb++)
            for (unsigned int idim=0 ; idim<ndim  ; idim++)
                posArray[ipb][idim] = probeParticles.position(idim,ipb);
        
        // Add array "positions" into the current HDF5 group
        H5::matrix_MPI(did, "positions", posArray[0][0], nPart_total, ndim, probesStart, nPart_local);
        
        probes.probesStart.push_back(probesStart);
        
        // Add arrays "p0", "p1", ... to the current group
        ostringstream pk;
        for (unsigned int iDimProbe=0; iDimProbe<=dimProbe; iDimProbe++) {
            pk.str("");
            pk << "p" << iDimProbe;
            H5::vector(did, pk.str(), allPos[iDimProbe][0], ndim);
        }
        
        // Add array "number" to the current group
        H5::vector(did, "number", vecNumber[0], dimProbe);
        
        // Add attribute every to the current group
        H5::attr(did, "every", every);
        // Add attribute "dimension" to the current group
        H5::attr(did, "dimension", dim);
        
        // Close current group
        H5Gclose(did);
        
    }
}

void Diagnostic::initPhases(PicParams& params, InputData &ifile, SmileiMPI *smpi) {
    
    //! create the particle structure
    phases.ndim=params.nDim_particle;    
    phases.my_part.pos.resize(params.nDim_particle);
    phases.my_part.mom.resize(3);
    
    bool ok;
    
    unsigned int numPhases=ifile.nComponents("DiagPhase");
    MESSAGE(1,"Initializing PhaseSpaces (" << numPhases << ")");
    for (unsigned int n_phase = 0; n_phase < numPhases; n_phase++) {
        
        phaseStructure my_phase;
        
        my_phase.every=0;
        ok=ifile.extract("every",my_phase.every,"DiagPhase",n_phase);
        if (!ok) {
            //            if (n_probephase>0) {
            //                my_phase.every=phases.vecDiagPhase.end()->every;
            //            } else {
            my_phase.every=params.global_every;
            //            }
        }
        
        vector<string> kind;
        ifile.extract("kind",kind,"DiagPhase",n_phase);        
        for (vector<string>::iterator it=kind.begin(); it!=kind.end();it++) {
            if (std::find(kind.begin(), it, *it) == it) {
                my_phase.kind.push_back(*it); 
            } else {
                WARNING("removed duplicate " << *it << " in \"DiagPhase\" " << n_phase);
            }
        }
        
        vector<double> time_range(2,0.);
        ok=ifile.extract("time_range",time_range,"DiagPhase",n_phase);        
        
        if (!ok) { 
            my_phase.tmin = 0.;
            my_phase.tmax = params.sim_time;
        }
        else {
            my_phase.tmin = time_range[0]*params.conv_fac;
            my_phase.tmax = time_range[1]*params.conv_fac;
        }
        
        
        ifile.extract("species",my_phase.species,"DiagPhase",n_phase);
        
        my_phase.deflate=0;
        ifile.extract("deflate",my_phase.deflate,"DiagPhase",n_phase);
        
        if (my_phase.species.size()==0) {
            WARNING("adding all species to the \"DiagPhase\" " << n_phase);
            for (unsigned int i=0;i<params.species_param.size(); i++) {
                my_phase.species.push_back(params.species_param[i].species_type);
            }
        }
        
        ifile.extract("pos_min",my_phase.pos_min,"DiagPhase",n_phase);
        ifile.extract("pos_max",my_phase.pos_max,"DiagPhase",n_phase);
        ifile.extract("pos_num",my_phase.pos_num,"DiagPhase",n_phase);
        for (unsigned int i=0; i<my_phase.pos_min.size(); i++) {
            my_phase.pos_min[i] *= params.conv_fac;
            my_phase.pos_max[i] *= params.conv_fac;
            if (my_phase.pos_min[i]==my_phase.pos_max[i]) {
                my_phase.pos_min[i] = 0.0;
                my_phase.pos_max[i] = params.sim_length[i];
            }
        }
        
        
        ifile.extract("mom_min",my_phase.mom_min,"DiagPhase",n_phase);
        ifile.extract("mom_max",my_phase.mom_max,"DiagPhase",n_phase);
        ifile.extract("mom_num",my_phase.mom_num,"DiagPhase",n_phase);
        
        ifile.extract("lor_min",my_phase.lor_min,"DiagPhase",n_phase);
        ifile.extract("lor_max",my_phase.lor_max,"DiagPhase",n_phase);
        ifile.extract("lor_num",my_phase.lor_num,"DiagPhase",n_phase);
        
        
        hid_t gidParent=0;
        if (n_phase == 0 && smpi->isMaster()) {
            ostringstream file_name("");
            file_name<<"PhaseSpace.h5";
            phases.fileId = H5Fcreate( file_name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            string ver(__VERSION);
            
            // write version
            hid_t aid3  = H5Screate(H5S_SCALAR);
            hid_t atype = H5Tcopy(H5T_C_S1);
            H5Tset_size(atype, ver.size());
            H5Tset_strpad(atype,H5T_STR_NULLTERM);
            hid_t attr3 = H5Acreate(phases.fileId, "Version", atype, aid3, H5P_DEFAULT, H5P_DEFAULT);
            
            H5Awrite(attr3, atype, ver.c_str());
            
            H5Aclose(attr3);
            H5Sclose(aid3);
            H5Tclose(atype);            
            
            ostringstream groupName("");
            groupName << "ps" << setw(4) << setfill('0') << n_phase;
            gidParent = H5Gcreate(phases.fileId, groupName.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
            
            hid_t sid = H5Screate(H5S_SCALAR);	
            hid_t aid = H5Acreate(gidParent, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(aid, H5T_NATIVE_UINT, &my_phase.every);
            H5Sclose(sid);
            H5Aclose(aid);
        }
        
        
        for (unsigned int ii=0 ; ii < my_phase.kind.size(); ii++) {
            DiagnosticPhase *diagPhase=NULL;
            
            // create DiagnosticPhase
            if (params.geometry == "1d3v") {
                if (my_phase.kind[ii] == "xpx") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,0);
                } else if (my_phase.kind[ii] == "xpy") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,1);
                } else if (my_phase.kind[ii] == "xpz") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,2);
                } else if (my_phase.kind[ii] == "xlor") {
                    diagPhase =  new DiagnosticPhasePosLor(my_phase,0);
                } else if (my_phase.kind[ii] == "pxpy") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,1);
                } else if (my_phase.kind[ii] == "pxpz") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,2);
                } else if (my_phase.kind[ii] == "pypz") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,1,2);
                } else {
                    ERROR("kind " << my_phase.kind[ii] << " not implemented for geometry " << params.geometry);
                }
            } else if (params.geometry == "2d3v") {
                if (my_phase.kind[ii] == "xpx") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,0);
                } else if (my_phase.kind[ii] == "xpy") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,1);
                } else if (my_phase.kind[ii] == "xpz") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,2);
                } else if (my_phase.kind[ii] == "ypx") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,0);
                } else if (my_phase.kind[ii] == "ypy") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,1);
                } else if (my_phase.kind[ii] == "ypz") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,2);
                } else if (my_phase.kind[ii] == "pxpy") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,1);
                } else if (my_phase.kind[ii] == "pxpz") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,2);
                } else if (my_phase.kind[ii] == "pypz") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,1,2);                    
                } else if (my_phase.kind[ii] == "xlor") {
                    diagPhase =  new DiagnosticPhasePosLor(my_phase,0);
                } else if (my_phase.kind[ii] == "ylor") {
                    diagPhase =  new DiagnosticPhasePosLor(my_phase,1);
                } else {
                    ERROR("kind " << my_phase.kind[ii] << " not implemented for geometry " << params.geometry);
                }
            } else if (params.geometry == "3d3v") {
                if (my_phase.kind[ii] == "xpx") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,0);
                } else if (my_phase.kind[ii] == "xpy") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,1);
                } else if (my_phase.kind[ii] == "xpz") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,0,2);
                } else if (my_phase.kind[ii] == "ypx") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,0);
                } else if (my_phase.kind[ii] == "ypy") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,1);
                } else if (my_phase.kind[ii] == "ypz") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,1,2);
                } else if (my_phase.kind[ii] == "zpx") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,2,0);
                } else if (my_phase.kind[ii] == "zpy") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,2,1);
                } else if (my_phase.kind[ii] == "zpz") {
                    diagPhase =  new DiagnosticPhasePosMom(my_phase,2,2);
                } else if (my_phase.kind[ii] == "pxpy") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,1);
                } else if (my_phase.kind[ii] == "pxpz") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,0,2);
                } else if (my_phase.kind[ii] == "pypz") {
                    diagPhase =  new DiagnosticPhaseMomMom(my_phase,1,2);                    
                } else if (my_phase.kind[ii] == "xlor") {
                    diagPhase =  new DiagnosticPhasePosLor(my_phase,0);
                } else if (my_phase.kind[ii] == "ylor") {
                    diagPhase =  new DiagnosticPhasePosLor(my_phase,1);
                } else if (my_phase.kind[ii] == "zlor") {
                    diagPhase =  new DiagnosticPhasePosLor(my_phase,2);
                } else {
                    ERROR("kind " << my_phase.kind[ii] << " not implemented for geometry " << params.geometry);
                }                
            } else {
                ERROR("DiagnosticPhase not implemented for geometry " << params.geometry);
            }
            if (diagPhase) {
                if (smpi->isMaster()) {
                    //! create a group for each species of this diag and keep track of its ID.
                    
                    hsize_t dims[3] = {0,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                    hsize_t max_dims[3] = {H5S_UNLIMITED,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                    hsize_t chunk_dims[3] = {1,diagPhase->my_data.dims()[0],diagPhase->my_data.dims()[1]};
                    
                    hid_t sid = H5Screate_simple (3, dims, max_dims);	
                    hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
                    H5Pset_layout(pid, H5D_CHUNKED);
                    H5Pset_chunk(pid, 3, chunk_dims);
                    
                    H5Pset_deflate (pid, std::min((unsigned int)9,my_phase.deflate));
                    
                    diagPhase->dataId = H5Dcreate (gidParent, my_phase.kind[ii].c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid,H5P_DEFAULT);
                    H5Pclose (pid);	
                    H5Sclose (sid);
                    
                    // write attribute of species present in the phaseSpace
                    string namediag;
                    for (unsigned int k=0; k<my_phase.species.size(); k++) {
                        namediag+=my_phase.species[k]+" ";
                    }
                    namediag=namediag.substr(0, namediag.size()-1);
                    sid = H5Screate(H5S_SCALAR);
                    hid_t tid = H5Tcopy(H5T_C_S1);
                    H5Tset_size(tid, namediag.size());
                    H5Tset_strpad(tid,H5T_STR_NULLTERM);
                    hid_t aid = H5Acreate(gidParent, "species", tid, sid, H5P_DEFAULT, H5P_DEFAULT);
                    H5Awrite(aid, tid, namediag.c_str());
                    
                    
                    // write attribute extent of the phaseSpace
                    hsize_t dimsPos[2] = {2,2};
                    sid = H5Screate_simple(2, dimsPos, NULL);
                    aid = H5Acreate (gidParent, "extents", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
                    double tmp[4] = {diagPhase->firstmin, diagPhase->firstmax, diagPhase->secondmin, diagPhase->secondmax};
                    H5Awrite(aid, H5T_NATIVE_DOUBLE, tmp);
                    H5Aclose(aid);
                    H5Sclose(sid);
                    
                }
                phases.vecDiagPhase.push_back(diagPhase);	
            }
            
        } 
        
        if (smpi->isMaster() ) {
            H5Gclose(gidParent);
        }
    }
}


void Diagnostic::initParticles(PicParams& params, InputData &ifile) {
    unsigned int every, time_average;
    string output;
    vector<string> species;
    vector<unsigned int> species_numbers;
    DiagnosticParticlesAxis  *tmpAxis;
    vector<DiagnosticParticlesAxis*> tmpAxes;
    DiagnosticParticles * tmpDiagParticles;
    vector<PyObject*> allAxes;
    
    bool ok;
    
    unsigned int numDiagParticles=ifile.nComponents("DiagParticles");
    MESSAGE(1,"Initializing Particles (" << numDiagParticles << ")");
    for (unsigned int n_diag_particles = 0; n_diag_particles < numDiagParticles; n_diag_particles++) {
        
        // get parameter "output" that determines the quantity to sum in the output array
        output = "";
        ok = ifile.extract("output",output,"DiagParticles",n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `output` required");
        
        // get parameter "every" which is the period (in timesteps) for getting the outputs
        every = 0;
        ok = ifile.extract("every",every,"DiagParticles",n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `every` required");
        
        // get parameter "time_average" that determines the number of timestep to average the outputs
        time_average = 1;
        ifile.extract("time_average",time_average,"DiagParticles",n_diag_particles);
        if (time_average > every)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": `time_average` cannot be larger than `every`");
        if (time_average < 1) time_average=1;
        
        // get parameter "species" that determines the species to use (can be a list of species)
        species.resize(0);
        ok = ifile.extract("species",species,"DiagParticles",n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `species` required");
        // verify that the species exist, remove duplicates and sort by number
        species_numbers = params.FindSpecies(species);
        
        
        // get parameter "axes" that adds axes to the diagnostic
        // Each axis should contain several items:
        //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
        allAxes=ifile.extract_pyVec("axes","DiagParticles",n_diag_particles);
        
        if (allAxes.size() == 0)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": axes must contain something");
        
        tmpAxes.resize(0);
        for (unsigned int iaxis=0; iaxis<allAxes.size(); iaxis++ ) {
            tmpAxis = new DiagnosticParticlesAxis();
            PyObject *oneAxis=allAxes[iaxis];
            if (PyTuple_Check(oneAxis) || PyList_Check(oneAxis)) {
                PyObject* seq = PySequence_Fast(oneAxis, "expected a sequence");
                unsigned int lenAxisArgs=PySequence_Size(seq);
                if (lenAxisArgs<4)
                    ERROR("Diagnotic Particles #" << n_diag_particles << ": axis #" << iaxis << " contain at least 4 arguments");
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 0),tmpAxis->type)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": First item must be a string (axis type)");
                } else {
                    if (   (tmpAxis->type == "z" && params.nDim_particle <3)
                        || (tmpAxis->type == "y" && params.nDim_particle <2) )
                        ERROR("Diagnotic Particles #" << n_diag_particles << ": axis " << tmpAxis->type << " cannot exist in " << params.nDim_particle << "D");
                }
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 1),tmpAxis->min)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": Second item must be a double (axis min)");
                }
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 2),tmpAxis->max)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": Third item must be a double (axis max)");
                }
                
                // If the axis is spatial, then we need to apply the conv_fac
                if (tmpAxis->type=="x" || tmpAxis->type=="y" || tmpAxis->type=="z") {
                    tmpAxis->min *= params.conv_fac;
                    tmpAxis->max *= params.conv_fac;
                }
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 3),tmpAxis->nbins)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": Fourth item must be an int (number of bins)");
                }
                
                // 5 - Check for  other keywords such as "logscale" and "edge_inclusive"
                tmpAxis->logscale = false;
                tmpAxis->edge_inclusive = false;
                for(unsigned int i=4; i<lenAxisArgs; i++) {
                    string my_str("");
                    PyTools::convert(PySequence_Fast_GET_ITEM(seq, i),my_str);
                    if(my_str=="logscale" ||  my_str=="log_scale" || my_str=="log")
                        tmpAxis->logscale = true;
                    else if(my_str=="edges" ||  my_str=="edge" ||  my_str=="edge_inclusive" ||  my_str=="edges_inclusive")
                        tmpAxis->edge_inclusive = true;
                    else
                        ERROR("Diagnotic Particles #" << n_diag_particles << ": keyword `" << my_str << "` not understood");
                }
                
                tmpAxes.push_back(tmpAxis);
                
                Py_DECREF(seq);
            }
            
            
            HEREIAM(iaxis);
        }
        
        
        
        // create new diagnostic object
        tmpDiagParticles = new DiagnosticParticles(n_diag_particles, output, every, time_average, species_numbers, tmpAxes);
        // add this object to the list
        vecDiagnosticParticles.push_back(tmpDiagParticles);
    }
}


