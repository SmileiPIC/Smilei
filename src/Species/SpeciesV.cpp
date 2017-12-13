#include "SpeciesV.h"

#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>

#include <omp.h>

// IDRIS
#include <cstring>
// IDRIS
#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "PartBoundCond.h"
#include "PartWall.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "Profile.h"

#include "Projector.h"
#include "ProjectorFactory.h"

#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
SpeciesV::SpeciesV(Params& params, Patch* patch) :
    Species(params, patch)
{
    initCluster( params );

}//END SpeciesV creator

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
SpeciesV::~SpeciesV()
{
}


void SpeciesV::initCluster(Params& params)
{
    //Temporary BECK
    int ncells = 1;
    for (int iDim=0 ; iDim<nDim_particle ; iDim++)
        ncells *= (params.n_space[iDim]+1);
    bmax.resize(ncells,0);
    bmin.resize(ncells,0);
    species_loc_bmax.resize(ncells,0);

    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)

    //Primal dimension of fields. 
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;

    b_dim.resize(params.nDim_field, 1);
    if (nDim_particle == 1){
        b_dim[0] =  (1 + clrw) + 2 * oversize[0];
        f_dim1 = 1;
        f_dim2 = 1;
    }
    if (nDim_particle == 2){
        b_dim[0] =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] =  f_dim1;
        f_dim2 = 1;
    }
    if (nDim_particle == 3){
        b_dim[0] =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] = f_dim1;
        b_dim[2] = f_dim2;
    }

    //Initialize specMPI
    MPIbuff.allocate(nDim_particle);

    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_lost = 0.;
    nrj_new_particles = 0.;

}//END initCluster


void SpeciesV::dynamics(double time_dual, unsigned int ispec,
                       ElectroMagn* EMfields, Interpolator* Interp,
                       Projector* Proj, Params &params, bool diag_flag,
                       PartWalls* partWalls,
                       Patch* patch, SmileiMPI* smpi,
                       RadiationTables & RadiationTables,
                       MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                       vector<Diagnostic*>& localDiags)
{
    int ithread;
    #ifdef _OPENMP
        ithread = omp_get_thread_num();
    #else
        ithread = 0;
    #endif

    unsigned int iPart;

    // Reset list of particles to exchange
    clearExchList();

    int tid(0);
    double ener_iPart(0.);
    std::vector<double> nrj_lost_per_thd(1, 0.);

    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if (time_dual>time_frozen) { // moving particle

        smpi->dynamics_resize(ithread, nDim_particle, bmax.back());

        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);


        // Interpolate the fields at the particle position
        for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
            (*Interp)(EMfields, *particles, smpi, &(bmin[scell]), &(bmax[scell]), ithread );

        for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin++) {

            // Ionization
            if (Ionize)
                (*Ionize)(particles, bmin[ibin], bmax[ibin], Epart, EMfields, Proj);

            // Radiation losses
            if (Radiate)
            {

                // Radiation process
                (*Radiate)(*particles, this->photon_species, smpi,
                         RadiationTables,
                         bmin[ibin], bmax[ibin], ithread );

                // Update scalar variable for diagnostics
                nrj_radiation += (*Radiate).getRadiatedEnergy();

                // Update the quantum parameter chi
                (*Radiate).compute_thread_chipa(*particles,
                                                smpi,
                                                bmin[ibin],
                                                bmax[ibin],
                                                ithread );
            }

            // Multiphoton Breit-Wheeler
            if (Multiphoton_Breit_Wheeler_process)
            {

                // Pair generation process
                (*Multiphoton_Breit_Wheeler_process)(*particles,
                         smpi,
                         MultiphotonBreitWheelerTables,
                         bmin[ibin], bmax[ibin], ithread );

                 // Update scalar variable for diagnostics
                 // We reuse nrj_radiation for the pairs
                 nrj_radiation += (*Multiphoton_Breit_Wheeler_process).getPairEnergy();

                 // Update the photon quantum parameter chi of all photons
                 (*Multiphoton_Breit_Wheeler_process).compute_thread_chiph(*particles,
                                                 smpi,
                                                 bmin[ibin],
                                                 bmax[ibin],
                                                 ithread );

                 // Suppression of the decayed photons into pairs
                 (*Multiphoton_Breit_Wheeler_process).decayed_photon_cleaning(
                                 *particles,ibin, bmin.size(), &bmin[0], &bmax[0]);

            }
        }

        // Push the particles and the photons
        (*Push)(*particles, smpi, 0, bmax[bmax.size()-1], ithread );
        //particles->test_move( bmin[ibin], bmax[ibin], params );

        //Prepare for sorting
        for (unsigned int i=0; i<species_loc_bmax.size(); i++)
            species_loc_bmax[i] = 0;

        for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin++) {
            // Apply wall and boundary conditions
            if (mass>0)
            {
                for(unsigned int iwall=0; iwall<partWalls->size(); iwall++) {
                    for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if ( !(*partWalls)[iwall]->apply(*particles, iPart, this, dtgf, ener_iPart)) {
                            nrj_lost_per_thd[tid] += mass * ener_iPart;
                        }
                    }
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {
                    if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += mass * ener_iPart;
                        (*particles).cell_keys[iPart] = -1;
                    }
                    else {
                        //First reduction of the count sort algorithm. Lost particles are not included.
                        species_loc_bmax[(*particles).cell_keys[iPart]] ++; //First reduction of the count sort algorithm. Lost particles are not included.
                    }
                 }


            } else if (mass==0) {
                for(unsigned int iwall=0; iwall<partWalls->size(); iwall++) {
                    for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if ( !(*partWalls)[iwall]->apply(*particles, iPart, this, dtgf, ener_iPart)) {
                                nrj_lost_per_thd[tid] += ener_iPart;
                        }
                    }
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {
                    if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += ener_iPart;
                        (*particles).cell_keys[iPart] = -1;
                    }
                    else {
                        //First reduction of the count sort algorithm. Lost particles are not included.
                        species_loc_bmax[(*particles).cell_keys[iPart]] ++; //First reduction of the count sort algorithm. Lost particles are not included.

                    }
                }

            }
        }
        //START EXCHANGE PARTICLES OF THE CURRENT BIN ?


        // Project currents if not a Test species and charges as well if a diag is needed.
        // Do not project if a photon
        if ((!particles->is_test) && (mass > 0))
            for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
                (*Proj)(EMfields, *particles, smpi, bmin[scell], bmax[scell], ithread, scell, clrw, diag_flag, b_dim, ispec );


        for (unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++)
            nrj_bc_lost += nrj_lost_per_thd[tid];

//        // Add the ionized electrons to the electron species
//        if (Ionize)
//            electron_species->importParticles( params, patch, Ionize->new_electrons, localDiags );
//
//        // Radiation losses
//        if (Radiate)
//        {
//            // If creation of macro-photon, we add them to photon_species
//            if (photon_species)
//            {
//                photon_species->importParticles(params,
//                                                patch,
//                                                Radiate->new_photons,
//                                                localDiags);
//            }
//        }
//
//        // Multiphoton Breit-Wheeler
//        if (Multiphoton_Breit_Wheeler_process)
//        {
//
//            // Addition of the electron-positron particles
//            for (int k=0; k<2; k++) {
//                mBW_pair_species[k]->importParticles(params,
//                                             patch,
//                                             Multiphoton_Breit_Wheeler_process->new_pair[k],
//                                             localDiags);
//            }
//        }

    }
    else { // immobile particle (at the moment only project density)
        if ( diag_flag &&(!particles->is_test)){
            double* b_rho=nullptr;
            for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin ++) { //Loop for projection on buffer_proj

                if (nDim_field==2)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(ibin*clrw*f_dim1) : &(*EMfields->rho_)(ibin*clrw*f_dim1) ;
                if (nDim_field==3)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(ibin*clrw*f_dim1*f_dim2) : &(*EMfields->rho_)(ibin*clrw*f_dim1*f_dim2) ;
                else if (nDim_field==1)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(ibin*clrw) : &(*EMfields->rho_)(ibin*clrw) ;
                for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {
                    (*Proj)(b_rho, (*particles), iPart, ibin*clrw, b_dim);
                } //End loop on particles
            }//End loop on bins

        }
    }//END if time vs. time_frozen

}//END dynamic


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - increment the charge (projection)
//   - used at initialisation for Poisson (and diags if required, not for now dynamics )
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesV::computeCharge(unsigned int ispec, ElectroMagn* EMfields, Projector* Proj)
{
    // -------------------------------
    // calculate the particle charge
    // -------------------------------
    if ( (!particles->is_test) ) {
        double* b_rho=&(*EMfields->rho_)(0);

        for (unsigned int iPart=bmin[0] ; (int)iPart<bmax[bmax.size()-1]; iPart++ )
            (*Proj)(b_rho, (*particles), iPart, 0, b_dim);

    }

}//END computeCharge


// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesV::sort_part(Params &params)
{
    unsigned int ip, ip0, npart, ixy, ncell, token;
    int ii, ip_dest;
    int IX;
    double X;
    unsigned int m, length[3];
    vector<int> indices, buf_cell_keys[3][2];
    std::vector<unsigned int> cycle;
    unsigned int ip_src, icell;

    //split cell into smaller sub_cells for refined sorting
    ncell = (params.n_space[0]+1);
    for ( unsigned int i=1; i < params.nDim_field; i++) ncell *= (params.n_space[i]+1);

    //particles is a pointer pointing either particles_sorted[0] or [1].
    //token determines which one is currently pointed at and point at the other.
    token = (particles == &particles_sorted[0]);

    indices.resize(ncell);
    npart = (*particles).size(); //Number of particles before exchange

    length[0]=0;
    length[1]=params.n_space[1]+1;
    length[2]=params.n_space[2]+1;

    //species_loc_bmax stores the # of particles in a given cell quadrant.

    //Loop over just arrived particles
    for (unsigned int idim=0; idim < nDim_particle ; idim++){
        for (unsigned int ineighbor=0 ; ineighbor < 2 ; ineighbor++){
            buf_cell_keys[idim][ineighbor].resize( MPIbuff.part_index_recv_sz[idim][ineighbor]);
            #pragma omp simd 
            for (ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++){
                for (unsigned int ipos=0; ipos < nDim_particle ; ipos++) {
                    X = MPIbuff.partRecv[idim][ineighbor].position(ipos,ip)-min_loc_vec[ipos];
                    IX = round(X * dx_inv_[ipos] );
                    buf_cell_keys[idim][ineighbor][ip] = buf_cell_keys[idim][ineighbor][ip] * length[ipos] + IX;
                }


            }
            //Can we vectorize this reduction ?
            for (ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++){
                species_loc_bmax[buf_cell_keys[idim][ineighbor][ip]] ++;
            }
        }
    }


    // second loop convert the count array in cumulative sum
    indices[0]=0;
    bmin[0]=0;
    ii = 0;
    for (unsigned int ic=1; ic < ncell; ic++)
        {
            indices[ic] = indices[ic-1] + species_loc_bmax[ic-1];
            bmin[ic] = indices[ic];
            ii += abs(bmin[ic]-bmax[ic-1]);//Measure movement with respect to previous time step
            bmax[ic-1]= indices[ic];
        }
    bmax[ncell-1] = bmax[ncell-2] + species_loc_bmax.back() ; //New total number of particles is stored as last element of bmax
    //bmax[0] = bmax.back();           //For legacy, total number of particle is also stored in bmax[0], assuming a single cluster per patch. 



    if ( ii > 0.8*npart ) { //Full copy on a new buffer 
        particles_sorted[token].initialize( bmax.back(), *particles );
        // Copy particles in the new particle buffer (former version of the count sort)
        ip0 = 0;
        addPartInExchList(npart);//In order to stop loop at npart
        // last loop puts the particles and update the count array
        for (unsigned int j=0; j < indexes_of_particles_to_exchange.size(); j++){
            for (ip=ip0; ip < indexes_of_particles_to_exchange[j] ; ip++){
                ixy = (*particles).cell_keys[ip];
                (*particles).overwrite_part(ip, particles_sorted[token] , indices[ixy]);
                indices[ixy]++;
            }
            ip0 = indexes_of_particles_to_exchange[j] + 1;
        }
        indexes_of_particles_to_exchange.pop_back();
        for (unsigned int idim=0; idim < nDim_particle ; idim++){
            for (unsigned int ineighbor=0 ; ineighbor < 2 ; ineighbor++){
                for (ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++){
                    ixy = buf_cell_keys[idim][ineighbor][ip];
                    MPIbuff.partRecv[idim][ineighbor].overwrite_part(ip, particles_sorted[token] , indices[ixy]);
                    indices[ixy] ++;
                }
            }
        }
        //Point towards the newly created and sorted vector.
        particles = &particles_sorted[token] ;


    } else { //Copy all particles in the real particles array (ip < bmax) then performs a compact in place sort.
                    if (MPIbuff.partRecv[0][0].size() == 0) MPIbuff.partRecv[0][0].initialize(0, *particles); //Is this correct ?

                    // Resize the particle vector 
                    if (bmax.back() > npart){
                        (*particles).resize(bmax.back(), nDim_particle);
                        (*particles).cell_keys.resize(bmax.back(),-1); // Merge this in particles.resize(..) ?
                        for (unsigned int ipart = npart; ipart < bmax.back(); ipart ++) addPartInExchList(ipart);
                    }

                    //Copy all particles from MPI buffers back to the writable particles (marked by the index to exchange)
                    ip0 =0;
                    indexes_of_particles_to_exchange.push_back(bmax.back());//While loop stopper
                    for (unsigned int idim=0; idim < nDim_particle ; idim++){
                        for (unsigned int ineighbor=0 ; ineighbor < 2 ; ineighbor++){
                            for (ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++){
                                MPIbuff.partRecv[idim][ineighbor].overwrite_part(ip, *particles , indexes_of_particles_to_exchange[ip0]);
                                (*particles).cell_keys[indexes_of_particles_to_exchange[ip0]] = buf_cell_keys[idim][ineighbor][ip];
                                ip0++;
                            }
                        }
                    }
                    //Copy valid particles siting over bmax.back() back into the real particles array (happens when more particles are lost than received)
                    unsigned int real_part = bmax.back();
                    while (indexes_of_particles_to_exchange[ip0] < bmax.back()) { //This index is increasing. If at some point it is > bmax it means that all the other will be too.
                        while ((*particles).cell_keys[real_part] == -1) real_part++;
                        (*particles).overwrite_part(real_part, indexes_of_particles_to_exchange[ip0]);
                        (*particles).cell_keys[indexes_of_particles_to_exchange[ip0]] = (*particles).cell_keys[real_part];
                        real_part ++;
                        ip0++;
                    }

                    // Resize the particle vector 
                    if (bmax.back() < npart){
                        (*particles).resize(bmax.back(), nDim_particle);
                        (*particles).cell_keys.resize(bmax.back()); // Merge this in particles.resize(..) ?
                    }

                    //Generates indices at which particles must be copied
                    for (unsigned int iicell=0; iicell < ncell; iicell ++)
                        indices[iicell] = bmin[iicell];

                    icell = 0;//Stores the current number of the cell
                    //Loop over all cells
                    for (icell = 0 ; icell < ncell; icell++){
                        for (ip=indices[icell]; ip < bmax[icell] ; ip++){
                            //update value of current cell 'icell' if necessary 
                            //if particle changes cell, build a cycle of exchange as long as possible. Treats all particles
                            if ((*particles).cell_keys[ip] != icell ){
                                cycle.resize(1);
                                cycle[0] = ip;
                                ip_src = ip;
                                //While the destination particle is not going out of the patch or back to the initial cell, keep building the cycle.
                                while ( (*particles).cell_keys[ip_src] != icell) {
                                   //Scan the next cell destination
                                    ip_dest = indices[(*particles).cell_keys[ip_src]];
                                    while ( (*particles).cell_keys[ip_dest] == (*particles).cell_keys[ip_src] ) ip_dest++;
                                    //In the destination cell, if a particle is going out of this cell, add it to the cycle.
                                    indices[(*particles).cell_keys[ip_src]] = ip_dest + 1 ;
                                    cycle.push_back(ip_dest);
                                    ip_src = ip_dest; //Destination becomes source for the next iteration
                                }
                                    //swap parts
                                    (*particles).swap_parts(cycle);
                                    //Must also swap the cell_keys to mark particles as treated
                                    for (unsigned int ipart = cycle.size()-1; ipart > 0; ipart--) (*particles).cell_keys[cycle[ipart]] = (*particles).cell_keys[cycle[ipart-1]];
                                    (*particles).cell_keys[ip] = icell;

                            }
                        }
                    } //end loop on cells
    }

}


void SpeciesV::compute_part_cell_keys(Params &params)
{
    //Compute part_cell_keys at patch creation. This operation is normally done in the pusher to avoid additional particles pass.  

    unsigned int ip, npart, ixy;
    int IX;
    double X;
    unsigned int length[3];

    npart = (*particles).size(); //Number of particles before exchange
    length[0]=0;
    length[1]=params.n_space[1]+1;
    length[2]=params.n_space[2]+1;

    #pragma omp simd 
    for (ip=0; ip < npart ; ip++){
    // Counts the # of particles in each cell (or sub_cell) and store it in sbmax.
        for (unsigned int ipos=0; ipos < nDim_particle ; ipos++) {
            X = (*particles).position(ipos,ip)-min_loc_vec[ipos];
            IX = round(X * dx_inv_[ipos] );
            (*particles).cell_keys[ip] = (*particles).cell_keys[ip] * length[ipos] + IX;
        }
    }
    for (ip=0; ip < npart ; ip++)
        species_loc_bmax[(*particles).cell_keys[ip]] ++ ;

}
