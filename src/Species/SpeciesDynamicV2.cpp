#include "SpeciesDynamicV2.h"

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

#include "DiagnosticTrack.h"

#include "SpeciesMetrics.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
SpeciesDynamicV2::SpeciesDynamicV2(Params& params, Patch* patch) :
    Species(params, patch)
{
    initCluster( params );

}//END SpeciesDynamicV2 creator

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
SpeciesDynamicV2::~SpeciesDynamicV2()
{
}


void SpeciesDynamicV2::initCluster(Params& params)
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

void SpeciesDynamicV2::resizeCluster(Params& params)
{

}// end resizeCluster

void SpeciesDynamicV2::dynamics(double time_dual, unsigned int ispec,
                       ElectroMagn* EMfields, Params &params, bool diag_flag,
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

        //smpi->dynamics_resize(ithread, nDim_particle, bmax.back());

        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);

        int npack    =  f_dim0-2*oversize[0];
        int packsize = (f_dim1-2*oversize[1]);
        if (nDim_particle == 3)
            packsize *= (f_dim2-2*oversize[2]);

        //Prepare for sorting
        for (unsigned int i=0; i<species_loc_bmax.size(); i++)
            species_loc_bmax[i] = 0;

        // Resize Cell_keys
        // Need to check if this is necessary here
        // (*particles).cell_keys.resize((*particles).size());

        for ( int ipack = 0 ; ipack < npack ; ipack++ ) {

            int nparts_in_pack = bmax[ (ipack+1) * packsize-1 ];
            smpi->dynamics_resize(ithread, nDim_particle, nparts_in_pack );

            // Interpolate the fields at the particle position
            //for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
            //    (*Interp)(EMfields, *particles, smpi, &(bmin[scell]), &(bmax[scell]), ithread );
            for (unsigned int scell = 0 ; scell < packsize ; scell++)
                (*Interp)(EMfields, *particles, smpi, &(bmin[ipack*packsize+scell]), &(bmax[ipack*packsize+scell]), ithread, bmin[ipack*packsize] );


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
            //(*Push)(*particles, smpi, 0, bmax[bmax.size()-1], ithread );
            (*Push)(*particles, smpi, bmin[ipack*packsize], bmax[ipack*packsize+packsize-1], ithread, bmin[ipack*packsize] );

            // Computation of the particle cell keys for all particles
            // this->compute_bin_cell_keys(params, bmin[ipack*packsize], bmax[ipack*packsize+packsize-1]);

            //for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin++) {
            for (unsigned int ibin = 0 ; ibin < packsize ; ibin++) {
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
                        //for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {

                        for (iPart=bmin[ipack*packsize+ibin] ; (int)iPart<bmax[ipack*packsize+ibin]; iPart++ ) {
                            if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                                addPartInExchList( iPart );
                                nrj_lost_per_thd[tid] += mass * ener_iPart;
                                (*particles).cell_keys[iPart] = -1;
                            }
                            else {
                                //Compute cell_keys of remaining particles
                                for ( int i = 0 ; i<nDim_particle; i++ ){
                                    (*particles).cell_keys[iPart] *= this->length[i];
                                    (*particles).cell_keys[iPart] += round( ((*particles).position(i,iPart)-min_loc_vec[i]+0.00000000000001) * dx_inv_[i] );
                                }
                                //First reduction of the count sort algorithm. Lost particles are not included.
                                species_loc_bmax[(*particles).cell_keys[iPart]] ++;
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
                            //Compute cell_keys of remaining particles
                            for ( int i = 0 ; i<nDim_particle; i++ ){
                                (*particles).cell_keys[iPart] *= this->length[i];
                                (*particles).cell_keys[iPart] += round( ((*particles).position(i,iPart)-min_loc_vec[i]+0.00000000000001) * dx_inv_[i] );
                            }
                            //First reduction of the count sort algorithm. Lost particles are not included.
                            species_loc_bmax[(*particles).cell_keys[iPart]] ++;

                        }
                    }

                }
            }
            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?


            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if ((!particles->is_test) && (mass > 0))
                //for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
                //    (*Proj)(EMfields, *particles, smpi, bmin[scell], bmax[scell], ithread, scell, clrw, diag_flag, params.is_spectral, b_dim, ispec );
                for (unsigned int scell = 0 ; scell < packsize ; scell++)
                    (*Proj)(EMfields, *particles, smpi, bmin[ipack*packsize+scell], bmax[ipack*packsize+scell], ithread, ipack*packsize+scell, clrw, diag_flag, params.is_spectral, b_dim, ispec, bmin[ipack*packsize] );

            for (unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++)
                nrj_bc_lost += nrj_lost_per_thd[tid];

        }

    }
    else { // immobile particle (at the moment only project density)
        if ( diag_flag &&(!particles->is_test)){
            double* b_rho=nullptr;
            for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin ++) { //Loop for projection on buffer_proj

                if (nDim_field==2)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(0*clrw*f_dim1) : &(*EMfields->rho_)(0*clrw*f_dim1) ;
                if (nDim_field==3)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(0*clrw*f_dim1*f_dim2) : &(*EMfields->rho_)(0*clrw*f_dim1*f_dim2) ;
                else if (nDim_field==1)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(0*clrw) : &(*EMfields->rho_)(0*clrw) ;
                for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {
                    (*Proj)(b_rho, (*particles), iPart, 0*clrw, b_dim);
                } //End loop on particles
            }//End loop on bins

        }
    }//END if time vs. time_frozen

}//END dynamic


void SpeciesDynamicV2::scalar_dynamics(double time_dual, unsigned int ispec,
                       ElectroMagn* EMfields, Params &params, bool diag_flag,
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
    if (time_dual>time_frozen)
    { // moving particle

        smpi->dynamics_resize(ithread, nDim_particle, bmax.back());

        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &(smpi->dynamics_Epart[ithread]);

        //Prepare for sorting
        for (unsigned int i=0; i<species_loc_bmax.size(); i++)
            species_loc_bmax[i] = 0;

        // Resize Cell_keys
        // Need to check if this is necessary here
        // (*particles).cell_keys.resize((*particles).size());

        // Interpolate the fields at the particle position
        //for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
        //    (*Interp)(EMfields, *particles, smpi, &(bmin[scell]), &(bmax[scell]), ithread );
        for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
        {
            (*Interp)(EMfields, *particles, smpi, &(bmin[scell]), &(bmax[scell]), ithread);

            // Ionization
            if (Ionize)
                (*Ionize)(particles, bmin[scell], bmax[scell], Epart, EMfields, Proj);

            // Radiation losses
            if (Radiate)
            {

                // Radiation process
                (*Radiate)(*particles, this->photon_species, smpi,
                           RadiationTables,
                           bmin[scell], bmax[scell], ithread );

                // Update scalar variable for diagnostics
                nrj_radiation += (*Radiate).getRadiatedEnergy();

                // Update the quantum parameter chi
                (*Radiate).compute_thread_chipa(*particles,
                                                smpi,
                                                bmin[scell],
                                                bmax[scell],
                                                ithread );
            }

            // Multiphoton Breit-Wheeler
            if (Multiphoton_Breit_Wheeler_process)
            {

                // Pair generation process
                (*Multiphoton_Breit_Wheeler_process)(*particles,
                                                     smpi,
                                                     MultiphotonBreitWheelerTables,
                                                     bmin[scell], bmax[scell], ithread );

                // Update scalar variable for diagnostics
                // We reuse nrj_radiation for the pairs
                nrj_radiation += (*Multiphoton_Breit_Wheeler_process).getPairEnergy();

                // Update the photon quantum parameter chi of all photons
                (*Multiphoton_Breit_Wheeler_process).compute_thread_chiph(*particles,
                                                                          smpi,
                                                                          bmin[scell],
                                                                          bmax[scell],
                                                                          ithread );

                // Suppression of the decayed photons into pairs
                (*Multiphoton_Breit_Wheeler_process).decayed_photon_cleaning(
                                                                             *particles,scell, bmin.size(), &bmin[0], &bmax[0]);

            }
        }

        // Push the particles and the photons
        (*Push)(*particles, smpi, 0, bmax.back(), ithread, 0.);

        // Computation of the particle cell keys for all particles
        // this->compute_bin_cell_keys(params,0, bmax.back());

        for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
        {
            // Apply wall and boundary conditions
            if (mass>0)
            {
                for(unsigned int iwall=0; iwall<partWalls->size(); iwall++) {
                    for (iPart=bmin[scell] ; (int)iPart<bmax[scell]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if ( !(*partWalls)[iwall]->apply(*particles, iPart, this, dtgf, ener_iPart)) {
                            nrj_lost_per_thd[tid] += mass * ener_iPart;
                        }
                    }
                }

                for (iPart=bmin[scell] ; (int)iPart<bmax[scell]; iPart++ ) {
                    if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += mass * ener_iPart;
                        (*particles).cell_keys[iPart] = -1;
                    }
                    else {
                        //Compute cell_keys of remaining particles
                        for ( int i = 0 ; i<nDim_particle; i++ ){
                            (*particles).cell_keys[iPart] *= this->length[i];
                            (*particles).cell_keys[iPart] += round( ((*particles).position(i,iPart)-min_loc_vec[i]+0.00000000000001) * dx_inv_[i] );
                        }
                        //First reduction of the count sort algorithm. Lost particles are not included.
                        species_loc_bmax[(*particles).cell_keys[iPart]] ++;
                    }
                }

            } else if (mass==0) {
                for(unsigned int iwall=0; iwall<partWalls->size(); iwall++) {
                    for (iPart=bmin[scell] ; (int)iPart<bmax[scell]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if ( !(*partWalls)[iwall]->apply(*particles, iPart, this, dtgf, ener_iPart)) {
                            nrj_lost_per_thd[tid] += ener_iPart;
                        }
                    }
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for (iPart=bmin[scell] ; (int)iPart<bmax[scell]; iPart++ ) {
                    if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += ener_iPart;
                        (*particles).cell_keys[iPart] = -1;
                    }
                    else {
                        //Compute cell_keys of remaining particles
                        for ( int i = 0 ; i<nDim_particle; i++ ){
                            (*particles).cell_keys[iPart] *= this->length[i];
                            (*particles).cell_keys[iPart] += round( ((*particles).position(i,iPart)-min_loc_vec[i]+0.00000000000001) * dx_inv_[i] );
                        }
                        //First reduction of the count sort algorithm. Lost particles are not included.
                        species_loc_bmax[(*particles).cell_keys[iPart]] ++;
                    }
                }
            } // end if mass > 0

            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if ((!particles->is_test) && (mass > 0))
                (*Proj)(EMfields, *particles, smpi, bmin[scell],
                                                    bmax[scell],
                                                    ithread, scell,
                                                    0, diag_flag,
                                                    params.is_spectral,
                                                    b_dim, ispec);

        } // end loop on cells

        for (unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++)
            nrj_bc_lost += nrj_lost_per_thd[tid];

    }
    else { // immobile particle (at the moment only project density)
        if ( diag_flag &&(!particles->is_test)){
            double* b_rho=nullptr;
            for (unsigned int scell = 0 ; scell < bmin.size() ; scell ++) { //Loop for projection on buffer_proj

                if (nDim_field==2)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(0*clrw*f_dim1) : &(*EMfields->rho_)(0*clrw*f_dim1) ;
                if (nDim_field==3)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(0*clrw*f_dim1*f_dim2) : &(*EMfields->rho_)(0*clrw*f_dim1*f_dim2) ;
                else if (nDim_field==1)
                    b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(0*clrw) : &(*EMfields->rho_)(0*clrw) ;
                for (iPart=bmin[scell] ; (int)iPart<bmax[scell]; iPart++ ) {
                    (*Proj)(b_rho, (*particles), iPart, 0*clrw, b_dim);
                } //End loop on particles
            }//End loop on bins

        }
    }//END if time vs. time_frozen

}//END scalar_dynamics


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - increment the charge (projection)
//   - used at initialisation for Poisson (and diags if required, not for now dynamics )
// ---------------------------------------------------------------------------------------------------------------------
void SpeciesDynamicV2::computeCharge(unsigned int ispec, ElectroMagn* EMfields)
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
void SpeciesDynamicV2::sort_part(Params &params)
{

    unsigned int npart, ncell;
    int ip_dest, cell_target;
    unsigned int length[3];
    vector<int> buf_cell_keys[3][2];
    std::vector<unsigned int> cycle;
    unsigned int ip_src;

    //split cell into smaller sub_cells for refined sorting
    ncell = (params.n_space[0]+1);
    for ( unsigned int i=1; i < params.nDim_field; i++) ncell *= (params.n_space[i]+1);

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
            for (unsigned int ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++){
                for (unsigned int ipos=0; ipos < nDim_particle ; ipos++) {
                    double X = MPIbuff.partRecv[idim][ineighbor].position(ipos,ip)-min_loc_vec[ipos]+0.00000000000001;
                    int IX = round(X * dx_inv_[ipos] );
                    buf_cell_keys[idim][ineighbor][ip] = buf_cell_keys[idim][ineighbor][ip] * length[ipos] + IX;
                }
            }
            //Can we vectorize this reduction ?
            for (unsigned int ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++){
                species_loc_bmax[buf_cell_keys[idim][ineighbor][ip]] ++;
            }
        }
    }

    // second loop convert the count array in cumulative sum
    bmin[0]=0;
    for (unsigned int ic=1; ic < ncell; ic++)
        {
            bmin[ic] = bmin[ic-1] + species_loc_bmax[ic-1];
            bmax[ic-1]= bmin[ic];
        }
    bmax[ncell-1] = bmax[ncell-2] + species_loc_bmax.back() ; //New total number of particles is stored as last element of bmax

    //Now proceed to the cycle sort

    if (MPIbuff.partRecv[0][0].size() == 0) MPIbuff.partRecv[0][0].initialize(0, *particles); //Is this correct ?

    // Resize the particle vector
    if (bmax.back() > npart){
        (*particles).resize(bmax.back(), nDim_particle);
        (*particles).cell_keys.resize(bmax.back(),-1); // Merge this in particles.resize(..) ?
        for (unsigned int ipart = npart; ipart < bmax.back(); ipart ++) addPartInExchList(ipart);
    }

    //Copy all particles from MPI buffers back to the writable particles via cycle sort pass.
    for (unsigned int idim=0; idim < nDim_particle ; idim++){
        for (unsigned int ineighbor=0 ; ineighbor < 2 ; ineighbor++){
            for (unsigned int ip=0; ip < MPIbuff.part_index_recv_sz[idim][ineighbor]; ip++){
                cycle.resize(1);
                cell_target = buf_cell_keys[idim][ineighbor][ip];
                ip_dest = bmin[cell_target];
                while ( (*particles).cell_keys[ip_dest] == cell_target ) ip_dest++;
                bmin[cell_target] = ip_dest + 1 ;
                cycle[0] = ip_dest;
                cell_target = (*particles).cell_keys[ip_dest];
                //As long as the particle is not erased, we can build up the cycle.
                while (cell_target != -1){
                    ip_dest = bmin[cell_target];
                    while ( (*particles).cell_keys[ip_dest] == cell_target ) ip_dest++;
                    bmin[cell_target] = ip_dest + 1 ;
                    cycle.push_back(ip_dest);
                    cell_target = (*particles).cell_keys[ip_dest];
                }
                //Last target_cell is -1, the particle must be erased:
                (*particles).translate_parts(cycle);
                //Eventually copy particle from the MPI buffer into the particle vector.
                MPIbuff.partRecv[idim][ineighbor].overwrite_part(ip, *particles , cycle[0]);
            }
        }
    }

    //Copy valid particles siting over bmax.back() back into the real particles array (happens when more particles are lost than received)
    for (unsigned int ip=bmax.back(); ip < npart; ip++){
        cell_target = (*particles).cell_keys[ip];
        if(cell_target == -1) continue;
        cycle.resize(0);
        cycle.push_back(ip);
        //As long as the particle is not erased, we can build up the cycle.
        while (cell_target != -1){
            ip_dest = bmin[cell_target];
            while ( (*particles).cell_keys[ip_dest] == cell_target ) ip_dest++;
            bmin[cell_target] = ip_dest + 1 ;
            cycle.push_back(ip_dest);
            cell_target = (*particles).cell_keys[ip_dest];
        }
        //Last target_cell is -1, the particle must be erased:
        (*particles).translate_parts(cycle);
    }

    // Resize the particle vector
    if (bmax.back() < npart){
        (*particles).resize(bmax.back(), nDim_particle);
        (*particles).cell_keys.resize(bmax.back()); // Merge this in particles.resize(..) ?
    }


    //Loop over all cells
    for (unsigned int icell = 0 ; icell < ncell; icell++){
        for (unsigned int ip=bmin[icell]; ip < bmax[icell] ; ip++){
            //update value of current cell 'icell' if necessary
            //if particle changes cell, build a cycle of exchange as long as possible. Treats all particles
            if ((*particles).cell_keys[ip] != icell ){
                cycle.resize(1);
                cycle[0] = ip;
                ip_src = ip;
                //While the destination particle is not going out of the patch or back to the initial cell, keep building the cycle.
                while ( (*particles).cell_keys[ip_src] != icell) {
                   //Scan the next cell destination
                    ip_dest = bmin[(*particles).cell_keys[ip_src]];
                    while ( (*particles).cell_keys[ip_dest] == (*particles).cell_keys[ip_src] ) ip_dest++;
                    //In the destination cell, if a particle is going out of this cell, add it to the cycle.
                    bmin[(*particles).cell_keys[ip_src]] = ip_dest + 1 ;
                    cycle.push_back(ip_dest);
                    ip_src = ip_dest; //Destination becomes source for the next iteration
                }
                    //swap parts
                    (*particles).swap_parts(cycle);
            }
        }
    } //end loop on cells
    // Restore bmin initial value
    bmin[0]=0;
    for (unsigned int ic=1; ic < ncell; ic++)
        bmin[ic] = bmax[ic-1];

}

// -----------------------------------------------------------------------------
//! Compute part_cell_keys at patch creation.
//! This operation is normally done in the pusher to avoid additional particles pass.
// -----------------------------------------------------------------------------
void SpeciesDynamicV2::compute_part_cell_keys(Params &params)
{

    unsigned int ip, nparts, ixy;
    int IX;
    double X;
    unsigned int length[3];

    //Number of particles before exchange
    nparts = (*particles).size();

    // Cell_keys is resized at the current number of particles
    (*particles).cell_keys.resize(nparts);

    length[0]=0;
    length[1]=params.n_space[1]+1;
    length[2]=params.n_space[2]+1;

    #pragma omp simd
    for (ip=0; ip < nparts ; ip++){
    // Counts the # of particles in each cell (or sub_cell) and store it in sbmax.
        for (unsigned int ipos=0; ipos < nDim_particle ; ipos++) {
            X = (*particles).position(ipos,ip)-min_loc_vec[ipos]+0.00000000000001;
            IX = round(X * dx_inv_[ipos] );
            (*particles).cell_keys[ip] = (*particles).cell_keys[ip] * length[ipos] + IX;
        }
    }

    // Reduction of the number of particles per cell in species_loc_bmax
    for (ip=0; ip < nparts ; ip++)
        species_loc_bmax[(*particles).cell_keys[ip]] ++ ;

}

// -----------------------------------------------------------------------------
//! Compute cell_keys for the specified bin boundaries.
//! params object that contains the global parameters
//! istart first bin index
//! iend last bin index
// -----------------------------------------------------------------------------
void SpeciesDynamicV2::compute_bin_cell_keys(Params &params, int istart, int iend)
{
    // Resize of cell_keys seems necessary here
    (*particles).cell_keys.resize((*particles).size());

    #pragma omp simd
    for (int ip=istart; ip < iend; ip++){
    // Counts the # of particles in each cell (or sub_cell) and store it in sbmax.
        for (unsigned int ipos=0; ipos < nDim_particle ; ipos++) {
            (*particles).cell_keys[ip] *= this->length[ipos];
            (*particles).cell_keys[ip] += round( ((*particles).position(ipos,ip)-min_loc_vec[ipos]+0.00000000000001) * dx_inv_[ipos] );
        }
    }
}

void SpeciesDynamicV2::importParticles( Params& params, Patch* patch, Particles& source_particles, vector<Diagnostic*>& localDiags )
{
    unsigned int npart = source_particles.size(), ibin, ii, nbin=bmin.size();
    double inv_cell_length = 1./ params.cell_length[0];

    // If this species is tracked, set the particle IDs
    if( particles->tracked )
        dynamic_cast<DiagnosticTrack*>(localDiags[tracking_diagnostic])->setIDs( source_particles );

    unsigned int length[3];
    length[0]=0;
    length[1]=params.n_space[1]+1;
    length[2]=params.n_space[2]+1;

    int IX;
    double X;

    // Move particles
    for( unsigned int i=0; i<npart; i++ ) {

        ibin = 0;
        for (unsigned int ipos=0; ipos < nDim_particle ; ipos++) {
            X = source_particles.position(ipos,i)-min_loc_vec[ipos]+0.00000000000001;
            IX = round(X * dx_inv_[ipos] );
            ibin = ibin * length[ipos] + IX;
        }

        // Copy particle to the correct bin
        source_particles.cp_particle(i, *particles, bmax[ibin] );

        // Update the bin counts
        bmax[ibin]++;
        for (ii=ibin+1; ii<nbin; ii++) {
            bmin[ii]++;
            bmax[ii]++;
        }

        particles->cell_keys.insert( particles->cell_keys.begin() + bmin[ibin] + species_loc_bmax[ibin], ibin);
        species_loc_bmax[ibin] ++ ;

    }

    source_particles.clear();
}

// -----------------------------------------------------------------------------
//! This function reconfigures the type of species according
//! to the vectorization mode
// -----------------------------------------------------------------------------
void SpeciesDynamicV2::reconfiguration(Params &params, Patch * patch)
{

    //unsigned int ncell;
    bool reasign_operators = false;
    //float ratio_number_of_vecto_cells;
    float vecto_time = 0.;
    float scalar_time = 0.;

    //split cell into smaller sub_cells for refined sorting
    // cell = (params.n_space[0]+1);
    //for ( unsigned int i=1; i < params.nDim_field; i++) ncell *= (params.n_space[i]+1);

    // --------------------------------------------------------------------
    // Metrics 1 - based on the ratio of vectorized cells
    // Compute the number of cells that contain more than 8 particles
    //ratio_number_of_vecto_cells = SpeciesMetrics::get_ratio_number_of_vecto_cells(species_loc_bmax,8);

    // Test metrics, if necessary we reasign operators
    //if ( (ratio_number_of_vecto_cells > 0.5 && this->vectorized_operators == false)
    //  || (ratio_number_of_vecto_cells < 0.5 && this->vectorized_operators == true))
    //{
    //    reasign_operators = true;
    //}
    // --------------------------------------------------------------------

    // --------------------------------------------------------------------
    // Metrics 2 - based on the evaluation of the computational time
    SpeciesMetrics::get_computation_time(species_loc_bmax,
                                        vecto_time,
                                        scalar_time);

    //std::cout << "vecto_time " << vecto_time << " " << scalar_time << '\n';

    if ( (vecto_time <= scalar_time && this->vectorized_operators == false)
      || (vecto_time > scalar_time && this->vectorized_operators == true))
    {
        reasign_operators = true;
    }
    // --------------------------------------------------------------------

    /*std::cout << "Vectorized_operators: " << this->vectorized_operators
              << " ratio_number_of_vecto_cells: " << this->ratio_number_of_vecto_cells
              << " number_of_vecto_cells: " << number_of_vecto_cells
              << " number_of_non_zero_cells: " << number_of_non_zero_cells
              << " ncells: " << ncell << "\n";*/

    // Operator reasignment if required by the metrics
    if (reasign_operators)
    {

        // The type of operator is changed
        this->vectorized_operators = !this->vectorized_operators;

        /*MESSAGE(1,"> Species " << this->name << " reconfiguration (" << this->vectorized_operators
                  << ") in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ")"
                  << " of MPI process "<< patch->MPI_me_);*/

        this->reconfigure_operators(params, patch);

    }

    /*std::cout << " bin number: " << bmin.size()
              << " nb particles: " << bmax[bmax.size()-1]
              << '\n';*/

}

// -----------------------------------------------------------------------------
//! This function reconfigures the type of species according
//! to the vectorization mode
// -----------------------------------------------------------------------------
void SpeciesDynamicV2::configuration(Params &params, Patch * patch)
{
    //float ratio_number_of_vecto_cells;
    float vecto_time = 0.;
    float scalar_time = 0.;

    //split cell into smaller sub_cells for refined sorting
    // cell = (params.n_space[0]+1);
    //for ( unsigned int i=1; i < params.nDim_field; i++) ncell *= (params.n_space[i]+1);

    // --------------------------------------------------------------------
    // Metrics 1 - based on the ratio of vectorized cells
    // Compute the number of cells that contain more than 8 particles
    //ratio_number_of_vecto_cells = SpeciesMetrics::get_ratio_number_of_vecto_cells(species_loc_bmax,8);


    // --------------------------------------------------------------------

    // --------------------------------------------------------------------
    // Metrics 2 - based on the evaluation of the computational time
    SpeciesMetrics::get_computation_time(species_loc_bmax,
                                        vecto_time,
                                        scalar_time);

    if (vecto_time <= scalar_time )
    {
        this->vectorized_operators = true;
    }
    else if (vecto_time > scalar_time)
    {
        this->vectorized_operators = false;
    }
    // --------------------------------------------------------------------

    /*std::cout << "Vectorized_operators: " << this->vectorized_operators
              << " ratio_number_of_vecto_cells: " << this->ratio_number_of_vecto_cells
              << " number_of_vecto_cells: " << number_of_vecto_cells
              << " number_of_non_zero_cells: " << number_of_non_zero_cells
              << " ncells: " << ncell << "\n";*/

    this->reconfigure_operators(params, patch);

}

// -----------------------------------------------------------------------------
//! This function reconfigures the operators
// -----------------------------------------------------------------------------
void SpeciesDynamicV2::reconfigure_operators(Params &params, Patch * patch)
{
    // Destroy current operators
    delete Interp;
    //delete Push;
    delete Proj;

    // Reassign the correct Interpolator
    Interp = InterpolatorFactory::create(params, patch, this->vectorized_operators);
    // Reassign the correct Pusher to Push
    //Push = PusherFactory::create(params, this);
    // Reassign the correct Projector
    Proj = ProjectorFactory::create(params, patch, this->vectorized_operators);
}
