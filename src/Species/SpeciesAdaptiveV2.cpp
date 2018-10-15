#include "SpeciesAdaptiveV2.h"

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
SpeciesAdaptiveV2::SpeciesAdaptiveV2(Params& params, Patch* patch) :
    SpeciesV(params, patch)
{
    initCluster( params );
    npack_ = 0 ;
    packsize_ = 0;
}//END SpeciesAdaptiveV2 creator

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
SpeciesAdaptiveV2::~SpeciesAdaptiveV2()
{
}

//! Method calculating the Particle dynamics (interpolation, pusher, projection)
//! without vectorized operators but with the cell sorting algorithm
void SpeciesAdaptiveV2::scalar_dynamics(double time_dual, unsigned int ispec,
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

#ifdef  __DETAILED_TIMERS
    double timer;
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

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif

        // Interpolate the fields at the particle position
        (*Interp)(EMfields, *particles, smpi, &(bmin[0]), &(bmax[bmax.size()-1]), ithread, bmin[0]);

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[0] += MPI_Wtime() - timer;
#endif

        // Interpolate the fields at the particle position
        //for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
        //    (*Interp)(EMfields, *particles, smpi, &(bmin[scell]), &(bmax[scell]), ithread );
        for (unsigned int scell = 0 ; scell < bmin.size() ; scell++)
        {

            // Ionization
            if (Ionize)
            {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
                (*Ionize)(particles, bmin[scell], bmax[scell], Epart, patch, Proj);
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[4] += MPI_Wtime() - timer;
#endif
            }


            // Radiation losses
            if (Radiate)
            {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
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
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[5] += MPI_Wtime() - timer;
#endif
            }

            // Multiphoton Breit-Wheeler
            if (Multiphoton_Breit_Wheeler_process)
            {
#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif
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
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[6] += MPI_Wtime() - timer;
#endif
            }
        }

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
        // Push the particles and the photons
        (*Push)(*particles, smpi, 0, bmax.back(), ithread, 0.);
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[1] += MPI_Wtime() - timer;
        timer = MPI_Wtime();
#endif

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
                        for ( unsigned int i = 0 ; i<nDim_particle; i++ ){
                            (*particles).cell_keys[iPart] *= this->length[i];
                            (*particles).cell_keys[iPart] += round( ((*particles).position(i,iPart)-min_loc_vec[i]) * dx_inv_[i] );
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
                        for ( unsigned int i = 0 ; i<nDim_particle; i++ ){
                            (*particles).cell_keys[iPart] *= this->length[i];
                            (*particles).cell_keys[iPart] += round( ((*particles).position(i,iPart)-min_loc_vec[i]) * dx_inv_[i] );
                        }
                        //First reduction of the count sort algorithm. Lost particles are not included.
                        species_loc_bmax[(*particles).cell_keys[iPart]] ++;
                    }
                }
            } // end if mass > 0
        } // end loop on cells

#ifdef  __DETAILED_TIMERS
        patch->patch_timers[3] += MPI_Wtime() - timer;
#endif

        // Project currents if not a Test species and charges as well if a diag is needed.
        // Do not project if a photon
        if ((!particles->is_test) && (mass > 0))
        {

#ifdef  __DETAILED_TIMERS
        timer = MPI_Wtime();
#endif
            (*Proj)(EMfields, *particles, smpi, bmin[0],
                                                bmax.back(),
                                                ithread, 0,
                                                0, diag_flag,
                                                params.is_spectral,
                                                b_dim, ispec);
#ifdef  __DETAILED_TIMERS
        patch->patch_timers[2] += MPI_Wtime() - timer;
#endif

        }

        for (unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++)
            nrj_bc_lost += nrj_lost_per_thd[tid];

    }
    else { // immobile particle (at the moment only project density)
        if ( diag_flag &&(!particles->is_test)){
            double* b_rho=nullptr;

            for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin ++) { //Loop for projection on buffer_proj

                b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(0) : &(*EMfields->rho_)(0) ;

                for (iPart=bmin[ibin] ; (int)iPart<bmax[ibin]; iPart++ ) {
                    (*Proj)(b_rho, (*particles), iPart, 0, b_dim);
                } //End loop on particles
            }//End loop on bins

        }
    }//END if time vs. time_frozen

}//END scalar_dynamics


// -----------------------------------------------------------------------------
//! Compute part_cell_keys at patch creation.
//! This operation is normally done in the pusher to avoid additional particles pass.
// -----------------------------------------------------------------------------
/*void SpeciesAdaptiveV2::compute_part_cell_keys(Params &params)
{

    unsigned int ip, nparts;
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
            X = (*particles).position(ipos,ip)-min_loc_vec[ipos];
            IX = round(X * dx_inv_[ipos] );
            (*particles).cell_keys[ip] = (*particles).cell_keys[ip] * length[ipos] + IX;
        }
    }

    // Reduction of the number of particles per cell in species_loc_bmax
    for (ip=0; ip < nparts ; ip++)
        species_loc_bmax[(*particles).cell_keys[ip]] ++ ;
}*/


// -----------------------------------------------------------------------------
//! This function reconfigures the type of species according
//! to the vectorization mode
// -----------------------------------------------------------------------------
void SpeciesAdaptiveV2::reconfiguration(Params &params, Patch * patch)
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
void SpeciesAdaptiveV2::configuration(Params &params, Patch * patch)
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

    if (vecto_time < scalar_time )
    {
        this->vectorized_operators = true;
    }
    else if (vecto_time > scalar_time)
    {
        this->vectorized_operators = false;
    }
    // Default mode where there is no particles
    else
    {
        this->vectorized_operators = (params.dynamic_default_mode == "on");
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
void SpeciesAdaptiveV2::reconfigure_operators(Params &params, Patch * patch)
{
    // Destroy current operators
    delete Interp;
    //delete Push;
    delete Proj;

    // Reassign the correct Interpolator
    this->Interp = InterpolatorFactory::create(params, patch, this->vectorized_operators);
    // Reassign the correct Pusher to Push
    //Push = PusherFactory::create(params, this);
    // Reassign the correct Projector
    this->Proj = ProjectorFactory::create(params, patch, this->vectorized_operators);
}
