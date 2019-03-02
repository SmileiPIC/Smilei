// ----------------------------------------------------------------------------
//! \file MerginngVranic.cpp
//
//! \brief Functions of the class MergingVranic
//! Particle merging with the method of Vranic et al.
//! Vranic CPC 191 65-73 (2015)
//
//! Creation - 01/2019 - Mathieu Lobet
//
// ----------------------------------------------------------------------------

#include "MergingVranic.h"

#include <cmath>

// -----------------------------------------------------------------------------
//! Constructor for RadiationNLandauLifshitz
//! Inherited from Radiation
// -----------------------------------------------------------------------------
MergingVranic::MergingVranic(Params& params,
                             Species * species)
      : Merging(params, species)
{
}

// -----------------------------------------------------------------------------
//! Destructor for MergingVranic
// -----------------------------------------------------------------------------
MergingVranic::~MergingVranic()
{
}


// ---------------------------------------------------------------------
//! Overloading of () operator: perform the Vranic particle merging
//! \param particles   particle object containing the particle
//!                    properties
//! \param smpi        MPI properties
//! \param istart      Index of the first particle
//! \param iend        Index of the last particle
//! \param ithread     Thread index
// ---------------------------------------------------------------------
void MergingVranic::operator() (
        Particles &particles,
        SmileiMPI* smpi,
        int istart,
        int iend,
        int ithread,
        int ipart_ref)
{

    // First of all, we check that there is enought particles per cell
    // to process the merging.
    if ((unsigned int)(iend - istart) > merging_ppc_min_threshold_) {

        // Minima
        double mx_min;
        double theta_min;
        double phi_min;

        // Maxima
        double mx_max;
        double theta_max;
        double phi_max;

        // Delta
        double mx_delta;
        double theta_delta;
        double phi_delta;

        // Index in each direction
        unsigned int mx_i;
        unsigned int theta_i;
        unsigned int phi_i;

        // Local particle index
        unsigned int ic;
        unsigned int ipack;
        int ip;
        unsigned ipart;

        // Dimensions
        dimensions_[0] = 5;
        dimensions_[1] = 5;
        dimensions_[2] = 5;
        unsigned int momentum_cells = dimensions_[0]
                                    * dimensions_[1]
                                    * dimensions_[2];

        // Total weight for merging process
        double total_weight;
        double total_momentum[3];
        double total_energy;

        // Momentum shortcut
        double* momentum[3];
        for ( int i = 0 ; i<3 ; i++ )
            momentum[i] =  &( particles.momentum(i,0) );

        // Weight shortcut
        double *weight = &( particles.weight( 0 ) );

        // Id shortcut
        uint64_t *id = &( particles.id( 0 ) );

        // Norm of the momentum
        double momentum_norm;

        // Local vector to store the momentum index in the momentum discretization
        std::vector <unsigned int> momentum_i(iend-istart,0);

        // Sorted array of particle index
        std::vector <unsigned int> sorted_particles(iend-istart,0);

        // Array containing the number of particles per momentum cells
        std::vector <unsigned int> particles_per_momentum_cells(momentum_cells,0);

        // Array containing the first particle index of each momentum cell
        // in the sorted particle array
        std::vector <unsigned int> momentum_cells_index(momentum_cells,0);

        // Local vector to store the momentum angles in the spherical base
        std::vector <double> phi(iend-istart,0);
        std::vector <double> theta(iend-istart,0);

        // Cell center coordinates in the spherical base
        std::vector <double> cell_center_mx(dimensions_[0],0);
        std::vector <double> cell_center_theta(dimensions_[1],0);
        std::vector <double> cell_center_phi(dimensions_[2],0);

        //std::cerr << iend-istart << std::endl;

        // ________________________________________________
        // First step: Computation of the maxima and minima

        momentum_norm = sqrt(momentum[0][istart]*momentum[0][istart]
                      + momentum[1][istart]*momentum[1][istart]
                      + momentum[2][istart]*momentum[2][istart]);

        mx_min = momentum[0][istart];
        mx_max = mx_min;

        theta[0] = atan2(momentum[1][istart],momentum[0][istart]);
        phi[0]   = asin(momentum[2][istart] / momentum_norm);

        theta_min = theta[0];
        phi_min   = phi[0];

        theta_max = theta_min;
        phi_max   = phi_min;

        for (ipart=istart+1 ; ipart<iend; ipart++ ) {

            // Local array index
            ip = ipart - istart;

            momentum_norm = sqrt(momentum[0][ipart]*momentum[0][ipart]
                          + momentum[1][ipart]*momentum[1][ipart]
                          + momentum[2][ipart]*momentum[2][ipart]);

            mx_min = fmin(mx_min,momentum[0][ipart]);
            mx_max = fmax(mx_max,momentum[0][ipart]);

            phi[ip]   = asin(momentum[2][ipart] / momentum_norm);
            theta[ip] = atan2(momentum[1][ipart] , momentum[0][ipart]);

            theta_min = fmin(theta_min,theta[ip]);
            theta_max = fmax(theta_max,theta[ip]);

            phi_min = fmin(phi_min,phi[ip]);
            phi_max = fmax(phi_max,phi[ip]);

        }

        // Extra to include the max in the discretization
        mx_max += (mx_max - mx_min)*0.01;
        theta_max += (theta_max - theta_min)*0.01;
        phi_max += (phi_max - phi_min)*0.01;

        // __________________________________________________________
        // Second step : Computation of the discretization and steps

        // Computation of the deltas (discretization steps)
        mx_delta = (mx_max - mx_min) / dimensions_[0];
        theta_delta = (theta_max - theta_min) / dimensions_[1];
        phi_delta = (phi_max - phi_min) / dimensions_[2];

        // Check if min and max are very close
        if (mx_delta < 1e-10) {
            mx_delta = 0.;
            dimensions_[0] = 1;
        }
        if (theta_delta < 1e-10) {
            theta_delta = 0.;
            dimensions_[1] = 1;
        }
        if (phi_delta < 1e-10) {
            phi_delta = 0.;
            dimensions_[2] = 1;
        }

        // Computation of the cell centers
        for (mx_i = 0 ; mx_i < dimensions_[0] ; mx_i ++) {
            cell_center_mx[mx_i] = mx_min + (mx_i + 0.5) * mx_delta;
        }
        for (theta_i = 0 ; theta_i < dimensions_[1] ; theta_i ++) {
            cell_center_theta[theta_i] = theta_min + (theta_i + 0.5) * theta_delta;
        }
        for (phi_i = 0 ; phi_i < dimensions_[2] ; phi_i ++) {
            cell_center_phi[phi_i] = phi_min + (phi_i + 0.5) * phi_delta;
        }

        // ___________________________________________________________________
        // Third step: for each particle, momentum indexes are computed in the
        // requested discretization. The number of particles per momentum bin
        // is also determined.

        #pragma omp simd
        for (int ipart=istart ; ipart<iend; ipart++ ) {

            ip = ipart - istart;

            // 3d indexes in the momentum discretization
            mx_i    = (unsigned int)( (momentum[0][ipart] - mx_min)/ mx_delta - 0.5);
            theta_i = (unsigned int)( (theta[ip] - theta_min)       / theta_delta - 0.5);
            phi_i   = (unsigned int)( (phi[ip] - phi_min)           / phi_delta - 0.5);

            // 1d Index in the momentum discretization
            momentum_i[ip] = mx_i    * dimensions_[1]*dimensions_[2]
                          + theta_i * dimensions_[2] + phi_i;

            // Number of particles per momentum cells
            particles_per_momentum_cells[momentum_i[ip]] += 1;

        }

        // Computation of the cell index in the sorted array of particles
        //std::cerr << "momentum_cells_index building" << std::endl;
        for (ic = 1 ; ic < momentum_cells ; ic++) {
            momentum_cells_index[ic]  = momentum_cells_index[ic-1] + particles_per_momentum_cells[ic-1];
            // std::cerr << "ic: " << ic
            //           << " momentum_cells_index[ic]: " << momentum_cells_index[ic]
            //           << " particles_per_momentum_cells[ic]: " << particles_per_momentum_cells[ic]
            //           << " Particles: " << iend-istart
            //           << std::endl;
            particles_per_momentum_cells[ic] = 0;
        }

        // ___________________________________________________________________
        // Fourth step: sort particles in correct bins according to their
        // momentum properties

        for (int ip=0 ; ip<iend-istart; ip++ ) {

            // Momentum cell for this particle
            ic = momentum_i[ip];

            // std::cerr << "Momentum index:" << momentum_i[ip]
            //           << " / " << momentum_cells
            //           << ", mom cell index: " << momentum_cells_index[ic]
            //           << " / " << iend-istart-1
            //           << ", ppmc: " << particles_per_momentum_cells[ic]
            //           << std::endl;

            sorted_particles[momentum_cells_index[ic]
            + particles_per_momentum_cells[ic]] = ip;

            particles_per_momentum_cells[ic] += 1;
        }


        // ___________________________________________________________________
        // Fifth step: for each momentum bin, merge packet of particles composed of
        // at least `min_packet_size_` and `max_packet_size_`

        // Loop over the the momentum cells that have enough particules
        for (ic=0 ; ic<momentum_cells; ic++ ) {
            if (particles_per_momentum_cells[ic] >= 4 )
            {
                /*std::cerr << "ic: " << ic
                          << ", ppmc: " << particles_per_momentum_cells[ic]
                          << std::endl;*/

                // Loop over the packets of particles that can be merged
                for (ipack = 0 ; ipack < particles_per_momentum_cells[ic] ; ipack += 4) {

                    total_weight = 0;
                    total_momentum[0] = 0;
                    total_momentum[1] = 0;
                    total_momentum[2] = 0;
                    total_energy = 0;

                    std::cerr << "Merging start: " << std::endl;

                    // Compute total weight, total momentum and total energy
                    for (ip = ipack ; ip < ipack + 4 ; ip ++) {

                        ipart = sorted_particles[momentum_cells_index[ic] + ip];

                        std::cerr << " ipart: " << ipart << std::endl;
                        std::cerr << " w: "  << weight[ipart]
                                  << " mx: " << momentum[0][ipart] << std::endl;

                        //total_weight += weight[ipart];

                        // total momentum vector (pt)
                        //total_momentum[0] += momentum[0][ipart]*weight[ipart];
                        //total_momentum[1] += momentum[1][ipart]*weight[ipart];
                        //total_momentum[2] += momentum[2][ipart]*weight[ipart];

                        // total energy
                        /*total_energy += sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                                                 + momentum[1][ipart]*momentum[1][ipart]
                                                 + momentum[2][ipart]*momentum[2][ipart]);*/
                    }

                }
            }
        }

    }

}
