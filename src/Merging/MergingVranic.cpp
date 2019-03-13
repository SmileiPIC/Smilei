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
    // Momentum cell discretization
    dimensions_[0] = (unsigned int)(species->merge_momentum_cell_size_[0]);
    dimensions_[1] = (unsigned int)(species->merge_momentum_cell_size_[1]);
    dimensions_[2] = (unsigned int)(species->merge_momentum_cell_size_[2]);

    // Min and max particle per cell number
    min_packet_size_ = species->merge_min_packet_size_;
    max_packet_size_ = species->merge_max_packet_size_;
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
//! \param count       Final number of particles
// ---------------------------------------------------------------------
void MergingVranic::operator() (
        Particles &particles,
        SmileiMPI* smpi,
        int istart,
        int iend,
        int & count)
        //unsigned int &remaining_particles,
        //unsigned int &merged_particles)
{

    unsigned int number_of_particles = (unsigned int)(iend - istart);

    // First of all, we check that there is enought particles per cell
    // to process the merging.
    if (number_of_particles > ppc_min_threshold_) {

        // Minima
        double mr_min;
        double theta_min;
        double phi_min;

        // Maxima
        double mr_max;
        double theta_max;
        double phi_max;

        // Delta
        double mr_delta;
        double theta_delta;
        double phi_delta;

        // Angles
        double phi;
        double theta;
        double omega;
        double cos_omega;
        double sin_omega;

        // Index in each direction
        unsigned int mr_i;
        unsigned int theta_i;
        unsigned int phi_i;

        // Local particle index
        unsigned int ic, icc;
        unsigned int ipack;
        unsigned int npack;
        unsigned int ip, ip_min, ip_max;
        unsigned int ipart;

        // Total number of momentum cells
        unsigned int momentum_cells = dimensions_[0]
                                    * dimensions_[1]
                                    * dimensions_[2];

        // Total number of angular momentum cells
        unsigned int momentum_angular_cells = dimensions_[1]
                                            * dimensions_[2];

        // Total weight for merging process
        double total_weight;
        double total_momentum[3];
        double total_momentum_norm;
        double total_energy;

        // New particle properties
        double new_energy;
        double new_momentum_norm;
        double e1_x,e1_y,e1_z;
        double e2_x,e2_y,e2_z;
        double e2_norm;

        // Momentum shortcut
        double* momentum[3];
        for ( int i = 0 ; i<3 ; i++ )
            momentum[i] =  &( particles.momentum(i,0) );

        // Weight shortcut
        double *weight = &( particles.weight( 0 ) );

        // CEll keys shortcut
        int *cell_keys = &( particles.cell_keys[0] );

        // Norm of the momentum
        double momentum_norm;

        // Local vector to store the momentum index in the momentum discretization
        std::vector <unsigned int> momentum_cell_index(iend-istart,0);

        // Sorted array of particle index
        std::vector <unsigned int> sorted_particles(iend-istart,0);

        // Array containing the number of particles per momentum cells
        std::vector <unsigned int> particles_per_momentum_cells(momentum_cells,0);

        // Array containing the first particle index of each momentum cell
        // in the sorted particle array
        std::vector <unsigned int> momentum_cell_particle_index(momentum_cells,0);

        // Local vector to store the momentum angles in the spherical base
        std::vector <double> particles_phi(iend-istart,0);
        std::vector <double> particles_theta(iend-istart,0);

        // Cell direction unit vector in the spherical base
        std::vector <double> cell_vec_x(dimensions_[1]*dimensions_[2],0);
        std::vector <double> cell_vec_y(dimensions_[1]*dimensions_[2],0);
        std::vector <double> cell_vec_z(dimensions_[1]*dimensions_[2],0);

        //std::cerr << iend-istart << std::endl;

        // ________________________________________________
        // Computation of the maxima and minima for each direction

        momentum_norm = sqrt(momentum[0][istart]*momentum[0][istart]
                      + momentum[1][istart]*momentum[1][istart]
                      + momentum[2][istart]*momentum[2][istart]);

        mr_min = momentum_norm;
        mr_max = mr_min;

        particles_theta[0] = atan2(momentum[1][istart],momentum[0][istart]);
        particles_phi[0]   = asin(momentum[2][istart] / momentum_norm);

        theta_min = particles_theta[0];
        phi_min   = particles_phi[0];

        theta_max = theta_min;
        phi_max   = phi_min;

        for (ipart=istart+1 ; ipart<iend; ipart++ ) {

            // Local array index
            ip = ipart - istart;

            momentum_norm = sqrt(momentum[0][ipart]*momentum[0][ipart]
                          + momentum[1][ipart]*momentum[1][ipart]
                          + momentum[2][ipart]*momentum[2][ipart]);

            mr_min = fmin(mr_min,momentum_norm);
            mr_max = fmax(mr_max,momentum_norm);

            particles_phi[ip]   = asin(momentum[2][ipart] / momentum_norm);
            particles_theta[ip] = atan2(momentum[1][ipart] , momentum[0][ipart]);

            theta_min = fmin(theta_min,particles_theta[ip]);
            theta_max = fmax(theta_max,particles_theta[ip]);

            phi_min = fmin(phi_min,particles_phi[ip]);
            phi_max = fmax(phi_max,particles_phi[ip]);

        }

        // Extra to include the max in the discretization
        mr_max += (mr_max - mr_min)*0.01;
        theta_max += (theta_max - theta_min)*0.01;
        phi_max += (phi_max - phi_min)*0.01;

        // __________________________________________________________
        // Second step : Computation of the discretization and steps

        // Computation of the deltas (discretization steps)
        mr_delta = (mr_max - mr_min) / dimensions_[0];
        theta_delta = (theta_max - theta_min) / dimensions_[1];
        phi_delta = (phi_max - phi_min) / dimensions_[2];

        // Check if min and max are very close
        if (mr_delta < 1e-10) {
            mr_delta = 0.;
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

        // Computation of the cell direction unit vector (vector d in Vranic et al.)
        #pragma omp simd collapse(2)
        for (theta_i = 0 ; theta_i < dimensions_[1] ; theta_i ++) {
            for (phi_i = 0 ; phi_i < dimensions_[2] ; phi_i ++) {

                icc = theta_i * dimensions_[2] + phi_i;

                theta = theta_min + (theta_i + 0.5) * theta_delta;
                phi = phi_min + (phi_i + 0.5) * phi_delta;

                cell_vec_x[icc] = cos(phi)*cos(theta);
                cell_vec_y[icc] = cos(phi)*sin(theta);
                cell_vec_z[icc] = sin(phi);
            }
        }

        /*std::cerr << " Position cell: "
                  << " particles: " << iend - istart
                  << " mr_min: " << mr_min << " mr_max: " << mr_max << " dmx: " << mr_delta
                  << " theta_min: " << theta_min << " theta_max: " << theta_max
                  << " phi_min: " << phi_min << " phi_max: " << phi_max
                  << std::endl;*/

        // For each particle, momentum cell indexes are computed in the
        // requested discretization.
        // This loop can be efficiently vectorized
        #pragma omp simd
        for (ipart= (unsigned int)(istart) ; ipart<(unsigned int)(iend); ipart++ ) {

            ip = ipart - istart;

            momentum_norm = sqrt(momentum[0][ipart]*momentum[0][ipart]
                          + momentum[1][ipart]*momentum[1][ipart]
                          + momentum[2][ipart]*momentum[2][ipart]);

            // 3d indexes in the momentum discretization
            mr_i    = (unsigned int)( (momentum_norm - mr_min) / mr_delta);
            theta_i = (unsigned int)( (particles_theta[ip] - theta_min)  / theta_delta);
            phi_i   = (unsigned int)( (particles_phi[ip] - phi_min)      / phi_delta);

            // 1D Index in the momentum discretization
            momentum_cell_index[ip] = phi_i    * dimensions_[0]*dimensions_[1]
                          + theta_i * dimensions_[0] + mr_i;
        }

        // The number of particles per momentum cells
        // is then determined.
        // No vectorization because of random memory accesses
        for (ip=0; ip<number_of_particles; ip++ ) {
            // Number of particles per momentum cells
            particles_per_momentum_cells[momentum_cell_index[ip]] += 1;

        }

        // Computation of the cell index in the sorted array of particles
        //std::cerr << "momentum_cell_particle_index building" << std::endl;
        for (ic = 1 ; ic < momentum_cells ; ic++) {
            momentum_cell_particle_index[ic]  = momentum_cell_particle_index[ic-1] + particles_per_momentum_cells[ic-1];
            // std::cerr << "ic: " << ic
            //           << " momentum_cell_particle_index[ic]: " << momentum_cell_particle_index[ic]
            //           << " particles_per_momentum_cells[ic]: " << particles_per_momentum_cells[ic]
            //           << " Particles: " << iend-istart
            //           << std::endl;
            particles_per_momentum_cells[ic-1] = 0;
        }
        particles_per_momentum_cells[momentum_cells-1] = 0;

        // sort particles in correct bins according to their
        // momentum properties
        for (ip=0 ; ip<number_of_particles; ip++ ) {

            // Momentum cell index for this particle
            ic = momentum_cell_index[ip];

            // std::cerr << "Momentum index:" << momentum_cell_index[ip]
            //           << " / " << momentum_cells
            //           << ", mom cell index: " << momentum_cell_particle_index[ic]
            //           << " / " << iend-istart-1
            //           << ", ppmc: " << particles_per_momentum_cells[ic]
            //           << std::endl;

            sorted_particles[momentum_cell_particle_index[ic]
            + particles_per_momentum_cells[ic]] = istart + ip;

            particles_per_momentum_cells[ic] += 1;
        }

        // Debugging
        /*for (mr_i=0 ; mr_i< dimensions_[0]; mr_i++ ) {
            for (theta_i=0 ; theta_i< dimensions_[1]; theta_i++ ) {
                for (phi_i=0 ; phi_i< dimensions_[2]; phi_i++ ) {

                    ic = mr_i * dimensions_[1]*dimensions_[2]
                       + theta_i * dimensions_[2] + phi_i;

                    std::cerr << " Momentum cell: " << ic
                              << "   Momentum cell p index: " << momentum_cell_particle_index[ic]
                              << "   Particles: " << particles_per_momentum_cells[ic]
                              << "  mx in [" << mr_i * mr_delta + mr_min  << ", "
                                            << (mr_i+1) * mr_delta + mr_min  << "]"
                              << "  theta in [" << theta_i * theta_delta + theta_min  << ", "
                                            << (theta_i+1) * theta_delta + theta_min  << "]"
                              << std::endl;
                    for (ip = 0 ; ip < particles_per_momentum_cells[ic] ; ip ++) {
                        ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];

                        momentum_norm = sqrt(momentum[0][ipart]*momentum[0][ipart]
                                      + momentum[1][ipart]*momentum[1][ipart]
                                      + momentum[2][ipart]*momentum[2][ipart]);

                        std::cerr << "   - Id: " << ipart - istart
                                  << "     id2: " << ipart
                                  << "     mx: " << momentum_norm
                                  << std::endl;
                    }
                }
            }
        }*/


        // For each momentum bin, merge packet of particles composed of
        // at least `min_packet_size_` and `max_packet_size_`

        // Loop over the the momentum cells that have enough particules
        for (phi_i=0 ; phi_i< dimensions_[2]; phi_i++ ) {
            for (theta_i=0 ; theta_i< dimensions_[1]; theta_i++ ) {

                // 1D cell direction index
                icc = theta_i + phi_i* dimensions_[1] ;

                for (mr_i=0 ; mr_i< dimensions_[0]; mr_i++ ) {

                    // 1D cell index
                    ic = mr_i + icc*dimensions_[0];

                    if (particles_per_momentum_cells[ic] >= min_packet_size_ ) {
                        /*std::cerr << "ic: " << ic
                                  << ", ppmc: " << particles_per_momentum_cells[ic]
                                  << std::endl;*/

                        /*std::cerr << " - Momentum cell: " << ic
                                << "  mx in [" << mr_i * mr_delta + mr_min  << ", "
                                              << (mr_i+1) * mr_delta + mr_min  << "]"
                                << "  theta in [" << theta_i * theta_delta + theta_min  << ", "
                                              << (theta_i+1) * theta_delta + theta_min  << "]"
                                << std::endl;*/

                        // Computation of the number of particle packets to merge
                        npack = particles_per_momentum_cells[ic]/max_packet_size_;

                        // Check if the rest is sufficient to add an additional smaller packet
                        if (particles_per_momentum_cells[ic]%max_packet_size_ >= min_packet_size_) {
                            npack += 1;
                        }

                        // Loop over the packets of particles that can be merged
                        for (ipack = 0 ; ipack < npack ; ipack += 1) {

                            total_weight = 0;
                            total_momentum[0] = 0;
                            total_momentum[1] = 0;
                            total_momentum[2] = 0;
                            total_energy = 0;

                            /*std::cerr << "   Merging start: " << std::endl;*/
                            ip_min = ipack*max_packet_size_;
                            ip_max = std::min((ipack+1)*max_packet_size_,number_of_particles);

                            // Compute total weight, total momentum and total energy
                            for (ip = ip_min ; ip < ip_max ; ip ++) {

                                ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];

                                /*std::cerr << "   ipart: " << ipart << std::endl;
                                std::cerr << "   w: "  << weight[ipart]
                                          << "   mx: " << fabs(momentum[0][ipart]) << std::endl;*/

                                // Total weight (wt)
                                total_weight += weight[ipart];

                                // total momentum vector (pt)
                                total_momentum[0] += momentum[0][ipart]*weight[ipart];
                                total_momentum[1] += momentum[1][ipart]*weight[ipart];
                                total_momentum[2] += momentum[2][ipart]*weight[ipart];

                                // total energy (\varespilon_t)
                                total_energy += weight[ipart]
                                                         * sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                                                         + momentum[1][ipart]*momentum[1][ipart]
                                                         + momentum[2][ipart]*momentum[2][ipart]);
                            }

                            // \varepsilon_a in Vranic et al
                            new_energy = total_energy / total_weight;

                            // pa in Vranic et al.
                            new_momentum_norm = sqrt(new_energy*new_energy - 1.0);

                            total_momentum_norm = sqrt(total_momentum[0]*total_momentum[0]
                                                +      total_momentum[1]*total_momentum[1]
                                                +      total_momentum[2]*total_momentum[2]);

                            // Angle between pa and pt, pb and pt in Vranic et al.
                            omega = acos(total_momentum_norm / (total_weight*new_momentum_norm));
                            sin_omega = sin(omega);
                            cos_omega = cos(omega);

                            // Computation of e1 unit vector
                            e1_x = total_momentum[0] / total_momentum_norm;
                            e1_y = total_momentum[1] / total_momentum_norm;
                            e1_z = total_momentum[2] / total_momentum_norm;

                            // Computation of e2  = e1 x e3 unit vector
                            // e3 = e1 x cell_vec
                            e2_x = e1_y*e1_y*cell_vec_x[icc]
                                 - e1_x * (e1_y*cell_vec_y[icc] + e1_z*cell_vec_z[icc])
                                 + e1_z*e1_z*cell_vec_x[icc];
                            e2_y = e1_z*e1_z*cell_vec_y[icc]
                                 - e1_y * (e1_z*cell_vec_z[icc] + e1_x*cell_vec_x[icc])
                                 + e1_x*e1_x*cell_vec_y[icc];
                            e2_z = e1_x*e1_x*cell_vec_z[icc]
                                 - e1_z * (e1_x*cell_vec_x[icc] + e1_y*cell_vec_y[icc])
                                 + e1_y*e1_y*cell_vec_z[icc];

                            e2_norm = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z);

                            // e2 is normalized to be a unit vector
                            e2_x = e2_x / e2_norm;
                            e2_y = e2_y / e2_norm;
                            e2_z = e2_z / e2_norm;

                            // The first 2 particles of the list will
                            // be the merged particles.

                            // Update momentum of the first particle
                            ipart = sorted_particles[momentum_cell_particle_index[ic] + ipack*max_packet_size_];
                            momentum[0][ipart] = new_momentum_norm*(cos_omega*e1_x + sin_omega*e2_x);
                            momentum[1][ipart] = new_momentum_norm*(cos_omega*e1_y + sin_omega*e2_y);
                            momentum[2][ipart] = new_momentum_norm*(cos_omega*e1_z + sin_omega*e2_z);
                            weight[ipart] = 0.5*total_weight;

                            // Update momentum of the second particle
                            ipart = sorted_particles[momentum_cell_particle_index[ic] + ipack*max_packet_size_ + 1];
                            momentum[0][ipart] = new_momentum_norm*(cos_omega*e1_x - sin_omega*e2_x);
                            momentum[1][ipart] = new_momentum_norm*(cos_omega*e1_y - sin_omega*e2_y);
                            momentum[2][ipart] = new_momentum_norm*(cos_omega*e1_z - sin_omega*e2_z);
                            weight[ipart] = 0.5*total_weight;

                            // std::cerr << " Total energy: " << sqrt(1.0
                            //                                   + momentum[0][ipart]*momentum[0][ipart]
                            //                                   + momentum[1][ipart]*momentum[1][ipart]
                            //                                   + momentum[2][ipart]*momentum[2][ipart])
                            //           << " New energy: " << new_energy
                            //           << " Momentum norm: " << sqrt(momentum[0][ipart]*momentum[0][ipart]
                            //           + momentum[1][ipart]*momentum[1][ipart]
                            //           + momentum[2][ipart]*momentum[2][ipart])
                            //           << " new momentum: " << new_momentum_norm
                            //           << std::endl;

                            // Other particles are tagged to be removed after
                            for (ip = ip_min + 2; ip < ip_max ; ip ++) {
                                ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];
                                cell_keys[ipart] = -1;
                                count--;
                            }

                        }
                    }

                }
            }
        }

    }

}
