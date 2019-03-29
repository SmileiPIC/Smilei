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
        double mass,
        Particles &particles,
        std::vector <int> &mask,
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
    if (number_of_particles > min_particles_per_cell) {

        // Momentum discretization
        unsigned int dim[3];
        for (unsigned int i = 0; i < 3 ; i++) {
            dim[i] = dimensions_[i];
        }

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

        // Inverse Delta
        double inv_mr_delta;
        double inv_theta_delta;
        double inv_phi_delta;

        // Angles
        double phi;
        double theta;
        double mr;
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

        // Total weight for merging process
        double total_weight;
        double total_momentum_x;
        double total_momentum_y;
        double total_momentum_z;
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

        // Cell keys shortcut
        // int *cell_keys = &( particles.cell_keys[0] );

        // Norm of the momentum
        // std::vector <double> momentum_norm(number_of_particles,0);;
        double  * momentum_norm = (double*) aligned_alloc(64, number_of_particles*sizeof(double));

        // Local vector to store the momentum index in the momentum discretization
        // std::vector <unsigned int> momentum_cell_index(number_of_particles,0);
        unsigned int  * momentum_cell_index = (unsigned int*) aligned_alloc(64, number_of_particles*sizeof(unsigned int));

        // Sorted array of particle index
        // std::vector <unsigned int> sorted_particles(number_of_particles,0);
        unsigned int  * sorted_particles = (unsigned int*) aligned_alloc(64, number_of_particles*sizeof(unsigned int));

        // Local vector to store the momentum angles in the spherical base
        // std::vector <double> particles_phi(number_of_particles,0);
        // std::vector <double> particles_theta(number_of_particles,0);
        double  * particles_phi = (double*) aligned_alloc(64, number_of_particles*sizeof(double));
        double  * particles_theta = (double*) aligned_alloc(64, number_of_particles*sizeof(double));

        // Arrays for inefficient vectorization strategy
        // double  * momentum_x_loc = (double*) aligned_alloc(64, max_packet_size_*sizeof(double));
        // double  * momentum_y_loc = (double*) aligned_alloc(64, max_packet_size_*sizeof(double));
        // double  * momentum_z_loc = (double*) aligned_alloc(64, max_packet_size_*sizeof(double));
        // double  * weight_loc = (double*)  aligned_alloc(64, max_packet_size_*sizeof(double));
        // unsigned int ipp;

        // Computation of the particle momentum properties
        #pragma omp simd aligned(momentum_norm, particles_theta, particles_phi: 64) private(ip)
        for (ipart=(unsigned int)(istart) ; ipart<(unsigned int) (iend); ipart++ ) {

            // Local array index
            ip = ipart - istart;

            momentum_norm[ip] = sqrt(momentum[0][ipart]*momentum[0][ipart]
                          + momentum[1][ipart]*momentum[1][ipart]
                          + momentum[2][ipart]*momentum[2][ipart]);

            particles_phi[ip]   = asin(momentum[2][ipart] / momentum_norm[ip]);
            particles_theta[ip] = atan2(momentum[1][ipart] , momentum[0][ipart]);
        }

        // Computation of the maxima and minima for each direction
        mr_min = momentum_norm[0];
        mr_max = mr_min;

        theta_min = particles_theta[0];
        phi_min   = particles_phi[0];

        theta_max = theta_min;
        phi_max   = phi_min;

        #pragma omp simd \
        reduction(min:mr_min) reduction(min:theta_min) reduction(min:phi_min) \
        reduction(max:mr_max) reduction(max:theta_max) reduction(max:phi_max) \
        aligned(momentum_norm, particles_theta, particles_phi: 64)
        for (ip=1 ; ip < number_of_particles; ip++ ) {
            mr_min = fmin(mr_min,momentum_norm[ip]);
            mr_max = fmax(mr_max,momentum_norm[ip]);

            theta_min = fmin(theta_min,particles_theta[ip]);
            theta_max = fmax(theta_max,particles_theta[ip]);

            phi_min = fmin(phi_min,particles_phi[ip]);
            phi_max = fmax(phi_max,particles_phi[ip]);

        }

        // Extra to include the max in the discretization
        // mr_max += (mr_max - mr_min)*0.01;
        // theta_max += (theta_max - theta_min)*0.01;
        // phi_max += (phi_max - phi_min)*0.01;

        // Computation of the deltas (discretization steps)
        // Check if min and max boundaries are very close
        if (abs((mr_max - mr_min)) < 1e-10) {
            mr_delta = 0.1;
            inv_mr_delta = 0;
            dim[0] = 1;
        } else {
            // mr_delta = (mr_max - mr_min) / dim[0];
            dim[0] += 1;
            mr_delta = (mr_max - mr_min) / (dim[0]-1);
            // A bit of chaos to kill the accumulation effect
            mr_min -= 0.99*mr_delta*Rand::uniform();
            inv_mr_delta = 1./mr_delta;
        }
        if (abs((theta_max - theta_min)) < 1e-10) {
            theta_delta = 0.1*M_PI;
            inv_theta_delta = 0;
            dim[1] = 1;
        } else {
            //theta_delta = (theta_max - theta_min) / dim[1];
            dim[1] += 1;
            theta_delta = (theta_max - theta_min) / (dim[1]-1);
            // A bit of chaos to kill the accumulation effect
            theta_min -= 0.99*theta_delta*Rand::uniform();
            inv_theta_delta = 1./theta_delta;
        }
        if (abs((phi_max - phi_min)) < 1e-10) {
            phi_delta = 0.1*M_PI;
            inv_phi_delta = 0;
            dim[2] = 1;
        } else {
            //phi_delta = (phi_max - phi_min) / dim[2];
            dim[2] += 1;
            phi_delta = (phi_max - phi_min) / (dim[2]-1);
            // A bit of chaos to kill the accumulation effect
            phi_min -= 0.99*phi_delta*Rand::uniform();
            inv_phi_delta = 1./phi_delta;
        }

        // std::cerr << " phi min: " << phi_min << " max: " << phi_max
        //           << " phi error: " << abs((phi_max - phi_min))
        //           << " phi_delta: " << phi_delta << std::endl;

        // Total number of momentum cells
        unsigned int momentum_cells = dim[0]
                                    * dim[1]
                                    * dim[2];

        // Total number of angular momentum cells
        unsigned int momentum_angular_cells = dim[1]
                                            * dim[2];

        // Array containing the number of particles per momentum cells
        // std::vector <unsigned int> particles_per_momentum_cells(momentum_cells,0);
        unsigned int  * particles_per_momentum_cells = (unsigned int*) aligned_alloc(64, momentum_cells*sizeof(unsigned int));

        // Array containing the first particle index of each momentum cell
        // in the sorted particle array
        // std::vector <unsigned int> momentum_cell_particle_index(momentum_cells,0);
        unsigned int  * momentum_cell_particle_index = (unsigned int*) aligned_alloc(64, momentum_cells*sizeof(unsigned int));

        // Initialization when using aligned_alloc
        for (ic = 0 ; ic < momentum_cells ; ic++) {
            momentum_cell_particle_index[ic] = 0;
            particles_per_momentum_cells[ic] = 0;
        }

        // Computation of the cell direction unit vector (vector d in Vranic et al.)
        // Only necessary for mass particles

        double  * cell_vec_x ;
        double  * cell_vec_y ;
        double  * cell_vec_z ;

        if (mass > 0) {

            // Cell direction unit vector in the spherical base
            // std::vector <double> cell_vec_x(momentum_angular_cells,0);
            // std::vector <double> cell_vec_y(momentum_angular_cells,0);
            // std::vector <double> cell_vec_z(momentum_angular_cells,0);
            cell_vec_x = (double*) aligned_alloc(64, momentum_angular_cells*sizeof(double));
            cell_vec_y = (double*) aligned_alloc(64, momentum_angular_cells*sizeof(double));
            cell_vec_z = (double*) aligned_alloc(64, momentum_angular_cells*sizeof(double));

            for (phi_i = 0 ; phi_i < dim[2] ; phi_i ++) {
                #pragma omp simd private(theta, phi, icc) \
                aligned(cell_vec_x, cell_vec_y, cell_vec_z: 64)
                for (theta_i = 0 ; theta_i < dim[1] ; theta_i ++) {

                    icc = theta_i + phi_i * dim[1];

                    theta = theta_min + (theta_i + 0.5) * theta_delta;
                    phi = phi_min + (phi_i + 0.5) * phi_delta;

                    cell_vec_x[icc] = cos(phi)*cos(theta);
                    cell_vec_y[icc] = cos(phi)*sin(theta);
                    cell_vec_z[icc] = sin(phi);

                    // std::cerr << " theta: " << theta
                    //           << " phi: " << phi
                    //           << std::endl;

                }
            }
        }

        // For each particle, momentum cell indexes are computed in the
        // requested discretization.
        // This loop can be efficiently vectorized
        #pragma omp simd \
        aligned(momentum_norm, particles_theta, particles_phi: 64) \
        aligned(momentum_cell_index: 64)
        for (ip= 0; ip < number_of_particles ; ip++ ) {

            // 3d indexes in the momentum discretization
            mr_i    = (unsigned int) floor( (momentum_norm[ip] - mr_min) * inv_mr_delta);
            theta_i = (unsigned int) floor( (particles_theta[ip] - theta_min)  * inv_theta_delta);
            phi_i   = (unsigned int) floor( (particles_phi[ip] - phi_min)      * inv_phi_delta);

            // 1D Index in the momentum discretization
            momentum_cell_index[ip] = phi_i * dim[0]*dim[1]
                          + theta_i * dim[0] + mr_i;

                // std::cerr << "momentum_cell_index: " << momentum_cell_index[ip]
                //           << " momentum_norm: " << momentum_norm[ip]
                //           << " particles_theta: " << particles_theta[ip]
                //           << " particles_phi: " << particles_phi[ip]
                //           << " phi_i: " << phi_i
                //           << " theta_i: " << theta_i
                //           << " mr_i: " << mr_i
                //           << std::endl;

        }

        // The number of particles per momentum cells
        // is then determined.
        // No vectorization because of random memory accesses
        for (ip=0; ip<number_of_particles; ip++ ) {
            // Number of particles per momentum cells
            particles_per_momentum_cells[momentum_cell_index[ip]] += 1;
        }

        // Computation of the cell index in the sorted array of particles
        for (ic = 1 ; ic < momentum_cells ; ic++) {
            momentum_cell_particle_index[ic]  = momentum_cell_particle_index[ic-1] + particles_per_momentum_cells[ic-1];
            particles_per_momentum_cells[ic-1] = 0;
        }
        particles_per_momentum_cells[momentum_cells-1] = 0;

        // sort particles in correct bins according to their
        // momentum properties
        // Not vectorizable
        for (ip=0 ; ip<number_of_particles; ip++ ) {

            // Momentum cell index for this particle
            ic = momentum_cell_index[ip];

            sorted_particles[momentum_cell_particle_index[ic]
            + particles_per_momentum_cells[ic]] = istart + ip;

            particles_per_momentum_cells[ic] += 1;
        }

        // Debugging
        /*for (mr_i=0 ; mr_i< dim[0]; mr_i++ ) {
            for (theta_i=0 ; theta_i< dim[1]; theta_i++ ) {
                for (phi_i=0 ; phi_i< dim[2]; phi_i++ ) {

                    ic = mr_i * dim[1]*dim[2]
                       + theta_i * dim[2] + phi_i;

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

                        std::cerr << "   - Id: " << ipart - istart
                                  << "     id2: " << ipart
                                  << "     mx: " << momentum_norm[ip]
                                  << std::endl;
                    }
                }
            }
        }*/


        // For each momentum bin, merge packet of particles composed of
        // at least `min_packet_size_` and `max_packet_size_`

        // Loop over the the momentum cells that have enough particules
        for (phi_i=0 ; phi_i< dim[2]; phi_i++ ) {
            for (theta_i=0 ; theta_i< dim[1]; theta_i++ ) {

                // 1D cell direction index
                icc = theta_i + phi_i* dim[1] ;

                for (mr_i=0 ; mr_i< dim[0]; mr_i++ ) {

                    // 1D cell index
                    ic = mr_i + icc*dim[0];

                    // Check if there is enought particles in the momentum
                    // cell to trigger the merging procecc
                    if (particles_per_momentum_cells[ic] >= min_packet_size_ ) {

                        // Computation of the number of particle packets to merge
                        npack = particles_per_momentum_cells[ic]/max_packet_size_;

                        // Check if the rest is sufficient to add an additional smaller packet
                        if (particles_per_momentum_cells[ic]%max_packet_size_ >= min_packet_size_) {
                            npack += 1;
                        }

                        // Loop over the packets of particles that can be merged
                        for (ipack = 0 ; ipack < npack ; ipack += 1) {

                            total_weight = 0;
                            // total_momentum_x = 0;
                            // total_momentum_y = 0;
                            // total_momentum_z = 0;
                            total_energy = 0;
                            mr = 0;
                            theta=0;
                            phi = 0;

                            // First index of the packet
                            ip_min = ipack*max_packet_size_;
                            // last index of the packet
                            ip_max = std::min((ipack+1)*max_packet_size_,particles_per_momentum_cells[ic]);

                            // _______________________________________________________________
                            // This simd optimization is not efficient
                            // kept for memory

                            // Store in an aligned vector the momentum and weight of
                            // the particles of the packet
                            // Vectorizable but important Random Memory Accesses
                            // for (ip = ip_min ; ip < ip_max ; ip ++) {
                            //    ipp = ip - ip_min;
                            //    ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];
                            //    momentum_x_loc[ipp] = momentum[0][ipart];
                            //    momentum_y_loc[ipp] = momentum[1][ipart];
                            //    momentum_z_loc[ipp] = momentum[2][ipart];
                            //    weight_loc[ipp] = weight[ipart];
                            //}

                            // Compute total weight, total momentum and total energy
                            // #pragma omp simd reduction(+:total_weight) \
                            // reduction(+:total_momentum_x) reduction(+:total_momentum_y) \
                            // reduction(+:total_momentum_z) reduction(+:total_energy) \
                            // aligned(weight_loc, momentum_x_loc, momentum_y_loc, momentum_z_loc: 64)
                            // for (ipp = 0 ; ipp < ip_max - ip_min ; ipp ++) {

                            //     // Total weight (wt)
                            //     total_weight += weight_loc[ipp];

                            //     // total momentum vector (pt)
                            //     total_momentum_x += momentum_x_loc[ipp]*weight_loc[ipp];
                            //     total_momentum_y += momentum_y_loc[ipp]*weight_loc[ipp];
                            //     total_momentum_z += momentum_z_loc[ipp]*weight_loc[ipp];

                                // total energy (\varespilon_t)
                            //     total_energy += weight_loc[ipp]
                            //                              * sqrt(1.0 + momentum_x_loc[ipp]*momentum_x_loc[ipp]
                            //                              + momentum_y_loc[ipp]*momentum_y_loc[ipp]
                            //                              + momentum_z_loc[ipp]*momentum_z_loc[ipp]);
                            // }

                            // _______________________________________________________________

                            // Compute total weight, total momentum and total energy
                            // Photons
                            if (mass == 0) {

                                for (ip = ip_min ; ip < ip_max ; ip ++) {

                                    // Particle index in Particles
                                    ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];

                                    // Total weight (wt)
                                    total_weight += weight[ipart];

                                    mr += weight[ipart]*momentum_norm[ipart - istart];
                                    theta += weight[ipart]*particles_theta[ipart - istart];
                                    phi += weight[ipart]*particles_phi[ipart - istart];

                                    // total momentum vector (pt)
                                    // total_momentum_x += momentum[0][ipart]*weight[ipart];
                                    // total_momentum_y += momentum[1][ipart]*weight[ipart];
                                    // total_momentum_z += momentum[2][ipart]*weight[ipart];

                                    // total energy (\varespilon_t)
                                    // total_energy += weight[ipart]*momentum_norm[ipart - istart];
                                    // std::cerr << " E: " << sqrt(momentum[0][ipart]*momentum[0][ipart]
                                    //                  + momentum[1][ipart]*momentum[1][ipart]
                                    //                  + momentum[2][ipart]*momentum[2][ipart])
                                    //         << " w: " << weight[ipart]
                                    //         << " mr: " << momentum_norm[ipart - istart]
                                    //         << " phi: " << particles_phi[ipart - istart]
                                    //         << " theta: " << particles_theta[ipart - istart]
                                    //         << std::endl;
                                }

                                // here for photon, this parameter is the inverse of total_weight
                                new_momentum_norm = 1 / total_weight;

                                mr = mr * new_momentum_norm;
                                theta = theta * new_momentum_norm;
                                phi = phi * new_momentum_norm;

                                total_momentum_x = mr * cos(theta)*cos(phi);
                                total_momentum_y = mr * sin(theta)*cos(phi);
                                total_momentum_z = mr * sin(phi);

                                // \varepsilon_a in Vranic et al
                                // new_energy = total_energy / total_weight;

                                // total_momentum_norm = 1./sqrt(total_momentum_x*total_momentum_x
                                //                     +      total_momentum_y*total_momentum_y
                                //                     +      total_momentum_z*total_momentum_z);

                                // Update momentum of the first photon
                                ipart = sorted_particles[momentum_cell_particle_index[ic] + ip_min];
                                // momentum[0][ipart] = new_momentum_norm*total_momentum_x*total_momentum_norm;
                                // momentum[1][ipart] = new_momentum_norm*total_momentum_y*total_momentum_norm;
                                // momentum[2][ipart] = new_momentum_norm*total_momentum_z*total_momentum_norm;
                                momentum[0][ipart] = total_momentum_x;
                                momentum[1][ipart] = total_momentum_y;
                                momentum[2][ipart] = total_momentum_z;
                                weight[ipart] = total_weight;

                                // std::cerr << " Energy final: " << sqrt(momentum[0][ipart]*momentum[0][ipart]
                                //                  + momentum[1][ipart]*momentum[1][ipart]
                                //                  + momentum[2][ipart]*momentum[2][ipart])
                                //                  << " phi: " << asin(momentum[2][ipart] / momentum_norm[ip])
                                //                  << " theta: " << atan2(momentum[1][ipart] , momentum[0][ipart])
                                //                  << std::endl;

                                // Other photons are tagged to be removed after
                                for (ip = ip_min + 1; ip < ip_max ; ip ++) {
                                    ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];
                                    mask[ipart] = -1;
                                    count--;
                                }

                            // Mass particles
                            } else {
                                for (ip = ip_min ; ip < ip_max ; ip ++) {

                                    // Particle index in Particles
                                    ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];

                                    // Total weight (wt)
                                    total_weight += weight[ipart];

                                    // total momentum vector (pt)
                                    // total_momentum_x += momentum[0][ipart]*weight[ipart];
                                    // total_momentum_y += momentum[1][ipart]*weight[ipart];
                                    // total_momentum_z += momentum[2][ipart]*weight[ipart];

                                    mr += weight[ipart]*momentum_norm[ipart - istart];
                                    theta += weight[ipart]*particles_theta[ipart - istart];
                                    phi += weight[ipart]*particles_phi[ipart - istart];

                                    // total energy (\varespilon_t)
                                    total_energy += weight[ipart]
                                                             * sqrt(1.0 + momentum[0][ipart]*momentum[0][ipart]
                                                             + momentum[1][ipart]*momentum[1][ipart]
                                                             + momentum[2][ipart]*momentum[2][ipart]);
                                    // total_energy += weight[ipart]
                                    //     * sqrt(1.0 + momentum_norm[ipart - istart]*momentum_norm[ipart - istart]);
                                }

                                // here for photon, this parameter is the inverse of total_weight
                                new_momentum_norm = 1 / total_weight;

                                //mr = mr * new_momentum_norm;
                                theta = theta * new_momentum_norm;
                                phi = phi * new_momentum_norm;

                                total_momentum_x = mr * cos(theta)*cos(phi);
                                total_momentum_y = mr * sin(theta)*cos(phi);
                                total_momentum_z = mr * sin(phi);

                                // \varepsilon_a in Vranic et al
                                new_energy = total_energy * new_momentum_norm;

                                // pa in Vranic et al.
                                new_momentum_norm = sqrt(new_energy*new_energy - 1.0);

                                total_momentum_norm = sqrt(total_momentum_x*total_momentum_x
                                                    +      total_momentum_y*total_momentum_y
                                                    +      total_momentum_z*total_momentum_z);

                                // Angle between pa and pt, pb and pt in Vranic et al.
                                omega = std::acos(total_momentum_norm / (total_weight*new_momentum_norm));
                                sin_omega = sin(omega);
                                cos_omega = cos(omega);

                                // std::cerr << "new_energy: " << new_energy
                                //           << " new_momentum_norm: " << new_momentum_norm
                                //           << " total_momentum_norm: " << total_momentum_norm
                                //           << " omega: " << omega
                                //           << std::endl;

                                // Computation of e1 unit vector
                                e1_x = total_momentum_x / total_momentum_norm;
                                e1_y = total_momentum_y / total_momentum_norm;
                                e1_z = total_momentum_z / total_momentum_norm;

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
                                ipart = sorted_particles[momentum_cell_particle_index[ic] + ip_min];
                                momentum[0][ipart] = new_momentum_norm*(cos_omega*e1_x + sin_omega*e2_x);
                                momentum[1][ipart] = new_momentum_norm*(cos_omega*e1_y + sin_omega*e2_y);
                                momentum[2][ipart] = new_momentum_norm*(cos_omega*e1_z + sin_omega*e2_z);
                                weight[ipart] = 0.5*total_weight;

                                // if (isnan(momentum[0][ipart])) {
                                    // std::cerr << " mx: " << momentum[0][ipart]
                                    //           << " my: " << momentum[1][ipart]
                                    //           << " mz: " << momentum[2][ipart]
                                    //           << " new_momentum_norm: " << new_momentum_norm
                                    //           << " total_weight: " << total_weight
                                    //           << " total_momentum_norm: " << total_momentum_norm
                                    //           << " weight[ipart]: " << weight[ipart]
                                    //           << " omega: " << omega
                                    //           << " " << total_momentum_norm / (total_weight*new_momentum_norm)
                                    //           << " cos_omega: " << cos_omega << " " << sin_omega
                                    //           << " cell_vec: " << cell_vec_x[icc] << " " << cell_vec_y[icc] << " " << cell_vec_z[icc]
                                    //           << " e1: " << e1_x << " " << e1_y << " " << e1_z
                                    //           << " e2: " << e2_x << " " << e2_y << " " << e2_z
                                    //           << std::endl;
                                // }

                                // Update momentum of the second particle
                                ipart = sorted_particles[momentum_cell_particle_index[ic] + ip_min + 1];
                                momentum[0][ipart] = new_momentum_norm*(cos_omega*e1_x - sin_omega*e2_x);
                                momentum[1][ipart] = new_momentum_norm*(cos_omega*e1_y - sin_omega*e2_y);
                                momentum[2][ipart] = new_momentum_norm*(cos_omega*e1_z - sin_omega*e2_z);
                                weight[ipart] = 0.5*total_weight;


                                // Other particles are tagged to be removed after
                                for (ip = ip_min + 2; ip < ip_max ; ip ++) {
                                    ipart = sorted_particles[momentum_cell_particle_index[ic] + ip];
                                    mask[ipart] = -1;
                                    count--;
                                }

                            }

                        }
                    }

                }
            }
        }

        // Free aligned arrays
        free(momentum_norm);
        free(momentum_cell_index);
        if (mass > 0) {
            free(cell_vec_x);
            free(cell_vec_y);
            free(cell_vec_z);
        }
        free(particles_phi);
        free(particles_theta);
        free(sorted_particles);
        free(particles_per_momentum_cells);
        free(momentum_cell_particle_index);
    }

}
