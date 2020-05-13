// ----------------------------------------------------------------------------
//! \file MergingVranicCartesian.cpp
//
//! \brief Functions of the class MergingVranicCartesian
//! Particle merging with the method of Vranic et al. using
//! a Cartesian geometry
//! Vranic CPC 191 65-73 (2015)
//
//! Creation - 01/2019 - Mathieu Lobet
//
// ----------------------------------------------------------------------------

#include "MergingVranicCartesian.h"

#include <cmath>
#include <cstdlib>

// -----------------------------------------------------------------------------
//! Constructor for RadiationNLandauLifshitz
//! Inherited from Radiation
// -----------------------------------------------------------------------------
MergingVranicCartesian::MergingVranicCartesian(Params& params,
                             Species * species, Random * rand)
      : Merging(params, species, rand)
{
    // Momentum cell discretization
    dimensions_[0] = (unsigned int)(species->merge_momentum_cell_size_[0]);
    dimensions_[1] = (unsigned int)(species->merge_momentum_cell_size_[1]);
    dimensions_[2] = (unsigned int)(species->merge_momentum_cell_size_[2]);

    // Min and max particle per cell number
    min_packet_size_ = species->merge_min_packet_size_;
    max_packet_size_ = species->merge_max_packet_size_;

    // Minimum momentum cell length
    min_momentum_cell_length_[0] = species->merge_min_momentum_cell_length_[0];
    min_momentum_cell_length_[1] = species->merge_min_momentum_cell_length_[1];
    min_momentum_cell_length_[2] = species->merge_min_momentum_cell_length_[2];
    
    // Accumulation correction
    accumulation_correction_ = species->merge_accumulation_correction_;
    
    // Discretization scale
    log_scale_ = species->merge_log_scale_;
    
    // Minimum momentum value in log scale
    min_momentum_log_scale_ = species->merge_min_momentum_log_scale_;
    
}

// -----------------------------------------------------------------------------
//! Destructor for MergingVranicCartesian
// -----------------------------------------------------------------------------
MergingVranicCartesian::~MergingVranicCartesian()
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
void MergingVranicCartesian::operator() (
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
    if (number_of_particles > min_particles_per_cell_) {

        // Momentum discretization
        unsigned int dim[3];
        for (unsigned int i = 0; i < 3 ; i++) {
            dim[i] = dimensions_[i];
        }

        // Minima
        double momentum_min[3];

        // Maxima
        double momentum_max[3];

        // Delta
        double momentum_delta[3];

        // Inverse Delta
        double inv_momentum_delta[3];

        // Angles
        // double omega;
        double cos_omega;
        double sin_omega;

        // Index in each direction
        unsigned int mx_i;
        unsigned int my_i;
        unsigned int mz_i;

        // Adjustement of the discretization
        double nb_delta;

        // Local particle index
        unsigned int ic;
        unsigned int icc;
        unsigned int ipack;
        unsigned int npack;
        unsigned int ip;
        unsigned int ipr, ipr_min, ipr_max;

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
        double e3_x,e3_y,e3_z;
        double e2_x,e2_y,e2_z;
        double e2_norm;

        // Momentum cell directions
        double cell_vec_x;
        double cell_vec_y;
        double cell_vec_z;

        // Momentum shortcut
        double* momentum[3];
        for ( int i = 0 ; i<3 ; i++ )
            momentum[i] =  &( particles.momentum(i,0) );

        // Weight shortcut
        double *weight = &( particles.weight( 0 ) );

        // Cell keys shortcut
        // int *cell_keys = &( particles.cell_keys[0] );

        // Local vector to store the momentum index in the momentum discretization
        std::vector <unsigned int> momentum_cell_index(number_of_particles,0);
//         unsigned int  * momentum_cell_index = (unsigned int*) aligned_alloc(64, number_of_particles*sizeof(unsigned int));

        // Sorted array of particle index
        std::vector <unsigned int> sorted_particles(number_of_particles,0);
//         unsigned int  * sorted_particles = (unsigned int*) aligned_alloc(64, number_of_particles*sizeof(unsigned int));

        // Particle gamma factor
        std::vector <double> gamma(number_of_particles,0);
//         double  * gamma = (double*) aligned_alloc(64, number_of_particles*sizeof(double));

        // Computation of the particle gamma factor
        if (mass == 0) {
            #pragma omp simd private(ipr)
//             aligned(gamma : 64)
            for (ip=(unsigned int)(istart) ; ip<(unsigned int) (iend); ip++ ) {

                // Local (relative) array index
                ipr = ip - istart;

                gamma[ipr] = sqrt(momentum[0][ip]*momentum[0][ip]
                              + momentum[1][ip]*momentum[1][ip]
                              + momentum[2][ip]*momentum[2][ip]);

            }
        } else {
            #pragma omp simd private(ipr)
//             aligned(gamma : 64)
            for (ip=(unsigned int)(istart) ; ip<(unsigned int) (iend); ip++ ) {

                // Local (relative) array index
                ipr = ip - istart;

                gamma[ipr] = sqrt(1.0 + momentum[0][ip]*momentum[0][ip]
                              + momentum[1][ip]*momentum[1][ip]
                              + momentum[2][ip]*momentum[2][ip]);

            }
        }

        // Computation of the maxima and minima for each direction
        momentum_min[0] = momentum[0][istart];
        momentum_max[0] = momentum[0][istart];

        momentum_min[1] = momentum[1][istart];
        momentum_max[1] = momentum[1][istart];

        momentum_min[2] = momentum[2][istart];
        momentum_max[2] = momentum[2][istart];

#if __INTEL_COMPILER > 18000
        #pragma omp simd \
        reduction(min:momentum_min)  \
        reduction(max:momentum_max)
#endif
        for (ip=(unsigned int) (istart) ; ip < (unsigned int) (iend); ip++ ) {
            momentum_min[0] = std::min(momentum_min[0],momentum[0][ip]);
            momentum_max[0] = std::max(momentum_max[0],momentum[0][ip]);

            momentum_min[1] = std::min(momentum_min[1],momentum[1][ip]);
            momentum_max[1] = std::max(momentum_max[1],momentum[1][ip]);

            momentum_min[2] = std::min(momentum_min[2],momentum[2][ip]);
            momentum_max[2] = std::max(momentum_max[2],momentum[2][ip]);
        }

        // std::cerr << " momentum_min[0]: " << momentum_min[0]
        //           << " momentum_max[0]: " << momentum_max[0]
        //           << " momentum_min[1]: " << momentum_min[1]
        //           << " momentum_max[1]: " << momentum_max[1]
        //           << " momentum_min[2]: " << momentum_min[2]
        //           << " momentum_max[2]: " << momentum_max[2]
        //           << std::endl;

        // Extra to include the max in the discretization
        // mr_max += (mr_max - mr_min)*0.01;
        // theta_max += (theta_max - theta_min)*0.01;
        // phi_max += (phi_max - phi_min)*0.01;

        // Computation of the deltas (discretization steps)
        // Check if min and max boundaries are very close in each momentum direction
        for (ip = 0 ; ip < 3 ; ip++) {
            if (fabs((momentum_max[ip] - momentum_min[ip])) < min_momentum_cell_length_[ip]) {
                // If momentum_min[ip] and momentum_max[ip] have the same sign
                if (momentum_max[ip] <= 0 || momentum_min[ip] >= 0) {
                    momentum_delta[ip] = min_momentum_cell_length_[ip];
                    momentum_min[ip] = (momentum_max[ip] + momentum_min[ip] - momentum_delta[ip])*0.5;
                    momentum_max[ip] = (momentum_max[ip] + momentum_min[ip] + momentum_delta[ip])*0.5;
                    inv_momentum_delta[ip] = 0;
                    dim[ip] = 1;
                // Else momentum_min[ip] and momentum_max[ip] have opposite signs,
                } else {
                    momentum_max[ip] = std::max(fabs((momentum_max[ip] + momentum_min[ip] + min_momentum_cell_length_[ip])*0.5),fabs((momentum_max[ip] + momentum_min[ip] - min_momentum_cell_length_[ip])*0.5));
                    momentum_min[ip] = -momentum_max[ip];
                    momentum_delta[ip] = momentum_max[ip];
                    inv_momentum_delta[ip] = 1.0/momentum_delta[ip];
                    dim[ip] = 2;
                }
            } else {
                // If the user request a single cell in this direction
                if (dim[ip] == 1) {
                    momentum_max[ip] += (momentum_max[ip] - momentum_min[ip])*0.01;
                    momentum_delta[ip] = (momentum_max[ip] - momentum_min[ip]);
                    inv_momentum_delta[ip] = 1.0/momentum_delta[ip];
                // If momentum_min[ip] and momentum_max[ip] have the same sign
                } else if (momentum_max[ip] <= 0 || momentum_min[ip] >= 0) {
                    momentum_max[ip] += (momentum_max[ip] - momentum_min[ip])*0.01;
                    if (accumulation_correction_) {
                        momentum_delta[ip] = (momentum_max[ip] - momentum_min[ip]) / (dim[ip]-1);
                        //momentum_min[ip] -= 0.99*momentum_delta[ip]*Rand::uniform();
                        momentum_min[ip] -= 0.99*momentum_delta[ip]*rand_->uniform();
                    } else {
                        momentum_delta[ip] = (momentum_max[ip] - momentum_min[ip]) / (dim[ip]);
                    }
                    inv_momentum_delta[ip] = 1.0/momentum_delta[ip];
                // Else momentum_min[ip] and momentum_max[ip] have opposite signs,
                // The 0 value is at the boundary between 2 cells
                } else {
                    if (accumulation_correction_) {
                        //dim[ip] = int(dim[ip]*(1+Rand::uniform()));
                        dim[ip] = int(dim[ip]*(1+rand_->uniform()));
                    }
                    momentum_delta[ip] = fabs(momentum_max[ip] - momentum_min[ip]) / dim[ip];
                    inv_momentum_delta[ip] = 1.0/momentum_delta[ip];
                    nb_delta = ceil(fabs(momentum_min[ip]) * inv_momentum_delta[ip]);
                    momentum_min[ip] = -nb_delta*momentum_delta[ip];
                    dim[ip] += 1;
                    momentum_max[ip] = momentum_min[ip] + dim[ip] * momentum_delta[ip];
                }
            }
        }

        // if (fabs(momentum_max[1] - momentum_min[1]) < min_momentum_cell_length_[1]) {
        //     if (momentum_max[1] <= 0 || momentum_min[1] >= 0) {
        //         momentum_delta[1] = min_momentum_cell_length_[1];
        //         momentum_min[1] = (momentum_max[1] + momentum_min[1] - momentum_delta[1])*0.5;
        //         momentum_max[1] = (momentum_max[1] + momentum_min[1] + momentum_delta[1])*0.5;
        //         inv_momentum_delta[1] = 0;
        //         dim[1] = 1;
        //     } else {
        //         momentum_max[1] = std::max(fabs((momentum_max[1] + momentum_min[1] + min_momentum_cell_length_[1])*0.5),fabs((momentum_max[1] + momentum_min[1] - min_momentum_cell_length_[1])*0.5));
        //         momentum_min[1] = -momentum_max[1];
        //         momentum_delta[1] = momentum_max[1];
        //         inv_momentum_delta[1] = 1.0/momentum_delta[1];
        //         dim[1] = 2;
        //     }
        // } else {
        //     if (dim[1] == 1) {
        //         momentum_max[1] += (momentum_max[1] - momentum_min[1])*0.01;
        //         momentum_delta[1] = (momentum_max[1] - momentum_min[1]);
        //         inv_momentum_delta[1] = 1.0/momentum_delta[1];
        //     // If momentum_min[1] and momentum_max[1] have the same sign
        //     } else if (momentum_max[1] <= 0 || momentum_min[1] >= 0) {
        //         momentum_max[1] += (momentum_max[1] - momentum_min[1])*0.01;
        //         momentum_delta[1] = (momentum_max[1] - momentum_min[1]) / (dim[1]-1);
        //         momentum_min[1] -= 0.99*momentum_delta[1]*Rand::uniform();
        //         inv_momentum_delta[1] = 1.0/momentum_delta[1];
        //     // else, discretization centerd in 0
        //     } else {
        //         dim[1] = int(dim[1]*(1+Rand::uniform()));
        //         momentum_delta[1] = fabs(momentum_max[1] - momentum_min[1]) / dim[1];
        //         inv_momentum_delta[1] = 1.0/momentum_delta[1];
        //         nb_delta = ceil(fabs(momentum_min[1]) * inv_momentum_delta[1]);
        //         momentum_min[1] = -nb_delta*momentum_delta[1];
        //         dim[1] += 1;
        //         momentum_max[1] = momentum_min[1] + dim[1] * momentum_delta[1];
        //     }
        // }
        
        // std::cerr << std::scientific << std::setprecision(15)
        //           << " My centering: "
        //           << " my interval: " << fabs(momentum_max[1] - momentum_min[1])
        //           << " min_momentum_cell_length_[1]: " << min_momentum_cell_length_[1]
        //           << " momentum_min[1]: " << momentum_min[1]
        //           << " momentum_max[1]: " << momentum_max[1]
        //           << " momentum_delta[1]: " << momentum_delta[1]
        //           << " my_dim: " << dim[1]
        //           << " nb_delta: " << nb_delta
        //           << std::endl;
        


        // std::cerr << std::scientific << std::setprecision(15)
        //           << " Mz centering: "
        //           << " mz interval: " << fabs(momentum_max[2] - momentum_min[2])
        //           << " min_momentum_cell_length_[0]: " << min_momentum_cell_length_[0]
        //           << " min_momentum_cell_length_[1]: " << min_momentum_cell_length_[1]
        //           << " min_momentum_cell_length_[2]: " << min_momentum_cell_length_[2]
        //           << " momentum_min[2]: " << momentum_min[2]
        //           << " momentum_max[2]: " << momentum_max[2]
        //           << " momentum_delta[2]: " << momentum_delta[2]
        //           << " mz_dim: " << dim[2]
        //           << " nb_delta: " << nb_delta
        //           << std::endl;

        // Total number of momentum cells
        unsigned int momentum_cells = dim[0]
                                    * dim[1]
                                    * dim[2];

        // Array containing the number of particles per momentum cells
        std::vector <unsigned int> particles_per_momentum_cells(momentum_cells,0);
//         unsigned int  * particles_per_momentum_cells = (unsigned int*) aligned_alloc(64, momentum_cells*sizeof(unsigned int));

        // Array containing the first particle index of each momentum cell
        // in the sorted particle array
        std::vector <unsigned int> momentum_cell_particle_index(momentum_cells,0);
//         unsigned int  * momentum_cell_particle_index = (unsigned int*) aligned_alloc(64, momentum_cells*sizeof(unsigned int));
//
//         // Initialization when using aligned_alloc
//         for (ic = 0 ; ic < momentum_cells ; ic++) {
//             momentum_cell_particle_index[ic] = 0;
//             particles_per_momentum_cells[ic] = 0;
//         }

        // std::cerr << "Cell index" << std::endl;

        // For each particle, momentum cell indexes are computed in the
        // requested discretization.
        // This loop can be efficiently vectorized
        #pragma omp simd \
        private(ipr,mx_i,my_i,mz_i)
//         aligned(momentum_cell_index: 64)
        for (ip=(unsigned int) (istart) ; ip < (unsigned int) (iend); ip++ ) {

            // Relative particle array index
            ipr = ip - istart;

            // 3d indexes in the momentum discretization
            mx_i = (unsigned int) floor( (momentum[0][ip] - momentum_min[0]) * inv_momentum_delta[0]);
            my_i = (unsigned int) floor( (momentum[1][ip] - momentum_min[1]) * inv_momentum_delta[1]);
            mz_i = (unsigned int) floor( (momentum[2][ip] - momentum_min[2]) * inv_momentum_delta[2]);

            // 1D Index in the momentum discretization
            momentum_cell_index[ipr] = mz_i * dim[0]*dim[1]
                          + my_i * dim[0] + mx_i;

                // std::cerr << "momentum_cell_index: " << momentum_cell_index[ipr]
                //           << " / " << dim[0]*dim[1]*dim[2]
                //           << " mx_i: " << mx_i
                //           << " my_i: " << my_i
                //           << " mz_i: " << mz_i
                //           << " / " << dim[2]
                //           << " mz: " << momentum[2][ip]
                //           << " momentum_min[2]: " << momentum_min[2]
                //           << " momentum_max[2]: " << momentum_max[2]
                //           << " momentum_delta[2]: " << momentum_delta[2]
                //           << std::endl;

        }

        // std::cerr << "number of particles per momentum cells" << std::endl;

        // The number of particles per momentum cells
        // is then determined.
        // No vectorization because of random memory accesses
        for (ipr=0; ipr<number_of_particles; ipr++ ) {
            // Number of particles per momentum cells
            particles_per_momentum_cells[momentum_cell_index[ipr]] += 1;
        }

        // std::cerr << "Computation of the cell index" << std::endl;

        // Computation of the cell index in the sorted array of particles
        for (ic = 1 ; ic < momentum_cells ; ic++) {
            momentum_cell_particle_index[ic]  = momentum_cell_particle_index[ic-1] + particles_per_momentum_cells[ic-1];
            particles_per_momentum_cells[ic-1] = 0;
        }
        particles_per_momentum_cells[momentum_cells-1] = 0;

        // std::cerr << "sort particles in correct bins" << std::endl;

        // sort particles in correct bins according to their
        // momentum properties
        // Not vectorizable
        for (ipr=0 ; ipr<number_of_particles; ipr++ ) {

            // Momentum cell index for this particle
            ic = momentum_cell_index[ipr];

            // std::cerr << " ic: " << ic
            //           << " ipr: " << ipr
            //           << " momentum_cell_particle_index[ic]: " << momentum_cell_particle_index[ic]
            //           << std::endl;

            sorted_particles[momentum_cell_particle_index[ic]
            + particles_per_momentum_cells[ic]] = istart + ipr;

            particles_per_momentum_cells[ic] += 1;
        }


        // For each momentum bin, merge packet of particles composed of
        // at least `min_packet_size_` and `max_packet_size_`

        // Loop over the the momentum cells that have enough particules
        for (mz_i=0 ; mz_i< dim[2]; mz_i++ ) {

            cell_vec_z = momentum_min[2] + (mz_i+0.5)*momentum_delta[2];

            for (my_i=0 ; my_i< dim[1]; my_i++ ) {

                cell_vec_y = momentum_min[1] + (my_i+0.5)*momentum_delta[1];

                // 1D cell direction index
                icc = my_i + mz_i* dim[1] ;

                for (mx_i=0 ; mx_i< dim[0]; mx_i++ ) {

                    cell_vec_x = momentum_min[0] + (mx_i+0.5)*momentum_delta[0];

                    // 1D cell index
                    ic = mx_i + icc*dim[0];

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
                            total_momentum_x = 0;
                            total_momentum_y = 0;
                            total_momentum_z = 0;
                            total_energy = 0;

                            // First index of the packet
                            ipr_min = ipack*max_packet_size_;
                            // last index of the packet
                            ipr_max = std::min((ipack+1)*max_packet_size_,particles_per_momentum_cells[ic]);

                            // _______________________________________________________________

                            // Compute total weight, total momentum and total energy

                            for (ipr = ipr_min ; ipr < ipr_max ; ipr ++) {

                                // Particle index in Particles
                                ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];

                                // Total weight (wt)
                                total_weight += weight[ip];

                                // total momentum vector (pt)
                                total_momentum_x += momentum[0][ip]*weight[ip];
                                total_momentum_y += momentum[1][ip]*weight[ip];
                                total_momentum_z += momentum[2][ip]*weight[ip];

                                // total energy
                                total_energy += weight[ip]*gamma[ip - istart];

                            }

                            // \varepsilon_a in Vranic et al
                            new_energy = total_energy / total_weight;

                            // pa in Vranic et al.
                            // For photons
                            if (mass == 0) {
                                new_momentum_norm = new_energy;
                            // For mass particles
                            } else {
                                new_momentum_norm = sqrt(new_energy*new_energy - 1.0);
                            }

                            // Total momentum norm
                            total_momentum_norm = sqrt(total_momentum_x*total_momentum_x
                                                +      total_momentum_y*total_momentum_y
                                                +      total_momentum_z*total_momentum_z);

                            // Angle between pa and pt, pb and pt in Vranic et al.
                            cos_omega = std::min(total_momentum_norm / (total_weight*new_momentum_norm),1.0);
                            sin_omega = sqrt(1 - cos_omega*cos_omega);

                            // Now, represents the inverse to avoid useless division
                            total_momentum_norm = 1/total_momentum_norm;

                            // Computation of e1 unit vector
                            e1_x = total_momentum_x*total_momentum_norm;
                            e1_y = total_momentum_y*total_momentum_norm;
                            e1_z = total_momentum_z*total_momentum_norm;

                            // e3 = e1 x cell_vec
                            e3_x = e1_y*cell_vec_z - e1_z*cell_vec_y;
                            e3_y = e1_z*cell_vec_x - e1_x*cell_vec_z;
                            e3_z = e1_x*cell_vec_y - e1_y*cell_vec_x;

                            // All particle momenta are not collinear
                            if (fabs(e3_x*e3_x + e3_y*e3_y + e3_z*e3_z) > 0)
                            {
                                
                                // Computation of e2  = e1 x e3 unit vector
                                // e2_x = e1_y*e1_y*cell_vec_x
                                //      - e1_x * (e1_y*cell_vec_y + e1_z*cell_vec_z)
                                //      + e1_z*e1_z*cell_vec_x;
                                // e2_y = e1_z*e1_z*cell_vec_y
                                //      - e1_y * (e1_z*cell_vec_z + e1_x*cell_vec_x)
                                //      + e1_x*e1_x*cell_vec_y;
                                // e2_z = e1_x*e1_x*cell_vec_z
                                //      - e1_z * (e1_x*cell_vec_x + e1_y*cell_vec_y)
                                //      + e1_y*e1_y*cell_vec_z;

                                e2_x = e1_y*e3_z - e1_z*e3_y;
                                e2_y = e1_z*e3_x - e1_x*e3_z;
                                e2_z = e1_x*e3_y - e1_y*e3_x;
                            
                                e2_norm = 1./sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z);

                                // e2 is normalized to be a unit vector
                                e2_x = e2_x * e2_norm;
                                e2_y = e2_y * e2_norm;
                                e2_z = e2_z * e2_norm;

                                // double mx1 = new_momentum_norm*(cos_omega*e1_x + sin_omega*e2_x);
                                // double mx2 = new_momentum_norm*(cos_omega*e1_x - sin_omega*e2_x);
                                // total_momentum_y = new_momentum_norm*(cos_omega*e1_y + sin_omega*e2_y);
                                // double my2 = new_momentum_norm*(cos_omega*e1_y - sin_omega*e2_y);
                                // total_momentum_z = new_momentum_norm*(cos_omega*e1_z + sin_omega*e2_z);
                                // double mz2 = new_momentum_norm*(cos_omega*e1_z - sin_omega*e2_z);
                                // if (isnan(total_momentum_x)
                                //     || isnan(total_momentum_y)
                                //     || isnan(total_momentum_z)
                                //     // || mx1 < momentum_min[0] + mx_i*momentum_delta[0]
                                //     // || total_momentum_y < momentum_min[1] + my_i*momentum_delta[1]
                                //     // || total_momentum_z < momentum_min[2] + mz_i*momentum_delta[2]
                                //     // || mx1 > momentum_min[0] + (mx_i+1)*momentum_delta[0]
                                //     // || total_momentum_y > momentum_min[1] + (my_i+1)*momentum_delta[1]
                                //     // || total_momentum_z > momentum_min[2] + (mz_i+1)*momentum_delta[2]
                                //     // || fabs(new_energy - sqrt(mx1*mx1
                                //     //                   + total_momentum_y*total_momentum_y
                                //     //                   + total_momentum_z*total_momentum_z)) > 1e-7
                                //     ) {
                                //
                                //     for (ipr = ipr_min; ipr < ipr_max ; ipr ++) {
                                //         ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                //         std::cerr
                                //         << " mx: " << momentum[0][ip]
                                //         << " my: " << momentum[1][ip]
                                //         << " mz: " << momentum[2][ip]
                                //         << std::endl;
                                //     }
                                //
                                //     ip = sorted_particles[momentum_cell_particle_index[ic] + ipr_min];
                                //
                                //     //std::cerr <<
                                //     ERROR(std::scientific << std::setprecision(15)
                                //       << " dim: " << dim[0] << " " << dim[1] << " " << dim[2] << "\n"
                                //       << " mx1: " << mx1
                                //       << " mx2: " << mx2
                                //       << " mx[i]: " << momentum_min[0] + mx_i*momentum_delta[0]
                                //       << " mx[i+1]: " << momentum_min[0] + (mx_i+1)*momentum_delta[0]
                                //       << " momentum_delta[0]: " << momentum_delta[0]
                                //       << "\n"
                                //       << " my1: " << total_momentum_y
                                //       << " my2: " << my2
                                //       << " my[i]: " << momentum_min[1] + my_i*momentum_delta[1]
                                //       << " my[i+1]: " << momentum_min[1] + (my_i+1)*momentum_delta[1] << "\n"
                                //       << " mz: " << total_momentum_z
                                //       << " mz2: " << mz2
                                //       << " mz[i]: " << momentum_min[2] + mz_i*momentum_delta[2]
                                //       << " mz[i+1]: " << momentum_min[2] + (mz_i+1)*momentum_delta[2] << "\n"
                                //       << " total_weight: " << 0.5*total_weight
                                //       << " energy: " << total_energy / total_weight
                                //       << " energy: " << sqrt(mx1*mx1
                                //                         + total_momentum_y*total_momentum_y
                                //                         + total_momentum_z*total_momentum_z)
                                //       << " " << fabs(new_energy - sqrt(total_momentum_x*total_momentum_x
                                //                         + total_momentum_y*total_momentum_y
                                //                         + total_momentum_z*total_momentum_z)) << "\n"
                                //       << " total energy: " << 0.5*total_energy
                                //       << " total energy: " << 0.5*total_weight*sqrt(mx1*mx1
                                //                         + total_momentum_y*total_momentum_y
                                //                         + total_momentum_z*total_momentum_z)
                                //       << " " << fabs(total_energy - total_weight*sqrt(mx1*mx1
                                //                         + total_momentum_y*total_momentum_y
                                //                         + total_momentum_z*total_momentum_z)) << "\n"
                                //     << " omega: " << std::acos(std::min(total_momentum_norm / (total_weight*new_momentum_norm),1.0))
                                //     << " " << total_momentum_norm / (total_weight*new_momentum_norm)
                                //     << " cos_omega: " << cos_omega << " sin_omega: " << sin_omega << "\n"
                                //     << " cell_vec: " << cell_vec_x << " " << cell_vec_y << " " << cell_vec_z << "\n"
                                //     << " e1: " << e1_x << " " << e1_y << " " << e1_z << "\n"
                                //     << " e3: " << e3_x << " " << e3_y << " " << e3_z << "\n"
                                //     << " e2: " << e2_x << " " << e2_y << " " << e2_z << "\n"
                                //     << " e1.e2: " << e1_x*e2_x + e1_y*e2_y + e1_z*e2_z
                                //     )
                                //     //<< std::endl;
                                // }
                                
                                // Create the merged particles
                                // --------------------------------------------
                                
                                // Method 1: determinist - use the position of
                                // the first particles of the list
                                
                                // Update momentum of the first photon
                                ip = sorted_particles[momentum_cell_particle_index[ic] + ipr_min];
                                
                                momentum[0][ip] = new_momentum_norm*(cos_omega*e1_x + sin_omega*e2_x);
                                momentum[1][ip] = new_momentum_norm*(cos_omega*e1_y + sin_omega*e2_y);
                                momentum[2][ip] = new_momentum_norm*(cos_omega*e1_z + sin_omega*e2_z);
                                weight[ip] = 0.5*total_weight;
                                
                                // Update momentum of the second particle
                                ip = sorted_particles[momentum_cell_particle_index[ic] + ipr_min + 1];
                                momentum[0][ip] = new_momentum_norm*(cos_omega*e1_x - sin_omega*e2_x);
                                momentum[1][ip] = new_momentum_norm*(cos_omega*e1_y - sin_omega*e2_y);
                                momentum[2][ip] = new_momentum_norm*(cos_omega*e1_z - sin_omega*e2_z);
                                weight[ip] = 0.5*total_weight;
                                
                                // Other photons are tagged to be removed after
                                for (ipr = ipr_min + 2; ipr < ipr_max ; ipr ++) {
                                    ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                    mask[ip] = -1;
                                    count--;
                                }
                                
                                // Method 2: random - pick up randomly old particle positions
                                
                                // unsigned int ipr1 = ipr_min + int(Rand::uniform()*(ipr_max - ipr_min));
                                // unsigned int ipr2 = ipr_min + int(Rand::uniform()*(ipr_max - ipr_min));
                                // while (ipr1 == ipr2)
                                // {
                                //     ipr2 = ipr_min + int(Rand::uniform()*(ipr_max - ipr_min));
                                // }
                                //
                                // for (ipr = ipr_min; ipr < ipr_max ; ipr ++) {
                                //     if (ipr == ipr1) {
                                //         ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                //         momentum[0][ip] = new_momentum_norm*(cos_omega*e1_x + sin_omega*e2_x);
                                //         momentum[1][ip] = new_momentum_norm*(cos_omega*e1_y + sin_omega*e2_y);
                                //         momentum[2][ip] = new_momentum_norm*(cos_omega*e1_z + sin_omega*e2_z);
                                //         weight[ip] = 0.5*total_weight;
                                //     } else if (ipr == ipr2) {
                                //         ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                //         momentum[0][ip] = new_momentum_norm*(cos_omega*e1_x - sin_omega*e2_x);
                                //         momentum[1][ip] = new_momentum_norm*(cos_omega*e1_y - sin_omega*e2_y);
                                //         momentum[2][ip] = new_momentum_norm*(cos_omega*e1_z - sin_omega*e2_z);
                                //         weight[ip] = 0.5*total_weight;
                                //     } else {
                                //         ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                //         mask[ip] = -1;
                                //         count--;
                                //     }
                                // }
                                
                            // Special treatment for collinear photons
                            // Collinear particles are merged
                            } else {
                                
                                if (mass == 0) {
                                    // Method 1: determinist - use the position of
                                    // the first particles of the list
                                    
                                    // Update momentum of the first photon
                                    ip = sorted_particles[momentum_cell_particle_index[ic] + ipr_min];
                                    
                                    momentum[0][ip] = new_momentum_norm*e1_x;
                                    momentum[1][ip] = new_momentum_norm*e1_y;
                                    momentum[2][ip] = new_momentum_norm*e1_z;
                                    weight[ip] = total_weight;
                                    
                                    // Other photons are tagged to be removed after
                                    for (ipr = ipr_min + 1; ipr < ipr_max ; ipr ++) {
                                        ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                        mask[ip] = -1;
                                        count--;
                                    }
                                    
                                    // Method 2: random - pick up randomly old particle positions
                                    
                                    // unsigned int ipr1 = ipr_min + int(Rand::uniform()*(ipr_max - ipr_min));
                                    // for (ipr = ipr_min; ipr < ipr_max ; ipr ++) {
                                    //     if (ipr == ipr1) {
                                    //         ip = sorted_particles[momentum_cell_particle_index[ic] + ipr1];
                                    //         momentum[0][ip] = new_momentum_norm*e1_x;
                                    //         momentum[1][ip] = new_momentum_norm*e1_y;
                                    //         momentum[2][ip] = new_momentum_norm*e1_z;
                                    //         weight[ip] = total_weight;
                                    //     } else {
                                    //         ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                    //         mask[ip] = -1;
                                    //         count--;
                                    //     }
                                    // }
                                }
                                
                            }// end check collinear momenta

                        }
                    }

                }
            }
        }

        // Free aligned arrays
//         free(gamma);
//         free(momentum_cell_index);
//         free(sorted_particles);
//         free(particles_per_momentum_cells);
//         free(momentum_cell_particle_index);
    }

}
