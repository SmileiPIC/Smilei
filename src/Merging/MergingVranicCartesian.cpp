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

// -----------------------------------------------------------------------------
//! Constructor for RadiationNLandauLifshitz
//! Inherited from Radiation
// -----------------------------------------------------------------------------
MergingVranicCartesian::MergingVranicCartesian(Params& params,
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

    // Minimum momentum cell length
    min_momentum_cell_length_[0] = species->merge_min_momentum_cell_length_[0];
    min_momentum_cell_length_[1] = species->merge_min_momentum_cell_length_[1];
    min_momentum_cell_length_[2] = species->merge_min_momentum_cell_length_[2];
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
    if (number_of_particles > min_particles_per_cell) {

        // Momentum discretization
        unsigned int dim[3];
        for (unsigned int i = 0; i < 3 ; i++) {
            dim[i] = dimensions_[i];
        }

        // Minima
        double mx_min;
        double my_min;
        double mz_min;

        // Maxima
        double mx_max;
        double my_max;
        double mz_max;

        // Delta
        double mx_delta;
        double my_delta;
        double mz_delta;

        // Inverse Delta
        double inv_mx_delta;
        double inv_my_delta;
        double inv_mz_delta;

        // Angles
        double omega;
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
        // std::vector <unsigned int> momentum_cell_index(number_of_particles,0);
        unsigned int  * momentum_cell_index = (unsigned int*) aligned_alloc(64, number_of_particles*sizeof(unsigned int));

        // Sorted array of particle index
        // std::vector <unsigned int> sorted_particles(number_of_particles,0);
        unsigned int  * sorted_particles = (unsigned int*) aligned_alloc(64, number_of_particles*sizeof(unsigned int));

        // Particle gamma factor
        // std::vector <double> gamma(number_of_particles,0);
        double  * gamma = (double*) aligned_alloc(64, number_of_particles*sizeof(double));

        // Computation of the particle gamma factor
        if (mass > 0) {

            #pragma omp simd aligned(gamma : 64) private(ipr)
            for (ip=(unsigned int)(istart) ; ip<(unsigned int) (iend); ip++ ) {

                // Local (relative) array index
                ipr = ip - istart;

                gamma[ipr] = sqrt(1.0 + momentum[0][ip]*momentum[0][ip]
                              + momentum[1][ip]*momentum[1][ip]
                              + momentum[2][ip]*momentum[2][ip]);

            }
        }

        // Computation of the maxima and minima for each direction
        mx_min = momentum[0][istart];
        mx_max = momentum[0][istart];

        my_min = momentum[1][istart];
        my_max = momentum[1][istart];

        mz_min = momentum[2][istart];
        mz_max = momentum[2][istart];

        #pragma omp simd \
        reduction(min:mx_min) reduction(min:my_min) reduction(min:mz_min) \
        reduction(max:mx_max) reduction(max:my_max) reduction(max:mz_max)
        for (ip=(unsigned int) (istart) ; ip < (unsigned int) (iend); ip++ ) {
            mx_min = fmin(mx_min,momentum[0][ip]);
            mx_max = fmax(mx_max,momentum[0][ip]);

            my_min = fmin(my_min,momentum[1][ip]);
            my_max = fmax(my_max,momentum[1][ip]);

            mz_min = fmin(mz_min,momentum[2][ip]);
            mz_max = fmax(mz_max,momentum[2][ip]);
        }

        // std::cerr << " mx_min: " << mx_min
        //           << " mx_max: " << mx_max
        //           << " my_min: " << my_min
        //           << " my_max: " << my_max
        //           << " mz_min: " << mz_min
        //           << " mz_max: " << mz_max
        //           << std::endl;

        // Extra to include the max in the discretization
        // mr_max += (mr_max - mr_min)*0.01;
        // theta_max += (theta_max - theta_min)*0.01;
        // phi_max += (phi_max - phi_min)*0.01;

        // Computation of the deltas (discretization steps)
        // Check if min and max boundaries are very close
        if (fabs((mx_max - mx_min)) < min_momentum_cell_length_[0]) {
            if (mx_max <= 0 || mx_min >= 0) {
                mx_delta = min_momentum_cell_length_[0];
                mx_min = (mx_max + mx_min - mx_delta)*0.5;
                mx_max = (mx_max + mx_min + mx_delta)*0.5;
                inv_mx_delta = 0;
                dim[0] = 1;
            } else {
                mx_max = std::max(fabs((mx_max + mx_min + min_momentum_cell_length_[0])*0.5),fabs((mx_max + mx_min - min_momentum_cell_length_[0])*0.5));
                mx_min = -mx_max;
                mx_delta = mx_max;
                inv_mx_delta = 1.0/mx_delta;
                dim[0] = 2;
            }
        } else {
            // 1 element in this direction
            if (dim[0] == 1) {
                mx_max += (mx_max - mx_min)*0.01;
                mx_delta = (mx_max - mx_min);
                inv_mx_delta = 1.0/mx_delta;
            // If mx_min and mx_max have the same sign
            } else if (mx_max <= 0 || mx_min >= 0) {
                mx_max += (mx_max - mx_min)*0.01;
                mx_delta = (mx_max - mx_min) / (dim[0]-1);
                mx_min -= 0.99*mx_delta*Rand::uniform();
                inv_mx_delta = 1.0/mx_delta;
            // else mx_min and mx_max have different signs,
            // discretization centerd in 0
            } else {
                dim[0] = int(dim[0]*(1+Rand::uniform()));
                mx_delta = fabs(mx_max - mx_min) / dim[0];
                inv_mx_delta = 1.0/mx_delta;
                nb_delta = ceil(fabs(mx_min) * inv_mx_delta);
                mx_min = -nb_delta*mx_delta;
                dim[0] += 1;
                mx_max = mx_min + dim[0] * mx_delta;
            }
        }

        if (fabs(my_max - my_min) < min_momentum_cell_length_[1]) {
            if (my_max <= 0 || my_min >= 0) {
                my_delta = min_momentum_cell_length_[1];
                my_min = (my_max + my_min - my_delta)*0.5;
                my_max = (my_max + my_min + my_delta)*0.5;
                inv_my_delta = 0;
                dim[1] = 1;
            } else {
                my_max = std::max(fabs((my_max + my_min + min_momentum_cell_length_[1])*0.5),fabs((my_max + my_min - min_momentum_cell_length_[1])*0.5));
                my_min = -my_max;
                my_delta = my_max;
                inv_my_delta = 1.0/my_delta;
                dim[1] = 2;
            }
        } else {
            if (dim[1] == 1) {
                my_max += (my_max - my_min)*0.01;
                my_delta = (my_max - my_min);
                inv_my_delta = 1.0/my_delta;
            // If my_min and my_max have the same sign
            } else if (my_max <= 0 || my_min >= 0) {
                my_max += (my_max - my_min)*0.01;
                my_delta = (my_max - my_min) / (dim[1]-1);
                my_min -= 0.99*my_delta*Rand::uniform();
                inv_my_delta = 1.0/my_delta;
            // else, discretization centerd in 0
            } else {
                dim[1] = int(dim[1]*(1+Rand::uniform()));
                my_delta = fabs(my_max - my_min) / dim[1];
                inv_my_delta = 1.0/my_delta;
                nb_delta = ceil(fabs(my_min) * inv_my_delta);
                my_min = -nb_delta*my_delta;
                dim[1] += 1;
                my_max = my_min + dim[1] * my_delta;
            }
        }
        
        // std::cerr << std::scientific << std::setprecision(15)
        //           << " My centering: "
        //           << " my interval: " << fabs(my_max - my_min)
        //           << " min_momentum_cell_length_[1]: " << min_momentum_cell_length_[1]
        //           << " my_min: " << my_min
        //           << " my_max: " << my_max
        //           << " my_delta: " << my_delta
        //           << " my_dim: " << dim[1]
        //           << " nb_delta: " << nb_delta
        //           << std::endl;
        
        // Momentum z direction
        if (fabs(mz_max - mz_min) < min_momentum_cell_length_[2]) {
            // If mz_min and mz_max have the same sign
            if (mz_max <= 0 || mz_min >= 0) {
                mz_delta = min_momentum_cell_length_[2];
                mz_min = (mz_max + mz_min - mz_delta)*0.5;
                mz_max = (mz_max + mz_min + mz_delta)*0.5;
                inv_mz_delta = 0;
                dim[2] = 1;
            // else if mz_min and mz_max does not have the same sign,
            // discretization centerd in 0 and dim = 2 instead of 1
            } else {
                mz_max = std::max(fabs((mz_max + mz_min + min_momentum_cell_length_[1])*0.5),fabs((mz_max + mz_min - min_momentum_cell_length_[2])*0.5));
                mz_min = -mz_max;
                mz_delta = mz_max;
                inv_mz_delta = 1.0/mz_delta;
                dim[2] = 2;
            }
        } else {
            if (dim[2] == 1) {
                mz_max += (mz_max - mz_min)*0.01;
                mz_delta = (mz_max - mz_min);
                inv_mz_delta = 1.0/mz_delta;
            // If mz_min and mz_max have the same sign
            } else if (mz_max <= 0 || mz_min >= 0) {
                mz_max += (mz_max - mz_min)*0.01;
                mz_delta = (mz_max - mz_min) / (dim[2]-1);
                mz_min -= 0.99*mz_delta*Rand::uniform();
                inv_mz_delta = 1.0/mz_delta;
            // else if mz_min and mz_max does not have the same sign,
            // discretization centerd in 0
            } else {
                dim[2] = int(dim[2]*(1+Rand::uniform()));
                mz_delta = fabs(mz_max - mz_min) / dim[2];
                inv_mz_delta = 1.0/mz_delta;
                nb_delta = ceil(fabs(mz_min) * inv_mz_delta);
                mz_min = -nb_delta*mz_delta;
                dim[2] += 1;
                mz_max = mz_min + dim[2] * mz_delta;
            }
        }

        // std::cerr << std::scientific << std::setprecision(15)
        //           << " Mz centering: "
        //           << " mz interval: " << fabs(mz_max - mz_min)
        //           << " min_momentum_cell_length_[0]: " << min_momentum_cell_length_[0]
        //           << " min_momentum_cell_length_[1]: " << min_momentum_cell_length_[1]
        //           << " min_momentum_cell_length_[2]: " << min_momentum_cell_length_[2]
        //           << " mz_min: " << mz_min
        //           << " mz_max: " << mz_max
        //           << " mz_delta: " << mz_delta
        //           << " mz_dim: " << dim[2]
        //           << " nb_delta: " << nb_delta
        //           << std::endl;

        // Total number of momentum cells
        unsigned int momentum_cells = dim[0]
                                    * dim[1]
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

        // std::cerr << "Cell index" << std::endl;

        // For each particle, momentum cell indexes are computed in the
        // requested discretization.
        // This loop can be efficiently vectorized
        #pragma omp simd \
        private(ipr,mx_i,my_i,mz_i)     \
        aligned(momentum_cell_index: 64)
        for (ip=(unsigned int) (istart) ; ip < (unsigned int) (iend); ip++ ) {

            // Relative particle array index
            ipr = ip - istart;

            // 3d indexes in the momentum discretization
            mx_i = (unsigned int) floor( (momentum[0][ip] - mx_min) * inv_mx_delta);
            my_i = (unsigned int) floor( (momentum[1][ip] - my_min) * inv_my_delta);
            mz_i = (unsigned int) floor( (momentum[2][ip] - mz_min) * inv_mz_delta);

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
                //           << " mz_min: " << mz_min
                //           << " mz_max: " << mz_max
                //           << " mz_delta: " << mz_delta
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

            cell_vec_z = mz_min + (mz_i+0.5)*mz_delta;

            for (my_i=0 ; my_i< dim[1]; my_i++ ) {

                cell_vec_y = my_min + (my_i+0.5)*my_delta;

                // 1D cell direction index
                icc = my_i + mz_i* dim[1] ;

                for (mx_i=0 ; mx_i< dim[0]; mx_i++ ) {

                    cell_vec_x = mx_min + (mx_i+0.5)*mx_delta;

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
                            // Photons
                            if (mass == 0) {

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
                                    total_energy += weight[ip]*sqrt(momentum[0][ip]*momentum[0][ip]
                                                      + momentum[1][ip]*momentum[1][ip]
                                                      + momentum[2][ip]*momentum[2][ip]);

                                }

                                // \varepsilon_a in Vranic et al
                                new_energy = total_energy / total_weight;

                                // pa in Vranic et al.
                                new_momentum_norm = new_energy;

                                // Total momentum norm
                                total_momentum_norm = sqrt(total_momentum_x*total_momentum_x
                                                    +      total_momentum_y*total_momentum_y
                                                    +      total_momentum_z*total_momentum_z);

                                // Angle between pa and pt, pb and pt in Vranic et al.
                                omega = std::acos(std::min(total_momentum_norm / (total_weight*new_momentum_norm),1.0));
                                cos_omega = std::min(total_momentum_norm / (total_weight*new_momentum_norm),1.0);
                                sin_omega = sqrt(1 - cos_omega*cos_omega);

                                // Now, represents the inverse to avoid useless division
                                total_momentum_norm = 1/total_momentum_norm;

                                // Computation of e1 unit vector
                                e1_x = total_momentum_x*total_momentum_norm;
                                e1_y = total_momentum_y*total_momentum_norm;
                                e1_z = total_momentum_z*total_momentum_norm;

                                e3_x = e1_y*cell_vec_z - e1_z*cell_vec_y;
                                e3_y = e1_z*cell_vec_x - e1_x*cell_vec_z;
                                e3_z = e1_x*cell_vec_y - e1_y*cell_vec_x;

                                // Computation of e2  = e1 x e3 unit vector
                                // e3 = e1 x cell_vec
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

                                // momentum[0][ip] = total_momentum_x * new_momentum_norm;
                                // momentum[1][ip] = total_momentum_y * new_momentum_norm;
                                // momentum[2][ip] = total_momentum_z * new_momentum_norm;
                                // weight[ip] = total_weight;

                                // std::cerr << " Energy final: " << sqrt(momentum[0][ip]*momentum[0][ip]
                                //                  + momentum[1][ip]*momentum[1][ip]
                                //                  + momentum[2][ip]*momentum[2][ip])
                                //                  << std::endl;
                                double mx1 = new_momentum_norm*(cos_omega*e1_x + sin_omega*e2_x);
                                double mx2 = new_momentum_norm*(cos_omega*e1_x - sin_omega*e2_x);
                                total_momentum_y = new_momentum_norm*(cos_omega*e1_y + sin_omega*e2_y);
                                double my2 = new_momentum_norm*(cos_omega*e1_y - sin_omega*e2_y);
                                total_momentum_z = new_momentum_norm*(cos_omega*e1_z + sin_omega*e2_z);
                                double mz2 = new_momentum_norm*(cos_omega*e1_z - sin_omega*e2_z);
                                if (isnan(total_momentum_x)
                                    || isnan(total_momentum_y)
                                    || isnan(total_momentum_z)
                                    // || mx1 < mx_min + mx_i*mx_delta
                                    // || total_momentum_y < my_min + my_i*my_delta
                                    // || total_momentum_z < mz_min + mz_i*mz_delta
                                    // || mx1 > mx_min + (mx_i+1)*mx_delta
                                    // || total_momentum_y > my_min + (my_i+1)*my_delta
                                    // || total_momentum_z > mz_min + (mz_i+1)*mz_delta
                                    || fabs(new_energy - sqrt(mx1*mx1
                                                      + total_momentum_y*total_momentum_y
                                                      + total_momentum_z*total_momentum_z)) > 1e-7
                                ) {
                                    
                                    for (ipr = ipr_min; ipr < ipr_max ; ipr ++) {
                                        ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                        std::cerr
                                        << " mx: " << momentum[0][ip]
                                        << " my: " << momentum[1][ip]
                                        << " mz: " << momentum[2][ip]
                                        << std::endl;
                                    }
                                    
                                    ip = sorted_particles[momentum_cell_particle_index[ic] + ipr_min];
                                    
                                    //std::cerr <<
                                    ERROR(std::scientific << std::setprecision(15)
                                              << " dim: " << dim[0] << " " << dim[1] << " " << dim[2] << "\n"
                                              << " mx1: " << mx1
                                              << " mx2: " << mx2
                                              << " mx[i]: " << mx_min + mx_i*mx_delta
                                              << " mx[i+1]: " << mx_min + (mx_i+1)*mx_delta
                                              << " mx_delta: " << mx_delta
                                              << "\n"
                                              << " my1: " << total_momentum_y
                                              << " my2: " << my2
                                              << " my[i]: " << my_min + my_i*my_delta
                                              << " my[i+1]: " << my_min + (my_i+1)*my_delta << "\n"
                                              << " mz: " << total_momentum_z
                                              << " mz2: " << mz2
                                              << " mz[i]: " << mz_min + mz_i*mz_delta
                                              << " mz[i+1]: " << mz_min + (mz_i+1)*mz_delta << "\n"
                                              << " total_weight: " << 0.5*total_weight
                                              << " energy: " << total_energy / total_weight
                                              << " energy: " << sqrt(mx1*mx1
                                                                + total_momentum_y*total_momentum_y
                                                                + total_momentum_z*total_momentum_z)
                                              << " " << fabs(new_energy - sqrt(total_momentum_x*total_momentum_x
                                                                + total_momentum_y*total_momentum_y
                                                                + total_momentum_z*total_momentum_z)) << "\n"
                                              << " total energy: " << 0.5*total_energy
                                              << " total energy: " << 0.5*total_weight*sqrt(mx1*mx1
                                                                + total_momentum_y*total_momentum_y
                                                                + total_momentum_z*total_momentum_z)
                                              << " " << fabs(total_energy - total_weight*sqrt(mx1*mx1
                                                                + total_momentum_y*total_momentum_y
                                                                + total_momentum_z*total_momentum_z)) << "\n"
                                            << " omega: " << omega
                                            << " " << total_momentum_norm / (total_weight*new_momentum_norm)
                                            << " cos_omega: " << cos_omega << " sin_omega: " << sin_omega << "\n"
                                            << " cell_vec: " << cell_vec_x << " " << cell_vec_y << " " << cell_vec_z << "\n"
                                            << " e1: " << e1_x << " " << e1_y << " " << e1_z << "\n"
                                            << " e3: " << e3_x << " " << e3_y << " " << e3_z << "\n"
                                            << " e2: " << e2_x << " " << e2_y << " " << e2_z << "\n"
                                            << " e1.e2: " << e1_x*e2_x + e1_y*e2_y + e1_z*e2_z
                                          )
                                    //<< std::endl;

                                }

                                
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

                            // Mass particles
                            } else {
                                for (ipr = ipr_min ; ipr < ipr_max ; ipr ++) {

                                    // Particle index in Particles
                                    ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];

                                    // Total weight (wt)
                                    total_weight += weight[ip];

                                    // total momentum vector (pt)
                                    total_momentum_x += momentum[0][ip]*weight[ip];
                                    total_momentum_y += momentum[1][ip]*weight[ip];
                                    total_momentum_z += momentum[2][ip]*weight[ip];

                                    // total energy (\varespilon_t)
                                    total_energy += weight[ip]//*gamma[ip-istart];
                                                              * sqrt(1.0 + momentum[0][ip]*momentum[0][ip]
                                                              + momentum[1][ip]*momentum[1][ip]
                                                              + momentum[2][ip]*momentum[2][ip]);

                                }

                                // \varepsilon_a in Vranic et al
                                new_energy = total_energy / total_weight;

                                // pa in Vranic et al.
                                new_momentum_norm = sqrt(new_energy*new_energy - 1.0);

                                // Total momentum norm
                                total_momentum_norm = sqrt(total_momentum_x*total_momentum_x
                                                    +      total_momentum_y*total_momentum_y
                                                    +      total_momentum_z*total_momentum_z);

                                // Angle between pa and pt, pb and pt in Vranic et al.
                                omega = std::acos(std::min(total_momentum_norm / (total_weight*new_momentum_norm),1.0));
                                sin_omega = sin(omega);
                                cos_omega = cos(omega);

                                // std::cerr << "new_energy: " << new_energy
                                //           << " new_momentum_norm: " << new_momentum_norm
                                //           << " total_momentum_norm: " << total_momentum_norm
                                //           << " omega: " << omega
                                //           << std::endl;

                                // Now, represents the inverse to avoid useless division
                                total_momentum_norm = 1/total_momentum_norm;

                                // Computation of e1 unit vector
                                e1_x = total_momentum_x*total_momentum_norm;
                                e1_y = total_momentum_y*total_momentum_norm;
                                e1_z = total_momentum_z*total_momentum_norm;

                                // Computation of e2  = e1 x e3 unit vector
                                // e3 = e1 x cell_vec
                                e2_x = e1_y*e1_y*cell_vec_x
                                     - e1_x * (e1_y*cell_vec_y + e1_z*cell_vec_z)
                                     + e1_z*e1_z*cell_vec_x;
                                e2_y = e1_z*e1_z*cell_vec_y
                                     - e1_y * (e1_z*cell_vec_z + e1_x*cell_vec_x)
                                     + e1_x*e1_x*cell_vec_y;
                                e2_z = e1_x*e1_x*cell_vec_z
                                     - e1_z * (e1_x*cell_vec_x + e1_y*cell_vec_y)
                                     + e1_y*e1_y*cell_vec_z;

                                e2_norm = 1./sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z);

                                // e2 is normalized to be a unit vector
                                e2_x = e2_x * e2_norm;
                                e2_y = e2_y * e2_norm;
                                e2_z = e2_z * e2_norm;

                                // The first 2 particles of the list will
                                // be the merged particles.

                                // Update momentum of the first particle
                                ip = sorted_particles[momentum_cell_particle_index[ic] + ipr_min];
                                momentum[0][ip] = new_momentum_norm*(cos_omega*e1_x + sin_omega*e2_x);
                                momentum[1][ip] = new_momentum_norm*(cos_omega*e1_y + sin_omega*e2_y);
                                momentum[2][ip] = new_momentum_norm*(cos_omega*e1_z + sin_omega*e2_z);
                                weight[ip] = 0.5*total_weight;

                                if (isnan(momentum[0][ip])
                                    || isnan(omega)
                                    || isnan(momentum[1][ip])
                                    || isnan(momentum[2][ip])
                                    || momentum[0][ip] < mx_min + mx_i*mx_delta
                                    || momentum[1][ip] < my_min + my_i*my_delta
                                    || momentum[2][ip] < mz_min + mz_i*mz_delta
                                    || momentum[0][ip] > mx_min + (mx_i+1)*mx_delta
                                    || momentum[1][ip] > my_min + (my_i+1)*my_delta
                                    || momentum[2][ip] > mz_min + (mz_i+1)*mz_delta
                                ) {
                                    //std::cerr <<
                                    ERROR(
                                                 " dim: " << dim[0] << " " << dim[1] << " " << dim[2]
                                              << std::scientific
                                              << " mx: " << momentum[0][ip]
                                              << " my: " << momentum[1][ip]
                                              << " mz: " << momentum[2][ip]
                                              << " new_momentum_norm: " << new_momentum_norm
                                              << " total_weight: " << total_weight
                                              << " total_momentum_norm: " << total_momentum_norm
                                              << " weight[ip]: " << weight[ip]
                                              << " omega: " << omega
                                              << " " << total_momentum_norm / (total_weight*new_momentum_norm)
                                              << " cos_omega: " << cos_omega << " sin_omega" << sin_omega
                                              << " cell_vec: " << cell_vec_x << " " << cell_vec_y << " " << cell_vec_z
                                              << " e1: " << e1_x << " " << e1_y << " " << e1_z
                                              << " e2: " << e2_x << " " << e2_y << " " << e2_z
                                          )
                                    //<< std::endl;
                                
                                }

                                // Update momentum of the second particle
                                ip = sorted_particles[momentum_cell_particle_index[ic] + ipr_min + 1];
                                momentum[0][ip] = new_momentum_norm*(cos_omega*e1_x - sin_omega*e2_x);
                                momentum[1][ip] = new_momentum_norm*(cos_omega*e1_y - sin_omega*e2_y);
                                momentum[2][ip] = new_momentum_norm*(cos_omega*e1_z - sin_omega*e2_z);
                                weight[ip] = 0.5*total_weight;


                                // Other particles are tagged to be removed after
                                for (ipr = ipr_min + 2; ipr < ipr_max ; ipr ++) {
                                    ip = sorted_particles[momentum_cell_particle_index[ic] + ipr];
                                    mask[ip] = -1;
                                    count--;
                                }

                            }

                        }
                    }

                }
            }
        }

        // Free aligned arrays
        free(gamma);
        free(momentum_cell_index);
        free(sorted_particles);
        free(particles_per_momentum_cells);
        free(momentum_cell_particle_index);
    }

}
