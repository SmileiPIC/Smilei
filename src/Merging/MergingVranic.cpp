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
MergingVranic::MergingVranic( Params &params,
                              Species *species )
    : Merging( params, species )
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
void MergingVranic::operator()(
    Particles &particles,
    SmileiMPI *smpi,
    int istart,
    int iend,
    int ithread,
    int ipart_ref )
{

    // First of all, we check that there is enought particles per cell
    // to process the merging.
    if( iend - istart > merging_ppc_min_threshold_ ) {
    
        // Minima
        double mx_min;
        float theta_min;
        float phi_min;
        
        // Maxima
        double mx_max;
        float theta_max;
        float phi_max;
        
        // Delta
        double mx_delta;
        float theta_delta;
        float phi_delta;
        
        // Index in each direction
        unsigned int mx_i;
        unsigned int theta_i;
        unsigned int phi_i;
        
        // Local particle index
        int k;
        
        // Momentum shortcut
        double *momentum[3];
        for( int i = 0 ; i<3 ; i++ ) {
            momentum[i] =  &( particles.momentum( i, 0 ) );
        }
        
        // Norm of the momentum
        double momentum_norm;
        
        // Local vector to store the momentum index in the momentum discretization
        std::vector <unsigned int> momentum_i( iend-istart, 0 );
        
        // Local vector to store the momentum angles in the spherical base
        std::vector <float> phi( iend-istart, 0 );
        std::vector <float> theta( iend-istart, 0 );
        
        //std::cerr << iend-istart << std::endl;
        
        // ________________________________________________
        // First step: Computation of the maxima and minima
        
        momentum_norm = sqrt( momentum[0][istart]*momentum[0][istart]
                              + momentum[1][istart]*momentum[1][istart]
                              + momentum[2][istart]*momentum[2][istart] );
                              
        mx_min = momentum[0][istart];
        mx_max = mx_min;
        
        theta[0] = atan2( momentum[1][istart], momentum[0][istart] );
        phi[0]   = asin( momentum[2][istart] / momentum_norm );
        
        theta_min = theta[0];
        phi_min   = phi[0];
        
        theta_max = theta_min;
        phi_max   = phi_min;
        
        for( int ipart=istart+1 ; ipart<iend; ipart++ ) {
        
            // Local array index
            k = ipart - istart;
            
            momentum_norm = sqrt( momentum[0][ipart]*momentum[0][ipart]
                                  + momentum[1][ipart]*momentum[1][ipart]
                                  + momentum[2][ipart]*momentum[2][ipart] );
                                  
            mx_min = fmin( mx_min, momentum[0][ipart] );
            mx_max = fmax( mx_max, momentum[0][ipart] );
            
            phi[k]   = asin( momentum[2][ipart] / momentum_norm );
            theta[k] = atan2( momentum[1][ipart], momentum[0][ipart] );
            
            theta_min = fmin( theta_min, theta[k] );
            theta_max = fmax( theta_max, theta[k] );
            
            phi_min = fmin( phi_min, phi[k] );
            phi_max = fmax( phi_max, phi[k] );
            
        }
        
        // __________________________________________________________
        // Second step : Computation of the discretization and steps
        
        // Dimensions
        dimensions_[0] = 50;
        dimensions_[1] = 10;
        dimensions_[2] = 10;
        
        // Computation of the deltas (discretization steps)
        mx_delta = ( mx_max - mx_min ) / dimensions_[0];
        theta_delta = ( theta_max - theta_min ) / dimensions_[1];
        phi_delta = ( phi_max - phi_min ) / dimensions_[2];
        
        // ___________________________________________________________________
        // Third step: for each particle, momentum indexes are computed in the
        // requested discretization. The number of particles per momentum bin
        // is also determined.
        
        for( int ipart=istart ; ipart<iend; ipart++ ) {
        
            k = ipart - istart;
            
            // 3d indexes in the momentum discretization
            mx_i    = ( momentum[0][ipart] - mx_min )/ mx_delta;
            theta_i = ( theta[k] - theta_min )       / theta_delta;
            phi_i   = ( phi[k] - phi_min )           / phi_delta;
            
            // 1d Index in the momentum discretization
            momentum_i[k] = mx_i*dimensions_[1]*dimensions_[2]
                            + theta_i*dimensions_[2] + phi_i;
                            
        }
        
        // ___________________________________________________________________
        // Fourth step: sort particles in correct bins according to their
        // momentum properties
        
        // ___________________________________________________________________
        // Fifth step: for each bin, merge packet of particles composed of
        // at least `min_packet_size_` and `max_packet_size_`
        
    }
    
}
