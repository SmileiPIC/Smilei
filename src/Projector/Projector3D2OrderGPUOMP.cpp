#include "Projector3D2OrderGPUOMP.h"

#include "ElectroMagn.h"
#include "Patch.h"
#include "gpu.h"

Projector3D2OrderGPU::Projector3D2OrderGPU( Params &parameters, Patch *a_patch )
    : Projector3D{ parameters, a_patch }
{
    // Shouldn't Projector2D's state initialization be done in Projector2D's
    // constructor ?
    Projector3D::dx_inv_         = 1.0 / parameters.cell_length[0];
    Projector3D::dy_inv_         = 1.0 / parameters.cell_length[1];
    Projector3D::dz_inv_         = 1.0 / parameters.cell_length[2];
    Projector3D::dx_ov_dt_       = parameters.cell_length[0] / parameters.timestep;
    Projector3D::dy_ov_dt_       = parameters.cell_length[1] / parameters.timestep;
    Projector3D::dz_ov_dt_       = parameters.cell_length[2] / parameters.timestep;
    Projector3D::i_domain_begin_ = a_patch->getCellStartingGlobalIndex( 0 );
    Projector3D::j_domain_begin_ = a_patch->getCellStartingGlobalIndex( 1 );
    Projector3D::k_domain_begin_ = a_patch->getCellStartingGlobalIndex( 2 );
    Projector3D::nprimy          = parameters.n_space[1] + 2 * parameters.oversize[1] + 1;
    Projector3D::nprimz          = parameters.n_space[2] + 2 * parameters.oversize[2] + 1;

    // Due to the initialization order (Projector2D's constructor does not
    // initialize it's member variable) we better initialize
    // Projector2D2OrderGPU's member variable after explicitly initializing
    // Projector2D.
    not_spectral  = !parameters.is_pxr;
    dt   = parameters.timestep;
    dts2 = dt / 2.0;
    dts4 = dts2 / 2.0;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    // When sorting is disabled, these values are invalid (-1) and the HIP
    // implementation can't be used.
    x_dimension_bin_count_ = parameters.getGPUBinCount( 1 );
    y_dimension_bin_count_ = parameters.getGPUBinCount( 2 );
    z_dimension_bin_count_ = parameters.getGPUBinCount( 3 );
#else
    ERROR( "Only usable in GPU mode! " );
#endif
}

Projector3D2OrderGPU::~Projector3D2OrderGPU()
{
    // EMPTY
}

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
extern "C" void
currentDepositionKernel( double *__restrict__ Jx,
                         double *__restrict__ Jy,
                         double *__restrict__ Jz,
                         int Jx_size,
                         int Jy_size,
                         int Jz_size,
                         const double *__restrict__ particle_position_x,
                         const double *__restrict__ particle_position_y,
                         const double *__restrict__ particle_position_z,
                         const short *__restrict__ particle_charge,
                         const double *__restrict__ particle_weight,
                         const int *__restrict__ host_bin_index,
                         unsigned int x_dimension_bin_count,
                         unsigned int y_dimension_bin_count,
                         unsigned int z_dimension_bin_count,
                         const double *__restrict__ invgf_,
                         const int *__restrict__ iold_,
                         const double *__restrict__ deltaold_,
                         double inv_cell_volume,
                         double dx_inv,
                         double dy_inv,
                         double dz_inv,
                         double dx_ov_dt,
                         double dy_ov_dt,
                         double dz_ov_dt,
                         int    i_domain_begin,
                         int    j_domain_begin,
                         int    k_domain_begin,
                         int    nprimy,
                         int    nprimz,
                         int    not_spectral );

extern "C" void
currentAndDensityDepositionKernel( double *__restrict__ Jx,
                         double *__restrict__ Jy,
                         double *__restrict__ Jz,
                         double *__restrict__ rho,
                         int Jx_size,
                         int Jy_size,
                         int Jz_size,
                         int rho_size,
                         const double *__restrict__ particle_position_x,
                         const double *__restrict__ particle_position_y,
                         const double *__restrict__ particle_position_z,
                         const short *__restrict__ particle_charge,
                         const double *__restrict__ particle_weight,
                         const int *__restrict__ host_bin_index,
                         unsigned int x_dimension_bin_count,
                         unsigned int y_dimension_bin_count,
                         unsigned int z_dimension_bin_count,
                         const double *__restrict__ invgf_,
                         const int *__restrict__ iold_,
                         const double *__restrict__ deltaold_,
                         double inv_cell_volume,
                         double dx_inv,
                         double dy_inv,
                         double dz_inv,
                         double dx_ov_dt,
                         double dy_ov_dt,
                         double dz_ov_dt,
                         int    i_domain_begin,
                         int    j_domain_begin,
                         int    k_domain_begin,
                         int    nprimy, 
                         int    nprimz,
                         int    not_spectral );



#endif

namespace { // Unnamed namespace == static == internal linkage == no exported symbols

    /// Project global current densities (EMfields->Jx_/Jy_/Jz_)
    ///
    /* inline */ void
    currents( double *__restrict__ Jx,
              double *__restrict__ Jy,
              double *__restrict__ Jz,
              int          Jx_size,
              int          Jy_size,
              int          Jz_size,
              Particles   &particles,
              unsigned int x_dimension_bin_count,
              unsigned int y_dimension_bin_count,
              unsigned int z_dimension_bin_count,
              const double *__restrict__ invgf_,
              const int *__restrict__ iold_,
              const double *__restrict__ deltaold_,
              double inv_cell_volume,
              double dx_inv,
              double dy_inv,
              double dz_inv,
              double dx_ov_dt,
              double dy_ov_dt,
              double dz_ov_dt,
              int    i_domain_begin,
              int    j_domain_begin,
              int    k_domain_begin,
              int    nprimy, 
              int    nprimz,
              double,
              int not_spectral )
    {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
        currentDepositionKernel( Jx,
                                 Jy,
                                 Jz,
                                 Jx_size,
                                 Jy_size,
                                 Jz_size,
                                 particles.getPtrPosition( 0 ),
                                 particles.getPtrPosition( 1 ),
                                 particles.getPtrPosition( 2 ),
                                 particles.getPtrCharge(),
                                 particles.getPtrWeight(),
                                 particles.last_index.data(),
                                 x_dimension_bin_count,
                                 y_dimension_bin_count,
                                 z_dimension_bin_count,
                                 invgf_,
                                 iold_,
                                 deltaold_,
                                 inv_cell_volume,
                                 dx_inv,
                                 dy_inv,
                                 dz_inv,
                                 dx_ov_dt,
                                 dy_ov_dt,
                                 dz_ov_dt,
                                 i_domain_begin,
                                 j_domain_begin,
                                 k_domain_begin,
                                 nprimy, nprimz,
                                 not_spectral );
#else
        SMILEI_ASSERT( false );
#endif
    }

    /// Like currents(), project the particle current on the grid (Jx_/Jy_/Jz_)
    /// but also compute global current densities rho used for diagFields timestep
    ///
    /* inline */ void
    currentsAndDensity( double *__restrict__ Jx,
                        double *__restrict__ Jy,
                        double *__restrict__ Jz,
                        double *__restrict__ rho,
                        int Jx_size,
                        int Jy_size,
                        int Jz_size,
                        int rho_size,
                        Particles   &particles,
                        unsigned int x_dimension_bin_count,
                        unsigned int y_dimension_bin_count,
                        unsigned int z_dimension_bin_count,
                        const double *__restrict__ invgf_,
                        const int *__restrict__ iold_,
                        const double *__restrict__ deltaold_,
                        double inv_cell_volume,
                        double dx_inv,
                        double dy_inv,
                        double dz_inv,
                        double dx_ov_dt,
                        double dy_ov_dt,
                        double dz_ov_dt,
                        int    i_domain_begin,
                        int    j_domain_begin,
                        int    k_domain_begin,
                        int    nprimy,
                        int    nprimz,
                        double,
                        int not_spectral )
    {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
        currentAndDensityDepositionKernel( Jx,
                                 Jy,
                                 Jz,
                                 rho,
                                 Jx_size,
                                 Jy_size,
                                 Jz_size,
                                 rho_size,
                                 particles.getPtrPosition( 0 ),
                                 particles.getPtrPosition( 1 ),
                                 particles.getPtrPosition( 2 ),
                                 particles.getPtrCharge(),
                                 particles.getPtrWeight(),
                                 particles.last_index.data(),
                                 x_dimension_bin_count,
                                 y_dimension_bin_count,
                                 z_dimension_bin_count,
                                 invgf_,
                                 iold_,
                                 deltaold_,
                                 inv_cell_volume,
                                 dx_inv,
                                 dy_inv,
                                 dz_inv,
                                 dx_ov_dt,
                                 dy_ov_dt,
                                 dz_ov_dt,
                                 i_domain_begin,
                                 j_domain_begin,
                                 k_domain_begin,
                                 nprimy, nprimz,
                                 not_spectral );
#else
        SMILEI_ASSERT( false );
#endif
    }

} // namespace

void Projector3D2OrderGPU::basic( double      *rhoj,
                                  Particles   &particles,
                                  unsigned int ipart,
                                  unsigned int type,
                                  int bin_shift )
{


    // Warning : this function is used for frozen species only. It is assumed that position = position_old !!!

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------

    int iloc, ny( nprimy );
    // (x,y,z) components of the current density for the macro-particle
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );

    if( type > 0 ) {
        charge_weight *= 1./std::sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
                                  
        if( type == 1 ) {
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) {
            charge_weight *= particles.momentum( 1, ipart );
            ny++;
        } else {
            charge_weight *= particles.momentum( 2, ipart );
        }
    }

    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2;
    double Sx1[5], Sy1[5], Sz1[5]; // arrays used for the Esirkepov projection method

    // Initialize all current-related arrays to zero
    for( unsigned int i=0; i<5; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }

    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------

    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip        = std::round( xpn + 0.5 * ( type==1 ) ); // index of the central node
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sx1[1] = 0.5 * ( delta2-delta+0.25 );
    Sx1[2] = 0.75-delta2;
    Sx1[3] = 0.5 * ( delta2+delta+0.25 );

    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = std::round( ypn + 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sy1[1] = 0.5 * ( delta2-delta+0.25 );
    Sy1[2] = 0.75-delta2;
    Sy1[3] = 0.5 * ( delta2+delta+0.25 );

    zpn = particles.position( 1, ipart ) * dz_inv_;
    int kp = std::round( zpn + 0.5*( type==2 ) );
    delta  = zpn - ( double )kp;
    delta2 = delta*delta;
    Sz1[1] = 0.5 * ( delta2-delta+0.25 );
    Sz1[2] = 0.75-delta2;
    Sz1[3] = 0.5 * ( delta2+delta+0.25 );

    // ---------------------------
    // Calculate the total current
    // ---------------------------
    ip -= i_domain_begin_ + 2;
    jp -= j_domain_begin_ + 2;
    kp -= k_domain_begin_ + 2;

    for( unsigned int i=0 ; i<5 ; i++ ) {
        for( unsigned int j=0 ; j<5 ; j++ ) {
            iloc = (( i+ip )*ny+jp+j)*nz + kp;
            for( unsigned int k=0 ; k<5 ; k++ ) {
                rhoj[iloc+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[j];
            }
        }
    }
}

void Projector3D2OrderGPU::ionizationCurrents( Field      *Jx,
                                               Field      *Jy,
                                               Field      *Jz,
                                               Particles  &particles,
                                               int         ipart,
                                               LocalFields Jion )
{
    ERROR( "Projector3D2OrderGPU::ionizationCurrents(): Not implemented !" );
}

void Projector3D2OrderGPU::currentsAndDensityWrapper( ElectroMagn *EMfields,
                                                      Particles   &particles,
                                                      SmileiMPI   *smpi,
                                                      int,
                                                      int,
                                                      int  ithread,
                                                      bool diag_flag,
                                                      bool is_spectral,
                                                      int  ispec,
                                                      int  icell,
                                                      int  ipart_ref )
{
    std::vector<int>    &iold  = smpi->dynamics_iold[ithread];
    std::vector<double> &delta = smpi->dynamics_deltaold[ithread];
    std::vector<double> &invgf = smpi->dynamics_invgf[ithread];
    Jx_                        = EMfields->Jx_->data();
    Jy_                        = EMfields->Jy_->data();
    Jz_                        = EMfields->Jz_->data();
    rho_                       = EMfields->rho_->data();

    if( diag_flag ) {
        // TODO(Etienne M): DIAGS. Find a way to get rho. We could:
        // - Pull everything we need from the GPU and compute on the host
        // - Implement currentsAndDensity on GPU which means:
        //      - The J<x>_s/Rho_s, if required by the diags must be on the GPU
        //          -
        //      - Rho_ must be on the GPU if species specific charge/current densities are not diagnosed
        //

        // double *const b_Jx  = EMfields->Jx_s[ispec] ? &( *EMfields->Jx_s[ispec] )( 0 ) : Jx_;
        // double *const b_Jy  = EMfields->Jy_s[ispec] ? &( *EMfields->Jy_s[ispec] )( 0 ) : Jy_;
        // double *const b_Jz  = EMfields->Jz_s[ispec] ? &( *EMfields->Jz_s[ispec] )( 0 ) : Jz_;
        // double *const b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : rho_;

        double *const __restrict__ b_Jx  = EMfields->Jx_s[ispec] ? EMfields->Jx_s[ispec]->data() : EMfields->Jx_->data();
        unsigned int Jx_size             = EMfields->Jx_s[ispec] ? EMfields->Jx_s[ispec]->size() : EMfields->Jx_->size();

        double *const __restrict__ b_Jy  = EMfields->Jy_s[ispec] ? EMfields->Jy_s[ispec]->data() : EMfields->Jy_->data();
        unsigned int Jy_size             = EMfields->Jy_s[ispec] ? EMfields->Jy_s[ispec]->size() : EMfields->Jy_->size();

        double *const __restrict__ b_Jz  = EMfields->Jz_s[ispec] ? EMfields->Jz_s[ispec]->data() : EMfields->Jz_->data();
        unsigned int Jz_size             = EMfields->Jz_s[ispec] ? EMfields->Jz_s[ispec]->size() : EMfields->Jz_->size();

        double *const __restrict__ b_rho  = EMfields->rho_s[ispec] ? EMfields->rho_s[ispec]->data() : EMfields->rho_->data();
        unsigned int rho_size             = EMfields->rho_s[ispec] ? EMfields->rho_s[ispec]->size() : EMfields->rho_->size();

        // for( int ipart = istart; ipart < iend; ipart++ ) {
        //     currentsAndDensity( b_Jx, b_Jy, b_Jz, b_rho,
        //                         particles, ipart,
        //                         invgf[ipart], &iold[ipart], &delta[ipart],
        //                         inv_cell_volume,
        //                         dx_inv_, dy_inv_,
        //                         dx_ov_dt_, dy_ov_dt_,
        //                         i_domain_begin_, j_domain_begin_,
        //                         nprimy,
        //                         one_third,
        //                         pxr );
        // }

        // Does not compute Rho !

        currentsAndDensity( b_Jx, b_Jy, b_Jz, b_rho,
                  Jx_size, Jy_size, Jz_size, rho_size,
                  particles, x_dimension_bin_count_, y_dimension_bin_count_, z_dimension_bin_count_,
                  invgf.data(), iold.data(), delta.data(),
                  inv_cell_volume,
                  dx_inv_, dy_inv_, dz_inv_,
                  dx_ov_dt_, dy_ov_dt_, dz_ov_dt_,
                  i_domain_begin_, j_domain_begin_, k_domain_begin_,
                  nprimy, nprimz,
                  one_third,
                  not_spectral );

    } else {
        // If no field diagnostics this timestep, then the projection is done directly on the total arrays
        if( is_spectral ) {
            ERROR( "Not implemented on GPU" );
            // for( int ipart = istart; ipart < iend; ipart++ ) {
            //     currentsAndDensity( Jx_, Jy_, Jz_, rho_,
            //                         particles, ipart,
            //                         invgf[ipart], &iold[ipart], &delta[ipart],
            //                         inv_cell_volume,
            //                         dx_inv_, dy_inv_,
            //                         dx_ov_dt_, dy_ov_dt_,
            //                         i_domain_begin_, j_domain_begin_,
            //                         nprimy,
            //                         one_third,
            //                         pxr );
            // }
        } else {
            currents( Jx_, Jy_, Jz_,
                      EMfields->Jx_->globalDims_, EMfields->Jy_->globalDims_, EMfields->Jz_->globalDims_,
                      particles, x_dimension_bin_count_, y_dimension_bin_count_, z_dimension_bin_count_,
                      invgf.data(), iold.data(), delta.data(),
                      inv_cell_volume,
                      dx_inv_, dy_inv_, dz_inv_,
                      dx_ov_dt_, dy_ov_dt_, dz_ov_dt_,
                      i_domain_begin_, j_domain_begin_, k_domain_begin_,
                      nprimy, nprimz,
                      one_third,
                      not_spectral );
        }
    }
}

void Projector3D2OrderGPU::susceptibility( ElectroMagn *EMfields,
                                           Particles   &particles,
                                           double       species_mass,
                                           SmileiMPI   *smpi,
                                           int          istart,
                                           int          iend,
                                           int          ithread,
                                           int          icell,
                                           int          ipart_ref )
{
    ERROR( "Projector3D2OrderGPU::susceptibility(): Not implemented !" );
}
