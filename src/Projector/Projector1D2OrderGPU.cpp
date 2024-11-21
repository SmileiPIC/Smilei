

#if defined( SMILEI_ACCELERATOR_GPU )
#include "Projector1D2OrderGPUKernelCUDAHIP.h"
#include <cmath>
#include "Tools.h"
#endif

#include "Projector1D2OrderGPU.h"

#include "ElectroMagn.h"
#include "Patch.h"
#include "gpu.h"


Projector1D2OrderGPU::Projector1D2OrderGPU( Params &parameters, Patch *a_patch )
    : Projector1D{ parameters, a_patch }
{
    Projector1D::dx_inv_         = 1.0 / parameters.cell_length[0];
    Projector1D::dx_ov_dt_       = parameters.cell_length[0] / parameters.timestep;
    Projector1D::i_domain_begin_ = a_patch->getCellStartingGlobalIndex( 0 );

    not_spectral_  = !parameters.is_pxr;
    dts2_ = parameters.timestep / 2.0;
    dts4_ = dts2_ / 2.0;
#if defined( SMILEI_ACCELERATOR_GPU ) 
    x_dimension_bin_count_ = parameters.getGPUBinCount( 1 );
#else
    ERROR( "Only usable in GPU mode! " );
#endif
}

Projector1D2OrderGPU::~Projector1D2OrderGPU()
{
}
#if defined( SMILEI_ACCELERATOR_GPU )


//! Project global current densities (EMfields->Jx_/Jy_/Jz_)
extern "C" void
currentDepositionKernel1DOnDevice( double *__restrict__ host_Jx,
                         double *__restrict__ host_Jy,
                         double *__restrict__ host_Jz,
                         int Jx_size,
                         int Jy_size,
                         int Jz_size,
                         const double *__restrict__ device_particle_position_x,
                         const double *__restrict__ device_particle_momentum_y,
                         const double *__restrict__ device_particle_momentum_z,
                         const short *__restrict__ device_particle_charge,
                         const double *__restrict__ device_particle_weight,
                         const int *__restrict__ host_bin_index,
                         unsigned int x_dimension_bin_count_,
                         const double *__restrict__ host_invgf_,
                         const int *__restrict__ host_iold_,
                         const double *__restrict__ host_deltaold_,
                         double inv_cell_volume,
                         double dx_inv_,
                         double dx_ov_dt_,
                         int    i_domain_begin_,
                         int    not_spectral_ )
{
    cudahip1d::currentDepositionKernel1D( host_Jx, host_Jy, host_Jz,
                                 Jx_size, Jy_size, Jz_size,
                                 device_particle_position_x, device_particle_momentum_y,
                                 device_particle_momentum_z,
                                 device_particle_charge,
                                 device_particle_weight,
                                 host_bin_index,
                                 x_dimension_bin_count_,
                                 host_invgf_,
                                 host_iold_, host_deltaold_,
                                 inv_cell_volume,
                                 dx_inv_,
                                 dx_ov_dt_,
                                 i_domain_begin_,
                                 not_spectral_ );
}


//! Project global current and charge densities (EMfields->Jx_/Jy_/Jz_/rho_)
//!
extern "C" void
currentAndDensityDepositionKernel1DOnDevice( double *__restrict__ host_Jx,
                                   double *__restrict__ host_Jy,
                                   double *__restrict__ host_Jz,
                                   double *__restrict__ host_rho,
                                   int Jx_size,
                                   int Jy_size,
                                   int Jz_size,
                                   int rho_size,
                                   const double *__restrict__ device_particle_position_x,
                                   const double *__restrict__ device_particle_momentum_y,
                                   const double *__restrict__ device_particle_momentum_z,
                                   const short *__restrict__ device_particle_charge,
                                   const double *__restrict__ device_particle_weight,
                                   const int *__restrict__ host_bin_index,
                                   unsigned int x_dimension_bin_count_,
                                   const double *__restrict__ host_invgf_,
                                   const int *__restrict__ host_iold_,
                                   const double *__restrict__ host_deltaold_,
                                   double inv_cell_volume,
                                   double dx_inv_,
                                   double dx_ov_dt_,
                                   int    i_domain_begin_,
                                   int    not_spectral_ )
{
    cudahip1d::currentAndDensityDepositionKernel1D( host_Jx, host_Jy, host_Jz, host_rho,
                                           Jx_size, Jy_size, Jz_size, rho_size,
                                           device_particle_position_x, device_particle_momentum_y,
                                           device_particle_momentum_z,
                                           device_particle_charge,
                                           device_particle_weight,
                                           host_bin_index,
                                           x_dimension_bin_count_,
                                           host_invgf_,
                                           host_iold_, host_deltaold_,
                                           inv_cell_volume,
                                           dx_inv_,
                                           dx_ov_dt_,
                                           i_domain_begin_,
                                           not_spectral_ );
}
#endif

// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep
// ---------------------------------------------------------------------------------------------------------------------
void Projector1D2OrderGPU::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int bin_shift )
{

    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3
    
    // The variable bin received is  number of bin * cluster width.
    // Declare local variables
    int ip;
    double xjn, xj_m_xip, xj_m_xip2;
    double S1[5];            // arrays used for the Esirkepov projection method
    
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    if( type > 0 ) {
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
                                  
        if( type == 1 ) {
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) {
            charge_weight *= particles.momentum( 1, ipart );
        } else {
            charge_weight *= particles.momentum( 2, ipart );
        }
    }
    
    // Initialize variables
    for( unsigned int i=0; i<5; i++ ) {
        S1[i]=0.;
    }//i
    
    // Locate particle new position on the primal grid
    xjn       = particles.position( 0, ipart ) * dx_inv_;
    ip        = round( xjn + 0.5 * ( type==1 ) );     // index of the central node
    xj_m_xip  = xjn - ( double )ip;                   // normalized distance to the nearest grid point
    xj_m_xip2 = xj_m_xip * xj_m_xip;                  // square of the normalized distance to the nearest grid point
    
    // coefficients 2nd order interpolation on 3 nodes
    //ip_m_ipo = ip-ipo;
    S1[1] = 0.5 * ( xj_m_xip2 - xj_m_xip + 0.25 );
    S1[2] = ( 0.75 - xj_m_xip2 );
    S1[3] = 0.5 * ( xj_m_xip2 + xj_m_xip + 0.25 );
    
    ip -= i_domain_begin_ + 2 + bin_shift;
    
    // 2nd order projection for charge density
    // At the 2nd order, oversize = 2.
    for( unsigned int i=0; i<5; i++ ) {
        rhoj[i + ip ] += charge_weight * S1[i];
    }
    
}


void Projector1D2OrderGPU::currentsAndDensityWrapper( ElectroMagn *EMfields,
                                                      Particles   &particles,
                                                      SmileiMPI   *smpi,
                                                      int,
                                                      int,
                                                      int  ithread,
                                                      bool diag_flag,
                                                      bool is_spectral,
                                                      int  ispec,
                                                      int  icell, 
                                                      int  ipart_ref ) // icell and ipart_ref unused
{
    std::vector<int>    &iold  = smpi->dynamics_iold[ithread];
    std::vector<double> &delta = smpi->dynamics_deltaold[ithread];
    std::vector<double> &invgf = smpi->dynamics_invgf[ithread];

    if( diag_flag ) {

        double *const __restrict__ b_Jx = EMfields->Jx_s[ispec] ? EMfields->Jx_s[ispec]->data() : EMfields->Jx_->data();
        unsigned int Jx_size            = EMfields->Jx_s[ispec] ? EMfields->Jx_s[ispec]->size() : EMfields->Jx_->size();

        double *const __restrict__ b_Jy = EMfields->Jy_s[ispec] ? EMfields->Jy_s[ispec]->data() : EMfields->Jy_->data();
        unsigned int Jy_size            = EMfields->Jy_s[ispec] ? EMfields->Jy_s[ispec]->size() : EMfields->Jy_->size();

        double *const __restrict__ b_Jz = EMfields->Jz_s[ispec] ? EMfields->Jz_s[ispec]->data() : EMfields->Jz_->data();
        unsigned int Jz_size            = EMfields->Jz_s[ispec] ? EMfields->Jz_s[ispec]->size() : EMfields->Jz_->size();

        double *const __restrict__ b_rho = EMfields->rho_s[ispec] ? EMfields->rho_s[ispec]->data() : EMfields->rho_->data();
        unsigned int rho_size            = EMfields->rho_s[ispec] ? EMfields->rho_s[ispec]->size() : EMfields->rho_->size();

        // Does not compute Rho !

#if defined( SMILEI_ACCELERATOR_GPU )

        currentAndDensityDepositionKernel1DOnDevice( b_Jx,b_Jy,b_Jz,b_rho,
                            Jx_size, Jy_size, Jz_size, rho_size,
                            particles.getPtrPosition( 0 ),
                            particles.getPtrMomentum( 1 ),
                            particles.getPtrMomentum( 2 ),
                            particles.getPtrCharge(),
                            particles.getPtrWeight(),
                            particles.last_index.data(),
                            x_dimension_bin_count_,
                            invgf.data(),
                            iold.data(),
                            delta.data(),
                            inv_cell_volume,
                            dx_inv_,
                            dx_ov_dt_,
                            i_domain_begin_,
                            not_spectral_ );

#else
        SMILEI_ASSERT( false );
#endif
    } else {
        if( is_spectral ) {
            ERROR( "Not implemented on GPU" );
        }
        else{

#if defined( SMILEI_ACCELERATOR_GPU )
            currentDepositionKernel1DOnDevice(EMfields->Jx_->data(), EMfields->Jz_->data(), EMfields->Jz_->data(),
                    EMfields->Jx_->size(), EMfields->Jy_->size(), EMfields->Jz_->size(),
                    particles.getPtrPosition( 0 ),
                    particles.getPtrMomentum( 1 ),
                    particles.getPtrMomentum( 2 ),
                    particles.getPtrCharge(),
                    particles.getPtrWeight(),
                    particles.last_index.data(),
                    x_dimension_bin_count_,
                    invgf.data(),
                    iold.data(),
                    delta.data(),
                    inv_cell_volume,
                    dx_inv_,
                    dx_ov_dt_,
                    i_domain_begin_,
                    not_spectral_ );
#else
        SMILEI_ASSERT( false );
#endif
        }
    }
}

void Projector1D2OrderGPU::ionizationCurrents( Field      *Jx,
                                               Field      *Jy,
                                               Field      *Jz,
                                               Particles  &particles,
                                               int         ipart,
                                               LocalFields Jion )
{
    ERROR( "Projector1D2OrderGPU::ionizationCurrents(): Not implemented !" );
}

void Projector1D2OrderGPU::susceptibility( ElectroMagn *EMfields,
                                           Particles   &particles,
                                           double       species_mass,
                                           SmileiMPI   *smpi,
                                           int          istart,
                                           int          iend,
                                           int          ithread,
                                           int          icell,
                                           int          ipart_ref )
{
    ERROR( "Projector1D2OrderGPU::susceptibility(): Not implemented !" );
}
