
#include "MA_Solver3D_Friedman.h"

#include "ElectroMagn.h"
#include "Field3D.h"

MA_Solver3D_Friedman::MA_Solver3D_Friedman( Params &params )
    : Solver3D( params )
{
    ftheta = params.Friedman_theta;
    alpha  = 1.-0.5*ftheta+0.5*ftheta*ftheta;
    beta   = ftheta*( 1.-0.5*ftheta );
    delta  = 0.5*ftheta*( 1.-ftheta )*( 1.-ftheta );
}

MA_Solver3D_Friedman::~MA_Solver3D_Friedman()
{
}

void MA_Solver3D_Friedman::operator()( ElectroMagn *fields )
{

    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    const unsigned int nz_p = fields->dimPrim[2];
    const unsigned int nz_d = fields->dimDual[2];
    // Static-cast of the fields
    double *Ex3D  = &(fields->Ex_->data_[0]);
    double *Ey3D  = &(fields->Ey_->data_[0]);
    double *Ez3D  = &(fields->Ez_->data_[0]);
    double *Bx3D  = &(fields->Bx_->data_[0]);
    double *By3D  = &(fields->By_->data_[0]);
    double *Bz3D  = &(fields->Bz_->data_[0]);
    double *Jx3D  = &(fields->Jx_->data_[0]);
    double *Jy3D  = &(fields->Jy_->data_[0]);
    double *Jz3D  = &(fields->Jz_->data_[0]);

    double *Ex_f  = &( fields->filter_->Ex_[0]->data_[0] );
    double *Ey_f  = &( fields->filter_->Ex_[1]->data_[0] );
    double *Ez_f  = &( fields->filter_->Ex_[2]->data_[0] );
    double *Ex_m1 = &( fields->filter_->Ex_[1]->data_[0] );
    double *Ey_m1 = &( fields->filter_->Ey_[1]->data_[0] );
    double *Ez_m1 = &( fields->filter_->Ez_[1]->data_[0] );
    double *Ex_m2 = &( fields->filter_->Ex_[2]->data_[0] );
    double *Ey_m2 = &( fields->filter_->Ey_[2]->data_[0] );
    double *Ez_m2 = &( fields->filter_->Ez_[2]->data_[0] );

    double adv = 0.;
    int index_i_j_k;
    // Electric field Ex^(d,p,p)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {

              adv             = -dt*Jx3D[ i*(ny_p*nz_p) + j*(nz_p) + k ]
                  +                 dt_ov_dy * ( Bz3D[ i*(ny_d*nz_p) + (j+1)*(nz_p) + k   ] - Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] )
                  -                 dt_ov_dz * ( By3D[ i*(ny_p*nz_d) +  j   *(nz_d) + k+1 ] - By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] );

              index_i_j_k = i*(ny_p*nz_p) + j*(nz_p) + k;
              // advance electric field
              Ex3D [ index_i_j_k ]   += adv;
              // compute the time-filtered field
              Ex_f [ index_i_j_k ]    = alpha*(Ex3D[ index_i_j_k ])
                                      + beta*adv
                                      + delta*( Ex_m1[ index_i_j_k ])
                                      + ftheta*( Ex_m2 [ index_i_j_k ] );
              // update Ex_m2 and Ex_m1
              Ex_m2[ index_i_j_k ]    = Ex_m1[ index_i_j_k ]
                                      - ftheta* (Ex_m2[ index_i_j_k ]);
              Ex_m1[ index_i_j_k ]    = Ex3D [ index_i_j_k ]  - adv;

            }
        }
    }

    // Electric field Ey^(p,d,p)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {

                adv             = -dt*Jy3D[ i*(ny_d*nz_p) + j*(nz_p) + k ]
                -                  dt_ov_dx * ( Bz3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + k   ] - Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] )
                +                  dt_ov_dz * ( Bx3D[  i   *(ny_d*nz_d) + j*(nz_d) + k+1 ] - Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] );

                index_i_j_k = i*(ny_d*nz_p) + j*(nz_p) + k;
                // advance electric field
                Ey3D [ index_i_j_k ]   += adv;
                // compute the time-filtered field
                Ey_f [ index_i_j_k ]    = alpha*(Ey3D[ index_i_j_k ])
                                        + beta*adv
                                        + delta*( Ey_m1[ index_i_j_k ])
                                        + ftheta*( Ey_m2 [ index_i_j_k ] );
                // update Ex_m2 and Ex_m1
                Ey_m2[ index_i_j_k ]    = Ey_m1[ index_i_j_k ]
                                        - ftheta* (Ey_m2[ index_i_j_k ]);
                Ey_m1[ index_i_j_k ]    = Ey3D [ index_i_j_k ]  - adv;

            }
        }
    }

    // Electric field Ez^(p,p,d)
    for( unsigned int i=0 ;  i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=0 ; k<nz_d ; k++ ) {

                adv             = -dt*Jz3D[ i*(ny_p*nz_d) + j*(nz_d) + k ]
                +                  dt_ov_dx * ( By3D[ (i+1)*(ny_p*nz_d) +  j   *(nz_d) + k ] - By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] )
                -                  dt_ov_dy * ( Bx3D[  i   *(ny_d*nz_d) + (j+1)*(nz_d) + k ] - Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] );

                index_i_j_k = i*(ny_p*nz_d) + j*(nz_d) + k;
                // advance electric field
                Ez3D [ index_i_j_k ]   += adv;
                // compute the time-filtered field
                Ez_f [ index_i_j_k ]    = alpha*(Ez3D[ index_i_j_k ])
                                        + beta*adv
                                        + delta*( Ez_m1[ index_i_j_k ])
                                        + ftheta*( Ez_m2 [ index_i_j_k ] );
                // update Ex_m2 and Ex_m1
                Ez_m2[ index_i_j_k ]    = Ez_m1[ index_i_j_k ]
                                        - ftheta* (Ez_m2[ index_i_j_k ]);
                Ez_m1[ index_i_j_k ]    = Ez3D [ index_i_j_k ]  - adv;

            }
        }
    }

}
