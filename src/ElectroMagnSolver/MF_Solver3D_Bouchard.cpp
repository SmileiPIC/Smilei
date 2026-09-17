#include "MF_Solver3D_Bouchard.h"

#include "ElectroMagn.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"

#include <algorithm>

MF_Solver3D_Bouchard::MF_Solver3D_Bouchard( Params &params )
    : Solver3D( params )
{
    //ERROR("Under development, not yet working");
    double dt = params.timestep;
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    double dx_ov_dt  = dx/dt;
    double dy_ov_dt  = dy/dt;
    double dz_ov_dt  = dz/dt;
    double dt_ov_dx  = dt/dx;
    double dt_ov_dy  = dt/dy;
    double dt_ov_dz  = dt/dz;
    //Not necessary to have dx=dy=dz, but dispersion law are modify
    //In particular if dz >> dx,dy then solver become like the 2d solver
    if( (dx!=dy)||(dx!=dz)||(dy!=dz) ) {
        WARNING( "Bouchard solver works best with identical cell-lengths in all directions" );
    }
    if( dx_ov_dt!=2 ) {
        WARNING( "Bouchard solver requires dx/dt = 2 (Magic Timestep)" );
    }

    double delta = 0.1222*(1-2*2)/4. ; // 0.1222*(1-pow(2.,2))/4. ;
    double beta = -0.1727*(1-0.5*2*2-4.*delta)/4. ; // -0.1727*(1-0.5*pow(2.,2)-4.*delta)/4. ;
    double alpha = 1-4.*beta-3.*delta ;

    delta_x = delta ;
    delta_y = delta ;
    delta_z = delta ;
    beta_xy = beta ;
    beta_yx = beta ;
    beta_xz = beta ;
    beta_zx = beta ;
    beta_yz = beta ;
    beta_zy = beta ;
    alpha_x = alpha ;
    alpha_y = alpha ;
    alpha_z = alpha ;

    Ax  = alpha_x*dt/dx;
    Ay  = alpha_y*dt/dy;
    Az  = alpha_z*dt/dz;
    Bxy = beta_xy*dt/dx;
    Byx = beta_yx*dt/dy;
    Bxz = beta_xz*dt/dx;
    Bzx = beta_zx*dt/dz;
    Byz = beta_yz*dt/dy;
    Bzy = beta_zy*dt/dz;
    Dx  = delta_x*dt/dx;
    Dy  = delta_y*dt/dy;
    Dz  = delta_z*dt/dz;

    isEFilterApplied = params.Friedman_filter;
}

MF_Solver3D_Bouchard::~MF_Solver3D_Bouchard()
{
}

void MF_Solver3D_Bouchard::operator()( ElectroMagn* fields )
{
    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    const unsigned int nz_p = fields->dimPrim[2];
    const unsigned int nz_d = fields->dimDual[2];

    // Static-cast of the fields

    /** Legacy
    Field3D *Ex3D = static_cast<Field3D *>( fields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( fields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( fields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( fields->Bx_ );
    Field3D *By3D = static_cast<Field3D *>( fields->By_ );
    Field3D *Bz3D = static_cast<Field3D *>( fields->Bz_ );
    **/

    const double *const __restrict__ Ex3D = isEFilterApplied ? fields->filter_->Ex_[0]->data() : fields->Ex_->data(); // [x * (ny_p*nz_p) + y*nz_p + z] : dual in x primal in y,z (dpp)
    const double *const __restrict__ Ey3D = isEFilterApplied ? fields->filter_->Ey_[0]->data() : fields->Ey_->data(); // [x * (ny_d*nz_p) + y*nz_p + z] : dual in y primal in x,z (pdp)
    const double *const __restrict__ Ez3D = isEFilterApplied ? fields->filter_->Ez_[0]->data() : fields->Ez_->data(); // [x * (ny_p*nz_d) + y*nz_d + z] : dual in z primal in x,y (ppd)
    double *const __restrict__ Bx3D       = fields->Bx_->data();                    // [x * (ny_d*nz_d) + y*nz_d + z] : dual in y,z primal in x
    double *const __restrict__ By3D       = fields->By_->data();                    // [x * (ny_p*nz_d) + y*nz_d + z] : dual in x,z primal in y
    double *const __restrict__ Bz3D       = fields->Bz_->data();                    // [x * (ny_d*nz_p) + y*nz_p + z] : dual in x,y primal in z


    // Magnetic field Bx^(p,d,d) Electric field Ey^(p,d,p) Ez^(p,p,d) 
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    const int sizeofEx = fields->Ex_->number_of_points_;
    const int sizeofEy = fields->Ey_->number_of_points_;
    const int sizeofEz = fields->Ez_->number_of_points_;
    const int sizeofBx = fields->Bx_->number_of_points_;
    const int sizeofBy = fields->By_->number_of_points_;
    const int sizeofBz = fields->Bz_->number_of_points_;

    #pragma acc parallel present( Bx3D[0:sizeofBx], Ey3D[0:sizeofEy], Ez3D[0:sizeofEz] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=1 ; i<nx_p-1;  i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc loop worker
#endif
        for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {

                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += Az * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j*(nz_p) + k-1 ] )
                                                     + Bzx * ( Ey3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i+1)*(ny_d*nz_p) +  j*(nz_p) + k-1 ] + Ey3D[ (i-1)*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j*(nz_p) + k-1 ] )
                                                     + Bzy * ( Ey3D[ i*(ny_d*nz_p) + (j+1)*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  (j+1)*(nz_p) + k-1 ] + Ey3D[ i*(ny_d*nz_p) + (j-1)*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  (j-1)*(nz_p) + k-1 ] )
                                                     +  Dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k+1 ] - Ey3D[ i*(ny_d*nz_p) +  j*(nz_p) + k-2 ] )
                                                     -  Ay * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) +  (j-1)*(nz_d) + k ] )
                                                     - Byx * ( Ez3D[ (i+1)*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i+1)*(ny_p*nz_d) +  (j-1)*(nz_d) + k ] + Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) +  (j-1)*(nz_d) + k ] )
                                                     - Byz * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + (k+1) ] - Ez3D[ i*(ny_p*nz_d) +  (j-1)*(nz_d) + k+1 ] + Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k-1 ] - Ez3D[ i*(ny_p*nz_d) +  (j-1)*(nz_d) + k-1 ] )
                                                     -  Dy * ( Ez3D[ i*(ny_p*nz_d) + (j+1)*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) +  (j-2)*(nz_d) + k ] );

                /** Legacy
                ( *Bx3D )( i, j, k ) += Az * ( ( *Ey3D )( i, j, k )-( *Ey3D )( i, j, k-1 ) )
                                     + Bzx * ( ( *Ey3D )( i+1, j, k ) - ( *Ey3D )( i+1, j, k-1 ) + ( *Ey3D )( i-1, j, k )-( *Ey3D )( i-1, j, k-1 ) )
                                     + Bzy * ( ( *Ey3D )( i, j+1, k ) - ( *Ey3D )( i, j+1, k-1 ) + ( *Ey3D )( i, j-1, k )-( *Ey3D )( i, j-1, k-1 ) )
                                     +  Dz * ( ( *Ey3D )( i, j, k+1 ) - ( *Ey3D )( i, j, k-2) )
                                     -  Ay * ( ( *Ez3D )( i,  j, k )  - ( *Ez3D )( i,  j-1, k ) )
                                     - Byx * ( ( *Ez3D )( i+1, j, k )  - ( *Ez3D )( i+1, j-1, k ) + ( *Ez3D )( i-1, j, k )-( *Ez3D )( i-1, j-1, k ) )
                                     - Byz * ( ( *Ez3D )( i, j, k+1 )  - ( *Ez3D )( i, j-1, k+1 ) + ( *Ez3D )( i, j, k-1 )-( *Ez3D )( i, j-1, k-1 ) )
                                     -  Dy * ( ( *Ez3D )( i, j+1, k )-( *Ez3D )( i, j-2, k ) );
                **/
            }
        }
    }

    // Magnetic field By^(d,p,d) Electric field Ex^(d,p,p) Ez^(p,p,d)
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( By3D[0:sizeofBy], Ex3D[0:sizeofEx], Ez3D[0:sizeofEz] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc loop worker
#endif
        for( unsigned int j=1 ; j<ny_p-1 ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {

                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += Ax * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k ] )
                                                     + Bxy * ( Ez3D[ i*(ny_p*nz_d) + (j+1)*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + (j+1)*(nz_d) + k ] + Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + (j-1)*(nz_d) + k ] )
                                                     + Bxz * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k+1 ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k+1 ] + Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k-1 ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k-1 ] )
                                                     +  Dx * ( Ez3D[ (i+1)*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-2)*(ny_p*nz_d) + j*(nz_d) + k ] )
                                                     -  Az * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     - Bzx * ( Ex3D[ (i+1)*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[ (i+1)*(ny_p*nz_p) + j*(nz_p) + k-1 ] + Ex3D[ (i-1)*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[ (i-1)*(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     - Bzy * ( Ex3D[ i*(ny_p*nz_p) + (j+1)*(nz_p) + k ] - Ex3D[ i*(ny_p*nz_p) + (j+1)*(nz_p) + k-1 ] + Ex3D[ i*(ny_p*nz_p) + (j-1)*(nz_p) + k ] - Ex3D[ i*(ny_p*nz_p) + (j-1)*(nz_p) + k-1 ] )
                                                     -  Dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k+1 ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-2 ] ) ;
                
                /** Legacy
                ( *By3D )( i, j, k ) += Ax * ( ( *Ez3D )( i,  j, k ) - ( *Ez3D )( i-1, j, k ) )
                                     + Bxy * ( ( *Ez3D )( i,  j+1, k ) - ( *Ez3D )( i-1, j+1, k ) + ( *Ez3D )( i, j-1, k )-( *Ez3D )( i-1, j-1, k ) )
                                     + Bxz * ( ( *Ez3D )( i, j, k+1 ) - ( *Ez3D )( i-1, j, k+1 ) + ( *Ez3D )( i, j, k-1 )-( *Ez3D )( i-1, j, k-1 ) )
                                     +  Dx * ( ( *Ez3D )( i+1, j, k ) - ( *Ez3D )( i-2, j, k ) )
                                     -  Az * ( ( *Ex3D )( i, j, k )-( *Ex3D )( i, j, k-1 ) )
                                     - Bzy * ( ( *Ex3D )( i, j+1, k )-( *Ex3D )( i, j+1, k-1 ) + ( *Ex3D )( i, j-1, k )-( *Ex3D )( i, j-1, k-1 ) )
                                     - Bzx * ( ( *Ex3D )( i+1, j, k )-( *Ex3D )( i+1, j, k-1 ) + ( *Ex3D )( i-1, j, k )-( *Ex3D )( i-1, j, k-1 ) )
                                     -  Dz * ( ( *Ex3D )( i, j, k+1 ) - ( *Ex3D )( i, j, k-2) ) ;
                **/
            }
        }
    }

    // Magnetic field Bz^(d,d,p) Electric field Ex^(d,p,p) Ey^(p,d,p)
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( Bz3D[0:sizeofBz], Ex3D[0:sizeofEx], Ey3D[0:sizeofEy] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc loop worker
#endif
        for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( unsigned int k=1 ; k<nz_p-1 ; k++ ) {

                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += Ay * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] )
                                                     + Byx * ( Ex3D[ (i+1)*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  (i+1)   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] + Ex3D[ (i-1)*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  (i-1)   *(ny_p*nz_p) + (j-1)*(nz_p) + k ])
                                                     + Byz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k+1 ] - Ex3D[  i*(ny_p*nz_p) + (j-1)*(nz_p) + k+1 ] + Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k-1 ] - Ex3D[  i*(ny_p*nz_p) + (j-1)*(nz_p) + k-1 ] )
                                                     + Dy  * ( Ex3D[ i*(ny_p*nz_p) + (j+1)*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-2)*(nz_p) + k ] )
                                                     -  Ax * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     - Bxy * ( Ey3D[ i*(ny_d*nz_p) + (j+1)*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  (j+1)   *(nz_p) + k ] + Ey3D[ i*(ny_d*nz_p) + (j-1)*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  (j-1)   *(nz_p) + k ] )
                                                     - Bxz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k+1 ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k+1 ] + Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k-1 ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] )
                                                     - Dx  * ( Ey3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-2)*(ny_d*nz_p) +  j   *(nz_p) + k ] ) ;

                /** Legacy
                ( *Bz3D )( i, j, k ) += Ay * ( ( *Ex3D )( i, j, k )-( *Ex3D )( i, j-1, k ) )
                                     + Byx * ( ( *Ex3D )( i+1, j, k )-( *Ex3D )( i+1, j-1, k ) + ( *Ex3D )( i-1, j, k )-( *Ex3D )( i-1, j-1, k ))
                                     + Byz * ( ( *Ex3D )( i, j, k+1 )-( *Ex3D )( i, j-1, k+1 ) + ( *Ex3D )( i, j, k-1 )-( *Ex3D )( i, j-1, k-1 ))
                                     + Dy  * ( ( *Ex3D )( i, j+1, k )-( *Ex3D )( i, j-2, k ) )
                                     -  Ax * ( ( *Ey3D )( i, j, k )-( *Ey3D )( i-1, j, k ) )
                                     - Bxz * ( ( *Ey3D )( i, j, k+1 )-( *Ey3D )( i-1, j, k+1 ) + ( *Ey3D )( i, j, k-1 )-( *Ey3D )( i-1, j, k-1 ))
                                     - Bxy * ( ( *Ey3D )( i, j+1, k )-( *Ey3D )( i-1, j+1, k ) + ( *Ey3D )( i, j-1, k )-( *Ey3D )( i-1, j-1, k ))
                                     - Dx  * ( ( *Ey3D )( i+1, j, k )-( *Ey3D )( i-2, j, k ) ) ;
                **/
            }
        }
    }

   //Additional boundaries treatment on the primal direction of each B field

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( By3D[0:sizeofBy], Ex3D[0:sizeofEx], Ez3D[0:sizeofEz] ) // modify the field here
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 2 )
#endif
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc loop vector
#endif
            for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
                // Magnetic field By^(d,p,d) Electric field Ex^(d,p,p) Ez^(p,p,d)
                unsigned int i=1 ;
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt_ov_dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     +   dt_ov_dx * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k   ] );
                // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
                // Magnetic field By^(d,p,d) Electric field Ex^(d,p,p) Ez^(p,p,d)
                i=nx_d-2 ;
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt_ov_dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     +   dt_ov_dx * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k   ] );
            }
        }

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( Bz3D[0:sizeofBz], Ey3D[0:sizeofEy], Ex3D[0:sizeofEx] ) // modify the field here
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 2 )
#endif
        for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc loop vector
#endif
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
                // Magnetic field Bz^(d,d,p) Electric field Ex^(d,p,p) Ey^(p,d,p)
                unsigned int i=1 ;
                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt_ov_dx * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     +   dt_ov_dy * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] );
                // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
                // Magnetic field Bz^(d,d,p) Electric field Ex^(d,p,p) Ey^(p,d,p)
                i=nx_d-2 ;
                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt_ov_dx * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     +   dt_ov_dy * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] );
            }
        }

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( Bx3D[0:sizeofBx], Ez3D[0:sizeofEz], Ey3D[0:sizeofEy] ) // modify the field here
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 2 )
#endif
        for( unsigned int i=0 ; i<nx_p;  i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc loop vector
#endif
            for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                //At Ymin
                //Additional boundaries treatment for j=1 and j=nx_d-2 for Bx and Bz
                // Magnetic field Bx^(p,d,d) Electric field Ey^(p,d,p) Ez^(p,p,d)
                unsigned int j=1 ;
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += -dt_ov_dy * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k   ] )
                                                     +   dt_ov_dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] );
                //At Ymax
                //Additional boundaries treatment for j=1 and j=nx_d-2 for Bx and Bz
                // Magnetic field Bx^(p,d,d) Electric field Ey^(p,d,p) Ez^(p,p,d)
                j=ny_d-2 ;
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += -dt_ov_dy * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k   ] )
                                                     +   dt_ov_dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] );
            }
        }

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( Bz3D[0:sizeofBz], Ey3D[0:sizeofEy], Ex3D[0:sizeofEx] ) // modify the field here
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 2 )
#endif
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc loop vector
#endif
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                unsigned int j=1 ;
                //At Ymin
                //Additional boundaries treatment for j=1 and j=nx_d-2 for Bx and Bz
                // Magnetic field Bz^(d,d,p) Electric field Ex^(d,p,p) Ey^(p,d,p)
                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt_ov_dx * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     +   dt_ov_dy * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] );
                //At Ymax
                //Additional boundaries treatment for j=1 and j=nx_d-2 for Bx and Bz
                // Magnetic field Bz^(d,d,p) Electric field Ex^(d,p,p) Ey^(p,d,p)
                j=ny_d-2 ;
                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt_ov_dx * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     +   dt_ov_dy * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] );
            }
        }

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( Bx3D[0:sizeofBx], Ez3D[0:sizeofEz], Ey3D[0:sizeofEy] ) // modify the field here
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 2 )
#endif
        for( unsigned int i=0 ; i<nx_p;  i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc loop vector
#endif
            for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
                //At Zmin
                //Additional boundaries treatment for k=1 and k=nx_d-2 for Bx and By
                // Magnetic field Bx^(p,d,d) Electric field Ey^(p,d,p) Ez^(p,p,d)
                unsigned int k=1 ;
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += -dt_ov_dy * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k   ] )
                                                     +   dt_ov_dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] );
                //At Zmax
                //Additional boundaries treatment for k=1 and k=nx_d-2 for Bx and By
                // Magnetic field Bx^(p,d,d) Electric field Ey^(p,d,p) Ez^(p,p,d)
                k=nz_d-2 ;
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += -dt_ov_dy * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k   ] )
                                                     +   dt_ov_dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] );
            }
        }

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( By3D[0:sizeofBy], Ex3D[0:sizeofEx], Ez3D[0:sizeofEz] ) // modify the field here
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 2 )
#endif
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc loop vector
#endif
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                //At Zmin
                //Additional boundaries treatment for k=1 and k=nx_d-2 for Bx and By
                // Magnetic field By^(d,p,d) Electric field Ex^(d,p,p) Ez^(p,p,d)
                unsigned int k=1 ;
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt_ov_dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     +   dt_ov_dx * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k   ] );
                //At Zmax
                //Additional boundaries treatment for k=1 and k=nx_d-2 for Bx and By
                // Magnetic field By^(d,p,d) Electric field Ex^(d,p,p) Ez^(p,p,d)                                                  
                k=nz_d-2 ;
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt_ov_dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     +   dt_ov_dx * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k   ] );
            }
        }

}//END solveMaxwellFaraday
