
#include "MA_Solver1D_norm.h"

#include "ElectroMagn.h"
#include "Field1D.h"

MA_Solver1D_norm::MA_Solver1D_norm( Params &params )
    : Solver1D( params )
{
}

MA_Solver1D_norm::~MA_Solver1D_norm()
{
}

void MA_Solver1D_norm::operator()( ElectroMagn *fields )
{
    {
    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    /*Field1D *Ex1D = static_cast<Field1D *>( fields->Ex_ );
    Field1D *Ey1D = static_cast<Field1D *>( fields->Ey_ );
    Field1D *Ez1D = static_cast<Field1D *>( fields->Ez_ );
    Field1D *By1D = static_cast<Field1D *>( fields->By_ );
    Field1D *Bz1D = static_cast<Field1D *>( fields->Bz_ );
    Field1D *Jx1D = static_cast<Field1D *>( fields->Jx_ );
    Field1D *Jy1D = static_cast<Field1D *>( fields->Jy_ );
    Field1D *Jz1D = static_cast<Field1D *>( fields->Jz_ );*/

    double *const __restrict__ Ex1D       = fields->Ex_->data(); // [x] : dual in x   primal in y,z
    double *const __restrict__ Ey1D       = fields->Ey_->data(); // [x] : dual in y   primal in x,z
    double *const __restrict__ Ez1D       = fields->Ez_->data(); // [x] : dual in z   primal in x,y
    //const double *const __restrict__ Bx1D = fields->Bx_->data(); // [x] : dual in y,z primal in x
    const double *const __restrict__ By1D = fields->By_->data(); // [x] : dual in x,z primal in y
    const double *const __restrict__ Bz1D = fields->Bz_->data(); // [x] : dual in x,y primal in z
    const double *const __restrict__ Jx1D = fields->Jx_->data(); // [x] : dual in x   primal in y,z
    const double *const __restrict__ Jy1D = fields->Jy_->data(); // [x] : dual in y   primal in x,z
    const double *const __restrict__ Jz1D = fields->Jz_->data(); // [x] : dual in z   primal in x,y 

    {
        fields->Ex_->copyFromDeviceToHost();
        fields->Ey_->copyFromDeviceToHost();
        fields->Ez_->copyFromDeviceToHost();
        fields->Jx_->copyFromDeviceToHost();
        fields->Jy_->copyFromDeviceToHost();
        fields->Jz_->copyFromDeviceToHost();
    }
    std::cout<< "printing before in MA solver ex, ey and ez for nx_d="<<nx_d<< "then jx,jy,jz" <<std::endl;
    for( unsigned int ix=0 ; ix<std::min(nx_d,nx_p) ; ++ix ) {
        std::cout<< std::setprecision (15)<<Ex1D[ix] << " " <<Ey1D[ix] << " "<<Ez1D[ix] << " " << Jx1D[ix] << " " <<Jy1D[ix] << " "<<Jz1D[ix]<<std::endl;
    }
    // --------------------
    // Solve Maxwell-Ampere
    // --------------------
    // Calculate the electrostatic field ex on the dual grid
#if defined( SMILEI_OPENACC_MODE )                                                                                                     
    const int sizeofEx = fields->Ex_->number_of_points_;                                                                               
    const int sizeofEy = fields->Ey_->number_of_points_;                                                                               
    const int sizeofEz = fields->Ez_->number_of_points_;                                                                               
    //const int sizeofBx = fields->Bx_->number_of_points_;                                                                               
    const int sizeofBy = fields->By_->number_of_points_;                                                                               
    const int sizeofBz = fields->Bz_->number_of_points_;                   
    #pragma acc parallel present( Ex1D[0:sizeofEx], Jx1D[0:sizeofEx] )
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for
#endif
#if !defined( SMILEI_ACCELERATOR_MODE )
        #pragma omp simd
#endif
    for( unsigned int ix=0 ; ix<nx_d ; ++ix ) {
        //( *Ex1D )( ix )= ( *Ex1D )( ix ) - dt * ( *Jx1D )( ix ) ;
        Ex1D[ix] -= dt * Jx1D[ix];
    }
    // Transverse fields ey, ez  are defined on the primal grid    #pragma acc parallel present( Ex1D[0:sizeofEx], Jx1D[0:sizeofEx], Bx1D[0:sizeofBz],Ey1D[0:sizeofEx], Jy1D[0:sizeofEx], By1D[0:sizeofBz],Ez1D[0:sizeofEx], Jz1D[0:sizeofEx], Bz1D[0:sizeofBz]  )                             

#if defined( SMILEI_OPENACC_MODE )                    
    #pragma acc parallel present(Ey1D[0:sizeofEy], Jy1D[0:sizeofEy], By1D[0:sizeofBy],Ez1D[0:sizeofEz], Jz1D[0:sizeofEz], Bz1D[0:sizeofBz])
    #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for
#endif
#if !defined( SMILEI_ACCELERATOR_MODE )
        #pragma omp simd
#endif
    for( unsigned int ix=0 ; ix<nx_p ; ++ix ) {
        Ey1D[ix] -= dt_ov_dx * Bz1D[ix+1] - Bz1D[ix] - dt * Jy1D[ix];
        Ez1D[ix] += dt_ov_dx * By1D[ix+1] - By1D[ix] - dt * Jz1D[ix];
        //( *Ey1D )( ix )= ( *Ey1D )( ix ) - dt_ov_dx * ( ( *Bz1D )( ix+1 ) - ( *Bz1D )( ix ) ) - dt * ( *Jy1D )( ix ) ;
        //( *Ez1D )( ix )= ( *Ez1D )( ix ) + dt_ov_dx * ( ( *By1D )( ix+1 ) - ( *By1D )( ix ) ) - dt * ( *Jz1D )( ix ) ;
    }

    {
        fields->Ex_->copyFromDeviceToHost();
        fields->Ey_->copyFromDeviceToHost();
        fields->Ez_->copyFromDeviceToHost();
    }
    }
    // to be deleted
    {
        const unsigned int nx_p = fields->dimPrim[0];
        const unsigned int nx_d = fields->dimDual[0];
        double *const __restrict__ Ex1D       = fields->Ex_->data(); // [x] : dual in x   primal in y,z
        double *const __restrict__ Ey1D       = fields->Ey_->data(); // [x] : dual in y   primal in x,z
        double *const __restrict__ Ez1D       = fields->Ez_->data(); // [x] : dual in z   primal in x,y

        std::cout<< "printing after in MA solver ex, ey and ez for nx_d="<<nx_d<< "nx_p = "<< nx_p<<std::endl;
        for( unsigned int ix=0 ; ix<std::min(nx_d,nx_p) ; ++ix ) {
            std::cout<< std::setprecision (15)<<Ex1D[ix] << " " <<Ey1D[ix] << " "<<Ez1D[ix] << " " <<std::endl;
        }
    }
}

