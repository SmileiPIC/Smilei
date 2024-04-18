
#include "MF_Solver1D_Yee.h"

#include "ElectroMagn.h"
#include "Field1D.h"

MF_Solver1D_Yee::MF_Solver1D_Yee( Params &params )
    : Solver1D( params )
{
    isEFilterApplied = params.Friedman_filter;
}

MF_Solver1D_Yee::~MF_Solver1D_Yee()
{
}

void MF_Solver1D_Yee::operator()( ElectroMagn *fields )
{
    // const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    
    // Static-cast of the fields
    /*Field1D* Ey1D;
    Field1D* Ez1D;
    if (isEFilterApplied) {
        Ey1D = static_cast<Field1D*>(fields->filter_->Ey_[0]);
        Ez1D = static_cast<Field1D*>(fields->filter_->Ez_[0]);
    } else {
        Ey1D = static_cast<Field1D*>(fields->Ey_);
        Ez1D = static_cast<Field1D*>(fields->Ez_);
    }*/
    const double *const __restrict__ Ey1D = isEFilterApplied ? fields->filter_->Ey_[0]->data() :
                                                               fields->Ey_->data(); // [ix] : dual in y   primal in x,z
    const double *const __restrict__ Ez1D = isEFilterApplied ? fields->filter_->Ez_[0]->data() :
                                                               fields->Ez_->data();// [ix] : dual in z   primal in x,y
    
    //Field1D *By1D   = static_cast<Field1D *>( fields->By_ );
    //Field1D *Bz1D   = static_cast<Field1D *>( fields->Bz_ );
    double *const __restrict__ By1D       = fields->By_->data();// [ix] : dual in x,z primal in y
    double *const __restrict__ Bz1D       = fields->Bz_->data();// [ix] : dual in x,y primal in z
    
    // to be deleted
    /*std::cout<< "printing before in FM solver by and bz for nx_d-1="<<nx_d-1<<std::endl;
    for( unsigned int ix=1 ; ix<nx_d-1 ; ++ix ) {
        std::cout<<By1D[ix] << " "<<Bz1D[ix] <<std::endl;
    }*/
    // ---------------------
    // Solve Maxwell-Faraday
    // ---------------------
    // NB: bx is given in 1d and defined when initializing the fields (here put to 0)
    // Transverse fields  by & bz are defined on the dual grid
#if defined( SMILEI_OPENACC_MODE )                                                                                                    
    const int sizeofEy = fields->Ey_->number_of_points_;
    const int sizeofEz = fields->Ez_->number_of_points_;
    const int sizeofBy = fields->By_->number_of_points_;
    const int sizeofBz = fields->Bz_->number_of_points_;
    #pragma acc parallel present( By1D[0:sizeofBy], Bz1D[0:sizeofBz],Ey1D[0:sizeofEy],Ez1D[0:sizeofEz] )                                               
    #pragma acc loop gang vector             
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for
#endif
#if !defined( SMILEI_ACCELERATOR_MODE )
    #pragma omp simd
#endif
    for( unsigned int ix=1 ; ix<nx_d-1 ; ix++ ) {
        By1D[ix] = By1D[ix] + dt_ov_dx * ( Ez1D[ix] - Ez1D[ix-1] );
        Bz1D[ix] = Bz1D[ix] - dt_ov_dx * ( Ey1D[ix] - Ey1D[ix-1] );
        //( *By1D )( ix )= ( *By1D )( ix ) + dt_ov_dx * ( ( *Ez1D )( ix ) - ( *Ez1D )( ix-1 ) ) ;
        //( *Bz1D )( ix )= ( *Bz1D )( ix ) - dt_ov_dx * ( ( *Ey1D )( ix ) - ( *Ey1D )( ix-1 ) ) ;
    }



    /*{
        fields->By_->copyFromDeviceToHost();
        fields->Bz_->copyFromDeviceToHost();
    }
    std::cout<< "printing after in FM solver by and bz for nx_d-1="<<nx_d-1<<std::endl;
    for( unsigned int ix=1 ; ix<nx_d-1 ; ++ix ) {
        std::cout<<By1D[ix] << " "<<Bz1D[ix] <<std::endl;
    }*/
}
