
#include "SyncVectorPatch.h"

#include <vector>
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #include <openacc.h>
#endif
#include "Params.h"
#include "SmileiMPI.h"
#include "VectorPatch.h"
#include "gpu.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------       PARTICLES         ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------


template void SyncVectorPatch::exchangeAlongAllDirections<double,Field>( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
template void SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
template void SyncVectorPatch::exchangeAlongAllDirectionsNoOMP<double,Field>( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
template void SyncVectorPatch::exchangeAlongAllDirectionsNoOMP<complex<double>,cField>( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );

void SyncVectorPatch::initExchParticles( VectorPatch &vecPatches, int ispec, Params &params, SmileiMPI *smpi )
{
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->copyExchParticlesToBuffers( ispec, params );
    }
    
    // Start exchange along dimension 0 only
    SyncVectorPatch::initExchParticlesAlongDimension( vecPatches, ispec, 0, params, smpi );
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function performs:
//! - the exchange of particles for each direction using the diagonal trick.
//! - the importation of the new particles in the particle property arrays
//! - the sorting of particles
// ---------------------------------------------------------------------------------------------------------------------
void SyncVectorPatch::finalizeExchParticlesAndSort( VectorPatch &vecPatches, int ispec, Params &params, SmileiMPI *smpi )
{
    // finish exchange along dimension 0 only
    SyncVectorPatch::finalizeExchParticlesAlongDimension( vecPatches, ispec, 0, params, smpi );
    
    // Other directions
    for( unsigned int iDim=1 ; iDim<params.nDim_field ; iDim++ ) {
        SyncVectorPatch::initExchParticlesAlongDimension( vecPatches, ispec, iDim, params, smpi );
        SyncVectorPatch::finalizeExchParticlesAlongDimension( vecPatches, ispec, iDim, params, smpi );
    }
    
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->importAndSortParticles( ispec, params );
    }


    /*
    // Debugging
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
    {
        unsigned int npart = vecPatches(ipatch)->vecSpecies[ispec]->particles->size();
        for( unsigned int i=0; i<npart; i++ ) {
            if (vecPatches(ipatch)->vecSpecies[ispec]->particles->position(0,i)< vecPatches(ipatch)->getDomainLocalMin(0)
             || vecPatches(ipatch)->vecSpecies[ispec]->particles->position(0,i) > vecPatches(ipatch)->getDomainLocalMax(0)
             || vecPatches(ipatch)->vecSpecies[ispec]->particles->position(1,i) < vecPatches(ipatch)->getDomainLocalMin(1)
             || vecPatches(ipatch)->vecSpecies[ispec]->particles->position(1,i) > vecPatches(ipatch)->getDomainLocalMax(1))
             {
            cerr << setprecision(12)
                 << " Patch: " << ipatch << "/" << vecPatches.size()-1
                 << " Species: " << ispec
                 << " ipart: " << i
                 << " " << vecPatches(ipatch)->vecSpecies[ispec]->particles->weight(i)
                 << " " << vecPatches(ipatch)->vecSpecies[ispec]->particles->charge(0)
                 << " " << vecPatches(ipatch)->getDomainLocalMin(0)
                 << "<" << vecPatches(ipatch)->vecSpecies[ispec]->particles->position(0,i)
                 << "<" << vecPatches(ipatch)->getDomainLocalMax(0)
                 << " " << vecPatches(ipatch)->getDomainLocalMin(1)
                 << "<" << vecPatches(ipatch)->vecSpecies[ispec]->particles->position(1,i)
                 << "<" << vecPatches(ipatch)->getDomainLocalMax(1)
                 << " " << vecPatches(ipatch)->vecSpecies[ispec]->particles->momentum(0,i)
                 << " " << vecPatches(ipatch)->vecSpecies[ispec]->particles->momentum(1,i)
                 << std::endl;
            }
        }
    }*/

}

void SyncVectorPatch::initExchParticlesAlongDimension( VectorPatch &vecPatches, int ispec, int iDim, Params &params, SmileiMPI *smpi )
{
    // Exchange numbers of particles in direction 0 only
#ifndef _NO_MPI_TM
    #pragma omp for schedule(runtime)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->exchNbrOfParticles( smpi, ispec, params, iDim, &vecPatches );
    }
}

void SyncVectorPatch::finalizeExchParticlesAlongDimension( VectorPatch &vecPatches, int ispec, int iDim, Params &params, SmileiMPI *smpi )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(runtime)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->endNbrOfParticles( ispec, iDim );
    }

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->prepareParticles( smpi, ispec, params, iDim, &vecPatches );
    }

#ifndef _NO_MPI_TM
    #pragma omp for schedule(runtime)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->exchParticles( smpi, ispec, params, iDim, &vecPatches );
    }

#ifndef _NO_MPI_TM
    #pragma omp for schedule(runtime)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->waitExchParticles( ispec, iDim );
    }

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->cornersParticles( ispec, params, iDim );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------       DENSITIES         ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

void SyncVectorPatch::sumRhoJ( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    // Sum Jx, Jy and Jz
    SyncVectorPatch::sumAllComponents( vecPatches.densities, vecPatches, smpi );
    // Sum rho
    if( ( vecPatches.diag_flag ) || ( params.is_spectral ) ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listrho_, vecPatches, smpi );
    }
}

void SyncVectorPatch::sumEnvChi( Params &, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    // Sum Env_Chi
    SyncVectorPatch::sum<double,Field>( vecPatches.listEnv_Chi_, vecPatches, smpi );
}

//sumRhoJ for AM geometry
void SyncVectorPatch::sumRhoJ( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
{
    SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJl_[imode], vecPatches, smpi );
    SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJr_[imode], vecPatches, smpi );
    SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJt_[imode], vecPatches, smpi );
    if( ( vecPatches.diag_flag ) || ( params.is_spectral ) ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listrho_AM_[imode], vecPatches, smpi );
        if (params.is_spectral)
            SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listrho_old_AM_[imode], vecPatches, smpi );
    }
}

void SyncVectorPatch::sumRhoJs( Params &, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    //  Sum Jx_s(ispec), Jy_s(ispec), Jz_s(ispec)
    if( vecPatches.listJxs_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listJxs_, vecPatches, smpi );
    }
    if( vecPatches.listJys_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listJys_, vecPatches, smpi );
    }
    if( vecPatches.listJzs_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listJzs_, vecPatches, smpi );
    }
    // Sum rho_s(ispec)
    if( vecPatches.listrhos_.size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listrhos_, vecPatches, smpi );
    }
}

void SyncVectorPatch::sumEnvChis( Params &, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    // Sum EnvChi_s(ispec)
    if( vecPatches.listEnv_Chis_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listEnv_Chis_, vecPatches, smpi );
    }

}
void SyncVectorPatch::sumRhoJs( Params &, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
{
    // Sum Jx_s(ispec), Jy_s(ispec) and Jz_s(ispec)
    if( vecPatches.listJls_[imode].size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJls_[imode], vecPatches, smpi );
    }
    if( vecPatches.listJrs_[imode].size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJrs_[imode], vecPatches, smpi );
    }
    if( vecPatches.listJts_[imode] .size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJts_[imode], vecPatches, smpi );
    }
    // Sum rho_s(ispec)
    if( vecPatches.listrhos_AM_[imode].size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listrhos_AM_[imode], vecPatches, smpi );
    }
}

// The idea is to minimize the number of implicit barriers and maximize the workload between barriers
// fields : contains all (Jx then Jy then Jz) components of a field for all patches of vecPatches
//     - fields is not directly used in the exchange process, just to find local neighbor's field
//     - sums are operated on sublists, members of vecPatches, which contains :
//         - densitiesMPIx   : fields which have MPI   neighbor along X
//         - densitiesLocalx : fields which have local neighbor along X (a same field can be adressed by both)
//         - ... for Y and Z
//     - These fields are identified with lists of index MPIxIdx and LocalxIdx (... for Y and Z)
void SyncVectorPatch::sumAllComponents( std::vector<Field *> &fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned int h0, oversize[3], size[3];
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    int nPatches( vecPatches.size() );

    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    size[0] = vecPatches( 0 )->EMfields->size_[0];
    size[1] = vecPatches( 0 )->EMfields->size_[1];
    size[2] = vecPatches( 0 )->EMfields->size_[2];

    int nDim = vecPatches( 0 )->EMfields->Jx_->dims_.size();

    // -----------------
    // Sum per direction :

    // iDim = 0, initialize comms : Isend/Irecv
    unsigned int nPatchMPIx = vecPatches.MPIxIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nPatchMPIx ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIxIdx[ifield];
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, iNeighbor ) ) {
                vecPatches.densitiesMPIx[ifield             ]->create_sub_fields ( 0, iNeighbor, 2*oversize[0]+1+1 ); // +1, Jx dual in X
                vecPatches.densitiesMPIx[ifield+nPatchMPIx  ]->create_sub_fields ( 0, iNeighbor, 2*oversize[0]+1+0 ); // +0, Jy prim in X
                vecPatches.densitiesMPIx[ifield+2*nPatchMPIx]->create_sub_fields ( 0, iNeighbor, 2*oversize[0]+1+0 ); // +0, Jz prim in X
                vecPatches.densitiesMPIx[ifield             ]->extract_fields_sum( 0, iNeighbor, oversize[0] );
                vecPatches.densitiesMPIx[ifield+nPatchMPIx  ]->extract_fields_sum( 0, iNeighbor, oversize[0] );
                vecPatches.densitiesMPIx[ifield+2*nPatchMPIx]->extract_fields_sum( 0, iNeighbor, oversize[0] );
// #ifdef SMILEI_ACCELERATOR_GPU_OACC
//                 Field* field = vecPatches.densitiesMPIx[ifield      ];
//                 double* Jx   = field->sendFields_[iNeighbor]->data_;
//                 int sizeofJx = field->sendFields_[iNeighbor]->size();
//                 field = vecPatches.densitiesMPIx[ifield+nPatchMPIx  ];
//                 double* Jy   = field->sendFields_[iNeighbor]->data_;
//                 int sizeofJy = field->sendFields_[iNeighbor]->size();
//                 field = vecPatches.densitiesMPIx[ifield+2*nPatchMPIx];
//                 double*   Jz = field->sendFields_[iNeighbor]->data_;
//                 int sizeofJz = field->sendFields_[iNeighbor]->size();
//                 //#pragma acc update host( Jx[0:sizeofJx], Jy[0:sizeofJy], Jz[0:sizeofJz] )
// #endif
            }
        }
        vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIx[ifield             ], 0, smpi, true ); // Jx
        vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIx[ifield+  nPatchMPIx], 0, smpi, true ); // Jy
        vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIx[ifield+2*nPatchMPIx], 0, smpi, true ); // Jz
    }

    // iDim = 0, local
    const int nFieldLocalx = vecPatches.densitiesLocalx.size() / 3;

#if defined( SMILEI_ACCELERATOR_GPU )
    // At initialization, we may get a CPU buffer than needs to be handled on the host.
    const bool is_memory_on_device = vecPatches.densitiesLocalx.size() > 0 &&
                                     smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( vecPatches.densitiesLocalx[0]->data() );
#endif

    for( int icomp=0 ; icomp<3 ; icomp++ ) {
        if( nFieldLocalx==0 ) {
            continue;
        }

        unsigned int gsp[3];
        //unsigned int nx_ =  vecPatches.densitiesLocalx[icomp*nFieldLocalx]->dims_[0];
        unsigned int ny_ = 1;
        unsigned int nz_ = 1;
        if( nDim>1 ) {
            ny_ = vecPatches.densitiesLocalx[icomp*nFieldLocalx]->dims_[1];
            if( nDim>2 ) {
                nz_ = vecPatches.densitiesLocalx[icomp*nFieldLocalx]->dims_[2];
            }
        }
        gsp[0] = 1+2*oversize[0]+vecPatches.densitiesLocalx[icomp*nFieldLocalx]->isDual_[0]; //Ghost size primal

        unsigned int istart =  icomp   *nFieldLocalx;
        unsigned int iend    = ( icomp+1 )*nFieldLocalx;
        #pragma omp for schedule(static) private(pt1,pt2)
        for( unsigned int ifield=istart ; ifield<iend ; ifield++ ) {
            int ipatch = vecPatches.LocalxIdx[ ifield-icomp*nFieldLocalx ];
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
                pt1 = &( fields[ vecPatches( ipatch )->neighbor_[0][0]-h0+icomp*nPatches ]->data_[size[0]*ny_*nz_] );
                pt2 = &( vecPatches.densitiesLocalx[ifield]->data_[0] );
                //Sum 2 ==> 1

                const unsigned int last = gsp[0] * ny_ * nz_;

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
                int ptsize = vecPatches.densitiesLocalx[ifield]->size();
                int nspace0 = size[0];
                #pragma acc parallel if ( is_memory_on_device) present(pt1[0-nspace0*ny_*nz_:ptsize],pt2[0:ptsize])
                #pragma acc loop worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                #pragma omp target if( is_memory_on_device )
                #pragma omp teams distribute parallel for
#endif
                for( unsigned int i = 0; i < last; i++ ) {
                    pt1[i] += pt2[i];
                    pt2[i]  = pt1[i];
                }
                //Copy back the results to 2
                //memcpy( pt2, pt1, gsp[0]*ny_*nz_*sizeof( double ) );
            }
        }
    }

    // iDim = 0, finalize (waitall)
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nPatchMPIx ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIxIdx[ifield];
        vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIx[ifield             ], 0 ); // Jx
        vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIx[ifield+nPatchMPIx  ], 0 ); // Jy
        vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIx[ifield+2*nPatchMPIx], 0 ); // Jz
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, ( iNeighbor+1 )%2 ) ) {
// #ifdef SMILEI_ACCELERATOR_GPU_OACC
//                 Field* field = vecPatches.densitiesMPIx[ifield      ];
//                 double* Jx   = field->recvFields_[(iNeighbor+1)%2]->data_;
//                 int sizeofJx = field->recvFields_[(iNeighbor+1)%2]->size();
//                 field = vecPatches.densitiesMPIx[ifield+nPatchMPIx  ];
//                 double* Jy   = field->recvFields_[(iNeighbor+1)%2]->data_;
//                 int sizeofJy = field->recvFields_[(iNeighbor+1)%2]->size();
//                 field = vecPatches.densitiesMPIx[ifield+2*nPatchMPIx];
//                 double*   Jz = field->recvFields_[(iNeighbor+1)%2]->data_;
//                 int sizeofJz = field->recvFields_[(iNeighbor+1)%2]->size();
//                 //#pragma acc update device( Jx[0:sizeofJx], Jy[0:sizeofJy], Jz[0:sizeofJz] )
// #endif
                vecPatches.densitiesMPIx[ifield             ]->inject_fields_sum( 0, iNeighbor, oversize[0] );
                vecPatches.densitiesMPIx[ifield+nPatchMPIx  ]->inject_fields_sum( 0, iNeighbor, oversize[0] );
                vecPatches.densitiesMPIx[ifield+2*nPatchMPIx]->inject_fields_sum( 0, iNeighbor, oversize[0] );
            }
        }


   }
    // END iDim = 0 sync
    // -----------------

    if( nDim>1 ) {
        // -----------------
        // Sum per direction :

        // iDim = 1, initialize comms : Isend/Irecv
        unsigned int nPatchMPIy = vecPatches.MPIyIdx.size();
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ifield=0 ; ifield<nPatchMPIy ; ifield++ ) {
            unsigned int ipatch = vecPatches.MPIyIdx[ifield];
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, iNeighbor ) ) {
                    vecPatches.densitiesMPIy[ifield             ]->create_sub_fields ( 1, iNeighbor, 2*oversize[1]+1+0 ); // +0, Jx dual in Y
                    vecPatches.densitiesMPIy[ifield+nPatchMPIy  ]->create_sub_fields ( 1, iNeighbor, 2*oversize[1]+1+1 ); // +1, Jy prim in Y
                    vecPatches.densitiesMPIy[ifield+2*nPatchMPIy]->create_sub_fields ( 1, iNeighbor, 2*oversize[1]+1+0 ); // +0, Jz prim in Y
                    vecPatches.densitiesMPIy[ifield             ]->extract_fields_sum( 1, iNeighbor, oversize[1] );
                    vecPatches.densitiesMPIy[ifield+nPatchMPIy  ]->extract_fields_sum( 1, iNeighbor, oversize[1] );
                    vecPatches.densitiesMPIy[ifield+2*nPatchMPIy]->extract_fields_sum( 1, iNeighbor, oversize[1] );
// #ifdef SMILEI_ACCELERATOR_GPU_OACC
//                     Field* field = vecPatches.densitiesMPIy[ifield      ];
//                     double* Jx   = field->sendFields_[iNeighbor+2]->data_;
//                     int sizeofJx = field->sendFields_[iNeighbor+2]->size();
//                     field = vecPatches.densitiesMPIy[ifield+nPatchMPIy  ];
//                     double* Jy   = field->sendFields_[iNeighbor+2]->data_;
//                     int sizeofJy = field->sendFields_[iNeighbor+2]->size();
//                     field = vecPatches.densitiesMPIy[ifield+2*nPatchMPIy];
//                     double*   Jz = field->sendFields_[iNeighbor+2]->data_;
//                     int sizeofJz = field->sendFields_[iNeighbor+2]->size();
//                     //#pragma acc update host( Jx[0:sizeofJx], Jy[0:sizeofJy], Jz[0:sizeofJz] )
// #endif
                }
            }
            vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIy[ifield             ], 1, smpi, true ); // Jx
            vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIy[ifield+nPatchMPIy  ], 1, smpi, true ); // Jy
            vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIy[ifield+2*nPatchMPIy], 1, smpi, true ); // Jz
        }

        // iDim = 1,
        const int nFieldLocaly = vecPatches.densitiesLocaly.size() / 3;

#if defined( SMILEI_ACCELERATOR_GPU )
        const bool is_memory_on_device = vecPatches.densitiesLocaly.size() > 0 &&
                                         smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( vecPatches.densitiesLocaly[0]->data() );
#endif

        for( int icomp=0 ; icomp<3 ; icomp++ ) {
            if( nFieldLocaly==0 ) {
                continue;
            }

            unsigned int gsp[3];
            unsigned int nx_ =  vecPatches.densitiesLocaly[icomp*nFieldLocaly]->dims_[0];
            unsigned int ny_ = 1;
            unsigned int nz_ = 1;
            if( nDim>1 ) {
                ny_ = vecPatches.densitiesLocaly[icomp*nFieldLocaly]->dims_[1];
                if( nDim>2 ) {
                    nz_ = vecPatches.densitiesLocaly[icomp*nFieldLocaly]->dims_[2];
                }
            }
            gsp[0] = 1+2*oversize[0]+vecPatches.densitiesLocaly[icomp*nFieldLocaly]->isDual_[0]; //Ghost size primal
            gsp[1] = 1+2*oversize[1]+vecPatches.densitiesLocaly[icomp*nFieldLocaly]->isDual_[1]; //Ghost size primal

            unsigned int istart =  icomp   *nFieldLocaly;
            unsigned int iend    = ( icomp+1 )*nFieldLocaly;
            #pragma omp for schedule(static) private(pt1,pt2)
            for( unsigned int ifield=istart ; ifield<iend ; ifield++ ) {
                int ipatch = vecPatches.LocalyIdx[ ifield-icomp*nFieldLocaly ];
                if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
                    //The patch to the south belongs to the same MPI process than I.
                    pt1 = &( fields[vecPatches( ipatch )->neighbor_[1][0]-h0+icomp*nPatches]->data_[size[1]*nz_] );
                    pt2 = &( vecPatches.densitiesLocaly[ifield]->data_[0] );

                    const unsigned int outer_last   = nx_ * ny_ * nz_;
                    const unsigned int outer_stride = ny_ * nz_;
                    const unsigned int inner_last   = gsp[1] * nz_;

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
                    int ptsize = vecPatches.densitiesLocaly[ifield]->size();
                    int blabla = size[1];
                    #pragma acc parallel if (is_memory_on_device) present(pt1[0-blabla*nz_:ptsize],pt2[0:ptsize])
                    #pragma acc loop worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                    #pragma omp target if( is_memory_on_device )
                    #pragma omp teams distribute parallel for collapse(2)
#endif
                    for( unsigned int j = 0; j < outer_last; j += outer_stride ) {
                        for( unsigned int i = 0; i < inner_last; i++ ) {
                            pt1[i+j] += pt2[i+j];
                            pt2[i+j]  = pt1[i+j];
                        }
                        //memcpy( pt2, pt1, gsp[1]*nz_*sizeof( double ) );
                        //pt1 += ny_*nz_;
                        //pt2 += ny_*nz_;
                    }
                }
            }
        }

        // iDim = 1, finalize (waitall)
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ifield=0 ; ifield<nPatchMPIy ; ifield=ifield+1 ) {
            unsigned int ipatch = vecPatches.MPIyIdx[ifield];
            vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIy[ifield             ], 1 ); // Jx
            vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIy[ifield+nPatchMPIy  ], 1 ); // Jy
            vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIy[ifield+2*nPatchMPIy], 1 ); // Jz
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, ( iNeighbor+1 )%2 ) ) {
// #ifdef SMILEI_ACCELERATOR_GPU_OACC
//                     Field* field = vecPatches.densitiesMPIy[ifield      ];
//                     double* Jx   = field->recvFields_[(iNeighbor+1)%2+2]->data_;
//                     int sizeofJx = field->recvFields_[(iNeighbor+1)%2+2]->size();
//                     field = vecPatches.densitiesMPIy[ifield+nPatchMPIy  ];
//                     double* Jy   = field->recvFields_[(iNeighbor+1)%2+2]->data_;
//                     int sizeofJy = field->recvFields_[(iNeighbor+1)%2+2]->size();
//                     field = vecPatches.densitiesMPIy[ifield+2*nPatchMPIy];
//                     double*   Jz = field->recvFields_[(iNeighbor+1)%2+2]->data_;
//                     int sizeofJz = field->recvFields_[(iNeighbor+1)%2+2]->size();
//                     //#pragma acc update device( Jx[0:sizeofJx], Jy[0:sizeofJy], Jz[0:sizeofJz] )
// #endif
                    vecPatches.densitiesMPIy[ifield             ]->inject_fields_sum( 1, iNeighbor, oversize[1] );
                    vecPatches.densitiesMPIy[ifield+nPatchMPIy  ]->inject_fields_sum( 1, iNeighbor, oversize[1] );
                    vecPatches.densitiesMPIy[ifield+2*nPatchMPIy]->inject_fields_sum( 1, iNeighbor, oversize[1] );
                }
            }
        }
        // END iDim = 1 sync
        // -----------------

        if( nDim>2 ) {
            // -----------------
            // Sum per direction :

            // iDim = 2, initialize comms : Isend/Irecv
            unsigned int nPatchMPIz = vecPatches.MPIzIdx.size();
#ifndef _NO_MPI_TM
            #pragma omp for schedule(static)
#else
            #pragma omp single
#endif
            for( unsigned int ifield=0 ; ifield<nPatchMPIz ; ifield++ ) {
                unsigned int ipatch = vecPatches.MPIzIdx[ifield];
                for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                    if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, iNeighbor ) ) {
                        vecPatches.densitiesMPIz[ifield             ]->create_sub_fields ( 2, iNeighbor, 2*oversize[2]+1+0 ); // +0, Jx dual in Z
                        vecPatches.densitiesMPIz[ifield+nPatchMPIz  ]->create_sub_fields ( 2, iNeighbor, 2*oversize[2]+1+0 ); // +0, Jy prim in Z
                        vecPatches.densitiesMPIz[ifield+2*nPatchMPIz]->create_sub_fields ( 2, iNeighbor, 2*oversize[2]+1+1 ); // +1, Jz prim in Z
                        vecPatches.densitiesMPIz[ifield             ]->extract_fields_sum( 2, iNeighbor, oversize[2] );
                        vecPatches.densitiesMPIz[ifield+nPatchMPIz  ]->extract_fields_sum( 2, iNeighbor, oversize[2] );
                        vecPatches.densitiesMPIz[ifield+2*nPatchMPIz]->extract_fields_sum( 2, iNeighbor, oversize[2] );
// #ifdef SMILEI_ACCELERATOR_GPU_OACC
//                         Field* field = vecPatches.densitiesMPIz[ifield      ];
//                         double* Jx   = field->sendFields_[iNeighbor+4]->data_;
//                         int sizeofJx = field->sendFields_[iNeighbor+4]->size();
//                         field = vecPatches.densitiesMPIz[ifield+nPatchMPIz  ];
//                         double* Jy   = field->sendFields_[iNeighbor+4]->data_;
//                         int sizeofJy = field->sendFields_[iNeighbor+4]->size();
//                         field = vecPatches.densitiesMPIz[ifield+2*nPatchMPIz];
//                         double*   Jz = field->sendFields_[iNeighbor+4]->data_;
//                         int sizeofJz = field->sendFields_[iNeighbor+4]->size();
//                         //#pragma acc update host( Jx[0:sizeofJx], Jy[0:sizeofJy], Jz[0:sizeofJz] )
// #endif
                    }
                }
                vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIz[ifield             ], 2, smpi, true ); // Jx
                vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIz[ifield+nPatchMPIz  ], 2, smpi, true ); // Jy
                vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIz[ifield+2*nPatchMPIz], 2, smpi, true ); // Jz
            }

            // iDim = 2 local
            const int nFieldLocalz = vecPatches.densitiesLocalz.size() / 3;

#if defined( SMILEI_ACCELERATOR_GPU )
            const bool is_memory_on_device = vecPatches.densitiesLocalz.size() > 0 &&
                                             smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( vecPatches.densitiesLocalz[0]->data() );
#endif

            for( int icomp=0 ; icomp<3 ; icomp++ ) {
                if( nFieldLocalz==0 ) {
                    continue;
                }

                unsigned int gsp[3];
                unsigned int nx_ =  vecPatches.densitiesLocalz[icomp*nFieldLocalz]->dims_[0];
                unsigned int ny_ = 1;
                unsigned int nz_ = 1;
                if( nDim>1 ) {
                    ny_ = vecPatches.densitiesLocalz[icomp*nFieldLocalz]->dims_[1];
                    if( nDim>2 ) {
                        nz_ = vecPatches.densitiesLocalz[icomp*nFieldLocalz]->dims_[2];
                    }
                }
                gsp[0] = 1+2*oversize[0]+vecPatches.densitiesLocalz[icomp*nFieldLocalz]->isDual_[0]; //Ghost size primal
                gsp[1] = 1+2*oversize[1]+vecPatches.densitiesLocalz[icomp*nFieldLocalz]->isDual_[1]; //Ghost size primal
                gsp[2] = 1+2*oversize[2]+vecPatches.densitiesLocalz[icomp*nFieldLocalz]->isDual_[2]; //Ghost size primal

                unsigned int istart  =  icomp   *nFieldLocalz;
                unsigned int iend    = ( icomp+1 )*nFieldLocalz;
                #pragma omp for schedule(static) private(pt1,pt2)
                for( unsigned int ifield=istart ; ifield<iend ; ifield++ ) {
                    int ipatch = vecPatches.LocalzIdx[ ifield-icomp*nFieldLocalz ];
                    if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
                        //The patch below me belongs to the same MPI process than I.
                        pt1 = &( fields[vecPatches( ipatch )->neighbor_[2][0]-h0+icomp*nPatches]->data_[size[2]] );
                        pt2 = &( vecPatches.densitiesLocalz[ifield]->data_[0] );

                        const unsigned int outer_last   = nx_ * ny_ * nz_;
                        const unsigned int outer_stride = nz_;
                        const unsigned int inner_last   = gsp[2];

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
                        int ptsize = vecPatches.densitiesLocalz[ifield]->size();
                        int blabla = size[2];
                        #pragma acc parallel if (is_memory_on_device) present(pt1[0-blabla:ptsize],pt2[0:ptsize])
                        #pragma acc loop worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                        #pragma omp target if( is_memory_on_device )
                         #pragma omp teams distribute parallel for collapse( 2 )
#endif
                        for( unsigned int j = 0; j < outer_last; j += outer_stride ) {
                            for( unsigned int i = 0; i < inner_last; i++ ) {
                                pt1[i+j] += pt2[i+j];
                                pt2[i+j] =  pt1[i+j];
                            }
                        }
                    }
                }
            }

            // iDim = 2, complete non local sync through MPIfinalize (waitall)
#ifndef _NO_MPI_TM
            #pragma omp for schedule(static)
#else
            #pragma omp single
#endif
            for( unsigned int ifield=0 ; ifield<nPatchMPIz ; ifield=ifield+1 ) {
                unsigned int ipatch = vecPatches.MPIzIdx[ifield];
                vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIz[ifield             ], 2 ); // Jx
                vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIz[ifield+nPatchMPIz  ], 2 ); // Jy
                vecPatches( ipatch )->finalizeSumField( vecPatches.densitiesMPIz[ifield+2*nPatchMPIz], 2 ); // Jz
                for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                    if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, ( iNeighbor+1 )%2 ) ) {
// #ifdef SMILEI_ACCELERATOR_GPU_OACC
//                         Field* field = vecPatches.densitiesMPIz[ifield      ];
//                         double* Jx   = field->recvFields_[(iNeighbor+1)%2+4]->data_;
//                         int sizeofJx = field->recvFields_[(iNeighbor+1)%2+4]->size();
//                         field = vecPatches.densitiesMPIz[ifield+nPatchMPIz  ];
//                         double* Jy   = field->recvFields_[(iNeighbor+1)%2+4]->data_;
//                         int sizeofJy = field->recvFields_[(iNeighbor+1)%2+4]->size();
//                         field = vecPatches.densitiesMPIz[ifield+2*nPatchMPIz ];
//                         double*   Jz = field->recvFields_[(iNeighbor+1)%2+4]->data_;
//                         int sizeofJz = field->recvFields_[(iNeighbor+1)%2+4]->size();
//                         //#pragma acc update device( Jx[0:sizeofJx], Jy[0:sizeofJy], Jz[0:sizeofJz] )
// #endif
                        vecPatches.densitiesMPIz[ifield             ]->inject_fields_sum( 2, iNeighbor, oversize[2] );
                        vecPatches.densitiesMPIz[ifield+nPatchMPIz  ]->inject_fields_sum( 2, iNeighbor, oversize[2] );
                        vecPatches.densitiesMPIz[ifield+2*nPatchMPIz]->inject_fields_sum( 2, iNeighbor, oversize[2] );
                    }
                }
            }
            // END iDim = 2 sync
            // -----------------

        } // End if dims_.size()>2

    } // End if dims_.size()>1

}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------         FIELDS          ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

void SyncVectorPatch::exchangeE( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    // full_B_exchange is true if (Buneman BC, Lehe, Bouchard or spectral solvers)
    // E is exchange if spectral solver and/or at the end of initialisation of non-neutral plasma

    if( !params.full_B_exchange ) {
        SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEx_, vecPatches, smpi );
        SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEy_, vecPatches, smpi );
        SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEz_, vecPatches, smpi );
    } else {
        SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listEx_, vecPatches, smpi );
        SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listEy_, vecPatches, smpi );
        SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listEz_, vecPatches, smpi );
    }

}

void SyncVectorPatch::finalizeexchangeE( Params &params, VectorPatch &vecPatches )
{
    // full_B_exchange is true if (Buneman BC, Lehe, Bouchard or spectral solvers)
    // E is exchange if spectral solver and/or at the end of initialisation of non-neutral plasma

    if( !params.full_B_exchange ) {
        SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEx_, vecPatches );
        SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEy_, vecPatches );
        SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEz_, vecPatches );
    }
    //else
    //    done in exchangeSynchronizedPerDirection
}

void SyncVectorPatch::exchangeB( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    // full_B_exchange is true if (Buneman BC, Lehe, Bouchard or spectral solvers)

    if( vecPatches.listBx_[0]->dims_.size()==1 ) {
        // Exchange Bs0 : By_ and Bz_ (dual in X)
        SyncVectorPatch::exchangeAllComponentsAlongX( vecPatches.Bs0, vecPatches, smpi );
    } else {
        if( params.full_B_exchange ) {
            // Exchange Bx_ in Y then X
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listBx_, vecPatches, smpi );
            // Exchange By_ in Y then X
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listBy_, vecPatches, smpi );
            // Exchange Bz_ in Y then X
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listBz_, vecPatches, smpi );

        } else {
            if( vecPatches.listBx_[0]->dims_.size()==2 ) {
                // Exchange Bs0 : By_ and Bz_ (dual in X)
                SyncVectorPatch::exchangeAllComponentsAlongX( vecPatches.Bs0, vecPatches, smpi );
                // Exchange Bs1 : Bx_ and Bz_ (dual in Y)
                SyncVectorPatch::exchangeAllComponentsAlongY( vecPatches.Bs1, vecPatches, smpi );
            } else if( vecPatches.listBx_[0]->dims_.size()==3 ) {
                // Exchange Bs0 : By_ and Bz_ (dual in X)
                SyncVectorPatch::exchangeAllComponentsAlongX( vecPatches.Bs0, vecPatches, smpi );
                // Exchange Bs1 : Bx_ and Bz_ (dual in Y)
                SyncVectorPatch::exchangeAllComponentsAlongY( vecPatches.Bs1, vecPatches, smpi );
                // Exchange Bs2 : Bx_ and By_ (dual in Z)
                SyncVectorPatch::exchangeAllComponentsAlongZ( vecPatches.Bs2, vecPatches, smpi );
            }
        }
    }
}

void SyncVectorPatch::finalizeexchangeB( Params &params, VectorPatch &vecPatches )
{
    // full_B_exchange is true if (Buneman BC, Lehe, Bouchard or spectral solvers)

    if( vecPatches.listBx_[0]->dims_.size()==1 ) {
        // Finalize exchange Bs0 : By_ and Bz_ (dual in X)
        SyncVectorPatch::finalizeExchangeAllComponentsAlongX( vecPatches );
    } else if( vecPatches.listBx_[0]->dims_.size()==2 ) {
        if( !params.full_B_exchange ) {
            // Finalize exchange Bs0 : By_ and Bz_ (dual in X)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongX( vecPatches );
            // Finalize exchange Bs1 : Bx_ and Bz_ (dual in Y)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongY( vecPatches );
        }
        //else
        //    done in exchangeSynchronizedPerDirection
    } else if( vecPatches.listBx_[0]->dims_.size()==3 ) {
        if( !params.full_B_exchange ) {
            // Finalize exchange Bs0 : By_ and Bz_ (dual in X)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongX( vecPatches );
            // Finalize exchange Bs1 : Bx_ and Bz_ (dual in Y)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongY( vecPatches );
            // Finalize exchange Bs2 : Bx_ and By_ (dual in Z)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongZ( vecPatches );
        }
        //else
        //    done in exchangeSynchronizedPerDirection
    }

}

void SyncVectorPatch::exchangeJ( Params &, VectorPatch &vecPatches, SmileiMPI *smpi )
{

    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listJx_, vecPatches, smpi );
    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listJy_, vecPatches, smpi );
    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listJz_, vecPatches, smpi );
}

void SyncVectorPatch::finalizeexchangeJ( Params &, VectorPatch &vecPatches )
{

    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listJx_, vecPatches );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listJy_, vecPatches );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listJz_, vecPatches );
}


void SyncVectorPatch::exchangeB( Params &, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
{
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBl_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBl_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBr_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBr_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBt_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBt_[imode], vecPatches );
}

void SyncVectorPatch::exchangeE( Params &, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
{
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listEl_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEl_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listEr_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEr_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listEt_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEt_[imode], vecPatches );
}

void SyncVectorPatch::exchangeBmBTIS3( Params &/*params*/, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
{
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBr_mBTIS3[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBr_mBTIS3[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBt_mBTIS3[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBt_mBTIS3[imode], vecPatches );
}

// void SyncVectorPatch::finalizeexchangeB( Params &, VectorPatch &, int )
// {
// }
// 
// void SyncVectorPatch::finalizeexchangeE( Params &, VectorPatch &, int )
// {
// }


//void SyncVectorPatch::exchangeB( Params &, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
//{
//
//    SyncVectorPatch::exchangeComplex( vecPatches.listBl_[imode], vecPatches, smpi );
//    SyncVectorPatch::exchangeComplex( vecPatches.listBr_[imode], vecPatches, smpi );
//    SyncVectorPatch::exchangeComplex( vecPatches.listBt_[imode], vecPatches, smpi );
//}

//void SyncVectorPatch::finalizeexchangeB( Params &, VectorPatch &vecPatches, int imode )
//{
//
//    SyncVectorPatch::finalizeexchangeComplex( vecPatches.listBl_[imode], vecPatches );
//    SyncVectorPatch::finalizeexchangeComplex( vecPatches.listBr_[imode], vecPatches );
//    SyncVectorPatch::finalizeexchangeComplex( vecPatches.listBt_[imode], vecPatches );
//}

void SyncVectorPatch::exchangeA( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    if( !params.full_Envelope_exchange ) {
        // current envelope value
        SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listA_, vecPatches, smpi );
        SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listA_, vecPatches );
        // value of envelope at previous timestep
        SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listA0_, vecPatches, smpi );
        SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listA0_, vecPatches );
    } else {
        // current envelope value
        SyncVectorPatch::exchangeSynchronizedPerDirection<complex<double>,cField>( vecPatches.listA_, vecPatches, smpi );
        // value of envelope at previous timestep
        SyncVectorPatch::exchangeSynchronizedPerDirection<complex<double>,cField>( vecPatches.listA0_, vecPatches, smpi );
    }
}

// void SyncVectorPatch::finalizeexchangeA( Params &, VectorPatch &vecPatches )
// {
//    // current envelope value
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listA_, vecPatches );
//    // current envelope value
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listA0_, vecPatches );
// }

// void SyncVectorPatch::exchangeEnvEEnvA( Params &, VectorPatch &vecPatches, SmileiMPI *smpi )
// {
//     // current envelope |E| value
//     SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEnvE_, vecPatches, smpi );
//     SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEnvE_, vecPatches );
//     // current envelope |A| value
//     SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEnvA_, vecPatches, smpi );
//     SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEnvA_, vecPatches );
// }
// 
// void SyncVectorPatch::finalizeexchangeEnvEEnvA( Params &, VectorPatch &vecPatches )
// {
//
// }

void SyncVectorPatch::exchangeEnvEx( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    // current envelope |Ex| value
    if( !params.full_Envelope_exchange ) {
        SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEnvEx_, vecPatches, smpi );
        SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEnvEx_, vecPatches );
    } else {
        SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listEnvEx_, vecPatches, smpi );
    }
}

void SyncVectorPatch::exchangeBmBTIS3( Params &/*params*/, VectorPatch &vecPatches, SmileiMPI *smpi )
{   // exchange BmBTIS3 in Cartesian geometries

    // exchange ByBTIS3 
    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listBy_mBTIS3, vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBy_mBTIS3, vecPatches );

    // exchange BzBTIS3 
    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listBz_mBTIS3, vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBz_mBTIS3, vecPatches );
}

// void SyncVectorPatch::finalizeexchangeBmBTIS3( Params &params, VectorPatch &vecPatches )
// {
// 
// }

// void SyncVectorPatch::finalizeexchangeEnvEx( Params &, VectorPatch &vecPatches )
// {
// 
// }

// void SyncVectorPatch::exchangePhi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
// {
//
//     if( !params.full_Envelope_exchange ) {
//         // current ponderomotive potential
//         SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listPhi_, vecPatches, smpi );
//         SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listPhi_, vecPatches );
//         // value of ponderomotive potential at previous timestep
//         SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listPhi0_, vecPatches, smpi );
//         SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listPhi0_, vecPatches );
//     } else {
//         // current ponderomotive potential
//         SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listPhi_, vecPatches, smpi );
//         // value of ponderomotive potential at previous timestep
//         SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listPhi0_, vecPatches, smpi );
//     }
//
// }
// 
// void SyncVectorPatch::finalizeexchangePhi( Params &, VectorPatch &vecPatches )
// {
//
// }




void SyncVectorPatch::exchangeGradPhi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    if (  params.geometry != "AMcylindrical" ) {
        if( !params.full_Envelope_exchange ) {
            // current Gradient value
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhix_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhix_, vecPatches );
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhiy_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiy_, vecPatches );
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhiz_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiz_, vecPatches );

            // value of Gradient at previous timestep
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhix0_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhix0_, vecPatches );
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhiy0_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiy0_, vecPatches );
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhiz0_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiz0_, vecPatches );
        } else {
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhix_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhiy_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhiz_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhix0_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhiy0_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhiz0_, vecPatches, smpi );
        }
    } else {
        if( !params.full_Envelope_exchange ) {
            // current Gradient value
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhil_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhil_, vecPatches );
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhir_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhir_, vecPatches );

            // value of Gradient at previous timestep
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhil0_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhil0_, vecPatches );
            SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listGradPhir0_, vecPatches, smpi );
            SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhir0_, vecPatches );
        } else {
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhil_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhil0_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhir_, vecPatches, smpi );
            SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listGradPhir0_, vecPatches, smpi );
        }
    }
}

// void SyncVectorPatch::finalizeexchangeGradPhi( Params &, VectorPatch &vecPatches )
// {
//    // current Gradient value
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhix_, vecPatches );
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiy_, vecPatches );
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiz_, vecPatches );
// }

void SyncVectorPatch::exchangeEnvChi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    if( !params.full_Envelope_exchange ) {
        // susceptibility
        SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEnv_Chi_, vecPatches, smpi );
        SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEnv_Chi_, vecPatches );
        
    } else {
        // susceptibility
        SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( vecPatches.listEnv_Chi_, vecPatches, smpi );
    }
}


void SyncVectorPatch::templateGenerator()
{
    SmileiMPI* smpi = NULL;
    VectorPatch patches;
    SyncVectorPatch::exchangeAlongAllDirectionsNoOMP<double         ,Field >( patches.listEx_, patches, smpi );
    SyncVectorPatch::exchangeAlongAllDirectionsNoOMP<complex<double>,cField>( patches.listEx_, patches, smpi );
    SyncVectorPatch::exchangeAlongAllDirectionsNoOMP<double         ,Field >( patches.listEx_, patches, smpi );
    SyncVectorPatch::exchangeAlongAllDirectionsNoOMP<double         ,Field >( patches.listEx_, patches, smpi );
}

// fields : contains a single field component (X, Y or Z) for all patches of vecPatches
template<typename T, typename F>
void SyncVectorPatch::exchangeAlongAllDirections( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned int  oversize[3];
    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                    fields[ipatch]->create_sub_fields  ( iDim, iNeighbor, oversize[iDim] );
                    fields[ipatch]->extract_fields_exch( iDim, iNeighbor, oversize[iDim] );
                }
            }
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange       ( fields[ipatch], iDim, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], iDim, smpi );
        }
    } // End for iDim

    unsigned int nx_, ny_( 1 ), nz_( 1 ), h0, size[3], gsp[3];
    T *pt1, *pt2;
    F *field1, *field2;
    h0 = vecPatches( 0 )->hindex;

    size[0] = vecPatches( 0 )->EMfields->size_[0];
    size[1] = vecPatches( 0 )->EMfields->size_[1];
    size[2] = vecPatches( 0 )->EMfields->size_[2];

    nx_ = fields[0]->dims_[0];
    if( fields[0]->dims_.size()>1 ) {
        ny_ = fields[0]->dims_[1];
        if( fields[0]->dims_.size()>2 ) {
            nz_ = fields[0]->dims_[2];
        }
    }


    gsp[0] = ( oversize[0] + 1 + fields[0]->isDual_[0] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
            field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[0][0]-h0] );
            field2 = static_cast<F *>( fields[ipatch] );
            pt1 = &( *field1 )( ( size[0] )*ny_*nz_ );
            pt2 = &( *field2 )( 0 );
            memcpy( pt2, pt1, oversize[0]*ny_*nz_*sizeof( T ) );
            memcpy( pt1+gsp[0]*ny_*nz_, pt2+gsp[0]*ny_*nz_, oversize[0]*ny_*nz_*sizeof( T ) );
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] )

        if( fields[0]->dims_.size()>1 ) {
            gsp[1] = ( oversize[1] + 1 + fields[0]->isDual_[1] ); //Ghost size primal
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
                field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[1][0]-h0] );
                field2 = static_cast<F *>( fields[ipatch] );
                pt1 = &( *field1 )( size[1]*nz_ );
                pt2 = &( *field2 )( 0 );
                for( unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_ ) {
                    for( unsigned int j = 0 ; j < oversize[1]*nz_ ; j++ ) {
                        // Rewrite with memcpy ?
                        pt2[i+j] = pt1[i+j] ;
                        pt1[i+j+gsp[1]*nz_] = pt2[i+j+gsp[1]*nz_] ;
                    }
                }
            } // End if ( MPI_me_ == MPI_neighbor_[1][0] )

            if( fields[0]->dims_.size()>2 ) {
                gsp[2] = ( oversize[2] + 1 + fields[0]->isDual_[2] ); //Ghost size primal
                if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
                    field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[2][0]-h0] );
                    field2 = static_cast<F *>( fields[ipatch] );
                    pt1 = &( *field1 )( size[2] );
                    pt2 = &( *field2 )( 0 );
                    for( unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_ ) {
                        for( unsigned int j = 0 ; j < ny_*nz_ ; j += nz_ ) {
                            for( unsigned int k = 0 ; k < oversize[2] ; k++ ) {
                                pt2[i+j+k] = pt1[i+j+k] ;
                                pt1[i+j+k+gsp[2]] = pt2[i+j+k+gsp[2]] ;
                            }
                        }
                    }
                }// End if ( MPI_me_ == MPI_neighbor_[2][0] )
            }// End if dims_.size()>2
        } // End if dims_.size()>1
    } // End for( ipatch )

}

// MPI_Wait for all communications initialised in exchangeAlongAllDirections
void SyncVectorPatch::finalizeExchangeAlongAllDirections( std::vector<Field *> fields, VectorPatch &vecPatches )
{
    unsigned oversize[3];
    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], iDim );

            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
                    fields[ipatch]->inject_fields_exch( iDim, iNeighbor, oversize[iDim] );
                }
            }
        }
    } // End for iDim

}


// fields : contains a single field component (X, Y or Z) for all patches of vecPatches
template<typename T, typename F>
void SyncVectorPatch::exchangeAlongAllDirectionsNoOMP( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned oversize[3];
    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                    fields[ipatch]->create_sub_fields  ( iDim, iNeighbor, oversize[iDim] );
                    fields[ipatch]->extract_fields_exch( iDim, iNeighbor, oversize[iDim] );
                }
            }
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange       ( fields[ipatch], iDim, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], iDim, smpi );
        }
    } // End for iDim


    unsigned int nx_, ny_( 1 ), nz_( 1 ), h0, size[3], gsp[3];
    T *pt1, *pt2;
    F* field1;
    F* field2;
    h0 = vecPatches( 0 )->hindex;

    size[0] = vecPatches( 0 )->EMfields->size_[0];
    size[1] = vecPatches( 0 )->EMfields->size_[1];
    size[2] = vecPatches( 0 )->EMfields->size_[2];

    nx_ = fields[0]->dims_[0];
    if( fields[0]->dims_.size()>1 ) {
        ny_ = fields[0]->dims_[1];
        if( fields[0]->dims_.size()>2 ) {
            nz_ = fields[0]->dims_[2];
        }
    }

    gsp[0] = ( oversize[0] + 1 + fields[0]->isDual_[0] ); //Ghost size primal

    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
            field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[0][0]-h0] );
            field2 = static_cast<F *>( fields[ipatch] );
            pt1 = &( *field1 )( size[0]*ny_*nz_ );
            pt2 = &( *field2 )( 0 );
            memcpy( pt2, pt1, oversize[0]*ny_*nz_*sizeof( T ) );
            memcpy( pt1+gsp[0]*ny_*nz_, pt2+gsp[0]*ny_*nz_, oversize[0]*ny_*nz_*sizeof( T ) );
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] )

        if( fields[0]->dims_.size()>1 ) {
            gsp[1] = ( oversize[1] + 1 + fields[0]->isDual_[1] ); //Ghost size primal
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
                field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[1][0]-h0] );
                field2 = static_cast<F *>( fields[ipatch] );
                pt1 = &( *field1 )( size[1]*nz_ );
                pt2 = &( *field2 )( 0 );
                for( unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_ ) {
                    for( unsigned int j = 0 ; j < oversize[1]*nz_ ; j++ ) {
                        // Rewrite with memcpy ?
                        pt2[i+j] = pt1[i+j] ;
                        pt1[i+j+gsp[1]*nz_] = pt2[i+j+gsp[1]*nz_] ;
                    }
                }
            } // End if ( MPI_me_ == MPI_neighbor_[1][0] )

            if( fields[0]->dims_.size()>2 ) {
                gsp[2] = ( oversize[2] + 1 + fields[0]->isDual_[2] ); //Ghost size primal
                if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
                    field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[2][0]-h0] );
                    field2 = static_cast<F *>( fields[ipatch] );
                    pt1 = &( *field1 )( size[2] );
                    pt2 = &( *field2 )( 0 );
                    for( unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_ ) {
                        for( unsigned int j = 0 ; j < ny_*nz_ ; j += nz_ ) {
                            for( unsigned int k = 0 ; k < oversize[2] ; k++ ) {
                                pt2[i+j+k] = pt1[i+j+k] ;
                                pt1[i+j+k+gsp[2]] = pt2[i+j+k+gsp[2]] ;
                            }
                        }
                    }
                }// End if ( MPI_me_ == MPI_neighbor_[2][0] )
            }// End if dims_.size()>2
        } // End if dims_.size()>1
    } // End for( ipatch )

}


// MPI_Wait for all communications initialised in exchangeAlongAllDirections
void SyncVectorPatch::finalizeExchangeAlongAllDirectionsNoOMP( std::vector<Field *> fields, VectorPatch &vecPatches )
{
    unsigned oversize[3];
    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], iDim );

            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
                    fields[ipatch]->inject_fields_exch( iDim, iNeighbor, oversize[iDim] );
                }
            }

        }
    } // End for iDim

}


//Proceed to the synchronization of field including corner ghost cells.
//This is done by exchanging one dimension at a time
template void SyncVectorPatch::exchangeSynchronizedPerDirection<double,Field>( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
template void SyncVectorPatch::exchangeSynchronizedPerDirection<complex<double>,cField>( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );

template<typename T, typename F>
void SyncVectorPatch::exchangeSynchronizedPerDirection( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{

    unsigned int nx_, ny_( 1 ), nz_( 1 ), h0, oversize[3], size[3], gsp[3];
    T *pt1, *pt2;
    F* field1;
    F* field2;
    h0 = vecPatches( 0 )->hindex;

    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    size[0] = vecPatches( 0 )->EMfields->size_[0];
    size[1] = vecPatches( 0 )->EMfields->size_[1];
    size[2] = vecPatches( 0 )->EMfields->size_[2];

    nx_ = fields[0]->dims_[0];
    if( fields[0]->dims_.size()>1 ) {
        ny_ = fields[0]->dims_[1];
        if( fields[0]->dims_.size()>2 ) {
            nz_ = fields[0]->dims_[2];
        }
    }

    if( fields[0]->dims_.size()>2 ) {

        // Dimension 2
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, iNeighbor ) ) {
                    fields[ipatch]->create_sub_fields  ( 2, iNeighbor, oversize[2] );
                    fields[ipatch]->extract_fields_exch( 2, iNeighbor, oversize[2] );
                }
            }
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange( fields[ipatch], 2, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], 2, smpi );
        }
        

#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], 2 );

            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, ( iNeighbor+1 )%2 ) ) {
                    fields[ipatch]->inject_fields_exch( 2, iNeighbor, oversize[2] );
                }
            }
        }

        #pragma omp for schedule(static) private(pt1,pt2)
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

            gsp[2] = ( oversize[2] + 1 + fields[0]->isDual_[2] ); //Ghost size primal
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
                field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[2][0]-h0]  );
                field2 = static_cast<F *>( fields[ipatch] );
                pt1 = &( *field1 )( size[2] );
                pt2 = &( *field2 )( 0 );
                //for (unsigned int in = oversize[0] ; in < nx_-oversize[0]; in ++){
                for( unsigned int in = 0 ; in < nx_ ; in ++ ) {
                    unsigned int i = in * ny_*nz_;
                    //for (unsigned int jn = oversize[1] ; jn < ny_-oversize[1] ; jn ++){
                    for( unsigned int jn = 0 ; jn < ny_ ; jn ++ ) {
                        unsigned int j = jn *nz_;
                        for( unsigned int k = 0 ; k < oversize[2] ; k++ ) {
                            pt2[i+j+k] = pt1[i+j+k] ;
                            pt1[i+j+k+gsp[2]] = pt2[i+j+k+gsp[2]] ;
                        }
                    }
                }
            }// End if ( MPI_me_ == MPI_neighbor_[2][0] )

        } // End for( ipatch )
    }

    // Dimension 1
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, iNeighbor ) ) {
                fields[ipatch]->create_sub_fields  ( 1, iNeighbor, oversize[1] );
                fields[ipatch]->extract_fields_exch( 1, iNeighbor, oversize[1] );
            }
        }
        if ( !dynamic_cast<cField*>( fields[ipatch] ) )
            vecPatches( ipatch )->initExchange( fields[ipatch], 1, smpi );
        else
            vecPatches( ipatch )->initExchangeComplex( fields[ipatch], 1, smpi );
    }

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->finalizeExchange( fields[ipatch], 1 );

        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, ( iNeighbor+1 )%2 ) ) {
                fields[ipatch]->inject_fields_exch( 1, iNeighbor, oversize[1] );
            }
        }
    }

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        gsp[1] = ( oversize[1] + 1 + fields[0]->isDual_[1] ); //Ghost size primal
        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
            field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[1][0]-h0]  );
            field2 = static_cast<F *>( fields[ipatch] );
            pt1 = &( *field1 )( size[1]*nz_ );
            pt2 = &( *field2 )( 0 );
            for( unsigned int in = 0 ; in < nx_ ; in ++ ) {
                //for (unsigned int in = oversize[0] ; in < nx_-oversize[0] ; in ++){ // <== This doesn't work. Why ??
                unsigned int i = in * ny_*nz_;
                for( unsigned int j = 0 ; j < oversize[1]*nz_ ; j++ ) {
                    // Rewrite with memcpy ?
                    pt2[i+j] = pt1[i+j] ;
                    pt1[i+j+gsp[1]*nz_] = pt2[i+j+gsp[1]*nz_] ;
                }
            }
        } // End if ( MPI_me_ == MPI_neighbor_[1][0] )

    } // End for( ipatch )

    // Dimension 0
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, iNeighbor ) ) {
                fields[ipatch]->create_sub_fields  ( 0, iNeighbor, oversize[0] );
                fields[ipatch]->extract_fields_exch( 0, iNeighbor, oversize[0] );
            }
        }
        if ( !dynamic_cast<cField*>( fields[ipatch] ) )
            vecPatches( ipatch )->initExchange( fields[ipatch], 0, smpi );
        else
            vecPatches( ipatch )->initExchangeComplex( fields[ipatch], 0, smpi );
    }

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->finalizeExchange( fields[ipatch], 0 );

        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, ( iNeighbor+1 )%2 ) ) {
                fields[ipatch]->inject_fields_exch( 0, iNeighbor, oversize[0] );
            }
        }
    }



    gsp[0] = ( oversize[0] + 1 + fields[0]->isDual_[0] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
            field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[0][0]-h0]  );
            field2 = static_cast<F *>( fields[ipatch] );
            pt1 = &( *field1 )( ( size[0] )*ny_*nz_ );
            pt2 = &( *field2 )( 0 );
            memcpy( pt2, pt1, oversize[0]*ny_*nz_*sizeof( T ) );
            memcpy( pt1+gsp[0]*ny_*nz_, pt2+gsp[0]*ny_*nz_, oversize[0]*ny_*nz_*sizeof( T ) );
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] )

    } // End for( ipatch )

}


// The idea is to minimize the number of implicit barriers and maximize the workload between barriers
// fields : contains components of B which required exchange in X (By and Bz)
//     - fields is not directly used in the exchange process, just to find local neighbor's field
//     - sexchanges are operated on sublists, members of vecPatches, which contains :
//         - B_MPIx   : fields which have MPI   neighbor along X
//         - B_Localx : fields which have local neighbor along X (a same field can be adressed by both)
//     - These fields are identified with lists of index MPIxIdx and LocalxIdx
void SyncVectorPatch::exchangeAllComponentsAlongX( std::vector<Field *> &fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[0];

    unsigned int nMPIx = vecPatches.MPIxIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIx ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIxIdx[ifield];
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, iNeighbor ) ) {
                vecPatches.B_MPIx[ifield      ]->create_sub_fields  ( 0, iNeighbor, oversize );
                vecPatches.B_MPIx[ifield      ]->extract_fields_exch( 0, iNeighbor, oversize );
                vecPatches.B_MPIx[ifield+nMPIx]->create_sub_fields  ( 0, iNeighbor, oversize );
                vecPatches.B_MPIx[ifield+nMPIx]->extract_fields_exch( 0, iNeighbor, oversize );
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                Field* field = vecPatches.B_MPIx[ifield      ];
                double* By   = field->sendFields_[iNeighbor]->data_;
                int sizeofBy = field->sendFields_[iNeighbor]->size();
                field        = vecPatches.B_MPIx[ifield+nMPIx];
                double* Bz   = field->sendFields_[iNeighbor]->data_;
                int sizeofBz = field->sendFields_[iNeighbor]->size();
                //#pragma acc update host(By[0:sizeofBy],Bz[0:sizeofBz])
#endif
            }
        }
        vecPatches( ipatch )->initExchange( vecPatches.B_MPIx[ifield      ], 0, smpi, true ); // By
        vecPatches( ipatch )->initExchange( vecPatches.B_MPIx[ifield+nMPIx], 0, smpi, true ); // Bz
    }

    unsigned int h0, size;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    size = vecPatches( 0 )->EMfields->size_[0];

    int nPatches( vecPatches.size() );
    int nDim = vecPatches( 0 )->EMfields->Bx_->dims_.size();

    // TODO(Etienne M): Can we somehow get a CPU pointer when GPU mode is enabled ? If not, remove the
    // is_memory_on_device check.
    const bool is_memory_on_device = vecPatches.B_localx.size() > 0 &&
                                     smilei::tools::gpu::HostDeviceMemoryManagement::IsHostPointerMappedOnDevice( &( vecPatches.B_localx[0]->data_[0] ) );

    int nFieldLocalx = vecPatches.B_localx.size()/2;
    for( int icomp=0 ; icomp<2 ; icomp++ ) {
        if( nFieldLocalx==0 ) {
            continue;
        }

        unsigned int ny_( 1 ), nz_( 1 ), gsp;
        if( nDim>1 ) {
            ny_ = vecPatches.B_localx[icomp*nFieldLocalx]->dims_[1];
            if( nDim>2 ) {
                nz_ = vecPatches.B_localx[icomp*nFieldLocalx]->dims_[2];
            }
        }
        gsp = ( oversize + 1 + vecPatches.B_localx[icomp*nFieldLocalx]->isDual_[0] ); //Ghost size primal

        unsigned int istart =  icomp   *nFieldLocalx;
        unsigned int iend    = ( icomp+1 )*nFieldLocalx;
        #pragma omp for schedule(static) private(pt1,pt2)
        for( unsigned int ifield=istart ; ifield<iend ; ifield++ ) {
            int ipatch = vecPatches.LocalxIdx[ ifield-icomp*nFieldLocalx ];

            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
                pt1 = &( fields[vecPatches( ipatch )->neighbor_[0][0]-h0+icomp*nPatches]->data_[size*ny_*nz_] );
                pt2 = &( vecPatches.B_localx[ifield]->data_[0] );
                //for filter

                if( is_memory_on_device ) {
                    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceMemoryCopy( smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( pt2 ),
                                                                                      smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( pt1 ),
                                                                                      oversize * ny_ * nz_ );

                    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceMemoryCopy( smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( pt1 ) + gsp * ny_ * nz_,
                                                                                      smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( pt2 ) + gsp * ny_ * nz_,
                                                                                      oversize * ny_ * nz_ );
                } else {
                    // If we have GPU support enabled and for some reason we have to handle a CPU buffer,
                    // IsHostPointerMappedOnDevice would prevent us from using the GPU memcpy function.
                    std::memcpy( pt2, pt1, oversize * ny_ * nz_ * sizeof( double ) );
                    std::memcpy( pt1 + gsp * ny_ * nz_, pt2 + gsp * ny_ * nz_, oversize * ny_ * nz_ * sizeof( double ) );
                }
            } // End if ( MPI_me_ == MPI_neighbor_[0][0] )

        } // End for( ipatch )
    }

}

// MPI_Wait for all communications initialised in exchangeAllComponentsAlongX
void SyncVectorPatch::finalizeExchangeAllComponentsAlongX( VectorPatch &vecPatches )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[0];

    unsigned int nMPIx = vecPatches.MPIxIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIx ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIxIdx[ifield];
        vecPatches( ipatch )->finalizeExchange( vecPatches.B_MPIx[ifield      ], 0 ); // By
        vecPatches( ipatch )->finalizeExchange( vecPatches.B_MPIx[ifield+nMPIx], 0 ); // Bz
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, ( iNeighbor+1 )%2 ) ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                Field* field = vecPatches.B_MPIx[ifield      ];
                double* By   = field->recvFields_[(iNeighbor+1)%2]->data_;
                int sizeofBy = field->recvFields_[(iNeighbor+1)%2]->size();
                field        = vecPatches.B_MPIx[ifield+nMPIx];
                double* Bz   = field->recvFields_[(iNeighbor+1)%2]->data_;
                int sizeofBz = field->recvFields_[(iNeighbor+1)%2]->size();
                //#pragma acc update device(By[0:sizeofBy],Bz[0:sizeofBz])
#endif
                vecPatches.B_MPIx[ifield      ]->inject_fields_exch( 0, iNeighbor, oversize );
                vecPatches.B_MPIx[ifield+nMPIx]->inject_fields_exch( 0, iNeighbor, oversize );
            }
        }
    }

}


// The idea is to minimize the number of implicit barriers and maximize the workload between barriers
// fields : contains components of B which required exchange in Y (Bx and Bz)
//     - fields is not directly used in the exchange process, just to find local neighbor's field
//     - sexchanges are operated on sublists, members of vecPatches, which contains :
//         - B_MPIy   : fields which have MPI   neighbor along Y
//         - B_Localy : fields which have local neighbor along Y (a same field can be adressed by both)
//     - These fields are identified with lists of index MPIyIdx and LocalyIdx
void SyncVectorPatch::exchangeAllComponentsAlongY( std::vector<Field *> &fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[1];

    unsigned int nMPIy = vecPatches.MPIyIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIy ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIyIdx[ifield];
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, iNeighbor ) ) {
                vecPatches.B1_MPIy[ifield      ]->create_sub_fields  ( 1, iNeighbor, oversize );
                vecPatches.B1_MPIy[ifield      ]->extract_fields_exch( 1, iNeighbor, oversize );
                vecPatches.B1_MPIy[ifield+nMPIy]->create_sub_fields  ( 1, iNeighbor, oversize );
                vecPatches.B1_MPIy[ifield+nMPIy]->extract_fields_exch( 1, iNeighbor, oversize );
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                Field* field = vecPatches.B1_MPIy[ifield      ];
                double* Bx   = field->sendFields_[iNeighbor+2]->data_;
                int sizeofBx = field->sendFields_[iNeighbor+2]->size();
                field        = vecPatches.B1_MPIy[ifield+nMPIy];
                double* Bz   = field->sendFields_[iNeighbor+2]->data_;
                int sizeofBz = field->sendFields_[iNeighbor+2]->size();
                //#pragma acc update host(Bx[0:sizeofBx],Bz[0:sizeofBz])
#endif
            }
        }
        vecPatches( ipatch )->initExchange( vecPatches.B1_MPIy[ifield      ], 1, smpi, true ); // Bx
        vecPatches( ipatch )->initExchange( vecPatches.B1_MPIy[ifield+nMPIy], 1, smpi, true ); // Bz
    }

    unsigned int h0, size;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    size = vecPatches( 0 )->EMfields->size_[1];

    int nPatches( vecPatches.size() );
    int nDim = vecPatches( 0 )->EMfields->Bx_->dims_.size();

    int nFieldLocaly = vecPatches.B1_localy.size()/2;
    for( int icomp=0 ; icomp<2 ; icomp++ ) {
        if( nFieldLocaly==0 ) {
            continue;
        }

        unsigned int nx_, ny_, nz_( 1 ), gsp;
        nx_ = vecPatches.B1_localy[icomp*nFieldLocaly]->dims_[0];
        ny_ = vecPatches.B1_localy[icomp*nFieldLocaly]->dims_[1];
        if( nDim>2 ) {
            nz_ = vecPatches.B1_localy[icomp*nFieldLocaly]->dims_[2];
        }
        //for filter
        gsp = ( oversize + 1 + vecPatches.B1_localy[icomp*nFieldLocaly]->isDual_[1] ); //Ghost size primal

        unsigned int istart =  icomp   *nFieldLocaly;
        unsigned int iend    = ( icomp+1 )*nFieldLocaly;
        #pragma omp for schedule(static) private(pt1,pt2)
        for( unsigned int ifield=istart ; ifield<iend ; ifield++ ) {

            int ipatch = vecPatches.LocalyIdx[ ifield-icomp*nFieldLocaly ];
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
                pt1 = &( fields[vecPatches( ipatch )->neighbor_[1][0]-h0+icomp*nPatches]->data_[size*nz_] );
                pt2 = &( vecPatches.B1_localy[ifield]->data_[0] );
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                int ptsize = vecPatches.B1_localy[ifield]->size();
                #pragma acc parallel present(pt1[0-size*nz_:ptsize],pt2[0:ptsize])
                #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                const int ptsize = ( nx_ * ny_ * nz_ ) - ( ny_ * nz_ ) + oversize * nz_ + gsp * nz_;
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 2 )
#endif
                for( unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_ ) {
                    // for filter
                    for( unsigned int j = 0 ; j < oversize*nz_ ; j++ ) {
                        pt2[i+j] = pt1[i+j] ;
                        pt1[i+j+gsp*nz_] = pt2[i+j+gsp*nz_] ;
                    } // mempy to do
                }
            } // End if ( MPI_me_ == MPI_neighbor_[1][0] )

        } // End for( ipatch )
    }
}


// MPI_Wait for all communications initialised in exchangeAllComponentsAlongY
void SyncVectorPatch::finalizeExchangeAllComponentsAlongY( VectorPatch &vecPatches )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[1];

    unsigned int nMPIy = vecPatches.MPIyIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIy ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIyIdx[ifield];
        vecPatches( ipatch )->finalizeExchange( vecPatches.B1_MPIy[ifield      ], 1 ); // By
        vecPatches( ipatch )->finalizeExchange( vecPatches.B1_MPIy[ifield+nMPIy], 1 ); // Bz
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, ( iNeighbor+1 )%2 ) ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                Field* field = vecPatches.B1_MPIy[ifield      ];
                double* Bx   = field->recvFields_[(iNeighbor+1)%2+2]->data_;
                int sizeofBx = field->recvFields_[(iNeighbor+1)%2+2]->size();
                field        = vecPatches.B1_MPIy[ifield+nMPIy];
                double* Bz   = field->recvFields_[(iNeighbor+1)%2+2]->data_;
                int sizeofBz = field->recvFields_[(iNeighbor+1)%2+2]->size();
                //#pragma acc update device(Bx[0:sizeofBx],Bz[0:sizeofBz])
#endif
                vecPatches.B1_MPIy[ifield      ]->inject_fields_exch( 1, iNeighbor, oversize );
                vecPatches.B1_MPIy[ifield+nMPIy]->inject_fields_exch( 1, iNeighbor, oversize );
            }
        }
    }

}


// The idea is to minimize the number of implicit barriers and maximize the workload between barriers
// fields : contains components of B which required exchange in Z (Bx and By)
//     - fields is not directly used in the exchange process, just to find local neighbor's field
//     - sexchanges are operated on sublists, members of vecPatches, which contains :
//         - B_MPIz   : fields which have MPI   neighbor along Z
//         - B_Localz : fields which have local neighbor along Z (a same field can be adressed by both)
//     - These fields are identified with lists of index MPIzIdx and LocalzIdx
void SyncVectorPatch::exchangeAllComponentsAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[2];

    unsigned int nMPIz = vecPatches.MPIzIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIz ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIzIdx[ifield];
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, iNeighbor ) ) {
                vecPatches.B2_MPIz[ifield      ]->create_sub_fields  ( 2, iNeighbor, oversize );
                vecPatches.B2_MPIz[ifield      ]->extract_fields_exch( 2, iNeighbor, oversize );
                vecPatches.B2_MPIz[ifield+nMPIz]->create_sub_fields  ( 2, iNeighbor, oversize );
                vecPatches.B2_MPIz[ifield+nMPIz]->extract_fields_exch( 2, iNeighbor, oversize );
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                Field* field = vecPatches.B2_MPIz[ifield      ];
                double* Bx   = field->sendFields_[iNeighbor+4]->data_;
                int sizeofBx = field->sendFields_[iNeighbor+4]->size();
                field        = vecPatches.B2_MPIz[ifield+nMPIz];
                double* By   = field->sendFields_[iNeighbor+4]->data_;
                int sizeofBy = field->sendFields_[iNeighbor+4]->size();
                //#pragma acc update host(Bx[0:sizeofBx],By[0:sizeofBy])
#endif
            }
        }
        vecPatches( ipatch )->initExchange( vecPatches.B2_MPIz[ifield],       2, smpi, true ); // Bx
        vecPatches( ipatch )->initExchange( vecPatches.B2_MPIz[ifield+nMPIz], 2, smpi, true ); // By
    }

    unsigned int h0, size;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    size = vecPatches( 0 )->EMfields->size_[2];

    int nPatches( vecPatches.size() );

    int nFieldLocalz = vecPatches.B2_localz.size()/2;
    for( int icomp=0 ; icomp<2 ; icomp++ ) {
        if( nFieldLocalz==0 ) {
            continue;
        }

        unsigned int nx_, ny_, nz_, gsp;
        nx_ = vecPatches.B2_localz[icomp*nFieldLocalz]->dims_[0];
        ny_ = vecPatches.B2_localz[icomp*nFieldLocalz]->dims_[1];
        nz_ = vecPatches.B2_localz[icomp*nFieldLocalz]->dims_[2];
        //for filter
        gsp = ( oversize + 1 + vecPatches.B2_localz[icomp*nFieldLocalz]->isDual_[2] ); //Ghost size primal

        unsigned int istart  =  icomp   *nFieldLocalz;
        unsigned int iend    = ( icomp+1 )*nFieldLocalz;
        #pragma omp for schedule(static) private(pt1,pt2)
        for( unsigned int ifield=istart ; ifield<iend ; ifield++ ) {

            int ipatch = vecPatches.LocalzIdx[ ifield-icomp*nFieldLocalz ];
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
                pt1 = &( fields[vecPatches( ipatch )->neighbor_[2][0]-h0+icomp*nPatches]->data_[size] );
                pt2 = &( vecPatches.B2_localz[ifield]->data_[0] );
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                int ptsize = vecPatches.B2_localz[ifield]->size();
                #pragma acc parallel present(pt1[0-size:ptsize],pt2[0:ptsize])
                #pragma acc loop gang worker vector
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                const int ptsize = ( nx_ * ny_ * nz_ ) - ( ny_ * nz_ ) + ( ny_ * nz_ ) - nz_ + oversize;
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 3 )
#endif
                for( unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_ ) {
                    for( unsigned int j = 0 ; j < ny_*nz_ ; j += nz_ ) {
                        for( unsigned int k = 0 ; k < oversize ; k++ ) {
                            pt2[i+j+k] = pt1[i+j+k] ;
                            pt1[i+j+k+gsp] = pt2[i+j+k+gsp] ;
                        }
                    }
                }
            } // End if ( MPI_me_ == MPI_neighbor_[2][0] )

        } // End for( ipatch )
    }
}

// MPI_Wait for all communications initialised in exchangeAllComponentsAlongZ
void SyncVectorPatch::finalizeExchangeAllComponentsAlongZ( VectorPatch &vecPatches )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[2];

    unsigned int nMPIz = vecPatches.MPIzIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIz ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIzIdx[ifield];
        vecPatches( ipatch )->finalizeExchange( vecPatches.B2_MPIz[ifield      ], 2 ); // Bx
        vecPatches( ipatch )->finalizeExchange( vecPatches.B2_MPIz[ifield+nMPIz], 2 ); // By
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, ( iNeighbor+1 )%2 ) ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
                Field* field = vecPatches.B2_MPIz[ifield      ];
                double* Bx   = field->recvFields_[(iNeighbor+1)%2+4]->data_;
                int sizeofBx = field->recvFields_[(iNeighbor+1)%2+4]->size();
                field        = vecPatches.B2_MPIz[ifield+nMPIz];
                double* By   = field->recvFields_[(iNeighbor+1)%2+4]->data_;
                int sizeofBy = field->recvFields_[(iNeighbor+1)%2+4]->size();
                //#pragma acc update device(Bx[0:sizeofBx],By[0:sizeofBy])
#endif
                vecPatches.B2_MPIz[ifield      ]->inject_fields_exch( 2, iNeighbor, oversize );
                vecPatches.B2_MPIz[ifield+nMPIz]->inject_fields_exch( 2, iNeighbor, oversize );
            }
        }
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// -------------------------------------------  DEPRECATED FIELD FUNCTIONS  --------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

template<typename T, typename F>
void SyncVectorPatch::exchangeAlongX( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[0];

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        if (fields[ipatch]!=NULL) {
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, iNeighbor ) ) {
                    fields[ipatch]->create_sub_fields  ( 0, iNeighbor, oversize );
                    fields[ipatch]->extract_fields_exch( 0, iNeighbor, oversize );
                }
            }
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange( fields[ipatch], 0, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], 0, smpi );
        }
    }

    unsigned int ny_( 1 ), nz_( 1 ), h0, neighbour_n_space, gsp;
    T *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    // fields[0] can be NULL, look for the 1st non null field
    int IPATCH = 0;
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ )
        if (fields[ipatch]!=NULL) {
            IPATCH = ipatch;
            break;
        }

    if( fields[IPATCH]==NULL ) return;

    gsp = ( oversize + 1 + fields[IPATCH]->isDual_[0] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
            if (fields[ipatch]!=NULL) {
                Field* f = fields[vecPatches( ipatch )->neighbor_[0][0]-h0];
                neighbour_n_space = f->dims_[0] - 1 - 2*oversize - f->isDual_[0]; // define here to be applied on usual or PML fields
                if( fields[IPATCH]->dims_.size()>1 ) {
                    ny_ = f->dims_[1]; // Some can have ny, other ny+npml, but neighboors along X have the same number
                    if( fields[IPATCH]->dims_.size()>2 ) {
                        nz_ = f->dims_[2]; // Some can have nz, other nz+npml, but neighboors along X have the same number
                    }
                }
                pt1 = &( *static_cast<F*>( f              ) )( neighbour_n_space*ny_*nz_ ); //my Xmin neighbour
                pt2 = &( *static_cast<F*>( fields[ipatch] ) )( 0 ); // me

                memcpy( pt2, pt1, oversize*ny_*nz_*sizeof( T ) ); //me <= my neighbour
                memcpy( pt1+gsp*ny_*nz_, pt2+gsp*ny_*nz_, oversize*ny_*nz_*sizeof( T ) );// my neighbour <= me


            } // End if ( MPI_me_ == MPI_neighbor_[0][0] )
        }
    } // End for( ipatch )

}

void SyncVectorPatch::finalizeExchangeAlongX( std::vector<Field *> fields, VectorPatch &vecPatches )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[0];

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        if (fields[ipatch]!=NULL) {
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], 0 );
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 0, ( iNeighbor+1 )%2 ) ) {
                    fields[ipatch]->inject_fields_exch( 0, iNeighbor, oversize );
                }
            }
        }
    }

}

template<typename T, typename F>
void SyncVectorPatch::exchangeAlongY( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[1];

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        if (fields[ipatch]!=NULL) {
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, iNeighbor ) ) {
                    fields[ipatch]->create_sub_fields  ( 1, iNeighbor, oversize );
                    fields[ipatch]->extract_fields_exch( 1, iNeighbor, oversize );
                }
            }
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange( fields[ipatch], 1, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], 1, smpi );
        }
    }

    unsigned int nx_, ny_, nz_( 1 ), h0, gsp, neighbour_n_space, my_ny;
    T *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    // fields[0] can be NULL, look for the 1st non null field
    int IPATCH = 0;
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ )
        if (fields[ipatch]!=NULL) {
            IPATCH = ipatch;
            break;
        }

    if( fields[IPATCH]==NULL ) return;

    gsp = ( oversize + 1 + fields[IPATCH]->isDual_[1] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
            if (fields[ipatch]!=NULL) {
                Field* f = fields[vecPatches( ipatch )->neighbor_[1][0]-h0];
                // In case of PML, ny can be different from e and my neighbour
                nx_ = f->dims_[0]; // Some can have nx, other nx+npml, but neighbours along Y have the same number
                ny_ =                       f->dims_[1]; // Neighbour ny_
                my_ny        = fields[ipatch]->dims_[1];
                neighbour_n_space = ny_ - 1 - 2*oversize - f->isDual_[1]; // define here to be applied on usual or PML fields
                if( f->dims_.size()>2 ) {
                    nz_ = f->dims_[2]; // Some can have nz, other nz+npml, but neighboors along Y have the same number
                }
                pt1 = &( *static_cast<F*>( f              ) )( neighbour_n_space * nz_ );//my Ymin neighbour
                pt2 = &( *static_cast<F*>( fields[ipatch] ) )( 0 );//me
                for( unsigned int i = 0 ; i < nx_ ; i ++ ) {
                    for( unsigned int j = 0 ; j < oversize*nz_ ; j++ ) {
                        pt2[i*my_ny*nz_ + j] = pt1[i*ny_*nz_ + j] ; // me <= my_neighbour
                        pt1[i*ny_*nz_ + j + gsp*nz_] = pt2[i*my_ny*nz_ + j + gsp*nz_] ; //my_neighbour <= me
                    }
                }
            }
        }
    }

}

void SyncVectorPatch::finalizeExchangeAlongY( std::vector<Field *> fields, VectorPatch &vecPatches )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[1];

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        if (fields[ipatch]!=NULL) {
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], 1 );
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 1, ( iNeighbor+1 )%2 ) ) {
                    fields[ipatch]->inject_fields_exch( 1, iNeighbor, oversize );
                }
            }
        }
    }

}

template<typename T, typename F>
void SyncVectorPatch::exchangeAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[2];

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        if (fields[ipatch]!=NULL) {
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, iNeighbor ) ) {
                    fields[ipatch]->create_sub_fields  ( 2, iNeighbor, oversize );
                    fields[ipatch]->extract_fields_exch( 2, iNeighbor, oversize );
                }
            }
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange( fields[ipatch], 2, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], 2, smpi );
        }
    }
    
    unsigned int nx_, ny_, nz_, h0, gsp;
    T *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    // fields[0] can be NULL, look for the 1st non null field
    int IPATCH = 0;
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ )
        if (fields[ipatch]!=NULL) {
            IPATCH = ipatch;
            break;
        }

    if( fields[IPATCH]==NULL ) return;

    gsp = ( oversize + 1 + fields[IPATCH]->isDual_[2] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
            if (fields[ipatch]!=NULL) {
                Field* f = fields[vecPatches( ipatch )->neighbor_[2][0]-h0];
                nx_ = f->dims_[0]; // Some can have nx, other nx+npml, but neighboors along Z have the same number
                ny_ = f->dims_[1]; // Some can have ny, other ny+npml, but neighboors along Z have the same number
                nz_ = f->dims_[2]; // nz is same for all PML
                pt1 = &( *static_cast<F*>( f              ) )( nz_ - 1 - 2*oversize - f->isDual_[2] ); //my neighbour
                pt2 = &( *static_cast<F*>( fields[ipatch] ) )( 0 ); //me
                for( unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_ ) {
                    for( unsigned int j = 0 ; j < ny_*nz_ ; j += nz_ ) {
                        for( unsigned int k = 0 ; k < oversize ; k++ ) {
                            pt2[i+j+k] = pt1[i+j+k] ; //me <= my neighbour
                            pt1[i+j+k+gsp] = pt2[i+j+k+gsp] ; // my_neighbour <= me
                        }
                    }
                }
            }
        }
    }

}

void SyncVectorPatch::finalizeExchangeAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches )
{
    unsigned oversize = vecPatches( 0 )->EMfields->oversize[2];

#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        if (fields[ipatch]!=NULL) {
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], 2 );
            for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
                if ( vecPatches( ipatch )->is_a_MPI_neighbor( 2, ( iNeighbor+1 )%2 ) ) {
                    fields[ipatch]->inject_fields_exch( 2, iNeighbor, oversize );
                }
            }
        }
    }

}

void SyncVectorPatch::exchangeForPML( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    if ( ( params.geometry != "AMcylindrical" ) ) {
        if (params.nDim_field==1) {
            return;
        }
        else {
            // Testing implementation of distributed PML on Xmin and Xmax
            int iDim = 0;
            if ( params.EM_BCs[iDim][0] == "PML" || params.EM_BCs[iDim][1] == "PML" ) { // If a PML along X
                for ( int min_max=0 ; min_max<2 ; min_max++ ) {
                    #pragma omp single
                    vecPatches.buildPMLList( "Bx", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "By", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Bz", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hx", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hy", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hz", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    if( params.Laser_Envelope_model) {
                        if( params.Env_BCs[iDim][min_max] == "PML" /* params.Laser_Envelope_model */ ) {
                            #pragma omp single
                            vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi  );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi  );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi  );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );
                        }
                    }
                }
                if (params.nDim_field>2) {
                    // In 3D, distributed PML on Xmin and Xmax require synchronization along Z
                    for ( int min_max=0 ; min_max<2 ; min_max++ ) {
                        #pragma omp single
                        vecPatches.buildPMLList( "Bx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "By", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Bz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hy", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        if( params.Laser_Envelope_model) {
                            if( params.Env_BCs[iDim][min_max] == "PML" /* params.Laser_Envelope_model */ ) {
                                #pragma omp single
                                vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongZ<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongZ<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongZ<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );
                            }
                        }
                    }
                }
            }
            // Testing implementation of distributed PML on Ymin and Ymax
            iDim = 1;
            if ( params.EM_BCs[iDim][0] == "PML" || params.EM_BCs[iDim][1] == "PML" ) { // If a PML along Y
                for ( int min_max=0 ; min_max<2 ; min_max++ ) {
                    #pragma omp single
                    vecPatches.buildPMLList( "Bx", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "By", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Bz", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hx", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hy", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hz", iDim, min_max, smpi  );

                    SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                    if( params.Laser_Envelope_model) {
                        if( params.Env_BCs[iDim][min_max] == "PML" /* params.Laser_Envelope_model */ ) {
                            #pragma omp single
                            vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi  );

                            SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi  );

                            SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi  );

                            SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        }
                    }
                } 
                if (params.nDim_field>2) {
                    // In 3D, distributed PML on Ymin and Ymax require synchronization along Z
                    for ( int min_max=0 ; min_max<2 ; min_max++ ) {
                        #pragma omp single
                        vecPatches.buildPMLList( "Bx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "By", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Bz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hy", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongZ<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                        if( params.Laser_Envelope_model) {
                            if( params.Env_BCs[iDim][min_max] == "PML" /* params.Laser_Envelope_model */ ) {
                                #pragma omp single
                                vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongZ<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongZ<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongZ<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongZ( vecPatches.listForPML_, vecPatches );
                            }
                        }
                    }
                }
            } 
            if (params.nDim_field>2) {
                int iDim = 2;
                if ( params.EM_BCs[iDim][0] == "PML" || params.EM_BCs[iDim][1] == "PML" ) { // If a PML along Z
                    for ( int min_max=0 ; min_max<2 ; min_max++ ) {
                        #pragma omp single
                        vecPatches.buildPMLList( "Bx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "By", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Bz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hy", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongX<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Bx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "By", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Bz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hx", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hy", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "Hz", iDim, min_max, smpi  );

                        SyncVectorPatch::exchangeAlongY<double,Field>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                        if( params.Laser_Envelope_model) {
                            if( params.Env_BCs[iDim][min_max] == "PML" /* params.Laser_Envelope_model */ ) {
                                #pragma omp single
                                vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                                #pragma omp single
                                vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi  );

                                SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                                SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );
                            }
                        }
                    }
                }
            } // End if (ndim_field>2)
        } // End if (ndim_field>1)
    } // End if (cartesian)
    else { // AM
        for (unsigned int imode = 0 ; imode < params.nmodes ; imode++) {
            // Testing implementation of distributed PML on Xmin and Xmax
            int iDim = 0;
            if ( params.EM_BCs[iDim][0] == "PML" || params.EM_BCs[iDim][1] == "PML" ) { // If a PML along X
                for ( int min_max=0 ; min_max<2 ; min_max++ ) {
                    #pragma omp single
                    vecPatches.buildPMLList( "Bl", iDim, min_max, smpi, imode );

                    SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Br", iDim, min_max, smpi, imode );

                    SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Bt", iDim, min_max, smpi, imode );

                    SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hl", iDim, min_max, smpi, imode );

                    SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Hr", iDim, min_max, smpi, imode );

                    SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    #pragma omp single
                    vecPatches.buildPMLList( "Ht", iDim, min_max, smpi, imode );

                    SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                    SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                    if (imode == 0 && params.Laser_Envelope_model ){
                        if( params.Env_BCs[iDim][min_max] == "PML" /* params.Laser_Envelope_model */ ) {
                            #pragma omp single
                            vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi, imode );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi, imode );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi, imode );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "G_np1_pml", iDim, min_max, smpi, imode );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "G_nm1_pml", iDim, min_max, smpi, imode );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );

                            #pragma omp single
                            vecPatches.buildPMLList( "G_n_pml", iDim, min_max, smpi, imode );

                            SyncVectorPatch::exchangeAlongY<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                            SyncVectorPatch::finalizeExchangeAlongY( vecPatches.listForPML_, vecPatches );
                        }
                    }
                }
            }
            // Testing implementation of distributed PML on Ymin and Ymax
            iDim = 1;
            if (params.EM_BCs[iDim][1] == "PML" ) { // If a PML along R
                int min_max=1 ; //No PML on axis
                #pragma omp single
                vecPatches.buildPMLList( "Bl", iDim, min_max, smpi, imode );

                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                #pragma omp single
                vecPatches.buildPMLList( "Br", iDim, min_max, smpi, imode );

                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                #pragma omp single
                vecPatches.buildPMLList( "Bt", iDim, min_max, smpi, imode );

                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                #pragma omp single
                vecPatches.buildPMLList( "Hl", iDim, min_max, smpi, imode );

                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                #pragma omp single
                vecPatches.buildPMLList( "Hr", iDim, min_max, smpi, imode );

                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                #pragma omp single
                vecPatches.buildPMLList( "Ht", iDim, min_max, smpi, imode );

                SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                if (imode == 0 && params.Laser_Envelope_model ){
                    if( params.Env_BCs[iDim][min_max] == "PML" /* params.Laser_Envelope_model */ ) {
                        #pragma omp single
                        vecPatches.buildPMLList( "A_np1_pml", iDim, min_max, smpi, imode );

                        SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "A_nm1_pml", iDim, min_max, smpi, imode );

                        SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "A_n_pml", iDim, min_max, smpi, imode );

                        SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "G_np1_pml", iDim, min_max, smpi, imode );

                        SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "G_nm1_pml", iDim, min_max, smpi, imode );

                        SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );

                        #pragma omp single
                        vecPatches.buildPMLList( "G_n_pml", iDim, min_max, smpi, imode );

                        SyncVectorPatch::exchangeAlongX<complex<double>,cField>( vecPatches.listForPML_, vecPatches, smpi );
                        SyncVectorPatch::finalizeExchangeAlongX( vecPatches.listForPML_, vecPatches );
                    }
                }
            }
        } // End for( imode )
    } // End else( AM )
}
