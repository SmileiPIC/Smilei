
#include "SyncVectorPatch.h"

#include <vector>

#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

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

void SyncVectorPatch::exchangeParticles( VectorPatch &vecPatches, int ispec, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->initExchParticles( smpi, ispec, params );
    }

    // Init comm in direction 0
#ifndef _NO_MPI_TM
    #pragma omp for schedule(runtime)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->exchNbrOfParticles( smpi, ispec, params, 0, &vecPatches );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! This function performs:
//! - the exhcange of particles for each direction using the diagonal trick.
//! - the importation of the new particles in the particle property arrays
//! - the sorting of particles
// ---------------------------------------------------------------------------------------------------------------------
void SyncVectorPatch::finalizeAndSortParticles( VectorPatch &vecPatches, int ispec, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
    SyncVectorPatch::finalizeExchangeParticles( vecPatches, ispec, 0, params, smpi, timers, itime );

    // Per direction
    for( unsigned int iDim=1 ; iDim<params.nDim_field ; iDim++ ) {
#ifndef _NO_MPI_TM
        #pragma omp for schedule(runtime)
#else
        #pragma omp single
#endif
        for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
            vecPatches( ipatch )->exchNbrOfParticles( smpi, ispec, params, iDim, &vecPatches );
        }

        SyncVectorPatch::finalizeExchangeParticles( vecPatches, ispec, iDim, params, smpi, timers, itime );
    }

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->importAndSortParticles( smpi, ispec, params, &vecPatches );
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


void SyncVectorPatch::finalizeExchangeParticles( VectorPatch &vecPatches, int ispec, int iDim, Params &params, SmileiMPI *smpi, Timers &timers, int itime )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(runtime)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->endNbrOfParticles( smpi, ispec, params, iDim, &vecPatches );
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
        vecPatches( ipatch )->finalizeExchParticles( smpi, ispec, params, iDim, &vecPatches );
    }

    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        vecPatches( ipatch )->cornersParticles( smpi, ispec, params, iDim, &vecPatches );
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------       DENSITIES         ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

void SyncVectorPatch::sumRhoJ( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi, Timers &timers, int itime )
{
    // Sum Jx, Jy and Jz
    SyncVectorPatch::sumAllComponents( vecPatches.densities, vecPatches, smpi, timers, itime );
    // Sum rho
    if( ( vecPatches.diag_flag ) || ( params.is_spectral ) ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listrho_, vecPatches, smpi, timers, itime );
    }
}

void SyncVectorPatch::sumEnvChi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi, Timers &timers, int itime )
{
    // Sum Env_Chi
    SyncVectorPatch::sum<double,Field>( vecPatches.listEnv_Chi_, vecPatches, smpi, timers, itime );
}

void SyncVectorPatch::sumRhoJ( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi, Timers &timers, int itime )
{
    SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJl_[imode], vecPatches, smpi, timers, itime );
    SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJr_[imode], vecPatches, smpi, timers, itime );
    SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJt_[imode], vecPatches, smpi, timers, itime );
    if( ( vecPatches.diag_flag ) || ( params.is_spectral ) ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listrho_AM_[imode], vecPatches, smpi, timers, itime );
    }
}

void SyncVectorPatch::sumRhoJs( Params &params, VectorPatch &vecPatches, int ispec, SmileiMPI *smpi, Timers &timers, int itime )
{
    // Sum Jx_s(ispec), Jy_s(ispec) and Jz_s(ispec)
    if( vecPatches.listJxs_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listJxs_, vecPatches, smpi, timers, itime );
    }
    if( vecPatches.listJys_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listJys_, vecPatches, smpi, timers, itime );
    }
    if( vecPatches.listJzs_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listJzs_, vecPatches, smpi, timers, itime );
    }
    // Sum rho_s(ispec)
    if( vecPatches.listrhos_.size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listrhos_, vecPatches, smpi, timers, itime );
    }
}

void SyncVectorPatch::sumEnvChis( Params &params, VectorPatch &vecPatches, int ispec, SmileiMPI *smpi, Timers &timers, int itime )
{
    // Sum EnvChi_s(ispec)
    if( vecPatches.listEnv_Chis_ .size()>0 ) {
        SyncVectorPatch::sum<double,Field>( vecPatches.listEnv_Chis_, vecPatches, smpi, timers, itime );
    }

}
void SyncVectorPatch::sumRhoJs( Params &params, VectorPatch &vecPatches, int imode, int ispec, SmileiMPI *smpi, Timers &timers, int itime )
{
    // Sum Jx_s(ispec), Jy_s(ispec) and Jz_s(ispec)
    if( vecPatches.listJls_[imode].size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJls_[imode], vecPatches, smpi, timers, itime );
    }
    if( vecPatches.listJrs_[imode].size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJrs_[imode], vecPatches, smpi, timers, itime );
    }
    if( vecPatches.listJts_[imode] .size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listJts_[imode], vecPatches, smpi, timers, itime );
    }
    // Sum rho_s(ispec)
    if( vecPatches.listrhos_AM_[imode].size()>0 ) {
        SyncVectorPatch::sum<complex<double>,cField>( vecPatches.listrhos_AM_[imode], vecPatches, smpi, timers, itime );
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
// timers and itime were here introduced for debugging
void SyncVectorPatch::sumAllComponents( std::vector<Field *> &fields, VectorPatch &vecPatches, SmileiMPI *smpi, Timers &timers, int itime )
{
    unsigned int h0, oversize[3], n_space[3];
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    int nPatches( vecPatches.size() );

    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    n_space[0] = vecPatches( 0 )->EMfields->n_space[0];
    n_space[1] = vecPatches( 0 )->EMfields->n_space[1];
    n_space[2] = vecPatches( 0 )->EMfields->n_space[2];

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
        vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIx[ifield             ], 0, smpi ); // Jx
        vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIx[ifield+  nPatchMPIx], 0, smpi ); // Jy
        vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIx[ifield+2*nPatchMPIx], 0, smpi ); // Jz
    }
    // iDim = 0, local
    int nFieldLocalx = vecPatches.densitiesLocalx.size()/3;
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
                pt1 = &( fields[ vecPatches( ipatch )->neighbor_[0][0]-h0+icomp*nPatches ]->data_[n_space[0]*ny_*nz_] );
                pt2 = &( vecPatches.densitiesLocalx[ifield]->data_[0] );
                //Sum 2 ==> 1
                for( unsigned int i = 0; i < gsp[0]* ny_*nz_ ; i++ ) {
                    pt1[i] += pt2[i];
                }
                //Copy back the results to 2
                memcpy( pt2, pt1, gsp[0]*ny_*nz_*sizeof( double ) );
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
            vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIy[ifield             ], 1, smpi ); // Jx
            vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIy[ifield+nPatchMPIy  ], 1, smpi ); // Jy
            vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIy[ifield+2*nPatchMPIy], 1, smpi ); // Jz
        }

        // iDim = 1,
        int nFieldLocaly = vecPatches.densitiesLocaly.size()/3;
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
                    pt1 = &( fields[vecPatches( ipatch )->neighbor_[1][0]-h0+icomp*nPatches]->data_[n_space[1]*nz_] );
                    pt2 = &( vecPatches.densitiesLocaly[ifield]->data_[0] );
                    for( unsigned int j = 0; j < nx_ ; j++ ) {
                        for( unsigned int i = 0; i < gsp[1]*nz_ ; i++ ) {
                            pt1[i] += pt2[i];
                        }
                        memcpy( pt2, pt1, gsp[1]*nz_*sizeof( double ) );
                        pt1 += ny_*nz_;
                        pt2 += ny_*nz_;
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
                vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIz[ifield             ], 2, smpi ); // Jx
                vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIz[ifield+nPatchMPIz  ], 2, smpi ); // Jy
                vecPatches( ipatch )->initSumField( vecPatches.densitiesMPIz[ifield+2*nPatchMPIz], 2, smpi ); // Jz
            }

            // iDim = 2 local
            int nFieldLocalz = vecPatches.densitiesLocalz.size()/3;
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
                        pt1 = &( fields[vecPatches( ipatch )->neighbor_[2][0]-h0+icomp*nPatches]->data_[n_space[2]] );
                        pt2 = &( vecPatches.densitiesLocalz[ifield]->data_[0] );
                        for( unsigned int j = 0; j < nx_*ny_ ; j++ ) {
                            for( unsigned int i = 0; i < gsp[2] ; i++ ) {
                                pt1[i] += pt2[i];
                                pt2[i] =  pt1[i];
                            }
                            pt1 += nz_;
                            pt2 += nz_;
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
        SyncVectorPatch::finalizeExchangeAllComponentsAlongX( vecPatches.Bs0, vecPatches );
    } else if( vecPatches.listBx_[0]->dims_.size()==2 ) {
        if( !params.full_B_exchange ) {
            // Finalize exchange Bs0 : By_ and Bz_ (dual in X)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongX( vecPatches.Bs0, vecPatches );
            // Finalize exchange Bs1 : Bx_ and Bz_ (dual in Y)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongY( vecPatches.Bs1, vecPatches );
        }
        //else
        //    done in exchangeSynchronizedPerDirection
    } else if( vecPatches.listBx_[0]->dims_.size()==3 ) {
        if( !params.full_B_exchange ) {
            // Finalize exchange Bs0 : By_ and Bz_ (dual in X)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongX( vecPatches.Bs0, vecPatches );
            // Finalize exchange Bs1 : Bx_ and Bz_ (dual in Y)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongY( vecPatches.Bs1, vecPatches );
            // Finalize exchange Bs2 : Bx_ and By_ (dual in Z)
            SyncVectorPatch::finalizeExchangeAllComponentsAlongZ( vecPatches.Bs2, vecPatches );
        }
        //else
        //    done in exchangeSynchronizedPerDirection
    }

}

void SyncVectorPatch::exchangeJ( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{

    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listJx_, vecPatches, smpi );
    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listJy_, vecPatches, smpi );
    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listJz_, vecPatches, smpi );
}

void SyncVectorPatch::finalizeexchangeJ( Params &params, VectorPatch &vecPatches )
{

    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listJx_, vecPatches );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listJy_, vecPatches );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listJz_, vecPatches );
}


void SyncVectorPatch::exchangeB( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
{
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBl_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBl_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBr_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBr_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listBt_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listBt_[imode], vecPatches );
}

void SyncVectorPatch::exchangeE( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
{
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listEl_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEl_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listEr_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEr_[imode], vecPatches );
    SyncVectorPatch::exchangeAlongAllDirections<complex<double>,cField>( vecPatches.listEt_[imode], vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEt_[imode], vecPatches );
}

void SyncVectorPatch::finalizeexchangeB( Params &params, VectorPatch &vecPatches, int imode )
{
}

void SyncVectorPatch::finalizeexchangeE( Params &params, VectorPatch &vecPatches, int imode )
{
}


//void SyncVectorPatch::exchangeB( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi )
//{
//
//    SyncVectorPatch::exchangeComplex( vecPatches.listBl_[imode], vecPatches, smpi );
//    SyncVectorPatch::exchangeComplex( vecPatches.listBr_[imode], vecPatches, smpi );
//    SyncVectorPatch::exchangeComplex( vecPatches.listBt_[imode], vecPatches, smpi );
//}

//void SyncVectorPatch::finalizeexchangeB( Params &params, VectorPatch &vecPatches, int imode )
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

void SyncVectorPatch::finalizeexchangeA( Params &params, VectorPatch &vecPatches )
{
//    // current envelope value
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listA_, vecPatches );
//    // current envelope value
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listA0_, vecPatches );
}

// void SyncVectorPatch::exchangeEnvEEnvA( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
// {
//     // current envelope |E| value
//     SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEnvE_, vecPatches, smpi );
//     SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEnvE_, vecPatches );
//     // current envelope |A| value
//     SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEnvA_, vecPatches, smpi );
//     SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEnvA_, vecPatches );
// }
// 
// void SyncVectorPatch::finalizeexchangeEnvEEnvA( Params &params, VectorPatch &vecPatches )
// {
// 
// }

void SyncVectorPatch::exchangeEnvEx( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    // current envelope |Ex| value
    SyncVectorPatch::exchangeAlongAllDirections<double,Field>( vecPatches.listEnvEx_, vecPatches, smpi );
    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listEnvEx_, vecPatches );
}

void SyncVectorPatch::finalizeexchangeEnvEx( Params &params, VectorPatch &vecPatches )
{

}

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
// void SyncVectorPatch::finalizeexchangePhi( Params &params, VectorPatch &vecPatches )
// {
// 
// }




void SyncVectorPatch::exchangeGradPhi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    if (  params.geometry != "AMcylindrical" ) {
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
    }
}

void SyncVectorPatch::finalizeexchangeGradPhi( Params &params, VectorPatch &vecPatches )
{
    // current Gradient value
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhix_, vecPatches );
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiy_, vecPatches );
//    SyncVectorPatch::finalizeExchangeAlongAllDirections( vecPatches.listGradPhiz_, vecPatches );
}

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
// timers and itime were here introduced for debugging
template<typename T, typename F>
void SyncVectorPatch::exchangeAlongAllDirections( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange       ( fields[ipatch], iDim, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], iDim, smpi );
        }
    } // End for iDim


    unsigned int nx_, ny_( 1 ), nz_( 1 ), h0, oversize[3], n_space[3], gsp[3];
    T *pt1, *pt2;
    F *field1, *field2;
    h0 = vecPatches( 0 )->hindex;

    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    n_space[0] = vecPatches( 0 )->EMfields->n_space[0];
    n_space[1] = vecPatches( 0 )->EMfields->n_space[1];
    n_space[2] = vecPatches( 0 )->EMfields->n_space[2];

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
            pt1 = &( *field1 )( ( n_space[0] )*ny_*nz_ );
            pt2 = &( *field2 )( 0 );
            memcpy( pt2, pt1, oversize[0]*ny_*nz_*sizeof( T ) );
            memcpy( pt1+gsp[0]*ny_*nz_, pt2+gsp[0]*ny_*nz_, oversize[0]*ny_*nz_*sizeof( T ) );
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] )

        if( fields[0]->dims_.size()>1 ) {
            gsp[1] = ( oversize[1] + 1 + fields[0]->isDual_[1] ); //Ghost size primal
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
                field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[1][0]-h0] );
                field2 = static_cast<F *>( fields[ipatch] );
                pt1 = &( *field1 )( n_space[1]*nz_ );
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
                    pt1 = &( *field1 )( n_space[2] );
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
    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
#ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
#else
        #pragma omp single
#endif
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->finalizeExchange       ( fields[ipatch], iDim );
            else
                vecPatches( ipatch )->finalizeExchangeComplex( fields[ipatch], iDim );
        }
    } // End for iDim

}


// fields : contains a single field component (X, Y or Z) for all patches of vecPatches
// timers and itime were here introduced for debugging
template<typename T, typename F>
void SyncVectorPatch::exchangeAlongAllDirectionsNoOMP( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initExchange       ( fields[ipatch], iDim, smpi );
            else
                vecPatches( ipatch )->initExchangeComplex( fields[ipatch], iDim, smpi );
        }
    } // End for iDim


    unsigned int nx_, ny_( 1 ), nz_( 1 ), h0, oversize[3], n_space[3], gsp[3];
    T *pt1, *pt2;
    F* field1;
    F* field2;
    h0 = vecPatches( 0 )->hindex;

    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    n_space[0] = vecPatches( 0 )->EMfields->n_space[0];
    n_space[1] = vecPatches( 0 )->EMfields->n_space[1];
    n_space[2] = vecPatches( 0 )->EMfields->n_space[2];

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
            pt1 = &( *field1 )( n_space[0]*ny_*nz_ );
            pt2 = &( *field2 )( 0 );
            memcpy( pt2, pt1, oversize[0]*ny_*nz_*sizeof( T ) );
            memcpy( pt1+gsp[0]*ny_*nz_, pt2+gsp[0]*ny_*nz_, oversize[0]*ny_*nz_*sizeof( T ) );
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] )

        if( fields[0]->dims_.size()>1 ) {
            gsp[1] = ( oversize[1] + 1 + fields[0]->isDual_[1] ); //Ghost size primal
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
                field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[1][0]-h0] );
                field2 = static_cast<F *>( fields[ipatch] );
                pt1 = &( *field1 )( n_space[1]*nz_ );
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
                    pt1 = &( *field1 )( n_space[2] );
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
    for( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->finalizeExchange       ( fields[ipatch], iDim );
            else
                vecPatches( ipatch )->finalizeExchangeComplex( fields[ipatch], iDim );
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

    unsigned int nx_, ny_( 1 ), nz_( 1 ), h0, oversize[3], n_space[3], gsp[3];
    T *pt1, *pt2;
    F* field1;
    F* field2;
    h0 = vecPatches( 0 )->hindex;

    oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
    oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
    oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

    n_space[0] = vecPatches( 0 )->EMfields->n_space[0];
    n_space[1] = vecPatches( 0 )->EMfields->n_space[1];
    n_space[2] = vecPatches( 0 )->EMfields->n_space[2];

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
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->finalizeExchange( fields[ipatch], 2 );
            else
                vecPatches( ipatch )->finalizeExchangeComplex( fields[ipatch], 2 );
        }

        #pragma omp for schedule(static) private(pt1,pt2)
        for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

            gsp[2] = ( oversize[2] + 1 + fields[0]->isDual_[2] ); //Ghost size primal
            if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
                field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[2][0]-h0]  );
                field2 = static_cast<F *>( fields[ipatch] );
                pt1 = &( *field1 )( n_space[2] );
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
        if ( !dynamic_cast<cField*>( fields[ipatch] ) )
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], 1 );
        else
            vecPatches( ipatch )->finalizeExchangeComplex( fields[ipatch], 1 );
    }

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        gsp[1] = ( oversize[1] + 1 + fields[0]->isDual_[1] ); //Ghost size primal
        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
            field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[1][0]-h0]  );
            field2 = static_cast<F *>( fields[ipatch] );
            pt1 = &( *field1 )( n_space[1]*nz_ );
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
        if ( !dynamic_cast<cField*>( fields[ipatch] ) )
            vecPatches( ipatch )->finalizeExchange( fields[ipatch], 0 );
        else
            vecPatches( ipatch )->finalizeExchangeComplex( fields[ipatch], 0 );
    }



    gsp[0] = ( oversize[0] + 1 + fields[0]->isDual_[0] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
            field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[0][0]-h0]  );
            field2 = static_cast<F *>( fields[ipatch] );
            pt1 = &( *field1 )( ( n_space[0] )*ny_*nz_ );
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
    unsigned int nMPIx = vecPatches.MPIxIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIx ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIxIdx[ifield];
        vecPatches( ipatch )->initExchange( vecPatches.B_MPIx[ifield      ], 0, smpi ); // By
        vecPatches( ipatch )->initExchange( vecPatches.B_MPIx[ifield+nMPIx], 0, smpi ); // Bz
    }


    unsigned int h0, oversize, n_space;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    oversize = vecPatches( 0 )->EMfields->oversize[0];

    n_space = vecPatches( 0 )->EMfields->n_space[0];

    int nPatches( vecPatches.size() );
    int nDim = vecPatches( 0 )->EMfields->Bx_->dims_.size();

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
                pt1 = &( fields[vecPatches( ipatch )->neighbor_[0][0]-h0+icomp*nPatches]->data_[n_space*ny_*nz_] );
                pt2 = &( vecPatches.B_localx[ifield]->data_[0] );
                //for filter
                memcpy( pt2, pt1, oversize*ny_*nz_*sizeof( double ) );
                memcpy( pt1+gsp*ny_*nz_, pt2+gsp*ny_*nz_, oversize*ny_*nz_*sizeof( double ) );
            } // End if ( MPI_me_ == MPI_neighbor_[0][0] )

        } // End for( ipatch )
    }

}

// MPI_Wait for all communications initialised in exchangeAllComponentsAlongX
void SyncVectorPatch::finalizeExchangeAllComponentsAlongX( std::vector<Field *> &fields, VectorPatch &vecPatches )
{
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
    unsigned int nMPIy = vecPatches.MPIyIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIy ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIyIdx[ifield];
        vecPatches( ipatch )->initExchange( vecPatches.B1_MPIy[ifield      ], 1, smpi ); // Bx
        vecPatches( ipatch )->initExchange( vecPatches.B1_MPIy[ifield+nMPIy], 1, smpi ); // Bz
    }

    unsigned int h0, oversize, n_space;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    oversize = vecPatches( 0 )->EMfields->oversize[1];
    n_space = vecPatches( 0 )->EMfields->n_space[1];

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
                pt1 = &( fields[vecPatches( ipatch )->neighbor_[1][0]-h0+icomp*nPatches]->data_[n_space*nz_] );
                pt2 = &( vecPatches.B1_localy[ifield]->data_[0] );
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
void SyncVectorPatch::finalizeExchangeAllComponentsAlongY( std::vector<Field *> &fields, VectorPatch &vecPatches )
{
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
    unsigned int nMPIz = vecPatches.MPIzIdx.size();
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ifield=0 ; ifield<nMPIz ; ifield++ ) {
        unsigned int ipatch = vecPatches.MPIzIdx[ifield];
        vecPatches( ipatch )->initExchange( vecPatches.B2_MPIz[ifield],       2, smpi ); // Bx
        vecPatches( ipatch )->initExchange( vecPatches.B2_MPIz[ifield+nMPIz], 2, smpi ); // By
    }

    unsigned int h0, oversize, n_space;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    oversize = vecPatches( 0 )->EMfields->oversize[2];
    n_space = vecPatches( 0 )->EMfields->n_space[2];

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
                pt1 = &( fields[vecPatches( ipatch )->neighbor_[2][0]-h0+icomp*nPatches]->data_[n_space] );
                pt2 = &( vecPatches.B2_localz[ifield]->data_[0] );
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
void SyncVectorPatch::finalizeExchangeAllComponentsAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches )
{
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
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// -------------------------------------------  DEPRECATED FIELD FUNCTIONS  --------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

void SyncVectorPatch::exchangeAlongX( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->initExchange( fields[ipatch], 0, smpi );
    }

    unsigned int ny_( 1 ), nz_( 1 ), h0, oversize, n_space, gsp;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    oversize = vecPatches( 0 )->EMfields->oversize[0];

    n_space = vecPatches( 0 )->EMfields->n_space[0];

    if( fields[0]->dims_.size()>1 ) {
        ny_ = fields[0]->dims_[1];
        if( fields[0]->dims_.size()>2 ) {
            nz_ = fields[0]->dims_[2];
        }
    }

    //gsp[0] = 2*oversize[0]+fields[0]->isDual_[0]; //Ghost size primal
    //for filter
    gsp = ( oversize + 1 + fields[0]->isDual_[0] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
            pt1 = &( *fields[vecPatches( ipatch )->neighbor_[0][0]-h0] )( n_space*ny_*nz_ );
            pt2 = &( *fields[ipatch] )( 0 );
            //memcpy( pt2, pt1, ny_*sizeof(double));
            //memcpy( pt1+gsp[0]*ny_, pt2+gsp[0]*ny_, ny_*sizeof(double));
            //for filter
            memcpy( pt2, pt1, oversize*ny_*nz_*sizeof( double ) );
            memcpy( pt1+gsp*ny_*nz_, pt2+gsp*ny_*nz_, oversize*ny_*nz_*sizeof( double ) );
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] )


    } // End for( ipatch )

}

void SyncVectorPatch::finalizeExchangeAlongX( std::vector<Field *> fields, VectorPatch &vecPatches )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->finalizeExchange( fields[ipatch], 0 );
    }

}

void SyncVectorPatch::exchangeAlongY( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->initExchange( fields[ipatch], 1, smpi );
    }

    unsigned int nx_, ny_, nz_( 1 ), h0, oversize, n_space, gsp;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    oversize = vecPatches( 0 )->EMfields->oversize[1];
    n_space = vecPatches( 0 )->EMfields->n_space[1];

    nx_ = fields[0]->dims_[0];
    ny_ = fields[0]->dims_[1];
    if( fields[0]->dims_.size()>2 ) {
        nz_ = fields[0]->dims_[2];
    }

    //gsp = 2*oversize[1]+fields[0]->isDual_[1]; //Ghost size primal
    //for filter
    gsp = ( oversize + 1 + fields[0]->isDual_[1] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
            pt1 = &( *fields[vecPatches( ipatch )->neighbor_[1][0]-h0] )( n_space*nz_ );
            pt2 = &( *fields[ipatch] )( 0 );
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

void SyncVectorPatch::finalizeExchangeAlongY( std::vector<Field *> fields, VectorPatch &vecPatches )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->finalizeExchange( fields[ipatch], 1 );
    }

}

void SyncVectorPatch::exchangeAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->initExchange( fields[ipatch], 2, smpi );
    }

    unsigned int nx_, ny_, nz_, h0, oversize, n_space, gsp;
    double *pt1, *pt2;
    h0 = vecPatches( 0 )->hindex;

    oversize = vecPatches( 0 )->EMfields->oversize[2];
    n_space = vecPatches( 0 )->EMfields->n_space[2];

    nx_ = fields[0]->dims_[0];
    ny_ = fields[0]->dims_[1];
    nz_ = fields[0]->dims_[2];

    //gsp = 2*oversize[1]+fields[0]->isDual_[1]; //Ghost size primal
    //for filter
    gsp = ( oversize + 1 + fields[0]->isDual_[2] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {

        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
            pt1 = &( *fields[vecPatches( ipatch )->neighbor_[2][0]-h0] )( n_space );
            pt2 = &( *fields[ipatch] )( 0 );
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

void SyncVectorPatch::finalizeExchangeAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches )
{
#ifndef _NO_MPI_TM
    #pragma omp for schedule(static)
#else
    #pragma omp single
#endif
    for( unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++ ) {
        vecPatches( ipatch )->finalizeExchange( fields[ipatch], 2 );
    }

}
