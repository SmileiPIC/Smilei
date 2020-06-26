
#ifndef SYNCVECTORPATCH_H
#define SYNCVECTORPATCH_H

#include <vector>

#include "VectorPatch.h"

class Params;
class SmileiMPI;
class Field;
class cField;
class Timers;

class SyncVectorPatch
{
public :

    //! Particles synchronization
    static void exchangeParticles( VectorPatch &vecPatches, int ispec, Params &params, SmileiMPI *smpi, Timers &timers, int itime );
    static void finalizeAndSortParticles( VectorPatch &vecPatches, int ispec, Params &params, SmileiMPI *smpi, Timers &timers, int itime );
    static void finalizeExchangeParticles( VectorPatch &vecPatches, int ispec, int iDim, Params &params, SmileiMPI *smpi, Timers &timers, int itime );

    //! Densities synchronization
    static void sumRhoJ( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi, Timers &timers, int itime );
    //! Densities synchronization per mode
    static void sumRhoJ( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi, Timers &timers, int itime );
    //! Densities synchronization per species
    static void sumRhoJs( Params &params, VectorPatch &vecPatches, int ispec, SmileiMPI *smpi, Timers &timers, int itime );
    //! Densities synchronization per species per mode
    static void sumRhoJs( Params &params, VectorPatch &vecPatches, int imode, int ispec, SmileiMPI *smpi, Timers &timers, int itime );
    //! Densities synchronization, including envelope
    static void sumEnvChi( Params &params, VectorPatch &vecPatches, SmileiMPI *smp, Timers &timers, int itime );
    static void sumEnvChis( Params &params, VectorPatch &vecPatches, int ispec, SmileiMPI *smp, Timers &timers, int itime );

    // fields : contains a single field component for all patches of vecPatches
    // timers and itime were here introduced for debugging
    template<typename T, typename F> static
    void sum( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi, Timers &timers, int itime )
    {
        unsigned int nx_, ny_, nz_, h0, oversize[3], n_space[3], gsp[3];
        T *pt1, *pt2;
        F* field1;
        F* field2;
        h0 = vecPatches( 0 )->hindex;

        int nPatches( vecPatches.size() );

        oversize[0] = vecPatches( 0 )->EMfields->oversize[0];
        oversize[1] = vecPatches( 0 )->EMfields->oversize[1];
        oversize[2] = vecPatches( 0 )->EMfields->oversize[2];

        n_space[0] = vecPatches( 0 )->EMfields->n_space[0];
        n_space[1] = vecPatches( 0 )->EMfields->n_space[1];
        n_space[2] = vecPatches( 0 )->EMfields->n_space[2];

        unsigned int nComp = fields.size()/nPatches;

        // -----------------
        // Sum per direction :

        // iDim = 0, initialize comms : Isend/Irecv
    #ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
    #else
        #pragma omp single
    #endif
        for( unsigned int ifield=0 ; ifield<fields.size() ; ifield++ ) {
            unsigned int ipatch = ifield%nPatches;
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->initSumField( fields[ifield], 0, smpi );
            else
                vecPatches( ipatch )->initSumFieldComplex( fields[ifield], 0, smpi );
        }

        // iDim = 0, local
        for( unsigned int icomp=0 ; icomp<nComp ; icomp++ ) {
            nx_ = fields[icomp*nPatches]->dims_[0];
            ny_ = 1;
            nz_ = 1;
            if( fields[icomp*nPatches]->dims_.size()>1 ) {
                ny_ = fields[icomp*nPatches]->dims_[1];
                if( fields[icomp*nPatches]->dims_.size()>2 ) {
                    nz_ = fields[icomp*nPatches]->dims_[2];
                }
            }
            gsp[0] = 1+2*oversize[0]+fields[icomp*nPatches]->isDual_[0]; //Ghost size primal
            #pragma omp for schedule(static) private(pt1,pt2)
            for( unsigned int ifield=icomp*nPatches ; ifield<( icomp+1 )*nPatches ; ifield++ ) {
                unsigned int ipatch = ifield%nPatches;
                if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[0][0] ) {
                    //The patch to the west belongs to the same MPI process than I.
                    field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[0][0]-h0+icomp*nPatches] );
                    field2 = static_cast<F *>( fields[ifield] );
                    pt1 = &( *field1 )( n_space[0]*ny_*nz_ );
                    pt2 = &( *field2 )( 0 );
                    //Sum 2 ==> 1
                    for( unsigned int i = 0; i < gsp[0]* ny_*nz_ ; i++ ) {
                        pt1[i] += pt2[i];
                    }
                    //Copy back the results to 2
                    memcpy( pt2, pt1, gsp[0]*ny_*nz_*sizeof( T ) );
                }
            }
        }

        // iDim = 0, finalize (waitall)
    #ifndef _NO_MPI_TM
        #pragma omp for schedule(static)
    #else
        #pragma omp single
    #endif
        for( unsigned int ifield=0 ; ifield<fields.size() ; ifield++ ) {
            unsigned int ipatch = ifield%nPatches;
            if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                vecPatches( ipatch )->finalizeSumField( fields[ifield], 0 );
            else
                vecPatches( ipatch )->finalizeSumFieldComplex( fields[ifield], 0 );
        }
        // END iDim = 0 sync
        // -----------------

        if( fields[0]->dims_.size()>1 ) {
            // -----------------
            // Sum per direction :

            // iDim = 1, initialize comms : Isend/Irecv
    #ifndef _NO_MPI_TM
            #pragma omp for schedule(static)
    #else
            #pragma omp single
    #endif
            for( unsigned int ifield=0 ; ifield<fields.size() ; ifield++ ) {
                unsigned int ipatch = ifield%nPatches;
                if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                    vecPatches( ipatch )->initSumField( fields[ifield], 1, smpi );
                else
                    vecPatches( ipatch )->initSumFieldComplex( fields[ifield], 1, smpi );
            }

            // iDim = 1, local
            for( unsigned int icomp=0 ; icomp<nComp ; icomp++ ) {
                nx_ = fields[icomp*nPatches]->dims_[0];
                ny_ = 1;
                nz_ = 1;
                if( fields[icomp*nPatches]->dims_.size()>1 ) {
                    ny_ = fields[icomp*nPatches]->dims_[1];
                    if( fields[icomp*nPatches]->dims_.size()>2 ) {
                        nz_ = fields[icomp*nPatches]->dims_[2];
                    }
                }
                gsp[0] = 1+2*oversize[0]+fields[icomp*nPatches]->isDual_[0]; //Ghost size primal
                gsp[1] = 1+2*oversize[1]+fields[icomp*nPatches]->isDual_[1]; //Ghost size primal

                #pragma omp for schedule(static) private(pt1,pt2)
                for( unsigned int ifield=icomp*nPatches ; ifield<( icomp+1 )*nPatches ; ifield++ ) {
                    unsigned int ipatch = ifield%nPatches;
                    if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[1][0] ) {
                        //The patch to the south belongs to the same MPI process than I.
                        field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[1][0]-h0+icomp*nPatches] );
                        field2 = static_cast<F *>( fields[ifield] );
                        pt1 = &( *field1 )( n_space[1]*nz_ );
                        pt2 = &( *field2 )( 0 );
                        for( unsigned int j = 0; j < nx_ ; j++ ) {
                            for( unsigned int i = 0; i < gsp[1]*nz_ ; i++ ) {
                                pt1[i] += pt2[i];
                            }
                            memcpy( pt2, pt1, gsp[1]*nz_*sizeof( T ) );
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
            for( unsigned int ifield=0 ; ifield<fields.size() ; ifield++ ) {
                unsigned int ipatch = ifield%nPatches;
                if ( !dynamic_cast<cField*>( fields[ipatch] ) )
                    vecPatches( ipatch )->finalizeSumField( fields[ifield], 1 );
                else
                    vecPatches( ipatch )->finalizeSumFieldComplex( fields[ifield], 1 );
            }
            // END iDim = 1 sync
            // -----------------

            if( fields[0]->dims_.size()>2 ) {
                // -----------------
                // Sum per direction :

                // iDim = 2, initialize comms : Isend/Irecv
    #ifndef _NO_MPI_TM
                #pragma omp for schedule(static)
    #else
                #pragma omp single
    #endif
                for( unsigned int ifield=0 ; ifield<fields.size() ; ifield++ ) {
                    unsigned int ipatch = ifield%nPatches;
                    vecPatches( ipatch )->initSumField( fields[ifield], 2, smpi );
                }

                // iDim = 2 local
                for( unsigned int icomp=0 ; icomp<nComp ; icomp++ ) {
                    nx_ = fields[icomp*nPatches]->dims_[0];
                    ny_ = 1;
                    nz_ = 1;
                    if( fields[icomp*nPatches]->dims_.size()>1 ) {
                        ny_ = fields[icomp*nPatches]->dims_[1];
                        if( fields[icomp*nPatches]->dims_.size()>2 ) {
                            nz_ = fields[icomp*nPatches]->dims_[2];
                        }
                    }
                    gsp[0] = 1+2*oversize[0]+fields[icomp*nPatches]->isDual_[0]; //Ghost size primal
                    gsp[1] = 1+2*oversize[1]+fields[icomp*nPatches]->isDual_[1]; //Ghost size primal
                    gsp[2] = 1+2*oversize[2]+fields[icomp*nPatches]->isDual_[2]; //Ghost size primal
                    #pragma omp for schedule(static) private(pt1,pt2)
                    for( unsigned int ifield=icomp*nPatches ; ifield<( icomp+1 )*nPatches ; ifield++ ) {
                        unsigned int ipatch = ifield%nPatches;
                        if( vecPatches( ipatch )->MPI_me_ == vecPatches( ipatch )->MPI_neighbor_[2][0] ) {
                            //The patch below me belongs to the same MPI process than I.
                            field1 = static_cast<F *>( fields[vecPatches( ipatch )->neighbor_[2][0]-h0+icomp*nPatches] );
                            field2 = static_cast<F *>( fields[ifield] );
                            pt1 = &( *field1 )( n_space[2] );
                            pt2 = &( *field2 )( 0 );
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
                for( unsigned int ifield=0 ; ifield<fields.size() ; ifield++ ) {
                    unsigned int ipatch = ifield%nPatches;
                    vecPatches( ipatch )->finalizeSumField( fields[ifield], 2 );
                }
                // END iDim = 2 sync
                // -----------------

            } // End if dims_.size()>2

        } // End if dims_.size()>1

    }

    static void sumAllComponents( std::vector<Field *> &fields, VectorPatch &vecPatches, SmileiMPI *smpi, Timers &timers, int itime );

    void templateGenerator();

    //! Fields synchronization
    static void exchangeE( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeexchangeE( Params &params, VectorPatch &vecPatches );
    static void exchangeB( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeexchangeB( Params &params, VectorPatch &vecPatches );

    static void exchangeE( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi );
    static void finalizeexchangeE( Params &params, VectorPatch &vecPatches, int imode );   
    static void exchangeB( Params &params, VectorPatch &vecPatches, int imode, SmileiMPI *smpi );
    static void finalizeexchangeB( Params &params, VectorPatch &vecPatches, int imode );

    static void exchangeJ( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeexchangeJ( Params &params, VectorPatch &vecPatches );

    static void exchangeA( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeexchangeA( Params &params, VectorPatch &vecPatches );
    // static void exchangeEnvEEnvA( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    // static void finalizeexchangeEnvEEnvA( Params &params, VectorPatch &vecPatches );
    // static void exchangePhi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    // static void finalizeexchangePhi( Params &params, VectorPatch &vecPatches );
    static void exchangeEnvEx( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeexchangeEnvEx( Params &params, VectorPatch &vecPatches );
    static void exchangeGradPhi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeexchangeGradPhi( Params &params, VectorPatch &vecPatches );
    static void exchangeEnvChi( Params &params, VectorPatch &vecPatches, SmileiMPI *smpi );

    template<typename T, typename MT> static void exchangeAlongAllDirections( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAlongAllDirections( std::vector<Field *> fields, VectorPatch &vecPatches );

    template<typename T, typename MT> static void exchangeAlongAllDirectionsNoOMP( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAlongAllDirectionsNoOMP( std::vector<Field *> fields, VectorPatch &vecPatches );

    template<typename T, typename MT> static void exchangeSynchronizedPerDirection( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void exchangeSynchronizedPerDirection( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );

    static void exchangeAllComponentsAlongX( std::vector<Field *> &fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAllComponentsAlongX( std::vector<Field *> &fields, VectorPatch &vecPatches );
    static void exchangeAllComponentsAlongY( std::vector<Field *> &fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAllComponentsAlongY( std::vector<Field *> &fields, VectorPatch &vecPatches );
    static void exchangeAllComponentsAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAllComponentsAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches );

    //! Deprecated field functions
    static void exchangeAlongX( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAlongX( std::vector<Field *> fields, VectorPatch &vecPatches );
    static void exchangeAlongY( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAlongY( std::vector<Field *> fields, VectorPatch &vecPatches );
    static void exchangeAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches, SmileiMPI *smpi );
    static void finalizeExchangeAlongZ( std::vector<Field *> fields, VectorPatch &vecPatches );

};

#endif
