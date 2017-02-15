
#include "SyncVectorPatch.h"

#include <vector>

#include "VectorPatch.h"
#include "Params.h"
#include "SmileiMPI.h"

using namespace std;

void SyncVectorPatch::exchangeParticles(VectorPatch& vecPatches, int ispec, Params &params, SmileiMPI* smpi, Timers &timers, int itime)
{
    //timers.syncPart.restart();
    #pragma omp for schedule(runtime)
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->initExchParticles(smpi, ispec, params);
    }

    // Per direction
    for (unsigned int iDim=0 ; iDim<1 ; iDim++) {
        #pragma omp for schedule(runtime)
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            vecPatches(ipatch)->initCommParticles(smpi, ispec, params, iDim, &vecPatches);
        }

//        #pragma omp for schedule(runtime)
//        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
//            vecPatches(ipatch)->CommParticles(smpi, ispec, params, iDim, &vecPatches);
//        }
//        #pragma omp for schedule(runtime)
//        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
//            vecPatches(ipatch)->finalizeCommParticles(smpi, ispec, params, iDim, &vecPatches);
//        }
    }

//    #pragma omp for schedule(runtime)
//    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
//        vecPatches(ipatch)->vecSpecies[ispec]->sort_part();
    //timers.syncPart.update( params.printNow( itime ) );
}


void SyncVectorPatch::finalize_and_sort_parts(VectorPatch& vecPatches, int ispec, Params &params, SmileiMPI* smpi, Timers &timers, int itime)
{
    //timers.syncPart.restart();
    #pragma omp for schedule(runtime)
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->CommParticles(smpi, ispec, params, 0, &vecPatches);
    }
    #pragma omp for schedule(runtime)
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->finalizeCommParticles(smpi, ispec, params, 0, &vecPatches);
    }

    // Per direction
    for (unsigned int iDim=1 ; iDim<params.nDim_particle ; iDim++) {
        #pragma omp for schedule(runtime)
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            vecPatches(ipatch)->initCommParticles(smpi, ispec, params, iDim, &vecPatches);
        }

        #pragma omp for schedule(runtime)
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            vecPatches(ipatch)->CommParticles(smpi, ispec, params, iDim, &vecPatches);
        }
        #pragma omp for schedule(runtime)
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            vecPatches(ipatch)->finalizeCommParticles(smpi, ispec, params, iDim, &vecPatches);
        }
    }
    //timers.syncPart.update( params.printNow( itime ) );

    //timers.injectPart.restart();
//    #pragma omp for schedule(runtime)
//    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
//        vecPatches(ipatch)->injectParticles(smpi, ispec, params, params.nDim_particle-1, &vecPatches); // wait


    #pragma omp for schedule(runtime)
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
        vecPatches(ipatch)->vecSpecies[ispec]->sort_part();
    //timers.injectPart.update( params.printNow( itime ) );
}


void SyncVectorPatch::sumRhoJ(VectorPatch& vecPatches, Timers &timers, int itime)
{
    SyncVectorPatch::sum( vecPatches.densities , vecPatches, timers, itime );
    if(vecPatches.diag_flag) SyncVectorPatch::sum( vecPatches.listrho_, vecPatches, timers, itime );
}

void SyncVectorPatch::sumRhoJs(VectorPatch& vecPatches, int ispec , Timers &timers, int itime)
{

    if(vecPatches.listJxs_ .size()>0) SyncVectorPatch::sum( vecPatches.listJxs_ , vecPatches, timers, itime  );
    if(vecPatches.listJys_ .size()>0) SyncVectorPatch::sum( vecPatches.listJys_ , vecPatches, timers, itime  );
    if(vecPatches.listJzs_ .size()>0) SyncVectorPatch::sum( vecPatches.listJzs_ , vecPatches, timers, itime  );
    if(vecPatches.listrhos_.size()>0) SyncVectorPatch::sum( vecPatches.listrhos_, vecPatches, timers, itime  );
}

void SyncVectorPatch::exchangeE( VectorPatch& vecPatches )
{

    SyncVectorPatch::exchange( vecPatches.listEx_, vecPatches );
    SyncVectorPatch::exchange( vecPatches.listEy_, vecPatches );
    SyncVectorPatch::exchange( vecPatches.listEz_, vecPatches );
}

void SyncVectorPatch::finalizeexchangeE( VectorPatch& vecPatches )
{

    SyncVectorPatch::finalizeexchange( vecPatches.listEx_, vecPatches );
    SyncVectorPatch::finalizeexchange( vecPatches.listEy_, vecPatches );
    SyncVectorPatch::finalizeexchange( vecPatches.listEz_, vecPatches );
}

void SyncVectorPatch::exchangeB( VectorPatch& vecPatches )
{

    if (vecPatches.listBx_[0]->dims_.size()==1) {
        SyncVectorPatch::exchange0( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::exchange0( vecPatches.listBz_, vecPatches );
    }
    else if ( vecPatches.listBx_[0]->dims_.size()==2 ) {
        SyncVectorPatch::exchange1( vecPatches.listBx_, vecPatches );
        SyncVectorPatch::exchange0( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::exchange ( vecPatches.listBz_, vecPatches );
    }
    else if ( vecPatches.listBx_[0]->dims_.size()==3 ) {
        SyncVectorPatch::exchange1( vecPatches.listBx_, vecPatches );
        SyncVectorPatch::exchange2( vecPatches.listBx_, vecPatches );
        SyncVectorPatch::exchange0( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::exchange2( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::exchange0( vecPatches.listBz_, vecPatches );
        SyncVectorPatch::exchange1( vecPatches.listBz_, vecPatches );
    }

}

void SyncVectorPatch::finalizeexchangeB( VectorPatch& vecPatches )
{
    if (vecPatches.listBx_[0]->dims_.size()==1) {
        SyncVectorPatch::finalizeexchange0( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::finalizeexchange0( vecPatches.listBz_, vecPatches );
    }
    else if ( vecPatches.listBx_[0]->dims_.size()==2 ) {
        SyncVectorPatch::finalizeexchange1( vecPatches.listBx_, vecPatches );
        SyncVectorPatch::finalizeexchange0( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::finalizeexchange ( vecPatches.listBz_, vecPatches );
    }
    else if ( vecPatches.listBx_[0]->dims_.size()==3 ) {
        SyncVectorPatch::finalizeexchange1( vecPatches.listBx_, vecPatches );
        SyncVectorPatch::finalizeexchange2( vecPatches.listBx_, vecPatches );
        SyncVectorPatch::finalizeexchange0( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::finalizeexchange2( vecPatches.listBy_, vecPatches );
        SyncVectorPatch::finalizeexchange0( vecPatches.listBz_, vecPatches );
        SyncVectorPatch::finalizeexchange1( vecPatches.listBz_, vecPatches );
    }

}


void SyncVectorPatch::sum( std::vector<Field*> fields, VectorPatch& vecPatches, Timers &timers, int itime )
{
    unsigned int nx_, ny_, nz_, h0, oversize[3], n_space[3], gsp[3];
    double *pt1,*pt2;
    h0 = vecPatches(0)->hindex;

    int nPatches( vecPatches.size() );

    oversize[0] = vecPatches(0)->EMfields->oversize[0];
    oversize[1] = vecPatches(0)->EMfields->oversize[1];
    oversize[2] = vecPatches(0)->EMfields->oversize[2];
    
    n_space[0] = vecPatches(0)->EMfields->n_space[0];
    n_space[1] = vecPatches(0)->EMfields->n_space[1];
    n_space[2] = vecPatches(0)->EMfields->n_space[2];

    int nComp = fields.size()/nPatches;

    // -----------------
    // Sum per direction :

    // iDim = 0, initialize comms : Isend/Irecv
    #pragma omp for schedule(static) 
    for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++) {
        unsigned int ipatch = ifield%nPatches;
        vecPatches(ipatch)->initSumField( fields[ifield], 0 );
    }

//    #pragma omp for schedule(static) 
//    for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++) {
//        unsigned int ipatch = ifield%nPatches;
//        vecPatches(ipatch)->testSumField( fields[ifield], 0 );
//    }

    // iDim = 0, local
    timers.syncDens.restart();
    for (unsigned int icomp=0 ; icomp<nComp ; icomp++) {
        nx_ = fields[icomp*nPatches]->dims_[0];
        ny_ = 1;
        nz_ = 1;
        if (fields[icomp*nPatches]->dims_.size()>1) {
            ny_ = fields[icomp*nPatches]->dims_[1];
            if (fields[icomp*nPatches]->dims_.size()>2) 
                nz_ = fields[icomp*nPatches]->dims_[2];
        }
        gsp[0] = 1+2*oversize[0]+fields[icomp*nPatches]->isDual_[0]; //Ghost size primal
        #pragma omp for schedule(static) private(pt1,pt2)
        for (unsigned int ifield=icomp*nPatches ; ifield<(icomp+1)*nPatches ; ifield++) {
            unsigned int ipatch = ifield%nPatches;
            if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[0][0]){
                //The patch to the west belongs to the same MPI process than I.
                pt1 = &(*fields[vecPatches(ipatch)->neighbor_[0][0]-h0+icomp*nPatches])(n_space[0]*ny_*nz_);
                pt2 = &(*fields[ifield])(0);
                //Sum 2 ==> 1
                for (unsigned int i = 0; i < gsp[0]* ny_*nz_ ; i++) pt1[i] += pt2[i];
                //Copy back the results to 2
                memcpy( pt2, pt1, gsp[0]*ny_*nz_*sizeof(double)); 
            }
        }
    }
    timers.syncDens.update();

    // iDim = 0, finalize (waitall)
    #pragma omp for schedule(static)
    for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++){
        unsigned int ipatch = ifield%nPatches;
        vecPatches(ipatch)->finalizeSumField( fields[ifield], 0 );
    }
    // END iDim = 0 sync
    // -----------------

    if (fields[0]->dims_.size()>1) {
        // -----------------
        // Sum per direction :
        
        // iDim = 1, initialize comms : Isend/Irecv
        #pragma omp for schedule(static)
        for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++) {
            unsigned int ipatch = ifield%nPatches;
            vecPatches(ipatch)->initSumField( fields[ifield], 1 );
        }

//        #pragma omp for schedule(static) 
//        for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++) {
//            unsigned int ipatch = ifield%nPatches;
//            vecPatches(ipatch)->testSumField( fields[ifield], 1 );
//        }

        // iDim = 1, local
        timers.syncDensY.restart();
        for (unsigned int icomp=0 ; icomp<nComp ; icomp++) {
            nx_ = fields[icomp*nPatches]->dims_[0];
            ny_ = 1;
            nz_ = 1;
            if (fields[icomp*nPatches]->dims_.size()>1) {
                ny_ = fields[icomp*nPatches]->dims_[1];
                if (fields[icomp*nPatches]->dims_.size()>2) 
                    nz_ = fields[icomp*nPatches]->dims_[2];
            }
            gsp[0] = 1+2*oversize[0]+fields[icomp*nPatches]->isDual_[0]; //Ghost size primal
            gsp[1] = 1+2*oversize[1]+fields[icomp*nPatches]->isDual_[1]; //Ghost size primal

            #pragma omp for schedule(static) private(pt1,pt2)
            for (unsigned int ifield=icomp*nPatches ; ifield<(icomp+1)*nPatches ; ifield++) {
                unsigned int ipatch = ifield%nPatches;
                if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[1][0]){
                    //The patch to the south belongs to the same MPI process than I.
                    pt1 = &(*fields[vecPatches(ipatch)->neighbor_[1][0]-h0+icomp*nPatches])(n_space[1]*nz_);
                    pt2 = &(*fields[ifield])(0);
                    for (unsigned int j = 0; j < nx_ ; j++){
                        for (unsigned int i = 0; i < gsp[1]*nz_ ; i++) pt1[i] += pt2[i];
                        memcpy( pt2, pt1, gsp[1]*nz_*sizeof(double)); 
                        pt1 += ny_*nz_;
                        pt2 += ny_*nz_;
                    }
                }
            }
        }
        timers.syncDensY.update();

        // iDim = 1, finalize (waitall)
        #pragma omp for schedule(static)
        for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++){
            unsigned int ipatch = ifield%nPatches;
            vecPatches(ipatch)->finalizeSumField( fields[ifield], 1 );
        }
        // END iDim = 1 sync
        // -----------------        

        if (fields[0]->dims_.size()>2) {
            // -----------------
            // Sum per direction :

            // iDim = 2, initialize comms : Isend/Irecv
            #pragma omp for schedule(static)
            for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++) {
                unsigned int ipatch = ifield%nPatches;
                vecPatches(ipatch)->initSumField( fields[ifield], 2 );
            }

            // iDim = 2 local
            for (unsigned int icomp=0 ; icomp<nComp ; icomp++) {
                nx_ = fields[icomp*nPatches]->dims_[0];
                ny_ = 1;
                nz_ = 1;
                if (fields[icomp*nPatches]->dims_.size()>1) {
                    ny_ = fields[icomp*nPatches]->dims_[1];
                    if (fields[icomp*nPatches]->dims_.size()>2) 
                        nz_ = fields[icomp*nPatches]->dims_[2];
                }
                gsp[0] = 1+2*oversize[0]+fields[icomp*nPatches]->isDual_[0]; //Ghost size primal
                gsp[1] = 1+2*oversize[1]+fields[icomp*nPatches]->isDual_[1]; //Ghost size primal
                gsp[2] = 1+2*oversize[2]+fields[icomp*nPatches]->isDual_[2]; //Ghost size primal
                #pragma omp for schedule(static) private(pt1,pt2)
                for (unsigned int ifield=icomp*nPatches ; ifield<(icomp+1)*nPatches ; ifield++) {
                    unsigned int ipatch = ifield%nPatches;
                    if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[2][0]){
                        //The patch below me belongs to the same MPI process than I.
                        pt1 = &(*fields[vecPatches(ipatch)->neighbor_[2][0]-h0+icomp*nPatches])(n_space[2]);
                        pt2 = &(*fields[ifield])(0);
                        for (unsigned int j = 0; j < nx_*ny_ ; j++){
                            for (unsigned int i = 0; i < gsp[2] ; i++){
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
            #pragma omp for schedule(static)
            for (unsigned int ifield=0 ; ifield<fields.size() ; ifield++){
                unsigned int ipatch = ifield%nPatches;
                vecPatches(ipatch)->finalizeSumField( fields[ifield], 2 );
            }
            // END iDim = 2 sync
            // -----------------

        } // End if dims_.size()>2

    } // End if dims_.size()>1

}


void SyncVectorPatch::exchange( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    for ( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
        #pragma omp for schedule(static)
        for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
            vecPatches(ipatch)->initExchange( fields[ipatch], iDim );
    } // End for iDim


    unsigned int nx_, ny_(1), nz_(1), h0, oversize[3], n_space[3], gsp[3];
    double *pt1,*pt2;
    h0 = vecPatches(0)->hindex;

    oversize[0] = vecPatches(0)->EMfields->oversize[0];
    oversize[1] = vecPatches(0)->EMfields->oversize[1];
    oversize[2] = vecPatches(0)->EMfields->oversize[2];

    n_space[0] = vecPatches(0)->EMfields->n_space[0];
    n_space[1] = vecPatches(0)->EMfields->n_space[1];
    n_space[2] = vecPatches(0)->EMfields->n_space[2];

    nx_ = fields[0]->dims_[0];
    if (fields[0]->dims_.size()>1) {
        ny_ = fields[0]->dims_[1];
        if (fields[0]->dims_.size()>2) 
            nz_ = fields[0]->dims_[2];
    }


    gsp[0] = ( oversize[0] + 1 + fields[0]->isDual_[0] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++) {

        if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[0][0]){
            pt1 = &(*fields[vecPatches(ipatch)->neighbor_[0][0]-h0])((n_space[0])*ny_*nz_);
            pt2 = &(*fields[ipatch])(0);
            memcpy( pt2, pt1, oversize[0]*ny_*nz_*sizeof(double)); 
            memcpy( pt1+gsp[0]*ny_*nz_, pt2+gsp[0]*ny_*nz_, oversize[0]*ny_*nz_*sizeof(double)); 
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] ) 

        if (fields[0]->dims_.size()>1) {
            gsp[1] = ( oversize[1] + 1 + fields[0]->isDual_[1] ); //Ghost size primal
            if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[1][0]){
                pt1 = &(*fields[vecPatches(ipatch)->neighbor_[1][0]-h0])(n_space[1]*nz_);
                pt2 = &(*fields[ipatch])(0);
                for (unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_){
                    for (unsigned int j = 0 ; j < oversize[1]*nz_ ; j++ ){
                        // Rewrite with memcpy ?
                        pt2[i+j] = pt1[i+j] ;
                        pt1[i+j+gsp[1]*nz_] = pt2[i+j+gsp[1]*nz_] ;
                    } 
                } 
            } // End if ( MPI_me_ == MPI_neighbor_[1][0] ) 

            if (fields[0]->dims_.size()>2) {
                gsp[2] = ( oversize[2] + 1 + fields[0]->isDual_[2] ); //Ghost size primal
                if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[2][0]){
                    pt1 = &(*fields[vecPatches(ipatch)->neighbor_[2][0]-h0])(n_space[2]);
                    pt2 = &(*fields[ipatch])(0);
                    for (unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_){
                        for (unsigned int j = 0 ; j < ny_*nz_ ; j += nz_){
                            for (unsigned int k = 0 ; k < oversize[2] ; k++ ){
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


void SyncVectorPatch::finalizeexchange( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    for ( unsigned int iDim=0 ; iDim<fields[0]->dims_.size() ; iDim++ ) {
        #pragma omp for schedule(static)
        for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
            vecPatches(ipatch)->finalizeExchange( fields[ipatch], iDim );
    } // End for iDim

}


void SyncVectorPatch::exchange0( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
        vecPatches(ipatch)->initExchange( fields[ipatch], 0 );

    unsigned int ny_(1), nz_(1), h0, oversize, n_space, gsp;
    double *pt1,*pt2;
    h0 = vecPatches(0)->hindex;

    oversize = vecPatches(0)->EMfields->oversize[0];

    n_space = vecPatches(0)->EMfields->n_space[0];

    if (fields[0]->dims_.size()>1) {
        ny_ = fields[0]->dims_[1];
        if (fields[0]->dims_.size()>2) 
            nz_ = fields[0]->dims_[2];
    }

    //gsp[0] = 2*oversize[0]+fields[0]->isDual_[0]; //Ghost size primal
    //for filter
    gsp = ( oversize + 1 + fields[0]->isDual_[0] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++) {

        if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[0][0]){
            pt1 = &(*fields[vecPatches(ipatch)->neighbor_[0][0]-h0])(n_space*ny_*nz_);
            pt2 = &(*fields[ipatch])(0);
            //memcpy( pt2, pt1, ny_*sizeof(double)); 
            //memcpy( pt1+gsp[0]*ny_, pt2+gsp[0]*ny_, ny_*sizeof(double)); 
            //for filter
            memcpy( pt2, pt1, oversize*ny_*nz_*sizeof(double)); 
            memcpy( pt1+gsp*ny_*nz_, pt2+gsp*ny_*nz_, oversize*ny_*nz_*sizeof(double)); 
        } // End if ( MPI_me_ == MPI_neighbor_[0][0] ) 


    } // End for( ipatch )

}


void SyncVectorPatch::finalizeexchange0( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
        vecPatches(ipatch)->finalizeExchange( fields[ipatch], 0 );

}

void SyncVectorPatch::exchange1( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
        vecPatches(ipatch)->initExchange( fields[ipatch], 1 );

    unsigned int nx_, ny_, nz_(1), h0, oversize, n_space, gsp;
    double *pt1,*pt2;
    h0 = vecPatches(0)->hindex;

    oversize = vecPatches(0)->EMfields->oversize[1];
    n_space = vecPatches(0)->EMfields->n_space[1];

    nx_ = fields[0]->dims_[0];
    ny_ = fields[0]->dims_[1];
    if (fields[0]->dims_.size()>2) 
        nz_ = fields[0]->dims_[2];

    //gsp = 2*oversize[1]+fields[0]->isDual_[1]; //Ghost size primal
    //for filter
    gsp = ( oversize + 1 + fields[0]->isDual_[1] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++) {

        if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[1][0]){
            pt1 = &(*fields[vecPatches(ipatch)->neighbor_[1][0]-h0])(n_space*nz_);
            pt2 = &(*fields[ipatch])(0);
            for (unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_){
                // for filter
                for (unsigned int j = 0 ; j < oversize*nz_ ; j++ ){
                    pt2[i+j] = pt1[i+j] ;
                    pt1[i+j+gsp*nz_] = pt2[i+j+gsp*nz_] ;
                } // mempy to do
            } 
        } // End if ( MPI_me_ == MPI_neighbor_[1][0] ) 

    } // End for( ipatch )

}


void SyncVectorPatch::finalizeexchange1( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
        vecPatches(ipatch)->finalizeExchange( fields[ipatch], 1 );

}



void SyncVectorPatch::exchange2( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
        vecPatches(ipatch)->initExchange( fields[ipatch], 2 );

    unsigned int nx_, ny_, nz_, h0, oversize, n_space, gsp;
    double *pt1,*pt2;
    h0 = vecPatches(0)->hindex;

    oversize = vecPatches(0)->EMfields->oversize[2];
    n_space = vecPatches(0)->EMfields->n_space[2];

    nx_ = fields[0]->dims_[0];
    ny_ = fields[0]->dims_[1];
    nz_ = fields[0]->dims_[2];

    //gsp = 2*oversize[1]+fields[0]->isDual_[1]; //Ghost size primal
    //for filter
    gsp = ( oversize + 1 + fields[0]->isDual_[2] ); //Ghost size primal

    #pragma omp for schedule(static) private(pt1,pt2)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++) {

        if (vecPatches(ipatch)->MPI_me_ == vecPatches(ipatch)->MPI_neighbor_[2][0]){
           pt1 = &(*fields[vecPatches(ipatch)->neighbor_[2][0]-h0])(n_space);
           pt2 = &(*fields[ipatch])(0);
           for (unsigned int i = 0 ; i < nx_*ny_*nz_ ; i += ny_*nz_){
               for (unsigned int j = 0 ; j < ny_*nz_ ; j += nz_){
                   for (unsigned int k = 0 ; k < oversize ; k++ ){
                       pt2[i+j+k] = pt1[i+j+k] ;
                       pt1[i+j+k+gsp] = pt2[i+j+k+gsp] ;
                   } 
               }
           } 
        } // End if ( MPI_me_ == MPI_neighbor_[2][0] ) 

    } // End for( ipatch )
}


void SyncVectorPatch::finalizeexchange2( std::vector<Field*> fields, VectorPatch& vecPatches )
{
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<fields.size() ; ipatch++)
        vecPatches(ipatch)->finalizeExchange( fields[ipatch], 2 );

}
