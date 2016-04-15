#include "VectorPatch.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"
#include "SmileiIOFactory.h"
#include "SimWindow.h"
#include "SolverFactory.h"
#include "DiagnosticFactory.h"

#include "SyncVectorPatch.h"
#include "DiagsVectorPatch.h"

#include "Timer.h"

#include <cstring>
//#include <string>

using namespace std;


VectorPatch::VectorPatch()
{
}


VectorPatch::~VectorPatch()
{
}

void VectorPatch::close(SmileiMPI * smpiData)
{
    closeAllDiags( smpiData );
    
    for (unsigned int idiag=0 ; idiag<globalDiags.size(); idiag++)
        delete globalDiags[idiag];
    globalDiags.clear();
    
    patches_[0]->EMfields->clean(); // This destructs the laser profiles only once
    for (unsigned int ipatch=0 ; ipatch<size(); ipatch++)
        delete patches_[ipatch];
    
    patches_.clear();
}

void VectorPatch::createGlobalDiags(Params& params, SmileiMPI* smpi)
{
    globalDiags = DiagnosticFactory::createGlobalDiagnostics(params, smpi, (*this)(0) );

}


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------       INTERFACES        ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// For all patch, move particles (restartRhoJ(s), dynamics and exchangeParticles)
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::dynamics(Params& params, SmileiMPI* smpi, SimWindow* simWindow,
                           int* diag_flag, double time_dual, vector<Timer>& timer)
{
    timer[1].restart();
    
    #pragma omp for schedule(runtime)
    for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
        (*this)(ipatch)->EMfields->restartRhoJ();
        for (unsigned int ispec=0 ; ispec<(*this)(ipatch)->vecSpecies.size() ; ispec++) {
            if ( (*this)(ipatch)->vecSpecies[ispec]->isProj(time_dual, simWindow) || diag_flag  ) {
                species(ipatch, ispec)->dynamics(time_dual, ispec,
                                                 emfields(ipatch), interp(ipatch), proj(ipatch),
                                                 params, *diag_flag, partwalls(ipatch),
                                                 (*this)(ipatch), smpi);
            }
        }
    
    }
    timer[1].update();
    
    timer[8].restart();
    for (unsigned int ispec=0 ; ispec<(*this)(0)->vecSpecies.size(); ispec++) {
        if ( (*this)(0)->vecSpecies[ispec]->isProj(time_dual, simWindow) ){
            SyncVectorPatch::exchangeParticles((*this), ispec, params, smpi ); // Included sort_part
        }
    }
    timer[8].update();

} // END dynamics


// ---------------------------------------------------------------------------------------------------------------------
// For all patch, sum densities on ghost cells (sum per species if needed, sync per patch and MPI sync)
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::sumDensities( int* diag_flag, vector<Timer>& timer )
{
    timer[4].restart();
    if  (*diag_flag){
        #pragma omp for schedule(static)
        for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
             // Per species in global, Attention if output -> Sync / per species fields
            (*this)(ipatch)->EMfields->computeTotalRhoJ();
        }
    }
    timer[4].update();
    
    timer[9].restart();
    SyncVectorPatch::sumRhoJ( (*this), *diag_flag ); // MPI
    
    if(*diag_flag){
        for (unsigned int ispec=0 ; ispec<(*this)(0)->vecSpecies.size(); ispec++) {
            update_field_list(ispec);
            SyncVectorPatch::sumRhoJs( (*this), ispec ); // MPI
        }
    }
    timer[9].update();
    
} // End sumDensities


// ---------------------------------------------------------------------------------------------------------------------
// For all patch, update E and B (Ampere, Faraday, boundary conditions, exchange B and center B)
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::solveMaxwell(Params& params, SimWindow* simWindow, int itime, double time_dual, vector<Timer>& timer)
{
    timer[2].restart();
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++){
        // Saving magnetic fields (to compute centered fields used in the particle pusher)
        // Stores B at time n in B_m.
        (*this)(ipatch)->EMfields->saveMagneticFields();
        // Computes Ex_, Ey_, Ez_ on all points.
        // E is already synchronized because J has been synchronized before.
        (*this)(ipatch)->EMfields->solveMaxwellAmpere();
    }
    //(*this).exchangeE();
    
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++){
        // Computes Bx_, By_, Bz_ at time n+1 on interior points.
        // (*this)(ipatch)->EMfields->solveMaxwellFaraday();
        (*(*this)(ipatch)->EMfields->MaxwellFaradaySolver_)((*this)(ipatch)->EMfields);
        // Applies boundary conditions on B
        (*this)(ipatch)->EMfields->boundaryConditions(itime, time_dual, (*this)(ipatch), params, simWindow);
    }
    //Synchronize B fields between patches.
    timer[2].update();
    
    timer[9].restart();
    SyncVectorPatch::exchangeB( (*this) );
    timer[9].update();
    
    timer[2].restart();
    // Computes B at time n+1/2 using B and B_m.
    #pragma omp for schedule(static)
    for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++)
        (*this)(ipatch)->EMfields->centerMagneticFields();
    timer[2].update();

} // END solveMaxwell



void VectorPatch::initAllDiags(Params& params, SmileiMPI* smpi)
{
    // globalDiags : scalars + particles
    for (unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++) {
        // For all diags createFile (if global, juste by mpi master, test inside)
        if ( smpi->isMaster() )
            globalDiags[idiag]->openFile( params, smpi, *this, true );
    
    } // END for globalDiags
    
    
    // Initialize scalars (globalDiags[0], especially energy balance)
    for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
        globalDiags[0]->run( (*this)(ipatch), 0 );
    }
    smpi->computeGlobalDiags( globalDiags[0], 0 );
    
    // localDiags : probes, track & fields
    for (unsigned int idiag = 0 ; idiag < (*this)(0)->localDiags.size() ; idiag++) {
        
        // track, compute global number of particles + compute global Idx
        if ( (*this)(0)->localDiags[idiag]->type_ == "Track" ) {
            DiagnosticTrack* diagTrack0 = static_cast<DiagnosticTrack*>( (*this)(0)->localDiags[idiag] );
            diagTrack0->setFileSplitting( params, smpi, *this );
            for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
                DiagnosticTrack* diagTrack  = static_cast<DiagnosticTrack*>( (*this)(ipatch)->localDiags[idiag] );
                diagTrack->setFile( diagTrack0->getFileId() );
            }
        }
        
        // For all diags createFile
        //(*this)(0)->localDiags[idiag]->openFile( params, smpi, *this, true );
        
        if ( (*this)(0)->localDiags[idiag]->type_ == "Probes" ) {
            (*this)(0)->localDiags[idiag]->openFile( params, smpi, *this, true );
            DiagnosticProbes* diagProbes0 = static_cast<DiagnosticProbes*>( (*this)(0)->localDiags[idiag] );
            diagProbes0->setFileSplitting( params, smpi, *this );
            for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
                DiagnosticProbes* diagProbes = static_cast<DiagnosticProbes*>( (*this)(ipatch)->localDiags[idiag] );
                diagProbes->setFile( diagProbes0->getFileId() );
                diagProbes->writePositionIn( params );
            }// END  ipatch
        
        } // END if Probes
    
    } // END for localDiags
    
} // END initAllDiags


void VectorPatch::closeAllDiags(SmileiMPI* smpi)
{
    for (unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++) {
        if ( smpi->isMaster() )
            globalDiags[idiag]->closeFile();
    }
    for (unsigned int idiag = 0 ; idiag < (*this)(0)->localDiags.size() ; idiag++) {
        (*this)(0)->localDiags[idiag]->closeFile();
    }

}


void VectorPatch::openAllDiags(Params& params,SmileiMPI* smpi)
{
    for (unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++) {
        // For all diags createFile (if global, juste by mpi master, test inside)
        if ( smpi->isMaster() )
            globalDiags[idiag]->openFile( params, smpi, *this, false );
    
    } // END for globalDiags
    
    
    for (unsigned int idiag = 0 ; idiag < (*this)(0)->localDiags.size() ; idiag++) {
        (*this)(0)->localDiags[idiag]->openFile( params, smpi, *this, false );
        
        if ( (*this)(0)->localDiags[idiag]->type_ == "Track" ) {
            DiagnosticTrack* diagTrack0 = static_cast<DiagnosticTrack*>( (*this)(0)->localDiags[idiag] );
            // Spli
            for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
                DiagnosticTrack* diagTrack  = static_cast<DiagnosticTrack*>( (*this)(ipatch)->localDiags[idiag] );
                diagTrack->setFile( diagTrack0->getFileId() );
            }
        }
        
        if ( (*this)(0)->localDiags[idiag]->type_ == "Probes" ) {
            DiagnosticProbes* diagProbes0 = static_cast<DiagnosticProbes*>( (*this)(0)->localDiags[idiag] );
            //diagProbes0->setFileSplitting( params, smpi, *this );
            for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
                DiagnosticProbes* diagProbes = static_cast<DiagnosticProbes*>( (*this)(ipatch)->localDiags[idiag] );
                diagProbes->setFile( diagProbes0->getFileId() );
            }
        }
    }
    
    /*
    // Should be OK with the implementation below :
    //   abstract splitting too even if completely non sense for global diags ?
    
    for (unsigned int idiag = 0 ; idiag < (*this)(0)->localDiags.size() ; idiag++) {
        (*this)(0)->localDiags[idiag]->openFile( params, smpi, *this, false );
        for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++)
            (*this)(ipatch)->localDiags[idiag]->setFile( (*this)(0)->localDiags[idiag] );
        if ( (*this)(0)->localDiags[idiag]->type_ == "Probes" )
            static_cast<DiagnosticProbes*>( (*this)(0)->localDiags[idiag] )->setFileSplitting( params, smpi, *this );
    }
    */


}


// ---------------------------------------------------------------------------------------------------------------------
// For all patch, Compute and Write all diags
//   - Scalars, Probes, Phases, TrackParticles, Fields, Average fields
//   - set diag_flag to 0 after write
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::runAllDiags(Params& params, SmileiMPI* smpi, int* diag_flag, int itime, vector<Timer>& timer)
{
    // Dump Fields
    // -------------------------------------------
    timer[6].restart();
    // diag_flag = 1 if  :
    //   vecPatches.Diags->field_timeSelection->theTimeIsNow(itime)
    if  (*diag_flag){
        for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
            
            // Write EM fields dump in Fields.h5
            if (ipatch==0) (*this)(ipatch)->sio->createTimeStepInSingleFileTime( itime );
            (*this)(ipatch)->sio->writeAllFieldsSingleFileTime( (*this)(ipatch)->EMfields->allFields, itime, 0 );
            
            // Check the dedicated fields output write frequency 
            if( (*this)(ipatch)->sio->dumpAvgFields_ && avgFieldTimeIsNow(itime) ) {
                // Write EM average fields dump in Fields_avg.h5
                (*this)(ipatch)->sio->writeAllFieldsSingleFileTime( (*this)(ipatch)->EMfields->allFields_avg, itime, 1 );
            }
            // Re-init rho, Jxyz per species for next diag timestep
            (*this)(ipatch)->EMfields->restartRhoJs();
        }
        *diag_flag = 0 ;
    }
    timer[6].update();
    
    
    // Diagnostics : compute locally
    //    Parallel write for Probes, TrackParticles
    // -------------------------------------------
    timer[3].restart();
    // globalDiags : scalars + particles
    static_cast<DiagnosticScalar*>( globalDiags[0] )->reset( itime );
    for (unsigned int idiag = 0 ; idiag < globalDiags.size() ; idiag++) {
        for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
            globalDiags[idiag]->run( (*this)(ipatch), itime );
        }
        // particles
        smpi->computeGlobalDiags( globalDiags[idiag], itime);
        
    } // END for globalDiags
    
    // move only scalars write call from SmileiMPI::computeGlobalDiags()
    // Diagnostic Particles includes dedicated cleaning 
    globalDiags[0]->write(itime);
    
    
    // localDiags : probes, track & fields
    for (unsigned int idiag = 0 ; idiag < (*this)(0)->localDiags.size() ; idiag++) {
        (*this)(0)->localDiags[idiag]->prepare( (*this)(0), itime );
        for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
                (*this)(ipatch)->localDiags[idiag]->run( (*this)(ipatch), itime );
            }
    }
    for (unsigned int idiag = 0 ; idiag < (*this)(0)->localDiags.size() ; idiag++) {
        for (unsigned int ipatch=0 ; ipatch<(*this).size() ; ipatch++) {
              (*this)(ipatch)->localDiags[idiag]->write( itime );
            }
    }
    /*
        // track, compute global number of particles + compute global Idx
        if ( (*this)(0)->localDiags[idiag]->type_ == "Track" ) {
            DiagnosticTrack* diagTrack = static_cast<DiagnosticTrack*>( (*this)(0)->localDiags[idiag] );
            diagTrack->setFileSize( params, smpi, *this );
        }
        
    } // END for localDiags
    */
    timer[3].update();   

} // END runAllDiags


// ---------------------------------------------------------------------------------------------------------------------
// Check if rho is null (MPI & patch sync)
// ---------------------------------------------------------------------------------------------------------------------
bool VectorPatch::isRhoNull( SmileiMPI* smpi )
{
    double norm2(0.);
    double locnorm2(0.);
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
        locnorm2 += (*this)(ipatch)->EMfields->computeRhoNorm2();
    
    MPI_Allreduce(&locnorm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    return (norm2<=0.);
} // END isRhoNull


// ---------------------------------------------------------------------------------------------------------------------
// Solve Poisson to initialize E
//   - all steps are done locally, sync per patch, sync per MPI process 
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::solvePoisson( Params &params, SmileiMPI* smpi )
{
    unsigned int nx_p2_global = (params.n_space_global[0]+1) * (params.n_space_global[1]+1);
    
    unsigned int iteration_max = 50000;
    double       error_max     = 1.e-14;
    unsigned int iteration=0;
    
    // Init & Store internal data (phi, r, p, Ap) per patch
    double rnew_dot_rnew_local(0.);
    double rnew_dot_rnew(0.);    
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
        (*this)(ipatch)->EMfields->initPoisson( (*this)(ipatch) );
        rnew_dot_rnew_local += (*this)(ipatch)->EMfields->compute_r();
    }
    MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    std::vector<Field*> Ex_;
    std::vector<Field*> Ap_;
    
    for (int ipatch=0 ; ipatch<this->size() ; ipatch++) {
        Ex_.push_back( (*this)(ipatch)->EMfields->Ex_ );
        Ap_.push_back( (*this)(ipatch)->EMfields->Ap_ );
    }
    
    // compute control parameter
    double ctrl = rnew_dot_rnew / (double)(nx_p2_global);
    
    // ---------------------------------------------------------
    // Starting iterative loop for the conjugate gradient method
    // ---------------------------------------------------------
    if (smpi->isMaster()) DEBUG("Starting iterative loop for CG method");
    while ( (ctrl > error_max) && (iteration<iteration_max) ) {
        
        iteration++;
        if (smpi->isMaster()) DEBUG("iteration " << iteration << " started with control parameter ctrl = " << ctrl*1.e14 << " x 1e-14");
        
        // scalar product of the residual
        double r_dot_r = rnew_dot_rnew;
        
        for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) 
            (*this)(ipatch)->EMfields->compute_Ap( (*this)(ipatch) );
        
        // Exchange Ap_ (intra & extra MPI)
        SyncVectorPatch::exchange( Ap_, *this );
        
       // scalar product p.Ap
        double p_dot_Ap       = 0.0;
        double p_dot_Ap_local = 0.0;
        for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
            p_dot_Ap_local += (*this)(ipatch)->EMfields->compute_pAp();
        }
        MPI_Allreduce(&p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
        // compute new potential and residual
        for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
            (*this)(ipatch)->EMfields->update_pand_r( r_dot_r, p_dot_Ap );
        }
        
        // compute new residual norm
        rnew_dot_rnew       = 0.0;
        rnew_dot_rnew_local = 0.0;
        for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
            rnew_dot_rnew_local += (*this)(ipatch)->EMfields->compute_r();
        }
        MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (smpi->isMaster()) DEBUG("new residual norm: rnew_dot_rnew = " << rnew_dot_rnew);
        
        // compute new directio
        for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
            (*this)(ipatch)->EMfields->update_p( rnew_dot_rnew, r_dot_r );
        }
        
        // compute control parameter
        ctrl = rnew_dot_rnew / (double)(nx_p2_global);
        if (smpi->isMaster()) DEBUG("iteration " << iteration << " done, exiting with control parameter ctrl = " << ctrl);
    
    }//End of the iterative loop
    
    
    // --------------------------------
    // Status of the solver convergence
    // --------------------------------
    if (iteration == iteration_max) {
        if (smpi->isMaster())
            WARNING("Poisson solver did not converge: reached maximum iteration number: " << iteration
                    << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    else {
        if (smpi->isMaster()) 
            MESSAGE(1,"Poisson solver converged at iteration: " << iteration
                    << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    
    // ------------------------------------------
    // Compute the electrostatic fields Ex and Ey
    // ------------------------------------------
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
        (*this)(ipatch)->EMfields->initE( (*this)(ipatch) );
    
    SyncVectorPatch::exchangeE( *this );    
    
    // Centering of the electrostatic fields
    // -------------------------------------
    
    vector<double> E_Add(Ex_[0]->dims_.size(),0.);
    if ( Ex_[0]->dims_.size()>1 ) {
        double Ex_WestNorth = 0.0;
        double Ey_WestNorth = 0.0;
        double Ex_EastSouth = 0.0;
        double Ey_EastSouth = 0.0;
        
        //The NorthWest patch has Patch coordinates X=0, Y=2^m1-1= number_of_patches[1]-1.
        //Its hindex is
        int patch_NorthWest = generalhilbertindex(params.mi[0], params.mi[1], 0,  params.number_of_patches[1]-1);
        //The MPI rank owning it is
        int rank_WestNorth = smpi->hrank(patch_NorthWest);
        //The SouthEast patch has Patch coordinates X=2^m0-1= number_of_patches[0]-1, Y=0.
        //Its hindex is
        int patch_SouthEast = generalhilbertindex(params.mi[0], params.mi[1], params.number_of_patches[0]-1, 0);
        //The MPI rank owning it is
        int rank_EastSouth = smpi->hrank(patch_SouthEast);
        
        
        //cout << params.mi[0] << " " << params.mi[1] << " " << params.number_of_patches[0] << " " << params.number_of_patches[1] << endl;
        //cout << patch_NorthWest << " " << rank_WestNorth << " " << patch_SouthEast << " " << rank_EastSouth << endl;
        
        if ( smpi->getRank() == rank_WestNorth ) {
            Ex_WestNorth = (*this)(patch_NorthWest-((*this).refHindex_))->EMfields->getEx_WestNorth();
            Ey_WestNorth = (*this)(patch_NorthWest-((*this).refHindex_))->EMfields->getEy_WestNorth();
        }
        
        // East-South corner
        if ( smpi->getRank() == rank_EastSouth ) {
            Ex_EastSouth = (*this)(patch_SouthEast-((*this).refHindex_))->EMfields->getEx_EastSouth();
            Ey_EastSouth = (*this)(patch_SouthEast-((*this).refHindex_))->EMfields->getEy_EastSouth();
        }
        
        MPI_Bcast(&Ex_WestNorth, 1, MPI_DOUBLE, rank_WestNorth, MPI_COMM_WORLD);
        MPI_Bcast(&Ey_WestNorth, 1, MPI_DOUBLE, rank_WestNorth, MPI_COMM_WORLD);
        
        MPI_Bcast(&Ex_EastSouth, 1, MPI_DOUBLE, rank_EastSouth, MPI_COMM_WORLD);
        MPI_Bcast(&Ey_EastSouth, 1, MPI_DOUBLE, rank_EastSouth, MPI_COMM_WORLD);
        
        //This correction is always done, independantly of the periodicity. Is this correct ?
        E_Add[0] = -0.5*(Ex_WestNorth+Ex_EastSouth);
        E_Add[1] = -0.5*(Ey_WestNorth+Ey_EastSouth);
    
    }
    else if( Ex_[0]->dims_.size()==1 ) {
        double Ex_West = 0.0;
        double Ex_East = 0.0;
        
        unsigned int rankWest = 0;
        if ( smpi->getRank() == 0 ) {
            //Ex_West = (*Ex1D)(index_bc_min[0]);
            Ex_West = (*this)( (0)-((*this).refHindex_))->EMfields->getEx_West();
        }
        MPI_Bcast(&Ex_West, 1, MPI_DOUBLE, rankWest, MPI_COMM_WORLD);
        
        unsigned int rankEast = smpi->getSize()-1;
        if ( smpi->getRank() == smpi->getSize()-1 ) {
            //Ex_East = (*Ex1D)(index_bc_max[0]);
            Ex_East = (*this)( (params.number_of_patches[0]-1)-((*this).refHindex_))->EMfields->getEx_East();
        }
        MPI_Bcast(&Ex_East, 1, MPI_DOUBLE, rankEast, MPI_COMM_WORLD);
        E_Add[0] = -0.5*(Ex_West+Ex_East);
        
    }
    
    // Centering electrostatic fields
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
        (*this)(ipatch)->EMfields->centeringE( E_Add );
    
    
    // Compute error on the Poisson equation
    double deltaPoisson_max = 0.0;
    int i_deltaPoisson_max  = -1;
    
#ifdef _A_FINALISER
    for (unsigned int i=0; i<nx_p; i++) {
        double deltaPoisson = abs( ((*Ex1D)(i+1)-(*Ex1D)(i))/dx - (*rho1D)(i) );
        if (deltaPoisson > deltaPoisson_max) {
            deltaPoisson_max   = deltaPoisson;
            i_deltaPoisson_max = i;
        }
    }
#endif
    
    //!\todo Reduce to find global max
    if (smpi->isMaster())
        MESSAGE(1,"Poisson equation solved. Maximum error = " << deltaPoisson_max << " at i= " << i_deltaPoisson_max);

} // END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------    BALANCING METHODS    ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Explicits patch movement regarding new patch distribution stored in smpi->patch_count
//   - compute send_patch_id_
//   - compute recv_patch_id_
//   - create empty (not really, created like at t0) new patch in recv_patches_
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::createPatches(Params& params, SmileiMPI* smpi, SimWindow* simWindow)
{
    unsigned int n_moved(0), nPatches_now;
    recv_patches_.resize(0);
    
    // Set Index of the 1st patch of the vector yet on current MPI rank
    // Is this really necessary ? It should be done already ...
    refHindex_ = (*this)(0)->Hindex();
    
    //When going to openMP, these two vectors must be stored by patch and not by vectorPatch.
    recv_patch_id_.clear();
    send_patch_id_.clear();
    
    
    // define recv_patches_ parsing patch_count
    // Go to 1st patch to recv (maybe yet on current CPU)
    // istart = Index of the futur 1st patch
    // recv : store real Hindex
    int istart( 0 );
    for (int irk=0 ; irk<smpi->getRank() ; irk++) istart += smpi->patch_count[irk];
    //recv_patch_id stores all the hindex this process must own at the end of the exchange.
    for (int ipatch=0 ; ipatch<smpi->patch_count[smpi->getRank()] ; ipatch++)
        recv_patch_id_.push_back( istart+ipatch );
    
    
    // define send_patches_ parsing patch_count
    // send_patch_id_ stores indices from 0 to current npatch(before exchange)
    //for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
    //    send_patch_id_.push_back( ipatch );
    //}
    //Current number of patch
    nPatches_now = this->size() ;
    
    
    //std::vector<int> tmp(0);
    //Loop on current patches...
    //for (unsigned int ipatch=0 ; ipatch<send_patch_id_.size() ; ipatch++)
    for (unsigned int ipatch=0 ; ipatch < nPatches_now ; ipatch++)
      //if        current hindex        <  future refHindex  OR current hindex > future last hindex...
        if ( ( refHindex_+ipatch < recv_patch_id_[0] ) || ( refHindex_+ipatch > recv_patch_id_.back() ) )
      //    put this patch in tmp. We will have to send it away.
            //tmp.push_back( ipatch );
            send_patch_id_.push_back( ipatch );
    
    //  nPatches <- future number of patches owned.
    // Backward loop on future patches...
    for ( int ipatch=recv_patch_id_.size()-1 ; ipatch>=0 ; ipatch--) {
      //if      future patch hindex  >= current refHindex             AND    future patch hindex <= current last hindex
        //if ( ( recv_patch_id_[ipatch]>=refHindex_+send_patch_id_[0] ) && ( recv_patch_id_[ipatch]<=refHindex_+send_patch_id_[send_patch_id_.size()-1] ) ) {
        //                                          send_patch_id_[0] should be equal to 0 ??
        if ( ( recv_patch_id_[ipatch]>=refHindex_ ) && ( recv_patch_id_[ipatch] <= refHindex_ + nPatches_now - 1 ) ) {
            //Remove this patch from the receive list because I already own it.
            recv_patch_id_.erase( recv_patch_id_.begin()+ipatch );
        }
    }
    
    //send_patch_id_ = tmp;
    
    if (simWindow) n_moved = simWindow->getNmoved(); 
    // Store in local vector future patches
    // Loop on the patches I have to receive and do not already own.
    for (unsigned int ipatch=0 ; ipatch < recv_patch_id_.size() ; ipatch++) {
        // density profile is initializes as if t = 0 !
        // Species will be cleared when, nbr of particles will be known
        //Creation of a new patch, ready to receive its content from MPI neighbours.
        Patch* newPatch = PatchesFactory::create(params, smpi, recv_patch_id_[ipatch], n_moved );
        //Store pointers to newly created patch in recv_patches_.
        recv_patches_.push_back( newPatch );
    }

} // END createPatches


// ---------------------------------------------------------------------------------------------------------------------
// Exchange patches, based on createPatches initialization
//   take care of reinitialize patch master and diag file managment
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::exchangePatches(SmileiMPI* smpi, Params& params)
{
    (*this).closeAllDiags(smpi);
    
    hid_t globalFile    = (*this)(0)->sio->global_file_id_;
    hid_t globalFileAvg = (*this)(0)->sio->global_file_id_avg;
    
    int nSpecies( (*this)(0)->vecSpecies.size() );
    int newMPIrank, oldMPIrank;
    //int newMPIrankbis, oldMPIrankbis, tmp;
    newMPIrank = smpi->getRank() -1;
    oldMPIrank = smpi->getRank() -1;
    int istart( 0 );
    int nmessage = 2*nSpecies+10;
    
    
    for (int irk=0 ; irk<smpi->getRank() ; irk++) istart += smpi->patch_count[irk];
    // Send part
    // Send particles
    for (unsigned int ipatch=0 ; ipatch < send_patch_id_.size() ; ipatch++) {
        // locate rank which will own send_patch_id_[ipatch]
        // We assume patches are only exchanged with neighbours.
        // Once all patches supposed to be sent to the left are done, we send the rest to the right.
      //if   hindex of patch to be sent              >  future hindex of the first patch owned by this process 
        if(send_patch_id_[ipatch]+refHindex_ > istart ) newMPIrank = smpi->getRank() + 1;
        //cout << "Rank " << smpi->getRank() << " sending patch " << send_patch_id_[ipatch]+refHindex_ << " to " << newMPIrank << endl; 
        //newMPIrankbis = 0 ;
        //tmp = smpi->patch_count[newMPIrankbis];
        //while ( tmp <= send_patch_id_[ipatch]+refHindex_ ) {
        //    newMPIrankbis++;
        //    tmp += smpi->patch_count[newMPIrankbis];
        //}
        
        //if (newMPIrank != newMPIrankbis){
        //    cout << "newMIPrank problem ! " << newMPIrank << endl;
        //    newMPIrank = newMPIrankbis ;
        //}
        
        smpi->isend( (*this)(send_patch_id_[ipatch]), newMPIrank, (refHindex_+send_patch_id_[ipatch])*nmessage );
    }
    
    for (unsigned int ipatch=0 ; ipatch < recv_patch_id_.size() ; ipatch++) {
      //if   hindex of patch to be received > first hindex actually owned, that means it comes from the next MPI process and not from the previous anymore. 
        if(recv_patch_id_[ipatch] > refHindex_ ) oldMPIrank = smpi->getRank() + 1;
        //cout << "Rank " << smpi->getRank() << " receiving patch " << recv_patch_id_[ipatch] << " from " << oldMPIrank << endl; 
        //oldMPIrankbis = 0 ; // Comparing recv_patch_id_[ipatch] to 1st yet on current MPI rank
        //if ( recv_patch_id_[ipatch] > refHindex_ )
        //    oldMPIrankbis = smpi->getRank()+1;
        //else
        //    oldMPIrankbis = smpi->getRank()-1;
        
        //if (oldMPIrank != oldMPIrankbis){
        //    cout << "oldMIPrank problem ! " << oldMPIrank << endl;
        //    oldMPIrank = oldMPIrankbis ;
        //}
        smpi->recv( recv_patches_[ipatch], oldMPIrank, recv_patch_id_[ipatch]*nmessage, params );
    }
    
    smpi->barrier();
    //Delete sent patches
    int nPatchSend(send_patch_id_.size());
    for (int ipatch=nPatchSend-1 ; ipatch>=0 ; ipatch--) {
        //Ok while at least 1 old patch stay inon current CPU
        (*this)(send_patch_id_[ipatch])->sio->setFiles(0,0);
        delete (*this)(send_patch_id_[ipatch]);
        patches_[ send_patch_id_[ipatch] ] = NULL;
        patches_.erase( patches_.begin() + send_patch_id_[ipatch] );
        
    }
    //Put received patches in the global vecPatches
    for (unsigned int ipatch=0 ; ipatch<recv_patch_id_.size() ; ipatch++) {
        if ( recv_patch_id_[ipatch] > refHindex_ )
            patches_.push_back( recv_patches_[ipatch] );
        else
            patches_.insert( patches_.begin()+ipatch, recv_patches_[ipatch] );
    }
    recv_patches_.clear();
    
    for (int ipatch=0 ; ipatch<patches_.size() ; ipatch++ ) { 
        (*this)(ipatch)->updateMPIenv(smpi);
    }
    
    (*this).openAllDiags(params,smpi);
    
    //definePatchDiagsMaster();
    DiagsVectorPatch::definePatchDiagsMaster( *this, globalFile, globalFileAvg );
    DiagsVectorPatch::updatePatchFieldDump( *this, params );
    
    (*this).set_refHindex() ;
    
    update_field_list() ;    
    
} // END exchangePatches

// ---------------------------------------------------------------------------------------------------------------------
// Write in a file patches communications
//   - Send/Recv MPI rank
//   - Send/Recv patch Id
// ---------------------------------------------------------------------------------------------------------------------
void VectorPatch::output_exchanges(SmileiMPI* smpi)
{
    ofstream output_file;
    ostringstream name("");
    name << "debug_output"<<smpi->getRank()<<".txt" ;
    output_file.open(name.str().c_str(), std::ofstream::out | std::ofstream::app);
    int newMPIrank, oldMPIrank;
    newMPIrank = smpi->getRank() -1;
    oldMPIrank = smpi->getRank() -1;
    int istart( 0 );
    for (int irk=0 ; irk<smpi->getRank() ; irk++) istart += smpi->patch_count[irk];
    for (unsigned int ipatch=0 ; ipatch < send_patch_id_.size() ; ipatch++) {
        if(send_patch_id_[ipatch]+refHindex_ > istart ) newMPIrank = smpi->getRank() + 1;
        output_file << "Rank " << smpi->getRank() << " sending patch " << send_patch_id_[ipatch]+refHindex_ << " to " << newMPIrank << endl; 
    }
    for (unsigned int ipatch=0 ; ipatch < recv_patch_id_.size() ; ipatch++) {
        if(recv_patch_id_[ipatch] > refHindex_ ) oldMPIrank = smpi->getRank() + 1;
        output_file << "Rank " << smpi->getRank() << " receiving patch " << recv_patch_id_[ipatch] << " from " << oldMPIrank << endl; 
    }
    output_file << "NEXT" << endl;
    output_file.close();
} // END output_exchanges

    //! Resize vector of field*
void VectorPatch::update_field_list()
{
    listJx_.resize( size() ) ;
    listJy_.resize( size() ) ;
    listJz_.resize( size() ) ;
    listrho_.resize( size() ) ;
    listEx_.resize( size() ) ;
    listEy_.resize( size() ) ;
    listEz_.resize( size() ) ;
    listBx_.resize( size() ) ;
    listBy_.resize( size() ) ;
    listBz_.resize( size() ) ;
    
    for (int ipatch=0 ; ipatch < size() ; ipatch++) {
        listJx_[ipatch] = patches_[ipatch]->EMfields->Jx_ ;
        listJy_[ipatch] = patches_[ipatch]->EMfields->Jy_ ;
        listJz_[ipatch] = patches_[ipatch]->EMfields->Jz_ ;
        listrho_[ipatch] =patches_[ipatch]->EMfields->rho_;
        listEx_[ipatch] = patches_[ipatch]->EMfields->Ex_ ;
        listEy_[ipatch] = patches_[ipatch]->EMfields->Ey_ ;
        listEz_[ipatch] = patches_[ipatch]->EMfields->Ez_ ;
        listBx_[ipatch] = patches_[ipatch]->EMfields->Bx_ ;
        listBy_[ipatch] = patches_[ipatch]->EMfields->By_ ;
        listBz_[ipatch] = patches_[ipatch]->EMfields->Bz_ ;
    }
}
void VectorPatch::update_field_list(int ispec)
{
    #pragma omp single
    {
        listJxs_.resize( size() ) ;
        listJys_.resize( size() ) ;
        listJzs_.resize( size() ) ;
        listrhos_.resize( size() ) ;
    }
    
    #pragma omp for schedule(static)
    for (int ipatch=0 ; ipatch < size() ; ipatch++) {
        listJxs_[ipatch] = patches_[ipatch]->EMfields->Jx_s[ispec] ;
        listJys_[ipatch] = patches_[ipatch]->EMfields->Jy_s[ispec] ;
        listJzs_[ipatch] = patches_[ipatch]->EMfields->Jz_s[ispec] ;
        listrhos_[ipatch] =patches_[ipatch]->EMfields->rho_s[ispec];
    }
}


