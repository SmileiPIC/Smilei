
#include "Region.h"

#include "PatchesFactory.h"
#include "DomainDecompositionFactory.h"
#include "DiagnosticFields1D.h"
#include "ElectroMagn.h"
#include "Solver.h"

#include "Params.h"
#include "OpenPMDparams.h"
#include "SimWindow.h"
#include "Timers.h"
#include "DoubleGrids.h"

using namespace std;

Region::Region( Params &params ) :
    vecPatch_( params ),
    decomposition_( NULL ),
    patch_( NULL ),
    diag_( NULL ),
    fake_patch( NULL ),
    coupled_( false )
{
}

void Region::build( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, OpenPMDparams &openPMD, bool global_region )
{
    // New_DD
    int rk(0);
    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
    vecPatch_.patches_.clear();
    //vecPatch_.refHindex_ = rk;

    if (global_region) 
        decomposition_ = NULL;
    else
        decomposition_ = DomainDecompositionFactory::createGlobal( params );

    if ( (global_region) && (rk) )
        return;

    patch_ = PatchesFactory::create( params, smpi, decomposition_, vecPatch_.refHindex_ );
    patch_->setLocationAndAllocateFields( params, decomposition_, vecPatches );
    vecPatch_.patches_.push_back( patch_ );
    
    //vecPatch_.refHindex_ = vecPatches.refHindex_ / vecPatches.size();
    vecPatch_.updateFieldList( smpi );
    
    //vecPatch_.updateFieldList(0, smpi);
    vecPatch_.patches_[0]->finalizeMPIenvironment( params );
    vecPatch_.nrequests = vecPatches( 0 )->requests_.size();
    vecPatch_.nAntennas = vecPatch_( 0 )->EMfields->antennas.size();
    vecPatch_.initExternals( params );
    if (!params.apply_rotational_cleaning)
        vecPatch_.applyExternalFields();
    
    fake_patch = PatchesFactory::clone(vecPatches(0), params, smpi, vecPatches.domain_decomposition_, 0, 0, false);
    if (params.is_spectral)
        patch_->EMfields->saveMagneticFields( true );
        
}


void Region::coupling( Params &params, bool global_region )
{
    vecPatch_( 0 )->EMfields->MaxwellAmpereSolver_->coupling( params, vecPatch_( 0 )->EMfields, global_region );
    if ( ( params.geometry == "AMcylindrical" ) && ( global_region ) )
        vecPatch_( 0 )->EMfields->MaxwellAmpereSolver_->rotational_cleaning( vecPatch_( 0 )->EMfields );
    coupled_ = true;
}

Region::~Region()
{
}

void Region::clean()
{
    if ( (vecPatch_.patches_.size()) && (coupled_) )
        vecPatch_( 0 )->EMfields->MaxwellAmpereSolver_->uncoupling();

    if( diag_ !=NULL ) {
        diag_->closeFile();
        delete diag_;
    }
    if( patch_!=NULL ) {
        delete patch_;
    }
    if( decomposition_ !=NULL ) {
        delete decomposition_;
    }
    
    if( fake_patch!=NULL ) {
        delete fake_patch;
    }
    
}

void Region::solveMaxwell( Params &params, SimWindow *simWindow, int itime, double time_dual, Timers &timers, SmileiMPI *smpi )
{
    vecPatch_.solveMaxwell( params, simWindow, itime, time_dual, timers, smpi );
    
    // current no more used for now, reinitialized for next timestep
    //vecPatch_.resetRhoJ();

}

void Region::identify_additional_patches(SmileiMPI* smpi, VectorPatch& vecPatches, Params& params, SimWindow* simWindow)
{
    double delta_moving_win( simWindow->getNmoved()*params.cell_length[0] );
    
    //cout << "patch_->getDomainLocalMin(0) = " << patch_->getDomainLocalMin(0) << endl;
    //cout << "patch_->getDomainLocalMax(0) = " << patch_->getDomainLocalMax(0) << endl;
    //cout << "patch_->getDomainLocalMin(1) = " << patch_->getDomainLocalMin(1) << endl;
    //cout << "patch_->getDomainLocalMax(1) = " << patch_->getDomainLocalMax(1) << endl;
    //cout << "patch_->getDomainLocalMin(2) = " << patch_->getDomainLocalMin(2) << endl;
    //cout << "patch_->getDomainLocalMax(2) = " << patch_->getDomainLocalMax(2) << endl;

    // vecPatch_.patches_
    double domain_min( 10 );
    double domain_max( 0 );

    for ( unsigned int ipatch = 0 ; ipatch < vecPatches.size() ; ipatch++ ) {
        
        bool patch_is_in( true );
        for ( unsigned int iDim = 0 ; iDim < params.nDim_field ; iDim++ ) {
            if ( vecPatch_.patches_.size() ) {
                domain_min = patch_->getDomainLocalMin(iDim);
                domain_max = patch_->getDomainLocalMax(iDim);
            }
            double center = ( vecPatches(ipatch)->getDomainLocalMin(iDim) + vecPatches(ipatch)->getDomainLocalMax(iDim) - 2. * delta_moving_win *(iDim==0)  ) /2.;
            //cout << iDim << " " << domain_min << " " << center << " " <<domain_max << endl;
            if ( ( center < domain_min ) || ( center > domain_max ) )
                patch_is_in = false;
        }
        if (!patch_is_in) {
            additional_patches_.push_back( vecPatches(ipatch)->hindex );
            additional_patches_ranks.push_back( hrank_global_region( vecPatches(ipatch)->hindex, params, vecPatches ) );
        }
    }
    
    //cout << smpi->getRank() << " - additional : ";
    //for (int i = 0 ; i < additional_patches_.size() ; i++)
    //    cout << additional_patches_[i] << " " ;
    //cout << endl;
    //cout << smpi->getRank() << " - additional on : ";
    //for (int i = 0 ; i < additional_patches_.size() ; i++)
    //    cout << additional_patches_ranks[i] << " " ;
    //cout << endl;

}

int Region::hrank_global_region( int hindex, Params& params, VectorPatch& vecPatches )
{
    int rank(0);

    if (decomposition_!=NULL) { // real double decomposition
        std::vector<unsigned int> patch_coordinates = vecPatches.domain_decomposition_->getDomainCoordinates( hindex );
        std::vector<int> rank_coordinates;
        //rank_coordinates.resize( params.nDim_field );
        rank_coordinates.resize( 3, 0 );

        for ( unsigned int iDim = 0 ; iDim < params.nDim_field ; iDim++ ) {
            int min =  patch_coordinates[iDim]    * params.n_space[iDim];
            int max = (patch_coordinates[iDim]+1) * params.n_space[iDim];
            int center = (min+max)/2;

            int idomain(0);
            while ( ( center > params.offset_map[iDim][idomain] ) && ( idomain < (int)params.number_of_region[iDim] ) )
                idomain++;

            rank_coordinates[iDim] = idomain-1;

        }
        rank = params.map_rank[ rank_coordinates[0] ][ rank_coordinates[1] ][ rank_coordinates[2] ];
    }
    else { // Gather everything on MPI 0 
        rank = 0;
    }

    return rank;
}

void Region::identify_missing_patches(SmileiMPI* smpi, VectorPatch& vecPatches, Params& params)
{
    if ( (decomposition_==NULL) && (smpi->getRank()) )
        return;
    //missing_patches_.push_back()

    // Loop on theroritical patches

    std::vector<int> patch_min_coord(patch_->getDomainLocalMin().size());
    std::vector<int> patch_max_coord(patch_->getDomainLocalMin().size());
    int npatch_domain(1);

    if (decomposition_!=NULL) { // SDMD
        for ( unsigned int iDim = 0 ; iDim < patch_->getDomainLocalMin().size() ; iDim++ ) {
            patch_min_coord[iDim] = params.offset_map[iDim][ patch_->Pcoordinates[iDim] ] / params.n_space[iDim];
            if ( patch_->Pcoordinates[iDim] < params.number_of_region[iDim]-1 )
                patch_max_coord[iDim] = params.offset_map[iDim][ patch_->Pcoordinates[iDim]+1 ] / params.n_space[iDim] - 1;
            else
                patch_max_coord[iDim] = params.n_space_global[iDim] / params.n_space[iDim] - 1;

            npatch_domain *= params.n_space_region[iDim] / params.n_space[iDim];
        }
    }
    else { // Global mode for rotational cleaning
        for ( unsigned int iDim = 0 ; iDim < patch_->getDomainLocalMin().size() ; iDim++ ) {
            patch_min_coord[iDim] = (int)( patch_->getDomainLocalMin(iDim) / params.cell_length[iDim] / (double)params.n_space[iDim] );
            patch_max_coord[iDim] = (int)( patch_->getDomainLocalMax(iDim) / params.cell_length[iDim] / (double)params.n_space[iDim] ) - 1;
            npatch_domain *= params.n_space_region[iDim] / params.n_space[iDim];
        }
    }

    //cout << patch_min_coord[0] << " " << patch_min_coord[1] << endl;
    //cout << patch_max_coord[0] << " " << patch_max_coord[1] << endl;
    
    //MPI_Alltoall or gather/bcast (npatch_domain);

    int patch_min_id = vecPatches.domain_decomposition_->getDomainId( patch_min_coord );
    int patch_max_id = vecPatches.domain_decomposition_->getDomainId( patch_max_coord );
    //cout << patch_min_id << " " << patch_max_id << "\t npatche theoritical : " << npatch_domain << endl; 

    int hmin = min (patch_min_id,patch_max_id); 
    int hmax = max (patch_min_id,patch_max_id); 


    if (params.nDim_field==1) {
        ERROR("Not mimplemented");
    }
    else if (params.nDim_field==2) {
        for ( int ix = patch_min_coord[0] ; ix <= patch_max_coord[0] ; ix++ ) {
            for ( int iy = patch_min_coord[1] ; iy <= patch_max_coord[1] ; iy++ ) {
                std::vector<int> coords(patch_->getDomainLocalMin().size());
                coords[0] = ix;
                coords[1] = iy;
                int hindex = vecPatches.domain_decomposition_->getDomainId( coords );
                if (hindex < hmin)
                    hmin = hindex;
                if (hindex > hmax)
                    hmax = hindex;
                
            }
        }
    }
    else if (params.nDim_field==3) {
        for ( int ix = patch_min_coord[0] ; ix <= patch_max_coord[0] ; ix++ ) {
            for ( int iy = patch_min_coord[1] ; iy <= patch_max_coord[1] ; iy++ ) {
                for ( int iz = patch_min_coord[2] ; iz <= patch_max_coord[2] ; iz++ ) {
                    std::vector<int> coords(patch_->getDomainLocalMin().size());
                    coords[0] = ix;
                    coords[1] = iy;
                    coords[2] = iz;
                    int hindex = vecPatches.domain_decomposition_->getDomainId( coords );
                    if (hindex < hmin)
                        hmin = hindex;
                    if (hindex > hmax)
                        hmax = hindex;
                }
            }
        }
    }

    for ( int idx = hmin ; idx <=hmax ; idx++ ) {
        
        // test if (idx is in domain)
        bool patch_is_in( true );

        vector<unsigned int> Pcoordinates = vecPatches.domain_decomposition_->getDomainCoordinates( idx );
        for ( unsigned int iDim = 0 ; iDim < patch_->getDomainLocalMin().size() ; iDim++ ) {
            double center = ( 2. * (double)Pcoordinates[iDim] + 1. )* (double)params.n_space[iDim]*params.cell_length[iDim] / 2.;
            
            if ( ( center < patch_->getDomainLocalMin(iDim) ) || ( center > patch_->getDomainLocalMax(iDim) ) ) {
                patch_is_in = false;
            }
        }
        if (!patch_is_in)
            continue;
        
        if ( (idx < (int)vecPatches(0)->hindex) || (idx > (int)vecPatches(vecPatches.size()-1)->hindex) ) {
            missing_patches_.push_back( idx );
            missing_patches_ranks.push_back(  smpi->hrank( idx ) );
        }
        else
            local_patches_.push_back( idx );
    }

    //cout << "Theoritical : " << hmin << " " << hmax << endl;
    //cout << smpi->getRank() << " - missing : " ;
    //for (int i = 0 ; i < missing_patches_.size() ; i++)
    //    cout << missing_patches_[i] << " " ;
    //cout << endl;
    //cout << smpi->getRank() << " - missing on : " ;
    //for (int i = 0 ; i < missing_patches_.size() ; i++)
    //    cout << missing_patches_ranks[i] << " " ;
    //cout << endl;
    // 
    //cout << smpi->getRank() << " - local : " ;
    //for (int i = 0 ; i < local_patches_.size() ; i++)
    //    cout << local_patches_[i] << " " ;
    //cout << endl;

}


void Region::reset_fitting(SmileiMPI* smpi, Params& params)
{
    // Reorder domain to minimize costs
    int max  = -1;
    int ref  = -1;
    int ref2 = -1;

    if (local_patches_.size()>0) {
        max = local_patches_.size();
        ref = smpi->getRank();
    }

    // Quel domaine je récupère en fonction de mes patchs présents sur le process mais non locaux au sens décomposition ?
    //    a qui appartiennent les patchs surnuméraires ?
    //    un proc possède t il plus de matchs que moi ?
    for (int irk=0;irk<smpi->getSize();irk++) {
        int count = std::count( additional_patches_ranks.begin(), additional_patches_ranks.end(), irk );
        if (count==0) continue;
        if (count>max ) {
            ref2 = ref;
            ref  = irk;
            max = count;
        }
        else if (count==max) {
            if ( irk%2 == smpi->getRank()%2 ) {
                ref2 = ref;
                ref = irk;
            }
            else {
                ref2 = irk;
            }
                
        }
    }
    //cout << "Target : " << ref << " or " << ref2 << "\n";

    int local_map[2];
    local_map[0] = ref;
    local_map[1] = ref2;

    int global_map[2*smpi->getSize()];

    MPI_Allgather(local_map, 2, MPI_INT,
               global_map, 2, MPI_INT,
               MPI_COMM_WORLD);

    int target_map[smpi->getSize()];
    for (int i=0 ; i< smpi->getSize() ;i++)
        target_map[i] = -1;

     //cout << "GLOBAL MAP" <<endl;
     //for (int i=0 ; i< smpi->getSize() ;i++)
     //    cout << global_map[2*i] <<" " ;
     //cout << endl;
     //for (int i=0 ; i< smpi->getSize() ;i++)
     //    cout << global_map[2*i+1] <<" " ;
     //cout << endl;

    //local
    for (int i=0 ; i< smpi->getSize() ;i++) {
        if (global_map[2*i] == i) {
            target_map[i] = global_map[2*i];
            for (int j=0 ; j< smpi->getSize() ;j++) {
                if (j==i)continue;
                if ( global_map[2*j] == target_map[i] )
                    global_map[2*j] = -1;
                else if ( global_map[2*j+1] == target_map[i] )
                    global_map[2*j+1] = -1;
            }
        }
        else if  (global_map[2*i+1] == i) {
            target_map[i] = global_map[2*i+1];
            for (int j=0 ; j< smpi->getSize() ;j++) {
                if (j==i)continue;
                if ( global_map[2*j] == target_map[i] )
                    global_map[2*j] = -1;
                else if ( global_map[2*j+1] == target_map[i] )
                    global_map[2*j+1] = -1;
            }
        }
    }

     //cout << "TARGET MAP - STEP 1 (local) " << endl;
     //for (int i=0 ; i< smpi->getSize() ;i++)
     //    cout << target_map[i] <<" " ;
     //cout << endl;


    // single
    for (int i=0 ; i< smpi->getSize() ;i++) {
        if (global_map[2*i] == -1) {
            target_map[i] = global_map[2*i+1];
            for (int j=0 ; j< smpi->getSize() ;j++) {
                if (j==i)continue;
                if ( global_map[2*j] == target_map[i] )
                    global_map[2*j] = -1;
                else if ( global_map[2*j+1] == target_map[i] )
                    global_map[2*j+1] = -1;
            }
        }
        else if  (global_map[2*i+1] == -1) {
            target_map[i] = global_map[2*i];
            for (int j=0 ; j< smpi->getSize() ;j++) {
                if (j==i)continue;
                if ( global_map[2*j] == target_map[i] )
                    global_map[2*j] = -1;
                else if ( global_map[2*j+1] == target_map[i] )
                    global_map[2*j+1] = -1;
            }
        }
    }

     //cout << "TARGET MAP - STEP 2 (single) " << endl;
     //for (int i=0 ; i< smpi->getSize() ;i++)
     //    cout << target_map[i] <<" " ;
     //cout << endl;

     // facor first
     for (int i=0 ; i< smpi->getSize() ;i++) {
         if ( (target_map[i]!=-1) && (global_map[2*i]==-1) ) continue;
         target_map[i] = global_map[2*i];
         global_map[2*i+1] = -1;
         for (int j=0 ; j< smpi->getSize() ;j++) {
             if (j==i)continue;
             if ( global_map[2*j] == target_map[i] )
                 global_map[2*j] = -1;
             else if ( global_map[2*j+1] == target_map[i] )
                 global_map[2*j+1] = -1;
         }
     }
     //cout << "TARGET MAP - STEP 2.5 (favor first) " << endl;
     //for (int i=0 ; i< smpi->getSize() ;i++)
     //    cout << target_map[i] <<" " ;
     //cout << endl;


    // Unattribuated
    for (int i=0 ; i< smpi->getSize() ;i++) {
        //int todo = -1;
        bool att(false);
        for (int j=0 ; j< smpi->getSize() ;j++) { //0 1 2 -1 4 -1
            if ( target_map[j] == i ) {
                att = true;
            }
        }
        if (att == false) {
            for (int j=0 ; j< smpi->getSize() ;j++) {
                if (target_map[j] == -1 ){                    
                    target_map[j] = i;
                    break;
                }
            }
        }
    }

    //cout << "TARGET MAP - STEP 3 (unattribuated) " << endl;
    //for (int i=0 ; i< smpi->getSize() ;i++)
    //    cout << target_map[i] <<" " ;
    //cout << endl;


    int mpi_map[smpi->getSize()];
    for (int i=0 ; i< smpi->getSize() ;i++) {

        for (int j=0 ; j< smpi->getSize() ;j++) {
            if ( target_map[j] == i )
                mpi_map[i] = j;
        }

        
    }

    //cout << "TARGET MAP " << endl;
    //for (int i=0 ; i< smpi->getSize() ;i++)
    //    cout << target_map[i] <<" " ;
    //cout << endl;

    // Revert the map for future usage :
    //   - the computed one contains the domain id per MPI process
    //   - the necessary one will contain the MPI rank associated to the domain distributer through a linearized arrangement
     for (int i=0 ; i< smpi->getSize() ;i++)
         target_map[i] = mpi_map[i];


    std::vector< std::vector< std::vector<int> > > new_map_rank;
    new_map_rank.resize( params.number_of_region[0] );
    for ( unsigned int iDim = 0 ; iDim < params.number_of_region[0] ; iDim++ ) {
        new_map_rank[iDim].resize( params.number_of_region[1] );
        for ( unsigned int jDim = 0 ; jDim < params.number_of_region[1] ; jDim++ ) {
            new_map_rank[iDim][jDim].resize( params.number_of_region[2], -1 );
        }
    }

    //for ( int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ ) {
    //    for ( int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ ) {
    //        cout << params.map_rank[xDom][yDom][0] << " ";
    //    }
    //    cout << endl;
    //}


    //update params.map_rank[xDom][yDom][zDom]
    for ( unsigned int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ )
        for ( unsigned int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ ) {
            for ( unsigned int zDom = 0 ; zDom < params.number_of_region[2] ; zDom++ ) {
                for (int i=0 ; i< smpi->getSize() ;i++)
                    if ( params.map_rank[xDom][yDom][zDom] == i ) {
                        new_map_rank[xDom][yDom][zDom] = target_map[i];
                    }
            }
        }


    //cout << "Befor : " << vecPatch_.refHindex_ << " -" << params.coordinates[0] << " " << params.coordinates[1] << " " << params.coordinates[2] << endl;
    //cout << "Befor : " << vecPatch_.refHindex_ << " -" << params.coordinates[0] << " " << params.coordinates[1] << " " << endl;


    // Define the domain id of the domain hosted per current MPI rank
    int targeted_rk = -1;
    for (int i=0 ; i< smpi->getSize() ;i++)
        if (target_map[i]==smpi->getRank())
            targeted_rk = i;

    if ( targeted_rk==-1 ) {
        MPI_Abort(MPI_COMM_WORLD,-1);
    }


    // Compute coordinates of current patch in 3D
    for ( unsigned int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ )
        for ( unsigned int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ )
            for ( unsigned int zDom = 0 ; zDom < params.number_of_region[2] ; zDom++ ) {

                if ( params.map_rank[xDom][yDom][zDom] ==  targeted_rk ) {
                    params.coordinates[0] = xDom;
                    params.coordinates[1] = yDom;
                    params.coordinates[2] = zDom;
                }
            }

    int count_ = 0;
    for ( unsigned int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ )
        for ( unsigned int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ ) {
            for ( unsigned int zDom = 0 ; zDom < params.number_of_region[2] ; zDom++ ) {
                params.map_rank[xDom][yDom][zDom] = new_map_rank[xDom][yDom][zDom];

                if ( new_map_rank[xDom][yDom][zDom] == smpi->getRank() )
                    vecPatch_.refHindex_ = count_; //domain Id are distributed linearly
                count_++;

            }
        }



    //for ( int xDom = 0 ; xDom < params.number_of_region[0] ; xDom++ ) {
    //    for ( int yDom = 0 ; yDom < params.number_of_region[1] ; yDom++ ) {
    //        cout << params.map_rank[xDom][yDom][0] << " ";
    //    }
    //    cout << endl;
    //}


    // Compute size of local domain
    for ( unsigned int iDim = 0 ; iDim < params.nDim_field ; iDim++ ) {
        if ( params.coordinates[iDim] != (int)params.number_of_region[iDim]-1 ) {
            params.n_space_region[iDim] = params.offset_map[iDim][params.coordinates[iDim]+1] - params.offset_map[iDim][params.coordinates[iDim]];
        }
        else {
            params.n_space_region[iDim] = params.n_space_global[iDim] - params.offset_map[iDim][params.coordinates[iDim]];
        }
    }

    //cout << "After : " << vecPatch_.refHindex_ << " -" << params.coordinates[0] << " " << params.coordinates[1] << endl;


}

void Region::solveEnvelope( Params &params, SimWindow *simWindow, int itime, double time_dual, Timers &timers, SmileiMPI *smpi )
{
    vecPatch_.solveEnvelope( params, simWindow, itime, time_dual, timers, smpi );
    
}

void Region::reset_mapping()
{
    additional_patches_.clear();
    additional_patches_ranks.clear();
    local_patches_.clear();
    missing_patches_.clear();
    missing_patches_ranks.clear();
}
