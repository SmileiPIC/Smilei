
#include "SmileiMPI.h"

#include <cmath>
#include <cstring>

#include <iostream>
#include <sstream>
#include <fstream>

#include "Params.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "Field.h"

#include "Species.h"
#include "Hilbert_functions.h"
#include "VectorPatch.h"

#include "Diagnostic.h"
#include "DiagnosticScalar.h"
#include "DiagnosticParticles.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI constructor :
//     - Call MPI_Init_thread, MPI_THREAD_MULTIPLE required
//     - Set MPI env
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI::SmileiMPI( int* argc, char*** argv )
{    
    int mpi_provided;

#ifdef _OPENMP
    MPI_Init_thread( argc, argv, MPI_THREAD_MULTIPLE, &mpi_provided );
    if (mpi_provided != MPI_THREAD_MULTIPLE){
        ERROR("MPI_THREAD_MULTIPLE not supported. Compile your MPI ibrary with THREAD_MULTIPLE support.");
    }
#else
    MPI_Init( argc, argv );
#endif

    SMILEI_COMM_WORLD = MPI_COMM_WORLD;
    MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
    MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );

} // END SmileiMPI::SmileiMPI


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI destructor :
//     - Call MPI_Finalize
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI::~SmileiMPI()
{
    delete[]periods_;

    MPI_Finalize();

} // END SmileiMPI::~SmileiMPI


// ---------------------------------------------------------------------------------------------------------------------
// Broadcast namelist in SMILEI_COMM_WORLD
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::bcast( string& val )
{
    int charSize=0;
    if (isMaster()) charSize = val.size()+1;
    MPI_Bcast(&charSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);

    char tmp[charSize];
    if (isMaster()) strcpy(tmp, val.c_str());
    MPI_Bcast(tmp, charSize, MPI_CHAR, 0, SMILEI_COMM_WORLD);

    if (!isMaster()) val=tmp;

} // END bcast( string )


// ---------------------------------------------------------------------------------------------------------------------
// Broadcast namelist
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::bcast( int& val )
{
    MPI_Bcast(&val, 1, MPI_INT, 0, SMILEI_COMM_WORLD);

} // END bcast( int ) in SMILEI_COMM_WORLD


// ---------------------------------------------------------------------------------------------------------------------
// Initialize MPI (per process) environment
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::init( Params& params )
{
    // Initialize patch environment 
    patch_count.resize(smilei_sz, 0);
    target_patch_count.resize(smilei_sz, 0);
    capabilities.resize(smilei_sz, 1);
    Tcapabilities = smilei_sz;

    // Initialize patch distribution
    init_patch_count(params);

    // Initialize buffers for particles push vectorization
    //     - 1 thread push particles for a unique patch at a given time
    //     - so 1 buffer per thread
#ifdef _OPENMP
    dynamics_Epart.resize(omp_get_max_threads());
    dynamics_Bpart.resize(omp_get_max_threads());
    dynamics_gf.resize(omp_get_max_threads());
    dynamics_iold.resize(omp_get_max_threads());
    dynamics_deltaold.resize(omp_get_max_threads());
#else
    dynamics_Epart.resize(1);
    dynamics_Bpart.resize(1);
    dynamics_gf.resize(1);
    dynamics_iold.resize(1);
    dynamics_deltaold.resize(1);
#endif

    // Set periodicity of the simulated problem
    periods_  = new int[params.nDim_field];
    for (unsigned int i=0 ; i<params.nDim_field ; i++) periods_[i] = 0;
    // Geometry periodic in x
    if (params.bc_em_type_x[0]=="periodic") {
        periods_[0] = 1;
        MESSAGE(1,"applied topology for periodic BCs in x-direction");
    }
    if (params.nDim_field>1) {
        // Geometry periodic in y
        if (params.bc_em_type_y[0]=="periodic") {
            periods_[1] = 1;
            MESSAGE(2,"applied topology for periodic BCs in y-direction");
        }
    }
    if (params.nDim_field>2) {
        // Geometry periodic in y
        if (params.bc_em_type_z[0]=="periodic") {
            periods_[2] = 1;
            MESSAGE(2,"applied topology for periodic BCs in z-direction");
        }
    }
} // END init


// ---------------------------------------------------------------------------------------------------------------------
//  Initialize patch distribution
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::init_patch_count( Params& params)
{

//#ifndef _NOTBALANCED
//    bool use_load_balancing(true);
//    if (!use_load_balancing) {
//        int Npatches = params.number_of_patches[0];
//        for (unsigned int i = 1; i < params.nDim_field; i++)
//            Npatches *=  params.number_of_patches[i]; // Total number of patches.
//        if (Npatches!=smilei_sz) ERROR("number of patches abd MPI processes");
//        for (unsigned int i=0; i<smilei_sz; i++) patch_count[i] = 1;
//        return;
//    }
//#endif
    
    unsigned int Npatches, r, Ncur, Pcoordinates[3], ncells_perpatch;
    double Tload,Tcur, Lcur, total_load, local_load, local_load_temp, above_target, below_target;
    
    unsigned int tot_species_number = PyTools::nComponents("Species");
    
    // Define capabilities here if not default.              
    //Capabilities of devices hosting the different mpi processes. All capabilities are assumed to be equal for the moment.
    //Compute total capability: Tcapabilities. Uncomment if cpability != 1 per MPI rank
    //Tcapabilities = 0;
    //for (unsigned int i = 0; i < smilei_sz; i++)
    //    Tcapabilities += capabilities[i];
    
    //Compute target load: Tload = Total load * local capability / Total capability.
    
    // Some initialization of the box parameters
    Npatches = params.tot_number_of_patches;
    ncells_perpatch = 1;
    vector<double> cell_xmin(3,0.), cell_xmax(3,1.), cell_dx(3,2.), x_cell(3,0);
    for (unsigned int i = 0; i < params.nDim_field; i++) {
        ncells_perpatch *= params.n_space[i]+2*params.oversize[i];
        if (params.cell_length[i]!=0.) cell_dx[i] = params.cell_length[i];
    }
    
    // First, distribute all patches evenly
    unsigned int Npatches_local = Npatches / smilei_sz, FirstPatch_local;
    unsigned int remainder = Npatches % smilei_sz;
    if( smilei_rk < remainder ) {
        Npatches_local++;
        FirstPatch_local = Npatches_local * smilei_rk;
    } else {
        FirstPatch_local = Npatches_local * smilei_rk + remainder;
    }
//    // Test
//    int tot, loc=Npatches_local;
//    MPI_Allreduce( &loc, &tot, 1, MPI_INT, MPI_SUM, SMILEI_COMM_WORLD );
//    if( tot != Npatches ) ERROR("Npatches should be "<<Npatches<<" but it is "<<tot);
    
    // Second, prepare the profiles for each species
    vector<Profile*> densityProfiles(0), ppcProfiles(0);
    for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++){
        std::string species_type("");
        PyTools::extract("species_type",species_type,"Species",ispecies);
        PyObject *profile1;
        std::string densityProfileType("");
        bool ok1 = PyTools::extract_pyProfile("nb_density"    , profile1, "Species", ispecies);
        bool ok2 = PyTools::extract_pyProfile("charge_density", profile1, "Species", ispecies);
        if( ok1 ) densityProfileType = "nb";
        if( ok2 ) densityProfileType = "charge";
        densityProfiles.push_back(new Profile(profile1, params.nDim_particle, densityProfileType+"_density "+species_type));
        PyTools::extract_pyProfile("n_part_per_cell", profile1, "Species", ispecies);
        ppcProfiles.push_back(new Profile(profile1, params.nDim_particle, "n_part_per_cell "+species_type));
    }
    
    // Third, loop over local patches to obtain their approximate load
    vector<double> PatchLoad (Npatches_local, 0.);
    for(unsigned int ipatch=0; ipatch<Npatches_local; ipatch++){
        // Get patch coordinates
        unsigned int hindex = FirstPatch_local + ipatch;
        generalhilbertindexinv(params.mi[0], params.mi[1], params.mi[2], &Pcoordinates[0], &Pcoordinates[1], &Pcoordinates[2], hindex);
        for (unsigned int i=0 ; i<params.nDim_field ; i++) {
            Pcoordinates[i] *= params.n_space[i];
            if (params.cell_length[i]!=0.) {
                cell_xmin[i] = (Pcoordinates[i]+0.5)*params.cell_length[i];
                cell_xmax[i] = (Pcoordinates[i]+0.5+params.n_space[i])*params.cell_length[i];
            }
        }
        //Accumulate particles load of the current patch
        for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++){
            local_load = 0.;
            
            // This commented block loops through all cells of the current patch to calculate the load
            //for (x_cell[0]=cell_xmin[0]; x_cell[0]<cell_xmax[0]; x_cell[0]+=cell_dx[0]) {
            //    for (x_cell[1]=cell_xmin[1]; x_cell[1]<cell_xmax[1]; x_cell[1]+=cell_dx[1]) {
            //        for (x_cell[2]=cell_xmin[2]; x_cell[2]<cell_xmax[2]; x_cell[2]+=cell_dx[2]) {
            //            double n_part_in_cell = floor(ppcProfiles[ispecies]->valueAt(x_cell));
            //            if( n_part_in_cell<=0. ) continue;
            //            else if( densityProfiles[ispecies]->valueAt(x_cell)==0. ) continue;
            //            local_load += n_part_in_cell;
            //        }
            //    }
            //}
            // Instead of looping all cells, the following takes only the central point (much faster)
            for (unsigned int i=0 ; i<params.nDim_field ; i++) {
                if (params.cell_length[i]==0.) x_cell[i] = 0.;
                else x_cell[i] = 0.5*(cell_xmin[i]+cell_xmax[i]);
            }
            double n_part_in_cell = floor(ppcProfiles[ispecies]->valueAt(x_cell));
            if( n_part_in_cell && densityProfiles[ispecies]->valueAt(x_cell)!=0.)
                local_load += n_part_in_cell * ncells_perpatch;
            
            // Consider whether this species is frozen
            double time_frozen(0.);
            PyTools::extract("time_frozen",time_frozen ,"Species",ispecies);
            if(time_frozen > 0.) local_load *= params.coef_frozen;
            // Add the load of the species to the current patch load
            PatchLoad[ipatch] += local_load;
        }
        //Add grid contribution to the load.
        PatchLoad[ipatch] += ncells_perpatch*params.coef_cell;
    }
    densityProfiles.resize(0); densityProfiles.clear();
    ppcProfiles.resize(0); ppcProfiles.clear();
    
    // Fourth, the arrangement of patches is balanced
    
    // Initialize loads
    total_load = Npatches*ncells_perpatch*params.coef_cell ; // We assume the load of one cell to be equal to coef_cell and account for ghost cells.
    Tload = total_load/Tcapabilities; //Target load for each mpi process.
    Tcur = Tload * capabilities[0];  //Init.
    r = 0;  //Start by finding work for rank 0.
    Ncur = 0; // Number of patches assigned to current rank r.
    Lcur = 0.; //Load assigned to current rank r.
    
    // MPI master loops patches and figures the best arrangement
    if( smilei_rk==0 ) {
        unsigned int rk = 0;
        MPI_Status status;
        while( true ) { // loop cpu ranks
            unsigned int hindex = 0;
            for(unsigned int ipatch=0; ipatch < Npatches_local; ipatch++){
                local_load = PatchLoad[ipatch];
                Lcur += local_load; //Add grid contribution to the load.
                Ncur++; // Try to assign current patch to rank r.
                
                //if (isMaster()) cout <<"h= " << hindex << " Tcur = " << Tcur << " Lcur = " << Lcur <<" Ncur = " << Ncur <<" r= " << r << endl;
                if (r < (unsigned int)smilei_sz-1){
                    
                    if ( Lcur > Tcur || smilei_sz-r >= Npatches-hindex){ //Load target is exceeded or we have as many patches as procs left.
                        above_target = Lcur - Tcur;  //Including current patch, we exceed target by that much.
                        below_target = Tcur - (Lcur-local_load); // Excluding current patch, we mis the target by that much.
                        if((above_target > below_target) && (Ncur!=1)) { // If we're closer to target without the current patch...
                            patch_count[r] = Ncur-1;      // ... include patches up to current one.
                            Ncur = 1;
                            //Lcur = local_load;
                        } else {                          //Else ...
                            patch_count[r] = Ncur;        //...assign patches including the current one.
                            Ncur = 0;
                            //Lcur = 0.;
                        }
                        r++; //Move on to the next rank.
                        //Tcur = Tload * capabilities[r];  //Target load for current rank r.
                        Tcur += Tload * capabilities[r];  //Target load for current rank r.
                    }
                }// End if on r.
                hindex++;
            }// End loop on patches for rank rk
            patch_count[smilei_sz-1] = Ncur; // the last MPI process takes what's left.
            
            // Go to next rank
            rk++;
            if( rk >= smilei_sz ) break;
            
            // Get the load of patches pre-calculated by the next rank
            if( rk == remainder ) {
                Npatches_local--;
                PatchLoad.resize(Npatches_local);
            }
            MPI_Recv(&PatchLoad[0], Npatches_local, MPI_DOUBLE, rk, rk, SMILEI_COMM_WORLD, &status);
        }
        
        // The master cpu also writes the patch count to the file
        ofstream fout;
        fout.open ("patch_load.txt");
        fout << "Total load = " << total_load << endl;
        for (rk=0; rk<smilei_sz; rk++)
            fout << "patch count = " << patch_count[rk]<<endl;
        fout.close();
        
    // The other MPIs send their pre-calculated information
    } else {
        MPI_Send(&PatchLoad[0], Npatches_local, MPI_DOUBLE, 0, smilei_rk, SMILEI_COMM_WORLD);
    }
    
    // Lastly, the patch count is broadcast to all ranks
    MPI_Bcast( &patch_count[0], smilei_sz, MPI_INT, 0, SMILEI_COMM_WORLD);
    
} // END init_patch_count


// ---------------------------------------------------------------------------------------------------------------------
//  Recompute patch distribution
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::recompute_patch_count( Params& params, VectorPatch& vecpatches, double time_dual )
{

    unsigned int Npatches, r,Ncur,ncells_perpatch, Lmin, Lmin1, Lmin2, Lmax;
    double Tload,Tcur, Lcur, above_target, below_target, cells_load;
    unsigned int npatchmin =1;
    //Load of a cell = coef_cell*load of a particle.
    //Load of a frozen particle = coef_frozen*load of a particle.
    std::vector<double> Lp,Lp_global;
    int recv_counts[smilei_sz];
    ofstream fout;

    if (isMaster()) {
        fout.open ("patch_load.txt", std::ofstream::out | std::ofstream::app);
    }
    
    Npatches = params.tot_number_of_patches;
    
    ncells_perpatch = params.n_space[0]+2*params.oversize[0]; //Initialization
    for (unsigned int idim = 1; idim < params.nDim_field; idim++)
        ncells_perpatch *= params.n_space[idim]+2*params.oversize[idim];
 
    unsigned int tot_species_number = vecpatches(0)->vecSpecies.size();
    cells_load = ncells_perpatch*params.coef_cell ;

    Lp.resize(patch_count[smilei_rk], cells_load);
    Lp_global.resize(Npatches,0.);

    Tload = 0.;
    r = 0;  //Start by finding work for rank 0.
    Ncur = 0; // Number of patches assigned to current rank r.
    Lcur = 0.; //Load assigned to current rank r.

    //Compute Local Loads of each Patch (Lp)
    for(unsigned int ipatch=0; ipatch < (unsigned int)patch_count[smilei_rk]; ipatch++){
        for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++) {
            Lp[ipatch] += vecpatches(ipatch)->vecSpecies[ispecies]->getNbrOfParticles()*(1+(params.coef_frozen-1)*(time_dual < vecpatches(ipatch)->vecSpecies[ispecies]->time_frozen)) ;
        }
    }

    //Allgatherv loads of all patches in Lp_global
  
    recv_counts[0] = 0;
    for(int i=1; i < smilei_sz ; i++) recv_counts[i] = recv_counts[i-1]+patch_count[i-1];

    MPI_Allgatherv(&Lp[0],patch_count[smilei_rk],MPI_DOUBLE,&Lp_global[0], &patch_count[0], recv_counts, MPI_DOUBLE,MPI_COMM_WORLD);

    //Compute total loads
    for(unsigned int ipatch=0; ipatch < Npatches; ipatch++) Tload += Lp_global[ipatch];
    Tload /= Tcapabilities; //Target load for each mpi process.
    Tcur = Tload * capabilities[0];  //Init.

    //Loop over all patches to determine target_patch_count.
    for(unsigned int ipatch=0; ipatch < Npatches; ipatch++){

        Lcur += Lp_global[ipatch]; 
        Ncur++; // Try to assign current patch to rank r.

        if (r < (unsigned int)smilei_sz-1){

            if ( Lcur > Tcur || smilei_sz-r >= Npatches-ipatch){ //Load target is exceeded or we have as many patches as procs left.
                above_target = Lcur - Tcur;  //Including current patch, we exceed target by that much.
                below_target = Tcur - (Lcur-Lp_global[ipatch]); // Excluding current patch, we mis the target by that much.
                if((above_target > below_target) && (Ncur!=1)) { // If we're closer to target without the current patch...
                    target_patch_count[r] = Ncur-1;      // ... include patches up to current one.
                    Ncur = 1;
                } else {                          //Else ...
                    target_patch_count[r] = Ncur;        //...assign patches including the current one.
                    Ncur = 0;
                }
                r++; //Move on to the next rank.
                Tcur += Tload * capabilities[r];  //Target load for current rank r.
            } 
        }// End if on r.
        if (ipatch == Npatches-1){
            target_patch_count[smilei_sz-1] = Ncur; //When we reach the last patch, the last MPI process takes what's left.
        }
    }// End loop on patches.


    //Make sure the new patch_count is not too different from the previous one.
    // First patch
    Ncur = 0;                           //Sold
    Lmin1 = npatchmin;
    Lmin2 = 1;
    Lmin = std::max(Lmin1, Lmin2);      //Pmin
    Lmax = patch_count[0] - npatchmin;  //Pmax
    Tcur = 0;                           //Plast 

    //Loop
    for(unsigned int i=0; i< (unsigned int)smilei_sz-1; i++){

        Lmin2 += patch_count[i];         // futur Pmin
        Tcur += target_patch_count[i];   //Plast
        Lmax += patch_count[i+1];        //Pmax
 
        if (Tcur < Lmin){                      
            patch_count[i] = Lmin - Ncur;
        } else if (Tcur > Lmax ){                      
            patch_count[i] = Lmax - Ncur;
        } else {
            patch_count[i] = Tcur-Ncur;
        }
        Ncur += patch_count[i];           //new Sold
        Lmin1 = Ncur + npatchmin;
        Lmin = std::max(Lmin1, Lmin2);    //new Pmin

    }

    //Last patch
    patch_count[smilei_sz-1] = Npatches-Ncur;

    //Write patch_load.txt
    if (smilei_rk==0) {
        fout << "\tt = " << time_dual << endl;
        for (int irk=0;irk<smilei_sz;irk++)
            fout << " patch_count[" << irk << "] = " << patch_count[irk] << " target patch_count = "<< target_patch_count[irk] << endl;
        fout.close();
    }

    return;

} // END recompute_patch_count


// ----------------------------------------------------------------------
// Returns the rank of the MPI process currently owning patch h.
// ----------------------------------------------------------------------
int SmileiMPI::hrank(int h)
{
    if (h == MPI_PROC_NULL) return MPI_PROC_NULL;

    int patch_counter,rank;
    rank=0;
    patch_counter = patch_count[0];
    while (h >= patch_counter) {
        rank++;
        patch_counter += patch_count[rank];
    }
    return rank;
} // END hrank


// ----------------------------------------------------------------------
// Create MPI type to exchange all particles properties of particles
// ----------------------------------------------------------------------
MPI_Datatype SmileiMPI::createMPIparticles( Particles* particles )
{
    int nbrOfProp = particles->double_prop.size() + particles->short_prop.size() + particles->uint_prop.size();

    MPI_Aint address[nbrOfProp];
    for ( unsigned int iprop=0 ; iprop<particles->double_prop.size() ; iprop++ )
        MPI_Get_address( &( (*(particles->double_prop[iprop]))[0] ), &(address[iprop]) );
    for ( unsigned int iprop=0 ; iprop<particles->short_prop.size() ; iprop++ )
        MPI_Get_address( &( (*(particles->short_prop[iprop]))[0] ), &(address[particles->double_prop.size()+iprop]) );
    for ( unsigned int iprop=0 ; iprop<particles->uint_prop.size() ; iprop++ )
        MPI_Get_address( &( (*(particles->uint_prop[iprop]))[0] ), &(address[particles->double_prop.size()+particles->short_prop.size()+iprop]) );

    int nbr_parts[nbrOfProp];
    // number of elements per property
    for (int i=0 ; i<nbrOfProp ; i++)
        nbr_parts[i] = particles->size();

    MPI_Aint disp[nbrOfProp];
    // displacement between 2 properties
    disp[0] = 0;
    for (int i=1 ; i<nbrOfProp ; i++)
        disp[i] = address[i] - address[0];

    MPI_Datatype partDataType[nbrOfProp];
    // define MPI type of each property, default is DOUBLE
    for ( unsigned int i=0 ; i<particles->double_prop.size() ; i++)
        partDataType[i] = MPI_DOUBLE;
    for ( unsigned int iprop=0 ; iprop<particles->short_prop.size() ; iprop++ )
        partDataType[ particles->double_prop.size()+iprop] = MPI_SHORT;
    for ( unsigned int iprop=0 ; iprop<particles->uint_prop.size() ; iprop++ )
        partDataType[ particles->double_prop.size()+particles->short_prop.size()+iprop] = MPI_UNSIGNED;

    MPI_Datatype typeParticlesMPI;
    MPI_Type_create_struct( nbrOfProp, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
    MPI_Type_commit( &typeParticlesMPI );
    
    return typeParticlesMPI;

} // END createMPIparticles


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// -----------------------------------------       PATCH SEND / RECV METHODS        ------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::isend(Patch* patch, int to, int tag, Params& params)
{
    //MPI_Request request;

    for (int ispec=0 ; ispec<(int)patch->vecSpecies.size() ; ispec++){
        isend( &(patch->vecSpecies[ispec]->bmax), to, tag+2*ispec+1 );
        if ( patch->vecSpecies[ispec]->getNbrOfParticles() > 0 ){
            patch->vecSpecies[ispec]->typePartSend = createMPIparticles( patch->vecSpecies[ispec]->particles );
            isend( patch->vecSpecies[ispec]->particles, to, tag+2*ispec, patch->vecSpecies[ispec]->typePartSend );
            MPI_Type_free( &(patch->vecSpecies[ispec]->typePartSend) );
        }
    }

    //! \todo Removed the following block because the probe particles are not exchanged
    //
    //// Send probes' particles
    //for ( int iprobe = 0 ; iprobe < (int)patch->probes.size() ; iprobe++ ) {
    //    isend( patch->probes[iprobe], to, tag+2*patch->vecSpecies.size()+iprobe, params.nDim_particle );
    //}

    // Count number max of comms :
    int maxtag = 2 * patch->vecSpecies.size() + (2+params.nDim_particle) * patch->probes.size();
    
    isend( patch->EMfields, to, maxtag );
    
} // END isend( Patch )


void SmileiMPI::recv(Patch* patch, int from, int tag, Params& params)
{
    int nbrOfPartsRecv;

    for (int ispec=0 ; ispec<(int)patch->vecSpecies.size() ; ispec++){
        //Receive bmax
        recv( &patch->vecSpecies[ispec]->bmax, from, tag+2*ispec+1 );
        //Reconstruct bmin from bmax
        memcpy(&(patch->vecSpecies[ispec]->bmin[1]), &(patch->vecSpecies[ispec]->bmax[0]), (patch->vecSpecies[ispec]->bmax.size()-1)*sizeof(int) );
        patch->vecSpecies[ispec]->bmin[0]=0;
        //Prepare patch for receiving particles
        nbrOfPartsRecv = patch->vecSpecies[ispec]->bmax.back(); 
        //cout << smilei_rk << " recv " << nbrOfPartsRecv << endl;
        patch->vecSpecies[ispec]->particles->initialize( nbrOfPartsRecv, params.nDim_particle );
        //Receive particles
        if ( nbrOfPartsRecv > 0 ) {
            patch->vecSpecies[ispec]->typePartSend = createMPIparticles( patch->vecSpecies[ispec]->particles );
            recv( patch->vecSpecies[ispec]->particles, from, tag+2*ispec, patch->vecSpecies[ispec]->typePartSend );
            MPI_Type_free( &(patch->vecSpecies[ispec]->typePartSend) );
        }
    }
    
    // Removed next block because the probe particles are not exchanged.
    //
    //// Receive probes' particles
    //for ( int iprobe = 0 ; iprobe < (int)patch->probes.size() ; iprobe++ ) {
    //    recv( patch->probes[iprobe], from, tag+2*patch->vecSpecies.size()+iprobe, params.nDim_particle );
    //}

    // Count number max of comms :
    int maxtag = 2 * patch->vecSpecies.size() + (2+params.nDim_particle) * patch->probes.size();

    patch->EMfields->initAntennas(patch);
    recv( patch->EMfields, from, maxtag );


} // END recv ( Patch )


void SmileiMPI::isend(Particles* particles, int to, int tag, MPI_Datatype typePartSend)
{
    MPI_Request request ;
    MPI_Isend( &(particles->position(0,0)), 1, typePartSend, to, tag, MPI_COMM_WORLD, &request );

} // END isend( Particles )


void SmileiMPI::recv(Particles* particles, int to, int tag, MPI_Datatype typePartRecv)
{
    MPI_Status status;
    MPI_Recv( &(particles->position(0,0)), 1, typePartRecv, to, tag, MPI_COMM_WORLD, &status );

} // END recv( Particles )


// Assuming vec.size() is known (number of species). Asynchronous.
void SmileiMPI::isend(std::vector<int>* vec, int to, int tag)
{
    MPI_Request request; 
    MPI_Isend( &((*vec)[0]), (*vec).size(), MPI_INT, to, tag, MPI_COMM_WORLD, &request );

} // End isend ( bmax )


void SmileiMPI::recv(std::vector<int> *vec, int from, int tag)
{
    MPI_Status status;
    MPI_Recv( &((*vec)[0]), vec->size(), MPI_INT, from, tag, MPI_COMM_WORLD, &status );

} // End recv ( bmax )


void SmileiMPI::isend(ElectroMagn* EM, int to, int tag)
{
    isend( EM->Ex_, to, tag+0);
    isend( EM->Ey_, to, tag+1);
    isend( EM->Ez_, to, tag+2);
    isend( EM->Bx_, to, tag+3);
    isend( EM->By_, to, tag+4);
    isend( EM->Bz_, to, tag+5);
    isend( EM->Bx_m, to, tag+6);
    isend( EM->By_m, to, tag+7);
    isend( EM->Bz_m, to, tag+8);
    
    for (int antennaId=0 ; antennaId<(int)EM->antennas.size() ; antennaId++) {
        isend( EM->antennas[antennaId].field, to, tag+9+antennaId );
    }
    
    tag += 10 + EM->antennas.size();
    
    for (unsigned int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if(! EM->emBoundCond[bcId]) continue;
        
        for (unsigned int laserId=0 ; laserId < EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
            
            Laser * laser = EM->emBoundCond[bcId]->vecLaser[laserId];
            if( !(laser->spacetime[0]) && !(laser->spacetime[1]) ){
                LaserProfileSeparable* profile;
                profile = static_cast<LaserProfileSeparable*> ( laser->profiles[0] );
                if( ! profile->space_envelope ) continue;
                isend( profile->space_envelope, to , tag );
                isend( profile->phase, to, tag + 1);
                profile = static_cast<LaserProfileSeparable*> ( laser->profiles[1] );
                isend( profile->space_envelope, to , tag + 2 );
                isend( profile->phase, to, tag + 3);
                tag = tag + 4;
            }
        }
    }
} // End isend ( ElectroMagn )


void SmileiMPI::recv(ElectroMagn* EM, int from, int tag)
{
    recv( EM->Ex_, from, tag+0 );
    recv( EM->Ey_, from, tag+1 );
    recv( EM->Ez_, from, tag+2 );
    recv( EM->Bx_, from, tag+3 );
    recv( EM->By_, from, tag+4 );
    recv( EM->Bz_, from, tag+5 );
    recv( EM->Bx_m, from, tag+6 );
    recv( EM->By_m, from, tag+7 );
    recv( EM->Bz_m, from, tag+8 );
    
    for (int antennaId=0 ; antennaId<(int)EM->antennas.size() ; antennaId++) {
        recv( EM->antennas[antennaId].field, from, tag+9+antennaId );
    }
    
    tag += 10 + EM->antennas.size();
    
    for (unsigned int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if(! EM->emBoundCond[bcId]) continue;
        
        for (unsigned int laserId=0 ; laserId<EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
            Laser * laser = EM->emBoundCond[bcId]->vecLaser[laserId];
            if( !(laser->spacetime[0]) && !(laser->spacetime[1]) ){
                LaserProfileSeparable* profile;
                profile = static_cast<LaserProfileSeparable*> ( laser->profiles[0] );
                if( ! profile->space_envelope ) continue;
                recv( profile->space_envelope, from , tag );
                recv( profile->phase, from, tag + 1);
                profile = static_cast<LaserProfileSeparable*> ( laser->profiles[1] );
                recv( profile->space_envelope, from , tag + 2 );
                recv( profile->phase, from, tag + 3);
                tag = tag + 4;
            }
        }
    }
} // End recv ( ElectroMagn )


void SmileiMPI::isend(Field* field, int to, int hindex)
{
    MPI_Request request;
    MPI_Isend( &((*field)(0)),field->globalDims_, MPI_DOUBLE, to, hindex, MPI_COMM_WORLD, &request );

} // End isend ( Field )


void SmileiMPI::recv(Field* field, int from, int hindex)
{
    MPI_Status status;
    MPI_Recv( &((*field)(0)),field->globalDims_, MPI_DOUBLE, from, hindex, MPI_COMM_WORLD, &status );

} // End recv ( Field )

void SmileiMPI::isend( ProbeParticles* probe, int to, int tag, unsigned int nDim_particles )
{
    MPI_Request request; 
    // send offset
    MPI_Isend( &(probe->offset_in_file), 1, MPI_INT, to, tag, MPI_COMM_WORLD, &request );
    // send number of particles
    int nPart = probe->particles.size();
    MPI_Isend( &nPart, 1, MPI_INT, to, tag+1, MPI_COMM_WORLD, &request );
    // send particles
    if( nPart>0 )
        for( unsigned int i=0; i<nDim_particles; i++)
            MPI_Isend( &(probe->particles.Position[i][0]), nPart, MPI_DOUBLE, to, tag+1+i, MPI_COMM_WORLD, &request );

} // End isend ( probes )


void SmileiMPI::recv( ProbeParticles* probe, int from, int tag, unsigned int nDim_particles )
{
    MPI_Status status;
    // receive offset
    MPI_Recv( &(probe->offset_in_file), 1, MPI_INT, from, tag, MPI_COMM_WORLD, &status );
    // receive number of particles
    int nPart;
    MPI_Recv( &nPart, 1, MPI_INT, from, tag+1, MPI_COMM_WORLD, &status );
    // Resize particles
    probe->particles.initialize(nPart, nDim_particles);
    // receive particles
    if( nPart>0 )
        for( unsigned int i=0; i<nDim_particles; i++)
            MPI_Recv( &(probe->particles.Position[i][0]), nPart, MPI_DOUBLE, from, tag+1+i, MPI_COMM_WORLD, &status );

} // End recv ( probes )


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ------------------------------------------      DIAGS MPI SYNC     --------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Wrapper of MPI synchronization of all computing diags
//   - concerns    : scalars, phasespace, particles
//   - not concern : probes, fields, track particles (each patch write its own data)
//   - called in VectorPatch::runAllDiags(...)
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags(Diagnostic* diag, int timestep)
{
    if ( DiagnosticScalar* scalar = dynamic_cast<DiagnosticScalar*>( diag ) ) {
        computeGlobalDiags(scalar, timestep);
    } else if (DiagnosticParticles* particles = dynamic_cast<DiagnosticParticles*>( diag )) {
        computeGlobalDiags(particles, timestep);
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// MPI synchronization of scalars diags
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags(DiagnosticScalar* scalars, int timestep)
{
    
    if ( !(scalars->printNow(timestep))
      && !(scalars->timeSelection->theTimeIsNow(timestep)) ) return;
    
    vector<string>::iterator iterKey = scalars->out_key.begin();
    for(vector<double>::iterator iter = scalars->out_value.begin(); iter !=scalars->out_value.end(); iter++) {
        if ( ( (*iterKey).find("Min") == std::string::npos ) && ( (*iterKey).find("Max") == std::string::npos ) ) {
            MPI_Reduce(isMaster()?MPI_IN_PLACE:&((*iter)), &((*iter)), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        else if ( (*iterKey).find("MinCell") != std::string::npos ) {
            vector<double>::iterator iterVal = iter-1;
            val_index minVal;
            minVal.val   = (*iterVal);
            minVal.index = (*iter);
            MPI_Reduce(isMaster()?MPI_IN_PLACE:&minVal, &minVal, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
            if (isMaster()) {
                (*iterVal) = minVal.val;
                (*iter)    = minVal.index;
            }
        }
        else if ( (*iterKey).find("MaxCell") != std::string::npos ) {
            vector<double>::iterator iterVal = iter-1;
            val_index maxVal;
            maxVal.val   = (*iterVal);
            maxVal.index = (*iter);
            MPI_Reduce(isMaster()?MPI_IN_PLACE:&maxVal, &maxVal, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
            if (isMaster()) {
                (*iterVal) = maxVal.val;
                (*iter)    = maxVal.index;
            }
        }
        iterKey++;
    }

    if (isMaster()) {

        double Ukin = scalars->getScalar("Ukin");
        double Uelm = scalars->getScalar("Uelm");

        // added & lost energies due to the moving window
        double Ukin_out_mvw = scalars->getScalar("Ukin_out_mvw");
        double Ukin_inj_mvw = scalars->getScalar("Ukin_inj_mvw");
        double Uelm_out_mvw = scalars->getScalar("Uelm_out_mvw");
        double Uelm_inj_mvw = scalars->getScalar("Uelm_inj_mvw");
        
        // added & lost energies at the boundaries
        double Ukin_bnd = scalars->getScalar("Ukin_bnd");
        double Uelm_bnd = scalars->getScalar("Uelm_bnd");

        // total energy in the simulation
        double Utot = Ukin + Uelm;
        
        // total energy at time 0
        if (timestep==0) {
            scalars->Energy_time_zero  = Utot;
        }
        
        // the normalized energy balanced is normalized with respect to the current energy
        scalars->EnergyUsedForNorm = Utot;
        
        // expected total energy
        double Uexp = scalars->Energy_time_zero + Uelm_bnd + Ukin_inj_mvw + Uelm_inj_mvw
            -           ( Ukin_bnd + Ukin_out_mvw + Uelm_out_mvw );
        
        // energy balance
        double Ubal = Utot - Uexp;
        
        // normalized energy balance
        double Ubal_norm(0.);
        if (scalars->EnergyUsedForNorm>0.)
            Ubal_norm = Ubal / scalars->EnergyUsedForNorm;

        scalars->setScalar("Ubal_norm",Ubal_norm);
        scalars->setScalar("Ubal",Ubal);
        scalars->setScalar("Uexp",Uexp);
        scalars->setScalar("Utot",Utot);

    }
} // END computeGlobalDiags(DiagnosticScalar& scalars ...)


// ---------------------------------------------------------------------------------------------------------------------
// MPI synchronization of diags particles
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::computeGlobalDiags(DiagnosticParticles* diagParticles, int timestep)
{
    if (timestep - diagParticles->timeSelection->previousTime() == diagParticles->time_average-1) {
        MPI_Reduce(diagParticles->filename.size()?MPI_IN_PLACE:&diagParticles->data_sum[0], &diagParticles->data_sum[0], diagParticles->output_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if( !isMaster() ) diagParticles->clear();
    }
} // END computeGlobalDiags(DiagnosticParticles* diagParticles ...)
