
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

    MPI_Init_thread( argc, argv, MPI_THREAD_MULTIPLE, &mpi_provided );
    if (mpi_provided != MPI_THREAD_MULTIPLE){
        ERROR("MPI_THREAD_MULTIPLE not supported. Compile your MPI ibrary with THREAD_MULTIPLE support.");
    }

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
    MPI_Bcast(&tmp, charSize, MPI_CHAR, 0, SMILEI_COMM_WORLD);

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
    for (int i=0 ; i<params.nDim_field ; i++) periods_[i] = 0;
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

    unsigned int Npatches, r,Ncur,Pcoordinates[3],ncells_perpatch;
    double Tload,Tcur, Lcur, local_load, local_load_temp, above_target, below_target;
    std::vector<unsigned int> mincell,maxcell; //Min and max values of non empty cells for each species and in each dimension.
    vector<double> density_length(params.nDim_field,0.);
    //Load of a cell = coef_cell*load of a particle.
    //Load of a frozen particle = coef_frozen*load of a particle.
    ofstream fout;

    unsigned int tot_species_number = PyTools::nComponents("Species");
    mincell.resize(tot_species_number*3);
    maxcell.resize(tot_species_number*3);
       
    // Define capabilities here if not default.              
    //Capabilities of devices hosting the different mpi processes. All capabilities are assumed to be equal for the moment.
    //Compute total capability: Tcapabilities. Uncomment if cpability != 1 per MPI rank
    //Tcapabilities = 0;
    //for (unsigned int i = 0; i < smilei_sz; i++)
    //    Tcapabilities += capabilities[i];

    //Compute target load: Tload = Total load * local capability / Total capability.
    
    Tload = 0.;
    Npatches = params.tot_number_of_patches;
    
    ncells_perpatch = params.n_space[0]+2*params.oversize[0]; //Initialization
    for (unsigned int idim = 1; idim < params.nDim_field; idim++)
        ncells_perpatch *= params.n_space[idim]+2*params.oversize[idim];

    r = 0;  //Start by finding work for rank 0.
    Ncur = 0; // Number of patches assigned to current rank r.
    Lcur = 0.; //Load assigned to current rank r.


    for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++){

        //Needs to be updated when dens_lenth is a vector in params.

        // Build profile
        std::string species_type("");
        PyTools::extract("species_type",species_type,"Species",ispecies);

        PyObject *profile1;
        PyTools::extract_pyProfile("nb_density"    , profile1, "Species", ispecies);
        PyTools::extract_pyProfile("charge_density", profile1, "Species", ispecies);
        PyTools::extract_pyProfile("n_part_per_cell", profile1, "Species", ispecies);
        Profile *ppcProfile = new Profile(profile1, params.nDim_particle, "n_part_per_cell "+species_type);

        local_load = 0;
        // Count global number of particles, 
        for (unsigned int i=0; i<params.n_space_global[0]; i++) {
            for (unsigned int j=0; j<params.n_space_global[1]; j++) {
                for (unsigned int k=0; k<params.n_space_global[2]; k++) {
                    vector<double> x_cell(3,0);
                    x_cell[0] = (i+0.5)*params.cell_length[0];
                    x_cell[1] = (j+0.5)*params.cell_length[1];
                    x_cell[2] = (k+0.5)*params.cell_length[2];

                    int n_part_in_cell = round(ppcProfile->valueAt(x_cell));
                    if ( n_part_in_cell<=0. )
                        continue;
                    else
                        local_load += n_part_in_cell;
                }
            }
        }

        double time_frozen(0.);
        PyTools::extract("time_frozen",time_frozen ,"Species",ispecies);
        if(time_frozen > 0.) local_load *= params.coef_frozen;
        Tload += local_load;

        delete ppcProfile;

    } // End for ispecies

    Tload += Npatches*ncells_perpatch*params.coef_cell ; // We assume the load of one cell to be equal to coef_cell and account for ghost cells.
    if (isMaster()) {
        fout.open ("patch_load.txt");
        fout << "Total load = " << Tload << endl;
    }
    Tload /= Tcapabilities; //Target load for each mpi process.
    Tcur = Tload * capabilities[0];  //Init.
    //Tcur = 0;  //Init.

    //Loop over all patches
    for(unsigned int hindex=0; hindex < Npatches; hindex++){
        generalhilbertindexinv(params.mi[0], params.mi[1], &Pcoordinates[0], &Pcoordinates[1], hindex);
        for (unsigned int idim = 0; idim < params.nDim_field; idim++) {
            Pcoordinates[idim] *= params.n_space[idim]; //Compute patch cells coordinates
        }
        local_load = 0.; //Accumulate load of the current patch
        for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++){

            // Build profile
            std::string species_type("");
            PyTools::extract("species_type",species_type,"Species",ispecies);

            PyObject *profile1;
            PyTools::extract_pyProfile("nb_density"    , profile1, "Species", ispecies);
            PyTools::extract_pyProfile("charge_density", profile1, "Species", ispecies);
            PyTools::extract_pyProfile("n_part_per_cell", profile1, "Species", ispecies);
            Profile *ppcProfile = new Profile(profile1, params.nDim_particle, "n_part_per_cell "+species_type);

            vector<double> cell_index(3,0);
            for (unsigned int i=0 ; i<params.nDim_field ; i++) {
                if (params.cell_length[i]!=0)
                    cell_index[i] = Pcoordinates[i]*params.cell_length[i];
            }
            local_load_temp = 0; 
            // Count global number of particles, 
            for (unsigned int i=0; i<params.n_space[0]; i++) {
                for (unsigned int j=0; j<params.n_space[1]; j++) {
                    for (unsigned int k=0; k<params.n_space[2]; k++) {
                        vector<double> x_cell(3,0);
                        x_cell[0] = cell_index[0] + (i+0.5)*params.cell_length[0];
                        x_cell[1] = cell_index[1] + (j+0.5)*params.cell_length[1];
                        x_cell[2] = cell_index[2] + (k+0.5)*params.cell_length[2];

                        int n_part_in_cell = round(ppcProfile->valueAt(x_cell));
                        if ( n_part_in_cell<=0. )
                            continue;
                        else
                            local_load_temp += n_part_in_cell;
                    }
                }
            }
            delete ppcProfile;
            double time_frozen(0.);
            PyTools::extract("time_frozen",time_frozen ,"Species",ispecies);
            if(time_frozen > 0.) local_load_temp *= params.coef_frozen;

            local_load += local_load_temp; // Accumulate species contribution to the load.
        } // End for ispecies

        local_load += ncells_perpatch*params.coef_cell; //Add grid contribution to the load.
        Lcur += local_load; //Add grid contribution to the load.
        Ncur++; // Try to assign current patch to rank r.

        //if (isMaster()) cout <<"h= " << hindex << " Tcur = " << Tcur << " Lcur = " << Lcur <<" Ncur = " << Ncur <<" r= " << r << endl;
        if (r < smilei_sz-1){

            if ( Lcur > Tcur || smilei_sz-r >= Npatches-hindex){ //Load target is exceeded or we have as many patches as procs left.
                above_target = Lcur - Tcur;  //Including current patch, we exceed target by that much.
                below_target = Tcur - (Lcur-local_load); // Excluding current patch, we mis the target by that much.
                if(above_target > below_target) { // If we're closer to target without the current patch...
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
        if (hindex == Npatches-1){
            patch_count[smilei_sz-1] = Ncur; //When we reach the last patch, the last MPI process takes what's left.
        }
    }// End loop on patches.
    if (isMaster()) {
        for (unsigned int i=0; i<smilei_sz; i++) fout << "patch count = " << patch_count[i]<<endl;
        fout.close();
    }
    
} // END init_patch_count


// ---------------------------------------------------------------------------------------------------------------------
//  Recompute patch distribution
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::recompute_patch_count( Params& params, VectorPatch& vecpatches, double time_dual )
{

    unsigned int Npatches, r,Ncur,Pcoordinates[3],ncells_perpatch, Lmin, Lmintemp;
    double Tload,Tcur, Lcur, above_target, below_target;
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
 
    unsigned int tot_species_number = PyTools::nComponents("Species");

    Lp.resize(patch_count[smilei_rk],0.);
    Lp_global.resize(Npatches,0.);

    Tload = 0.;
    r = 0;  //Start by finding work for rank 0.
    Ncur = 0; // Number of patches assigned to current rank r.
    Lcur = 0.; //Load assigned to current rank r.

    //Compute Local Loads of each Patch (Lp)
    for(unsigned int hindex=0; hindex < patch_count[smilei_rk]; hindex++){
        for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++) {
            Lp[hindex] += vecpatches(hindex)->vecSpecies[ispecies]->getNbrOfParticles()*(1+(params.coef_frozen-1)*(time_dual > vecpatches(hindex)->vecSpecies[ispecies]->time_frozen)) ;
        }
        Lp[hindex] += ncells_perpatch*params.coef_cell ;
    }

    //Allgatherv loads of all patches
  
    recv_counts[0] = 0;
    for(unsigned int i=1; i < smilei_sz ; i++) recv_counts[i] = recv_counts[i-1]+patch_count[i-1];

    MPI_Allgatherv(&Lp[0],patch_count[smilei_rk],MPI_DOUBLE,&Lp_global[0], &patch_count[0], recv_counts, MPI_DOUBLE,MPI_COMM_WORLD);

    //Compute total loads
    for(unsigned int hindex=0; hindex < Npatches; hindex++) Tload += Lp_global[hindex];
    Tload /= Tcapabilities; //Target load for each mpi process.
    Tcur = Tload * capabilities[0];  //Init.

    //Loop over all patches
    for(unsigned int hindex=0; hindex < Npatches; hindex++){

        Lcur += Lp_global[hindex]; 
        Ncur++; // Try to assign current patch to rank r.

        if (r < smilei_sz-1){

            if ( Lcur > Tcur || smilei_sz-r >= Npatches-hindex){ //Load target is exceeded or we have as many patches as procs left.
                above_target = Lcur - Tcur;  //Including current patch, we exceed target by that much.
                below_target = Tcur - (Lcur-Lp_global[hindex]); // Excluding current patch, we mis the target by that much.
                if(above_target > below_target) { // If we're closer to target without the current patch...
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
        if (hindex == Npatches-1){
            target_patch_count[smilei_sz-1] = Ncur; //When we reach the last patch, the last MPI process takes what's left.
        }
    }// End loop on patches.


    //Make sure the new patch_count is not too different from the previous one.
    // First patch
    Lcur = patch_count[0]+ patch_count[1]-1;
    Ncur = 0;
    Tcur = target_patch_count[0];  
    Lmintemp = patch_count[0]+1;
    if (Tcur > Lcur ){
        patch_count[0] = Lcur - Ncur;
    } else {
        patch_count[0] = target_patch_count[0] - Ncur;
    }

    //Loop
    for(unsigned int i=1; i< smilei_sz-1; i++){
        Ncur += patch_count[i-1];                  //Ncur = sum n=0..i-1 new patch_count[n]
        Tcur += target_patch_count[i];            //Tcur = sum n=0..i target_patch_count[n]
        Lcur += patch_count[i+1];                  //Lcur = sum n=0..i+1 initial patch_count[n] -1
        Lmin = Lmintemp;
        Lmintemp += patch_count[i];
 
        if (Tcur < Lmin){                      
            patch_count[i] = Lmin - Ncur;
        } else if (Tcur > Lcur ){                      
            patch_count[i] = Lcur - Ncur;
        } else {
            patch_count[i] = Tcur-Ncur;
        }
    }

    //Last patch
    Ncur += patch_count[smilei_sz-2];                  //Ncur = sum n=0..i-1 new patch_count[n]
    patch_count[smilei_sz-1] = Npatches-Ncur;

    if (smilei_rk==0) {
        //fout << "\tt = " << scientific << setprecision(3) << time_dual << endl;
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

    unsigned int patch_counter,rank;
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
    for ( int iprop=0 ; iprop<particles->double_prop.size() ; iprop++ )
        MPI_Get_address( &( (*(particles->double_prop[iprop]))[0] ), &(address[iprop]) );
    for ( int iprop=0 ; iprop<particles->short_prop.size() ; iprop++ )
        MPI_Get_address( &( (*(particles->short_prop[iprop]))[0] ), &(address[particles->double_prop.size()+iprop]) );
    for ( int iprop=0 ; iprop<particles->uint_prop.size() ; iprop++ )
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
    for (int i=0 ; i<particles->double_prop.size() ; i++)
        partDataType[i] = MPI_DOUBLE;
    for ( int iprop=0 ; iprop<particles->short_prop.size() ; iprop++ )
        partDataType[ particles->double_prop.size()+iprop] = MPI_SHORT;
    for ( int iprop=0 ; iprop<particles->uint_prop.size() ; iprop++ )
        partDataType[ particles->double_prop.size()+particles->short_prop.size()+iprop] = MPI_UNSIGNED;

    MPI_Datatype typeParticlesMPI;
    MPI_Type_struct( nbrOfProp, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
    MPI_Type_commit( &typeParticlesMPI );
    
    return typeParticlesMPI;

} // END createMPIparticles


// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// -----------------------------------------       PATCH SEND / RECV METHODS        ------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI::isend(Patch* patch, int to, int tag)
{
    //MPI_Request request;

    for (int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++){
        isend( &(patch->vecSpecies[ispec]->bmax), to, tag+2*ispec+1 );
        if ( patch->vecSpecies[ispec]->getNbrOfParticles() > 0 ){
            patch->vecSpecies[ispec]->typePartSend = createMPIparticles( patch->vecSpecies[ispec]->particles );
            isend( patch->vecSpecies[ispec]->particles, to, tag+2*ispec, patch->vecSpecies[ispec]->typePartSend );
            MPI_Type_free( &(patch->vecSpecies[ispec]->typePartSend) );
        }
    }
    isend( patch->EMfields, to, tag+2*patch->vecSpecies.size() );

    for ( int idiag = 0 ; idiag < patch->localDiags.size() ; idiag++ ) {
	// just probes (track data = species, managed above, meta-data : H5S_select)
	if ( patch->localDiags[idiag]->type_ == "Probes" )
	    isend( static_cast<DiagnosticProbes*>(patch->localDiags[idiag]), to, tag+2*patch->vecSpecies.size()+9+idiag );
    }

} // END isend( Patch )


void SmileiMPI::recv(Patch* patch, int from, int tag, Params& params)
{
    int nbrOfPartsRecv;

    for (int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++){
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
   
    recv( patch->EMfields, from, tag+2*patch->vecSpecies.size() );

    for ( int idiag = 0 ; idiag < patch->localDiags.size() ; idiag++ ) {
	// just probes (track data = species, managed above, meta-data : H5S_select)
	if ( patch->localDiags[idiag]->type_ == "Probes" )
	    recv( static_cast<DiagnosticProbes*>(patch->localDiags[idiag]), from, tag+2*patch->vecSpecies.size()+9+idiag );
    }

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
    
    for (int antennaId=0 ; antennaId<EM->antennas.size() ; antennaId++)
        isend( EM->antennas[antennaId].field, to, tag+9+antennaId );
    
    tag += 10 + EM->antennas.size();
    
    for (int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if(! EM->emBoundCond[bcId]) continue;
        
        for (int laserId=0 ; laserId < EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
            
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
    
    for (int antennaId=0 ; antennaId<EM->antennas.size() ; antennaId++)
        recv( EM->antennas[antennaId].field, from, tag+9+antennaId );
    
    tag += 10 + EM->antennas.size();
    
    for (int bcId=0 ; bcId<EM->emBoundCond.size() ; bcId++ ) {
        if(! EM->emBoundCond[bcId]) continue;
        
        for (int laserId=0 ; laserId<EM->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
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


void SmileiMPI::isend( DiagnosticProbes* diags, int to, int tag )
{
    MPI_Request request; 
    MPI_Isend( &(diags->probesStart), 1, MPI_INT, to, tag, MPI_COMM_WORLD, &request );

} // End isend ( Diagnostics )


void SmileiMPI::recv( DiagnosticProbes* diags, int from, int tag )
{
    MPI_Status status;
    MPI_Recv( &(diags->probesStart), 1, MPI_INT, from, tag, MPI_COMM_WORLD, &status );


} // End recv ( Diagnostics )


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
    if ( diag->type_ == "Scalar" ) {
	DiagnosticScalar* scalar = static_cast<DiagnosticScalar*>( diag );
	computeGlobalDiags(scalar, timestep);
    }
    else if ( diag->type_ == "Particles" ) {
	DiagnosticParticles* particles = static_cast<DiagnosticParticles*>( diag );
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
    
    int nscalars(0);

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
        }
        else if ( (*iterKey).find("MaxCell") != std::string::npos ) {
            vector<double>::iterator iterVal = iter-1;
            val_index maxVal;
            maxVal.val   = (*iterVal);
            maxVal.index = (*iter);
            MPI_Reduce(isMaster()?MPI_IN_PLACE:&maxVal, &maxVal, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
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
        
        if (timestep==0) {
            scalars->Energy_time_zero  = Utot;
            scalars->EnergyUsedForNorm = scalars->Energy_time_zero;
        }
        
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
