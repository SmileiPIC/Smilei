// CLOBAL COORDINATES: 
//                             MPI_minGlobal                                                                        MPI_maxGlobal
//                      --------<===================================== gs ===================================>------------
//     GLOBAL INDICES:          0                                  .                                        nspace_global
//                           ix+oversize                                                                  ix+oversize
//                      ------------------------------------       .              ------------------------------------
//                      |   |   |     ...          |   |   |       .              |   |   |   |   ...    |   |   |   |
//                      |   |   |     ...          |   |   |       .              |   |   |   |   ...    |   |   |   |
//                      ------------------------------------       .              ------------------------------------
//                            MPI_minLocal      MPI_maxLocal       .               MPI_minLocal          MPI_maxLocal
//                                                 ----------------------------------------                 
//                                                 |   |   |       .              |   |   |
//                                                 |   |   |       .              |   |   |
//                                                 ----------------------------------------
// LOCAL COORDINATES:                             x(0) rlb        x(ix)             rub  x(nspace)
//                                                 ----<============= length =========>----
//     LOCAL INDICES:                              0   lb                            ub   nspace

#include "SmileiMPI.h"

#include <cmath>
#include <cstring>

#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "DiagParams.h"
#include "Diagnostic.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "Field.h"

#include "Species.h"
#include "Hilbert_functions.h"
#include "Patch.h"

using namespace std;

SmileiMPI::SmileiMPI( int* argc, char*** argv )
{    
    int mpi_provided;

    MPI_Init_thread( argc, argv, MPI_THREAD_FUNNELED, &mpi_provided );
    if (mpi_provided == MPI_THREAD_SINGLE){
        MESSAGE("openMP not supported");
    }

    SMILEI_COMM_WORLD = MPI_COMM_WORLD;
    MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
    MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );

}

SmileiMPI::SmileiMPI( SmileiMPI *smpi )
{
    SMILEI_COMM_WORLD = smpi->SMILEI_COMM_WORLD;
    MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
    MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );

    oversize = smpi->oversize;
    cell_starting_global_index = smpi->cell_starting_global_index;
    min_local = smpi->min_local;
    max_local = smpi->max_local;

    n_space_global = smpi->n_space_global;
    patch_count.resize(smilei_sz, 0);

}

SmileiMPI::~SmileiMPI()
{
    delete[]periods_;

    int status = 0;
    MPI_Finalized( &status );
    if (!status) MPI_Finalize();

}

void SmileiMPI::bcast( InputData& input_data )
{
    DEBUG(10,"broadcast namelist");
    bcast(input_data.namelist);    

    input_data.parseStream();
    
    // Randomization
    unsigned long seedTime=0;
    input_data.extract("random_seed",seedTime);
    srand(seedTime+getRank());
    
}

void SmileiMPI::bcast( string& val )
{
    int charSize=0;
    if (isMaster()) charSize = val.size()+1;
    MPI_Bcast(&charSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);

    char tmp[charSize];
    strcpy(tmp, val.c_str());
    MPI_Bcast(&tmp, charSize, MPI_CHAR, 0, SMILEI_COMM_WORLD);

    if (!isMaster()) val=tmp;

}

void SmileiMPI::init( PicParams& params )
{


    oversize.resize(params.nDim_field, 0);
    cell_starting_global_index.resize(params.nDim_field, 0);
    min_local.resize(params.nDim_field, 0.);
    max_local.resize(params.nDim_field, 0.);
    n_space_global.resize(params.nDim_field, 0);
    patch_count.resize(smilei_sz, 0);
    target_patch_count.resize(smilei_sz, 0);

    interParticles.initialize(0,params.nDim_particle); 
 
    init_patch_count(params);

    periods_  = new int[params.nDim_field];
    for (int i=0 ; i<params.nDim_field ; i++) periods_[i] = 0;
    // Geometry periodic in x
    if (params.bc_em_type_long=="periodic") {
        periods_[0] = 1;
        MESSAGE(1,"applied topology for periodic BCs in x-direction");
    }
    if (params.nDim_field>1) {
	// Geometry periodic in y
	if (params.bc_em_type_trans=="periodic") {
	    periods_[1] = 1;
	    MESSAGE(2,"applied topology for periodic BCs in y-direction");
	}
    }
}

void SmileiMPI::init_patch_count( PicParams& params)
{
    unsigned int Npatches, r,Ncur,Pcoordinates[3],ncells_perpatch, Tcapabilities;
    double Tload,Tcur, Lcur, local_load, local_load_temp, above_target, below_target;
    std::vector<unsigned int> mincell,maxcell,capabilities; //Min and max values of non empty cells for each species and in each dimension.
    double density_length[3];
    //Load of a cell = coef_cell*load of a particle.
    //Load of a frozen particle = coef_frozen*load of a particle.
    double coef_cell, coef_frozen; 

    coef_cell = 0.1;
    coef_frozen = 0.1;
 
    mincell.resize(params.n_species*3);
    maxcell.resize(params.n_species*3);
    capabilities.resize(smilei_sz, 1); //Capabilities of devices hosting the different mpi processes. All capabilities are assumed to be equal for the moment.
    Tcapabilities = 0;
    for (unsigned int i = 0; i < smilei_sz; i++)
        Tcapabilities += capabilities[i];

    //Compute total Load
    
    Tload = 0.;
    Npatches = params.number_of_patches[0];
    for (unsigned int i = 1; i < params.nDim_field; i++)
        Npatches *=  params.number_of_patches[i]; // Total number of patches.

    ncells_perpatch = params.n_space[0]+2*params.oversize[0]; //Initialization
    for (unsigned int idim = 1; idim < params.nDim_field; idim++)
        ncells_perpatch *= params.n_space[idim]+2*params.oversize[idim];

    r = 0;  //Start by finding work for rank 0.
    Ncur = 0; // Number of patches assigned to current rank r.
    Lcur = 0.; //Load assigned to current rank r.


    for (unsigned int ispecies = 0; ispecies < params.n_species; ispecies++){
        density_length[0] = params.species_param[ispecies].dens_length_x[0];
        density_length[1] = params.species_param[ispecies].dens_length_y[0];
        //Needs to be updated when dens_lenth is a vector in params.
        //density_length[2] = params.species_param[ispecies].dens_length_z[0];

        local_load = params.species_param[ispecies].n_part_per_cell ;
        if(params.species_param[ispecies].time_frozen > 0.) local_load *= coef_frozen; // Assumes time_dual = 0. Might be false at restart !
        for (unsigned int idim = 0; idim < params.nDim_field; idim++){
            mincell[ispecies*3+idim] = params.species_param[ispecies].vacuum_length[idim]/params.cell_length[idim];
            // This strange way of writing this below prevents rounding errors.
            maxcell[ispecies*3+idim] = min (params.sim_length[idim]/params.cell_length[idim],  mincell[ispecies*3+idim] + density_length[idim]/params.cell_length[idim]);
            local_load *= (maxcell[ispecies*3+idim]-mincell[ispecies*3+idim]);
        }
        Tload += local_load; //Particle contribution to the load
    }
    Tload += Npatches*ncells_perpatch*coef_cell ; // We assume the load of one cell to be equal to be coef_cell and account for ghost cells.
    if (isMaster()) cout << "Total load = " << Tload << endl;
    Tload /= Tcapabilities; //Target load for each mpi process.
    Tcur = Tload * capabilities[0];  //Init.

    //Loop over all patches
    for(unsigned int hindex=0; hindex < Npatches; hindex++){
        generalhilbertindexinv(params.mi[0], params.mi[1], &Pcoordinates[0], &Pcoordinates[1], hindex);
        for (unsigned int idim = 0; idim < params.nDim_field; idim++) {
            Pcoordinates[idim] *= params.n_space[idim]; //Compute patch cells coordinates
        }
        local_load = 0.; //Accumulate load of the current patch
        for (unsigned int ispecies = 0; ispecies < params.n_species; ispecies++){
            local_load_temp = params.species_param[ispecies].n_part_per_cell; //Accumulate load of the current species.
            if(params.species_param[ispecies].time_frozen > 0.) local_load_temp *= coef_frozen;
            for (unsigned int idim = 0; idim < params.nDim_field; idim++){
                local_load_temp *= min (min( (int)(maxcell[ispecies*3+idim]-Pcoordinates[idim]), (int)params.n_space[idim]), min((int)(Pcoordinates[idim]+params.n_space[idim]-mincell[ispecies*3+idim]), (int)params.n_space[idim]));
                if (local_load_temp < 0.) local_load_temp = 0.;
            } 
            local_load += local_load_temp; // Accumulate species contribution to the load.
        }

        local_load += ncells_perpatch*coef_cell; //Add grid contribution to the load.
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
                    Lcur = local_load;
                } else {                          //Else ...
                    patch_count[r] = Ncur;        //...assign patches including the current one.
                    Ncur = 0;
                    Lcur = 0.;
                }
                r++; //Move on to the next rank.
                Tcur = Tload * capabilities[r];  //Target load for current rank r.
            } 
        }// End if on r.
        if (hindex == Npatches-1){
            patch_count[smilei_sz-1] = Ncur; //When we reach the last patch, the last MPI process takes what's left.
        }
    }// End loop on patches.
    if (isMaster()) for (unsigned int i=0; i<smilei_sz; i++) cout << "patch count = " << patch_count[i]<<endl;

}

void SmileiMPI::recompute_patch_count( PicParams& params, VectorPatch& vecpatches, double time_dual )
{
    unsigned int Npatches, r,Ncur,Pcoordinates[3],ncells_perpatch, Tcapabilities, Lmin, Lmintemp;
    double Tload,Tcur, Lcur, above_target, below_target;
    std::vector<unsigned int> mincell,maxcell,capabilities; //Min and max values of non empty cells for each species and in each dimension.
    //Load of a cell = coef_cell*load of a particle.
    //Load of a frozen particle = coef_frozen*load of a particle.
    double coef_cell, coef_frozen; 
    std::vector<double> Lp,Lp_global;
    int recv_counts[smilei_sz];

    coef_cell = 0.1;
    coef_frozen = 0.1;

    Npatches = params.number_of_patches[0];
    for (unsigned int i = 1; i < params.nDim_field; i++) 
        Npatches *=  params.number_of_patches[i]; // Total number of patches.

    ncells_perpatch = params.n_space[0]+2*params.oversize[0]; //Initialization
    for (unsigned int idim = 1; idim < params.nDim_field; idim++)
        ncells_perpatch *= params.n_space[idim]+2*params.oversize[idim];
 
    mincell.resize(params.n_species*3);
    maxcell.resize(params.n_species*3);
    capabilities.resize(smilei_sz, 1); //Capabilities of devices hosting the different mpi processes. All capabilities are assumed to be equal for the moment.
    Tcapabilities = 0;
    for (unsigned int i = 0; i < smilei_sz; i++)
        Tcapabilities += capabilities[i];
    Lp.resize(patch_count[smilei_rk],0.);
    Lp_global.resize(Npatches,0.);

    Tload = 0.;
    r = 0;  //Start by finding work for rank 0.
    Ncur = 0; // Number of patches assigned to current rank r.
    Lcur = 0.; //Load assigned to current rank r.

    //Compute Local Loads of each Patch (Lp)
    for(unsigned int hindex=0; hindex < patch_count[smilei_rk]; hindex++){
        for (unsigned int ispecies = 0; ispecies < params.n_species; ispecies++)
            Lp[hindex] += params.species_param[ispecies].n_part_per_cell*vecpatches(hindex)->vecSpecies[ispecies]->getNbrOfParticles()*(1+(coef_frozen-1)*(time_dual > params.species_param[ispecies].time_frozen)) ;
        Lp[hindex] += ncells_perpatch*coef_cell ;
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
                    Lcur = Lp_global[hindex];
                } else {                          //Else ...
                    target_patch_count[r] = Ncur;        //...assign patches including the current one.
                    Ncur = 0;
                    Lcur = 0;
                }
                r++; //Move on to the next rank.
                Tcur = Tload * capabilities[r];  //Target load for current rank r.
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

    return;
}



void SmileiMPI::sumRho( ElectroMagn* EMfields )
{
    sumField( EMfields->rho_ );

}

void SmileiMPI::sumRhoJ( ElectroMagn* EMfields )
{
    // sum total charge density and currents
    //sumField( EMfields->rho_ );
    sumField( EMfields->Jx_ );
    sumField( EMfields->Jy_ );
    sumField( EMfields->Jz_ );
   
}
void SmileiMPI::sumRhoJs( ElectroMagn* EMfields, int ispec, bool currents )
{
   // sum density and currents for all species
   sumField( EMfields->rho_s[ispec] );
   if(currents){
       sumField( EMfields->Jx_s[ispec] );
       sumField( EMfields->Jy_s[ispec] );
       sumField( EMfields->Jz_s[ispec] );
   }
}

void SmileiMPI::exchangeE( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Ex_ );
    exchangeField( EMfields->Ey_ );
    exchangeField( EMfields->Ez_ );

}
void SmileiMPI::exchangeE( ElectroMagn* EMfields, unsigned int clrw )
{
    exchangeField_movewin( EMfields->Ex_, clrw );
    exchangeField_movewin( EMfields->Ey_, clrw );
    exchangeField_movewin( EMfields->Ez_, clrw );

}

void SmileiMPI::exchangeB( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Bx_ );
    exchangeField( EMfields->By_ );
    exchangeField( EMfields->Bz_ );

}
void SmileiMPI::exchangeB( ElectroMagn* EMfields, unsigned int clrw )
{
    exchangeField_movewin( EMfields->Bx_, clrw );
    exchangeField_movewin( EMfields->By_, clrw);
    exchangeField_movewin( EMfields->Bz_, clrw );

}

void SmileiMPI::exchangeBm( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Bx_m );
    exchangeField( EMfields->By_m );
    exchangeField( EMfields->Bz_m );

}
void SmileiMPI::exchangeBm( ElectroMagn* EMfields, unsigned int clrw )
{
    exchangeField_movewin( EMfields->Bx_m, clrw );
    exchangeField_movewin( EMfields->By_m, clrw );
    exchangeField_movewin( EMfields->Bz_m, clrw );

}

void SmileiMPI::exchangeAvg( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Ex_avg );
    exchangeField( EMfields->Ey_avg );
    exchangeField( EMfields->Ez_avg );
    exchangeField( EMfields->Bx_avg );
    exchangeField( EMfields->By_avg );
    exchangeField( EMfields->Bz_avg );
}

// ! Attention 1D
MPI_Datatype SmileiMPI::createMPIparticles( Particles* particles, int nbrOfProp )
{
    if (particles->Position.size()!=2)
	cout << "1D Particles not managed" << endl;

    MPI_Datatype typeParticlesMPI;

    MPI_Aint address[nbrOfProp];
    MPI_Get_address( &(particles->position(0,0)), &(address[0]) );
    MPI_Get_address( &(particles->position(1,0)), &(address[1]) );
    //MPI_Get_address( &(particles->position_old(0,0)), &(address[2]) );
    //MPI_Get_address( &(particles->position_old(1,0)), &(address[3]) );
    MPI_Get_address( &(particles->momentum(0,0)), &(address[2]) );
    MPI_Get_address( &(particles->momentum(1,0)), &(address[3]) );
    MPI_Get_address( &(particles->momentum(2,0)), &(address[4]) );
    MPI_Get_address( &(particles->weight(0)),     &(address[5]) );
    MPI_Get_address( &(particles->charge(0)),     &(address[6]) );

    int nbr_parts[nbrOfProp];
    MPI_Aint disp[nbrOfProp];
    MPI_Datatype partDataType[nbrOfProp];

    for (int i=0 ; i<nbrOfProp ; i++)
	nbr_parts[i] = particles->size();
    disp[0] = 0;
    for (int i=1 ; i<nbrOfProp ; i++)
	disp[i] = address[i] - address[0];
    for (int i=0 ; i<nbrOfProp ; i++)
	partDataType[i] = MPI_DOUBLE;
    partDataType[nbrOfProp-1] = MPI_SHORT;

    MPI_Type_struct( nbrOfProp, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
    MPI_Type_commit( &typeParticlesMPI );

    return typeParticlesMPI;
} // END createMPIparticles


// Returns the rank of the MPI process currently owning patch h.
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
}

void SmileiMPI::computeGlobalDiags(Diagnostic* diags, int timestep)
{
    if (timestep % diags->scalars.every == 0) computeGlobalDiags(diags->scalars, timestep);
    //computeGlobalDiags(probes); // HDF5 write done per patch in DiagProbes::*
    //computeGlobalDiags(phases);
}

void SmileiMPI::computeGlobalDiags(DiagnosticScalar& scalars, int timestep)
{
    int nscalars(0);
    for(vector<pair<string,double> >::iterator iter = scalars.out_list.begin(); iter !=scalars.out_list.end(); iter++) {
	if ( ( (iter->first).find("Min") == std::string::npos ) && ( (iter->first).find("Max") == std::string::npos ) ) {
	    MPI_Reduce(isMaster()?MPI_IN_PLACE:&((*iter).second), &((*iter).second), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else if ( (iter->first).find("MinCell") != std::string::npos ) {
	    vector<pair<string,double> >::iterator iterVal = iter-1;
	    val_index minVal;
	    minVal.val   = (*iterVal).second;
	    minVal.index = (*iter).second;
	    MPI_Reduce(isMaster()?MPI_IN_PLACE:&minVal, &minVal, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
	}
	else if ( (iter->first).find("MaxCell") != std::string::npos ) {
	    vector<pair<string,double> >::iterator iterVal = iter-1;
	    val_index maxVal;
	    maxVal.val   = (*iterVal).second;
	    maxVal.index = (*iter).second;
	    MPI_Reduce(isMaster()?MPI_IN_PLACE:&maxVal, &maxVal, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
	}
    }

    if (isMaster()) {

	double Etot_part = scalars.getScalar("Eparticles");
	double Etot_fields = scalars.getScalar("EFields");
	
	double Energy_time_zero = 0.;//scalars.Energy_time_zero;
	double poyTot = scalars.getScalar("Poynting");

	double Elost_part = scalars.getScalar("Elost");
	double Emw_lost  = scalars.getScalar("Emw_lost");
	double Emw_lost_fields = scalars.getScalar("Emw_lost_fields");

	double Total_Energy=Etot_part+Etot_fields;

	double Energy_Balance=Total_Energy-(Energy_time_zero+poyTot)+Elost_part+Emw_lost+Emw_lost_fields;
	//double Energy_Balance=Total_Energy-(Energy_time_zero);
	double Energy_Bal_norm(0.);
	if (scalars.EnergyUsedForNorm>0.)
	    Energy_Bal_norm=Energy_Balance/scalars.EnergyUsedForNorm;
	scalars.EnergyUsedForNorm = Total_Energy;

	scalars.setScalar("Etot",Total_Energy);
	scalars.setScalar("Ebalance",Energy_Balance);
	scalars.setScalar("Ebal_norm",Energy_Bal_norm);

	if (timestep==0) {
	    scalars.Energy_time_zero  = Total_Energy;
	    scalars.EnergyUsedForNorm = Energy_time_zero;
	}


    }
	

    scalars.write(timestep);

}

void SmileiMPI::send(Patch* patch, int to, int hindex)
{
    for (int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++)
	send( patch->vecSpecies[ispec], to, hindex );

    send( patch->EMfields, to, hindex );

}

void SmileiMPI::recv(Patch* patch, int from, int hindex)
{
    for (int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++)
	recv( patch->vecSpecies[ispec], from, hindex );

    recv( patch->EMfields, from, hindex );


}

void SmileiMPI::send(Species* species, int to, int hindex)
{
    send( species->particles, to, hindex );
    send( species->bmin, to, hindex );
    send( species->bmax, to, hindex );
}

void SmileiMPI::recv(Species* species, int from, int hindex)
{
    recv( species->particles, from, hindex );
    recv( &species->bmin, from, hindex );
    recv( &species->bmax, from, hindex );
}

void SmileiMPI::send(Particles* particles, int to, int hindex)
{
    // Number of properties per particles = nDim_Particles + 3 + 1 + 1
    int nbrOfProp( 7 );
    MPI_Datatype typePartSend = createMPIparticles( particles, nbrOfProp );
    MPI_Send( &(particles->position(0,0)), 1, typePartSend, to, hindex, MPI_COMM_WORLD );
    MPI_Type_free( &typePartSend );

}

void SmileiMPI::recv(Particles* particles, int from, int hindex)
{
    MPI_Status status;

    // Number of properties per particles = nDim_Particles + 3 + 1 + 1
    int nbrOfProp( 7 );
    MPI_Datatype typePartRecv = createMPIparticles( particles, nbrOfProp );
    MPI_Recv( &(particles->position(0,0)), 1, typePartRecv, from, hindex, MPI_COMM_WORLD, &status );
    MPI_Type_free( &typePartRecv );

}

// Assuming vec.size() is known (number of species)
void SmileiMPI::send(std::vector<int> vec, int to, int hindex)
{
    MPI_Send( &(vec[0]), vec.size(), MPI_INT, to, hindex, MPI_COMM_WORLD );
    barrier();

}

void SmileiMPI::recv(std::vector<int> *vec, int from, int hindex)
{
    MPI_Status status;
    MPI_Recv( &((*vec)[0]), vec->size(), MPI_INT, from, hindex, MPI_COMM_WORLD, &status );
    barrier();

}

void SmileiMPI::send(ElectroMagn* fields, int to, int hindex)
{
    send( fields->Ex_, to, hindex );
    send( fields->Ey_, to, hindex );
    send( fields->Ez_, to, hindex );
    send( fields->Bx_, to, hindex );
    send( fields->By_, to, hindex );
    send( fields->Bz_, to, hindex );
}

void SmileiMPI::recv(ElectroMagn* fields, int from, int hindex)
{
    recv( fields->Ex_, from, hindex );
    recv( fields->Ey_, from, hindex );
    recv( fields->Ez_, from, hindex );
    recv( fields->Bx_, from, hindex );
    recv( fields->By_, from, hindex );
    recv( fields->Bz_, from, hindex );
}

void SmileiMPI::send(Field* field, int to, int hindex)
{
    MPI_Send( &((*field)(0)),field->globalDims_, MPI_DOUBLE, to, hindex, MPI_COMM_WORLD );
}

void SmileiMPI::recv(Field* field, int from, int hindex)
{
    MPI_Status status;
    MPI_Recv( &((*field)(0)),field->globalDims_, MPI_DOUBLE, from, hindex, MPI_COMM_WORLD, &status );

}
