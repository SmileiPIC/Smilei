
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

#include "Params.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "Field.h"

#include "Species.h"
#include "Hilbert_functions.h"
#include "VectorPatch.h"

using namespace std;

SmileiMPI::SmileiMPI( int* argc, char*** argv )
{    
    int mpi_provided;

    MPI_Init_thread( argc, argv, MPI_THREAD_MULTIPLE, &mpi_provided );
    if (mpi_provided != MPI_THREAD_MULTIPLE){
        MESSAGE("MPI_THREAD_MULTIPLE not supported");
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

void SmileiMPI::bcast( string& val )
{
    int charSize=0;
    if (isMaster()) charSize = val.size()+1;
    MPI_Bcast(&charSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);

    char tmp[charSize];
    if (isMaster()) strcpy(tmp, val.c_str());
    MPI_Bcast(&tmp, charSize, MPI_CHAR, 0, SMILEI_COMM_WORLD);

    if (!isMaster()) val=tmp;

}

void SmileiMPI::bcast( int& val )
{
    MPI_Bcast(&val, 1, MPI_INT, 0, SMILEI_COMM_WORLD);
}

void SmileiMPI::init( Params& params )
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
}

void SmileiMPI::init_patch_count( Params& params)
{
//#ifndef _NOTBALANCED
//    bool use_load_balancing(true);
//    if (!use_load_balancing) {
//	int Npatches = params.number_of_patches[0];
//	for (unsigned int i = 1; i < params.nDim_field; i++)
//	    Npatches *=  params.number_of_patches[i]; // Total number of patches.
//	if (Npatches!=smilei_sz) ERROR("number of patches abd MPI processes");
//	for (unsigned int i=0; i<smilei_sz; i++) patch_count[i] = 1;
//	return;
//    }
//#endif

    unsigned int Npatches, r,Ncur,Pcoordinates[3],ncells_perpatch, Tcapabilities;
    double Tload,Tcur, Lcur, local_load, local_load_temp, above_target, below_target;
    std::vector<unsigned int> mincell,maxcell,capabilities; //Min and max values of non empty cells for each species and in each dimension.
    vector<double> density_length(params.nDim_field,0.);
    //Load of a cell = coef_cell*load of a particle.
    //Load of a frozen particle = coef_frozen*load of a particle.
    double coef_cell, coef_frozen; 

    coef_cell = 50;
    coef_frozen = 0.1;
 
    unsigned int tot_species_number = PyTools::nComponents("Species");
    mincell.resize(tot_species_number*3);
    maxcell.resize(tot_species_number*3);
                     
    capabilities.resize(smilei_sz, 1); //Capabilities of devices hosting the different mpi processes. All capabilities are assumed to be equal for the moment.
    //Compute total capability: Tcapabilities
    Tcapabilities = 0;
    for (unsigned int i = 0; i < smilei_sz; i++)
        Tcapabilities += capabilities[i];

    //Compute target load: Tload = Total load * local capability / Total capability.
    
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


    for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++){

        //Needs to be updated when dens_lenth is a vector in params.

#ifdef _BEFORE_PY
        density_length[0] = params.species_param[ispecies].dens_length_x[0];
	if (params.nDim_field>1) {
	    density_length[1] = params.species_param[ispecies].dens_length_y[0];
	    if (params.nDim_field>2)
		density_length[2] = params.species_param[ispecies].dens_length_z[0];
	}
	

        local_load = params.species_param[ispecies].n_part_per_cell ;
        if(params.species_param[ispecies].time_frozen > 0.) local_load *= coef_frozen; // Assumes time_dual = 0. Might be false at restart !
        for (unsigned int idim = 0; idim < params.nDim_field; idim++){
            mincell[ispecies*3+idim] = params.species_param[ispecies].vacuum_length[idim]/params.cell_length[idim];
            // This strange way of writing this below prevents rounding errors.
            maxcell[ispecies*3+idim] = min (params.sim_length[idim]/params.cell_length[idim],  mincell[ispecies*3+idim] + density_length[idim]/params.cell_length[idim]);
            local_load *= (maxcell[ispecies*3+idim]-mincell[ispecies*3+idim]);
        }
        Tload += local_load; //Particle contribution to the load
#else
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
	if(time_frozen > 0.) local_load *= coef_frozen;
	Tload += local_load;

	delete ppcProfile;
#endif
    } // End for ispecies

    Tload += Npatches*ncells_perpatch*coef_cell ; // We assume the load of one cell to be equal to coef_cell and account for ghost cells.
    if (isMaster()) cout << "Total load = " << Tload << endl;
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

#ifdef _BEFORE_PY
            local_load_temp = params.species_param[ispecies].n_part_per_cell; //Accumulate load of the current species.
            if(params.species_param[ispecies].time_frozen > 0.) local_load_temp *= coef_frozen;
            for (unsigned int idim = 0; idim < params.nDim_field; idim++){
                local_load_temp *= min (min( (int)(maxcell[ispecies*3+idim]-Pcoordinates[idim]), (int)params.n_space[idim]), min((int)(Pcoordinates[idim]+params.n_space[idim]-mincell[ispecies*3+idim]), (int)params.n_space[idim]));
                if (local_load_temp < 0.) local_load_temp = 0.;
            } 
#else
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
	    if(time_frozen > 0.) local_load_temp *= coef_frozen;

#endif
            local_load += local_load_temp; // Accumulate species contribution to the load.
        } // End for ispecies

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
    if (isMaster()) for (unsigned int i=0; i<smilei_sz; i++) cout << "patch count = " << patch_count[i]<<endl;

}

void SmileiMPI::recompute_patch_count( Params& params, VectorPatch& vecpatches, double time_dual )
{

#ifdef _PY_PATCH_COUNT
    MESSAGE( "Plasma not defined, may used python" );
#else

    unsigned int Npatches, r,Ncur,Pcoordinates[3],ncells_perpatch, Tcapabilities, Lmin, Lmintemp;
    double Tload,Tcur, Lcur, above_target, below_target;
    std::vector<unsigned int> mincell,maxcell,capabilities; //Min and max values of non empty cells for each species and in each dimension.
    //Load of a cell = coef_cell*load of a particle.
    //Load of a frozen particle = coef_frozen*load of a particle.
    double coef_cell, coef_frozen; 
    std::vector<double> Lp,Lp_global;
    int recv_counts[smilei_sz];

    coef_cell = 50;
    coef_frozen = 0.1;

    Npatches = params.number_of_patches[0];
    for (unsigned int i = 1; i < params.nDim_field; i++) 
        Npatches *=  params.number_of_patches[i]; // Total number of patches.

    ncells_perpatch = params.n_space[0]+2*params.oversize[0]; //Initialization
    for (unsigned int idim = 1; idim < params.nDim_field; idim++)
        ncells_perpatch *= params.n_space[idim]+2*params.oversize[idim];
 
    unsigned int tot_species_number = PyTools::nComponents("Species");

    mincell.resize(tot_species_number*3);
    maxcell.resize(tot_species_number*3);
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
        for (unsigned int ispecies = 0; ispecies < tot_species_number; ispecies++) {
	    //Lp[hindex] += params.species_param[ispecies].n_part_per_cell*vecpatches(hindex)->vecSpecies[ispecies]->getNbrOfParticles()*(1+(coef_frozen-1)*(time_dual > params.species_param[ispecies].time_frozen)) ;
	    //Lp[hindex] += params.species_param[ispecies].ppc_profile*vecpatches(hindex)->vecSpecies[ispecies]->getNbrOfParticles()*(1+(coef_frozen-1)*(time_dual > vecpatches(hindex)->vecSpecies[ispecies]->time_frozen)) ;
	    Lp[hindex] += vecpatches(hindex)->vecSpecies[ispecies]->getNbrOfParticles()*(1+(coef_frozen-1)*(time_dual > vecpatches(hindex)->vecSpecies[ispecies]->time_frozen)) ;
	}
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
    //Tcur = 0;  //Init.


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
                    //Lcur = Lp_global[hindex];
                } else {                          //Else ...
                    target_patch_count[r] = Ncur;        //...assign patches including the current one.
                    Ncur = 0;
                    //Lcur = 0;
                }
                r++; //Move on to the next rank.
                Tcur += Tload * capabilities[r];  //Target load for current rank r.
                //Tcur = Tload * capabilities[r];  //Target load for current rank r.
            } 
        }// End if on r.
        if (hindex == Npatches-1){
            target_patch_count[smilei_sz-1] = Ncur; //When we reach the last patch, the last MPI process takes what's left.
        }
    }// End loop on patches.


    // ------------- target patch count
//#ifdef _TESTPATCHEXCH
//    for (int irk=0;irk<smilei_sz;irk++) {
//	patch_count[irk] = 2;
//	//if (irk==4) patch_count[irk] = 29;
//    }
//#endif
    // ------------- target patch count

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

	if (smilei_rk==0)
	    for (int irk=0;irk<smilei_sz;irk++)
		cout << " patch_count[" << irk << "] = " << patch_count[irk] << " target patch_count = "<< target_patch_count[irk] << endl;

    return;
#endif
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

MPI_Datatype SmileiMPI::createMPIparticles( Particles* particles, int nbrOfProp )
{
    int nbrOfProp2 = particles->double_prop.size() + particles->short_prop.size() + particles->uint_prop.size();

    MPI_Aint address[nbrOfProp2];
    for ( int iprop=0 ; iprop<particles->double_prop.size() ; iprop++ )
	MPI_Get_address( &( (*(particles->double_prop[iprop]))[0] ), &(address[iprop]) );
    for ( int iprop=0 ; iprop<particles->short_prop.size() ; iprop++ )
        MPI_Get_address( &( (*(particles->short_prop[iprop]))[0] ), &(address[particles->double_prop.size()+iprop]) );
    for ( int iprop=0 ; iprop<particles->uint_prop.size() ; iprop++ )
        MPI_Get_address( &( (*(particles->uint_prop[iprop]))[0] ), &(address[particles->double_prop.size()+particles->short_prop.size()+iprop]) );

    int nbr_parts[nbrOfProp2];
    // number of elements per property
    for (int i=0 ; i<nbrOfProp2 ; i++)
        nbr_parts[i] = particles->size();

    MPI_Aint disp[nbrOfProp2];
    // displacement between 2 properties
    disp[0] = 0;
    for (int i=1 ; i<nbrOfProp2 ; i++)
        disp[i] = address[i] - address[0];

    MPI_Datatype partDataType[nbrOfProp2];
    // define MPI type of each property, default is DOUBLE
    for (int i=0 ; i<particles->double_prop.size() ; i++)
        partDataType[i] = MPI_DOUBLE;
    for ( int iprop=0 ; iprop<particles->short_prop.size() ; iprop++ )
        partDataType[ particles->double_prop.size()+iprop] = MPI_SHORT;
    for ( int iprop=0 ; iprop<particles->uint_prop.size() ; iprop++ )
        partDataType[ particles->double_prop.size()+particles->short_prop.size()+iprop] = MPI_UNSIGNED;

    MPI_Datatype typeParticlesMPI;
    MPI_Type_struct( nbrOfProp2, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
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
    computeGlobalDiags(diags->phases, timestep);
    for (unsigned int i=0; i<diags->vecDiagnosticParticles.size(); i++)
	computeGlobalDiags(diags->vecDiagnosticParticles[i], timestep);
}

void SmileiMPI::computeGlobalDiags(DiagnosticScalar& scalars, int timestep)
{
    int nscalars(0);

    vector<string>::iterator iterKey = scalars.out_key.begin();
    for(vector<double>::iterator iter = scalars.out_value.begin(); iter !=scalars.out_value.end(); iter++) {
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

	double Ukin = scalars.getScalar("Ukin");
	double Uelm = scalars.getScalar("Uelm");

	// added & lost energies due to the moving window
	double Ukin_out_mvw = scalars.getScalar("Ukin_out_mvw");
	double Ukin_inj_mvw = scalars.getScalar("Ukin_inj_mvw");
	double Uelm_out_mvw = scalars.getScalar("Uelm_out_mvw");
	double Uelm_inj_mvw = scalars.getScalar("Uelm_inj_mvw");
        
	// added & lost energies at the boundaries
	double Ukin_bnd = scalars.getScalar("Ukin_bnd");
	double Uelm_bnd = scalars.getScalar("Uelm_bnd");

	// total energy in the simulation
	double Utot = Ukin + Uelm;

	// expected total energy
	double Uexp = scalars.Energy_time_zero + Uelm_bnd + Ukin_inj_mvw + Uelm_inj_mvw
	    -           ( Ukin_bnd + Ukin_out_mvw + Uelm_out_mvw );
        
	// energy balance
	double Ubal = Utot - Uexp;
        
	// energy used for normalization
	scalars.EnergyUsedForNorm = Utot;
        
	// normalized energy balance
	double Ubal_norm(0.);
	if (scalars.EnergyUsedForNorm>0.)
	    Ubal_norm = Ubal / scalars.EnergyUsedForNorm;

	scalars.setScalar("Ubal_norm",Ubal_norm);
	scalars.setScalar("Ubal",Ubal);
	scalars.setScalar("Uexp",Uexp);
	scalars.setScalar("Utot",Utot);

	if (timestep==0) {
	    scalars.Energy_time_zero  = Utot;
	    scalars.EnergyUsedForNorm = scalars.Energy_time_zero;
	}


    }
	

    scalars.write(timestep);

}

void SmileiMPI::computeGlobalDiags(DiagnosticPhaseSpace& phases, int timestep)
{
    int nDiags( phases.vecDiagPhaseToRun.size() );
    for (vector<DiagnosticPhase*>::const_iterator diag=phases.vecDiagPhaseToRun.begin() ; diag != phases.vecDiagPhaseToRun.end(); diag++) 
	(*diag)->writeData();
    phases.vecDiagPhaseToRun.clear();
}


//for (unsigned int i=0; i<vecDiagnosticParticles.size(); i++)
//    computeGlobalDiags(diags->vecDiagnosticParticles[i], timestep);
void SmileiMPI::computeGlobalDiags(DiagnosticParticles* diagParticles, int timestep)
{
    if (timestep % diagParticles->every == diagParticles->time_average-1) {
	MPI_Reduce(diagParticles->filename.size()?MPI_IN_PLACE:&diagParticles->data_sum[0], &diagParticles->data_sum[0], diagParticles->output_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	
	if (smilei_rk==0)
	    diagParticles->write( timestep );
    }
    if (diagParticles->time_average == 1)
	diagParticles->clean();
}

void SmileiMPI::send(Patch* patch, int to, int tag)
{
    for (int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++)
	send( patch->vecSpecies[ispec], to, tag+2*ispec );

    send( patch->EMfields, to, tag+2*patch->vecSpecies.size() );

    send( patch->Diags, to, tag+2*patch->vecSpecies.size()+9 );

}
void SmileiMPI::isend(Patch* patch, int to, int tag)
{
    int nbrOfProp = 7;
    //MPI_Request request;

    for (int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++){
        isend( &(patch->vecSpecies[ispec]->bmax), to, tag+2*ispec+1 );
	//cout << smilei_rk << " sedn " << patch->vecSpecies[ispec]->getNbrOfParticles() << "(" << patch->vecSpecies[ispec]->bmax[patch->vecSpecies[ispec]->bmax.size()-1] << ")" << endl;
        if ( patch->vecSpecies[ispec]->getNbrOfParticles() > 0 ){
            patch->vecSpecies[ispec]->typePartSend = createMPIparticles( patch->vecSpecies[ispec]->particles, nbrOfProp );
            isend( patch->vecSpecies[ispec]->particles, to, tag+2*ispec, patch->vecSpecies[ispec]->typePartSend );
	    MPI_Type_free( &(patch->vecSpecies[ispec]->typePartSend) );
        }
    }
    isend( patch->EMfields, to, tag+2*patch->vecSpecies.size() );
    isend( patch->Diags, to, tag+2*patch->vecSpecies.size()+9 );
}

void SmileiMPI::recv(Patch* patch, int from, int tag)
{
    for (int ispec=0 ; ispec<patch->vecSpecies.size() ; ispec++)
	recv( patch->vecSpecies[ispec], from, tag+2*ispec );

    recv( patch->EMfields, from, tag+2*patch->vecSpecies.size() );

    recv( patch->Diags, from, tag+2*patch->vecSpecies.size()+9 );

}
void SmileiMPI::new_recv(Patch* patch, int from, int tag, Params& params)
{
    int nbrOfProp = 7;
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
	    patch->vecSpecies[ispec]->typePartSend = createMPIparticles( patch->vecSpecies[ispec]->particles, nbrOfProp );
    	    new_recv( patch->vecSpecies[ispec]->particles, from, tag+2*ispec, patch->vecSpecies[ispec]->typePartSend );
	    MPI_Type_free( &(patch->vecSpecies[ispec]->typePartSend) );
	}
    }
   
    recv( patch->EMfields, from, tag+2*patch->vecSpecies.size() );
    recv( patch->Diags, from, tag+2*patch->vecSpecies.size()+9 );

}

void SmileiMPI::send(Species* species, int to, int tag)
{
    if ( species->getNbrOfParticles() )
	send( species->particles, to, tag );
    send( species->bmin, to, tag+1 );
    send( species->bmax, to, tag+2 );
}
//void SmileiMPI::new_send(Species* species, int to, int tag)
//{
//    int nbrOfProp = 7;
//    send( species->bmin, to, tag+1 );
//    send( species->bmax, to, tag+2 );
//    if ( species->getNbrOfParticles() > 0 ){
//        species->typePartSend = createMPIparticles( species->particles, nbrOfProp );
//	new_send( species->particles, to, tag, species->typePartSend );
//    }
//}


void SmileiMPI::recv(Species* species, int from, int tag)
{
    if ( species->getNbrOfParticles() )
	recv( species->particles, from, tag );
    recv( &species->bmin, from, tag+1 );
    recv( &species->bmax, from, tag+2 );
}

void SmileiMPI::send(Particles* particles, int to, int tag)
{
    // Number of properties per particles = nDim_Particles + 3 + 1 + 1
    int nbrOfProp( 7 );
    MPI_Datatype typePartSend = createMPIparticles( particles, nbrOfProp );
    MPI_Send( &(particles->position(0,0)), 1, typePartSend, to, tag, MPI_COMM_WORLD );
    MPI_Type_free( &typePartSend );

}
void SmileiMPI::isend(Particles* particles, int to, int tag, MPI_Datatype typePartSend)
{
    MPI_Request request ;
    MPI_Isend( &(particles->position(0,0)), 1, typePartSend, to, tag, MPI_COMM_WORLD, &request );
}

void SmileiMPI::recv(Particles* particles, int from, int tag)
{
    MPI_Status status;

    // Number of properties per particles = nDim_Particles + 3 + 1 + 1
    int nbrOfProp( 7 );
    MPI_Datatype typePartRecv = createMPIparticles( particles, nbrOfProp );
    MPI_Recv( &(particles->position(0,0)), 1, typePartRecv, from, tag, MPI_COMM_WORLD, &status );
    MPI_Type_free( &typePartRecv );

}
void SmileiMPI::new_recv(Particles* particles, int to, int tag, MPI_Datatype typePartRecv)
{
    MPI_Status status;
    MPI_Recv( &(particles->position(0,0)), 1, typePartRecv, to, tag, MPI_COMM_WORLD, &status );
}

// Assuming vec.size() is known (number of species)
void SmileiMPI::send(std::vector<int> vec, int to, int hindex)
{
    MPI_Send( &(vec[0]), vec.size(), MPI_INT, to, hindex, MPI_COMM_WORLD );

}

// Assuming vec.size() is known (number of species). Asynchronous.
void SmileiMPI::isend(std::vector<int>* vec, int to, int tag)
{
    MPI_Request request; 
    MPI_Isend( &((*vec)[0]), (*vec).size(), MPI_INT, to, tag, MPI_COMM_WORLD, &request );

}

void SmileiMPI::recv(std::vector<int> *vec, int from, int tag)
{
    MPI_Status status;
    MPI_Recv( &((*vec)[0]), vec->size(), MPI_INT, from, tag, MPI_COMM_WORLD, &status );

}

void SmileiMPI::send(ElectroMagn* fields, int to, int tag)
{
    send( fields->Ex_, to, tag+0 );
    send( fields->Ey_, to, tag+1 );
    send( fields->Ez_, to, tag+2 );
    send( fields->Bx_, to, tag+3 );
    send( fields->By_, to, tag+4 );
    send( fields->Bz_, to, tag+5 );
    send( fields->Bx_m, to, tag+6);
    send( fields->By_m, to, tag+7);
    send( fields->Bz_m, to, tag+8);}

void SmileiMPI::isend(ElectroMagn* fields, int to, int tag)
{

    isend( fields->Ex_, to, tag+0);
    isend( fields->Ey_, to, tag+1);
    isend( fields->Ez_, to, tag+2);
    isend( fields->Bx_, to, tag+3);
    isend( fields->By_, to, tag+4);
    isend( fields->Bz_, to, tag+5);
    isend( fields->Bx_m, to, tag+6);
    isend( fields->By_m, to, tag+7);
    isend( fields->Bz_m, to, tag+8);
}

void SmileiMPI::recv(ElectroMagn* fields, int from, int tag)
{
    recv( fields->Ex_, from, tag+0 );
    recv( fields->Ey_, from, tag+1 );
    recv( fields->Ez_, from, tag+2 );
    recv( fields->Bx_, from, tag+3 );
    recv( fields->By_, from, tag+4 );
    recv( fields->Bz_, from, tag+5 );
    recv( fields->Bx_m, from, tag+6 );
    recv( fields->By_m, from, tag+7 );
    recv( fields->Bz_m, from, tag+8 );
}

void SmileiMPI::send(Field* field, int to, int hindex)
{
    MPI_Send( &((*field)(0)),field->globalDims_, MPI_DOUBLE, to, hindex, MPI_COMM_WORLD );
}

void SmileiMPI::isend(Field* field, int to, int hindex)
{
    MPI_Request request;
    MPI_Isend( &((*field)(0)),field->globalDims_, MPI_DOUBLE, to, hindex, MPI_COMM_WORLD, &request );
}

void SmileiMPI::recv(Field* field, int from, int hindex)
{
    MPI_Status status;
    MPI_Recv( &((*field)(0)),field->globalDims_, MPI_DOUBLE, from, hindex, MPI_COMM_WORLD, &status );

}

void SmileiMPI::send( Diagnostic* diags, int to, int tag )
{
    send( diags->probes.probesStart, to, tag );
}

void SmileiMPI::isend( Diagnostic* diags, int to, int tag )
{
    isend( &(diags->probes.probesStart), to, tag );
}

void SmileiMPI::recv( Diagnostic* diags, int from, int tag )
{
    recv( &(diags->probes.probesStart), from, tag );
}
