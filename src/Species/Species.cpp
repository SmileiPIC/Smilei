#include "Species.h"

#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>

#include <omp.h>

#include "PusherFactory.h"
#include "IonizationFactory.h"

#include "PartBoundCond.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"

#include "DensityFactory.h"

#include "Projector.h"

#include "SmileiMPI.h"
#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

using namespace std;

//Obsolete
Species::Species(PicParams& params, int ispec, SmileiMPI* smpi) :
densityProfile(DensityFactory::create(params, ispec)),
speciesNumber(ispec),
cell_length(params.cell_length),
oversize(params.oversize),
ndim(params.nDim_particle),
min_loc(smpi->getDomainLocalMin(0)),
min_loc_vec(smpi->getDomainLocalMin()),
clrw(params.clrw),
species_param(params.species_param[ispec]),
particles(&particles_sorted[0])
//i_domain_begin( smpi->getCellStartingGlobalIndex(0) ),
//j_domain_begin( smpi->getCellStartingGlobalIndex(1) )
{
    specMPI.init();
    
    initSpecies(params);
    initCluster(params);

    if (!params.restart) {
        unsigned int npart_effective=0;
        
        // Create particles in a space starting at cell_index
        vector<double> cell_index(3,0);
        for (unsigned int i=0 ; i<params.nDim_field ; i++) {
            if (cell_length[i]!=0)
	        cell_index[i] = smpi->getDomainLocalMin(i);
        }
        
        int starting_bin_idx = 0;
        // does a loop over all cells in the simulation
        // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
        npart_effective = createParticles(params.n_space, cell_index, starting_bin_idx, params );
        
        //PMESSAGE( 1, smpi->getRank(),"Species "<< speciesNumber <<" # part "<< npart_effective );
    }
	
    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond( params, ispec, smpi );

}

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species::Species(PicParams& params, int ispec, SmileiMPI* smpi, Patch* patch) :
densityProfile(DensityFactory::create(params, ispec)),
speciesNumber(ispec),
cell_length(params.cell_length),
oversize(params.oversize),
ndim(params.nDim_particle),
min_loc(patch->getDomainLocalMin(0)),
min_loc_vec(patch->getDomainLocalMin()),
clrw(params.clrw),
species_param(params.species_param[ispec]),
particles(&particles_sorted[0])

//i_domain_begin( patch->getCellStartingGlobalIndex(0) ),
//j_domain_begin( patch->getCellStartingGlobalIndex(1) )

{
	
	
    initSpecies(params);
    initCluster(params);

    if (!params.restart) {
        unsigned int npart_effective=0;
        
        // Create particles in a space starting at cell_index
        vector<double> cell_index(3,0);
        for (unsigned int i=0 ; i<params.nDim_field ; i++) {
            if (cell_length[i]!=0)
                cell_index[i] = patch->getDomainLocalMin(i);
        }
        
        int starting_bin_idx = 0;
        // does a loop over all cells in the simulation
        // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
	//cout << "Patch start idx = " << cell_index[0] << "\t" << cell_index[1] << endl;
        npart_effective = createParticles(params.n_space, cell_index, starting_bin_idx, params );
        
        //PMESSAGE( 1, smpi->getRank(),"Species "<< speciesNumber <<" # part "<< npart_effective );
    }
    /*double part_min[2], part_max[2];
    for ( int iDim = 0 ; iDim < 2 ; iDim++ ) {
	part_min[iDim] = 1000000.;
	part_max[iDim] = -1.;
	for ( int iPart = 0; iPart < getNbrOfParticles() ; iPart++ ) {
	    if ( particles->Position[iDim][iPart] < part_min[iDim] )
		part_min[iDim] = particles->Position[iDim][iPart];
	    if ( particles->Position[iDim][iPart] > part_max[iDim] )
		part_max[iDim] = particles->Position[iDim][iPart];
	}
	if (getNbrOfParticles()>0)
	    cout << " iDim = " << iDim << "\t" << part_min[iDim] << "\t" << part_max[iDim] << endl;
    }*/

	
	
    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond( params, ispec, smpi, patch);



}//END Species creator


void Species::initSpecies(PicParams& params)
{
    // -------------------
    // Variable definition
    // -------------------
    PI2 = 2.0 * M_PI;
    dx_inv_ = 1./cell_length[0];
    dy_inv_ = 1./cell_length[1];
    DEBUG(species_param.species_type);

    electron_species = NULL;
	
    // atomic number

    // assign the correct Pusher to Push
    Push = PusherFactory::create( params, speciesNumber );

    // assign the Ionization model (if needed) to Ionize
    Ionize = IonizationFactory::create( params, speciesNumber );
    if (Ionize) DEBUG("Species " << speciesNumber << " can be ionized!");

    unsigned int nthds(1);
//#pragma omp parallel shared(nthds) 
//    {
//#ifdef _OMP
//        nthds = omp_get_num_threads();	  
//#endif
//    }
    indexes_of_particles_to_exchange_per_thd.resize(nthds);
    
    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_lost = 0.;
    nrj_new_particles = 0.;

}

void Species::initCluster(PicParams& params)
{

    // Width of clusters:
    // Clusters must all have the same size:
#ifdef _BLABLA
    if (params.n_space[0]%clrw != 0)
        ERROR("Cluster width (clrw) = " << clrw << "should divide n_space[0] = " << params.n_space[0] );
    //Testing if clusters width (clrw) is large enough for the spliting technique:
    if (clrw < 2*oversize[0]+2){
        ERROR("Cluster width (clrw) = "<< clrw << " must be greater than 2*oversize[0]+1 = " << 2*oversize[0]+1 );
    }
#endif

    // Arrays of the min and max indices of the particle bins
    bmin.resize(params.n_space[0]/clrw);
    bmax.resize(params.n_space[0]/clrw);
    if (ndim == 3){
        bmin.resize(params.n_space[0]/clrw*params.n_space[1]);
        bmax.resize(params.n_space[0]/clrw*params.n_space[1]);
    }
	
    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)
    
    //Primal dimension of fields. 
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;
    
    if (ndim == 1){
        b_dim0 =  (1 + clrw) + 2 * oversize[0];
        b_dim1 =  1;
        b_dim2 =  1;
        b_lastdim = b_dim0;
    }
    if (ndim == 2){
        b_dim0 =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim1 =  f_dim1;
        b_dim2 =  1;
        b_lastdim = b_dim1;
    }
    if (ndim == 3){
        b_dim0 =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim1 = f_dim1;
        b_dim2 = f_dim2;
        b_lastdim = b_dim2;
    }
    
    size_proj_buffer = b_dim0*b_dim1*b_dim2; //primal size of a single bufefr.

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
    delete Push;
    if (Ionize) delete Ionize;
    if (partBoundCond) delete partBoundCond;
    if (densityProfile) delete densityProfile;
    
    DEBUG(10,"Species deleted ");
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its numerical weight (equivalent to a number density)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initWeight(PicParams* params, unsigned int ispec, unsigned int iPart, double density)
{
    for (unsigned  p= iPart; p<iPart+params->species_param[ispec].n_part_per_cell; p++) {
        (*particles).weight(p) = density / params->species_param[ispec].n_part_per_cell;
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its charge state
// ---------------------------------------------------------------------------------------------------------------------
void Species::initCharge(PicParams* params, unsigned int ispec, unsigned int iPart, double density)
{
    for (unsigned  p= iPart; p<iPart+params->species_param[ispec].n_part_per_cell; p++) {
        (*particles).charge(p) = params->species_param[ispec].charge;
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their position
//   - either using regular distribution in the mesh (regular)
//   - or using uniform random distribution (for cold and maxwell-juettner distribution)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initPosition(unsigned int np, unsigned int iPart, double *indexes, unsigned int ndim,
                           std::vector<double> cell_length, string initialization_type)
{
    for (unsigned  p= iPart; p<iPart+np; p++)
	{
	    for (unsigned  i=0; i<ndim ; i++)
		{
		    if (initialization_type == "regular") {
                            (*particles).position(i,p)=indexes[i]+(p-iPart+0.5)*cell_length[i]/np;
		    } else if (initialization_type == "cold" || initialization_type == "maxwell-juettner") {
                        (*particles).position(i,p)=indexes[i]+(((double)rand() / RAND_MAX))*cell_length[i];
		    }
		    (*particles).position_old(i,p) = (*particles).position(i,p);
		}// i
	}// p
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their momentum
//   - at zero if regular or cold
//   - using random distribution if maxwell-juettner
// ---------------------------------------------------------------------------------------------------------------------
void Species::initMomentum(unsigned int np, unsigned int iPart, double *temp, double *vel, string initialization_type,
                           vector<double>& max_jutt_cumul)
{
	
    // average mean-momentum (used to center the distribution)
    double pMean[3]= {0.0,0.0,0.0};
	
	
    if (initialization_type == "regular" || initialization_type == "cold")
	{
	    // initialize momentum at 0 for regular or cold initialization type
	    for (unsigned int p= iPart; p<iPart+np; p++) {
            for (unsigned int i=0; i<3 ; i++) {
                (*particles).momentum(i,p) = 0.0;
            }
			
			
	    }
		
	} else if (initialization_type == "maxwell-juettner")
	{
	    // initialize using the Maxwell-Juettner distribution function
	    for (unsigned int p= iPart; p<iPart+np; p++)
		{
		    double Renergy=(double)rand() / RAND_MAX;
		    double phi=acos(1.0-2.0*(double)rand() / RAND_MAX);
		    double theta=2.0*M_PI*(double)rand() / RAND_MAX;
			
		    int il=0;
		    int ir=max_jutt_cumul.size();
		    while (ir > il+1)  {
                int im=(il+ir)/2;
                if (Renergy > max_jutt_cumul[im]) {
                    il=im;
                } else {
                    ir=im;
                }
		    }
		    double right_w=(Renergy-max_jutt_cumul[il])/(max_jutt_cumul[il+1]);
		    double left_w=1-right_w;
			
		    double Ener=left_w*il*dE +right_w*(il+1)*dE;
		    double psm = sqrt(pow(1.0+Ener,2)-1.0);
			
		    (*particles).momentum(0,p) = psm*cos(theta)*sin(phi);
		    (*particles).momentum(1,p) = psm*sin(theta)*sin(phi);
		    (*particles).momentum(2,p) = psm*cos(phi);
		    for (unsigned int i=0; i<3 ; i++)
			{
			    pMean[i] += (*particles).momentum(i,p);
			}
		}//p
		
	    // center the distribution function around pMean
	    // \todo{Allow for non-zero mean-velocity (MG)}
	    for (unsigned int p= iPart; p<iPart+np; p++)
		{
		    for (unsigned int i=0; i<3 ; i++) {
                (*particles).momentum(i,p) -= pMean[i]/np;
		    }
		}
		
    }//END if initialization_type
    
    
    // Adding the mean velocity (using relativistic composition)
    if ( (vel[0]!=0.0) || (vel[1]!=0.0) || (vel[2]!=0.0) ){
        
        double vx, vy, vz, v2, g, gm1, Lxx, Lyy, Lzz, Lxy, Lxz, Lyz, gp, px, py, pz;
        
        // mean-velocity
        vx  = -vel[0];
        vy  = -vel[1];
        vz  = -vel[2];
        v2  = vx*vx + vy*vy + vz*vz;
        g   = 1.0/sqrt(1.0-v2);
        gm1 = g - 1.0;
        
        // compute the different component of the Matrix block of the Lorentz transformation
        Lxx = 1.0 + gm1 * vx*vx/v2;
        Lyy = 1.0 + gm1 * vy*vy/v2;
        Lzz = 1.0 + gm1 * vz*vz/v2;
        Lxy = gm1 * vx*vy/v2;
        Lxz = gm1 * vx*vz/v2;
        Lyz = gm1 * vy*vz/v2;
        
        // Lorentz transformation of the momentum
        for (unsigned int p=iPart; p<iPart+np; p++)
	    {
            gp = sqrt(1.0 + pow((*particles).momentum(0,p),2) + pow((*particles).momentum(1,p),2) + pow((*particles).momentum(2,p),2));
            px = -gp*g*vx + Lxx * (*particles).momentum(0,p) + Lxy * (*particles).momentum(1,p) + Lxz * (*particles).momentum(2,p);
            py = -gp*g*vy + Lxy * (*particles).momentum(0,p) + Lyy * (*particles).momentum(1,p) + Lyz * (*particles).momentum(2,p);
            pz = -gp*g*vz + Lxz * (*particles).momentum(0,p) + Lyz * (*particles).momentum(1,p) + Lzz * (*particles).momentum(2,p);
            (*particles).momentum(0,p) = px;
            (*particles).momentum(1,p) = py;
            (*particles).momentum(2,p) = pz;
	    }
        
    }//ENDif vel != 0
	
    
	
}//END initMomentum


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamics(double time_dual, unsigned int ispec, ElectroMagn* EMfields, Interpolator* Interp,
                       Projector* Proj, PicParams &params, int diag_flag)
{
    //Interpolator* LocInterp = InterpolatorFactory::create(params, smpi, NULL);
    
    // Electric field at the particle position
    LocalFields Epart;
    // Magnetic field at the particle position
    LocalFields Bpart;
    // Ionization current
    LocalFields Jion;
	
    int iloc,jloc;
    unsigned int i,j,ibin,iPart;
    
    //! buffers for currents and charge
    double *b_Jx,*b_Jy,*b_Jz,*b_rho;

    // Reset list of particles to exchange
    int tid(0);
    double gf = 1.0;
    std::vector<double> nrj_lost_per_thd(1, 0.);
//#ifdef _OMP
//    tid = omp_get_thread_num();
//    int nthds = omp_get_num_threads();
//    nrj_lost_per_thd.resize(nthds, 0.);
//#endif
    clearExchList(tid);
    	
    //ener_tot  = 0.;
    //ener_lost = 0.;
    double ener_iPart(0.);
    //bool contribute(true);
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if (time_dual>species_param.time_frozen) { // moving particle

	//Allocate buffer for projection  *****************************
	// *4 accounts for Jy, Jz and rho. * nthds accounts for each thread.
	//b_Jx = (double *) malloc(4 * size_proj_buffer * sizeof(double));
	//Point buffers of each thread to the correct position
	//b_Jy = b_Jx + size_proj_buffer ;
	//b_Jz = b_Jy + size_proj_buffer ;
	//b_rho = b_Jz + size_proj_buffer ;

        for (ibin = 0 ; ibin < bmin.size() ; ibin++) {

            //memset( &(b_Jx[0]), 0, 4*size_proj_buffer*sizeof(double)); 
            if (diag_flag == 0){
	        b_Jx =  &(*EMfields->Jx_ )(ibin*clrw*f_dim1);
	        b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(f_dim1+1));
	        b_Jz =  &(*EMfields->Jz_ )(ibin*clrw*f_dim1);

            } else { 
	        b_Jx =  &(*EMfields->Jx_s[ispec] )(ibin*clrw*f_dim1);
	        b_Jy =  &(*EMfields->Jy_s[ispec] )(ibin*clrw*(f_dim1+1));
	        b_Jz =  &(*EMfields->Jz_s[ispec] )(ibin*clrw*f_dim1);
	        b_rho = &(*EMfields->rho_s[ispec])(ibin*clrw*f_dim1);
            }

            for (iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) {
				
                // Interpolate the fields at the particle position
                //(*LocInterp)(EMfields, *particles, iPart, &Epart, &Bpart);
                (*Interp)(EMfields, *particles, iPart, &Epart, &Bpart);
				
                // Do the ionization
                if (Ionize && (*particles).charge(iPart) < (int) species_param.atomic_number) {
                    //!\todo Check if it is necessary to put to 0 or if LocalFields ensures it
                    Jion.x=0.0;
                    Jion.y=0.0;
                    Jion.z=0.0;
                    (*Ionize)(*particles, iPart, Epart, Jion);
                    (*Proj)(EMfields->Jx_, EMfields->Jy_, EMfields->Jz_, *particles, iPart, Jion);
                }
				
				
                // Push the particle
                (*Push)(*particles, iPart, Epart, Bpart, gf);
				
                // Apply boundary condition on the particles
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is no more in the domain local
                //	if omp, create a list per thread
                if ( !partBoundCond->apply( *particles, iPart, params.species_param[ispec], ener_iPart ) ) {
                    addPartInExchList( tid, iPart );
		    nrj_lost_per_thd[tid] += ener_iPart;
                }

		if (diag_flag == 0){ 
		    (*Proj)(b_Jx , b_Jy , b_Jz , *particles,  iPart, gf, ibin*clrw, b_lastdim);
                } else {
		    (*Proj)(b_Jx , b_Jy , b_Jz ,b_rho, *particles,  iPart, gf, ibin*clrw, b_lastdim);
                }

            }//iPart

                //for (i = 0; i < b_dim0 ; i++) {
                //    iloc = ibin*clrw + i ;
                //    //! \todo Here b_dim0 is the dual size. Make sure no problems arise when i == b_dim0-1 for primal arrays.
                //    for (j = 0; j < b_dim1 ; j++) {
                //        (*EMfields->Jx_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jx[i*b_dim1+j];   //  primal along y
                //        (*EMfields->Jy_s[ispec]) (iloc*(f_dim1+1)+j) +=  b_Jy[i*b_dim1+j];   //+1 because dual along y
                //        (*EMfields->Jz_s[ispec]) (iloc*(f_dim1  )+j) +=  b_Jz[i*b_dim1+j];   // primal along y
                //        (*EMfields->rho_s[ispec])(iloc*(f_dim1  )+j) += b_rho[i*b_dim1+j];   // primal along y
                //    }
                //}
        }// ibin
        //free(b_Jx);

	for (int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++)
	    nrj_bc_lost += nrj_lost_per_thd[tid];
       

        // Needs to be reviewed 
        //if (Ionize && electron_species) {
        //    for (unsigned int i=0; i < Ionize->new_electrons.size(); i++) {
        //        // electron_species->particles.push_back(Ionize->new_electrons[i]);
	//			
        //        int ibin = (int) ((Ionize->new_electrons).position(0,i) / cell_length[0]) - ( smpi->getCellStartingGlobalIndex(0) + oversize[0] );
        //        DEBUG("here " << ibin << " " << (Ionize->new_electrons).position(0,i)/(2*M_PI));
        //        // Copy Ionize->new_electrons(i) in electron_species->particles at position electron_species->bmin[ibin]
        //        Ionize->new_electrons.cp_particle(i, (*electron_species->particles), electron_species->bmin[ibin] );
	//			
        //        // Update bins status
        //        // (ugly update, memory is allocated anywhere, OK with vectors per particles parameters)
        //        electron_species->bmax[ibin]++;
        //        for (int i=ibin+1; i<bmin.size(); i++) {
        //            electron_species->bmin[i]++;
        //            electron_species->bmax[i]++;
        //        }
        //        DEBUG("here");
        //    }
	//		
        //    // if (Ionize->new_electrons.size())
        //    //      DEBUG("number of electrons " << electron_species->particles.size() << " " << );
        //    Ionize->new_electrons.clear();
        //}
    }
    else { // immobile particle (at the moment only project density)
        if (diag_flag == 1){
	    //b_rho = (double *) malloc(size_proj_buffer * sizeof(double));

            for (ibin = 0 ; ibin < bmin.size() ; ibin ++) { //Loop for projection on buffer_proj

		b_rho = &(*EMfields->rho_s[ispec])(ibin*clrw*f_dim1);    
		//memset( &(b_rho[0]), 0, size_proj_buffer*sizeof(double)); 
                for (iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) {
                    //Update position_old because it is required for the Projection.
                    for ( int i = 0 ; i<ndim ; i++ ) {
                        (*particles).position_old(i, iPart)  = (*particles).position(i, iPart);
                    }

                    (*Proj)(b_rho, (*particles), iPart, ibin*clrw, b_lastdim);
                } //End loop on particles
            }//End loop on bins
	    //free(b_rho);
        }
    }//END if time vs. time_frozen
    //delete LocInterp;
}//END dynamic

// ---------------------------------------------------------------------------------------------------------------------
// Dump for the Species
//! \todo{Check if one wants to keep it like this (MG)}
// ---------------------------------------------------------------------------------------------------------------------
void Species::dump(std::ofstream& ofile)
{
    for (unsigned int i=0; i<(*particles).size(); i++ )
	{
	    ofile << i ;
	    for (unsigned int m=0; m<ndim; m++) ofile << "\t" << (*particles).position(m,i);
	    for (unsigned int m=0; m<3; m++)    ofile << "\t" << (*particles).momentum(m,i);
	    ofile << "\t" << (*particles).weight(i); //<< "\t" << Push->getMass() << "\t" << Push->getCharge();
	    ofile << endl;
	}
    ofile << endl;
}

// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::sort_part()
{
    //The width of one bin is cell_length[0] * clrw.
	
    int p1,p2,bmin_init;
    unsigned int bin;
    double limit;

	
    //Backward pass
    for (bin=0; bin<bmin.size()-1; bin++) { //Loop on the bins. 
        limit = min_loc + (bin+1)*cell_length[0]*clrw;
        p1 = bmax[bin]-1;
        //If first particles change bin, they do not need to be swapped.
        while (p1 == bmax[bin]-1 && p1 >= bmin[bin]) {
            if ((*particles).position(0,p1) >= limit ) {
                bmax[bin]--;
            }
            p1--;
        }
        //         Now particles have to be swapped
        for( p2 = p1 ; p2 >= bmin[bin] ; p2-- ) { //Loop on the bin's particles.
            if ((*particles).position(0,p2) >= limit ) {
                //This particle goes up one bin.
                (*particles).swap_part(p2,bmax[bin]-1);
                bmax[bin]--;
            }
        }
    }
    //Forward pass + Rebracketting
    for (bin=1; bin<bmin.size(); bin++) { //Loop on the bins. 
        limit = min_loc + bin*cell_length[0]*clrw;
        bmin_init = bmin[bin];
        p1 = bmin[bin];
        while (p1 == bmin[bin] && p1 < bmax[bin]) {
            if ((*particles).position(0,p1) < limit ) {
                bmin[bin]++;
            }
            p1++;
        }
        for( p2 = p1 ; p2 < bmax[bin] ; p2++ ) { //Loop on the bin's particles.
            if ((*particles).position(0,p2) < limit ) {
                //This particle goes down one bin.
                (*particles).swap_part(p2,bmin[bin]);
                bmin[bin]++;
            }
        }
		
        //Rebracketting
        //Number of particles from bin going down is: bmin[bin]-bmin_init.
        //Number of particles from bin-1 going up is: bmin_init-bmax[bin-1].
        //Total number of particles we need to swap is the min of both.
        p2 = min(bmin[bin]-bmin_init,bmin_init-bmax[bin-1]);
        if (p2 >0) (*particles).swap_part(bmax[bin-1],bmin[bin]-p2,p2);
        bmax[bin-1] += bmin[bin] - bmin_init;
        bmin[bin] = bmax[bin-1];
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::count_sort_part(PicParams &params)
{

unsigned int ip, npart, ndim,ixy,tot, oc, nxy, bin, token;
int ix,iy;
double x,y;

nxy = params.n_space[0]*params.n_space[1];
token = (particles == &particles_sorted[0]);

int indices[nxy];
npart = (*particles).size();
//particles_sorted = particles ;
particles_sorted[token].initialize(npart, params.nDim_particle);

for (unsigned int i=0; i < nxy ; i++) indices[i] = 0 ;

// first loop counts the # of particles in each cell
for (ip=0; ip < npart; ip++)
{
    x = (*particles).position(0,ip)-min_loc;
    y = (*particles).position(1,ip)-min_loc_vec[1];

    ix = floor(x * dx_inv_) ;
    iy = floor(y * dy_inv_) ;

    ixy = iy + ix*params.n_space[1];


    indices[ixy] ++;
}

// second loop convert the count array in cumulative sum
tot=0;
for (unsigned int ixy=0; ixy < nxy; ixy++)
{
    oc = indices[ixy];
    indices[ixy] = tot;
    tot += oc;
}

//Bookmarking is not needed if normal sort is called before.
//bmin[0] = 0;    
//for (bin=0; bin<bmin.size()-1; bin++) { //Loop on the bins. 
//
//    bmin[bin+1] = indices[(bin+1)*params.n_space[1]*clrw] ;   
//    bmax[bin] = bmin[bin+1];   
//}
//bin = bmin.size()-1 ;
//bmax[bin] = npart;   

// last loop puts the particles and update the count array
for (ip=0; ip < npart; ip++) {
    x = (*particles).position(0,ip)-min_loc;
    y = (*particles).position(1,ip)-min_loc_vec[1];

    ix = floor(x * dx_inv_) ;
    iy = floor(y * dy_inv_) ;

    ixy = iy + ix*params.n_space[1];

    (*particles).overwrite_part2D(ip, particles_sorted[token] , indices[ixy]);
    indices[ixy]++;
}

particles = &particles_sorted[token] ;

}


//Obsolete
void Species::movingWindow_x(unsigned int shift, SmileiMPI *smpi, PicParams& params)
{
    //i_domain_begin+=shift;
    // Update BC positions
    partBoundCond->moveWindow_x( shift*cell_length[0], smpi );
    // Set for bin managment
    //min_loc += shift*cell_length[0];
    
    // Send particles of first bin on process rank-1
    // If no rank-1 -> particles deleted
    //clearExchList(0);
    //for (unsigned int ibin = 0 ; ibin < 1 ; ibin++)
    //    for (unsigned int iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) {
    //        addPartInExchList( 0, iPart );
    //        nrj_mw_lost += (*particles).weight(iPart)*((*particles).lor_fac(iPart)-1.0);
    //    }
    //
    //// bin 0 empty
    //// Shifts all the bins by 1. 
    //bmin.erase( bmin.begin() );
    //bmax.erase( bmax.begin() );
    //// Create new bin at the end
    //// Update last values of bmin and bmax to handle correctly received particles
    //bmin.push_back( bmax[bmax.size()-1] );
    //bmax.push_back( bmax[bmax.size()-1] );
    //bmin[0] = 0;

    //int iDim(0);    
    //smpi->exchangeParticles( this, speciesNumber,params, 0, iDim );
    //
    //// Create new particles
    //if (smpi->isEastern() ) {
    //    defineNewCells(shift, smpi, params);
    //}
    
}

void Species::defineNewCells(unsigned int shift, SmileiMPI *smpi, PicParams& params)
{
    // does a loop over all cells in the simulation
    // considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]
    vector<double> cell_index(3,0);
    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        if (cell_length[i]!=0) {
            cell_index[i] = smpi->getDomainLocalMin(i);
        }
    }
    // cell_index[0] goes to end of the domain minus cell to create
    cell_index[0] = smpi->getDomainLocalMax(0) - shift*cell_length[0];
    
    // Next bin to create
    int new_bin_idx = bmin.size() - 1;
    
    vector<unsigned int> n_space_created(3,0);
    n_space_created[0] = shift;
    n_space_created[1] = params.n_space[1];
    n_space_created[2] = params.n_space[2];
    
    unsigned int npart_effective = createParticles(n_space_created, cell_index, new_bin_idx, params );
}


int Species::createParticles(vector<unsigned int> n_space_to_create, vector<double> cell_index, int new_bin_idx, PicParams& params  )
{
    // ---------------------------------------------------------
    // Calculate density and number of particles for the species
    // ---------------------------------------------------------
    
    // field containing the density distribution (always 3d)
    Field3D density(n_space_to_create);
	
    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D temperature[3];
	
    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D velocity[3];
	
    for (unsigned int i=0; i<3; i++) {
        velocity[i].allocateDims(n_space_to_create);
        temperature[i].allocateDims(n_space_to_create);
    }
    
    int npart_effective = 0;
    for (unsigned int i=0; i<n_space_to_create[0]; i++) {
        for (unsigned int j=0; j<n_space_to_create[1]; j++) {
            for (unsigned int k=0; k<n_space_to_create[2]; k++) {
                
                vector<double> x_cell(3,0);
                x_cell[0] = cell_index[0] + (i+0.5)*cell_length[0];
                x_cell[1] = cell_index[1] + (j+0.5)*cell_length[1];
                x_cell[2] = cell_index[2] + (k+0.5)*cell_length[2];
                
                // assign density its correct value in the cell
                density(i,j,k) = species_param.density
                *                (*densityProfile)(x_cell);
                
                // for non-zero density define temperature & mean-velocity and increment the nb of particles
                if (density(i,j,k)!=0.0) {
                    
                    // assign the temperature & mean-velocity their correct value in the cell
                    for (unsigned int m=0; m<3; m++)	{
                        temperature[m](i,j,k) = species_param.temperature[m];
                        velocity[m](i,j,k) = species_param.mean_velocity[m];
                    }
                    
                    // increment the effective number of particle by n_part_per_cell
                    // for each cell with as non-zero density
                    npart_effective += species_param.n_part_per_cell;
                    //DEBUG(10,"Specie "<< speciesNumber <<" # part "<<npart_effective<<" "<<i<<" "<<j<<" "<<k);
                    
                }//ENDif non-zero density
                
            }//i
        }//j
    }//k end the loop on all cells
    
    // defines npart_effective for the Species & create the corresponding particles
    // -----------------------------------------------------------------------
    
    // if moving_win
    //     particles.create_particles(npart_effective);
    // else {
    //    // reserve included in initialize if particles emty
    //    particles.reserve(round( params->species_param[speciesNumber].c_part_max * npart_effective ), ndim);
    //    particles.initialize(n_existing_particles+npart_effective, params_->nDim_particle);
    // }
    
    int n_existing_particles = (*particles).size();
    (*particles).initialize(n_existing_particles+npart_effective, params.nDim_particle);
    
    
    // define Maxwell-Juettner related quantities
    // ------------------------------------------
    
    // Maxwell-Juettner cumulative function (array)
    std::vector<double> max_jutt_cumul;
    
    if (species_param.initialization_type=="maxwell-juettner") {
        //! \todo{Pass this parameters in a code constants class (MG)}
        nE     = 20000;
        muEmax = 20.0;
        
        max_jutt_cumul.resize(nE);
        double mu=species_param.mass/species_param.temperature[0];
        double Emax=muEmax/mu;
        dE=Emax/nE;
        
        double fl=0;
        double fr=0;
        max_jutt_cumul[0]=0.0;
        for (unsigned  i=1; i<nE; i++ ) {
            //! \todo{this is just the isotropic case, generalise to non-isotropic (MG)}
            fr=(1+i*dE)*sqrt(pow(1.0+i*dE,2)-1.0) * exp(-mu*i*dE);
            max_jutt_cumul[i]=max_jutt_cumul[i-1] + 0.5*dE*(fr+fl);
            fl=fr;
        }
        for (unsigned int i=0; i<nE; i++) max_jutt_cumul[i]/=max_jutt_cumul[nE-1];
        
    }
    
    
    // Initialization of the particles properties
    // ------------------------------------------
    unsigned int iPart=n_existing_particles;
    double *indexes=new double[params.nDim_particle];
    double *temp=new double[3];
    double *vel=new double[3];
    
    // start a loop on all cells
    
    //bmin[bin] point to begining of bin (first particle)
    //bmax[bin] point to end of bin (= bmin[bin+1])
    //if bmax = bmin, bin is empty of particle.
    
    for (unsigned int i=0; i<n_space_to_create[0]; i++) {
        if (i%clrw == 0) bmin[new_bin_idx+i/clrw] = iPart;
        for (unsigned int j=0; j<n_space_to_create[1]; j++) {
            for (unsigned int k=0; k<n_space_to_create[2]; k++) {
                // initialize particles in meshes where the density is non-zero
                if (density(i,j,k)>0) {
                    
                    temp[0] = temperature[0](i,j,k);
                    vel[0]  = velocity[0](i,j,k);
                    temp[1] = temperature[1](i,j,k);
                    vel[1]  = velocity[1](i,j,k);
                    temp[2] = temperature[2](i,j,k);
                    vel[2]  = velocity[2](i,j,k);
                    
                    indexes[0]=i*cell_length[0]+cell_index[0];
                    if (ndim > 1) {
                        indexes[1]=j*cell_length[1]+cell_index[1];
                        if (ndim > 2) {
                            indexes[2]=k*cell_length[2]+cell_index[2];
                        }//ndim > 2
                    }//ndim > 1
                    
                    initPosition(species_param.n_part_per_cell,iPart, indexes, params.nDim_particle,
                                 cell_length, species_param.initialization_type);
                    initMomentum(species_param.n_part_per_cell,iPart, temp, vel,
                                 species_param.initialization_type, max_jutt_cumul);
                    initWeight(&params, speciesNumber, iPart, density(i,j,k));
                    initCharge(&params, speciesNumber, iPart, density(i,j,k));
                    
                    //calculate new iPart (jump to next cell)
                    iPart+=species_param.n_part_per_cell;
                }//END if density > 0
            }//k end the loop on all cells
        }//j
        if (i%clrw == clrw -1) bmax[new_bin_idx+i/clrw] = iPart;
    }//i
    
    delete [] indexes;
    delete [] temp;
    delete [] vel;
    
    for (unsigned int iPart=n_existing_particles; iPart<n_existing_particles+npart_effective; iPart++) {
	nrj_new_particles += (*particles).weight(iPart)*((*particles).lor_fac(iPart)-1.0);
    }

    return npart_effective;
    
}

// Obsolete ...
void Species::updateMvWinLimits(double x_moved) {
    partBoundCond->updateMvWinLimits(x_moved);
    min_loc += x_moved;
}
//Do we have to project this species ?
 bool Species::isProj(double time_dual, SimWindow* simWindow) {
    bool isproj;
    //Recompute frozen particles density if
    //moving window is activated, actually moving at this time step, and we are not in a density slope.
    isproj =(time_dual > species_param.time_frozen  ||
                 (simWindow && simWindow->isMoving(time_dual) &&
                     (species_param.species_geometry == "gaussian" ||
                         (species_param.species_geometry == "trapezoidal" &&
                            //Before end of density ramp up.
                            (simWindow->getXmoved() < species_param.vacuum_length[0] + species_param.dens_length_x[1] + oversize[0]*cell_length[0] || 
                            //After begining of density ramp down. 
                            simWindow->getXmoved() +  simWindow->getNspace_win_x()*cell_length[0] > species_param.vacuum_length[0] + species_param.dens_length_x[1]+ species_param.dens_length_x[0]
                            )
                        )
                    ) 
                )
            );
    return isproj;
}
