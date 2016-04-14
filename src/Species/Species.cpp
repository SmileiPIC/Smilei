#include "Species.h"

#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>

#include <omp.h>

// IDRIS
#include <cstring>
// IDRIS
#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "PartBoundCond.h"
#include "PartWall.h"
//#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "Profile.h"

#include "Projector.h"

#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species::Species(Params& params, Patch* patch) :
c_part_max(1),
dynamics_type("norm"), 
time_frozen(0), 
radiating(false), 
ionization_model("none"),
particles(&particles_sorted[0]),
clrw(params.clrw),  
oversize(params.oversize), 
cell_length(params.cell_length), 
velocityProfile(3,NULL),
temperatureProfile(3,NULL),
electron_species(NULL),
nDim_particle(params.nDim_particle),
min_loc_vec(patch->getDomainLocalMin()), 
partBoundCond(NULL),
min_loc(patch->getDomainLocalMin(0)) 
{
    DEBUG(species_type);
    
    PI2 = 2.0 * M_PI;
    
    dx_inv_ = 1./cell_length[0];
    dy_inv_ = 1./cell_length[1];
    
    initCluster(params);

}//END Species creator

void Species::initCluster(Params& params)
{

    // Clusters must all have the same size:
#ifdef _BLABLA
    // -------------------
    // Variable definition
    // -------------------

    // Width of clusters:
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
    if (nDim_particle == 3){
        bmin.resize(params.n_space[0]/clrw*params.n_space[1]);
        bmax.resize(params.n_space[0]/clrw*params.n_space[1]);
    }
    
    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)
    
    //Primal dimension of fields. 
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;
    
    if (nDim_particle == 1){
        b_dim0 =  (1 + clrw) + 2 * oversize[0];
        b_dim1 =  1;
        b_dim2 =  1;
        b_lastdim = b_dim0;
    }
    if (nDim_particle == 2){
        b_dim0 =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim1 =  f_dim1;
        b_dim2 =  1;
        b_lastdim = b_dim1;
    }
    if (nDim_particle == 3){
        b_dim0 =  (1 + clrw) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim1 = f_dim1;
        b_dim2 = f_dim2;
        b_lastdim = b_dim2;
    }
    
    size_proj_buffer = b_dim0*b_dim1*b_dim2; //primal size of a single bufefr.
    
    //Initialize specMPI
    //specMPI.init();
       
    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_lost = 0.;
    nrj_new_particles = 0.;
   
}//END initCluster


// Initialize the operators (Push, Ionize, PartBoundCond)
// This must be separate from the parameters because the Species cloning copies
// the parameters but not the operators.
void Species::initOperators(Params& params, Patch* patch)
{
    // assign the correct Pusher to Push
    Push = PusherFactory::create(params, this);
    
    // Assign the Ionization model (if needed) to Ionize
    //  Needs to be placed after createParticles() because requires the knowledge of max_charge
    // \todo pay attention to restart
    Ionize = IonizationFactory::create(params, this);
    if (Ionize) {
        DEBUG("Species " << species_type << " can be ionized!");
    }
    
    if (Ionize && species_type=="electron") {
        ERROR("Species " << species_type << " can be ionized but species_type='electron'");
    }
    
    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond(params, this, patch);
}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
    delete Push;
    if (Ionize) delete Ionize;
    if (partBoundCond) delete partBoundCond;
    if (chargeProfile) delete chargeProfile;
    if (densityProfile) delete densityProfile;
    for (unsigned int i=0; i<velocityProfile.size(); i++)
        delete velocityProfile[i];
    for (unsigned int i=0; i<temperatureProfile.size(); i++)
        delete temperatureProfile[i];
    
    DEBUG("Species deleted");
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its numerical weight (equivalent to a number density)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initWeight(unsigned int nPart, unsigned int iPart, double density)
{
    for (unsigned  p= iPart; p<iPart+nPart; p++) {
        (*particles).weight(p) = density / nPart;
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its charge state
// ---------------------------------------------------------------------------------------------------------------------
void Species::initCharge(unsigned int nPart, unsigned int iPart, double q)
{
    short Z = (short)q;
    double r = q-(double)Z;
    
    // if charge is integer, then all particles have the same charge
    if ( r == 0. ) {
        for (unsigned int p = iPart; p<iPart+nPart; p++)
            (*particles).charge(p) = Z;
    // if charge is not integer, then particles can have two different charges
    } else {
        int tot = 0, Nm, Np;
        double rr=r/(1-r), diff;
        Np = (int)round(r*(double)nPart);
        Nm = (int)nPart - Np;
        for (unsigned int p = iPart; p<iPart+nPart; p++) {
            if (Np > rr*Nm) {
                (*particles).charge(p) = Z+1;
                Np--;
            } else {
                (*particles).charge(p) = Z;
                Nm--;
            }
            tot += (*particles).charge(p);
        }
        diff = ((double)nPart)*q - (double)tot; // missing charge
        if (diff != 0.) {
            WARNING("Could not match exactly charge="<<q<<" for species "<< species_type <<" (difference of "<<diff<<"). Try to add particles.");
        }
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their position
//   - either using regular distribution in the mesh (initPosition_type = regular)
//   - or using uniform random distribution (initPosition_type = random)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initPosition(unsigned int nPart, unsigned int iPart, double *indexes)
{
    for (unsigned  p= iPart; p<iPart+nPart; p++) {
        for (unsigned  i=0; i<nDim_particle ; i++) {
            
            // define new position (either regular or random)
            if (initPosition_type == "regular") {
                (*particles).position(i,p)=indexes[i]+(p-iPart+0.5)*cell_length[i]/nPart;
            } else if (initPosition_type == "random") {
                (*particles).position(i,p)=indexes[i]+(((double)rand() / RAND_MAX))*cell_length[i];
            }
            //(*particles).position_old(i,p) = (*particles).position(i,p);
        }// i
    }// p
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their momentum
//   - at zero (init_momentum_type = cold)
//   - using random distribution (init_momentum_type = maxwell-juettner)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initMomentum(unsigned int nPart, unsigned int iPart, double *temp, double *vel, vector<double>& max_jutt_cumul)
{
    
    // average mean-momentum (used to center the distribution)
    double pMean[3]= {0.0,0.0,0.0};
    
    if (initMomentum_type == "cold") {
        
        for (unsigned int p= iPart; p<iPart+nPart; p++) {
            for (unsigned int i=0; i<3 ; i++) {
                (*particles).momentum(i,p) = 0.0;
            }
        }
        
    } else if (initMomentum_type == "maxwell-juettner")
    {
        // initialize using the Maxwell-Juettner distribution function
        
        for (unsigned int p= iPart; p<iPart+nPart; p++)
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
        for (unsigned int p= iPart; p<iPart+nPart; p++)
            {
                for (unsigned int i=0; i<3 ; i++) {
                    (*particles).momentum(i,p) -= pMean[i]/nPart;
                }
            }
        
        for (unsigned int p= iPart; p<iPart+nPart; p++) {
            (*particles).momentum(1,p) *= sqrt(temp[1]/temp[0]);
            (*particles).momentum(2,p) *= sqrt(temp[2]/temp[0]);
        }

        
    // Rectangular distribution
    } else if (initMomentum_type == "rectangular") {
        
        for (unsigned int p= iPart; p<iPart+nPart; p++) {
            (*particles).momentum(0,p) = (2.*(double)rand() / RAND_MAX - 1.) * sqrt(temp[0]/mass);
            (*particles).momentum(1,p) = (2.*(double)rand() / RAND_MAX - 1.) * sqrt(temp[1]/mass);
            (*particles).momentum(2,p) = (2.*(double)rand() / RAND_MAX - 1.) * sqrt(temp[2]/mass);
        }
    }//END if initMomentum_type
    
    // Adding the mean velocity (using relativistic composition)
    // ---------------------------------------------------------
    double vx, vy, vz, v2, g, gm1, Lxx, Lyy, Lzz, Lxy, Lxz, Lyz, gp, px, py, pz;
    // mean-velocity
    vx  = -vel[0];
    vy  = -vel[1];
    vz  = -vel[2];
    v2  = vx*vx + vy*vy + vz*vz;
    if ( v2>0. ){
        
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
        for (unsigned int p=iPart; p<iPart+nPart; p++)
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
                       Projector* Proj, Params &params, int diag_flag, PartWalls* partWalls, Patch* patch, SmileiMPI* smpi)
{
    int ithread;
    #ifdef _OPENMP
        ithread = omp_get_thread_num();
    #else
        ithread = 0;
    #endif

    // Ionization current
    LocalFields Jion;
    
    unsigned int iPart;
    
    // Reset list of particles to exchange
    clearExchList();

    int tid(0);
    double ener_iPart(0.);
    std::vector<double> nrj_lost_per_thd(1, 0.);
            
    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if (time_dual>time_frozen) { // moving particle
    
        smpi->dynamics_resize(ithread, nDim_particle, bmax.back());

        //Point to local thread dedicated buffers
        //Still needed for ionization
        std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);

        for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin++) {

            // Interpolate the fields at the particle position
            (*Interp)(EMfields, *particles, smpi, bmin[ibin], bmax[ibin], ithread );

            //Ionization
            if (Ionize){                                
                for (iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) {
                    // Do the ionization (!for testParticles)
                    if ( (*particles).charge(iPart) < (int) atomic_number) {
                        //!\todo Check if it is necessary to put to 0 or if LocalFields ensures it
                        Jion.x=0.0;
                        Jion.y=0.0;
                        Jion.z=0.0;
                        (*Ionize)(*particles, iPart, (*Epart)[iPart], Jion);
                        (*Proj)(EMfields->Jx_, EMfields->Jy_, EMfields->Jz_, *particles, iPart, Jion);
                    }
                }
            }    
                
            // Push the particles
            (*Push)(*particles, smpi, bmin[ibin], bmax[ibin], ithread );
            //for (iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) 
            //    (*Push)(*particles, iPart, (*Epart)[iPart], (*Bpart)[iPart] , (*gf)[iPart]);

            //particles->test_move( bmin[ibin], bmax[ibin], params );


            // Apply wall and boundary conditions
            for (iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) {
                for(unsigned int iwall=0; iwall<partWalls->size(); iwall++) {
                    if ( !(*partWalls)[iwall]->apply(*particles, iPart, this, ener_iPart)) {
                        nrj_lost_per_thd[tid] += mass * ener_iPart;
                    }
                }
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                    addPartInExchList( iPart );
                    //nrj_lost_per_thd[tid] += ener_iPart;
                    nrj_lost_per_thd[tid] += mass * ener_iPart;
                }
             }

            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?

             // Project currents if not a Test species and charges as well if a diag is needed. 
             if (!(*particles).isTest)
                 (*Proj)(EMfields, *particles, smpi, bmin[ibin], bmax[ibin], ithread, ibin, clrw, diag_flag, b_lastdim, ispec );

        }// ibin

        for (int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++)
            nrj_bc_lost += nrj_lost_per_thd[tid];

        // Needs to be reviewed 
        if (Ionize && electron_species) {
            for (unsigned int i=0; i < Ionize->new_electrons.size(); i++) {
                // electron_species->(*particles).push_back(Ionize->new_electrons[i]);
                                
                int ibin = (int) ((Ionize->new_electrons).position(0,i) / cell_length[0]) - ( patch->getCellStartingGlobalIndex(0) + oversize[0] );
                DEBUG("here " << ibin << " " << (Ionize->new_electrons).position(0,i)/(2*M_PI));

                // Copy Ionize->new_electrons(i) in electron_species->particles at position electron_species->bmin[ibin]
                Ionize->new_electrons.cp_particle(i, (*electron_species->particles), electron_species->bmin[ibin] );
                                
                // Update bins status
                // (ugly update, memory is allocated anywhere, OK with vectors per particles parameters)
                electron_species->bmax[ibin]++;
                DEBUG("e- " << i << " to bin " << ibin << " (" <<bmin.size() << "," <<bmax.size()<<")" );
                for (unsigned int ii=ibin+1; ii<bmin.size(); ii++) {
                    electron_species->bmin[ii]++;
                    electron_species->bmax[ii]++;
                }
            }
                        
            // if (Ionize->new_electrons.size())
            //      DEBUG("number of electrons " << electron_species->(*particles).size() << " " << );
            Ionize->new_electrons.clear();
        }
    }
    else { // immobile particle (at the moment only project density)
        if ((diag_flag == 1)&&(!(*particles).isTest)){
            double* b_rho;
            for (unsigned int ibin = 0 ; ibin < bmin.size() ; ibin ++) { //Loop for projection on buffer_proj

                if (params.nDim_field==2)
                    b_rho = &(*EMfields->rho_s[ispec])(ibin*clrw*f_dim1);    
                else if (params.nDim_field==1)
                    b_rho = &(*EMfields->rho_s[ispec])(ibin*clrw);    
                for (iPart=bmin[ibin] ; iPart<bmax[ibin]; iPart++ ) {
                    (*Proj)(b_rho, (*particles), iPart, ibin*clrw, b_lastdim);
                } //End loop on particles
            }//End loop on bins
            
        }
    }//END if time vs. time_frozen

}//END dynamic


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
void Species::count_sort_part(Params &params)
{
    unsigned int ip, npart, ixy,tot, oc, nxy, token;
    int ix,iy;
    double x,y;

    nxy = params.n_space[0]*params.n_space[1];
    token = (particles == &particles_sorted[0]);

    int indices[nxy];
    npart = (*particles).size();
    //particles_sorted = particles ;
    particles_sorted[token].initialize(npart, *particles);

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

        (*particles).overwrite_part(ip, particles_sorted[token] , indices[ixy]);
        indices[ixy]++;
    }

    particles = &particles_sorted[token] ;

}


int Species::createParticles(vector<unsigned int> n_space_to_create, vector<double> cell_index, int new_bin_idx)
{
    // ---------------------------------------------------------
    // Calculate density and number of particles for the species
    // ---------------------------------------------------------
    
    // field containing the charge distribution (always 3d)
    Field3D charge(n_space_to_create);
    max_charge = 0.;
    
    // field containing the density distribution (always 3d)
    Field3D density(n_space_to_create);
    
    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D temperature[3];
    
    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D velocity[3];
    
    // field containing the number of particles in each cell
    Field3D n_part_in_cell(n_space_to_create);

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
                
                n_part_in_cell(i,j,k) = round(ppcProfile->valueAt(x_cell));
                if( n_part_in_cell(i,j,k)<=0. ) {
                    n_part_in_cell(i,j,k) = 0.;
                    density(i,j,k) = 0.;
                    continue;
                }
                
                // assign charge its correct value in the cell
                charge(i,j,k) = chargeProfile->valueAt(x_cell);
                if( charge(i,j,k)>max_charge ) max_charge=charge(i,j,k);
                // assign density its correct value in the cell
                density(i,j,k) = densityProfile->valueAt(x_cell);
                if(density(i,j,k)!=0. && densityProfileType=="charge") {
                    if(charge(i,j,k)==0.) ERROR("Encountered non-zero charge density and zero charge at the same location");
                    density(i,j,k) /= charge(i,j,k);
                }
                density(i,j,k) = abs(density(i,j,k));
                
                // for non-zero density define temperature & mean-velocity and increment the nb of particles
                if (density(i,j,k)!=0.0) {
                    
                    // assign the temperature & mean-velocity their correct value in the cell
                    for (unsigned int m=0; m<3; m++) {
                        temperature[m](i,j,k) = temperatureProfile[m]->valueAt(x_cell);
                        //MESSAGE("temp 1 :" <<  temperature[m](i,j,k))
                        velocity[m](i,j,k) = velocityProfile[m]->valueAt(x_cell);
                    }
                    
                    // increment the effective number of particle by n_part_in_cell(i,j,k)
                    // for each cell with as non-zero density
                    npart_effective += n_part_in_cell(i,j,k);
                    
                }//ENDif non-zero density
                
            }//i
        }//j
    }//k end the loop on all cells
    
    // defines npart_effective for the Species & create the corresponding particles
    // -----------------------------------------------------------------------
    
    // if moving_win
    //     (*particles).create_particles(npart_effective);
    // else {
    //    // reserve included in initialize if particles emty
    //    (*particles).reserve(round( params->species_param[speciesNumber].c_part_max * npart_effective ), ndim);
    //    (*particles).initialize(n_existing_particles+npart_effective, params_->nDim_particle);
    // }
    
    int n_existing_particles = (*particles).size();
    (*particles).initialize(n_existing_particles+npart_effective, nDim_particle);
    
    // define Maxwell-Juettner related quantities
    // ------------------------------------------
    
    // Maxwell-Juettner cumulative function (array)
    std::vector<double> max_jutt_cumul;
    
    /*
    if (initMomentum_type=="maxwell-juettner") {
        //! \todo{Pass this parameters in a code constants class (MG)}
        nE     = 20000;
        muEmax = 20.0;
        
        max_jutt_cumul.resize(nE);
        double mu=mass/temperature[0];
        //double mu=mass/temperature[m](i,j,k);
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
    */
    
    // Initialization of the particles properties
    // ------------------------------------------
    unsigned int nPart;
    unsigned int iPart=n_existing_particles;
    double *indexes=new double[nDim_particle];
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
                    
                    if (initMomentum_type=="maxwell-juettner") {
                        //! \todo{Pass this parameters in a code constants class (MG)}
                        nE     = 20000;
                        muEmax = 20.0;
                        
                        max_jutt_cumul.resize(nE);
                        double mu=mass/temperature[0](i,j,k); // For Temperature profile
                        double Emax=muEmax/mu;
                        dE=Emax/nE;
                        
                        double fl=0;
                        double fr=0;
                        max_jutt_cumul[0]=0.0;
                        for (unsigned int l=1; l<nE; l++ ) {
                            //! \todo{this is just the isotropic case, generalise to non-isotropic (MG)}
                            fr=(1.+l*dE)*sqrt(pow(1.0+l*dE,2)-1.0) * exp(-mu*l*dE);
                            max_jutt_cumul[l]=max_jutt_cumul[l-1] + 0.5*dE*(fr+fl);
                            fl=fr;
                        }
                        for (unsigned int l=0; l<nE; l++) max_jutt_cumul[l]/=max_jutt_cumul[nE-1];
                    }
                    
                    temp[0] = temperature[0](i,j,k);
                    vel[0]  = velocity[0](i,j,k);
                    temp[1] = temperature[1](i,j,k);
                    vel[1]  = velocity[1](i,j,k);
                    temp[2] = temperature[2](i,j,k);
                    vel[2]  = velocity[2](i,j,k);
                    nPart = n_part_in_cell(i,j,k);
                    
                    indexes[0]=i*cell_length[0]+cell_index[0];
                    if (nDim_particle > 1) {
                        indexes[1]=j*cell_length[1]+cell_index[1];
                        if (nDim_particle > 2) {
                            indexes[2]=k*cell_length[2]+cell_index[2];
                        }//nDim_particle > 2
                    }//nDim_particle > 1
                    
                    initPosition(nPart, iPart, indexes);
                    
                    initMomentum(nPart,iPart, temp, vel, max_jutt_cumul);
                    
                    initWeight(nPart, iPart, density(i,j,k));
                    initCharge(nPart, iPart, charge(i,j,k));
                    
                    //calculate new iPart (jump to next cell)
                    iPart+=nPart;
                }//END if density > 0
            }//k end the loop on all cells
        }//j
        if (i%clrw == clrw -1) bmax[new_bin_idx+i/clrw] = iPart;
    }//i
    
    
    delete [] indexes;
    delete [] temp;
    delete [] vel;
    
    // Recalculate former position using the particle velocity
    // (necessary to calculate currents at time t=0 using the Esirkepov projection scheme)
    for (int iPart=n_existing_particles; iPart<n_existing_particles+npart_effective; iPart++) {
        /*897 for (int i=0; i<(int)nDim_particle; i++) {
            (*particles).position_old(i,iPart) -= (*particles).momentum(i,iPart)/(*particles).lor_fac(iPart) * params.timestep;
        }897*/
        nrj_new_particles += (*particles).weight(iPart)*((*particles).lor_fac(iPart)-1.0);
    }
    
    if ((*particles).tracked)
        (*particles).setIds();

    return npart_effective;
    
} // End createParticles


// ------------------------------------------------
// Set position when using restart & moving window
// patch are initialized with t0 position
// ------------------------------------------------
void Species::updateMvWinLimits(double x_moved)
{
    partBoundCond->updateMvWinLimits(x_moved);
    min_loc += x_moved;

} // End updateMvWinLimits


//Do we have to project this species ?
bool Species::isProj(double time_dual, SimWindow* simWindow) {

    return time_dual > time_frozen  || (simWindow && simWindow->isMoving(time_dual)) ;
  
    //Recompute frozen particles density if
    //moving window is activated, actually moving at this time step, and we are not in a density slope.
    /*    bool isproj =(time_dual > species_param.time_frozen  ||
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
            return isproj;*/
    //return time_dual > species_param.time_frozen  || (simWindow && simWindow->isMoving(time_dual)) ;
}
