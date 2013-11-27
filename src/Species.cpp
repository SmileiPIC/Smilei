
#include "Species.h"

#include "ParticleFactory.h"
#include "PusherFactory.h"
#include "IonizationFactory.h"

#include "PartBoundCond.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "Projector.h"

#include "SmileiMPI.h"

#include "Field3D.h"
#include "Tools.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species::Species(PicParams* params, int ispec, SmileiMPI* smpi) {

    // -------------------
	// Variable definition
	// -------------------
    
	name_str=params->species_param[ispec].species_type;
    
	DEBUG(name_str);
	
    electron_species = NULL;
    
    // species index
    speciesNumber = ispec;
    
	// number of spatial dimensions for the particles
	ndim = params->nDim_particle;
	
	// time over which particles remain frozen
	time_frozen = params->species_param[ispec].time_frozen;
    
    // atomic number
    atomic_number = params->species_param[ispec].atomic_number;
    
	// field containing the density distribution (always 3d)
	Field3D density(params->n_space);
    
	// field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
	Field3D temperature[3];
    
	// field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
	Field3D velocity[3];
    
	for (unsigned int i=0;i<3;i++) {
		velocity[i].allocateDims(params->n_space);
		temperature[i].allocateDims(params->n_space);
	}

        // Arrays of the min and max indices of the particle bins
        bmin.resize(params->n_space[ndim-1]);
        bmax.resize(params->n_space[ndim-1]);
    
        //number of cell per bin
        int ncell_per_bin;

        switch (ndim){
            case 1:
                ncell_per_bin = 1;
            break;
            case 2:
                ncell_per_bin = params->n_space[0];
            break;
            case 3:
                ncell_per_bin = params->n_space[0]*params->n_space[1];
            break;
        } 
        
        //initialisation of the indices
        for (unsigned int i=0;i<bmin.size();i++) {

            bmin[i]=i*ncell_per_bin*params->species_param[ispec].n_part_per_cell;
            bmax[i]=(i+1)*ncell_per_bin*params->species_param[ispec].n_part_per_cell-1;
        } 
	
    // ---------------------------------------------------------
	// Calculate density and number of particles for the species
	// ---------------------------------------------------------
    
	npart_effective=0;

	// does a loop over all cells in the simulation
	// considering a 3d volume with size n_space[0]*n_space[1]*n_space[2]

	vector<int> cell_index(3,0);
	for (unsigned int i=0 ; i<params->nDim_field ; i++) {
		if (params->cell_length[i]!=0)
			cell_index[i] = round (smpi->getDomainLocalMin(i)/params->cell_length[i]);
	}

	for (unsigned int k=0; k<params->n_space[2]; k++) {
		for (unsigned int j=0; j<params->n_space[1]; j++) {
			for (unsigned int i=0; i<params->n_space[0]; i++) {
                
				// ------------------------
				// Constant density profile
				// ------------------------
				if (params->plasma_geometry=="constant") {

					if (((params->cell_length[0]==0.0) || (
						 (cell_index[0]+i+0.5)*params->cell_length[0] > params->vacuum_length[0] &&
						 (cell_index[0]+i+0.5)*params->cell_length[0] < params->vacuum_length[0]+params->plasma_length[0]
						 )) &&
						((params->cell_length[1]==0.0) || (
						 (cell_index[1]+j+0.5)*params->cell_length[1] > params->vacuum_length[1] &&
						 (cell_index[1]+j+0.5)*params->cell_length[1] < params->vacuum_length[1]+params->plasma_length[1]
						 )) &&
						((params->cell_length[2]==0.0) || (
						 (cell_index[2]+k+0.5)*params->cell_length[2] > params->vacuum_length[2] &&
						 (cell_index[2]+k+0.5)*params->cell_length[2] < params->vacuum_length[2]+params->plasma_length[2])
						)) {

						// assign density its correct value in the cell
						density(i,j,k) = params->species_param[ispec].density;

                        
						// assign the temperature & mean-velocity their correct value in the cell
						for (unsigned int m=0; m<3;m++)	{
							temperature[m](i,j,k) = params->species_param[ispec].temperature[m];
							velocity[m](i,j,k) = params->species_param[ispec].mean_velocity[m];
						}
                        
						// increment the effective number of particle by n_part_per_cell
						// for each cell with as non-zero density
						npart_effective += params->species_param[ispec].n_part_per_cell;
						//DEBUG(10,"Specie "<< ispec <<" # part "<<npart_effective<<" "<<i<<" "<<j<<" "<<k);
						
					} else {
						density(i,j,k)=0.0;
					}
					
				} else { // plasma_geometry
					
					// ---------------------------
					// Not defined density profile
					// ---------------------------
					ERROR("geometry not implemented");
					
				}//END if plasma_geometry
			}//i
		}//j
	}//k end the loop on all cells

    // defines npart_effective for the Species & create the corresponding particles
    // -----------------------------------------------------------------------
    params->species_param[ispec].n_part_max = round( params->species_param[ispec].c_part_max * npart_effective );
	particles.reserve(params->species_param[ispec].n_part_max);
    particles = ParticleFactory::createVector(params, ispec, npart_effective);
    
    // define Maxwell-Juettner related quantities
    // ------------------------------------------
	
    // Maxwell-Juettner cumulative function (array)
    std::vector<double> max_jutt_cumul;
    
    if (params->species_param[ispec].initialization_type=="maxwell-juettner") {
	    //! \todo{Pass this parameters in a code constants class (MG)}
	    nE     = 20000;
	    muEmax = 20.0;

	    max_jutt_cumul.resize(nE);
	    double mu=params->species_param[ispec].mass/params->species_param[ispec].temperature[0];
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
    unsigned int iPart=0;
    unsigned int *indexes=new unsigned int[ndim];
    double *temp=new double[ndim];
    double *vel=new double[ndim];

    // start a loop on all cells

    // Rappel :
    // int cell_index[0] = process_coord_x*(params->n_space[0]);
    // int cell_index[1] = process_coord_y*(params->n_space[1]);
    // int cell_index[2] = process_coord_z*(params->n_space[2]);

    for (unsigned int k=0; k<params->n_space[2]; k++) {
	    for (unsigned int j=0; j<params->n_space[1]; j++) {
		    for (unsigned int i=0; i<params->n_space[0]; i++) {
			    // initialize particles in meshes where the density is non-zero
				if (density(i,j,k)>0) {
				  //DEBUG(0,i);
					indexes[0]=i+cell_index[0];
					temp[0]=temperature[0](i,j,k);
					vel[0]=velocity[0](i,j,k);
					if (ndim > 1) {
						indexes[1]=j+cell_index[1];
						temp[1]=temperature[1](i,j,k);
						vel[1]=velocity[1](i,j,k);
						if (ndim > 2) {
							indexes[2]=k+cell_index[2];
							temp[2]=temperature[2](i,j,k);
							vel[2]=velocity[2](i,j,k);
						}//ndim > 2
					}//ndim > 1
                    
					initPosition(params->species_param[ispec].n_part_per_cell,iPart, indexes, ndim, 
						     params->cell_length, params->species_param[ispec].initialization_type);
					initMomentum(params->species_param[ispec].n_part_per_cell,iPart, temp, vel, 
						     params->species_param[ispec].initialization_type, max_jutt_cumul);
					initWeight(params, ispec, iPart, density(i,j,k));
                    initCharge(params, ispec, iPart, density(i,j,k));
                    
					//calculate new iPart (jump to next cell)
					iPart+=params->species_param[ispec].n_part_per_cell;
                    
				}//END if density > 0
		    }//i
	    }//j
    }//k end the loop on all cells
    
    delete indexes;
    delete temp;
    delete vel;
   
    // assign the correct Pusher to Push
    Push = PusherFactory::create( params, ispec );
	  
	// assign the Ionization model (if needed) to Ionize
	Ionize = IonizationFactory::create( params, ispec );
    if (Ionize) DEBUG("------------------------------------------------------------ " <<ispec);
    
    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond( params, ispec, smpi);

    PMESSAGE( 1, smpi->getRank(),"Species "<< ispec <<" # part "<< npart_effective );
    		
    // define a template particle for particles sorting
    swapPart = ParticleFactory::create(params, ispec);

}//END Species creator


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
	int nParticles = particles.size();
	for (int iPart=0 ; iPart<nParticles; iPart++ ) {
		delete particles[iPart];
	}
	delete Push;
	if (Ionize) delete Ionize;
	if (partBoundCond) delete partBoundCond;
	DEBUG(10,"Species deleted ");
	delete swapPart;
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its numerical weight (equivalent to a number density)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initWeight(PicParams* params, unsigned int ispec, unsigned int iPart, double density)
{
	for (unsigned  p= iPart; p<iPart+params->species_param[ispec].n_part_per_cell; p++) {
        particles[p]->weight() = density / params->species_param[ispec].n_part_per_cell;
/*
		particles[p]->weight() = density * params->species_param[ispec].charge
       /                        params->species_param[ispec].n_part_per_cell;
*/
	}
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize its charge state
// ---------------------------------------------------------------------------------------------------------------------
void Species::initCharge(PicParams* params, unsigned int ispec, unsigned int iPart, double density)
{
	for (unsigned  p= iPart; p<iPart+params->species_param[ispec].n_part_per_cell; p++) {
		particles[p]->charge() = params->species_param[ispec].charge;
	}
}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their position
//   - either using regular distribution in the mesh (regular)
//   - or using uniform random distribution (for cold and maxwell-juettner distribution)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initPosition(unsigned int np, unsigned int iPart, unsigned int *indexes, unsigned int ndim,
                           std::vector<double> cell_length, string initialization_type)
{
	for (unsigned  p= iPart; p<iPart+np; p++)
    {
		for (unsigned  i=0; i<ndim ; i++)
        {
			if (initialization_type == "regular") {
				particles[p]->position(i)=indexes[i]*cell_length[i]+(p-iPart)*cell_length[i]/np;
			} else if (initialization_type == "cold" || initialization_type == "maxwell-juettner") {
				particles[p]->position(i)=(indexes[i]+((double)rand() / RAND_MAX))*cell_length[i];
			}
			particles[p]->position_old(i) = particles[p]->position(i);
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
	double pMean[3]={0.0,0.0,0.0};
    
    
	if (initialization_type == "regular" || initialization_type == "cold")
    {
        // initialize momentum at 0 for regular or cold initialization type
		for (unsigned int p= iPart; p<iPart+np; p++) {
			for (unsigned int i=0; i<3 ; i++) {
				particles[p]->momentum(i) = 0.0;
			}
		}
        
	} else if (initialization_type == "maxwell-juettner")
    {
        // initialize using the Maxwell-Juettner distribution function
        //! \todo{Generalize to non-isotrop temperature (MG)}
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
			
			particles[p]->momentum(0) = psm*cos(theta)*sin(phi);
			particles[p]->momentum(1) = psm*sin(theta)*sin(phi);
			particles[p]->momentum(2) = psm*cos(phi);
			for (unsigned int i=0; i<3 ; i++)
            {
				pMean[i] += particles[p]->momentum(i);
			}
		}//p
        
        // center the distribution function around pMean
        // \todo{Allow for non-zero mean-velocity (MG)}
		for (unsigned int p= iPart; p<iPart+np; p++)
        {
			for (unsigned int i=0; i<3 ; i++) {
				particles[p]->momentum(i) -= pMean[i]/np;
			}
		}
        
	}//END if initialization_type

}//END initMomentum


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamic(double time_dual, ElectroMagn* Champs, Interpolator* Interp, Projector* Proj, SmileiMPI* smpi)
{
 	//! \todo{benefit from locate to create list of particles to send (JD)}
  	// SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

	// electric field at the particle position
	LocalFields Epart;
	// magnetic field at the particle position
	LocalFields Bpart;
	
	// number of particles for this Species
	int unsigned nParticles = getNbrOfParticles();
	smpi->clearExchList();

	// -------------------------------
	// calculate the particle dynamics
	// -------------------------------
	if (time_dual>time_frozen) {
		// moving particle
		//int locate;
		double gf = 1.0;

		// for all particles of the Species
		//#pragma omp parallel for shared (Champs)
		for (unsigned int iPart=0 ; iPart<nParticles; iPart++ ) {
			// Interpolate the fields at the particle position
			(*Interp)(Champs, particles[iPart], &Epart, &Bpart);
			
			// Do the ionization
			if (Ionize && particles[iPart]->charge() < (int) atomic_number) { //AND
				(*Ionize)(particles[iPart], Epart);
			}

			// Push the particle
			//(*Push)(particles[iPart], Epart, Bpart, gf);

			// Apply boundary condition on the particles
			// Boundary Condition may be physical or due to domain decomposition
			// apply returns 0 if iPart is no more in the domain local
			//	if omp then critical on smpi->addPartInExchList, may be applied after // loop
			//partBoundCond->apply( particles[iPart] );
			if ( !partBoundCond->apply( particles[iPart] ) ) smpi->addPartInExchList( iPart );

			(*Proj)(Champs, particles[iPart], gf);

		}// iPart
//		for (unsigned int iPart=0 ; iPart<nParticles; iPart++ ) {
//			if ( !partBoundCond->apply( particles[iPart] ) ) smpi->addPartInExchList( iPart );
//		}
        if (Ionize && electron_species) {
//            if (Ionize->new_electrons.size()) DEBUG("passing ionized electrons " << Ionize->new_electrons.size());
            for (unsigned int i=0; i < Ionize->new_electrons.size(); i++) {
                electron_species->particles.push_back(Ionize->new_electrons[i]);
            }
//			cerr << "****************** " << speciesNumber << " " << Ionize->new_electrons.size() << " " << electron_species->particles.size() << endl;
            Ionize->new_electrons.clear();
        }
	}
	else {
		// immobile particle (at the moment only project density)
                
		//! \todo{Implement Esirkepov method for the longitudinal currents (MG)}
		//#pragma omp parallel for shared (Champs)
		for (unsigned int iPart=0 ; iPart<nParticles; iPart++ ) {
			(*Proj)(Champs->rho_, particles[iPart]);
		}
	}//END if time vs. time_frozen
	

}//END dynamic



// ---------------------------------------------------------------------------------------------------------------------
// Dump for the Species
//! \todo{Check if one wants to keep it like this (MG)}
// ---------------------------------------------------------------------------------------------------------------------
void Species::dump(std::ofstream& ofile)
{
	for (unsigned int i=0; i<npart_effective; i++ )
    {
		ofile << i ;
		for (unsigned int m=0; m<ndim; m++) ofile << "\t" << particles[i]->position(m);
		for (unsigned int m=0; m<3; m++)    ofile << "\t" << particles[i]->momentum(m);
		ofile << "\t" << particles[i]->weight(); //<< "\t" << Push->getMass() << "\t" << Push->getCharge();
		ofile << endl;
	}
	ofile << endl;
}
// It computes the method on the singole specie. You can add here your parameter for a new diagnostic.
void Species::computeScalar(){
	double charge_tot=0.0;
	if (getNbrOfParticles()>0) {
		for (unsigned int iPart=0 ; iPart<getNbrOfParticles(); iPart++ ) {
			charge_tot+=(double)particles[iPart]->charge();
		}
	}
	scalars["charge_tot"]=charge_tot;
	scalars["part_number"]=getNbrOfParticles();
	
//	scalars["toto"]=max(0.0,charge_tot);
}

// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::sort_part()
{
	swap_part(particles[0], particles[1]);
}

// ---------------------------------------------------------------------------------------------------------------------
// Exchange particles part1 & part2 memory location
// ---------------------------------------------------------------------------------------------------------------------
void Species::swap_part(Particle* part1, Particle* part2)
{
	// hard compute of part_mem_size must be replaced
	int part_mem_size=(2*ndim+3+1)*sizeof(double)+sizeof(short);
	memcpy( &(swapPart->buf), &(part1->buf), part_mem_size);
	memcpy( &(part1->buf), &(part2->buf), part_mem_size);
	memcpy( &(part2->buf), &(swapPart->buf), part_mem_size);

}
