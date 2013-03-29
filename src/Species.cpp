#include "Tools.h"
#include "Species.h"
#include "Particle.h"
#include "ParticleRad.h"

#include "Interpolator.h"
#include "Interpolator1D2Order.h"
#include "Projector.h"
#include "Pusher.h"
#include "PusherBoris.h"

#include "Field3D.h"

#include "ElectroMagn1D.h"

#include "PartBoundCond.h"
#include "BoundaryConditionType.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

Species::Species(PicParams* params, int ispec) {
	ndim=params->nDim_particle;
	
	time_frozen = params->species_param[ispec].time_frozen;
	Field3D density(params->n_space);
	Field3D temperature[3];
	Field3D velocity[3];
	for (size_t i=0;i<3;i++) {
		velocity[i].allocateDims(params->n_space);
		temperature[i].allocateDims(params->n_space);
	}
	
	// calculate density and number of particles for species ispec
	npart_effective=0;
	for (unsigned int k=0; k<params->n_space[2]; k++) {
		for (unsigned int j=0; j<params->n_space[1]; j++) {
			for (unsigned int i=0; i<params->n_space[0]; i++) {
			//	int icell = i + j*params->n_space[0] +k*params->n_space[0]*params->n_space[1];
				if (params->plasma_geometry=="constant") {
		//			DEBUG(i << " " << j << " " << k << " ");
					if (((params->cell_length[0]==0.0) || (
						 (i+0.5)*params->cell_length[0] > params->vacuum_length[0] && 
						 (i+0.5)*params->cell_length[0] < params->vacuum_length[0]+params->plasma_length[0]
						 )) &&
						((params->cell_length[1]==0.0) || (
						 (j+0.5)*params->cell_length[1] > params->vacuum_length[1] && 
						 (j+0.5)*params->cell_length[1] < params->vacuum_length[1]+params->plasma_length[1]
						 )) &&
						((params->cell_length[2]==0.0) || (
						 (k+0.5)*params->cell_length[2] > params->vacuum_length[2] && 
						 (k+0.5)*params->cell_length[2] < params->vacuum_length[2]+params->plasma_length[2])
						)) {

						density(i,j,k) = params->species_param[ispec].density;
						for (size_t m=0; m<3;m++) {
							temperature[m](i,j,k) = params->species_param[ispec].temperature[m];
							velocity[m](i,j,k) = params->species_param[ispec].mean_velocity[m];
						}					
						npart_effective += params->species_param[ispec].n_part_per_cell;
						DEBUG(10,"Specie "<< ispec <<" # part "<<npart_effective<<" "<<i<<" "<<j<<" "<<k);

					} else {
						density(i,j,k)=0.0;
					}
				} else {
					ERROR("geometry not implemented");
				}
			}
		}
	}


    params->species_param[ispec].n_part_max = round( params->species_param[ispec].c_part_max*npart_effective );
	
/*	if (npart_effective > params->species_param[ispec].n_part_max) {
		WARNING( "Number of particles for species " << ispec << " redefined to " << npart_effective );
		params->species_param[ispec].n_part_max=npart_effective;
	}
*/
    
	particles.resize(params->species_param[ispec].n_part_max);
	for (unsigned int k=0; k<npart_effective; k++) {
		if (params->species_param[ispec].radiating) {
			particles[k]=new ParticleRad(ndim);
		} else {
			particles[k]=new Particle(ndim);
		}
	}
	
	
	//! Maxwell-Juettner array
	std::vector<double> max_jutt_cumul;
	if (params->species_param[ispec].initialization_type=="maxwell-juettner") {
		//! todo{control this}
		nE=20000;
		muEmax=20.0;

		max_jutt_cumul.resize(nE);
		double mu=params->species_param[ispec].mass/params->species_param[ispec].temperature[0];
		double Emax=muEmax/mu;
		dE=Emax/nE;
	
		double fl=0;
		double fr=0;
		max_jutt_cumul[0]=0.0;
		for (unsigned  i=1; i<nE; i++ ) {
			//! \todo{this is just the isotropic case}
			fr=(1+i*dE)*sqrt(pow(1.0+i*dE,2)-1.0) * exp(-mu*i*dE);
			max_jutt_cumul[i]=max_jutt_cumul[i-1] + 0.5*dE*(fr+fl);
			fl=fr;
		}
		for (size_t i=0; i<nE; i++) max_jutt_cumul[i]/=max_jutt_cumul[nE-1];
		
	}
	
	size_t iPart=0;
	size_t *indexes=new size_t[ndim];
	double *temp=new double[ndim];
	double *vel=new double[ndim];
	
	for (unsigned int k=0; k<params->n_space[2]; k++) {
		for (unsigned int j=0; j<params->n_space[1]; j++) {
			for (unsigned int i=0; i<params->n_space[0]; i++) {
				if (density(i,j,k)>0) {
					DEBUG(0,i);
					indexes[0]=i;
					temp[0]=temperature[0](i,j,k);
					vel[0]=velocity[0](i,j,k);
					if (ndim > 1) {
						indexes[1]=j;
						temp[1]=temperature[1](i,j,k);
						vel[1]=velocity[1](i,j,k);
						if (ndim > 2) {
							indexes[2]=k;
							temp[2]=temperature[2](i,j,k);
							vel[2]=velocity[2](i,j,k);
						}
					}
					initPosition(params->species_param[ispec].n_part_per_cell,iPart, indexes, ndim, 
								 params->cell_length, params->species_param[ispec].initialization_type);
					initMomentum(params->species_param[ispec].n_part_per_cell,iPart, temp, vel, 
								 params->species_param[ispec].initialization_type, max_jutt_cumul);
					initMassChargeWeight(params, ispec, iPart, density(i,j,k));
					iPart+=params->species_param[ispec].n_part_per_cell;
				}
			}
		}
	}
	delete indexes;
	delete temp;
	delete vel;

	if ( params->species_param[ispec].dynamics_type == "norm" )
		Push = new PusherBoris( params, ispec );
	else
		ERROR( "Unknwon dynamics : " << params->species_param[ispec].dynamics_type );
	
	//! \todo{other dimensions to store in class members, n_ord_proj_max to define as input}
	partBoundCond = new PartBoundCond(  params, ispec );

	MESSAGE(1,"Specie "<< ispec <<" # part "<< npart_effective);
}

Species::~Species()
{
	int nParticles = particles.size();
	for (int iPart=0 ; iPart<nParticles; iPart++ ) {
		delete particles[iPart];
	}
	DEBUG(10,"Species deleted ");

	delete Push;
}

void Species::initMassChargeWeight(PicParams* params, int ispec, int iPart, double density){
	for (unsigned  p= iPart; p<iPart+params->species_param[ispec].n_part_per_cell; p++) {
//		particles[p]->weight_ = density / (params->species_param[ispec].n_part_per_cell * params->cell_volume);
		particles[p]->chargeDensity() = density / params->species_param[ispec].n_part_per_cell * params->species_param[ispec].charge;
	}
}

void Species::initPosition(unsigned int np, unsigned int iPart, size_t *indexes, unsigned int ndim, std::vector<double> cell_length, string initialization_type) {
	for (unsigned  p= iPart; p<iPart+np; p++) {
		for (unsigned  i=0; i<ndim ; i++) {
			if (initialization_type == "regular") {
				particles[p]->position(0)=indexes[i]*cell_length[i]+(p-iPart)*cell_length[i]/np;
			} else if (initialization_type == "cold" || initialization_type == "maxwell-juettner") {
				particles[p]->position(0)=(indexes[i]+((double)rand() / RAND_MAX))*cell_length[i];
			}
		}
	}
}


void Species::initMomentum(unsigned int np, unsigned int iPart, double *temp, double *vel, string initialization_type, vector<double>& max_jutt_cumul) {
	double pMean[3]={0.0,0.0,0.0};
	if (initialization_type == "regular" || initialization_type == "cold") {
		for (size_t p= iPart; p<iPart+np; p++) {
			for (size_t i=0; i<3 ; i++) {
				particles[p]->moments(i) = 0.0;
			}
		}
	} else if (initialization_type == "maxwell-juettner") {
		for (size_t p= iPart; p<iPart+np; p++) {
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
			
			particles[p]->moments(0) = psm*cos(theta)*sin(phi);
			particles[p]->moments(1) = psm*sin(theta)*sin(phi);
			particles[p]->moments(2) = psm*cos(phi);
			for (size_t i=0; i<3 ; i++) {
				pMean[i] += particles[p]->moments(i);
			}
		}
		for (size_t p= iPart; p<iPart+np; p++) {
			for (size_t i=0; i<3 ; i++) {
				particles[p]->moments(i) -= pMean[i]/np;
			}
		}
	}

}

void Species::dynamic(double time_dual, ElectroMagn* Champs, Interpolator* Interp, Projector* Proj)
{
	chLocaux Epart;
	chLocaux Bpart;
	
	int nParticles = getNbrOfParticles();

	if (time_dual>time_frozen) {
		int locate;
		double gf = 1.0;
		for (int iPart=0 ; iPart<nParticles; iPart++ ) {
			DEBUG(5,"ipart= "<<iPart);
			(*Interp)(Champs, particles[iPart], &Epart, &Bpart);
			gf = 1.0;
			(*Push)(particles[iPart], Epart, Bpart, gf);

			//! \todo{locate : define 0 to ...}
			locate = partBoundCond->locateParticle( particles[iPart] );
			if ( locate == 0 )	{}
			else if ( locate == 1 )	(*partBoundCond->bc_west)( particles[iPart], 2.*partBoundCond->x_min );
			else if ( locate == 2 )	(*partBoundCond->bc_east)( particles[iPart], 2.*partBoundCond->x_max );
			// else ...

			(*Proj)(Champs, particles[iPart], gf);
		}
	}
	else {
		for (int iPart=0 ; iPart<nParticles; iPart++ ) {
			(*Proj)(Champs->rho_, particles[iPart]);
		}
	}
	
}

void Species::dump(std::ofstream& ofile) {
	for (size_t i=0; i<npart_effective; i++ ) {
		ofile << i ;
		for (size_t m=0; m<ndim; m++) ofile << "\t" << particles[i]->position(m);
		for (size_t m=0; m<3; m++)    ofile << "\t" << particles[i]->moments(m);
		ofile << "\t" << particles[i]->chargeDensity() << "\t" << Push->getMass() << "\t" << Push->getCharge();
		ofile << endl;
	}
	ofile << endl;
}


