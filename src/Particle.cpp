
#include "Particle.h"

#include "PicParams.h"

#include <iostream>

using namespace std;

Particle::Particle(int nDim)
{
	Pos_ = new double[nDim];
//	DEBUG("Particle created "<<nDim<<"D");
}

Particle::~Particle()
{
	delete [] Pos_;
}

//void Particle::Initialize(PicParams* params, int ispec, int np, double x0_pos, double x1_pos, double x2_pos, 
//						  double density, double temperature, double velocity_x, double velocity_y, double velocity_z)
//{	
//	weight_ = density / params->species_param[ispec].n_part_per_cell; 
//
//	if (params->species_param[ispec].initialization_type=="regular") {
//		Psm_[0]=0.0;
//		Psm_[1]=0.0;
//		Psm_[2]=0.0;
//	}
//	
//		
//	/*do j=ndeb_s(is)+(ix-1)*nppm_s(is),ndeb_s(is)+ix*nppm_s(is)-1
//     ! RANDOMLY DISTRIBUTED PARTICLES
//     particle(posx,j)     = (real(nv+ix-1,rprec)+ran2(seed_ran2))*dx
//     particle(pxsm,j)     = 0.d0
//     particle(pysm,j)     = 0.d0
//     particle(pzsm,j)     = 0.d0
//     particle(weight,j)   = dens_loc/nppm_s(is)
//     particle(charge,j)   = q_s(is)
//     particle(mass,j)     = m_s(is)
//	 enddo!j*/
//	
//}

void Particle::Print(PicParams* params)
{
	for (unsigned int i=0 ; i<params->nDim_field ; i++ ) cout << Pos_[i] << " ";
	for (unsigned int i=0 ; i<3 ; i++ )                  cout << Psm_[i] << " ";
	cout <<  charge_density_ << " " << endl;
}

