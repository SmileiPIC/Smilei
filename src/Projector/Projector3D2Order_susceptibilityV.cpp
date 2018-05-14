#include "Projector3D2Order_susceptibilityV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector3D2Order_susceptibilityV
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2Order_susceptibilityV::Projector3D2Order_susceptibilityV (Params& params, Patch* patch) : Projector3D(params, patch)
{
    dx_inv_   = 1.0/params.cell_length[0];
    dy_inv_   = 1.0/params.cell_length[1];
    dz_inv_   = 1.0/params.cell_length[2];
  
    one_third = 1.0/3.0;

    i_domain_begin = patch->getCellStartingGlobalIndex(0);
    j_domain_begin = patch->getCellStartingGlobalIndex(1);
    k_domain_begin = patch->getCellStartingGlobalIndex(2);

    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;

    DEBUG("cell_length "<< params.cell_length[0]);

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector3D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2Order_susceptibilityV::~Projector3D2Order_susceptibilityV()
{
}


// ---------------------------------------------------------------------------------------------------------------------
//! Project local currents (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibilityV::operator() (double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
} // END Project local current densities (Jx, Jy, Jz, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local current densities (sort)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibilityV::operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold)
{
} // END Project local densities (Jx, Jy, Jz, rho, sort)


// ---------------------------------------------------------------------------------------------------------------------
//! Project local densities only (Frozen species)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibilityV::operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim)
{
} // END Project local current densities (Frozen species)

// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities (ionize)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2Order_susceptibilityV::operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion)
{
} // END Project global current densities (ionize)

//Wrapper for projection
void Projector3D2Order_susceptibilityV::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref)
{
    // std::vector<int> *iold = &(smpi->dynamics_iold[ithread]);
    // std::vector<double> *delta = &(smpi->dynamics_deltaold[ithread]);
    // std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);
    // 
    // int dim1 = EMfields->dimPrim[1];
    // int dim2 = EMfields->dimPrim[2];
    // 
    // // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    // if (!diag_flag){ 
    //     if (!is_spectral) {
    //         double* b_Jx =  &(*EMfields->Jx_ )(ibin*clrw* dim1   * dim2   );
    //         double* b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)* dim2   );
    //         double* b_Jz =  &(*EMfields->Jz_ )(ibin*clrw* dim1   *(dim2+1));
    //         for ( int ipart=istart ; ipart<iend; ipart++ )
    //             (*this)(b_Jx , b_Jy , b_Jz , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
    //     }
    //     else {
    //         double* b_Jx =  &(*EMfields->Jx_ )(ibin*clrw* dim1   * dim2   );
    //         double* b_Jy =  &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)* dim2   );
    //         double* b_Jz =  &(*EMfields->Jz_ )(ibin*clrw* dim1   *(dim2+1));
    //         double* b_rho=  &(*EMfields->rho_)(ibin*clrw* dim1   * dim2   );
    //         for ( int ipart=istart ; ipart<iend; ipart++ )
    //             (*this)(b_Jx , b_Jy , b_Jz , b_rho , particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
    //     }
    // // Otherwise, the projection may apply to the species-specific arrays
    // } else {
    //     double* b_Jx  = EMfields->Jx_s [ispec] ? &(*EMfields->Jx_s [ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->Jx_ )(ibin*clrw* dim1   *dim2) ;
    //     double* b_Jy  = EMfields->Jy_s [ispec] ? &(*EMfields->Jy_s [ispec])(ibin*clrw*(dim1+1)*dim2) : &(*EMfields->Jy_ )(ibin*clrw*(dim1+1)*dim2) ;
    //     double* b_Jz  = EMfields->Jz_s [ispec] ? &(*EMfields->Jz_s [ispec])(ibin*clrw*dim1*(dim2+1)) : &(*EMfields->Jz_ )(ibin*clrw*dim1*(dim2+1)) ;
    //     double* b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(ibin*clrw* dim1   *dim2) : &(*EMfields->rho_)(ibin*clrw* dim1   *dim2) ;
    //     for ( int ipart=istart ; ipart<iend; ipart++ )
    //         (*this)(b_Jx , b_Jy , b_Jz ,b_rho, particles,  ipart, (*invgf)[ipart], ibin*clrw, b_dim, &(*iold)[ipart], &(*delta)[ipart]);
    // }

}


// Wrapper for projector of susceptibility




// Projector for susceptibility used as source term in envelope equation
                                                     
void Projector3D2Order_susceptibilityV::project_susceptibility(double* Chi_envelope, Particles &particles, int istart, int iend, unsigned int bin, std::vector<unsigned int> &b_dim, SmileiMPI* smpi, int ithread, double species_mass, int* iold, int ipart_ref)
{
    std::vector<double> *Epart       = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Phipart     = &(smpi->dynamics_PHIpart[ithread]);
    std::vector<double> *GradPhipart = &(smpi->dynamics_GradPHIpart[ithread]);

    int nparts = smpi->dynamics_invgf[ithread].size();
    double* Ex       = &( (*Epart)[0*nparts] );
    double* Ey       = &( (*Epart)[1*nparts] );
    double* Ez       = &( (*Epart)[2*nparts] );
    double* Phi      = &( (*Phipart)[0*nparts] );
    double* GradPhix = &( (*GradPhipart)[0*nparts] );
    double* GradPhiy = &( (*GradPhipart)[1*nparts] );
    double* GradPhiz = &( (*GradPhipart)[2*nparts] );

    int vecSize = 8;
    int bsize = 5*5*5*vecSize; // primal grid, particles did not yet move (3x3x3 enough)
    double bChi[bsize] __attribute__((aligned(64)));
    
    double Sx1[40] __attribute__((aligned(64)));
    double Sy1[40] __attribute__((aligned(64)));
    double Sz1[40] __attribute__((aligned(64)));
    double charge_weight[8] __attribute__((aligned(64)));

    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( (int)iend-(int)istart );
    int nbVec = ( iend-istart+(cell_nparts-1)-((iend-istart-1)&(cell_nparts-1)) ) / vecSize;
    if (nbVec*vecSize != cell_nparts)
        nbVec++;

    #pragma omp simd
    for (unsigned int j=0; j<bsize; j++)
        bChi[j] = 0.;
    
    for (int ivect=0 ; ivect < cell_nparts; ivect += vecSize ){

        int np_computed(min(cell_nparts-ivect,vecSize));
        int istart0 = (int)istart + ivect;
        
        //for (unsigned int ipart=istart ; (int)ipart<iend; ipart++ ) {
        //#pragma omp simd
        for (int ipart=0 ; ipart<np_computed; ipart++ ) {
            
            int iloc,jloc;
     
            double momentum[3];
        
            double inv_gamma_ponderomotive,inv_gamma0;
            double charge_over_mass_dts2,charge_sq_over_mass_dts4,charge_sq_over_mass_sq;
            double pxsm, pysm, pzsm;
            double one_over_mass=1./species_mass;

            charge_over_mass_dts2    = particles.charge(istart+ipart)*dts2*one_over_mass;
            // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
            charge_sq_over_mass_dts4 = particles.charge(istart+ipart)*dts4*one_over_mass;      
            // (charge over mass)^2
            charge_sq_over_mass_sq   = particles.charge(istart+ipart)*particles.charge(istart+ipart)*one_over_mass*one_over_mass;

            for ( int i = 0 ; i<3 ; i++ )
                momentum[i] = particles.momentum(i,istart+ipart);
 
            // compute initial ponderomotive gamma (more precisely, its inverse) 
            inv_gamma0 = 1./sqrt( 1. + momentum[0]*momentum[0]+ momentum[1]*momentum[1] + momentum[2]*momentum[2] + *(Phi+istart-ipart_ref+ipart)*charge_sq_over_mass_sq );
         
            // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
            pxsm = inv_gamma0 * (charge_over_mass_dts2*(*(Ex+istart-ipart_ref+ipart)) - charge_sq_over_mass_dts4*(*(GradPhix+istart-ipart_ref+ipart)) * inv_gamma0 ) * momentum[0];
            pysm = inv_gamma0 * (charge_over_mass_dts2*(*(Ey+istart-ipart_ref+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiy+istart-ipart_ref+ipart)) * inv_gamma0 ) * momentum[1];
            pzsm = inv_gamma0 * (charge_over_mass_dts2*(*(Ez+istart-ipart_ref+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiz+istart-ipart_ref+ipart)) * inv_gamma0 ) * momentum[2];
         
            // update of gamma ponderomotive (more precisely, the inverse)
            inv_gamma_ponderomotive = 1./( 1./inv_gamma0 + (pxsm+pysm+pzsm)*0.5 );
 
            // (x,y,z) components of the current density for the macro-particle
            charge_weight[ipart] = (double)(particles.charge(istart+ipart))*(double)(particles.charge(istart+ipart))*particles.weight(istart+ipart)*inv_gamma_ponderomotive*one_over_mass; 
 
            // variable declaration
            double xpn, ypn, zpn;
            double delta, delta2;
 
            // Initialize all current-related arrays to zero
            for (unsigned int i=0; i<5; i++) {
                Sx1[i*vecSize+ipart] = 0.;
                Sy1[i*vecSize+ipart] = 0.;
                Sz1[i*vecSize+ipart] = 0.;
            }

            // --------------------------------------------------------
            // Locate particles & Calculate Esirkepov coef. S, DS and W
            // --------------------------------------------------------
 
            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            xpn = particles.position(0, istart+ipart) * dx_inv_;
            int ip = round(xpn);
            delta  = xpn - (double)ip;
            delta2 = delta*delta;
            Sx1[1*vecSize+ipart] = 0.5 * (delta2-delta+0.25);
            Sx1[2*vecSize+ipart] = 0.75-delta2;
            Sx1[3*vecSize+ipart] = 0.5 * (delta2+delta+0.25);
 
            ypn = particles.position(1, istart+ipart) * dy_inv_;
            int jp = round(ypn);
            delta  = ypn - (double)jp;
            delta2 = delta*delta;
            Sy1[1*vecSize+ipart] = 0.5 * (delta2-delta+0.25);
            Sy1[2*vecSize+ipart] = 0.75-delta2;
            Sy1[3*vecSize+ipart] = 0.5 * (delta2+delta+0.25);
 
            zpn = particles.position(2, istart+ipart) * dz_inv_;
            int kp = round(zpn);
            delta  = zpn - (double)kp;
            delta2 = delta*delta;
            Sz1[1*vecSize+ipart] = 0.5 * (delta2-delta+0.25);
            Sz1[2*vecSize+ipart] = 0.75-delta2;
            Sz1[3*vecSize+ipart] = 0.5 * (delta2+delta+0.25);
 
        } // end ipart loop

        for (int ipart=0 ; ipart<np_computed; ipart++ ) {
            for (unsigned int i=0 ; i<5 ; i++) {
                for (unsigned int j=0 ; j<5 ; j++) {
                    int index( ( i*25 + j*5 )*vecSize+ipart );
                    for (unsigned int k=0 ; k<5 ; k++) {
                        bChi [ index+k*vecSize ] +=  charge_weight[ipart] * Sx1[i*vecSize+ipart]*Sy1[j*vecSize+ipart]*Sz1[k*vecSize+ipart];
                    }
                }
            }//i

        } // end ipart loop
        
    } // end ivect

    // ---------------------------
    // Calculate the total charge
    // ---------------------------    
    int ipo = iold[0];
    int jpo = iold[1];
    int kpo = iold[2];
    int ipom2 = ipo-2;
    int jpom2 = jpo-2;
    int kpom2 = kpo-2;

    int iloc0 = ipom2*b_dim[1]*b_dim[2]+jpom2*b_dim[2]+kpom2;
    int iloc = iloc0;
    for (unsigned int i=0 ; i<5 ; i++) {
        for (unsigned int j=0 ; j<5 ; j++) {
            #pragma omp simd
            for (unsigned int k=0 ; k<5 ; k++) {
                double tmpChi = 0.;
                int ilocal = ((i)*25+j*5+k)*vecSize;
                #pragma unroll(8)
                for (int ipart=0 ; ipart<8; ipart++ ){
                    tmpChi +=  bChi[ilocal+ipart];
                }
                 Chi_envelope [iloc + (j)*(b_dim[2]) + k] +=  tmpChi;
            }
        }
        iloc += b_dim[1]*(b_dim[2]);
    }

    
    
}
