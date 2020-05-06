#include "Projector3D2OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector3D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2OrderV::Projector3D2OrderV( Params &params, Patch *patch ) : Projector3D( params, patch )
{
    dx_inv_   = 1.0/params.cell_length[0];
    dx_ov_dt  = params.cell_length[0] / params.timestep;
    dy_inv_   = 1.0/params.cell_length[1];
    dy_ov_dt  = params.cell_length[1] / params.timestep;
    dz_inv_   = 1.0/params.cell_length[2];
    dz_ov_dt  = params.cell_length[2] / params.timestep;
    
    i_domain_begin = patch->getCellStartingGlobalIndex( 0 );
    j_domain_begin = patch->getCellStartingGlobalIndex( 1 );
    k_domain_begin = patch->getCellStartingGlobalIndex( 2 );
    
    nscelly = params.n_space[1] + 1;
    nscellz = params.n_space[2] + 1;
    oversize[0] = params.oversize[0];
    oversize[1] = params.oversize[1];
    oversize[2] = params.oversize[2];
    nprimy = nscelly + 2*oversize[1];
    nprimz = nscellz + 2*oversize[2];
    dq_inv[0] = dx_inv_;
    dq_inv[1] = dy_inv_;
    dq_inv[2] = dz_inv_;
    
    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;
    
    DEBUG( "cell_length "<< params.cell_length[0] );
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Projector3D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Projector3D2OrderV::~Projector3D2OrderV()
{
}

// ---------------------------------------------------------------------------------------------------------------------
//!  Project current densities & charge : diagFields timstep (not vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, int *iold, double *deltaold, int ipart_ref )
{

    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int npart_total = invgf->size();
    int ipo = iold[0];
    int jpo = iold[1];
    int kpo = iold[2];
    int ipom2 = ipo-2;
    int jpom2 = jpo-2;
    int kpom2 = kpo-2;
    
    int vecSize = 8;
    unsigned int bsize = 5*5*5*vecSize;
    
    double bJx[bsize] __attribute__( ( aligned( 64 ) ) );
    
    double DSx[40] __attribute__( ( aligned( 64 ) ) );
    double DSy[40] __attribute__( ( aligned( 64 ) ) );
    double DSz[40] __attribute__( ( aligned( 64 ) ) );
    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );
    
    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );
    int nbVec = ( iend-istart+( cell_nparts-1 )-( ( iend-istart-1 )&( cell_nparts-1 ) ) ) / vecSize;
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }
    
    
    // Jx, Jy, Jz
    currents( Jx, Jy, Jz, particles, istart, iend, invgf, iold, deltaold, ipart_ref );
    
    
    // rho^(p,p,d)
    cell_nparts = ( int )iend-( int )istart;
    #pragma omp simd
    for( unsigned int j=0; j<1000; j++ ) {
        bJx[j] = 0.;
    }
    
    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {
    
        int np_computed( min( cell_nparts-ivect, vecSize ) );
        int istart0 = ( int )istart + ivect;
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            compute_distances( particles, npart_total, ipart, istart0, ipart_ref, deltaold, iold, DSx, DSy, DSz );
            charge_weight[ipart] = inv_cell_volume * ( double )( particles.charge( istart0+ipart ) )*particles.weight( istart0+ipart );
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            for( unsigned int i=0 ; i<5 ; i++ ) {
                for( unsigned int j=0 ; j<5 ; j++ ) {
                    int index( ( i*25 + j*5 )*vecSize+ipart );
                    for( unsigned int k=0 ; k<5 ; k++ ) {
                        bJx [ index+k*vecSize ] +=  charge_weight[ipart] * DSx[i*vecSize+ipart]*DSy[j*vecSize+ipart]*DSz[k*vecSize+ipart];
                    }
                }
            }//i
            
            
        } // END ipart (compute coeffs)
        
    }
    
    int iloc0 = ipom2*nprimy*nprimz+jpom2*nprimz+kpom2;
    int iloc = iloc0;
    for( unsigned int i=0 ; i<5 ; i++ ) {
        for( unsigned int j=0 ; j<5 ; j++ ) {
            #pragma omp simd
            for( unsigned int k=0 ; k<5 ; k++ ) {
                double tmpRho = 0.;
                int ilocal = ( ( i )*25+j*5+k )*vecSize;
#pragma unroll(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpRho +=  bJx[ilocal+ipart];
                }
                rho [iloc + ( j )*( nprimz ) + k] +=  tmpRho;
            }
        }
        iloc += nprimy*( nprimz );
    }
    
} // END Project local current densities at diagFields timestep.


// ---------------------------------------------------------------------------------------------------------------------
//! Project charge : frozen & diagFields timstep (not vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type )
{
    //Warning : this function is used for frozen species or initialization only and doesn't use the standard scheme.
    //rho type = 0
    //Jx type = 1
    //Jy type = 2
    //Jz type = 3
    
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int iloc, jloc;
    int ny( nprimy ), nz( nprimz ), nyz;
    // (x,y,z) components of the current density for the macro-particle
    
    // variable declaration
    double xpn, ypn, zpn;
    double delta, delta2;
    double Sx1[5], Sy1[5], Sz1[5]; // arrays used for the Esirkepov projection method
    
    double charge_weight = inv_cell_volume * ( double )( particles.charge( ipart ) )*particles.weight( ipart );
    if( type > 0 ) {
        charge_weight *= 1./sqrt( 1.0 + particles.momentum( 0, ipart )*particles.momentum( 0, ipart )
                                  + particles.momentum( 1, ipart )*particles.momentum( 1, ipart )
                                  + particles.momentum( 2, ipart )*particles.momentum( 2, ipart ) );
                                  
        if( type == 1 ) {
            charge_weight *= particles.momentum( 0, ipart );
        } else if( type == 2 ) {
            charge_weight *= particles.momentum( 1, ipart );
            ny ++;
        } else {
            charge_weight *= particles.momentum( 2, ipart );
            nz ++;
        }
    }
    nyz = ny*nz;
    
    // Initialize all current-related arrays to zero
    for( unsigned int i=0; i<5; i++ ) {
        Sx1[i] = 0.;
        Sy1[i] = 0.;
        Sz1[i] = 0.;
    }
    
    // --------------------------------------------------------
    // Locate particles & Calculate Esirkepov coef. S, DS and W
    // --------------------------------------------------------
    
    // locate the particle on the primal grid at current time-step & calculate coeff. S1
    xpn = particles.position( 0, ipart ) * dx_inv_;
    int ip = round( xpn + 0.5*( type==1 ) );
    delta  = xpn - ( double )ip;
    delta2 = delta*delta;
    Sx1[1] = 0.5 * ( delta2-delta+0.25 );
    Sx1[2] = 0.75-delta2;
    Sx1[3] = 0.5 * ( delta2+delta+0.25 );
    
    ypn = particles.position( 1, ipart ) * dy_inv_;
    int jp = round( ypn + 0.5*( type==2 ) );
    delta  = ypn - ( double )jp;
    delta2 = delta*delta;
    Sy1[1] = 0.5 * ( delta2-delta+0.25 );
    Sy1[2] = 0.75-delta2;
    Sy1[3] = 0.5 * ( delta2+delta+0.25 );
    
    zpn = particles.position( 2, ipart ) * dz_inv_;
    int kp = round( zpn + 0.5*( type==3 ) );
    delta  = zpn - ( double )kp;
    delta2 = delta*delta;
    Sz1[1] = 0.5 * ( delta2-delta+0.25 );
    Sz1[2] = 0.75-delta2;
    Sz1[3] = 0.5 * ( delta2+delta+0.25 );
    
    // ---------------------------
    // Calculate the total charge
    // ---------------------------
    ip -= i_domain_begin + 2;
    jp -= j_domain_begin + 2;
    kp -= k_domain_begin + 2;
    
    for( unsigned int i=0 ; i<5 ; i++ ) {
        iloc = ( i+ip ) * nyz;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            jloc = ( jp+j ) * nz;
            for( unsigned int k=0 ; k<5 ; k++ ) {
                rhoj[iloc+jloc+kp+k] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
            }
        }
    }//i
    
} // END Project local current densities (Frozen species)


// ---------------------------------------------------------------------------------------------------------------------
//! Project global current densities : ionization (WARNING: Not Vectorized)
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion )
{
    Field3D *Jx3D  = static_cast<Field3D *>( Jx );
    Field3D *Jy3D  = static_cast<Field3D *>( Jy );
    Field3D *Jz3D  = static_cast<Field3D *>( Jz );
    
    
    //Declaration of local variables
    int ip, id, jp, jd, kp, kd;
    double xpn, xpmxip, xpmxip2, xpmxid, xpmxid2;
    double ypn, ypmyjp, ypmyjp2, ypmyjd, ypmyjd2;
    double zpn, zpmzkp, zpmzkp2, zpmzkd, zpmzkd2;
    double Sxp[3], Sxd[3], Syp[3], Syd[3], Szp[3], Szd[3];
    
    // weighted currents
    double weight = inv_cell_volume * particles.weight( ipart );
    double Jx_ion = Jion.x * weight;
    double Jy_ion = Jion.y * weight;
    double Jz_ion = Jion.z * weight;
    
    //Locate particle on the grid
    xpn    = particles.position( 0, ipart ) * dx_inv_; // normalized distance to the first node
    ypn    = particles.position( 1, ipart ) * dy_inv_; // normalized distance to the first node
    zpn    = particles.position( 1, ipart ) * dz_inv_; // normalized distance to the first node
    
    // x-primal index
    ip      = round( xpn );                  // x-index of the central node
    xpmxip  = xpn - ( double )ip;            // normalized distance to the nearest grid point
    xpmxip2 = xpmxip*xpmxip;                 // square of the normalized distance to the nearest grid point
    
    // x-dual index
    id      = round( xpn+0.5 );              // x-index of the central node
    xpmxid  = xpn - ( double )id + 0.5;      // normalized distance to the nearest grid point
    xpmxid2 = xpmxid*xpmxid;                 // square of the normalized distance to the nearest grid point
    
    // y-primal index
    jp      = round( ypn );                  // y-index of the central node
    ypmyjp  = ypn - ( double )jp;            // normalized distance to the nearest grid point
    ypmyjp2 = ypmyjp*ypmyjp;                 // square of the normalized distance to the nearest grid point
    
    // y-dual index
    jd      = round( ypn+0.5 );              // y-index of the central node
    ypmyjd  = ypn - ( double )jd + 0.5;      // normalized distance to the nearest grid point
    ypmyjd2 = ypmyjd*ypmyjd;                 // square of the normalized distance to the nearest grid point
    
    // z-primal index
    kp      = round( zpn );                  // z-index of the central node
    zpmzkp  = zpn - ( double )kp;            // normalized distance to the nearest grid point
    zpmzkp2 = zpmzkp*zpmzkp;                 // square of the normalized distance to the nearest grid point
    
    // z-dual index
    kd      = round( zpn+0.5 );              // z-index of the central node
    zpmzkd  = zpn - ( double )kd + 0.5;      // normalized distance to the nearest grid point
    zpmzkd2 = zpmzkd*zpmzkd;                 // square of the normalized distance to the nearest grid point
    
    Sxp[0] = 0.5 * ( xpmxip2-xpmxip+0.25 );
    Sxp[1] = ( 0.75-xpmxip2 );
    Sxp[2] = 0.5 * ( xpmxip2+xpmxip+0.25 );
    
    Sxd[0] = 0.5 * ( xpmxid2-xpmxid+0.25 );
    Sxd[1] = ( 0.75-xpmxid2 );
    Sxd[2] = 0.5 * ( xpmxid2+xpmxid+0.25 );
    
    Syp[0] = 0.5 * ( ypmyjp2-ypmyjp+0.25 );
    Syp[1] = ( 0.75-ypmyjp2 );
    Syp[2] = 0.5 * ( ypmyjp2+ypmyjp+0.25 );
    
    Syd[0] = 0.5 * ( ypmyjd2-ypmyjd+0.25 );
    Syd[1] = ( 0.75-ypmyjd2 );
    Syd[2] = 0.5 * ( ypmyjd2+ypmyjd+0.25 );
    
    Szp[0] = 0.5 * ( zpmzkp2-zpmzkp+0.25 );
    Szp[1] = ( 0.75-zpmzkp2 );
    Szp[2] = 0.5 * ( zpmzkp2+zpmzkp+0.25 );
    
    Szd[0] = 0.5 * ( zpmzkd2-zpmzkd+0.25 );
    Szd[1] = ( 0.75-zpmzkd2 );
    Szd[2] = 0.5 * ( zpmzkd2+zpmzkd+0.25 );
    
    ip  -= i_domain_begin;
    id  -= i_domain_begin;
    jp  -= j_domain_begin;
    jd  -= j_domain_begin;
    kp  -= k_domain_begin;
    kd  -= k_domain_begin;
    
    for( unsigned int i=0 ; i<3 ; i++ ) {
        int iploc=ip+i-1;
        int idloc=id+i-1;
        for( unsigned int j=0 ; j<3 ; j++ ) {
            int jploc=jp+j-1;
            int jdloc=jd+j-1;
            for( unsigned int k=0 ; k<3 ; k++ ) {
                int kploc=kp+k-1;
                int kdloc=kd+k-1;
                // Jx^(d,p,p)
                ( *Jx3D )( idloc, jploc, kploc ) += Jx_ion * Sxd[i]*Syp[j]*Szp[k];
                // Jy^(p,d,p)
                ( *Jy3D )( iploc, jdloc, kploc ) += Jy_ion * Sxp[i]*Syd[j]*Szp[k];
                // Jz^(p,p,d)
                ( *Jz3D )( iploc, jploc, kdloc ) += Jz_ion * Sxp[i]*Syp[j]*Szd[k];
            }//k
        }//j
    }//i
    
} // END Project global current densities (ionize)


// ---------------------------------------------------------------------------------------------------------------------
//! Project current densities : main projector vectorized
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::currents( double *Jx, double *Jy, double *Jz, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, int *iold, double *deltaold, int ipart_ref )
{
    // -------------------------------------
    // Variable declaration & initialization
    // -------------------------------------
    
    int npart_total = invgf->size();
    int ipo = iold[0];
    int jpo = iold[1];
    int kpo = iold[2];
    int ipom2 = ipo-2;
    int jpom2 = jpo-2;
    int kpom2 = kpo-2;
    int nyz = nprimy*nprimz;
    
    int vecSize = 8;
    unsigned int bsize = 5*5*5*vecSize;
    
    double bJx[bsize] __attribute__( ( aligned( 64 ) ) );
    
    double Sx0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    double Sy0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    double Sz0_buff_vect[32] __attribute__( ( aligned( 64 ) ) );
    double DSx[40] __attribute__( ( aligned( 64 ) ) );
    double DSy[40] __attribute__( ( aligned( 64 ) ) );
    double DSz[40] __attribute__( ( aligned( 64 ) ) );
    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );
    
    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );
    int nbVec = ( iend-istart+( cell_nparts-1 )-( ( iend-istart-1 )&( cell_nparts-1 ) ) ) / vecSize;
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }
    
    // Jx^(d,p,p)
    #pragma omp simd
    for( unsigned int j=0; j<1000; j++ ) {
        bJx[j] = 0.;
    }
    
    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {
    
        int np_computed( min( cell_nparts-ivect, vecSize ) );
        int istart0 = ( int )istart + ivect;
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            compute_distances( particles, npart_total, ipart, istart0, ipart_ref, deltaold, iold, Sx0_buff_vect, Sy0_buff_vect, Sz0_buff_vect, DSx, DSy, DSz );
            charge_weight[ipart] = inv_cell_volume * ( double )( particles.charge( istart0+ipart ) )*particles.weight( istart0+ipart );
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            computeJ( ipart, charge_weight, DSx, DSy, DSz, Sy0_buff_vect, Sz0_buff_vect, bJx, dx_ov_dt, 25, 5, 1 );
        } // END ipart (compute coeffs)
        
    } // END ivect
    
    int iglobal0 = ipom2*nyz+jpom2*nprimz+kpom2;
    
    int iglobal  = iglobal0;
    for( unsigned int i=1 ; i<5 ; i++ ) {
        iglobal += nyz;
        for( unsigned int j=0 ; j<5 ; j++ ) {
            #pragma omp simd
            for( unsigned int k=0 ; k<5 ; k++ ) {
                double tmpJx = 0.;
                int ilocal = ( ( i )*25+j*5+k )*vecSize;
#pragma unroll(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpJx += bJx [ilocal+ipart];
                }
                Jx[iglobal+j*nprimz+k]         += tmpJx;
            }
        }
    }
    
    
    // Jy^(p,d,p)
    #pragma omp simd
    for( unsigned int j=0; j<1000; j++ ) {
        bJx[j] = 0.;
    }
    
    
    cell_nparts = ( int )iend-( int )istart;
    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {
    
        int np_computed( min( cell_nparts-ivect, vecSize ) );
        int istart0 = ( int )istart + ivect;
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            compute_distances( particles, npart_total, ipart, istart0, ipart_ref, deltaold, iold, Sx0_buff_vect, Sy0_buff_vect, Sz0_buff_vect, DSx, DSy, DSz );
            charge_weight[ipart] = inv_cell_volume * ( double )( particles.charge( istart0+ipart ) )*particles.weight( istart0+ipart );
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            computeJ( ipart, charge_weight, DSy, DSx, DSz, Sx0_buff_vect, Sz0_buff_vect, bJx, dy_ov_dt, 5, 25, 1 );
        } // END ipart (compute coeffs)
    }
    
    iglobal = iglobal0+ipom2*nprimz;
    for( unsigned int i=0 ; i<5 ; i++ ) {
        for( unsigned int j=1 ; j<5 ; j++ ) {
            #pragma omp simd
            for( unsigned int k=0 ; k<5 ; k++ ) {
                double tmpJy = 0.;
                int ilocal = ( ( i )*25+j*5+k )*vecSize;
#pragma unroll(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpJy += bJx [ilocal+ipart];
                }
                Jy[iglobal+j*nprimz+k] += tmpJy;
            }
        }
        iglobal += ( nprimy+1 )*nprimz;
    }
    
    
    // Jz^(p,p,d)
    cell_nparts = ( int )iend-( int )istart;
    #pragma omp simd
    for( unsigned int j=0; j<1000; j++ ) {
        bJx[j] = 0.;
    }
    
    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {
    
        int np_computed( min( cell_nparts-ivect, vecSize ) );
        int istart0 = ( int )istart + ivect;
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            compute_distances( particles, npart_total, ipart, istart0, ipart_ref, deltaold, iold, Sx0_buff_vect, Sy0_buff_vect, Sz0_buff_vect, DSx, DSy, DSz );
            charge_weight[ipart] = inv_cell_volume * ( double )( particles.charge( istart0+ipart ) )*particles.weight( istart0+ipart );
        }
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            computeJ( ipart, charge_weight, DSz, DSx, DSy, Sx0_buff_vect, Sy0_buff_vect, bJx, dz_ov_dt, 1, 25, 5 );
        } // END ipart (compute coeffs)
        
    }
    
    iglobal = iglobal0  + jpom2 +ipom2*nprimy;
    for( unsigned int i=0 ; i<5 ; i++ ) {
        for( unsigned int j=0 ; j<5 ; j++ ) {
            #pragma omp simd
            for( unsigned int k=1 ; k<5 ; k++ ) {
                double tmpJz = 0.;
                int ilocal = ( ( i )*25+j*5+k )*vecSize;
#pragma unroll(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpJz +=  bJx[ilocal+ipart];
                }
                Jz [iglobal + ( j )*( nprimz+1 ) + k] +=  tmpJz;
            }
        }
        iglobal += nprimy*( nprimz+1 );
    }
    
    
} // END Project vectorized


// ---------------------------------------------------------------------------------------------------------------------
//! Wrapper for projection
// ---------------------------------------------------------------------------------------------------------------------
void Projector3D2OrderV::currentsAndDensityWrapper( ElectroMagn *EMfields,
        Particles &particles,
        SmileiMPI *smpi,
        int istart, int iend,
        int ithread,
        bool diag_flag,
        bool is_spectral,
        int ispec, int scell, int ipart_ref )
{
    if( istart == iend ) {
        return;    //Don't treat empty cells.
    }
    
    //Independent of cell. Should not be here
    //{
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    //}
    int iold[3];
    
    iold[0] = scell/( nscelly*nscellz )+oversize[0];
    
    iold[1] = ( ( scell%( nscelly*nscellz ) ) / nscellz )+oversize[1];
    iold[2] = ( ( scell%( nscelly*nscellz ) ) % nscellz )+oversize[2];
    
    
    // If no field diagnostics this timestep, then the projection is done directly on the total arrays
    if( !diag_flag ) {
        if( !is_spectral ) {
            double *b_Jx =  &( *EMfields->Jx_ )( 0 );
            double *b_Jy =  &( *EMfields->Jy_ )( 0 );
            double *b_Jz =  &( *EMfields->Jz_ )( 0 );
            currents( b_Jx, b_Jy, b_Jz, particles,  istart, iend, invgf, iold, &( *delta )[0], ipart_ref );
        } else {
            ERROR( "TO DO with rho" );
        }
        
        // Otherwise, the projection may apply to the species-specific arrays
    } else {
        double *b_Jx  = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
        double *b_Jy  = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
        double *b_Jz  = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;
        double *b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
        currentsAndDensity( b_Jx, b_Jy, b_Jz, b_rho, particles,  istart, iend, invgf, iold, &( *delta )[0], ipart_ref );
    }
}


void Projector3D2OrderV::susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int scell, int ipart_ref )
{

    double *Chi_envelope = &( *EMfields->Env_Chi_ )( 0 ) ;
    
    int iold[3];
    
    iold[0] = scell/( nscelly*nscellz )+oversize[0];
    iold[1] = ( ( scell%( nscelly*nscellz ) ) / nscellz )+oversize[1];
    iold[2] = ( ( scell%( nscelly*nscellz ) ) % nscellz )+oversize[2];
    
    
    std::vector<double> *Epart       = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Phipart     = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPhipart = &( smpi->dynamics_GradPHIpart[ithread] );
    std::vector<double> *inv_gamma_ponderomotive = &( smpi->dynamics_inv_gamma_ponderomotive[ithread] );
    
    int nparts = smpi->dynamics_invgf[ithread].size();
    double *Ex       = &( ( *Epart )[0*nparts] );
    double *Ey       = &( ( *Epart )[1*nparts] );
    double *Ez       = &( ( *Epart )[2*nparts] );
    double *Phi      = &( ( *Phipart )[0*nparts] );
    double *GradPhix = &( ( *GradPhipart )[0*nparts] );
    double *GradPhiy = &( ( *GradPhipart )[1*nparts] );
    double *GradPhiz = &( ( *GradPhipart )[2*nparts] );
    
    int vecSize = 8;
    unsigned int bsize = 3*3*3*vecSize; // primal grid, particles did not yet move (3x3x3 enough)
    double bChi[bsize] __attribute__( ( aligned( 64 ) ) );
    
    double Sx1[24] __attribute__( ( aligned( 64 ) ) );
    double Sy1[24] __attribute__( ( aligned( 64 ) ) );
    double Sz1[24] __attribute__( ( aligned( 64 ) ) );
    double charge_weight[8] __attribute__( ( aligned( 64 ) ) );
    
    // Closest multiple of 8 higher or equal than npart = iend-istart.
    int cell_nparts( ( int )iend-( int )istart );
    int nbVec = ( iend-istart+( cell_nparts-1 )-( ( iend-istart-1 )&( cell_nparts-1 ) ) ) / vecSize;
    if( nbVec*vecSize != cell_nparts ) {
        nbVec++;
    }
    double one_over_mass=1./species_mass;
    
    #pragma omp simd
    for( unsigned int j=0; j<bsize; j++ ) {
        bChi[j] = 0.;
    }
    
    for( int ivect=0 ; ivect < cell_nparts; ivect += vecSize ) {
    
        int np_computed( min( cell_nparts-ivect, vecSize ) );
        int istart0 = ( int )istart + ivect;
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
        
            double momentum[3];
            
            double gamma_ponderomotive, gamma0, gamma0_sq;
            double charge_over_mass_dts2, charge_sq_over_mass_sq_dts4, charge_sq_over_mass_sq;
            double pxsm, pysm, pzsm;
            
            double c = particles.charge( istart0+ipart );
            
            charge_over_mass_dts2       = c *dts2*one_over_mass;
            // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
            charge_sq_over_mass_sq_dts4 = c*c*dts4*one_over_mass*one_over_mass;
            // (charge over mass)^2
            charge_sq_over_mass_sq      = c*c*one_over_mass*one_over_mass;
            
            for( int i = 0 ; i<3 ; i++ ) {
                momentum[i] = particles.momentum( i, istart0+ipart );
            }
            
            // compute initial ponderomotive gamma
            gamma0_sq = 1. + momentum[0]*momentum[0]+ momentum[1]*momentum[1] + momentum[2]*momentum[2] + *( Phi+istart0-ipart_ref+ipart )*charge_sq_over_mass_sq ;
            gamma0    = sqrt( gamma0_sq ) ;
            
            // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
            pxsm = ( gamma0 * charge_over_mass_dts2*( *( Ex+istart-ipart_ref+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhix+istart-ipart_ref+ipart ) ) ) * momentum[0] / gamma0_sq;
            pysm = ( gamma0 * charge_over_mass_dts2*( *( Ey+istart-ipart_ref+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhiy+istart-ipart_ref+ipart ) ) ) * momentum[1] / gamma0_sq;
            pzsm = ( gamma0 * charge_over_mass_dts2*( *( Ez+istart-ipart_ref+ipart ) ) - charge_sq_over_mass_sq_dts4*( *( GradPhiz+istart-ipart_ref+ipart ) ) ) * momentum[2] / gamma0_sq;
            
            // update of gamma ponderomotive
            gamma_ponderomotive = gamma0 + ( pxsm+pysm+pzsm )*0.5 ;
            // buffer inverse of ponderomotive gamma to use it in ponderomotive momentum pusher
            ( *inv_gamma_ponderomotive )[ipart-ipart_ref] = 1./gamma_ponderomotive;
            
            // susceptibility for the macro-particle
            charge_weight[ipart] = c*c*inv_cell_volume * particles.weight( istart0+ipart )*one_over_mass/gamma_ponderomotive ;
            
            // variable declaration
            double xpn, ypn, zpn;
            double delta, delta2;
            
            // Initialize all current-related arrays to zero
            for( unsigned int i=0; i<3; i++ ) {
                Sx1[i*vecSize+ipart] = 0.;
                Sy1[i*vecSize+ipart] = 0.;
                Sz1[i*vecSize+ipart] = 0.;
            }
            
            // --------------------------------------------------------
            // Locate particles & Calculate Esirkepov coef. S, DS and W
            // --------------------------------------------------------
            
            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            xpn = particles.position( 0, istart0+ipart ) * dx_inv_;
            int ip = round( xpn );
            delta  = xpn - ( double )ip;
            delta2 = delta*delta;
            Sx1[0*vecSize+ipart] = 0.5 * ( delta2-delta+0.25 );
            Sx1[1*vecSize+ipart] = 0.75-delta2;
            Sx1[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
            
            ypn = particles.position( 1, istart0+ipart ) * dy_inv_;
            int jp = round( ypn );
            delta  = ypn - ( double )jp;
            delta2 = delta*delta;
            Sy1[0*vecSize+ipart] = 0.5 * ( delta2-delta+0.25 );
            Sy1[1*vecSize+ipart] = 0.75-delta2;
            Sy1[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
            
            zpn = particles.position( 2, istart0+ipart ) * dz_inv_;
            int kp = round( zpn );
            delta  = zpn - ( double )kp;
            delta2 = delta*delta;
            Sz1[0*vecSize+ipart] = 0.5 * ( delta2-delta+0.25 );
            Sz1[1*vecSize+ipart] = 0.75-delta2;
            Sz1[2*vecSize+ipart] = 0.5 * ( delta2+delta+0.25 );
            
        } // end ipart loop
        
        #pragma omp simd
        for( int ipart=0 ; ipart<np_computed; ipart++ ) {
            for( unsigned int i=0 ; i<3 ; i++ ) {
                for( unsigned int j=0 ; j<3 ; j++ ) {
                    int index( ( i*9 + j*3 )*vecSize+ipart );
                    for( unsigned int k=0 ; k<3 ; k++ ) {
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
    int nyz = nprimy*nprimz;
    
    int iglobal = ( ipo-1 )*nyz+( jpo-1 )*nprimz+kpo-1;
    
    for( unsigned int i=0 ; i<3 ; i++ ) {
        for( unsigned int j=0 ; j<3 ; j++ ) {
            #pragma omp simd
            for( unsigned int k=0 ; k<3 ; k++ ) {
                double tmpChi = 0.;
                int ilocal = ( i*9+j*3+k )*vecSize;
#pragma unroll(8)
                for( int ipart=0 ; ipart<8; ipart++ ) {
                    tmpChi +=  bChi[ilocal+ipart];
                }
                Chi_envelope [iglobal + j*nprimz + k] +=  tmpChi;
            }
        }
        iglobal += nyz;
    }
    
    
    
}
