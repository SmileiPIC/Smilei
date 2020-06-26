#ifndef INTERPOLATOR2D2ORDER_H
#define INTERPOLATOR2D2ORDER_H


#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D2Order : public Interpolator2D
{

public:
    Interpolator2D2Order( Params &, Patch * );
    ~Interpolator2D2Order() override final {};
    
    inline void fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc );
    void fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc ) override final ;
    void fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final ;
    void fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection ) override final;
    void oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1=NULL, double *l2=NULL, double *l3=NULL ) override final;
    
    inline double compute( double *coeffx, double *coeffy, Field2D *f, int idx, int idy )
    {
        double interp_res( 0. );
        //unroll ?
        for( int iloc=-1 ; iloc<2 ; iloc++ ) {
            for( int jloc=-1 ; jloc<2 ; jloc++ ) {
                interp_res += *( coeffx+iloc ) * *( coeffy+jloc ) * ( *f )( idx+iloc, idy+jloc );
            }
        }
        return interp_res;
    };
    
    void fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc ) override final;
    void envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;    

private:
    inline void coeffs( double xpn, double ypn )
    {
        // Indexes of the central nodes
        ip_ = round( xpn );
        id_ = round( xpn+0.5 );
        jp_ = round( ypn );
        jd_ = round( ypn+0.5 );
        
        // Declaration and calculation of the coefficient for interpolation
        double delta2;
        
        deltax   = xpn - ( double )id_ + 0.5;
        delta2  = deltax*deltax;
        coeffxd_[0] = 0.5 * ( delta2-deltax+0.25 );
        coeffxd_[1] = 0.75 - delta2;
        coeffxd_[2] = 0.5 * ( delta2+deltax+0.25 );
        
        deltax   = xpn - ( double )ip_;
        delta2  = deltax*deltax;
        coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
        coeffxp_[1] = 0.75 - delta2;
        coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );
        
        deltay   = ypn - ( double )jd_ + 0.5;
        delta2  = deltay*deltay;
        coeffyd_[0] = 0.5 * ( delta2-deltay+0.25 );
        coeffyd_[1] = 0.75 - delta2;
        coeffyd_[2] = 0.5 * ( delta2+deltay+0.25 );
        
        deltay   = ypn - ( double )jp_;
        delta2  = deltay*deltay;
        coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
        coeffyp_[1] = 0.75 - delta2;
        coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );
        
        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        id_ = id_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        jd_ = jd_ - j_domain_begin;
    }
    
    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax, deltay;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3];
    
    
};//END class

#endif
