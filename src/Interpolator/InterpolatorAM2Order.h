#ifndef INTERPOLATORRZ2ORDER_H
#define INTERPOLATORRZ2ORDER_H


#include "InterpolatorAM.h"
#include "cField2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class InterpolatorAM2Order : public InterpolatorAM
{

public:
    InterpolatorAM2Order(Params&, Patch*);
    ~InterpolatorAM2Order() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;


    inline std::complex<double> compute( double* coeffx, double* coeffy, cField2D* f, int idx, int idy) {
        std::complex<double> interp_res(0.);
        for (int iloc=-1 ; iloc<2 ; iloc++) {
            for (int jloc=-1 ; jloc<2 ; jloc++) {
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc)) ;

                //std::cout<<"f "<<std::fixed << std::setprecision(3)<<(*f)(idx+iloc,idy+jloc)<<std::endl;
            }
        }
        //std::cout<<"interp res "<< interp_res <<std::endl;
        return interp_res;
    };


    inline std::complex<double> compute_p( double* coeffx, double* coeffy, cField2D* f, int idx, int idy) {
	std::complex<double> interp_res(0.);
        //unroll ?
	for (int iloc=-1 ; iloc<2 ; iloc++) {
	    for (int jloc=-1 ; jloc<2 ; jloc++) {
                #ifdef _TODO_RZ
                #endif
                //std::cout<< "idy+jloc "<< idy+jloc << std::endl;
                //interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ( real( (*f)(idx+iloc,idy+jloc) ) + imag( (*f)(idx+iloc,idy+jloc) ) );
                if (idy+jloc+j_domain_begin ==0){
                    interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc))*6./dr ;
                }
                else{
                    interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc))/((idy+jloc+j_domain_begin)*dr);           
                }
	    }
	}
	return interp_res;
    };
  
    inline std::complex<double> compute_d( double* coeffx, double* coeffy, cField2D* f, int idx, int idy) {
        std::complex<double> interp_res(0.);
        for (int iloc=-1 ; iloc<1 ; iloc++) {
            for (int jloc=-1 ; jloc<1 ; jloc++) {
                if (jloc+idy+j_domain_begin==0) {
                    //std::cout<<"jloc+idy are on axes"<<std::endl;
                }
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc))/((idy+jloc+j_domain_begin+0.5)*dr) ;
            }
        }
        return interp_res;
    };

    inline std::complex<double> compute_0_T( double* coeffx, double* coeffy, cField2D* f, int idx, int idy) {
        std::complex<double> interp_res(0.);
        for (int iloc=-1 ; iloc<2 ; iloc++) {
            for (int jloc=-1 ; jloc<2 ; jloc++) {
                if (jloc+idy+j_domain_begin==0) {
                    interp_res -= *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc)) ;
                }else{
                    interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc)) ;
                }
            }
        }
        return interp_res;
    };
    inline std::complex<double> compute_0_L( double* coeffx, double* coeffy, cField2D* f, int idx, int idy) {
        std::complex<double> interp_res(0.);
        for (int iloc=-1 ; iloc<2 ; iloc++) {
            for (int jloc=-1 ; jloc<2 ; jloc++) {
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc)) ;
            }
        } 
        return interp_res;
    };  
    inline std::complex<double> compute_1_T( double* coeffx, double* coeffy, cField2D* f, int idx, int idy, std::complex<double>* exptheta ) {
        std::complex<double> interp_res(0.);
        for (int iloc=-1 ; iloc<2 ; iloc++) {
            for (int jloc=-1 ; jloc<2 ; jloc++) {
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc))*(*exptheta) ;
            }
        } 
        return interp_res;
    };
    inline std::complex<double> compute_1_L( double* coeffx, double* coeffy, cField2D* f, int idx, int idy, std::complex<double>* exptheta ) {
        std::complex<double> interp_res(0.);
        for (int iloc=-1 ; iloc<2 ; iloc++) {
            for (int jloc=-1 ; jloc<2 ; jloc++) {
                if (jloc+idy+j_domain_begin==0) {
                    interp_res -= *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc)) ;
                }else{
                    interp_res += *(coeffx+iloc) * *(coeffy+jloc) * ((*f)(idx+iloc,idy+jloc))*(*exptheta);
                }
            }
        }
        return interp_res;
    }; 
private:
    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax, deltar ;
    // exp m theta
    std::complex<double> exp_m_theta;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3];
    //! Number of modes;
    unsigned int nmodes;


};//END class

#endif
