#ifndef ELECTROMAGN3DRZ_H
#define ELECTROMAGN3DRZ_H

#include "ElectroMagn.h"
#include "Field.h"
#include "Field2D.h"
#include "cField2D.h"

class Params;

//! class ElectroMagn3D containing all information on the electromagnetic fields & currents for 3d3v simulations
class ElectroMagn3DRZ : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn3DRZ
    ElectroMagn3DRZ(Params &params, DomainDecomposition* domain_decomposition, std::vector<Species*>& vecSpecies, Patch* patch);
    ElectroMagn3DRZ( ElectroMagn3DRZ* emFields, Params &params, Patch* patch );

    //! Destructor for ElectroMagn3DRZ
    ~ElectroMagn3DRZ();

    std::vector<cField2D*> El_;
    std::vector<cField2D*> Er_;
    std::vector<cField2D*> Et_;
    std::vector<cField2D*> Bl_;
    std::vector<cField2D*> Br_;
    std::vector<cField2D*> Bt_;
    std::vector<cField2D*> Bl_m;
    std::vector<cField2D*> Br_m;
    std::vector<cField2D*> Bt_m;
    std::vector<cField2D*> Jl_;
    std::vector<cField2D*> Jr_;
    std::vector<cField2D*> Jt_;
    std::vector<cField2D*> rho_RZ_;
    //define a vector of vectors
    std::vector<cField2D*> Jl_s;
    std::vector<cField2D*> Jr_s;
    std::vector<cField2D*> Jt_s;
    std::vector<cField2D*> rho_RZ_s;
    void restartRhoJ() override;
    void restartRhoJs() override;

    void initPoisson(Patch *patch);
    double compute_r();
    void compute_Ap(Patch *patch);
    void compute_Ap_relativistic_Poisson(Patch* patch, double gamma_mean) {;}
    //Access to Ap
    double compute_pAp();
    void update_pand_r(double r_dot_r, double p_dot_Ap);
    void update_p(double rnew_dot_rnew, double r_dot_r);
    void initE(Patch *patch);
    void initE_relativistic_Poisson(Patch *patch, double gamma_mean) {;}
    void initB_relativistic_Poisson(Patch *patch, double gamma_mean) {;}
    void center_fields_from_relativistic_Poisson(Patch *patch) {;}
    void initRelativisticPoissonFields(Patch *patch) {;}
    void sum_rel_fields_to_em_fields(Patch *patch) {;}
    void centeringE( std::vector<double> E_Add );
    void centeringErel( std::vector<double> E_Add ) {;}

    double getEx_Xmin() { return 0.; }
    double getEx_Xmax() { return 0.; }

    double getExrel_Xmin() { return 0.; }
    double getExrel_Xmax() { return 0.; }

    double getEx_XminYmax() { return 0.; }
    double getEy_XminYmax() { return 0.; }
    double getEx_XmaxYmin() { return 0.; }
    double getEy_XmaxYmin() { return 0.; }
  
    double getExrel_XminYmax() { return 0.; }
    double getEyrel_XminYmax() { return 0.; }
    double getExrel_XmaxYmin() { return 0.; }
    double getEyrel_XmaxYmin() { return 0.; }

    //! Total number of modes in Fourier poloidal decomposition.
    const unsigned int nmodes = 2;

    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields(bool);

    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields();
    
    //! Method used to apply a single-pass binomial filter on currents
    void binomialCurrentFilter();
    
    //! Creates a new field with the right characteristics, depending on the name
    Field * createField(std::string fieldname);
    
    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ();
    void addToGlobalRho(int ispec, unsigned int clrw);
    void computeTotalRhoJs(unsigned int clrw);
 
    //! Method used to compute the total susceptibility by summing over all species
    void computeTotalEnvChi();

    //! Method used to gather species densities and currents on a single array
    void synchronizePatch(unsigned int clrw);
    void finalizePatch(unsigned int clrw);

    //! \todo Create properties the laser time-profile (MG & TV)

    //! Number of nodes on the primal grid in the x-direction
    unsigned int nl_p;

    //! Number of nodes on the dual grid in the x-direction
    unsigned int nl_d;

    //! Number of nodes on the primal grid in the y-direction
    unsigned int nr_p;

    //! Number of nodes on the dual grid in the y-direction
    unsigned int nr_d;

    //! Spatial step dl for 3D3V cartesian simulations
    double dl;

    //! Spatial step dr for 3D3V cartesian simulations
    double dr;

    //! Ratio of the time-step by the spatial-step dt/dl for 3D3V cartesian simulations
    double dt_ov_dl;

    //! Ratio of the time-step by the spatial-step dt/dr for 3D3V cartesian simulations
    double dt_ov_dr;

    //! Ratio of the spatial-step by the time-step dl/dt for 3D3V cartesian simulations
    double dl_ov_dt;

    //! Ratio of the spatial-step by the time-step dr/dt for 3D3V cartesian simulations
    double dr_ov_dt;
    //! Minimum radius in the current patch
    int j_glob_;

    //! compute Poynting on borders
    void computePoynting();
    //! Method used to impose external fields
    void applyExternalFields(Patch* patch) override;
    //! Method used to impose external fields
    void applyExternalField(Field*, Profile*, Patch*);
    
    void initAntennas(Patch* patch);
    
    //! Compute local square norm of charge denisty is not null
    double computeRhoNorm2() override {
        double norm2(0);
        for (unsigned int imode = 0 ; imode<nmodes ; imode++ ) 
            rho_RZ_[imode]->norm2(istart, bufsize);
        return norm2;
    }

    //! Fold EM fields modes correctly around axis
    void fold_fields(bool diag_flag);

    
    //! from smpi is ymax
    const bool isYmin;
    
    //! from smpi is ymin
    const bool isYmax;
    
    //! Initialize quantities needed in the creators of ElectroMagn3DRZ
    void initElectroMagn3DRZQuantities(Params &params, Patch* patch);

    void finishInitialization(int nspecies, Patch* patch) override final;

};

#endif
