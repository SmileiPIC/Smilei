#ifndef ELECTROMAGNAM_H
#define ELECTROMAGNAM_H

#include "ElectroMagn.h"
#include "Field.h"
#include "Field2D.h"
#include "cField2D.h"

class Params;

//! class ElectroMagn3D containing all information on the electromagnetic fields & currents for 3d3v simulations
class ElectroMagnAM : public ElectroMagn
{
public:
    //! Constructor for ElectroMagnAM
    ElectroMagnAM( Params &params, DomainDecomposition *domain_decomposition, std::vector<Species *> &vecSpecies, Patch *patch );
    ElectroMagnAM( ElectroMagnAM *emFields, Params &params, Patch *patch );
    
    //! Destructor for ElectroMagnAM
    ~ElectroMagnAM();
    
    std::vector<cField2D *> El_;
    std::vector<cField2D *> Er_;
    std::vector<cField2D *> Et_;
    std::vector<cField2D *> Bl_;
    std::vector<cField2D *> Br_;
    std::vector<cField2D *> Bt_;
    std::vector<cField2D *> Bl_m;
    std::vector<cField2D *> Br_m;
    std::vector<cField2D *> Bt_m;
    std::vector<cField2D *> Jl_;
    std::vector<cField2D *> Jr_;
    std::vector<cField2D *> Jt_;
    std::vector<cField2D *> rho_AM_;
    //define a vector of vectors
    std::vector<cField2D *> Jl_s;
    std::vector<cField2D *> Jr_s;
    std::vector<cField2D *> Jt_s;
    std::vector<cField2D *> rho_AM_s;
    void restartRhoJ() override;
    void restartRhoJs() override;
    
    void initPoisson( Patch *patch ) override;
    double compute_r() override;
    void compute_Ap( Patch *patch ) override;
    void compute_Ap_relativistic_Poisson( Patch *patch, double gamma_mean ) override {;}
    //Access to Ap
    double compute_pAp() override;
    void update_pand_r( double r_dot_r, double p_dot_Ap ) override;
    void update_p( double rnew_dot_rnew, double r_dot_r ) override;
    void initE( Patch *patch ) override;
    void initE_relativistic_Poisson( Patch *patch, double gamma_mean ) override {;}
    void initB_relativistic_Poisson( Patch *patch, double gamma_mean ) override {;}
    void center_fields_from_relativistic_Poisson( Patch *patch ) override {;}
    void initRelativisticPoissonFields( Patch *patch ) override {;}
    void sum_rel_fields_to_em_fields( Patch *patch ) override {;}
    void centeringE( std::vector<double> E_Add ) override;
    void centeringErel( std::vector<double> E_Add ) override {;}
    
    double getEx_Xmin() override
    {
        return 0.;
    }
    double getEx_Xmax() override
    {
        return 0.;
    }
    
    double getExrel_Xmin() override
    {
        return 0.;
    }
    double getExrel_Xmax() override
    {
        return 0.;
    }
    
    double getEx_XminYmax() override
    {
        return 0.;
    }
    double getEy_XminYmax() override
    {
        return 0.;
    }
    double getEx_XmaxYmin() override
    {
        return 0.;
    }
    double getEy_XmaxYmin() override
    {
        return 0.;
    }
    
    double getExrel_XminYmax() override
    {
        return 0.;
    }
    double getEyrel_XminYmax() override
    {
        return 0.;
    }
    double getExrel_XmaxYmin() override
    {
        return 0.;
    }
    double getEyrel_XmaxYmin() override
    {
        return 0.;
    }
    
    //! Total number of modes in Fourier poloidal decomposition.
    unsigned int nmodes;
    
    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields( bool ) override;
    
    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields() override;
    
    //! Method used to apply a single-pass binomial filter on currents
    void binomialCurrentFilter() override;
    
    //! Creates a new field with the right characteristics, depending on the name
    Field *createField( std::string fieldname ) override;
    
    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ() override;
    void addToGlobalRho( int ispec, unsigned int clrw );
    void computeTotalRhoJs( unsigned int clrw );
    
    //! Method used to compute the total susceptibility by summing over all species
    void computeTotalEnvChi() override;
    
    //! Method used to gather species densities and currents on a single array
    void synchronizePatch( unsigned int clrw );
    void finalizePatch( unsigned int clrw );
    
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

    //! Inverse radius
    //double *ptr_invR, *ptr_invRd;
    
    //! compute Poynting on borders
    void computePoynting() override;
    //! Method used to impose external fields
    void applyExternalFields( Patch *patch ) override;
    //! Method used to impose external fields
    void applyExternalField( Field *, Profile *, Patch * ) override;
    
    void initAntennas( Patch *patch ) override;
    
    //! Compute local square norm of charge denisty is not null
    double computeRhoNorm2() override
    {
        double norm2( 0 );
        for( unsigned int imode = 0 ; imode<nmodes ; imode++ ) {
            rho_AM_[imode]->norm2( istart, bufsize );
        }
        return norm2;
    }
    
    void on_axis_J( bool diag_flag );
    //! from smpi is ymax
    const bool isYmin;
    
    //! from smpi is ymin
    const bool isYmax;
    
    //! Initialize quantities needed in the creators of ElectroMagnAM
    void initElectroMagnAMQuantities( Params &params, Patch *patch );
    
    void finishInitialization( int nspecies, Patch *patch ) override final;
    
};

#endif
