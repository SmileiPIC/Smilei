#ifndef ELECTROMAGNAM_H
#define ELECTROMAGNAM_H

#include "ElectroMagn.h"
#include "Field.h"
#include "cField.h"
#include "Field2D.h"
#include "cField2D.h"

class Params;

//! class ElectroMagn3D containing all information on the electromagnetic fields & currents for 3d3v simulations
class ElectroMagnAM : public ElectroMagn
{
public:
    //! Constructor for ElectroMagnAM
    ElectroMagnAM( Params &params, std::vector<Species *> &vecSpecies, Patch *patch );
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
    std::vector<cField2D *> rho_old_AM_;
    //define a vector of vectors
    std::vector<cField2D *> Jl_s;
    std::vector<cField2D *> Jr_s;
    std::vector<cField2D *> Jt_s;
    std::vector<cField2D *> rho_AM_s;
    void restartRhoJ() override;
    void restartRhoold() ;
    void restartRhos() ;
    void restartRhoJs() override;
    
    // fields for Poisson solver
    cField2D *El_Poisson_;
    cField2D *Er_Poisson_;
    cField2D *Et_Poisson_;

    void initPoisson( Patch *patch ) override;
    double compute_r() override;
    void compute_Ap( Patch *patch ) override;
    void compute_Ap_relativistic_Poisson( Patch *, double ) override {};
    void compute_Ap_relativistic_Poisson_AM( Patch *patch, double gamma_mean, unsigned int imode );
    void compute_Ap_Poisson_AM( Patch *patch, unsigned int imode );
    //Access to Ap
    double compute_pAp() override { return 0.; };
    std::complex<double> compute_pAp_AM();
    void update_pand_r( double, double ) override {};
    void update_p( double rnew_dot_rnew, double r_dot_r ) override;
    void update_pand_r_AM( double r_dot_r, std::complex<double> p_dot_Ap );
    void initE( Patch *patch ) override;
    void delete_phi_r_p_Ap();
    void delete_relativistic_fields();
    void delete_Poisson_fields();
    void initE_relativistic_Poisson( Patch *, double ) override {};
    void initE_relativistic_Poisson_AM( double gamma_mean, unsigned int imode );
    void initE_Poisson_AM( unsigned int imode );
    void initB_relativistic_Poisson( double ) override {};
    void initB_relativistic_Poisson_AM( double gamma_mean );
    void center_fields_from_relativistic_Poisson() override {};
    void center_fields_from_relativistic_Poisson_AM();
    void initRelativisticPoissonFields() override;
    void initPoissonFields();
    void initPoisson_init_phi_r_p_Ap( unsigned int imode );
    void sum_rel_fields_to_em_fields() override {};
    void sum_rel_fields_to_em_fields_AM( Params &params, unsigned int imode );
    void sum_Poisson_fields_to_em_fields_AM( unsigned int imode );
    void centeringE( std::vector<double> E_Add ) override;
    void centeringErel( std::vector<double> ) override {};
    
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
    void binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes) override;

    //! Method used to apply a single-pass custom FIR based filter on currents
    void customFIRCurrentFilter(unsigned int, std::vector<unsigned int>, std::vector<double> ) override {};
 
    //! Creates a new field with the right characteristics, depending on the name
    Field *createField( std::string fieldname, Params& params ) override;
    
    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ() override;
    void computeTotalRhoold();
    void addToGlobalRho( int ispec, unsigned int clrw );
    
    //! Method used to compute the total susceptibility by summing over all species
    void computeTotalEnvChi() override;
    
    //! Method used to gather species densities and currents on a single array
    void synchronizePatch( unsigned int clrw );
    void finalizePatch( unsigned int clrw );
    
    //! Spatial step dl for 2D3V cylindrical simulations
    double dl;
    
    //! Spatial step dr for 2D3V cylindrical simulations
    double dr;
    
    //! Ratio of the time-step by the spatial-step dt/dl for 2D3V cylindrical simulations
    double dt_ov_dl;
    
    //! Ratio of the time-step by the spatial-step dt/dr for 2D3V cylindrical simulations
    double dt_ov_dr;
    
    //! Ratio of the spatial-step by the time-step dl/dt for 2D3V cylindrical simulations
    double dl_ov_dt;
    
    //! Ratio of the spatial-step by the time-step dr/dt for 2D3V cylindrical simulations
    double dr_ov_dt;
    //! Minimum radius in the current patch
    int j_glob_;
    
    //! compute Poynting on borders
    void computePoynting( unsigned int axis, unsigned int side ) override;
    //! Method used to impose external fields
    void applyExternalFields( Patch *patch ) override;
    //! Method used to impose external time fields
    void applyPrescribedFields( Patch *patch, double time ) override;
    //! Method used to impose one external field
    void applyExternalField( Field *, Profile *, Patch * ) override;
    //! Method used to impose one external time field
    void applyPrescribedField( Field *, Profile *, Patch *, double time ) override;
    
    void initAntennas( Patch *patch, Params& params ) override;
    
    //! Compute local square norm of charge denisty is not null
    double computeRhoNorm2() override
    {
        double norm2( 0 );
        for( unsigned int imode = 0 ; imode<nmodes ; imode++ ) {
            norm2 += rho_AM_[imode]->norm2( istart, bufsize );
        }
        return norm2;
    }
    
    //! Compute total energy in fields
    double computeEnergy() override
    {
        double nrj( 0. );
        
        for( unsigned int imode=0; imode<nmodes; imode++ ) {
            nrj += El_[imode]->norm2_cylindrical( istart, bufsize, j_glob_ );
            nrj += Er_[imode]->norm2_cylindrical( istart, bufsize, j_glob_ );
            nrj += Et_[imode]->norm2_cylindrical( istart, bufsize, j_glob_ );
            
            nrj += Bl_m[imode]->norm2_cylindrical( istart, bufsize, j_glob_ );
            nrj += Br_m[imode]->norm2_cylindrical( istart, bufsize, j_glob_ );
            nrj += Bt_m[imode]->norm2_cylindrical( istart, bufsize, j_glob_ );
            
            if( imode == 0 ) {
                nrj *= 2;
            }
        }
        
        return 0.5 * nrj * dr; // multiplication by cell_volume is done in scalars
    }
    
    //! from smpi is ymax
    const bool isYmin;
    
    //! from smpi is ymin
    const bool isYmax;
    
    //! Initialize quantities needed in the creators of ElectroMagnAM
    void initElectroMagnAMQuantities( Params &params, Patch *patch );
    
    void finishInitialization( int nspecies, Patch *patch ) override final;
    
    //!Pointers toward R inverse values stored in patch
    double *invR, *invRd;

    // copy currents projected on sub-buffers to global currents
    void copyInLocalDensities(int, int, 
                              double*, double*, double*, double*, 
                              std::vector<unsigned int>, bool ) override final {};

    // copy currents projected on sub-buffers to global currents
    void copyInLocalAMDensities(int ispec, int ibin, 
                              std::complex<double> *b_Jl, std::complex<double> *b_Jr, 
                              std::complex<double> *b_Jt, std::complex<double> *b_rhoAM, 
                              std::vector<unsigned int> b_dim, bool diag_flag) override final;  

    // copy susceptibility projected on sub-buffers to global susceptibility
    void copyInLocalSusceptibility(int ispec, int ibin, 
                                   double* b_Chi, std::vector<unsigned int> b_dim, bool diag_flag) override final;

};

#endif
