#ifndef ELECTROMAGN3D_H
#define ELECTROMAGN3D_H

#include "ElectroMagn.h"
#include "Field.h"
#include "Field3D.h"

class Params;

//! class ElectroMagn3D containing all information on the electromagnetic fields & currents for 3Dcartesian simulations
class ElectroMagn3D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn3D
    ElectroMagn3D( Params &params, std::vector<Species *> &vecSpecies, Patch *patch );
    ElectroMagn3D( ElectroMagn3D *emFields, Params &params, Patch *patch );

    //! Destructor for ElectroMagn3D
    ~ElectroMagn3D();

    void initPoisson( Patch *patch ) override;
    double compute_r() override;
    void compute_Ap( Patch *patch ) override;
    void compute_Ap_relativistic_Poisson( Patch *patch, double gamma_mean ) override;
    //Access to Ap
    double compute_pAp() override;
    void update_pand_r( double r_dot_r, double p_dot_Ap ) override;
    void update_p( double rnew_dot_rnew, double r_dot_r ) override;
    void initE( Patch *patch ) override;
    void initE_relativistic_Poisson( Patch *patch, double gamma_mean ) override;
    void initB_relativistic_Poisson( double gamma_mean ) override;
    void center_fields_from_relativistic_Poisson() override;
    void initRelativisticPoissonFields() override;
    void sum_rel_fields_to_em_fields() override;
    void centeringE( std::vector<double> E_Add ) override;
    void centeringErel( std::vector<double> E_Add ) override;

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

//    //! Method used to solve Maxwell-Ampere equation
//    void solveMaxwellAmpere();

    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields( bool ) override;

    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields() override;

    //! Method used to apply a single-pass binomial filter on currents
    void binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes) override;

    //! Method used to apply a single-pass custom FIR based filter on currents
    void customFIRCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes, std::vector<double> filtering_coeff) override;

    //! Creates a new field with the right characteristics, depending on the name
    Field *createField( std::string fieldname, Params& params ) override;

    //! Method used to compute the total charge density and currents by summing over all species on CPU (Host)
    void computeTotalRhoJ() override;

// #if defined( SMILEI_ACCELERATOR_GPU )
//     //! Method used to compute the total charge density and currents by summing over all species on Device
//     void computeTotalRhoJOnDevice() override;
// #endif


    void addToGlobalRho( int ispec, unsigned int clrw );

    //! Method used to compute the total susceptibility by summing over all species
    void computeTotalEnvChi() override;
    //void addToGlobalEnvChi(int ispec, unsigned int clrw);
    //void computeTotalEnvChis(unsigned int clrw);

    //! Method used to gather species densities and currents on a single array
    void synchronizePatch( unsigned int clrw );
    void finalizePatch( unsigned int clrw );

    //! Spatial step dx for 3D3V cartesian simulations
    double dx;

    //! Spatial step dy for 3D3V cartesian simulations
    double dy;

    //! Spatial step dz for 3D3V cartesian simulations
    double dz;

    //! Ratio of the time-step by the spatial-step dt/dx for 3D3V cartesian simulations
    double dt_ov_dx;

    //! Ratio of the time-step by the spatial-step dt/dy for 3D3V cartesian simulations
    double dt_ov_dy;

    //! Ratio of the time-step by the spatial-step dt/dz for 3D3V cartesian simulations
    double dt_ov_dz;

    //! Ratio of the spatial-step by the time-step dx/dt for 3D3V cartesian simulations
    double dx_ov_dt;

    //! Ratio of the spatial-step by the time-step dy/dt for 3D3V cartesian simulations
    double dy_ov_dt;

    //! Ratio of the spatial-step by the time-step dz/dt for 3D3V cartesian simulations
    double dz_ov_dt;

    //! compute Poynting on borders
    void computePoynting( unsigned int axis, unsigned int side ) override;

    //! Method used to impose external fields
    void applyExternalField( Field *, Profile *, Patch * ) override;

    //! Method used to impose external fields
    void applyPrescribedField( Field *, Profile *, Patch *, double time ) override;

    void initAntennas( Patch *patch );

    void initAntennas( Patch* patch, Params& params ) override;

    const bool isYmin, isYmax, isZmin, isZmax;

    // copy currents projected on sub-buffers to global currents
    void copyInLocalDensities(int ispec, int ibin,
                              double* b_Jx, double* b_Jy, double* b_Jz, double* b_rho,
                              std::vector<unsigned int> b_dim, bool diag_flag) override final;

    // copy susceptibility projected on sub-buffers to global susceptibility
    void copyInLocalSusceptibility(int ispec, int ibin,
                                   double* b_Chi, std::vector<unsigned int> b_dim, bool diag_flag) override final;

private:


    //! Initialize quantities needed in the creators of ElectroMagn3D
    void initElectroMagn3DQuantities( Params &params, Patch *patch );
};

#endif
