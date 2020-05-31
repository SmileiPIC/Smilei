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
    ElectroMagn3D( Params &params, DomainDecomposition *domain_decomposition, std::vector<Species *> &vecSpecies, Patch *patch );
    ElectroMagn3D( ElectroMagn3D *emFields, Params &params, Patch *patch );
    
    //! Destructor for ElectroMagn3D
    ~ElectroMagn3D();
    
    void initPoisson( Patch *patch );
    double compute_r();
    void compute_Ap( Patch *patch );
    void compute_Ap_relativistic_Poisson( Patch *patch, double gamma_mean );
    //Access to Ap
    double compute_pAp();
    void update_pand_r( double r_dot_r, double p_dot_Ap );
    void update_p( double rnew_dot_rnew, double r_dot_r );
    void initE( Patch *patch );
    void initE_relativistic_Poisson( Patch *patch, double gamma_mean );
    void initB_relativistic_Poisson( Patch *patch, double gamma_mean );
    void center_fields_from_relativistic_Poisson( Patch *patch );
    void initRelativisticPoissonFields( Patch *patch );
    void sum_rel_fields_to_em_fields( Patch *patch );
    void centeringE( std::vector<double> E_Add );
    void centeringErel( std::vector<double> E_Add );
    
    double getEx_Xmin()
    {
        return 0.;
    }
    double getEx_Xmax()
    {
        return 0.;
    }
    
    double getExrel_Xmin()
    {
        return 0.;
    }
    double getExrel_Xmax()
    {
        return 0.;
    }
    
    double getEx_XminYmax()
    {
        return 0.;
    }
    double getEy_XminYmax()
    {
        return 0.;
    }
    double getEx_XmaxYmin()
    {
        return 0.;
    }
    double getEy_XmaxYmin()
    {
        return 0.;
    }
    
    double getExrel_XminYmax()
    {
        return 0.;
    }
    double getEyrel_XminYmax()
    {
        return 0.;
    }
    double getExrel_XmaxYmin()
    {
        return 0.;
    }
    double getEyrel_XmaxYmin()
    {
        return 0.;
    }
    
//    //! Method used to solve Maxwell-Ampere equation
//    void solveMaxwellAmpere();

    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields( bool );
    
    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields();
    
    //! Method used to apply a single-pass binomial filter on currents
    void binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes);

    //! Method used to apply a single-pass custom FIR based filter on currents
    void customFIRCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes, std::vector<double> filtering_coeff);
 
    //! Creates a new field with the right characteristics, depending on the name
    Field *createField( std::string fieldname, Params& params );
    
    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ();
    void addToGlobalRho( int ispec, unsigned int clrw );
    
    //! Method used to compute the total susceptibility by summing over all species
    void computeTotalEnvChi();
    //void addToGlobalEnvChi(int ispec, unsigned int clrw);
    //void computeTotalEnvChis(unsigned int clrw);
    
    //! Method used to gather species densities and currents on a single array
    void synchronizePatch( unsigned int clrw );
    void finalizePatch( unsigned int clrw );
    
    //! \todo Create properties the laser time-profile (MG & TV)
    
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;
    
    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;
    
    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;
    
    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;
    
    //! Number of nodes on the primal grid in the z-direction
    unsigned int nz_p;
    
    //! Number of nodes on the dual grid in the z-direction
    unsigned int nz_d;
    
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
    void computePoynting();
    
    //! Method used to impose external fields
    void applyExternalField( Field *, Profile *, Patch * );
    
    //! Method used to impose external fields
    void applyPrescribedField( Field *, Profile *, Patch *, double time );
    
    void initAntennas( Patch *patch );
    
    void initAntennas( Patch* patch, Params& params );
    
    //! from smpi is ymax
    const bool isYmin;
    
    //! from smpi is ymin
    const bool isYmax;
    
    //! from smpi is zmax
    const bool isZmax;
    
    //! from smpi is zmin
    const bool isZmin;
    
private:


    //! Initialize quantities needed in the creators of ElectroMagn3D
    void initElectroMagn3DQuantities( Params &params, Patch *patch );
};

#endif
