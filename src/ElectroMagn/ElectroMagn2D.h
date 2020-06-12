#ifndef ELECTROMAGN2D_H
#define ELECTROMAGN2D_H

#include "ElectroMagn.h"
#include "Field.h"
#include "Field2D.h"

class Params;

//! class ElectroMagn2D containing all information on the electromagnetic fields & currents for 2Dcartesian simulations
class ElectroMagn2D : public ElectroMagn
{
public:
    //! Constructor for ElectroMagn2D
    ElectroMagn2D( Params &params, DomainDecomposition *domain_decomposition, std::vector<Species *> &vecSpecies, Patch *patch );
    ElectroMagn2D( ElectroMagn2D *emFields, Params &params, Patch *patch );
    
    //! Destructor for ElectroMagn2D
    ~ElectroMagn2D();
    
    // --------------------------------------
    //  --------- PATCH IN PROGRESS ---------
    // --------------------------------------
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
        return ( *Ex_ )( 0, ny_p-1 );
    }
    double getEy_XminYmax()
    {
        return ( *Ey_ )( 0, ny_d-1 );
    }
    double getEx_XmaxYmin()
    {
        return ( *Ex_ )( nx_d-1, 0 );
    }
    double getEy_XmaxYmin()
    {
        return ( *Ey_ )( nx_p-1, 0 );
    }
    
    double getExrel_XminYmax()
    {
        return ( *Ex_rel_ )( 0, ny_p-1 );
    }
    double getEyrel_XminYmax()
    {
        return ( *Ey_rel_ )( 0, ny_d-1 );
    }
    double getExrel_XmaxYmin()
    {
        return ( *Ex_rel_ )( nx_d-1, 0 );
    }
    double getEyrel_XmaxYmin()
    {
        return ( *Ey_rel_ )( nx_p-1, 0 );
    }
    
    // --------------------------------------
    //  --------- PATCH IN PROGRESS ---------
    // --------------------------------------
    
//    //! Method used to solve Maxwell-Ampere equation
//    void solveMaxwellAmpere();

    //! Method used to save the Magnetic fields (used to center them)
    void saveMagneticFields( bool );
    
    //! Method used to center the Magnetic fields (used to push the particles)
    void centerMagneticFields();
    
    //! Method used to apply a single-pass binomial filter on currents
    void binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes);

    //! Method used to apply a single-pass custom FIR based filter on currents
    void customFIRCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes,std::vector<double> filtering_coeff);
    
    //! Creates a new field with the right characteristics, depending on the name
    Field *createField( std::string fieldname, Params& params );
    
    //! Method used to compute the total charge density and currents by summing over all species
    void computeTotalRhoJ();
    void addToGlobalRho( int ispec, unsigned int clrw );
    
    //! Method used to compute the total susceptibility by summing over all species
    void computeTotalEnvChi();
    
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
    
    //! Spatial step dx for 2D3V cartesian simulations
    double dx;
    
    //! Spatial step dy for 2D3V cartesian simulations
    double dy;
    
    //! Ratio of the time-step by the spatial-step dt/dx for 2D3V cartesian simulations
    double dt_ov_dx;
    
    //! Ratio of the time-step by the spatial-step dt/dy for 2D3V cartesian simulations
    double dt_ov_dy;
    
    //! Ratio of the spatial-step by the time-step dx/dt for 2D3V cartesian simulations
    double dx_ov_dt;
    
    //! Ratio of the spatial-step by the time-step dy/dt for 2D3V cartesian simulations
    double dy_ov_dt;
    
    //! compute Poynting on borders
    void computePoynting();
    
    //! Method used to impose external fields
    void applyExternalField( Field *, Profile *, Patch * );
    
    void initAntennas( Patch *patch, Params& params );
    //! Method used to impose external time fields
    void applyPrescribedField( Field *, Profile *, Patch *, double time );
    
    
    //! from smpi is ymax
    const bool isYmin;
    
    //! from smpi is ymin
    const bool isYmax;
    
private:


    //! Initialize quantities needed in the creators of ElectroMagn2D
    void initElectroMagn2DQuantities( Params &params, Patch *patch );
};

#endif
