#ifndef ELECTROMAGN_H
#define ELECTROMAGN_H

#include <vector>
#include <string>
#include <map>

#include "Field.h"
#include "Tools.h"
#include "Profile.h"
#include "Species.h"


class Params;
class Projector;
class Laser;
class ElectroMagnBC;
class SimWindow;
class Patch;
class Solver;
class DomainDecomposition;
class LaserEnvelope;


inline std::string LowerCase( std::string in )
{
    std::string out=in;
    std::transform( out.begin(), out.end(), out.begin(), ::tolower );
    return out;
}

// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each ExtField
// ---------------------------------------------------------------------------------------------------------------------
struct ExtField {
    //! field to which apply the external field
    std::string field;
    
    Profile *profile;
    
    unsigned int index;
};

// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each PrescribedField
// ---------------------------------------------------------------------------------------------------------------------
struct PrescribedField { 
   
    Profile *profile;
    
    Field *savedField;
    
    unsigned int index;
};

// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each Antenna
// ---------------------------------------------------------------------------------------------------------------------
struct Antenna {
    //! Jx, Jy or Jz
    std::string fieldName;
    
    Profile *time_profile;
    Profile *space_profile;
    Profile *space_time_profile;
    
    Field *field;
    
    unsigned int index;
    
    bool spacetime;
};


//! Structure for the Friedman field filters
struct FriedmanFields {
    ~FriedmanFields() {
        for( unsigned int i=0; i<Ex_.size(); i++ ) {
            delete Ex_[i];
        }
        for( unsigned int i=0; i<Ey_.size(); i++ ) {
            delete Ey_[i];
        }
        for( unsigned int i=0; i<Ez_.size(); i++ ) {
            delete Ez_[i];
        }
    }
    std::vector<Field *> Ex_;
    std::vector<Field *> Ey_;
    std::vector<Field *> Ez_;
};


//! class ElectroMagn: generic class containing all information on the electromagnetic fields and currents
class ElectroMagn
{
    friend class SmileiMPI;
    friend class SimWindow;
    friend class Checkpoint;
    friend class DiagnosticScalar;
    
public:
    //! Constructor for Electromagn
    ElectroMagn( Params &params, std::vector<Species *> &vecSpecies, Patch *patch );
    ElectroMagn( ElectroMagn *emFields, Params &params, Patch *patch );
    void initElectroMagnQuantities();
    //! Extra initialization. Used in ElectroMagnFactory
    virtual void finishInitialization( int nspecies, Patch *patch );

    // copy currents projected on sub-buffers to global currents
    virtual void copyInLocalDensities(int ispec, int ibin, 
                                      double* b_Jx, double* b_Jy, double* b_Jz, double* b_rho, 
                                      std::vector<unsigned int> b_dim, bool diag_flag) = 0;

    // copy currents projected on sub-buffers to global currents
    virtual void copyInLocalAMDensities(int, int, 
                                        std::complex<double> *, std::complex<double> *, 
                                        std::complex<double> *, std::complex<double> *, 
                                        std::vector<unsigned int>, bool ){}; 

    // copy susceptibility projected on sub-buffers to global susceptibility
    virtual void copyInLocalSusceptibility(int ispec, int ibin, 
                                           double* b_Chi, std::vector<unsigned int> b_dim, bool diag_flag) = 0;

    //! Destructor for Electromagn
    virtual ~ElectroMagn();
    
    void updateGridSize( Params &params, Patch *patch );
    
    void clean();
    
    std::vector<unsigned int> dimPrim;
    std::vector<unsigned int> dimDual;
    
    std::vector<unsigned int> index_bc_min;
    std::vector<unsigned int> index_bc_max;
    
    //! time-step (from Params)
    const double timestep;
    
    //! cell length (from Params)
    const std::vector<double> cell_length;
    
    //! x-component of the electric field
    Field *Ex_;
    
    //! y-component of the electric field
    Field *Ey_;
    
    //! z-component of the electric field
    Field *Ez_;
    
    //! x-component of the magnetic field
    Field *Bx_;
    
    //! y-component of the magnetic field
    Field *By_;
    
    //! z-component of the magnetic field
    Field *Bz_;
    
    //! x-component of the time-centered magnetic field
    Field *Bx_m;
    
    //! y-component of the time-centered magnetic field
    Field *By_m;
    
    //! z-component of the time-centered magnetic field
    Field *Bz_m;
    
    //! x-component of the total charge current
    Field *Jx_;
    
    //! y-component of the total charge current
    Field *Jy_;
    
    //! z-component of the total charge current
    Field *Jz_;
    
    //! Total charge density
    Field *rho_;
    Field *rhoold_;
    
    // Fields for relativistic Initialization
    Field *Ex_rel_;
    Field *Ey_rel_;
    Field *Ez_rel_;
    Field *Bx_rel_;
    Field *By_rel_;
    Field *Bz_rel_;
    Field *Bx_rel_t_minus_halfdt_;
    Field *By_rel_t_minus_halfdt_;
    Field *Bz_rel_t_minus_halfdt_;
    Field *Bx_rel_t_plus_halfdt_;
    Field *By_rel_t_plus_halfdt_;
    Field *Bz_rel_t_plus_halfdt_;

    // Fields for relativistic Initialization in AM
    cField *El_rel_;
    cField *Er_rel_;
    cField *Et_rel_;
    cField *Bl_rel_;
    cField *Br_rel_;
    cField *Bt_rel_;
    cField *Bl_rel_t_minus_halfdt_;
    cField *Br_rel_t_minus_halfdt_;
    cField *Bt_rel_t_minus_halfdt_;
    cField *Bl_rel_t_plus_halfdt_;
    cField *Br_rel_t_plus_halfdt_;
    cField *Bt_rel_t_plus_halfdt_;
    
    //PXR quantities:
    Field *Ex_pxr;
    Field *Ey_pxr;
    Field *Ez_pxr;
    Field *Bx_pxr;
    Field *By_pxr;
    Field *Bz_pxr;
    Field *Jx_pxr;
    Field *Jy_pxr;
    Field *Jz_pxr;
    Field *rho_pxr;
    Field *rhoold_pxr;
    
    //! Laser envelope
    LaserEnvelope *envelope;
    
    //! In linear polarization, this is the absolute value of the envelope of the vector potential in the polarization direction
    //! In circular polarization, this is the absolute value of the envelope of the vector potential along one direction
    //! This choice means that if the laser has a given a0 at waist in vacuum, 
    //! in the same point max(Env_A_abs_) = a0 in linear polarization and max(Env_A_abs_) = a0/sqrt(2) in circular polarization. 
    //! However in the same point the ponderomotive potential will always have max(Phi) = a0^2/2.
    Field *Env_A_abs_;
    
    //! Chi field (i.e. susceptibility) for envelope equation
    Field *Env_Chi_;
    
    //! Envelope of laser electric field along the same direction of Env_A_abs, absolute value
    Field *Env_E_abs_;
    //! Envelope of laser electric field along the longitudinal direction, absolute value
    Field *Env_Ex_abs_;
    
    FriedmanFields * filter_;
    
    //! all Fields in electromagn
    std::vector<Field *> allFields;
    
    //! all Fields averages required in diagnostic Fields
    std::vector<std::vector<Field *> > allFields_avg;
    
    //! Vector of charge density and currents for each species
    const unsigned int n_species;
    std::vector<Field *> Jx_s;
    std::vector<Field *> Jy_s;
    std::vector<Field *> Jz_s;
    std::vector<Field *> rho_s;
    
    // vector of susceptibility for each species
    std::vector<Field *> Env_Chi_s;
    
    //! Index where each species starts in allFields (concerns only species-related fields)
    std::vector<unsigned int> species_starts;
    
    //! Creates a new field with the right characteristics, depending on the name
    virtual Field *createField( std::string fieldname, Params& params ) = 0;
    
    //! nDim_field (from params)
    const unsigned int nDim_field;
    
    //! Volume of the single cell (from params)
    const double cell_volume;
    
    //! size_ (from params) always 3D
    std::vector<unsigned int> size_;
    
    //! Index of starting elements in arrays without duplicated borders
    //! By constuction 1 element is shared in primal field, 2 in dual
    //! 3 : Number of direction (=1, if dim not defined)
    //! 2 : isPrim/isDual
    unsigned int istart[3][2];
    //! Number of elements in arrays without duplicated borders
    unsigned int bufsize[3][2];
    
    //!\todo should this be just an integer???
    //! Oversize domain to exchange less particles (from params)
    std::vector<unsigned int> oversize;
    
    //! Constructor for Electromagn
    ElectroMagn( Params &params, Patch *patch );
    
    //! Method used to initialize the total charge currents and densities
    virtual void restartRhoJ();
    //! Method used to initialize the total charge currents and densities of species
    virtual void restartRhoJs();
    
    //! Method used to initialize the total susceptibility
    virtual void restartEnvChi();
    //! Method used to initialize the total susceptibility of species
    virtual void restartEnvChis();
    
    
    //! Method used to sum all species densities and currents to compute the total charge density and currents
    virtual void computeTotalRhoJ() = 0;
    
    //! Method used to sum all species susceptibility to compute the total susceptibility
    virtual void computeTotalEnvChi() = 0;
    
    virtual void initPoisson( Patch *patch ) = 0;
    virtual double compute_r() = 0;
    virtual void compute_Ap( Patch *patch ) = 0;
    virtual void compute_Ap_relativistic_Poisson( Patch *patch, double gamma_mean ) = 0;
    //Access to Ap
    virtual double compute_pAp() = 0;
    virtual void update_pand_r( double r_dot_r, double p_dot_Ap ) = 0;
    virtual void update_p( double rnew_dot_rnew, double r_dot_r ) = 0;
    virtual void initE( Patch *patch ) = 0;
    virtual void initE_relativistic_Poisson( Patch *patch, double gamma_mean ) = 0;
    virtual void initB_relativistic_Poisson( double gamma_mean ) = 0;
    virtual void center_fields_from_relativistic_Poisson() = 0; // centers in Yee cells the fields
    virtual void sum_rel_fields_to_em_fields() = 0;
    virtual void initRelativisticPoissonFields() = 0;
    virtual void centeringE( std::vector<double> E_Add ) = 0;
    virtual void centeringErel( std::vector<double> E_Add ) = 0;
    
    virtual double getEx_Xmin() = 0; // 2D !!!
    virtual double getEx_Xmax() = 0; // 2D !!!
    
    virtual double getExrel_Xmin() = 0; // 2D !!!
    virtual double getExrel_Xmax() = 0; // 2D !!!
    
    virtual double getEx_XminYmax() = 0; // 1D !!!
    virtual double getEy_XminYmax() = 0; // 1D !!!
    virtual double getEx_XmaxYmin() = 0; // 1D !!!
    virtual double getEy_XmaxYmin() = 0; // 1D !!!
    
    virtual double getExrel_XminYmax() = 0; // 1D !!!
    virtual double getEyrel_XminYmax() = 0; // 1D !!!
    virtual double getExrel_XmaxYmin() = 0; // 1D !!!
    virtual double getEyrel_XmaxYmin() = 0; // 1D !!!
    
    std::vector<unsigned int> index_min_p_;
    std::vector<unsigned int> index_max_p_;
    Field *phi_;
    Field *r_;
    Field *p_;
    Field *Ap_;

    cField *phi_AM_;
    cField *r_AM_;
    cField *p_AM_;
    cField *Ap_AM_;
    
    //! \todo check time_dual or time_prim (MG)
//    //! method used to solve Maxwell's equation (takes current time and time-step as input parameter)
//    virtual void solveMaxwellAmpere() = 0;
    //! Maxwell Ampere Solver
    Solver *MaxwellAmpereSolver_;
    //! Maxwell Faraday Solver
    Solver *MaxwellFaradaySolver_;
    virtual void saveMagneticFields( bool ) = 0;
    virtual void centerMagneticFields() = 0;
    virtual void binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes ) = 0;
    virtual void customFIRCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes, std::vector<double> filtering_coeff) = 0;
    
    void boundaryConditions( double time_dual, Patch *patch, SimWindow *simWindow );
    
    void laserDisabled();
    
    void incrementAvgField( Field *field, Field *field_avg );
    
    //! compute Poynting on borders
    virtual void computePoynting( unsigned int axis, unsigned int side ) = 0;
    
    //! pointing vector on borders
    //! 1D: poynting[0][0]=left , poynting[1][0]=right
    //! 2D: poynting[0][0]=xmin , poynting[1][0]=xmax
    //!     poynting[1][0]=ymin, poynting[1][0]=ymax
    std::vector<double> poynting[2];
    
    //same as above but instantaneous
    std::vector<double> poynting_inst[2];
    
    //! Compute local square norm of charge denisty is not null
    virtual double computeRhoNorm2()
    {
        return rho_->norm2( istart, bufsize );
    }
    
    //! Compute local sum of Ex
    inline double computeExSum()
    {
        return Ex_->sum( istart, bufsize );
    }
    //! Compute local sum of Ey
    inline double computeEySum()
    {
        return Ey_->sum( istart, bufsize );
    }
    //! Compute local sum of Ez
    inline double computeEzSum()
    {
        return Ez_->sum( istart, bufsize );
    }
    
    //! Compute local sum of Ex_rel_
    inline double computeExrelSum()
    {
        return Ex_rel_->sum( istart, bufsize );
    }
    //! Compute local sum of Ey_rel_
    inline double computeEyrelSum()
    {
        return Ey_rel_->sum( istart, bufsize );
    }
    //! Compute local sum of Ez_rel_
    inline double computeEzrelSum()
    {
        return Ez_rel_->sum( istart, bufsize );
    }
    
    //! external fields parameters the key string is the name of the field and the value is a vector of ExtField
    std::vector<ExtField> extFields;
    
    //! external time fields parameters the key string is the name of the field and the value is a vector of ExtField
    std::vector<PrescribedField> prescribedFields;
    
    //! Method used to impose external fields (apply to all Fields)
    virtual void applyExternalFields( Patch * );

    //! Method used to copy real field and impose an external time fields (apply to all Fields)
    virtual void applyPrescribedFields( Patch *, double time );

	//! Method use to reset the real value of all fields on which we imposed an external time field
	virtual void resetPrescribedFields();

    void saveExternalFields( Patch * );
    
    //! Method used to impose external fields (apply to a given Field)
    virtual void applyExternalField( Field *, Profile *, Patch * ) = 0 ;
    
    //! Method used to impose a prescribed fields (apply to a given Field)
    virtual void applyPrescribedField( Field *, Profile *, Patch *, double time) = 0 ;
    
    //! Antenna
    std::vector<Antenna> antennas;
    
    //! Method used to impose external currents (aka antennas)
    void applyAntenna( unsigned int iAntenna, double intensity );
    
    //! Method that fills the initial spatial profile of the antenna
    virtual void initAntennas( Patch *, Params& ) {};
    
    virtual double computeEnergy();
    
    inline int getMemFootPrint()
    {
    
        int emSize = 9+4; // 3 x (E, B, Bm) + 3 x J, rho
        
        if( Env_Chi_ ) {
            emSize += 4;    //Env_Chi, Env_A_abs, Env_E_abs, Env_Ex_abs;
        }
        
        for( unsigned int ispec=0 ; ispec<Jx_s.size() ; ispec++ ) {
            if( Jx_s [ispec] ) {
                emSize++;
            }
            if( Jy_s [ispec] ) {
                emSize++;
            }
            if( Jz_s [ispec] ) {
                emSize++;
            }
            if( rho_s [ispec] ) {
                emSize++;
            }
            if( Env_Chi_s [ispec] ) {
                emSize++;
            }
        }
        
        for( unsigned int idiag = 0 ; idiag < allFields_avg.size() ; idiag++ ) {
            emSize += allFields_avg[idiag].size() ;
        }
        
        
        for( size_t i=0 ; i<nDim_field ; i++ ) {
            emSize *= dimPrim[i];
        }
        
        emSize *= sizeof( double );
        return emSize;
    }
    
    //! Vector of boundary-condition per side for the fields
    std::vector<ElectroMagnBC *> emBoundCond;
    
    bool isXmin, isXmax;
    
    //! Corners coefficient for BC
    std::vector<double> beta_edge;
    std::vector<std::vector<double>> S_edge;
    
protected :
    bool is_pxr;
    
    //! Accumulate nrj lost / gained with moving window
    double nrj_mw_out;
    
    //! Accumulate nrj added with new fields
    double nrj_mw_inj;
    
private:

    
};

#endif
