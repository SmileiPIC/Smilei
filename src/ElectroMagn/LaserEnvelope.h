
#ifndef LaserENVELOPE_H
#define LaserENVELOPE_H

#include <vector>
#include <complex>

#include "Params.h"
#include "Patch.h"

class Field;
class ElectroMagn;
class Profile;
class EnvelopeBC;
class SimWindow;

// Class for envelope
class LaserEnvelope
{
public:
    LaserEnvelope( Params &params, Patch *patch, ElectroMagn *EMfields ); // Main constructor
    LaserEnvelope( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved ); // Cloning constructor
    virtual void initEnvelope( Patch *patch, ElectroMagn *EMfields ) = 0;
    virtual ~LaserEnvelope();
    virtual void updateEnvelope( ElectroMagn *EMfields ) = 0;
    virtual void updateEnvelopeReducedDispersion( ElectroMagn *EMfields ) = 0;
    virtual void computePhiEnvAEnvE( ElectroMagn *EMfields ) = 0;
    virtual void computeGradientPhi( ElectroMagn *EMfields ) = 0;
    void boundaryConditions( int itime, double time_dual, Patch *patch, Params &params, SimWindow *simWindow, ElectroMagn *EMfields );
    virtual void savePhiAndGradPhi() = 0;
    virtual void centerPhiAndGradPhi() = 0;
    
    Profile *profile_;
    const std::vector<double> cell_length;
    const double timestep;

    double polarization_phi;
    double ellipticity; // 0 for linear polarization, 1 for circular polarization
    double ellipticity_factor; // 1 for linear polarization, 2 for circular polarization. 
    // This coefficient is used for the ponderomotive potential Phi = ellipticity_factor*|A|^2/2.

    std:: string envelope_solver  = "explicit"; // default value
    
    Field *A_;         // envelope value at timestep n
    Field *A0_;        // envelope value at timestep n-1
    
    Field *Phi_;       // ponderomotive potential Phi=|A|^2/2 value at timestep n
    Field *GradPhix_;  // x component of the gradient of Phi at timestep n
    Field *GradPhiy_;  // y component of the gradient of Phi at timestep n
    Field *GradPhiz_;  // z component of the gradient of Phi at timestep n
    Field *GradPhil_;  // l component of the gradient of Phi at timestep n
    Field *GradPhir_;  // r component of the gradient of Phi at timestep n
    
    // Correspondent time-centered quantities
    Field *Phi_m;
    Field *GradPhix_m;
    Field *GradPhiy_m;
    Field *GradPhiz_m;
    Field *GradPhil_m;
    Field *GradPhir_m;
    
    //! Vector of boundary-condition per side for the envelope field
    std::vector<EnvelopeBC *> EnvBoundCond;
    //EnvBoundCond = EnvelopeBC_Factory::create(params, patch);

    ////// auxiliary quantities for the solver, computation of gradients etc
    double one_ov_2dt,dt_sq;         // 1/(2dt), dt^2
    double one_ov_dx_sq,one_ov_2dx;  // 1/(2dx), 1/(dx^2)     // for all geometries
    double one_ov_dy_sq,one_ov_2dy;  // 1/(2dy), 1/(dy^2)     // for 2D, 3D geometries
    double one_ov_dz_sq,one_ov_2dz;  // 1/(2dz), 1/(dz^2)     // for 3D geometry
    double one_ov_dl_sq,one_ov_2dl;  // 1/(2dl), 1/(dl^2)     // for AMcylindrical geometry
    double one_ov_dr_sq,one_ov_2dr;  // 1/(2dr), 1/(dr^2), dr // for AMcylindrical geometry
    double dr;

    // imaginary unit and quantities using it
    std::complex<double> i1 = std::complex<double>( 0., 1 );
    std::complex<double> i1_2k0_over_2dx;
    std::complex<double> i1_2k0_over_2dl;
    std::complex<double> one_plus_ik0dt;
    std::complex<double> one_plus_ik0dt_ov_one_plus_k0sq_dtsq;

    // coefficient necessary for the reduced dispersion envelope solver  
    double delta; 
};

// Class for envelope
class LaserEnvelope1D : public LaserEnvelope
{
public:
    LaserEnvelope1D( Params &params, Patch *patch, ElectroMagn *EMfields );
    LaserEnvelope1D( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved );
    void initEnvelope( Patch *patch, ElectroMagn *EMfields ) override final;
    ~LaserEnvelope1D();
    void updateEnvelope( ElectroMagn *EMfields ) override final;
    void updateEnvelopeReducedDispersion( ElectroMagn *EMfields ) override final;
    void computePhiEnvAEnvE( ElectroMagn *EMfields ) override final;
    void computeGradientPhi( ElectroMagn *EMfields ) override final;
    void savePhiAndGradPhi() override final;
    void centerPhiAndGradPhi() override final;
};

// Class for envelope
class LaserEnvelope2D : public LaserEnvelope
{
public:
    LaserEnvelope2D( Params &params, Patch *patch, ElectroMagn *EMfields );
    LaserEnvelope2D( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved );
    void initEnvelope( Patch *patch, ElectroMagn *EMfields ) override final;
    ~LaserEnvelope2D();
    void updateEnvelope( ElectroMagn *EMfields ) override final;
    void updateEnvelopeReducedDispersion( ElectroMagn *EMfields ) override final;
    void computePhiEnvAEnvE( ElectroMagn *EMfields ) override final;
    void computeGradientPhi( ElectroMagn *EMfields ) override final;
    void savePhiAndGradPhi() override final;
    void centerPhiAndGradPhi() override final;
};

// Class for envelope
class LaserEnvelope3D : public LaserEnvelope
{
public:
    LaserEnvelope3D( Params &params, Patch *patch, ElectroMagn *EMfields );
    LaserEnvelope3D( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved );
    void initEnvelope( Patch *patch, ElectroMagn *EMfields ) override final;
    ~LaserEnvelope3D();
    void updateEnvelope( ElectroMagn *EMfields ) override final;
    void updateEnvelopeReducedDispersion( ElectroMagn *EMfields ) override final;
    void computePhiEnvAEnvE( ElectroMagn *EMfields ) override final;
    void computeGradientPhi( ElectroMagn *EMfields ) override final;
    void savePhiAndGradPhi() override final;
    void centerPhiAndGradPhi() override final;
};

// Class for envelope with cylindrical symmetry
class LaserEnvelopeAM : public LaserEnvelope
{
public:
    LaserEnvelopeAM( Params &params, Patch *patch, ElectroMagn *EMfields );
    LaserEnvelopeAM( LaserEnvelope *envelope, Patch *patch, ElectroMagn *EMfields, Params &params, unsigned int n_moved );
    void initEnvelope( Patch *patch, ElectroMagn *EMfields ) override final;
    ~LaserEnvelopeAM();
    void updateEnvelope( ElectroMagn *EMfields ) override final;
    void updateEnvelopeReducedDispersion( ElectroMagn *EMfields ) override final;
    void computePhiEnvAEnvE( ElectroMagn *EMfields ) override final;
    void computeGradientPhi( ElectroMagn *EMfields ) override final;
    void savePhiAndGradPhi() override final;
    void centerPhiAndGradPhi() override final;
};



#endif
