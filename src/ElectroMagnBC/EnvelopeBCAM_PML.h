#ifndef ENVELOPEBCAM_PML_H
#define ENVELOPEBCAM_PML_H


#include <vector>
#include <complex>
#include "Tools.h"
#include "EnvelopeBC.h"
#include "LaserEnvelope.h"
#include "Field2D.h"
#include "cField2D.h"
#include "ElectroMagn.h"

class Params;
class ElectroMagn;
class LaserEnvelope;
class PML_SolverAM_Envelope;

class EnvelopeBCAM_PML : public EnvelopeBC
{
public:

    EnvelopeBCAM_PML( Params &params, Patch *patch, unsigned int i_boundary );
    ~EnvelopeBCAM_PML();

    virtual void apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch ) override;

    //void save_fields( Field *, Patch *patch ) override;
    //void disableExternalFields() override;

    cField2D* A2D_np1_ = NULL;
    cField2D* A2D_n_ = NULL;
    cField2D* A2D_nm1_ = NULL;

    cField2D* G2D_np1_ = NULL;
    cField2D* G2D_n_ = NULL;
    cField2D* G2D_nm1_ = NULL;

    cField2D* u1_np1_l_ = NULL;
    cField2D* u2_np1_l_ = NULL;
    cField2D* u3_np1_l_ = NULL;
    cField2D* u1_nm1_l_ = NULL;
    cField2D* u2_nm1_l_ = NULL;
    cField2D* u3_nm1_l_ = NULL;

    cField2D* u1_np1_r_ = NULL;
    cField2D* u2_np1_r_ = NULL;
    cField2D* u3_np1_r_ = NULL;
    cField2D* u1_nm1_r_ = NULL;
    cField2D* u2_nm1_r_ = NULL;
    cField2D* u3_nm1_r_ = NULL;

    Field2D* Phi2D_ = NULL;

    cField* getA2Dnp1PML() override { return A2D_np1_; };
    cField* getA2DnPML() override { return A2D_n_; };
    cField* getA2Dnm1PML() override { return A2D_nm1_; };

    cField* getG2Dnp1PML() override { return A2D_np1_; };
    cField* getG2DnPML() override { return A2D_n_; };
    cField* getG2Dnm1PML() override { return A2D_nm1_; };

    cField* getu1np1lPML() override { return u1_np1_l_; };
    cField* getu2np1lPML() override { return u2_np1_l_; };
    cField* getu3np1lPML() override { return u3_np1_l_; };
    cField* getu1nm1lPML() override { return u1_nm1_l_; };
    cField* getu2nm1lPML() override { return u2_nm1_l_; };
    cField* getu3nm1lPML() override { return u3_nm1_l_; };

    cField* getu1np1rPML() override { return u1_np1_r_; };
    cField* getu2np1rPML() override { return u2_np1_r_; };
    cField* getu3np1rPML() override { return u3_np1_r_; };
    cField* getu1nm1rPML() override { return u1_nm1_r_; };
    cField* getu2nm1rPML() override { return u2_nm1_r_; };
    cField* getu3nm1rPML() override { return u3_nm1_r_; };

    Field* getPhi2DPML() override { return Phi2D_; };

    int domain_oversize_l;
    int domain_oversize_r;

    int ncells_pml_lmin ;
    int ncells_pml_lmax ;
    int ncells_pml_rmin ;
    int ncells_pml_rmax ;

    int ncells_pml_domain_lmin;
    int ncells_pml_domain_lmax;
    int ncells_pml_domain_rmin;
    int ncells_pml_domain_rmax;

    int ncells_pml;
    int ncells_pml_domain;

    int nsolver;

    int min2exchange;
    int max2exchange;
    int solvermin;
    int solvermax;

    int rpml_size_in_l;
    int startpml;

    int j_glob_pml;
};

#endif
