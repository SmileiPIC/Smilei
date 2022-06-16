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
class PML_SolverAM_EnvelopeReducedDispersion;

class EnvelopeBCAM_PML : public EnvelopeBC
{
public:

    EnvelopeBCAM_PML( Params &params, Patch *patch, unsigned int i_boundary );
    ~EnvelopeBCAM_PML();

    virtual void apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch ) override;

    //void save_fields( Field *, Patch *patch ) override;
    //void disableExternalFields() override;

    cField2D* A_np1_ = NULL;
    cField2D* A_n_ = NULL;
    cField2D* A_nm1_ = NULL;

    cField2D* G_np1_ = NULL;
    cField2D* G_n_ = NULL;
    cField2D* G_nm1_ = NULL;

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

    Field2D* Phi_ = NULL;
    Field2D* Chi_ = NULL;

    cField* getAnp1PML() override { return A_np1_; };
    cField* getAnPML() override { return A_n_; };
    cField* getAnm1PML() override { return A_nm1_; };

    cField* getGnp1PML() override { return G_np1_; };
    cField* getGnPML() override { return G_n_; };
    cField* getGnm1PML() override { return G_nm1_; };

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

    Field* getPhiPML() override { return Phi_; };
    Field* getChiPML() override { return Chi_; };

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
