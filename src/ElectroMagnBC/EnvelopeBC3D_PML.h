#ifndef ENVELOPEBC3D_PML_H
#define ENVELOPEBC3D_PML_H


#include <vector>
#include <complex>
#include "Tools.h"
#include "EnvelopeBC.h"
#include "LaserEnvelope.h"
#include "Field3D.h"
#include "cField3D.h"
#include "ElectroMagn.h"

class Params;
class ElectroMagn;
class LaserEnvelope;
class PML_Solver3D_Envelope;

class EnvelopeBC3D_PML : public EnvelopeBC
{
public:

    EnvelopeBC3D_PML( Params &params, Patch *patch, unsigned int i_boundary );
    ~EnvelopeBC3D_PML();

    virtual void apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch ) override;

    //void save_fields( Field *, Patch *patch ) override;
    //void disableExternalFields() override;

    cField3D* A_np1_ = NULL;
    cField3D* A_n_ = NULL;
    cField3D* A_nm1_ = NULL;

    cField3D* u1_np1_x_ = NULL;
    cField3D* u2_np1_x_ = NULL;
    cField3D* u3_np1_x_ = NULL;
    cField3D* u1_nm1_x_ = NULL;
    cField3D* u2_nm1_x_ = NULL;
    cField3D* u3_nm1_x_ = NULL;

    cField3D* u1_np1_y_ = NULL;
    cField3D* u2_np1_y_ = NULL;
    cField3D* u3_np1_y_ = NULL;
    cField3D* u1_nm1_y_ = NULL;
    cField3D* u2_nm1_y_ = NULL;
    cField3D* u3_nm1_y_ = NULL;

    cField3D* u1_np1_z_ = NULL;
    cField3D* u2_np1_z_ = NULL;
    cField3D* u3_np1_z_ = NULL;
    cField3D* u1_nm1_z_ = NULL;
    cField3D* u2_nm1_z_ = NULL;
    cField3D* u3_nm1_z_ = NULL;

    Field3D* Phi_ = NULL;
    Field3D* Chi_ = NULL;

    cField* getAnp1PML() override { return A_np1_; };
    cField* getAnPML() override { return A_n_; };
    cField* getAnm1PML() override { return A_nm1_; };

    cField* getu1np1xPML() override { return u1_np1_x_; };
    cField* getu2np1xPML() override { return u2_np1_x_; };
    cField* getu3np1xPML() override { return u3_np1_x_; };
    cField* getu1nm1xPML() override { return u1_nm1_x_; };
    cField* getu2nm1xPML() override { return u2_nm1_x_; };
    cField* getu3nm1xPML() override { return u3_nm1_x_; };

    cField* getu1np1yPML() override { return u1_np1_y_; };
    cField* getu2np1yPML() override { return u2_np1_y_; };
    cField* getu3np1yPML() override { return u3_np1_y_; };
    cField* getu1nm1yPML() override { return u1_nm1_y_; };
    cField* getu2nm1yPML() override { return u2_nm1_y_; };
    cField* getu3nm1yPML() override { return u3_nm1_y_; };

    cField* getu1np1zPML() override { return u1_np1_z_; };
    cField* getu2np1zPML() override { return u2_np1_z_; };
    cField* getu3np1zPML() override { return u3_np1_z_; };
    cField* getu1nm1zPML() override { return u1_nm1_z_; };
    cField* getu2nm1zPML() override { return u2_nm1_z_; };
    cField* getu3nm1zPML() override { return u3_nm1_z_; };

    Field* getPhiPML() override { return Phi_; };
    Field* getChiPML() override { return Chi_; };

    int domain_oversize_x;
    int domain_oversize_y;
    int domain_oversize_z;

    int ncells_pml_xmin ;
    int ncells_pml_xmax ;
    int ncells_pml_ymin ;
    int ncells_pml_ymax ;
    int ncells_pml_zmin ;
    int ncells_pml_zmax ;

    int ncells_pml_domain_xmin;
    int ncells_pml_domain_xmax;
    int ncells_pml_domain_ymin;
    int ncells_pml_domain_ymax;
    int ncells_pml_domain_zmin;
    int ncells_pml_domain_zmax;

    int ncells_pml;
    int ncells_pml_domain;

    int nsolver;

    int min2exchange;
    int max2exchange;
    int solvermin;
    int solvermax;

    int ypml_size_in_x;
    int zpml_size_in_x;
    int zpml_size_in_y;
    int startpml;
};

#endif
