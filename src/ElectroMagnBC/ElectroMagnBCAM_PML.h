#ifndef ELECTROMAGNBCAM_PML_H
#define ELECTROMAGNBCAM_PML_H


#include <vector>
#include <complex>
#include "Tools.h"
#include "ElectroMagnBCAM.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "dcomplex.h"

class Params;
class ElectroMagn;
class Field;
class PML_SolverAM;

class ElectroMagnBCAM_PML : public ElectroMagnBCAM
{
public:

    ElectroMagnBCAM_PML( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBCAM_PML();

    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;

    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;

    std::vector<cField2D *> El_ ;//= NULL;
    std::vector<cField2D *> Er_ ;//= NULL;
    std::vector<cField2D *> Et_ ;//= NULL;
    std::vector<cField2D *> Bl_ ;//= NULL;
    std::vector<cField2D *> Br_ ;//= NULL;
    std::vector<cField2D *> Bt_ ;//= NULL;
    std::vector<cField2D *> Dl_ ;//= NULL;
    std::vector<cField2D *> Dr_ ;//= NULL;
    std::vector<cField2D *> Dt_ ;//= NULL;
    std::vector<cField2D *> Hl_ ;//= NULL;
    std::vector<cField2D *> Hr_ ;//= NULL;
    std::vector<cField2D *> Ht_ ;//= NULL;

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

    double factor_laser_space_time ;
    double factor_laser_angle_W ;
    double factor_laser_angle_E ;

};

#endif
