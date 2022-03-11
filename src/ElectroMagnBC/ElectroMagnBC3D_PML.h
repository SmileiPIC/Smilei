#ifndef ELECTROMAGNBC3D_PML_H
#define ELECTROMAGNBC3D_PML_H

#include <vector>
#include "Tools.h"
#include "ElectroMagnBC3D.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"

class Params;
class ElectroMagn;
class Field;
class PML_Solver3D_Yee;
class PML_Solver3D_Bouchard;

class ElectroMagnBC3D_PML : public ElectroMagnBC3D
{
public:

    ElectroMagnBC3D_PML( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBC3D_PML();

    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;

    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;

    Field3D* Ex_ = NULL;
    Field3D* Ey_ = NULL;
    Field3D* Ez_ = NULL;
    Field3D* Bx_ = NULL;
    Field3D* By_ = NULL;
    Field3D* Bz_ = NULL;
    Field3D* Dx_ = NULL;
    Field3D* Dy_ = NULL;
    Field3D* Dz_ = NULL;
    Field3D* Hx_ = NULL;
    Field3D* Hy_ = NULL;
    Field3D* Hz_ = NULL;

    Field* getExPML() override { return Ex_; };
    Field* getEyPML() override { return Ey_; };
    Field* getEzPML() override { return Ez_; };
    Field* getBxPML() override { return Bx_; };
    Field* getByPML() override { return By_; };
    Field* getBzPML() override { return Bz_; };
    Field* getDxPML() override { return Dx_; };
    Field* getDyPML() override { return Dy_; };
    Field* getDzPML() override { return Dz_; };
    Field* getHxPML() override { return Hx_; };
    Field* getHyPML() override { return Hy_; };
    Field* getHzPML() override { return Hz_; };

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

    double factor_laser_space_time ;
    double factor_laser_angle_W ;
    double factor_laser_angle_E ;
    double factor_laser_angle_S ;
    double factor_laser_angle_N ;
    double factor_laser_angle_B ;
    double factor_laser_angle_T ;
};

#endif
