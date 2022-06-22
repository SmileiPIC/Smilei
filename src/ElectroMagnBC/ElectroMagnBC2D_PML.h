
#ifndef ELECTROMAGNBC2D_PML_H
#define ELECTROMAGNBC2D_PML_H


#include <vector>
#include "Tools.h"
#include "ElectroMagnBC2D.h"
#include "ElectroMagn2D.h"
#include "Field2D.h"

class Params;
class ElectroMagn;
class Field;
class PML_Solver2D_Yee;
class PML_Solver2D_Bouchard;

class ElectroMagnBC2D_PML : public ElectroMagnBC2D
{
public:

    ElectroMagnBC2D_PML( Params &params, Patch *patch, unsigned int _min_max );
    ~ElectroMagnBC2D_PML();
    
    virtual void apply( ElectroMagn *EMfields, double time_dual, Patch *patch ) override;
    
    void save_fields( Field *, Patch *patch ) override;
    void disableExternalFields() override;

    Field2D* Ex_ = NULL;
    Field2D* Ey_ = NULL;
    Field2D* Ez_ = NULL;
    Field2D* Bx_ = NULL;
    Field2D* By_ = NULL;
    Field2D* Bz_ = NULL;
    Field2D* Dx_ = NULL;
    Field2D* Dy_ = NULL;
    Field2D* Dz_ = NULL;
    Field2D* Hx_ = NULL;
    Field2D* Hy_ = NULL;
    Field2D* Hz_ = NULL;

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

    int ncells_pml_xmin ;
    int ncells_pml_xmax ;
    int ncells_pml_ymin ;
    int ncells_pml_ymax ;
        
    int ncells_pml_domain_xmin;
    int ncells_pml_domain_xmax;
    int ncells_pml_domain_ymin;
    int ncells_pml_domain_ymax; 

    int ncells_pml;
    int ncells_pml_domain;

    int nsolver;

    int min2exchange;
    int max2exchange;
    unsigned int solvermin;
    unsigned int solvermax;   

    int ypml_size_in_x;
    int startpml;

    double factor_laser_space_time ;
    double factor_laser_angle_W ;
    double factor_laser_angle_E ;
    double factor_laser_angle_S ;
    double factor_laser_angle_N ;
    
};

#endif

