
#ifndef FIELDSBC2D_Trans_DAMPING_H
#define FIELDSBC2D_Trans_DAMPING_H

#include "FieldsBC.h" 

class PicParams;
class ElectroMagn;

class FieldsBC2D_Trans_Damping : public FieldsBC {
public:
    FieldsBC2D_Trans_Damping( PicParams &params );
    ~FieldsBC2D_Trans_Damping();

    virtual void apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi);

 private:
    //! Number of nodes on the primal grid in the x-direction
    unsigned int nx_p;

    //! Number of nodes on the dual grid in the x-direction
    unsigned int nx_d;

    //! Number of nodes on the primal grid in the y-direction
    unsigned int ny_p;

    //! Number of nodes on the dual grid in the y-direction
    unsigned int ny_d;

    
    // number of dumping layers
    int ny_l;
    // Damping coefficient
    double cdamp; 
    // array of coefficient per layer
    double* coeff;

    
};

#endif

