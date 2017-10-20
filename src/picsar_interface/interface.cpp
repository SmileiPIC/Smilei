#include "interface.h"



void copy_field(Field3D* out, Field3D * in)
{
        unsigned int n1,n2,n3;
        n1 = ((in->dims_[0]) < (out->dims_[0]) ? (in->dims_[0]) : (out->dims_[0]));
        n2 = ((in->dims_[1]) < (out->dims_[1]) ? (in->dims_[1]) : (out->dims_[1]));
        n3 = ((in->dims_[2]) < (out->dims_[2]) ? (in->dims_[2]) : (out->dims_[2]));
//	for(unsigned int i =0;i<(out->dims_[0]);i++){
//	for(unsigned int j =0;j<(out->dims_[1]);j++){
//	for(unsigned int k =0;k<(out->dims_[2]);k++)	
//		(*out)(i,j,k) = 0.0;
//		}
//
//	}
        for(unsigned int i=0;i<n1;i++){
                for(unsigned int j=0;j<n2;j++){
                        for(unsigned int k=0;k<n3;k++){
                                (*out)(i,j,k) =  (*in)(i,j,k);
                        }
                }
        }
}
void duplicate_field_into_pxr(ElectroMagn * fields,Params& params)
{
Field3D* Ex3D_pxr = static_cast<Field3D*>(fields->Ex_pxr);
Field3D* Ey3D_pxr = static_cast<Field3D*>(fields->Ey_pxr);
Field3D* Ez3D_pxr = static_cast<Field3D*>(fields->Ez_pxr);
Field3D* Bx3D_pxr = static_cast<Field3D*>(fields->Bx_pxr);
Field3D* By3D_pxr = static_cast<Field3D*>(fields->By_pxr);
Field3D* Bz3D_pxr = static_cast<Field3D*>(fields->Bz_pxr);
Field3D* Jx3D_pxr = static_cast<Field3D*>(fields->Jx_pxr);
Field3D* Jy3D_pxr = static_cast<Field3D*>(fields->Jy_pxr);
Field3D* Jz3D_pxr = static_cast<Field3D*>(fields->Jz_pxr);
Field3D* rho3D_pxr = static_cast<Field3D*>(fields->rho_pxr);
Field3D* rhoold3D_pxr = static_cast<Field3D*>(fields->rhoold_pxr);

Field3D* Ex3D = static_cast<Field3D*>(fields->Ex_);
Field3D* Ey3D = static_cast<Field3D*>(fields->Ey_);
Field3D* Ez3D = static_cast<Field3D*>(fields->Ez_);
Field3D* Bx3D = static_cast<Field3D*>(fields->Bx_);
Field3D* By3D = static_cast<Field3D*>(fields->By_);
Field3D* Bz3D = static_cast<Field3D*>(fields->Bz_);
Field3D* Jx3D = static_cast<Field3D*>(fields->Jx_);
Field3D* Jy3D = static_cast<Field3D*>(fields->Jy_);
Field3D* Jz3D = static_cast<Field3D*>(fields->Jz_);
Field3D* rho3D = static_cast<Field3D*>(fields->rho_);
Field3D* rhoold3D = static_cast<Field3D*>(fields->rhoold_);
copy_field(Ex3D_pxr,Ex3D);
copy_field(Ey3D_pxr,Ey3D);
copy_field(Ez3D_pxr,Ez3D);
copy_field(Bx3D_pxr,Bx3D);
copy_field(By3D_pxr,By3D);
copy_field(Bz3D_pxr,Bz3D);
copy_field(Jx3D_pxr,Jx3D);
copy_field(Jy3D_pxr,Jy3D);
copy_field(Jz3D_pxr,Jz3D);
copy_field(rho3D_pxr,rho3D);
copy_field(rhoold3D_pxr,rhoold3D);
}

void duplicate_field_into_smilei(ElectroMagn * fields,Params& params)
{

Field3D* Ex3D_pxr = static_cast<Field3D*>(fields->Ex_pxr);
Field3D* Ey3D_pxr = static_cast<Field3D*>(fields->Ey_pxr);
Field3D* Ez3D_pxr = static_cast<Field3D*>(fields->Ez_pxr);
Field3D* Bx3D_pxr = static_cast<Field3D*>(fields->Bx_pxr);
Field3D* By3D_pxr = static_cast<Field3D*>(fields->By_pxr);
Field3D* Bz3D_pxr = static_cast<Field3D*>(fields->Bz_pxr);

//Field3D* Jx3D_pxr = static_cast<Field3D*>(fields->Jx_pxr);
//Field3D* Jy3D_pxr = static_cast<Field3D*>(fields->Jy_pxr);
//Field3D* Jz3D_pxr = static_cast<Field3D*>(fields->Jz_pxr);
//Field3D* rho3D_pxr = static_cast<Field3D*>(fields->rho_pxr);
//Field3D* rhoold3D_pxr = static_cast<Field3D*>(fields->rhoold_pxr);

Field3D* Ex3D = static_cast<Field3D*>(fields->Ex_);
Field3D* Ey3D = static_cast<Field3D*>(fields->Ey_);
Field3D* Ez3D = static_cast<Field3D*>(fields->Ez_);
Field3D* Bx3D = static_cast<Field3D*>(fields->Bx_);
Field3D* By3D = static_cast<Field3D*>(fields->By_);
Field3D* Bz3D = static_cast<Field3D*>(fields->Bz_);

//Field3D* Jx3D = static_cast<Field3D*>(fields->Jx_);
//Field3D* Jy3D = static_cast<Field3D*>(fields->Jy_);
//Field3D* Jz3D = static_cast<Field3D*>(fields->Jz_);
//Field3D* rho3D = static_cast<Field3D*>(fields->rho_);
//Field3D* rhoold3D = static_cast<Field3D*>(fields->rhoold_);

copy_field(Ex3D,Ex3D_pxr);
copy_field(Ey3D,Ey3D_pxr);
copy_field(Ez3D,Ez3D_pxr);
copy_field(Bx3D,Bx3D_pxr);
copy_field(By3D,By3D_pxr);
copy_field(Bz3D,Bz3D_pxr);

//copy_field(Jx3D,Jx3D_pxr);
//copy_field(Jy3D,Jy3D_pxr);
//copy_field(Jz3D,Jz3D_pxr);
//copy_field(rho3D,rho3D_pxr);
//copy_field(rhoold3D,rhoold3D_pxr);
}

