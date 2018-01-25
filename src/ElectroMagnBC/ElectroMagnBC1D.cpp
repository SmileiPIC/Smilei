
#include "ElectroMagnBC1D.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field1D.h"


ElectroMagnBC1D::ElectroMagnBC1D( Params &params, Patch* patch, unsigned int _min_max )
  : ElectroMagnBC( params, patch, _min_max )
{
    
}

ElectroMagnBC1D::~ElectroMagnBC1D()
{
}


void ElectroMagnBC1D::applyBConEdges( ElectroMagn* EMfields, Patch* patch )
{
    // Static cast of the fields
   Field1D* Ex1D = static_cast<Field1D*>(EMfields->Ex_);
   Field1D* Ey1D = static_cast<Field1D*>(EMfields->Ey_);
   Field1D* Ez1D = static_cast<Field1D*>(EMfields->Ez_);
   Field1D* Bx1D = static_cast<Field1D*>(EMfields->Bx_);
   Field1D* By1D = static_cast<Field1D*>(EMfields->By_);
   Field1D* Bz1D = static_cast<Field1D*>(EMfields->Bz_);

   double one_ov_dbeta ;
  
   // Do we chose to let extremal points here ?
   if (patch->isXmin()){

       unsigned int i = 0;

   } //End series of Xmin edges

   if (patch->isXmax()){

       unsigned int i = nx_p - 1;

   } //End series of Xmax edges

}
