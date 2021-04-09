
#include "ElectroMagnBC2D.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field2D.h"


ElectroMagnBC2D::ElectroMagnBC2D( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC( params, patch, i_boundary )
{
}

ElectroMagnBC2D::~ElectroMagnBC2D()
{
}


void ElectroMagnBC2D::applyBConEdges( ElectroMagn *EMfields, Patch *patch )
{
//    // Static cast of the fields
//   Field2D* Ex2D = static_cast<Field2D*>(EMfields->Ex_);
//   Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
//   Field2D* Ez2D = static_cast<Field2D*>(EMfields->Ez_);
//   Field2D* Bx2D = static_cast<Field2D*>(EMfields->Bx_);
//   Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
//   Field2D* Bz2D = static_cast<Field2D*>(EMfields->Bz_);
//
//   double one_ov_dbeta ;
//
//   if (patch->isXmin()){
//
//       unsigned int i = 0;
//       if (patch->isYmin()){
//           // Xmin/Ymin
//           // edge 0 : By[0,0,k] + beta(-x)Bx[0,0,k] = S(-x)
//           // edge 8 : Bx[0,0,k] + beta(-y)By[0,0,k] = S(-y)
//       }//End Xmin Ymin edge
//
//       if (patch->isYmax()){
//           // Xmin/Ymax
//           //edge 1 : By[0,n_p[1]-1,k] + beta(-x)Bx[0,n_p[1],k] = S(-x)
//           //edge 12 : Bx[0,n_p[1],k] + beta(-y)By[0,n_p[1]-1,k] = S(-y)
//       }// End Xmin Ymax edge
//
//   } //End series of Xmin edges
//
//   if (patch->isXmax()){
//
//       unsigned int i = n_p[0] - 1;
//       if (patch->isYmin()){
//           // Xmax/Ymin
//           // edge 4 : By[n_p[0],0,k] + beta(+x)Bx[n_p[0]-1,0,k] = S(+x)
//           // edge 9 : Bx[n_p[0]-1,0,k] + beta(-y)By[n_p[0]-1,0,k] = S(-y)
//       }//End Xmax Ymin edge
//
//       if (patch->isYmax()){
//           // Xmax/Ymax
//           //edge 5 :  By[n_p[0],n_p[1]-1,k] + beta(+x)Bx[n_p[0]-1,n_p[1],k] = S(+x)
//           //edge 13 : Bx[n_p[0]-1,n_p[1],k] + beta(-y)By[n_p[0],n_p[1]-1,k] = S(-y)
//       }// End Xmax Ymax edge
//
//   } //End series of Xmax edges
//
}
