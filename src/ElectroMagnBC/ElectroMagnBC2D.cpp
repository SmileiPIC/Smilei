
#include "ElectroMagnBC2D.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field2D.h"


ElectroMagnBC2D::ElectroMagnBC2D( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC( params, patch, _min_max )
{
    std::vector<unsigned int> n_space(params.n_space);
    if (params.uncoupled_grids)
        n_space = params.n_space_region;

    std::vector<unsigned int> oversize(params.oversize);
    if (params.uncoupled_grids)
        oversize = params.region_oversize;

    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = nx_p+1-params.is_pxr;
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = ny_p+1-params.is_pxr;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = params.cell_length[0];
    dt_ov_dx = dt/dx;
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = params.cell_length[1];
    dt_ov_dy = dt/dy;
    dy_ov_dt = 1.0/dt_ov_dy;
    
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
//           //edge 1 : By[0,ny_p-1,k] + beta(-x)Bx[0,ny_p,k] = S(-x)
//           //edge 12 : Bx[0,ny_p,k] + beta(-y)By[0,ny_p-1,k] = S(-y)
//       }// End Xmin Ymax edge
//
//   } //End series of Xmin edges
//
//   if (patch->isXmax()){
//
//       unsigned int i = nx_p - 1;
//       if (patch->isYmin()){
//           // Xmax/Ymin
//           // edge 4 : By[nx_p,0,k] + beta(+x)Bx[nx_p-1,0,k] = S(+x)
//           // edge 9 : Bx[nx_p-1,0,k] + beta(-y)By[nx_p-1,0,k] = S(-y)
//       }//End Xmax Ymin edge
//
//       if (patch->isYmax()){
//           // Xmax/Ymax
//           //edge 5 :  By[nx_p,ny_p-1,k] + beta(+x)Bx[nx_p-1,ny_p,k] = S(+x)
//           //edge 13 : Bx[nx_p-1,ny_p,k] + beta(-y)By[nx_p,ny_p-1,k] = S(-y)
//       }// End Xmax Ymax edge
//
//   } //End series of Xmax edges
//
}
