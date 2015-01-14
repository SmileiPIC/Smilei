#include "ElectroMagn2D.h"

#include <cmath>

#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "Field2D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart2D.h"

#include "ExtFieldProfile2D.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::ElectroMagn2D(PicParams &params, LaserParams &laser_params, SmileiMPI* smpi) : 
ElectroMagn(params, laser_params, smpi),
isWestern(smpi->isWestern()),
isEastern(smpi->isEastern()),
isNorthern(smpi->isNorthern()),
isSouthern(smpi->isSouthern())
{
    // local dt to store
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    int process_coord_x = smpi2D->getProcCoord(0);
    
    
    // --------------------------------------------------
    // Calculate quantities related to the simulation box
    // --------------------------------------------------
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = cell_length[0];
    dt_ov_dx = timestep/dx;
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = cell_length[1];
    dt_ov_dy = timestep/dy;
    dy_ov_dt = 1.0/dt_ov_dy;
    
    // ----------------------
    // Electromagnetic fields
    // ----------------------
    //! \todo Homogenize 1D/2D dimPrim/dimDual or nx_p/nx_d/ny_p/ny_d
    
    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    
    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = n_space[0]+2+2*oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = n_space[1]+2+2*oversize[1];
    
    // Allocation of the EM fields
    Ex_  = new Field2D(dimPrim, 0, false, "Ex");
    Ey_  = new Field2D(dimPrim, 1, false, "Ey");
    Ez_  = new Field2D(dimPrim, 2, false, "Ez");
    Bx_  = new Field2D(dimPrim, 0, true,  "Bx");
    By_  = new Field2D(dimPrim, 1, true,  "By");
    Bz_  = new Field2D(dimPrim, 2, true,  "Bz");
    Bx_m = new Field2D(dimPrim, 0, true,  "Bx_m");
    By_m = new Field2D(dimPrim, 1, true,  "By_m");
    Bz_m = new Field2D(dimPrim, 2, true,  "Bz_m");
    
    // Total charge currents and densities
    Jx_   = new Field2D(dimPrim, 0, false, "Jx");
    Jy_   = new Field2D(dimPrim, 1, false, "Jy");
    Jz_   = new Field2D(dimPrim, 2, false, "Jz");
    rho_  = new Field2D(dimPrim, "Rho" );
    
    // Allocation of the time-averaged EM fields
    Ex_avg  = new Field2D(dimPrim, 0, false, "Ex_avg");
    Ey_avg  = new Field2D(dimPrim, 1, false, "Ey_avg");
    Ez_avg  = new Field2D(dimPrim, 2, false, "Ez_avg");
    Bx_avg  = new Field2D(dimPrim, 0, true,  "Bx_avg");
    By_avg  = new Field2D(dimPrim, 1, true,  "By_avg");
    Bz_avg  = new Field2D(dimPrim, 2, true,  "Bz_avg");
    
    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field2D(dimPrim, 0, false, ("Jx_"+params.species_param[ispec].species_type).c_str());
        Jy_s[ispec]  = new Field2D(dimPrim, 1, false, ("Jy_"+params.species_param[ispec].species_type).c_str());
        Jz_s[ispec]  = new Field2D(dimPrim, 2, false, ("Jz_"+params.species_param[ispec].species_type).c_str());
        rho_s[ispec] = new Field2D(dimPrim, ("Rho_"+params.species_param[ispec].species_type).c_str());
    }

    
//    ostringstream file_name("");
//    for (unsigned int ispec=0; ispec<n_species; ispec++) {
//        file_name.str("");
//        file_name << "Jx_s" << ispec;
//        Jx_s[ispec]  = new Field2D(dimPrim, 0, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "Jy_s" << ispec;
//        Jy_s[ispec]  = new Field2D(dimPrim, 1, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "Jz_s" << ispec;
//        Jz_s[ispec]  = new Field2D(dimPrim, 2, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "rho_s" << ispec;
//        rho_s[ispec] = new Field2D(dimPrim, file_name.str().c_str());
//    }
    
    
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }
    /*
     MESSAGE("index_bc_min / index_bc_max / nx_p / nx_d" << index_bc_min[0]
            << " " << index_bc_max[0] << " " << nx_p<< " " << nx_d);
     */
    
    
    // Define limits of non duplicated elements
    // (by construction 1 (prim) or 2 (dual) elements shared between per MPI process)
    // istart 
    for (unsigned int i=0 ; i<3 ; i++)
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
            istart[i][isDual] = 0;
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        for (unsigned int isDual=0 ; isDual<2 ; isDual++) {
            istart[i][isDual] = oversize[i];
            if (smpi2D->getProcCoord(i)!=0) istart[i][isDual]+=1;
        }
    }
    
    // bufsize = nelements
    for (unsigned int i=0 ; i<3 ; i++) 
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
            bufsize[i][isDual] = 1;
    
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        for (int isDual=0 ; isDual<2 ; isDual++)
            bufsize[i][isDual] = n_space[i] + 1;
        
        for (int isDual=0 ; isDual<2 ; isDual++) {
            bufsize[i][isDual] += isDual; 
            if ( smpi2D->getNbrOfProcs(i)!=1 ) {
                
                if ( ( !isDual ) && (smpi2D->getProcCoord(i)!=0) )
                    bufsize[i][isDual]--;
                else if  (isDual) {
                    bufsize[i][isDual]--;
                    if ( (smpi2D->getProcCoord(i)!=0) && (smpi2D->getProcCoord(i)!=smpi2D->getNbrOfProcs(i)-1) ) 
                        bufsize[i][isDual]--;
                }
                
            } // if ( smpi2D->getNbrOfProcs(i)!=1 )
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field 
    
    
}//END constructor Electromagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::~ElectroMagn2D()
{
    
}//END ElectroMagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Solve Poisson
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solvePoisson(SmileiMPI* smpi)
{
    
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    
    unsigned int iteration_max = 50000;
    double       error_max     = 1.e-14;
    
    Field2D* Ex2D  = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D  = static_cast<Field2D*>(Ey_);
    Field2D* rho2D = static_cast<Field2D*>(rho_);
    
    double one_ov_dx_sq       = 1.0/(dx*dx);
    double one_ov_dy_sq       = 1.0/(dy*dy);
    double two_ov_dx2dy2      = 2.0*(1.0/(dx*dx)+1.0/(dy*dy));
    
    unsigned int nx_p2_global = (smpi2D->n_space_global[0]+1) * (smpi2D->n_space_global[1]+1);
    unsigned int smilei_sz    = smpi2D->smilei_sz;
    unsigned int smilei_rk    = smpi2D->smilei_rk;
    
    
    // Boundary condition for the direction vector (all put to 0)
    vector<double> pSouth(nx_p);
    for (unsigned int i=0; i<nx_p; i++) pSouth[i]=0.0;
    vector<double> pNorth(nx_p);
    for (unsigned int i=0; i<nx_p; i++) pNorth[i]=0.0;
    vector<double> pWest(ny_p);
    for (unsigned int j=0; j<ny_p; j++) pWest[j]=0.0;
    vector<double> pEast(ny_p);
    for (unsigned int j=0; j<ny_p; j++) pEast[j]=0.0;
    
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // ----------------------------------------------------------------------------------
    vector<unsigned int> index_min_p(2);
    vector<unsigned int> index_min_d(2);
    vector<unsigned int> index_max_p(2);
    vector<unsigned int> index_max_d(2);
    index_min_p[0] = oversize[0];
    index_min_p[1] = oversize[1];
    index_min_d[0] = oversize[0];
    index_min_d[1] = oversize[1];
    index_max_p[0] = nx_p - 2 - oversize[0];
    index_max_p[1] = ny_p - 2 - oversize[1];
    index_max_d[0] = nx_d - 2 - oversize[0];
    index_max_d[1] = ny_d - 2 - oversize[1];
    if (smpi2D->isWestern()) {
        index_min_p[0] = 0;
        index_min_d[0] = 0;
    }
    if (smpi2D->isEastern()) {
        index_max_p[0] = nx_p-1;
        index_max_d[0] = nx_d-1;
    }
    
    
    // Initialization of the variables
    // -------------------------------
    DEBUG(1,"Initialisation for the iterative CG method started");
    unsigned int iteration=0;
    
    Field2D phi(dimPrim);    // scalar potential
    Field2D r(dimPrim);      // residual vector
    Field2D p(dimPrim);      // direction vector
    Field2D Ap(dimPrim);     // A*p vector
    
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            phi(i,j)   = 0.0;
            r(i,j)     = -(*rho2D)(i,j);
            p(i,j)     = r(i,j);
        }//j
    }//i
    
    // norm of the residual
    double rnew_dot_rnew       = 0.0;
    double rnew_dot_rnew_local = 0.0;
    for (unsigned int i=index_min_p[0]; i<=index_max_p[0]; i++) {
        for (unsigned int j=index_min_p[1]; j<=index_max_p[1]; j++) {
            rnew_dot_rnew_local += r(i,j)*r(i,j);
        }
    }
    MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // compute control parameter
    double ctrl = rnew_dot_rnew / (double)(nx_p2_global);
    
    DEBUG(1,"Initialisation for the iterative CG method done");
    
    
    // ---------------------------------------------------------
    // Starting iterative loop for the conjugate gradient method
    // ---------------------------------------------------------
    DEBUG(1,"Starting iterative loop for CG method");
    while ( (ctrl > error_max) && (iteration<iteration_max) ) {
        
        iteration++;
        DEBUG(5,"iteration " << iteration << " started with control parameter ctrl = " << ctrl*1.e14 << " x 1e-14");
        
        // scalar product of the residual
        double r_dot_r = rnew_dot_rnew;
        
        
        // vector product Ap = A*p
        for (unsigned int i=1; i<nx_p-1; i++) {
            for (unsigned int j=1; j<ny_p-1; j++) {
                Ap(i,j) = one_ov_dx_sq*(p(i-1,j)+p(i+1,j)) + one_ov_dy_sq*(p(i,j-1)+p(i,j+1)) - two_ov_dx2dy2*p(i,j);
            }//j
        }//i
        
        
        // Western BC
        if ( smpi2D->isWestern() ) {
            for (unsigned int j=1; j<ny_p-1; j++) {
                Ap(0,j)      = one_ov_dx_sq*(pWest[j]+p(1,j))
                +              one_ov_dy_sq*(p(0,j-1)+p(0,j+1))
                -              two_ov_dx2dy2*p(0,j);
            }
            // at corners
            Ap(0,0)           = one_ov_dx_sq*(pWest[0]+p(1,0))               // West/South
            +                   one_ov_dy_sq*(pSouth[0]+p(0,1))
            -                   two_ov_dx2dy2*p(0,0);
            Ap(0,ny_p-1)      = one_ov_dx_sq*(pWest[ny_p-1]+p(1,ny_p-1))     // West/North
            +                   one_ov_dy_sq*(p(0,ny_p-2)+pNorth[0])
            -                   two_ov_dx2dy2*p(0,ny_p-1);
        }
        
        // Eastern BC
        if ( smpi2D->isEastern() ) {
            
            for (unsigned int j=1; j<ny_p-1; j++) {
                Ap(nx_p-1,j) = one_ov_dx_sq*(p(nx_p-2,j)+pEast[j])
                +              one_ov_dy_sq*(p(nx_p-1,j-1)+p(nx_p-1,j+1))
                -              two_ov_dx2dy2*p(nx_p-1,j);
            }
            // at corners
            Ap(nx_p-1,0)      = one_ov_dx_sq*(p(nx_p-2,0)+pEast[0])                 // East/South
            +                   one_ov_dy_sq*(pSouth[nx_p-1]+p(nx_p-1,1))
            -                   two_ov_dx2dy2*p(nx_p-1,0);
            Ap(nx_p-1,ny_p-1) = one_ov_dx_sq*(p(nx_p-2,ny_p-1)+pEast[ny_p-1])       // East/North
            +                   one_ov_dy_sq*(p(nx_p-1,ny_p-2)+pNorth[nx_p-1])
            -                   two_ov_dx2dy2*p(nx_p-1,ny_p-1);
        }
        
        // Periodic BC on the y-direction
        smpi2D->exchangeField(&Ap);
        
        // scalar product p.Ap
        double p_dot_Ap       = 0.0;
        double p_dot_Ap_local = 0.0;
        for (unsigned int i=index_min_p[0]; i<=index_max_p[0]; i++) {
            for (unsigned int j=index_min_p[1]; j<=index_max_p[1]; j++) {
                p_dot_Ap_local += p(i,j)*Ap(i,j);
            }
        }
        MPI_Allreduce(&p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // compute new potential and residual
        double alpha_k = r_dot_r/p_dot_Ap;
        for(unsigned int i=0; i<nx_p; i++) {
            for(unsigned int j=0; j<ny_p; j++) {
                phi(i,j) += alpha_k * p(i,j);
                r(i,j)   -= alpha_k * Ap(i,j);
            }
        }
        
        // compute new residual norm
        rnew_dot_rnew       = 0.0;
        rnew_dot_rnew_local = 0.0;
        for (unsigned int i=index_min_p[0]; i<=index_max_p[0]; i++) {
            for (unsigned int j=index_min_p[1]; j<=index_max_p[1]; j++) {
                rnew_dot_rnew_local += r(i,j)*r(i,j);
            }
        }
        MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        DEBUG(10,"new residual norm: rnew_dot_rnew = " << rnew_dot_rnew);
        
        // compute new direction
        double beta_k = rnew_dot_rnew/r_dot_r;
        for (unsigned int i=0; i<nx_p; i++) {
            for(unsigned int j=0; j<ny_p; j++) {
                p(i,j) = r(i,j) + beta_k * p(i,j);
            }
        }
        
        // compute control parameter
        ctrl = rnew_dot_rnew / (double)(nx_p2_global);
        DEBUG(10,"iteration " << iteration << " done, exiting with control parameter ctrl = " << ctrl);
        
    }//End of the iterative loop
    
    
    // --------------------------------
    // Status of the solver convergence
    // --------------------------------
    if (iteration == iteration_max) {
        if (smpi->isMaster())
            WARNING("Poisson solver did not converge: reached maximum iteration number: " << iteration
                    << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    else {
        if (smpi->isMaster()) 
            MESSAGE(1,"Poisson solver converged at iteration: " << iteration
                    << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    
    
    // ------------------------------------------
    // Compute the electrostatic fields Ex and Ey
    // ------------------------------------------
    
    // Ex
    DEBUG(2, "Computing Ex from scalar potential");
    for (unsigned int i=1; i<nx_d-1; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(i,j) = (phi(i-1,j)-phi(i,j))/dx;
        }
    }
    // Ey
    DEBUG(2, "Computing Ey from scalar potential");
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=1; j<ny_d-1; j++) {
            (*Ey2D)(i,j) = (phi(i,j-1)-phi(i,j))/dy;
        }
    }
    
    
    // Apply BC on Ex and Ey
    // ---------------------
    // Ex / West
    if (smpi2D->isWestern()) {
        DEBUG(2, "Computing Western BC on Ex");
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(0,j) = (*Ex2D)(1,j) + ((*Ey2D)(0,j+1)-(*Ey2D)(0,j))*dx/dy  - dx*(*rho2D)(0,j);
        }
    }
    // Ex / East
    if (smpi2D->isEastern()) {
        DEBUG(2, "Computing Eastern BC on Ex");
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(nx_d-1,j) = (*Ex2D)(nx_d-2,j) - ((*Ey2D)(nx_p-1,j+1)-(*Ey2D)(nx_p-1,j))*dx/dy + dx*(*rho2D)(nx_p-1,j);
        }
    }
    
    smpi2D->exchangeField(Ex2D);
    smpi2D->exchangeField(Ey2D);
    
    
    // Centering of the electrostatic fields
    // -------------------------------------
    
    double Ex_WestNorth = 0.0;
    double Ey_WestNorth = 0.0;
    double Ex_EastSouth = 0.0;
    double Ey_EastSouth = 0.0;
    
    // West-North corner
    int rank_WestNorth = smpi2D->extrem_ranks[0][1];
    if ( (smpi2D->isWestern()) && (smpi2D->isNorthern())) {
        if (smpi2D->smilei_rk != smpi2D->extrem_ranks[0][1]) ERROR("west-north process rank not well defined");
        Ex_WestNorth = (*Ex2D)(0,ny_p-1);
        Ey_WestNorth = (*Ey2D)(0,ny_d-1);
    }
    MPI_Bcast(&Ex_WestNorth, 1, MPI_DOUBLE, rank_WestNorth, MPI_COMM_WORLD);
    MPI_Bcast(&Ey_WestNorth, 1, MPI_DOUBLE, rank_WestNorth, MPI_COMM_WORLD);
    
    // East-South corner
    int rank_EastSouth = smpi2D->extrem_ranks[1][0];
    if ((smpi2D->isEastern()) && (smpi2D->isSouthern())) {
        if (smpi2D->smilei_rk != smpi2D->extrem_ranks[1][0]) ERROR("east-south process rank not well defined");
        Ex_EastSouth = (*Ex2D)(nx_d-1,0);
        Ey_EastSouth = (*Ey2D)(nx_p-1,0);
    }
    MPI_Bcast(&Ex_EastSouth, 1, MPI_DOUBLE, rank_EastSouth, MPI_COMM_WORLD);
    MPI_Bcast(&Ey_EastSouth, 1, MPI_DOUBLE, rank_EastSouth, MPI_COMM_WORLD);
    
    
    // Centering electrostatic fields
    double Ex_Add = -0.5*(Ex_WestNorth+Ex_EastSouth);
    double Ey_Add = -0.5*(Ey_WestNorth+Ey_EastSouth);
    
    for (unsigned int i=0; i<nx_d; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(i,j) += Ex_Add;
        }
    }
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_d; j++) {
            (*Ey2D)(i,j) += Ey_Add;
        }
    }
    
}//END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::saveMagneticFields()
{
    // Static cast of the fields
    Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
    Field2D* By2D   = static_cast<Field2D*>(By_);
    Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
    Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);
    
    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D_m)(i,j)=(*Bx2D)(i,j);
        }
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j)=(*By2D)(i,j);
        }
    }
    
    // Magnetic field Bz^(d,d)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_m)(i,j)=(*Bz2D)(i,j);
        }
    }
    
}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Solve the Maxwell-Ampere equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwellAmpere()
{
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(Bx_);
    Field2D* By2D = static_cast<Field2D*>(By_);
    Field2D* Bz2D = static_cast<Field2D*>(Bz_);
    Field2D* Jx2D = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D = static_cast<Field2D*>(Jz_);
    
    // Electric field Ex^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ex2D)(i,j) += -timestep*(*Jx2D)(i,j) + dt_ov_dy * ( (*Bz2D)(i,j+1) - (*Bz2D)(i,j) );
        }
    }
    
    // Electric field Ey^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Ey2D)(i,j) += -timestep*(*Jy2D)(i,j) - dt_ov_dx * ( (*Bz2D)(i+1,j) - (*Bz2D)(i,j) );
        }
    }
    
    // Electric field Ez^(p,p)
    for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez2D)(i,j) += -timestep*(*Jz2D)(i,j)
            +               dt_ov_dx * ( (*By2D)(i+1,j) - (*By2D)(i,j) )
            -               dt_ov_dy * ( (*Bx2D)(i,j+1) - (*Bx2D)(i,j) );
        }
    }
    
}//END solveMaxwellAmpere



// ---------------------------------------------------------------------------------------------------------------------
// Solve the Maxwell-Faraday equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwellFaraday()
{
    
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(Bx_);
    Field2D* By2D = static_cast<Field2D*>(By_);
    Field2D* Bz2D = static_cast<Field2D*>(Bz_);
    
    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p;  i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*Bx2D)(i,j) -= dt_ov_dy * ( (*Ez2D)(i,j) - (*Ez2D)(i,j-1) );
        }
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D)(i,j) += dt_ov_dx * ( (*Ez2D)(i,j) - (*Ez2D)(i-1,j) );
        }
    }
    
    // Magnetic field Bz^(d,d)
    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*Bz2D)(i,j) += dt_ov_dy * ( (*Ex2D)(i,j) - (*Ex2D)(i,j-1) )
            -               dt_ov_dx * ( (*Ey2D)(i,j) - (*Ey2D)(i-1,j) );
        }
    }
    
}//END solveMaxwellFaraday


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::centerMagneticFields()
{
    // Static cast of the fields
    Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
    Field2D* By2D   = static_cast<Field2D*>(By_);
    Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
    Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);
    
    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D_m)(i,j) = ( (*Bx2D)(i,j) + (*Bx2D_m)(i,j) )*0.5;
        }
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j) = ( (*By2D)(i,j) + (*By2D_m)(i,j) )*0.5;
        }
    }
    
    // Magnetic field Bz^(d,d)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_m)(i,j) = ( (*Bz2D)(i,j) + (*Bz2D_m)(i,j) )*0.5;
        }
    }
    
}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::incrementAvgFields(unsigned int time_step, unsigned int ntime_step_avg)
{
    // Static cast of the fields
    Field2D* Ex2D     = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D     = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D     = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D_m   = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m   = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m   = static_cast<Field2D*>(Bz_m);
    Field2D* Ex2D_avg = static_cast<Field2D*>(Ex_avg);
    Field2D* Ey2D_avg = static_cast<Field2D*>(Ey_avg);
    Field2D* Ez2D_avg = static_cast<Field2D*>(Ez_avg);
    Field2D* Bx2D_avg = static_cast<Field2D*>(Bx_avg);
    Field2D* By2D_avg = static_cast<Field2D*>(By_avg);
    Field2D* Bz2D_avg = static_cast<Field2D*>(Bz_avg);
    
    // reset the averaged fields for (time_step-1)%ntime_step_avg == 0
    if ( (time_step-1)%ntime_step_avg==0 ){
        Ex2D_avg->put_to(0.0);
        Ey2D_avg->put_to(0.0);
        Ez2D_avg->put_to(0.0);
        Bx2D_avg->put_to(0.0);
        By2D_avg->put_to(0.0);
        Bz2D_avg->put_to(0.0);
    }
    
    // increment the time-averaged fields
    
    // Electric field Ex^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ex2D_avg)(i,j) += (*Ex2D)(i,j);
        }
    }
    
    // Electric field Ey^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Ey2D_avg)(i,j) += (*Ey2D)(i,j);
        }
    }
    
    // Electric field Ez^(p,p)
    for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez2D_avg)(i,j) += (*Ez2D)(i,j);
        }
    }
    
    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D_avg)(i,j) += (*Bx2D_m)(i,j);
        }
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_avg)(i,j) += (*By2D_m)(i,j);
        }
    }
    
    // Magnetic field Bz^(d,d)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_avg)(i,j) += (*Bz2D_m)(i,j);
        }
    }
    
    
}//END incrementAvgFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge densities and currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::restartRhoJ()
{
    // --------------------------
    // Total currents and density
    // --------------------------
    
    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
    Field2D* rho2D   = static_cast<Field2D*>(rho_);
    
    // Charge density rho^(p,p) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*rho2D)(i,j) = 0.0;
        }
    }
    
    // Current Jx^(d,p) to 0
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Jx2D)(i,j) = 0.0;
        }
    }
    
    // Current Jy^(p,d) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Jy2D)(i,j) = 0.0;
        }
    }
    
    // Current Jz^(p,p) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Jz2D)(i,j) = 0.0;
        }
    }
}//END restartRhoJ
    
    
void ElectroMagn2D::restartRhoJs(int ispec, bool currents)
{
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
    Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
    Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
    Field2D* rho2D_s = static_cast<Field2D*>(rho_s[ispec]);
    
    // Charge density rho^(p,p) to 0
    #pragma omp for schedule(static)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*rho2D_s)(i,j) = 0.0;
        }
    }
    if (currents){ 
        // Current Jx^(d,p) to 0
        #pragma omp for schedule(static)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jx2D_s)(i,j) = 0.0;
            }
        }
        
        // Current Jy^(p,d) to 0
        #pragma omp for schedule(static)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Jy2D_s)(i,j) = 0.0;
            }
        }
        
        // Current Jz^(p,p) to 0
        #pragma omp for schedule(static)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jz2D_s)(i,j) = 0.0;
            }
        }
    }
}//END restartRhoJs
    



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::computeTotalRhoJ()
{
    
    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
    Field2D* rho2D   = static_cast<Field2D*>(rho_);
    
    
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
        Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
        Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
        Field2D* rho2D_s = static_cast<Field2D*>(rho_s[ispec]);
        
        // Charge density rho^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*rho2D)(i,j) += (*rho2D_s)(i,j);
            }
        }
        
        // Current Jx^(d,p) to 0
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jx2D)(i,j) += (*Jx2D_s)(i,j);
            }
        }
        
        // Current Jy^(p,d) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Jy2D)(i,j) += (*Jy2D_s)(i,j);
            }
        }
        
        // Current Jz^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jz2D)(i,j) += (*Jz2D_s)(i,j);
            }
        }
        
    }//END loop on species ispec
    
}//END computeTotalRhoJ


void ElectroMagn2D::computePoynting() {

    if (isWestern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)];
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBy=istart[0][By_m->isDual(0)];
        
        unsigned int jEy=istart[1][Ey_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];
        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBy=istart[1][By_m->isDual(1)];
        
        for (unsigned int j=0; j<=bufsize[1][Ez_->isDual(1)]; j++) {
            
            double Ey__ = 0.5*((*Ey_)(iEy,jEy+j) + (*Ey_)(iEy, jEy+j+1));
            double Bz__ = 0.25*((*Bz_m)(iBz,jBz+j)+(*Bz_m)(iBz+1,jBz+j)+(*Bz_m)(iBz,jBz+j+1)+(*Bz_m)(iBz+1,jBz+j+1));
            double Ez__ = (*Ez_)(iEz,jEz+j);
            double By__ = 0.5*((*By_m)(iBy,jBy+j) + (*By_m)(iBy+1, jBy+j));
            
            poynting_inst[0][0] = dy*timestep*(Ey__*Bz__ - Ez__*By__);
            poynting[0][0]+= poynting_inst[0][0];
        }
    }//if Western
    
    
    if (isEastern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)]  + bufsize[0][Ey_->isDual(0)] -1;
        unsigned int iBz=istart[0][Bz_m->isDual(0)] + bufsize[0][Bz_m->isDual(0)]-1;
        unsigned int iEz=istart[0][Ez_->isDual(0)]  + bufsize[0][Ez_->isDual(0)] -1;
        unsigned int iBy=istart[0][By_m->isDual(0)] + bufsize[0][By_m->isDual(0)]-1;
        
        unsigned int jEy=istart[1][Ey_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];
        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBy=istart[1][By_m->isDual(1)];
        
        for (unsigned int j=0; j<=bufsize[1][Ez_->isDual(1)]; j++) {
            
            double Ey__ = 0.5*((*Ey_)(iEy,jEy+j) + (*Ey_)(iEy, jEy+j+1));
            double Bz__ = 0.25*((*Bz_m)(iBz,jBz+j)+(*Bz_m)(iBz+1,jBz+j)+(*Bz_m)(iBz,jBz+j+1)+(*Bz_m)(iBz+1,jBz+j+1));
            double Ez__ = (*Ez_)(iEz,jEz+j);
            double By__ = 0.5*((*By_m)(iBy,jBy+j) + (*By_m)(iBy+1, jBy+j));
            
            poynting_inst[1][0] = dy*timestep*(Ey__*Bz__ - Ez__*By__);
            poynting[1][0] -= poynting_inst[1][0];
        }
    }//if Easter
    
    if (isSouthern) {
        
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)]; 
        
        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBx=istart[1][Bx_m->isDual(1)];
        unsigned int jEx=istart[1][Ex_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];
        
        for (unsigned int i=0; i<=bufsize[0][Ez_->isDual(0)]; i++) {
            double Ez__ = (*Ez_)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx_m)(iBx+i,jBx) + (*Bx_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex_)(iEx+i,jEx) + (*Ex_)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz_m)(iBz+i,jBz)+(*Bz_m)(iBz+i+1,jBz)+(*Bz_m)(iBz+i,jBz+1)+(*Bz_m)(iBz+i+1,jBz+1));
            
            poynting_inst[0][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[0][1] += poynting_inst[0][1];
        }
    }// if South
    
    if (isNorthern) {
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)];
        
        unsigned int jEz=istart[1][Ez_->isDual(1)]  + bufsize[1][Ez_->isDual(1)] -1;
        unsigned int jBx=istart[1][Bx_m->isDual(1)] + bufsize[1][Bx_m->isDual(1)]-1;
        unsigned int jEx=istart[1][Ex_->isDual(1)]  + bufsize[1][Ex_->isDual(1)] -1;
        unsigned int jBz=istart[1][Bz_m->isDual(1)] + bufsize[1][Bz_m->isDual(1)]-1;
        
        for (unsigned int i=0; i<=bufsize[0][Ez_->isDual(0)]; i++) {
            double Ez__ = (*Ez_)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx_m)(iBx+i,jBx) + (*Bx_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex_)(iEx+i,jEx) + (*Ex_)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz_m)(iBz+i,jBz)+(*Bz_m)(iBz+i+1,jBz)+(*Bz_m)(iBz+i,jBz+1)+(*Bz_m)(iBz+i+1,jBz+1));
            
            poynting_inst[1][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[1][1] -= poynting_inst[1][1];
        }
    }//if North

}

void ElectroMagn2D::applyExternalField(Field* my_field,  ExtFieldProfile *my_profile, SmileiMPI* smpi) {
    
    Field2D* field2D=static_cast<Field2D*>(my_field);
    ExtFieldProfile2D* profile=static_cast<ExtFieldProfile2D*> (my_profile);
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

    vector<double> pos(2,0);
    
    for (int i=0 ; i<field2D->dims()[0] ; i++) {
        pos[0] = ( (double)(smpi2D->getCellStartingGlobalIndex(0)+i +(field2D->isDual(0)?-0.5:0)) )*dx;
        
        for (int j=0 ; j<field2D->dims()[1] ; j++) {
            
            pos[1] = ( (double)(smpi2D->getCellStartingGlobalIndex(1)+i +(field2D->isDual(1)?-0.5:0)) )*dy;
            
            (*field2D)(i,j) = (*field2D)(i,j) + (*profile)(pos);
        }
    }
}
