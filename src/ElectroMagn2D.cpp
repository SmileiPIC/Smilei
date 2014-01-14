#include "ElectroMagn2D.h"
#include "PicParams.h"
#include "Field2D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart2D.h"

#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::ElectroMagn2D(PicParams* params, SmileiMPI* smpi)
: ElectroMagn(params, smpi)
{
    // local dt to store
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    int process_coord_x = smpi2D->getProcCoord(0);
    
    
    // --------------------------------------------------
    // Calculate quantities related to the simulation box
    // --------------------------------------------------
    
    // time-step
    dt       = params->timestep;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = params->cell_length[0];
    dt_ov_dx = dt/dx;
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = params->cell_length[1];
    dt_ov_dy = dt/dy;
    dy_ov_dt = 1.0/dt_ov_dy;
    
    //!\todo Check if oversize really needs to be a vector
    oversize_ = params->oversize[0];
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    // West boundary
    double theta  = 0.0; //! \todo Introduce in parameters for Boundary cond., e.g., params->EMBoundary->theta_W
    double factor = 1.0 / (cos(theta) + dt_ov_dx);
    Alpha_SM_W    = 2.0                     * factor;
    Beta_SM_W     = - (cos(theta)-dt_ov_dx) * factor;
    Gamma_SM_W    = 4.0 * cos(theta)        * factor;
    Delta_SM_W    = - (sin(theta)+dt_ov_dy) * factor;
    Epsilon_SM_W  = - (sin(theta)-dt_ov_dy) * factor;
    MESSAGE("WEST : " << Alpha_SM_W << Beta_SM_W << Gamma_SM_W);
    
    // East boundary
    theta         = M_PI;
    factor        = 1.0 / (cos(theta) - dt_ov_dx);
    Alpha_SM_E    = 2.0                      * factor;
    Beta_SM_E     = - (cos(theta)+dt_ov_dx)  * factor;
    Gamma_SM_E    = 4.0 * cos(theta)         * factor;
    Delta_SM_E    = - (sin(theta)+dt_ov_dy)  * factor;
    Epsilon_SM_E  = - (sin(theta)-dt_ov_dy)  * factor;
    MESSAGE("EAST : " << Alpha_SM_E << Beta_SM_E << Gamma_SM_E);
    
    // ----------------------
    // Electromagnetic fields
    // ----------------------
    //! \todo{Homogenize 1D/2D dimPrim/dimDual or nx_p/nx_d/ny_p/ny_d}
    
    dimPrim.resize( params->nDim_field );
    dimDual.resize( params->nDim_field );
    
    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<params->nDim_field ; i++) {
		// Standard scheme
		dimPrim[i] = params->n_space[i]+1;
		dimDual[i] = params->n_space[i]+2;
		// + Ghost domain
		dimPrim[i] += 2*params->oversize[i];
		dimDual[i] += 2*params->oversize[i];
    }
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params->n_space[0]+1+2*params->oversize[0];
    nx_d = params->n_space[0]+2+2*params->oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params->n_space[1]+1+2*params->oversize[1];
    ny_d = params->n_space[1]+2+2*params->oversize[1];
	
	// Allocation of the EM fields
    Ex_ = new Field2D( dimPrim, 0, false, "Ex" );
	Ey_ = new Field2D( dimPrim, 1, false, "Ey");
	Ez_ = new Field2D( dimPrim, 2, false, "Ez");
	Bx_ = new Field2D( dimPrim, 0, true, "Bx");
	By_ = new Field2D( dimPrim, 1, true, "By");
	Bz_ = new Field2D( dimPrim, 2, true, "Bz");
	Bx_m = new Field2D(dimPrim, 0, true, "Bx_m");
	By_m = new Field2D(dimPrim, 1, true, "By_m");
	Bz_m = new Field2D(dimPrim, 2, true, "Bz_m");
	
	// Total charge currents and densities
	Jx_ = new Field2D(dimPrim, 0, false, "Jx");
	Jy_ = new Field2D(dimPrim, 1, false, "Jy");
	Jz_ = new Field2D(dimPrim, 2, false, "Jz");
	rho_ = new Field2D(dimPrim, "Rho" );
	rho_o = new Field2D(dimPrim, "Rho_old" );
    
    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++){
        Jx_s[ispec]  = new Field2D(dimPrim, 0, false);
        Jy_s[ispec]  = new Field2D(dimPrim, 1, false);
        Jz_s[ispec]  = new Field2D(dimPrim, 2, false);
        rho_s[ispec] = new Field2D(dimPrim);
    }
	
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( params->nDim_field, 0 );
    index_bc_max.resize( params->nDim_field, 0 );
    for (unsigned int i=0 ; i<params->nDim_field ; i++) {
        index_bc_min[i] = params->oversize[i];
        index_bc_max[i] = dimDual[i]-params->oversize[i]-1;
    }
    MESSAGE("index_bc_min / index_bc_max / nx_p / nx_d" << index_bc_min[0]
            << " " << index_bc_max[0] << " " << nx_p<< " " << nx_d);
	
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
    DEBUG("Entering Poisson Solver");
    
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
    vector<double> pSouth(nx_p); for (unsigned int i=0; i<nx_p; i++) pSouth[i]=0.0;
    vector<double> pNorth(nx_p); for (unsigned int i=0; i<nx_p; i++) pNorth[i]=0.0;
    vector<double> pWest(ny_p);  for (unsigned int j=0; j<ny_p; j++) pWest[j]=0.0;
    vector<double> pEast(ny_p);  for (unsigned int j=0; j<ny_p; j++) pEast[j]=0.0;
    
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // ----------------------------------------------------------------------------------
    vector<unsigned int> index_min_p(2);
    vector<unsigned int> index_min_d(2);
    vector<unsigned int> index_max_p(2);
    vector<unsigned int> index_max_d(2);
    index_min_p[0] = oversize_;
    index_min_p[1] = oversize_;
    index_min_d[0] = oversize_;
    index_min_d[1] = oversize_;
    index_max_p[0] = nx_p - 2 - oversize_;
    index_max_p[1] = ny_p - 2 - oversize_;
    index_max_d[0] = nx_d - 2 - oversize_;
    index_max_d[1] = ny_d - 2 - oversize_;
    if (smpi2D->isWester()) {
        index_min_p[0] = 0;
        index_min_d[0] = 0;
    }
    if (smpi2D->isEaster()) {
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
    
    for (unsigned int i=0; i<nx_p; i++){
        for (unsigned int j=0; j<ny_p; j++){
            phi(i,j)   = 0.0;
            r(i,j)     = -(*rho2D)(i,j);
            p(i,j)     = r(i,j);
        }//j
    }//i
    
    // norm of the residual
    double rnew_dot_rnew       = 0.0;
    double rnew_dot_rnew_local = 0.0;
    for (unsigned int i=index_min_p[0]; i<=index_max_p[0]; i++){
        for (unsigned int j=index_min_p[1]; j<=index_max_p[1]; j++){
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
        for (unsigned int i=1; i<nx_p-1; i++){
            for (unsigned int j=1; j<ny_p-1; j++){
                Ap(i,j) = one_ov_dx_sq*(p(i-1,j)+p(i+1,j)) + one_ov_dy_sq*(p(i,j-1)+p(i,j+1)) - two_ov_dx2dy2*p(i,j);
            }//j
        }//i
        
        
        // Western BC
        if ( smpi2D->isWester() ) {
            for (unsigned int j=1; j<ny_p-1; j++){
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
        if ( smpi2D->isEaster() ) {
            
            for (unsigned int j=1; j<ny_p-1; j++){
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
        for (unsigned int i=index_min_p[0]; i<=index_max_p[0]; i++){
            for (unsigned int j=index_min_p[1]; j<=index_max_p[1]; j++){
                p_dot_Ap_local += p(i,j)*Ap(i,j);
            }
        }
        MPI_Allreduce(&p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // compute new potential and residual
        double alpha_k = r_dot_r/p_dot_Ap;
        for(unsigned int i=0; i<nx_p; i++){
            for(unsigned int j=0; j<ny_p; j++){
                phi(i,j) += alpha_k * p(i,j);
                r(i,j)   -= alpha_k * Ap(i,j);
            }
        }
        
        // compute new residual norm
        rnew_dot_rnew       = 0.0;
        rnew_dot_rnew_local = 0.0;
        for (unsigned int i=index_min_p[0]; i<=index_max_p[0]; i++){
            for (unsigned int j=index_min_p[1]; j<=index_max_p[1]; j++){
                rnew_dot_rnew_local += r(i,j)*r(i,j);
            }
        }
        MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        DEBUG(10,"new residual norm: rnew_dot_rnew = " << rnew_dot_rnew);
        
        // compute new direction
        double beta_k = rnew_dot_rnew/r_dot_r;
        for (unsigned int i=0; i<nx_p; i++){
            for(unsigned int j=0; j<ny_p; j++){
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
    if (iteration == iteration_max){
        WARNING("Poisson solver did not converge: reached maximum iteration number: " << iteration
                << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    else{
        MESSAGE("Poisson solver converged at iteration: " << iteration
                << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    
    
    // ------------------------------------------
    // Compute the electrostatic fields Ex and Ey
    // ------------------------------------------
    
    // Ex
    DEBUG(2, "Computing Ex from scalar potential");
    for (unsigned int i=1; i<nx_d-1; i++){
        for (unsigned int j=0; j<ny_p; j++){
            (*Ex2D)(i,j) = (phi(i-1,j)-phi(i,j))/dx;
        }
    }
    // Ey
    DEBUG(2, "Computing Ey from scalar potential");
    for (unsigned int i=0; i<nx_p; i++){
        for (unsigned int j=1; j<ny_d-1; j++){
            (*Ey2D)(i,j) = (phi(i,j-1)-phi(i,j))/dy;
        }
    }
    
    
    // Apply BC on Ex and Ey
    // ---------------------
    // Ex / West
    if (smpi2D->isWester()){
        DEBUG(2, "Computing Western BC on Ex");
        for (unsigned int j=0; j<ny_p; j++){
            (*Ex2D)(0,j) = (*Ex2D)(1,j) + ((*Ey2D)(0,j+1)-(*Ey2D)(0,j))*dx/dy  - dx*(*rho2D)(0,j);
        }
    }
    // Ex / East
    if (smpi2D->isEaster()){
        DEBUG(2, "Computing Eastern BC on Ex");
        for (unsigned int j=0; j<ny_p; j++){
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
    if ( (smpi2D->isWester()) && (smpi2D->isNorthern())) {
        if (smpi2D->smilei_rk != smpi2D->extrem_ranks[0][1]) ERROR("west-north process rank not well defined");
        Ex_WestNorth = (*Ex2D)(0,ny_p-1);
        Ey_WestNorth = (*Ey2D)(0,ny_d-1);
    }
    MPI_Bcast(&Ex_WestNorth, 1, MPI_DOUBLE, rank_WestNorth, MPI_COMM_WORLD);
    MPI_Bcast(&Ey_WestNorth, 1, MPI_DOUBLE, rank_WestNorth, MPI_COMM_WORLD);
    
    // East-South corner
    int rank_EastSouth = smpi2D->extrem_ranks[1][0];
    if ((smpi2D->isEaster()) && (smpi2D->isSouthern())) {
        if (smpi2D->smilei_rk != smpi2D->extrem_ranks[1][0]) ERROR("east-south process rank not well defined");
        Ex_EastSouth = (*Ex2D)(nx_d-1,0);
        Ey_EastSouth = (*Ey2D)(nx_p-1,0);
    }
    MPI_Bcast(&Ex_EastSouth, 1, MPI_DOUBLE, rank_EastSouth, MPI_COMM_WORLD);
    MPI_Bcast(&Ey_EastSouth, 1, MPI_DOUBLE, rank_EastSouth, MPI_COMM_WORLD);
    
    cerr << "Ex_WestNorth = " << Ex_WestNorth << "  -  Ey_WestNorth = " << Ey_WestNorth << endl;
    cerr << "Ex_EastSouth = " << Ex_EastSouth << "  -  Ey_EastSouth = " << Ey_EastSouth << endl;
    
    // Centering electrostatic fields
    double Ex_Add = -0.5*(Ex_WestNorth+Ex_EastSouth);
    double Ey_Add = -0.5*(Ey_WestNorth+Ey_EastSouth);
    
    for (unsigned int i=0; i<nx_d; i++){
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(i,j) += Ex_Add;
        }
    }
    for (unsigned int i=0; i<nx_p; i++){
        for (unsigned int j=0; j<ny_d; j++) {
            (*Ey2D)(i,j) += Ey_Add;
        }
    }
    
}//END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwell(double time_dual, SmileiMPI* smpi)
{
    //! \todo All this is generic (does not depend on geometry) move to ElectroMagn.cpp ??? (MG to JD)
    saveMagneticFields();
    solveMaxwellAmpere();
    solveMaxwellFaraday();
    smpi->exchangeB( this );
    applyEMBoundaryConditions(time_dual, smpi);
    centerMagneticFields();
    
}//END ElectroMagn2D



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
			(*Ex2D)(i,j) += -dt*(*Jx2D)(i,j) + dt_ov_dy * ( (*Bz2D)(i,j+1) - (*Bz2D)(i,j) );
		}
    }
    
    // Electric field Ey^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
		for (unsigned int j=0 ; j<ny_d ; j++) {
			(*Ey2D)(i,j) += -dt*(*Jy2D)(i,j) - dt_ov_dx * ( (*Bz2D)(i+1,j) - (*Bz2D)(i,j) );
		}
    }
    
    // Electric field Ez^(p,p)
    for (unsigned int i=0 ;  i<nx_p ; i++) {
		for (unsigned int j=0 ; j<ny_p ; j++) {
			(*Ez2D)(i,j) += -dt*(*Jz2D)(i,j)
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
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::applyEMBoundaryConditions(double time_dual, SmileiMPI* smpi)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    
    // Static cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D = static_cast<Field2D*>(Bx_);
    Field2D* By2D = static_cast<Field2D*>(By_);
    Field2D* Bz2D = static_cast<Field2D*>(Bz_);
    
    
    // ! \todo Transverse profile & incidence angle is not yet introduced (MG)
    // -----------------------------------------
    // Laser temporal profile
    // -----------------------------------------
    double byW=0.0, bzW=0.0, byE=0.0, bzE=0.0;
    
    for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
        
		if (laser_[ilaser]->laser_struct.angle == 0){
            // Incident field (west boundary)
            byW += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);
            bzW += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
        } else if (laser_[ilaser]->laser_struct.angle == 180){
            // Incident field (east boundary)
            byE += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);
            bzE += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
        } else {
            ERROR("Angle not yet implemented for laser " << ilaser);
        }

    }//ilaser
    
	
    // If !periodic
    //   BC : Bx(i=0...nx_p, 0) & Bx(i=0...nx_p, ny_d-1)
    //   BC : Bz(i=0...nx_d-1, 0) & Bz(i=0...nx_d-1, ny_d-1)
    // else
    //   Ez(i,-1)/Ez(i,ny_d-1) defined -> OK
    //   Ex(i,-1)/Ex(i,ny_d-1) defined -> OK
	
    // BC : By(0, j=0...ny_p  ) & By(nx_d-1, j=0...ny_p  )	    -> TO DO
    // BC : Bz(O, j=0...ny_d-1) & Bz(nx_d-1, j=0...ny_d-1)		-> TO DO
	
	
    // -----------------------------------------
    // Silver-Mueller boundary conditions (West)
    // -----------------------------------------
    if ( smpi2D->isWester() ) {
		//MESSAGE( smpi->getRank() << " is wester" );
		// for By^(d,p)
		for (unsigned int j=0 ; j<ny_p ; j++) {
			(*By2D)(0,j) = Alpha_SM_W   * (*Ez2D)(0,j)
			+              Beta_SM_W    * (*By2D)(1,j)
			+              Gamma_SM_W   * byW
			+              Delta_SM_W   * (*Bx2D)(0,j+1)
			+              Epsilon_SM_W * (*Bx2D)(0,j);
		}
		// for Bz^(d,d)
		for (unsigned int j=0 ; j<ny_d ; j++) {
			(*Bz2D)(0,j) = -Alpha_SM_W * (*Ey2D)(0,j)
			+               Beta_SM_W  * (*Bz2D)(1,j)
			+               Gamma_SM_W * bzW;
		}
		
		/*		// Correction on unused extreme ghost cells : put the fields to 0
		 // --------------------------------------------------------------
		 //! \todo{Alloc Fields only if necessary to not doing this correction}
		 for (unsigned int i=0 ; i<index_bc_min[0] ; i++) {
		 // for Ey^(p,d), Bx^(p,d), Bz^(d,d)
		 for (unsigned int j=0 ; j<ny_d ; j++) {
		 (*Ey2D)(i,j)=0.0;
		 (*Bx2D)(i,j)=0.0;
		 (*Bz2D)(i,j)=0.0;
		 }
		 // for Ex^(d,p), Ez^(p,p), By^(d,p)
		 for (unsigned int j=0 ; j<ny_p ; j++) {
		 (*Ex2D)(i,j)=0.0;
		 (*Ez2D)(i,j)=0.0;
		 (*By2D)(i,j)=0.0;
		 }
		 }
		 */        
        //MESSAGE( "WEST BC ==> By_inc, By " << byW << " " << (*By2D)(index_bc_min[0],250) );
        
    }//if West
    
    // -----------------------------------------
    // Silver-Mueller boundary conditions (East)
    // -----------------------------------------
    if ( smpi2D->isEaster() ) {
		//MESSAGE( smpi->getRank() << " is easter" );
        //MESSAGE("Here once");
		// for By^(d,p)
		for (unsigned int j=0 ; j<ny_p ; j++) {
			(*By2D)(nx_d-1,j) = Alpha_SM_E   * (*Ez2D)(nx_p-1,j)
			+                   Beta_SM_E    * (*By2D)(nx_d-2,j)
			+                   Gamma_SM_E   * byE
			+                   Delta_SM_E   * (*Bx2D)(nx_p-1,j+1) // Check x-index
			+                   Epsilon_SM_E * (*Bx2D)(nx_p-1,j);
		}
		// for Bz^(d,d)
		for (unsigned int j=0 ; j<ny_d ; j++) {
			(*Bz2D)(nx_d-1,j) = -Alpha_SM_E * (*Ey2D)(nx_p-1,j)
			+                    Beta_SM_E  * (*Bz2D)(nx_d-2,j)
			+                    Gamma_SM_E * bzE;
		}
		
		/*		// Correction on unused extreme ghost cells : put the fields to 0
		 // --------------------------------------------------------------
		 // for Ex^(d,p), By^(d,p), Bz^(d,d)
		 for (unsigned int i=index_bc_max[0]+1 ; i<nx_d ; i++) {
		 for (unsigned int j=0 ; j<ny_p ; j++) {
		 (*Ex2D)(i,j) = 0.0;
		 (*By2D)(i,j) = 0.0;
		 }
		 for (unsigned int j=0 ; j<ny_d ; j++) {
		 (*Bz2D)(i,j) = 0.0;
		 }
		 }
		 // for Ey^(p,d), Ez^(p,p), Bx^(p,d)
		 for (unsigned int i=index_bc_max[0] ; i<nx_p ; i++) {
		 for (unsigned int j=0 ; j<ny_d ; j++) {
		 (*Ey2D)(i,j)=0.0;
		 (*Bx2D)(i,j)=0.0;
		 }
		 }
		 // for Ez^(p,p)
		 for (unsigned int i=index_bc_max[0] ; i<nx_p ; i++) {
		 for (unsigned int j=0 ; j<ny_p ; j++) {
		 (*Ez2D)(i,j)=0.0;
		 }
		 }
		 */
        //MESSAGE( "EAST BC ==> By_inc, By " << byE << " " << (*By2D)(index_bc_max[0],250) );
		
    }//if East
    
}// END applyEMBoundaryConditions



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
    Field2D* rho2D_o = static_cast<Field2D*>(rho_o);
	
    // Save current charge density as old charge density
    for (unsigned int i=0 ; i<nx_p ; i++) {
		for (unsigned int j=0 ; j<ny_p ; j++) {
			(*rho2D_o)(i,j) = (*rho2D)(i,j);
		}
    }
    
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
    
    
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++){
        Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
        Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
        Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
        Field2D* rho2D_s = static_cast<Field2D*>(rho_s[ispec]);
        
        // Charge density rho^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*rho2D_s)(i,j) = 0.0;
            }
        }
        
        // Current Jx^(d,p) to 0
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jx2D_s)(i,j) = 0.0;
            }
        }
        
        // Current Jy^(p,d) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Jy2D_s)(i,j) = 0.0;
            }
        }
        
        // Current Jz^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jz2D_s)(i,j) = 0.0;
            }
        }
        
    }//END loop on species ispec
    
}//END restartRhoJ
