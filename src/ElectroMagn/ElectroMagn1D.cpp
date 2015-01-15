#include "ElectroMagn1D.h"

#include <cmath>

#include <sstream>
#include <string>
#include <iostream>

#include "PicParams.h"
#include "Field1D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::ElectroMagn1D(PicParams &params,  LaserParams &laser_params, SmileiMPI* smpi)
: ElectroMagn(params, laser_params, smpi),
isWestern(smpi->isWestern()),
isEastern(smpi->isEastern())
{
    // local dt to store
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
    int process_coord_x = smpi1D->getProcCoord(0);
    
    oversize_ = oversize[0];
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step
    dx       = cell_length[0];
    dt_ov_dx = timestep/cell_length[0];
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // Electromagnetic fields
    // ----------------------
    // number of nodes of the primal-grid
    nx_p = n_space[0]+1 + 2*oversize[0];
    // number of nodes of the dual-grid
    nx_d = n_space[0]+2 + 2*oversize[0];
    // dimPrim/dimDual = nx_p/nx_d
    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    for (size_t i=0 ; i<nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    
    // Allocation of the EM fields
    Ex_  = new Field1D(dimPrim, 0, false, "Ex");
    Ey_  = new Field1D(dimPrim, 1, false, "Ey");
    Ez_  = new Field1D(dimPrim, 2, false, "Ez");
    Bx_  = new Field1D(dimPrim, 0, true,  "Bx");
    By_  = new Field1D(dimPrim, 1, true,  "By");
    Bz_  = new Field1D(dimPrim, 2, true,  "Bz");
    Bx_m = new Field1D(dimPrim, 0, true,  "Bx_m");
    By_m = new Field1D(dimPrim, 1, true,  "By_m");
    Bz_m = new Field1D(dimPrim, 2, true,  "Bz_m");
    
    // Allocation of time-averaged EM fields
    Ex_avg  = new Field1D(dimPrim, 0, false, "Ex_avg");
    Ey_avg  = new Field1D(dimPrim, 1, false, "Ey_avg");
    Ez_avg  = new Field1D(dimPrim, 2, false, "Ez_avg");
    Bx_avg  = new Field1D(dimPrim, 0, true,  "Bx_avg");
    By_avg  = new Field1D(dimPrim, 1, true,  "By_avg");
    Bz_avg  = new Field1D(dimPrim, 2, true,  "Bz_avg");
    
    // Total charge currents and densities
    Jx_   = new Field1D(dimPrim, 0, false, "Jx");
    Jy_   = new Field1D(dimPrim, 1, false, "Jy");
    Jz_   = new Field1D(dimPrim, 2, false, "Jz");
    rho_  = new Field1D(dimPrim, "Rho" );
    
    // Charge currents currents and density for each species

    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field1D(dimPrim, 0, false, ("Jx_"+params.species_param[ispec].species_type).c_str());
        Jy_s[ispec]  = new Field1D(dimPrim, 1, false, ("Jy_"+params.species_param[ispec].species_type).c_str());
        Jz_s[ispec]  = new Field1D(dimPrim, 2, false, ("Jz_"+params.species_param[ispec].species_type).c_str());
        rho_s[ispec] = new Field1D(dimPrim, ("Rho_"+params.species_param[ispec].species_type).c_str());
        Jx_s[n_species+ispec]  = new Field1D(dimPrim, 0, false, ("Jx_"+params.species_param[ispec].species_type).c_str());
        Jy_s[n_species+ispec]  = new Field1D(dimPrim, 1, false, ("Jy_"+params.species_param[ispec].species_type).c_str());
        Jz_s[n_species+ispec]  = new Field1D(dimPrim, 2, false, ("Jz_"+params.species_param[ispec].species_type).c_str());
        rho_s[n_species+ispec] = new Field1D(dimPrim, ("Rho_"+params.species_param[ispec].species_type).c_str());
    }
    
//    ostringstream file_name("");
//    for (unsigned int ispec=0; ispec<n_species; ispec++) {
//        file_name.str("");
//        file_name << "Jx_s" << ispec;
//        Jx_s[ispec]  = new Field1D(dimPrim, 0, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "Jy_s" << ispec;
//        Jy_s[ispec]  = new Field1D(dimPrim, 1, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "Jz_s" << ispec;
//        Jz_s[ispec]  = new Field1D(dimPrim, 2, false, file_name.str().c_str());
//        file_name.str("");
//        file_name << "rho_s" << ispec;
//        rho_s[ispec] = new Field1D(dimPrim, file_name.str().c_str());
//    }
    
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (size_t i=0 ; i<nDim_field ; i++) {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }
    
    
    // Define limits of non duplicated elements
    // (by construction 1 (prim) or 2 (dual) elements shared between per MPI process)
    // istart
    for (unsigned int i=0 ; i<3 ; i++)
        for (unsigned int isDual=0 ; isDual<2 ; isDual++)
            istart[i][isDual] = 0;
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        for (unsigned int isDual=0 ; isDual<2 ; isDual++) {
            istart[i][isDual] = oversize[i];
            if (smpi1D->getProcCoord(i)!=0) istart[i][isDual]+=1;
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
            if ( smpi1D->getNbrOfProcs(i)!=1 ) {
                
                if ( ( !isDual ) && (smpi1D->getProcCoord(i)!=0) )
                    bufsize[i][isDual]--;
                else if  (isDual) {
                    bufsize[i][isDual]--;
                    if ( (smpi1D->getProcCoord(i)!=0) && (smpi1D->getProcCoord(i)!=smpi1D->getNbrOfProcs(i)-1) )
                        bufsize[i][isDual]--;
                }
                
            } // if ( smpi1D->getNbrOfProcs(i)!=1 )
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field
    
    
}//END constructor Electromagn1D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::~ElectroMagn1D()
{
}


// ---------------------------------------------------------------------------------------------------------------------
// Solve Poisson
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solvePoisson(SmileiMPI* smpi)
{
    
    SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
    //int process_coord_x = smpi1D->getProcCoord(0);
    
    unsigned int iteration_max  = 100000;
    double       error_max      = 1.e-14; //!\todo Check what should be used to relate to machine precision
    
    double       dx_sq          = dx*dx;
    unsigned int nx_p_global    = smpi1D->n_space_global[0] + 1;
    unsigned int smilei_sz      = smpi1D->smilei_sz;
    unsigned int smilei_rk      = smpi1D->smilei_rk;
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // ----------------------------------------------------------------------------------
    vector<unsigned int> index_min_p(1);
    vector<unsigned int> index_min_d(1);
    vector<unsigned int> index_max_p(1);
    vector<unsigned int> index_max_d(1);
    index_min_p[0] = oversize_;
    index_min_d[0] = oversize_;
    index_max_p[0] = nx_p - 2 - oversize_;
    index_max_d[0] = nx_d - 2 - oversize_;
    if (smpi1D->isWestern()) {
        index_min_p[0] = 0;
        index_min_d[0] = 0;
    }
    if (smpi1D->isEastern()) {
        index_max_p[0] = nx_p-1;
        index_max_d[0] = nx_d-1;
    }
    
    Field1D* Ex1D  = static_cast<Field1D*>(Ex_);
    Field1D* rho1D = static_cast<Field1D*>(rho_);
    
    double pW=0.0;
    double pE=0.0;
    
    Field1D phi(dimPrim);    // scalar potential
    Field1D r(dimPrim);      // residual vector
    Field1D p(dimPrim);      // direction vector
    Field1D Ap(dimPrim);     // A*p vector
    
    
    // -------------------------------
    // Initialization of the variables
    // -------------------------------
    DEBUG(1,"Initialize variables");
    
    unsigned int iteration = 0;
    
    // phi: scalar potential, r: residual and p: direction
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        phi(i)   = 0.0;
        r(i)     = -dx_sq * (*rho1D)(i);
        p(i)     = r(i);
    }
    
    // scalar product of the residual (first local then reduction over all procs)
    //!\todo Generalize parallel computation of scalar product as a method for FieldND
    double rnew_dot_rnew       = 0.0;
    double rnew_dot_rnew_local = 0.0;
    //for (unsigned int i=index_bc_min[0] ; i<index_bc_max[0]-1 ; i++) rnew_dot_rnew_local += r(i)*r(i);
    for (unsigned int i=index_min_p[0] ; i<=index_max_p[0] ; i++) rnew_dot_rnew_local += r(i)*r(i);
    MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // control parameter
    //!\todo Check control parameter for // case
    double ctrl = rnew_dot_rnew/(double)(nx_p_global);
    
    
    // ---------------------------------------------------------
    // Starting iterative loop for the conjugate gradient method
    // ---------------------------------------------------------
    DEBUG(1,"Starting iterative loop");
    while ( (iteration<5) || ((ctrl > error_max) && (iteration<iteration_max)) ) {
        //  NB: forced at least 5 iterations to avoid problem when initial fields are small
        
        iteration++;
        DEBUG(5,"iteration " << iteration << " started with control parameter ctrl = " << ctrl*1.e14 << " x 1e-14");
        
        double r_dot_r = rnew_dot_rnew;
        
        // vector product Ap = A*p
        for (unsigned int i=1 ; i<dimPrim[0]-1 ; i++) Ap(i) = p(i-1) - 2.0*p(i) + p(i+1);
        
        // apply BC on Ap
        if (smpi1D->isWestern()) Ap(0)      = pW        - 2.0*p(0)      + p(1);
        if (smpi1D->isEastern()) Ap(nx_p-1) = p(nx_p-2) - 2.0*p(nx_p-1) + pE;
        smpi1D->exchangeField(&Ap);
        
        // scalar product p.Ap
        double p_dot_Ap = 0.0;
        double p_dot_Ap_local = 0.0;
        //for (unsigned int i=index_bc_min[0] ; i<index_bc_max[0]-1 ; i++) p_dot_Ap_local += p(i)*Ap(i);
        for (unsigned int i=index_min_p[0] ; i<=index_max_p[0] ; i++) p_dot_Ap_local += p(i)*Ap(i);
        MPI_Allreduce(&p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // compute new potential and residual
        double alpha_k = r_dot_r/p_dot_Ap;
        for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
            phi(i) += alpha_k * p(i);
            r(i)   -= alpha_k * Ap(i);
        }
        
        // compute new residual norm
        rnew_dot_rnew       = 0.0;
        rnew_dot_rnew_local = 0.0;
        //for (unsigned int i=index_bc_min[0] ; i<index_bc_max[0]-1 ; i++) rnew_dot_rnew_local += r(i)*r(i);
        for (unsigned int i=index_min_p[0] ; i<=index_max_p[0] ; i++) rnew_dot_rnew_local += r(i)*r(i);
        MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // compute new direction
        double beta_k = rnew_dot_rnew/r_dot_r;
        for (unsigned int i=0 ; i<dimPrim[0] ; i++)  p(i) = r(i) + beta_k * p(i);
        
        // compute control parameter
        ctrl = rnew_dot_rnew/(double)(nx_p_global);
        
    }//End of the iterative loop
    
    
    // Status of the solver convergence
    if (iteration == iteration_max) {
        if (smpi1D->isMaster())
            WARNING("Poisson solver did not converge: reached maximum iteration number: " << iteration
                    << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    else {
        if (smpi1D->isMaster())
            MESSAGE(1,"Poisson solver converged at iteration: " << iteration
                    << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    
    // ----------------------------------
    // Compute the electrostatic field Ex
    // ----------------------------------
    
    for (unsigned int i=1; i<nx_d-1; i++) (*Ex1D)(i) = (phi(i-1)-phi(i))/dx;
    
    
    // BC on Ex
    smpi1D->exchangeField(Ex1D);
    if (smpi1D->isWestern()) (*Ex1D)(0)      = (*Ex1D)(1)      - dx*(*rho1D)(0);
    if (smpi1D->isEastern()) (*Ex1D)(nx_d-1) = (*Ex1D)(nx_d-2) + dx*(*rho1D)(nx_p-1);
    
    
    // Find field to be added to ensure BC: Ex_West = -Ex_East
    
    double Ex_West = 0.0;
    double Ex_East = 0.0;
    
    unsigned int rankWest = smpi1D->extrem_ranks[0][0];
    if (smpi1D->isWestern()) {
        if (smilei_rk != smpi1D->extrem_ranks[0][0]) ERROR("western process not well defined");
        Ex_West = (*Ex1D)(index_bc_min[0]);
    }
    MPI_Bcast(&Ex_West, 1, MPI_DOUBLE, rankWest, MPI_COMM_WORLD);
    
    unsigned int rankEast = smpi1D->extrem_ranks[0][1];
    if (smpi1D->isEastern()) {
        if (smilei_rk != smpi1D->extrem_ranks[0][1]) ERROR("eastern process not well defined");
        Ex_East = (*Ex1D)(index_bc_max[0]);
    }
    MPI_Bcast(&Ex_East, 1, MPI_DOUBLE, rankEast, MPI_COMM_WORLD);
    
    
    double Ex_Add = -0.5*(Ex_West+Ex_East);
    
    // Center the electrostatic field
    for (unsigned int i=0; i<nx_d; i++) (*Ex1D)(i) += Ex_Add;
    
    
    // Compute error on the Poisson equation
    double deltaPoisson_max = 0.0;
    int i_deltaPoisson_max  = -1;
    for (unsigned int i=0; i<nx_p; i++) {
        double deltaPoisson = abs( ((*Ex1D)(i+1)-(*Ex1D)(i))/dx - (*rho1D)(i) );
        if (deltaPoisson > deltaPoisson_max) {
            deltaPoisson_max   = deltaPoisson;
            i_deltaPoisson_max = i;
        }
    }
    
    //!\todo Reduce to find global max
    if (smpi1D->isMaster())
        MESSAGE(1,"Poisson equation solved. Maximum error = " << deltaPoisson_max << " at i= " << i_deltaPoisson_max);
    
    
    
}//END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::saveMagneticFields()
{
    // Static cast of the fields
    Field1D* Bx1D   = static_cast<Field1D*>(Bx_);
    Field1D* By1D   = static_cast<Field1D*>(By_);
    Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
    Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
    Field1D* By1D_m = static_cast<Field1D*>(By_m);
    Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);
    
    // for Bx^(p)
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*Bx1D_m)(i)=(*Bx1D)(i);
    }
    //for By^(d) & Bz^(d)
    for (unsigned int i=0 ; i<dimDual[0] ; i++) {
        (*By1D_m)(i) = (*By1D)(i);
        (*Bz1D_m)(i) = (*Bz1D)(i);
    }
    
}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwellAmpere()
{
    
    Field1D* Ex1D = static_cast<Field1D*>(Ex_);
    Field1D* Ey1D = static_cast<Field1D*>(Ey_);
    Field1D* Ez1D = static_cast<Field1D*>(Ez_);
    Field1D* By1D = static_cast<Field1D*>(By_);
    Field1D* Bz1D = static_cast<Field1D*>(Bz_);
    Field1D* Jx1D = static_cast<Field1D*>(Jx_);
    Field1D* Jy1D = static_cast<Field1D*>(Jy_);
    Field1D* Jz1D = static_cast<Field1D*>(Jz_);
    
    // --------------------
    // Solve Maxwell-Ampere
    // --------------------
    // Calculate the electrostatic field ex on the dual grid
    //for (unsigned int ix=0 ; ix<nx_d ; ix++){
    for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
        (*Ex1D)(ix)= (*Ex1D)(ix) - timestep* (*Jx1D)(ix) ;
    }
    // Transverse fields ey, ez  are defined on the primal grid
    //for (unsigned int ix=0 ; ix<nx_p ; ix++) {
    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
        (*Ey1D)(ix)= (*Ey1D)(ix) - dt_ov_dx * ( (*Bz1D)(ix+1) - (*Bz1D)(ix)) - timestep * (*Jy1D)(ix) ;
        (*Ez1D)(ix)= (*Ez1D)(ix) + dt_ov_dx * ( (*By1D)(ix+1) - (*By1D)(ix)) - timestep * (*Jz1D)(ix) ;
    }
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwellFaraday()
{
    Field1D* Ey1D   = static_cast<Field1D*>(Ey_);
    Field1D* Ez1D   = static_cast<Field1D*>(Ez_);
    Field1D* By1D   = static_cast<Field1D*>(By_);
    Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
    
    // ---------------------
    // Solve Maxwell-Faraday
    // ---------------------
    // NB: bx is given in 1d and defined when initializing the fields (here put to 0)
    // Transverse fields  by & bz are defined on the dual grid
    //for (unsigned int ix=1 ; ix<nx_p ; ix++) {
    for (unsigned int ix=1 ; ix<dimDual[0]-1 ; ix++) {
        (*By1D)(ix)= (*By1D)(ix) + dt_ov_dx * ( (*Ez1D)(ix) - (*Ez1D)(ix-1)) ;
        (*Bz1D)(ix)= (*Bz1D)(ix) - dt_ov_dx * ( (*Ey1D)(ix) - (*Ey1D)(ix-1)) ;
    }
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::centerMagneticFields()
{
    // Static cast of the fields
    Field1D* Bx1D   = static_cast<Field1D*>(Bx_);
    Field1D* By1D   = static_cast<Field1D*>(By_);
    Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
    Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
    Field1D* By1D_m = static_cast<Field1D*>(By_m);
    Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);
    
    // for Bx^(p)
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*Bx1D_m)(i) = ( (*Bx1D)(i)+ (*Bx1D_m)(i))*0.5 ;
    }
    
    // for By^(d) & Bz^(d)
    for (unsigned int i=0 ; i<dimDual[0] ; i++) {
        (*By1D_m)(i)= ((*By1D)(i)+(*By1D_m)(i))*0.5 ;
        (*Bz1D_m)(i)= ((*Bz1D)(i)+(*Bz1D_m)(i))*0.5 ;
    }
    
}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::incrementAvgFields(unsigned int time_step, unsigned int ntime_step_avg)
{
    // Static cast of the fields
    Field1D* Ex1D     = static_cast<Field1D*>(Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(Bz_m);
    Field1D* Ex1D_avg = static_cast<Field1D*>(Ex_avg);
    Field1D* Ey1D_avg = static_cast<Field1D*>(Ey_avg);
    Field1D* Ez1D_avg = static_cast<Field1D*>(Ez_avg);
    Field1D* Bx1D_avg = static_cast<Field1D*>(Bx_avg);
    Field1D* By1D_avg = static_cast<Field1D*>(By_avg);
    Field1D* Bz1D_avg = static_cast<Field1D*>(Bz_avg);
    
    // reset the averaged fields for (time_step-1)%ntime_step_avg == 0
    if ( (time_step-1)%ntime_step_avg==0 ){
        Ex1D_avg->put_to(0.0);
        Ey1D_avg->put_to(0.0);
        Ez1D_avg->put_to(0.0);
        Bx1D_avg->put_to(0.0);
        By1D_avg->put_to(0.0);
        Bz1D_avg->put_to(0.0);
    }
    
    // for Ey^(p), Ez^(p) & Bx^(p)
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*Ey1D_avg)(i) += (*Ey1D)(i);
        (*Ez1D_avg)(i) += (*Ez1D)(i);
        (*Bx1D_avg)(i) += (*Bx1D_m)(i);
    }
    
    // for Ex^(d), By^(d) & Bz^(d)
    for (unsigned int i=0 ; i<dimDual[0] ; i++) {
        (*Ex1D_avg)(i) += (*Ex1D)(i);
        (*By1D_avg)(i) += (*By1D_m)(i);
        (*Bz1D_avg)(i) += (*Bz1D_m)(i);
    }
    
}//END incrementAvgFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge density and transverse currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::restartRhoJ()
{
    Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
    Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
    Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
    Field1D* rho1D   = static_cast<Field1D*>(rho_);
    
    // --------------------------
    // Total currents and density
    // --------------------------
    
    // put longitudinal current to zero on the dual grid
    for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
        (*Jx1D)(ix)    = 0.0;
    }
    
    // all fields are defined on the primal grid
    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
        (*rho1D)(ix)   = 0.0;
        (*Jy1D)(ix)    = 0.0;
        (*Jz1D)(ix)    = 0.0;
    }
}    
void ElectroMagn1D::restartRhoJs(unsigned int ispec, bool currents)
{
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
    Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
    Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
    Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);
    Field1D* Jx1D_s2  = static_cast<Field1D*>(Jx_s[n_species+ispec]);
    Field1D* Jy1D_s2  = static_cast<Field1D*>(Jy_s[n_species+ispec]);
    Field1D* Jz1D_s2  = static_cast<Field1D*>(Jz_s[n_species+ispec]);
    Field1D* rho1D_s2 = static_cast<Field1D*>(rho_s[n_species+ispec]);

    #pragma omp for schedule(static) 
    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
        (*rho1D_s)(ix) = 0.0;
        (*rho1D_s2)(ix) = 0.0;
    }
    if (currents){
        // put longitudinal current to zero on the dual grid
        #pragma omp for schedule(static) 
        for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
            (*Jx1D_s)(ix)  = 0.0;
            (*Jx1D_s2)(ix)  = 0.0;
        }
        #pragma omp for schedule(static) 
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
        // all fields are defined on the primal grid
            (*Jy1D_s)(ix)  = 0.0;
            (*Jz1D_s)(ix)  = 0.0;
            (*Jy1D_s2)(ix)  = 0.0;
            (*Jz1D_s2)(ix)  = 0.0;
        }
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::computeTotalRhoJ()
{
    Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
    Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
    Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
    Field1D* rho1D   = static_cast<Field1D*>(rho_);
    
    for (unsigned int ispec=0; ispec<n_species*2; ispec++) {
        Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
        Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
        Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
        Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);
        
        #pragma omp for schedule(static) nowait
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
            (*Jx1D)(ix)  += (*Jx1D_s)(ix);
            (*Jy1D)(ix)  += (*Jy1D_s)(ix);
            (*Jz1D)(ix)  += (*Jz1D_s)(ix);
            (*rho1D)(ix) += (*rho1D_s)(ix);
        }

        #pragma omp single
        {
            (*Jx1D)(dimPrim[0])  += (*Jx1D_s)(dimPrim[0]);
        }
    }//END loop on species ispec
}
void ElectroMagn1D::addToGlobalRho(int ispec, int clrw)
{
 }
void ElectroMagn1D::computeTotalRhoJs(int clrw)
{
 }

// ---------------------------------------------------------------------------------------------------------------------
// Gather the total density and currents for species on a single array instead of twin arrays.
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::sumtwins()
{
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
        Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
        Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
        Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);
        Field1D* Jx1D_s2  = static_cast<Field1D*>(Jx_s[n_species+ispec]);
        Field1D* Jy1D_s2  = static_cast<Field1D*>(Jy_s[n_species+ispec]);
        Field1D* Jz1D_s2  = static_cast<Field1D*>(Jz_s[n_species+ispec]);
        Field1D* rho1D_s2 = static_cast<Field1D*>(rho_s[n_species+ispec]);
        
        // Charge density rho^(p) 
        for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
            (*rho1D_s)(i) += (*rho1D_s2)(i);
            (*Jy1D_s)(i) += (*Jy1D_s2)(i);
            (*Jz1D_s)(i) += (*Jz1D_s2)(i);
        }
        
        // Current Jx^(d)
        for (unsigned int i=0 ; i<dimDual[0] ; i++) {
            (*Jx1D_s)(i) += (*Jx1D_s2)(i);
        }
        
    }//END loop on species ispec
    
}//END sumtwins

void ElectroMagn1D::computePoynting() {
    if (isWestern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)];
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBy=istart[0][By_m->isDual(0)];
        
        poynting_inst[0][0]=0.5*timestep*((*Ey_)(iEy) * ((*Bz_m)(iBz) + (*Bz_m)(iBz+1)) -
                                          (*Ez_)(iEz) * ((*By_m)(iBy) + (*By_m)(iBy+1)));
        poynting[0][0] += poynting_inst[0][0];
    } 
    if (isEastern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)]  + bufsize[0][Ey_->isDual(0)]-1;
        unsigned int iBz=istart[0][Bz_m->isDual(0)] + bufsize[0][Bz_m->isDual(0)]-1;
        unsigned int iEz=istart[0][Ez_->isDual(0)]  + bufsize[0][Ez_->isDual(0)]-1;
        unsigned int iBy=istart[0][By_m->isDual(0)] + bufsize[0][By_m->isDual(0)]-1;
        
        poynting_inst[1][0]=0.5*timestep*((*Ey_)(iEy) * ((*Bz_m)(iBz-1) + (*Bz_m)(iBz)) - 
                                          (*Ez_)(iEz) * ((*By_m)(iBy-1) + (*By_m)(iBy)));
        poynting[1][0] -= poynting_inst[1][0];
        
    }    
}
