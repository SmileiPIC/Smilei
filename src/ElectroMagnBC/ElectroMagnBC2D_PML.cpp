#include "ElectroMagnBC2D_PML.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field2D.h"
#include "Tools.h"
#include "Laser.h"

#include "SolverFactory.h"

using namespace std;

ElectroMagnBC2D_PML::ElectroMagnBC2D_PML( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC2D( params, patch, i_boundary )
{

    pml_solver_ = SolverFactory::createPML( params );
    if (params.maxwell_sol=="Yee"){
        nsolver=2;
        //MESSAGE("FDTD scheme in PML region : Yee.");
    }
    else if (params.maxwell_sol=="Bouchard"){
        nsolver=4;
        //MESSAGE("FDTD scheme in PML region : Bouchard.");
    }
    else {
        WARNING("The solver you use in the main domain is not the same as in the PML region. FDTD scheme in PML region : Yee.");
        nsolver=2;
    }

    if( patch->isBoundary( i_boundary_ ) ) {

        int iDim = i_boundary / 2;
        int min_or_max = i_boundary_ % 2;

        domain_oversize_x =  patch->oversize[0] ;
        domain_oversize_y =  patch->oversize[1] ;

        if (patch->isXmin() ) {//&& i_boundary_ == 0 ) {
            ncells_pml_xmin = params.number_of_pml_cells[0][0];
            ncells_pml_domain_xmin = ncells_pml_xmin + 1*patch->oversize[0] + nsolver/2;
            domain_oversize_x = patch->oversize[0] ;
        }
        else {
            ncells_pml_xmin = 0;
            ncells_pml_domain_xmin = 0;
        }
        if (patch->isXmax() ) {//&& i_boundary_ == 1 ) {
            ncells_pml_xmax = params.number_of_pml_cells[0][1];
            ncells_pml_domain_xmax = ncells_pml_xmax + 1*patch->oversize[0] + nsolver/2;
            domain_oversize_x = patch->oversize[0] ;
        }
        else {
            ncells_pml_xmax = 0;
            ncells_pml_domain_xmax = 0;
        }
        if (patch->isYmin() ) {//&& i_boundary_ == 2 ) {
            ncells_pml_ymin = params.number_of_pml_cells[1][0];
            ncells_pml_domain_ymin = ncells_pml_ymin + 1*patch->oversize[1] + nsolver/2;
            domain_oversize_y = patch->oversize[1] ;
        }
        else {
            ncells_pml_ymin = 0;
            ncells_pml_domain_ymin = 0;
        }
        if (patch->isYmax() ) {//&& i_boundary_ == 3 ) {
            ncells_pml_ymax = params.number_of_pml_cells[1][1];
            ncells_pml_domain_ymax = ncells_pml_ymax + 1*patch->oversize[1] + nsolver/2;
            domain_oversize_y = patch->oversize[1] ;
        }
        else {
            ncells_pml_ymax = 0;
            ncells_pml_domain_ymax = 0;
        }

        ncells_pml = params.number_of_pml_cells[iDim][min_or_max];
        ncells_pml_domain = ncells_pml+1*patch->oversize[iDim] + nsolver/2;

        // Define min and max idx to exchange
        // the good data f(solver,oversize)
        if (min_or_max==0){
            // if min border : Exchange of data (for domain to pml-domain)
            // min2exchange <= i < max2exchange
            min2exchange = 1*nsolver/2 ;
            max2exchange = 2*nsolver/2 ;
            // Solver
            solvermin = nsolver/2 ;
            solvermax = ncells_pml_domain - patch->oversize[iDim] ;
        }
        else if (min_or_max==1){
            // if max border : Exchange of data (for domain to pml-domain)
            // min2exchange <= i < max2exchange
            min2exchange = 1*nsolver/2 ;
            max2exchange = 2*nsolver/2 ;
            // Solver
            solvermin = patch->oversize[iDim] + nsolver/2 - nsolver/2 + 1 ;
            solvermax = ncells_pml_domain-nsolver/2 ;
        }

        if (ncells_pml==0){
            ERROR("PML domain have to be >0 cells in thickness");
        }
        
        if( iDim == 0 ) {
            dimPrim = { (unsigned int)ncells_pml_domain, patch->size_[1]+1+2*patch->oversize[1] };
        } else {
            dimPrim = { patch->size_[0]+1+2*patch->oversize[0] + ncells_pml_xmin + ncells_pml_xmax, (unsigned int)ncells_pml_domain };
            ypml_size_in_x = dimPrim[0];
        }

        startpml = patch->oversize[iDim]+nsolver/2;

        int ncells_pml_min[1];
        ncells_pml_min[0] = ncells_pml_xmin;
        int ncells_pml_max[1];
        ncells_pml_max[0] = ncells_pml_xmax;

        pml_solver_->setDomainSizeAndCoefficients( iDim, min_or_max, dimPrim, ncells_pml_domain, startpml, ncells_pml_min, ncells_pml_max, patch );

        std::string si_boundary = std::to_string(i_boundary_);

        Ex_ = new Field2D( dimPrim, 0, false, "Ex_pml"+si_boundary );
        Ey_ = new Field2D( dimPrim, 1, false, "Ey_pml"+si_boundary );
        Ez_ = new Field2D( dimPrim, 2, false, "Ez_pml"+si_boundary );
        Bx_ = new Field2D( dimPrim, 0, true , "Bx_pml"+si_boundary );
        By_ = new Field2D( dimPrim, 1, true , "By_pml"+si_boundary );
        Bz_ = new Field2D( dimPrim, 2, true , "Bz_pml"+si_boundary );
        Dx_ = new Field2D( dimPrim, 0, false, "Dx_pml"+si_boundary );
        Dy_ = new Field2D( dimPrim, 1, false, "Dy_pml"+si_boundary );
        Dz_ = new Field2D( dimPrim, 2, false, "Dz_pml"+si_boundary );
        Hx_ = new Field2D( dimPrim, 0, true , "Hx_pml"+si_boundary );
        Hy_ = new Field2D( dimPrim, 1, true , "Hy_pml"+si_boundary );
        Hz_ = new Field2D( dimPrim, 2, true , "Hz_pml"+si_boundary );

        //Laser parameter
        double pyKx, pyKy; //, pyKz;
        double kx, ky; //, kz;
        double Knorm;
        double omega = 1. ;

        factor_laser_space_time = 2.*dt_ov_d[0] ;

        // Xmin boundary
        pyKx = params.EM_BCs_k[0][0];
        pyKy = params.EM_BCs_k[0][1];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;

        factor_laser_angle_W = kx/Knorm;

        // Xmax boundary
        pyKx = params.EM_BCs_k[1][0];
        pyKy = params.EM_BCs_k[1][1];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;

        factor_laser_angle_E = kx/Knorm;

        // Ymin boundary
        pyKx = params.EM_BCs_k[2][0];
        pyKy = params.EM_BCs_k[2][1];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;

        factor_laser_angle_S = ky/Knorm;

        // Ymax boundary
        pyKx = params.EM_BCs_k[3][0];
        pyKy = params.EM_BCs_k[3][1];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;

        factor_laser_angle_N = ky/Knorm;
    }
}


ElectroMagnBC2D_PML::~ElectroMagnBC2D_PML()
{
    delete Ex_;
    delete Dx_;
    delete Hx_;
    delete Bx_;
    delete Ey_;
    delete Dy_;
    delete Hy_;
    delete By_;
    delete Ez_;
    delete Dz_;
    delete Hz_;
    delete Bz_;

    if (pml_solver_!=NULL) {
         delete pml_solver_;
    }
}


void ElectroMagnBC2D_PML::save_fields( Field *, Patch * )
{
}


void ElectroMagnBC2D_PML::disableExternalFields()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------

void ElectroMagnBC2D_PML::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    int iDim = i_boundary_ / 2;
    int min_or_max = (i_boundary_)%2;

    Field2D *Ex_domain = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey_domain = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez_domain = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx_domain = static_cast<Field2D *>( EMfields->Bx_ );
    Field2D *By_domain = static_cast<Field2D *>( EMfields->By_ );
    Field2D *Bz_domain = static_cast<Field2D *>( EMfields->Bz_ );
    
    if( ! patch->isBoundary( i_boundary_ ) ) return;
    
    if( i_boundary_ == 0 ) {

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        

        // 2. Exchange field PML <- Domain
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            // MESSAGE("Copy PML < Domain");
            // MESSAGE(ncells_pml_domain-domain_oversize_x-nsolver/2+i<<"<"<<i);
            for( unsigned int j = 0 ; j<n_d[1] ; j++ ) {
                (*Ey_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Ey_domain)(i,j);
                (*Dy_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Ey_domain)(i,j);
                (*Hx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Bx_domain)(i,j);
                (*Bx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Bx_domain)(i,j);
                (*Hz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Bz_domain)(i,j);
                (*Bz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Bz_domain)(i,j);
            }
            for( unsigned int j = 0 ; j<n_p[1] ; j++ ) {
                (*Hy_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*By_domain)(i,j);
                (*By_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*By_domain)(i,j);
                (*Ex_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Ex_domain)(i,j);
                (*Dx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Ex_domain)(i,j);
                (*Ez_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Ez_domain)(i,j);
                (*Dz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j) = (*Ez_domain)(i,j);
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);

        //Injecting a laser
        vector<double> yp( 1 );
        for( unsigned int j=patch->isYmin() ; j<n_p[1]-patch->isYmax() ; j++ ) {

            double byW = 0.;
            yp[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*d[1];

            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                byW += vecLaser[ilaser]->getAmplitude0( yp, time_dual, j, 0 );
            }
            (*Hy_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_W*byW ;
            (*By_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_W*byW ;
        }

        vector<double> yd( 1 );
        for( unsigned int j=patch->isYmin() ; j<n_d[1]-patch->isYmax() ; j++ ) {

            double bzW = 0.;
            yd[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*d[1];

            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bzW += vecLaser[ilaser]->getAmplitude1( yd, time_dual, j, 0 );
            }
            (*Hz_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_W*bzW ;
            (*Bz_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_W*bzW ;
        }

        // 4. Exchange PML -> Domain
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for( unsigned int j = 0 ; j<n_p[1] ; j++ ) {
                (*Ez_domain)(i,j) = (*Ez_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j);
            }
            for( unsigned int j = 0 ; j<n_d[1] ; j++ ) {
                (*Ey_domain)(i,j) = (*Ey_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j);
                (*Bx_domain)(i,j) = (*Hx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j);
            }
        }
        // Duals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for( unsigned int j = 0 ; j<n_p[1] ; j++ ) {
                (*Ex_domain)(i,j) = (*Ex_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j);
                (*By_domain)(i,j) = (*Hy_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j);
            }
            for( unsigned int j = 0 ; j<n_d[1] ; j++ ) {
                (*Bz_domain)(i,j) = (*Hz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j);
            }
        }
    }
    else if( i_boundary_ == 1 ) {

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 2. Exchange field Domain -> PML
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            for( unsigned int j = 0 ; j<n_d[1] ; j++ ) {
                (*Ey_)(domain_oversize_x+nsolver/2-i,j) = (*Ey_domain)(n_p[0]-i,j);
                (*Dy_)(domain_oversize_x+nsolver/2-i,j) = (*Ey_domain)(n_p[0]-i,j);
                (*Hx_)(domain_oversize_x+nsolver/2-i,j) = (*Bx_domain)(n_p[0]-i,j);
                (*Bx_)(domain_oversize_x+nsolver/2-i,j) = (*Bx_domain)(n_p[0]-i,j);
                (*Hz_)(domain_oversize_x+nsolver/2-i,j) = (*Bz_domain)(n_p[0]-i,j);
                (*Bz_)(domain_oversize_x+nsolver/2-i,j) = (*Bz_domain)(n_p[0]-i,j);
            }
            for( unsigned int j = 0 ; j<n_p[1] ; j++ ) {
                (*Hy_)(domain_oversize_x+nsolver/2-i,j) = (*By_domain)(n_p[0]-i,j);
                (*By_)(domain_oversize_x+nsolver/2-i,j) = (*By_domain)(n_p[0]-i,j);
                (*Ex_)(domain_oversize_x+nsolver/2-i,j) = (*Ex_domain)(n_p[0]-i,j);
                (*Dx_)(domain_oversize_x+nsolver/2-i,j) = (*Ex_domain)(n_p[0]-i,j);
                (*Ez_)(domain_oversize_x+nsolver/2-i,j) = (*Ez_domain)(n_p[0]-i,j);
                (*Dz_)(domain_oversize_x+nsolver/2-i,j) = (*Ez_domain)(n_p[0]-i,j);
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);

        //Injecting a laser
        vector<double> yp( 1 );
        for( unsigned int j=patch->isYmin() ; j<n_p[1]-patch->isYmax() ; j++ ) {

            double byE = 0.;
            yp[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*d[1];

            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                byE += vecLaser[ilaser]->getAmplitude0( yp, time_dual, j, 0 );
            }
            (*Hy_)( domain_oversize_x+nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_E*byE ;
            (*By_)( domain_oversize_x+nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_E*byE ;
        }

        vector<double> yd( 1 );
        for( unsigned int j=patch->isYmin() ; j<n_d[1]-patch->isYmax() ; j++ ) {

            double bzE = 0.;
            yd[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*d[1];

            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bzE += vecLaser[ilaser]->getAmplitude1( yd, time_dual, j, 0 );
            }
            (*Hz_)( domain_oversize_x+nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_E*bzE ;
            (*Bz_)( domain_oversize_x+nsolver/2, j ) += factor_laser_space_time*factor_laser_angle_E*bzE ;
        }

        // 4. Exchange Domain -> PML
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2-1 ; i++){
            for( unsigned int j = 0 ; j<n_p[1] ; j++ ) {
                (*Ez_domain)(n_p[0]-1-i,j) = (*Ez_)(domain_oversize_x+nsolver/2-1-i,j);
            }
            for( unsigned int j = 0 ; j<n_d[1] ; j++ ) {
                (*Ey_domain)(n_p[0]-1-i,j) = (*Ey_)(domain_oversize_x+nsolver/2-1-i,j);
                (*Bx_domain)(n_p[0]-1-i,j) = (*Hx_)(domain_oversize_x+nsolver/2-1-i,j);
            }
        }
        // Duals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for( unsigned int j = 0 ; j<n_p[1] ; j++ ) {
                (*Ex_domain)(n_d[0]-1-i,j) = (*Ex_)(domain_oversize_x+nsolver/2-i,j);
                (*By_domain)(n_d[0]-1-i,j) = (*Hy_)(domain_oversize_x+nsolver/2-i,j);
            }
            for( unsigned int j = 0 ; j<n_d[1] ; j++ ) {
                (*Bz_domain)(n_d[0]-1-i,j) = (*Hz_)(domain_oversize_x+nsolver/2-i,j);
            }
        }
    }
    else if( i_boundary_ == 2 ) {

        ElectroMagnBC2D_PML* pml_fields_xmin = NULL;
        ElectroMagnBC2D_PML* pml_fields_xmax = NULL;

        if(ncells_pml_xmin != 0){
            pml_fields_xmin = static_cast<ElectroMagnBC2D_PML*>( EMfields->emBoundCond[0] );
        }
        if(ncells_pml_xmax != 0){
            pml_fields_xmax = static_cast<ElectroMagnBC2D_PML*>( EMfields->emBoundCond[1] );
        }

        Field2D* Ex_pml_xmin = NULL;
        Field2D* Ey_pml_xmin = NULL;
        Field2D* Ez_pml_xmin = NULL;
        Field2D* Hx_pml_xmin = NULL;
        Field2D* Hy_pml_xmin = NULL;
        Field2D* Hz_pml_xmin = NULL;
        Field2D* Dx_pml_xmin = NULL;
        Field2D* Dy_pml_xmin = NULL;
        Field2D* Dz_pml_xmin = NULL;
        Field2D* Bx_pml_xmin = NULL;
        Field2D* By_pml_xmin = NULL;
        Field2D* Bz_pml_xmin = NULL;

        Field2D* Ex_pml_xmax = NULL;
        Field2D* Ey_pml_xmax = NULL;
        Field2D* Ez_pml_xmax = NULL;
        Field2D* Hx_pml_xmax = NULL;
        Field2D* Hy_pml_xmax = NULL;
        Field2D* Hz_pml_xmax = NULL;
        Field2D* Dx_pml_xmax = NULL;
        Field2D* Dy_pml_xmax = NULL;
        Field2D* Dz_pml_xmax = NULL;
        Field2D* Bx_pml_xmax = NULL;
        Field2D* By_pml_xmax = NULL;
        Field2D* Bz_pml_xmax = NULL;

        if(ncells_pml_xmin != 0){
            Ex_pml_xmin = pml_fields_xmin->Ex_;
            Ey_pml_xmin = pml_fields_xmin->Ey_;
            Ez_pml_xmin = pml_fields_xmin->Ez_;
            Hx_pml_xmin = pml_fields_xmin->Hx_;
            Hy_pml_xmin = pml_fields_xmin->Hy_;
            Hz_pml_xmin = pml_fields_xmin->Hz_;
            Dx_pml_xmin = pml_fields_xmin->Dx_;
            Dy_pml_xmin = pml_fields_xmin->Dy_;
            Dz_pml_xmin = pml_fields_xmin->Dz_;
            Bx_pml_xmin = pml_fields_xmin->Bx_;
            By_pml_xmin = pml_fields_xmin->By_;
            Bz_pml_xmin = pml_fields_xmin->Bz_;
        }

        if(ncells_pml_xmax != 0){
            Ex_pml_xmax = pml_fields_xmax->Ex_;
            Ey_pml_xmax = pml_fields_xmax->Ey_;
            Ez_pml_xmax = pml_fields_xmax->Ez_;
            Hx_pml_xmax = pml_fields_xmax->Hx_;
            Hy_pml_xmax = pml_fields_xmax->Hy_;
            Hz_pml_xmax = pml_fields_xmax->Hz_;
            Dx_pml_xmax = pml_fields_xmax->Dx_;
            Dy_pml_xmax = pml_fields_xmax->Dy_;
            Dz_pml_xmax = pml_fields_xmax->Dz_;
            Bx_pml_xmax = pml_fields_xmax->Bx_;
            By_pml_xmax = pml_fields_xmax->By_;
            Bz_pml_xmax = pml_fields_xmax->Bz_;
        }

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        int idx_start = 0;
                        // Les qtes Primals
                        (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bx_pml_xmin)(i,j);
                        (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Hx_pml_xmin)(i,j);
                        (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ey_pml_xmin)(i,j);
                        (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Dy_pml_xmin)(i,j);
                        (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ez_pml_xmin)(i,j);
                        (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Dz_pml_xmin)(i,j);
                        // Les qtes Duals
                        (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ex_pml_xmin)(i,j);
                        (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Dx_pml_xmin)(i,j);
                        (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Hy_pml_xmin)(i,j);
                        (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*By_pml_xmin)(i,j);
                        (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Hz_pml_xmin)(i,j);
                        (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bz_pml_xmin)(i,j);
                    }
                }
            }
            for( unsigned int i = 0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ex_domain)(i,j);
                (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ex_domain)(i,j);
                (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*By_domain)(i,j);
                (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*By_domain)(i,j);
                (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bz_domain)(i,j);
                (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bz_domain)(i,j);
            }
            for( unsigned int i = 0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bx_domain)(i,j);
                (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bx_domain)(i,j);
                (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ey_domain)(i,j);
                (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ey_domain)(i,j);
                (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ez_domain)(i,j);
                (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ez_domain)(i,j);
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        // Les qtes Primals commencent a (ypml_size_in_x+1)-ncells_pml_xmax
                        (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bx_pml_xmax)(domain_oversize_x+nsolver/2+i,j);
                        (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Hx_pml_xmax)(domain_oversize_x+nsolver/2+i,j);
                        (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ey_pml_xmax)(domain_oversize_x+nsolver/2+i,j);
                        (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Dy_pml_xmax)(domain_oversize_x+nsolver/2+i,j);
                        (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ez_pml_xmax)(domain_oversize_x+nsolver/2+i,j);
                        (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Dz_pml_xmax)(domain_oversize_x+nsolver/2+i,j);
                        // Toutes les qtes Duals commence a +1
                        (*Ex_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Ex_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j);
                        (*Dx_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Dx_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j);
                        (*Hy_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Hy_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j);
                        (*By_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*By_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j);
                        (*Hz_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Hz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j);
                        (*Bz_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j) = (*Bz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j);
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);

        //Injecting a laser
        vector<double> xp( 1 );
        for( unsigned int i=patch->isXmin() ; i<n_p[0]-patch->isXmax() ; i++ ) {

            double bxS = 0.;
            xp[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*d[0];

            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bxS += vecLaser[ilaser]->getAmplitude0( xp, time_dual, i, 0 );
            }
            (*Hx_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2 ) += factor_laser_space_time*factor_laser_angle_S*bxS ;
            (*Bx_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2 ) += factor_laser_space_time*factor_laser_angle_S*bxS ;
        }

        vector<double> xd( 1 );
        for( unsigned int i=patch->isXmin() ; i<n_d[0]-patch->isXmax() ; i++ ) {
            
            double bzS = 0.;
            xd[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - 0.5 - ( int )EMfields->oversize[0] )*d[0];
            
            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bzS += vecLaser[ilaser]->getAmplitude1( xd, time_dual, i, 0 );
            }
            (*Hz_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2 ) += factor_laser_space_time*factor_laser_angle_S*bzS ;
            (*Bz_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2 ) += factor_laser_space_time*factor_laser_angle_S*bzS ;
        }

        // 4. Exchange PML -> Domain
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for( unsigned int i = 0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ez_domain)(i,j) = (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
            }
            for( unsigned int i = 0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ex_domain)(i,j) = (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                (*By_domain)(i,j) = (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
            }
        }
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for( unsigned int i = 0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ey_domain)(i,j) = (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                (*Bx_domain)(i,j) = (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
            }
            for( unsigned int i = 0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Bz_domain)(i,j) = (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primal in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        // Primals
                        (*Ez_pml_xmin)(i,j) = (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Dz_pml_xmin)(i,j) = (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        // Duals
                        (*Ex_pml_xmin)(i,j) = (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Dx_pml_xmin)(i,j) = (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Hy_pml_xmin)(i,j) = (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*By_pml_xmin)(i,j) = (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                    }
                }
            }
        }
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        // Primals
                        (*Ey_pml_xmin)(i,j) = (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Dy_pml_xmin)(i,j) = (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Hx_pml_xmin)(i,j) = (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Bx_pml_xmin)(i,j) = (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        // Duals
                        (*Hz_pml_xmin)(i,j) = (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Bz_pml_xmin)(i,j) = (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                    }
                }
            }
        }

        // 6. Exchange PML y -> PML x MAX
        // Primal in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_domain_xmax ;
                        // Primals
                        (*Ez_pml_xmax)(i,j) = (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Dz_pml_xmax)(i,j) = (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        // Dual
                        (*Ex_pml_xmax)(i,j) = (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Dx_pml_xmax)(i,j) = (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Hy_pml_xmax)(i,j) = (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*By_pml_xmax)(i,j) = (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                    }
                }
            }
        }
        // Dual in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_domain_xmax ;
                        // Primals
                        (*Ey_pml_xmax)(i,j) = (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Dy_pml_xmax)(i,j) = (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Hx_pml_xmax)(i,j) = (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Bx_pml_xmax)(i,j) = (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        // Dual
                        (*Hz_pml_xmax)(i,j) = (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                        (*Bz_pml_xmax)(i,j) = (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j);
                    }
                }
            }
        }
    }
    else if( i_boundary_ == 3 ) {

        ElectroMagnBC2D_PML* pml_fields_xmin = NULL ;
        ElectroMagnBC2D_PML* pml_fields_xmax = NULL ;

        if(ncells_pml_xmin != 0){
            pml_fields_xmin = static_cast<ElectroMagnBC2D_PML*>( EMfields->emBoundCond[0] );
        }
        if(ncells_pml_xmax != 0){        
            pml_fields_xmax = static_cast<ElectroMagnBC2D_PML*>( EMfields->emBoundCond[1] );
        }

        Field2D* Ex_pml_xmin = NULL;
        Field2D* Ey_pml_xmin = NULL;
        Field2D* Ez_pml_xmin = NULL;
        Field2D* Hx_pml_xmin = NULL;
        Field2D* Hy_pml_xmin = NULL;
        Field2D* Hz_pml_xmin = NULL;
        Field2D* Dx_pml_xmin = NULL;
        Field2D* Dy_pml_xmin = NULL;
        Field2D* Dz_pml_xmin = NULL;
        Field2D* Bx_pml_xmin = NULL;
        Field2D* By_pml_xmin = NULL;
        Field2D* Bz_pml_xmin = NULL;

        Field2D* Ex_pml_xmax = NULL;
        Field2D* Ey_pml_xmax = NULL;
        Field2D* Ez_pml_xmax = NULL;
        Field2D* Hx_pml_xmax = NULL;
        Field2D* Hy_pml_xmax = NULL;
        Field2D* Hz_pml_xmax = NULL;
        Field2D* Dx_pml_xmax = NULL;
        Field2D* Dy_pml_xmax = NULL;
        Field2D* Dz_pml_xmax = NULL;
        Field2D* Bx_pml_xmax = NULL;
        Field2D* By_pml_xmax = NULL;
        Field2D* Bz_pml_xmax = NULL;

        if(ncells_pml_xmin != 0){
            Ex_pml_xmin = pml_fields_xmin->Ex_;
            Ey_pml_xmin = pml_fields_xmin->Ey_;
            Ez_pml_xmin = pml_fields_xmin->Ez_;
            Hx_pml_xmin = pml_fields_xmin->Hx_;
            Hy_pml_xmin = pml_fields_xmin->Hy_;
            Hz_pml_xmin = pml_fields_xmin->Hz_;
            Dx_pml_xmin = pml_fields_xmin->Dx_;
            Dy_pml_xmin = pml_fields_xmin->Dy_;
            Dz_pml_xmin = pml_fields_xmin->Dz_;
            Bx_pml_xmin = pml_fields_xmin->Bx_;
            By_pml_xmin = pml_fields_xmin->By_;
            Bz_pml_xmin = pml_fields_xmin->Bz_;
        }

        if(ncells_pml_xmax != 0){
            Ex_pml_xmax = pml_fields_xmax->Ex_;
            Ey_pml_xmax = pml_fields_xmax->Ey_;
            Ez_pml_xmax = pml_fields_xmax->Ez_;
            Hx_pml_xmax = pml_fields_xmax->Hx_;
            Hy_pml_xmax = pml_fields_xmax->Hy_;
            Hz_pml_xmax = pml_fields_xmax->Hz_;
            Dx_pml_xmax = pml_fields_xmax->Dx_;
            Dy_pml_xmax = pml_fields_xmax->Dy_;
            Dz_pml_xmax = pml_fields_xmax->Dz_;
            Bx_pml_xmax = pml_fields_xmax->Bx_;
            By_pml_xmax = pml_fields_xmax->By_;
            Bz_pml_xmax = pml_fields_xmax->Bz_;
        }

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        int idx_start = 0;
                        // Les qtes Primals
                        (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Bx_pml_xmin)(i,n_p[1]-j);
                        (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Hx_pml_xmin)(i,n_p[1]-j);
                        (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ey_pml_xmin)(i,n_p[1]-j);
                        (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Dy_pml_xmin)(i,n_p[1]-j);
                        (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ez_pml_xmin)(i,n_p[1]-j);
                        (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Dz_pml_xmin)(i,n_p[1]-j);
                        // Les qtes Duals
                        (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ex_pml_xmin)(i,n_p[1]-j);
                        (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Dx_pml_xmin)(i,n_p[1]-j);
                        (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Hy_pml_xmin)(i,n_p[1]-j);
                        (*By_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*By_pml_xmin)(i,n_p[1]-j);
                        (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Hz_pml_xmin)(i,n_p[1]-j);
                        (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Bz_pml_xmin)(i,n_p[1]-j);
                    }
                }
            }
            for( unsigned int i = 0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ex_domain)(i,n_p[1]-j);
                (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ex_domain)(i,n_p[1]-j);
                (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*By_domain)(i,n_p[1]-j);
                (*By_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*By_domain)(i,n_p[1]-j);
                (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Bz_domain)(i,n_p[1]-j);
                (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Bz_domain)(i,n_p[1]-j);
            }
            for( unsigned int i = 0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Bx_domain)(i,n_p[1]-j);
                (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Bx_domain)(i,n_p[1]-j);
                (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ey_domain)(i,n_p[1]-j);
                (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ey_domain)(i,n_p[1]-j);
                (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ez_domain)(i,n_p[1]-j);
                (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ez_domain)(i,n_p[1]-j);
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        // Les qtes Primals commencent a (ypml_size_in_x+1)-ncells_pml_xmax
                        (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Bx_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j);
                        (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Hx_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j);
                        (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ey_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j);
                        (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Dy_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j);
                        (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Ez_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j);
                        (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-j) = (*Dz_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j);
                        // Toutes les qtes Duals commence a +1
                        (*Ex_)(idx_start+1+i,domain_oversize_y+nsolver/2-j) = (*Ex_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j);
                        (*Dx_)(idx_start+1+i,domain_oversize_y+nsolver/2-j) = (*Dx_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j);
                        (*Hy_)(idx_start+1+i,domain_oversize_y+nsolver/2-j) = (*Hy_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j);
                        (*By_)(idx_start+1+i,domain_oversize_y+nsolver/2-j) = (*By_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j);
                        (*Hz_)(idx_start+1+i,domain_oversize_y+nsolver/2-j) = (*Hz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j);
                        (*Bz_)(idx_start+1+i,domain_oversize_y+nsolver/2-j) = (*Bz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j);
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, dimPrim, solvermin, solvermax);

        //Injecting a laser
        vector<double> xp( 1 );
        for( unsigned int i=patch->isXmin() ; i<n_p[0]-patch->isXmax() ; i++ ) {

            double bxN = 0.;
            xp[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*d[0];

            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bxN += vecLaser[ilaser]->getAmplitude0( xp, time_dual, i, 0 );
            }
            (*Hx_)( ncells_pml_xmin+i , domain_oversize_y+nsolver/2 ) += factor_laser_space_time*factor_laser_angle_N*bxN ;
            (*Bx_)( ncells_pml_xmin+i , domain_oversize_y+nsolver/2 ) += factor_laser_space_time*factor_laser_angle_N*bxN ;
        }

        vector<double> xd( 1 );
        for( unsigned int i=patch->isXmin() ; i<n_d[0]-patch->isXmax() ; i++ ) {

            double bzN = 0.;
            xd[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - 0.5 - ( int )EMfields->oversize[0] )*d[0];

            // Lasers
            for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                bzN += vecLaser[ilaser]->getAmplitude1( xd, time_dual, i, 0 );
            }
            (*Hz_)( ncells_pml_xmin+i, domain_oversize_y+nsolver/2 ) += factor_laser_space_time*factor_laser_angle_N*bzN ;
            (*Bz_)( ncells_pml_xmin+i, domain_oversize_y+nsolver/2 ) += factor_laser_space_time*factor_laser_angle_N*bzN ;
        }

        // 4. Exchange PML -> Domain
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2-1 ; j++){
            for( unsigned int i = 0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ez_domain)(i,n_p[1]-1-j) = (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
            }
            for( unsigned int i = 0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ex_domain)(i,n_p[1]-1-j) = (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                (*By_domain)(i,n_p[1]-1-j) = (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
            }
        }
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for( unsigned int i = 0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Ey_domain)(i,n_d[1]-1-j) = (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                (*Bx_domain)(i,n_d[1]-1-j) = (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j);
            }
            for( unsigned int i = 0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                (*Bz_domain)(i,n_d[1]-1-j) = (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j);
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2-1 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        // Primals
                        (*Ez_pml_xmin)(i,n_p[1]-1-j) = (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*Dz_pml_xmin)(i,n_p[1]-1-j) = (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        // Duals
                        (*Ex_pml_xmin)(i,n_p[1]-1-j) = (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*Dx_pml_xmin)(i,n_p[1]-1-j) = (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*Hy_pml_xmin)(i,n_p[1]-1-j) = (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*By_pml_xmin)(i,n_p[1]-1-j) = (*By_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                    }
                }
            }
        }
        // Dual in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        // Primals
                        (*Ey_pml_xmin)(i,n_d[1]-1-j) = (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Dy_pml_xmin)(i,n_d[1]-1-j) = (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Hx_pml_xmin)(i,n_d[1]-1-j) = (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Bx_pml_xmin)(i,n_d[1]-1-j) = (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        // Duals
                        (*Hz_pml_xmin)(i,n_d[1]-1-j) = (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Bz_pml_xmin)(i,n_d[1]-1-j) = (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                    }
                }
            }
        }

        // 6. Exchange PML y -> PML x MAX
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2-1 ; j++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_domain_xmax ;
                        // Primals
                        (*Ez_pml_xmax)(i,n_p[1]-1-j) = (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*Dz_pml_xmax)(i,n_p[1]-1-j) = (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        // Duals
                        (*Ex_pml_xmax)(i,n_p[1]-1-j) = (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*Dx_pml_xmax)(i,n_p[1]-1-j) = (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*Hy_pml_xmax)(i,n_p[1]-1-j) = (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                        (*By_pml_xmax)(i,n_p[1]-1-j) = (*By_)(idx_start+i,domain_oversize_y+nsolver/2-1-j);
                    }
                }
            }
        }
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_domain_xmax ;
                        // Les qtes Primals
                        (*Ey_pml_xmax)(i,n_d[1]-1-j) = (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Dy_pml_xmax)(i,n_d[1]-1-j) = (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Hx_pml_xmax)(i,n_d[1]-1-j) = (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Bx_pml_xmax)(i,n_d[1]-1-j) = (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        // Duals
                        (*Hz_pml_xmax)(i,n_d[1]-1-j) = (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                        (*Bz_pml_xmax)(i,n_d[1]-1-j) = (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j);
                    }
                }
            }
        }
    }
}

