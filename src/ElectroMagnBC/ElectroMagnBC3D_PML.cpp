#include "ElectroMagnBC3D_PML.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field3D.h"
#include "Tools.h"
#include "Laser.h"

#include "SolverFactory.h"

using namespace std;

ElectroMagnBC3D_PML::ElectroMagnBC3D_PML( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC3D( params, patch, i_boundary )
{

    std::vector<unsigned int> n_space(params.n_space);
    std::vector<unsigned int> oversize(params.oversize);
    if( params.multiple_decomposition ) {
        n_space = params.n_space_region;
        oversize = params.region_oversize;
    }

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

    if ( ( i_boundary_ == 0 && patch->isXmin() )
         || ( i_boundary_ == 1 && patch->isXmax() )
         || ( i_boundary_ == 2 && patch->isYmin() )
         || ( i_boundary_ == 3 && patch->isYmax() )
         || ( i_boundary_ == 4 && patch->isZmin() )
         || ( i_boundary_ == 5 && patch->isZmax() ) ) {

        int iDim = 0*((i_boundary_==0)||(i_boundary_==1))+1*((i_boundary_==2)||(i_boundary_==3))+2*((i_boundary_==4)||(i_boundary_==5));
        int min_or_max = (i_boundary_)%2;

        domain_oversize_x =  oversize[0] ;
        domain_oversize_y =  oversize[1] ;
        domain_oversize_z =  oversize[2] ;

        if (patch->isXmin() ) {//&& i_boundary_ == 0 ) {
            ncells_pml_xmin = params.number_of_pml_cells[0][0];
            ncells_pml_domain_xmin = ncells_pml_xmin + 1*oversize[0] + nsolver/2;
            domain_oversize_x = oversize[0] ;
        }
        else {
            ncells_pml_xmin = 0;
            ncells_pml_domain_xmin = 0;
        }
        if (patch->isXmax() ) {//&& i_boundary_ == 1 ) {
            ncells_pml_xmax = params.number_of_pml_cells[0][1];
            ncells_pml_domain_xmax = ncells_pml_xmax + 1*oversize[0] + nsolver/2;
            domain_oversize_x = oversize[0] ;
        }
        else {
            ncells_pml_xmax = 0;
            ncells_pml_domain_xmax = 0;
        }
        if (patch->isYmin() ) {//&& i_boundary_ == 2 ) {
            ncells_pml_ymin = params.number_of_pml_cells[1][0];
            ncells_pml_domain_ymin = ncells_pml_ymin + 1*oversize[1] + nsolver/2;
            domain_oversize_y = oversize[1] ;
        }
        else {
            ncells_pml_ymin = 0;
            ncells_pml_domain_ymin = 0;
        }
        if (patch->isYmax() ) {//&& i_boundary_ == 3 ) {
            ncells_pml_ymax = params.number_of_pml_cells[1][1];
            ncells_pml_domain_ymax = ncells_pml_ymax + 1*oversize[1] + nsolver/2;
            domain_oversize_y = oversize[1] ;
        }
        else {
            ncells_pml_ymax = 0;
            ncells_pml_domain_ymax = 0;
        }
        if (patch->isZmin() ) {//&& i_boundary_ == 2 ) {
            ncells_pml_zmin = params.number_of_pml_cells[2][0];
            ncells_pml_domain_zmin = ncells_pml_zmin + 1*oversize[2] + nsolver/2;
            domain_oversize_z = oversize[2] ;
        }
        else {
            ncells_pml_zmin = 0;
            ncells_pml_domain_zmin = 0;
        }
        if (patch->isZmax() ) {//&& i_boundary_ == 3 ) {
            ncells_pml_zmax = params.number_of_pml_cells[2][1];
            ncells_pml_domain_zmax = ncells_pml_zmax + 1*oversize[2] + nsolver/2;
            domain_oversize_z = oversize[2] ;
        }
        else {
            ncells_pml_zmax = 0;
            ncells_pml_domain_zmax = 0;
        }

        ncells_pml = params.number_of_pml_cells[iDim][min_or_max];
        ncells_pml_domain = ncells_pml+1*oversize[iDim] + nsolver/2;

        // Define min and max idx to exchange
        // the good data f(solver,oversize)
        if (min_or_max==0){
            // if min border : Exchange of data (for domain to pml-domain)
            // min2exchange <= i < max2exchange
            min2exchange = 1*nsolver/2 ;
            max2exchange = 2*nsolver/2 ;
            // Solver
            solvermin = nsolver/2 ;
            solvermax = ncells_pml_domain - oversize[iDim] ;
        }
        else if (min_or_max==1){
            // if max border : Exchange of data (for domain to pml-domain)
            // min2exchange <= i < max2exchange
            min2exchange = 1*nsolver/2 ;
            max2exchange = 2*nsolver/2 ;
            // Solver
            solvermin = oversize[iDim] + nsolver/2 - nsolver/2 + 1 ;
            solvermax = ncells_pml_domain-nsolver/2 ;
        }

        if (ncells_pml==0){
            ERROR("PML domain have to be >0 cells in thickness");
        }

        std::vector<unsigned int> dimPrim( params.nDim_field );
        for( int i=0 ; i<params.nDim_field ; i++ ) {
            dimPrim[i] = n_space[i]+1+2*oversize[i];
        }
        dimPrim[iDim] = ncells_pml_domain;
        // Redefine the size of the PMLx in x, PMLy in y and PMLz in z (thickness)
        // -----------------------------------------------------------
        // ncells_pml_domain = ncells_pml+1*oversize[iDim] + nsolver/2;
        // -----------------------------------------------------------
        //  n_space    ->  ncells_pml
        // +1          -> +nsolver/2
        // +2*oversize -> +1*oversize
        if ( iDim==1 ){
            // If the PML domain in Y is in Xmin or Xmax too, add cell orthogonally
            dimPrim[iDim-1] += ncells_pml_xmin + ncells_pml_xmax ;
            ypml_size_in_x = dimPrim[iDim-1] ;
        }
        if ( iDim==2 ){
            dimPrim[0] += ncells_pml_xmin + ncells_pml_xmax ;
            dimPrim[1] += ncells_pml_ymin + ncells_pml_ymax ;
            zpml_size_in_x = dimPrim[0] ;
            zpml_size_in_y = dimPrim[1] ;
        }

        startpml = oversize[iDim]+nsolver/2;

        int ncells_pml_min[2];
        ncells_pml_min[0] = ncells_pml_xmin;
        ncells_pml_min[1] = ncells_pml_ymin;
        int ncells_pml_max[2];
        ncells_pml_max[0] = ncells_pml_xmax;
        ncells_pml_max[1] = ncells_pml_ymax;

        pml_solver_->setDomainSizeAndCoefficients( iDim, min_or_max, ncells_pml_domain, startpml, ncells_pml_min, ncells_pml_max, patch );

        Ex_ = new Field3D( dimPrim, 0, false, "Ex_pml" );
        Ey_ = new Field3D( dimPrim, 1, false, "Ey_pml" );
        Ez_ = new Field3D( dimPrim, 2, false, "Ez_pml" );
        Bx_ = new Field3D( dimPrim, 0, true, "Bx_pml" );
        By_ = new Field3D( dimPrim, 1, true, "By_pml" );
        Bz_ = new Field3D( dimPrim, 2, true, "Bz_pml" );
        Dx_ = new Field3D( dimPrim, 0, false, "Dx_pml" );
        Dy_ = new Field3D( dimPrim, 1, false, "Dy_pml" );
        Dz_ = new Field3D( dimPrim, 2, false, "Dz_pml" );
        Hx_ = new Field3D( dimPrim, 0, true, "Hx_pml" );
        Hy_ = new Field3D( dimPrim, 1, true, "Hy_pml" );
        Hz_ = new Field3D( dimPrim, 2, true, "Hz_pml" );

        //Laser parameter
        double pyKx, pyKy, pyKz;
        double kx, ky, kz;
        double Knorm;
        double omega = 1. ;

        // Xmin boundary
        pyKx = params.EM_BCs_k[0][0];
        pyKy = params.EM_BCs_k[0][1];
        pyKz = params.EM_BCs_k[0][2];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;
        kz = omega*pyKz/Knorm;

        factor_laser_space_time = 2.*dt_ov_d[0] ;
        factor_laser_angle_W = factor_laser_space_time*kx/Knorm;

        // Xmax boundary
        pyKx = params.EM_BCs_k[1][0];
        pyKy = params.EM_BCs_k[1][1];
        pyKz = params.EM_BCs_k[1][2];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;
        kz = omega*pyKz/Knorm;

        factor_laser_space_time = 2.*dt_ov_d[0] ;
        factor_laser_angle_E = factor_laser_space_time*kx/Knorm;

        // Ymin boundary
        pyKx = params.EM_BCs_k[2][0];
        pyKy = params.EM_BCs_k[2][1];
        pyKz = params.EM_BCs_k[2][2];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;
        kz = omega*pyKz/Knorm;

        factor_laser_space_time = 2.*dt_ov_d[1] ;
        factor_laser_angle_S = factor_laser_space_time*ky/Knorm;

        // Ymax boundary
        pyKx = params.EM_BCs_k[3][0];
        pyKy = params.EM_BCs_k[3][1];
        pyKz = params.EM_BCs_k[3][2];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;
        kz = omega*pyKz/Knorm;

        factor_laser_space_time = 2.*dt_ov_d[1] ;
        factor_laser_angle_N = factor_laser_space_time*ky/Knorm;

        // Zmin boundary
        pyKx = params.EM_BCs_k[4][0];
        pyKy = params.EM_BCs_k[4][1];
        pyKz = params.EM_BCs_k[4][2];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;
        kz = omega*pyKz/Knorm;

        factor_laser_space_time = 2.*dt_ov_d[2] ;
        factor_laser_angle_B = factor_laser_space_time*kz/Knorm;

        // Zmax boundary
        pyKx = params.EM_BCs_k[5][0];
        pyKy = params.EM_BCs_k[5][1];
        pyKz = params.EM_BCs_k[5][2];
        Knorm = sqrt( pyKx*pyKx + pyKy*pyKy + pyKz*pyKz ) ;
        kx = omega*pyKx/Knorm;
        ky = omega*pyKy/Knorm;
        kz = omega*pyKz/Knorm;

        factor_laser_space_time = 2.*dt_ov_d[2] ;
        factor_laser_angle_T = factor_laser_space_time*kz/Knorm;


    }
}


ElectroMagnBC3D_PML::~ElectroMagnBC3D_PML()
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


void ElectroMagnBC3D_PML::save_fields( Field *my_field, Patch *patch )
{
}


void ElectroMagnBC3D_PML::disableExternalFields()
{
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC3D_PML::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{

    int iDim = 0*((i_boundary_==0)||(i_boundary_==1))+1*((i_boundary_==2)||(i_boundary_==3))+2*((i_boundary_==4)||(i_boundary_==5));
    int min_or_max = (i_boundary_)%2;

    Field3D *Ex_domain = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey_domain = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez_domain = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx_domain = static_cast<Field3D *>( EMfields->Bx_ );
    Field3D *By_domain = static_cast<Field3D *>( EMfields->By_ );
    Field3D *Bz_domain = static_cast<Field3D *>( EMfields->Bz_ );

    vector<double> pos( 2 );

    if( i_boundary_ == 0 && patch->isXmin() ) {

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        // 2. Exchange field PML <- Domain
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            // MESSAGE("Copy PML < Domain");
            // MESSAGE(ncells_pml_domain-domain_oversize_x-nsolver/2+i<<"<"<<i);
            for ( int j=0 ; j<n_d[1] ; j++ ) {
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Hx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Bx_domain)(i,j,k);
                    (*Bx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Bx_domain)(i,j,k);
                }
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Ey_domain)(i,j,k);
                    (*Dy_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Ey_domain)(i,j,k);
                    (*Hz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Bz_domain)(i,j,k);
                    (*Bz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Bz_domain)(i,j,k);
                }
            }
            for ( int j=0 ; j<n_p[1] ; j++ ) {
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Ez_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Ez_domain)(i,j,k);
                    (*Dz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Ez_domain)(i,j,k);
                    (*Hy_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*By_domain)(i,j,k);
                    (*By_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*By_domain)(i,j,k);
                }
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Ex_domain)(i,j,k);
                    (*Dx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*Ex_domain)(i,j,k);
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        //Injecting a laser
        vector<double> by( n_p[1]*n_d[2], 0. );
        for( unsigned int j=patch->isYmin() ; j<n_p[1]-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*d[1];
            for( unsigned int k=patch->isZmin() ; k<n_d[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - 0.5 - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    by[ j*n_d[2]+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                }
                (*Hy_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j , k ) += factor_laser_angle_W*by[j*n_d[2]+k] ;
                (*By_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j , k ) += factor_laser_angle_W*by[j*n_d[2]+k] ;
            }
        }

        vector<double> bz( n_d[1]*n_p[2], 0. );
        for( unsigned int j=patch->isYmin() ; j<n_d[1]-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*d[1];
            for( unsigned int k=patch->isZmin() ; k<n_p[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ j*n_p[2]+k] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, k );
                }
                (*Hz_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j , k) += factor_laser_angle_W*bz[ j*n_p[2]+k ] ;
                (*Bz_)( ncells_pml_domain-domain_oversize_x-nsolver/2, j , k) += factor_laser_angle_W*bz[ j*n_p[2]+k ] ;
            }
        }

        // 4. Exchange PML -> Domain
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( int j=0 ; j<n_p[1] ; j++ ) {
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Ez_domain)(i,j,k) = (*Ez_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                }
            }
            for ( int j=0 ; j<n_d[1] ; j++ ) {
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_domain)(i,j,k) = (*Ey_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Bx_domain)(i,j,k) = (*Hx_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                }
            }
        }
        // Duals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( int j=0 ; j<n_p[1] ; j++ ) {
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_domain)(i,j,k) = (*Ex_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*By_domain)(i,j,k) = (*Hy_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                }
            }
            for ( int j=0 ; j<n_d[1] ; j++ ) {
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Bz_domain)(i,j,k) = (*Hz_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                }
            }
        }
    }

    else if( i_boundary_ == 1 && patch->isXmax() ) {

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        // 2. Exchange field Domain -> PML
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            for ( int j=0 ; j<n_d[1] ; j++ ) {
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Hx_)(domain_oversize_x+nsolver/2-i,j,k) = (*Bx_domain)(n_p[0]-i,j,k);
                    (*Bx_)(domain_oversize_x+nsolver/2-i,j,k) = (*Bx_domain)(n_p[0]-i,j,k);
                }
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_)(domain_oversize_x+nsolver/2-i,j,k) = (*Ey_domain)(n_p[0]-i,j,k);
                    (*Dy_)(domain_oversize_x+nsolver/2-i,j,k) = (*Ey_domain)(n_p[0]-i,j,k);
                    (*Hz_)(domain_oversize_x+nsolver/2-i,j,k) = (*Bz_domain)(n_p[0]-i,j,k);
                    (*Bz_)(domain_oversize_x+nsolver/2-i,j,k) = (*Bz_domain)(n_p[0]-i,j,k);
                }
            }
            for ( int j=0 ; j<n_p[1] ; j++ ) {
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Ez_)(domain_oversize_x+nsolver/2-i,j,k) = (*Ez_domain)(n_p[0]-i,j,k);
                    (*Dz_)(domain_oversize_x+nsolver/2-i,j,k) = (*Ez_domain)(n_p[0]-i,j,k);
                    (*Hy_)(domain_oversize_x+nsolver/2-i,j,k) = (*By_domain)(n_p[0]-i,j,k);
                    (*By_)(domain_oversize_x+nsolver/2-i,j,k) = (*By_domain)(n_p[0]-i,j,k);
                }
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_)(domain_oversize_x+nsolver/2-i,j,k) = (*Ex_domain)(n_p[0]-i,j,k);
                    (*Dx_)(domain_oversize_x+nsolver/2-i,j,k) = (*Ex_domain)(n_p[0]-i,j,k);
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        //Injecting a laser
        vector<double> by( n_p[1]*n_d[2], 0. );
        for( unsigned int j=patch->isYmin() ; j<n_p[1]-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*d[1];
            for( unsigned int k=patch->isZmin() ; k<n_d[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - 0.5 - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    by[ j*n_d[2]+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, k );
                }
                (*Hy_)( domain_oversize_x+nsolver/2, j, k ) += factor_laser_angle_E*by[j*n_d[2]+k] ;
                (*By_)( domain_oversize_x+nsolver/2, j, k ) += factor_laser_angle_E*by[j*n_d[2]+k] ;
            }
        }

        vector<double> bz( n_d[1]*n_p[2], 0. );
        for( unsigned int j=patch->isYmin() ; j<n_d[1]-patch->isYmax() ; j++ ) {
            pos[0] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*d[1];
            for( unsigned int k=patch->isZmin() ; k<n_p[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ j*n_p[2]+k] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, k );
                }
                (*Hz_)( domain_oversize_x+nsolver/2, j, k ) += factor_laser_angle_E*bz[ j*n_p[2]+k ] ;
                (*Bz_)( domain_oversize_x+nsolver/2, j, k ) += factor_laser_angle_E*bz[ j*n_p[2]+k ] ;
            }
        } 

        // 4. Exchange Domain -> PML
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2-1 ; i++){
            for ( int j=0 ; j<n_p[1] ; j++ ) {
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Ez_domain)(n_p[0]-1-i,j,k) = (*Ez_)(domain_oversize_x+nsolver/2-1-i,j,k);
                }
            }
            for ( int j=0 ; j<n_d[1] ; j++ ) {
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_domain)(n_p[0]-1-i,j,k) = (*Ey_)(domain_oversize_x+nsolver/2-1-i,j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Bx_domain)(n_p[0]-1-i,j,k) = (*Hx_)(domain_oversize_x+nsolver/2-1-i,j,k);
                }
            }
        }
        // Duals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( int j=0 ; j<n_p[1] ; j++ ) {
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_domain)(n_d[0]-1-i,j,k) = (*Ex_)(domain_oversize_x+nsolver/2-i,j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*By_domain)(n_d[0]-1-i,j,k) = (*Hy_)(domain_oversize_x+nsolver/2-i,j,k);
                }
            }
            for ( int j=0 ; j<n_d[1] ; j++ ) {
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Bz_domain)(n_d[0]-1-i,j,k) = (*Hz_)(domain_oversize_x+nsolver/2-i,j,k);
                }
            }
        }
    }
    else if( i_boundary_ == 2 && patch->isYmin() ) {

        ElectroMagnBC3D_PML* pml_fields_xmin = NULL ;
        ElectroMagnBC3D_PML* pml_fields_xmax = NULL ;

        if(ncells_pml_xmin != 0){
            pml_fields_xmin = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[0] );
        }
        if(ncells_pml_xmax != 0){
            pml_fields_xmax = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[1] );
        }

        Field3D* Ex_pml_xmin = NULL;
        Field3D* Ey_pml_xmin = NULL;
        Field3D* Ez_pml_xmin = NULL;
        Field3D* Hx_pml_xmin = NULL;
        Field3D* Hy_pml_xmin = NULL;
        Field3D* Hz_pml_xmin = NULL;
        Field3D* Dx_pml_xmin = NULL;
        Field3D* Dy_pml_xmin = NULL;
        Field3D* Dz_pml_xmin = NULL;
        Field3D* Bx_pml_xmin = NULL;
        Field3D* By_pml_xmin = NULL;

        Field3D* Bz_pml_xmin = NULL;
        Field3D* Ex_pml_xmax = NULL;
        Field3D* Ey_pml_xmax = NULL;
        Field3D* Ez_pml_xmax = NULL;
        Field3D* Hx_pml_xmax = NULL;
        Field3D* Hy_pml_xmax = NULL;
        Field3D* Hz_pml_xmax = NULL;
        Field3D* Dx_pml_xmax = NULL;
        Field3D* Dy_pml_xmax = NULL;
        Field3D* Dz_pml_xmax = NULL;
        Field3D* Bx_pml_xmax = NULL;
        Field3D* By_pml_xmax = NULL;
        Field3D* Bz_pml_xmax = NULL;

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
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        int idx_start = 0;
                        // Les qtes i-Primals commencent a 0
                        // Toutes les qtes i-Duals 0
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ey_pml_xmin)(i,j,k);
                            (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Dy_pml_xmin)(i,j,k);
                            // i-Duals
                            (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ex_pml_xmin)(i,j,k);
                            (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Dx_pml_xmin)(i,j,k);
                            (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Hz_pml_xmin)(i,j,k);
                            (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bz_pml_xmin)(i,j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ez_pml_xmin)(i,j,k);
                            (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Dz_pml_xmin)(i,j,k);
                            (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bx_pml_xmin)(i,j,k);
                            (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Hx_pml_xmin)(i,j,k);
                            // i-Duals
                            (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Hy_pml_xmin)(i,j,k);
                            (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*By_pml_xmin)(i,j,k);
                        }
                    }
                }
            }
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ey_domain)(i,j,k);
                    (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ey_domain)(i,j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bx_domain)(i,j,k);
                    (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bx_domain)(i,j,k);
                    (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ez_domain)(i,j,k);
                    (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ez_domain)(i,j,k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ex_domain)(i,j,k);
                    (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ex_domain)(i,j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*By_domain)(i,j,k);
                    (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*By_domain)(i,j,k);
                    (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bz_domain)(i,j,k);
                    (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bz_domain)(i,j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        // Les qtes i-Primals commencent a (ypml_size_in_x+1)-ncells_pml_xmax
                        // Toutes les qtes i-Duals commence a idx_start+1
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ey_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Dy_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            // i-Duals
                            (*Ex_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ex_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            (*Dx_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Dx_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            (*Hz_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Hz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            (*Bz_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Ez_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Dz_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Bx_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Hx_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            // i-Duals
                            (*Hy_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*Hy_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            (*By_)(idx_start+1+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*By_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        //Injecting a laser
        vector<double> bx( n_p[0]*n_d[2], 0. );
        for( unsigned int i=patch->isXmin() ; i<n_p[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int k=patch->isZmin() ; k<n_d[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k -0.5 - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bx[ i*n_d[2]+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, i, k );
                }
                (*Hx_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2, k ) += factor_laser_angle_S*bx[ i*n_d[2]+k ] ;
                (*Bx_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2, k ) += factor_laser_angle_S*bx[ i*n_d[2]+k ] ;
            }
        }

        vector<double> bz( n_d[0]*n_p[2], 0. );
        for( unsigned int i=patch->isXmin() ; i<n_d[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i -0.5 - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int k=patch->isZmin() ; k<n_p[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ i*n_p[2]+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, i, k );
                }
                (*Hz_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2, k ) += factor_laser_angle_S*bz[ i*n_p[2]+k ] ;
                (*Bz_)( ncells_pml_xmin+i, ncells_pml_domain-domain_oversize_y-nsolver/2, k ) += factor_laser_angle_S*bz[ i*n_p[2]+k ] ;
            }
        }

        // 4. Exchange PML -> Domain
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Ez_domain)(i,j,k) = (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_domain)(i,j,k) = (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*By_domain)(i,j,k) = (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                }
            }
        }
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_domain)(i,j,k) = (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Bx_domain)(i,j,k) = (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Bz_domain)(i,j,k) = (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                }
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primal in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Duals
                            (*Ex_pml_xmin)(i,j,k) = (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Dx_pml_xmin)(i,j,k) = (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Ez_pml_xmin)(i,j,k) = (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Dz_pml_xmin)(i,j,k) = (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            // i-Duals
                            (*Hy_pml_xmin)(i,j,k) = (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*By_pml_xmin)(i,j,k) = (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
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
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            (*Ey_pml_xmin)(i,j,k) = (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Dy_pml_xmin)(i,j,k) = (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            // i-Duals
                            (*Hz_pml_xmin)(i,j,k) = (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Bz_pml_xmin)(i,j,k) = (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Hx_pml_xmin)(i,j,k) = (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Bx_pml_xmin)(i,j,k) = (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
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
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            // Nothing
                            // i-Dual
                            (*Ex_pml_xmax)(i,j,k) = (*Ex_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Dx_pml_xmax)(i,j,k) = (*Dx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Ez_pml_xmax)(i,j,k) = (*Ez_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Dz_pml_xmax)(i,j,k) = (*Dz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            // i-Dual
                            (*Hy_pml_xmax)(i,j,k) = (*Hy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*By_pml_xmax)(i,j,k) = (*By_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
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
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            (*Ey_pml_xmax)(i,j,k) = (*Ey_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Dy_pml_xmax)(i,j,k) = (*Dy_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            // i-Dual
                            (*Hz_pml_xmax)(i,j,k) = (*Hz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Bz_pml_xmax)(i,j,k) = (*Bz_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Hx_pml_xmax)(i,j,k) = (*Hx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                            (*Bx_pml_xmax)(i,j,k) = (*Bx_)(idx_start+i,ncells_pml_domain-domain_oversize_x-nsolver/2+j,k);
                        }
                    }
                }
            }
        }
    }
    else if( i_boundary_ == 3 && patch->isYmax() ) {

        ElectroMagnBC3D_PML* pml_fields_xmin = NULL ;
        ElectroMagnBC3D_PML* pml_fields_xmax = NULL ;

        if(ncells_pml_xmin != 0){
            pml_fields_xmin = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[0] );
        }
        if(ncells_pml_xmax != 0){
            pml_fields_xmax = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[1] );
        }

        Field3D* Ex_pml_xmin = NULL;
        Field3D* Ey_pml_xmin = NULL;
        Field3D* Ez_pml_xmin = NULL;
        Field3D* Hx_pml_xmin = NULL;
        Field3D* Hy_pml_xmin = NULL;
        Field3D* Hz_pml_xmin = NULL;
        Field3D* Dx_pml_xmin = NULL;
        Field3D* Dy_pml_xmin = NULL;
        Field3D* Dz_pml_xmin = NULL;
        Field3D* Bx_pml_xmin = NULL;
        Field3D* By_pml_xmin = NULL;
        Field3D* Bz_pml_xmin = NULL;

        Field3D* Ex_pml_xmax = NULL;
        Field3D* Ey_pml_xmax = NULL;
        Field3D* Ez_pml_xmax = NULL;
        Field3D* Hx_pml_xmax = NULL;
        Field3D* Hy_pml_xmax = NULL;
        Field3D* Hz_pml_xmax = NULL;
        Field3D* Dx_pml_xmax = NULL;
        Field3D* Dy_pml_xmax = NULL;
        Field3D* Dz_pml_xmax = NULL;
        Field3D* Bx_pml_xmax = NULL;
        Field3D* By_pml_xmax = NULL;
        Field3D* Bz_pml_xmax = NULL;

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
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        int idx_start = 0;
                        // Les qtes i-Primals commencent a 0
                        // Toutes les qtes i-Duals 0
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // Les qtes i-Primals
                            (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ey_pml_xmin)(i,n_p[1]-j,k);
                            (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Dy_pml_xmin)(i,n_p[1]-j,k);
                            // Les qtes i-Duals
                            (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ex_pml_xmin)(i,n_p[1]-j,k);
                            (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Dx_pml_xmin)(i,n_p[1]-j,k);
                            (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Hz_pml_xmin)(i,n_p[1]-j,k);
                            (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Bz_pml_xmin)(i,n_p[1]-j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // Les qtes i-Primals
                            (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ez_pml_xmin)(i,n_p[1]-j,k);
                            (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Dz_pml_xmin)(i,n_p[1]-j,k);
                            (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Bx_pml_xmin)(i,n_p[1]-j,k);
                            (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Hx_pml_xmin)(i,n_p[1]-j,k);
                            // Les qtes i-Duals
                            (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Hy_pml_xmin)(i,n_p[1]-j,k);
                            (*By_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*By_pml_xmin)(i,n_p[1]-j,k);
                        }
                    }
                }
            }
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ey_domain)(i,n_p[1]-j,k);
                    (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ey_domain)(i,n_p[1]-j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Bx_domain)(i,n_p[1]-j,k);
                    (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Bx_domain)(i,n_p[1]-j,k);
                    (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ez_domain)(i,n_p[1]-j,k);
                    (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ez_domain)(i,n_p[1]-j,k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ex_domain)(i,n_p[1]-j,k);
                    (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ex_domain)(i,n_p[1]-j,k);
                    (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Bz_domain)(i,n_p[1]-j,k);
                    (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Bz_domain)(i,n_p[1]-j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*By_domain)(i,n_p[1]-j,k);
                    (*By_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*By_domain)(i,n_p[1]-j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        // Les qtes i-Primals commencent a (ypml_size_in_x+1)-ncells_pml_xmax
                        // Toutes les qtes i-Duals commence a idx_start + 1
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ey_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j,k);
                            (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Dy_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j,k);
                            // i-Duals
                            (*Ex_)(idx_start+1+i,domain_oversize_y+nsolver/2-j,k) = (*Ex_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j,k);
                            (*Dx_)(idx_start+1+i,domain_oversize_y+nsolver/2-j,k) = (*Dx_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j,k);
                            (*Hz_)(idx_start+1+i,domain_oversize_y+nsolver/2-j,k) = (*Hz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j,k);
                            (*Bz_)(idx_start+1+i,domain_oversize_y+nsolver/2-j,k) = (*Bz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Bx_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j,k);
                            (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Hx_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j,k);
                            (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Ez_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j,k);
                            (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*Dz_pml_xmax)(domain_oversize_x+nsolver/2+i,n_p[1]-j,k);
                            // i-Duals
                            (*Hy_)(idx_start+1+i,domain_oversize_y+nsolver/2-j,k) = (*Hy_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j,k);
                            (*By_)(idx_start+1+i,domain_oversize_y+nsolver/2-j,k) = (*By_pml_xmax)(domain_oversize_x+nsolver/2+1+i,n_p[1]-j,k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        //Injecting a laser
        vector<double> bx( n_p[0]*n_d[2], 0. );
        for( unsigned int i=patch->isXmin() ; i<n_p[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int k=patch->isZmin() ; k<n_d[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k -0.5 - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bx[ i*n_d[2]+k ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, i, k );
                }
                (*Hx_)( ncells_pml_xmin+i , domain_oversize_y+nsolver/2, k ) += factor_laser_angle_N*bx[ i*n_d[2]+k ] ;
                (*Bx_)( ncells_pml_xmin+i , domain_oversize_y+nsolver/2, k ) += factor_laser_angle_N*bx[ i*n_d[2]+k ] ;
            }
        }

        vector<double> bz( n_d[0]*n_p[2], 0. );
        for( unsigned int i=patch->isXmin() ; i<n_d[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i -0.5 - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int k=patch->isZmin() ; k<n_p[2]-patch->isZmax() ; k++ ) {
                pos[1] = patch->getDomainLocalMin( 2 ) + ( ( int )k - ( int )EMfields->oversize[2] )*d[2];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bz[ i*n_p[2]+k ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, i, k );
                }
                (*Hz_)( ncells_pml_xmin+i, domain_oversize_y+nsolver/2, k ) += factor_laser_angle_N*bz[ i*n_p[2]+k ] ;
                (*Bz_)( ncells_pml_xmin+i, domain_oversize_y+nsolver/2, k ) += factor_laser_angle_N*bz[ i*n_p[2]+k ] ;
            }
        }

        // 4. Exchange PML -> Domain
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2-1 ; j++){
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Ez_domain)(i,n_p[1]-1-j,k) = (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ex_domain)(i,n_p[1]-1-j,k) = (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*By_domain)(i,n_p[1]-1-j,k) = (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                }
            }
        }
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Ey_domain)(i,n_d[1]-1-j,k) = (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                }
                for ( int k=0 ; k<n_d[2] ; k++ ) {
                    (*Bx_domain)(i,n_d[1]-1-j,k) = (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( int k=0 ; k<n_p[2] ; k++ ) {
                    (*Bz_domain)(i,n_d[1]-1-j,k) = (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                }
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2-1 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            // Nothing
                            // i-Duals
                            (*Ex_pml_xmin)(i,n_p[1]-1-j,k) = (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*Dx_pml_xmin)(i,n_p[1]-1-j,k) = (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Ez_pml_xmin)(i,n_p[1]-1-j,k) = (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*Dz_pml_xmin)(i,n_p[1]-1-j,k) = (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            // i-Duals
                            (*Hy_pml_xmin)(i,n_p[1]-1-j,k) = (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*By_pml_xmin)(i,n_p[1]-1-j,k) = (*By_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                        }
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
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            (*Ey_pml_xmin)(i,n_d[1]-1-j,k) = (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*Dy_pml_xmin)(i,n_d[1]-1-j,k) = (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            // i-Duals
                            (*Hz_pml_xmin)(i,n_d[1]-1-j,k) = (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*Bz_pml_xmin)(i,n_d[1]-1-j,k) = (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Hx_pml_xmin)(i,n_d[1]-1-j,k) = (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*Bx_pml_xmin)(i,n_d[1]-1-j,k) = (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            // i-Duals
                            // Nothing
                        }
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
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            // Nothing
                            // i-Duals
                            (*Ex_pml_xmax)(i,n_p[1]-1-j,k) = (*Ex_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*Dx_pml_xmax)(i,n_p[1]-1-j,k) = (*Dx_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Ez_pml_xmax)(i,n_p[1]-1-j,k) = (*Ez_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*Dz_pml_xmax)(i,n_p[1]-1-j,k) = (*Dz_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            // i-Duals
                            (*Hy_pml_xmax)(i,n_p[1]-1-j,k) = (*Hy_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*By_pml_xmax)(i,n_p[1]-1-j,k) = (*By_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                        }
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
                        for ( int k=0 ; k<n_p[2] ; k++ ) {
                            // i-Primals
                            (*Ey_pml_xmax)(i,n_d[1]-1-j,k) = (*Ey_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*Dy_pml_xmax)(i,n_d[1]-1-j,k) = (*Dy_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            // i-Duals
                            (*Hz_pml_xmax)(i,n_d[1]-1-j,k) = (*Hz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*Bz_pml_xmax)(i,n_d[1]-1-j,k) = (*Bz_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                        }
                        for ( int k=0 ; k<n_d[2] ; k++ ) {
                            // i-Primals
                            (*Hx_pml_xmax)(i,n_d[1]-1-j,k) = (*Hx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*Bx_pml_xmax)(i,n_d[1]-1-j,k) = (*Bx_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            // i-Duals
                            // Nothing
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 4 && patch->isZmin() ) {

        ElectroMagnBC3D_PML* pml_fields_xmin = NULL ;
        ElectroMagnBC3D_PML* pml_fields_xmax = NULL ;
        ElectroMagnBC3D_PML* pml_fields_ymin = NULL ;
        ElectroMagnBC3D_PML* pml_fields_ymax = NULL ;

        if(ncells_pml_xmin != 0){
            pml_fields_xmin = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[0] );
        }
        if(ncells_pml_xmax != 0){
            pml_fields_xmax = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[1] );
        }
        if(ncells_pml_ymin != 0){
            pml_fields_ymin = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[2] );
        }
        if(ncells_pml_ymax != 0){
            pml_fields_ymax = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[3] );
        }

        Field3D* Ex_pml_xmin = NULL;
        Field3D* Ey_pml_xmin = NULL;
        Field3D* Ez_pml_xmin = NULL;
        Field3D* Hx_pml_xmin = NULL;
        Field3D* Hy_pml_xmin = NULL;
        Field3D* Hz_pml_xmin = NULL;
        Field3D* Dx_pml_xmin = NULL;
        Field3D* Dy_pml_xmin = NULL;
        Field3D* Dz_pml_xmin = NULL;
        Field3D* Bx_pml_xmin = NULL;
        Field3D* By_pml_xmin = NULL;
        Field3D* Bz_pml_xmin = NULL;
        Field3D* Ex_pml_xmax = NULL;
        Field3D* Ey_pml_xmax = NULL;
        Field3D* Ez_pml_xmax = NULL;
        Field3D* Hx_pml_xmax = NULL;
        Field3D* Hy_pml_xmax = NULL;
        Field3D* Hz_pml_xmax = NULL;
        Field3D* Dx_pml_xmax = NULL;
        Field3D* Dy_pml_xmax = NULL;
        Field3D* Dz_pml_xmax = NULL;
        Field3D* Bx_pml_xmax = NULL;
        Field3D* By_pml_xmax = NULL;
        Field3D* Bz_pml_xmax = NULL;

        Field3D* Ex_pml_ymin = NULL;
        Field3D* Ey_pml_ymin = NULL;
        Field3D* Ez_pml_ymin = NULL;
        Field3D* Hx_pml_ymin = NULL;
        Field3D* Hy_pml_ymin = NULL;
        Field3D* Hz_pml_ymin = NULL;
        Field3D* Dx_pml_ymin = NULL;
        Field3D* Dy_pml_ymin = NULL;
        Field3D* Dz_pml_ymin = NULL;
        Field3D* Bx_pml_ymin = NULL;
        Field3D* By_pml_ymin = NULL;
        Field3D* Bz_pml_ymin = NULL;
        Field3D* Ex_pml_ymax = NULL;
        Field3D* Ey_pml_ymax = NULL;
        Field3D* Ez_pml_ymax = NULL;
        Field3D* Hx_pml_ymax = NULL;
        Field3D* Hy_pml_ymax = NULL;
        Field3D* Hz_pml_ymax = NULL;
        Field3D* Dx_pml_ymax = NULL;
        Field3D* Dy_pml_ymax = NULL;
        Field3D* Dz_pml_ymax = NULL;
        Field3D* Bx_pml_ymax = NULL;
        Field3D* By_pml_ymax = NULL;
        Field3D* Bz_pml_ymax = NULL;

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

        if(ncells_pml_ymin != 0){
            Ex_pml_ymin = pml_fields_ymin->Ex_;
            Ey_pml_ymin = pml_fields_ymin->Ey_;
            Ez_pml_ymin = pml_fields_ymin->Ez_;
            Hx_pml_ymin = pml_fields_ymin->Hx_;
            Hy_pml_ymin = pml_fields_ymin->Hy_;
            Hz_pml_ymin = pml_fields_ymin->Hz_;
            Dx_pml_ymin = pml_fields_ymin->Dx_;
            Dy_pml_ymin = pml_fields_ymin->Dy_;
            Dz_pml_ymin = pml_fields_ymin->Dz_;
            Bx_pml_ymin = pml_fields_ymin->Bx_;
            By_pml_ymin = pml_fields_ymin->By_;
            Bz_pml_ymin = pml_fields_ymin->Bz_;
        }

        if(ncells_pml_ymax != 0){
            Ex_pml_ymax = pml_fields_ymax->Ex_;
            Ey_pml_ymax = pml_fields_ymax->Ey_;
            Ez_pml_ymax = pml_fields_ymax->Ez_;
            Hx_pml_ymax = pml_fields_ymax->Hx_;
            Hy_pml_ymax = pml_fields_ymax->Hy_;
            Hz_pml_ymax = pml_fields_ymax->Hz_;
            Dx_pml_ymax = pml_fields_ymax->Dx_;
            Dy_pml_ymax = pml_fields_ymax->Dy_;
            Dz_pml_ymax = pml_fields_ymax->Dz_;
            Bx_pml_ymax = pml_fields_ymax->Bx_;
            By_pml_ymax = pml_fields_ymax->By_;
            Bz_pml_ymax = pml_fields_ymax->Bz_;
        }

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        /*
        Dans le bloc du dessous (2.)
        Peut être faire passer les echanges XMAX et YMAX apres l'echange du domaine
        */

        // 2. Exchange field PML <- Domain
        for ( int k=min2exchange ; k<max2exchange ; k++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    // Les qtes i-Primals commencent a 0
                    // Toutes les qtes i-Duals 0
                    for ( int j=0 ; j<n_p[1] ; j++ ) {
                        for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                            //Ex
                            (*Ex_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ex_pml_xmin)(i,j,k);
                            (*Dx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dx_pml_xmin)(i,j,k);
                            //Ez
                            (*Ez_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ez_pml_xmin)(i,j,k);
                            (*Dz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dz_pml_xmin)(i,j,k);
                            //By
                            (*By_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*By_pml_xmin)(i,j,k);
                            (*Hy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hy_pml_xmin)(i,j,k);
                        }
                    }
                    for ( int j=0 ; j<n_d[1] ; j++ ) {
                        for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                            //Ey
                            (*Ey_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ey_pml_xmin)(i,j,k);
                            (*Dy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dy_pml_xmin)(i,j,k);
                            //Bx
                            (*Bx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bx_pml_xmin)(i,j,k);
                            (*Hx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hx_pml_xmin)(i,j,k);
                            //Bz
                            (*Bz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bz_pml_xmin)(i,j,k);
                            (*Hz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hz_pml_xmin)(i,j,k);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    // Les qtes j-Primals commencent a 0
                    // Toutes les qtes j-Duals 0
                    for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            // Ey
                            (*Ey_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ey_pml_ymin)(i,j,k);
                            (*Dy_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dy_pml_ymin)(i,j,k);
                            // Ez
                            (*Ez_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ez_pml_ymin)(i,j,k);
                            (*Dz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dz_pml_ymin)(i,j,k);
                            // Bx
                            (*Bx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bx_pml_ymin)(i,j,k);
                            (*Hx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hx_pml_ymin)(i,j,k);
                        }
                    }
                    for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            // Ex
                            (*Ex_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ex_pml_ymin)(i,j,k);
                            (*Dx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dx_pml_ymin)(i,j,k);
                            // By
                            (*By_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*By_pml_ymin)(i,j,k);
                            (*Hy_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hy_pml_ymin)(i,j,k);
                            // Bz
                            (*Bz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bz_pml_ymin)(i,j,k);
                            (*Hz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hz_pml_ymin)(i,j,k);
                        }
                    }
                }
            }
            // Au milieu du domain
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            // Si on est sur les bord min du domaine, au fait un petit shift
            // Les donnees ont deja ete echange au dessus
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    // Ez
                    (*Ez_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ez_domain)(i,j,k);
                    (*Dz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ez_domain)(i,j,k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    // Ey
                    (*Ey_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ey_domain)(i,j,k);
                    (*Dy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ey_domain)(i,j,k);
                    // Bx
                    (*Bx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bx_domain)(i,j,k);
                    (*Hx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bx_domain)(i,j,k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    // Ex
                    (*Ex_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ex_domain)(i,j,k);
                    (*Dx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ex_domain)(i,j,k);
                    // By
                    (*By_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*By_domain)(i,j,k);
                    (*Hy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*By_domain)(i,j,k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    // Bz
                    (*Bz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bz_domain)(i,j,k);
                    (*Hz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bz_domain)(i,j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = (zpml_size_in_x-1)-(ncells_pml_xmax-1) ;
                    int jdx_start = ncells_pml_ymin;
                    // Les qtes i-Primals commencent a (zpml_size_in_x+1)-ncells_pml_xmax
                    // Toutes les qtes i-Duals commence a idx_start+1
                    for ( int j=0 ; j<n_p[1] ; j++ ) {
                        for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                            //Ex
                            (*Ex_)(idx_start+1+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ex_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            (*Dx_)(idx_start+1+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dx_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            //Ez
                            (*Ez_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ez_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            (*Dz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dz_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            //By
                            (*By_)(idx_start+1+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*By_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            (*Hy_)(idx_start+1+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hy_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                        }
                    }
                    for ( int j=0 ; j<n_d[1] ; j++ ) {
                        for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                            //Ey
                            (*Ey_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ey_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            (*Dy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dy_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            //Bx
                            (*Bx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bx_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            (*Hx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hx_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                            //Bz
                            (*Bz_)(idx_start+1+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                            (*Hz_)(idx_start+1+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,k);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = (zpml_size_in_y-1)-(ncells_pml_ymax-1) ;
                    // Les qtes j-Primals commencent a (zpml_size_in_y+1)-ncells_pml_ymax
                    // Toutes les qtes j-Duals commence a jdx_start+1
                    for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            // Ey
                            (*Ey_)(i,jdx_start+j+1,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ey_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,k);
                            (*Dy_)(i,jdx_start+j+1,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dy_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,k);
                            // Ez
                            (*Ez_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ez_pml_ymax)(i,domain_oversize_y+nsolver/2+j,k);
                            (*Dz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dz_pml_ymax)(i,domain_oversize_y+nsolver/2+j,k);
                            // Bx
                            (*Bx_)(i,jdx_start+j+1,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bx_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,k);
                            (*Hx_)(i,jdx_start+j+1,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hx_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,k);
                        }
                    }
                    for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            // Ex
                            (*Ex_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Ex_pml_ymax)(i,domain_oversize_y+nsolver/2+j,k);
                            (*Dx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Dx_pml_ymax)(i,domain_oversize_y+nsolver/2+j,k);
                            // By
                            (*By_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*By_pml_ymax)(i,domain_oversize_y+nsolver/2+j,k);
                            (*Hy_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hy_pml_ymax)(i,domain_oversize_y+nsolver/2+j,k);
                            // Bz
                            (*Bz_)(i,jdx_start+1+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Bz_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,k);
                            (*Hz_)(i,jdx_start+1+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*Hz_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,k);
                        }
                    }
                }
            }
        }

        // /*
        // Dans le bloc du dessus (2.)
        // Peut être faire passer les echanges XMAX et YMAX apres l'echange du domaine
        // */

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        //Injecting a laser
        vector<double> bx( n_p[0]*n_d[1], 0. ); // Bx(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<n_p[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int j=patch->isYmin() ; j<n_d[1]-patch->isYmax() ; j++ ) {
                pos[1] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*d[1];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bx[ i*n_d[1]+j ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, i, j );
                }
                (*Hx_)( ncells_pml_xmin+i, ncells_pml_ymin+j, ncells_pml_domain-domain_oversize_z-nsolver/2 ) += factor_laser_angle_B*bx[ i*n_d[1]+j ] ;
                (*Bx_)( ncells_pml_xmin+i, ncells_pml_ymin+j, ncells_pml_domain-domain_oversize_z-nsolver/2 ) += factor_laser_angle_B*bx[ i*n_d[1]+j ] ;
            }
        }

        vector<double> by( n_d[0]*n_p[1], 0. ); // By(d,p,d)
        for( unsigned int i=patch->isXmin() ; i<n_d[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i -0.5 - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int j=patch->isYmin() ; j<n_p[1]-patch->isYmax() ; j++ ) {
                pos[1] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*d[1];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    by[ i*n_p[1]+j ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, i, j );
                }
                (*Hy_)( ncells_pml_xmin+i, ncells_pml_ymin+j, ncells_pml_domain-domain_oversize_z-nsolver/2 ) += factor_laser_angle_B*by[ i*n_p[1]+j ] ;
                (*By_)( ncells_pml_xmin+i, ncells_pml_ymin+j, ncells_pml_domain-domain_oversize_z-nsolver/2 ) += factor_laser_angle_B*by[ i*n_p[1]+j ] ;
            }
        }

        // 4. Exchange PML -> Domain
        // Ici il faut a priori remplacer y<->z et j<->k
        // Primals in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                // for ( int j=0 ; j<n_p[1] ; j++ ) {
                //     No field
                // }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Ey_domain)(i,j,k) = (*Ey_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_x-nsolver/2+k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Ex_domain)(i,j,k) = (*Ex_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_x-nsolver/2+k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Bz_domain)(i,j,k) = (*Hz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_x-nsolver/2+k);
                }
            }
        }
        // Duals in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Ez_domain)(i,j,k) = (*Ez_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_x-nsolver/2+k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Bx_domain)(i,j,k) = (*Hx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_x-nsolver/2+k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*By_domain)(i,j,k) = (*Hy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_x-nsolver/2+k);
                }
                // for ( int j=0 ; j<n_d[1] ; j++ ) {
                //     No field
                // }
            }
        }

        // 5. Exchange PML z -> PML xMIN
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Primals
                            // i-Duals
                            (*Ex_pml_xmin)(i,j,k) = (*Ex_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dx_pml_xmin)(i,j,k) = (*Dx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Ey_pml_xmin)(i,j,k) = (*Ey_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dy_pml_xmin)(i,j,k) = (*Dy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // i-Duals
                            (*Hz_pml_xmin)(i,j,k) = (*Hz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bz_pml_xmin)(i,j,k) = (*Bz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }
        // Duals in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Primals
                            (*Ez_pml_xmin)(i,j,k) = (*Ez_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dz_pml_xmin)(i,j,k) = (*Dz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // i-Duals
                            (*Hy_pml_xmin)(i,j,k) = (*Hy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*By_pml_xmin)(i,j,k) = (*By_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Hx_pml_xmin)(i,j,k) = (*Hx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bx_pml_xmin)(i,j,k) = (*Bx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // i-Duals
                        }
                    }
                }
            }
        }

        // 7. Exchange PML z -> PML xMAX
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = zpml_size_in_x-ncells_pml_domain_xmax ;
                        int jdx_start = ncells_pml_ymin;
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Dual
                            (*Ex_pml_xmax)(i,j,k) = (*Ex_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dx_pml_xmax)(i,j,k) = (*Dx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Ey_pml_xmax)(i,j,k) = (*Ey_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dy_pml_xmax)(i,j,k) = (*Dy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // i-Dual
                            (*Hz_pml_xmax)(i,j,k) = (*Hz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bz_pml_xmax)(i,j,k) = (*Bz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }
        // Dual in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = zpml_size_in_x-ncells_pml_domain_xmax ;
                        int jdx_start = ncells_pml_ymin;
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Primals
                            (*Ez_pml_xmax)(i,j,k) = (*Ez_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dz_pml_xmax)(i,j,k) = (*Dz_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // i-Dual
                            (*Hy_pml_xmax)(i,j,k) = (*Hy_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*By_pml_xmax)(i,j,k) = (*By_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Hx_pml_xmax)(i,j,k) = (*Hx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bx_pml_xmax)(i,j,k) = (*Bx_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }

        // 6. Exchange PML z -> PML yMin
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int j=0 ; j<ncells_pml_domain_ymin ; j++ ) {
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            // j-Duals
                            (*Ey_pml_ymin)(i,j,k) = (*Ey_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dy_pml_ymin)(i,j,k) = (*Dy_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ex_pml_ymin)(i,j,k) = (*Ex_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dx_pml_ymin)(i,j,k) = (*Dx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // j-Duals
                            (*Hz_pml_ymin)(i,j,k) = (*Hz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bz_pml_ymin)(i,j,k) = (*Bz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }
        // Dual in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int j=0 ; j<ncells_pml_domain_ymin ; j++ ) {
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ez_pml_ymin)(i,j,k) = (*Ez_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dz_pml_ymin)(i,j,k) = (*Dz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // j-Duals
                            (*Hx_pml_ymin)(i,j,k) = (*Hx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bx_pml_ymin)(i,j,k) = (*Bx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Hy_pml_ymin)(i,j,k) = (*Hy_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*By_pml_ymin)(i,j,k) = (*By_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // j-Duals
                        }
                    }
                }
            }
        }

        // 8. Exchange PML z -> PML yMax
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    for ( int j=0 ; j<ncells_pml_domain_ymax ; j++ ) {
                        int jdx_start = zpml_size_in_y-ncells_pml_domain_ymax ;
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            // j-Dual
                            (*Ey_pml_ymax)(i,j,k) = (*Ey_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dy_pml_ymax)(i,j,k) = (*Dy_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ex_pml_ymax)(i,j,k) = (*Ex_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dx_pml_ymax)(i,j,k) = (*Dx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // j-Duals
                            (*Hz_pml_ymax)(i,j,k) = (*Hz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bz_pml_ymax)(i,j,k) = (*Bz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }
        // Dual in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    for ( int j=0 ; j<ncells_pml_domain_ymax ; j++ ) {
                        int jdx_start = zpml_size_in_y-ncells_pml_domain_ymax ;
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ez_pml_ymax)(i,j,k) = (*Ez_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Dz_pml_ymax)(i,j,k) = (*Dz_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // j-Duals
                            (*Hx_pml_ymax)(i,j,k) = (*Hx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*Bx_pml_ymax)(i,j,k) = (*Bx_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Hy_pml_ymax)(i,j,k) = (*Hy_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*By_pml_ymax)(i,j,k) = (*By_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            // j-Duals
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 5 && patch->isZmax() ) {

        ElectroMagnBC3D_PML* pml_fields_xmin = NULL ;
        ElectroMagnBC3D_PML* pml_fields_xmax = NULL ;
        ElectroMagnBC3D_PML* pml_fields_ymin = NULL ;
        ElectroMagnBC3D_PML* pml_fields_ymax = NULL ;

        if(ncells_pml_xmin != 0){
            pml_fields_xmin = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[0] );
        }
        if(ncells_pml_xmax != 0){
            pml_fields_xmax = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[1] );
        }
        if(ncells_pml_ymin != 0){
            pml_fields_ymin = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[2] );
        }
        if(ncells_pml_ymax != 0){
            pml_fields_ymax = static_cast<ElectroMagnBC3D_PML*>( EMfields->emBoundCond[3] );
        }

        Field3D* Ex_pml_xmin = NULL;
        Field3D* Ey_pml_xmin = NULL;
        Field3D* Ez_pml_xmin = NULL;
        Field3D* Hx_pml_xmin = NULL;
        Field3D* Hy_pml_xmin = NULL;
        Field3D* Hz_pml_xmin = NULL;
        Field3D* Dx_pml_xmin = NULL;
        Field3D* Dy_pml_xmin = NULL;
        Field3D* Dz_pml_xmin = NULL;
        Field3D* Bx_pml_xmin = NULL;
        Field3D* By_pml_xmin = NULL;
        Field3D* Bz_pml_xmin = NULL;
        Field3D* Ex_pml_xmax = NULL;
        Field3D* Ey_pml_xmax = NULL;
        Field3D* Ez_pml_xmax = NULL;
        Field3D* Hx_pml_xmax = NULL;
        Field3D* Hy_pml_xmax = NULL;
        Field3D* Hz_pml_xmax = NULL;
        Field3D* Dx_pml_xmax = NULL;
        Field3D* Dy_pml_xmax = NULL;
        Field3D* Dz_pml_xmax = NULL;
        Field3D* Bx_pml_xmax = NULL;
        Field3D* By_pml_xmax = NULL;
        Field3D* Bz_pml_xmax = NULL;

        Field3D* Ex_pml_ymin = NULL;
        Field3D* Ey_pml_ymin = NULL;
        Field3D* Ez_pml_ymin = NULL;
        Field3D* Hx_pml_ymin = NULL;
        Field3D* Hy_pml_ymin = NULL;
        Field3D* Hz_pml_ymin = NULL;
        Field3D* Dx_pml_ymin = NULL;
        Field3D* Dy_pml_ymin = NULL;
        Field3D* Dz_pml_ymin = NULL;
        Field3D* Bx_pml_ymin = NULL;
        Field3D* By_pml_ymin = NULL;
        Field3D* Bz_pml_ymin = NULL;
        Field3D* Ex_pml_ymax = NULL;
        Field3D* Ey_pml_ymax = NULL;
        Field3D* Ez_pml_ymax = NULL;
        Field3D* Hx_pml_ymax = NULL;
        Field3D* Hy_pml_ymax = NULL;
        Field3D* Hz_pml_ymax = NULL;
        Field3D* Dx_pml_ymax = NULL;
        Field3D* Dy_pml_ymax = NULL;
        Field3D* Dz_pml_ymax = NULL;
        Field3D* Bx_pml_ymax = NULL;
        Field3D* By_pml_ymax = NULL;
        Field3D* Bz_pml_ymax = NULL;

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

        if(ncells_pml_ymin != 0){
            Ex_pml_ymin = pml_fields_ymin->Ex_;
            Ey_pml_ymin = pml_fields_ymin->Ey_;
            Ez_pml_ymin = pml_fields_ymin->Ez_;
            Hx_pml_ymin = pml_fields_ymin->Hx_;
            Hy_pml_ymin = pml_fields_ymin->Hy_;
            Hz_pml_ymin = pml_fields_ymin->Hz_;
            Dx_pml_ymin = pml_fields_ymin->Dx_;
            Dy_pml_ymin = pml_fields_ymin->Dy_;
            Dz_pml_ymin = pml_fields_ymin->Dz_;
            Bx_pml_ymin = pml_fields_ymin->Bx_;
            By_pml_ymin = pml_fields_ymin->By_;
            Bz_pml_ymin = pml_fields_ymin->Bz_;
        }

        if(ncells_pml_ymax != 0){
            Ex_pml_ymax = pml_fields_ymax->Ex_;
            Ey_pml_ymax = pml_fields_ymax->Ey_;
            Ez_pml_ymax = pml_fields_ymax->Ez_;
            Hx_pml_ymax = pml_fields_ymax->Hx_;
            Hy_pml_ymax = pml_fields_ymax->Hy_;
            Hz_pml_ymax = pml_fields_ymax->Hz_;
            Dx_pml_ymax = pml_fields_ymax->Dx_;
            Dy_pml_ymax = pml_fields_ymax->Dy_;
            Dz_pml_ymax = pml_fields_ymax->Dz_;
            Bx_pml_ymax = pml_fields_ymax->Bx_;
            By_pml_ymax = pml_fields_ymax->By_;
            Bz_pml_ymax = pml_fields_ymax->Bz_;
        }

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        /*
        Dans le bloc du dessous (2.)
        Peut être faire passer les echanges XMAX et YMAX apres l'echange du domaine
        */

        // 2. Exchange field PML <- Domain
        for ( int k=min2exchange ; k<max2exchange ; k++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            //Ex
                            (*Ex_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ex_pml_xmin)(i,j,n_p[2]-k);
                            (*Dx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dx_pml_xmin)(i,j,n_p[2]-k);
                            //Ez
                            (*Ez_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ez_pml_xmin)(i,j,n_p[2]-k);
                            (*Dz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dz_pml_xmin)(i,j,n_p[2]-k);
                            //By
                            (*By_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*By_pml_xmin)(i,j,n_p[2]-k);
                            (*Hy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hy_pml_xmin)(i,j,n_p[2]-k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            //Ey
                            (*Ey_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ey_pml_xmin)(i,j,n_p[2]-k);
                            (*Dy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dy_pml_xmin)(i,j,n_p[2]-k);
                            //Bx
                            (*Bx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bx_pml_xmin)(i,j,n_p[2]-k);
                            (*Hx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hx_pml_xmin)(i,j,n_p[2]-k);
                            //Bz
                            (*Bz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bz_pml_xmin)(i,j,n_p[2]-k);
                            (*Hz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hz_pml_xmin)(i,j,n_p[2]-k);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            // Ey
                            (*Ey_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ey_pml_ymin)(i,j,n_p[2]-k);
                            (*Dy_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dy_pml_ymin)(i,j,n_p[2]-k);
                            // Ez
                            (*Ez_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ez_pml_ymin)(i,j,n_p[2]-k);
                            (*Dz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dz_pml_ymin)(i,j,n_p[2]-k);
                            // Bx
                            (*Bx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bx_pml_ymin)(i,j,n_p[2]-k);
                            (*Hx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hx_pml_ymin)(i,j,n_p[2]-k);
                        }
                    }
                    for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            // Ex
                            (*Ex_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ex_pml_ymin)(i,j,n_p[2]-k);
                            (*Dx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dx_pml_ymin)(i,j,n_p[2]-k);
                            // By
                            (*By_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*By_pml_ymin)(i,j,n_p[2]-k);
                            (*Hy_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hy_pml_ymin)(i,j,n_p[2]-k);
                            // Bz
                            (*Bz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bz_pml_ymin)(i,j,n_p[2]-k);
                            (*Hz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hz_pml_ymin)(i,j,n_p[2]-k);
                        }
                    }
                }
            }
            // Au milieu du domain
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    // Ez
                    (*Ez_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ez_domain)(i,j,n_p[2]-k);
                    (*Dz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ez_domain)(i,j,n_p[2]-k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    // Ey
                    (*Ey_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ey_domain)(i,j,n_p[2]-k);
                    (*Dy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ey_domain)(i,j,n_p[2]-k);
                    // Bx
                    (*Bx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bx_domain)(i,j,n_p[2]-k);
                    (*Hx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bx_domain)(i,j,n_p[2]-k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    // Ex
                    (*Ex_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ex_domain)(i,j,n_p[2]-k);
                    (*Dx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ex_domain)(i,j,n_p[2]-k);
                    // By
                    (*By_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*By_domain)(i,j,n_p[2]-k);
                    (*Hy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*By_domain)(i,j,n_p[2]-k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    // Bz
                    (*Bz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bz_domain)(i,j,n_p[2]-k);
                    (*Hz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bz_domain)(i,j,n_p[2]-k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = (zpml_size_in_x-1)-(ncells_pml_xmax-1) ;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            //Ex
                            (*Ex_)(idx_start+1+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ex_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,n_p[2]-k);
                            (*Dx_)(idx_start+1+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dx_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,n_p[2]-k);
                            //Ez
                            (*Ez_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ez_pml_xmax)(domain_oversize_x+nsolver/2+i,j,n_p[2]-k);
                            (*Dz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dz_pml_xmax)(domain_oversize_x+nsolver/2+i,j,n_p[2]-k);
                            //By
                            (*By_)(idx_start+1+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*By_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,n_p[2]-k);
                            (*Hy_)(idx_start+1+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hy_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,n_p[2]-k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            //Ey
                            (*Ey_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ey_pml_xmax)(domain_oversize_x+nsolver/2+i,j,n_p[2]-k);
                            (*Dy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dy_pml_xmax)(domain_oversize_x+nsolver/2+i,j,n_p[2]-k);
                            //Bx
                            (*Bx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bx_pml_xmax)(domain_oversize_x+nsolver/2+i,j,n_p[2]-k);
                            (*Hx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hx_pml_xmax)(domain_oversize_x+nsolver/2+i,j,n_p[2]-k);
                            //Bz
                            (*Bz_)(idx_start+1+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Bz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,n_p[2]-k);
                            (*Hz_)(idx_start+1+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hz_pml_xmax)(domain_oversize_x+nsolver/2+1+i,j,n_p[2]-k);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = (zpml_size_in_y-1)-(ncells_pml_ymax-1) ;
                    // Les qtes j-Primals commencent a (zpml_size_in_y+1)-ncells_pml_ymax
                    // Toutes les qtes j-Duals commence a jdx_start+1
                    for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            // Ey
                            (*Ey_)(i,jdx_start+j+1,domain_oversize_z+nsolver/2-k) = (*Ey_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,n_p[2]-k);
                            (*Dy_)(i,jdx_start+j+1,domain_oversize_z+nsolver/2-k) = (*Dy_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,n_p[2]-k);
                            // Ez
                            (*Ez_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ez_pml_ymax)(i,domain_oversize_y+nsolver/2+j,n_p[2]-k);
                            (*Dz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dz_pml_ymax)(i,domain_oversize_y+nsolver/2+j,n_p[2]-k);
                            // Bx
                            (*Bx_)(i,jdx_start+j+1,domain_oversize_z+nsolver/2-k) = (*Bx_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,n_p[2]-k);
                            (*Hx_)(i,jdx_start+j+1,domain_oversize_z+nsolver/2-k) = (*Hx_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,n_p[2]-k);
                        }
                    }
                    for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            // Ex
                            (*Ex_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Ex_pml_ymax)(i,domain_oversize_y+nsolver/2+j,n_p[2]-k);
                            (*Dx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Dx_pml_ymax)(i,domain_oversize_y+nsolver/2+j,n_p[2]-k);
                            // By
                            (*By_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*By_pml_ymax)(i,domain_oversize_y+nsolver/2+j,n_p[2]-k);
                            (*Hy_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*Hy_pml_ymax)(i,domain_oversize_y+nsolver/2+j,n_p[2]-k);
                            // Bz
                            (*Bz_)(i,jdx_start+1+j,domain_oversize_z+nsolver/2-k) = (*Bz_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,n_p[2]-k);
                            (*Hz_)(i,jdx_start+1+j,domain_oversize_z+nsolver/2-k) = (*Hz_pml_ymax)(i,domain_oversize_y+nsolver/2+1+j,n_p[2]-k);
                        }
                    }
                }
            }
        }

        /*
        Dans le bloc du dessus (2.)
        Peut être faire passer les echanges XMAX et YMAX apres l'echange du domaine
        */

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        //Injecting a laser
        vector<double> bx( n_p[0]*n_d[1], 0. ); // Bx(p,d,d)
        for( unsigned int i=patch->isXmin() ; i<n_p[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int j=patch->isYmin() ; j<n_d[1]-patch->isYmax() ; j++ ) {
                pos[1] = patch->getDomainLocalMin( 1 ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[1] )*d[1];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    bx[ i*n_d[1]+j ] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, i, j );
                }
                (*Hx_)( ncells_pml_xmin+i, ncells_pml_ymin+j, domain_oversize_z+nsolver/2 ) += factor_laser_angle_T*bx[ i*n_d[1]+j ] ;
                (*Bx_)( ncells_pml_xmin+i, ncells_pml_ymin+j, domain_oversize_z+nsolver/2 ) += factor_laser_angle_T*bx[ i*n_d[1]+j ] ;
            }
        }

        vector<double> by( n_d[0]*n_p[1], 0. ); // By(d,p,d)
        for( unsigned int i=patch->isXmin() ; i<n_d[0]-patch->isXmax() ; i++ ) {
            pos[0] = patch->getDomainLocalMin( 0 ) + ( ( int )i -0.5 - ( int )EMfields->oversize[0] )*d[0];
            for( unsigned int j=patch->isYmin() ; j<n_p[1]-patch->isYmax() ; j++ ) {
                pos[1] = patch->getDomainLocalMin( 1 ) + ( ( int )j - ( int )EMfields->oversize[1] )*d[1];
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    by[ i*n_p[1]+j ] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, i, j );
                }
                (*Hy_)( ncells_pml_xmin+i, ncells_pml_ymin+j, domain_oversize_z+nsolver/2 ) += factor_laser_angle_T*by[ i*n_p[1]+j ] ;
                (*By_)( ncells_pml_xmin+i, ncells_pml_ymin+j, domain_oversize_z+nsolver/2 ) += factor_laser_angle_T*by[ i*n_p[1]+j ] ;
            }
        }

        // 4. Exchange PML -> Domain
        // Ici il faut a priori remplacer y<->z et j<->k
        // Primals in z-direction
        for (int k=0 ; k < nsolver/2-1 ; k++){
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                // for ( int j=0 ; j<n_p[1] ; j++ ) {
                //     No field
                // }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Ey_domain)(i,j,n_p[2]-1-k) = (*Ey_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Ex_domain)(i,j,n_p[2]-1-k) = (*Ex_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Bz_domain)(i,j,n_p[2]-1-k) = (*Hz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                }
            }
        }
        // Duals in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( int i=0 ; i<n_p[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Ez_domain)(i,j,n_d[2]-1-k) = (*Ez_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Bx_domain)(i,j,n_d[2]-1-k) = (*Hx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                }
            }
            for ( int i=0 ; i<n_d[0] ; i++ ) {
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*By_domain)(i,j,n_d[2]-1-k) = (*Hy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                }
                // for ( int j=0 ; j<n_d[1] ; j++ ) {
                //     No field
                // }
            }
        }

        // Here Change y<->z for PML xMIN
        // 5. Exchange PML z -> PML xMIN
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2-1 ; k++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Primals
                            // i-Duals
                            (*Ex_pml_xmin)(i,j,n_p[2]-1-k) = (*Ex_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dx_pml_xmin)(i,j,n_p[2]-1-k) = (*Dx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Ey_pml_xmin)(i,j,n_p[2]-1-k) = (*Ey_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dy_pml_xmin)(i,j,n_p[2]-1-k) = (*Dy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            // i-Duals
                            (*Hz_pml_xmin)(i,j,n_p[2]-1-k) = (*Hz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Bz_pml_xmin)(i,j,n_p[2]-1-k) = (*Bz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                    }
                }
            }
        }
        // Duals in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Primals
                            (*Ez_pml_xmin)(i,j,n_d[2]-1-k) = (*Ez_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Dz_pml_xmin)(i,j,n_d[2]-1-k) = (*Dz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // i-Duals
                            (*Hy_pml_xmin)(i,j,n_d[2]-1-k) = (*Hy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*By_pml_xmin)(i,j,n_d[2]-1-k) = (*By_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Hx_pml_xmin)(i,j,n_d[2]-1-k) = (*Hx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Bx_pml_xmin)(i,j,n_d[2]-1-k) = (*Bx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // i-Duals
                        }
                    }
                }
            }
        }

        // Here Change y<->z for PML xMax
        // 6. Exchange PML z -> PML xMAX
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2-1 ; k++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = zpml_size_in_x-ncells_pml_domain_xmax ;
                        int jdx_start = ncells_pml_ymin;
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Dual
                            (*Ex_pml_xmax)(i,j,n_p[2]-1-k) = (*Ex_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dx_pml_xmax)(i,j,n_p[2]-1-k) = (*Dx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Ey_pml_xmax)(i,j,n_p[2]-1-k) = (*Ey_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dy_pml_xmax)(i,j,n_p[2]-1-k) = (*Dy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // i-Dual
                            (*Hz_pml_xmax)(i,j,n_p[2]-1-k) = (*Hz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Bz_pml_xmax)(i,j,n_p[2]-1-k) = (*Bz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                    }
                }
            }
        }
        // Dual in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = zpml_size_in_x-ncells_pml_domain_xmax ;
                        int jdx_start = ncells_pml_ymin;
                        for ( int j=0 ; j<n_p[1] ; j++ ) {
                            // i-Primals
                            (*Ez_pml_xmax)(i,j,n_d[2]-1-k) = (*Ez_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Dz_pml_xmax)(i,j,n_d[2]-1-k) = (*Dz_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // i-Dual
                            (*Hy_pml_xmax)(i,j,n_d[2]-1-k) = (*Hy_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*By_pml_xmax)(i,j,n_d[2]-1-k) = (*By_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                        for ( int j=0 ; j<n_d[1] ; j++ ) {
                            // i-Primals
                            (*Hx_pml_xmax)(i,j,n_d[2]-1-k) = (*Hx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Bx_pml_xmax)(i,j,n_d[2]-1-k) = (*Bx_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                    }
                }
            }
        }

        // 7. Exchange PML z -> PML yMin
        // Here Change x<->y with respect to previous block 5.for PML xMin -> PML yMin
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2-1 ; k++){
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int j=0 ; j<ncells_pml_domain_ymin ; j++ ) {
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            // j-Duals
                            (*Ey_pml_ymin)(i,j,n_p[2]-1-k) = (*Ey_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dy_pml_ymin)(i,j,n_p[2]-1-k) = (*Dy_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ex_pml_ymin)(i,j,n_p[2]-1-k) = (*Ex_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dx_pml_ymin)(i,j,n_p[2]-1-k) = (*Dx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // j-Duals
                            (*Hz_pml_ymin)(i,j,n_p[2]-1-k) = (*Hz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Bz_pml_ymin)(i,j,n_p[2]-1-k) = (*Bz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                    }
                }
            }
        }
        // Dual in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int j=0 ; j<ncells_pml_domain_ymin ; j++ ) {
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ez_pml_ymin)(i,j,n_d[2]-1-k) = (*Ez_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Dz_pml_ymin)(i,j,n_d[2]-1-k) = (*Dz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // j-Duals
                            (*Hx_pml_ymin)(i,j,n_d[2]-1-k) = (*Hx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Bx_pml_ymin)(i,j,n_d[2]-1-k) = (*Bx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Hy_pml_ymin)(i,j,n_d[2]-1-k) = (*Hy_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*By_pml_ymin)(i,j,n_d[2]-1-k) = (*By_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // j-Duals
                        }
                    }
                }
            }
        }

        // 8. Exchange PML z -> PML yMax
        // Here Change x<->y with respect to previous block 6.for PML xMax -> PML yMax
        // Primal in z-direction
        for (int k=0 ; k < nsolver/2-1 ; k++){
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    for ( int j=0 ; j<ncells_pml_domain_ymax ; j++ ) {
                        int jdx_start = zpml_size_in_y-ncells_pml_domain_ymax ;
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            // j-Dual
                            (*Ey_pml_ymax)(i,j,n_p[2]-1-k) = (*Ey_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dy_pml_ymax)(i,j,n_p[2]-1-k) = (*Dy_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ex_pml_ymax)(i,j,n_p[2]-1-k) = (*Ex_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Dx_pml_ymax)(i,j,n_p[2]-1-k) = (*Dx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            // j-Duals
                            (*Hz_pml_ymax)(i,j,n_p[2]-1-k) = (*Hz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*Bz_pml_ymax)(i,j,n_p[2]-1-k) = (*Bz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                    }
                }
            }
        }
        // Dual in z-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    for ( int j=0 ; j<ncells_pml_domain_ymax ; j++ ) {
                        int jdx_start = zpml_size_in_y-ncells_pml_domain_ymax ;
                        for ( int i=0 ; i<n_p[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Ez_pml_ymax)(i,j,n_d[2]-1-k) = (*Ez_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Dz_pml_ymax)(i,j,n_d[2]-1-k) = (*Dz_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // j-Duals
                            (*Hx_pml_ymax)(i,j,n_d[2]-1-k) = (*Hx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*Bx_pml_ymax)(i,j,n_d[2]-1-k) = (*Bx_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                        for ( int i=0 ; i<n_d[0]+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            // j-Primals
                            (*Hy_pml_ymax)(i,j,n_d[2]-1-k) = (*Hy_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*By_pml_ymax)(i,j,n_d[2]-1-k) = (*By_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            // j-Duals
                        }
                    }
                }
            }
        }
    }
}
