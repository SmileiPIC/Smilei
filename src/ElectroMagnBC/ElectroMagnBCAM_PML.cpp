#include "ElectroMagnBCAM_PML.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "FieldFactory.h"
#include "cField2D.h"
#include "Tools.h"
#include "Laser.h"
#include <complex>
#include "dcomplex.h"

#include "SolverFactory.h"

using namespace std;

ElectroMagnBCAM_PML::ElectroMagnBCAM_PML( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBCAM( params, patch, i_boundary )
{
    std::vector<unsigned int> n_space(params.n_space);
    std::vector<unsigned int> oversize(params.oversize);
    if( params.multiple_decomposition ) {
        n_space = params.n_space_region;
        oversize = params.region_oversize;
    }

    //Number of modes
    Nmode = params.nmodes;

    pml_solver_ = SolverFactory::createPML( params );
    if (params.maxwell_sol=="Yee"){
        nsolver=2;
        //MESSAGE("FDTD scheme in PML region : Yee.");
    }
    else {
        //WARNING("The solver you use in the main domain is not the same as in the PML region. FDTD scheme in PML region : Yee.");
        nsolver=2;
    }

    El_.resize( Nmode );
    Dl_.resize( Nmode );
    Hl_.resize( Nmode );
    Bl_.resize( Nmode );
    Er_.resize( Nmode );
    Dr_.resize( Nmode );
    Hr_.resize( Nmode );
    Br_.resize( Nmode );
    Et_.resize( Nmode );
    Dt_.resize( Nmode );
    Ht_.resize( Nmode );
    Bt_.resize( Nmode );

    if ( ( i_boundary_ == 0 && patch->isXmin() )
         || ( i_boundary_ == 1 && patch->isXmax() )
         || ( i_boundary_ == 2 && patch->isYmin() )
         || ( i_boundary_ == 3 && patch->isYmax() ) ) {

        int iDim = 0*((i_boundary_==0)||(i_boundary_==1))+1*((i_boundary_==2)||(i_boundary_==3));
        int min_or_max = (i_boundary_)%2;

        domain_oversize_l =  oversize[0] ;
        domain_oversize_r =  oversize[1] ;

        if (patch->isXmin() ) {//&& i_boundary_ == 0 ) {
            ncells_pml_lmin = params.number_of_pml_cells[0][0];
            ncells_pml_domain_lmin = ncells_pml_lmin + 1*oversize[0] + nsolver/2;
            domain_oversize_l = oversize[0] ;
        }
        else {
            ncells_pml_lmin = 0;
            ncells_pml_domain_lmin = 0;
        }
        if (patch->isXmax() ) {//&& i_boundary_ == 1 ) {
            ncells_pml_lmax = params.number_of_pml_cells[0][1];
            ncells_pml_domain_lmax = ncells_pml_lmax + 1*oversize[0] + nsolver/2;
            domain_oversize_l = oversize[0] ;
        }
        else {
            ncells_pml_lmax = 0;
            ncells_pml_domain_lmax = 0;
        }
        if (patch->isYmin() ) {//&& i_boundary_ == 2 ) {
            ncells_pml_rmin = params.number_of_pml_cells[1][0];
            ncells_pml_domain_rmin = ncells_pml_rmin + 1*oversize[1] + nsolver/2;
            domain_oversize_r = oversize[1] ;
        }
        else {
            ncells_pml_rmin = 0;
            ncells_pml_domain_rmin = 0;
        }
        if (patch->isYmax() ) {//&& i_boundary_ == 3 ) {
            ncells_pml_rmax = params.number_of_pml_cells[1][1];
            ncells_pml_domain_rmax = ncells_pml_rmax + 1*oversize[1] + nsolver/2;
            domain_oversize_r = oversize[1] ;
        }
        else {
            ncells_pml_rmax = 0;
            ncells_pml_domain_rmax = 0;
        }

        ncells_pml = params.number_of_pml_cells[iDim][min_or_max];
        ncells_pml_domain = ncells_pml+1*oversize[iDim]+ nsolver/2;

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
        if ( iDim==1 ){
            dimPrim[iDim-1] += ncells_pml_lmin + ncells_pml_lmax ;
            rpml_size_in_l = dimPrim[iDim-1] ;
        }

        int startpml = oversize[iDim]+nsolver/2;

        int ncells_pml_min[1];
        ncells_pml_min[0] = ncells_pml_lmin;
        int ncells_pml_max[1];
        ncells_pml_max[0] = ncells_pml_lmax;

        pml_solver_->setDomainSizeAndCoefficients( iDim, min_or_max, ncells_pml_domain, startpml, ncells_pml_min, ncells_pml_max, patch );

        for ( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            ostringstream mode_id( "" );
            mode_id << "_mode_" << imode;
            El_[imode] = FieldFactory::createComplex( dimPrim, 0, false, ( "El_pml_"+mode_id.str() ).c_str(), params );
            Dl_[imode] = FieldFactory::createComplex( dimPrim, 0, false, ( "Dl_pml_"+mode_id.str() ).c_str(), params );
            Hl_[imode] = FieldFactory::createComplex( dimPrim, 0, true,  ( "Hl_pml_"+mode_id.str() ).c_str(), params );
            Bl_[imode] = FieldFactory::createComplex( dimPrim, 0, true,  ( "Bl_pml_"+mode_id.str() ).c_str(), params );
            Er_[imode] = FieldFactory::createComplex( dimPrim, 1, false, ( "Er_pml_"+mode_id.str() ).c_str(), params );
            Dr_[imode] = FieldFactory::createComplex( dimPrim, 1, false, ( "Dr_pml_"+mode_id.str() ).c_str(), params );
            Hr_[imode] = FieldFactory::createComplex( dimPrim, 1, true,  ( "Hr_pml_"+mode_id.str() ).c_str(), params );
            Br_[imode] = FieldFactory::createComplex( dimPrim, 1, true,  ( "Br_pml_"+mode_id.str() ).c_str(), params );
            Et_[imode] = FieldFactory::createComplex( dimPrim, 2, false, ( "Et_pml_"+mode_id.str() ).c_str(), params );
            Dt_[imode] = FieldFactory::createComplex( dimPrim, 2, false, ( "Dt_pml_"+mode_id.str() ).c_str(), params );
            Ht_[imode] = FieldFactory::createComplex( dimPrim, 2, true,  ( "Ht_pml_"+mode_id.str() ).c_str(), params );
            Bt_[imode] = FieldFactory::createComplex( dimPrim, 2, true,  ( "Bt_pml_"+mode_id.str() ).c_str(), params );
        }

    //Laser parameter
    double pyKx, pyKy; //, pyKz;
    double kl, kr; //, kz;
    double Knorm;
    double omega = 1. ;

    factor_laser_space_time = 2.*dt_ov_d[0] ;

    // Xmin boundary
    pyKx = +1.;
    pyKy = 0.;
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
    kl = omega*pyKx/Knorm;
    kr = omega*pyKy/Knorm;

    factor_laser_angle_W = kl/Knorm;

    // Xmax boundary
    pyKx = -1.;
    pyKy = 0.;
    Knorm = sqrt( pyKx*pyKx + pyKy*pyKy ) ;
    kl = omega*pyKx/Knorm;
    kr = omega*pyKy/Knorm;

    factor_laser_angle_E = kl/Knorm;

    }
}


ElectroMagnBCAM_PML::~ElectroMagnBCAM_PML()
{
    for ( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
        delete El_[imode];
        delete Dl_[imode];
        delete Hl_[imode];
        delete Bl_[imode];
        delete Er_[imode];
        delete Dr_[imode];
        delete Hr_[imode];
        delete Br_[imode];
        delete Et_[imode];
        delete Dt_[imode];
        delete Ht_[imode];
        delete Bt_[imode];
    }
    if (pml_solver_!=NULL) {
         delete pml_solver_;
    }
}


void ElectroMagnBCAM_PML::save_fields( Field *my_field, Patch *patch )
{
}


void ElectroMagnBCAM_PML::disableExternalFields()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCAM_PML::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    int iDim = 0*((i_boundary_==0)||(i_boundary_==1))+1*((i_boundary_==2)||(i_boundary_==3));
    int min_or_max = (i_boundary_)%2;

    if( i_boundary_ == 0 && patch->isXmin() ) {

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];

            // 2. Exchange field PML <- Domain
            for ( int i=min2exchange ; i<max2exchange ; i++ ) {
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Er_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Er_domain)(i,j);
                    (*Dr_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Er_domain)(i,j);
                    (*Hl_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Bl_domain)(i,j);
                    (*Bl_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Bl_domain)(i,j);
                    (*Ht_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Bt_domain)(i,j);
                    (*Bt_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Bt_domain)(i,j);
                }
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Br_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Br_domain)(i,j);
                    (*Hr_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Br_domain)(i,j);
                    (*El_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*El_domain)(i,j);
                    (*Dl_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*El_domain)(i,j);
                    (*Et_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Et_domain)(i,j);
                    (*Dt_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*Et_domain)(i,j);
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            //cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            //cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            //cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            //cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            //cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            //cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
            // for Br^(d,p)
            vector<double> yp( 1 );
            for( unsigned int j=3*( patch->isYmin() ) ; j<n_p[1]-patch->isYmax() ; j++ ) {
                
                std::complex<double> byW = 0.;
                yp[0] = patch->getDomainLocalMin( 1 ) +( (double)j - (double)EMfields->oversize[1] )*d[1];
                
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    if (vecLaser[ilaser]->spacetime.size() > 2){
                        byW += vecLaser[ilaser]->getAmplitudecomplexN(yp, time_dual, 0, 0, 2*imode);
                    } else {
                        if( imode==1 ) {
                            byW +=   vecLaser[ilaser]->getAmplitude0( yp, time_dual, 1+2*j, 0 )
                                + Icpx * vecLaser[ilaser]->getAmplitude1( yp, time_dual, 1+2*j, 0 );
                        }
                    }
                }
                // unsigned int i=120;
                // ( *Br_domain )( i, j ) += factor_laser_space_time*factor_laser_angle_W*byW;
                unsigned int i=ncells_pml_domain-domain_oversize_l-nsolver/2;
                ( *Hr_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_W*byW;
                ( *Br_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_W*byW;
            }
            // for Bt^(d,d)
            vector<double> yd( 1 );
            for( unsigned int j=3*( patch->isYmin() ); j<n_d[1]-patch->isYmax() ; j++ ) {

                std::complex<double> bzW = 0.;
                yd[0] = patch->getDomainLocalMin( 1 ) + ( (double)j - 0.5 - (double)EMfields->oversize[1] )*d[1];
                
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    if (vecLaser[ilaser]->spacetime.size() > 2){
                        bzW += vecLaser[ilaser]->getAmplitudecomplexN(yd, time_dual, 0, 0, 2*imode+1);
                    } else {
                        if( imode==1 ) {
                            bzW +=   vecLaser[ilaser]->getAmplitude1( yd, time_dual, 2*j, 0 )
                                   - Icpx * vecLaser[ilaser]->getAmplitude0( yd, time_dual, 2*j, 0 );
                        }
                    }
                }
                // unsigned int i=120;
                // ( *Bt_domain )( i, j ) += factor_laser_space_time*factor_laser_angle_W*bzW;
                unsigned int i=ncells_pml_domain-domain_oversize_l-nsolver/2;
                ( *Ht_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_W*bzW;
                ( *Bt_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_W*bzW;
            }
            //Redo condition on axis for Bt because it was modified
            if( patch->isYmin() && imode != 1 ) {
                unsigned int i=ncells_pml_domain-domain_oversize_l-nsolver/2;
                ( *Bt_[imode] )( i, 2 ) = -( *Bt_[imode] )( i, 3 );
                ( *Ht_[imode] )( i, 2 ) = -( *Ht_[imode] )( i, 3 );
            }
            if( patch->isYmin() && imode == 1 ) {
                unsigned int i=ncells_pml_domain-domain_oversize_l-nsolver/2;
                ( *Bt_[imode] )( i, 2 )= -2.*Icpx*( *Br_[imode] )( i, 2 )-( *Bt_[imode] )( i, 3 );
                ( *Ht_[imode] )( i, 2 )= -2.*Icpx*( *Hr_[imode] )( i, 2 )-( *Ht_[imode] )( i, 3 );
            }
        }

        // 4. Exchange PML -> Domain
        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
            // Primals in x-direction
            for (int i=0 ; i < nsolver/2 ; i++){
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Et_domain)(i,j) = (*Et_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Er_domain)(i,j) = (*Er_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                    (*Bl_domain)(i,j) = (*Hl_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                }
            }
            // Duals in x-direction
            for (int i=0 ; i < nsolver/2 ; i++){
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*El_domain)(i,j) = (*El_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                    (*Br_domain)(i,j) = (*Hr_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Bt_domain)(i,j) = (*Ht_[imode])(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                }
            }
        }
    }
    else if( i_boundary_ == 1 && patch->isXmax() ) {

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];

            for ( int i=min2exchange ; i<max2exchange ; i++ ) {
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Er_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Er_domain)(n_p[0]-i,j);
                    (*Dr_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Er_domain)(n_p[0]-i,j);
                    (*Hl_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Bl_domain)(n_p[0]-i,j);
                    (*Bl_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Bl_domain)(n_p[0]-i,j);
                    (*Ht_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Bt_domain)(n_p[0]-i,j);
                    (*Bt_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Bt_domain)(n_p[0]-i,j);
                }
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Br_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Br_domain)(n_p[0]-i,j);
                    (*Hr_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Br_domain)(n_p[0]-i,j);
                    (*El_[imode])(domain_oversize_l+nsolver/2-i,j) = (*El_domain)(n_p[0]-i,j);
                    (*Dl_[imode])(domain_oversize_l+nsolver/2-i,j) = (*El_domain)(n_p[0]-i,j);
                    (*Et_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Et_domain)(n_p[0]-i,j);
                    (*Dt_[imode])(domain_oversize_l+nsolver/2-i,j) = (*Et_domain)(n_p[0]-i,j);
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        //Injecting a laser
        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            //cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            //cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            //cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            //cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            //cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            //cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
            // for Br^(d,p)
            vector<double> yp( 1 );
            for( unsigned int j=3*( patch->isYmin() ) ; j<n_p[1] ; j++ ) {
                
                std::complex<double> byE = 0.;
                yp[0] = patch->getDomainLocalMin( 1 ) +( (double)j - (double)EMfields->oversize[1] )*d[1];
                
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    if (vecLaser[ilaser]->spacetime.size() > 2){
                        byE += vecLaser[ilaser]->getAmplitudecomplexN(yp, time_dual, 0, 0, 2*imode);
                    } else {
                        if( imode==1 ) {
                            byE +=   vecLaser[ilaser]->getAmplitude0( yp, time_dual, 1+2*j, 0 )
                                + Icpx * vecLaser[ilaser]->getAmplitude1( yp, time_dual, 1+2*j, 0 );
                        }
                    }
                }
                // unsigned int i=120;
                // ( *Br_domain )( i, j ) += factor_laser_space_time*factor_laser_angle_E*byE;
                unsigned int i=domain_oversize_l+nsolver/2;
                ( *Hr_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_E*byE;
                ( *Br_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_E*byE;
            }
            // for Bt^(d,d)
            vector<double> yd( 1 );
            for( unsigned int j=3*( patch->isYmin() ); j<n_d[1] ; j++ ) {

                std::complex<double> bzE = 0.;
                yd[0] = patch->getDomainLocalMin( 1 ) + ( (double)j - 0.5 - (double)EMfields->oversize[1] )*d[1];
                
                // Lasers
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    if (vecLaser[ilaser]->spacetime.size() > 2){
                        bzE += vecLaser[ilaser]->getAmplitudecomplexN(yd, time_dual, 0, 0, 2*imode+1);
                    } else {
                        if( imode==1 ) {
                            bzE +=   vecLaser[ilaser]->getAmplitude1( yd, time_dual, 2*j, 0 )
                                   - Icpx * vecLaser[ilaser]->getAmplitude0( yd, time_dual, 2*j, 0 );
                        }
                    }
                }
                // unsigned int i=120;
                // ( *Bt_domain )( i, j ) += factor_laser_space_time*factor_laser_angle_E*bzE;
                unsigned int i=domain_oversize_l+nsolver/2;
                ( *Ht_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_E*bzE;
                ( *Bt_[imode] )( i, j ) += factor_laser_space_time*factor_laser_angle_E*bzE;
            }
            //Redo condition on axis for Bt because it was modified
            if( patch->isYmin() && imode != 1 ) {
                unsigned int i=domain_oversize_l+nsolver/2;
                ( *Bt_[imode] )( i, 2 ) = -( *Bt_[imode] )( i, 3 );
                ( *Ht_[imode] )( i, 2 ) = -( *Ht_[imode] )( i, 3 );
            }
            if( patch->isYmin() && imode == 1 ) {
                unsigned int i=domain_oversize_l+nsolver/2;
                ( *Bt_[imode] )( i, 2 )= -2.*Icpx*( *Br_[imode] )( i, 2 )-( *Bt_[imode] )( i, 3 );
                ( *Ht_[imode] )( i, 2 )= -2.*Icpx*( *Hr_[imode] )( i, 2 )-( *Ht_[imode] )( i, 3 );
            }
        }

        // 4. Exchange Domain -> PML
        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];
            // Primals in x-direction
            for (int i=0 ; i < nsolver/2-1 ; i++){
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*Et_domain)(n_p[0]-1-i,j) = (*Et_[imode])(domain_oversize_l+nsolver/2-1-i,j);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Er_domain)(n_p[0]-1-i,j) = (*Er_[imode])(domain_oversize_l+nsolver/2-1-i,j);
                    (*Bl_domain)(n_p[0]-1-i,j) = (*Hl_[imode])(domain_oversize_l+nsolver/2-1-i,j);
                }
            }
            // Duals in x-direction
            for (int i=0 ; i < nsolver/2 ; i++){
                for ( int j=0 ; j<n_p[1] ; j++ ) {
                    (*El_domain)(n_d[0]-1-i,j) = (*El_[imode])(domain_oversize_l+nsolver/2-i,j);
                    (*Br_domain)(n_d[0]-1-i,j) = (*Hr_[imode])(domain_oversize_l+nsolver/2-i,j);
                }
                for ( int j=0 ; j<n_d[1] ; j++ ) {
                    (*Bt_domain)(n_d[0]-1-i,j) = (*Ht_[imode])(domain_oversize_l+nsolver/2-i,j);
                }
            }
        }
    }
    else if( i_boundary_ == 2 && patch->isYmin() ) {
        // // ERROR("PML not allow on the symetric axis")

        // ElectroMagnBCAM_PML* pml_fields_lmin = NULL ;
        // ElectroMagnBCAM_PML* pml_fields_lmax = NULL ;

        // if(ncells_pml_lmin != 0){
        //     pml_fields_lmin = static_cast<ElectroMagnBCAM_PML*>( EMfields->emBoundCond[0] );
        // }
        // if(ncells_pml_lmax != 0){
        //     pml_fields_lmax = static_cast<ElectroMagnBCAM_PML*>( EMfields->emBoundCond[1] );
        // }

        // cField2D* El_pml_lmin = NULL;
        // cField2D* Er_pml_lmin = NULL;
        // cField2D* Et_pml_lmin = NULL;
        // cField2D* Hl_pml_lmin = NULL;
        // cField2D* Hr_pml_lmin = NULL;
        // cField2D* Ht_pml_lmin = NULL;
        // cField2D* Dl_pml_lmin = NULL;
        // cField2D* Dr_pml_lmin = NULL;
        // cField2D* Dt_pml_lmin = NULL;
        // cField2D* Bl_pml_lmin = NULL;
        // cField2D* Br_pml_lmin = NULL;
        // cField2D* Bt_pml_lmin = NULL;

        // cField2D* El_pml_lmax = NULL;
        // cField2D* Er_pml_lmax = NULL;
        // cField2D* Et_pml_lmax = NULL;
        // cField2D* Hl_pml_lmax = NULL;
        // cField2D* Hr_pml_lmax = NULL;
        // cField2D* Ht_pml_lmax = NULL;
        // cField2D* Dl_pml_lmax = NULL;
        // cField2D* Dr_pml_lmax = NULL;
        // cField2D* Dt_pml_lmax = NULL;
        // cField2D* Bl_pml_lmax = NULL;
        // cField2D* Br_pml_lmax = NULL;
        // cField2D* Bt_pml_lmax = NULL;

        // // 1. Solve Maxwell_PML for E-field :
        // // As if B-field isn't updated
        // pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        // //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        // for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
        //     cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        //     cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        //     cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        //     cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        //     cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        //     cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];

        //     if(ncells_pml_lmin != 0){
        //         El_pml_lmin = pml_fields_lmin->El_[imode];
        //         Er_pml_lmin = pml_fields_lmin->Er_[imode];
        //         Et_pml_lmin = pml_fields_lmin->Et_[imode];
        //         Hl_pml_lmin = pml_fields_lmin->Hl_[imode];
        //         Hr_pml_lmin = pml_fields_lmin->Hr_[imode];
        //         Ht_pml_lmin = pml_fields_lmin->Ht_[imode];
        //         Dl_pml_lmin = pml_fields_lmin->Dl_[imode];
        //         Dr_pml_lmin = pml_fields_lmin->Dr_[imode];
        //         Dt_pml_lmin = pml_fields_lmin->Dt_[imode];
        //         Bl_pml_lmin = pml_fields_lmin->Bl_[imode];
        //         Br_pml_lmin = pml_fields_lmin->Br_[imode];
        //         Bt_pml_lmin = pml_fields_lmin->Bt_[imode];
        //     }

        //     if(ncells_pml_lmax != 0){
        //         El_pml_lmax = pml_fields_lmax->El_[imode];
        //         Er_pml_lmax = pml_fields_lmax->Er_[imode];
        //         Et_pml_lmax = pml_fields_lmax->Et_[imode];
        //         Hl_pml_lmax = pml_fields_lmax->Hl_[imode];
        //         Hr_pml_lmax = pml_fields_lmax->Hr_[imode];
        //         Ht_pml_lmax = pml_fields_lmax->Ht_[imode];
        //         Dl_pml_lmax = pml_fields_lmax->Dl_[imode];
        //         Dr_pml_lmax = pml_fields_lmax->Dr_[imode];
        //         Dt_pml_lmax = pml_fields_lmax->Dt_[imode];
        //         Bl_pml_lmax = pml_fields_lmax->Bl_[imode];
        //         Br_pml_lmax = pml_fields_lmax->Br_[imode];
        //         Bt_pml_lmax = pml_fields_lmax->Bt_[imode];
        //     }

        //     // 2. Exchange field PML <- Domain
        //     for ( int j=min2exchange ; j<max2exchange ; j++ ) {
        //         if (patch->isXmin()) {
        //             if(ncells_pml_lmin != 0){
        //                 for ( int i=0 ; i<ncells_pml_lmin ; i++ ) {
        //                     int idl_start = 0;
        //                     // Les qtes Primals
        //                     (*Bl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bl_pml_lmin)(i,j);
        //                     (*Hl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Hl_pml_lmin)(i,j);
        //                     (*Er_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Er_pml_lmin)(i,j);
        //                     (*Dr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Dr_pml_lmin)(i,j);
        //                     (*Et_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Et_pml_lmin)(i,j);
        //                     (*Dt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Dt_pml_lmin)(i,j);
        //                     // Les qtes Duals
        //                     (*El_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*El_pml_lmin)(i,j);
        //                     (*Dl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Dl_pml_lmin)(i,j);
        //                     (*Hr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Hr_pml_lmin)(i,j);
        //                     (*Br_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Br_pml_lmin)(i,j);
        //                     (*Ht_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Ht_pml_lmin)(i,j);
        //                     (*Bt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bt_pml_lmin)(i,j);
        //                 }
        //             }
        //         }
        //         for ( int i=0 ; i<n_d[0] ; i++ ) {
        //             int idl_start = ncells_pml_lmin;
        //             (*El_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*El_domain)(i,j);
        //             (*Dl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*El_domain)(i,j);
        //             (*Hr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Br_domain)(i,j);
        //             (*Br_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Br_domain)(i,j);
        //             (*Ht_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bt_domain)(i,j);
        //             (*Bt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bt_domain)(i,j);
        //         }
        //         for ( int i=0 ; i<n_p[0] ; i++ ) {
        //             int idl_start = ncells_pml_lmin;
        //             (*Hl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bl_domain)(i,j);
        //             (*Bl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bl_domain)(i,j);
        //             (*Er_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Er_domain)(i,j);
        //             (*Dr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Er_domain)(i,j);
        //             (*Et_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Et_domain)(i,j);
        //             (*Dt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Et_domain)(i,j);
        //         }
        //         if (patch->isXmax()) {
        //             if(ncells_pml_lmax != 0){
        //                 for ( int i=0 ; i<ncells_pml_lmax ; i++ ) {
        //                     int idl_start = (rpml_size_in_l-1)-(ncells_pml_lmax-1) ;
        //                     // Les qtes Primals commencent a (rpml_size_in_l+1)-(ncells_pml_lmax-1)
        //                     (*Bl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bl_pml_lmax)(domain_oversize_l+nsolver/2+i,j);
        //                     (*Hl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Hl_pml_lmax)(domain_oversize_l+nsolver/2+i,j);
        //                     (*Er_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Er_pml_lmax)(domain_oversize_l+nsolver/2+i,j);
        //                     (*Dr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Dr_pml_lmax)(domain_oversize_l+nsolver/2+i,j);
        //                     (*Et_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Et_pml_lmax)(domain_oversize_l+nsolver/2+i,j);
        //                     (*Dt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Dt_pml_lmax)(domain_oversize_l+nsolver/2+i,j);
        //                     // Toutes les qtes Duals commence a +1
        //                     (*El_[imode])(idl_start+1+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*El_pml_lmax)(domain_oversize_l+nsolver/2+1+i,j);
        //                     (*Dl_[imode])(idl_start+1+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Dl_pml_lmax)(domain_oversize_l+nsolver/2+1+i,j);
        //                     (*Hr_[imode])(idl_start+1+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Hr_pml_lmax)(domain_oversize_l+nsolver/2+1+i,j);
        //                     (*Br_[imode])(idl_start+1+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Br_pml_lmax)(domain_oversize_l+nsolver/2+1+i,j);
        //                     (*Ht_[imode])(idl_start+1+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Ht_pml_lmax)(domain_oversize_l+nsolver/2+1+i,j);
        //                     (*Bt_[imode])(idl_start+1+i,ncells_pml_domain-domain_oversize_r-nsolver/2+j) = (*Bt_pml_lmax)(domain_oversize_l+nsolver/2+1+i,j);
        //                 }
        //             }
        //         }
        //     }
        // }

        // // 3. Solve Maxwell_PML for B-field :
        // //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        // pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        // for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
        //     cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        //     cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        //     cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        //     cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
        //     cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
        //     cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];

        //     if (patch->isXmin()) {
        //         if(ncells_pml_lmin != 0){
        //             El_pml_lmin = pml_fields_lmin->El_[imode];
        //             Er_pml_lmin = pml_fields_lmin->Er_[imode];
        //             Et_pml_lmin = pml_fields_lmin->Et_[imode];
        //             Hl_pml_lmin = pml_fields_lmin->Hl_[imode];
        //             Hr_pml_lmin = pml_fields_lmin->Hr_[imode];
        //             Ht_pml_lmin = pml_fields_lmin->Ht_[imode];
        //             Dl_pml_lmin = pml_fields_lmin->Dl_[imode];
        //             Dr_pml_lmin = pml_fields_lmin->Dr_[imode];
        //             Dt_pml_lmin = pml_fields_lmin->Dt_[imode];
        //             Bl_pml_lmin = pml_fields_lmin->Bl_[imode];
        //             Br_pml_lmin = pml_fields_lmin->Br_[imode];
        //             Bt_pml_lmin = pml_fields_lmin->Bt_[imode];
        //         }
        //     }

        //     if (patch->isXmax()) {
        //         if(ncells_pml_lmax != 0){
        //             El_pml_lmax = pml_fields_lmax->El_[imode];
        //             Er_pml_lmax = pml_fields_lmax->Er_[imode];
        //             Et_pml_lmax = pml_fields_lmax->Et_[imode];
        //             Hl_pml_lmax = pml_fields_lmax->Hl_[imode];
        //             Hr_pml_lmax = pml_fields_lmax->Hr_[imode];
        //             Ht_pml_lmax = pml_fields_lmax->Ht_[imode];
        //             Dl_pml_lmax = pml_fields_lmax->Dl_[imode];
        //             Dr_pml_lmax = pml_fields_lmax->Dr_[imode];
        //             Dt_pml_lmax = pml_fields_lmax->Dt_[imode];
        //             Bl_pml_lmax = pml_fields_lmax->Bl_[imode];
        //             Br_pml_lmax = pml_fields_lmax->Br_[imode];
        //             Bt_pml_lmax = pml_fields_lmax->Bt_[imode];
        //         }
        //     }

        //     // 4. Exchange PML -> Domain
        //     // Primals in y-direction
        //     for (int j=0 ; j < nsolver/2 ; j++){
        //         for ( int i=0 ; i<n_p[0] ; i++ ) {
        //             int idl_start = ncells_pml_lmin;
        //             (*Et_domain)(i,j) = (*Et_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //         }
        //         for ( int i=0 ; i<n_d[0] ; i++ ) {
        //             int idl_start = ncells_pml_lmin;
        //             (*El_domain)(i,j) = (*El_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //             (*Br_domain)(i,j) = (*Hr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //         }
        //     }
        //     // Duals in y-direction
        //     for (int j=0 ; j < nsolver/2 ; j++){
        //         for ( int i=0 ; i<n_p[0] ; i++ ) {
        //             int idl_start = ncells_pml_lmin;
        //             (*Er_domain)(i,j) = (*Er_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //             (*Bl_domain)(i,j) = (*Hl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //         }
        //         for ( int i=0 ; i<n_d[0] ; i++ ) {
        //             int idl_start = ncells_pml_lmin;
        //             (*Bt_domain)(i,j) = (*Ht_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //         }
        //     }

        //     // 5. Exchange PML y -> PML x MIN
        //     // Primal in y-direction
        //     for (int j=0 ; j < nsolver/2 ; j++){
        //         if (patch->isXmin()) {
        //             if(ncells_pml_lmin != 0){
        //                 for ( int i=0 ; i<ncells_pml_domain_lmin ; i++ ) {
        //                     int idl_start = 0;
        //                     // Primals
        //                     (*Et_pml_lmin)(i,j) = (*Et_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Dt_pml_lmin)(i,j) = (*Dt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     // Duals
        //                     (*El_pml_lmin)(i,j) = (*El_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Dl_pml_lmin)(i,j) = (*Dl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Hr_pml_lmin)(i,j) = (*Hr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Br_pml_lmin)(i,j) = (*Br_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                 }
        //             }
        //         }
        //     }
        //     // Duals in y-direction
        //     for (int j=0 ; j < nsolver/2 ; j++){
        //         if (patch->isXmin()) {
        //             if(ncells_pml_lmin != 0){
        //                 for ( int i=0 ; i<ncells_pml_domain_lmin ; i++ ) {
        //                     int idl_start = 0;
        //                     // Primals
        //                     (*Er_pml_lmin)(i,j) = (*Er_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Dr_pml_lmin)(i,j) = (*Dr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Hl_pml_lmin)(i,j) = (*Hl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Bl_pml_lmin)(i,j) = (*Bl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     // Duals
        //                     (*Ht_pml_lmin)(i,j) = (*Ht_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Bt_pml_lmin)(i,j) = (*Bt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                 }
        //             }
        //         }
        //     }

        //     // 6. Exchange PML y -> PML x MAX
        //     // Primal in y-direction
        //     for (int j=0 ; j < nsolver/2 ; j++){
        //         if (patch->isXmax()) {
        //             if(ncells_pml_lmax != 0){
        //                 for ( int i=0 ; i<ncells_pml_domain_lmax ; i++ ) {
        //                     int idl_start = rpml_size_in_l-ncells_pml_domain_lmax ;
        //                     // Primals
        //                     (*Et_pml_lmax)(i,j) = (*Et_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Dt_pml_lmax)(i,j) = (*Dt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     // Dual
        //                     (*El_pml_lmax)(i,j) = (*El_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Dl_pml_lmax)(i,j) = (*Dl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Hr_pml_lmax)(i,j) = (*Hr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Br_pml_lmax)(i,j) = (*Br_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                 }
        //             }
        //         }
        //     }
        //     // Dual in y-direction
        //     for (int j=0 ; j < nsolver/2 ; j++){
        //         if (patch->isXmax()) {
        //             if(ncells_pml_lmax != 0){
        //                 for ( int i=0 ; i<ncells_pml_domain_lmax ; i++ ) {
        //                     int idl_start = rpml_size_in_l-ncells_pml_domain_lmax ;
        //                     // Primals
        //                     (*Er_pml_lmax)(i,j) = (*Er_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Dr_pml_lmax)(i,j) = (*Dr_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Hl_pml_lmax)(i,j) = (*Hl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Bl_pml_lmax)(i,j) = (*Bl_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     // Dual
        //                     (*Ht_pml_lmax)(i,j) = (*Ht_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                     (*Bt_pml_lmax)(i,j) = (*Bt_[imode])(idl_start+i,ncells_pml_domain-domain_oversize_l-nsolver/2+j);
        //                 }
        //             }
        //         }
        //     }
        // }
    }
    else if( i_boundary_ == 3 && patch->isYmax() ) {

        ElectroMagnBCAM_PML* pml_fields_lmin = NULL ;
        ElectroMagnBCAM_PML* pml_fields_lmax = NULL ;

        if(ncells_pml_lmin != 0){
            pml_fields_lmin = static_cast<ElectroMagnBCAM_PML*>( EMfields->emBoundCond[0] );
        }
        if(ncells_pml_lmax != 0){
            pml_fields_lmax = static_cast<ElectroMagnBCAM_PML*>( EMfields->emBoundCond[1] );
        }

        cField2D* El_pml_lmin = NULL;
        cField2D* Er_pml_lmin = NULL;
        cField2D* Et_pml_lmin = NULL;
        cField2D* Hl_pml_lmin = NULL;
        cField2D* Hr_pml_lmin = NULL;
        cField2D* Ht_pml_lmin = NULL;
        cField2D* Dl_pml_lmin = NULL;
        cField2D* Dr_pml_lmin = NULL;
        cField2D* Dt_pml_lmin = NULL;
        cField2D* Bl_pml_lmin = NULL;
        cField2D* Br_pml_lmin = NULL;
        cField2D* Bt_pml_lmin = NULL;

        cField2D* El_pml_lmax = NULL;
        cField2D* Er_pml_lmax = NULL;
        cField2D* Et_pml_lmax = NULL;
        cField2D* Hl_pml_lmax = NULL;
        cField2D* Hr_pml_lmax = NULL;
        cField2D* Ht_pml_lmax = NULL;
        cField2D* Dl_pml_lmax = NULL;
        cField2D* Dr_pml_lmax = NULL;
        cField2D* Dt_pml_lmax = NULL;
        cField2D* Bl_pml_lmax = NULL;
        cField2D* Br_pml_lmax = NULL;
        cField2D* Bt_pml_lmax = NULL;

        // 1. Solve Maxwell_PML for E-field :
        // As if B-field isn't updated
        pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        //pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];

            if(ncells_pml_lmin != 0){
                El_pml_lmin = pml_fields_lmin->El_[imode];
                Er_pml_lmin = pml_fields_lmin->Er_[imode];
                Et_pml_lmin = pml_fields_lmin->Et_[imode];
                Hl_pml_lmin = pml_fields_lmin->Hl_[imode];
                Hr_pml_lmin = pml_fields_lmin->Hr_[imode];
                Ht_pml_lmin = pml_fields_lmin->Ht_[imode];
                Dl_pml_lmin = pml_fields_lmin->Dl_[imode];
                Dr_pml_lmin = pml_fields_lmin->Dr_[imode];
                Dt_pml_lmin = pml_fields_lmin->Dt_[imode];
                Bl_pml_lmin = pml_fields_lmin->Bl_[imode];
                Br_pml_lmin = pml_fields_lmin->Br_[imode];
                Bt_pml_lmin = pml_fields_lmin->Bt_[imode];
            }

            if(ncells_pml_lmax != 0){
                El_pml_lmax = pml_fields_lmax->El_[imode];
                Er_pml_lmax = pml_fields_lmax->Er_[imode];
                Et_pml_lmax = pml_fields_lmax->Et_[imode];
                Hl_pml_lmax = pml_fields_lmax->Hl_[imode];
                Hr_pml_lmax = pml_fields_lmax->Hr_[imode];
                Ht_pml_lmax = pml_fields_lmax->Ht_[imode];
                Dl_pml_lmax = pml_fields_lmax->Dl_[imode];
                Dr_pml_lmax = pml_fields_lmax->Dr_[imode];
                Dt_pml_lmax = pml_fields_lmax->Dt_[imode];
                Bl_pml_lmax = pml_fields_lmax->Bl_[imode];
                Br_pml_lmax = pml_fields_lmax->Br_[imode];
                Bt_pml_lmax = pml_fields_lmax->Bt_[imode];
            }

            for ( int j=min2exchange ; j<max2exchange ; j++ ) {
                if (patch->isXmin()) {
                    if(ncells_pml_lmin != 0){
                        for ( int i=0 ; i<ncells_pml_lmin ; i++ ) {
                            int idl_start = 0;
                            // Les qtes Primals
                            (*Bl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Bl_pml_lmin)(i,n_p[1]-j);
                            (*Hl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Hl_pml_lmin)(i,n_p[1]-j);
                            (*Er_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Er_pml_lmin)(i,n_p[1]-j);
                            (*Dr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Dr_pml_lmin)(i,n_p[1]-j);
                            (*Et_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Et_pml_lmin)(i,n_p[1]-j);
                            (*Dt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Dt_pml_lmin)(i,n_p[1]-j);
                            // Les qtes Duals
                            (*El_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*El_pml_lmin)(i,n_p[1]-j);
                            (*Dl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Dl_pml_lmin)(i,n_p[1]-j);
                            (*Hr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Hr_pml_lmin)(i,n_p[1]-j);
                            (*Br_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Br_pml_lmin)(i,n_p[1]-j);
                            (*Ht_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Ht_pml_lmin)(i,n_p[1]-j);
                            (*Bt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Bt_pml_lmin)(i,n_p[1]-j);
                        }
                    }
                }
                for ( int i=0 ; i<n_d[0] ; i++ ) {
                    int idl_start = ncells_pml_lmin;
                    (*El_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*El_domain)(i,n_p[1]-j);
                    (*Dl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*El_domain)(i,n_p[1]-j);
                    (*Hr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Br_domain)(i,n_p[1]-j);
                    (*Br_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Br_domain)(i,n_p[1]-j);
                    (*Ht_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Bt_domain)(i,n_p[1]-j);
                    (*Bt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Bt_domain)(i,n_p[1]-j);
                }
                for ( int i=0 ; i<n_p[0] ; i++ ) {
                    int idl_start = ncells_pml_lmin;
                    (*Hl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Bl_domain)(i,n_p[1]-j);
                    (*Bl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Bl_domain)(i,n_p[1]-j);
                    (*Er_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Er_domain)(i,n_p[1]-j);
                    (*Dr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Er_domain)(i,n_p[1]-j);
                    (*Et_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Et_domain)(i,n_p[1]-j);
                    (*Dt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Et_domain)(i,n_p[1]-j);
                }
                if (patch->isXmax()) {
                    if(ncells_pml_lmax != 0){
                        for ( int i=0 ; i<ncells_pml_lmax ; i++ ) {
                            int idl_start = (rpml_size_in_l-1)-(ncells_pml_lmax-1) ;
                            // Les qtes Primals commencent a (rpml_size_in_l+1)-(ncells_pml_lmax-1)
                            (*Bl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Bl_pml_lmax)(domain_oversize_l+nsolver/2+i,n_p[1]-j);
                            (*Hl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Hl_pml_lmax)(domain_oversize_l+nsolver/2+i,n_p[1]-j);
                            (*Er_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Er_pml_lmax)(domain_oversize_l+nsolver/2+i,n_p[1]-j);
                            (*Dr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Dr_pml_lmax)(domain_oversize_l+nsolver/2+i,n_p[1]-j);
                            (*Et_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Et_pml_lmax)(domain_oversize_l+nsolver/2+i,n_p[1]-j);
                            (*Dt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j) = (*Dt_pml_lmax)(domain_oversize_l+nsolver/2+i,n_p[1]-j);
                            // Toutes les qtes Duals commence a +1
                            (*El_[imode])(idl_start+1+i,domain_oversize_r+nsolver/2-j) = (*El_pml_lmax)(domain_oversize_l+nsolver/2+1+i,n_p[1]-j);
                            (*Dl_[imode])(idl_start+1+i,domain_oversize_r+nsolver/2-j) = (*Dl_pml_lmax)(domain_oversize_l+nsolver/2+1+i,n_p[1]-j);
                            (*Hr_[imode])(idl_start+1+i,domain_oversize_r+nsolver/2-j) = (*Hr_pml_lmax)(domain_oversize_l+nsolver/2+1+i,n_p[1]-j);
                            (*Br_[imode])(idl_start+1+i,domain_oversize_r+nsolver/2-j) = (*Br_pml_lmax)(domain_oversize_l+nsolver/2+1+i,n_p[1]-j);
                            (*Ht_[imode])(idl_start+1+i,domain_oversize_r+nsolver/2-j) = (*Ht_pml_lmax)(domain_oversize_l+nsolver/2+1+i,n_p[1]-j);
                            (*Bt_[imode])(idl_start+1+i,domain_oversize_r+nsolver/2-j) = (*Bt_pml_lmax)(domain_oversize_l+nsolver/2+1+i,n_p[1]-j);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for B-field :
        //pml_solver_->compute_E_from_D( EMfields, iDim, min_or_max, solvermin, solvermax);
        pml_solver_->compute_H_from_B( EMfields, iDim, min_or_max, solvermin, solvermax);

        for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
            cField2D *El_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            cField2D *Er_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            cField2D *Et_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            cField2D *Bl_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_[imode];
            cField2D *Br_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_[imode];
            cField2D *Bt_domain = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_[imode];

            if (patch->isXmin()) {
                if(ncells_pml_lmin != 0){
                    El_pml_lmin = pml_fields_lmin->El_[imode];
                    Er_pml_lmin = pml_fields_lmin->Er_[imode];
                    Et_pml_lmin = pml_fields_lmin->Et_[imode];
                    Hl_pml_lmin = pml_fields_lmin->Hl_[imode];
                    Hr_pml_lmin = pml_fields_lmin->Hr_[imode];
                    Ht_pml_lmin = pml_fields_lmin->Ht_[imode];
                    Dl_pml_lmin = pml_fields_lmin->Dl_[imode];
                    Dr_pml_lmin = pml_fields_lmin->Dr_[imode];
                    Dt_pml_lmin = pml_fields_lmin->Dt_[imode];
                    Bl_pml_lmin = pml_fields_lmin->Bl_[imode];
                    Br_pml_lmin = pml_fields_lmin->Br_[imode];
                    Bt_pml_lmin = pml_fields_lmin->Bt_[imode];
                }
            }

            if (patch->isXmax()) {
                if(ncells_pml_lmax != 0){
                    El_pml_lmax = pml_fields_lmax->El_[imode];
                    Er_pml_lmax = pml_fields_lmax->Er_[imode];
                    Et_pml_lmax = pml_fields_lmax->Et_[imode];
                    Hl_pml_lmax = pml_fields_lmax->Hl_[imode];
                    Hr_pml_lmax = pml_fields_lmax->Hr_[imode];
                    Ht_pml_lmax = pml_fields_lmax->Ht_[imode];
                    Dl_pml_lmax = pml_fields_lmax->Dl_[imode];
                    Dr_pml_lmax = pml_fields_lmax->Dr_[imode];
                    Dt_pml_lmax = pml_fields_lmax->Dt_[imode];
                    Bl_pml_lmax = pml_fields_lmax->Bl_[imode];
                    Br_pml_lmax = pml_fields_lmax->Br_[imode];
                    Bt_pml_lmax = pml_fields_lmax->Bt_[imode];
                }
            }

            // 4. Exchange PML -> Domain
            // Primals in y-direction
            for (int j=0 ; j < nsolver/2-1 ; j++){
                for ( int i=0 ; i<n_p[0] ; i++ ) {
                    int idl_start = ncells_pml_lmin;
                    (*Et_domain)(i,n_p[1]-1-j) = (*Et_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                }
                for ( int i=0 ; i<n_d[0] ; i++ ) {
                    int idl_start = ncells_pml_lmin;
                    (*El_domain)(i,n_p[1]-1-j) = (*El_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                    (*Br_domain)(i,n_p[1]-1-j) = (*Hr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                }
            }
            // Duals in y-direction
            for (int j=0 ; j < nsolver/2 ; j++){
                for ( int i=0 ; i<n_p[0] ; i++ ) {
                    int idl_start = ncells_pml_lmin;
                    (*Er_domain)(i,n_d[1]-1-j) = (*Er_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                    (*Bl_domain)(i,n_d[1]-1-j) = (*Hl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                }
                for ( int i=0 ; i<n_d[0] ; i++ ) {
                    int idl_start = ncells_pml_lmin;
                    (*Bt_domain)(i,n_d[1]-1-j) = (*Ht_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                }
            }

            // 5. Exchange PML y -> PML x MIN
            // Primal in y-direction
            for (int j=0 ; j < nsolver/2-1 ; j++){
                if (patch->isXmin()) {
                    if(ncells_pml_lmin != 0){
                        for ( int i=0 ; i<ncells_pml_domain_lmin ; i++ ) {
                            int idl_start = 0;
                            // Primals
                            (*Et_pml_lmin)(i,n_p[1]-1-j) = (*Et_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Dt_pml_lmin)(i,n_p[1]-1-j) = (*Dt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            // Duals
                            (*El_pml_lmin)(i,n_p[1]-1-j) = (*El_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Dl_pml_lmin)(i,n_p[1]-1-j) = (*Dl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Hr_pml_lmin)(i,n_p[1]-1-j) = (*Hr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Br_pml_lmin)(i,n_p[1]-1-j) = (*Br_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                        }
                    }
                }
            }
            // Duals in y-direction
            for (int j=0 ; j < nsolver/2 ; j++){
                if (patch->isXmin()) {
                    if(ncells_pml_lmin != 0){
                        for ( int i=0 ; i<ncells_pml_domain_lmin ; i++ ) {
                            int idl_start = 0;
                            // Primals
                            (*Er_pml_lmin)(i,n_d[1]-1-j) = (*Er_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Dr_pml_lmin)(i,n_d[1]-1-j) = (*Dr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Hl_pml_lmin)(i,n_d[1]-1-j) = (*Hl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Bl_pml_lmin)(i,n_d[1]-1-j) = (*Bl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            // Duals
                            (*Ht_pml_lmin)(i,n_d[1]-1-j) = (*Ht_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Bt_pml_lmin)(i,n_d[1]-1-j) = (*Bt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                        }
                    }
                }
            }

            // 6. Exchange PML y -> PML x MAX
            // Primal in y-direction
            for (int j=0 ; j < nsolver/2-1 ; j++){
                if (patch->isXmax()) {
                    if(ncells_pml_lmax != 0){
                        for ( int i=0 ; i<ncells_pml_domain_lmax ; i++ ) {
                            int idl_start = rpml_size_in_l-ncells_pml_domain_lmax ;
                            // Primals
                            (*Et_pml_lmax)(i,n_p[1]-1-j) = (*Et_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Dt_pml_lmax)(i,n_p[1]-1-j) = (*Dt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            // Dual
                            (*El_pml_lmax)(i,n_p[1]-1-j) = (*El_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Dl_pml_lmax)(i,n_p[1]-1-j) = (*Dl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Hr_pml_lmax)(i,n_p[1]-1-j) = (*Hr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                            (*Br_pml_lmax)(i,n_p[1]-1-j) = (*Br_[imode])(idl_start+i,domain_oversize_r+nsolver/2-1-j);
                        }
                    }
                }
            }
            // Dual in y-direction
            for (int j=0 ; j < nsolver/2 ; j++){
                if (patch->isXmax()) {
                    if(ncells_pml_lmax != 0){
                        for ( int i=0 ; i<ncells_pml_domain_lmax ; i++ ) {
                            int idl_start = rpml_size_in_l-ncells_pml_domain_lmax ;
                            // Primals
                            (*Er_pml_lmax)(i,n_d[1]-1-j) = (*Er_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Dr_pml_lmax)(i,n_d[1]-1-j) = (*Dr_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Hl_pml_lmax)(i,n_d[1]-1-j) = (*Hl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Bl_pml_lmax)(i,n_d[1]-1-j) = (*Bl_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            // Dual
                            (*Ht_pml_lmax)(i,n_d[1]-1-j) = (*Ht_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                            (*Bt_pml_lmax)(i,n_d[1]-1-j) = (*Bt_[imode])(idl_start+i,domain_oversize_r+nsolver/2-j);
                        }
                    }
                }
            }
        }
    }
}
