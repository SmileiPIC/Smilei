#include "EnvelopeBC3D_PML.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "LaserEnvelope.h"
#include "Field3D.h"
#include "cField3D.h"
#include "Tools.h"

#include "SolverFactory.h"

using namespace std;

EnvelopeBC3D_PML::EnvelopeBC3D_PML( Params &params, Patch *patch, unsigned int i_boundary )
    : EnvelopeBC( params, patch, i_boundary )
{

    std::vector<unsigned int> n_space(params.n_space);
    std::vector<unsigned int> oversize(params.oversize);
    if (params.multiple_decomposition) {
        n_space = params.n_space_region;
        oversize = params.region_oversize;
    }

    pml_solver_envelope_ = SolverFactory::createPMLenvelope( params );
    if (params.envelope_solver == "explicit"){
        nsolver=2;
    }
    else if (params.envelope_solver == "explicit_reduced_dispersion"){
        nsolver=2;
    }
    else {
        WARNING("The solver you use in the main domain for envelope is not the same as in the PML region.");
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

        if (patch->isXmin() ) {//&& min_max == 0 ) {
            ncells_pml_xmin = params.number_of_pml_cells[0][0];
            ncells_pml_domain_xmin = ncells_pml_xmin + 1*oversize[0] + nsolver/2;
            domain_oversize_x = oversize[0] ;
        }
        else {
            ncells_pml_xmin = 0;
            ncells_pml_domain_xmin = 0;
        }
        if (patch->isXmax() ) {//&& min_max == 1 ) {
            ncells_pml_xmax = params.number_of_pml_cells[0][1];
            ncells_pml_domain_xmax = ncells_pml_xmax + 1*oversize[0] + nsolver/2;
            domain_oversize_x = oversize[0] ;
        }
        else {
            ncells_pml_xmax = 0;
            ncells_pml_domain_xmax = 0;
        }
        if (patch->isYmin() ) {//&& min_max == 2 ) {
            ncells_pml_ymin = params.number_of_pml_cells[1][0];
            ncells_pml_domain_ymin = ncells_pml_ymin + 1*oversize[1] + nsolver/2;
            domain_oversize_y = oversize[1] ;
        }
        else {
            ncells_pml_ymin = 0;
            ncells_pml_domain_ymin = 0;
        }
        if (patch->isYmax() ) {//&& min_max == 3 ) {
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
            solvermax = ncells_pml_domain - 1 - oversize[iDim] ;
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
            dimPrim[iDim-1] += (ncells_pml_xmin-1*(patch->isXmin())) + (ncells_pml_xmax-1*(patch->isXmax())) ;
            ypml_size_in_x = dimPrim[iDim-1] ;
            //std::cout << "size : " << ypml_size_in_x << ".";
        }
        if ( iDim==2 ){
            dimPrim[0] += (ncells_pml_xmin-1*(patch->isXmin())) + (ncells_pml_xmax-1*(patch->isXmax())) ;
            dimPrim[1] += (ncells_pml_ymin-1*(patch->isYmin())) + (ncells_pml_ymax-1*(patch->isYmax())) ;
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

        pml_solver_envelope_->setDomainSizeAndCoefficients( iDim, min_or_max, ncells_pml_domain, startpml, ncells_pml_min, ncells_pml_max, patch );

        // A-field
        A3D_np1_ = new cField3D( dimPrim, "A_np1_pml" );
        A3D_n_ = new cField3D( dimPrim, "A_n_pml" );
        A3D_nm1_ = new cField3D( dimPrim, "A_nm1_pml" );
        // Auxillary Variable
        u1_np1_x_ = new cField3D( dimPrim, "u1_np1_x_pml" );
        u2_np1_x_ = new cField3D( dimPrim, "u2_np1_x_pml" );
        u3_np1_x_ = new cField3D( dimPrim, "u3_np1_x_pml" );
        u1_nm1_x_ = new cField3D( dimPrim, "u1_nm1_x_pml" );
        u2_nm1_x_ = new cField3D( dimPrim, "u2_nm1_x_pml" );
        u3_nm1_x_ = new cField3D( dimPrim, "u3_nm1_x_pml" );
        // ----
        u1_np1_y_ = new cField3D( dimPrim, "u1_np1_y_pml" );
        u2_np1_y_ = new cField3D( dimPrim, "u2_np1_y_pml" );
        u3_np1_y_ = new cField3D( dimPrim, "u3_np1_y_pml" );
        u1_nm1_y_ = new cField3D( dimPrim, "u1_nm1_y_pml" );
        u2_nm1_y_ = new cField3D( dimPrim, "u2_nm1_y_pml" );
        u3_nm1_y_ = new cField3D( dimPrim, "u3_nm1_y_pml" );
        // ----
        u1_np1_z_ = new cField3D( dimPrim, "u1_np1_z_pml" );
        u2_np1_z_ = new cField3D( dimPrim, "u2_np1_z_pml" );
        u3_np1_z_ = new cField3D( dimPrim, "u3_np1_z_pml" );
        u1_nm1_z_ = new cField3D( dimPrim, "u1_nm1_z_pml" );
        u2_nm1_z_ = new cField3D( dimPrim, "u2_nm1_z_pml" );
        u3_nm1_z_ = new cField3D( dimPrim, "u3_nm1_z_pml" );

        // Ponderomoteur Potential
        Phi3D_ = new Field3D( dimPrim, "Phi_pml" );
    }
}


EnvelopeBC3D_PML::~EnvelopeBC3D_PML()
{
    delete A3D_np1_;
    delete A3D_n_;
    delete A3D_nm1_;
    // ----
    delete u1_np1_x_;
    delete u2_np1_x_;
    delete u3_np1_x_;
    delete u1_nm1_x_;
    delete u2_nm1_x_;
    delete u3_nm1_x_;
    // ----
    delete u1_np1_y_;
    delete u2_np1_y_;
    delete u3_np1_y_;
    delete u1_nm1_y_;
    delete u2_nm1_y_;
    delete u3_nm1_y_;
    // ----
    delete u1_np1_z_;
    delete u2_np1_z_;
    delete u3_np1_z_;
    delete u1_nm1_z_;
    delete u2_nm1_z_;
    delete u3_nm1_z_;
    // ----
    delete Phi3D_;

    if (pml_solver_envelope_!=NULL) {
        delete pml_solver_envelope_;
    }
}

/*
void EnvelopeBC3D_PML::save_fields( Field *my_field, Patch *patch )
{
}


void EnvelopeBC3D_PML::disableExternalFields()
{
}
*/

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------

void EnvelopeBC3D_PML::apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    int iDim = 0*((i_boundary_==0)||(i_boundary_==1))+1*((i_boundary_==2)||(i_boundary_==3))+2*((i_boundary_==4)||(i_boundary_==5));
    int min_or_max = (i_boundary_)%2;

    cField3D *A3D_np1_domain  = static_cast<cField3D *>( envelope->A_ );  // A_ is the envelope at timestep n BUT at this point A_ is already update so in fact its correspond to A3D_np1_domain
    cField3D *A3D_n_domain    = static_cast<cField3D *>( envelope->A0_ ); // A0_ is the envelope at timestep n-1 BUT at this point A0_ is already update so in fact its A3D_n_domain
    Field3D  *Phi3D_domain    = static_cast<Field3D *>( envelope->Phi_ ); // the ponderomotive potential Phi=|A|^2/2 at timestep n

    double ellipticity_factor = envelope->ellipticity_factor;

    if( i_boundary_ == 0 && patch->isXmin() ) {

        // 2. Exchange field PML <- Domain
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            for ( int j=1 ; j<ny_p-1 ; j++ ) {
                for ( int k=1 ; k<nz_p-1 ; k++ ) {
                    (*A3D_n_)(ncells_pml_domain-1-domain_oversize_x-nsolver/2+i,j,k) = (*A3D_n_domain)(i,j,k);
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( int j=0 ; j<ny_p ; j++ ) {
                for ( int k=0 ; k<nz_p ; k++ ) {
                    (*A3D_np1_domain)(i,j,k) = (*A3D_n_)(ncells_pml_domain-1-domain_oversize_x-nsolver/2+i,j,k);
                    (*A3D_n_domain)(i,j,k) = (*A3D_nm1_)(ncells_pml_domain-1-domain_oversize_x-nsolver/2+i,j,k);
                    (*Phi3D_domain)(i,j,k) = 0.5*ellipticity_factor*std::abs(  (*A3D_np1_domain)(i,j,k) )*std::abs(  (*A3D_np1_domain)(i,j,k) );
                }
            }
        }
    }

    else if( i_boundary_ == 1 && patch->isXmax() ) {

        // 2. Exchange field Domain -> PML
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            for ( int j=1 ; j<ny_p-1 ; j++ ) {
                for ( int k=1 ; k<nz_p-1 ; k++ ) {
                    (*A3D_n_)(domain_oversize_x+nsolver/2-i,j,k) = (*A3D_n_domain)(nx_p-1-i,j,k);
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, solvermin, solvermax);

        // 4. Exchange Domain -> PML
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( int j=0 ; j<ny_p ; j++ ) {
                for ( int k=0 ; k<nz_p ; k++ ) {
                    (*A3D_np1_domain)(nx_p-1-i,j,k) = (*A3D_n_)(domain_oversize_x+nsolver/2-i,j,k);
                    (*A3D_n_domain)(nx_p-1-i,j,k) = (*A3D_nm1_)(domain_oversize_x+nsolver/2-i,j,k);
                    (*Phi3D_domain)(nx_p-1-i,j,k) = 0.5*ellipticity_factor*std::abs( (*A3D_np1_domain)(nx_p-1-i,j,k) )*std::abs( (*A3D_np1_domain)(nx_p-1-i,j,k) ) ;
                }
            }
        }
    }

    else if( i_boundary_ == 2 && patch->isYmin() ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );

        cField3D* A3D_np1_pml_xmin = NULL;
        cField3D* A3D_n_pml_xmin   = NULL;
        Field3D*  Phi3D_pml_xmin   = NULL;

        cField3D* A3D_np1_pml_xmax = NULL;
        cField3D* A3D_n_pml_xmax   = NULL;
        Field3D*  Phi3D_pml_xmax   = NULL;

        if(ncells_pml_xmin != 0){
            A3D_np1_pml_xmin = pml_fields_xmin->A3D_n_;
            A3D_n_pml_xmin   = pml_fields_xmin->A3D_nm1_;
            Phi3D_pml_xmin   = pml_fields_xmin->Phi3D_;
        }

        if(ncells_pml_xmax != 0){
            A3D_np1_pml_xmax = pml_fields_xmax->A3D_n_;
            A3D_n_pml_xmax   = pml_fields_xmax->A3D_nm1_;
            Phi3D_pml_xmax   = pml_fields_xmax->Phi3D_;
        }

        // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_n_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k) = (*A3D_n_pml_xmin)(i,j,k);
                        }
                    }
                }
            }
            for ( int i=1 ; i<nx_p-1 ; i++ ) {
                for ( int k=1 ; k<nz_p-1 ; k++ ) {
                    int idx_start = ncells_pml_xmin-1*(patch->isXmin());
                    (*A3D_n_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k) = (*A3D_n_domain)(i,j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_n_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k) = (*A3D_n_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, solvermin, solvermax);

        // std::cout << "nx_p : " << nx_p << "."; // 1029

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for ( int i=1 ; i<nx_p-1 ; i++ ) { // From i=0 to i=1028 with ypml_size_in_x=1048
                int idx_start = ncells_pml_xmin-1*(patch->isXmin());
                for ( int k=1 ; k<nz_p-1 ; k++ ) {
                    (*A3D_np1_domain)(i,j,k) = (*A3D_n_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k);
                    (*A3D_n_domain)(i,j,k) = (*A3D_nm1_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k);
                    (*Phi3D_domain)(i,j,k) = 0.5*ellipticity_factor*std::abs( (*A3D_np1_domain)(i,j,k) )*std::abs( (*A3D_np1_domain)(i,j,k) );
                }
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        // Dual in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_np1_pml_xmin)(i,j,k) = (*A3D_n_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k);
                            (*A3D_n_pml_xmin)(i,j,k) = (*A3D_nm1_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k);
                        }
                    }
                }
            }
        }

        // 6. Exchange PML y -> PML x MAX
        // Primals in y-direction
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_domain_xmax ;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_np1_pml_xmax  )(i,j,k) = (*A3D_n_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k);
                            (*A3D_n_pml_xmax  )(i,j,k) = (*A3D_nm1_)(idx_start+i,ncells_pml_domain-1-domain_oversize_y-nsolver/2+j,k);
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 3 && patch->isYmax() ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );

        cField3D* A3D_np1_pml_xmin = NULL;
        cField3D* A3D_n_pml_xmin   = NULL;
        Field3D*  Phi3D_pml_xmin   = NULL;

        cField3D* A3D_np1_pml_xmax = NULL;
        cField3D* A3D_n_pml_xmax   = NULL;
        Field3D*  Phi3D_pml_xmax   = NULL;

        if(ncells_pml_xmin != 0){
            A3D_np1_pml_xmin = pml_fields_xmin->A3D_n_;
            A3D_n_pml_xmin   = pml_fields_xmin->A3D_nm1_;
            Phi3D_pml_xmin   = pml_fields_xmin->Phi3D_;
        }

        if(ncells_pml_xmax != 0){
            A3D_np1_pml_xmax = pml_fields_xmax->A3D_n_;
            A3D_n_pml_xmax   = pml_fields_xmax->A3D_nm1_;
            Phi3D_pml_xmax   = pml_fields_xmax->Phi3D_;
        }

        // Attention à la gestion des pas de temps
        // A priori pour calculer A_np1 on a juste besoin de A_n
        // Même les variables intermediaires sont inutiles car calculé à partir de An pour (i,j)

        // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*A3D_n_pml_xmin)(i,ny_p-1-j,k);
                        }
                    }
                }
            }
            for ( int i=1 ; i<nx_p-1 ; i++ ) {
                int idx_start = ncells_pml_xmin-1*(patch->isXmin());
                for ( int k=1 ; k<nz_p-1 ; k++ ) {
                    (*A3D_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*A3D_n_domain)(i,ny_p-1-j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*A3D_n_pml_xmax)(domain_oversize_x+nsolver/2+i,ny_p-1-j,k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for ( int i=1 ; i<nx_p-1 ; i++ ) {
                int idx_start = ncells_pml_xmin-1*(patch->isXmin());
                for ( int k=1 ; k<nz_p-1 ; k++ ) {
                    (*A3D_np1_domain)(i,ny_p-1-j,k) = (*A3D_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                    (*A3D_n_domain)(i,ny_p-1-j,k) = (*A3D_nm1_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                    (*Phi3D_domain)(i,ny_p-1-j,k) = 0.5*ellipticity_factor*std::abs( (*A3D_np1_domain)(i,ny_p-1-j,k) )*std::abs( (*A3D_np1_domain)(i,ny_p-1-j,k) );
                }
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        // Dual in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_np1_pml_xmin)(i,ny_p-1-j,k) = (*A3D_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*A3D_n_pml_xmin)(i,ny_p-1-j,k) = (*A3D_nm1_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                        }
                    }
                }
            }
        }

        // 6. Exchange PML y -> PML x MAX
        // Primals in y-direction
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_domain_xmax ;
                        for ( int k=0 ; k<nz_p ; k++ ) {
                            (*A3D_np1_pml_xmax  )(i,ny_p-1-j,k) = (*A3D_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                            (*A3D_n_pml_xmax  )(i,ny_p-1-j,k) = (*A3D_nm1_)(idx_start+i,domain_oversize_y+nsolver/2-j,k);
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 4 && patch->isZmin() ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );
        EnvelopeBC3D_PML* pml_fields_ymin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[2] );
        EnvelopeBC3D_PML* pml_fields_ymax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[3] );

        cField3D* A3D_np1_pml_xmin = NULL;
        cField3D* A3D_n_pml_xmin   = NULL;
        Field3D*  Phi3D_pml_xmin   = NULL;

        cField3D* A3D_np1_pml_xmax = NULL;
        cField3D* A3D_n_pml_xmax   = NULL;
        Field3D*  Phi3D_pml_xmax   = NULL;

        cField3D* A3D_np1_pml_ymin = NULL;
        cField3D* A3D_n_pml_ymin   = NULL;
        Field3D*  Phi3D_pml_ymin   = NULL;

        cField3D* A3D_np1_pml_ymax = NULL;
        cField3D* A3D_n_pml_ymax   = NULL;
        Field3D*  Phi3D_pml_ymax   = NULL;

        if(ncells_pml_xmin != 0){
            A3D_np1_pml_xmin = pml_fields_xmin->A3D_n_;
            A3D_n_pml_xmin   = pml_fields_xmin->A3D_nm1_;
            Phi3D_pml_xmin   = pml_fields_xmin->Phi3D_;
        }

        if(ncells_pml_xmax != 0){
            A3D_np1_pml_xmax = pml_fields_xmax->A3D_n_;
            A3D_n_pml_xmax   = pml_fields_xmax->A3D_nm1_;
            Phi3D_pml_xmax   = pml_fields_xmax->Phi3D_;
        }

        if(ncells_pml_ymin != 0){
            A3D_np1_pml_ymin = pml_fields_ymin->A3D_n_;
            A3D_n_pml_ymin   = pml_fields_ymin->A3D_nm1_;
            Phi3D_pml_ymin   = pml_fields_ymin->Phi3D_;
        }

        if(ncells_pml_ymax != 0){
            A3D_np1_pml_ymax = pml_fields_ymax->A3D_n_;
            A3D_n_pml_ymax   = pml_fields_ymax->A3D_nm1_;
            Phi3D_pml_ymax   = pml_fields_ymax->Phi3D_;
        }

        // 2. Exchange field PML <- Domain
        for ( int k=min2exchange ; k<max2exchange ; k++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        for ( int j=0 ; j<ny_p ; j++ ) {
                            (*A3D_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k) = (*A3D_n_pml_xmin)(i,j,k);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            (*A3D_n_)(i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k) = (*A3D_n_pml_xmin)(i,j,k);
                        }
                    }
                }
            }
            int idx_start = ncells_pml_xmin-1*(patch->isXmin());
            int jdx_start = ncells_pml_ymin-1*(patch->isYmin());
            for ( int i=1 ; i<nx_p-1 ; i++ ) {
                for ( int j=1 ; k<ny_p-1 ; j++ ) {
                    (*A3D_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k) = (*A3D_n_domain)(i,j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = (zpml_size_in_x-1)-(ncells_pml_xmax-1) ;
                    int jdx_start = ncells_pml_ymin ;
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ny_p ; j++ ) {
                            (*A3D_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k) = (*A3D_n_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = (zpml_size_in_y-1)-(ncells_pml_ymax-1) ;
                    for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            (*A3D_n_)(i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k) = (*A3D_n_pml_xmax)(i,domain_oversize_y+nsolver/2+j,k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            int idx_start = ncells_pml_xmin-1*(patch->isXmin());
            int jdx_start = ncells_pml_ymin-1*(patch->isYmin());
            for ( int i=1 ; i<nx_p-1 ; i++ ) {
                for ( int j=1 ; j<ny_p-1 ; j++ ) {
                    (*A3D_np1_domain)(i,j,k) = (*A3D_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                    (*A3D_n_domain)(i,j,k) = (*A3D_nm1_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                    (*Phi3D_domain)(i,j,k) = 0.5*ellipticity_factor*std::abs( (*A3D_np1_domain)(i,j,k) )*std::abs( (*A3D_np1_domain)(i,j,k) );
                }
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        // Dual in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        for ( int j=0 ; j<ny_p ; k++ ) {
                            (*A3D_np1_pml_xmin)(i,j,k) = (*A3D_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                            (*A3D_n_pml_xmin)(i,j,k) = (*A3D_nm1_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }

        // 6. Exchange PML y -> PML x MAX
        // Primals in y-direction
        // Duals in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = zpml_size_in_x-ncells_pml_domain_xmax ;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        for ( int j=0 ; j<ny_p ; j++ ) {
                            (*A3D_np1_pml_xmax  )(i,j,k) = (*A3D_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                            (*A3D_n_pml_xmax  )(i,j,k) = (*A3D_nm1_)(idx_start+i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }

        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int j=0 ; j<ncells_pml_domain_ymin ; j++ ) {
                        for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            (*A3D_np1_pml_ymin)(i,j,k) = (*A3D_n_)(i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                            (*A3D_n_pml_ymin)(i,j,k) = (*A3D_nm1_)(i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }

        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = zpml_size_in_y-ncells_pml_domain_ymax ;
                    for ( int j=0 ; j<ncells_pml_domain_ymax ; j++ ) {
                        for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            (*A3D_np1_pml_ymax)(i,j,k) = (*A3D_n_)(i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                            (*A3D_n_pml_ymax)(i,j,k) = (*A3D_nm1_)(i,jdx_start+j,ncells_pml_domain-1-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 5 && patch->isZmax() ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );
        EnvelopeBC3D_PML* pml_fields_ymin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[2] );
        EnvelopeBC3D_PML* pml_fields_ymax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[3] );

        cField3D* A3D_np1_pml_xmin = NULL;
        cField3D* A3D_n_pml_xmin   = NULL;
        Field3D*  Phi3D_pml_xmin   = NULL;

        cField3D* A3D_np1_pml_xmax = NULL;
        cField3D* A3D_n_pml_xmax   = NULL;
        Field3D*  Phi3D_pml_xmax   = NULL;

        cField3D* A3D_np1_pml_ymin = NULL;
        cField3D* A3D_n_pml_ymin   = NULL;
        Field3D*  Phi3D_pml_ymin   = NULL;

        cField3D* A3D_np1_pml_ymax = NULL;
        cField3D* A3D_n_pml_ymax   = NULL;
        Field3D*  Phi3D_pml_ymax   = NULL;

        if(ncells_pml_xmin != 0){
            A3D_np1_pml_xmin = pml_fields_xmin->A3D_n_;
            A3D_n_pml_xmin   = pml_fields_xmin->A3D_nm1_;
            Phi3D_pml_xmin   = pml_fields_xmin->Phi3D_;
        }

        if(ncells_pml_xmax != 0){
            A3D_np1_pml_xmax = pml_fields_xmax->A3D_n_;
            A3D_n_pml_xmax   = pml_fields_xmax->A3D_nm1_;
            Phi3D_pml_xmax   = pml_fields_xmax->Phi3D_;
        }

        if(ncells_pml_ymin != 0){
            A3D_np1_pml_ymin = pml_fields_ymin->A3D_n_;
            A3D_n_pml_ymin   = pml_fields_ymin->A3D_nm1_;
            Phi3D_pml_ymin   = pml_fields_ymin->Phi3D_;
        }

        if(ncells_pml_ymax != 0){
            A3D_np1_pml_ymax = pml_fields_ymax->A3D_n_;
            A3D_n_pml_ymax   = pml_fields_ymax->A3D_nm1_;
            Phi3D_pml_ymax   = pml_fields_ymax->Phi3D_;
        }

        // Attention à la gestion des pas de temps
        // A priori pour calculer A_np1 on a juste besoin de A_n
        // Même les variables intermediaires sont inutiles car calculé à partir de An pour (i,j)

        // 2. Exchange field PML <- Domain
        for ( int k=min2exchange ; k<max2exchange ; k++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        for ( int j=0 ; j<ny_p ; j++ ) {
                            (*A3D_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A3D_n_pml_xmin)(i,j,nz_p-1-k);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            (*A3D_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A3D_n_pml_xmin)(i,j,nz_p-1-k);
                        }
                    }
                }
            }
            int idx_start = ncells_pml_xmin-1*(patch->isXmin());
            int jdx_start = ncells_pml_ymin-1*(patch->isYmin());
            for ( int i=1 ; i<nx_p-1 ; i++ ) {
                for ( int j=1 ; k<ny_p-1 ; j++ ) {
                    (*A3D_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A3D_n_domain)(i,j,nz_p-1-k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = (zpml_size_in_x-1)-(ncells_pml_xmax-1) ;
                    int jdx_start = ncells_pml_ymin ;
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ny_p ; j++ ) {
                            (*A3D_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A3D_n_pml_xmax)(domain_oversize_x+nsolver/2+i,j,nz_p-1-k);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = (zpml_size_in_y-1)-(ncells_pml_ymax-1) ;
                    for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            (*A3D_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A3D_n_pml_xmax)(i,domain_oversize_y+nsolver/2+j,nz_p-1-k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            int idx_start = ncells_pml_xmin-1*(patch->isXmin());
            int jdx_start = ncells_pml_ymin-1*(patch->isYmin());
            for ( int i=1 ; i<nx_p-1 ; i++ ) {
                for ( int j=1 ; j<ny_p-1 ; j++ ) {
                    (*A3D_np1_domain)(i,j,nz_p-1-k) = (*A3D_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                    (*A3D_n_domain)(i,j,nz_p-1-k) = (*A3D_nm1_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                    (*Phi3D_domain)(i,j,nz_p-1-k) = 0.5*ellipticity_factor*std::abs( (*A3D_np1_domain)(i,j,nz_p-1-k) )*std::abs( (*A3D_np1_domain)(i,j,nz_p-1-k) );
                }
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        // Dual in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        for ( int j=0 ; j<ny_p ; k++ ) {
                            (*A3D_np1_pml_xmin)(i,j,nz_p-1-k) = (*A3D_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*A3D_n_pml_xmin)(i,j,nz_p-1-k) = (*A3D_nm1_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                    }
                }
            }
        }

        // 6. Exchange PML y -> PML x MAX
        // Primals in y-direction
        // Duals in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = zpml_size_in_x-ncells_pml_domain_xmax ;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_domain_xmax ; i++ ) {
                        for ( int j=0 ; j<ny_p ; j++ ) {
                            (*A3D_np1_pml_xmax  )(i,j,ny_p-1-k) = (*A3D_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*A3D_n_pml_xmax  )(i,j,ny_p-1-k) = (*A3D_nm1_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                    }
                }
            }
        }

        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( int j=0 ; j<ncells_pml_domain_ymin ; j++ ) {
                        for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            (*A3D_np1_pml_ymin)(i,j,nz_p-1-k) = (*A3D_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*A3D_n_pml_ymin)(i,j,nz_p-1-k) = (*A3D_nm1_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                    }
                }
            }
        }

        for (int k=0 ; k < nsolver/2 ; k++){
            if (patch->isYmax()) {
                int jdx_start = zpml_size_in_y-ncells_pml_domain_ymax ;
                if(ncells_pml_ymax != 0){
                    for ( int j=0 ; j<ncells_pml_domain_ymax ; j++ ) {
                        for ( int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            (*A3D_np1_pml_ymax)(i,j,nz_p-1-k) = (*A3D_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                            (*A3D_n_pml_ymax)(i,j,nz_p-1-k) = (*A3D_nm1_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k);
                        }
                    }
                }
            }
        }
    }
}
