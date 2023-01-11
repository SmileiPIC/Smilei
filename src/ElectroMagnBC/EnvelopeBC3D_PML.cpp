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
    pml_solver_envelope_ = SolverFactory::createPMLenvelope( params );
    if (params.envelope_solver == "explicit"){
        nsolver=4;
    }
    else if (params.envelope_solver == "explicit_reduced_dispersion"){
        nsolver=4;
    }
    else {
        WARNING("The solver you use in the main domain for envelope is not the same as in the PML region.");
        nsolver=4;
    }

    if( patch->isBoundary( i_boundary_ ) ) {

        int iDim = i_boundary_ / 2;
        int min_or_max = (i_boundary_)%2;

        domain_oversize_x =  patch->oversize[0] ;
        domain_oversize_y =  patch->oversize[1] ;
        domain_oversize_z =  patch->oversize[2] ;

        if (patch->isXmin() ) {//&& min_max == 0 ) {
            ncells_pml_xmin = params.number_of_pml_cells[0][0];
            ncells_pml_domain_xmin = ncells_pml_xmin + 1*patch->oversize[0] + nsolver/2;
            domain_oversize_x = patch->oversize[0] ;
        }
        else {
            ncells_pml_xmin = 0;
            ncells_pml_domain_xmin = 0;
        }
        if (patch->isXmax() ) {//&& min_max == 1 ) {
            ncells_pml_xmax = params.number_of_pml_cells[0][1];
            ncells_pml_domain_xmax = ncells_pml_xmax + 1*patch->oversize[0] + nsolver/2;
            domain_oversize_x = patch->oversize[0] ;
        }
        else {
            ncells_pml_xmax = 0;
            ncells_pml_domain_xmax = 0;
        }
        if (patch->isYmin() ) {//&& min_max == 2 ) {
            ncells_pml_ymin = params.number_of_pml_cells[1][0];
            ncells_pml_domain_ymin = ncells_pml_ymin + 1*patch->oversize[1] + nsolver/2;
            domain_oversize_y = patch->oversize[1] ;
        }
        else {
            ncells_pml_ymin = 0;
            ncells_pml_domain_ymin = 0;
        }
        if (patch->isYmax() ) {//&& min_max == 3 ) {
            ncells_pml_ymax = params.number_of_pml_cells[1][1];
            ncells_pml_domain_ymax = ncells_pml_ymax + 1*patch->oversize[1] + nsolver/2;
            domain_oversize_y = patch->oversize[1] ;
        }
        else {
            ncells_pml_ymax = 0;
            ncells_pml_domain_ymax = 0;
        }
        if (patch->isZmin() ) {//&& i_boundary_ == 2 ) {
            ncells_pml_zmin = params.number_of_pml_cells[2][0];
            ncells_pml_domain_zmin = ncells_pml_zmin + 1*patch->oversize[2] + nsolver/2;
            domain_oversize_z = patch->oversize[2] ;
        }
        else {
            ncells_pml_zmin = 0;
            ncells_pml_domain_zmin = 0;
        }
        if (patch->isZmax() ) {//&& i_boundary_ == 3 ) {
            ncells_pml_zmax = params.number_of_pml_cells[2][1];
            ncells_pml_domain_zmax = ncells_pml_zmax + 1*patch->oversize[2] + nsolver/2;
            domain_oversize_z = patch->oversize[2] ;
        }
        else {
            ncells_pml_zmax = 0;
            ncells_pml_domain_zmax = 0;
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
            min2exchange = 1*nsolver/2+1 ;
            max2exchange = 2*nsolver/2+1 ;
            // Solver
            solvermin = patch->oversize[iDim]  ;
            solvermax = ncells_pml_domain-nsolver/2 ;
        }

        if (ncells_pml==0){
            ERROR("PML domain have to be >0 cells in thickness");
        }

        if( iDim == 0 ) {
            dimPrim = { (unsigned int)ncells_pml_domain, patch->size_[1]+1+2*patch->oversize[1], patch->size_[2]+1+2*patch->oversize[2] };
        } else if( iDim == 1 ) {
            ypml_size_in_x = patch->size_[0]+1+2*patch->oversize[0] + ncells_pml_xmin + ncells_pml_xmax;
            dimPrim = { (unsigned int)ypml_size_in_x, (unsigned int)ncells_pml_domain, patch->size_[2]+1+2*patch->oversize[2] };
        } else {
            zpml_size_in_x = patch->size_[0]+1+2*patch->oversize[0] + ncells_pml_xmin + ncells_pml_xmax;
            zpml_size_in_y = patch->size_[1]+1+2*patch->oversize[1] + ncells_pml_ymin + ncells_pml_ymax;
            dimPrim = { (unsigned int)ypml_size_in_x, (unsigned int)zpml_size_in_y, (unsigned int)ncells_pml_domain };
        }
        
        startpml = patch->oversize[iDim]+nsolver/2;

        int ncells_pml_min[2];
        ncells_pml_min[0] = ncells_pml_xmin;
        ncells_pml_min[1] = ncells_pml_ymin;
        int ncells_pml_max[2];
        ncells_pml_max[0] = ncells_pml_xmax;
        ncells_pml_max[1] = ncells_pml_ymax;

        pml_solver_envelope_->setDomainSizeAndCoefficients( iDim, min_or_max, dimPrim, ncells_pml_domain, startpml, ncells_pml_min, ncells_pml_max, patch );

        std::string si_boundary = std::to_string(i_boundary_);
        // A-field
        A_np1_ = new cField3D( dimPrim, "A_np1_pml"+si_boundary );
        A_n_ = new cField3D( dimPrim, "A_n_pml"+si_boundary );
        A_nm1_ = new cField3D( dimPrim, "A_nm1_pml"+si_boundary );
        // Auxillary Variable
        u1_np1_x_ = new cField3D( dimPrim, "u1_np1_x_pml"+si_boundary );
        u2_np1_x_ = new cField3D( dimPrim, "u2_np1_x_pml"+si_boundary );
        u3_np1_x_ = new cField3D( dimPrim, "u3_np1_x_pml"+si_boundary );
        u1_nm1_x_ = new cField3D( dimPrim, "u1_nm1_x_pml"+si_boundary );
        u2_nm1_x_ = new cField3D( dimPrim, "u2_nm1_x_pml"+si_boundary );
        u3_nm1_x_ = new cField3D( dimPrim, "u3_nm1_x_pml"+si_boundary );
        // ----
        u1_np1_y_ = new cField3D( dimPrim, "u1_np1_y_pml"+si_boundary );
        u2_np1_y_ = new cField3D( dimPrim, "u2_np1_y_pml"+si_boundary );
        u3_np1_y_ = new cField3D( dimPrim, "u3_np1_y_pml"+si_boundary );
        u1_nm1_y_ = new cField3D( dimPrim, "u1_nm1_y_pml"+si_boundary );
        u2_nm1_y_ = new cField3D( dimPrim, "u2_nm1_y_pml"+si_boundary );
        u3_nm1_y_ = new cField3D( dimPrim, "u3_nm1_y_pml"+si_boundary );
        // ----
        u1_np1_z_ = new cField3D( dimPrim, "u1_np1_z_pml"+si_boundary );
        u2_np1_z_ = new cField3D( dimPrim, "u2_np1_z_pml"+si_boundary );
        u3_np1_z_ = new cField3D( dimPrim, "u3_np1_z_pml"+si_boundary );
        u1_nm1_z_ = new cField3D( dimPrim, "u1_nm1_z_pml"+si_boundary );
        u2_nm1_z_ = new cField3D( dimPrim, "u2_nm1_z_pml"+si_boundary );
        u3_nm1_z_ = new cField3D( dimPrim, "u3_nm1_z_pml"+si_boundary );

        // Ponderomoteur Potential
        Phi_ = new Field3D( dimPrim, "Phi_pml"+si_boundary );
        Chi_ = new Field3D( dimPrim, "Chi_pml"+si_boundary );
    }
}


EnvelopeBC3D_PML::~EnvelopeBC3D_PML()
{
    delete A_np1_;
    delete A_n_;
    delete A_nm1_;
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
    delete Phi_;
    delete Chi_;

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

void EnvelopeBC3D_PML::apply( LaserEnvelope *envelope, ElectroMagn *EMfields/*time_dual*/, Patch *patch )
{
    int iDim = i_boundary_ / 2;
    int min_or_max = (i_boundary_)%2;

    cField3D *A_np1_domain  = static_cast<cField3D *>( envelope->A_ );  // A_ is the envelope at timestep n BUT at this point A_ is already update so in fact its correspond to A_np1_domain
    cField3D *A_n_domain    = static_cast<cField3D *>( envelope->A0_ ); // A0_ is the envelope at timestep n-1 BUT at this point A0_ is already update so in fact its A_n_domain
    Field3D  *Phi_domain    = static_cast<Field3D *>( envelope->Phi_ ); // the ponderomotive potential Phi=|A|^2/2 at timestep n
    Field3D  *Chi_domain    = static_cast<Field3D *>( EMfields->Env_Chi_ );

    double ellipticity_factor = envelope->ellipticity_factor;

    if( ! patch->isBoundary( i_boundary_ ) ) return;

    if( i_boundary_ == 0 ) {

        // 2. Exchange field PML <- Domain
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*A_n_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k) = (*A_n_domain)(i,j,k);
                }
            }
        }

        // Here we extrude Chi in PML domain according Chi at the xmin border
        for ( int i=0 ; i<ncells_pml_domain ; i++ ) {
            for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*Chi_)(i,j,k) = (*Chi_domain)(domain_oversize_x+1,j,k);
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*A_np1_domain)(i,j,k) = (*A_n_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                    (*A_n_domain)(i,j,k) = (*A_nm1_)(ncells_pml_domain-domain_oversize_x-nsolver/2+i,j,k);
                    (*Phi_domain)(i,j,k) = 0.5*ellipticity_factor*std::abs(  (*A_np1_domain)(i,j,k) )*std::abs(  (*A_np1_domain)(i,j,k) );
                }
            }
        }
    }

    else if( i_boundary_ == 1 ) {

        // 2. Exchange field Domain -> PML
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*A_n_)(domain_oversize_x+nsolver/2-i,j,k) = (*A_n_domain)(nx_p-i,j,k);
                }
            }
        }

        // Here we extrude Chi in PML domain according Chi at the xmax border
        for ( int i=0 ; i<ncells_pml_domain ; i++ ) {
            for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*Chi_)(i,j,k) = (*Chi_domain)(nx_p-2-domain_oversize_x,j,k);
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange Domain -> PML
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*A_np1_domain)(nx_p-1-i,j,k) = (*A_n_)(domain_oversize_x+nsolver/2-1-i,j,k);
                    (*A_n_domain)(nx_p-1-i,j,k) = (*A_nm1_)(domain_oversize_x+nsolver/2-1-i,j,k);
                    (*Phi_domain)(nx_p-1-i,j,k) = 0.5*ellipticity_factor*std::abs( (*A_np1_domain)(nx_p-1-i,j,k) )*std::abs( (*A_np1_domain)(nx_p-1-i,j,k) ) ;
                }
            }
        }
    }

    else if( i_boundary_ == 2 ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );

        cField3D* A_np1_pml_xmin = NULL;
        cField3D* A_n_pml_xmin   = NULL;
        Field3D*  Chi_pml_xmin   = NULL;

        cField3D* A_np1_pml_xmax = NULL;
        cField3D* A_n_pml_xmax   = NULL;
        Field3D*  Chi_pml_xmax   = NULL;

        if(ncells_pml_xmin != 0){
            A_np1_pml_xmin = pml_fields_xmin->A_n_;
            A_n_pml_xmin   = pml_fields_xmin->A_nm1_;
            Chi_pml_xmin   = pml_fields_xmin->Chi_;
        }

        if(ncells_pml_xmax != 0){
            A_np1_pml_xmax = pml_fields_xmax->A_n_;
            A_n_pml_xmax   = pml_fields_xmax->A_nm1_;
            Chi_pml_xmax   = pml_fields_xmax->Chi_;
        }

        // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_n_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*A_n_pml_xmin)(i,j,k);
                        }
                    }
                }
            }
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    int idx_start = ncells_pml_xmin;
                    (*A_n_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*A_n_domain)(i,j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_n_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k) = (*A_n_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                        }
                    }
                }
            }
        }

        // Here we extrude Chi in PML domain according Chi at the ymin border
        for ( int j=0 ; j<ncells_pml_domain ; j++ ) {
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*Chi_)(idx_start+i,j,k) = (*Chi_domain)(i,domain_oversize_y+1,k);
                }
            }
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin+domain_oversize_x ; i++ ) {
                        int idx_start = 0 ;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*Chi_)(idx_start+i,j,k) = (*Chi_pml_xmin)(i,domain_oversize_y+1,k);
                        }
                    }
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax+domain_oversize_x ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_xmax-domain_oversize_x ;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*Chi_)(idx_start+i,j,k) = (*Chi_pml_xmax)(i,domain_oversize_y+1,k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // std::cout << "nx_p : " << nx_p << "."; // 1029

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for ( unsigned int i=0 ; i<nx_p ; i++ ) { // From i=0 to i=1028 with ypml_size_in_x=1048
                int idx_start = ncells_pml_xmin;
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*A_np1_domain)(i,j,k) = (*A_n_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k);
                    (*A_n_domain)(i,j,k) = (*A_nm1_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k);
                    (*Phi_domain)(i,j,k) = 0.5*ellipticity_factor*std::abs( (*A_np1_domain)(i,j,k) )*std::abs( (*A_np1_domain)(i,j,k) );
                }
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_xmin ; i++ ) {
                        int idx_start = 0;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_np1_pml_xmin)(i,j,k) = (*A_n_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k);
                            (*A_n_pml_xmin)(i,j,k) = (*A_nm1_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k);
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
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_np1_pml_xmax  )(i,j,k) = (*A_n_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k);
                            (*A_n_pml_xmax  )(i,j,k) = (*A_nm1_)(idx_start+i,ncells_pml_domain-domain_oversize_y-nsolver/2+j,k);
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 3 ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );

        cField3D* A_np1_pml_xmin = NULL;
        cField3D* A_n_pml_xmin   = NULL;
        Field3D*  Chi_pml_xmin   = NULL;

        cField3D* A_np1_pml_xmax = NULL;
        cField3D* A_n_pml_xmax   = NULL;
        Field3D*  Chi_pml_xmax   = NULL;

        if(ncells_pml_xmin != 0){
            A_np1_pml_xmin = pml_fields_xmin->A_n_;
            A_n_pml_xmin   = pml_fields_xmin->A_nm1_;
            Chi_pml_xmin   = pml_fields_xmin->Chi_;
        }

        if(ncells_pml_xmax != 0){
            A_np1_pml_xmax = pml_fields_xmax->A_n_;
            A_n_pml_xmax   = pml_fields_xmax->A_nm1_;
            Chi_pml_xmax   = pml_fields_xmax->Chi_;
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
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*A_n_pml_xmin)(i,ny_p-j,k);
                        }
                    }
                }
            }
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*A_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*A_n_domain)(i,ny_p-j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        int idx_start = (ypml_size_in_x-1)-(ncells_pml_xmax-1) ;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_n_)(idx_start+i,domain_oversize_y+nsolver/2-j,k) = (*A_n_pml_xmax)(domain_oversize_x+nsolver/2+i,ny_p-j,k);
                        }
                    }
                }
            }
        }

        // Here we extrude Chi in PML domain according Chi at the ymax border
        for ( int j=0 ; j<ncells_pml_domain ; j++ ) {
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*Chi_)(idx_start+i,j,k) = (*Chi_domain)(i,ny_p-2-domain_oversize_y,k);
                }
            }
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    for ( int i=0 ; i<ncells_pml_xmin+domain_oversize_x ; i++ ) {
                        int idx_start = 0 ;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*Chi_)(idx_start+i,j,k) = (*Chi_pml_xmin)(i,ny_p-2-domain_oversize_y,k);
                        }
                    }
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    for ( int i=0 ; i<ncells_pml_xmax+domain_oversize_x ; i++ ) {
                        int idx_start = ypml_size_in_x-ncells_pml_xmax-domain_oversize_x ;
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*Chi_)(idx_start+i,j,k) = (*Chi_pml_xmax)(i,ny_p-2-domain_oversize_y,k);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int j=0 ; j < nsolver/2 ; j++){
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                int idx_start = ncells_pml_xmin;
                for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                    (*A_np1_domain)(i,ny_p-1-j,k) = (*A_n_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                    (*A_n_domain)(i,ny_p-1-j,k) = (*A_nm1_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                    (*Phi_domain)(i,ny_p-1-j,k) = 0.5*ellipticity_factor*std::abs( (*A_np1_domain)(i,ny_p-1-j,k) )*std::abs( (*A_np1_domain)(i,ny_p-1-j,k) );
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
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_np1_pml_xmin)(i,ny_p-1-j,k) = (*A_n_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*A_n_pml_xmin)(i,ny_p-1-j,k) = (*A_nm1_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
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
                        for ( unsigned int k=0 ; k<nz_p ; k++ ) {
                            (*A_np1_pml_xmax  )(i,ny_p-1-j,k) = (*A_n_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                            (*A_n_pml_xmax  )(i,ny_p-1-j,k) = (*A_nm1_)(idx_start+i,domain_oversize_y+nsolver/2-1-j,k);
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 4 ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );
        EnvelopeBC3D_PML* pml_fields_ymin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[2] );
        EnvelopeBC3D_PML* pml_fields_ymax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[3] );

        cField3D* A_np1_pml_xmin = NULL;
        cField3D* A_n_pml_xmin   = NULL;
        Field3D*  Chi_pml_xmin   = NULL;

        cField3D* A_np1_pml_xmax = NULL;
        cField3D* A_n_pml_xmax   = NULL;
        Field3D*  Chi_pml_xmax   = NULL;

        cField3D* A_np1_pml_ymin = NULL;
        cField3D* A_n_pml_ymin   = NULL;
        Field3D*  Chi_pml_ymin   = NULL;

        cField3D* A_np1_pml_ymax = NULL;
        cField3D* A_n_pml_ymax   = NULL;
        Field3D*  Chi_pml_ymax   = NULL;

        if(ncells_pml_xmin != 0){
            A_np1_pml_xmin = pml_fields_xmin->A_n_;
            A_n_pml_xmin   = pml_fields_xmin->A_nm1_;
            Chi_pml_xmin   = pml_fields_xmin->Chi_;
        }

        if(ncells_pml_xmax != 0){
            A_np1_pml_xmax = pml_fields_xmax->A_n_;
            A_n_pml_xmax   = pml_fields_xmax->A_nm1_;
            Chi_pml_xmax   = pml_fields_xmax->Chi_;
        }

        if(ncells_pml_ymin != 0){
            A_np1_pml_ymin = pml_fields_ymin->A_n_;
            A_n_pml_ymin   = pml_fields_ymin->A_nm1_;
            Chi_pml_ymin   = pml_fields_ymin->Chi_;
        }

        if(ncells_pml_ymax != 0){
            A_np1_pml_ymax = pml_fields_ymax->A_n_;
            A_n_pml_ymax   = pml_fields_ymax->A_nm1_;
            Chi_pml_ymax   = pml_fields_ymax->Chi_;
        }

        // 2. Exchange field PML <- Domain
        for ( int k=min2exchange ; k<max2exchange ; k++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*A_n_pml_xmin)(i,j,k);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            (*A_n_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*A_n_pml_ymin)(i,j,k);
                        }
                    }
                }
            }
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                    (*A_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*A_n_domain)(i,j,k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = (zpml_size_in_x-1)-(ncells_pml_xmax-1) ;
                    int jdx_start = ncells_pml_ymin ;
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*A_n_pml_xmax)(domain_oversize_x+nsolver/2+i,j,k);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = (zpml_size_in_y-1)-(ncells_pml_ymax-1) ;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            (*A_n_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k) = (*A_n_pml_ymax)(i,domain_oversize_y+nsolver/2+j,k);
                        }
                    }
                }
            }
        }

        // Here we extrude Chi in PML domain according Chi at the zmin border
        for ( int k=0 ; k<ncells_pml_domain ; k++ ) {
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                    (*Chi_)(idx_start+i,jdx_start+j,k) = 1*(*Chi_domain)(i,j,domain_oversize_z+1);
                }
            }
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmin+domain_oversize_x ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*Chi_)(idx_start+i,jdx_start+j,k) = (*Chi_pml_xmin)(i,j,domain_oversize_z+1);
                        }
                    }
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = zpml_size_in_x-ncells_pml_xmax-domain_oversize_x ;
                    int jdx_start = ncells_pml_ymin ;
                    for ( int i=0 ; i<ncells_pml_xmax+domain_oversize_x ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*Chi_)(idx_start+i,jdx_start+j,k) = (*Chi_pml_xmax)(i,j,domain_oversize_z+1);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            (*Chi_)(i,jdx_start+j,k) = (*Chi_pml_ymin)(i,j,domain_oversize_z+1);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = zpml_size_in_y-ncells_pml_ymax-domain_oversize_y ;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax+domain_oversize_y ; j++ ) {
                            (*Chi_)(i,jdx_start+j,k) = (*Chi_pml_ymax)(i,j,domain_oversize_z+1);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                    (*A_np1_domain)(i,j,k) = (*A_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                    (*A_n_domain)(i,j,k) = (*A_nm1_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                    (*Phi_domain)(i,j,k) = 0.5*ellipticity_factor*std::abs( (*A_np1_domain)(i,j,k) )*std::abs( (*A_np1_domain)(i,j,k) );
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
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_np1_pml_xmin)(i,j,k) = (*A_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*A_n_pml_xmin)(i,j,k) = (*A_nm1_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
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
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_np1_pml_xmax  )(i,j,k) = (*A_n_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*A_n_pml_xmax  )(i,j,k) = (*A_nm1_)(idx_start+i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
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
                        for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            (*A_np1_pml_ymin)(i,j,k) = (*A_n_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*A_n_pml_ymin)(i,j,k) = (*A_nm1_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
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
                        for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            (*A_np1_pml_ymax)(i,j,k) = (*A_n_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                            (*A_n_pml_ymax)(i,j,k) = (*A_nm1_)(i,jdx_start+j,ncells_pml_domain-domain_oversize_z-nsolver/2+k);
                        }
                    }
                }
            }
        }
    }

    else if( i_boundary_ == 5 ) {

        EnvelopeBC3D_PML* pml_fields_xmin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[0] );
        EnvelopeBC3D_PML* pml_fields_xmax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[1] );
        EnvelopeBC3D_PML* pml_fields_ymin = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[2] );
        EnvelopeBC3D_PML* pml_fields_ymax = static_cast<EnvelopeBC3D_PML*>( envelope->EnvBoundCond[3] );

        cField3D* A_np1_pml_xmin = NULL;
        cField3D* A_n_pml_xmin   = NULL;
        Field3D*  Chi_pml_xmin   = NULL;

        cField3D* A_np1_pml_xmax = NULL;
        cField3D* A_n_pml_xmax   = NULL;
        Field3D*  Chi_pml_xmax   = NULL;

        cField3D* A_np1_pml_ymin = NULL;
        cField3D* A_n_pml_ymin   = NULL;
        Field3D*  Chi_pml_ymin   = NULL;

        cField3D* A_np1_pml_ymax = NULL;
        cField3D* A_n_pml_ymax   = NULL;
        Field3D*  Chi_pml_ymax   = NULL;

        if(ncells_pml_xmin != 0){
            A_np1_pml_xmin = pml_fields_xmin->A_n_;
            A_n_pml_xmin   = pml_fields_xmin->A_nm1_;
            Chi_pml_xmin   = pml_fields_xmin->Chi_;
        }

        if(ncells_pml_xmax != 0){
            A_np1_pml_xmax = pml_fields_xmax->A_n_;
            A_n_pml_xmax   = pml_fields_xmax->A_nm1_;
            Chi_pml_xmax   = pml_fields_xmax->Chi_;
        }

        if(ncells_pml_ymin != 0){
            A_np1_pml_ymin = pml_fields_ymin->A_n_;
            A_n_pml_ymin   = pml_fields_ymin->A_nm1_;
            Chi_pml_ymin   = pml_fields_ymin->Chi_;
        }

        if(ncells_pml_ymax != 0){
            A_np1_pml_ymax = pml_fields_ymax->A_n_;
            A_n_pml_ymax   = pml_fields_ymax->A_nm1_;
            Chi_pml_ymax   = pml_fields_ymax->Chi_;
        }

    //     // Attention à la gestion des pas de temps
    //     // A priori pour calculer A_np1 on a juste besoin de A_n
    //     // Même les variables intermediaires sont inutiles car calculé à partir de An pour (i,j)

    //     // 2. Exchange field PML <- Domain
        for ( int k=min2exchange ; k<max2exchange ; k++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmin ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A_n_pml_xmin)(i,j,nz_p-k);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            (*A_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A_n_pml_ymin)(i,j,nz_p-k);
                        }
                    }
                }
            }
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                    (*A_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A_n_domain)(i,j,nz_p-k);
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = (zpml_size_in_x-1)-(ncells_pml_xmax-1) ;
                    int jdx_start = ncells_pml_ymin ;
                    for ( int i=0 ; i<ncells_pml_xmax ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A_n_pml_xmax)(domain_oversize_x+nsolver/2+i,j,nz_p-k);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = (zpml_size_in_y-1)-(ncells_pml_ymax-1) ;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax ; j++ ) {
                            (*A_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-k) = (*A_n_pml_ymax)(i,domain_oversize_y+nsolver/2+j,nz_p-k);
                        }
                    }
                }
            }
        }

        // Here we extrude Chi in PML domain according Chi at the zmax border
        for ( int k=0 ; k<ncells_pml_domain ; k++ ) {
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                    (*Chi_)(idx_start+i,jdx_start+j,k) = 1*(*Chi_domain)(i,j,nz_p-2-domain_oversize_z);
                }
            }
            if (patch->isXmin()) {
                if(ncells_pml_xmin != 0){
                    int idx_start = 0;
                    int jdx_start = ncells_pml_ymin;
                    for ( int i=0 ; i<ncells_pml_xmin+domain_oversize_x ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*Chi_)(idx_start+i,jdx_start+j,k) = (*Chi_pml_xmin)(i,j,nz_p-2-domain_oversize_z);
                        }
                    }
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_xmax != 0){
                    int idx_start = zpml_size_in_x-ncells_pml_xmax-domain_oversize_x ;
                    int jdx_start = ncells_pml_ymin ;
                    for ( int i=0 ; i<ncells_pml_xmax+domain_oversize_x ; i++ ) {
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*Chi_)(idx_start+i,jdx_start+j,k) = (*Chi_pml_xmax)(i,j,nz_p-2-domain_oversize_z);
                        }
                    }
                }
            }
            if (patch->isYmin()) {
                if(ncells_pml_ymin != 0){
                    int jdx_start = 0;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymin ; j++ ) {
                            (*Chi_)(i,jdx_start+j,k) = (*Chi_pml_ymin)(i,j,nz_p-2-domain_oversize_z);
                        }
                    }
                }
            }
            if (patch->isYmax()) {
                if(ncells_pml_ymax != 0){
                    int jdx_start = zpml_size_in_y-ncells_pml_ymax-domain_oversize_y ;
                    for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                        for ( int j=0 ; j<ncells_pml_ymax+domain_oversize_y ; j++ ) {
                            (*Chi_)(i,jdx_start+j,k) = (*Chi_pml_ymax)(i,j,nz_p-2-domain_oversize_z);
                        }
                    }
                }
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope, iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        for (int k=0 ; k < nsolver/2 ; k++){
            //std::cout << k << std::endl ;
            int idx_start = ncells_pml_xmin;
            int jdx_start = ncells_pml_ymin;
            for ( unsigned int i=0 ; i<nx_p ; i++ ) {
                for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                    (*A_np1_domain)(i,j,nz_p-1-k) = (*A_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                    (*A_n_domain)(i,j,nz_p-1-k) = (*A_nm1_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                    (*Phi_domain)(i,j,nz_p-1-k) = 0.5*ellipticity_factor*std::abs( (*A_np1_domain)(i,j,nz_p-1-k) )*std::abs( (*A_np1_domain)(i,j,nz_p-1-k) );
                    if (k == 1) {
                    //std::cout << nz_p << std::endl ;
                    //std::cout << "(" << idx_start+i << "," << jdx_start+j << "," << domain_oversize_z+nsolver/2-1-k << ")" << "=>" << "(" << i << "," << j << "," << nz_p-1-k << ")" << std::endl ;
                    }
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
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_np1_pml_xmin)(i,j,nz_p-1-k) = (*A_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*A_n_pml_xmin)(i,j,nz_p-1-k) = (*A_nm1_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
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
                        for ( unsigned int j=0 ; j<ny_p ; j++ ) {
                            (*A_np1_pml_xmax  )(i,j,nz_p-1-k) = (*A_n_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*A_n_pml_xmax  )(i,j,nz_p-1-k) = (*A_nm1_)(idx_start+i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
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
                        for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax ; i++ ) {
                            (*A_np1_pml_ymin)(i,j,nz_p-1-k) = (*A_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*A_n_pml_ymin)(i,j,nz_p-1-k) = (*A_nm1_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
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
                        for ( unsigned int i=0 ; i<nx_p+ncells_pml_xmin+ncells_pml_xmax; i++ ) {
                            (*A_np1_pml_ymax)(i,j,nz_p-1-k) = (*A_n_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                            (*A_n_pml_ymax)(i,j,nz_p-1-k) = (*A_nm1_)(i,jdx_start+j,domain_oversize_z+nsolver/2-1-k);
                        }
                    }
                }
            }
        }
    }
}
