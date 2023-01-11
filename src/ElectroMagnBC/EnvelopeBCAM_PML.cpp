#include "EnvelopeBCAM_PML.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "LaserEnvelope.h"
#include "Field2D.h"
#include "cField2D.h"
#include "Tools.h"
#include <complex>
#include "dcomplex.h"

#include "SolverFactory.h"

using namespace std;

EnvelopeBCAM_PML::EnvelopeBCAM_PML( Params &params, Patch *patch, unsigned int i_boundary )
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
        nsolver=2;
    }

    if( patch->isBoundary( i_boundary_ ) ) {

        int iDim = i_boundary / 2;
        int min_or_max = i_boundary_ % 2;

        domain_oversize_l =  patch->oversize[0] ;
        domain_oversize_r =  patch->oversize[1] ;

        if (patch->isXmin() ) {//&& i_boundary_ == 0 ) {
            ncells_pml_lmin = params.number_of_pml_cells[0][0];
            ncells_pml_domain_lmin = ncells_pml_lmin + 1*patch->oversize[0] + nsolver/2;
            domain_oversize_l = patch->oversize[0] ;
        }
        else {
            ncells_pml_lmin = 0;
            ncells_pml_domain_lmin = 0;
        }
        if (patch->isXmax() ) {//&& i_boundary_ == 1 ) {
            ncells_pml_lmax = params.number_of_pml_cells[0][1];
            ncells_pml_domain_lmax = ncells_pml_lmax + 1*patch->oversize[0] + nsolver/2;
            domain_oversize_l = patch->oversize[0] ;
        }
        else {
            ncells_pml_lmax = 0;
            ncells_pml_domain_lmax = 0;
        }
        if (patch->isYmin() ) {//&& i_boundary_ == 2 ) {
            ncells_pml_rmin = params.number_of_pml_cells[1][0];
            ncells_pml_domain_rmin = ncells_pml_rmin + 1*patch->oversize[1] + nsolver/2;
            domain_oversize_r = patch->oversize[1] ;
        }
        else {
            ncells_pml_rmin = 0;
            ncells_pml_domain_rmin = 0;
        }
        if (patch->isYmax() ) {//&& i_boundary_ == 3 ) {
            ncells_pml_rmax = params.number_of_pml_cells[1][1];
            ncells_pml_domain_rmax = ncells_pml_rmax + 1*patch->oversize[1] + nsolver/2;
            domain_oversize_r = patch->oversize[1] ;
        }
        else {
            ncells_pml_rmax = 0;
            ncells_pml_domain_rmax = 0;
        }

        ncells_pml = params.number_of_pml_cells[iDim][min_or_max];
        ncells_pml_domain = ncells_pml+1*patch->oversize[iDim] + nsolver/2;

        // Define min and max idx to exchange
        // the good data f(solver,patch->oversize)
        if (min_or_max==0){
            // if min border : Exchange of data (for domain to pml-domain)
            // min2exchange <= i < max2exchange
            min2exchange = 1*nsolver/2 ;
            max2exchange = 2*nsolver/2 ;
            // Solver
            solvermin = nsolver/2 ;
            // solvermax = ncells_pml_domain - 1 - patch->oversize[iDim] ;
            solvermax = ncells_pml_domain-patch->oversize[iDim] ;
        }
        else if (min_or_max==1){
            // if max border : Exchange of data (for domain to pml-domain)
            // min2exchange <= i < max2exchange
            // min2exchange = 1*nsolver/2 ;
            // max2exchange = 2*nsolver/2 ;
            min2exchange = 1*nsolver/2+1 ;
            max2exchange = 2*nsolver/2+1 ;
            // Solver
            //solvermin = patch->oversize[iDim] + nsolver/2 - nsolver/2 + 1 ;
            solvermin = patch->oversize[iDim];
            solvermax = ncells_pml_domain-nsolver/2 ;
        }

        if (ncells_pml==0){
            ERROR("PML domain have to be >0 cells in thickness");
        }
        
        if( iDim == 0 ) {
            dimPrim = { (unsigned int)ncells_pml_domain, patch->size_[1]+1+2*patch->oversize[1] };
        } else {
            dimPrim = { patch->size_[0]+1+2*patch->oversize[0] + ncells_pml_lmin + ncells_pml_lmax, (unsigned int)ncells_pml_domain };
            rpml_size_in_l = dimPrim[0];
        }

        startpml = patch->oversize[iDim]+nsolver/2;

        int ncells_pml_min[1];
        ncells_pml_min[0] = ncells_pml_lmin;
        int ncells_pml_max[1];
        ncells_pml_max[0] = ncells_pml_lmax;

        pml_solver_envelope_->setDomainSizeAndCoefficients( iDim, min_or_max, dimPrim, ncells_pml_domain, startpml, ncells_pml_min, ncells_pml_max, patch );

        std::string si_boundary = std::to_string(i_boundary_);
        // A-field
        A_np1_ = new cField2D( dimPrim, "A_np1_pml"+si_boundary );
        A_n_ = new cField2D( dimPrim, "A_n_pml"+si_boundary );
        A_nm1_ = new cField2D( dimPrim, "A_nm1_pml"+si_boundary );
        // G-field
        G_np1_ = new cField2D( dimPrim, "G_np1_pml"+si_boundary );
        G_n_ = new cField2D( dimPrim, "G_n_pml"+si_boundary );
        G_nm1_ = new cField2D( dimPrim, "G_nm1_pml"+si_boundary );
        // Auxillary Variable
        u1_np1_l_ = new cField2D( dimPrim, "u1_np1_l_pml"+si_boundary );
        u2_np1_l_ = new cField2D( dimPrim, "u2_np1_l_pml"+si_boundary );
        u3_np1_l_ = new cField2D( dimPrim, "u3_np1_l_pml"+si_boundary );
        u1_nm1_l_ = new cField2D( dimPrim, "u1_nm1_l_pml"+si_boundary );
        u2_nm1_l_ = new cField2D( dimPrim, "u2_nm1_l_pml"+si_boundary );
        u3_nm1_l_ = new cField2D( dimPrim, "u3_nm1_l_pml"+si_boundary );
        // ----
        u1_np1_r_ = new cField2D( dimPrim, "u1_np1_r_pml"+si_boundary );
        u2_np1_r_ = new cField2D( dimPrim, "u2_np1_r_pml"+si_boundary );
        u3_np1_r_ = new cField2D( dimPrim, "u3_np1_r_pml"+si_boundary );
        u1_nm1_r_ = new cField2D( dimPrim, "u1_nm1_r_pml"+si_boundary );
        u2_nm1_r_ = new cField2D( dimPrim, "u2_nm1_r_pml"+si_boundary );
        u3_nm1_r_ = new cField2D( dimPrim, "u3_nm1_r_pml"+si_boundary );

        // Ponderomoteur Potential
        Phi_ = new Field2D( dimPrim, "Phi_pml"+si_boundary );
        Chi_ = new Field2D( dimPrim, "Chi_pml"+si_boundary );
    }

    j_glob_pml = patch->getCellStartingGlobalIndex( 1 );
    // To do for Ymax in order to be coherent with the solver
    // j_glob_pml = patch->getCellStartingGlobalIndex( 1 )+nr_p-patch->oversize[iDim]-3;
    // std::cout << j_glob_pml << std::endl;

}


EnvelopeBCAM_PML::~EnvelopeBCAM_PML()
{
    delete A_np1_;
    delete A_n_;
    delete A_nm1_;
    // ----
    delete G_np1_;
    delete G_n_;
    delete G_nm1_;
    // ----
    delete u1_np1_l_;
    delete u2_np1_l_;
    delete u3_np1_l_;
    delete u1_nm1_l_;
    delete u2_nm1_l_;
    delete u3_nm1_l_;
    // ----
    delete u1_np1_r_;
    delete u2_np1_r_;
    delete u3_np1_r_;
    delete u1_nm1_r_;
    delete u2_nm1_r_;
    delete u3_nm1_r_;
    // ----
    delete Phi_;
    delete Chi_;

    if (pml_solver_envelope_!=NULL) {
        delete pml_solver_envelope_;
    }
}

/*
void EnvelopeBCAM_PML::save_fields( Field *my_field, Patch *patch )
{
}


void EnvelopeBCAM_PML::disableExternalFields()
{
}
*/

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------

void EnvelopeBCAM_PML::apply( LaserEnvelope *envelope, ElectroMagn *EMfields, Patch *patch )
{
    int iDim = i_boundary_ / 2;
    int min_or_max = (i_boundary_)%2;

    cField2D *A_np1_domain  = static_cast<cField2D *>( envelope->A_ );  // A_ is the envelope at timestep n BUT at this point A_ is already update so in fact its correspond to A_np1_domain
    cField2D *A_n_domain    = static_cast<cField2D *>( envelope->A0_ ); // A0_ is the envelope at timestep n-1 BUT at this point A0_ is already update so in fact its A_n_domain
    Field2D  *Phi_domain    = static_cast<Field2D *>( envelope->Phi_ ); // the ponderomotive potential Phi=|A|^2/2 at timestep n
    Field2D  *Chi_domain    = static_cast<Field2D *>( EMfields->Env_Chi_ );

    double ellipticity_factor = envelope->ellipticity_factor;

    if( ! patch->isBoundary( i_boundary_ ) ) return;

    if( i_boundary_ == 0 ) {

        // 2. Exchange field PML <- Domain
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            // for ( int j=1 ; j<nr_p-1 ; j++ ) {
            for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                // (*A_n_)(ncells_pml_domain-1-domain_oversize_l-nsolver/2+i,j) = (*A_n_domain)(i,j);
                // (*G_n_)(ncells_pml_domain-1-domain_oversize_l-nsolver/2+i,j) = (*A_n_domain)(i,j)*( (double) (j_glob_pml+j)*dr );
                (*A_n_)(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*A_n_domain)(i,j);
                (*G_n_)(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j) = (*A_n_domain)(i,j)*( (double) (j_glob_pml+j)*dr );
            }
        }

        // Source Chi Update (1 cells after xmin)
        for ( int i=0 ; i<ncells_pml_domain ; i++ ) {
            for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                (*Chi_)(i,j) = (*Chi_domain)(domain_oversize_l+1,j);
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope,iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Primals in x-direction
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                // (*A_np1_domain)(i,j) = (*A_n_)(ncells_pml_domain-1-domain_oversize_l-nsolver/2+i,j);
                // (*A_n_domain)(i,j) = (*A_nm1_)(ncells_pml_domain-1-domain_oversize_l-nsolver/2+i,j);
                (*A_np1_domain)(i,j) = (*A_n_)(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                (*A_n_domain)(i,j) = (*A_nm1_)(ncells_pml_domain-domain_oversize_l-nsolver/2+i,j);
                (*Phi_domain)(i,j) = 0.5*ellipticity_factor*std::abs(  (*A_np1_domain)(i,j) )*std::abs(  (*A_np1_domain)(i,j) );
            }
        }
    }

    else if( i_boundary_ == 1 ) {

        // 2. Exchange field Domain -> PML
        for ( int i=min2exchange ; i<max2exchange ; i++ ) {
            // for ( int j=1 ; j<nr_p-1 ; j++ ) {
            for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                // (*A_n_)(domain_oversize_l+nsolver/2-i,j) = (*A_n_domain)(nl_p-1-i,j);
                // (*G_n_)(domain_oversize_l+nsolver/2-i,j) = (*A_n_domain)(nl_p-1-i,j)*( (double) (j_glob_pml+j)*dr );
                (*A_n_)(domain_oversize_l+nsolver/2-i,j) = (*A_n_domain)(nl_p-i,j);
                (*G_n_)(domain_oversize_l+nsolver/2-i,j) = (*A_n_domain)(nl_p-i,j)*( (double) (j_glob_pml+j)*dr );
                // std::cout << "jglob" << j_glob_pml << std::endl;
                // std::cout << "j" << j << std::endl;
                // std::cout << "tot" << j_glob_pml+j << std::endl;
                // std::cout << "qte" << (*A_n_domain)(nl_p-i,j)*( (double) (j_glob_pml+j)*dr ) << std::endl;
            }
        }
        // Source Chi Update (1 cells before xmax)
        for ( int i=0 ; i<ncells_pml_domain ; i++ ) {
            for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                (*Chi_)(i,j) = (*Chi_domain)(nl_p-2-domain_oversize_l,j);
            }
        }

        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope,iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange Domain -> PML
        // Primals in x-direction
        // for (int i=0 ; i < nsolver/2 ; i++){
        for (int i=0 ; i < nsolver/2 ; i++){
            for ( unsigned int j=0 ; j<nr_p ; j++ ) {
                // (*A_np1_domain)(nl_p-1-i,j) = (*A_n_)(domain_oversize_l+nsolver/2-i,j);
                // (*A_n_domain)(nl_p-1-i,j) = (*A_nm1_)(domain_oversize_l+nsolver/2-i,j);
                (*A_np1_domain)(nl_p-1-i,j) = (*A_n_)(domain_oversize_l+nsolver/2-1-i,j);
                (*A_n_domain)(nl_p-1-i,j) = (*A_nm1_)(domain_oversize_l+nsolver/2-1-i,j);
                (*Phi_domain)(nl_p-1-i,j) = 0.5*ellipticity_factor*std::abs( (*A_np1_domain)(nl_p-1-i,j) )*std::abs( (*A_np1_domain)(nl_p-1-i,j) ) ;
            }
        }
    }

    else if( i_boundary_ == 2 ) {

        // NO BC on axis here

    }

    else if( i_boundary_ == 3 ) {

        EnvelopeBCAM_PML* pml_fields_lmin = NULL;
        EnvelopeBCAM_PML* pml_fields_lmax = NULL;

        if(ncells_pml_lmin != 0){
            pml_fields_lmin = static_cast<EnvelopeBCAM_PML*>( envelope->EnvBoundCond[0] );
        }

        if(ncells_pml_lmax != 0){
            pml_fields_lmax = static_cast<EnvelopeBCAM_PML*>( envelope->EnvBoundCond[1] );
        }

        cField2D* A_np1_pml_lmin = NULL;
        cField2D* A_n_pml_lmin   = NULL;
        cField2D* G_np1_pml_lmin = NULL;
        cField2D* G_n_pml_lmin   = NULL;
        Field2D*  Chi_pml_lmin   = NULL;

        cField2D* A_np1_pml_lmax = NULL;
        cField2D* A_n_pml_lmax   = NULL;
        cField2D* G_np1_pml_lmax = NULL;
        cField2D* G_n_pml_lmax   = NULL;
        Field2D*  Chi_pml_lmax   = NULL;

        if(ncells_pml_lmin != 0){
            A_np1_pml_lmin = pml_fields_lmin->A_n_;
            A_n_pml_lmin   = pml_fields_lmin->A_nm1_;
            G_np1_pml_lmin = pml_fields_lmin->G_n_;
            G_n_pml_lmin   = pml_fields_lmin->G_nm1_;
            Chi_pml_lmin   = pml_fields_lmin->Chi_;
        }

        if(ncells_pml_lmax != 0){
            A_np1_pml_lmax = pml_fields_lmax->A_n_;
            A_n_pml_lmax   = pml_fields_lmax->A_nm1_;
            G_np1_pml_lmax = pml_fields_lmax->G_n_;
            G_n_pml_lmax   = pml_fields_lmax->G_nm1_;
            Chi_pml_lmax   = pml_fields_lmax->Chi_;
        }

        // Attention à la gestion des pas de temps
        // A priori pour calculer A_np1 on a juste besoin de A_n
        // Même les variables intermediaires sont inutiles car calculé à partir de An pour (i,j)

        // 2. Exchange field PML <- Domain
        for ( int j=min2exchange ; j<max2exchange ; j++ ) {
            if (patch->isXmin()) {
                if(ncells_pml_lmin != 0){
                    for ( int i=0 ; i<ncells_pml_lmin ; i++ ) {
                        int idx_start = 0;
                        // Les qtes Primals
                        // (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_pml_lmin)(i,nr_p-1-j);
                        // (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*G_n_pml_lmin)(i,nr_p-1-j);
                        (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_pml_lmin)(i,nr_p-j);
                        (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*G_n_pml_lmin)(i,nr_p-j);
                    }
                }
            }
            // for ( int i=1 ; i<nl_p-1 ; i++ ) {
            for ( unsigned int i=0 ; i<nl_p ; i++ ) {
                // int idx_start = ncells_pml_lmin-1*(patch->isXmin());
                int idx_start = ncells_pml_lmin;
                // (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_domain)(i,nr_p-1-j);
                // (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_domain)(i,nr_p-1-j)*( (double) (j_glob_pml+nr_p-1-j)*dr );
                (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_domain)(i,nr_p-j);
                (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_domain)(i,nr_p-j)*( (double) (j_glob_pml+nr_p-j)*dr );
            }
            if (patch->isXmax()) {
                if(ncells_pml_lmax != 0){
                    for ( int i=0 ; i<ncells_pml_lmax ; i++ ) {
                        int idx_start = (rpml_size_in_l-1)-(ncells_pml_lmax-1) ;
                        // (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_pml_lmax)(domain_oversize_l+nsolver/2+i,nr_p-1-j);
                        // (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*G_n_pml_lmax)(domain_oversize_l+nsolver/2+i,nr_p-1-j);
                        (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*A_n_pml_lmax)(domain_oversize_l+nsolver/2+i,nr_p-j);
                        (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-j) = (*G_n_pml_lmax)(domain_oversize_l+nsolver/2+i,nr_p-j);
                    }
                }
            }
        }
 
        // Source Chi Update (1 cells before rmax)
        for ( int j=0 ; j<ncells_pml_domain ; j++ ) {
            for ( unsigned int i=0 ; i<nl_p ; i++ ) {
                int idx_start = ncells_pml_lmin;
                (*Chi_)(idx_start+i,j) = 1*(*Chi_domain)(i,nr_p-2-domain_oversize_r);
            }
            if (patch->isXmin()) {
                if(ncells_pml_lmin != 0){
                    for ( int i=0 ; i<ncells_pml_lmin+domain_oversize_l ; i++ ) {
                        int idx_start = 0 ;
                        (*Chi_)(idx_start+i,j) = 1*(*Chi_pml_lmin)(i,nr_p-2-domain_oversize_r);
                    }
                }
            }
            if (patch->isXmax()) {
                if(ncells_pml_lmax != 0){
                    for ( int i=0 ; i<ncells_pml_lmax+domain_oversize_l ; i++ ) {
                        int idx_start = rpml_size_in_l-ncells_pml_lmax-domain_oversize_l ;
                        (*Chi_)(idx_start+i,j) = 1*(*Chi_pml_lmax)(i,nr_p-2-domain_oversize_r);
                    }
                }
            }
        }


        // 3. Solve Maxwell_PML for A-field :
        pml_solver_envelope_->compute_A_from_G( envelope,iDim, min_or_max, dimPrim, solvermin, solvermax);

        // 4. Exchange PML -> Domain
        // Duals in y-direction
        // for (int j=0 ; j < nsolver/2 ; j++){
        for (int j=0 ; j < nsolver/2 ; j++){
            // for ( int i=1 ; i<nl_p-1 ; i++ ) {
            for ( unsigned int i=0 ; i<nl_p ; i++ ) {
                // int idx_start = ncells_pml_lmin-1*(patch->isXmin());
                int idx_start = ncells_pml_lmin;
                // (*A_np1_domain)(i,nr_p-1-j) = (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-j);
                // (*A_n_domain)(i,nr_p-1-j) = (*A_nm1_)(idx_start+i,domain_oversize_r+nsolver/2-j);
                (*A_np1_domain)(i,nr_p-1-j) = (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                (*A_n_domain)(i,nr_p-1-j) = (*A_nm1_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                (*Phi_domain)(i,nr_p-1-j) = 0.5*ellipticity_factor*std::abs( (*A_np1_domain)(i,nr_p-1-j) )*std::abs( (*A_np1_domain)(i,nr_p-1-j) );
            }
        }

        // 5. Exchange PML y -> PML x MIN
        // Primals in y-direction
        // Dual in y-direction
        // for (int j=0 ; j < nsolver/2 ; j++){
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmin()) {
                if(ncells_pml_lmin != 0){
                    for ( int i=0 ; i<ncells_pml_domain_lmin ; i++ ) {
                        int idx_start = 0;
                        // Primals
                        (*A_np1_pml_lmin)(i,nr_p-1-j) = (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                        (*A_n_pml_lmin)(i,nr_p-1-j) = (*A_nm1_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                        (*G_np1_pml_lmin)(i,nr_p-1-j) = (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                        (*G_n_pml_lmin)(i,nr_p-1-j) = (*G_nm1_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                    }
                }
            }
        }

        // 6. Exchange PML y -> PML x MAX
        // Primals in y-direction
        // Duals in y-direction
        // for (int j=0 ; j < nsolver/2 ; j++){
        for (int j=0 ; j < nsolver/2 ; j++){
            if (patch->isXmax()) {
                if(ncells_pml_lmax != 0){
                    for ( int i=0 ; i<ncells_pml_domain_lmax ; i++ ) {
                        int idx_start = rpml_size_in_l-ncells_pml_domain_lmax ;
                        // Les qtes Primals
                        (*A_np1_pml_lmax  )(i,nr_p-1-j) = (*A_n_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                        (*A_n_pml_lmax  )(i,nr_p-1-j) = (*A_nm1_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                        (*G_np1_pml_lmax  )(i,nr_p-1-j) = (*G_n_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                        (*G_n_pml_lmax  )(i,nr_p-1-j) = (*G_nm1_)(idx_start+i,domain_oversize_r+nsolver/2-1-j);
                    }
                }
            }
        }
    }
}
