#include "PeekAtSpecies.h"
#include "Hilbert_functions.h"

using namespace std;


PeekAtSpecies::PeekAtSpecies( Params &p, unsigned int species_id ) :
    density_profile_( NULL ),
    particles_per_cell_profile_( NULL ),
    params( &p )
{
    // Obtain the profiles of the species
    std::string species_name( "" );
    PyTools::extractOrNone( "name", species_name, "Species", species_id );
    PyObject *profile1=nullptr;
    std::string density_profile_type( "" );
    std::string peek_position_initialization;
    PyObject *py_pos_init = PyTools::extract_py( "position_initialization", "Species", species_id );
    if( PyTools::py2scalar( py_pos_init, peek_position_initialization ) ) {
        if( peek_position_initialization == "regular"
         || peek_position_initialization == "random"
         || peek_position_initialization == "centered" ) {
            bool ok1 = PyTools::extract_pyProfile( "number_density", profile1, "Species", species_id );
            bool ok2 = PyTools::extract_pyProfile( "charge_density", profile1, "Species", species_id );
            if( ok1 ) {
                density_profile_type = "nb";
            }
            if( ok2 ) {
                density_profile_type = "charge";
            }
            if( ! ( ok1 ^ ok2 ) ) {
                ERROR( "Species `" << species_name << "`: Missing density profile (`charge_density` or `number_density` option)" );
            }
            density_profile_ = new Profile( profile1, params->nDim_field, Tools::merge( density_profile_type, "_density ", species_name ), *params, true, true );
            PyTools::extract_pyProfile( "particles_per_cell", profile1, "Species", species_id );
            particles_per_cell_profile_ = new Profile( profile1, params->nDim_field, Tools::merge( "particles_per_cell ", species_name ), *params, true, true );
        }
    }
    Py_DECREF( py_pos_init );
}


PeekAtSpecies::~PeekAtSpecies()
{
    delete density_profile_;
    delete particles_per_cell_profile_;
}


double PeekAtSpecies::numberOfParticlesInPatch( unsigned int hindex )
{
    if( particles_per_cell_profile_ ) {
        // Get the patch x_cell
        unsigned int Px_cell[3];
        generalhilbertindexinv( params->mi[0], params->mi[1], params->mi[2], &Px_cell[0], &Px_cell[1], &Px_cell[2], hindex );
        // Get the patch center
        vector<Field *> x_cell( params->nDim_field );
        vector<unsigned int> n = {1, 1, 1};
        for( unsigned int i=0 ; i<params->nDim_field ; i++ ) {
            x_cell[i] = new Field3D( n );
            (*x_cell[i])(0) = ( ( double )( Px_cell[i] )+0.5 ) * params->patch_dimensions[i];
        }
        // Evaluate profile
        vector<double> global_origin( 3, 0. );
        Field3D nppc( n ), dens( n );
        particles_per_cell_profile_->valuesAt( x_cell, global_origin, nppc );
        density_profile_->valuesAt( x_cell, global_origin, dens );
        double n_part_in_cell = floor( nppc( 0, 0, 0 ) );
        if( n_part_in_cell<=0. || dens( 0, 0, 0 )==0. ) {
            n_part_in_cell = 0.;
        }
        for( unsigned int i=0 ; i<params->nDim_field ; i++ ) {
            delete x_cell[i] ;
        }
        return n_part_in_cell * params->n_cell_per_patch;
    } else {
        return 0. ;  //Do not account for particles initialized from numpy array
    }
}


double PeekAtSpecies::numberOfParticlesInPatch( vector<double> x )
{
    // Evaluate the profile at that location
    if( particles_per_cell_profile_ ) {
        vector<Field *> x_cell( params->nDim_field );
        vector<unsigned int> n = {1, 1, 1};
        for( unsigned int i=0 ; i<params->nDim_field ; i++ ) {
            x_cell[i] = new Field3D( n );
            (*x_cell[i])(0) = x[i];
        }
        vector<double> global_origin( 3, 0. );
        Field3D nppc( n ), dens( n );
        particles_per_cell_profile_->valuesAt( x_cell, global_origin, nppc );
        density_profile_->valuesAt( x_cell, global_origin, dens );
        double n_part_in_cell = floor( nppc( 0, 0, 0 ) );
        if( n_part_in_cell<=0. || dens( 0, 0, 0 )==0. ) {
            n_part_in_cell = 0.;
        }
        for( unsigned int i=0 ; i<params->nDim_field ; i++ ) {
            delete x_cell[i] ;
        }
        return n_part_in_cell * params->n_cell_per_patch;
    } else {
        return 0. ;  //Do not account for particles initialized from numpy array
    }
}


double PeekAtSpecies::totalNumberofParticles()
{
    // Loop over the box to obtain an approximate number of particles
    vector<double> x_cell( params->nDim_field, 0. );
    if( particles_per_cell_profile_ ) {
        vector<unsigned int> i_cell( params->nDim_field, 0 );
        vector<Field *> x_cell( params->nDim_field );
        vector<unsigned int> n = {1, 1, 1};
        vector<double> global_origin( 3, 0. );
        Field3D nppc( n ), dens( n );
        for( unsigned int idim=0; idim<params->nDim_field; idim++ ) {
            i_cell[idim] = 0;
            x_cell[idim] = new Field3D( n );
            (*x_cell[idim])(0) = params->patch_dimensions[idim] * 0.5;
        }
        // Loop patches
        double npart_total = 0.;
        for( unsigned int k=0; k<params->tot_number_of_patches; k++ ) {
            // Find the approximate number of particles in this patch
            particles_per_cell_profile_->valuesAt( x_cell, global_origin, nppc );
            density_profile_->valuesAt( x_cell, global_origin, dens );
            double n_part_in_cell = floor( nppc( 0, 0, 0 ) );
            if( n_part_in_cell>0. && dens( 0, 0, 0 )!=0. ) {
                npart_total += n_part_in_cell * params->n_cell_per_patch;
            }
            // Find next patch position
            for( unsigned int idim=0; idim<params->nDim_field; idim++ ) {
                if( i_cell[idim] < params->number_of_patches[idim]-1 ) {
                    i_cell[idim]++;
                    (*x_cell[idim])(0) += params->patch_dimensions[idim];
                    break;
                } else if( idim < params->nDim_field-1 ) {
                    i_cell[idim] = 0;
                    (*x_cell[idim])(0) = params->patch_dimensions[idim] * 0.5;
                }
            }
        }
        for( unsigned int i=0 ; i<params->nDim_field ; i++ ) {
            delete x_cell[i] ;
        }
        
        return npart_total;
    } else {
        return 0. ;  //Do not account for particles initialized from numpy array
    }
}
