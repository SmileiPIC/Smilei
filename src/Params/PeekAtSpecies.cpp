#include "PeekAtSpecies.h"
#include "Hilbert_functions.h"

using namespace std;


PeekAtSpecies::PeekAtSpecies(Params& p, unsigned int species_id) :
    densityProfile(NULL),
    ppcProfile(NULL),
    params( &p )
{
    // Obtain the profiles of the species
    std::string species_name("");
    PyTools::extract("name",species_name,"Species",species_id);
    PyObject *profile1=nullptr;
    std::string densityProfileType("");
    std::string peek_position_initialization;
    PyObject *py_pos_init = PyTools::extract_py("position_initialization", "Species",species_id);
    if ( PyTools::convert(py_pos_init, peek_position_initialization) ){
        bool ok1 = PyTools::extract_pyProfile("number_density", profile1, "Species", species_id);
        bool ok2 = PyTools::extract_pyProfile("charge_density", profile1, "Species", species_id);
        if( ok1 ) densityProfileType = "nb";
        if( ok2 ) densityProfileType = "charge";
        densityProfile = new Profile(profile1, params->nDim_field, Tools::merge(densityProfileType,"_density ",species_name));
        PyTools::extract_pyProfile("particles_per_cell", profile1, "Species", species_id);
        ppcProfile = new Profile(profile1, params->nDim_field, Tools::merge("particles_per_cell ",species_name));
    }
    Py_DECREF(py_pos_init);
}


PeekAtSpecies::~PeekAtSpecies()
{
    delete densityProfile;
    delete ppcProfile;
}


double PeekAtSpecies::numberOfParticlesInPatch( unsigned int hindex )
{
    // Get the patch coordinates
    unsigned int Pcoordinates[3];
    vector<double> x_cell( params->nDim_field, 0. );
    generalhilbertindexinv(params->mi[0], params->mi[1], params->mi[2], &Pcoordinates[0], &Pcoordinates[1], &Pcoordinates[2], hindex);
    for (unsigned int i=0 ; i<params->nDim_field ; i++) {
        x_cell[i] = ((double)(Pcoordinates[i])+0.5) * params->patch_dimensions[i];
    }
    // Evaluate the profile at that location
    if (ppcProfile){
        double n_part_in_cell = floor(ppcProfile->valueAt(x_cell));
        if( n_part_in_cell<=0. || densityProfile->valueAt(x_cell)==0. )
            n_part_in_cell = 0.;
        return n_part_in_cell * params->n_cell_per_patch;
    } else {
        return 0. ;  //Do not account for particles initialized from numpy array
    }
}


double PeekAtSpecies::numberOfParticlesInPatch( vector<double> x_cell )
{
    // Evaluate the profile at that location
    if (ppcProfile){
        double n_part_in_cell = floor(ppcProfile->valueAt(x_cell));
        if( n_part_in_cell<=0. || densityProfile->valueAt(x_cell)==0. )
            n_part_in_cell = 0.;
        return n_part_in_cell * params->n_cell_per_patch;
    } else {
        return 0. ;  //Do not account for particles initialized from numpy array
    }
}


double PeekAtSpecies::totalNumberofParticles()
{
    // Loop over the box to obtain an approximate number of particles
    vector<unsigned int> i_cell( params->nDim_field, 0 );
    vector<double> x_cell( params->nDim_field, 0. );
    if (ppcProfile){
        for( unsigned int idim=0; idim<params->nDim_field; idim++ ) {
            i_cell[idim] = 0;
            x_cell[idim] = params->patch_dimensions[idim] * 0.5;
        }
        // Loop patches
        double npart_total = 0.;
        for( unsigned int k=0; k<params->tot_number_of_patches; k++ ) {
            // Find the approximate number of particles in this patch
            double n_part_in_cell = floor(ppcProfile->valueAt(x_cell));
            if( n_part_in_cell>0. && densityProfile->valueAt(x_cell)!=0.)
                npart_total += n_part_in_cell * params->n_cell_per_patch;
            // Find next patch position
            for( unsigned int idim=0; idim<params->nDim_field; idim++ ) {
                if( i_cell[idim] < params->number_of_patches[idim]-1 ) {
                    i_cell[idim]++;
                    x_cell[idim] += params->patch_dimensions[idim];
                    break;
                } else if( idim < params->nDim_field-1 ) {
                    i_cell[idim] = 0;
                    x_cell[idim] = params->patch_dimensions[idim] * 0.5;
                }
            }
        }
        
        return npart_total;
    } else {
        return 0. ;  //Do not account for particles initialized from numpy array
    }
}
