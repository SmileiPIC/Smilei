#include "PyTools.h"
#include <iomanip>

#include "DiagnosticScreen.h"
#include "HistogramFactory.h"


using namespace std;


// Constructor
DiagnosticScreen::DiagnosticScreen(
    Params &params,
    SmileiMPI *smpi,
    Patch *patch,
    int diagId
) :
    DiagnosticParticleBinningBase( params, smpi, patch, diagId, "Screen", true, nullptr, excludedAxes(diagId) )
{
    dt = params.timestep;
    
    ostringstream name( "" );
    name << "DiagScreen #" << diagId;
    string errorPrefix = name.str();
    
    // get parameter "shape" that determines if the screen is a plane or a sphere
    PyTools::extract( "shape", screen_shape, "DiagScreen", diagId );
    if( screen_shape == "plane" ) {
        screen_type = 0;
    } else if( screen_shape == "sphere" ) {
        screen_type = 1;
    } else {
        ERROR( errorPrefix << ": parameter `shape` must be 'plane' or 'sphere'" );
    }
    
    // get parameter "point" that determines the plane reference point, or the sphere center
    if( !PyTools::extractV( "point", screen_point, "DiagScreen", diagId ) ) {
        ERROR( errorPrefix << ": parameter `point` required" );
    }
    if( screen_point.size() != params.nDim_particle ) {
        ERROR( errorPrefix << ": parameter `point` must have "<<params.nDim_particle<<" elements in "<<params.nDim_particle<<"D" );
    }
    
    // get parameter "vector" that determines the plane normal, or the sphere radius
    if( !PyTools::extractV( "vector", screen_vector, "DiagScreen", diagId ) ) {
        if( params.nDim_particle == 1 ) {
            screen_vector.resize( 1, 1. );
        } else {
            ERROR( errorPrefix << ": parameter `vector` required" );
        }
    }
    if( screen_vector.size() != params.nDim_particle ) {
        ERROR( errorPrefix << ": parameter `vector` must have "<<params.nDim_particle<<" elements in "<<params.nDim_particle<<"D" );
    }
    // Calculate the unit vector
    screen_unitvector.resize( params.nDim_particle, 0. );
    screen_vectornorm = 0.;
    for( unsigned int i=0; i<params.nDim_particle; i++ ) {
        screen_vectornorm += screen_vector[i]*screen_vector[i];
    }
    if( screen_vectornorm == 0. ) {
        ERROR( errorPrefix << ": parameter `vector` must not be a null vector" );
    }
    screen_vectornorm = sqrt( screen_vectornorm );
    for( unsigned int i=0; i<params.nDim_particle; i++ ) {
        screen_unitvector[i] = screen_vector[i]/screen_vectornorm;
    }
    // Calculate other unit vectors to form an orthogonal base
    screen_vector_a.resize( params.nDim_particle, 0. );
    screen_vector_b.resize( params.nDim_particle, 0. );
    if( params.nDim_particle > 1 ) {
        screen_vector_a[0] = -screen_unitvector[1];
        screen_vector_a[1] =  screen_unitvector[0];
        double norm = sqrt( pow( screen_vector_a[0], 2 ) + pow( screen_vector_a[1], 2 ) );
        if( norm < 1.e-8 ) {
            screen_vector_a[0] = 0.;
            screen_vector_a[1] = 1.;
        } else {
            screen_vector_a[0] /= norm;
            screen_vector_a[1] /= norm;
        }
    }
    if( params.nDim_particle > 2 ) {
        screen_vector_b[0] = screen_unitvector[1]*screen_vector_a[2] - screen_unitvector[2]*screen_vector_a[1];
        screen_vector_b[1] = screen_unitvector[2]*screen_vector_a[0] - screen_unitvector[0]*screen_vector_a[2];
        screen_vector_b[2] = screen_unitvector[0]*screen_vector_a[1] - screen_unitvector[1]*screen_vector_a[0];
    }
    
    // get parameter "oriented", true if particles coming from the other side count negatively
    PyTools::extract( "direction", direction, "DiagScreen", diagId );
    if( direction=="both" ) {
        direction_type = 0;
    } else if( direction=="canceling" ) {
        direction_type = 1;
    } else if( direction=="forward" ) {
        direction_type = 2;
    } else if( direction=="backward" ) {
        direction_type = 3;
    } else {
        ERROR( errorPrefix << ": parameter `direction` not understood" );
    }
    
    // If axes are "a", "b", "theta" or "phi", they need some coefficients
    unsigned int idim;
    for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
        string type = histogram->axes[i]->type;
        if( type=="a" || type=="b" || type=="theta" || type=="phi" ) {
            vector<double> coefficients( params.nDim_particle * 2 );
            for( idim=0; idim<params.nDim_particle; idim++ ) {
                coefficients[idim] = screen_point[idim];
            }
            if( type == "a" ) {
                for( idim=0; idim<params.nDim_particle; idim++ ) {
                    coefficients[params.nDim_particle+idim] = screen_vector_a[idim];
                }
            } else if( type == "b" ) {
                for( idim=0; idim<params.nDim_particle; idim++ ) {
                    coefficients[params.nDim_particle+idim] = screen_vector_b[idim];
                }
            } else if( type == "theta" ) {
                for( idim=0; idim<params.nDim_particle; idim++ ) {
                    coefficients[params.nDim_particle+idim] = screen_vector[idim] / pow( screen_vectornorm, 2 );
                }
            } else if( type == "phi" ) {
                coefficients.resize( params.nDim_particle * 3 );
                for( idim=0; idim<params.nDim_particle; idim++ ) {
                    coefficients[params.nDim_particle+idim] = screen_vector_a[idim];
                }
                for( idim=0; idim<params.nDim_particle; idim++ ) {
                    coefficients[2*params.nDim_particle+idim] = screen_vector_b[idim];
                }
            }
            histogram->axes[i]->coefficients = coefficients;
        }
    }
    
    data_sum.resize( output_size, 0. );
    
} // END DiagnosticScreen::DiagnosticScreen


DiagnosticScreen::~DiagnosticScreen()
{
} // END DiagnosticScreen::~DiagnosticScreen


bool DiagnosticScreen::prepare( int timestep )
{

    // This diag always runs, but the output is not done at every timestep
    return true;
    
} // END prepare


// run one screen diagnostic
void DiagnosticScreen::run( Patch *patch, int timestep, SimWindow *simWindow )
{

    vector<int> int_buffer;
    vector<double> double_buffer;
    vector<bool> opposite;
    unsigned int npart, ndim = screen_point.size(), ipart, idim, nuseful;
    double side, side_old, dtg;
    
    // Verify that this patch is in a useful region for this diag
    if( screen_type == 0 ) { // plane
        double distance_to_plane = 0.;
        for( idim=0; idim<ndim; idim++ ) {
            distance_to_plane += ( patch->center[idim] - screen_point[idim] ) * screen_unitvector[idim];
        }
        if( abs( distance_to_plane ) > patch->radius ) {
            return;
        }
    } else { // sphere
        double distance_to_center = 0.;
        for( idim=0; idim<ndim; idim++ ) {
            distance_to_center += pow( patch->center[idim] - screen_point[idim], 2 );
        }
        distance_to_center = sqrt( distance_to_center );
        if( abs( screen_vectornorm - distance_to_center ) > patch->radius ) {
            return;
        }
    }
    
    // loop species
    for( unsigned int ispec=0 ; ispec < species.size() ; ispec++ ) {
    
        Species *s = patch->vecSpecies[species[ispec]];
        npart = s->particles->size();
        nuseful = 0;
        int_buffer   .resize( npart );
        double_buffer.resize( npart );
        opposite     .resize( npart, false );
        
        // Fill the int_buffer with -1 (not crossing screen) and 0 (crossing screen)
        if( screen_type == 0 ) { // plane
            for( ipart=0; ipart<npart; ipart++ ) {
                side = 0.;
                side_old = 0.;
                dtg = dt / s->particles->LorentzFactor( ipart );
                for( idim=0; idim<ndim; idim++ ) {
                    side += ( s->particles->Position[idim][ipart] - screen_point[idim] ) * screen_unitvector[idim];
                    side_old += ( s->particles->Position[idim][ipart] - dtg*( s->particles->Momentum[idim][ipart] ) - screen_point[idim] ) * screen_unitvector[idim];
                }
                if( side*side_old < 0. ) {
                    int_buffer[ipart] = 0;
                    nuseful++;
                    if( side < 0. ) {
                        opposite[ipart] = true;
                    }
                } else {
                    int_buffer[ipart] = -1;
                }
            }
        } else { // sphere
            for( ipart=0; ipart<npart; ipart++ ) {
                side = 0.;
                side_old = 0.;
                dtg = dt / s->particles->LorentzFactor( ipart );
                for( idim=0; idim<ndim; idim++ ) {
                    side += pow( s->particles->Position[idim][ipart] - screen_point[idim], 2 );
                    side_old += pow( s->particles->Position[idim][ipart] - dtg*( s->particles->Momentum[idim][ipart] ) - screen_point[idim], 2 );
                }
                side     = screen_vectornorm-sqrt( side );
                side_old = screen_vectornorm-sqrt( side_old );
                if( side*side_old < 0. ) {
                    int_buffer[ipart] = 0;
                    nuseful++;
                    if( side > 0. ) {
                        opposite[ipart] = true;
                    }
                } else {
                    int_buffer[ipart] = -1;
                }
            }
        }
        
        if( nuseful == 0 ) {
            continue;
        }
        
        histogram->digitize( s, double_buffer, int_buffer, simWindow );
        histogram->valuate( s, double_buffer, int_buffer );
        
        if( direction_type == 1 ) { // canceling
            for( ipart=0; ipart<npart; ipart++ )
                if( opposite[ipart] ) {
                    double_buffer[ipart] = -double_buffer[ipart];
                }
        } else if( direction_type == 2 ) { // forward
            for( ipart=0; ipart<npart; ipart++ )
                if( opposite[ipart] ) {
                    double_buffer[ipart] = 0.;
                }
        } else if( direction_type == 3 ) { // backward
            for( ipart=0; ipart<npart; ipart++ )
                if( int_buffer[ipart]>=0 && !opposite[ipart] ) {
                    double_buffer[ipart] = 0.;
                }
        }
        
        histogram->distribute( double_buffer, int_buffer, data_sum );
        
    }
    
} // END run

bool DiagnosticScreen::writeNow( int timestep ) {
    return timeSelection->theTimeIsNow( timestep );
}

//! Zero the array
void DiagnosticScreen::clear()
{
    fill( data_sum.begin(), data_sum.end(), 0. );
}
