#include "IonizationFromRate.h"

#include <cmath>

#include "Particles.h"
#include "ParticleData.h"
#include "Species.h"

using namespace std;



IonizationFromRate::IonizationFromRate( Params &params, Species *species ) : Ionization( params, species )
{

    DEBUG( "Creating the FromRate Ionizaton class" );
    
    maximum_charge_state_ = species->maximum_charge_state_;
    ionization_rate_ = species->ionization_rate_;
    
    DEBUG( "Finished Creating the FromRate Ionizaton class" );
    
}



void IonizationFromRate::operator()( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, vector<double> *, Patch *patch, Projector *, int )
{

    //unsigned int Z, Zp1, newZ, k_times;
    unsigned int Z, k_times;
    vector<double> rate;
    
    // Leave if nothing to do
    if( ipart_min >= ipart_max ) {
        return;
    }
    
#ifdef SMILEI_USE_NUMPY
    // Run python to evaluate the ionization rate for each particle
    unsigned int npart = ipart_max - ipart_min;
    SMILEI_PY_ACQUIRE_GIL
    {
        ParticleData particleData( npart );
        particleData.startAt( ipart_min );
        particleData.set( particles );
        PyArrayObject *ret = ( PyArrayObject * )PyObject_CallFunctionObjArgs( ionization_rate_, particleData.get(), NULL );
        PyTools::checkPyError();
        if( ret == NULL ) {
            ERROR( "ionization_rate profile has not provided a correct result" );
        }
        double *arr = ( double * ) PyArray_GETPTR1( ret, 0 );
        rate.resize( npart );
        // Loop the return value and store
        for( unsigned int i=0; i<npart; i++ ) {
            rate[i] = arr[i];
        }
        Py_DECREF( ret );
    }
    SMILEI_PY_RELEASE_GIL
#endif
    
    
    for( unsigned int ipart=ipart_min ; ipart<ipart_max; ipart++ ) {
    
        // Current charge state of the ion
        Z = ( unsigned int )( particles->charge( ipart ) );
        
        // If ion already fully ionized then skip
        if( Z==maximum_charge_state_ ) {
            continue;
        }
        
        // Start of the Monte-Carlo routine  (At the moment, only 1 ionization per timestep is possible)
        // k_times will give the nb of ionization events
        k_times = 0;
        double ran_p = patch->rand_->uniform();
        if( ran_p < 1.0 - exp( -rate[ipart-ipart_min]*dt ) ) {
            k_times        = 1;
        }
        
        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------
        if( k_times!=0 ) {
            new_electrons.createParticle();
            int idNew = new_electrons.size() - 1;
            for( unsigned int i=0; i<new_electrons.dimension(); i++ ) {
                new_electrons.position( i, idNew )=particles->position( i, ipart );
            }
            for( unsigned int i=0; i<3; i++ ) {
                new_electrons.momentum( i, idNew ) = particles->momentum( i, ipart )*ionized_species_invmass;
            }
            new_electrons.weight( idNew )=double( k_times )*particles->weight( ipart );
            new_electrons.charge( idNew )=-1;
            
            if( save_ion_charge_ ) {
                ion_charge_.push_back( particles->charge( ipart ) );
            }
            
            // Increase the charge of the particle
            particles->charge( ipart ) += k_times;
        }
        
        
    } // Loop on particles
}
