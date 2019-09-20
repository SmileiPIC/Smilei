
#include "Species.h"

#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>

#include <omp.h>

// IDRIS
#include <cstring>
// IDRIS
#include "PusherFactory.h"
#include "IonizationFactory.h"
#include "RadiationFactory.h"
#include "MultiphotonBreitWheelerFactory.h"
#include "MergingFactory.h"
#include "PartBoundCond.h"
#include "PartWall.h"
#include "BoundaryConditionType.h"

#include "ElectroMagn.h"
#include "Interpolator.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "Profile.h"
#include "ElectroMagnAM.h"
#include "Projector.h"
#include "ProjectorFactory.h"
#include "ParticleCreator.h"

#include "SimWindow.h"
#include "Patch.h"

// #include "Field.h"
#include "Field1D.h"
#include "Field2D.h"
#include "Field3D.h"
#include "Tools.h"

#include "DiagnosticTrack.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Species
// input: simulation parameters & Species index
// ---------------------------------------------------------------------------------------------------------------------
Species::Species( Params &params, Patch *patch ) :
    c_part_max_( 1 ),
    ionization_rate( Py_None ),
    pusher( "boris" ),
    radiation_model( "none" ),
    time_frozen( 0 ),
    radiating( false ),
    relativistic_field_initialization( false ),
    time_relativistic_initialization( 0 ),
    multiphoton_Breit_Wheeler( 2, "" ),
    ionization_model( "none" ),
    density_profile_type_( "none" ),
    chargeProfile( NULL ),
    density_profile_( NULL ),
    velocity_profile_( 3, NULL ),
    temperature_profile_( 3, NULL ),
    particles_per_cell_profile_( NULL ),
    max_charge( 0. ),
    particles( &particles_sorted[0] ),
    regular_number_array_(0),
    position_initialization_array_( NULL ),
    momentum_initialization_array_( NULL ),
    n_numpy_particles_( 0 ),
    position_initialization_on_species_( false ),
    position_initialization_on_species_index( -1 ),
    electron_species( NULL ),
    electron_species_index( -1 ),
    photon_species( NULL ),
    //photon_species_index(-1),
    radiation_photon_species( "" ),
    mBW_pair_creation_sampling( 2, 1 ),
    clrw( params.clrw ),
    oversize( params.oversize ),
    cell_length( params.cell_length ),
    min_loc_vec( patch->getDomainLocalMin() ),
    tracking_diagnostic( 10000 ),
    nDim_particle( params.nDim_particle ),
    partBoundCond( NULL ),
    min_loc( patch->getDomainLocalMin( 0 ) ),
    merging_method_( "none" ),
    merging_time_selection_( 0 )
{

    PI2 = 2.0 * M_PI;
    PI_ov_2 = 0.5*M_PI;

    dx_inv_[0] = 1./cell_length[0];
    dx_inv_[1] = 1./cell_length[1];
    dx_inv_[2] = 1./cell_length[2];

    initCluster( params );
    nDim_field = params.nDim_field;
    inv_nDim_particles = 1./( ( double )nDim_particle );

    length_[0]=0;
    length_[1]=params.n_space[1]+1;
    length_[2]=params.n_space[2]+1;

    merge_momentum_cell_size_.resize(3);

    merge_min_momentum_cell_length_.resize(3);

}//END Species creator

void Species::initCluster( Params &params )
{
    // Arrays of the min and max indices of the particle bins
    first_index.resize( params.n_space[0]/clrw );
    last_index.resize( params.n_space[0]/clrw );

    //Size in each dimension of the buffers on which each bin are projected
    //In 1D the particles of a given bin can be projected on 6 different nodes at the second order (oversize = 2)

    //Primal dimension of fields.
    f_dim0 =  params.n_space[0] + 2 * oversize[0] +1;
    f_dim1 =  params.n_space[1] + 2 * oversize[1] +1;
    f_dim2 =  params.n_space[2] + 2 * oversize[2] +1;

    b_dim.resize( params.nDim_field, 1 );
    if( nDim_particle == 1 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0];
        f_dim1 = 1;
        f_dim2 = 1;
    }
    if( nDim_particle == 2 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] =  f_dim1;
        f_dim2 = 1;
    }
    if( nDim_particle == 3 ) {
        b_dim[0] = ( 1 + clrw ) + 2 * oversize[0]; // There is a primal number of bins.
        b_dim[1] = f_dim1;
        b_dim[2] = f_dim2;
    }

    //Initialize specMPI
    MPIbuff.allocate( nDim_particle );

    //ener_tot = 0.;
    nrj_bc_lost = 0.;
    nrj_mw_lost = 0.;
    new_particles_energy_ = 0.;
    nrj_radiation = 0.;

}//END initCluster

// -----------------------------------------------------------------------------
//! This function enables to resize the number of bins
// -----------------------------------------------------------------------------
void Species::resizeCluster( Params &params )
{

    // We keep the current number of particles
    int npart = particles->size();
    int size = params.n_space[0]/clrw;

    // Arrays of the min and max indices of the particle bins
    first_index.resize( size );
    last_index.resize( size );

    // We redistribute the particles between the bins
    int quotient = npart / size; // Fixed part for all bin
    int remainder = npart - quotient*size; // To be distributed among the first bin

    for( int ibin=0 ; ibin < size ; ibin++ ) {
        if( ibin < remainder ) {
            first_index[ibin] = ibin*quotient + ibin;
            last_index[ibin] = first_index[ibin] + quotient + 1;
        } else {
            first_index[ibin] = ibin*quotient + remainder;
            last_index[ibin] = first_index[ibin] + quotient;
        }
    }

    //std::cout << "size: " << size << " " << npart << " " << first_index[0] << " " << last_index[0] << '\n';

    // Recommended: A sorting process may be needed for best porfermance after this step

}// end resizeCluster

// Initialize the operators (Push, Ionize, PartBoundCond)
// This must be separate from the parameters because the Species cloning copies
// the parameters but not the operators.
void Species::initOperators( Params &params, Patch *patch )
{

    // interpolation operator (virtual)
    Interp = InterpolatorFactory::create( params, patch, this->vectorized_operators && !params.cell_sorting ); // + patchId -> idx_domain_begin (now = ref smpi)
    
    // assign the correct Pusher to Push
    Push = PusherFactory::create( params, this );
    if( this->ponderomotive_dynamics ) {
        Push_ponderomotive_position = PusherFactory::create_ponderomotive_position_updater( params, this );
    }

    // projection operator (virtual)
    Proj = ProjectorFactory::create( params, patch, this->vectorized_operators && !params.cell_sorting );  // + patchId -> idx_domain_begin (now = ref smpi)
    
    // Assign the Ionization model (if needed) to Ionize
    //  Needs to be placed after ParticleCreator() because requires the knowledge of max_charge
    // \todo pay attention to restart
    Ionize = IonizationFactory::create( params, this );

    // Create the radiation model
    Radiate = RadiationFactory::create( params, this );

    // Create the multiphoton Breit-Wheeler model
    Multiphoton_Breit_Wheeler_process = MultiphotonBreitWheelerFactory::create( params, this );

    // assign the correct Merging method to Merge
    Merge = MergingFactory::create( params, this );

    // define limits for BC and functions applied and for domain decomposition
    partBoundCond = new PartBoundCond( params, this, patch );
    for( unsigned int iDim=0 ; iDim < nDim_particle ; iDim++ ) {
        for( unsigned int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++ ) {
            MPIbuff.partRecv[iDim][iNeighbor].initialize( 0, ( *particles ) );
            MPIbuff.partSend[iDim][iNeighbor].initialize( 0, ( *particles ) );
            MPIbuff.part_index_send[iDim][iNeighbor].resize( 0 );
            MPIbuff.part_index_recv_sz[iDim][iNeighbor] = 0;
            MPIbuff.part_index_send_sz[iDim][iNeighbor] = 0;
        }
    }
    typePartSend.resize( nDim_particle*2, MPI_DATATYPE_NULL );
    typePartRecv.resize( nDim_particle*2, MPI_DATATYPE_NULL );
    exchangePatch = MPI_DATATYPE_NULL;

}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species
// ---------------------------------------------------------------------------------------------------------------------
Species::~Species()
{
    delete Push;
    delete Interp;
    delete Proj;

    if( Merge ) {
        delete Merge;
    }

    if( Ionize ) {
        delete Ionize;
    }
    if( Radiate ) {
        delete Radiate;
    }
    if( Multiphoton_Breit_Wheeler_process ) {
        delete Multiphoton_Breit_Wheeler_process;
    }
    if( partBoundCond ) {
        delete partBoundCond;
    }
    if( particles_per_cell_profile_ ) {
        delete particles_per_cell_profile_;
    }
    if( chargeProfile ) {
        delete chargeProfile;
    }
    if( density_profile_ ) {
        delete density_profile_;
    }
    for( unsigned int i=0; i<velocity_profile_.size(); i++ ) {
        delete velocity_profile_[i];
    }
    for( unsigned int i=0; i<temperature_profile_.size(); i++ ) {
        delete temperature_profile_[i];
    }
    if( ionization_rate!=Py_None ) {
        Py_DECREF( ionization_rate );
    }

}



// ---------------------------------------------------------------------------------------------------------------------
// For all (np) particles in a mesh initialize their position
//   - either using regular distribution in the mesh (position_initialization_ = regular)
//   - or using uniform random distribution (position_initialization_ = random)
// ---------------------------------------------------------------------------------------------------------------------
void Species::initPosition( unsigned int nPart, unsigned int iPart, double *indexes, Params &params )
{
    if( position_initialization_ == "regular" ) {

        

        if( params.geometry != "AMcylindrical" ) {
            double inv_coeff_array[3];
            int    coeff_array[3];

            if ( regular_number_array_.size()==0){
                double coeff = pow( ( double )nPart, inv_nDim_particles );
                if( nPart != ( unsigned int ) pow( round( coeff ), ( double )nDim_particle ) ) {
                    ERROR( "Impossible to put "<<nPart<<" particles regularly spaced in one cell. Use a square number, or `position_initialization = 'random'`" );
                }
                for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
                    coeff_array[idim] = coeff;
                    inv_coeff_array[idim] = 1./coeff;
                }
            } else{
               int npart_check=1;
               for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
                   npart_check *= regular_number_array_[idim];
               }
               if( nPart != npart_check) {
                   ERROR( "The number of particles required per cell and per dimension is not coherent with the total number of particles per cell." );
               }
               for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
                   coeff_array[idim] = regular_number_array_[idim];
                   inv_coeff_array[idim] = 1./(double)coeff_array[idim];
               }

            }

            for( unsigned int  p=iPart; p<iPart+nPart; p++ ) {
                int i = ( int )( p-iPart );
                for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
                    particles->position( idim, p ) = indexes[idim] + cell_length[idim] * 0.975 * inv_coeff_array[idim] * ( 0.5 + i%coeff_array[idim] );
                    i /= coeff_array[idim]; // integer division
                }
            }
        } else {

            //Trick to derive number of particles per dimension from total number of particles per cell
            unsigned int Np_array[nDim_particle];
            int Np = nPart;
            int counter = 0;
            unsigned int prime = 2;
            double dx, dr, dtheta, theta_offset;
            for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
                Np_array[idim] = 1;
            }

            while( prime <= 23 && Np > 1 ) {
                if( Np%prime == 0 ) {
                    Np = Np/prime;
                    Np_array[counter%nDim_particle] *= prime;
                    counter++;
                } else {
                    prime++;
                }
            }
            Np_array[counter%nDim_particle] *= Np; //At that point, if Np is not equal to 1, it means that nPart has a prime divisor greater than 23.
            std::sort( Np_array, Np_array + nDim_particle ); //sort so that the largest number of particles per dimension is used along theta.

            dx = cell_length[0]/Np_array[0];
            dr = cell_length[1]/Np_array[1];
            dtheta = 2.*M_PI   /Np_array[2];

            for( unsigned int ix = 0 ; ix < Np_array[0]; ix++ ) {
                double qx = indexes[0] + dx*( ix+0.5 );
                int nx = ix*( Np_array[2]*Np_array[1] );
                for( unsigned int ir = 0 ; ir < Np_array[1]; ir++ ) {
                    double qr = indexes[1] + dr*( ir+0.5 );
                    int nr = ir*( Np_array[2] );
                    theta_offset = Rand::uniform()*2.*M_PI;
                    for( unsigned int itheta = 0 ; itheta < Np_array[2]; itheta++ ) {
                        int p = nx+nr+itheta+iPart;
                        double theta = theta_offset + itheta*dtheta;
                        particles->position( 0, p ) = qx ;
                        particles->position( 1, p ) = qr*cos( theta );
                        particles->position( 2, p ) = qr*sin( theta );
                    }
                }
            }
        }

    } else if( position_initialization_ == "random" ) {
        if( params.geometry=="AMcylindrical" ) {
            double particles_r, particles_theta;
            for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
                particles->position( 0, p )=indexes[0]+Rand::uniform()*cell_length[0];
                particles_r=sqrt( indexes[1]*indexes[1]+ 2.*Rand::uniform()*( indexes[1]+cell_length[1]*0.5 )*cell_length[1] );
                particles_theta=Rand::uniform()*2.*M_PI;
                particles->position( 2, p )=particles_r*sin( particles_theta );
                particles->position( 1, p )= particles_r*cos( particles_theta );
            }
        } else {
            for( unsigned int p= iPart; p<iPart+nPart; p++ ) {
                for( unsigned int i=0; i<nDim_particle ; i++ ) {
                    particles->position( i, p )=indexes[i]+Rand::uniform()*cell_length[i];
                }
            }
        }
    } else if( position_initialization_ == "centered" ) {

        for( unsigned int p=iPart; p<iPart+nPart; p++ )
            for( unsigned int i=0; i<nDim_particle ; i++ ) {
                particles->position( i, p )=indexes[i]+0.5*cell_length[i];
            }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - perform ionization
//   - perform the radiation reaction
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::dynamics( double time_dual, unsigned int ispec,
                        ElectroMagn *EMfields,
                        Params &params, bool diag_flag,
                        PartWalls *partWalls,
                        Patch *patch, SmileiMPI *smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                        vector<Diagnostic *> &localDiags )
{
    int ithread, tid( 0 );
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    unsigned int iPart;

    // Reset list of particles to exchange
    clearExchList();

    double ener_iPart( 0. );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle dynamics
    // -------------------------------
    if( time_dual>time_frozen || Ionize) { // moving particle
    
        smpi->dynamics_resize( ithread, nDim_field, last_index.back(), params.geometry=="AMcylindrical" );
        //Point to local thread dedicated buffers
        //Still needed for ionization
        vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );

        for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin++ ) {

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Interpolate the fields at the particle position
            Interp->fieldsWrapper( EMfields, *particles, smpi, &( first_index[ibin] ), &( last_index[ibin] ), ithread );

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[0] += MPI_Wtime() - timer;
#endif

            // Ionization
            if( Ionize ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                ( *Ionize )( particles, first_index[ibin], last_index[ibin], Epart, patch, Proj );

#ifdef  __DETAILED_TIMERS
                patch->patch_timers[4] += MPI_Wtime() - timer;
#endif
            }
            
            if( time_dual<=time_frozen ) continue; // Do not push nor project frozen particles

            // Radiation losses
            if( Radiate ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Radiation process
                ( *Radiate )( *particles, this->photon_species, smpi,
                              RadiationTables,
                              first_index[ibin], last_index[ibin], ithread );

                // Update scalar variable for diagnostics
                nrj_radiation += Radiate->getRadiatedEnergy();

                // Update the quantum parameter chi
                Radiate->computeParticlesChi( *particles,
                                              smpi,
                                              first_index[ibin],
                                              last_index[ibin],
                                              ithread );
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[5] += MPI_Wtime() - timer;
#endif

            }


            // Multiphoton Breit-Wheeler
            if( Multiphoton_Breit_Wheeler_process ) {

#ifdef  __DETAILED_TIMERS
                timer = MPI_Wtime();
#endif

                // Pair generation process
                ( *Multiphoton_Breit_Wheeler_process )( *particles,
                                                        smpi,
                                                        MultiphotonBreitWheelerTables,
                                                        first_index[ibin], last_index[ibin], ithread );

                // Update scalar variable for diagnostics
                // We reuse nrj_radiation for the pairs
                nrj_radiation += Multiphoton_Breit_Wheeler_process->getPairEnergy();

                // Update the photon quantum parameter chi of all photons
                Multiphoton_Breit_Wheeler_process->compute_thread_chiph( *particles,
                        smpi,
                        first_index[ibin],
                        last_index[ibin],
                        ithread );
                
                // Suppression of the decayed photons into pairs
                Multiphoton_Breit_Wheeler_process->decayed_photon_cleaning(
                    *particles, smpi, ibin, first_index.size(), &first_index[0], &last_index[0], ithread );
                    
#ifdef  __DETAILED_TIMERS
                patch->patch_timers[6] += MPI_Wtime() - timer;
#endif

            }

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Push the particles and the photons
            ( *Push )( *particles, smpi, first_index[ibin], last_index[ibin], ithread );
            //particles->test_move( first_index[ibin], last_index[ibin], params );

        } //ibin

        // Particle injection
        // if (0) {
        //     vector<unsigned int> init_space( 3, 1 );
        //     int previous_particle_number = getNbrOfParticles();
        //     init_space[0] = 1;
        //     init_space[1] = params.n_space[1];
        //     init_space[2] = params.n_space[2];
        //     if ( ( (patch->isXmin()) && (species_number_<2) ) ||
        //          ( (patch->isXmax()) && (species_number_>1) ) ) {
        //         int icell;
        //         if (patch->isXmin())
        //             icell = 0;
        //         if (patch->isXmax())
        //             icell = params.n_space[0]-1;
        //         this->ParticleCreator( init_space, params, patch, icell );
        //     }
        //     // Suppr not interesting parts ...
        //     int new_particle_number = getNbrOfParticles();
        //     // cerr << "Number of new particles: " << new_particle_number-previous_particle_number
        //     //      << endl;
        //     for ( int ip = new_particle_number-1 ; ip >=  previous_particle_number; ip-- ){
        //         if ( patch->isXmin() ) {
        //             particles->Position[0][ip] += ( params.timestep*particles->Momentum[0][ip]*
        //                                             particles->inv_lor_fac(ip)-params.cell_length[0]);
        //             if ( ( particles->Position[0][ip] < 0) ) {
        //                 particles->erase_particle(ip);
        //                 new_particle_number--;
        //             }
        //         } else if ( patch->isXmax() ) {
        //             particles->Position[0][ip] += ( params.timestep*particles->Momentum[0][ip]*
        //                                             particles->inv_lor_fac(ip)+params.cell_length[0] );
        //             if (  particles->Position[0][ip] > params.grid_length[0]) {
        //                 particles->erase_particle(ip);
        //                 new_particle_number--;
        //             }
        //         }
        //     }
        //     for ( int ip = new_particle_number-1 ; ip >=  previous_particle_number; ip-- ){
        //         particles->Position[1][ip] += ( params.timestep*particles->Momentum[1][ip]*
        //                                         particles->inv_lor_fac(ip));
        //         if ( patch->isYmin() ) {
        //             particles->Position[0][ip] -= params.cell_length[1];
        //             if ( ( particles->Position[1][ip] < 0 ) ) {
        //                 particles->erase_particle(ip);
        //                 new_particle_number--;
        //             }
        //         } else if ( patch->isYmax() ) {
        //             particles->Position[0][ip] += params.cell_length[1] ;
        //             if (  particles->Position[1][ip] > params.grid_length[1] ) {
        //                 particles->erase_particle(ip);
        //                 new_particle_number--;
        //             }
        //         }
        //     }
        //
        //     double xpn,ypn,inv_gamma;
        //     int new_particle_index;
        //     int ioldx, ioldy;
        //     // smpi->dynamics_resize( ithread, nDim_field, new_particle_number, params.geometry=="AMcylindrical");
        //     // for ( int ip = previous_particle_number ; ip <  new_particle_number; ip++ ){
        //     //     inv_gamma = particles->inv_lor_fac(ip);
        //     //     xpn = (particles->Position[0][ip]
        //     //         - particles->Momentum[0][ip]*inv_gamma*params.timestep)*dx_inv_[0];
        //     //     ypn = (particles->Position[1][ip]
        //     //         - particles->Momentum[1][ip]*inv_gamma*params.timestep)*dx_inv_[1];
        //     //     smpi->dynamics_invgf[ithread][ip] = inv_gamma;
        //     //
        //     //     smpi->dynamics_iold[ithread][ip] = round(xpn);
        //     //     smpi->dynamics_iold[ithread][ip + new_particle_number] = round(ypn);
        //     //
        //     //     smpi->dynamics_deltaold[ithread][ip] = xpn - ( double )(smpi->dynamics_iold[ithread][ip]);
        //     //     smpi->dynamics_deltaold[ithread][ip + new_particle_number] = ypn - ( double )(smpi->dynamics_iold[ithread][ip + new_particle_number]);
        //     //
        //     //     smpi->dynamics_iold[ithread][ip] -= patch->getCellStartingGlobalIndex( 0 );
        //     //     smpi->dynamics_iold[ithread][ip + new_particle_number] -= patch->getCellStartingGlobalIndex( 1 );
        //     //
        //     // }
        //
        //
        //     // cerr << " clrw: " << params.clrw
        //     //      << " params.grid_length[0]: " << params.grid_length[0]
        //     //      << endl;
        //     // Move interesting parts to their place
        //     for ( int ip = previous_particle_number ; ip < new_particle_number ; ip++ ){
        //         // cerr << " Particle x: " << particles->Position[0][ip]
        //         //      << " Particle y: " << particles->Position[1][ip]
        //         //      << " Ip: " << ip
        //         //      << endl;
        //
        //         inv_gamma = particles->inv_lor_fac(ip);
        //
        //         xpn = (particles->Position[0][ip]
        //             - particles->Momentum[0][ip]*inv_gamma*params.timestep)*dx_inv_[0];
        //         ypn = (particles->Position[1][ip]
        //             - particles->Momentum[1][ip]*inv_gamma*params.timestep)*dx_inv_[1];
        //
        //         if ( patch->isXmin() ) {
        //             if ( ( particles->Position[0][ip] >= 0. ) ) {
        //                 int new_cell_idx=0;
        //
        //                 new_particle_index = first_index[(new_cell_idx)/params.clrw];
        //
        //                 particles->mv_particles(ip,new_particle_index);
        //
        //                 smpi->dynamics_invgf[ithread].insert(smpi->dynamics_invgf[ithread].begin()+new_particle_index,inv_gamma);
        //
        //                 ioldx = round(xpn);
        //                 ioldy = round(ypn);
        //
        //                 smpi->dynamics_iold[ithread].insert(smpi->dynamics_iold[ithread].begin()
        //                                                 + new_particle_index,ioldx - patch->getCellStartingGlobalIndex( 0 ));
        //                 smpi->dynamics_iold[ithread].insert(smpi->dynamics_iold[ithread].begin()
        //                                                 + new_particle_index + new_particle_number,ioldy - patch->getCellStartingGlobalIndex( 1 ));
        //
        //                 smpi->dynamics_deltaold[ithread].insert(smpi->dynamics_deltaold[ithread].begin()
        //                                                 + new_particle_index,xpn - (double)ioldx);
        //                 smpi->dynamics_deltaold[ithread].insert(smpi->dynamics_deltaold[ithread].begin()
        //                                                 + new_particle_index + new_particle_number,ypn - (double)ioldy);
        //
        //                 last_index[(new_cell_idx)/params.clrw]++;
        //                 for ( int idx=(new_cell_idx)/params.clrw+1 ; idx<last_index.size() ; idx++ ) {
        //                     first_index[idx]++;
        //                     last_index[idx]++;
        //                 }
        //             }
        //         }
        //         else if ( patch->isXmax() ) {
        //             if (  particles->Position[0][ip] <= params.grid_length[0] ) {
        //                 int new_cell_idx=params.n_space[0]-1;
        //
        //                 new_particle_index = first_index[(new_cell_idx)/params.clrw];
        //
        //                 particles->mv_particles(ip,first_index[(new_cell_idx)/params.clrw]);
        //
        //                 smpi->dynamics_invgf[ithread].insert(smpi->dynamics_invgf[ithread].begin()+new_particle_index,inv_gamma);
        //
        //                 ioldx = round(xpn);
        //                 ioldy = round(ypn);
        //
        //                 smpi->dynamics_iold[ithread].insert(smpi->dynamics_iold[ithread].begin()
        //                                                 + new_particle_index,ioldx - patch->getCellStartingGlobalIndex( 0 ));
        //                 smpi->dynamics_iold[ithread].insert(smpi->dynamics_iold[ithread].begin()
        //                                                 + new_particle_index + new_particle_number,ioldy - patch->getCellStartingGlobalIndex( 1 ));
        //
        //                 smpi->dynamics_deltaold[ithread].insert(smpi->dynamics_deltaold[ithread].begin()
        //                                                 + new_particle_index,xpn - (double)ioldx);
        //                 smpi->dynamics_deltaold[ithread].insert(smpi->dynamics_deltaold[ithread].begin()
        //                                                 + new_particle_index + new_particle_number,ypn - (double)ioldy);
        //
        //                 last_index[(new_cell_idx)/params.clrw]++;
        //                 for ( int idx=(new_cell_idx)/params.clrw+1 ; idx<last_index.size() ; idx++ ) {
        //                     first_index[idx]++;
        //                     last_index[idx]++;
        //                 }
        //             }
        //         }
        //     }
        // }
        
        
        
        for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin++ ) {

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[1] += MPI_Wtime() - timer;
            timer = MPI_Wtime();
#endif

            // Apply wall and boundary conditions
            if( mass>0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                            nrj_lost_per_thd[tid] += mass * ener_iPart;
                        }
                    }
                }
                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                    if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += mass * ener_iPart;
                        //}
                        //else if ( partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        //std::cout<<"removed particle position"<< particles->position(0,iPart)<<" , "<<particles->position(1,iPart)<<" ,"<<particles->position(2,iPart)<<std::endl;
                    }
                }

            } else if( mass==0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                            nrj_lost_per_thd[tid] += ener_iPart;
                        }
                    }
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                    if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += ener_iPart;
                    }
                }

            }

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[3] += MPI_Wtime() - timer;
#endif

            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif

            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
            if( ( !particles->is_test ) && ( mass > 0 ) ) {
                Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, first_index[ibin], last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
            }

#ifdef  __DETAILED_TIMERS
            patch->patch_timers[2] += MPI_Wtime() - timer;
#endif

        }// ibin


        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

//        // Add the ionized electrons to the electron species
//        if (Ionize)
//            electron_species->importParticles( params, patch, Ionize->new_electrons, localDiags );
//
//        // Radiation losses
//        if (Radiate)
//        {
//            // If creation of macro-photon, we add them to photon_species
//            if (photon_species)
//            {
//                photon_species->importParticles(params,
//                                                patch,
//                                                Radiate->new_photons_,
//                                                localDiags);
//            }
//        }
//
//        // Multiphoton Breit-Wheeler
//        if (Multiphoton_Breit_Wheeler_process)
//        {
//
//            // Addition of the electron-positron particles
//            for (int k=0; k<2; k++) {
//                mBW_pair_species[k]->importParticles(params,
//                                             patch,
//                                             Multiphoton_Breit_Wheeler_process->new_pair[k],
//                                             localDiags);
//            }
//        }

    } //End if moving or ionized particles

    if(time_dual <= time_frozen && diag_flag &&( !particles->is_test ) ) { //immobile particle (at the moment only project density)
        if( params.geometry != "AMcylindrical" ) {
            double *b_rho=nullptr;
            for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                }
            }
        } else {
            int n_species = patch->vecSpecies.size();
            complex<double> *b_rho=nullptr;
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
                int ifield = imode*n_species+ispec;
                for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                    b_rho = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
                    for( int iPart=first_index[ibin] ; iPart<last_index[ibin]; iPart++ ) {
                        Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                    }
                }
            }
        }
    } // End projection for frozen particles
} //END dynamics


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - interpolate the fields at the particle position
//   - perform ionization
//   - perform the radiation reaction
//   - perform the multiphoton Breit-Wheeler
//   - calculate the new velocity
//   - calculate the new position
//   - apply the boundary conditions
//   - increment the currents (projection)
// ---------------------------------------------------------------------------------------------------------------------
void Species::scalar_dynamics( double time_dual, unsigned int ispec,
                               ElectroMagn *EMfields,
                               Params &params, bool diag_flag,
                               PartWalls *partWalls,
                               Patch *patch, SmileiMPI *smpi,
                               RadiationTables &RadiationTables,
                               MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                               vector<Diagnostic *> &localDiags )
{

}

void Species::projection_for_diags( double time_dual, unsigned int ispec,
                                    ElectroMagn *EMfields,
                                    Params &params, bool diag_flag,
                                    Patch *patch, SmileiMPI *smpi )
{
    if( diag_flag &&( !particles->is_test ) ) {

        if( params.geometry != "AMcylindrical" ) {
            double *buf[4];

            for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj

                buf[0] = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                buf[1] = EMfields->Jx_s [ispec] ? &( *EMfields->Jx_s [ispec] )( 0 ) : &( *EMfields->Jx_ )( 0 ) ;
                buf[2] = EMfields->Jy_s [ispec] ? &( *EMfields->Jy_s [ispec] )( 0 ) : &( *EMfields->Jy_ )( 0 ) ;
                buf[3] = EMfields->Jz_s [ispec] ? &( *EMfields->Jz_s [ispec] )( 0 ) : &( *EMfields->Jz_ )( 0 ) ;

                for( int iPart=first_index[ibin] ; iPart<last_index[ibin]; iPart++ ) {
                    for( unsigned int quantity=0; quantity < 4; quantity++ ) {
                        Proj->basic( buf[quantity], ( *particles ), iPart, quantity );
                    }
                } //End loop on particles
            }//End loop on bins
        } else { // AM case
            complex<double> *buf[4];
            ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
            int n_species = patch->vecSpecies.size();
            for( unsigned int imode = 0; imode<params.nmodes; imode++ ) {
                int ifield = imode*n_species+ispec;

                for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj

                    buf[0] = emAM->rho_AM_s[ifield] ? &( *emAM->rho_AM_s[ifield] )( 0 ) : &( *emAM->rho_AM_[imode] )( 0 ) ;
                    buf[1] = emAM->Jl_s [ifield] ? &( *emAM->Jl_s [ifield] )( 0 ) : &( *emAM->Jl_[imode] )( 0 ) ;
                    buf[2] = emAM->Jr_s [ifield] ? &( *emAM->Jr_s [ifield] )( 0 ) : &( *emAM->Jr_[imode] )( 0 ) ;
                    buf[3] = emAM->Jt_s [ifield] ? &( *emAM->Jt_s [ifield] )( 0 ) : &( *emAM->Jt_[imode] )( 0 ) ;

                    for( int iPart=first_index[ibin] ; iPart<last_index[ibin]; iPart++ ) {
                        for( unsigned int quantity=0; quantity < 4; quantity++ ) {
                            Proj->basicForComplex( buf[quantity], ( *particles ), iPart, quantity, imode );
                        }
                    } //End loop on particles
                }//End loop on bins
            } //End loop on modes
        }
        
    }
}

// -----------------------------------------------------------------------------
//! For all particles of the species, import the new particles generated
//! from these different physical processes:
//! - ionization
//! - radiation reaction
//! - multiphoton Breit-Wheeler
// -----------------------------------------------------------------------------
void Species::dynamics_import_particles( double time_dual, unsigned int ispec,
        Params &params,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags )
{
    // if moving particle
    if( time_dual>time_frozen ) { // moving particle

        // Add the ionized electrons to the electron species
        if( Ionize ) {
            electron_species->importParticles( params, patch, Ionize->new_electrons, localDiags );
        }

        // Radiation losses
        if( Radiate ) {
            // If creation of macro-photon, we add them to photon_species
            if( photon_species ) {
                photon_species->importParticles( params,
                                                 patch,
                                                 Radiate->new_photons_,
                                                 localDiags );
            }
        }

        // Multiphoton Breit-Wheeler
        if( Multiphoton_Breit_Wheeler_process ) {
            // Addition of the electron-positron particles
            for( int k=0; k<2; k++ ) {
                mBW_pair_species[k]->importParticles( params,
                                                      patch,
                                                      Multiphoton_Breit_Wheeler_process->new_pair[k],
                                                      localDiags );
            }
        }
    }//END if time vs. time_frozen
}




// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species
//   - increment the charge (projection)
//   - used at initialisation for Poisson (and diags if required, not for now dynamics )
// ---------------------------------------------------------------------------------------------------------------------
void Species::computeCharge( unsigned int ispec, ElectroMagn *EMfields )
{
    // -------------------------------
    // calculate the particle charge
    // -------------------------------
    if( ( !particles->is_test ) ) {
        for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
            // Not for now, else rho is incremented twice. Here and dynamics. Must add restartRhoJs and manage independantly diags output
            //b_rho = EMfields->rho_s[ispec] ? &(*EMfields->rho_s[ispec])(bin_start) : &(*EMfields->rho_)(bin_start);
            if( !dynamic_cast<ElectroMagnAM *>( EMfields ) ) {
                double *b_rho = &( *EMfields->rho_ )( 0 );

                for( unsigned int iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                }
            } else {
                ElectroMagnAM *emAM = static_cast<ElectroMagnAM *>( EMfields );
                unsigned int Nmode = emAM->rho_AM_.size();
                for( unsigned int imode=0; imode<Nmode; imode++ ) {
                    complex<double> *b_rho = &( *emAM->rho_AM_[imode] )( 0 );
                    for( unsigned int iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                        Proj->basicForComplex( b_rho, ( *particles ), iPart, 0, imode );
                    }
                }
            }
        }

    }
}//END computeCharge


// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::sort_part( Params &params )
{

    int ndim = params.nDim_field;
    int idim;
    //cleanup_sent_particles(ispec, indexes_of_particles_to_exchange);

    //We have stored in indexes_of_particles_to_exchange the list of all particles that needs to be removed.
    /********************************************************************************/
    // Delete Particles included in the index of particles to exchange. Assumes indexes are sorted.
    /********************************************************************************/
    int ii, iPart;

    // Push lost particles at the end of bins
    for( unsigned int ibin = 0 ; ibin < last_index.size() ; ibin++ ) {
        ii = indexes_of_particles_to_exchange.size()-1;
        if( ii >= 0 ) { // Push lost particles to the end of the bin
            iPart = indexes_of_particles_to_exchange[ii];
            while( iPart >= last_index[ibin] && ii > 0 ) {
                ii--;
                iPart = indexes_of_particles_to_exchange[ii];
            }
            while( iPart == last_index[ibin]-1 && iPart >= first_index[ibin] && ii > 0 ) {
                last_index[ibin]--;
                ii--;
                iPart = indexes_of_particles_to_exchange[ii];
            }
            while( iPart >= first_index[ibin] && ii > 0 ) {
                particles->overwrite_part( last_index[ibin]-1, iPart );
                last_index[ibin]--;
                ii--;
                iPart = indexes_of_particles_to_exchange[ii];
            }
            //On traite la dernière particule (qui peut aussi etre la premiere)
            if( iPart >= first_index[ibin] && iPart < last_index[ibin] ) {
                particles->overwrite_part( last_index[ibin]-1, iPart );
                last_index[ibin]--;
            }
        }
    }


    //Shift the bins in memory
    //Warning: this loop must be executed sequentially. Do not use openMP here.
    for( int unsigned ibin = 1 ; ibin < last_index.size() ; ibin++ ) { //First bin don't need to be shifted
        ii = first_index[ibin]-last_index[ibin-1]; // Shift the bin in memory by ii slots.
        iPart = min( ii, last_index[ibin]-first_index[ibin] ); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
        if( iPart > 0 ) {
            particles->overwrite_part( last_index[ibin]-iPart, last_index[ibin-1], iPart );
        }
        last_index[ibin] -= ii;
        first_index[ibin] = last_index[ibin-1];
    }



    int nmove, lmove; // local, OK
    int shift[last_index.size()+1];//how much we need to shift each bin in order to leave room for the new particle
    double dbin;

    dbin = params.cell_length[0]*params.clrw; //width of a bin.
    for( unsigned int j=0; j<last_index.size()+1 ; j++ ) {
        shift[j]=0;
    }


    int nbNeighbors_ = 2;
    int n_part_recv;

    indexes_of_particles_to_exchange.clear();
    particles->erase_particle_trail( last_index.back() );

    //Evaluation of the necessary shift of all bins.2
    for( unsigned int j=0; j<last_index.size()+1 ; j++ ) {
        shift[j]=0;
    }

    //idim=0
    shift[1] += MPIbuff.part_index_recv_sz[0][0];//Particles coming from xmin all go to bin 0 and shift all the other bins.
    shift[last_index.size()] += MPIbuff.part_index_recv_sz[0][1];//Used only to count the total number of particles arrived.
    //idim>0
    for( idim = 1; idim < ndim; idim++ ) {
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            n_part_recv = MPIbuff.part_index_recv_sz[idim][iNeighbor];
            for( unsigned int j=0; j<( unsigned int )n_part_recv ; j++ ) {
                //We first evaluate how many particles arrive in each bin.
                ii = int( ( MPIbuff.partRecv[idim][iNeighbor].position( 0, j )-min_loc )/dbin ); //bin in which the particle goes.
                shift[ii+1]++; // It makes the next bins shift.
            }
        }
    }


    //Must be done sequentially
    for( unsigned int j=1; j<last_index.size()+1; j++ ) { //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
        shift[j]+=shift[j-1];
    }
    //Make room for new particles
    if( shift[last_index.size()] ) {
        //! vecor::resize of Charge crashed ! Temporay solution : push_back / Particle
        //particles->initialize( particles->size()+shift[last_index.size()], particles->Position.size() );
        for( int inewpart=0 ; inewpart<shift[last_index.size()] ; inewpart++ ) {
            particles->create_particle();
        }
    }

    //Shift bins, must be done sequentially
    for( unsigned int j=last_index.size()-1; j>=1; j-- ) {
        int n_particles = last_index[j]-first_index[j]; //Nbr of particle in this bin
        nmove = min( n_particles, shift[j] ); //Nbr of particles to move
        lmove = max( n_particles, shift[j] ); //How far particles must be shifted
        if( nmove>0 ) {
            particles->overwrite_part( first_index[j], first_index[j]+lmove, nmove );
        }
        first_index[j] += shift[j];
        last_index[j] += shift[j];
    }

    //Space has been made now to write the arriving particles into the correct bins
    //idim == 0  is the easy case, when particles arrive either in first or last bin.
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        n_part_recv = MPIbuff.part_index_recv_sz[0][iNeighbor];
        //if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
        if( ( n_part_recv!=0 ) ) {
            ii = iNeighbor*( last_index.size()-1 ); //0 if iNeighbor=0(particles coming from Xmin) and last_index.size()-1 otherwise.
            MPIbuff.partRecv[0][iNeighbor].overwrite_part( 0, *particles, last_index[ii], n_part_recv );
            last_index[ii] += n_part_recv ;
        }
    }
    //idim > 0; this is the difficult case, when particles can arrive in any bin.
    for( idim = 1; idim < ndim; idim++ ) {
        //if (idim!=iDim) continue;
        for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
            n_part_recv = MPIbuff.part_index_recv_sz[idim][iNeighbor];
            //if ( (neighbor_[idim][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
            if( ( n_part_recv!=0 ) ) {
                for( unsigned int j=0; j<( unsigned int )n_part_recv; j++ ) {
                    ii = int( ( MPIbuff.partRecv[idim][iNeighbor].position( 0, j )-min_loc )/dbin ); //bin in which the particle goes.
                    MPIbuff.partRecv[idim][iNeighbor].overwrite_part( j, *particles, last_index[ii] );
                    last_index[ii] ++ ;
                }
            }
        }
    }


    //The width of one bin is cell_length[0] * clrw.

    int p1, p2, first_index_init;
    unsigned int bin;
    double limit;


    //Backward pass
    for( bin=0; bin<first_index.size()-1; bin++ ) { //Loop on the bins.
        limit = min_loc + ( bin+1 )*cell_length[0]*clrw;
        p1 = last_index[bin]-1;
        //If first particles change bin, they do not need to be swapped.
        while( p1 == last_index[bin]-1 && p1 >= first_index[bin] ) {
            if( particles->position( 0, p1 ) >= limit ) {
                last_index[bin]--;
            }
            p1--;
        }
        //         Now particles have to be swapped
        for( p2 = p1 ; p2 >= first_index[bin] ; p2-- ) { //Loop on the bin's particles.
            if( particles->position( 0, p2 ) >= limit ) {
                //This particle goes up one bin.
                particles->swap_part( p2, last_index[bin]-1 );
                last_index[bin]--;
            }
        }
    }
    //Forward pass + Rebracketting
    for( bin=1; bin<first_index.size(); bin++ ) { //Loop on the bins.
        limit = min_loc + bin*cell_length[0]*clrw;
        first_index_init = first_index[bin];
        p1 = first_index[bin];
        while( p1 == first_index[bin] && p1 < last_index[bin] ) {
            if( particles->position( 0, p1 ) < limit ) {
                first_index[bin]++;
            }
            p1++;
        }
        for( p2 = p1 ; p2 < last_index[bin] ; p2++ ) { //Loop on the bin's particles.
            if( particles->position( 0, p2 ) < limit ) {
                //This particle goes down one bin.
                particles->swap_part( p2, first_index[bin] );
                first_index[bin]++;
            }
        }

        //Rebracketting
        //Number of particles from bin going down is: first_index[bin]-first_index_init.
        //Number of particles from bin-1 going up is: first_index_init-last_index[bin-1].
        //Total number of particles we need to swap is the min of both.
        p2 = min( first_index[bin]-first_index_init, first_index_init-last_index[bin-1] );
        if( p2 >0 ) {
            particles->swap_part( last_index[bin-1], first_index[bin]-p2, p2 );
        }
        last_index[bin-1] += first_index[bin] - first_index_init;
        first_index[bin] = last_index[bin-1];
    }
}

void Species::initial_configuration( Params &param, Patch *patch )
{
}

void Species::configuration( Params &param, Patch *patch )
{
}

void Species::reconfiguration( Params &param, Patch *patch )
{
}

// ---------------------------------------------------------------------------------------------------------------------
// Sort particles
// ---------------------------------------------------------------------------------------------------------------------
void Species::count_sort_part( Params &params )
{
    unsigned int ip, npart, ixy, tot, oc, nxy, token;
    int ix, iy;
    double x, y;

    nxy = params.n_space[0]*params.n_space[1];
    token = ( particles == &particles_sorted[0] );

    int indices[nxy];
    npart = particles->size();
    //particles_sorted = particles ;
    particles_sorted[token].initialize( npart, *particles );

    for( unsigned int i=0; i < nxy ; i++ ) {
        indices[i] = 0 ;
    }

    // first loop counts the # of particles in each cell
    for( ip=0; ip < npart; ip++ ) {
        x = particles->position( 0, ip )-min_loc;
        y = particles->position( 1, ip )-min_loc_vec[1];

        ix = floor( x * dx_inv_[0] ) ;
        iy = floor( y * dx_inv_[1] ) ;

        ixy = iy + ix*params.n_space[1];


        indices[ixy] ++;
    }

    // second loop convert the count array in cumulative sum
    tot=0;
    for( ixy=0; ixy < nxy; ixy++ ) {
        oc = indices[ixy];
        indices[ixy] = tot;
        tot += oc;
    }
    
    // last loop puts the particles and update the count array
    for( ip=0; ip < npart; ip++ ) {
        x = particles->position( 0, ip )-min_loc;
        y = particles->position( 1, ip )-min_loc_vec[1];

        ix = floor( x * dx_inv_[1] ) ;
        iy = floor( y * dx_inv_[2] ) ;

        ixy = iy + ix*params.n_space[1];
        particles->overwrite_part( ip, particles_sorted[token], indices[ixy] );
        indices[ixy]++;
    }

    particles = &particles_sorted[token] ;

}


int Species::ParticleCreator( vector<unsigned int> n_space_to_create, Params &params, Patch *patch, int new_cell_idx )
{
    // n_space_to_create_generalized = n_space_to_create, + copy of 2nd direction data among 3rd direction
    // same for local Species::cell_length[2]
    vector<unsigned int> n_space_to_create_generalized( n_space_to_create );
    unsigned int nPart, i, j, k;
    unsigned int npart_effective = 0 ;
    double *momentum[nDim_particle], *position[nDim_particle], *weight_arr;
    std::vector<int> my_particles_indices;
    vector<Field *> xyz( nDim_field );

    // Create particles in a space starting at cell_position
    vector<double> cell_position( 3, 0 );
    vector<double> cell_index( 3, 0 );
    for( unsigned int idim=0 ; idim<nDim_field ; idim++ ) {
        //if (params.cell_length[idim]!=0) { // Useless, nDim_field defined for (params.cell_length[idim>=nDim_field]==0)
        cell_position[idim] = patch->getDomainLocalMin( idim );
        cell_index   [idim] = ( double ) patch->getCellStartingGlobalIndex( idim );
        xyz[idim] = new Field3D( n_space_to_create_generalized );
        //}
    }
    // Create the x,y,z maps where profiles will be evaluated
    vector<double> ijk( 3 );
    for( ijk[0]=0; ijk[0]<n_space_to_create_generalized[0]; ijk[0]++ ) {
        for( ijk[1]=0; ijk[1]<n_space_to_create_generalized[1]; ijk[1]++ ) {
            for( ijk[2]=0; ijk[2]<n_space_to_create_generalized[2]; ijk[2]++ ) {
                for( unsigned int idim=0 ; idim<nDim_field ; idim++ ) {
                    ( *xyz[idim] )( ijk[0], ijk[1], ijk[2] ) = cell_position[idim] + ( ijk[idim]+0.5 )*cell_length[idim];
                }
                ( *xyz[0] )( ijk[0], ijk[1], ijk[2] ) += new_cell_idx*cell_length[0];
            }
        }
    }
    // ---------------------------------------------------------
    // Calculate density and number of particles for the species
    // ---------------------------------------------------------


    // field containing the charge distribution (always 3d)
    Field3D charge( n_space_to_create_generalized );
    max_charge = 0.;

    // field containing the number of particles in each cell
    Field3D n_part_in_cell( n_space_to_create_generalized );

    // field containing the density distribution (always 3d)
    Field3D density( n_space_to_create_generalized );

    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D temperature[3];
    // field containing the temperature distribution along all 3 momentum coordinates (always 3d * 3)
    Field3D velocity[3];

    if( momentum_initialization_array_ != NULL ) {
        for( unsigned int idim = 0; idim < 3; idim++ ) {
            momentum[idim] = &( momentum_initialization_array_[idim*n_numpy_particles_] );
        }
    } else {
        //Initialize velocity and temperature profiles
        for( i=0; i<3; i++ ) {
            velocity[i].allocateDims( n_space_to_create_generalized );
            temperature[i].allocateDims( n_space_to_create_generalized );
        }
        // Evaluate profiles
        for( unsigned int m=0; m<3; m++ ) {
            if( temperature_profile_[m] ) {
                temperature_profile_[m]->valuesAt( xyz, temperature[m] );
            } else {
                temperature[m].put_to( 0.0000000001 ); // default value
            }

            if( velocity_profile_[m] ) {
                velocity_profile_[m]   ->valuesAt( xyz, velocity   [m] );
            } else {
                velocity[m].put_to( 0.0 ); //default value
            }
            
        }
    }
    // Initialize charge profile
    if( this->mass > 0 ) {
        chargeProfile ->valuesAt( xyz, charge );
    }
    if( position_initialization_array_ != NULL ) {
        for( unsigned int idim = 0; idim < nDim_particle; idim++ ) {
            position[idim] = &( position_initialization_array_[idim*n_numpy_particles_] );
        }
        weight_arr =         &( position_initialization_array_[nDim_particle*n_numpy_particles_] );
        //Idea to speed up selection, provides xmin, xmax of the bunch and check if there is an intersection with the patch instead of going through all particles for all patches.
        for( unsigned int ip = 0; ip < n_numpy_particles_; ip++ ) {
            //If the particle belongs to this patch
            if (params.geometry!="AMcylindrical") {
                if( position[0][ip] >= patch->getDomainLocalMin( 0 ) && position[0][ip] < patch->getDomainLocalMax( 0 )
                    && ( nDim_particle < 2  || ( position[1][ip] >= patch->getDomainLocalMin( 1 ) && position[1][ip] < patch->getDomainLocalMax( 1 ) ) )
                    && ( nDim_particle < 3  || ( position[2][ip] >= patch->getDomainLocalMin( 2 ) && position[2][ip] < patch->getDomainLocalMax( 2 ) ) ) ) {
                    my_particles_indices.push_back( ip ); //This vector stores particles initially sittinig in the current patch.
                }
            }
            else {
                double distance =  sqrt( position[1][ip]*position[1][ip]+position[2][ip]*position[2][ip] );
                if( position[0][ip] >= patch->getDomainLocalMin( 0 ) && position[0][ip] < patch->getDomainLocalMax( 0 )
                    && ( distance >= patch->getDomainLocalMin( 1 ) && distance < patch->getDomainLocalMax( 1 ) ) ) {
                    my_particles_indices.push_back( ip ); //This vector stores particles initially sittinig in the current patch.
                }
            }
        }
        npart_effective = my_particles_indices.size();
    } else {
        //Initialize density and ppc profiles
        density_profile_->valuesAt( xyz, density );
        particles_per_cell_profile_    ->valuesAt( xyz, n_part_in_cell );
        weight_arr = NULL;
        //Now compute number of particles per cell
        double remainder, nppc;
        for( i=0; i<n_space_to_create_generalized[0]; i++ ) {
            for( j=0; j<n_space_to_create_generalized[1]; j++ ) {
                for( k=0; k<n_space_to_create_generalized[2]; k++ ) {

                    // Obtain the number of particles per cell
                    nppc = n_part_in_cell( i, j, k );
                    n_part_in_cell( i, j, k ) = floor( nppc );
                    // If not a round number, then we need to decide how to round
                    double intpart;
                    if( modf( nppc, &intpart ) > 0 ) {
                        remainder = pow( nppc - floor( nppc ), -inv_nDim_particles );
                        if( fmod( cell_index[0]+( double )i, remainder ) < 1.
                                && fmod( cell_index[1]+( double )j, remainder ) < 1.
                                && fmod( cell_index[2]+( double )k, remainder ) < 1. ) {
                            n_part_in_cell( i, j, k )++;
                        }
                    }

                    // assign charge its correct value in the cell
                    if( this->mass > 0 ) {
                        if( charge( i, j, k )>max_charge ) {
                            max_charge=charge( i, j, k );
                        }
                    }

                    // If zero or less, zero particles
                    if( n_part_in_cell( i, j, k )<=0. || density( i, j, k )==0. ) {
                        n_part_in_cell( i, j, k ) = 0.;
                        density( i, j, k ) = 0.;
                        continue;
                    }

                    // assign density its correct value in the cell
                    if( density_profile_type_=="charge" ) {
                        if( charge( i, j, k )==0. ) {
                            ERROR( "Encountered non-zero charge density and zero charge at the same location" );
                        }
                        density( i, j, k ) /= charge( i, j, k );
                    }
                    
                    density( i, j, k ) = abs( density( i, j, k ) );
                    
                    // multiply by the cell volume
                    density( i, j, k ) *= params.cell_volume;
                    if( params.geometry=="AMcylindrical" && position_initialization_ != "regular" ) {
                        //Particles weight in regular is normalized later.
                        density( i, j, k ) *= ( *xyz[1] )( i, j, k );
                    }
                    // increment the effective number of particle by n_part_in_cell(i,j,k)
                    // for each cell with as non-zero density
                    npart_effective += ( unsigned int ) n_part_in_cell( i, j, k );

                }//i
            }//j
        }//k end the loop on all cells
    }

    // defines npart_effective for the Species & create the corresponding particles
    // -----------------------------------------------------------------------

    // if moving_win
    //     particles->create_particles(npart_effective);
    // else {
    //    // reserve included in initialize if particles emty
    //    particles->reserve(round( params->species_param[species_number_].c_part_max_ * npart_effective ), ndim);
    //    particles->initialize(n_existing_particles+npart_effective, params_->nDim_particle);
    // }

    unsigned int n_existing_particles = particles->size();
    //if (!n_existing_particles)
    particles->initialize( n_existing_particles+npart_effective, nDim_particle );

    // Initialization of the particles properties
    // ------------------------------------------
    unsigned int iPart=n_existing_particles;
    double *indexes=new double[nDim_particle];
    double *temp=new double[3];
    double *vel=new double[3];

    // start a loop on all cells

    //first_index[bin] point to begining of bin (first particle)
    //last_index[bin] point to end of bin (= first_index[bin+1])
    //if last_index = first_index, bin is empty of particle.

    if( position_initialization_array_ == NULL ) {
        for( i=0; i<n_space_to_create_generalized[0]; i++ ) {
            if((!n_existing_particles)&&( i%clrw == 0 )) {
                first_index[(new_cell_idx+i)/clrw] = iPart;
            }
            for( j=0; j<n_space_to_create_generalized[1]; j++ ) {
                for( k=0; k<n_space_to_create_generalized[2]; k++ ) {
                    // initialize particles in meshes where the density is non-zero
                    if( density( i, j, k )>0 ) {

                        vel[0]  = velocity[0]( i, j, k );
                        vel[1]  = velocity[1]( i, j, k );
                        vel[2]  = velocity[2]( i, j, k );
                        temp[0] = temperature[0]( i, j, k );
                        temp[1] = temperature[1]( i, j, k );
                        temp[2] = temperature[2]( i, j, k );
                        nPart = n_part_in_cell( i, j, k );
                        
                        //if (n_existing_particles) {
                        //    iPart = n_existing_particles;
                        //    iPart = first_index[(new_cell_idx+i)/clrw];
                        //    last_index[(new_cell_idx+i)/clrw] += nPart;
                        //    for ( int idx=(new_cell_idx+i)/clrw+1 ; idx<last_index.size() ; idx++ ) {
                        //        first_index[idx] += nPart;
                        //        last_index[idx] += nPart;
                        //    }
                        //    particles->create_particles( nPart, iPart );
                        //
                        //}

                        indexes[0]=i*cell_length[0]+cell_position[0] + new_cell_idx*cell_length[0];;
                        if( nDim_particle > 1 ) {
                            indexes[1]=j*cell_length[1]+cell_position[1];
                            if( nDim_particle > 2 ) {
                                indexes[2]=k*cell_length[2]+cell_position[2];
                            }
                        }
                        if( !position_initialization_on_species_ ) {
                            //initPosition( nPart, iPart, indexes, params );
                            ParticleCreator::createPosition( position_initialization_,
                                            particles,
                                            this,
                                            nPart,
                                            iPart,
                                            indexes,
                                            params );
                        }
                        // initMomentum( nPart, iPart, temp, vel );
                        ParticleCreator::createMomentum( momentum_initialization_,
                                                particles,
                                                this,
                                                nPart,
                                                iPart,
                                                temp,
                                                vel);
                        
                        //initWeight( nPart, iPart, density( i, j, k ) );
                        ParticleCreator::createWeight( position_initialization_, particles, nPart, iPart, density( i, j, k ), params );
                        
                        // initCharge( nPart, iPart, charge( i, j, k ) );
                        ParticleCreator::createCharge( particles, this,
                                                        nPart, iPart, charge( i, j, k ) );

                        //if (n_existing_particles) {
                        //    // operate filter
                        //    for ( int ip = nPart-1 ; ip >= 0 ; ip-- ){
                        //        double lf=1;
                        //        for (int iDim=0;iDim<3;iDim++)
                        //            lf += particles->Momentum[iDim][iPart+ip]*particles->Momentum[iDim][iPart+ip];
                        //        lf =sqrt( lf );
                        //        double X = particles->Position[0][iPart+ip] - cell_length[0]+params.timestep*particles->Momentum[0][iPart+ip]/lf;
                        //        if ( patch->isXmin() ) {
                        //            if ( ( X < 0. ) ) {
                        //                nPart--; // ne sert à rien ici
                        //                particles->erase_particle(iPart+ip);
                        //                last_index[(new_cell_idx+i)/clrw]--;
                        //                for ( int idx=(new_cell_idx+i)/clrw+1 ; idx<last_index.size() ; idx++ ) {
                        //                    first_index[idx]--;
                        //                    last_index[idx]--;
                        //                }
                        //            }
                        //            else
                        //                particles->Position[0][iPart+ip] = X;
                        //        }
                        //        else if ( patch->isXmax() ) {
                        //            X = particles->Position[0][iPart+ip] + cell_length[0]+params.timestep*particles->Momentum[0][iPart+ip]/lf;
                        //            if (  X > params.grid_length[0] ) {
                        //                //cout << "params.grid_length[0]+ cell_length[0]*vel[0] = " << params.grid_length[0]+ cell_length[0]*vel[0] << endl;
                        //                //cout << "params.grid_length[0]                        = " << params.grid_length[0] << endl;
                        //                nPart--; // ne sert à rien ici
                        //                particles->erase_particle(iPart+ip);
                        //                last_index[(new_cell_idx+i)/clrw]--;
                        //                for ( int idx=(new_cell_idx+i)/clrw+1 ; idx<last_index.size() ; idx++ ) {
                        //                    first_index[idx]--;
                        //                    last_index[idx]--;
                        //                }
                        //            }
                        //            else
                        //                particles->Position[0][iPart+ip] = X;
                        //        }
                        //
                        //    }
                        //}

                        
                        iPart+=nPart;
                    }//END if density > 0
                }//k end the loop on all cells
            }//j
            if((!n_existing_particles)&&( i%clrw == clrw -1 )) {
                last_index[(new_cell_idx+i)/clrw] = iPart;
            }

        }//i
    } else if( n_existing_particles == 0 ) {
        // Here position are created from a numpy array.
        // Do not recreate particles from numpy array again after initialization. Is this condition enough ?
        // Initializing particles from numpy array and based on a count sort to comply with initial sorting.
        int nbins = first_index.size();
        int indices[nbins];
        double one_ov_dbin = 1. / ( cell_length[0] * clrw ) ;

        for( int ibin=0; ibin < nbins ; ibin++ ) {
            indices[ibin] = 0 ;
        }

        ///Compute proper indices for particle susing a count sort
        for( unsigned int ipart = 0; ipart < npart_effective ; ipart++ ) {
            unsigned int ip = my_particles_indices[ipart];
            double x = position[0][ip]-min_loc ;
            int ix = int( x * one_ov_dbin ) ;
            indices[ix] ++;
        }
        unsigned int tot=0;
        for( int ibin=0; ibin < nbins; ibin++ ) {
            unsigned int oc = indices[ibin];
            indices[ibin] = tot;
            tot += oc;
        }
        for( int ibin=0; ibin < nbins   ; ibin++ ) {
            first_index[ibin] = indices[ibin] ;
        }
        for( int ibin=0; ibin < nbins-1 ; ibin++ ) {
            last_index[ibin] = first_index[ibin+1] ;
        }
        last_index[nbins-1] = npart_effective ;

        //Now initialize particles at their proper indices
        for( unsigned int ipart = 0; ipart < npart_effective ; ipart++ ) {
            unsigned int ippy = my_particles_indices[ipart];//Indice of the particle in the python array.
            double x = position[0][ippy]-min_loc ;
            unsigned int ibin = int( x * one_ov_dbin ) ;
            int ip = indices[ibin] ; //Indice of the position of the particle in the particles array.

            unsigned int int_ijk[3] = {0, 0, 0};
            if ( params.geometry != "AMcylindrical") {
                for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
                    particles->position( idim, ip ) = position[idim][ippy];
                    int_ijk[idim] = ( unsigned int )( ( particles->position( idim, ip ) - min_loc_vec[idim] )/cell_length[idim] );
                }
            }
            else {
                for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
                    particles->position( idim, ip ) = position[idim][ippy];
                }
                int_ijk[0] = ( unsigned int )( ( particles->position( 0, ip ) - min_loc_vec[0] )/cell_length[0] );
                int_ijk[1] = ( unsigned int )( ( sqrt( position[1][ippy]*position[1][ippy]+position[2][ippy]*position[2][ippy] ) - min_loc_vec[1] )/cell_length[1] );
            }
            if( !momentum_initialization_array_ ) {
                vel [0] = velocity   [0]( int_ijk[0], int_ijk[1], int_ijk[2] );
                vel [1] = velocity   [1]( int_ijk[0], int_ijk[1], int_ijk[2] );
                vel [2] = velocity   [2]( int_ijk[0], int_ijk[1], int_ijk[2] );
                temp[0] = temperature[0]( int_ijk[0], int_ijk[1], int_ijk[2] );
                temp[1] = temperature[1]( int_ijk[0], int_ijk[1], int_ijk[2] );
                temp[2] = temperature[2]( int_ijk[0], int_ijk[1], int_ijk[2] );
                //initMomentum( 1, ip, temp, vel );
                ParticleCreator::createMomentum( momentum_initialization_,
                                        particles,
                                        this,
                                        1,
                                        ip,
                                        temp,
                                        vel);
            } else {
                for( unsigned int idim=0; idim < 3; idim++ ) {
                    particles->momentum( idim, ip ) = momentum[idim][ippy]/this->mass ;
                }
            }

            particles->weight( ip ) = weight_arr[ippy] ;
            // initCharge( 1, ip, charge( int_ijk[0], int_ijk[1], int_ijk[2] ) );
            ParticleCreator::createCharge( particles, this,
                            1, ip, charge( int_ijk[0], int_ijk[1], int_ijk[2] ) );
            indices[ibin]++;
        }
    }

    // Delete map xyz.
    for( unsigned int idim=0 ; idim<nDim_field ; idim++ ) {
        delete xyz[idim];
    }

    delete [] indexes;
    delete [] temp;
    delete [] vel;
    // Recalculate former position using the particle velocity
    // (necessary to calculate currents at time t=0 using the Esirkepov projection scheme)
    if( patch->isXmax() ) {
        // Matter particle case
        if( mass > 0 ) {
            for( iPart=n_existing_particles; iPart<n_existing_particles+npart_effective; iPart++ ) {
                /*897 for (int i=0; i<(int)nDim_particle; i++) {
                  particles->position_old(i,iPart) -= particles->momentum(i,iPart)/particles->lor_fac(iPart) * params.timestep;
                  }897*/
                new_particles_energy_ += particles->weight( iPart )*( particles->lor_fac( iPart )-1.0 );
            }
        }
        // Photon case
        else if( mass == 0 ) {
            for( iPart=n_existing_particles; iPart<n_existing_particles+npart_effective; iPart++ ) {
                /*897 for (int i=0; i<(int)nDim_particle; i++) {
                  particles->position_old(i,iPart) -= particles->momentum(i,iPart)/particles->lor_fac(iPart) * params.timestep;
                  }897*/
                new_particles_energy_ += particles->weight( iPart )*( particles->momentum_norm( iPart ) );
            }
        }
    }
    if( particles->tracked ) {
        particles->resetIds();
    }
    return npart_effective;

} // End ParticleCreator


// Move all particles from another species to this one
void Species::importParticles( Params &params, Patch *patch, Particles &source_particles, vector<Diagnostic *> &localDiags )
{
    unsigned int npart = source_particles.size(), ibin, ii, nbin=first_index.size();
    double inv_cell_length = 1./ params.cell_length[0];

    // std::cerr << "Species::importParticles "
    //           << " for "<< this->name
    //           << " in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ") "
    //           << " mpi process " << patch->MPI_me_ << " - "
    //           << " mode: " << this->vectorized_operators << " - "
    //           << " nb bin: " << first_index.size() << " - "
    //           << " nbp: " << npart
    //           << std::endl;

    // If this species is tracked, set the particle IDs
    if( particles->tracked ) {
        dynamic_cast<DiagnosticTrack *>( localDiags[tracking_diagnostic] )->setIDs( source_particles );
    }

    // Move particles
    vector<int> src_bin_keys( npart, 0 );
    for( unsigned int i=0; i<npart; i++ ) {
        // Copy particle to the correct bin
        src_bin_keys[i] = source_particles.position( 0, i )*inv_cell_length - ( patch->getCellStartingGlobalIndex( 0 ) + params.oversize[0] );
        src_bin_keys[i] /= params.clrw;
    }

    vector<int> bin_count( nbin, 0 );
    for( unsigned int ip=0; ip < npart ; ip++ )
        bin_count[src_bin_keys[ip]] ++;

    // sort new parts par bins
    int istart = 0;
    int istop  = bin_count[0];

    for ( int ibin = 0 ; ibin < nbin ; ibin++ ) {
        if (bin_count[ibin]!=0) {
            for( unsigned int ip=istart; ip < istop ; ip++ ) {
                if ( src_bin_keys[ip] == ibin )
                    continue;
                else { // rearrange particles
                    int ip_swap = istop;
                    while (( src_bin_keys[ip_swap] != ibin ) && (ip_swap<npart))
                        ip_swap++;
                    source_particles.swap_part(ip, ip_swap);
                    int tmp = src_bin_keys[ip];
                    src_bin_keys[ip] = src_bin_keys[ip_swap];
                    src_bin_keys[ip_swap] = tmp;
                } // rearrange particles
            } // end loop on particles of a cell

            // inject in main data structure per cell
            source_particles.cp_particles( istart, bin_count[ibin],
                                        *particles,
                                        first_index[ibin] );
            last_index[ibin] += bin_count[ibin];
            for ( int idx=ibin+1 ; idx<last_index.size() ; idx++ ) {
                first_index[idx] += bin_count[ibin];
                last_index[idx]  += bin_count[ibin];
            }

        }
        // update istart/istop fot the next cell
        istart += bin_count[ibin];
        if ( ibin != nbin-1  )
            istop  += bin_count[ibin+1];
        else
            istop = npart;

    } // End cell loop

    source_particles.clear();
}


// ------------------------------------------------
// Set position when using restart & moving window
// patch are initialized with t0 position
// ------------------------------------------------
//void Species::updateMvWinLimits(double x_moved)
//{
//    partBoundCond->updateMvWinLimits(x_moved);
//    min_loc += x_moved;
//
//} // End updateMvWinLimits


//Do we have to project this species ?
bool Species::isProj( double time_dual, SimWindow *simWindow )
{

    return time_dual > time_frozen  || ( simWindow->isMoving( time_dual ) ) ;

    //Recompute frozen particles density if
    //moving window is activated, actually moving at this time step, and we are not in a density slope.
    /*    bool isproj =(time_dual > species_param.time_frozen  ||
                 (simWindow && simWindow->isMoving(time_dual) &&
                     (species_param.species_geometry == "gaussian" ||
                         (species_param.species_geometry == "trapezoidal" &&
                            //Before end of density ramp up.
                            (simWindow->getXmoved() < species_param.vacuum_length[0] + species_param.dens_length_x[1] + oversize[0]*cell_length[0] ||
                            //After begining of density ramp down.
                            simWindow->getXmoved() +  simWindow->getNspace_win_x()*cell_length[0] > species_param.vacuum_length[0] + species_param.dens_length_x[1]+ species_param.dens_length_x[0]
                            )
                        )
                    )
                )
            );
            return isproj;*/
    //return time_dual > species_param.time_frozen  || (simWindow && simWindow->isMoving(time_dual)) ;
}

void Species::disableXmax()
{
    partBoundCond->bc_xmax   = NULL;
}

void Species::setXminBoundaryCondition()
{
    partBoundCond->bc_xmin   = &remove_particle;
}

// ---------------------------------------------------------------------------------------------------------------------
// Particle merging cell by cell
// ---------------------------------------------------------------------------------------------------------------------
void Species::mergeParticles( double time_dual, unsigned int ispec,
                              Params &params,
                              Patch *patch, SmileiMPI *smpi,
                              std::vector<Diagnostic *> &localDiags ) {}

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the fields at the particle position
//   - deposit susceptibility
//   - calculate the new momentum
// ---------------------------------------------------------------------------------------------------------------------
void Species::ponderomotive_update_susceptibility_and_momentum( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags )
{

    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen ) { // moving particle

        smpi->dynamics_resize( ithread, nDim_field, last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin++ ) { // loop on ibin

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( first_index[ibin] ), &( last_index[ibin] ), ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[7] += MPI_Wtime() - timer;
#endif


            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Proj->susceptibility( EMfields, *particles, mass, smpi, first_index[ibin], last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[8] += MPI_Wtime() - timer;
#endif


#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push only the particle momenta
            ( *Push )( *particles, smpi, first_index[ibin], last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[9] += MPI_Wtime() - timer;
#endif

        } // end loop on ibin
    } else { // immobile particle
    } //END if time vs. time_frozen
} // ponderomotive_update_susceptibility_and_momentum

// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the fields at the particle position
//   - deposit susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void Species::ponderomotive_project_susceptibility( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags )
{

    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    // -------------------------------
    // calculate the particle updated momentum
    // -------------------------------
    if( time_dual>time_frozen ) { // moving particle

        smpi->dynamics_resize( ithread, nDim_particle, last_index.back(), false );

        for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin++ ) { // loop on ibin

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Interp->fieldsAndEnvelope( EMfields, *particles, smpi, &( first_index[ibin] ), &( last_index[ibin] ), ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[7] += MPI_Wtime() - timer;
#endif


            // Project susceptibility, the source term of envelope equation
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Proj->susceptibility( EMfields, *particles, mass, smpi, first_index[ibin], last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[8] += MPI_Wtime() - timer;
#endif


        } // end loop on ibin
    } else { // immobile particle
    } //END if time vs. time_frozen
} // ponderomotive_project_susceptibility


// ---------------------------------------------------------------------------------------------------------------------
// For all particles of the species reacting to laser envelope
//   - interpolate the ponderomotive potential and its gradient at the particle position, for present and previous timestep
//   - calculate the new particle position
//   - particles BC
//   - project charge and current density
// ---------------------------------------------------------------------------------------------------------------------
void Species::ponderomotive_update_position_and_currents( double time_dual, unsigned int ispec,
        ElectroMagn *EMfields,
        Params &params, bool diag_flag, PartWalls *partWalls,
        Patch *patch, SmileiMPI *smpi,
        vector<Diagnostic *> &localDiags )
{

    int ithread;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

#ifdef  __DETAILED_TIMERS
    double timer;
#endif

    unsigned int iPart;

    // Reset list of particles to exchange - WARNING Should it be reset?
    clearExchList();

    int tid( 0 );
    double ener_iPart( 0. );
    std::vector<double> nrj_lost_per_thd( 1, 0. );

    // -------------------------------
    // calculate the particle updated position
    // -------------------------------
    if( time_dual>time_frozen ) { // moving particle

        smpi->dynamics_resize( ithread, nDim_field, last_index.back(), params.geometry=="AMcylindrical" );

        for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin++ ) {

            // Interpolate the ponderomotive potential and its gradient at the particle position, present and previous timestep
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            Interp->timeCenteredEnvelope( EMfields, *particles, smpi, &( first_index[ibin] ), &( last_index[ibin] ), ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[10] += MPI_Wtime() - timer;
#endif

#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            // Push only the particle position
            ( *Push_ponderomotive_position )( *particles, smpi, first_index[ibin], last_index[ibin], ithread );
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[11] += MPI_Wtime() - timer;
#endif

            // Apply wall and boundary conditions
            if( mass>0 ) {
                for( unsigned int iwall=0; iwall<partWalls->size(); iwall++ ) {
                    for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                        double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                        if( !( *partWalls )[iwall]->apply( *particles, iPart, this, dtgf, ener_iPart ) ) {
                            nrj_lost_per_thd[tid] += mass * ener_iPart;
                        }
                    }
                }

                // Boundary Condition may be physical or due to domain decomposition
                // apply returns 0 if iPart is not in the local domain anymore
                //        if omp, create a list per thread
                for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                    if( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                        addPartInExchList( iPart );
                        nrj_lost_per_thd[tid] += mass * ener_iPart;
                    }
                }

            } else if( mass==0 ) {
                ERROR( "Particles with zero mass cannot interact with envelope" );
                // for(unsigned int iwall=0; iwall<partWalls->size(); iwall++) {
                //     for (iPart=first_index[ibin] ; (int)iPart<last_index[ibin]; iPart++ ) {
                //         double dtgf = params.timestep * smpi->dynamics_invgf[ithread][iPart];
                //         if ( !(*partWalls)[iwall]->apply(*particles, iPart, this, dtgf, ener_iPart)) {
                //                 nrj_lost_per_thd[tid] += ener_iPart;
                //         }
                //     }
                // }
                //
                // // Boundary Condition may be physical or due to domain decomposition
                // // apply returns 0 if iPart is not in the local domain anymore
                // //        if omp, create a list per thread
                // for (iPart=first_index[ibin] ; (int)iPart<last_index[ibin]; iPart++ ) {
                //     if ( !partBoundCond->apply( *particles, iPart, this, ener_iPart ) ) {
                //         addPartInExchList( iPart );
                //         nrj_lost_per_thd[tid] += ener_iPart;
                //     }
                //  }

            } // end mass = 0? condition

            //START EXCHANGE PARTICLES OF THE CURRENT BIN ?

            // Project currents if not a Test species and charges as well if a diag is needed.
            // Do not project if a photon
#ifdef  __DETAILED_TIMERS
            timer = MPI_Wtime();
#endif
            if( ( !particles->is_test ) && ( mass > 0 ) ) {
                Proj->currentsAndDensityWrapper( EMfields, *particles, smpi, first_index[ibin], last_index[ibin], ithread, diag_flag, params.is_spectral, ispec );
            }
#ifdef  __DETAILED_TIMERS
            patch->patch_timers[12] += MPI_Wtime() - timer;
#endif

        } // end ibin loop

        for( unsigned int ithd=0 ; ithd<nrj_lost_per_thd.size() ; ithd++ ) {
            nrj_bc_lost += nrj_lost_per_thd[tid];
        }

    } // end case of moving particle
    else { // immobile particle

        if( diag_flag &&( !particles->is_test ) ) {
            double *b_rho=nullptr;
            for( unsigned int ibin = 0 ; ibin < first_index.size() ; ibin ++ ) { //Loop for projection on buffer_proj
                // only 3D is implemented actually
                if( nDim_field==2 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                }
                if( nDim_field==3 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                } else if( nDim_field==1 ) {
                    b_rho = EMfields->rho_s[ispec] ? &( *EMfields->rho_s[ispec] )( 0 ) : &( *EMfields->rho_ )( 0 ) ;
                }
                for( iPart=first_index[ibin] ; ( int )iPart<last_index[ibin]; iPart++ ) {
                    Proj->basic( b_rho, ( *particles ), iPart, 0 );
                } //End loop on particles
            }//End loop on bins
        } // end condition on diag and not particle test

    }//END if time vs. time_frozen
} // End ponderomotive_position_update

void Species::check( Patch *patch, std::string title )
{
    double sum_x = 0;
    double sum_y = 0;
    double sum_z = 0;
    double sum_px = 0;
    double sum_py = 0;
    double sum_pz = 0;
    double sum_w = 0;
    unsigned int sum_ck = 0;
    for( unsigned int ip=0; ip < particles->size() ; ip++ ) {
        sum_x += particles->position( 0, ip );
        sum_y += particles->position( 1, ip );
        sum_z += particles->position( 2, ip );
        sum_px += particles->momentum( 0, ip );
        sum_py += particles->momentum( 1, ip );
        sum_pz += particles->momentum( 1, ip );
        sum_w += particles->weight( ip );
        sum_ck += particles->cell_keys[ip];
    }
    std::cerr << "Check sum at " << title
              << " for "<< this->name
              << " in patch (" << patch->Pcoordinates[0] << "," <<  patch->Pcoordinates[1] << "," <<  patch->Pcoordinates[2] << ") "
              << " mpi process " << patch->MPI_me_
              << " - mode: " << this->vectorized_operators
              << " - nb bin: " << first_index.size()
              << " - nbp: " << particles->size()
              << " - w: " << sum_w
              << " - x: " << sum_x
              << " - y: " << sum_y
              << " - z: " << sum_z
              << " - px: " << sum_px
              << " - py: " << sum_py
              << " - pz: " << sum_pz
              << " - ck: " << sum_ck
              << '\n';
};
