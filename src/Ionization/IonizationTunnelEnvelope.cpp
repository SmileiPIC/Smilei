#include "IonizationTunnelEnvelope.h"
#include "IonizationTables.h"

#include <cmath>

#include "Particles.h"
#include "Species.h"

using namespace std;



IonizationTunnelEnvelope::IonizationTunnelEnvelope( Params &params, Species *species ) : Ionization( params, species )
{
    DEBUG( "Creating the Tunnel Envelope Ionizaton class" );
    
    atomic_number_          = species->atomic_number;
    
    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    Potential.resize( atomic_number_ );
    Azimuthal_quantum_number.resize( atomic_number_ );
    for( int Zstar=0; Zstar<( int )atomic_number_; Zstar++ ) {
        Potential               [Zstar] = IonizationTables::ionization_energy( atomic_number_, Zstar ) * eV_to_au;
        Azimuthal_quantum_number[Zstar] = IonizationTables::azimuthal_atomic_number( atomic_number_, Zstar );
    }
    
    for( unsigned int i=0; i<atomic_number_; i++ ) {
        DEBUG( "ioniz: i " << i << " potential: " << Potential[i] << " Az.q.num: " << Azimuthal_quantum_number[i] );
    }
    
    one_third = 1.0/3.0;
    
    alpha_tunnel.resize( atomic_number_ );
    beta_tunnel.resize( atomic_number_ );
    gamma_tunnel.resize( atomic_number_ );
    //Ip_times2_to_minus3ov4.resize( atomic_number_ );
    
    for( unsigned int Z=0 ; Z<atomic_number_ ; Z++ ) {
        DEBUG( "Z : " << Z );
        double cst      = ( ( double )Z+1.0 ) * sqrt( 2.0/Potential[Z] );
        alpha_tunnel[Z] = cst-1.0; // 2(n^*)-1
        beta_tunnel[Z]  = pow( 2, alpha_tunnel[Z] ) * ( 8.*Azimuthal_quantum_number[Z]+4.0 ) / ( cst*tgamma( cst ) ) * Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow( 2.0*Potential[Z], 1.5 );   // 2*(2I_p)^{3/2}
        //Ip_times2_to_minus3ov4[Z] = pow(2.*Potential[Z],-0.75); // (2I_p)^{-3/4}
    }
    
    ellipticity         = params.envelope_ellipticity;
    cos_phi             = cos(params.envelope_polarization_phi);
    sin_phi             = sin(params.envelope_polarization_phi);
    ellipticity_factor1 = (1.-params.envelope_ellipticity)/3./params.envelope_ellipticity;
    ellipticity_factor2 = params.envelope_ellipticity*(1+params.envelope_ellipticity)/2.;
    ellipticity_factor2 = pow(ellipticity_factor2,-0.5);

    DEBUG( "Finished Creating the Tunnel Envelope Ionizaton class" );
    
}



void IonizationTunnelEnvelope::operator()( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref )
{}


void IonizationTunnelEnvelope::envelopeIonization( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, std::vector<double> *Epart, std::vector<double> *EnvEabs_part, std::vector<double> *Phipart, Patch *patch, Projector *Proj, int ipart_ref )
{
    unsigned int Z, Zp1, newZ, k_times;
    double E, E_sq, EnvE_sq, Aabs, invE, delta, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    //double rms_momentum_major_axis, rms_momentum_minor_axis, delta_momentum_spread, coeff_ellipticity;
    double momentum_major_axis, momentum_minor_axis; //, rand_gaussian;
    vector<double> IonizRate_tunnel_envelope( atomic_number_ ), Dnom_tunnel( atomic_number_ );
    double ran_p_times_2pi;
    
    
    int nparts = Epart->size()/3;
    double *Ex      = &( ( *Epart )[0*nparts] );
    double *Ey      = &( ( *Epart )[1*nparts] );
    double *Ez      = &( ( *Epart )[2*nparts] );
    double *E_env   = &( ( *EnvEabs_part )[0*nparts] );
    double *Phi_env = &( ( *Phipart )[0*nparts] );
    
    for( unsigned int ipart=ipart_min ; ipart<ipart_max; ipart++ ) {
    
        // Current charge state of the ion
        Z = ( unsigned int )( particles->charge( ipart ) );
    
        // If ion already fully ionized then skip
        if( Z==atomic_number_ ) {
            continue;
        }
    
    
        // Absolute value of the electric field normalized in atomic units
        E_sq    = pow(EC_to_au,2) * (pow( *( Ex+ipart-ipart_ref ), 2 )
                                    +pow( *( Ey+ipart-ipart_ref ), 2 )
                                    +pow( *( Ez+ipart-ipart_ref ), 2 ) );
        // Envelope field normalized in atomic units
        EnvE_sq = pow(EC_to_au,2)*( pow( *( E_env+ipart-ipart_ref ), 2 ) );
        // Absolute value of envelope, necessary for the computation of the momentum of new electrons
        Aabs    = sqrt(2. * (*(Phi_env+ipart-ipart_ref))  );  // because Phi = |A|^2/2

        // Effective electric field for ionization:
        // sqrt(E_plasma^2+E_envelope^2)
        E = sqrt(E_sq+EnvE_sq);
    
        if( E<1e-10 ) {
            continue;
        }
    
        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------
    
        invE = 1./E;
        delta      = gamma_tunnel[Z]*invE; // 2*(2I_p)^{3/2}/E
        ran_p = patch->xorshift32() * patch->xorshift32_invmax;
        IonizRate_tunnel_envelope[Z] = beta_tunnel[Z] * exp( -delta*one_third + alpha_tunnel[Z]*log( delta ) );
    
        // k_times will give the nb of ionization events
        k_times = 0;
        Zp1=Z+1;

        if( Zp1 == atomic_number_ ) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if( ran_p < 1.0 -exp( -IonizRate_tunnel_envelope[Z]*dt ) ) {
                k_times        = 1;
            }
    
        } else {
            // else : multiple ionization can occur in one time-step
            //        partial & final ionization are decoupled (see Nuter Phys. Plasmas)
            // -------------------------------------------------------------------------
    
            // initialization
            Mult = 1.0;
            Dnom_tunnel[0]=1.0;
            Pint_tunnel = exp( -IonizRate_tunnel_envelope[Z]*dt ); // cummulative prob.
    
            //multiple ionization loop while Pint_tunnel < ran_p and still partial ionization
            while( ( Pint_tunnel < ran_p ) and ( k_times < atomic_number_-Zp1 ) ) {
                newZ = Zp1+k_times;
                delta = gamma_tunnel[newZ]*invE;
                
                IonizRate_tunnel_envelope[newZ] =  beta_tunnel[newZ]
                                         *         exp( -delta*one_third+alpha_tunnel[newZ]*log( delta ) );

                D_sum = 0.0;
                P_sum = 0.0;
                Mult  *= IonizRate_tunnel_envelope[Z+k_times];
                for( unsigned int i=0; i<k_times+1; i++ ) {
                    Dnom_tunnel[i]=Dnom_tunnel[i]/( IonizRate_tunnel_envelope[newZ]-IonizRate_tunnel_envelope[Z+i] );
                    D_sum += Dnom_tunnel[i];
                    P_sum += exp( -IonizRate_tunnel_envelope[Z+i]*dt )*Dnom_tunnel[i];
                }
                Dnom_tunnel[k_times+1] -= D_sum;
                P_sum                   = P_sum + Dnom_tunnel[k_times+1]*exp( -IonizRate_tunnel_envelope[newZ]*dt );
                Pint_tunnel             = Pint_tunnel + P_sum*Mult;
    
                k_times++;
            }//END while
    
            // final ionization (of last electron)
            if( ( ( 1.0-Pint_tunnel )>ran_p ) && ( k_times==atomic_number_-Zp1 ) ) {
                k_times++;
            }
        }//END Multiple ionization routine
    
        // Ionization ion current cannot be computed with the envelope ionization model
        
        // Creation of the new electrons
      
        // each electron will acquire a transverse momentum equal to p_perp = eA 
        // where the vector A(t) is equal to A*(vec_major_axis*cos(omega*t) + vec_minor_axis*epsilon*sin(omega*t))
        // where vec_major_axis and vec_minor_axis are the major and minor axes of the polarization ellipse
        // epsilon is the polarization ellipticity
        // since the envelope does not contain information on the oscillation phase,
        // omega*t will be extracted randomly between 0 and 2*pi
          
        // draw from uniform distribution
        ran_p_times_2pi= 2. * M_PI * patch->xorshift32() * patch->xorshift32_invmax;
        
        momentum_major_axis =               Aabs * cos(ran_p_times_2pi);
        momentum_minor_axis = ellipticity * Aabs * sin(ran_p_times_2pi);
    
        if( k_times !=0 ) {
            new_electrons.create_particle();
            //new_electrons.initialize( new_electrons.size()+1, new_electrons.dimension() );
            int idNew = new_electrons.size() - 1;
            for( unsigned int i=0; i<new_electrons.dimension(); i++ ) {
                new_electrons.position( i, idNew )=particles->position( i, ipart );
            }
            for( unsigned int i=0; i<3; i++ ) {
                new_electrons.momentum( i, idNew ) = particles->momentum( i, ipart )*ionized_species_invmass;
            }
    
            // add the momentum p=eA, back-transformed to y-z coordinates
            // no changes are made on the x momentum, as the envelope propagates in that direction
            new_electrons.momentum( 1, idNew ) += momentum_major_axis*cos_phi-momentum_minor_axis*sin_phi;
            new_electrons.momentum( 2, idNew ) += momentum_major_axis*sin_phi+momentum_minor_axis*cos_phi;
    
            new_electrons.weight( idNew )=double( k_times )*particles->weight( ipart );
            new_electrons.charge( idNew )=-1;
    
            // Increase the charge of the particle
            particles->charge( ipart ) += k_times;
        }
    
    
    } // Loop on particles

}

