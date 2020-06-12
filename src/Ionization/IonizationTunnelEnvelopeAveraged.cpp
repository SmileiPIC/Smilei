#include "IonizationTunnelEnvelopeAveraged.h"
#include "IonizationTables.h"

#include <cmath>

#include "Particles.h"
#include "Species.h"

using namespace std;



IonizationTunnelEnvelopeAveraged::IonizationTunnelEnvelopeAveraged( Params &params, Species *species ) : Ionization( params, species )
{
    DEBUG( "Creating the Tunnel Envelope Ionizaton Averaged class" );
    
    atomic_number_          = species->atomic_number_;
    
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
    
    alpha_tunnel.resize( atomic_number_+1 );
    beta_tunnel.resize( atomic_number_+1 );
    gamma_tunnel.resize( atomic_number_+1 );
    Ip_times2_to_minus3ov4.resize( atomic_number_+1 );
    
    for( unsigned int Z=0 ; Z<atomic_number_ ; Z++ ) {
        DEBUG( "Z : " << Z );
        double cst      = ( ( double )Z+1.0 ) * sqrt( 2.0/Potential[Z] );
        alpha_tunnel[Z] = cst-1.0; // 2(n^*)-1
        beta_tunnel[Z]  = pow( 2, alpha_tunnel[Z] ) * ( 8.*Azimuthal_quantum_number[Z]+4.0 ) / ( cst*tgamma( cst ) ) * Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow( 2.0*Potential[Z], 1.5 );   // 2*(2I_p)^{3/2}
        Ip_times2_to_minus3ov4[Z] = pow(2.*Potential[Z],-0.75); // (2I_p)^{-3/4}
    }
    
    ellipticity         = params.envelope_ellipticity;
    
    cos_phi             = cos(params.envelope_polarization_phi);
    sin_phi             = sin(params.envelope_polarization_phi);
    
    DEBUG( "Finished Creating the Tunnel Envelope Ionizaton Averaged class" );
    
}



void IonizationTunnelEnvelopeAveraged::operator()( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref )
{}


void IonizationTunnelEnvelopeAveraged::envelopeIonization( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, std::vector<double> *Epart, std::vector<double> *EnvEabs_part, std::vector<double> *EnvExabs_part, std::vector<double> *Phipart, Patch *patch, Projector *Proj, int ipart_ref )
{
    unsigned int Z, Zp1, newZ, k_times;
    double E, E_sq, EnvE_sq, Aabs, invE, delta, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    double coeff_ellipticity_in_ionization_rate;
    double p_perp; 
    vector<double> IonizRate_tunnel_envelope( atomic_number_ ), Dnom_tunnel( atomic_number_ );
    
    
    int nparts = Epart->size()/3;
    double *Ex      = &( ( *Epart )[0*nparts] );
    double *Ey      = &( ( *Epart )[1*nparts] );
    double *Ez      = &( ( *Epart )[2*nparts] );
    double *E_env   = &( ( *EnvEabs_part )[0*nparts] );
    double *Ex_env  = &( ( *EnvExabs_part )[0*nparts] );
    double *Phi_env = &( ( *Phipart )[0*nparts] );
    
    for( unsigned int ipart=ipart_min ; ipart<ipart_max; ipart++ ) {
    
        // Current charge state of the ion
        Z = ( unsigned int )( particles->charge( ipart ) );
    
        // If ion already fully ionized then skip
        if( Z==atomic_number_ ) {
            continue;
        }
    
    
        // Absolute value of the electric field |E_plasma| (from the plasma) normalized in atomic units
        E_sq    = pow(EC_to_au,2) * (pow( *( Ex+ipart-ipart_ref ), 2 )
                                    +pow( *( Ey+ipart-ipart_ref ), 2 )
                                    +pow( *( Ez+ipart-ipart_ref ), 2 ) );
        // Laser envelope electric field normalized in atomic units, using both transverse and longitudinal components:
        // |E_envelope|^2 = |Env_E|^2 + |Env_Ex|^2

        EnvE_sq = pow(EC_to_au,2)*( pow( *( E_env+ipart-ipart_ref ), 2 ) ) + pow(EC_to_au,2)*( pow( *( Ex_env+ipart-ipart_ref ), 2 ) );

        // Effective electric field for ionization:
        // |E| = sqrt(|E_plasma|^2+|E_envelope|^2)
        E = sqrt(E_sq+EnvE_sq);
    
        if( E<1e-10 ) {
            continue;
        }
    
        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------
    
        invE = 1./E;
        delta      = gamma_tunnel[Z]*invE; // 2*(2I_p)^{3/2}/E
        ran_p = patch->rand_->uniform();
        IonizRate_tunnel_envelope[Z] = beta_tunnel[Z] * exp( -delta*one_third + alpha_tunnel[Z]*log( delta ) );
    
        // Corrections on averaged ionization rate given by the polarization ellipticity  
        if (ellipticity==0.){ // linear polarization
            coeff_ellipticity_in_ionization_rate = pow((3./M_PI)/delta*2.,0.5);
        } else if (ellipticity==1.){ // circular polarization
            coeff_ellipticity_in_ionization_rate = 1.; // for circular polarization, the ionization rate is unchanged
        }

        IonizRate_tunnel_envelope[Z] = coeff_ellipticity_in_ionization_rate * IonizRate_tunnel_envelope[Z];
        double Ip_times2_power_minus3ov4=0.;
    
        // k_times will give the nb of ionization events
        k_times = 0;
        Zp1=Z+1;

        if( Zp1 == atomic_number_ ) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if( ran_p < 1.0 -exp( -IonizRate_tunnel_envelope[Z]*dt ) ) {
                k_times        = 1;
                Ip_times2_power_minus3ov4 = Ip_times2_to_minus3ov4[Z];
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

                // Corrections on averaged ionization rate given by the polarization ellipticity  
                if (ellipticity==0.){ // linear polarization
                    coeff_ellipticity_in_ionization_rate = pow((3./M_PI)/(gamma_tunnel[newZ-1]*invE)*2.,0.5);
                } else if (ellipticity==1.){ // circular polarization
                    coeff_ellipticity_in_ionization_rate = 1.; // for circular polarization, the ionization rate is unchanged
                }

                IonizRate_tunnel_envelope[newZ] = coeff_ellipticity_in_ionization_rate * beta_tunnel[newZ]
                                                  * exp( -delta*one_third+alpha_tunnel[newZ]*log( delta ) );

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
           
                Ip_times2_power_minus3ov4 += Ip_times2_to_minus3ov4[newZ-1];
            }//END while
    
            // final ionization (of last electron)
            if( ( ( 1.0-Pint_tunnel )>ran_p ) && ( k_times==atomic_number_-Zp1 ) ) {
                k_times++;
                Ip_times2_power_minus3ov4 += Ip_times2_to_minus3ov4[atomic_number_-1];
            }
        }//END Multiple ionization routine
    
        // ---- Ionization ion current cannot be computed with the envelope ionization model
      
        // ---- Creation of the new electrons
        
        if( k_times !=0 ) {
            new_electrons.createParticle();
            //new_electrons.initialize( new_electrons.size()+1, new_electrons.dimension() );
            int idNew = new_electrons.size() - 1;

            // The new electron is in the same position of the atom where it originated from
            for( unsigned int i=0; i<new_electrons.dimension(); i++ ) {
                new_electrons.position( i, idNew )=particles->position( i, ipart );
            }
            for( unsigned int i=0; i<3; i++ ) {
                new_electrons.momentum( i, idNew ) = particles->momentum( i, ipart )*ionized_species_invmass;
            }

           
            // ----  Initialise the momentum, weight and charge of the new electron

            if (ellipticity==0.){ // linear polarization

                double rand_gaussian  = patch->rand_->normal();

                Aabs    = sqrt(2. * (*(Phi_env+ipart-ipart_ref))  ); // envelope of the laser vector potential component along the polarization direction
                
                // recreate gaussian distribution with rms momentum spread for linear polarization, estimated by C.B. Schroeder
                // C. B. Schroeder et al., Phys. Rev. ST Accel. Beams 17, 2014, first part of Eqs. 7,10 
                p_perp = rand_gaussian * Aabs * sqrt(1.5*E) * Ip_times2_power_minus3ov4;         

                // add the transverse momentum p_perp to obtain a gaussian distribution 
                // in the momentum in the polarization direction p_perp, following Schroeder's result
                new_electrons.momentum( 1, idNew ) += p_perp*cos_phi;
                new_electrons.momentum( 2, idNew ) += p_perp*sin_phi;

                // initialize px to take into account the average drift <px>=A^2/4 and the px=|p_perp|^2/2 relation
                // Note: the agreement in the phase space between envelope and standard laser simulation will be seen only after the passage of the ionizing laser
                new_electrons.momentum( 0, idNew ) += Aabs*Aabs/4. + p_perp*p_perp/2.;

            } else if (ellipticity==1.){ // circular polarization

                // extract a random angle between 0 and 2pi, and assign p_perp = eA
                double rand_times_2pi = patch->rand_->uniform_2pi(); // from uniform distribution between [0,2pi]
                
                Aabs    = sqrt(2. * (*(Phi_env+ipart-ipart_ref))  );                 

                p_perp = Aabs;   // in circular polarization it corresponds to a0/sqrt(2)
                new_electrons.momentum( 1, idNew ) += p_perp*cos(rand_times_2pi)/sqrt(2);
                new_electrons.momentum( 2, idNew ) += p_perp*sin(rand_times_2pi)/sqrt(2); 
     
                // initialize px to take into account the average drift <px>=A^2/4 and the px=|p_perp|^2/2 result
                // Note: the agreement in the phase space between envelope and standard laser simulation will be seen only after the passage of the ionizing laser
                new_electrons.momentum( 0, idNew ) += Aabs*Aabs/2.; 
            
            }

            // weight and charge of the new electron
            new_electrons.weight( idNew )=double( k_times )*particles->weight( ipart );
            new_electrons.charge( idNew )=-1;
    
            // Increase the charge of the ion particle
            particles->charge( ipart ) += k_times;
        }
    
    
    } // Loop on particles

}

