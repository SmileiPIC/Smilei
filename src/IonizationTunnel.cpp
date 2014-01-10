#include "IonizationTunnel.h"
#include "Particle.h"
#include <cmath>

using namespace std;



IonizationTunnel::IonizationTunnel(PicParams *params, int ispec) : Ionization(params, ispec) {
	DEBUG("Creating the Tunnel Ionizaton class");
    
    one_third = 1.0/3.0;
    
    alpha_tunnel.resize(atomic_number_);
    beta_tunnel.resize(atomic_number_);
    gamma_tunnel.resize(atomic_number_);
    
    for (unsigned int Z=0 ; Z<atomic_number_ ; Z++){
        DEBUG("Z : " << Z);
        double cst      = ((double)Z+1.0) * sqrt(2.0/Potential[Z]);
        alpha_tunnel[Z] = cst-1.0;
        beta_tunnel[Z]  = pow(2,cst) * (4.0*Azimuthal_quantum_number[Z]+2.0) / tgamma(cst+1.0) * Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0*Potential[Z],1.5);
    }
    
    DEBUG("Finished Creating the Tunnel Ionizaton class");
    
}



void IonizationTunnel::operator() (Particle* part, LocalFields Epart) {
    // Charge state of the ion (particle)
    unsigned int Z = (unsigned int)(part->charge());
    
    // Number of successive ionization in one time-step
    unsigned int k_times = 0;
    
    // Generate a random number between 0 and 1
    double ran_p = (double)rand() / RAND_MAX;
    
	// Absolute value of the electric field normalized in atomic units
    double E = EC_to_au * sqrt( Epart.x*Epart.x + Epart.y*Epart.y + Epart.z*Epart.z );
	 
    // Ionization rate in normalized (SMILEI) units
    double delta     = gamma_tunnel[Z] / E;
    vector<double> IonizRate_tunnel(atomic_number_);
    IonizRate_tunnel[Z] = beta_tunnel[Z] * pow(delta,alpha_tunnel[Z]) * exp(-delta*one_third);

    // --------------------------------
	// Start of the Monte-Carlo routine
	// --------------------------------
    if (E!=0 && Z!=atomic_number_) {

		// if ionization of the last electron: single ionization
		// -----------------------------------------------------
        if ( Z == atomic_number_-1){
        	if ( ran_p < 1.0 -exp(-IonizRate_tunnel[Z]*dt) ) {
            	k_times = 1;
        	}
        	
            // else : multiple ionization can occur in one time-step
            //        partial & final ionization are decoupled (see Nuter Phys. Plasmas)
            // -------------------------------------------------------------------------
        } else {
            // initialization
        	double Mult = 1.0;
            vector<double> Dnom_tunnel(atomic_number_);
            Dnom_tunnel[0]=1.0;
        	double Pint_tunnel = exp(-IonizRate_tunnel[Z]*dt); // cummulative prob.
        	
        	//multiple ionization loop while Pint_tunnel < ran_p and still partial ionization
       		while ((Pint_tunnel < ran_p) and (k_times < atomic_number_-Z-1)) {
                unsigned int newZ = Z+k_times+1;
            	double Prob  = 0.0;
                delta = gamma_tunnel[newZ] / E;
            	IonizRate_tunnel[newZ] = beta_tunnel[newZ]
                *                        pow(delta,alpha_tunnel[newZ])
                *                        exp(-delta*one_third);
	            double D_sum = 0.0;
	            double P_sum = 0.0;
	            Mult  *= IonizRate_tunnel[Z+k_times];
            	for (unsigned int i=0; i<k_times+1; i++) {
                	Dnom_tunnel[i]=Dnom_tunnel[i]/(IonizRate_tunnel[newZ]-IonizRate_tunnel[Z+i]);
                	D_sum += Dnom_tunnel[i];
                	P_sum += exp(-IonizRate_tunnel[Z+i]*dt)*Dnom_tunnel[i];
            	}
            	Dnom_tunnel[k_times+1] = -D_sum;
            	P_sum                  = P_sum + Dnom_tunnel[k_times+1]*exp(-IonizRate_tunnel[newZ]*dt);
            	Prob                   = P_sum * Mult;
            	Pint_tunnel            = Pint_tunnel+Prob;
            	
            	k_times++;
        	}//END while
        	
        	// final ionization
        	if ( ((1.0-Pint_tunnel)>ran_p) && (k_times==atomic_number_-Z-1) ) {
            	k_times++;
        	}
        }//END Multiple ionization routine
        
        
        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------
        if (k_times !=0) {
            Particle* newParticle = new Particle(nDim_particle);
            for (unsigned int i=0;i<nDim_particle; i++) {
                newParticle->position(i)=part->position(i);
                newParticle->position_old(i)=part->position(i);
            }
            for (unsigned int i=0; i<3; i++){
                newParticle->momentum(i) = part->momentum(i)/ionized_species_mass;
            }
            newParticle->weight()=double(k_times)*part->weight();
            newParticle->charge()=-1;
            new_electrons.push_back(newParticle);
            
            // Increase the charge of the particle
            part->charge() += k_times;
        }
        
    }//END test on electric field && Z_atomic
    
}