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
        double cst      = (double)Z * sqrt(2.0/Potential[Z]);
        alpha_tunnel[Z] = cst-1.0;
        if (Z==0) {
            beta_tunnel[Z]  = pow(2,cst) * (4.0*Azimuthal_quantum_number[Z]+2.0) * Potential[Z] * au_to_w0;
        }
        else {
            beta_tunnel[Z]  = pow(2,cst) * (4.0*Azimuthal_quantum_number[Z]+2.0) / cst / gamma(cst) * Potential[Z] * au_to_w0;
        }
        
        gamma_tunnel[Z] = 2.0 * pow(2.0*Potential[Z],1.5);
    }

    DEBUG("Finished Creating the Tunnel Ionizaton class");
    
}



void IonizationTunnel::operator() (Particle* part, LocalFields Epart) {

    
    // Charge state of the ion (particle)
    unsigned int Z = (unsigned int)(part->charge());
    
	// Absolute value of the electric field normalized in atomic units
    double E = EC_to_au * sqrt( Epart.x*Epart.x + Epart.y*Epart.y + Epart.z*Epart.z );
    //MESSAGE( " >>> Ex, Ey, Ez, E = " << Epart.x << " " << Epart.y << " " << Epart.z << " " << EC_to_au );
    
//    if (std::isfinite(E)) {
    if (E!=0) {
        // Calculation of the ionization rate in normalized (SMILEI) units
        double delta     = gamma_tunnel[Z] / E;
        double IonizRate = beta_tunnel[Z] * pow(delta,alpha_tunnel[Z]) * exp(-delta*one_third);
        //MESSAGE(" >>>> alpha , beta, delta = " << alpha_tunnel[Z] << " " << beta_tunnel[Z] << " " << delta);
        
        // Generate a random number between 0 and 1
        double ran_p = (double)rand() / RAND_MAX;
        //MESSAGE(" >>>> ran_p = " << ran_p << " cond. " << " Ionization rate = " << IonizRate << " dt = " << dt);
        
//		DEBUG(ran_p << " " << beta_tunnel[Z] << " " << alpha_tunnel[Z] << " " << gamma_tunnel[Z] << " " << E << " " << dt);
        
        // Calculate new charge state of the particle
        if ( ran_p < 1.0 -exp(-IonizRate*dt) ) {
            
            // Create the electron
            Particle* newParticle = new Particle(nDim_particle);
            for (unsigned int i=0;i<nDim_particle; i++) {
                newParticle->position(i)=part->position(i);
                newParticle->position_old(i)=part->position(i);
            }
            for (unsigned int i=0; i<3; i++){
                newParticle->momentum(i) = part->momentum(i)/ionized_species_mass;
            }
            newParticle->weight()=part->weight();
            //newParticle->charge()=part->charge();
            newParticle->charge()=1;
            new_electrons.push_back(newParticle);
            
            // Increase the charge of the particle
            part->charge() += 1;
        }// ENDIF ionization takes place
    }

}