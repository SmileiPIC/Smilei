#include "IonizationTunnel.h"

#include <cmath>

#include "Particles.h"
#include "Species.h"

using namespace std;



IonizationTunnel::IonizationTunnel(Params& params, Species * species) : Ionization(params, species) {
    DEBUG("Creating the Tunnel Ionizaton class");
    
    one_third = 1.0/3.0;
    
    alpha_tunnel.resize(atomic_number_);
    beta_tunnel.resize(atomic_number_);
    gamma_tunnel.resize(atomic_number_);
    
    for (unsigned int Z=0 ; Z<atomic_number_ ; Z++) {
        DEBUG("Z : " << Z);
        double cst      = ((double)Z+1.0) * sqrt(2.0/Potential[Z]);
        alpha_tunnel[Z] = cst-1.0;
        beta_tunnel[Z]  = pow(2,alpha_tunnel[Z]) * (8.*Azimuthal_quantum_number[Z]+4.0) / (cst*tgamma(cst)) * Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0*Potential[Z],1.5);
    }
    
    new_electrons.initialize(0, params.nDim_particle );
    DEBUG("Finished Creating the Tunnel Ionizaton class");

}



void IonizationTunnel::operator() (Particles &particles, int ipart, LocalFields Epart) {
    for (unsigned int i=0; i< particles.size(); i++) {
        particles.print(i);
    }
    
    // Charge state of the ion (particle)
    unsigned int Z = (unsigned int)(particles.charge(ipart));
    
    // Number of successive ionization in one time-step
    unsigned int k_times = 0;
    
    // Generate a random number between 0 and 1
    double ran_p = (double)rand() / RAND_MAX;
    
    // Absolute value of the electric field normalized in atomic units
    double E = EC_to_au * sqrt( Epart.x*Epart.x + Epart.y*Epart.y + Epart.z*Epart.z );
    
    // Ionization rate in normalized (SMILEI) units
    double delta = gamma_tunnel[Z] / E;
    vector<double> IonizRate_tunnel(atomic_number_);
    IonizRate_tunnel[Z] = beta_tunnel[Z] * pow(delta,alpha_tunnel[Z]) * exp(-delta*one_third);
    
    // --------------------------------
    // Start of the Monte-Carlo routine
    // --------------------------------
    if (E!=0 && Z!=atomic_number_) {
        
        // if ionization of the last electron: single ionization
        // -----------------------------------------------------
        if ( Z == atomic_number_-1) {
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
            new_electrons.create_particle();
            //new_electrons.initialize( new_electrons.size()+1, new_electrons.dimension() );
            int idNew = new_electrons.size() - 1;
            for (unsigned int i=0; i<new_electrons.dimension(); i++) {
                new_electrons.position(i,idNew)=particles.position(i, ipart);
                new_electrons.position_old(i,idNew)=particles.position(i, ipart);
            }
            for (unsigned int i=0; i<3; i++) {
                new_electrons.momentum(i,idNew) = particles.momentum(i, ipart)/ionized_species_mass;
            }
            new_electrons.weight(idNew)=double(k_times)*particles.weight(ipart);
            new_electrons.charge(idNew)=-1;
            
            // Increase the charge of the particle
            particles.charge(ipart) += k_times;
        }
        
    }//END test on electric field && Z_atomic

}



void IonizationTunnel::operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Jion) {
    
    // Total ionization potential (used to compute the ionization current)
    double TotalIonizPot = 0.0;
    
    // Charge state of the ion (particle)
    unsigned int Z = (unsigned int)(particles.charge(ipart));
    
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
    
        if ( Z == atomic_number_-1) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if ( ran_p < 1.0 -exp(-IonizRate_tunnel[Z]*dt) ) {
                k_times        = 1;
            }
            
        } else {
            // else : multiple ionization can occur in one time-step
            //        partial & final ionization are decoupled (see Nuter Phys. Plasmas)
            // -------------------------------------------------------------------------
            
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
            
            // final ionization (of last electron)
            if ( ((1.0-Pint_tunnel)>ran_p) && (k_times==atomic_number_-Z-1) ) {
                k_times++;
            }
        }//END Multiple ionization routine
        
        
        // Calculation of the ionization current
        // -------------------------------------
        for (unsigned int k=0; k<k_times; k++) {
            TotalIonizPot += Potential[Z+k];
        }//END for k
        
        double factorJion = TotalIonizPot/dt/E/E;
        Jion.x = factorJion * Epart.x;
        Jion.y = factorJion * Epart.y;
        Jion.z = factorJion * Epart.z;
        
        
        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------
        if (k_times !=0) {
            new_electrons.create_particle();
            //new_electrons.initialize( new_electrons.size()+1, new_electrons.dimension() );
            int idNew = new_electrons.size() - 1;
            for (unsigned int i=0; i<new_electrons.dimension(); i++) {
                new_electrons.position(i,idNew)=particles.position(i, ipart);
                new_electrons.position_old(i,idNew)=particles.position(i, ipart);
            }
            for (unsigned int i=0; i<3; i++) {
                new_electrons.momentum(i,idNew) = particles.momentum(i, ipart)/ionized_species_mass;
            }
            new_electrons.weight(idNew)=double(k_times)*particles.weight(ipart);
            new_electrons.charge(idNew)=-1;
            
            // Increase the charge of the particle
            particles.charge(ipart) += k_times;
        }
        
    }//END test on electric field && Z_atomic

}
