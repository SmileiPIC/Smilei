
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>
#include <fstream>

#include "CollisionsSingle.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"
#include "Patch.h"
#include "VectorPatch.h"

using namespace std;


// Calculates the collisions for a given Collisions object
void CollisionsSingle::collide(Params& params, Patch* patch, int itime, vector<Diagnostic*>& localDiags)
{

    vector<unsigned int> index1, index2;
    unsigned int npairs; // number of pairs of macro-particles
    unsigned int np1, np2; // numbers of macro-particles in each species
    double n1, n2, n12, n123, n223; // densities of particles
    unsigned int i1, i2, N2max;
    Species   *s1, *s2, *stmp;
    Particles *p1, *p2;
    double m1, m2, m12, W1, W2, qqm, qqm2, gamma1, gamma2, gamma12, gamma12_inv,
           COM_vx, COM_vy, COM_vz, COM_vsquare, COM_gamma,
           term1, term2, term3, term4, term5, term6, coeff3, coeff4,
           vcv1, vcv2, px_COM, py_COM, pz_COM, p2_COM, p_COM, gamma1_COM, gamma2_COM,
           logL, bmin, s, vrel, smax,
           cosX, sinX, phi, sinXcosPhi, sinXsinPhi, p_perp, inv_p_perp,
           newpx_COM, newpy_COM, newpz_COM, U, vcp, ncol;
    bool not_duplicated_particle;
    
    s1 = patch->vecSpecies[species_group1[0]];
    s2 = patch->vecSpecies[species_group2[0]];
    
    bool debug = (debug_every > 0 && itime % debug_every == 0); // debug only every N timesteps
    
    if( debug ) {
        ncol = 0.;
        smean       = 0.;
        logLmean    = 0.;
        //temperature = 0.;
    }
    
    // Loop bins of particles (typically, cells, but may also be clusters)
    unsigned int nbin = patch->vecSpecies[0]->bmin.size();
    for (unsigned int ibin = 0 ; ibin < nbin ; ibin++) {
        
        // get number of particles for all necessary species
        for (unsigned int i=0; i<2; i++) { // try twice to ensure species 1 has more macro-particles
            np1 = s1->bmax[ibin] - s1->bmin[ibin];
            np2 = s2->bmax[ibin] - s2->bmin[ibin];
            if (np2 <= np1) break; // ok if species1 has more macro-particles
            else { // otherwise, we exchange species and try again
                stmp = s1; s1 = s2; s2 = stmp;
            }
        }
        
        // skip to next bin if no particles
        if (np1==0 || np2==0) continue;
        
        // Shuffle particles to have random pairs
        //    (It does not really exchange them, it is just a temporary re-indexing)
        index1.resize(np1);
        for (unsigned int i=0; i<np1; i++) index1[i] = i; // first, we make an ordered array
        random_shuffle(index1.begin(), index1.end()); // shuffle the index array
        if (intra_collisions) { // In the case of collisions within one species
            npairs = (int) ceil(((double)np1)/2.); // half as many pairs as macro-particles
            index2.resize(npairs);
            for (unsigned int i=0; i<npairs; i++) index2[i] = index1[(i+npairs)%np1]; // index2 is second half
            index1.resize(npairs); // index1 is first half
            N2max = np1 - npairs; // number of not-repeated particles (in group 2 only)
        } else { // In the case of collisions between two species
            npairs = np1; // as many pairs as macro-particles in species 1 (most numerous)
            index2.resize(npairs);
            for (unsigned int i=0; i<np1; i++) index2[i] = i % np2;
            N2max = np2; // number of not-repeated particles (in species 2 only)
        }
        
        // Offset the index array according to the current bins
        for (unsigned int i=0; i<npairs; i++)
            index1[i] += s1->bmin[ibin];
        for (unsigned int i=0; i<npairs; i++)
            index2[i] += s2->bmin[ibin];
         
        // Prepare the ionization
        Ionization->prepare1(s1->atomic_number);
        
        // Calculate the densities
        n1  = 0.; // density of group 1
        n2  = 0.; // density of group 2
        n12 = 0.; // "hybrid" density
        p1 = s1->particles;
        p2 = s2->particles;
        for (unsigned int i=0; i<npairs; i++) { // for each pair of particles
            // sum weights
            n1 += p1->weight(index1[i]);
            not_duplicated_particle = (i<N2max);
            if( not_duplicated_particle ) n2 += p2->weight(index2[i]); // special case for group 2 to avoid repeated particles
            n12 += min( p1->weight(index1[i]),  p2->weight(index2[i]) );
            // Same for ionization
            Ionization->prepare2(p1, index1[i], p2, index2[i], not_duplicated_particle);
        }
        if( intra_collisions ) { n1 += n2; n2 = n1; }
        n1  *= n_patch_per_cell;
        n2  *= n_patch_per_cell;
        n12 *= n_patch_per_cell;
        
        // Pre-calculate some numbers before the big loop
        n123 = pow(n1,2./3.);
        n223 = pow(n2,2./3.);
        coeff3 = params.timestep * n1*n2/n12;
        coeff4 = pow( 3.*coeff2 , -1./3. ) * coeff3;
        coeff3 *= coeff2;
        m1 = s1->mass;
        m2 = s2->mass;
        m12  = m1 / m2; // mass ratio
        
        // Prepare the ionization
        Ionization->prepare3(params.timestep, n_patch_per_cell);
        
        // Now start the real loop on pairs of particles
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------

        for (unsigned int i=0; i<npairs; i++) {
            i1 = index1[i];
            i2 = index2[i];
            W1 = p1->weight(i1);
            W2 = p2->weight(i2);
            
            // Calculate stuff
            qqm  = p1->charge(i1) * p2->charge(i2) / m1;
            qqm2 = qqm * qqm;
            
            // Get momenta and calculate gammas
            gamma1 = sqrt(1. + pow(p1->momentum(0,i1),2) + pow(p1->momentum(1,i1),2) + pow(p1->momentum(2,i1),2));
            gamma2 = sqrt(1. + pow(p2->momentum(0,i2),2) + pow(p2->momentum(1,i2),2) + pow(p2->momentum(2,i2),2));
            gamma12 = m12 * gamma1 + gamma2;
            gamma12_inv = 1./gamma12;
            
            // Calculate the center-of-mass (COM) frame
            // Quantities starting with "COM" are those of the COM itself, expressed in the lab frame.
            // They are NOT quantities relative to the COM.
            COM_vx = ( m12 * (p1->momentum(0,i1)) + p2->momentum(0,i2) ) * gamma12_inv;
            COM_vy = ( m12 * (p1->momentum(1,i1)) + p2->momentum(1,i2) ) * gamma12_inv;
            COM_vz = ( m12 * (p1->momentum(2,i1)) + p2->momentum(2,i2) ) * gamma12_inv;
            COM_vsquare = COM_vx*COM_vx + COM_vy*COM_vy + COM_vz*COM_vz;
            
            // Change the momentum to the COM frame (we work only on particle 1)
            // Quantities ending with "COM" are quantities of the particle expressed in the COM frame.
            if( COM_vsquare != 0.) {
                COM_gamma = pow( 1.-COM_vsquare , -0.5);
                term1 = (COM_gamma - 1.) / COM_vsquare;
                vcv1  = (COM_vx*(p1->momentum(0,i1)) + COM_vy*(p1->momentum(1,i1)) + COM_vz*(p1->momentum(2,i1)))/gamma1;
                vcv2  = (COM_vx*(p2->momentum(0,i2)) + COM_vy*(p2->momentum(1,i2)) + COM_vz*(p2->momentum(2,i2)))/gamma2;
                term2 = (term1*vcv1 - COM_gamma) * gamma1;
                px_COM = (p1->momentum(0,i1)) + term2*COM_vx;
                py_COM = (p1->momentum(1,i1)) + term2*COM_vy;
                pz_COM = (p1->momentum(2,i1)) + term2*COM_vz;
                gamma1_COM = (1.-vcv1)*COM_gamma*gamma1;
                gamma2_COM = (1.-vcv2)*COM_gamma*gamma2;
            } else {
                COM_gamma = 1.;
                term1 = 0.5;
                term2 = gamma1;
                px_COM = (p1->momentum(0,i1));
                py_COM = (p1->momentum(1,i1));
                pz_COM = (p1->momentum(2,i1));
                gamma1_COM = gamma1;
                gamma2_COM = gamma2;
            }
            p2_COM = px_COM*px_COM + py_COM*py_COM + pz_COM*pz_COM;
            p_COM  = sqrt(p2_COM);
            
            // Calculate some intermediate quantities
            term3 = COM_gamma * gamma12_inv;
            term4 = gamma1_COM * gamma2_COM;
            term5 = term4/p2_COM + m12;
            
            // Calculate coulomb log if necessary
            logL = coulomb_log;
            if( logL <= 0. ) { // if auto-calculation requested
                bmin = max( coeff1/m1/p_COM , abs(coeff2*qqm*term3*term5) ); // min impact parameter
                logL = 0.5*log(1.+patch->debye_length_squared[ibin]/pow(bmin,2));
                if (logL < 2.) logL = 2.;
            }
            
            // Calculate the collision parameter s12 (similar to number of real collisions)
            s = coeff3 * logL * qqm2 * term3 * p_COM * term5*term5 / (gamma1*gamma2);
            
            // Low-temperature correction
            vrel = p_COM/term3/term4; // relative velocity
            smax = coeff4 * (m12+1.) * vrel / max(m12*n123,n223);
            if (s>smax) s = smax;
            
            // Pick the deflection angles according to Nanbu's theory
            cosX = cos_chi(s);
            sinX = sqrt( 1. - cosX*cosX );
            //!\todo make a faster rand by preallocating ??
            phi = twoPi * Rand::uniform();
            
            // Calculate combination of angles
            sinXcosPhi = sinX*cos(phi);
            sinXsinPhi = sinX*sin(phi);
            
            // Apply the deflection
            p_perp = sqrt( px_COM*px_COM + py_COM*py_COM );
            if( p_perp > 1.e-10*p_COM ) { // make sure p_perp is not too small
                inv_p_perp = 1./p_perp;
                newpx_COM = (px_COM * pz_COM * sinXcosPhi - py_COM * p_COM * sinXsinPhi) * inv_p_perp + px_COM * cosX;
                newpy_COM = (py_COM * pz_COM * sinXcosPhi + px_COM * p_COM * sinXsinPhi) * inv_p_perp + py_COM * cosX;
                newpz_COM = -p_perp * sinXcosPhi  +  pz_COM * cosX;
            } else { // if p_perp is too small, we use the limit px->0, py=0
                newpx_COM = p_COM * sinXcosPhi;
                newpy_COM = p_COM * sinXsinPhi;
                newpz_COM = p_COM * cosX;
            }
            
            // Random number to choose whether deflection actually applies.
            // This is to conserve energy in average when weights are not equal.
            //!\todo make a faster rand by preallocating ??
            U = Rand::uniform();
            
            // Go back to the lab frame and store the results in the particle array
            vcp = COM_vx * newpx_COM + COM_vy * newpy_COM + COM_vz * newpz_COM;
            if( U < W2/W1 ) { // deflect particle 1 only with some probability
                term6 = term1*vcp + gamma1_COM * COM_gamma;
                p1->momentum(0,i1) = newpx_COM + COM_vx * term6;
                p1->momentum(1,i1) = newpy_COM + COM_vy * term6;
                p1->momentum(2,i1) = newpz_COM + COM_vz * term6;
            }
            if( U < W1/W2 ) { // deflect particle 2 only with some probability
                term6 = -m12 * term1*vcp + gamma2_COM * COM_gamma;
                p2->momentum(0,i2) = -m12 * newpx_COM + COM_vx * term6;
                p2->momentum(1,i2) = -m12 * newpy_COM + COM_vy * term6;
                p2->momentum(2,i2) = -m12 * newpz_COM + COM_vz * term6;
            }
            
            // Handle ionization
            Ionization->apply(p1, i1, p2, i2);
            
            if( debug ) {
                ncol     += npairs;
                smean    += s;
                logLmean += logL;
                //temperature += m1 * (sqrt(1.+pow(p1->momentum(0,i1),2)+pow(p1->momentum(1,i1),2)+pow(p1->momentum(2,i1),2))-1.);
            }
            
        } // end loop on pairs of particles
        
        Ionization->finish(s1, s2, params, patch, localDiags);
    }
    
    if(debug && ncol>0. ) {
        smean    /= ncol;
        logLmean /= ncol;
        //temperature /= ncol;
    }
}




