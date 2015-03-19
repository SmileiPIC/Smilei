#include "Collisions.h"

#include <cmath>
#include <iomanip>
#include <ostream>

using namespace std;


// Constructor
Collisions::Collisions(PicParams& param, unsigned int ncol, vector<unsigned int> sgroup1, vector<unsigned int>sgroup2, double clog, bool intra)
{

    n_collisions     = ncol;
    species_group1   = sgroup1;
    species_group2   = sgroup2;
    coulomb_log      = clog;
    intra_collisions = intra;
    
}


// Declare static variables here
bool               Collisions::debye_length_required;
vector<double>     Collisions::debye_length_squared;



// Calculates the debye length squared in each cluster
// The formula for the inverse debye length squared is sumOverSpecies(density*charge^2/temperature)
void Collisions::calculate_debye_length(PicParams& params, vector<Species*>& vecSpecies)
{

    // get info on particle binning
    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    unsigned int nspec = vecSpecies.size(); // number of species
    unsigned int bmin, bmax;
    double p2, density, density_max, charge, temperature, rmin2;
    Species   * s;
    Particles * p;
    double coeff = params.wavelength_SI/(6.*M_PI*2.8179403267e-15); // normLength/(3*electronRadius) = wavelength/(6*pi*electronRadius)
    
    debye_length_squared.resize(nbins);
    
    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {
        
        density_max = 0.;
        debye_length_squared[ibin] = 0.;
        for (unsigned int ispec=0 ; ispec<nspec ; ispec++) { // loop all species
            // Calculation of particles density, mean charge, and temperature
            // Density is the sum of weights
            // Temperature basic definition is the average <v*p> divided by 3
            //    (instead of v*p, we use p^2/gamma)
            s  = vecSpecies[ispec];
            p  = &(s->particles);
            bmin = s->bmin[ibin];
            bmax = s->bmax[ibin];
            density     = 0.;
            charge      = 0.;
            temperature = 0.;
            // loop particles to calculate average quantities
            for (unsigned int iPart=bmin ; iPart<bmax ; iPart++ ) {
                p2 = pow(p->momentum(0,iPart),2)+pow(p->momentum(1,iPart),2)+pow(p->momentum(2,iPart),2);
                density     += p->weight(iPart);
                charge      += p->weight(iPart) * p->charge(iPart);
                temperature += p->weight(iPart) * p2/sqrt(1.+p2);
            }
            if (density <= 0.) continue;
            charge /= density; // average charge
            temperature *= (s->species_param.mass) / (3.*density); // Te in units of me*c^2
            density /= params.n_cell_per_cluster; // density in units of critical density
            // compute inverse debye length squared
            if (temperature>0.) debye_length_squared[ibin] += density*charge*charge/temperature;
            // compute maximum density of species
            if (density>density_max) density_max = density;
        }
        
        // if there were particles, 
        if (debye_length_squared[ibin] > 0.) {
            // compute debye length squared in code units
            debye_length_squared[ibin] = 1./(debye_length_squared[ibin]);
            // apply lower limit to the debye length (minimum interatomic distance)
            rmin2 = pow(coeff*density_max, -2./3.);
            if (debye_length_squared[ibin] < rmin2) debye_length_squared[ibin] = rmin2;
        }
        
    }
    
#ifdef  __DEBUG
    // calculate and print average debye length
    double mean_debye_length = 0.;
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++)
        mean_debye_length += sqrt(debye_length_squared[ibin]);
    mean_debye_length /= (double)nbins;
    //DEBUG("Mean Debye length in code length units = " << scientific << setprecision(3) << mean_debye_length);
    mean_debye_length *= params.wavelength_SI/(2.*M_PI); // switch to SI
    DEBUG("Mean Debye length in meters = " << scientific << setprecision(3) << mean_debye_length );
#endif

}

