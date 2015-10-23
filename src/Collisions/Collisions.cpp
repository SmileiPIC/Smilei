#include "Collisions.h"
#include "SmileiMPI.h"
#include "Field2D.h"
#include "H5.h"


#include <cmath>
#include <iomanip>
#include <algorithm>
#include <ostream>

using namespace std;


// Constructor
Collisions::Collisions(SmileiMPI* smpi,
                       unsigned int n_collisions,
                       vector<unsigned int> species_group1, 
                       vector<unsigned int> species_group2, 
                       double coulomb_log, 
                       bool intra_collisions,
                       int debug_every,
                       unsigned int nbins,
                       int Z ) :
n_collisions    (n_collisions    ),
species_group1  (species_group1  ),
species_group2  (species_group2  ),
coulomb_log     (coulomb_log     ),
intra_collisions(intra_collisions),
debug_every     (debug_every     ),
start           (0               ),
atomic_number   (Z               )
{
    
    // Calculate total number of bins
    totbins = nbins;
    MPI_Allreduce( smpi->isMaster()?MPI_IN_PLACE:&totbins, &totbins, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    
    // if debug requested, prepare hdf5 file
    fileId = 0;
    if( debug_every>0 ) {
        ostringstream mystream;
        mystream.str("");
        mystream << "Collisions" << n_collisions << ".h5";
        // Create the HDF5 file 
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId = H5Fcreate(mystream.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        // write all parameters as HDF5 attributes
        H5::attr(fileId, "Version", string(__VERSION));
        mystream.str("");
        mystream << species_group1[0];
        for(unsigned int i=1; i<species_group1.size(); i++) mystream << "," << species_group1[i];
        H5::attr(fileId, "species1" , mystream.str());
        mystream.str("");
        mystream << species_group2[0];
        for(unsigned int i=1; i<species_group2.size(); i++) mystream << "," << species_group2[i];
        H5::attr(fileId, "species2" , mystream.str());
        H5::attr(fileId, "coulomb_log" , coulomb_log);
        H5::attr(fileId, "debug_every"  , debug_every);
        
        // Find out where the proc will start writing in the overall array
        MPI_Status status;
        // Receive the location where to start from the previous node
        if (smpi->getRank()>0) MPI_Recv( &(start), 1, MPI_INTEGER, smpi->getRank()-1, 0, MPI_COMM_WORLD, &status );
        // Send the location where to end to the next node
        int end = start+nbins;
        if (smpi->getRank()!=smpi->getSize()-1) MPI_Send( &end, 1, MPI_INTEGER, smpi->getRank()+1, 0, MPI_COMM_WORLD );
    }
}

Collisions::~Collisions()
{
    if (fileId != 0) H5Fclose(fileId);
}

// Reads the input file and creates the Collisions objects accordingly
vector<Collisions*> Collisions::create(Params& params, vector<Species*>& vecSpecies, SmileiMPI* smpi)
{
    vector<Collisions*> vecCollisions;
    
    vector<string> sg1, sg2;
    vector<vector<unsigned int>> sgroup;
    double clog;
    bool intra, debye_length_required = false, ionizing;
    int debug_every, Z, Z0, Z1;
    ostringstream mystream;
    Species *s0, *s;
    
    // Needs wavelength_SI to be defined
    unsigned int numcollisions=PyTools::nComponents("Collisions");
    if (numcollisions > 0)
        if (params.wavelength_SI <= 0.)
            ERROR("The parameter `wavelength_SI` needs to be defined and positive in order to compute collisions");
    
    // Loop over each binary collisions group and parse info
    for (unsigned int n_collisions = 0; n_collisions < numcollisions; n_collisions++) {
        
        MESSAGE("Parameters for collisions #" << n_collisions << " :");
        
        // Read the input file by searching for the keywords "species1" and "species2"
        // which are the names of the two species that will collide
        sg1.resize(0);
        sg2.resize(0);
        PyTools::extract("species1",sg1,"Collisions",n_collisions);
        PyTools::extract("species2",sg2,"Collisions",n_collisions);
        
        // Obtain the lists of species numbers from the lists of species names.
        sgroup.resize(2);
        sgroup[0] = params.FindSpecies(vecSpecies, sg1);
        sgroup[1] = params.FindSpecies(vecSpecies, sg2);
        
        // Each group of species sgroup[0] and sgroup[1] must not be empty
        if (sgroup[0].size()==0) ERROR("In collisions #" << n_collisions << ": No valid `species1`");
        if (sgroup[1].size()==0) ERROR("In collisions #" << n_collisions << ": No valid `species2`");
        
        // sgroup[0] and sgroup[1] can be equal, but cannot have common species if they are not equal
        if (sgroup[0] != sgroup[1]) {
            for (unsigned int i0=0; i0<sgroup[0].size(); i0++) {
                for (unsigned int i1=0; i1<sgroup[1].size(); i1++) {
                    if (sgroup[0][i0] == sgroup[1][i1])
                        ERROR("In collisions #" << n_collisions << ": species #" << sgroup[0][i0]
                            << " cannot collide with itself");
                }
            }
            intra = false;
        } else {
            intra = true;
        }
        
        // Coulomb logarithm (if negative or unset, then automatically computed)
        clog = 0.; // default
        PyTools::extract("coulomb_log",clog,"Collisions",n_collisions);
        if (clog <= 0.) debye_length_required = true; // auto coulomb log requires debye length
        
        // Number of timesteps between each debug output (if 0 or unset, no debug)
        debug_every = 0; // default
        PyTools::extract("debug_every",debug_every,"Collisions",n_collisions);
        
        // Collisional ionization
        ionizing = false; Z = 0; // default
        PyTools::extract("ionizing",ionizing,"Collisions",n_collisions);
        if( ionizing ) {
            if( intra )
                ERROR("In collisions #" << n_collisions << ": cannot ionize with intra-collisions");
            
            for (int g=0; g<2; g++) { // do sgroup[0], then sgroup[1]
                s0 = vecSpecies[sgroup[g][0]]; // first species of this group
                for (unsigned int i=1; i<sgroup[g].size(); i++) { // loop other species of same group
                    s = vecSpecies[sgroup[g][i]]; // current species
                    if( s->mass != s0->mass )
                        ERROR("In collisions #" << n_collisions << ": species in group `species"
                            << g+1 << "` must all have same masses for ionization");
                    
                    if( s->atomic_number != s0->atomic_number ) {
                        if( s->atomic_number * s0->atomic_number ==0 ) {
                            ERROR("In collisions #" << n_collisions << ": species in group `species"
                            << g+1 << "` cannot be mixed electrons and ions for ionization");
                        } else {
                            ERROR("In collisions #" << n_collisions << ": species in group `species"
                            << g+1 << "` must all have same atomic_number for ionization");
                        }
                    }
                }
            }
            // atomic number
            Z0 = vecSpecies[sgroup[0][0]]->atomic_number;
            Z1 = vecSpecies[sgroup[1][0]]->atomic_number;
            Z = (int) ( Z0>Z1 ? Z0 : Z1 );
            if( Z0*Z1!=0 )
                ERROR("In collisions #" << n_collisions << ": ionization requires electrons (no or null atomic_number)");
            if( Z==0 )
                ERROR("In collisions #" << n_collisions << ": ionization requires ions (atomic_number>0)");
        }
        
        // Print collisions parameters
        mystream.str(""); // clear
        mystream << "(" << sgroup[0][0];
        for (unsigned int rs=1 ; rs<sgroup[0].size() ; rs++) mystream << " " << sgroup[0][rs];
        if( intra ) {
            MESSAGE(1,"Intra collisions within species " << mystream.str() << ")");
        } else {
            mystream << ") and (" << sgroup[1][0];
            for (unsigned int rs=1 ; rs<sgroup[1].size() ; rs++) mystream << " " << sgroup[1][rs];
            MESSAGE(1,"Collisions between species " << mystream.str() << ")");
        }
        MESSAGE(1,"Coulomb logarithm: " << clog);
        if( debug_every>0 ) MESSAGE(1,"Debug every " << debug_every << " timesteps");
        mystream.str(""); // clear
        if( ionizing>0 ) MESSAGE(1,"Collisional ionization with atomic number "<<Z);
        
        // Add new Collisions objects to vector
        vecCollisions.push_back( new Collisions(smpi,n_collisions,sgroup[0],sgroup[1],clog,intra,
            debug_every, vecSpecies[0]->bmin.size(), Z) );
        
        // Create the ionization object
        if( ionizing ) {
            vecCollisions[n_collisions]->Ionization = NULL;
            // But do not duplicate if ionization already existing with same Z
            for( int i=0; i<n_collisions; i++ ) {
                if( Z == vecCollisions[i]->atomic_number ) {
                    vecCollisions[n_collisions]->Ionization = vecCollisions[i]->Ionization;
                    break;
                }
            }
            // Otherwise, create new
            if( vecCollisions[n_collisions]->Ionization == NULL){
                vecCollisions[n_collisions]->Ionization = new CollisionalIonization(Z, params.wavelength_SI, smpi);
                vecCollisions[n_collisions]->Ionization->new_electrons.initialize(0, params.nDim_particle );// temporary TO BE REMOVED
            }
        } else {
            // If no ionization, create 'empty' ionization object
            vecCollisions[n_collisions]->Ionization = new CollisionalNoIonization();
        }
    }
    
    // pass the variable "debye_length_required" into the Collision class
    Collisions::debye_length_required = debye_length_required;
    
    return vecCollisions;
}


// Declare other static variables here
bool               Collisions::debye_length_required;
vector<double>     Collisions::debye_length_squared;

// Calculates the debye length squared in each cluster
// The formula for the inverse debye length squared is sumOverSpecies(density*charge^2/temperature)
void Collisions::calculate_debye_length(Params& params, vector<Species*>& vecSpecies)
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
    //! \todo Make OpenMP parallelization (MF & JD)
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
            temperature *= (s->mass) / (3.*density); // Te in units of me*c^2
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


// Calculates the collisions for a given Collisions object
void Collisions::collide(Params& params, vector<Species*>& vecSpecies, int itime)
{

    unsigned int nbins = vecSpecies[0]->bmin.size(); // number of bins
    vector<unsigned int> *sg1, *sg2, *sgtmp, bmin1, bmax1, bmin2, bmax2, index1, index2;
    unsigned int nspec1, nspec2; // numbers of species in each group
    unsigned int npart1, npart2; // numbers of macro-particles in each group
    unsigned int npairs; // number of pairs of macro-particles
    vector<unsigned int> np1, np2; // numbers of macro-particles in each species, in each group
    double n1, n2, n12, n123, n223; // densities of particles
    unsigned int i1, i2, ispec1, ispec2;
    Species   *s1, *s2;
    Particles *p1, *p2;
    double m1, m2, m12, W1, W2, qqm, qqm2, gamma1, gamma2, gamma12, gamma12_inv,
           COM_vx, COM_vy, COM_vz, COM_vsquare, COM_gamma,
           term1, term2, term3, term4, term5, term6, coeff1, coeff2, coeff3, coeff4, twoPi,
           vcv1, vcv2, px_COM, py_COM, pz_COM, p2_COM, p_COM, gamma1_COM, gamma2_COM,
           logL, bmin, s, vrel, smax,
           cosX, sinX, phi, sinXcosPhi, sinXsinPhi, p_perp, inv_p_perp, 
           newpx_COM, newpy_COM, newpz_COM, U, vcp;
    Field2D *smean=NULL, *logLmean=NULL, *ncol=NULL;//, *temperature
    ostringstream name;
    hid_t did(0);
    
    sg1 = &species_group1;
    sg2 = &species_group2;
    
    
    bool debug = (debug_every > 0 && itime % debug_every == 0); // debug only every N timesteps
    
    if( debug ) {
        // Create H5 group for the current timestep
        name.str("");
        name << "t" << setfill('0') << setw(8) << itime;
        did = H5::group(fileId, name.str());
        // Prepare storage arrays
        vector<unsigned int> outsize(2); outsize[0]=nbins; outsize[1]=1;
        smean       = new Field2D(outsize);
        logLmean    = new Field2D(outsize);
        //temperature = new Field2D(outsize);
        ncol        = new Field2D(outsize);
    }
    
    // Loop on bins
    for (unsigned int ibin=0 ; ibin<nbins ; ibin++) {
        
        // get bin start/end for all necessary species, and number of particles
        for (unsigned int i=0; i<2; i++) { // try twice to ensure group 1 has more macro-particles
            nspec1 = sg1->size();
            nspec2 = sg2->size();
            bmin1.resize(nspec1); // bin starting point, for each of species group 1
            bmax1.resize(nspec1); // bin  ending  point, for each of species group 1
            np1  .resize(nspec1); // number of particles in that bin
            bmin2.resize(nspec2); // bin starting point, for each of species group 2
            bmax2.resize(nspec2); // bin  ending  point, for each of species group 2
            np2  .resize(nspec2); // number of particles in that bin
            npart1 = 0; npart2 = 0;
            for (ispec1=0 ; ispec1<nspec1 ; ispec1++) {
                s1 = vecSpecies[(*sg1)[ispec1]];
                bmin1[ispec1] = s1->bmin[ibin];
                bmax1[ispec1] = s1->bmax[ibin];
                np1[ispec1] = bmax1[ispec1] - bmin1[ispec1];
                npart1 += np1[ispec1];
            }
            for (ispec2=0 ; ispec2<nspec2 ; ispec2++) {
                s2 = vecSpecies[(*sg2)[ispec2]];
                bmin2[ispec2] = s2->bmin[ibin];
                bmax2[ispec2] = s2->bmax[ibin];
                np2[ispec2] = bmax2[ispec2] - bmin2[ispec2];
                npart2 += np2[ispec2];
            }
            if (npart2 <= npart1) break; // ok if group1 has more macro-particles
            else { // otherwise, we exchange groups and try again
                sgtmp = sg1; sg1 = sg2; sg2 = sgtmp;
            }
        }
        // now group1 has more macro-particles than group2
        
        // skip to next bin if no particles
        if (npart1==0 or npart2==0) continue;
        
        // Shuffle particles to have random pairs
        //    (It does not really exchange them, it is just a temporary re-indexing)
        index1.resize(npart1);
        for (unsigned int i=0; i<npart1; i++) index1[i] = i; // first, we make an ordered array
        //! \todo benchmark and improve the shuffling method ?
        random_shuffle(index1.begin(), index1.end()); // shuffle the index array
        if (intra_collisions) { // In the case of collisions within one species
            npairs = (int) ceil(((double)npart1)/2.); // half as many pairs as macro-particles
            index2.resize(npairs);
            for (unsigned int i=0; i<npairs; i++) index2[i] = index1[i+npart1-npairs]; // index2 is second half
            index1.resize(npairs); // index1 is first half
        } else { // In the case of collisions between two species
            npairs = npart1; // as many pairs as macro-particles in group 1 (most numerous)
            index2.resize(npairs);
            for (unsigned int i=0; i<npart1; i++) index2[i] = i % npart2;
        }
        
        // Prepare the ionization
        Ionization->prepare1(vecSpecies[(*sg1)[0]]->atomic_number);
        
        // Calculate the densities
        n1  = 0.; // density of group 1
        n2  = 0.; // density of group 2
        n12 = 0.; // "hybrid" density
        for (unsigned int i=0; i<npairs; i++) { // for each pair of particles
            // find species and index i1 of particle "1"
            i1 = index1[i];
            for (ispec1=0 ; i1>=np1[ispec1]; ispec1++) i1 -= np1[ispec1];
            i1 += bmin1[ispec1];
            // find species and index i2 of particle "2"
            i2 = index2[i];
            for (ispec2=0 ; i2>=np2[ispec2]; ispec2++) i2 -= np2[ispec2];
            i2 += bmin2[ispec2];
            // Pointers to particles
            p1 = &(vecSpecies[(*sg1)[ispec1]]->particles);
            p2 = &(vecSpecies[(*sg2)[ispec2]]->particles);
            // sum weights
            n1  += p1->weight(i1);
            n2  += p2->weight(i2);
            n12 += min( p1->weight(i1), p2->weight(i2) );
            // Same for ionization
            Ionization->prepare2(p1, i1, p2, i2);
        }
        if( intra_collisions ) { n1 += n2; n2 = n1; }
        n1  /= params.n_cell_per_cluster;
        n2  /= params.n_cell_per_cluster;
        n12 /= params.n_cell_per_cluster;
        
        // Pre-calculate some numbers before the big loop
        n123 = pow(n1,2./3.); n223 = pow(n2,2./3.);
        twoPi = 2. * M_PI;
        coeff1 = M_PI*6.62606957e-34/(9.10938215e-31*299792458.*params.wavelength_SI); // h/(2*me*c*normLength) = pi*h/(me*c*wavelength)
        coeff2 = 2.*M_PI*2.817940327e-15/params.wavelength_SI; // re/normLength = 2*pi*re/wavelength
        coeff3 = coeff2 * params.timestep * n1*n2/n12;
        coeff4 = pow( 3.*coeff2 , -1./3. ) * params.timestep * n1*n2/n12;
        
        // Prepare the ionization
        Ionization->prepare3(params.timestep, params.n_cell_per_cluster, npart1, npart2);
        
        if( debug ) {
            smean      ->data_2D[ibin][0] = 0.;
            logLmean   ->data_2D[ibin][0] = 0.;
            //temperature->data_2D[ibin][0] = 0.;
            ncol       ->data_2D[ibin][0] = (double)npairs;
        }
        
        
        // Now start the real loop on pairs of particles
        // See equations in http://dx.doi.org/10.1063/1.4742167
        // ----------------------------------------------------
        for (unsigned int i=0; i<npairs; i++) {
        
            // find species and index i1 of particle "1"
            i1 = index1[i];
            for (ispec1=0 ; i1>=np1[ispec1]; ispec1++) i1 -= np1[ispec1];
            i1 += bmin1[ispec1];
            // find species and index i2 of particle "2"
            i2 = index2[i];
            for (ispec2=0 ; i2>=np2[ispec2]; ispec2++) i2 -= np2[ispec2];
            i2 += bmin2[ispec2];
            
            s1 = vecSpecies[(*sg1)[ispec1]]; s2 = vecSpecies[(*sg2)[ispec2]];
            p1 = &(s1->particles);           p2 = &(s2->particles);
            m1 = s1->mass;                   m2 = s2->mass;
            W1 = p1->weight(i1);             W2 = p2->weight(i2);
            
            // Calculate stuff
            m12  = m1 / m2; // mass ratio
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
            COM_gamma = pow( 1.-COM_vsquare , -0.5);
            
            // Change the momentum to the COM frame (we work only on particle 1)
            // Quantities ending with "COM" are quantities of the particle expressed in the COM frame.
            term1 = (COM_gamma - 1.) / COM_vsquare;
            vcv1  = (COM_vx*(p1->momentum(0,i1)) + COM_vy*(p1->momentum(1,i1)) + COM_vz*(p1->momentum(2,i1)))/gamma1;
            vcv2  = (COM_vx*(p2->momentum(0,i2)) + COM_vy*(p2->momentum(1,i2)) + COM_vz*(p2->momentum(2,i2)))/gamma2;
            term2 = (term1*vcv1 - COM_gamma) * gamma1;
            px_COM = (p1->momentum(0,i1)) + term2*COM_vx;
            py_COM = (p1->momentum(1,i1)) + term2*COM_vy;
            pz_COM = (p1->momentum(2,i1)) + term2*COM_vz;
            p2_COM = px_COM*px_COM + py_COM*py_COM + pz_COM*pz_COM;
            p_COM  = sqrt(p2_COM);
            gamma1_COM = (1.-vcv1)*COM_gamma*gamma1;
            gamma2_COM = (1.-vcv2)*COM_gamma*gamma2;
            
            // Calculate some intermediate quantities
            term3 = COM_gamma * gamma12_inv;
            term4 = gamma1_COM * gamma2_COM;
            term5 = term4/p2_COM + m12;
            
            // Calculate coulomb log if necessary
            logL = coulomb_log;
            if( logL <= 0. ) { // if auto-calculation requested
                bmin = max( coeff1/m1/p_COM , abs(coeff2*qqm*term3*term5) ); // min impact parameter
                logL = 0.5*log(1.+debye_length_squared[ibin]/pow(bmin,2));
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
            phi = twoPi * ((double)rand() / RAND_MAX);
            
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
            U = ((double)rand() / RAND_MAX);
            
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
            Ionization->apply(gamma1, gamma2, p1, i1, p2, i2);
            
            if( debug ) {
                smean      ->data_2D[ibin][0] += s;
                logLmean   ->data_2D[ibin][0] += logL;
                //temperature->data_2D[ibin][0] += m1 * (sqrt(1.+pow(p1->momentum(0,i1),2)+pow(p1->momentum(1,i1),2)+pow(p1->momentum(2,i1),2))-1.);
            }
            
        } // end loop on pairs of particles
        
        if( debug && ncol->data_2D[ibin][0]>0.) {
            smean      ->data_2D[ibin][0] /= ncol->data_2D[ibin][0];
            logLmean   ->data_2D[ibin][0] /= ncol->data_2D[ibin][0];
            //temperature->data_2D[ibin][0] /= ncol->data_2D[ibin][0];
        }
        
    } // end loop on bins
    
    // temporary to be removed
    Ionization->finish(vecSpecies[(*sg1)[0]], vecSpecies[(*sg2)[0]], params);
    
    if( debug ) {
        name.str("");
        name << "/t" << setfill('0') << setw(8) << itime << "/s";
        H5::matrix_MPI(did, name.str(), smean      ->data_2D[0][0], (int)totbins, 1, (int)start, (int)nbins);
        name.str("");
        name << "/t" << setfill('0') << setw(8) << itime << "/coulomb_log";
        H5::matrix_MPI(did, name.str(), logLmean   ->data_2D[0][0], (int)totbins, 1, (int)start, (int)nbins);
        //name.str("");
        //name << "/t" << setfill('0') << setw(8) << itime << "/temperature";
        //H5::matrix_MPI(did, name.str(), temperature->data_2D[0][0], (int)totbins, 1, (int)start, (int)nbins);
        if(debye_length_squared.size()>0) {
            // We reuse the smean array for the debye length
            for(unsigned int i=0; i<nbins; i++)
                smean->data_2D[i][0] = sqrt(debye_length_squared[i]) * params.wavelength_SI/(2.*M_PI);
            name.str("");
            name << "/t" << setfill('0') << setw(8) << itime << "/debyelength";
            H5::matrix_MPI(did, name.str(), smean->data_2D[0][0], (int)totbins, 1, (int)start, (int)nbins);
        }
        // Close the group
        H5Gclose(did);
    }

}


// Technique given by Nanbu in http://dx.doi.org/10.1103/PhysRevE.55.4642
//   to pick randomly the deflection angle cosine, in the center-of-mass frame.
// It involves the "s" parameter (~ collision frequency * deflection expectation)
//   and a random number "U".
// Technique slightly modified in http://dx.doi.org/10.1063/1.4742167
inline double Collisions::cos_chi(double s)
{
    
    double A, invA;
    //!\todo make a faster rand by preallocating ??
    double U = (double)rand() / RAND_MAX;
    
    if( s < 0.1 ) {
        if ( U<0.0001 ) U=0.0001; // ensures cos_chi > 0
        return 1. + s*log(U);
    }
    if( s < 3.  ) {
        // the polynomial has been modified from the article in order to have a better form
        invA = 0.00569578 +(0.95602 + (-0.508139 + (0.479139 + ( -0.12789 + 0.0238957*s )*s )*s )*s )*s;
        A = 1./invA;
        return  invA  * log( exp(-A) + 2.*U*sinh(A) );
    }
    if( s < 6.  ) {
        A = 3.*exp(-s);
        return (1./A) * log( exp(-A) + 2.*U*sinh(A) );
    }
    return 2.*U - 1.;
    
}




