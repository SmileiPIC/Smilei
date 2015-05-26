#include "CollisionsFactory.h"
#include "DiagParams.h"

#include "Tools.h"

using namespace std;


// Reads the input file and creates the Collisions objects accordingly
vector<Collisions*> CollisionsFactory::create(PicParams& params, InputData &ifile, vector<Species*>& vecSpecies)
{
    std::vector<Collisions*> vecCollisions;

    vector<string> sg1, sg2;
    vector<unsigned int> sgroup1, sgroup2;
    double clog;
    bool intra, debye_length_required = false;
    ostringstream mystream;
    
    // Loop over each binary collisions group and parse info
    int n_collisions = 0;
    while (ifile.existComponent("collisions",n_collisions)) {
        
        MESSAGE("Parameters for collisions #" << n_collisions << " :");
        
        // Read the input file by searching for the keywords "species1" and "species2"
        // which are the names of the two species that will collide
        sg1.resize(0);
        sg2.resize(0);
        ifile.extract("species1",sg1,"collisions",n_collisions);
        ifile.extract("species2",sg2,"collisions",n_collisions);
        
        // Obtain the lists of species numbers from the lists of species names.
        sgroup1 = DiagParams::FindSpecies(sg1,params);
        sgroup2 = DiagParams::FindSpecies(sg2,params);
        
        // Each group of species sgroup1 and sgroup2 must not be empty
        if (sgroup1.size()==0) ERROR("No valid `species1` requested in collisions #" << n_collisions);
        if (sgroup2.size()==0) ERROR("No valid `species2` requested in collisions #" << n_collisions);
        
        // sgroup1 and sgroup2 can be equal, but cannot have common species if they are not equal
        if (sgroup1 != sgroup2) {
            for (unsigned int i1=0; i1<sgroup1.size(); i1++) {
                for (unsigned int i2=0; i2<sgroup2.size(); i2++) {
                    if (sgroup1[i1] == sgroup2[i2])
                        ERROR("Unauthorized species (#" << sgroup1[i1]
                              << ") in collisions #" << n_collisions
                              << " (inter-collisions must not have a species colliding with itself)");
                }
            }
            intra = false;
        } else {
            intra = true;
        }
        
        // Coulomb logarithm (if negative or unset, then automatically computed)
        clog = 0.; // default
        ifile.extract("coulomb_log",clog,"collisions",n_collisions);
        if (clog <= 0.) debye_length_required = true; // auto coulomb log requires debye length
        
        // Print collisions parameters
        mystream.str(""); // clear
        for (unsigned int rs=0 ; rs<sgroup1.size() ; rs++) mystream << " #" << sgroup1[rs];
        MESSAGE(1,"First  group of species :" << mystream.str());
        mystream.str(""); // clear
        for (unsigned int rs=0 ; rs<sgroup2.size() ; rs++) mystream << " #" << sgroup2[rs];
        MESSAGE(1,"Second group of species :" << mystream.str());
        MESSAGE(1,"Coulomb logarithm       : " << clog);
        MESSAGE(1,"Intra collisions        : " << (intra?"True":"False"));
        
        // Add new Collisions objects to vector
        vecCollisions.push_back( new Collisions(params,n_collisions,sgroup1,sgroup2,clog,intra) );
        n_collisions++;
    }
    
    // Needs wavelength_SI to be defined
    if (n_collisions > 0)
        if (params.wavelength_SI <= 0.)
            ERROR("The parameter `wavelength_SI` needs to be defined and positive in order to compute collisions");
    
    // pass the variable "debye_length_required" into the Collision class
    Collisions::debye_length_required = debye_length_required;
    
    return vecCollisions;
}



