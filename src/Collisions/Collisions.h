#ifndef COLLISIONS_H
#define COLLISIONS_H

#include <vector>
#include <cmath>

#include "Tools.h"
#include "H5.h"
#include "Random.h"
#include "Params.h"
#include "BinaryProcessData.h"

class Patch;
class Species;
class VectorPatch;

class Collisions
{
public:
    //! Constructor for Collisions between two species
    Collisions(
        Params &params,
        double coulomb_log,
        double coulomb_log_factor
    );
    //! Cloning Constructor
    Collisions( Collisions * );
    //! destructor
    ~Collisions();
    
    void prepare();
    #pragma acc routine vector nohost
    void apply( Random *random, BinaryProcessData &D );
    void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool intra, std::vector<unsigned int> sg1, std::vector<unsigned int> sg2, int itime );
    std::string name() {
        std::ostringstream t;
        t << "Collisions with Coulomb logarithm: ";
        if( coulomb_log_ == 0. ) {
            t << "auto";
        } else {
            t << coulomb_log_;
        }
        if( coulomb_log_factor_ != 1. ) {
            t << " x " << coulomb_log_factor_;
        }
        return t.str();
    };
    
    unsigned int npairs_tot_;
    double smean_, logLmean_;
    
protected:
    
    //! Coulomb logarithm (zero or negative means automatic)
    const double coulomb_log_;
    
    //! Coulomb logarithm factor
    const double coulomb_log_factor_;
    
    const double twoPi = 2. * 3.14159265358979323846;
    const double coeff1_, coeff2_, coeff3_, coeff4_;
    
};


#endif
