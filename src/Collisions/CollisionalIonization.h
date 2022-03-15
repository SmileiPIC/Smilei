
#ifndef COLLISIONALIONIZATION_H
#define COLLISIONALIONIZATION_H

#include <vector>

#include "Tools.h"
#include "Species.h"
#include "Params.h"
#include "BinaryProcess.h"

class Patch;

class CollisionalIonization : public BinaryProcess
{

public:
    //! Constructor
    CollisionalIonization( int, Params*, int, Particles* );
    //! Cloning Constructor
    CollisionalIonization( CollisionalIonization * );
    //! destructor
    ~CollisionalIonization() {};
    
    void prepare() {};
    void apply( Random *random, BinaryProcessData &D );
    void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool intra, std::vector<unsigned int> sg1, std::vector<unsigned int> sg2, int itime );
    std::string name() {
        std:: ostringstream t;
        t << "Collisional ionization with atomic number "<<atomic_number<<" towards species #"<<ionization_electrons_;
        return t.str();
    };
    
    //! Initializes the arrays in the database and returns the index of these arrays in the DB
    virtual unsigned int createDatabase( double );
    //! Assigns the correct databases
    virtual void assignDatabase( unsigned int );
    
    //! Gets the k-th binding energy of any neutral or ionized atom with atomic number Z and charge Zstar
    double binding_energy( int Zstar, int k );
    
    //! Coefficients used for interpolating the energy over a given initial list
    static const double a1, a2, npointsm1;
    static const int npoints;
    
    //! Local table of integrated cross-section
    std::vector<std::vector<double> > *crossSection;
    //! Local table of average secondary electron energy
    std::vector<std::vector<double> > *transferredEnergy;
    //! Local table of average incident electron energy lost
    std::vector<std::vector<double> > *lostEnergy;
    
    //! New electrons temporary species
    Particles new_electrons;
    
    //! Index of the atomic number in the databases
    unsigned int dataBaseIndex;
    
private:
    
    //! Atomic number
    int atomic_number;
    
    //! Species where new electrons are sent
    int ionization_electrons_;
    
    //! Global table of atomic numbers
    static std::vector<int> DB_Z;
    //! Global table of integrated cross-section
    static std::vector<std::vector<std::vector<double> > > DB_crossSection;
    //! Global table of average secondary electron energy
    static std::vector<std::vector<std::vector<double> > > DB_transferredEnergy;
    //! Global table of average incident electron energy lost
    static std::vector<std::vector<std::vector<double> > > DB_lostEnergy;
    
    //! True if first group of species is the electron
    bool electronFirst;
    
    //! Current ionization rate array (one cell per number of ionization events)
    std::vector<double> rate;
    std::vector<double> irate;
    //! Current ionization probability array (one cell per number of ionization events)
    std::vector<double> prob;
    
    //! Method called by ::apply to calculate the ionization, being sure that electrons are the first species
    void calculate( double, double, double, Particles *pe, int ie, Particles *pi, int ii, double U1, double U2, double coeff );
    
};

#endif
