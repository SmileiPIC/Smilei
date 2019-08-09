#ifndef PARTICLEINJECTOR_H
#define PARTICLEINJECTOR_H

#include <vector>
#include <string>

#include "Params.h"

//! class ParticleInjector
class ParticleInjector
{
public:
    
    // -----------------------------------------------------------------------------
    //  1. Constructor/destructor

    //! ParticleInjector constructor
    ParticleInjector( Params &, Patch * );

    //! ParticleInjector destructor
    virtual ~ParticleInjector();

    
    // -----------------------------------------------------------------------------
    //  2. Injector parameters
    
    //! number of the injector
    unsigned int injector_number_;

    //! Name of this injector
    std::string name_;

    //! kind/name of the species associated to this injector
    std::string species_name_;

    //! number of the species associated to this injector
    unsigned int species_number_;

    //! box side to where inject particles
    std::string box_side_;
    
    //! position initialization of the injector
    std::string position_initialization_;

    // -----------------------------------------------------------------------------
    //  3. Methods

    //! Method returning the number of the associated species
    inline unsigned int getSpeciesNumber() const
    {
        return species_number_;
    }

    //! Return if the injector is from Xmin
    inline bool isXmin()
    {
        return (box_side_ == "xmin");
    }

    //! Return if the injector is from Xmin
    inline bool isXmax()
    {
        return (box_side_ == "xmax");
    }

};

#endif
