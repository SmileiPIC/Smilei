// ----------------------------------------------------------------------------
//! \file Merging.h
//
//! \brief Header for the generic class Merging
//! dedicated to the particle merging.
//
//! Mathieu Lobet
//! 01/2019
//
// ----------------------------------------------------------------------------

#ifndef MERGING_H
#define MERGING_H

//  ----------------------------------------------------------------------------
//! Class Merging
//  ----------------------------------------------------------------------------
class Merging
{
public:

    //! Creator for Merging
    Merging(Params& params, Species *species);
    virtual ~Merging();

    //! Overloading of () operator
    //! \param particles   particle object containing the particle
    //!                    properties of the current species
    //! \param smpi        MPI properties
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    virtual void operator() (
            Particles &particles,
            SmileiMPI* smpi,
            int istart,
            int iend,
            int ithread) = 0;


protected:
private:
};

#endif
