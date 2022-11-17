#ifndef PARTCOMPTIME1D2ORDER_H
#define PARTCOMPTIME1D2ORDER_H

#include "PartCompTime.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PartCompTime1D2Order
//! Evaluation of the particle time in 1D at Order 2
//  --------------------------------------------------------------------------------------------------------------------
class PartCompTime1D2Order final : public PartCompTime
{
public:
    PartCompTime1D2Order();
    ~PartCompTime1D2Order() override final {};
    
    // -------------------------------------------------------------------------
    //! Evaluate the time (simple precision) to compute all particles
    //! in the current patch with vectorized operators
    //! @param count the numer of particles per cell
    //! @aram vecto_time time in vector mode
    //! @aram scalar_time time in scalar mode
    // -------------------------------------------------------------------------
    virtual void operator() (   const std::vector<int> &count,
                        float &vecto_time,
                        float &scalar_time ) override final;
    
private:

};//END class PartCompTime1D2Order

#endif
