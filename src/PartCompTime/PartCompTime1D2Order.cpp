#include "PartCompTime1D2Order.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

PartCompTime1D2Order::PartCompTime1D2Order( ) : PartCompTime() {};

// -----------------------------------------------------------------------------
//! Evaluate the time (simple precision) to compute all particles
//! in the current patch with vectorized operators
//! @param count the numer of particles per cell
//! @param vecto_time time in vector mode
//! @param scalar_time time in scalar mode
// -----------------------------------------------------------------------------
void PartCompTime1D2Order::operator()(  const std::vector<int> &/*count*/,
                                float &vecto_time,
                                float &scalar_time  )
{
    scalar_time = 0;
    vecto_time = 0;
};
