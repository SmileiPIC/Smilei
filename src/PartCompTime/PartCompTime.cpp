#include "PartCompTime.h"

PartCompTime::PartCompTime( )
{
};

// -----------------------------------------------------------------------------
//! Evaluate the time (simple precision) to compute all particles
//! in the current patch with vectorized operators
//! @param count the numer of particles per cell
//! @aram vecto_time time in vector mode
//! @aram scalar_time time in scalar mode
// -----------------------------------------------------------------------------
void PartCompTime::operator()(  const std::vector<int> &count,
                                float &vecto_time,
                                float &scalar_time  )
{};
