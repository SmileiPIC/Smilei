#include "PartCompTime1D4Order.h"

PartCompTime1D4Order::PartCompTime1D4Order() : PartCompTime() {};

void PartCompTime1D4Order::operator()( const std::vector<int> &/*count*/,
                                    float &vecto_time,
                                    float &scalar_time )
{
    // As in second-order 1D, there is no distinct vectorized operator cost.
    scalar_time = 0;
    vecto_time = 0;
}
