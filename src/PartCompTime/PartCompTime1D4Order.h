#ifndef PARTCOMPTIME1D4ORDER_H
#define PARTCOMPTIME1D4ORDER_H

#include "PartCompTime.h"

//! Placeholder computation costs for fourth-order 1D operators.
//! Scalar and vectorized paths use the same operators in this geometry.
class PartCompTime1D4Order final : public PartCompTime
{
public:
    PartCompTime1D4Order();
    ~PartCompTime1D4Order() override final {};

    void operator()( const std::vector<int> &count,
                     float &vecto_time,
                     float &scalar_time ) override final;
};

#endif
