#ifndef RANDOM_H
#define RANDOM_H

#include <cstdlib>
#include <inttypes.h>
#include <cmath>
#include "userFunctions.h"

namespace Random_namespace // in order to use the random functions without having access to the class random
{

    inline uint32_t xorshift32(uint32_t xorshift32_state)
    {
        // Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" 
        xorshift32_state ^= xorshift32_state << 13;
        xorshift32_state ^= xorshift32_state >> 17;
        xorshift32_state ^= xorshift32_state << 5;
        return xorshift32_state;
    }

    inline double uniform1(uint32_t xorshift32_state) {
        constexpr double xorshift32_invmax1 = (1.-1e-11)/4294967296.;
        return xorshift32(xorshift32_state) * xorshift32_invmax1;
    }

    inline double perp_rand_dp(uint32_t xorshift32_state) {
        double a = userFunctions::erfinv_dp( uniform1(xorshift32_state) ); 
        // technically we could also use the erfinv() function fron cuda, it would require compiling with -cuda though ...
        // the study showed the gap in perf for BC thermal was not worth the added depend
        if( xorshift32(xorshift32_state) & 1  ) { // rand->cointoss()
            a *= -1.;
        }
        return a;
    }
}

class Random
{
public:
    Random( unsigned int seed ) {
        // Initialize the state of the random number generator
        xorshift32_state = seed;
        // zero is not acceptable for xorshift
        if( xorshift32_state==0 ) {
            xorshift32_state = 1073741824;
        }
    }

    //! random integer
    inline uint32_t integer() {
        return xorshift32();
    }
    //! Random true/false
    inline bool cointoss() {
        return xorshift32() & 1;
    }
    //! Uniform rand from xorshift32 generator, between 0 (excluded) and 1 (included)
    inline double uniform() {
        return xorshift32() * xorshift32_invmax;
    }
    //! Uniform rand from xorshift32 generator, between 0 (excluded) and 1-10^-11
    inline double uniform1() {
        return xorshift32() * xorshift32_invmax1;
    }
    //! Uniform rand from xorshift32 generator, between -1. (excluded) and 1. (included)
    inline double uniform2() {
        return xorshift32() * xorshift32_invmax2 - 1.;
    }
    //! Uniform rand from xorshift32 generator, between 0. (excluded) and 2 pi (included)
    inline double uniform_2pi() {
        return xorshift32() * xorshift32_invmax_2pi;
    }
    //! Normal rand from xorshift32 generator (std deviation = 1.)
    inline double normal() {
        static double spare;
        static bool has_spare = false;
        if( has_spare ) {
            has_spare = false;
            return spare;
        } else {
            double u, v, s;
            do {
                u = uniform2();
                v = uniform2();
                s = u*u + v*v;
            } while( s >= 1. );
            s = std::sqrt( -2. * std::log(s) / s );
            spare = v * s;
            has_spare = true;
            return u * s;
        }
    }

    //! State of the random number generator
    uint32_t xorshift32_state;

private:
    
    //! Random number generator
    inline uint32_t xorshift32()
    {
        /* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
        xorshift32_state ^= xorshift32_state << 13;
        xorshift32_state ^= xorshift32_state >> 17;
        xorshift32_state ^= xorshift32_state << 5;
        return xorshift32_state;
    }
    //! Inverse of the maximum value of the random number generator
    static constexpr double xorshift32_invmax = 1./4294967296.;
    //! Almost inverse of the maximum value of the random number generator
    static constexpr double xorshift32_invmax1 = (1.-1e-11)/4294967296.;
    //! Twice inverse of the maximum value of the random number generator
    static constexpr double xorshift32_invmax2 = 2./4294967296.;
     //! two pi * inverse of the maximum value of the random number generator
    static constexpr double xorshift32_invmax_2pi = 2.*M_PI/4294967296.;
    
};


#endif
