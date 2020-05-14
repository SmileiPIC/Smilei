#ifndef RANDOM_H
#define RANDOM_H

#include <cstdlib>
#include <inttypes.h>
#include <cmath>

class Random
{
public:
    Random( unsigned int seed ) {
        // Initialize the state of the random number generator
        xorshift32_state = seed;
        // Ensure that the random seed is different for each patch
        xorshift32_state += std::rand();
        // zero is not acceptable for xorshift
        if( xorshift32_state==0 ) {
            xorshift32_state = 1073741824;
        }
    };
    
    ~Random() {};
    
    //! random integer
    inline uint32_t integer() {
        return xorshift32();
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
