#ifndef RANDOMSHUFFLE_H
#define RANDOMSHUFFLE_H

#include "Random.h"

#define SMILEI_SHUFFLE_THRESHOLD 8

static constexpr size_t best_prime[] = {2, 2, 5, 11, 17, 37, 67, 131, 257, 521, 1031, 2053, 4099, 8209, 16411, 32771, 65537, 131101, 262147, 524309, 1048583, 2097169, 4194319, 8388617};

// A random shuffler which ensures a small memory occupancy (fixed)
// at the cost of a reasonable loss in performance
class RandomShuffle
{
public:

    RandomShuffle( Random &rand, size_t length )
    : length_( length ), mask_( 1 ), P_( 0 ), i_( 0 )
    {
        
        if( length_ < SMILEI_SHUFFLE_THRESHOLD ) {
            
            // Fischer-Yates method
            
            for( size_t i = 0; i < length; i++ ){
                random_array_[i] = i;
            }
            for( size_t i = length-1; i > 0; i-- ){
                size_t j = rand.integer() % (i+1);
                size_t a = random_array_[i], b = random_array_[j];
                random_array_[i] = b; random_array_[j] = a;
            }
            
        } else {
            
            // Method from https://crypto.stackexchange.com/questions/107291/very-small-domains-in-ff1-or-similar
            
            // Find number of bits set in length
            unsigned int nbits = 0, N = length;
            while( N > 0 ) {
                nbits++; N >>= 1;
            }
            // Make a bit mask
            mask_ = (1 << nbits) - 1;
            // Find the closest prime above 2^nbits
            P_ = best_prime[nbits];
            // Generate 7 random numbers that constitutes the key
            for( size_t i=0; i<7; i++ ){
                size_t r = rand.integer();
                // Prevent key = zero
                do {
                    random_array_[i] = r & mask_;
                    r >>= 1;
                } while( random_array_[i] == 0 );
            }
        }
        
    }
    
    //! Get the next shuffled position
    size_t next() {
        
        size_t next_position;
        
        if( length_ < SMILEI_SHUFFLE_THRESHOLD ){
            
            next_position = random_array_[i_];
            i_ = ( i_ + 1 ) % length_;
            
        } else {
            
            do{
                next_position = i_;
                i_ = ( i_ + 1 ) % mask_;
                next_position ^= random_array_[0];
                do {
                    next_position = (random_array_[3] * next_position + random_array_[4]) % P_;
                } while( next_position > mask_ );
                next_position ^= random_array_[1];
                do {
                    next_position = (random_array_[5] * next_position + random_array_[6]) % P_;
                } while( next_position > mask_ );
                next_position ^= random_array_[2];
            } while( next_position >= length_ );
            
        }
        
        return next_position;
    }
    
private:
    
    //! number of elements to permute
    size_t length_;
    
    //! A random array required initially for the algorithm
    size_t random_array_[SMILEI_SHUFFLE_THRESHOLD];
    
    //! Bitmask required for method 1
    size_t mask_;
    
    //! Prime number required for method 1
    size_t P_;
    
    //! A number related to the number of calls to next()
    size_t i_;
    
};


#endif
