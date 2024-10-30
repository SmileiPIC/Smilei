#ifndef RANDOMSHUFFLE_H
#define RANDOMSHUFFLE_H

#include "Random.h"

static constexpr uint32_t best_prime[] = {2, 2, 5, 11, 17, 37, 67, 131, 257, 521, 1031, 2053, 4099, 8209, 16411, 32771, 65537, 131101, 262147, 524309, 1048583, 2097169, 4194319, 8388617};

//! Maximum number of elements for method 0
static constexpr uint32_t shuffle_threshold = 8;

// A random shuffler which ensures a small memory occupancy (fixed)
// at the cost of a reasonable loss in performance
class RandomShuffle
{
public:

    RandomShuffle( Random &rand, uint32_t length )
    : length_( length ), mask_( 1 ), P_( 0 ), i_( 0 )
    {
        reinit( rand, length );
    }
    
    //! Reinitialize the shuffler
    SMILEI_ACCELERATOR_DECLARE_ROUTINE
    void reinit( Random &rand, uint32_t length ) {
        i_ = 0;
        length_ = length;
        
        if( length_ < shuffle_threshold ) {
            
            // Fischer-Yates method
            
            for( uint32_t i = 0; i < length; i++ ){
                random_array_[i] = i;
            }
            for( uint32_t i = length-1; i > 0; i-- ){
                uint32_t j = rand.integer() % (i+1);
                uint32_t a = random_array_[i], b = random_array_[j];
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
            for( uint32_t i=0; i<7; i++ ){
                uint32_t r = rand.integer();
                // Prevent key = zero
                do {
                    random_array_[i] = r & mask_;
                    r >>= 1;
                } while( random_array_[i] == 0 );
            }
        }
    }
    SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
    
    //! Get the next shuffled position
    uint32_t next() {
        if( length_ < shuffle_threshold ){
            return next0();
        } else {
            return next1();
        }
    }
    
    //! Get n next shuffled positions
    void next( uint32_t n, uint32_t * shuffled ) {
        if( length_ < shuffle_threshold ){
            SMILEI_ACCELERATOR_LOOP_SEQ
            for( uint32_t i = 0; i < n; i++ ) {
                shuffled[i] = next0();
            }
        } else {
            SMILEI_ACCELERATOR_LOOP_SEQ
            for( uint32_t i = 0; i < n; i++ ) {
                shuffled[i] = next1();
            }
        }
    }
    
    //! Get the next shuffled position for method 0
    uint32_t next0() {
        uint32_t next_position = random_array_[i_];
        i_ = ( i_ + 1 ) % length_;
        return next_position;
    }
    
    //! Get the next shuffled position for method 1
    uint32_t next1() {
        uint32_t next_position;
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
        return next_position;
    }
    
private:
    
    //! number of elements to permute
    uint32_t length_;
    
    //! A random array required initially for the algorithm
    uint32_t random_array_[shuffle_threshold];
    
    //! Bitmask required for method 1
    uint32_t mask_;
    
    //! Prime number required for method 1
    uint32_t P_;
    
    //! A number related to the number of calls to next()
    uint32_t i_;
    
};


#endif
