
#ifdef SMILEI_ACCELERATOR_GPU_OACC
#include <openacc.h>
#endif

#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

namespace userFunctions
{
// removed all static properties since we are using a namespace

    //! inverse error function is taken from M.B. Giles. 'Approximating the erfinv function'. In GPU Computing Gems, volume 2, Morgan Kaufmann, 2011.
    inline double erfinv_sp( double x )
    {
        double w, p;
        w = -log( ( 1.0 - x ) * ( 1.0 + x ) );
        
        if( w < 5.000000 ) {
            w = w - 2.500000;
            p = +2.81022636000e-08      ;
            p = +3.43273939000e-07 + p*w;
            p = -3.52338770000e-06 + p*w;
            p = -4.39150654000e-06 + p*w;
            p = +0.00021858087e+00 + p*w;
            p = -0.00125372503e+00 + p*w;
            p = -0.00417768164e+00 + p*w;
            p = +0.24664072700e+00 + p*w;
            p = +1.50140941000e+00 + p*w;
        } else {
            w = sqrt( w ) - 3.000000;
            p = -0.000200214257      ;
            p = +0.000100950558 + p*w;
            p = +0.001349343220 + p*w;
            p = -0.003673428440 + p*w;
            p = +0.005739507730 + p*w;
            p = -0.007622461300 + p*w;
            p = +0.009438870470 + p*w;
            p = +1.001674060000 + p*w;
            p = +2.832976820000 + p*w;
        }
        return p*x;
    }
    /**
    * copied from erfinv_DP_1.cu by Prof. Mike Giles.
    * https://people.maths.ox.ac.uk/gilesm/
    * https://people.maths.ox.ac.uk/gilesm/codes/erfinv/
    *
    * Original code is written for CUDA.
    * Mutsuo Saito modified original code for C++.
    */
    inline double erfinv_dp(double x)
    {
        double w, p;
        double sign;
        if (x > 0) {
            sign = 1.0;
        } else {
            sign = -1.0;
            x = abs(x);
        }
        w = - log( (1.0 - x) * (1.0 + x) );

        if ( w < 6.250000 ) {
            w = w - 3.125000;
            p =  -3.6444120640178196996e-21;
            p =   -1.685059138182016589e-19 + p*w;
            p =   1.2858480715256400167e-18 + p*w;
            p =    1.115787767802518096e-17 + p*w;
            p =   -1.333171662854620906e-16 + p*w;
            p =   2.0972767875968561637e-17 + p*w;
            p =   6.6376381343583238325e-15 + p*w;
            p =  -4.0545662729752068639e-14 + p*w;
            p =  -8.1519341976054721522e-14 + p*w;
            p =   2.6335093153082322977e-12 + p*w;
            p =  -1.2975133253453532498e-11 + p*w;
            p =  -5.4154120542946279317e-11 + p*w;
            p =    1.051212273321532285e-09 + p*w;
            p =  -4.1126339803469836976e-09 + p*w;
            p =  -2.9070369957882005086e-08 + p*w;
            p =   4.2347877827932403518e-07 + p*w;
            p =  -1.3654692000834678645e-06 + p*w;
            p =  -1.3882523362786468719e-05 + p*w;
            p =    0.0001867342080340571352 + p*w;
            p =  -0.00074070253416626697512 + p*w;
            p =   -0.0060336708714301490533 + p*w;
            p =      0.24015818242558961693 + p*w;
            p =       1.6536545626831027356 + p*w;
        }
        else if ( w < 16.000000 ) {
            w = sqrt(w) - 3.250000;
            p =   2.2137376921775787049e-09;
            p =   9.0756561938885390979e-08 + p*w;
            p =  -2.7517406297064545428e-07 + p*w;
            p =   1.8239629214389227755e-08 + p*w;
            p =   1.5027403968909827627e-06 + p*w;
            p =   -4.013867526981545969e-06 + p*w;
            p =   2.9234449089955446044e-06 + p*w;
            p =   1.2475304481671778723e-05 + p*w;
            p =  -4.7318229009055733981e-05 + p*w;
            p =   6.8284851459573175448e-05 + p*w;
            p =   2.4031110387097893999e-05 + p*w;
            p =   -0.0003550375203628474796 + p*w;
            p =   0.00095328937973738049703 + p*w;
            p =   -0.0016882755560235047313 + p*w;
            p =    0.0024914420961078508066 + p*w;
            p =   -0.0037512085075692412107 + p*w;
            p =     0.005370914553590063617 + p*w;
            p =       1.0052589676941592334 + p*w;
            p =       3.0838856104922207635 + p*w;
        }
        else {
            w = sqrt(w) - 5.000000;
            p =  -2.7109920616438573243e-11;
            p =  -2.5556418169965252055e-10 + p*w;
            p =   1.5076572693500548083e-09 + p*w;
            p =  -3.7894654401267369937e-09 + p*w;
            p =   7.6157012080783393804e-09 + p*w;
            p =  -1.4960026627149240478e-08 + p*w;
            p =   2.9147953450901080826e-08 + p*w;
            p =  -6.7711997758452339498e-08 + p*w;
            p =   2.2900482228026654717e-07 + p*w;
            p =  -9.9298272942317002539e-07 + p*w;
            p =   4.5260625972231537039e-06 + p*w;
            p =  -1.9681778105531670567e-05 + p*w;
            p =   7.5995277030017761139e-05 + p*w;
            p =  -0.00021503011930044477347 + p*w;
            p =  -0.00013871931833623122026 + p*w;
            p =       1.0103004648645343977 + p*w;
            p =       4.8499064014085844221 + p*w;
        }
        return sign * p * x;
    }
                                     
    //! Load repartition in 1d between MPI processes    
    // ----------------------------------------------------------------------------
    //! \brief Distribute equally the load into chunk of an array
    //! and return the number of elements for the specified chunk number.
    //
    //! \param chunk number
    //! \param nb_chunks Total number of MPI tasks
    //! \param nb_elems Total number of element to be distributed
    //! \param imin Index of the first element for chunk
    //! \param nb_loc_elems Number of element for chunk
    // ----------------------------------------------------------------------------
    inline void distributeArray( int chunk,
                                            int nb_chunks,
                                            int nb_elems,
                                            int &imin,
                                            int &nb_loc_elems )
    {
        // If more ranks than elements,
        // only a part of the processes will work
        if( nb_chunks >= nb_elems ) {
            if( chunk < nb_elems ) {
                imin = chunk;
                nb_loc_elems = 1;
            } else {
                imin = nb_elems;
                nb_loc_elems = 0;
            }
        } else {
        
            int quotient;
            int remainder;
            
            // Part of the load equally distributed
            quotient = nb_elems/nb_chunks;
            
            // Remaining load to be distributed after balanced repartition
            remainder = nb_elems%nb_chunks;
            
            if( chunk < remainder ) {
                imin =  chunk*quotient+chunk;
                nb_loc_elems = quotient + 1;
            } else {
                imin = remainder + chunk*quotient;
                nb_loc_elems = quotient;
            }
        }
    }
                                    
    //! Load repartition in 1d between MPI processes.
    // ----------------------------------------------------------------------------
    //! \brief Distribute equally 1D array into chunks
    //! This function returns tables of indexes and length for each chunk
    //
    //! \param nb_chunks Total number of chunks
    //! \param nb_elems Total number of element to be distributed
    //! \param imin_table Index of the first element for each chunk
    //! \param length_table Number of element for each chunk
    // ----------------------------------------------------------------------------
    inline void distributeArray(
        int nb_chunks,
        int nb_elems,
        int *imin_table,
        int *length_table )
    {

        // If more chunks than elements,
        // only a part of the processes will work
        if( nb_chunks >= nb_elems ) {
            #pragma omp simd
            for( int chunk = 0 ; chunk < nb_elems ; chunk ++ ) {
                imin_table[chunk] = chunk;
                length_table[chunk] = 1;
            }
            #pragma omp simd
            for( int chunk = nb_elems ; chunk < nb_chunks ; chunk ++ ) {
                imin_table[chunk] = nb_elems;
                length_table[chunk] = 0;
            }
        } else {
        
            int quotient;
            int remainder;
            
            // Part of the load equally distributed
            quotient = nb_elems/nb_chunks;
            
            // Remaining load to be distributed after balanced repartition
            remainder = nb_elems%nb_chunks;
            
            #pragma omp simd
            for( int chunk = 0 ; chunk < remainder ; chunk ++ ) {
                imin_table[chunk] =  chunk*quotient+chunk;
                length_table[chunk] = quotient + 1;
            }
            #pragma omp simd
            for( int chunk = remainder ; chunk < nb_chunks ; chunk ++ ) {
                imin_table[chunk] = remainder + chunk*quotient;
                length_table[chunk] = quotient;
            }
        }
    }
        
    //! \brief This function uses a bijection algorithm in a monotonic double array
    //! to find the corresponding index i so that elem is between array[i]
    //! and array[i+1].
    //
    //! \param array array in which to find the value
    //! \param elem element to be found
    //! \param nb_elem number of elements
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    template <class T>
    inline int searchValuesInMonotonicArray(  T * array,
                                     T elem,
                                     int nb_elems )
    {
        int imin = 0; // lower bound
        int imax = nb_elems-1; // upper bound
        int imid = 0;
        
        if( elem == array[0] ) {
            return 0;
        } else if( elem == array[nb_elems-1] ) {
            return nb_elems-2;
        } else {
            while( imax - imin > 1 ) {
                imid= ( imin + imax )/2;
                //imid= (imin + imax)>>1;
                if( elem >= array[imid] ) {
                    imin = imid;
                } else {
                    imax = imid;
                }
            }
            return imin;
        }
    }
}
#endif
