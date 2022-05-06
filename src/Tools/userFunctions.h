
#ifdef _GPU
#include <openacc.h>
#endif

#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

class userFunctions
{

public:

    static double erfinv( double x );
    static double erfinv2( double x );
                                     
    //! Load repartition in 1d between MPI processes
    static void distributeArray( int rank,
                                    int nb_ranks,
                                    int nb_elems,
                                    int &imin,
                                    int &nb_loc_elems );
                                    
    //! Load repartition in 1d between MPI processes.
    //! This function returns tables of indexes and length for all rank
    static void distributeArray(
        int nb_chunks,
        int nb_elems,
        int *imin_table,
        int *length_table );
        
    //! \brief This function uses a bijection algorithm in a monotonic double array
    //! to find the corresponding index i so that elem is between array[i]
    //! and array[i+1].
    //
    //! \param array array in which to find the value
    //! \param elem element to be found
    //! \param nb_elem number of elements
    
#ifdef _GPU
    #pragma acc routine seq
#endif
    template <class T>
    static int searchValuesInMonotonicArray(  T * array,
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
                                 
private:


};
#endif
