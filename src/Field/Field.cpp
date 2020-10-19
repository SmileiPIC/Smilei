#include "Field.h"

#include <iostream>
#ifdef __PGI
#include <openacc.h>
#endif

void Field::put_to( double val )
{
    if( data_ ) {
        double* hostptr = data_;
#ifdef __PGI
        void* dataptr = acc_deviceptr( hostptr );
        // Test if data exists on GPU, put_to can be used on CPU and GPU during a simulation
        #pragma acc  parallel present(hostptr[0:globalDims_]) if (dataptr!=NULL)
        #pragma acc loop gang worker vector
#endif
        for( unsigned int i=0; i<globalDims_; i++ ) {
            hostptr[i] = val;
        }
    }
}
