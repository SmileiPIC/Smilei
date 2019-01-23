#ifndef TABULATEDFUNCTIONS_H
#define TABULATEDFUNCTIONS_H

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Tools.h"



//! singleton class of tabulated functions

class erfinv
{
public:
    static erfinv &instance()
    {
        static erfinv one_and_only_instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return one_and_only_instance;
    }
    
    //! returns inverse error function of a double x
    double call( double x );
    
    //! needs to be called one time before using erfinv
    void prepare();
    
protected:
    // creator is private for singletons
    erfinv() {};
    erfinv( erfinv const & ); // avoid implementation of this
    void operator=( erfinv const & ); // avoid implementation of this
    
private:

    //! number of points used to sample the fct
    unsigned int erfinv_tabSize_;
    
    //! mininum value for x
    double erfinv_xmin_;
    
    //! maximum value for x
    double erfinv_xmax_;
    
    //! constant used to compute the different values of x used during the sampling
    double erfinv_alpha_;
    
    //! vector storing the values of x used to sample erfinv
    std::vector<double> erfinv_x_;
    
    //! vector storing the sampled values of erfinv
    std::vector<double> erfinv_tab_;
    
};


#endif




