#ifndef TABULATEDFUNCTIONS_H
#define TABULATEDFUNCTIONS_H

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>

//! class of tabulated functions
class tabulatedFunctions {
    
public:
    
    //! constructor of tabulatedFunctions
    tabulatedFunctions();
    
    //! destructor of tabulatedFunctions
    ~tabulatedFunctions();
    
    //! returns inverse error function of a double x
    double erfinv(double x);
    
    //! needs to be called one time before using erfinv
    void erfinv_loadTab();
    
private:
    
    //! boolean .true. as long as the table as never been loaded
    bool erfinv_needLoad_;
    
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
    
    
};//tabulatedFunctions

#endif




