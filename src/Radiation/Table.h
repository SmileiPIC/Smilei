// ----------------------------------------------------------------------------
//! \file Table.h
//
//! \brief Class Table definition (for QED)
//
// ----------------------------------------------------------------------------

#ifndef TABLE_H
#define TABLE_H

#include <cstring>
#include <string>
#include <vector>

#include "Tools.h"
#include "SmileiMPI.h"

//------------------------------------------------------------------------------
//! Table class: manage tabulated data
//------------------------------------------------------------------------------
class Table
{

public:

    //! Constructor for Table
    Table();

    //! Destructor for Table
    ~Table();

    // --------------------------------------------------------
    // Main functions
    
    //! Allocation of data_
    void allocate();
    
    //! Bcast data_ and metadata to all MPI processes
    void bcast(SmileiMPI *smpi);
    
    //! get value using linear interpolation
    double get(double x);

    //! Copy values from input_data to the table data
    double set(std::vector<double> input_data);

    // --------------------------------------------------------
    // Parameters

    // Main data array
    double * data_ = nullptr;
    
    // Size of the table
    unsigned int size_;
    
    //! Minimum values
    double min_;
    
    //! Log10 of the minimum values
    double log10_min_;
    
    //! Max values
    double max_;
    
    //! Delta between values (when constant)
    double delta_;
    
    //! Inverse delta
    double inv_delta_;
    
    //! table name
    std::string name_;
    
    //! axis names
    std::string axis_name_;
    
protected:
    
private:

};


#endif 