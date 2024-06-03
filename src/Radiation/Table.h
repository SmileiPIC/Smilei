// ----------------------------------------------------------------------------
//! \file Table.h
//
//! \brief  Class Table used to store and manage 1D database for physical processes
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
    virtual ~Table();

    // --------------------------------------------------------
    // Main functions
    
    //! Allocation of data_
    void allocate();
    
    //! Bcast data_ and metadata to all MPI processes
    void bcast(SmileiMPI *smpi);
    
    //! Set the size of the data
    void set_size(unsigned int * dim_size);
    
    //! Compute usefull parameters using inputs
    void compute_parameters();
    
    //! get value using linear interpolation at position x
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    #pragma acc routine seq
#endif
    double get(double x);

    //! Copy values from input_data to the table data
    //! params[in] std::vector<double> & input_data : vector to be used to initialize table data
    void set(std::vector<double> & input_data);
    
    //! Return pointer to the data_ array
    inline double * __attribute__((always_inline)) data()
    {
        return &data_[0];
    };
    
    virtual void set(std::vector<double> &, std::vector<double> &) {};

    // --------------------------------------------------------
    // Parameters

    // Main data array
    double * data_ = nullptr;
    
    // Min values for axis 1 (only relevant if dimension_ > 1)
    double * axis1_min_ = nullptr;
    
    // Table dimension
    unsigned int dimension_;
    
    // Size of the table
    unsigned int size_;
    
    // Size for each dimension
    unsigned int dim_size_[2];
    
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
    
    double inv_dim_size_minus_one_[2];
    
    //! table name
    std::string name_;
    
    //! axis names
    std::string axis_name_[2];
    
protected:
    
private:

};


#endif 