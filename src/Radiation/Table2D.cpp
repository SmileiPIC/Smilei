// ----------------------------------------------------------------------------
//! \file Table2D.cpp
//
//! \brief Class Table 2D used to store and manage 2D database for physical processes
//
// ----------------------------------------------------------------------------

#include "Table2D.h"

// -----------------------------------------------------------------------------
// Constructor for Table
// -----------------------------------------------------------------------------
Table2D::Table2D()
{
    dimension_ = 2;
}

// -----------------------------------------------------------------------------
// Destructor for Table
// -----------------------------------------------------------------------------
Table2D::~Table2D()
{
    if (axis1_min_) {
        delete [] axis1_min_;
        axis1_min_ = nullptr;
    }
    if (data_) {
        delete [] data_;
        data_ = nullptr;
        size_ = 0;
        for (unsigned int i = 0 ; i < dimension_ ; i++) {
            dim_size_[i] = 0;
        }
    }
}

// -----------------------------------------------------------------------------
//! Copy values from input_data to the table data
// -----------------------------------------------------------------------------
void Table2D::set (std::vector<double> & input_axis1_min, std::vector<double> & input_data) {

    if (input_data.size() != size_) {
        ERROR("Impossible to initialize data in Table " << name_ << " because sizes do not match.")
    }

    for (unsigned int i = 0 ; i < size_ ; i++) {
        data_[i] = input_data[i];
    }
    
    for (unsigned int i = 0 ; i < dim_size_[0] ; i++) {
        axis1_min_[i] = input_axis1_min[i];
    }
}