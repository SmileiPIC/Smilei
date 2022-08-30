// ----------------------------------------------------------------------------
//! \file Table2D.h
//
//! \brief Class Table 2D used to store and manage 2D database for physical processes
//
// ----------------------------------------------------------------------------

#ifndef TABLE2D_H
#define TABLE2D_H

#include <cstring>
#include <string>
#include <vector>

#include "Tools.h"
#include "SmileiMPI.h"
#include "Table.h"

//------------------------------------------------------------------------------
//! class Table2D: manage 2D tabulated data
//------------------------------------------------------------------------------
class Table2D final : public Table 
{
    
public:
    
    // --------------------------------------------------------
    // Main functions
    
    //! Constructor for Table2D
    Table2D();

    //! Destructor for Table2D
    ~Table2D() override final;
    
    // --------------------------------------------------------
    // Parameters
    
    void set(std::vector<double> & input_axis1_min, std::vector<double> & input_data) override final;
    
protected:
    
private:
    
};

#endif 