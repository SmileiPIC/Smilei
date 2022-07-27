// ----------------------------------------------------------------------------
//! \file RadiationTablesDefault.h
//
//! \brief This class contains default tables to initalize QED tables
//
// ----------------------------------------------------------------------------

#ifndef RADIATIONTABLESDEFAULT_H
#define RADIATIONTABLESDEFAULT_H

#include <vector>
#include <cmath>
#include "RadiationTables.h"

// ---------------------------------------------
// Default values (initialization)
// ---------------------------------------------    
class RadiationTablesDefault
{
public :
    static void setDefault( Table & niel, Table & integfochi, Table & xi );
};

#endif
