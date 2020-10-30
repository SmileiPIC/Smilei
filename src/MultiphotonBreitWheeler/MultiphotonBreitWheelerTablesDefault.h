
#ifndef MBWTABLESDEFAULT_H
#define MBWTABLESDEFAULT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include "MultiphotonBreitWheelerTables.h"

// ---------------------------------------------
// Default values (initialization)
// ---------------------------------------------
class MultiphotonBreitWheelerTablesDefault
{
public :
    static void setDefault( MultiphotonBreitWheelerTables::T& table, MultiphotonBreitWheelerTables::Xi& xi );
};

#endif
