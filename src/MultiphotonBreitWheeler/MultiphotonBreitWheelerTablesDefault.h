
#ifndef MBWTABLESDEFAULT_H
#define MBWTABLESDEFAULT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include "MultiphotonBreitWheelerTables.h"
#include "Table.h"

// ---------------------------------------------
// Default values (initialization)
// ---------------------------------------------
class MultiphotonBreitWheelerTablesDefault
{
public :
    static void setDefault( Table & T, MultiphotonBreitWheelerTables::Xi& xi );
};

#endif
