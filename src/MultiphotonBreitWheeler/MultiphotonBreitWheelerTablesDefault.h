
#ifndef MBWTABLESDEFAULT_H
#define MBWTABLESDEFAULT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include "MultiphotonBreitWheelerTables.h"
#include "Table.h"
#include "Table2D.h"

// ---------------------------------------------
// Default values (initialization)
// ---------------------------------------------
class MultiphotonBreitWheelerTablesDefault
{
public :
    static void setDefault( Table & T, Table & xi );
};

#endif
