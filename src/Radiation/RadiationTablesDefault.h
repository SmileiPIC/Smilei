
#ifndef RADIATIONTABLESDEFAULT_H
#define RADIATIONTABLESDEFAULT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>
#include "userFunctions.h"
#include "Params.h"
#include "RadiationTables.h"
#include "H5.h"
#include "Random.h"

// ---------------------------------------------
// Default values (initialization)
// ---------------------------------------------    
class RadiationTablesDefault
{
public :
    static void setDefault( RadiationTables::Niel& niel, RadiationTables::IntegrationFoverChi& integfochi, RadiationTables::Xi& xi );
};

#endif
