#include "PyTools.h"
#include <iomanip>

#include "DiagnosticParticleBinning.h"

using namespace std;


// Constructor
DiagnosticParticleBinning::DiagnosticParticleBinning(
    Params &params,
    SmileiMPI *smpi,
    Patch *patch,
    int diagId
) : DiagnosticParticleBinningBase( params, smpi, patch, diagId, "ParticleBinning", false, nullptr, excludedAxes() )
{
}

DiagnosticParticleBinning::~DiagnosticParticleBinning()
{
}

