#ifndef DIAGNOSTICPARTICLEBINNING_H
#define DIAGNOSTICPARTICLEBINNING_H

#include "DiagnosticParticleBinningBase.h"

class DiagnosticParticleBinning : public DiagnosticParticleBinningBase
{
    friend class SmileiMPI;
    
public :

    //! Default constructor
    DiagnosticParticleBinning(
        Params &params,
        SmileiMPI *smpi,
        Patch *patch,
        int diagId
    );
    //! Default destructor
    ~DiagnosticParticleBinning() override;
    
    static std::vector<std::string> excludedAxes() {
        std::vector<std::string> excluded_axes( 0 );
        excluded_axes.push_back( "a" );
        excluded_axes.push_back( "b" );
        excluded_axes.push_back( "theta" );
        excluded_axes.push_back( "phi" );
        return excluded_axes;
    }
    
protected:
    
};

#endif

