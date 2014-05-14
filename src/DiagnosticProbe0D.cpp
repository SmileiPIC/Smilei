#include "DiagnosticProbe0D.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Field.h"

using namespace std;

DiagnosticProbe0D::DiagnosticProbe0D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) : DiagnosticProbe(params,diagParams,smpi,0) {
    every.resize(diagParams->probe0DStruc.size());
    probeParticles.resize(diagParams->probe0DStruc.size());
    probeId.resize(diagParams->probe0DStruc.size());
        
    for (unsigned int np=0; np<diagParams->probe0DStruc.size(); np++) {
        
        every[np]=diagParams->probe0DStruc[np].every;

        unsigned int ndim=params->nDim_particle;
 
        vector<unsigned int> vecNumber(1);
        vecNumber[0]=1; // just one probe per group

        probeParticles[np].initialize(vecNumber[0], ndim);
        probeId[np].resize(vecNumber[0]);
        
        vector<double> partPos(ndim);
        
        int found=smpi->getRank();
        for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
            partPos[iDim]=diagParams->probe0DStruc[np].pos[iDim];
            probeParticles[np].position(iDim,0) = 2*M_PI*partPos[iDim];
            if(smpi->getDomainLocalMin(iDim) >  probeParticles[np].position(iDim,0) || smpi->getDomainLocalMax(iDim) <= probeParticles[np].position(iDim,0)) {
                found=-1;
            }
        }
        probeId[np][0] = found;
        
        addProbe(np, partPos, vecNumber);
    }
}

