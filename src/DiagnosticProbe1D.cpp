#include "DiagnosticProbe1D.h"

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

DiagnosticProbe1D::DiagnosticProbe1D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) : DiagnosticProbe(params,diagParams,smpi,1) {
    every.resize(diagParams->probe1DStruc.size());
    probeParticles.resize(diagParams->probe1DStruc.size());
    probeId.resize(diagParams->probe1DStruc.size());

    for (unsigned int np=0; np<diagParams->probe1DStruc.size(); np++) {
        
        every[np]=diagParams->probe1DStruc[np].every;
        unsigned int ndim=params->nDim_particle;

        vector<unsigned int> vecNumber=diagParams->probe1DStruc[np].number;
        
        probeParticles[np].initialize(vecNumber[0], ndim);
        probeId[np].resize(vecNumber[0]);

        vector<double> partPos(ndim*vecNumber[0]);
        for(unsigned int count=0; count!=vecNumber[0]; ++count) {
            int found=smpi->getRank();
            for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
                unsigned int k=iDim+count*ndim;
                partPos[k]=diagParams->probe1DStruc[np].pos[iDim]+count*(diagParams->probe1DStruc[np].posFirst[iDim]-diagParams->probe1DStruc[np].pos[iDim])/(vecNumber[0]-1);
                
                probeParticles[np].position(iDim,count) = 2*M_PI*partPos[k];
                
                if(smpi->getDomainLocalMin(iDim) >  probeParticles[np].position(iDim,count) || smpi->getDomainLocalMax(iDim) <= probeParticles[np].position(iDim,count)) {
                    found=-1;
                }
            }
            probeId[np][count] = found;
        }
        addProbe(np, partPos, vecNumber);
    }
}

