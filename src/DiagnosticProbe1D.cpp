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

    open("Probes1D.h5");

    for (unsigned int np=0; np<diagParams->probe1DStruc.size(); np++) {
        
        every[np]=diagParams->probe1DStruc[np].every;
        unsigned int ndim=params->nDim_particle;

        vector<unsigned int> vecNumber(1);
        vecNumber[0]=diagParams->probe1DStruc[np].number;
        
        probeParticles[np].initialize(vecNumber[0], ndim);
        probeId[np].resize(diagParams->probe1DStruc[np].number);

        vector<double> partPos(ndim*vecNumber[0]);
        for(unsigned int count=0; count!=vecNumber[0]; ++count) {
            int found=smpi->getRank();
            for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
                if (diagParams->probe1DStruc[np].number>1) {
                    partPos[iDim+count*ndim]=diagParams->probe1DStruc[np].posStart[iDim]+count*(diagParams->probe1DStruc[np].posEnd[iDim]-diagParams->probe1DStruc[np].posStart[iDim])/(diagParams->probe1DStruc[np].number-1);
                } else {
                    partPos[iDim+count*ndim]=0.5*(diagParams->probe1DStruc[np].posStart[iDim]+diagParams->probe1DStruc[np].posEnd[iDim]);
                }
                probeParticles[np].position(iDim,count) = 2*M_PI*partPos[iDim+count*ndim];

                if(smpi->getDomainLocalMin(iDim) >  probeParticles[np].position(iDim,count) || smpi->getDomainLocalMax(iDim) <= probeParticles[np].position(iDim,count)) {
                    found=-1;
                }
            }
            probeId[np][count] = found;
        }
        addProbe(np, partPos, vecNumber);
    }
}

