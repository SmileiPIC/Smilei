#include "DiagnosticProbe2D.h"

#include <iomanip>
#include <string>
#include <iomanip>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field.h"

using namespace std;

DiagnosticProbe2D::DiagnosticProbe2D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi) : DiagnosticProbe(params,diagParams,smpi,2) {
    every.resize(diagParams->probe2DStruc.size());
    probeParticles.resize(diagParams->probe2DStruc.size());
    probeId.resize(diagParams->probe2DStruc.size());

    for (unsigned int np=0; np<diagParams->probe2DStruc.size(); np++) {
        
        every[np]=diagParams->probe2DStruc[np].every;
        unsigned int ndim=params->nDim_particle;

        vector<unsigned int> vecNumber=diagParams->probe2DStruc[np].number;
        
        probeParticles[np].initialize(vecNumber[0]*vecNumber[1], ndim);
        probeId[np].resize(vecNumber[0]*vecNumber[1]);

        vector<double> partPos(ndim*vecNumber[0]*vecNumber[1]);
        for(unsigned int count1=0; count1!=vecNumber[0]; ++count1) {
            for(unsigned int count2=0; count2!=vecNumber[1]; ++count2) {
                unsigned int count=count2+count1*vecNumber[1];
                
                int found=smpi->getRank();
                for(unsigned int iDim=0; iDim!=ndim; ++iDim) {
                    unsigned int k=iDim+count*ndim;
                    
                    partPos[k] = diagParams->probe2DStruc[np].pos[iDim]+count1*(diagParams->probe2DStruc[np].posFirst[iDim]-diagParams->probe2DStruc[np].pos[iDim])/(vecNumber[0]-1) +
                                 diagParams->probe2DStruc[np].pos[iDim]+count2*(diagParams->probe2DStruc[np].posSecond[iDim]-diagParams->probe2DStruc[np].pos[iDim])/(vecNumber[1]-1);
                                        
                    probeParticles[np].position(iDim,count1) = 2*M_PI*partPos[k];
                    
                    if(smpi->getDomainLocalMin(iDim) >  probeParticles[np].position(iDim,count1) || smpi->getDomainLocalMax(iDim) <= probeParticles[np].position(iDim,count1)) {
                        found=-1;
                    }
                    cerr << smpi->getRank() << " " << iDim  << " " << count1 << " " << count2 << " " << count << " " << k << " " << partPos[k] << " --- ";
                }
                cerr << found << endl;
                probeId[np][count1] = found;
            }
        }
        addProbe(np, partPos, vecNumber);
    }
}

