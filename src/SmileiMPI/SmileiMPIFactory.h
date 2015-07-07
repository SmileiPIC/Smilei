#ifndef SMILEIMPIFACTORY_H
#define SMILEIMPIFACTORY_H

#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"
#include "SmileiMPI_Cart2D.h"

#include "PicParams.h"

#include "Tools.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPIFactory
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPIFactory {
public:
    //  --------------------------------------------------------------------------------------------------------------------
    //! Create appropriate MPI environment for the geometry 
    //! \param params : Parameters
    //! \param smpiData : Initial MPI environment (data broadcast)
    //  --------------------------------------------------------------------------------------------------------------------
    static SmileiMPI* create(PicParams& params, SmileiMPI* smpiData) {
        SmileiMPI* smpi = NULL;
        MESSAGE(1, "Geometry:" << params.geometry);
        if ( params.geometry == "1d3v" ) {
            smpi = new  SmileiMPI_Cart1D(smpiData);
        }
        else if ( params.geometry == "2d3v" ) {
            smpi = new  SmileiMPI_Cart2D(smpiData);
        }
        else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }

        // Creation of a cartesian topology
        smpi->createTopology(params);

	// Creation of derivated datatypes for grid borders
        if ( params.geometry == "2d3v" ) smpi->createType(params);

        return smpi;
    }

};

#endif

