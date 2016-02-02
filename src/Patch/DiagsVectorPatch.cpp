
#include "DiagsVectorPatch.h"
#include "VectorPatch.h"

#include <cstring>
#include <vector>

using namespace std;

void DiagsVectorPatch::computeGlobalDiags( VectorPatch& vecPatches, int timestep)
{
    
    DiagsVectorPatch::computeScalarsDiags  ( vecPatches, timestep );
    DiagsVectorPatch::computePhaseSpace    ( vecPatches );
    DiagsVectorPatch::computeParticlesDiags( vecPatches, timestep );
    
}

void DiagsVectorPatch::computeScalarsDiags( VectorPatch& vecPatches, int timestep)
{
    //cout << "In Global Compute Scalar Diags " << vecPatches(0)->Diags->scalars.every << " \t timestep = " << timestep << endl;
    int scalars_every( vecPatches(0)->Diags->scalars.every );
    if (timestep % scalars_every != 0) return;

    //cout << "In Global Compute Scalar Daigs\n";

    //std::vector<std::pair<std::string,double> > out_list;
    //std::vector<std::string> out_key;
    //std::vector<double>      out_value;
    //std::vector<unsigned int> out_width;
    //std::vector<std::pair<std::string,double> >::iterator itDiagScalar;


    int nDiags( vecPatches(0)->Diags->scalars.out_value.size() );
    // Initialize scalars iterator on 1st diag
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
	vecPatches(ipatch)->Diags->scalars.itDiagScalarName  =  vecPatches(ipatch)->Diags->scalars.out_key.begin();
	vecPatches(ipatch)->Diags->scalars.itDiagScalarValue =  vecPatches(ipatch)->Diags->scalars.out_value.begin();
    }


    for (int idiags = 0 ; idiags<nDiags ; idiags++) {
	string diagName( *vecPatches(0)->Diags->scalars.itDiagScalarName );

	if ( ( diagName.find("Min") == std::string::npos ) && ( diagName.find("Max") == std::string::npos ) ) {
	    double sum(0.);
	    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
		sum += *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue;
		if (ipatch) {
		    vecPatches(ipatch)->Diags->scalars.itDiagScalarName++;
		    vecPatches(ipatch)->Diags->scalars.itDiagScalarValue++;
		}
	    }
	    *vecPatches(0)->Diags->scalars.itDiagScalarValue = sum;
	    vecPatches(0)->Diags->scalars.itDiagScalarName++;
	    vecPatches(0)->Diags->scalars.itDiagScalarValue++;
	}
	else if ( diagName.find("MinCell") != std::string::npos ) {
	    vector<double>::iterator iterVal    = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
	    vector<double>::iterator iterValRef = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
	    double min( *iterValRef );

	    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
		if ( *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue < min ) {
		    min = *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue;
		    iterVal = vecPatches(ipatch)->Diags->scalars.itDiagScalarValue-1;
		}
		if (ipatch) {
		    vecPatches(ipatch)->Diags->scalars.itDiagScalarName++;
		    vecPatches(ipatch)->Diags->scalars.itDiagScalarValue++;
		}
	    }
	    *vecPatches(0)->Diags->scalars.itDiagScalarValue = min;
	    iterValRef = iterVal;

	    vecPatches(0)->Diags->scalars.itDiagScalarName++;	    
	    vecPatches(0)->Diags->scalars.itDiagScalarValue++;	    
	}
	else if ( diagName.find("MaxCell") != std::string::npos ) {
	    vector<double>::iterator iterVal    = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
	    vector<double>::iterator iterValRef = vecPatches(0)->Diags->scalars.itDiagScalarValue-1;
	    double max( *iterValRef );

	    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
		if ( *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue > max ) {
		    max = *vecPatches(ipatch)->Diags->scalars.itDiagScalarValue;
		    iterVal = vecPatches(ipatch)->Diags->scalars.itDiagScalarValue-1;
		}
		if (ipatch) {
		    vecPatches(ipatch)->Diags->scalars.itDiagScalarName++;
		    vecPatches(ipatch)->Diags->scalars.itDiagScalarValue++;
		}
	    }
	    *vecPatches(0)->Diags->scalars.itDiagScalarValue = max;
	    iterValRef = iterVal;

	    vecPatches(0)->Diags->scalars.itDiagScalarName++;	    
	    vecPatches(0)->Diags->scalars.itDiagScalarValue++;	    
	}

	// Go to next diag
    }

    // After MPI sync
    //vecPatches(0)->Diags->scalars.write(timestep);

}

void DiagsVectorPatch::computePhaseSpace( VectorPatch& vecPatches )
{
    // A définir : DiagPhaseSpace::itDiagPhase

    int nDiags( vecPatches(0)->Diags->phases.vecDiagPhaseToRun.size() );

    // Initialize scalars iterator on 1st diag
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
	vecPatches(ipatch)->Diags->phases.itDiagPhase =  vecPatches(ipatch)->Diags->phases.vecDiagPhaseToRun.begin();
    
    for (int idiags = 0 ; idiags<nDiags ; idiags++) {
	vector<unsigned int> diagSize = (*vecPatches(0)->Diags->phases.itDiagPhase)->my_data.dims_;
	for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
	    for (int i=0 ; i<diagSize[0] ; i++)
		for (int j=0 ; j<diagSize[1] ; j++)
		    (*vecPatches(0)->Diags->phases.itDiagPhase)->my_data(i,j) += (*vecPatches(ipatch)->Diags->phases.itDiagPhase)->my_data(i,j);
	    vecPatches(ipatch)->Diags->phases.itDiagPhase++;
	} // for ipatch
	vecPatches(0)->Diags->phases.itDiagPhase++;

    } // for idiags

    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++)
	vecPatches(ipatch)->Diags->phases.vecDiagPhaseToRun.clear();

}


void DiagsVectorPatch::computeParticlesDiags( VectorPatch& vecPatches, int timestep)
{
    int nDiags( vecPatches(0)->Diags->vecDiagnosticParticles.size() );

    for (int idiags = 0 ; idiags<nDiags ; idiags++) {
	if (timestep % vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->every != vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->time_average-1) continue;

	int output_size = vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->output_size;
	for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++) {
	    for (int i=0 ; i<output_size ; i++)
		vecPatches(0)->Diags->vecDiagnosticParticles[idiags]->data_sum[i] += vecPatches(ipatch)->Diags->vecDiagnosticParticles[idiags]->data_sum[i];
	} // for ipatch

    } // for idiags

    
    for (unsigned int ipatch=1 ; ipatch<vecPatches.size() ; ipatch++)
	for (unsigned int i=0; i<vecPatches(ipatch)->Diags->vecDiagnosticParticles.size(); i++)
	       if (vecPatches(ipatch)->Diags->vecDiagnosticParticles[i]->time_average == 1)
		   vecPatches(ipatch)->Diags->vecDiagnosticParticles[i]->clean();

}
