#include "DiagParams.h"
#include "Tools.h"
#include <cmath>

using namespace std;

DiagParams::DiagParams() {
}

void DiagParams::parseInputData(InputData &ifile, PicParams& params) {
	ifile.extract("every",scalar_every,"diagnostic scalar");
	ifile.extract("every",map_every,"diagnostic map");
	ifile.extract("every",probe0d_every,"diagnostic probe0d");
<<<<<<< HEAD
    
=======

    ps_coor.resize(params.nDim_field);
	ifile.extract("x",ps_coor[0],"diagnostic probe0d");
    if (params.nDim_field>1) ifile.extract("y",ps_coor[1],"diagnostic probe0d");
	if (params.nDim_field>2) ifile.extract("z",ps_coor[2],"diagnostic probe0d");
    
    for (unsigned int k=0; k<ps_coor.size(); k++) 
        for (unsigned int i=0; i<ps_coor[k].size(); i++) 
            DEBUG(k << " " << i << " " << ps_coor[k][i]);
>>>>>>> d9430c98fb22fc65390629334be5bfeebc514228
}

