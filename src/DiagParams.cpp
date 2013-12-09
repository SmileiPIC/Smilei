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

    ps_coor.resize(params.nDim_field);
	ifile.extract("x",ps_coor[0],"diagnostic probe0d");
    if (params.nDim_field>1) {
		ifile.extract("y",ps_coor[1],"diagnostic probe0d");
		if (ps_coor[0].size() != ps_coor[1].size()) {
			ERROR("diagnostic probe0d: different dimension of x and y");
		}
	}
	if (params.nDim_field>2) {
		ifile.extract("z",ps_coor[2],"diagnostic probe0d");
		if (ps_coor[0].size() != ps_coor[2].size()) {
			ERROR("diagnostic probe0d: different dimension of x and z");
		}
    }
	
    for (unsigned int k=0; k<ps_coor.size(); k++) {
        for (unsigned int i=0; i<ps_coor[k].size(); i++) {
			ps_coor[k][i]*=2*M_PI;
            DEBUG("new coordinates" << k << " " << i << " " << ps_coor[k][i]);
		}
	}
}

