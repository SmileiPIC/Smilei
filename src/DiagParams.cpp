#include "DiagParams.h"
#include "Tools.h"
#include <cmath>

using namespace std;

DiagParams::DiagParams(InputData &ifile, PicParams& params) {
	ifile.extract("every",scalar_every,"diagnostic scalar");
	ifile.extract("every",map_every,"diagnostic map");
	ifile.extract("every",probe0d_every,"diagnostic probe0d");
    

    ps_coord.resize(params.nDim_field);
	ifile.extract("x",ps_coord[0],"diagnostic probe0d");
    if (params.nDim_field>1) {
		ifile.extract("y",ps_coord[1],"diagnostic probe0d");
		if (ps_coord[0].size() != ps_coord[1].size()) {
			ERROR("diagnostic probe0d: different dimension of x and y");
		}
	}
	if (params.nDim_field>2) {
		ifile.extract("z",ps_coord[2],"diagnostic probe0d");
		if (ps_coord[0].size() != ps_coord[2].size()) {
			ERROR("diagnostic probe0d: different dimension of x and z");
		}
    }
	
    for (unsigned int k=0; k<ps_coord.size(); k++) {
        for (unsigned int i=0; i<ps_coord[k].size(); i++) {
            if (ps_coord[k][i]<0||ps_coord[k][i]>params.sim_length[k]) ERROR("diagnostic probe0d: probe outside the domain");
			ps_coord[k][i]*=2*M_PI;
            DEBUG(10, "new coordinates " << k << " " << i << " " << ps_coord[k][i]);
		}
	}
}

