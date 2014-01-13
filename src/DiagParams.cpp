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
    
    n_probe1d=0;
    while (ifile.existGroup("diagnostic probe1d",n_probe1d)) {
        vector<vector<double> > ps_1d_c;
        ifile.extract("every",probe1d_every,"diagnostic probe1d",0,n_probe1d);
        ifile.extract("space_res",probe1d_res,"diagnostic probe1d",0,n_probe1d);
        vector<vector<double> > probe1d_ext(params.nDim_field,vector<double>(2));
        double probe1d_length;
        double probe1d_dl = 2*M_PI/probe1d_res;
        
        
        ifile.extract("x",probe1d_ext[0],"diagnostic probe1d",0,n_probe1d);
        probe1d_length=sqrt(pow(probe1d_ext[0][1]-probe1d_ext[0][0],2));
        if (params.nDim_field>1) {
            ifile.extract("y",probe1d_ext[1],"diagnostic probe1d",0,n_probe1d);
            probe1d_length=sqrt(pow(probe1d_ext[0][1]-probe1d_ext[0][0],2)+pow(probe1d_ext[1][1]-probe1d_ext[1][0],2));
            if (probe1d_ext[0].size() != probe1d_ext[1].size()) {
                ERROR("diagnostic probe1d: different dimension of x and y");
            }
        }
        if (params.nDim_field>2) {
            ifile.extract("z",probe1d_ext[2],"diagnostic probe1d",0,n_probe1d);
            probe1d_length=sqrt(pow(probe1d_ext[0][1]-probe1d_ext[0][0],2)+pow(probe1d_ext[1][1]-probe1d_ext[1][0],2)+pow(probe1d_ext[2][1]-probe1d_ext[2][0],2));
            if (probe1d_ext[0].size() != probe1d_ext[1].size()) {
                ERROR("diagnostic probe1d: different dimension of x and z");
            }
        }
        
        unsigned int probe1d_n= (int)(probe1d_length/probe1d_dl)+1;
        probe1d_dl=probe1d_length/(probe1d_n-1);
        double probe1d_dcoor;
        
        
        ps_1d_c.resize(params.nDim_field);
        for(unsigned int k=0;k<params.nDim_field;++k){
            probe1d_dcoor=(probe1d_ext[k][1]-probe1d_ext[k][0])*probe1d_dl/probe1d_length;
            for(unsigned int i=0;i<probe1d_res;++i){
                ps_1d_c[k].push_back(probe1d_ext[k][0]+i*probe1d_dcoor);
            }
        }
        ps_1d_coord.push_back(ps_1d_c);
        ps_1d_every.push_back(probe1d_every);
        n_probe1d++; 
    }
}

