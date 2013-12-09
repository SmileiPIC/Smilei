#ifndef DiagnosticProbe0D_H
#define DiagnosticProbe0D_H

#include "Tools.h"


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include <math.h>
#include "Species.h"
#include "Interpolator.h"
#include "Particle.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;


class DiagnosticProbe0D{
    
public:
    
	DiagnosticProbe0D(PicParams* params, SmileiMPI* smpi,std::vector<std::vector<double> > ps_coord);

	inline void set_p_coor(unsigned int i, std::vector<double> values){ps_coor[i]=values;}//write it in a smarter way
	inline std::vector<std::vector<double> > get_ps_coor(){return ps_coor;}
	inline std::vector<double> get_ps_coor(unsigned int i){return ps_coor[i];}
	
	void set_proc();
    void set_file_name();
    
    void run(int timestep, ElectroMagn* EMfields, Interpolator* interp);
    void open_file();
 
private:
    PicParams* params_;
    SmileiMPI* smpi_;
    std::vector<Particle*> newParticle;
    std::vector<LocalFields> Eloc_fields;
    std::vector<LocalFields> Bloc_fields;
    unsigned int n_probe;
    unsigned int n_probe_loc;
    std::vector<std::vector<double> > ps_coor;
    std::vector<bool> here;

};
#endif 