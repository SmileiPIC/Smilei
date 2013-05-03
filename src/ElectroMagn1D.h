
#ifndef ELECTROMAGN1D_H
#define ELECTROMAGN1D_H

#include "ElectroMagn.h"

//class PicParams;

class ElectroMagn1D : public ElectroMagn {
public:
	ElectroMagn1D(PicParams* params, SmileiMPI* smpi);
	~ElectroMagn1D();
	
	double A_, B_, C_;
	
	double tau1_L_, tau2_L_, tau1_R_, tau2_R_;
	int   las_tordr_L_, las_tordr_R_;
	
	void initMaxwell(SmileiMPI* smpi);
	void chargeConserving(SmileiMPI* smpi);
	void initRhoJ();
	virtual void solveMaxwell(double time_dual, double dt, SmileiMPI* smpi);
	
private:
};

#endif

