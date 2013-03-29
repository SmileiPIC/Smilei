//**********************************************************************************************************************************
//! Header for the virtual class ElectroMagn
//**********************************************************************************************************************************
#ifndef ELECTROMAGN_H
#define ELECTROMAGN_H

#include <vector>
#include <string>
#include "Tools.h"

class PicParams;
class Species;
class Projector;

class Field;
class Laser;

class ElectroMagn {
public:
	std::vector<unsigned int> dimPrim;
	std::vector<unsigned int> dimDual;
	std::vector<double> dspace;
	std::vector<double> dspacesdt;
	std::vector<double> dtsdspace;
	double time_dual;
	
	
	Field* Ex_;
	Field* Ey_;
	Field* Ez_;
	
	Field* Bx_;
	Field* By_;
	Field* Bz_;
	
	Field* Bx_m;
	Field* By_m;
	Field* Bz_m;
	
	Field* Jx_;
	Field* Jy_;
	Field* Jz_;
	
	Field* rho_;
	Field* rho_o;
	
	std::vector<Laser*> laser_;
	
	
	ElectroMagn(PicParams* paramsdouble);
	virtual ~ElectroMagn();
	
	void dump();
	
	void initRho(std::vector<Species*> vecSpecies, Projector* Proj);
	
	virtual void initMaxwell() = 0;
	virtual void initRhoJ() = 0;
	virtual void chargeConserving() = 0;
	virtual void solveMaxwell(double time_dual, double dt) = 0;
	
private:
};

#endif

