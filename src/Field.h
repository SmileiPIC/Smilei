
#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <cstdlib> 
#include <string>
#include <iostream>
#include <fstream>
#include "Tools.h"

struct chLocaux {
	double x;
	double y;
	double z;
	
};

class Field {
public:
	std::ofstream fdata_;
	
	Field() {;} ; 
	Field( std::vector<unsigned int> dims ) {;} ; 
	Field( std::vector<unsigned int> dims, std::string name ) {;} ; 
	virtual ~Field() {;} ;
	virtual void allocateDims(std::vector<unsigned int> dims) = 0;
	virtual void dump(std::vector<unsigned int> dims) = 0;

protected:
	std::vector<unsigned int> dims_;
	
private:
	
};

#endif

