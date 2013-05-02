#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include <cstdlib> 
#include <string>
#include <iostream>
#include <fstream>
#include "Tools.h"

//! Structure containing the fields at a given position (e.g. at a Particle position)
struct LocalFields
{
    //! value of the field component along the x-direction
	double x;
	//! value of the field component along the y-direction
    double y;
    //! value of the field component along the z-direction
	double z;
	
};

//! Class Field: generic class allowing to define vectors
class Field
{
    
public:
    //! link to file for dump
	std::ofstream fdata_;
	
    //! Constructor for Field: with no input argument
	Field() {;};

    //! Constructor for Field: with the Field dimensions as input argument
	Field( std::vector<unsigned int> dims ) {;};

    //! Constructor for Field: with the Field dimensions and dump file name as input argument
	Field( std::vector<unsigned int> dims, std::string name ) {;} ;

    //! Destructor for Field
	virtual ~Field() {;} ;
    
    //! Virtual method used to allocate Field
	virtual void allocateDims(std::vector<unsigned int> dims) = 0;

	//! Virtual method used to make a dump of the Field data
    virtual void dump(std::vector<unsigned int> dims) = 0;

protected:
    //! vector containing the dimensions of the Field
	std::vector<unsigned int> dims_;
	
private:
	
};

#endif
