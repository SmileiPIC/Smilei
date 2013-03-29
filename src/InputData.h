/*! @file InputData.h 
 
 @brief InputData.h is the definition of the class InputData which interpretates a namelist-like structure
 
 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef InputData_H
#define InputData_H

#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include "Tools.h"

#include <vector>

class PicParams;

/*! \brief This is the text parser (similar to namelists). It reads once the datafile (at constructor time or later with parsefile) and stores the read values in a double map dictionary (map_)
*/
class InputData {
	
public:	
	//! parse file
	void parseFile(std::string=std::string());
	
	//! read input data parameters in PicParams
	void InputDataParam(PicParams&);

	//! write namelist onto file
	void write(std::string, std::string=std::string());

	//! print the namelist on stream
	void write(std::ostream& = std::cerr, std::string=std::string());
		
	//! generic template to InputData a single value
	template <class T> void extract(std::string data, T &val, std::string group=std::string("")) {
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		DEBUGEXEC(if(map_[group].find(data)==map_[group].end()) WARNING("Namelist: not found [" << group << "] "<< data));		
		std::stringstream iss(map_[group][data]);
		iss >> std::boolalpha >> val;
	}
	
	//! generic class to InputData a pointer (pay attention that T will be allocated and it's up to you to delete it!) return value is the number of values InputData
	template <class T> size_t extract(std::string data, T*&val, std::string group=std::string("")) {
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		DEBUGEXEC(if(map_[group].find(data)==map_[group].end()) WARNING("Namelist: not found [" << group << "] "<< data));		
		std::stringstream iss(map_[group][data]);	
		std::vector<std::string> strVec;
		copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), 
			 std::back_inserter<std::vector<std::string> >(strVec));
		size_t nums=strVec.size();
		val= new T[nums];
		for (size_t i=0;i<nums;i++) {
			std::stringstream(strVec[i]) >> val[i];
		}	
		return nums;
	}
	
	//! generic class to InputData a vector (pay attention that T will be allocated and it's up to you to delete it!) return value is the number of values InputData
	template <class T> void extract(std::string data, std::vector<T>&val, std::string group=std::string("")) {
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		DEBUGEXEC(if(map_[group].find(data)==map_[group].end()) WARNING("Namelist: not found [" << group << "] "<< data));		
		std::stringstream iss(map_[group][data]);
		std::vector<std::string> strVec;
		copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), 
			 std::back_inserter<std::vector<std::string> >(strVec));
		size_t nums=strVec.size();
		val.resize(nums);
		for (size_t i=0;i<nums;i++) {
			std::stringstream(strVec[i]) >> val[i];
		}	
	}
		
	//! return the list of groups found
	std::vector<std::string> getGroups();

private:
	//! this is a function that removes triling spaces and tabs from the beginning and the end of a string
	std::string cleanString(std::string);
	
	//! this is a map that contains a couple of key values that have been InputData from dataFile file
	//! the first string is the group the second the key and the third is the value(s) 
	std::map < std::string, std::map <std::string,std::string> > map_;

};

#endif

