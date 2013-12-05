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
	InputData(std::string);
	//! parse file
	void parseFile(std::string=std::string());
	
	//! read input data parameters in PicParams
	void InputDataParam(PicParams&);

	//! write namelist onto file
	void write(std::string, std::string=std::string());

	//! print the namelist on stream
	void write(std::ostream& = std::cerr, std::string=std::string());
		
	//! generic template to InputData a single value
	template <class T> bool extract(std::string data, T &val, std::string group=std::string(""), unsigned int occurrence=0) {
		//		transform(data.begin(), data.end(), data.begin(), ::tolower);
		//		DEBUGEXEC(if(map_[group].find(data)==map_[group].end()) WARNING("Namelist: not found [" << group << "] "<< data));		
		//		std::stringstream iss(map_[group][data]);
		//		iss >> std::boolalpha >> val;
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		bool found=false;
		unsigned int n_occur=0;
		for (unsigned int i=0; i<map_.size(); i++) {
			if (group == map_[i].first) {
				if (occurrence==n_occur) {
					std::stringstream iss(map_[i].second[data]);					
					iss >> std::boolalpha >> val;
					found=true;
				}
				n_occur++;
			}
		}
		return found;
	}

	
	// REMOVED BECAUSE UNSAFE!!!!
	//! generic class to InputData a pointer (pay attention that T will be allocated and it's up to you to delete it!) return value is the number of values InputData
//	template <class T> size_t extract(std::string data, T*&val, std::string group=std::string("")) {
//		transform(data.begin(), data.end(), data.begin(), ::tolower);
//		DEBUGEXEC(if(map_[group].find(data)==map_[group].end()) WARNING("Namelist: not found [" << group << "] "<< data));		
//		std::stringstream iss(map_[group][data]);	
//		std::vector<std::string> strVec;
//		copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), 
//			 std::back_inserter<std::vector<std::string> >(strVec));
//		size_t nums=strVec.size();
//		val= new T[nums];
//		for (size_t i=0;i<nums;i++) {
//			std::stringstream(strVec[i]) >> val[i];
//		}	
//		return nums;
//	}
	
	//! generic class to InputData a vector (pay attention that T will be allocated and it's up to you to delete it!) return value is the number of values InputData
	template <class T> bool extract(std::string data, std::vector<T>&val, std::string group=std::string(""), unsigned int occurrence=0) {
//		transform(data.begin(), data.end(), data.begin(), ::tolower);
//		DEBUGEXEC(if(map_[group].find(data)==map_[group].end()) WARNING("Namelist: not found [" << group << "] "<< data));		
//		std::stringstream iss(map_[group][data]);
//		std::vector<std::string> strVec;
//		copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), 
//			 std::back_inserter<std::vector<std::string> >(strVec));
//		size_t nums=strVec.size();
//		val.resize(nums);
//		for (size_t i=0;i<nums;i++) {
//			std::stringstream(strVec[i]) >> val[i];
//		}
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		bool found=false;
		unsigned int n_occur=0;
		for (unsigned int i=0; i<map_.size(); i++) {
			if (group == map_[i].first) {
				if (occurrence==n_occur) {
					std::stringstream iss(map_[i].second[data]);					
					std::vector<std::string> strVec;
					copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), 
						 std::back_inserter<std::vector<std::string> >(strVec));
					size_t nums=strVec.size();
					val.resize(nums);
					for (size_t i=0;i<nums;i++) {
						std::stringstream(strVec[i]) >> val[i];
					}	
					found=true;
				}
				n_occur++;
			}
		}
		
		return found;
	}
	
	//! return true if the nth group exists
	bool existGroup(std::string groupName, unsigned int occurrence=0);

	//! return true if the key exists in the nth group
	bool existKey(std::string key, std::string groupName="", unsigned int occurrence=0);

	//! return the list of groups found
	std::vector<std::string> getGroups();

private:
	//! this is a function that removes triling spaces and tabs from the beginning and the end of a string
	std::string cleanString(std::string);
	
	//! this is a vector of pairs string,map 
	//! the string is the name of the group and the map key contains the name of the veriable and the string is it's value
	std::vector<std::pair<std::string , std::map <std::string,std::string> > > map_;
};

#endif

