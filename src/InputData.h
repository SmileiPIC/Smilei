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

/*! \brief This is the text parser (similar to namelists). It reads once the datafile (at constructor time or later with parsefile) and stores the read values in a double map dictionary (allData)
*/
class InputData {
	
public:	
	InputData();
	//! parse file
	void parseFile(std::string=std::string());
	
	//! parse stringstream
	void parseStream();
	
	//! read input data parameters in PicParams
	void InputDataParam(PicParams&);

	//! write namelist onto file (or cerr if unable)
	void write(std::string);
	
	//! write namelist onto cerr
	void write(){write(std::cerr);};
	
	//! generic template to InputData a single value
	template <class T> bool extract(std::string data, T &val, std::string group=std::string(""), unsigned int occurrenceItem=0, unsigned int occurrenceGroup=0) {
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		transform(group.begin(), group.end(), group.begin(), ::tolower);
		unsigned int n_occur_group=0;
		for (unsigned int i=0; i<allData.size(); i++) {
			if (group == allData[i].first) {
				if (occurrenceGroup==n_occur_group) {
					unsigned int n_occur_item=0;
					for (unsigned int j=0; j<allData[i].second.size(); j++) {
						if (data == allData[i].second[j].first) {
							if (occurrenceItem==n_occur_item) {
								std::stringstream iss(allData[i].second[j].second);					
								iss >> std::boolalpha >> val;
								DEBUG(100,"scalar " << data << "[" << occurrenceItem << "] g:" << group << " [" << occurrenceGroup << "] = " << val );
								return true;
							}
							n_occur_item++;
						}
					}
				}
				n_occur_group++;
			}
		}
		DEBUG(10,"============================== NOT FOUND! searching for scalar \"" << data << "\" [" << occurrenceItem << "] in group \"" << group << "\" [" << occurrenceGroup << "]");
		return false;
	}

	
	//! generic class to InputData a vector (pay attention that T will be allocated and it's up to you to delete it!) return value is the number of values InputData
	template <class T> bool extract(std::string data, std::vector<T>&val, std::string group=std::string(""), unsigned int occurrenceItem=0, unsigned int occurrenceGroup=0) {
		transform(data.begin(), data.end(), data.begin(), ::tolower);
		transform(group.begin(), group.end(), group.begin(), ::tolower);
		bool found=false;
		unsigned int n_occur_group=0;
		for (unsigned int i=0; i<allData.size(); i++) {
			if (group == allData[i].first) {
				if (occurrenceGroup==n_occur_group) {
					unsigned int n_occur_item=0;
					for (unsigned int j=0; j<allData[i].second.size(); j++) {
						if (data == allData[i].second[j].first) {
							if (occurrenceItem==n_occur_item) {
								std::stringstream iss(allData[i].second[j].second);					
								std::vector<std::string> strVec;
								copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), 
									 std::back_inserter<std::vector<std::string> >(strVec));
								size_t nums=strVec.size();
								val.resize(nums);
								for (size_t i=0;i<nums;i++) {
									std::stringstream(strVec[i]) >> val[i];
									DEBUG(100,"vector " << data << "[" << occurrenceItem << "] g:" << group << " [" << occurrenceGroup << "]" << val[i]);
								}	
								return true;
							}
							n_occur_item++;
						}
					}
				}
				n_occur_group++;
			}
		}
		DEBUG(10,"============================== NOT FOUND! searching for vector \"" << data << "\" [" << occurrenceItem << "] in group \"" << group << "\" [" << occurrenceGroup << "]");
		return found;

	}
	
	//! return true if the nth group exists
	bool existGroup(std::string groupName, unsigned int occurrenceGroup=0);

	//! return true if the key exists in the nth group
	bool existKey(std::string key, std::string groupName="", unsigned int occurrenceItem=0, unsigned int occurrenceGroup=0);

//	//! return the list of groups found
//	std::vector<std::string> getGroups();
	
	//! string containing the whole clean namelist
	std::string namelist;

private:
	//! print the namelist on stream
	void write(std::ostream&); 
		
	//! this is a function that removes triling spaces and tabs from the beginning and the end of a string
	std::string cleanString(std::string);
	
	//! this is a vector of pairs string,map 
	//! the string is the name of the group and the map key contains the name of the veriable and the string is it's value
	//test new branch
	std::vector<std::pair<std::string , std::vector< std::pair <std::string,std::string> > > > allData;

};

#endif

