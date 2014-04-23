/*! @file InputData.h

 @brief InputData.h is the definition of the class InputData which interpretates a namelist-like structure

 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef INPUTDATA_H
#define INPUTDATA_H

#include <cstdlib>

#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <iterator>
#include <vector>

#include "Tools.h"

class PicParams;

/*! \brief This is the text parser (similar to namelists).
 It reads once the datafile (at constructor time or later with parsefile) and stores the whole read text
 (after being cleaned) in a string variable (namelist) then this variable is passed to all nodes and parsed by filling the structure (allData)
 then you can extract the values with the extract methos (2 templates: one for single variables and one for vectors).
 You can also query the structure with existGroup
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
    void write() {
        write(std::cerr);
    };

    //! generic template to InputData a single value
    template <class T> bool extract(std::string data, T &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0) {
        data=cleanString(data);
        group=cleanString(group);
        int n_occur_group=0;
        for (unsigned int i=0; i<allData.size(); i++) {
            if (group == allData[i].first) {
                if (occurrenceGroup==n_occur_group || occurrenceGroup < 0) {
                    int n_occur_item=0;
                    for (unsigned int j=0; j<allData[i].second.size(); j++) {
                        if (data == allData[i].second[j].first) {
                            if (occurrenceItem==n_occur_item || occurrenceItem < 0) {
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
        DEBUG(10,"NOT FOUND! searching for scalar \"" << data << "\" [" << occurrenceItem << "] in group \"" << group << "\" [" << occurrenceGroup << "]");
        return false;
    }


    //! generic class to InputData a vector (pay attention that T will be allocated and it's up to you to delete it!) return value is the number of values InputData
    template <class T> bool extract(std::string data, std::vector<T>&val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0) {
        data=cleanString(data);
        group=cleanString(group);
        bool found=false;
        int n_occur_group=0;
        for (unsigned int i=0; i<allData.size(); i++) {
            if (group == allData[i].first) {
                if (occurrenceGroup==n_occur_group || occurrenceGroup < 0) {
                    int n_occur_item=0;
                    for (unsigned int j=0; j<allData[i].second.size(); j++) {
                        if (data == allData[i].second[j].first) {
                            if (occurrenceItem==n_occur_item || occurrenceItem < 0) {
                                std::stringstream iss(allData[i].second[j].second);
                                std::vector<std::string> strVec;
                                copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                                     std::back_inserter<std::vector<std::string> >(strVec));
                                size_t nums=strVec.size();
                                val.resize(nums);
                                for (size_t i=0; i<nums; i++) {
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
        DEBUG(10,"NOT FOUND! searching for vector \"" << data << "\" [" << occurrenceItem << "] in group \"" << group << "\" [" << occurrenceGroup << "]");
        return found;

    }

	bool extract(std::string data, std::vector<std::string>&val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0) {
        data=cleanString(data);
        group=cleanString(group);
        bool found=false;
        int n_occur_group=0;
        for (unsigned int i=0; i<allData.size(); i++) {
            if (group == allData[i].first) {
                if (occurrenceGroup==n_occur_group || occurrenceGroup < 0) {
                    int n_occur_item=0;
                    for (unsigned int j=0; j<allData[i].second.size(); j++) {
                        if (data == allData[i].second[j].first) {
                            if (occurrenceItem==n_occur_item || occurrenceItem < 0) {
                                std::stringstream iss(allData[i].second[j].second);
								do {
									std::string sub;
									iss >> sub;
									sub = cleanString(sub);
									if (!sub.empty()) val.push_back(sub);
								} while (iss);								
                                return true;
                            }
                            n_occur_item++;
                        }
                    }
                }
                n_occur_group++;
            }
        }
        DEBUG(10,"NOT FOUND! searching for vector \"" << data << "\" [" << occurrenceItem << "] in group \"" << group << "\" [" << occurrenceGroup << "]");
        return found;
		
    }
	
    //! return true if the nth group exists
    bool existGroup(std::string groupName, unsigned int occurrenceGroup=0);

    //! string containing the whole clean namelist
    std::string namelist;

    template <class T> void addVar(std::string nameVar, T &valVar,std::string nameGroup="") {
        std::vector< std::pair <std::string,std::string> > myvec;
        std::ostringstream s;
        s << valVar;
        myvec.push_back(make_pair(nameVar,s.str()));
        allData.push_back(make_pair(nameGroup,myvec));
    };

private:
    //! print the namelist on stream
    void write(std::ostream&);

    //! this is a function that removes trailing spaces and tabs from the beginning and the end of a string (and transforms in lowercase)
    std::string cleanString(std::string);

    //! this is a vector of pairs string, vector of pairs string,string....
    //! the first string is the name of the group and the second vector of pairs contains the variable name and it's value (as string)
    std::vector<std::pair<std::string , std::vector< std::pair <std::string,std::string> > > > allData;

};

#endif

