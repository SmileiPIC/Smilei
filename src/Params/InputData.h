/*! @file InputData.h

 @brief InputData.h is the definition of the class InputData which interpretates a namelist-like structure

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
#include "PicParams.h"
#include "SmileiMPI.h"



/*! \brief This is the text parser (similar to namelists).
 It reads once the datafile (at constructor time or later with readFile) and stores the whole read text
 (after being cleaned) in a string variable (namelist) then this variable is passed to all nodes and parsed by filling the structure (allData)
 then you can extract the values with the extract methos (2 templates: one for single variables and one for vectors).
 You can also query the structure with existGroup
*/
class InputData {

public:
    InputData(SmileiMPI*, std::vector<std::string>);
    ~InputData();

    //! string containing the whole clean namelist
    std::string namelist;
        
    //! get bool from python
    bool extract(std::string name, bool &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! get short from python
    bool extract(std::string name, short int &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! get uint from python
    bool extract(std::string name, unsigned int &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! get int from python
    bool extract(std::string name, int &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! get double from python
    bool extract(std::string name, double &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! get string from python
    bool extract(std::string name, std::string &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);

    //! get get vector of uint from python
    bool extract(std::string name, std::vector<unsigned int> &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! get get vector of int from python
    bool extract(std::string name, std::vector<int> &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! get get vector of double from python
    bool extract(std::string name, std::vector<double> &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    	
    //! get vector of string from python
    bool extract(std::string name, std::vector<std::string> &val, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);

    PyObject* extract_py(std::string name, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);

    bool extract_py_profile(ProfileStructure &prof, std::string name, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);

    std::vector<PyObject*> extract_pyVvec(std::string name, std::string group=std::string(""), int occurrenceItem=0, int occurrenceGroup=0);
    
    //! return true if the nth group exists
    bool existGroup(std::string groupName, unsigned int occurrenceGroup=0);
    
    void pyRunScript(std::string, std::string);
    
private:
    // python object: main namelist
    PyObject* py_namelist;
    
    //! print the namelist on stream
    void write(std::ostream&);

};

#endif

