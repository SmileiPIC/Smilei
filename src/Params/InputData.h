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
#include "PyTools.h"
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
        
    PyObject* extract_py(std::string name, std::string component=std::string(""), int nComponent=0);
    
    std::vector<PyObject*> extract_pyVec(std::string name, std::string component=std::string(""), int nComponent=0);
    
    //! get T from python
    template< typename T>
    bool extract(std::string name, T &val, std::string component=std::string(""), int nComponent=0) {
        PyObject* py_val = extract_py(name,component,nComponent);
        return PyTools::convert(py_val,val);
    }
    
    template< typename T>
    bool extract(std::string name, std::vector<T> &val, std::string component=std::string(""), int nComponent=0) {
        std::vector<PyObject*> py_val = extract_pyVec(name,component,nComponent);
        if (py_val.size())
            return PyTools::convert(py_val,val);
        return false;
    }
    
    //! return true if the nth component exists
    bool existComponent(std::string componentName, unsigned int nComponent=0);
    
    void pyRunScript(std::string, std::string);
    
private:
    // python object: main namelist
    PyObject* py_namelist;
    
    //! print the namelist on stream
    void write(std::ostream&);

};

#endif

