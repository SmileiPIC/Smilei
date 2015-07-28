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



/*! \brief This is namelist helper.
 Namelist must be a valid python script.
*/
class InputData {

public:
    InputData(SmileiMPI*, std::vector<std::string>);
    ~InputData();
    
    //! string containing the whole clean namelist
    std::string namelist;
    
    //! extract python object from namelist
    PyObject* extract_py(std::string name, std::string component=std::string(""), int nComponent=0);
    
    //! extract pytohn vector of objects from namelist
    std::vector<PyObject*> extract_pyVec(std::string name, std::string component=std::string(""), int nComponent=0);
    
    //! get T from python
    template< typename T>
    bool extract(std::string name, T &val, std::string component=std::string(""), int nComponent=0) {
        PyObject* py_val = extract_py(name,component,nComponent);
        PyTools::checkPyError();        
        return PyTools::convert(py_val,val);
    }
    
    //! extract vectors
    template< typename T>
    bool extract(std::string name, std::vector<T> &val, std::string component=std::string(""), int nComponent=0) {
        std::vector<PyObject*> py_val = extract_pyVec(name,component,nComponent);
        if (py_val.size())
            return PyTools::convert(py_val,val);
        return false;
    }
    
    //! return the number of components (see pyinit.py)
    int nComponents(std::string componentName);
    
    //! check if python can be closed (e.g. there is no laser python profile)
    void cleanup();
    
private:
    //! passing named command to python
    void pyRunScript(std::string command, std::string name=std::string(""));
    
    // python object: main namelist
    PyObject* py_namelist;
    
    //! print the namelist on stream
    void write(std::ostream&);

};

#endif

