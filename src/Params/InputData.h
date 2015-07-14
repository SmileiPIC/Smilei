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
    
    //! call the python cleanup function and 
    //! check if python can be closed (e.g. there is no laser python profile)
    //! by calling the _keep_python_running python function (part of pycontrol.pyh)
    void cleanup();
    
private:
    //! passing named command to python
    void pyRunScript(std::string command, std::string name=std::string(""));
    
};

#endif

