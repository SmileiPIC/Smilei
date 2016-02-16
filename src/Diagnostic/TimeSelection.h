/*
 * TimeSelection.h
 *
 *  Created on: Feb 2016
 */

#ifndef TimeSelection_H
#define TimeSelection_H

#include "PyTools.h"
#include <vector>
#include <limits>
#include <string>

class Params;

//! Class for selecting ranges of times when outputs are done
class TimeSelection {

public:
    //! Constructor
    TimeSelection(PyObject*, std::string);
    
    //! Destructor
    ~TimeSelection(){};
    
    //! Tell whether the current timestep is within the selection
    bool theTimeIsNow(int timestep);
    
    //! Tell what is the next timestep within the selection
    int nextTime(int timestep);
    
    //! Tell what is the previous timestep within the selection
    int previousTime(int timestep);
    
private:
    //! Starting timestep
    int start;
    //! Ending timestep
    int end;
    //! Period between each group
    int period;
    //! Number of repeats inside each group
    int repeat;
    //! Spacing between each repeat
    int spacing;
    
    //! Width of each group
    int groupWidth;
    
    //! Maximum integer
    int maxint = std::numeric_limits<int>::max();
    //! Minimum integer
    int minint = std::numeric_limits<int>::min();
    
};

#endif
