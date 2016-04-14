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
    //! Default Constructor
    TimeSelection(PyObject*, std::string);
    
    //! Cloning Constructor
    TimeSelection(TimeSelection*);
    
    //! Default Destructor
    ~TimeSelection(){};
    
    //! Tell whether the current timestep is within the selection
    bool theTimeIsNow(int timestep);
    //! Get the last answer of theTimeIsNow(int timestep)
    inline bool theTimeIsNow() { return TheTimeIsNow; };
    
    //! Get the next timestep within the selection
    int nextTime(int timestep);
    //! Get the last answer of nextTime(int timestep)
    inline int nextTime() { return NextTime; };
    
    //! Get the previous timestep within the selection
    int previousTime(int timestep);
    //! Get the last answer of previousTime(int timestep)
    inline int previousTime() { return PreviousTime; };
    
    //! Get the smallest interval between two selected timesteps
    inline int smallestInterval() { return (repeat==1) ? period : spacing; };
    
    //! Tell whether the timestep is between start and end
    inline bool inProgress(int timestep) { return timestep>=start && timestep<=end; };
    
    //! Tell whether this is an empty selection (no timesteps)
    inline bool isEmpty() { return period ? false : true; };
    
    //! Get the number of times in the selection
    int numberOfEvents(int tmin, int tmax);
    
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
    
    //! Last answer of theTimeIsNow(int timestep)
    bool TheTimeIsNow;
    //! Last answer of nextTime(int timestep)
    int NextTime;
    //! Last answer of previousTime(int timestep)
    int PreviousTime;
    
};

#endif
