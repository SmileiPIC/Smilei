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
class TimeSelection
{

public:
    //! Default Constructor
    TimeSelection( PyObject *, std::string );
    
    //! Empty selection constructor
    TimeSelection();
    
    //! Basic selection constructor
    TimeSelection( int period );
    
    //! Cloning Constructor
    TimeSelection( TimeSelection * );
    
    //! Default Destructor
    ~TimeSelection() {};
    
    //! Tell whether the current timestep is within the selection
    bool theTimeIsNow( int timestep );
    //! Get the last answer of theTimeIsNow(int timestep)
    inline bool theTimeIsNow()
    {
        return TheTimeIsNow;
    };
    
    //! Get the next timestep within the selection
    int nextTime( int timestep );
    //! Get the last answer of nextTime(int timestep)
    inline int nextTime()
    {
        return NextTime;
    };
    
    //! Get the previous timestep within the selection
    int previousTime( int timestep );
    //! Get the last answer of previousTime(int timestep)
    inline int previousTime()
    {
        return PreviousTime;
    };
    
    //! Get the smallest interval between two selected timesteps
    inline int smallestInterval()
    {
        return SmallestInterval;
    };
    
    //! Get the number of occurrences before the given timestep
    int howManyTimesBefore( int timestep );
    
    //! Tell whether the timestep is between start and end
    inline bool inProgress( int timestep )
    {
        return timestep>=start && timestep<=end;
    };
    
    //! Tell whether this is an empty selection (no timesteps)
    inline bool isEmpty()
    {
        return period>0. ? false : true;
    };
    
    //! Set the parameters of the time selection
    void set( double start, double end, double period );
    
    //! Obtain some information about the time selection
    std::string info();
    
private:
    //! Starting timestep
    double start;
    //! Ending timestep
    double end;
    //! Period between each group
    double period;
    //! Number of repeats inside each group
    int repeat;
    //! Spacing between each repeat
    double spacing;
    
    //! Smallest interval between two selected timesteps
    int SmallestInterval;
    //! Width of each group
    double groupWidth;
    
    //! Last answer of theTimeIsNow(int timestep)
    bool TheTimeIsNow;
    //! Last answer of nextTime(int timestep)
    int NextTime;
    //! Last answer of previousTime(int timestep)
    int PreviousTime;
    
};

#endif
