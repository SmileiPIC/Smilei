#include "TimeSelection.h"
#include "PyTools.h"

using namespace std;


// Constructor
TimeSelection::TimeSelection(PyObject* timeSelection, string name)
{
    
    start   = 0;
    end     = maxint;
    period  = 1;
    repeat  = 1;
    spacing = 1;
    
    // If the selection is an int
    if( PyNumber_Check(timeSelection) ) {
        
        // Interpret as the period
        if( !PyTools::convert(timeSelection, period) )
            ERROR(name << ": time selection must be an integer or a list of integers");
        // Default to 1
        if( !period ) period = 1;
    }
    // If the selection is a list or tuple
    else if( PyTuple_Check(timeSelection) || PyList_Check(timeSelection) ) {
        
        PyObject* seq = PySequence_Fast(seq, "expected a sequence");
        int nitems = PySequence_Size(seq);
        
        if ( nitems<1 )
            ERROR(name << ": time selection must have at least one argument");
        if ( nitems>5 )
            ERROR(name << ": time selection must not have more than 5 arguments");
        
        // Extract all values
        vector<int> items (nitems);
        for( int i=0; i<nitems; i++ )
            if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, i), items[i]))
                ERROR(name << ": time selection must be a list of integers");
        
        // Depending on the number of values, they mean different things (see doc)
        if        ( nitems==1 ) {
            if(items[0]) period  = items[0];
        } else if ( nitems==2 ) {
            if(items[0]) start   = items[0];
            if(items[1]) period  = items[1];
        } else if ( nitems==3 ) {
            if(items[0]) start   = items[0];
            if(items[1]) end     = items[1];
            if(items[2]) period  = items[2];
        } else if ( nitems==4 ) {
            if(items[0]) start   = items[0];
            if(items[1]) end     = items[1];
            if(items[2]) period  = items[2];
            if(items[3]) repeat  = items[3];
        } else if ( nitems==5 ) {
            if(items[0]) start   = items[0];
            if(items[1]) end     = items[1];
            if(items[2]) period  = items[2];
            if(items[3]) repeat  = items[3];
            if(items[4]) spacing = items[4];
        }
        
        Py_DECREF(seq);
        
    } else {
        
        ERROR(name << ": time selection must be an integer or a list of integers");
        
    }
    
    // Calculate the group width
    groupWidth = (repeat-1)*spacing + 1;
    
    // Verify a few things
    if( period<1 )
        ERROR(name << ": time selection's period must be > 0");
    if( repeat<1 )
        ERROR(name << ": time selection's repeat must be > 0");
    if( spacing<1 )
        ERROR(name << ": time selection's spacing must be > 0");
    if( groupWidth >= period )
        ERROR(name << ": time selection must have repeat*spacing<period");
    
}


// Tell whether the current timestep is within the selection
bool TimeSelection::theTimeIsNow(int timestep)
{
    
    // Not in selection if outside the start/end bounds
    if( timestep<start || timestep>end ) return false;
    
    // Calculate the remainder to the period
    int r = (timestep - start) % period; 
    
    // The time is now if the remainder is a multiple of the spacing, but still inside the group
    if( r%spacing==0 && r<groupWidth ) return true;
    
    return false;
    
}


// Tell what is the next timestep within the selection
// Returns the same timestep if already within the selection
int TimeSelection::nextTime(int timestep)
{
    
    if( timestep<=start ) return start;
    if( timestep>end    ) return maxint;
    
    int t = timestep-start; // timestep with offset
    int p = t / period;     // current period
    int r = t % period;     // remainder to the current period
    
    // If inside a group
    if( r < groupWidth ) {
        if( r%spacing==0 ) { return timestep; } // return self if already good timestep
        else { return p * period + (r/spacing+1)*spacing; } // otherwise, return next good timestep
    // If after group, return next group's beginning
    } else {
        return (p+1) * period;
    }
    
}


// Tell what is the previous timestep within the selection
// Returns the same timestep if already within the selection
int TimeSelection::previousTime(int timestep)
{
    
    if( timestep<start ) return minint;
    if( timestep>=end  ) return end;
    
    int t = timestep-start; // timestep with offset
    int p = t / period;     // current period
    int r = t % period;     // remainder to the current period
    
    // If inside a group
    if( r < groupWidth ) {
        return p * period + (r/spacing)*spacing; // return previous good timestep
    // If after group, return next group's beginning
    } else {
        return p * period + groupWidth - 1;
    }
    
}

