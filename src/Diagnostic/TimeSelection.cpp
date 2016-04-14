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
        // If zero, no output, ever
        if( !period ) {
            start = maxint;
            return;
        }
        
    }
    // If the selection is a list or tuple
    else if( PyTuple_Check(timeSelection) || PyList_Check(timeSelection) ) {
        
        PyObject* seq = PySequence_Fast(timeSelection, "expected a sequence");
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
        
        ERROR(name << ": time selection (`every`) must be an integer or a list of integers");
        
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
    if( groupWidth > period )
        ERROR(name << ": time selection must have repeat*spacing<period");
    
}


// Cloning Constructor
TimeSelection::TimeSelection(TimeSelection * timeSelection)
{
    start        = timeSelection->start       ;
    end          = timeSelection->end         ;
    period       = timeSelection->period      ;
    repeat       = timeSelection->repeat      ;
    spacing      = timeSelection->spacing     ;
    groupWidth   = timeSelection->groupWidth  ;
    TheTimeIsNow = timeSelection->TheTimeIsNow;
    NextTime     = timeSelection->NextTime    ;
    PreviousTime = timeSelection->PreviousTime;
}


// Tell whether the current timestep is within the selection
bool TimeSelection::theTimeIsNow(int timestep)
{
    TheTimeIsNow = false;
    // In selection if inside the start/end bounds
    if( timestep>=start && timestep<=end ) {
        // Calculate the remainder to the period
        int r = (timestep - start) % period; 
        // The time is now if the remainder is a multiple of the spacing, but still inside the group
        if( r%spacing==0 && r<groupWidth ) TheTimeIsNow = true;
    }
    return TheTimeIsNow;
}


// Tell what is the next timestep within the selection
// Returns the same timestep if already within the selection
int TimeSelection::nextTime(int timestep)
{
    if( timestep<=start ) { 
        NextTime = start;
    } else if( timestep>end ) { 
        NextTime = maxint;
    } else {
        int t = timestep-start; // timestep with offset
        int p = t / period;     // current period
        int r = t % period;     // remainder to the current period
        
        // If inside a group
        if( r < groupWidth ) {
            if( r%spacing==0 ) { NextTime = timestep; } // return current timestep if good
            else { NextTime = start + p * period + (r/spacing+1)*spacing; } // otherwise, return next good timestep
        // If after group, return next group's beginning
        } else {
            NextTime = start + (p+1) * period;
        }
    }
    return NextTime;
}


// Tell what is the previous timestep within the selection
// Returns the same timestep if already within the selection
int TimeSelection::previousTime(int timestep)
{
    if( timestep<start ) {
        PreviousTime = minint;
    } else if( timestep>=end ) {
        PreviousTime = end;
    } else {
        int t = timestep-start; // timestep with offset
        int p = t / period;     // current period
        int r = t % period;     // remainder to the current period
        
        // If inside a group
        if( r < groupWidth ) {
            PreviousTime = start + p * period + (r/spacing)*spacing; // return previous good timestep
        // If after group, return end of that group
        } else {
            PreviousTime = start + p * period + groupWidth - 1;
        }
    }
    return PreviousTime;
}


// Get the number of selected times between min and max timesteps
int TimeSelection::numberOfEvents(int tmin, int tmax)
{
    if( tmin<start ) tmin = start;
    if( tmin>end   ) tmax = end;
    if( tmax<start ) tmax = start;
    if( tmax>end   ) tmax = end;
    
    if( tmax<tmin ) return 0;
    
    int N = 0;
    
    // Work out the period containing tmin
    tmin -= start; // min timestep with offset
    int pmin = tmin / period;   // current period
    int rmin = tmin % period;   // remainder to the current period
    if( rmin < groupWidth ) N += (groupWidth-rmin-1)/spacing + 1; // Add the times in the current group
    
    // Work out the period containing tmax
    tmax -= start; // timesteps with offset
    int pmax = tmax / period;   // current period
    int rmax = tmax % period;   // remainder to the current period
    if( rmax < groupWidth ) N += rmax/spacing + 1; // Add the times in the current group
    else                    N += repeat          ; // or add the full group
    
    // Add the other periods
    N += (pmax-pmin-1) * repeat;
    
    return N;
}

