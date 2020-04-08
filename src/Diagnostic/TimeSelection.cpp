#include "PyTools.h"
#include "TimeSelection.h"
#include <math.h>
#include <sstream>
#include <iomanip>

using namespace std;


// Constructor
TimeSelection::TimeSelection( PyObject *timeSelection, string name )
{

    start   = 0.;
    end     = std::numeric_limits<double>::max();
    period  = 1.;
    repeat  = 1;
    spacing = 1.;
    
    // If the selection is an int
    if( PyNumber_Check( timeSelection ) ) {
    
        // Interpret as the period
        if( !PyTools::py2scalar( timeSelection, period ) ) {
            ERROR( name << ": time selection must be a number or a list of numbers" );
        }
        // If zero, no output, ever
        if( !period ) {
            start = std::numeric_limits<double>::max();
            return;
        }
        
    }
    // If the selection is a list or tuple
    else if( PyTuple_Check( timeSelection ) || PyList_Check( timeSelection ) ) {
    
        PyObject *seq = PySequence_Fast( timeSelection, "expected a sequence" );
        int nitems = PySequence_Size( seq );
        
        if( nitems<1 ) {
            ERROR( name << ": time selection must have at least one argument" );
        }
        if( nitems>5 ) {
            ERROR( name << ": time selection must not have more than 5 arguments" );
        }
        
        // Extract all values
        vector<double> items( nitems );
        for( int i=0; i<nitems; i++ )
            if( !PyTools::py2scalar( PySequence_Fast_GET_ITEM( seq, i ), items[i] ) ) {
                ERROR( name << ": time selection must be a list of numbers" );
            }
            
        // Depending on the number of values, they mean different things (see doc)
        if( nitems==1 ) {
            if( items[0] ) {
                period  = items[0];
            }
        } else if( nitems==2 ) {
            if( items[0] ) {
                start   = items[0];
            }
            if( items[1] ) {
                period  = items[1];
            }
        } else if( nitems==3 ) {
            if( items[0] ) {
                start   = items[0];
            }
            if( items[1] ) {
                end     = items[1];
            }
            if( items[2] ) {
                period  = items[2];
            }
        } else if( nitems==4 ) {
            if( items[0] ) {
                start   = items[0];
            }
            if( items[1] ) {
                end     = items[1];
            }
            if( items[2] ) {
                period  = items[2];
            }
            if( items[3] ) {
                repeat  = ( int )round( items[3] );
            }
        } else if( nitems==5 ) {
            if( items[0] ) {
                start   = items[0];
            }
            if( items[1] ) {
                end     = items[1];
            }
            if( items[2] ) {
                period  = items[2];
            }
            if( items[3] ) {
                repeat  = ( int )round( items[3] );
            }
            if( items[4] ) {
                spacing = items[4];
            }
        }
        
        Py_DECREF( seq );
        
    } else {
    
        ERROR( name << ": time selection (`every`) must be a number or a list of numbers" );
        
    }
    
    // Calculate the smallest interval
    SmallestInterval = ( repeat<=1 ) ? ( ( int )period ) : ( ( int )spacing );
    // Calculate the group width
    groupWidth = ( repeat-1 )*spacing + 1.;
    
    // Verify a few things
    if( period<1. ) {
        ERROR( name << ": time selection's period must be >= 1." );
    }
    if( repeat<1 ) {
        ERROR( name << ": time selection's repeat must be >= 1" );
    }
    if( spacing<1. ) {
        ERROR( name << ": time selection's spacing must be >= 1." );
    }
    if( groupWidth > period ) {
        ERROR( name << ": time selection must have repeat*spacing<period" );
    }
    
    Py_DECREF( timeSelection );
}


// Empty time selection
TimeSelection::TimeSelection()
{
    start   = std::numeric_limits<double>::max();
    end     = std::numeric_limits<double>::max();
    period  = 0.;
    repeat  = 1;
    spacing = 1.;
    groupWidth   = 0.;
    TheTimeIsNow = false;
    NextTime     = std::numeric_limits<int>::max();
    PreviousTime = std::numeric_limits<int>::max();
}

// Basic time selection
TimeSelection::TimeSelection( int copy_period )
{
    start         = 0;
    end           = std::numeric_limits<double>::max();
    period        = copy_period;
    repeat        = 1;
    spacing       = 1.;
    groupWidth    = 0.;
    TheTimeIsNow  = false;
    NextTime      = std::numeric_limits<int>::max();
    PreviousTime  = 0;
}

// Cloning Constructor
TimeSelection::TimeSelection( TimeSelection *timeSelection )
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
bool TimeSelection::theTimeIsNow( int timestep )
{
    TheTimeIsNow = false;
    // In selection if inside the start/end bounds
    if( timestep>=round( start ) && timestep<=round( end ) ) {
        // Calculate the number of timesteps since the start
        double t = ( double )timestep - start + 0.4999;
        // Calculate the time from the previous period
        t -= floor( t/period )*period;
        // If within group, calculate the time to the previous repeat
        if( t < groupWidth ) {
            t -= floor( t/spacing )*spacing;
        }
        // The time is now if closest repeat it within 0.5
        if( t < 1. ) {
            TheTimeIsNow = true;
        }
    }
    return TheTimeIsNow;
}


// Tell what is the next timestep within the selection
// Returns the same timestep if already within the selection
int TimeSelection::nextTime( int timestep )
{
    if( timestep<=round( start ) ) {
        NextTime = start;
    } else if( timestep>round( end ) ) {
        NextTime = std::numeric_limits<int>::max();
    } else {
        double t = ( double )( timestep )-start + 0.4999; // number of timesteps since the start
        double p = floor( t/period )*period; // previous period
        double T = t - p; // time to the previous period
        
        // If within group
        if( T < groupWidth ) {
            double r = floor( T/spacing )*spacing; // previous repeat
            if( T-r < 1. ) {
                NextTime = timestep;    // return current timestep if good
            } else {
                NextTime = ( int ) round( start + p + r + spacing );    // otherwise, return next repeat
            }
            // If after group, return next group's beginning
        } else {
            NextTime = ( int ) round( start + p+period );
        }
    }
    return NextTime;
}


// Tell what is the previous timestep within the selection
// Returns the same timestep if already within the selection
int TimeSelection::previousTime( int timestep )
{
    if( timestep<round( start ) ) {
        PreviousTime = std::numeric_limits<int>::min();
    } else if( timestep>=round( end ) ) {
        PreviousTime = end;
    } else {
        double t = ( double )( timestep )-start + 0.4999; // number of timesteps since the start
        double p = floor( t/period )*period; // previous period
        double T = t - p; // time to the previous period
        
        // If within group
        if( T < groupWidth ) {
            PreviousTime = ( int ) round( start + p + floor( T/spacing )*spacing ); // return previous good timestep
            // If after group, return end of that group
        } else {
            PreviousTime = ( int ) round( start + p + groupWidth - 1 );
        }
    }
    return PreviousTime;
}


// Tell how many selected timesteps have occured before the timestep requested
int TimeSelection::howManyTimesBefore( int timestep )
{
    int nt = 0;
    if( timestep>=round( start ) ) {
        if( timestep>=round( end ) ) {
            timestep = round( end );
        }
        double t = ( double )( timestep )-start + 0.4999; // number of timesteps since the start
        double p = floor( t/period ); // previous period number
        nt = p * repeat; // number of occurences before this period
        double T = t - p*period; // time to the previous period
        
        // If within group
        if( T < groupWidth ) {
            nt += ( int ) floor( T/spacing ); // add number of repeats to account for
            // If after group
        } else {
            nt += repeat; // add all repeats
        }
    }
    return nt;
}

//! Set the parameters of the time selection
void TimeSelection::set( double new_start, double new_end, double new_period )
{
    start  = new_start;
    end    = new_end;
    period = new_period;
}

//! Obtain some information about the time selection
std::string TimeSelection::info()
{
    ostringstream t( "" );
    if( period == 0. ) {
        t << "never";
    } else {
        if( round( period ) == period ) {
            t << "every " << ( int ) period << " iterations";
        } else {
            t << "every " << fixed << setprecision( 1 ) << period << " iterations";
        }
        
        if( repeat > 1 ) {
            t << " (" << repeat << " repeats each";
            if( spacing > 1. ) {
                if( round( spacing ) == spacing ) {
                    t << ", spaced by " << ( int ) spacing;
                } else {
                    t << ", spaced by " << fixed << setprecision( 1 ) << spacing;
                }
            }
            t << ")";
        }
        
        if( start > 0. ) {
            if( round( start ) == start ) {
                t << " from " << ( int ) start;
            } else {
                t << " from " << fixed << setprecision( 1 ) << start;
            }
        }
        
        if( end < std::numeric_limits<double>::max() ) {
            if( round( end ) == end ) {
                t << " until " << ( int ) end;
            } else {
                t << " until " << fixed << setprecision( 1 ) << end;
            }
        }
    }
    
    return t.str();
}
