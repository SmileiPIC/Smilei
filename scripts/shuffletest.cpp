/*

This program tests the speed and accuracy of the shuffling method used for binary processes.

*/


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>

#include "../src/Tools/RandomShuffle.h"

using namespace std;

static double test_shuffle(int n, Random &R)
{
    int repeats = 3e6/n;
    
    vector<int> all( n*n, 0 );
    vector<int> result( n, 0 );
    
    const clock_t initial_time = clock();
    for( int r=0; r<repeats; r++) {
        RandomShuffle shuffler( R, n );
        for( int i=0; i<n; i++ ) {
            result[i] = shuffler.next();
            // cout << result[i] << " ";
            all[i*n + result[i]] ++;
        }
        // cout << endl;
    }
    clock_t time = clock()-initial_time;
    
    // cout << "n = " << n << endl;
    // cout << "Time per shuffle = " << ((double)(time))/repeats << endl;
    // cout << "Total time = " << time << endl;
    
    double max_diff = 0;
    for( int i=0; i<n; i++) {
        for( int j=0; j<n; j++) {
            double diff = fabs( ((double) all[i*n+j])/repeats - 1./((double)n) );
            if( diff > max_diff ) max_diff = diff;
        }
    }
    cout << "n = " << n << "   repeats = " << repeats << "     max diff = " << max_diff << endl;
    
    return ((double)(time))/repeats;
}


int main( int argc, char *argv[] )
{
    Random R( clock() );
    
    size_t nmax = 2000;
    vector<double> times( nmax, 0. );
    
    for(int n=2; n<nmax; n+=100) {
        times[n] = test_shuffle(n, R);
    }
    
    cout << "plot([";
    for(int n=2; n<nmax; n+=100) {
        cout << n << ",";
    }
    cout << "],[";
    for(int n=2; n<nmax; n+=100) {
        cout << times[n] << ",";
    }
    cout << "])" << endl;
    
    return 0;
}