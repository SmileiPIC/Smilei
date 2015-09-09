/**This scripts generates hindex and coordinates of patches just like in Smilei**/

#include <cstdlib>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Hilbert_functions.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
//                                                   MAIN CODE
// ---------------------------------------------------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    std::cout.setf( std::ios::fixed, std:: ios::floatfield ); // floatfield set to fixed

/******* USER DEFINED PARAMETERS*****/
   unsigned int m0 = 5;
   unsigned int m1 = 3;
/************************************/


   unsigned int X,Y;
   ofstream myfile;
   myfile.open ("data.txt");

   cout << (1<<m0) << "patches selon x, " << (1<<m1) << " patches selon y." << endl;

   int total_patch;
   total_patch = (1<<m0)*(1<<m1);
   for (unsigned int ipatch = 0; ipatch < total_patch; ipatch++){
       generalhilbertindexinv(m0, m1, &X, &Y, ipatch);
       myfile << ipatch << " " << X << " " << Y << endl;

   }
    
   myfile.close();
   return 0;
    
}//END MAIN



