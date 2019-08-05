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
   unsigned int m0 = 10;
   unsigned int m1 = 6;
   unsigned int m2 = 6;
/************************************/


   unsigned int X,Y,Z;
   ofstream myfile;
   myfile.open ("data.txt");

   cout << (1<<m0) << " patches selon x, " << (1<<m1) << " patches selon y."<< (1<<m2) << " patches selon z."  << endl;

   int total_patch;
   total_patch = (1<<m0)*(1<<m1)*(1<<m2);
   for (unsigned int ipatch = 0; ipatch < total_patch; ipatch++){
       generalhilbertindexinv(m0, m1, m2, &X, &Y, &Z, ipatch);
       myfile << ipatch << " " << X << " " << Y <<  " " << Z << endl;

   }
    
   myfile.close();
   return 0;
    
}//END MAIN



