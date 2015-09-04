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

   unsigned int m0 = 5;
   unsigned int m1 = 3;
   unsigned int X,Y;
   int nmpi = 15;
   ofstream myfile;
   myfile.open ("data.txt");

   


   cout << (1<<m0) << "patches selon x, " << (1<<m1) << " patches selon y." << endl;

   int step, total_patch;
   total_patch = (1<<m0)*(1<<m1);
   //step = total_patch/nmpi + 1;
   step = total_patch/nmpi ;
   for (unsigned int ipatch = 0; ipatch < total_patch; ipatch++){
       generalhilbertindexinv(m0, m1, &X, &Y, ipatch);
       myfile << ipatch << " " << X << " " << Y << endl;

   }
   int first_patch = 0;
   int last_patch=-1;


   for (unsigned int impi=0; impi < nmpi-1; impi++){
       first_patch = last_patch+1;
       last_patch = first_patch-1+step;
       cout << "rank "<<impi<<" first "<<first_patch<<" last "<<last_patch<< " npatches= " << last_patch+1-first_patch << endl;
   } 
       unsigned int impi = nmpi-1;
       first_patch = last_patch+1;
       last_patch = total_patch-1;
       cout << "rank "<<impi<<" first "<<first_patch<<" last "<<last_patch<< " npatches= " << last_patch+1-first_patch << endl;
    
   myfile.close();
   return 0;
    
}//END MAIN



