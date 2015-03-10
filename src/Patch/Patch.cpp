#include "Patch.h"

void Patch::dtoxy(int npatches, int *x, int *y, int ipatch){

int rx, ry, s, t;

//convert ipatch (number of the patch along the Hilbert curve) to (x,y), the coordinate of the patch in the Patch Grid.
    t = ipatch;
    *x = 0;
    *y = 0;
    for (s=1; s<npatches; s*=2) {
        rx = 1 & (t/2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

//rotate/flip a quadrant appropriately
void Patch::rot(int n, int *x, int *y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            *x = n-1 - *x;
            *y = n-1 - *y;
        }
 
        //Swap x and y
        int t  = *x;
        *x = *y;
        *y = t;
    }
}

