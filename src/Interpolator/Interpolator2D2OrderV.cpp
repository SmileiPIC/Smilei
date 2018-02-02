#include "Interpolator2D2OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D2OrderV
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D2OrderV::Interpolator2D2OrderV(Params &params, Patch *patch) : Interpolator2D(params, patch)
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    D_inv[0] = 1.0/params.cell_length[0];
    D_inv[1] = 1.0/params.cell_length[1];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd OrderV Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D2OrderV::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double* ELoc, double* BLoc)
{
}

void Interpolator2D2OrderV::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread)
{
    if ( istart[0] == iend[0] ) return; //Don't treat empty cells.

    int nparts( particles.size() );

    double *Epart[3], *Bpart[3];
    double E,E2;

    double *deltaO[2];
    deltaO[0] = &(smpi->dynamics_deltaold[ithread][0]);
    deltaO[1] = &(smpi->dynamics_deltaold[ithread][nparts]);

    for (unsigned int k=0; k<3;k++) {   
        Epart[k]= &(smpi->dynamics_Epart[ithread][k*nparts]);
        Bpart[k]= &(smpi->dynamics_Bpart[ithread][k*nparts]);
    }

    int idx[2], idxO[2];
    //Primal indices are constant over the all cell
    idx[0]  = round( particles.position(0,*istart) * D_inv[0] );
    idxO[0] = idx[0] - i_domain_begin -1 ;
    idx[1]  = round( particles.position(1,*istart) * D_inv[1] );
    idxO[1] = idx[1] - j_domain_begin -1 ;

    double **Egrid[3], **Bgrid[3];

    Egrid[0] = (static_cast<Field2D*>(EMfields->Ex_))->data_2D;
    Egrid[1] = (static_cast<Field2D*>(EMfields->Ey_))->data_2D;
    Egrid[2] = (static_cast<Field2D*>(EMfields->Ez_))->data_2D;
    Bgrid[0] = (static_cast<Field2D*>(EMfields->Bx_m))->data_2D;
    Bgrid[1] = (static_cast<Field2D*>(EMfields->By_m))->data_2D;
    Bgrid[2] = (static_cast<Field2D*>(EMfields->Bz_m))->data_2D;

    double coeff[2][2][3][32]; 
    int dual[2][32]; // Size ndim. Boolean indicating if the part has a dual indice equal to the primal one (dual=0) or if it is +1 (dual=1).

    int vecSize = 32;

    int cell_nparts( (int)iend[0]-(int)istart[0] );
    int nbVec = ( iend[0]-istart[0]+(cell_nparts-1)-((iend[0]-istart[0]-1)&(cell_nparts-1)) ) / vecSize;

    if (nbVec*vecSize != cell_nparts)
        nbVec++;

    for (int iivect=0 ; iivect<nbVec; iivect++ ){
        int ivect = vecSize*iivect;

        int np_computed(0);
        if (cell_nparts > vecSize ) {
            np_computed = vecSize;
            cell_nparts -= vecSize;
        }       
        else
            np_computed = cell_nparts;


        #pragma omp simd
        for (int ipart=0 ; ipart<np_computed; ipart++ ){

            double delta0, delta;
            double delta2;

            for (int i=0;i<2;i++) { // for X/Y
                delta0 = particles.position(i,ipart+ivect+istart[0])*D_inv[i];
                dual [i][ipart] = ( delta0 - (double)idx[i] >=0. );

                for (int j=0;j<2;j++) { // for dual

                    delta   = delta0 - (double)idx[i] + (double)j*(0.5-dual[i][ipart]);
                    delta2  = delta*delta;

                    coeff[i][j][0][ipart]    =  0.5 * (delta2-delta+0.25);
                    coeff[i][j][1][ipart]    =  (0.75 - delta2);
                    coeff[i][j][2][ipart]    =  0.5 * (delta2+delta+0.25);
    
                    if (j==0) deltaO[i][ipart+ivect+istart[0]] = delta;
    
                }
            }
        }

        #pragma omp simd
        for (int ipart=0 ; ipart<np_computed; ipart++ ){
           //Ex(dual, primal)
            E  =   coeff[1][0][0][ipart]* Egrid[0][idxO[0]][idxO[1]];
            E  +=  coeff[1][0][1][ipart]* Egrid[0][idxO[0]][idxO[1]+1];
            E  +=  coeff[1][0][2][ipart]* Egrid[0][idxO[0]][idxO[1]+2];
            Epart[0][ipart+ivect+istart[0]] =   (1-dual[0][ipart])*coeff[0][1][0][ipart]* E ;   
            E  =   coeff[1][0][0][ipart]* Egrid[0][idxO[0]+1][idxO[1]];
            E  +=  coeff[1][0][1][ipart]* Egrid[0][idxO[0]+1][idxO[1]+1];
            E  +=  coeff[1][0][2][ipart]* Egrid[0][idxO[0]+1][idxO[1]+2];
            Epart[0][ipart+ivect+istart[0]] +=  ( (1-dual[0][ipart])*coeff[0][1][1][ipart]+ dual[0][ipart]*coeff[0][1][0][ipart]) * E ;   
            E  =   coeff[1][0][0][ipart]* Egrid[0][idxO[0]+2][idxO[1]];
            E  +=  coeff[1][0][1][ipart]* Egrid[0][idxO[0]+2][idxO[1]+1];
            E  +=  coeff[1][0][2][ipart]* Egrid[0][idxO[0]+2][idxO[1]+2];
            Epart[0][ipart+ivect+istart[0]] +=  ( (1-dual[0][ipart])*coeff[0][1][2][ipart]+ dual[0][ipart]*coeff[0][1][1][ipart]) * E ;   
            E  =   coeff[1][0][0][ipart]* Egrid[0][idxO[0]+3][idxO[1]];
            E  +=  coeff[1][0][1][ipart]* Egrid[0][idxO[0]+3][idxO[1]+1];
            E  +=  coeff[1][0][2][ipart]* Egrid[0][idxO[0]+3][idxO[1]+2];
            Epart[0][ipart+ivect+istart[0]] +=   dual[0][ipart]*coeff[0][1][2][ipart]* E ;   


            //Ey(primal, dual)
            E =   coeff[1][1][0][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]][idxO[1]  ]  + dual[1][ipart]*Egrid[1][idxO[0]][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]][idxO[1]+1]  + dual[1][ipart]*Egrid[1][idxO[0]][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]][idxO[1]+2]  + dual[1][ipart]*Egrid[1][idxO[0]][idxO[1]+3] ) ;
            Epart[1][ipart+ivect+istart[0]] =   coeff[0][0][0][ipart]* E;
            E =   coeff[1][1][0][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]+1][idxO[1]  ]  + dual[1][ipart]*Egrid[1][idxO[0]+1][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]+1][idxO[1]+1]  + dual[1][ipart]*Egrid[1][idxO[0]+1][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]+1][idxO[1]+2]  + dual[1][ipart]*Egrid[1][idxO[0]+1][idxO[1]+3] ) ;
            Epart[1][ipart+ivect+istart[0]] +=   coeff[0][0][1][ipart]* E;
            E =   coeff[1][1][0][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]+2][idxO[1]  ]  + dual[1][ipart]*Egrid[1][idxO[0]+2][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]+2][idxO[1]+1]  + dual[1][ipart]*Egrid[1][idxO[0]+2][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart])*Egrid[1][idxO[0]+2][idxO[1]+2]  + dual[1][ipart]*Egrid[1][idxO[0]+2][idxO[1]+3] ) ;
            Epart[1][ipart+ivect+istart[0]] +=   coeff[0][0][2][ipart]* E;

            //Ez(primal, primal)
            E =  coeff[1][0][0][ipart]* Egrid[2][idxO[0]][idxO[1]];
            E +=  coeff[1][0][1][ipart]* Egrid[2][idxO[0]][idxO[1]+1];
            E +=  coeff[1][0][2][ipart]* Egrid[2][idxO[0]][idxO[1]+2];
            Epart[2][ipart+ivect+istart[0]] =   coeff[0][0][0][ipart]* E;
            E =  coeff[1][0][0][ipart]* Egrid[2][idxO[0]+1][idxO[1]];
            E +=  coeff[1][0][1][ipart]* Egrid[2][idxO[0]+1][idxO[1]+1];
            E +=  coeff[1][0][2][ipart]* Egrid[2][idxO[0]+1][idxO[1]+2];
            Epart[2][ipart+ivect+istart[0]] +=   coeff[0][0][1][ipart]* E;
            E =  coeff[1][0][0][ipart]* Egrid[2][idxO[0]+2][idxO[1]];
            E +=  coeff[1][0][1][ipart]* Egrid[2][idxO[0]+2][idxO[1]+1];
            E +=  coeff[1][0][2][ipart]* Egrid[2][idxO[0]+2][idxO[1]+2];
            Epart[2][ipart+ivect+istart[0]] +=   coeff[0][0][2][ipart]* E;

            //Bx(primal, dual  )
            E =   coeff[1][1][0][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]][idxO[1]  ]  + dual[1][ipart]*Bgrid[0][idxO[0]][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]][idxO[1]+1]  + dual[1][ipart]*Bgrid[0][idxO[0]][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]][idxO[1]+2]  + dual[1][ipart]*Bgrid[0][idxO[0]][idxO[1]+3] ) ;
            Bpart[0][ipart+ivect+istart[0]] =   coeff[0][0][0][ipart]* E;
            E =   coeff[1][1][0][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]+1][idxO[1]  ]  + dual[1][ipart]*Bgrid[0][idxO[0]+1][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]+1][idxO[1]+1]  + dual[1][ipart]*Bgrid[0][idxO[0]+1][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]+1][idxO[1]+2]  + dual[1][ipart]*Bgrid[0][idxO[0]+1][idxO[1]+3] ) ;
            Bpart[0][ipart+ivect+istart[0]] +=   coeff[0][0][1][ipart]* E ;
            E =   coeff[1][1][0][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]+2][idxO[1]  ]  + dual[1][ipart]*Bgrid[0][idxO[0]+2][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]+2][idxO[1]+1]  + dual[1][ipart]*Bgrid[0][idxO[0]+2][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart])*Bgrid[0][idxO[0]+2][idxO[1]+2]  + dual[1][ipart]*Bgrid[0][idxO[0]+2][idxO[1]+3] ) ;
            Bpart[0][ipart+ivect+istart[0]] +=   coeff[0][0][2][ipart]* E ;

            //By(dual, primal )
            E =  coeff[1][0][0][ipart] * Bgrid[1][idxO[0]][idxO[1]+0];
            E +=  coeff[1][0][1][ipart]* Bgrid[1][idxO[0]][idxO[1]+1];
            E +=  coeff[1][0][2][ipart]* Bgrid[1][idxO[0]][idxO[1]+2];
            Bpart[1][ipart+ivect+istart[0]] =   (1-dual[0][ipart])*coeff[0][1][0][ipart]* E;
            E =  coeff[1][0][0][ipart] * Bgrid[1][idxO[0]+1][idxO[1]+0];
            E +=  coeff[1][0][1][ipart]* Bgrid[1][idxO[0]+1][idxO[1]+1];
            E +=  coeff[1][0][2][ipart]* Bgrid[1][idxO[0]+1][idxO[1]+2];
            Bpart[1][ipart+ivect+istart[0]] +=  ( (1-dual[0][ipart])*coeff[0][1][1][ipart]+ dual[0][ipart]*coeff[0][1][0][ipart]) * E ;   
            E =  coeff[1][0][0][ipart] * Bgrid[1][idxO[0]+2][idxO[1]+0];
            E +=  coeff[1][0][1][ipart]* Bgrid[1][idxO[0]+2][idxO[1]+1];
            E +=  coeff[1][0][2][ipart]* Bgrid[1][idxO[0]+2][idxO[1]+2];
            Bpart[1][ipart+ivect+istart[0]] +=  ( (1-dual[0][ipart])*coeff[0][1][2][ipart]+ dual[0][ipart]*coeff[0][1][1][ipart]) * E ;   
            E =  coeff[1][0][0][ipart] * Bgrid[1][idxO[0]+3][idxO[1]+0];
            E +=  coeff[1][0][1][ipart]* Bgrid[1][idxO[0]+3][idxO[1]+1];
            E +=  coeff[1][0][2][ipart]* Bgrid[1][idxO[0]+3][idxO[1]+2];
            Bpart[1][ipart+ivect+istart[0]] +=  dual[0][ipart]*coeff[0][1][2][ipart] * E ;   

            //Bz(dual, dual )
            E =  coeff[1][1][0][ipart] * ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]][idxO[1]  ] + dual[1][ipart]*Bgrid[2][idxO[0]][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]][idxO[1]+1] + dual[1][ipart]*Bgrid[2][idxO[0]][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]][idxO[1]+2] + dual[1][ipart]*Bgrid[2][idxO[0]][idxO[1]+3] ) ;
            Bpart[2][ipart+ivect+istart[0]] =   (1-dual[0][ipart])*coeff[0][1][0][ipart]* E;
            E =  coeff[1][1][0][ipart] * ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+1][idxO[1]+0] + dual[1][ipart]*Bgrid[2][idxO[0]+1][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+1][idxO[1]+1] + dual[1][ipart]*Bgrid[2][idxO[0]+1][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+1][idxO[1]+2] + dual[1][ipart]*Bgrid[2][idxO[0]+1][idxO[1]+3] ) ;
            Bpart[2][ipart+ivect+istart[0]] +=  ( (1-dual[0][ipart])*coeff[0][1][1][ipart]+ dual[0][ipart]*coeff[0][1][0][ipart]) * E ;   
            E =  coeff[1][1][0][ipart] * ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+2][idxO[1]  ] + dual[1][ipart]*Bgrid[2][idxO[0]+2][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+2][idxO[1]+1] + dual[1][ipart]*Bgrid[2][idxO[0]+2][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+2][idxO[1]+2] + dual[1][ipart]*Bgrid[2][idxO[0]+2][idxO[1]+3] ) ;
            Bpart[2][ipart+ivect+istart[0]] +=  ( (1-dual[0][ipart])*coeff[0][1][2][ipart]+ dual[0][ipart]*coeff[0][1][1][ipart]) * E ;   
            E =  coeff[1][1][0][ipart] * ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+3][idxO[1]  ] + dual[1][ipart]*Bgrid[2][idxO[0]+3][idxO[1]+1] ) ;
            E +=  coeff[1][1][1][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+3][idxO[1]+1] + dual[1][ipart]*Bgrid[2][idxO[0]+3][idxO[1]+2] ) ;
            E +=  coeff[1][1][2][ipart]* ( (1-dual[1][ipart]) * Bgrid[2][idxO[0]+3][idxO[1]+2] + dual[1][ipart]*Bgrid[2][idxO[0]+3][idxO[1]+3] ) ;
            Bpart[2][ipart+ivect+istart[0]] +=  dual[0][ipart]*coeff[0][1][2][ipart] * E ;   
        }
    }

} // END Interpolator2D2OrderV

void Interpolator2D2OrderV::operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc)
{
}
