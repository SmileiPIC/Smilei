#include "ElectroMagnBC2D_Trans_Damping.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "ElectroMagn2D.h"
#include "Field2D.h"
#include "Tools.h"

using namespace std;

ElectroMagnBC2D_Trans_Damping::ElectroMagnBC2D_Trans_Damping( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC2D( params, patch, _min_max )
{
    // number of dumping layers
    ny_l = 50;// To be read in file.in
    cdamp = 1.l;// To be read in file.in
    
    coeff.resize( ny_l );
    coeff[0] = 0.;
    
    for( unsigned int j=1 ; j<ny_l ; j++ ) {
        coeff[j] = 1.-cdamp*pow( ( ny_l-( double )j )/ny_l, 2 );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_Trans_Damping::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    // Static cast of the fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
//    Field2D* Ey2D = static_cast<Field2D*>(EMfields->Ey_);
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_ );
//    Field2D* By2D = static_cast<Field2D*>(EMfields->By_);
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_ );
    
    
    //   BC : Bx(i=0...nx_p, 0) & Bx(i=0...nx_p, ny_d-1)
    if( patch->isYmin() ) {
        // for Bx^(p,d)
        for( unsigned int i=0 ; i<nx_p ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Bx2D )( i, j ) *= coeff[j];
            }
            
        // for Bz^(d,d)
        for( unsigned int i=0 ; i<nx_d ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Bz2D )( i, j ) *= coeff[j];
            }
            
        // for By^(d,p)
        //for (unsigned int i=0 ; i<nx_d ; i++)
        //    for (unsigned int j=0 ; j<ny_l ; j++)
        //        (*By2D)(i,j) *= coeff[j];
        
        
        // for Ex^(d,p)
        for( unsigned int i=0 ; i<nx_d ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Ex2D )( i, j ) *= coeff[j];
            }
            
        // for Ez^(p,p)
        for( unsigned int i=0 ; i<nx_p ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Ez2D )( i, j ) *= coeff[j];
            }
            
        // for Ey^(p,d)
        //for (unsigned int i=0 ; i<nx_p ; i++)
        //    for (unsigned int j=0 ; j<ny_l ; j++)
        //        (*Ey2D)(i,j) *= coeff[j];
        
    }
    
    //   BC : Bz(i=0...nx_d-1, 0) & Bz(i=0...nx_d-1, ny_d-1)
    if( patch->isYmax() ) {
        // for Bx^(p,d)
        for( unsigned int i=0 ; i<nx_p ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Bx2D )( i, ny_d-1-j ) *= coeff[j];
            }
            
        // for Bz^(d,d)
        for( unsigned int i=0 ; i<nx_d ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Bz2D )( i, ny_d-1-j ) *= coeff[j];
            }
            
        // for By^(d,p)
        //for (unsigned int i=0 ; i<nx_d ; i++)
        //    for (unsigned int j=0 ; j<ny_l ; j++)
        //        (*By2D)(i,ny_p-1-j) *= coeff[j];
        
        // for Ex^(d,p)
        for( unsigned int i=0 ; i<nx_d ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Ex2D )( i, ny_p-1-j ) *= coeff[j];
            }
            
        // for Ez^(p,p)
        for( unsigned int i=0 ; i<nx_p ; i++ )
            for( unsigned int j=0 ; j<ny_l ; j++ ) {
                ( *Ez2D )( i, ny_p-1-j ) *= coeff[j];
            }
            
        // for Ey^(p,d)
        //for (unsigned int i=0 ; i<nx_p ; i++)
        //    for (unsigned int j=0 ; j<ny_l ; j++)
        //        (*Ey2D)(i,ny_d-1-j) *= coeff[j];
    }
    
    
}// END apply

