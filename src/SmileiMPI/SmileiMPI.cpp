#include "SmileiMPI.h"

#include <cmath>
#include <cstring>

#include <iostream>
#include <sstream>

#include "PicParams.h"
#include "DiagParams.h"
#include "Tools.h"

#include "ElectroMagn.h"
#include "Field.h"

#include "Species.h"

using namespace std;

SmileiMPI::SmileiMPI( int* argc, char*** argv )
{    
    int mpi_provided;

    MPI_Init_thread( argc, argv, MPI_THREAD_FUNNELED, &mpi_provided );
    if (mpi_provided == MPI_THREAD_SINGLE){
        MESSAGE("openMP not supported");
    }

    SMILEI_COMM_WORLD = MPI_COMM_WORLD;
    MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
    MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );

}

SmileiMPI::SmileiMPI( SmileiMPI *smpi )
{
    SMILEI_COMM_WORLD = smpi->SMILEI_COMM_WORLD;
    MPI_Comm_size( SMILEI_COMM_WORLD, &smilei_sz );
    MPI_Comm_rank( SMILEI_COMM_WORLD, &smilei_rk );

    oversize = smpi->oversize;
    cell_starting_global_index = smpi->cell_starting_global_index;
    min_local = smpi->min_local;
    max_local = smpi->max_local;

    n_space_global = smpi->n_space_global;

}

SmileiMPI::~SmileiMPI()
{
    int status = 0;
    MPI_Finalized( &status );
    if (!status) MPI_Finalize();

}

void SmileiMPI::bcast( InputData& input_data )
{
    DEBUG(10,"broadcast namelist");
    bcast(input_data.namelist);    

    input_data.parseStream();
    
}

void SmileiMPI::bcast( string& val )
{
    int charSize=0;
    if (isMaster()) charSize = val.size()+1;
    MPI_Bcast(&charSize, 1, MPI_INT, 0, SMILEI_COMM_WORLD);

    char tmp[charSize];
    strcpy(tmp, val.c_str());
    MPI_Bcast(&tmp, charSize, MPI_CHAR, 0, SMILEI_COMM_WORLD);

    if (!isMaster()) val=tmp;

}

void SmileiMPI::init( PicParams& params )
{
    oversize.resize(params.nDim_field, 0);
    cell_starting_global_index.resize(params.nDim_field, 0);
    min_local.resize(params.nDim_field, 0.);
    max_local.resize(params.nDim_field, 0.);
    n_space_global.resize(params.nDim_field, 0);
}


void SmileiMPI::sumRho( ElectroMagn* EMfields )
{
    sumField( EMfields->rho_ );

}

void SmileiMPI::sumRhoJ( ElectroMagn* EMfields )
{
    // sum total charge density and currents
    sumField( EMfields->rho_ );
    sumField( EMfields->Jx_ );
    sumField( EMfields->Jy_ );
    sumField( EMfields->Jz_ );
   
}
void SmileiMPI::sumRhoJs( ElectroMagn* EMfields, int ispec, bool currents )
{
   // sum density and currents for all species
   sumField( EMfields->rho_s[ispec] );
   if(currents){
       sumField( EMfields->Jx_s[ispec] );
       sumField( EMfields->Jy_s[ispec] );
       sumField( EMfields->Jz_s[ispec] );
   }
}

void SmileiMPI::exchangeE( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Ex_ );
    exchangeField( EMfields->Ey_ );
    exchangeField( EMfields->Ez_ );

}
void SmileiMPI::exchangeE( ElectroMagn* EMfields, int clrw )
{
    exchangeField_movewin( EMfields->Ex_, clrw );
    exchangeField_movewin( EMfields->Ey_, clrw );
    exchangeField_movewin( EMfields->Ez_, clrw );

}

void SmileiMPI::exchangeB( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Bx_ );
    exchangeField( EMfields->By_ );
    exchangeField( EMfields->Bz_ );

}
void SmileiMPI::exchangeB( ElectroMagn* EMfields, int clrw )
{
    exchangeField_movewin( EMfields->Bx_, clrw );
    exchangeField_movewin( EMfields->By_, clrw);
    exchangeField_movewin( EMfields->Bz_, clrw );

}

void SmileiMPI::exchangeBm( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Bx_m );
    exchangeField( EMfields->By_m );
    exchangeField( EMfields->Bz_m );

}
void SmileiMPI::exchangeBm( ElectroMagn* EMfields, int clrw )
{
    exchangeField_movewin( EMfields->Bx_m, clrw );
    exchangeField_movewin( EMfields->By_m, clrw );
    exchangeField_movewin( EMfields->Bz_m, clrw );

}

void SmileiMPI::exchangeAvg( ElectroMagn* EMfields )
{
    exchangeField( EMfields->Ex_avg );
    exchangeField( EMfields->Ey_avg );
    exchangeField( EMfields->Ez_avg );
    exchangeField( EMfields->Bx_avg );
    exchangeField( EMfields->By_avg );
    exchangeField( EMfields->Bz_avg );
}
