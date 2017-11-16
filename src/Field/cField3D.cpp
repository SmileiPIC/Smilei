#include "cField3D.h"

#include <iostream>
#include <vector>
#include <cstring>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for cField3D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
cField3D::cField3D() : Field()
{
    data_1D=NULL;
}

// with the dimensions as input argument
cField3D::cField3D(vector<unsigned int> dims) : Field(dims)
{
    data_1D=NULL;
    allocateDims(dims);
}

// with the dimensions and output (dump) file name as input argument
cField3D::cField3D(vector<unsigned int> dims, string name) : Field(dims, name)
{
    data_1D=NULL;
    allocateDims(dims);
}

// with the dimensions as input argument
cField3D::cField3D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal) : Field(dims, mainDim, isPrimal)
{
    data_1D=NULL;
    allocateDims(dims, mainDim, isPrimal);
}

// with the dimensions and output (dump) file name as input argument
cField3D::cField3D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name) : Field(dims, mainDim, isPrimal, name)
{
    data_1D=NULL;
    allocateDims(dims, mainDim, isPrimal);
}

// without allocating
cField3D::cField3D(string name, vector<unsigned int> dims) : Field(dims, name)
{
    data_1D=NULL;
    dims_=dims;
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for cField3D
// ---------------------------------------------------------------------------------------------------------------------
cField3D::~cField3D()
{

    if (data_1D!=NULL) {
        delete [] data_1D;
        delete [] data_3D;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField3D
// ---------------------------------------------------------------------------------------------------------------------
void cField3D::allocateDims()
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    if (dims_.size()!=3) ERROR("Alloc error must be 3 : " << dims_.size());
    if (data_1D!=NULL) delete [] data_1D;

    isDual_.resize( dims_.size(), 0 );

    data_1D = new complex<double>[dims_[0]*dims_[1]*dims_[2]];
    //! \todo{check row major order!!! (JD)}

    data_3D= new complex<double>**[dims_[0]];
    for (unsigned int i=0; i<dims_[0]; i++) {
        data_3D[i]= new complex<double>*[dims_[1]];
        for (unsigned int j=0; j<dims_[1]; j++)
        {
            data_3D[i][j] = data_1D + i*dims_[1]*dims_[2] + j*dims_[2];
            for (unsigned int k=0; k<dims_[2]; k++) data_3D[i][j][k] = 0.0;
        }
    }
    globalDims_ = dims_[0]*dims_[1]*dims_[2];

}

void cField3D::deallocateDims()
{
    delete [] data_1D;
    data_1D = NULL;
    delete [] data_3D;
    data_3D = NULL;
        
}

void cField3D::allocateDims(unsigned int dims1, unsigned int dims2, unsigned int dims3)
{
    vector<unsigned int> dims(3);
    dims[0]=dims1;
    dims[1]=dims2;
    dims[2]=dims3;
    allocateDims(dims);
}

// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField3D
// ---------------------------------------------------------------------------------------------------------------------
void cField3D::allocateDims(unsigned int mainDim, bool isPrimal )
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    if (dims_.size()!=3) ERROR("Alloc error must be 3 : " << dims_.size());
    if (data_1D) delete [] data_1D;
    
    // isPrimal define if mainDim is Primal or Dual
    isDual_.resize( dims_.size(), 0 );
    for ( unsigned int j=0 ; j<dims_.size() ; j++ ) {
        if ( (j==mainDim) && (!isPrimal) )
            isDual_[j] = 1;
        else if ( (j!=mainDim) && (isPrimal) )
            isDual_[j] = 1;
    }
    
    for ( unsigned int j=0 ; j<dims_.size() ; j++ )
        dims_[j] += isDual_[j];
    
    data_1D = new complex<double>[dims_[0]*dims_[1]*dims_[2]];
    //! \todo{check row major order!!! (JD)}
    
    data_3D= new complex<double>**[dims_[0]];
    for (unsigned int i=0; i<dims_[0]; i++)  {
                data_3D[i]= new complex<double>*[dims_[1]];
        for (unsigned int j=0; j<dims_[1]; j++)
        {
            data_3D[i][j] = data_1D + i*dims_[1]*dims_[2] + j*dims_[2];
            for (unsigned int k=0; k<dims_[2]; k++) data_3D[i][j][k] = 0.0;
        }

    }
    
    globalDims_ = dims_[0]*dims_[1]*dims_[2];
    
}



// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void cField3D::shift_x(unsigned int delta)
{
    memmove( &(data_3D[0][0]), &(data_3D[delta][0]), (dims_[1]*dims_[0]-delta*dims_[1])*sizeof(complex<double>) );
    memset( &(data_3D[dims_[0]-delta][0]), 0, delta*dims_[1]*sizeof(complex<double>));
    
}

double cField3D::norm2(unsigned int istart[3][2], unsigned int bufsize[3][2]) {
    double nrj(0.);
    
    int idxlocalstart[3];
    int idxlocalend[3];
    for ( int i=0 ; i<3 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    
    for ( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        for ( int j=idxlocalstart[1] ; j<idxlocalend[1] ; j++ ) {
            for ( int k=idxlocalstart[2] ; k<idxlocalend[2] ; k++ ) {
                nrj += (data_3D[i][j][k]).real()*(data_3D[i][j][k]).real()+ (data_3D[i][j][k]).imag()*(data_3D[i][j][k]).imag();
            }
        }
    }
   
    return nrj;
}
