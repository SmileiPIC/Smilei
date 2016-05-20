
#ifndef DIAGNOSTICFIELDS2D_H
#define DIAGNOSTICFIELDS2D_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFields2D : public DiagnosticFields {
public:
    DiagnosticFields2D( Params &params, SmileiMPI* smpi, Patch* patch, int );
    DiagnosticFields2D( DiagnosticFields* diag, Params &params, Patch* patch );
    //! Destructor for DiagnosticFields
    ~DiagnosticFields2D() ;

    //! Build memory and file space for // HDF5 write/read
    void createPattern( Params& params, Patch* patch );
    void updatePattern( Params& params, Patch* patch )override;

    //! Basic write current field in specified group of the global file
    void writeField( Field* field, hid_t group_id ) override;

    //! Basic write field on its own file (debug)
    void write( Field* field );

private:
    //! memory space for // HDF5 write/read
    //! Size = 2 x 2 : 0 if prim, 1 if dual per direction
    hid_t memspace_ [2][2];
    //! file space for // HDF5 write/read
    //! Size = 2 x 2 : 0 if prim, 1 if dual per direction
    hid_t filespace_[2][2];

    //! \todo Define chunk size of output for interpolated output
    //hsize_t chunk_dims[2];

};

#endif
