
#ifndef DIAGNOSTICFIELDSAM_H
#define DIAGNOSTICFIELDSAM_H

#include <string>
#include <vector>

#include "DiagnosticFields.h"
#include "Tools.h"

class DiagnosticFieldsAM : public DiagnosticFields
{
public:
    DiagnosticFieldsAM( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int, OpenPMDparams & );
    ~DiagnosticFieldsAM();
    
    void setFileSplitting( SmileiMPI *smpi, VectorPatch &vecPatches ) override;
    
    //! Copy patch field to current "data" buffer
    void getField( Patch *patch, unsigned int ) override;
    template<typename T, typename F>  void getField( Patch *patch, unsigned int, F& out_data );
    
    H5Write writeField( H5Write*, std::string ) override;
    template<typename F> H5Write writeField( H5Write*, std::string, F& linearized_data );

private:
    std::vector<unsigned int> buffer_skip_x, buffer_skip_y;
    
    std::vector<std::complex<double> > idata;
    bool is_complex_;
};

#endif
