
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
    
    void writeField( hid_t, int ) override;
    template<typename F> void writeField( hid_t dset_id, int itime, F& linearized_data, F& read_data, F& final_data );

private:

    unsigned int rewrite_npatch, rewrite_xmin, rewrite_ymin, rewrite_npatchx, rewrite_npatchy;
    std::vector<std::vector<unsigned int> > rewrite_patch;
    unsigned int rewrite_size[2], rewrite_start_in_file[2];
    
    std::vector<std::complex<double>> idata_reread, idata_rewrite, idata;
    
    int factor_;

};

#endif
