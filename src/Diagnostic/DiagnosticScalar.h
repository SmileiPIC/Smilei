#ifndef DIAGNOSTICSCALAR_H
#define DIAGNOSTICSCALAR_H

#include "Diagnostic.h"

#include "Params.h"
#include "Patch.h"
#include "SmileiMPI.h"


//! double-int structure to communicate min/max and location trough MPI 
struct val_index
{
    //! min/max 
    double val;
    //! cell location index
    int index;
};


class DiagnosticScalar : public Diagnostic {
    friend class SmileiMPI;

public :
    //! Default constructor
    DiagnosticScalar( Params &params, SmileiMPI* smpi, Patch* patch );

    //! Default destructor
    ~DiagnosticScalar() override;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    void closeFile() override;
    
    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;
    
    bool prepare( int timestep ) override;
    
    void run( Patch* patch, int timestep ) override;
    
    void write(int timestep) override;
    
    virtual bool needsRhoJs(int timestep) override;
    
    //! get a particular scalar
    double getScalar(std::string name);
    
    //! every for the standard pic timeloop output
    unsigned int print_every;
    
    //! initial energy (kinetic + EM)
    double Energy_time_zero;
    
    //! energy used for the normalization of energy balance (always uses current energy for normalization)
    double EnergyUsedForNorm;
    
    
    //! Compute the various scalars when requested
    void compute( Patch* patch, int timestep );
    
    //! Tell whether a printout is needed now
    inline bool printNow( int timestep ) {
        return (timestep % print_every == 0.);
    }
    
private :
    
    //! Calculate the length of a string when output to the file
    unsigned int calculateWidth( std::string key);
    
    //! sets a scalar name in the list of scalars (for initialization)
    int setKey( std::string key, int &currentIndex );
    //! sets a scalar name in the list of scalars (for initialization), in the case of min scalar
    int setKey_MINLOC( std::string key, int &currentIndex );
    //! sets a scalar name in the list of scalars (for initialization), in the case of max scalar
    int setKey_MAXLOC( std::string key, int &currentIndex );
    
    //! check if key is allowed
    bool allowedKey(std::string);
    
    //! write precision
    unsigned int precision;
    
    //! list of keys for scalars to be written
    std::vector<std::string> vars;
    
    //! these are lists to keep variable names, values, and allowance
    std::vector<std::string>  out_key;
    std::vector<double>       out_value;
    std::vector<unsigned int> out_width;
    std::vector<bool>         allowed;
    //! these are lists to keep variable names, values, and allowance, for min scalars
    std::vector<std::string>  out_key_MINLOC;
    std::vector<val_index>    out_value_MINLOC;
    std::vector<unsigned int> out_width_MINLOC;
    std::vector<bool>         allowed_MINLOC;
    //! these are lists to keep variable names, values, and allowance, for max scalars
    std::vector<std::string>  out_key_MAXLOC;
    std::vector<val_index>    out_value_MAXLOC;
    std::vector<unsigned int> out_width_MAXLOC;
    std::vector<bool>         allowed_MAXLOC;
    //! width of each field
    
    //! copied from params
    double cell_volume;
    
    //! this is copied from params
    double res_time;
    
    double dt;
    
    //! output stream
    std::ofstream fout;
    
    //! Variables that hold the index of several scalars in the arrays
    int index_Utot, index_Uexp, index_Ubal, index_Ubal_norm;
    int index_Uelm, index_Ukin, index_Uelm_bnd, index_Ukin_bnd;
    int index_Ukin_out_mvw, index_Ukin_inj_mvw, index_Uelm_out_mvw, index_Uelm_inj_mvw;
    std::vector<int> index_sDens, index_sNtot, index_sZavg, index_sUkin, index_fieldUelm;
    std::vector<int> index_fieldMin, index_fieldMax;
    std::vector<int> index_poy, index_poyInst;

};

#endif

