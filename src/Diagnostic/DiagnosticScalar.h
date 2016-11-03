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


//! Class for containing info on one particular scalar data
class Scalar {
public:
    Scalar(std::string name, std::string secondname, unsigned int width, bool allowed, double * value, int * loc, double reset_value):
        name(name), secondname(secondname), width(width), allowed(allowed), value(value), loc(loc), reset_value(reset_value)
        {};
    ~Scalar() {};
    virtual inline double get() { return *value; };
    virtual inline int getloc() { return loc?(*loc):0; };
    virtual inline void reset() { *value = reset_value; if(loc) *loc = -1; };
    std::string name, secondname;
    unsigned int width;
    bool allowed;
    double * value, reset_value;
    int * loc;
};

//! Class for the diagnostic of scalars
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
        return (timestep % print_every == 0);
    }
    
private :
    
    //! Calculate the length of a string when output to the file
    unsigned int calculateWidth( std::string key);
    
    //! sets a scalar in the list of scalars (for initialization)
    Scalar* newScalar_SUM( std::string name );
    //! sets a scalar in the list of scalars (for initialization), in the case of min scalar
    Scalar* newScalar_MINLOC( std::string name );
    //! sets a scalar in the list of scalars (for initialization), in the case of max scalar
    Scalar* newScalar_MAXLOC( std::string name );
    
    //! check if key is allowed
    bool allowedKey(std::string);
    
    //! write precision
    unsigned int precision;
    
    //! list of keys for scalars to be written
    std::vector<std::string> vars;
    
    //! The list of scalars data
    std::vector<Scalar*> allScalars;
    //! List of scalar values to be summed by MPI
    std::vector<double> values_SUM;
    //! List of scalar values to be MINLOCed by MPI
    std::vector<val_index> values_MINLOC;
    //! List of scalar values to be MAXLOCed by MPI
    std::vector<val_index> values_MAXLOC;
    
    //! Volume of a cell (copied from params)
    double cell_volume;
    
    //! Time resolution (copied from params)
    double res_time;
    
    //! Timestep (copied from params)
    double dt;
    
    //! output stream
    std::ofstream fout;
    
    //! Pointers to the various scalars
    Scalar *Utot, *Uexp, *Ubal, *Ubal_norm;
    Scalar *Uelm, *Ukin, *Uelm_bnd, *Ukin_bnd;
    Scalar *Ukin_out_mvw, *Ukin_inj_mvw, *Uelm_out_mvw, *Uelm_inj_mvw;
    std::vector<Scalar *> sDens, sNtot, sZavg, sUkin, fieldUelm;
    std::vector<Scalar *> fieldMin, fieldMax;
    std::vector<Scalar *> poy, poyInst;
    
    //! Booleans to tell which scalars should be computed or not
    bool necessary_Ubal_norm, necessary_Ubal, necessary_Utot, necessary_Uexp;
    bool necessary_Ukin, necessary_Ukin_BC;
    bool necessary_Uelm, necessary_Uelm_BC;
    bool necessary_fieldMinMax_any;
    std::vector<bool> necessary_species, necessary_fieldUelm, necessary_fieldMinMax, necessary_poy;

};

#endif

