#ifndef DIAGNOSTICSCALAR_H
#define DIAGNOSTICSCALAR_H

#include <fstream>

#include "Diagnostic.h"

class Patch;
class Params;
class SmileiMPI;


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
    Scalar(std::string name, std::string secondname, unsigned int width, bool allowed):
        name_(name), secondname_(secondname), width_(width), allowed_(allowed)
        {};
    virtual ~Scalar() {};
    virtual inline operator double() const { return 0.; }
    virtual inline void reset() { };
    std::string name_, secondname_;
    unsigned int width_;
    bool allowed_;
};

//! Child class of `Scalar` specific for scalars which only have one value
class Scalar_value : public Scalar {
public:
    Scalar_value(std::string name, unsigned int width, bool allowed, std::vector<double>* values):
        Scalar(name, "", width, allowed), values_(values), index(values->size())
        {};
    ~Scalar_value() {};
    inline Scalar_value& operator= (double v) {
        (*values_)[index] = v;
        return *this;
    };
    inline Scalar_value& operator+= (double v) {
        #pragma omp atomic
        (*values_)[index]+=v;
        return *this;
    };
    inline operator double() const override { return (*values_)[index]; } ;
    inline void reset() override { (*values_)[index]=0.; };
    std::vector<double>* values_;
    unsigned int index;
};

//! Child class of `Scalar` specific for scalars which have a value and a location
class Scalar_value_location : public Scalar {
public:
    Scalar_value_location(std::string name, std::string secondname, unsigned int width, bool allowed, std::vector<val_index>* values, double reset_value):
        Scalar(name, secondname, width, allowed), values_(values), index(values->size()), reset_value_(reset_value)
        {};
    ~Scalar_value_location() {};
    inline Scalar_value_location& operator= (val_index v) {
        (*values_)[index] = v;
        return *this;
    };
    inline operator double() const override { return (*values_)[index].val; }
    inline operator int() const { return (*values_)[index].index; }
    inline void reset() override { (*values_)[index].val=reset_value_; (*values_)[index].index=-1; };
    std::vector<val_index>* values_;
    unsigned int index;
    double reset_value_;
};

//! Class for the diagnostic of scalars
class DiagnosticScalar : public Diagnostic {
    friend class SmileiMPI;
    friend class Checkpoint;
public :
    //! Default constructor
    DiagnosticScalar( Params &params, SmileiMPI* smpi, Patch* patch );
    
    //! Default destructor
    ~DiagnosticScalar() override;
    
    void openFile( Params& params, SmileiMPI* smpi, bool newfile ) override;
    
    void closeFile() override;
    
    void init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches) override;
    
    bool prepare( int timestep ) override;
    
    void run( Patch* patch, int timestep, SimWindow* simWindow ) override;
    
    void write(int timestep, SmileiMPI* smpi) override;
    
    virtual bool needsRhoJs(int timestep) override;
    
    //! get a particular scalar
    double getScalar(std::string name);
    
    //! initial energy (kinetic + EM)
    double Energy_time_zero;
    
    //! energy used for the normalization of energy balance (always uses current energy for normalization)
    double EnergyUsedForNorm;
    
    //! Compute the various scalars when requested
    void compute( Patch* patch, int timestep );
    
    //! Latest timestep dumped
    int latest_timestep;
    
    //! Get memory footprint of current diagnostic
    int getMemFootPrint() override {
        return 0;
    }
    
    //! Get disk footprint of current diagnostic
    uint64_t getDiskFootPrint(int istart, int istop, Patch* patch) override;
    
private :

    //! Calculate the length of a string when output to the file
    unsigned int calculateWidth( std::string key);
    
    //! sets a scalar in the list of scalars (for initialization)
    Scalar_value* newScalar_SUM( std::string name );
    //! sets a scalar in the list of scalars (for initialization), in the case of min scalar
    Scalar_value_location* newScalar_MINLOC( std::string name );
    //! sets a scalar in the list of scalars (for initialization), in the case of max scalar
    Scalar_value_location* newScalar_MAXLOC( std::string name );
    
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
    
    //! Number of cells in a patch (copied from params)
    std::vector<unsigned int> n_space;
    
    //! Overall number of cells (copied from params)
    std::vector<unsigned int> n_space_global;
    
    //! output stream
    std::ofstream fout;
    
    //! Pointers to the various scalars
    Scalar_value *Utot, *Uexp, *Ubal, *Ubal_norm;
    Scalar_value *Uelm, *Ukin, *Uelm_bnd, *Ukin_bnd;
    Scalar_value *Ukin_out_mvw, *Ukin_inj_mvw, *Uelm_out_mvw, *Uelm_inj_mvw;
    // For the radiated energy
    Scalar_value *Urad;
    // Energy of the pairs created via the multiphoton Breit-Wheeler process
    Scalar_value *UmBWpairs;
    
    std::vector<Scalar_value *> sDens, sNtot, sZavg, sUkin, fieldUelm;
    // For the radiated energy per species
    std::vector<Scalar_value *> sUrad;
    std::vector<Scalar_value_location *> fieldMin, fieldMax;
    std::vector<Scalar_value *> poy, poyInst;
    
    //! Booleans to tell which scalars should be computed or not
    bool necessary_Ubal_norm, necessary_Ubal, necessary_Utot, necessary_Uexp;
    bool necessary_Ukin, necessary_Ukin_BC;
    bool necessary_Uelm, necessary_Uelm_BC;
    // For the radiated energy
    bool necessary_Urad;
    // For the pair generation via the multiphoton Breit-Wheeler
    bool necessary_UmBWpairs;
    bool necessary_fieldMinMax_any;
    std::vector<bool> necessary_species, necessary_fieldUelm, necessary_fieldMinMax, necessary_poy;
};

#endif
