#ifndef BINARYPROCESSES_H
#define BINARYPROCESSES_H

#include <vector>

#include "H5.h"
#include "BinaryProcess.h"

class BinaryProcesses
{

public:
    BinaryProcesses(
        Params &params,
        std::vector<unsigned int> species_group1,
        std::vector<unsigned int> species_group2,
        bool intra,
        std::vector<BinaryProcess*> processes,
        int every,
        int debug_every,
        double time_frozen,
        std::string filename
    );
    
    BinaryProcesses( BinaryProcesses *BPs );
    
    ~BinaryProcesses();
    
    //! Method to calculate the Debye length in each bin
    static void calculate_debye_length( Params &, Patch * );
    
    //! True if any of the BinaryProcesses objects need automatically-computed coulomb log
    static bool debye_length_required_;
    
    //! Apply processes at each timestep
    void apply( Params &, Patch *, int, std::vector<Diagnostic *> & );
    
    //! Outputs the debug info if requested
    static void debug( Params &params, int itime, unsigned int icoll, VectorPatch &vecPatches );
    
    H5Write * debug_file_;
    
    //! List of processes
    std::vector<BinaryProcess*> processes_;
    
private:
    
    //! First group of species
    std::vector<unsigned int> species_group1_;
    
    //! Second group of species
    std::vector<unsigned int> species_group2_;
    
    //! Whether the two groups are the same
    bool intra_;
    
    //! Number of timesteps between each pairing
    unsigned int every_;
    
    //! Number of timesteps between each debugging file output
    unsigned int debug_every_;
    
    //! Time before which binary processes do not happen
    double timesteps_frozen_;
    
    //! Debugging file name
    std::string filename_;
    
};

#endif
