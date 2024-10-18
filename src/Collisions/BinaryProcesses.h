#ifndef BINARYPROCESSES_H
#define BINARYPROCESSES_H

#include <vector>

#include "H5.h"
#include "CollisionalNuclearReaction.h"
#include "Collisions.h"
#include "CollisionalIonization.h"

class CellVolumeCalculator
{
public:
// On GPU, particles are sorted in dual cells,
// but on CPU, they are sorted in primal cells
    CellVolumeCalculator( Params & params, Patch * patch ) :
        isAM_( params.geometry == "AMcylindrical" ),
        ndim_( params.nDim_field ),
#ifdef SMILEI_ACCELERATOR_GPU
        size_{ params.patch_size_[0], params.patch_size_[1], params.patch_size_[2] },
#else
        size_{ params.patch_size_[0] + 1, params.patch_size_[1] + 1, params.patch_size_[2] + 1 },
#endif
        dual_cell_volume_( params.cell_volume ),
        rmin_( 0 ),
        dr_( 0 ),
        nbin_( 1 )
    {
        for( int idim = 0; idim < ndim_; idim++ ) {
            nbin_ *= size_[idim]; 
        }
        
        if( isAM_ ) {
#ifdef SMILEI_ACCELERATOR_GPU
            rmin_ = patch->getDomainLocalMin( 1 );
#else
            rmin_ = patch->getDomainLocalMin( 1 ) + 0.5 * params.cell_length[1];
#endif
            dr_ = params.cell_length[1];
        }
    };
    
    #pragma acc routine seq
    double operator()( size_t ibin ) {
        
        double volume = dual_cell_volume_;
        size_t rem = ibin;
        
#ifdef SMILEI_ACCELERATOR_GPU
        if( isAM_ ) {
            size_t ir = rem % size_[1];
            volume *= rmin_ + ir * dr_;
        }
#else
        if( isAM_ ) {
            
            size_t ir = rem % size_[1];
            rem /= size_[1];
            if( ir == 0 || ir == size_[1] - 1 ) {
                volume *= 0.5;
            }
            size_t ix = rem % size_[0];
            if( ix == 0 || ix == size_[0] - 1 ) {
                volume *= 0.5;
            }
            volume *= rmin_ + ir * dr_;
            
        } else {
            
            for( int idim = ndim_ - 1; idim >= 0; idim-- ) {
                size_t i = rem % size_[idim];
                rem /= size_[idim];
                if( i == 0 || i == size_[idim] - 1 ) {
                    volume *= 0.5;
                }
            }
        }
#endif
        
        return volume;
    }
    
    size_t nbin_;
private:
    bool isAM_;
    int ndim_;
    size_t size_[3];
    double dual_cell_volume_;
    double rmin_;
    double dr_;
};

class BinaryProcesses
{

public:
    BinaryProcesses(
        Params &params,
        std::vector<unsigned int> species_group1,
        std::vector<unsigned int> species_group2,
        bool intra,
        int screening_group,
        CollisionalNuclearReaction * nuclear_reactions,
        double clog,
        double coulomb_log,
        CollisionalIonization * collisional_ionization,
        int every,
        int debug_every,
        double time_frozen,
        std::string filename
    );
    
    BinaryProcesses( BinaryProcesses *BPs );
    
    ~BinaryProcesses();
    
    //! Method to calculate the Debye length in each bin
    void calculate_debye_length( Params &, Patch * );
    
    //! True if any of the BinaryProcesses objects need automatically-computed coulomb log
    static bool debye_length_required_;
    
    //! Apply processes at each timestep
    void apply( Params &, Patch *, int, std::vector<Diagnostic *> & );
    
    //! Outputs the debug info if requested
    static void debug( Params &params, int itime, unsigned int icoll, VectorPatch &vecPatches );
    
    H5Write * debug_file_;
    
    //! Processes
    CollisionalNuclearReaction * nuclear_reactions_;
    Collisions collisions_;
    CollisionalIonization * collisional_ionization_;
    
private:
    
    //! First group of species
    std::vector<unsigned int> species_group1_;
    
    //! Second group of species
    std::vector<unsigned int> species_group2_;
    
    //! Whether the two groups are the same
    bool intra_;
    
    //! Group that makes atomic screening (e-i collisions)
    int screening_group_;
    
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
