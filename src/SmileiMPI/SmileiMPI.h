#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <mpi.h>

#include <string>
#include <vector>

#include "Field.h"
#include "Particles.h"
#include "Tools.h"
#include "gpu.h"

class Params;
class Species;
class VectorPatch;
class DomainDecomposition;

class ElectroMagn;
class ProbeParticles;

class Diagnostic;
class DiagnosticScalar;
class DiagnosticParticleBinning;
class DiagnosticScreen;
class DiagnosticRadiationSpectrum;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiMPI
//  --------------------------------------------------------------------------------------------------------------------
class SmileiMPI
{
    friend class Checkpoint;
    friend class PatchesFactory;
    friend class Patch;
    friend class VectorPatch;
    friend class SimWindow;
    friend class AsyncMPIbuffers;

public:
    SmileiMPI() {};
    //! Create intial MPI environment
    SmileiMPI( int *argc, char ***argv );
    //! Destructor for SmileiMPI
    virtual ~SmileiMPI();

    // Broadcast a string in current communicator
    void bcast( std::string &val );
    // Broadcast an int in current communicator
    void bcast( int &val );

    //! Initialize  MPI (per process) environment
    //! \param params Parameters
    virtual void init( Params &params, DomainDecomposition *domain_decomposition );

    // Initialize the patch_count vector. Patches are distributed in order to balance the load between MPI processes.
    virtual void init_patch_count( Params &params, DomainDecomposition *domain_decomposition );

    // Recompute the patch_count vector. Browse patches and redistribute them in order to balance the load between MPI processes.
    void recompute_patch_count( Params &params, VectorPatch &vecpatches, double time_dual );
    // Returns the rank of the MPI process currently owning patch h.
    int hrank( int h );

    // Create MPI type to exchange all particles properties of particles
    MPI_Datatype createMPIparticles( Particles *particles );


    // PATCH SEND / RECV METHODS
    //     - during load balancing process
    //     - during moving window
    // -----------------------------------
    void isend( Patch *patch, int to, int tag, Params &params, bool send_xmax_bc = true );
    void waitall( Patch *patch );
    void recv( Patch *patch, int from, int tag, Params &params, bool recv_xmin_bc = true);

    void isend_fields( Patch *patch, int to, int &irequest, int tag, Params &params, bool send_xmax_bc = true );
    void recv_fields( Patch *patch, int from, int &tag, Params &params, bool recv_xmin_bc = true );
    void isend_species( Patch *patch, int to, int &irequest, int tag, Params &params );
    void recv_species( Patch *patch, int from, int &tag, Params &params );

    void isend( Particles *particles, int to, int tag, MPI_Datatype datatype, MPI_Request &request );
    void recv( Particles *partictles, int from, int tag, MPI_Datatype datatype );
    void isend( std::vector<int> *vec, int to, int tag, MPI_Request &request );
    void recv( std::vector<int> *vec, int from, int tag );

    void isend( std::vector<double> *vec, int to, int tag, MPI_Request &request );
    void recv( std::vector<double> *vec, int from, int tag );

    //Sending and reveiving ElectroMagn and ElectroMagAM
    void isend( ElectroMagn *fields, int to, int &irequest, std::vector<MPI_Request> &requests, int tag, bool send_xmax_bc );
    void isend( ElectroMagn *fields, int to, int &irequest, std::vector<MPI_Request> &requests, int tag, unsigned int nmodes, bool send_xmax_bc );
    void recv( ElectroMagn *fields, int from, int &tag, bool recv_xmax_bc );
    void recv( ElectroMagn *fields, int from, int &tag, unsigned int nmodes, bool recv_xmax_bc );

    //Templates to send/receive PML for both 2D and 3D
    template <typename Tpml>
    int  recv_PML(ElectroMagn *EM, Tpml embc, int bcId, int from, int tag, bool recv_xmin_bc);
    template <typename Tpml>
    void  send_PML(ElectroMagn *EM, Tpml embc, int bcId, int to, int &irequest, std::vector<MPI_Request> &requests, int tag, bool send_xmax_bc);


    //Sending and reveiving Fields and cFields
    //! Sends the whole Field
    void isend( Field *field, int to, int tag, MPI_Request &request );
    //! Sends the whole Field Device to Device (assuming MPI enables it)
#if defined (SMILEI_ACCELERATOR_GPU)
    void isendOnDevice( Field *field, int to, int tag, MPI_Request &request );
#endif

    void isend( Field *field, int to, int tag, MPI_Request &request, int x_first );       // Sends the first "x_first" columns of the Field
    void isendComplex( Field *field, int to, int tag, MPI_Request &request );          // Sends the whole cField
    void isendComplex( Field *field, int to, int tag, MPI_Request &request, int x_first );// Sends only the first x_first columns of the cField

    //! Receives the whole Field
    void recv( Field *field, int from, int tag);     
    //! Receives the whole Field Device to Device (assuming MPI enables it)
#if defined (SMILEI_ACCELERATOR_GPU)
    void recvOnDevice( Field *field, int from, int tag);     
#endif

    void recvShifted( Field *field, int from, int tag, int xshift ); //Shifts the reception adress by xshift columns and reduces the reception buffer size
    void recvComplex( Field *field, int from, int tag);              //Receives the whole cField
    void recvComplexShifted( Field *field, int from, int tag, int xshift ); //Shifts the reception adress by xshift columns and reduces the reception buffer size

    void sendComplex( Field *field, int to, int tag );
    void irecvComplex( Field *field, int from, int tag, MPI_Request &request );

    void isend( ProbeParticles *probe, int to, int tag, unsigned int );
    void recv( ProbeParticles *probe, int from, int tag, unsigned int );

    void isend( int *integer, int to, int tag, MPI_Request &request );
    void recv( int *integer, int from, int tag );

    // Functions for double grid exchange
    void send( Field* field, int to  , int tag );
    void irecv( Field* field, int from, int tag, MPI_Request& request );

    // DIAGS MPI SYNC
    // --------------

    // Wrapper of MPI synchronization of all computing diags
    void computeGlobalDiags(Diagnostic*                  diag, int timestep);
    // MPI synchronization of scalars diags
    void computeGlobalDiags(DiagnosticScalar*            diag, int timestep);
    // MPI synchronization of diags particles
    void computeGlobalDiags(DiagnosticParticleBinning*   diag, int timestep);
    // MPI synchronization of screen diags
    void computeGlobalDiags(DiagnosticScreen*            diag, int timestep);
    // MPI synchronization of radiation spectrum diags
    void computeGlobalDiags(DiagnosticRadiationSpectrum* diag, int timestep);

    // MPI basic methods
    // -----------------

    //! Method to identify the rank 0 MPI process
    inline bool isMaster()
    {
        return ( smilei_rk==0 );
    }
    //! Method to synchronize MPI process in the current MPI communicator
    inline void barrier()
    {
        MPI_Barrier( world_ );
    }
    //! Return MPI_Comm_rank
    inline int getRank()
    {
        return smilei_rk;
    }
    //! Return MPI_Comm_size
    inline int getSize()
    {
        return smilei_sz;
    }

    //! Return MPI_Comm_world
    inline MPI_Comm& world()
    {
        return world_;
    }

    //! Return omp_max_threads
    inline int getOMPMaxThreads()
    {
        return smilei_omp_max_threads;
    }

    //! Return local number of cores
    inline int getNumCores()
    {
        return number_of_cores;
    }

    //! Return global number of cores
    inline int getGlobalNumCores()
    {
        return global_number_of_cores;
    }

    //! Return tag upper bound of this MPI implementation
    inline int getTagUB()
    {
        int flag;
        int* tag_ub_ptr;
        MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub_ptr, &flag);
        return *tag_ub_ptr;
    }

    // Global buffers for vectorization of Species::dynamics
    // -----------------------------------------------------

    //! value of the Efield
    std::vector<std::vector<double>> dynamics_Epart;
    //! value of the Bfield
    std::vector<std::vector<double>> dynamics_Bpart;
    //! gamma factor
    std::vector<std::vector<double>> dynamics_invgf;
    //! iold_pos
    std::vector<std::vector<int>> dynamics_iold;
    //! delta_old_pos
    std::vector<std::vector<double>> dynamics_deltaold;
    //! theta old
    std::vector<std::vector<std::complex<double>>> dynamics_eithetaold;
    //! value of the By field for BTIS3
    std::vector<std::vector<double>> dynamics_Bpart_yBTIS3;
    //! value of the Bz field for BTIS3
    std::vector<std::vector<double>> dynamics_Bpart_zBTIS3;

    //! value of the grad(AA*) at itime and itime-1
    std::vector<std::vector<double>> dynamics_GradPHIpart;
    std::vector<std::vector<double>> dynamics_GradPHI_mpart;
    //! value of the AA* at itime and itime-1
    std::vector<std::vector<double>> dynamics_PHIpart;
    std::vector<std::vector<double>> dynamics_PHI_mpart;
    //! inverse of the ponderomotive gamma, used in susceptibility and ponderomotive momentum Pusher
    std::vector<std::vector<double>> dynamics_inv_gamma_ponderomotive;
    //! value of the EnvEabs used for envelope ionization
    std::vector<std::vector<double>> dynamics_EnvEabs_part;
    //! value of the EnvEabs used for envelope ionization
    std::vector<std::vector<double>> dynamics_EnvExabs_part;

    //! Return buffer size in thread ithread
    inline int __attribute__((always_inline)) getBufferSize(const int ithread)
    {
        return dynamics_invgf[ithread].size();
    }

    //! Erase Particles from istart ot the end in the buffers of thread ithread
    void eraseBufferParticleTrail( const int ndim, const int istart, const int ithread, bool isAM = false );

#if defined( SMILEI_ACCELERATOR_GPU_OMP ) || defined( SMILEI_ACCELERATOR_GPU_OACC )
    //! Map CPU buffers onto the GPU to at least accommodate particle_count
    //! particles. This method tries to reduce the number of
    //! allocation/deallocation which produces a lot of fragmentation on some
    //! GPUs. In fact, it emulates the reserve behavior of an std::vector.
    //!
    //! If you change the capacity of dynamics_Epart, dynamics_Bpart,
    //! dynamics_invgf, dynamics_iold or dynamics_deltaold vector outside of
    //! this method, you expose yourself to memory leak or corruption (in
    //! GPU mode).
    //!
    //! TODO(Etienne M): FREE To avoid a leak at the end of the program or when
    //! a SmileiMPI is destroyed. We should free the device memory. This can be
    //! done using a resizeDeviceBuffers( ithread, ndim_field, 0).
    //!
    void resizeDeviceBuffers( unsigned int ithread,
                              unsigned int ndim_field,
                              unsigned int particle_count,
                              float        growth_factor = 1.3F );
#endif

    inline void resizeBuffers( int ithread, int ndim_field, int npart, bool isAM = false )
    {
        dynamics_Epart[ithread].resize( 3*npart );
        dynamics_Bpart[ithread].resize( 3*npart );
        dynamics_invgf[ithread].resize( npart );
        dynamics_iold[ithread].resize( ndim_field*npart );
        dynamics_deltaold[ithread].resize( ndim_field*npart );
        if(use_BTIS3){
            dynamics_Bpart_yBTIS3[ithread].resize( npart );
            dynamics_Bpart_zBTIS3[ithread].resize( npart );
        }
        if( isAM ) {
            dynamics_eithetaold[ithread].resize( npart );
        }

        if( dynamics_GradPHIpart.size() > 0 ) {
            dynamics_GradPHIpart[ithread].resize( 3*npart );
            dynamics_GradPHI_mpart[ithread].resize( 3*npart );
            dynamics_PHIpart[ithread].resize( npart );
            dynamics_PHI_mpart[ithread].resize( npart );
            dynamics_inv_gamma_ponderomotive[ithread].resize( npart );
            if ( dynamics_EnvEabs_part.size() > 0 ){
                dynamics_EnvEabs_part[ithread].resize( npart );
                dynamics_EnvExabs_part[ithread].resize( npart );
            }
        }
    }


        // Resize buffers vector for a given number of buffers
    inline void resizeBuffers( int n_buffers, bool isAM = false)
    {
        dynamics_Epart.resize( n_buffers );
        dynamics_Bpart.resize( n_buffers );
        dynamics_invgf.resize( n_buffers );
        dynamics_iold.resize( n_buffers );
        dynamics_deltaold.resize( n_buffers );
        if(use_BTIS3){
            dynamics_Bpart_yBTIS3.resize( n_buffers );
            dynamics_Bpart_zBTIS3.resize( n_buffers );
        }
        if( isAM ) {
            dynamics_eithetaold.resize( n_buffers );
        }

        if( dynamics_GradPHIpart.size() > 0 ) {
            dynamics_GradPHIpart.resize( n_buffers );
            dynamics_GradPHI_mpart.resize( n_buffers );
            dynamics_PHIpart.resize( n_buffers );
            dynamics_PHI_mpart.resize( n_buffers );
            dynamics_inv_gamma_ponderomotive.resize( n_buffers );
            if ( dynamics_EnvEabs_part.size() > 0 ){
                dynamics_EnvEabs_part.resize( n_buffers );
                dynamics_EnvExabs_part.resize( n_buffers );
            }
        }
    }

    // Resize buffers to avoid memory leak with tasks
    inline void reduceDynamicsBufferSize( int buffer_id, bool isAM = false )
    {
        dynamics_Epart[buffer_id].resize( 1 );
        dynamics_Bpart[buffer_id].resize( 1 );
        dynamics_invgf[buffer_id].resize( 1 );
        dynamics_iold[buffer_id].resize( 1 );
        dynamics_deltaold[buffer_id].resize( 1 );
        if(use_BTIS3){
            dynamics_Bpart_yBTIS3[buffer_id].resize( 1 );
            dynamics_Bpart_zBTIS3[buffer_id].resize( 1 );
        }
        if( isAM ) {
            dynamics_eithetaold[buffer_id].resize( 1 );
        }

        if( dynamics_GradPHIpart.size() > 0 ) {
            dynamics_GradPHIpart[buffer_id].resize( 1 );
            dynamics_GradPHI_mpart[buffer_id].resize( 1 );
            dynamics_PHIpart[buffer_id].resize( 1 );
            dynamics_PHI_mpart[buffer_id].resize( 1 );
            dynamics_inv_gamma_ponderomotive[buffer_id].resize( 1 );
            if ( dynamics_EnvEabs_part.size() > 0 ){
                dynamics_EnvEabs_part[buffer_id].resize( 1 );
                dynamics_EnvExabs_part[buffer_id].resize( 1 );
            }
        }
    }

    // Resize buffers for old properties only
    inline void resizeOldPropertiesBuffer( int ithread, int ndim_field, int npart, bool isAM = false )
    {
        dynamics_iold[ithread].resize( ndim_field*npart );
        dynamics_deltaold[ithread].resize( ndim_field*npart );
        if( isAM ) {
            dynamics_eithetaold[ithread].resize( npart );
        }
    }

    bool test_mode;

    // Task tracing diag
    std::vector<std::vector<double>> particle_event_tracing_event_time_;
    std::vector<std::vector<unsigned int>> particle_event_tracing_start_or_end_;
    std::vector<std::vector<int>> particle_event_tracing_event_name_;
    int iter_frequency_particle_event_tracing_;
    double reference_time_;

    // determine if "task" tracing is performed at this iteration
    bool diagPartEventTracing(double time_dual, double timestep ){
        bool diagTracing = false;
        if (int((time_dual-0.5*timestep)/timestep)%(iter_frequency_particle_event_tracing_)==0){
            diagTracing = true;
        }
        return diagTracing;
    }
    // trace event or "task"
    void trace_event(int thread, double event_time,unsigned int event_start_or_end, int event_name)
    {
        particle_event_tracing_event_time_[thread].push_back(event_time);           // write time
        particle_event_tracing_start_or_end_[thread].push_back(event_start_or_end); // write Start/End
        particle_event_tracing_event_name_[thread].push_back(event_name);           // write Event Name
    };

    // If particle event tracing diagnostic is activated, trace event
#ifdef _PARTEVENTTRACING
    void traceEventIfDiagTracing( bool, int, unsigned int, int ) {};
#else
    void traceEventIfDiagTracing( bool diag_PartEventTracing, int thread,
                                  unsigned int event_start_or_end, int event_name )
    {
        if( diag_PartEventTracing ) trace_event( thread, (MPI_Wtime()-reference_time_), event_start_or_end, event_name );
    };
#endif

    bool use_BTIS3;

protected:
    //! Global MPI Communicator
    MPI_Comm world_;

    //! Number of MPI process in the current communicator
    int smilei_sz;
    //! MPI process Id in the current communicator
    int smilei_rk;
    //! OMP max number of threads in one MPI
    int smilei_omp_max_threads;
    //! OMP available cores in one MPI
    int number_of_cores;
    //! Global number of cores
    int global_number_of_cores;

    // Store periodicity (0/1) per direction
    // Should move in Params : last parameters of this type in this class
    int *periods_;

    //! For patch decomposition
    //Number of patches owned by each mpi process.
    std::vector<int>  patch_count, capabilities, patch_refHindexes;
    int Tcapabilities; //Default = smilei_sz (1 per MPI rank)
};


#endif
