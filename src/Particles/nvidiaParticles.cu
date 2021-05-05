
#include <iostream>

#include <nvidiaParticles.h>

//! Structure with specific function count_if_out for thrust::tuple operator
//! Return True if the entry is -1 as in the cell keys vector for instance
struct count_if_out
{
    __host__ __device__
  bool operator()(const int x)
  {
    return (x  == -1);
  }
};

//! Structure with specific function thrust::tuple operator
struct remove_if_out
{
    typedef thrust::tuple<double,double,double,double,double,double,double,short,int> Tuple;
    __host__ __device__
    bool operator()(const Tuple&  t)
    {
        // thrust::get<8>(t) = cell_keys
        return (thrust::get<8>(t) < 0);

    }
};

nvidiaParticles::nvidiaParticles() : Particles()
{
    gpu_nparts_ = 0;

}

void nvidiaParticles::initializeDataOnDevice()
{
    if (!Position.size()) {
        std::cout << "Set CPU first" << std::endl;
        return ;
    }

    double res_fac = 1.;
    nvidia_position_.resize( Position.size() );
    for (int idim=0;idim<Position.size();idim++) {
        nvidia_position_[idim].reserve( res_fac*Position[idim].size() );
        nvidia_position_[idim].resize( Position[idim].size() );
    }
    nvidia_momentum_.resize( Momentum.size() );
    for (int idim=0;idim<Momentum.size();idim++) {
        nvidia_momentum_[idim].reserve( res_fac*Momentum[idim].size() );
        nvidia_momentum_[idim].resize( Momentum[idim].size() );
    }
    nvidia_weight_.reserve( res_fac*Weight.size() );
    nvidia_weight_.resize( Weight.size() );
    nvidia_charge_.reserve( res_fac*Charge.size() );
    nvidia_charge_.resize( Charge.size() );

    if( isQuantumParameter ) {
        nvidia_chi_.reserve( res_fac*Chi.size() );
        nvidia_chi_.resize( Chi.size() );
    }

    nvidia_cell_keys_.reserve( res_fac*Charge.size() );
    nvidia_cell_keys_.resize( Charge.size() );
    gpu_nparts_ = Charge.size();

    if (gpu_nparts_!=0)
        syncGPU();
    else {
        for (int idim=0;idim<Position.size();idim++)
            nvidia_position_[idim].reserve( 100 );
        for (int idim=0;idim<Momentum.size();idim++)
            nvidia_momentum_[idim].reserve( 100 );
        nvidia_weight_.reserve( 100 );
        nvidia_charge_.reserve( 100 );
        if( isQuantumParameter ) {
            nvidia_chi_.reserve( 100 );
        }
        nvidia_cell_keys_.reserve( 100 );
    }

}

void nvidiaParticles::syncCPU()
{
    for (int idim=0;idim<Position.size();idim++) {
        Position[idim].resize( gpu_nparts_ );
        thrust::copy((nvidia_position_[idim]).begin(), (nvidia_position_[idim]).begin()+gpu_nparts_, (Position[idim]).begin());
    }
    for (int idim=0;idim<Momentum.size();idim++) {
        Momentum[idim].resize( gpu_nparts_ );
        thrust::copy((nvidia_momentum_[idim]).begin(), (nvidia_momentum_[idim]).begin()+gpu_nparts_, (Momentum[idim]).begin());
    }
    Weight.resize( gpu_nparts_ );
    thrust::copy((nvidia_weight_).begin(), (nvidia_weight_).begin()+gpu_nparts_, (Weight).begin());
    Charge.resize( gpu_nparts_ );
    thrust::copy((nvidia_charge_).begin(), (nvidia_charge_).begin()+gpu_nparts_, (Charge).begin());
    if (isQuantumParameter) {
        Chi.resize( gpu_nparts_ );
        thrust::copy((nvidia_chi_).begin(), (nvidia_chi_).begin()+gpu_nparts_, (Chi).begin());
    }
}

void nvidiaParticles::syncGPU()
{
    for (int idim=0;idim<Position.size();idim++) {
        nvidia_position_[idim].resize( Position[idim].size() );
        thrust::copy((Position[idim]).begin(), (Position[idim]).end(), (nvidia_position_[idim]).begin());
    }
    for (int idim=0;idim<Momentum.size();idim++) {
        nvidia_momentum_[idim].resize( Momentum[idim].size() );
        thrust::copy((Momentum[idim]).begin(), (Momentum[idim]).end(), (nvidia_momentum_[idim]).begin());
    }
    nvidia_weight_.resize( Weight.size() );
    thrust::copy((Weight).begin(), (Weight).end(), (nvidia_weight_).begin());
    nvidia_charge_.resize( Charge.size() );
    thrust::copy((Charge).begin(), (Charge).end(), (nvidia_charge_).begin());
    if (isQuantumParameter) {
        nvidia_chi_.resize( Chi.size() );
        thrust::copy((Chi).begin(), (Chi).end(), (nvidia_chi_).begin());
    }
    gpu_nparts_ = Charge.size();
}

void nvidiaParticles::extractParticles( Particles* particles_to_move )
{
    int nparts = gpu_nparts_;

    // Iterator of the main data structure
    ZipIterParts iter(thrust::make_tuple(nvidia_position_[0].begin(),
                                         nvidia_position_[1].begin(),
                                         nvidia_position_[2].begin(),
                                         nvidia_momentum_[0].begin(),
                                         nvidia_momentum_[1].begin(),
                                         nvidia_momentum_[2].begin(),
                                         nvidia_weight_.begin(),
                                         nvidia_charge_.begin()));

    nparts_to_move_ = thrust::count(thrust::device, nvidia_cell_keys_.begin(), nvidia_cell_keys_.begin()+nparts, -1);

    // Manage the send data structure
    nvidiaParticles* cp_parts = static_cast<nvidiaParticles*>( particles_to_move );
    // Resize it, if too small (copy_if do not resize)
    if ( nparts_to_move_ > cp_parts->nvidia_weight_.size() ) {
        cp_parts->nvidia_position_[0].resize( nparts_to_move_ );
        cp_parts->nvidia_position_[1].resize( nparts_to_move_ );
        cp_parts->nvidia_position_[2].resize( nparts_to_move_ );
        cp_parts->nvidia_momentum_[0].resize( nparts_to_move_ );
        cp_parts->nvidia_momentum_[1].resize( nparts_to_move_ );
        cp_parts->nvidia_momentum_[2].resize( nparts_to_move_ );
        cp_parts->nvidia_weight_.resize( nparts_to_move_ );
        cp_parts->nvidia_charge_.resize( nparts_to_move_ );
        if (isQuantumParameter) {
            cp_parts->nvidia_chi_.resize( nparts_to_move_ );
        }
        cp_parts->nvidia_cell_keys_.resize( nparts_to_move_ );
    }
    
    
    
    // Iterator of the send data structure (once it has been resized)
    ZipIterParts iter_copy(thrust::make_tuple(cp_parts->nvidia_position_[0].begin(),
                                              cp_parts->nvidia_position_[1].begin(),
                                              cp_parts->nvidia_position_[2].begin(),
                                              cp_parts->nvidia_momentum_[0].begin(),
                                              cp_parts->nvidia_momentum_[1].begin(),
                                              cp_parts->nvidia_momentum_[2].begin(),
                                              cp_parts->nvidia_weight_.begin(),
                                              cp_parts->nvidia_charge_.begin() ) );

    // Copy send particles in dedicated data structure if nvidia_cell_keys_=0 (currently = 1 if keeped, new PartBoundCond::apply(...))
    thrust::copy_if(thrust::device, iter, iter+nparts, nvidia_cell_keys_.begin(), iter_copy, count_if_out());
    // Special treatment for chi if radiation emission
    if (isQuantumParameter) {
        thrust::copy_if(thrust::device, nvidia_chi_.begin(), nvidia_chi_.begin()+nparts, nvidia_cell_keys_.begin(), cp_parts->nvidia_chi_.begin(), count_if_out());
    }
    
    cp_parts->gpu_nparts_ = nparts_to_move_;

    particles_to_move->syncCPU();
}


int nvidiaParticles::injectParticles( Particles* particles_to_move )
{
    int nparts = gpu_nparts_;
    // Remove particles which leaves current patch
    thrust::remove_if(thrust::device,
                      thrust::make_zip_iterator(thrust::make_tuple(
                                                                   nvidia_position_[0].begin(),
                                                                   nvidia_position_[1].begin(),
                                                                   nvidia_position_[2].begin(),
                                                                   nvidia_momentum_[0].begin(),
                                                                   nvidia_momentum_[1].begin(),
                                                                   nvidia_momentum_[2].begin(),
                                                                   nvidia_weight_.begin(),
                                                                   nvidia_charge_.begin()
                                                                  // , nvidia_cell_keys_.begin()
                                                                   )
                                                ),
                      thrust::make_zip_iterator(thrust::make_tuple(
                                                                   nvidia_position_[0].begin()+nparts,
                                                                   nvidia_position_[1].begin()+nparts,
                                                                   nvidia_position_[2].begin()+nparts,
                                                                   nvidia_momentum_[0].begin()+nparts,
                                                                   nvidia_momentum_[1].begin()+nparts,
                                                                   nvidia_momentum_[2].begin()+nparts,
                                                                   nvidia_weight_.begin()+nparts,
                                                                   nvidia_charge_.begin()+nparts
                                                                //,   nvidia_cell_keys_.begin()+nparts
                                                                   )
                                                ),
                      nvidia_cell_keys_.begin(),
                      count_if_out()
                      );
                      
    if (isQuantumParameter) {
        thrust::remove_if(thrust::device,
                          nvidia_chi_.begin(),
                          nvidia_chi_.begin()+nparts,
                          nvidia_cell_keys_.begin(),
                          count_if_out()
        );
    }
    
    // Update current number of particles
    gpu_nparts_ -= nparts_to_move_;
    nparts = gpu_nparts_;

    // Just resize cell keys, no need to remove
    // nvidia_cell_keys_.resize(gpu_nparts_);

    // Manage the recv data structure
    nvidiaParticles* cp_parts = static_cast<nvidiaParticles*>( particles_to_move );

    ZipIterParts iter_copy(thrust::make_tuple(cp_parts->nvidia_position_[0].begin(),
                                              cp_parts->nvidia_position_[1].begin(),
                                              cp_parts->nvidia_position_[2].begin(),
                                              cp_parts->nvidia_momentum_[0].begin(),
                                              cp_parts->nvidia_momentum_[1].begin(),
                                              cp_parts->nvidia_momentum_[2].begin(),
                                              cp_parts->nvidia_weight_.begin(),
                                              cp_parts->nvidia_charge_.begin() ) );
    int nparts_add = cp_parts->gpu_nparts_;
    int tot_parts = nparts + nparts_add;
 
    // Resize main data structure, if too small (copy_n do not resize)
    if ( tot_parts > nvidia_weight_.size() ) {
        nvidia_position_[0].resize( tot_parts, 0. );
        nvidia_position_[1].resize( tot_parts, 0. );
        nvidia_position_[2].resize( tot_parts, 0. );
        nvidia_momentum_[0].resize( tot_parts, 0. );
        nvidia_momentum_[1].resize( tot_parts, 0. );
        nvidia_momentum_[2].resize( tot_parts, 0. );
        nvidia_weight_.resize( tot_parts, 0. );
        nvidia_charge_.resize( tot_parts, 0 );
        nvidia_cell_keys_.resize( tot_parts, 0 );
        if (isQuantumParameter) {
            nvidia_chi_.resize(  tot_parts, 0 );
        }
    }
    // Iterator of the main data structure (once it has been resized)
    ZipIterParts iter(thrust::make_tuple(nvidia_position_[0].begin(),
                                         nvidia_position_[1].begin(),
                                         nvidia_position_[2].begin(),
                                         nvidia_momentum_[0].begin(),
                                         nvidia_momentum_[1].begin(),
                                         nvidia_momentum_[2].begin(),
                                         nvidia_weight_.begin(),
                                         nvidia_charge_.begin() ) );
 
 
    // Copy recv particles in main data structure
    thrust::copy_n(thrust::device, iter_copy, nparts_add, iter+nparts);
    
    if (isQuantumParameter) {
        thrust::copy_n(thrust::device,
                       cp_parts->nvidia_chi_.begin(),
                       nparts_add,
                       nvidia_chi_.begin()+nparts);
    }
    
    gpu_nparts_ += nparts_add;

    // Resize below useless : nvidia_cell_keys_ resized if necessary above, cell_keys not used on cpu
    //nvidia_cell_keys_.resize( gpu_nparts_ );
    //cell_keys.resize       ( gpu_nparts_ );
 
    cp_parts->gpu_nparts_ = 0;
    nparts_add -= nparts_to_move_; // update return value to update last_index
    nparts_to_move_ = 0;

    return nparts_add;
}
