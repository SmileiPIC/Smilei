
#include <iostream>

#include <nvidiaParticles.h>

struct count_if_out
{
    __host__ __device__
  bool operator()(const int x)
  {
    return (x  == -1);
  }
};

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

void nvidiaParticles::initGPU()
{
    if (!Position.size()) {
        std::cout << "Set CPU first" << std::endl;
        return ;
    }

    double res_fac = 1.;
    nvidiaPosition.resize( Position.size() );
    for (int idim=0;idim<Position.size();idim++) {
        nvidiaPosition[idim].reserve( res_fac*Position[idim].size() );
        nvidiaPosition[idim].resize( Position[idim].size() );
    }
    nvidiaMomentum.resize( Momentum.size() );
    for (int idim=0;idim<Momentum.size();idim++) {
        nvidiaMomentum[idim].reserve( res_fac*Momentum[idim].size() );
        nvidiaMomentum[idim].resize( Momentum[idim].size() );
    }
    nvidiaWeight.reserve( res_fac*Weight.size() );
    nvidiaWeight.resize( Weight.size() );
    nvidiaCharge.reserve( res_fac*Charge.size() );
    nvidiaCharge.resize( Charge.size() );

    if( isQuantumParameter ) {
        nvidia_chi_.reserve( res_fac*Chi.size() );
        nvidia_chi_.resize( Chi.size() );
    }

    nvidia_cell_keys.reserve( res_fac*Charge.size() );
    nvidia_cell_keys.resize( Charge.size() );
    gpu_nparts_ = Charge.size();

    if (gpu_nparts_!=0)
        syncGPU();
    else {
        for (int idim=0;idim<Position.size();idim++)
            nvidiaPosition[idim].reserve( 100 );
        for (int idim=0;idim<Momentum.size();idim++)
            nvidiaMomentum[idim].reserve( 100 );
        nvidiaWeight.reserve( 100 );
        nvidiaCharge.reserve( 100 );
        if( isQuantumParameter ) {
            nvidia_chi_.reserve( 100 );
        }
        nvidia_cell_keys.reserve( 100 );
    }

}

void nvidiaParticles::syncCPU()
{
    for (int idim=0;idim<Position.size();idim++) {
        Position[idim].resize( gpu_nparts_ );
        thrust::copy((nvidiaPosition[idim]).begin(), (nvidiaPosition[idim]).begin()+gpu_nparts_, (Position[idim]).begin());
    }
    for (int idim=0;idim<Momentum.size();idim++) {
        Momentum[idim].resize( gpu_nparts_ );
        thrust::copy((nvidiaMomentum[idim]).begin(), (nvidiaMomentum[idim]).begin()+gpu_nparts_, (Momentum[idim]).begin());
    }
    Weight.resize( gpu_nparts_ );
    thrust::copy((nvidiaWeight).begin(), (nvidiaWeight).begin()+gpu_nparts_, (Weight).begin());
    Charge.resize( gpu_nparts_ );
    thrust::copy((nvidiaCharge).begin(), (nvidiaCharge).begin()+gpu_nparts_, (Charge).begin());
    if (isQuantumParameter) {
        Chi.resize( gpu_nparts_ );
        thrust::copy((nvidia_chi_).begin(), (nvidia_chi_).begin()+gpu_nparts_, (Chi).begin());
    }
}

void nvidiaParticles::syncGPU()
{
    for (int idim=0;idim<Position.size();idim++) {
        nvidiaPosition[idim].resize( Position[idim].size() );
        thrust::copy((Position[idim]).begin(), (Position[idim]).end(), (nvidiaPosition[idim]).begin());
    }
    for (int idim=0;idim<Momentum.size();idim++) {
        nvidiaMomentum[idim].resize( Momentum[idim].size() );
        thrust::copy((Momentum[idim]).begin(), (Momentum[idim]).end(), (nvidiaMomentum[idim]).begin());
    }
    nvidiaWeight.resize( Weight.size() );
    thrust::copy((Weight).begin(), (Weight).end(), (nvidiaWeight).begin());
    nvidiaCharge.resize( Charge.size() );
    thrust::copy((Charge).begin(), (Charge).end(), (nvidiaCharge).begin());
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
    ZipIterParts iter(thrust::make_tuple(nvidiaPosition[0].begin(),
                                         nvidiaPosition[1].begin(),
                                         nvidiaPosition[2].begin(),
                                         nvidiaMomentum[0].begin(),
                                         nvidiaMomentum[1].begin(),
                                         nvidiaMomentum[2].begin(),
                                         nvidiaWeight.begin(),
                                         nvidiaCharge.begin()));

    nparts_to_move_ = thrust::count(thrust::device, nvidia_cell_keys.begin(), nvidia_cell_keys.begin()+nparts, -1);

    // Manage the send data structure
    nvidiaParticles* cp_parts = static_cast<nvidiaParticles*>( particles_to_move );
    // Resize it, if too small (copy_if do not resize)
    if ( nparts_to_move_ > cp_parts->nvidiaWeight.size() ) {
        cp_parts->nvidiaPosition[0].resize( nparts_to_move_ );
        cp_parts->nvidiaPosition[1].resize( nparts_to_move_ );
        cp_parts->nvidiaPosition[2].resize( nparts_to_move_ );
        cp_parts->nvidiaMomentum[0].resize( nparts_to_move_ );
        cp_parts->nvidiaMomentum[1].resize( nparts_to_move_ );
        cp_parts->nvidiaMomentum[2].resize( nparts_to_move_ );
        cp_parts->nvidiaWeight.resize( nparts_to_move_ );
        cp_parts->nvidiaCharge.resize( nparts_to_move_ );
        if (isQuantumParameter) {
            cp_parts->nvidia_chi_.resize( nparts_to_move_ );
        }
        cp_parts->nvidia_cell_keys.resize( nparts_to_move_ );
    }
    
    
    
    // Iterator of the send data structure (once it has been resized)
    ZipIterParts iter_copy(thrust::make_tuple(cp_parts->nvidiaPosition[0].begin(),
                                              cp_parts->nvidiaPosition[1].begin(),
                                              cp_parts->nvidiaPosition[2].begin(),
                                              cp_parts->nvidiaMomentum[0].begin(),
                                              cp_parts->nvidiaMomentum[1].begin(),
                                              cp_parts->nvidiaMomentum[2].begin(),
                                              cp_parts->nvidiaWeight.begin(),
                                              cp_parts->nvidiaCharge.begin() ) );

    // Copy send particles in dedicated data structure if nvidia_cell_keys=0 (currently = 1 if keeped, new PartBoundCond::apply(...))
    thrust::copy_if(thrust::device, iter, iter+nparts, nvidia_cell_keys.begin(), iter_copy, count_if_out());
    // Special treatment for chi if radiation emission
    if (isQuantumParameter) {
        thrust::copy_if(thrust::device, nvidia_chi_.begin(), nvidia_chi_.begin()+nparts, nvidia_cell_keys.begin(), cp_parts->nvidia_chi_.begin(), count_if_out());
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
                                                                   nvidiaPosition[0].begin(), nvidiaPosition[1].begin(), nvidiaPosition[2].begin(),
                                                                   nvidiaMomentum[0].begin(), nvidiaMomentum[1].begin(), nvidiaMomentum[2].begin(),
                                                                   nvidiaWeight.begin(), nvidiaCharge.begin(), nvidia_cell_keys.begin()
                                                                   )
                                                ),
                      thrust::make_zip_iterator(thrust::make_tuple(
                                                                   nvidiaPosition[0].begin()+nparts, nvidiaPosition[1].begin()+nparts, nvidiaPosition[2].begin()+nparts,
                                                                   nvidiaMomentum[0].begin()+nparts, nvidiaMomentum[1].begin()+nparts, nvidiaMomentum[2].begin()+nparts,
                                                                   nvidiaWeight.begin()+nparts, nvidiaCharge.begin()+nparts, nvidia_cell_keys.begin()+nparts
                                                                   )
                                                ),
                      remove_if_out()
                      );
    // Update current nmuber of particles
    gpu_nparts_ -= nparts_to_move_;
    nparts = gpu_nparts_;


    // Manage the recv data structure
    nvidiaParticles* cp_parts = static_cast<nvidiaParticles*>( particles_to_move );

    ZipIterParts iter_copy(thrust::make_tuple(cp_parts->nvidiaPosition[0].begin(), cp_parts->nvidiaPosition[1].begin(), cp_parts->nvidiaPosition[2].begin(),
                                              cp_parts->nvidiaMomentum[0].begin(), cp_parts->nvidiaMomentum[1].begin(), cp_parts->nvidiaMomentum[2].begin(),
                                              cp_parts->nvidiaWeight.begin(), cp_parts->nvidiaCharge.begin() ) );
    int nparts_add = cp_parts->gpu_nparts_;
    int tot_parts = nparts + nparts_add;
 
    // Resize main data structure, if too small (copy_n do not resize)
    if ( tot_parts > nvidiaWeight.size() ) {
        nvidiaPosition[0].resize( tot_parts, 0. );
        nvidiaPosition[1].resize( tot_parts, 0. );
        nvidiaPosition[2].resize( tot_parts, 0. );
        nvidiaMomentum[0].resize( tot_parts, 0. );
        nvidiaMomentum[1].resize( tot_parts, 0. );
        nvidiaMomentum[2].resize( tot_parts, 0. );
        nvidiaWeight.resize( tot_parts, 0. );
        nvidiaCharge.resize( tot_parts, 0 );
        nvidia_cell_keys.resize( tot_parts, 0 );
    }
    // Iterator of the main data structure (once it has been resized)
    ZipIterParts iter(thrust::make_tuple(nvidiaPosition[0].begin(), nvidiaPosition[1].begin(), nvidiaPosition[2].begin(),
                                         nvidiaMomentum[0].begin(), nvidiaMomentum[1].begin(), nvidiaMomentum[2].begin(),
                                         nvidiaWeight.begin(), nvidiaCharge.begin() ) );
 
 
    // Copy recv particles in main data structure
    thrust::copy_n(thrust::device, iter_copy, nparts_add, iter+nparts);
    gpu_nparts_ += nparts_add;

    // Resize below useless : nvidia_cell_keys resized if necessary above, cell_keys not used on cpu
    //nvidia_cell_keys.resize( gpu_nparts_ );
    //cell_keys.resize       ( gpu_nparts_ );
 
    cp_parts->gpu_nparts_ = 0;
    nparts_add -= nparts_to_move_; // update return value to update last_index
    nparts_to_move_ = 0;

    return nparts_add;
}
