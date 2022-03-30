#include "ParticlesFactory.h"

#if defined(_GPU) || defined(SMILEI_ACCELERATOR_GPU_OMP)
// TODO(Etienne M): find a way to put that in the function body, under the macro
// guards, where CreateGPUParticles is called.
extern "C" void* CreateGPUParticles();
#endif

Particles* ParticlesFactory::create(const Params& params) {
    Particles* particles = nullptr;

    if(params.gpu_computing) {

        // Because we can potentially use 2 compilers, one for the cuda/hip
        // sources (thrust), and an other for the host code, we dont want to
        // include nvidiaParticles.cu in this translation unit as it may not
        // generate the same datastructures. For instance the rocThrust lib,
        // when used with the cray compiler (a fork of clang), will chose to
        // use a cuda backend and we cant (except by using non documented
        // macros #risky) change this behavior. When using hipcc, rocThrust
        // lib will choose the hip backend. This discrepancy can lead to
        // ABI problems due to having 2 different definition of a given
        // datastructure.
        // The potential solutions are:
        // - Put the ParticlesFactory::create implementation in a separate
        // .cpp file and use the makefile to either compile using hipcc/nvcc
        // or SMILEICXX.
        // - Compile everything using the same, cuda/hip compatible
        // compiler. This is what's done for the nvidia gpu code as of
        // 16 march 22.
        // - Use an extern, stable C API, define the function in the
        // nvidiaParticles.cu implementatino file. Let the linker do it's
        // job.
        //
        // TLDR: We may use 2 compiler, one for the host code, one for the
        // thrust/gpu code. This could mean 2 different interpretation of
        // the nvidiaParticles.h header (and it's thrust content). This is
        // dangerous, we should hide everything that is related to an other
        // compiler behind a stable C API.
        //
#if defined(_GPU) || defined(SMILEI_ACCELERATOR_GPU_OMP)
        particles = static_cast<Particles*>(CreateGPUParticles());
#else
        ERROR("Unreachable state reached, Smilei was not built with GPU support!")
#endif
    } else {
        particles = new Particles();
    }

    if(particles == nullptr) {
        ERROR("particles could not be allocated!")
    }

    return particles;
}
