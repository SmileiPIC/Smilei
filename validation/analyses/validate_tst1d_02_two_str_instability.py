import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

tracked_x = S.TrackParticles.ion(axes=["x"], timesteps=0).getData()["x"][0]
Validate("Regularly spaced particles", tracked_x, 1e-7)

tracked_px = S.TrackParticles.ion(axes=["px"], timesteps=5000).getData()["px"][0]
Validate("Tracked particles px", tracked_px, 1.)

tracked_Ex = S.TrackParticles("ion",axes=["Ex"], timesteps=5000).getData()["Ex"][0]
Validate("Tracked particles Ex", tracked_Ex, 1e-7)

itimes = S.ParticleBinning.Diag0().getAvailableTimesteps()
Validate("Timesteps in particle diagnostic", itimes )

px = S.ParticleBinning.Diag0(sum={"x":"all"}, timesteps=0   ).getData()[0]
Validate("Initial electron momentum distribution", px, 0.1 )

px = S.ParticleBinning.Diag0(sum={"x":"all"}, timesteps=2000).getData()[0]
Validate("Final electron momentum distribution", px, 0.1)

Validate("List of fields in Probe", S.Probe(0).getFields() )

Validate("Species probe", S.Probe(0, 'Rho_eon1', timesteps=2000).getData()[0], 1e-4 )

# CONSISTENCY OF DIAGPERFORMANCES
#   check available quantities
Validate("Performances available quantities", sorted(S.Performances().getAvailableQuantities()))
#   check number of MPI processes
with h5py.File("./restart000/Performances.h5") as f:
	MPI_SIZE = f.attrs["MPI_SIZE"]
performances_size = len(S.Performances(raw="hindex").getData(timestep=0)[0])
Validate("Performances MPI_SIZE", MPI_SIZE == performances_size)
#   check total numbers
number_of_cells = np.sum(S.Performances(raw="number_of_cells").getData(), axis=1)
Validate("Performances number_of_cells", number_of_cells)
number_of_particles = np.sum(S.Performances(raw="number_of_particles").getData(), axis=1)
Validate("Performances number_of_particles", number_of_particles)
number_of_frozen_particles = np.sum(S.Performances(raw="number_of_frozen_particles").getData(), axis=1)
Validate("Performances number_of_frozen_particles", number_of_frozen_particles)
