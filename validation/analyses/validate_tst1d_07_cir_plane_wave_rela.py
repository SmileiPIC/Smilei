# ____________________________________________________________________________
#
# This script validates the relativist pushers in 1d by analyzing 
# the trajectory of a particle in a circular Gaussian plane wave.
#
# _____________________________________________________________________________


import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# Step represents the step between trajectory points what we consider
# This enables to reduce the size of the array
step = 2

# List of relativistic pushers
pusher_list = ["boris","vay","higueracary"]

# We load successively the particle track associated to each pusher
for pusher in pusher_list:
  
  # Data from the Track diagnostic
  Track = S.TrackParticles("electron_" + pusher, axes=["x","px","py","pz"], timesteps=[0,6000]).get()

  # We extract x,px,py,pz from the first particle
  x = np.array(Track['x'][::step,0])
  px = np.array(Track['px'][::step,0])
  py = np.array(Track['py'][::step,0])
  pz = np.array(Track['pz'][::step,0])

  # We determine p_\perp, p_\parallel and gamma
  p_perp = np.sqrt(py**2 + pz**2)
  gamma = np.sqrt(1. + px**2 +  py**2 + pz**2)

  # Validation of p_\perp, p_\parallel and gamma
  Validate("Electron p_perp for pusher: " + pusher, p_perp, 1e-7 )
  Validate("Electron p_x for pusher: " + pusher, px, 1e-7 )
  Validate("Electron gamma for pusher: " + pusher, gamma, 1e-7 )


# Ey laser field after reflection on the rhs boundary
Ey = S.Field.Field0("Ey").getData(timestep=10000)[0]
Ey_energy = (Ey**2).sum()
Validate("Final Ey field energy: ", Ey_energy, 50.)