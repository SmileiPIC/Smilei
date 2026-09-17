import os, re, numpy as np, h5py
import happi

#S = happi.Open(["./restart*"], verbose=False)
S = happi.Open(verbose=False)

Validate("List of fields in Probe0", S.probeInfo(0)["fields"] )

timesteps = list(S.Probe.Probe0("Ez").getAvailableTimesteps())
Validate("List of timesteps in Field0", timesteps )

Uelm = S.Scalar("Uelm").getData()
Validate("Uelm in the simulation domain", Uelm,0.05)

Ukin = S.Scalar("Ukin").getData()
Validate("Ukin in the simulation domain", Ukin,0.05)

#Ey = np.array(S.Probe.Probe0("Ey", timesteps=timesteps[-1],subset={"axis0":S.namelist.Lz/2}).getData())[0,::16,::8]
#Erms = np.mean(Ey**2)**0.5
#Validate("Ey profile in Field0 at Tsim", Ey, Erms*1e-7 )

#Bz = np.array(S.Probe.Probe0("Bz", timesteps=timesteps[-1],subset={"axis0":S.namelist.Lz/2}).getData())[0,::16,::8]
#Brms = np.mean(Bz**2)**0.5
#Validate("Bz profile in Field0 at Tsim", Bz, Brms*1e-7 )
