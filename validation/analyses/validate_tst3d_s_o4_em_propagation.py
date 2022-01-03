import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(average={"z":S.namelist.Main.grid_length[2]*0.3}, timesteps=240).getData()[0]
Validate("Ey field at iteration 240", Ey, 0.01)

# SUBGRID OF FIELD DIAG
Ey  = S.Field.Field0.Ey(timesteps=240).getData()[0]
Ey1 = S.Field.Field1.Ey(timesteps=240).getData()[0]
subgrid = S.namelist.DiagFields[1].subgrid
Validate("Field subgrid works", (Ey1==Ey[subgrid]).all())

# TIME-AVERAGED FIELD DIAG
Ey = S.Field.Field2("Ey", subset={"x":30,"y":[20,50,4],"z":[0,30,4]}, timesteps=209).getData()[0]
Validate("Ey profile in Field2", Ey, 1e-7 )

# 0-D PROBE IN 3D
Ey = S.Probe.Probe0.Ey().getData()
Validate("0-D probe Ey vs time", Ey, 0.01)
# 1-D PROBE IN 3D
Ey = S.Probe.Probe1.Ey(timesteps=240).getData()[0]
Validate("1-D probe Ey at last iteration", Ey, 0.01)
# 2-D PROBE IN 3D
Ey = S.Probe.Probe2.Ey(timesteps=240).getData()[0]
Validate("2-D probe Ey at last iteration", Ey, 0.01)
# 3-D PROBE IN 3D
Ey = S.Probe.Probe3.Ey(timesteps=240).getData()[0]
Validate("3-D probe Ey at last iteration", Ey, 0.01)

# TEST THAT Ubal_norm STAYS OK
max_ubal_norm = np.max( np.abs(S.Scalar.Ubal_norm().getData()) )
Validate("Max Ubal_norm is below 10%", max_ubal_norm<.1 )

# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5", "r") as f:
	dt = f["data/0000000000"].attrs["dt"]
	dx = f["data/0000000000/Ex"].attrs["gridSpacing"]
	patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)

# COMPARE THE FIELDS OF TRACKED PARTICLES
#attributes = ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]
#data = S.TrackParticles.eon(axes=attributes, select="any(t==0, Id==100)").getData()
#for f in attributes:
#	Validate("Field "+f+" of tracked electrons", data[f].flatten(), 5e-4)

#attributes = ["x", "y", "z"]
#data = S.TrackParticles.eon(axes=attributes, select="any(t==0, Id==100)").getData()
#for f in attributes:
#	Validate("trajectory in "+f+" of tracked electrons", data[f].flatten(), 5e-4)

#attributes = ["px", "py", "pz"]
#data = S.TrackParticles.eon(axes=attributes, select="any(t==0, Id==100)").getData()
#for f in attributes:
#	Validate("momentum in "+f+" of tracked electrons", data[f].flatten(), 5e-4)

axes=["x","y","z","px","py","pz","Ex",'Ey','Ez','Bx','By','Bz']

track_data = S.TrackParticles.eon(axes=axes,sort=True).getData()
for ip in [1,3,5,9]:
    for key, values in track_data.items():
        if (key in axes):
            array = track_data[key].T[ip]
            #values = np.array(track_data[key].T[ip])
            for i in range(len(array)):
                if np.isnan(array[i]):
                    array[i] = 0.
            Validate("Particle {} - attribute {}".format(ip,key),array, 5e-4)
