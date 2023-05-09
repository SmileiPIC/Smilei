import os, re, numpy as np, math 
from scipy.ndimage.filters import median_filter
import happi

S = happi.Open(["./restart*"], verbose=False)



ncel = [int(i) for i in S._ncels]

# SCALARS RELATED TO MIN/MAX OF FIELDS
for field in ["Ex", "Ey", "Ez", "Bx_m", "By_m", "Bz_m", "Jx", "Jy", "Jz", "Rho"]:
	if field in ["Ex", "Rho"]:
		precision = 0.25
	else:
		precision = 0.025
	Validate("Minimum of scalar "+field, S.Scalar(field+"Min").getData()[-1], precision)
	Validate("Maximum of scalar "+field, S.Scalar(field+"Max").getData()[-1], precision)
	
	MinLoc = np.unravel_index(int(S.Scalar(field+"MinCell").getData()[-1]), ncel)
	MaxLoc = np.unravel_index(int(S.Scalar(field+"MaxCell").getData()[-1]), ncel)
	Validate("Location of minimum of scalar "+field, MinLoc, 10)
	Validate("Location of maximum of scalar "+field, MaxLoc, 10)

# FIELD DIAGNOSTICS
fields     = ["Ex","Ey" ,"Ez" ,"Bx" ,"By" ,"Bz" ,"Bx_m","By_m","Bz_m","Jx" ,"Jy" ,"Jz" ,"Rho","Jx_test0","Jy_test0","Jz_test0","Rho_test0"]
precisions = [ 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 , 0.01 , 0.01 , 0.01, 0.01, 0.01, 0.1 , 3e-4     , 3e-4     , 0.001    , 0.1       ]
for field, precision in zip(fields, precisions) :
	F = S.Field(0, field, average={"z":"all"}, timesteps=40).getData()[-1]
	F = median_filter(F, 3)[::4,::4]
	Validate("Field "+field, F, precision)

# PROBE DIAGNOSTICS
fields     = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Jx", "Jy", "Jz", "Rho"]
precisions = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 , 0.3  ]
for field, precision in zip(fields, precisions):
	Validate("Probe field "+field, S.Probe(0, field).getData()[-1], precision)

# PARTICLE BINNING DIAGNOSTICS
precision = [0.01, 0.01, 0.01, 100., 100., 100., 100., 1., 5e7, 100., 100., 100., 100., 5e7, 2e-3, 0.01, 0.01]
for i, axis in enumerate(S.namelist.axes):
	name = axis[0] if type(axis[0]) is str else "user_function"
	Validate("Particle binning axis "+name, S.ParticleBinning(i, timesteps=40).getData()[-1], precision[i])
precision = [2e-3, 0.01, 0.01, 1e-7, 1e-7, 1e-7, 2e-12, 2e-7, 1e-7, 1e-7, 1e-7, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 0.01]
for j, deposited_quantity in enumerate(S.namelist.quantities):
	name = deposited_quantity if type(deposited_quantity) is str else "user_function"
	Validate("Particle binning deposited_quantity "+name, S.ParticleBinning(i+j+1, timesteps=40).getData()[-1], precision[j])

# CHECK THE RESULT OF keep_interpolated_fields in ParticleBinning
Validate("keep_interpolated_fields Ex", S.ParticleBinning(34, timesteps=45).getData()[0], 0.001)
Validate("keep_interpolated_fields Ey", S.ParticleBinning(35, timesteps=45).getData()[0], 0.001)
Validate("keep_interpolated_fields Ez", S.ParticleBinning(36, timesteps=45).getData()[0], 0.001)
Validate("keep_interpolated_fields Wx", S.ParticleBinning(37, timesteps=45).getData()[0], 1e-8)
Validate("keep_interpolated_fields Wy", S.ParticleBinning(38, timesteps=45).getData()[0], 1e-8)
Validate("keep_interpolated_fields Wz", S.ParticleBinning(39, timesteps=45).getData()[0], 1e-8)

# CHECK THE RESULT OF keep_interpolated_fields in ParticleBinning
d = S.TrackParticles("test2", axes=["Ex", "Ey", "Ez", "Wx", "Wy", "Wz"], timesteps=45).getData()
Validate("Tracked interpolated fields Ex", d["Ex"].squeeze(), 0.001)
Validate("Tracked interpolated fields Ey", d["Ey"].squeeze(), 0.001)
Validate("Tracked interpolated fields Ez", d["Ez"].squeeze(), 0.001)
Validate("Tracked interpolated fields Wx", d["Wx"].squeeze(), 1e-8)
Validate("Tracked interpolated fields Wy", d["Wy"].squeeze(), 1e-8)
Validate("Tracked interpolated fields Wz", d["Wz"].squeeze(), 1e-8)
