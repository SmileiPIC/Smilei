import os, re, numpy as np, math, h5py
import happi

# Read the radiation tables
f = h5py.File('radiation.h5','r')

table_h = f["h"]
table_intefochi = f["intefochi"]
table_chiphmin = f["chiphmin"]
table_xip = f["xip"]

Validate("Table H for Niel", table_h.value, 0.01)
Validate("Table integration f / chi (intefochi)", table_intefochi.value, 0.01)
Validate("Table photon chi min (chiph min)", table_chiphmin.value, 0.01)
Validate("Table integration f / chi (intefochi)", table_xip.value, 0.01)

# Read the multiphoton Breit-Wheeler tables
f = h5py.File('multiphoton_breit_wheeler_tables.h5','r')




# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5") as f:
	dt = f["data/0000000000"].attrs["dt"]
	dx = f["data/0000000000/Ex"].attrs["gridSpacing"]
	patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)
