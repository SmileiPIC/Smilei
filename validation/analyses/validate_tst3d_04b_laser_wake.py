import os, re, numpy as np, math 
import happi

S = happi.Open(["./restart*"], verbose=False)


# Ey FIELD ALONG CENTRAL AXIS
Ey = S.Probe(0, "Ey")
Validate("Field Ey on central axis timestep 300" , Ey.getData(timestep= 300)[0][::4], 0.01)
Validate("Field Ey on central axis timestep 1000", Ey.getData(timestep=1000)[0][::4], 0.01)
