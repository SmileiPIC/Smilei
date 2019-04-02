import os, re, numpy as np, math, h5py
import happi

# Read the radiation tables
f = h5py.File('./restart000/radiation_tables.h5','r')

table_h = f["./h"]
table_h_attrs = [table_h.attrs["chipa_min"],table_h.attrs["chipa_max"],table_h.attrs["chipa_dim"]]
table_integfochi = f["./integfochi"]
table_xip_chiphmin = f["./xip_chiphmin"]
table_xip = f["./xip"]

Validate("Table H for Niel", table_h.value, 0.01)
Validate("Table H attributes", table_h_attrs, 0.01)
Validate("Table integration f / chi (integfochi)", table_integfochi.value, 0.01)
Validate("Table photon chi min (xip_chiphmin)", table_xip_chiphmin.value, 0.01)
Validate("Table xip", table_xip.value, 0.01)

# Read the multiphoton Breit-Wheeler tables
f_bw = h5py.File('./restart000/multiphoton_Breit_Wheeler_tables.h5','r')

table_h_bw = f_bw["./h"]
table_xip_chipamin = f_bw["./xip_chipamin"]
table_xip_bw = f_bw["./xip"]

Validate("Table H for for BW", table_h_bw.value, 0.01)
Validate("Table particle chi min (xip_chipamin)", table_xip_chipamin.value, 0.01)
Validate("Table xip for BW", table_xip_bw.value, 0.01)
