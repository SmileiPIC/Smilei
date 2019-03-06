import os, re, numpy as np, math, glob
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=1600).getData()[0][::10,:]
Validate("Ey field at iteration 1600", Ey, 0.1)

# VERIFY SPECIES PROBE
Validate("Probe Rho_electron == Rho",
	all(S.Probe(0,'Rho_electron',timesteps=1800).getData()[0] == S.Probe(0,'Rho',timesteps=1800).getData()[0])
)

# CHECK THE LOAD BALANCING
txt = ""
restarts = glob.glob("restart*")
for folder in restarts:
	with open(folder+"/patch_load.txt") as f:
		txt += f.read()
patch_count0 = re.findall(r"patch_count\[0\] = (\d+)",txt)
patch_count1 = re.findall(r"patch_count\[1\] = (\d+)",txt)
initial_balance = [int(patch_count0[0] ), int(patch_count1[0] )]
final_balance   = [int(patch_count0[-1]), int(patch_count1[-1])]
Validate("Initial load balance", initial_balance, 1)
Validate("Final load balance", final_balance, 1)

# SCALARS RELATED TO BOUNDARIES AND MOVING WINDOW
Validate("Scalar Ukin_bnd"    , S.Scalar.Ukin_bnd    ().getData(), 0.0001)
Validate("Scalar Uelm_bnd"    , S.Scalar.Uelm_bnd    ().getData(), 1000. )
Validate("Scalar Ukin_out_mvw", S.Scalar.Ukin_out_mvw().getData(), 1.    )
Validate("Scalar Ukin_inj_mvw", S.Scalar.Ukin_inj_mvw().getData(), 2.e-7 )
Validate("Scalar Uelm_out_mvw", S.Scalar.Uelm_out_mvw().getData(), 1.    )
Validate("Scalar Uelm_inj_mvw", S.Scalar.Uelm_inj_mvw().getData(), 1.    )
