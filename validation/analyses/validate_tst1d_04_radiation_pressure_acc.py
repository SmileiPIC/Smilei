import os, re, numpy as np
import happi

S = happi.Open(["./restart*"], verbose=False)




# PHYSICS OF ION ACCELERATION
ion_profiles = S.Field.Field0.Rho_eon().getData()
ion_front = []
for p in ion_profiles: ion_front += [np.flatnonzero(p<-10)[0]]
ion_front = np.array(ion_front)
Validate("Ion front position vs time", ion_front, 100.)


# SCALARS RELATED TO GLOBAL ENERGIES
Validate("Scalar Utot"     , S.Scalar.Utot     ().getData(), 200.)
Validate("Scalar Ukin"     , S.Scalar.Ukin     ().getData(), 100.)
Validate("Scalar Uelm"     , S.Scalar.Uelm     ().getData(), 100.)
Validate("Scalar Uexp"     , S.Scalar.Uexp     ().getData(), 200.)
Validate("Scalar Ubal"     , S.Scalar.Ubal     ().getData(), 1.  )
Validate("Scalar Ubal_norm", S.Scalar.Ubal_norm().getData(), 0.005)
Validate("Scalar Uelm_Ex"  , S.Scalar.Uelm_Ex  ().getData(), 100.)
Validate("Scalar Uelm_Ey"  , S.Scalar.Uelm_Ey  ().getData(), 100.)
Validate("Scalar Uelm_Ez"  , S.Scalar.Uelm_Ez  ().getData(), 100.)
Validate("Scalar Uelm_Bx_m", S.Scalar.Uelm_Bx_m().getData(), 100.)
Validate("Scalar Uelm_By_m", S.Scalar.Uelm_By_m().getData(), 100.)
Validate("Scalar Uelm_Bz_m", S.Scalar.Uelm_Bz_m().getData(), 100.)

# SCALARS RELATED TO SPECIES
Validate("Scalar Ukin_ion" , S.Scalar.Ukin_ion ().getData(), 100.)

# OPERATION ON SCALARS
Validate("Scalar Ukin+Uelm", S.Scalar("Ukin+Uelm").getData(), 100.)

# 1D SCREEN DIAGS
for i,d in enumerate(S.namelist.DiagScreen):
	last_data = S.Screen(i, timesteps=21000).getData()[-1]
	if d.direction=="backward":
		precision = 0.4
	elif d.direction=="canceling":
		precision = 1.
	else:
		precision = 1.5
	Validate("Screen "+d.shape+" diag with "+d.direction+" direction", last_data, precision)
