import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE Bz FIELD
Bz = S.Field.Field0.Bz(average={"z":S.namelist.Main.grid_length[2]*0.5}, subset={"x":[20,50,2], "y":[10,35,2]}, timesteps=350).getData()[0]
Validate("Bz field at iteration 350", Bz, 0.01)
