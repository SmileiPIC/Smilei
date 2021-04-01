import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE Bz FIELD
Bz = S.Field.Field0.Bz(timesteps=1800, subset={"x":[50,180,8], "y":[0,100,8]}).getData()[0]
Validate("Bz field at iteration 1800", Bz, 0.01)
