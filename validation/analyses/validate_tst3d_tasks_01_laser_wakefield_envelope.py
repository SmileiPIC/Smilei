import os, re, numpy as np, math, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

# COMPARE THE Envelope FIELD, absolute value
Env_E_abs = S.Field.Field0.Env_E_abs(subset={"z":S.namelist.Main.grid_length[2]/2 }, timesteps=800.).getData()[0]
Validate("Env_E_abs field at iteration 800", Env_E_abs, 0.01)

# COMPARE THE Ex FIELD
Ex = S.Field.Field0.Ex(subset={"z":S.namelist.Main.grid_length[2]/2 }, timesteps=800.).getData()[0]
Validate("Ex field at iteration 800", Ex, 0.01)

# COMPARE THE Rho FIELD
Rho = S.Field.Field0.Rho(subset={"z":S.namelist.Main.grid_length[2]/2 }, timesteps=800.).getData()[0]
Validate("Rho field at iteration 800", Rho, 0.01)


# TEST THE GRID PARAMETERS
with h5py.File("./restart000/Fields0.h5", "r") as f:
        dt = f["data/0000000000"].attrs["dt"]
        dx = f["data/0000000000/Ex"].attrs["gridSpacing"]
        patchSize = f["data/0000000000"].attrs["patchSize"]
Validate("Value of the timestep" , dt, 1e-6)
Validate("Value of the grid step", dx, 1e-6)
Validate("Patch size", patchSize)                                            
