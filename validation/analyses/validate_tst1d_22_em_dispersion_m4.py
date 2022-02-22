import os, re, numpy as np, h5py
import happi

S = happi.Open(["./restart*"], verbose=False)

Validate("List of fields in Field0", S.fieldInfo(0)["fields"] )

timesteps = list(S.Field.Field0("Ez").getAvailableTimesteps())
Validate("List of timesteps in Field0", timesteps )

Ex = S.Field.Field0("Ex", timesteps=timesteps[-1]).getData()[0]
Ey = S.Field.Field0("Ey", timesteps=timesteps[-1]).getData()[0]
Ez = S.Field.Field0("Ez", timesteps=timesteps[-1]).getData()[0]
Erms = np.mean(Ex**2+Ey**2+Ez**2)**0.5 
Validate("Last Ex profile in Field0", Ex, Erms*1e-7 )
Validate("Last Ey profile in Field0", Ey, Erms*1e-7 )
Validate("Last Ez profile in Field0", Ez, Erms*1e-7 )

Bx = S.Field.Field0("Bx", timesteps=timesteps[-1]).getData()[0]
By = S.Field.Field0("By", timesteps=timesteps[-1]).getData()[0]
Bz = S.Field.Field0("Bz", timesteps=timesteps[-1]).getData()[0]
Brms = np.mean(Bx**2+By**2+Bz**2)**0.5 
Validate("Last Bx profile in Field0", Bx, Brms*1e-7 )
Validate("Last By profile in Field0", By, Brms*1e-7 )
Validate("Last Bz profile in Field0", Bz, Brms*1e-7 )
