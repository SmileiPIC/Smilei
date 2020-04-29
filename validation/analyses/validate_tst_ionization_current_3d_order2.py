import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)

# check the time dependance of the electromagnetic energy
Uelm = np.array( S.Scalar('Uelm').getData() )

Uelm_max = np.max(Uelm)
Uelm_end = Uelm[-1]
print("- ratio of energy lost in ionizing the hydrogen plasma: ",(Uelm_end-Uelm_max)/Uelm_max)

Validate("Uelm_max", Uelm_max, 1.e-6*Uelm_max)
Validate("Uelm_end", Uelm_max, 1.e-6*Uelm_end)
Validate("Uelm",     Uelm,     1.e-2*Uelm_max)
