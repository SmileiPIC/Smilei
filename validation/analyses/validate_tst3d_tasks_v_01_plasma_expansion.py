import os, re, numpy as np
import happi

S = happi.Open(["./restart*"], verbose=False)

max_ubal = np.max( np.abs(S.Scalar.Ubal().getData()) )
Validate("Max Ubal is below 2%", max_ubal<0.02 )


