# ______________________________________________________________________________
#
# Validation script for the radiation in the collision of
# a GeV electron bunch with a counter-propagatin circularly polarized wave
#
# In this tests case, an electron bunch is initialized per radiation
# loss models at the same positions with an energy of 1 GeV near the right
# boundary of the box. They propagate to the left of the box where a circularly
# polarized laser plane wave is injected. This wave has an hyper-guassian
# profile of wavelength \lambda.
#
# Validation:
# - Discontinuous radiation loss
# - Continuous radiation loss
# ______________________________________________________________________________

import os, re, numpy as np, h5py
from Smilei import *

S = Smilei(".", verbose=False)

# List of relativistic pushers
radiation_list = ["cont","disc"]

ukin_dict = {}
urad_dict = {}

# We load successively the particle track associated
# to each radiation algorithm
for radiation in radiation_list:

  # Read scalar diagnostics
  ScalarUkinDiag = S.Scalar("Ukin_electron_"+radiation).get()
  ukin = np.array(ScalarUkinDiag["data"])
  times = np.array(ScalarUkinDiag["times"])

  ScalarUradDiag = S.Scalar("Urad_electron_"+radiation).get()
  urad = np.array(ScalarUradDiag["data"])

  utot = ukin+urad

  print ' Electron_'+radiation+':'
  print ' Final kinetic energy: ', ukin[-1]
  print ' Maximum radiated energy: ', urad.max()
  print

  # Validation of the kinetic energy
  Validate("Kinetic energy evolution: ", ukin/ukin[0], 0.5e-1 )

  # Validation of the radiated energy
  Validate("Radiated energy evolution: " , urad/ukin[0], 0.5e-1 )

  # Validation of the total energy
  Validate("Total energy error (max - min)/uref: " ,
           (utot.max() - utot.min())/utot[0], 1e-2)

  ukin_dict[radiation] = ukin
  urad_dict[radiation] = urad

# ______________________________________________________________________________
# Comparision continuous and discontinuous methods

urad_rel_err = abs(urad_dict["disc"] - urad_dict["cont"]) / urad_dict["cont"].max()
ukin_rel_err = abs(ukin_dict["disc"] - ukin_dict["cont"]) / ukin_dict["cont"][0]

print ' Comparision continuous/discontinuous methods'
print ' Maximum relative error kinetic energy',ukin_rel_err.max()
print ' Maximum relative error radiative energy',urad_rel_err.max()

# Validation difference between continuous and discontinuous methods
Validate("Relative error on the kinetic energy / ukin at t=0: " , ukin_rel_err.max(), 0.03 )
Validate("Relative error on the radiative energy / urad max " , urad_rel_err.max(), 0.03 )
