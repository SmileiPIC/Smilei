import os, re, numpy as np, math
from scipy.interpolate import interp1d as interp
import happi

S = happi.Open(["./restart*"], verbose=False)

species_list = ["electron_cartesian","proton_cartesian",
				"electron_spherical_lin","proton_spherical_lin",
				"electron_spherical_log","proton_spherical_log"]

relative_error = 0.05

Scalar = {}

for ispecies,species in enumerate(species_list):
    name = "Ntot_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

    name = "Ukin_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

    name = "Dens_{}".format(species)
    Scalar[name] = np.array(S.Scalar(name).getData())
    for index,value in enumerate(Scalar[name]):
        Validate("Scalar {}[{}]".format(name,index) , value, value*relative_error)

# Energy _________________________________________________________________

relative_error = 0.05

cases = ["cartesian","spherical_lin","spherical_log"]

for icase,case in enumerate(cases):

    spectrum = S.ParticleBinning(diagNumber=icase,timesteps=50)
    data = np.array(spectrum.getData()[0])
    gamma = np.array(spectrum.get()["gamma"])
    integration = np.sum(data*gamma)*(gamma[1] - gamma[0])

    Validate("Integrated positron gamma spectrum for {}: ".format(case),integration,integration*relative_error)
    print("Integrated positron gamma spectrum for {}: {}.".format(case,integration))

    # for index,value in enumerate(data):
    #    Validate("Proton spectrum for {} at {}: ".format(case,index) , value, value*relative_error)

for icase,case in enumerate(cases):

    spectrum = S.ParticleBinning(diagNumber=icase+3,timesteps=50)
    data = np.array(spectrum.getData()[0])
    gamma = np.array(spectrum.get()["gamma"])
    integration = np.sum(data*gamma)*(gamma[1] - gamma[0])

    Validate("Integrated electron gamma spectrum for {}: ".format(case),integration,integration*relative_error)
    print("Integrated electron gamma spectrum for {}: {}.".format(case,integration))

    # for index,value in enumerate(data):
    #    Validate("Electron spectrum for {} at {}: ".format(case,index) , value, value*relative_error)
