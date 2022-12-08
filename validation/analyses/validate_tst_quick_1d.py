import os, re, numpy as np, math
import happi

S = happi.Open(["./restart*"], verbose=False)

# Get some Smilei parameters
timestep = S.namelist.Main.timestep
simulation_time = S.namelist.Main.simulation_time

# ------------------------------------------------------
# Scalar diags
# ------------------------------------------------------

print(" ------------------------------")
print(" Scalars")
print(" ------------------------------")

species_list = ["eon1", "pon1", "eon2", "pon2"]

relative_error = 0.05
scalar_list = ["Ukin", "Utot", "Uelm"]

for species in species_list:
    scalar_list.append("Ntot_{}".format(species))
    scalar_list.append("Ukin_{}".format(species))
    scalar_list.append("Dens_{}".format(species))


Scalar = {}

for scalar_name in scalar_list:
    Scalar[scalar_name] = np.array(S.Scalar(scalar_name).getData())
    print(" Validate {}".format(scalar_name))
    for index,value in enumerate(Scalar[scalar_name]):
        Validate("Scalar {}[{}]".format(scalar_name, index) , value, value*relative_error)

print("")

# ------------------------------------------------------
# Binning diags
# ------------------------------------------------------

print(" ------------------------------")
print(" Diags Binning")
print(" ------------------------------")

for i in range(4):
    
    particle_binning = S.ParticleBinning(diagNumber=i,timesteps=10)
    data = np.array(particle_binning.getData()[0])
    sum = np.sum(np.abs(data))

    print(" Valide sum of binning {}: {}".format(i,sum))

    Validate("Sum of diag {}".format(i) , sum, sum*relative_error)

print("")

# ------------------------------------------------------
# Field diags
# ------------------------------------------------------

print(" ------------------------------")
print(" Diags Field")
print(" ------------------------------")

fields = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Rho", "Jx", "Jy", "Jz"]

for species in species_list:
    fields.append("Jx_{}".format(species))	
    fields.append("Jy_{}".format(species))	
    fields.append("Jz_{}".format(species))	
    fields.append("Rho_{}".format(species))	

for field_name in fields:

    field = np.array(S.Field(0,field_name,timesteps=10).getData()[0])

    sum = np.sum(np.abs(field))

    print(" Valide sum of {}: {}".format(field_name,sum))

    Validate("Sum for field {}".format(field_name) , sum, sum*relative_error)

print("")

# ------------------------------------------------------
# Probe diags
# ------------------------------------------------------

print(" ------------------------------")
print(" Diags Probe")
print(" ------------------------------")

probe_list = ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "Rho", "Jx", "Jy", "Jz"]

for species in species_list:
    probe_list.append("Jx_{}".format(species))	
    probe_list.append("Jy_{}".format(species))	
    probe_list.append("Jz_{}".format(species))	
    probe_list.append("Rho_{}".format(species))	

for probe_name in probe_list:

    probe = np.array(S.Probe(0, probe_name,timesteps=10).getData()[0])

    sum = np.sum(np.abs(probe))

    print(" Valide sum of {}: {}".format(probe_name,sum))

    Validate("Sum for probe {}".format(probe_name) , sum, sum*relative_error)


print("")