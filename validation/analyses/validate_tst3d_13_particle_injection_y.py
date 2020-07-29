import os, re, numpy as np, math
import happi
from matplotlib.pyplot import *
from matplotlib.colors import LogNorm
import yt

S = happi.Open(["./restart*"], verbose=False)

# Get some Smilei parameters
timestep = S.namelist.Main.timestep
simulation_time = S.namelist.Main.simulation_time

diag_every = int(simulation_time / timestep)

species_list = ["eon1", "pon1", "eon2", "pon2"]

Scalar = {}

relative_error = 0.05

for species in species_list:
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

Scalar["time"] = np.array(S.Scalar("Ntot_eon1").get()["times"])*timestep

# Density binning _________________________________________________________________

print(" 2) Density particle binning diagnostic")

x_axis = {}
y_axis = {}
z_axis = {}
density_binning_final = {}

for i,species in enumerate(species_list):
    
    x_axis[species] = np.array(S.ParticleBinning(diagNumber=i,timesteps=0).get()["x"])
    y_axis[species] = np.array(S.ParticleBinning(diagNumber=i,timesteps=0).get()["y"])
    z_axis[species] = np.array(S.ParticleBinning(diagNumber=i,timesteps=0).get()["z"])
    density_binning_final[species] = np.array(S.ParticleBinning(diagNumber=i,timesteps=diag_every).getData()[0])

    average_density = np.mean(density_binning_final[species][density_binning_final[species]>0])
    
    Validate("Average density for {}".format(species_list[i]),average_density,0.01)
    
    print(" Average density for {}: {}".format(species_list[i], average_density))

# Energy _________________________________________________________________

print(" 3) Gamma spectrum particle binning diagnostic")

initial_gamma_spectrum = {}
final_gamma_spectrum = {}
gamma_bins = {}

sum_initial_gamma_spectrum = {}
sum_final_gamma_spectrum = {}

for i,species in enumerate(species_list):
    
    gamma_bins[species] = np.array(S.ParticleBinning(diagNumber=i+4,timesteps=0).get()["gamma"])
    initial_gamma_spectrum[species] = np.array(S.ParticleBinning(diagNumber=i+4,timesteps=0).getData()[0])
    final_gamma_spectrum[species] = np.array(S.ParticleBinning(diagNumber=i+4,timesteps=diag_every).getData()[0])
    
    sum_final_gamma_spectrum[species] = np.sum(final_gamma_spectrum[species])

    Validate("Sum of the gamma spectrum for {}".format(species_list[i]) , sum_final_gamma_spectrum[species], sum_final_gamma_spectrum[species]*relative_error)
    
    sum_initial_gamma_spectrum[species] = np.sum(initial_gamma_spectrum[species])
    
    error = (np.abs(np.subtract(initial_gamma_spectrum[species] / sum_initial_gamma_spectrum[species] , initial_gamma_spectrum[species] / sum_final_gamma_spectrum[species])))
    
    print(' Gamma spectrum max error for {}: {}'.format(species, np.max(error)))
    
    Validate("Gamma spectrum max error for {}".format(species) , np.max(error), relative_error)

# Figures __________________________________________________________________

if False:

    fig0 = figure(figsize=(10, 8))
    gs = GridSpec(3, 2)
    ax0 = subplot(gs[0,:])
    ax1 = subplot(gs[1,:])
    ax2 = subplot(gs[2,:])

    for species in species_list:
        name = "Ntot_{}".format(species)
        ax0.plot(Scalar["time"],Scalar[name],label=species)

    ax0.legend(loc="best")

    for species in species_list:
        name = "Ukin_{}".format(species)
        ax1.plot(Scalar["time"],Scalar[name])

    ax1.set_yscale("log")

    for species in species_list:
        name = "Dens_{}".format(species)
        ax2.plot(Scalar["time"],Scalar[name])

    fig1 = figure(figsize=(16, 8))
    gs = GridSpec(4,4)
    ax = {}
    ax["eon1"] = subplot(gs[0:2,0:2])
    ax["pon1"] = subplot(gs[0:2,2:4])
    ax["eon2"] = subplot(gs[2:4,0:2])
    ax["pon2"] = subplot(gs[2:4,2:4])

    for species in species_list:
        ax[species].plot(gamma_bins[species], initial_gamma_spectrum[species]/sum_initial_gamma_spectrum[species], color="C0")

    for species in species_list:
        ax[species].plot(gamma_bins[species], final_gamma_spectrum[species]/sum_final_gamma_spectrum[species], color="C1")

    fig2 = figure(figsize=(16, 8))
    gs = GridSpec(4,4)
    ax = {}
    ax["eon1"] = subplot(gs[0:2,0:2])
    ax["pon1"] = subplot(gs[0:2,2:4])
    ax["eon2"] = subplot(gs[2:4,0:2])
    ax["pon2"] = subplot(gs[2:4,2:4])

    for species in species_list:
        im = ax[species].pcolormesh(x_axis[species],y_axis[species],np.mean(density_binning_final[species].T,axis=0),
                        cmap=get_cmap('jet'),
                        shading='none')

        im.set_norm(LogNorm())
        cb = colorbar(im,ax=ax[species])
        ax[species].set_xlabel(r'$x$')
        ax[species].set_ylabel(r'$y$')
        t = ax[species].set_title(r'{} ({})'.format(species,diag_every))
        t.set_y(1.02)
        
    fig2.tight_layout()
    
    # data = dict(density = (density_binning_final[species], "g/cm**3"))
    # shape = density_binning_final[species].shape
    # bbox = np.array([[x_axis[species][0], x_axis[species][-1]], [-1.5, 1.5], [-1.5, 1.5]])
    # ds = yt.load_uniform_grid(data, shape, length_unit="Mpc", bbox=bbox)
    # yt.interactive_render(ds)
    # sc = yt.create_scene(ds)
    #
    # sc.camera.set_width(ds.quan(20, 'kpc'))
    # source = sc.sources['source_00']
    #
    # tf = yt.ColorTransferFunction((-28, -24))
    # tf.add_layers(4, w=0.01)
    #
    # source.set_transfer_function(tf)
    #
    # sc.show()
