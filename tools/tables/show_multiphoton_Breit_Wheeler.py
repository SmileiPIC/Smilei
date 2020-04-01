# ______________________________________________________________________________
#
# This script reads and plots the tables for synchroton radiation
# ______________________________________________________________________________


# ______________________________________________________________________________
# Importation

from matplotlib.pyplot import *
from matplotlib.colors import LogNorm
import h5py as h5py
import numpy as np
import sys

# ______________________________________________________________________________
# RCparams

rcParams['font.size'] = 20
rcParams['figure.facecolor'] = 'w'
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['axes.labelsize'] = 20

rcParams['xtick.major.size'] = 10
rcParams['ytick.major.size'] = 10

rcParams['xtick.minor.size'] = 5
rcParams['ytick.minor.size'] = 5

rcParams['axes.linewidth'] = 1.5

rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5

# ______________________________________________________________________________
# Checks

try:
    path = sys.argv[1]
except:
    print("\n Please, provide a path to the tables.\n")
    raise


# ______________________________________________________________________________
# Functions

def if_file_exist(filename):
  """
  """
  flag = True
  try:
    file = open(filename)
  except IOError:
    flag = False
    print ' No file : ', filename
  return flag
  

# ______________________________________________________________________________
# Read the table integfochi

if if_file_exist(path):
    
    f = h5py.File(path, "r")

    if 'integration_dt_dchi' in f:
    
        fig = figure(figsize=(8, 6))
        gs = GridSpec(2, 2)
        ax0 = subplot(gs[:,:])

        dataset = f['integration_dt_dchi']

        size_photon_chi = dataset.attrs["size_photon_chi"]
        chiph_min = dataset.attrs["min_photon_chi"]
        chiph_max = dataset.attrs["max_photon_chi"]

        print("Table T")
        print("Dimension: {}".format(size_photon_chi))
        print("Min particle chi: {}".format(chiph_min))
        print("Max particle chi: {}".format(chiph_max))

        chiph = np.logspace(np.log10(chiph_min),np.log10(chiph_max),size_photon_chi)

        ax0.plot(chiph,dataset,lw=2)

        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.set_xlabel(r'$\chi_\gamma$')
        ax0.set_ylabel(r'$T(\chi_\gamma)$')
        
        fig.tight_layout()

# ______________________________________________________________________________
# Read the table xip_chipamin
    
    if "min_particle_chi_for_xi" in f:
    
        fig = figure(figsize=(8, 6))
        gs = GridSpec(2, 2)
        ax0 = subplot(gs[:,:])
    
        xip_chipamin = f['min_particle_chi_for_xi']
        
        size_photon_chi = xip_chipamin.attrs["size_photon_chi"]
        chiph_min = xip_chipamin.attrs["min_photon_chi"]
        chiph_max = xip_chipamin.attrs["max_photon_chi"]

        print("")
        print("Table xi_chipamin")
        print("Dimension: {}".format(size_photon_chi))
        print("Min particle chi: {}".format(chiph_min))
        print("Max particle chi: {}".format(chiph_max))

        chiph = np.logspace(np.log10(chiph_min),np.log10(chiph_max),size_photon_chi)

        ax0.plot(chiph, np.power(10,xip_chipamin))

        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.set_xlabel(r'$\chi_\gamma$')
        ax0.set_ylabel(r'$\chi_{\pm,min}$')
        
        fig.tight_layout()


# ______________________________________________________________________________
# Read the table xip

    if "xi" in f:
        
        fig = figure(figsize=(8, 6))
        gs = GridSpec(2, 2)
        ax0 = subplot(gs[:,:])
        
        xip = np.array(f['xi']).T
        
        size_photon_chi = f['xi'].attrs["size_photon_chi"]
        chiph_min = f['xi'].attrs["min_photon_chi"]
        chiph_max = f['xi'].attrs["max_photon_chi"]
        chipa_dim = f['xi'].attrs["size_particle_chi"]

        print("")
        print("Table xi")
        print("Dimension: {}".format(size_photon_chi))
        print("Min particle chi: {}".format(chiph_min))
        print("Max particle chi: {}".format(chiph_max))

        chipa = np.linspace(1,chipa_dim+1,chipa_dim+1)
        chiph = np.logspace(np.log10(chiph_min),np.log10(chiph_max),size_photon_chi)
        
        im = ax0.pcolormesh(chiph,chipa,xip,
                        cmap=get_cmap('plasma'),
                        shading='none')

        ax0.set_xscale('log')
        #ax2.set_yscale('log')
        ax0.set_xlabel(r'$\chi_\gamma$')
        ax0.set_ylabel(r'index from $\chi_{\pm,min}$ to $\chi_\gamma$')
        ax0.set_ylim([1,chipa[-1]])

        im.set_norm(LogNorm())
        cb = colorbar(im,ax=ax0)

        fig.tight_layout()

show()
