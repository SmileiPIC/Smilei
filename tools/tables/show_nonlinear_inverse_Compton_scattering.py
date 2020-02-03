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

if if_file_exist('./radiation_tables.h5'):
    
    f = h5py.File('./radiation_tables.h5', "r")

    if 'integfochi' in f:

        fig = figure(figsize=(10, 6))
        gs = GridSpec(2, 2)
        ax0 = subplot(gs[:,:])

        dataset = f['integfochi']

        size_particle_chi = dataset.attrs["size_particle_chi"]
        chipa_min = dataset.attrs["min_particle_chi"]
        chipa_max = dataset.attrs["max_particle_chi"]

        print("")
        print("Table 'integfochi'")
        print("Size: {}".format(size_particle_chi))
        print("Min particle chi: {}".format(chipa_min))
        print("Max particle chi: {}".format(chipa_max))

        chi = np.logspace(np.log10(chipa_min),np.log10(chipa_max),size_particle_chi)

        ax0.plot(chi,dataset,lw=2,label=r'$\int{F/\chi d\chi}$')

        ax0.set_xscale('log')
        ax0.set_xlabel(r'$\chi_\pm$')
        ax0.set_ylabel(r'$\int{F/\chi d\chi}$')

        fig.tight_layout()

# ______________________________________________________________________________
# Read the table min photon chi for xi

if if_file_exist('./radiation_tables.h5'):

    f = h5py.File('./radiation_tables.h5', "r")
    
    if 'min_photon_chi_for_xi' in f:
        
        fig = figure(figsize=(10, 6))
        gs = GridSpec(2, 2)
        ax0 = subplot(gs[:,:])
        
        dataset = f['min_photon_chi_for_xi']
        
        size_particle_chi = dataset.attrs["size_particle_chi"]
        min_particle_chi = dataset.attrs["min_particle_chi"]
        max_particle_chi = dataset.attrs["max_particle_chi"]

        print("")
        print("Table 'min_photon_chi_for_xi'")
        print("Size: {}".format(size_particle_chi))
        print("Min particle chi: {}".format(min_particle_chi))
        print("Max particle chi: {}".format(max_particle_chi))

        chipa = np.logspace(np.log10(min_particle_chi),np.log10(max_particle_chi),size_particle_chi)

        ax0.plot(chipa, np.power(10,dataset))

        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.set_xlabel(r'$\chi_\pm$')
        ax0.set_ylabel(r'$\chi_{\gamma,min}$')

        fig.tight_layout()

# ______________________________________________________________________________
# Read the table xip

if if_file_exist('./radiation_tables.h5'):

    f = h5py.File('./radiation_tables.h5', "r")

    if 'xi' in f:
        
        fig = figure(figsize=(8, 6))
        gs = GridSpec(2, 2)
        ax0 = subplot(gs[:,:])
        
        xi = np.array(f['xi']).T
        
        # chipa_dim = f['xip'].attrs["chipa_dim"]
        # chipa_min = f['xip'].attrs["chipa_min"]
        # chipa_max = f['xip'].attrs["chipa_max"]
        # chiph_dim = f['xip'].attrs["chiph_dim"]
        size_particle_chi = f['xi'].attrs["size_particle_chi"]
        chipa_min = f['xi'].attrs["min_particle_chi"]
        chipa_max = f['xi'].attrs["max_particle_chi"]
        chiph_dim = f['xi'].attrs["size_photon_chi"]

        print("")
        print("Table 'xi'")
        print("Size particle chi axis: {}".format(size_particle_chi))
        print("Size photon chi axis: {}".format(chiph_dim))
        print("Min particle chi: {}".format(chipa_min))
        print("Max particle chi: {}".format(chipa_max))

        chiph = np.linspace(1,chiph_dim+1,chiph_dim+1)
        chipa = np.logspace(np.log10(chipa_min),np.log10(chipa_max),size_particle_chi)

        im = ax0.pcolormesh(chipa,chiph,xi,
                        cmap=get_cmap('plasma'),
                        shading='none')

        ax0.set_xscale('log')
        #ax2.set_yscale('log')
        ax0.set_xlabel(r'$\chi_\pm$')
        ax0.set_ylabel(r'index from $\chi_{\gamma,min}$ to $\chi_\pm$')
        ax0.set_ylim([1,chiph[-1]])

        im.set_norm(LogNorm())
        cb = colorbar(im,ax=ax0)

        fig.tight_layout()

# ______________________________________________________________________________
# Read the table h

if if_file_exist('./radiation_tables.h5'):
    
    f = h5py.File('./radiation_tables.h5', "r")

    if 'h' in f:

        fig = figure(figsize=(8, 6))
        gs = GridSpec(2, 2)
        ax0 = subplot(gs[:,:])

        dataset = f['h']

        size_particle_chi = dataset.attrs["size_particle_chi"]
        chipa_min = dataset.attrs["min_particle_chi"]
        chipa_max = dataset.attrs["max_particle_chi"]

        print("")
        print("Table 'h'")
        print("Size: {}".format(size_particle_chi))
        print("Min particle chi: {}".format(chipa_min))
        print("Max particle chi: {}".format(chipa_max))

        chi = np.logspace(np.log10(chipa_min),np.log10(chipa_max),size_particle_chi)

        ax0.plot(chi,dataset,lw=2,label=r'$H$')

        ax0.set_xscale('log')
        ax0.set_yscale('log')
        ax0.set_xlabel(r'$\chi_\pm$')
        ax0.set_ylabel(r'$H$')
        
        fig.tight_layout()

show()
