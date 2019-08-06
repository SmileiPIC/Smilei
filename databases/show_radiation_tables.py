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
# Figures

fig = figure(figsize=(16, 6))
gs = GridSpec(2, 6)
ax0 = subplot(gs[:,0:2])
ax1 = subplot(gs[:,2:4])
ax2 = subplot(gs[:,4:6])

# ______________________________________________________________________________
# Read the table integfochi

if if_file_exist('./radiation_tables.h5'):
    
    f = h5py.File('./radiation_tables.h5', "r")

    dataset = f['integfochi']

    chipa_dim = dataset.attrs["chipa_dim"]
    chipa_min = dataset.attrs["chipa_min"]
    chipa_max = dataset.attrs["chipa_max"]

    print("Table integfochi")
    print("Dimension: {}".format(chipa_dim))
    print("Min particle chi: {}".format(chipa_min))
    print("Max particle chi: {}".format(chipa_max))

    chi = np.logspace(np.log10(chipa_min),np.log10(chipa_max),chipa_dim)

    ax0.plot(chi,dataset,lw=2,label=r'$\int{F/\chi d\chi}$')

    ax0.set_xscale('log')
    ax0.set_xlabel(r'$\chi_\pm$')
    ax0.set_ylabel(r'$\int{F/\chi d\chi}$')

# ______________________________________________________________________________
# Read the table xip_chiphmin

if if_file_exist('./radiation_tables.h5'):

    f = h5py.File('./radiation_tables.h5', "r")
    
    xip_chiphmin = f['xip_chiphmin']
    
    chipa_dim = f['xip_chiphmin'].attrs["chipa_dim"]
    chipa_min = f['xip_chiphmin'].attrs["chipa_min"]
    chipa_max = f['xip_chiphmin'].attrs["chipa_max"]

    print("")
    print("Table 'xip_chiphmin")
    print("Dimension: {}".format(chipa_dim))
    print("Min particle chi: {}".format(chipa_min))
    print("Max particle chi: {}".format(chipa_max))

    chipa = np.logspace(np.log10(chipa_min),np.log10(chipa_max),chipa_dim)

    ax1.plot(chipa, np.power(10,xip_chiphmin))

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$\chi_\pm$')
    ax1.set_ylabel(r'$\chi_{\gamma,min}$')

# ______________________________________________________________________________
# Read the table xip

if if_file_exist('./radiation_tables.h5'):

    f = h5py.File('./radiation_tables.h5', "r")

    xip = np.array(f['xip']).T
    
    chipa_dim = f['xip'].attrs["chipa_dim"]
    chipa_min = f['xip'].attrs["chipa_min"]
    chipa_max = f['xip'].attrs["chipa_max"]
    chiph_dim = f['xip'].attrs["chiph_dim"]

    print("")
    print("Table xi")
    print("Dimension particle chi: {}".format(chipa_dim))
    print("Min particle chi: {}".format(chipa_min))
    print("Max particle chi: {}".format(chipa_max))

    chiph = np.linspace(1,chiph_dim+1,chiph_dim+1)
    chipa = np.logspace(np.log10(chipa_min),np.log10(chipa_max),chipa_dim)

    im = ax2.pcolormesh(chipa,chiph,xip,
                    cmap=get_cmap('plasma'),
                    shading='none')

    ax2.set_xscale('log')
    #ax2.set_yscale('log')
    ax2.set_xlabel(r'$\chi_\pm$')
    ax2.set_ylabel(r'index from $\chi_{\gamma,min}$ to $\chi_\pm$')
    ax2.set_ylim([1,128])

    im.set_norm(LogNorm())
    cb = colorbar(im,ax=ax2)

# ______________________________________________________________________________

fig.tight_layout()

show()
