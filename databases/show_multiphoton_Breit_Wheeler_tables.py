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

if if_file_exist('./multiphoton_Breit_Wheeler_tables.h5'):
    
    f = h5py.File('./multiphoton_Breit_Wheeler_tables.h5', "r")

    dataset = f['h']

    chiph_dim = dataset.attrs["chiph_dim"]
    chiph_min = dataset.attrs["chiph_min"]
    chiph_max = dataset.attrs["chiph_max"]

    print("Table dN/dt")
    print("Dimension: {}".format(chiph_dim))
    print("Min particle chi: {}".format(chiph_min))
    print("Max particle chi: {}".format(chiph_max))

    chiph = np.logspace(np.log10(chiph_min),np.log10(chiph_max),chiph_dim)

    ax0.plot(chiph,dataset,lw=2)

    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xlabel(r'$\chi_\gamma$')
    ax0.set_ylabel(r'$dN_\pm/dt$')

# ______________________________________________________________________________
# Read the table xip_chipamin

if if_file_exist('./multiphoton_Breit_Wheeler_tables.h5'):

    f = h5py.File('./multiphoton_Breit_Wheeler_tables.h5', "r")
    
    xip_chipamin = f['xip_chipamin']
    
    chiph_dim = xip_chipamin.attrs["chiph_dim"]
    chiph_min = xip_chipamin.attrs["chiph_min"]
    chiph_max = xip_chipamin.attrs["chiph_max"]

    print("")
    print("Table xi_chipamin")
    print("Dimension: {}".format(chiph_dim))
    print("Min particle chi: {}".format(chiph_min))
    print("Max particle chi: {}".format(chiph_max))

    chiph = np.logspace(np.log10(chiph_min),np.log10(chiph_max),chiph_dim)

    ax1.plot(chiph, np.power(10,xip_chipamin))

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$\chi_\gamma$')
    ax1.set_ylabel(r'$\chi_{\pm,min}$')

# ______________________________________________________________________________
# Read the table xip

if if_file_exist('./multiphoton_Breit_Wheeler_tables.h5'):

    f = h5py.File('./multiphoton_Breit_Wheeler_tables.h5', "r")
    
    xip = np.array(f['xip']).T
    
    chiph_dim = f['xip'].attrs["chiph_dim"]
    chiph_min = f['xip'].attrs["chiph_min"]
    chiph_max = f['xip'].attrs["chiph_max"]
    chipa_dim = f['xip'].attrs["chipa_dim"]

    print("")
    print("Table xi")
    print("Dimension: {}".format(chiph_dim))
    print("Min particle chi: {}".format(chiph_min))
    print("Max particle chi: {}".format(chiph_max))

    chipa = np.linspace(1,chipa_dim+1,chipa_dim+1)
    chiph = np.logspace(np.log10(chiph_min),np.log10(chiph_max),chiph_dim)
    
    im = ax2.pcolormesh(chiph,chipa,xip,
                    cmap=get_cmap('plasma'),
                    shading='none')

    ax2.set_xscale('log')
    #ax2.set_yscale('log')
    ax2.set_xlabel(r'$\chi_\gamma$')
    ax2.set_ylabel(r'index from $\chi_{\pm,min}$ to $\chi_\gamma$')
    ax2.set_ylim([0,127])

    im.set_norm(LogNorm())
    cb = colorbar(im,ax=ax2)

# ______________________________________________________________________________


fig.tight_layout()

show()
