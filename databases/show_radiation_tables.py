# ______________________________________________________________________________
#
# This script reads and plots the tables for synchroton radiation
# ______________________________________________________________________________


# ______________________________________________________________________________
# Importation

from matplotlib.pyplot import *
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

fig = figure(figsize=(12, 6))
gs = GridSpec(2, 2)
ax0 = subplot(gs[:,0:1])
ax1 = subplot(gs[:,1:2])

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
# Read the table xip

if if_file_exist('./radiation_tables.h5'):

    f = h5py.File('./radiation_tables.h5', "r")
    
    dataset = f['xip_chiphmin']
    
    chipa_dim = dataset.attrs["chipa_dim"]
    chipa_min = dataset.attrs["chipa_min"]
    chipa_max = dataset.attrs["chipa_max"]

fig.tight_layout()

show()
