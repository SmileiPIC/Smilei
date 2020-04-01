# ______________________________________________________________________________
#
# This script dumps the tables for the multiphoton Breit-Wheeler process
# ______________________________________________________________________________

import h5py as h5py
import numpy as np

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
# Read the table integration_dt_dchi

if if_file_exist('./multiphoton_Breit_Wheeler_tables.h5'):
    
    f = h5py.File('./multiphoton_Breit_Wheeler_tables.h5', "r")

    if 'integration_dt_dchi' in f:

        dataset = f['integration_dt_dchi']
        
        size_photon_chi = dataset.attrs["size_photon_chi"]
        
        print(" Table integration_dt_dchi: ")
        
        message = ""
        
        for i in range(0,size_photon_chi,10):
            for j in range(i,min(i+10,size_photon_chi)):
                v = dataset[i]
                message += "{:}".format(v) + ", "
            message += "\n"
                
        print(message)

# ______________________________________________________________________________
# Read the table min particle chi for xi

if if_file_exist('./multiphoton_Breit_Wheeler_tables.h5'):

    f = h5py.File('./multiphoton_Breit_Wheeler_tables.h5', "r")
    
    if 'min_particle_chi_for_xi' in f:
        
        dataset = f['min_particle_chi_for_xi']
        
        size_photon_chi = dataset.attrs["size_photon_chi"]
        
        print(" Table min_particle_chi_for_xi: ")
        
        message = ""
        
        for i in range(0,size_photon_chi,10):
            for j in range(i,min(i+10,size_photon_chi)):
                v = dataset[i]
                message += "{:}".format(v) + ", "
            message += "\n"
                
        print(message)

# ______________________________________________________________________________
# Read the table xi

if if_file_exist('./multiphoton_Breit_Wheeler_tables.h5'):

    f = h5py.File('./multiphoton_Breit_Wheeler_tables.h5', "r")

    if 'xi' in f:
        
        xi = np.array(f['xi'])
        
        size_particle_chi = f['xi'].attrs["size_particle_chi"]
        size_photon_chi = f['xi'].attrs["size_photon_chi"]

        print(" Table xi: ")
        
        message = ""
        
        for i in range(size_photon_chi):
            for j in range(size_particle_chi):
                v = xi[i,j]
                message += "{:}".format(v) + ", "
            message += "\n"
                
        print(message)
