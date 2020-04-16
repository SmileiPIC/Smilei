# ______________________________________________________________________________
#
# This script dumps the tables for synchroton radiation
# ______________________________________________________________________________

import h5py as h5py
import numpy as np
import sys
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
    print(' No file: {}'.format(filename))
  return flag

# ______________________________________________________________________________
# Dump the table integfochi

if if_file_exist(path):
    
    f = h5py.File(path, "r")

    if 'integfochi' in f:

        dataset = f['integfochi']
        
        size_particle_chi = dataset.attrs["size_particle_chi"]
        
        print(" Table integfochi: ")
        
        message = ""
        
        for i in range(0,size_particle_chi,10):
            for j in range(i,min(i+10,size_particle_chi)):
                v = dataset[i]
                message += "{:}".format(v) + ", "
            message += "\n"
                
        print(message)

# ______________________________________________________________________________
# Read the table h

if if_file_exist(path):
    
    f = h5py.File(path, "r")

    if 'h' in f:

        dataset = f['h']
        
        size_particle_chi = dataset.attrs["size_particle_chi"]
        
        print(" Table h for Niel: ")
        
        message = ""
        
        for i in range(0,size_particle_chi,10):
            for j in range(i,min(i+10,size_particle_chi)):
                v = dataset[i]
                message += "{:}".format(v) + ", "
            message += "\n"
                
        print(message)

# ______________________________________________________________________________
# Read the table min photon chi for xi

if if_file_exist(path):

    f = h5py.File(path, "r")
    
    if 'min_photon_chi_for_xi' in f:
        
        dataset = f['min_photon_chi_for_xi']
        
        size_particle_chi = dataset.attrs["size_particle_chi"]
        
        print(" Table min_photon_chi_for_xi: ")
        
        message = ""
        
        for i in range(0,size_particle_chi,10):
            for j in range(i,min(i+10,size_particle_chi)):
                v = dataset[i]
                message += "{:}".format(v) + ", "
            message += "\n"
                
        print(message)

# ______________________________________________________________________________
# Read the table xip

if if_file_exist(path):

    f = h5py.File(path, "r")

    if 'xi' in f:
        
        xi = np.array(f['xi'])
        
        size_particle_chi = f['xi'].attrs["size_particle_chi"]
        size_photon_chi = f['xi'].attrs["size_photon_chi"]

        print(" Table xi: ")
        
        message = ""
        
        for i in range(size_particle_chi):
            message += "   "
            for j in range(size_photon_chi):
                v = xi[i,j]
                message += "{:}".format(v) + ", "
            message += "\n"
                
        print(message)
