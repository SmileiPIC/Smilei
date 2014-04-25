#! /usr/bin/env python

# -----------
# All imports
# -----------
import tables as tbl
import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import sys
import argparse


from mpi4py import MPI
from pylab import *
# ---------------------------------------------
# Initialize MPI environment and folder
# ---------------------------------------------
comm=MPI.COMM_WORLD
my_rank=comm.Get_rank()
size=comm.Get_size()

# ---------------------------------------------
# Parsing of the arguments from the command line
# ---------------------------------------------
parser = argparse.ArgumentParser(description=' this is a parallel python code to plot the information from the file scalars.txt. You need to install py27-mpi4py @1.3_3 or other version. To run the code :  openmpirun -np ncpus python scalars.py -a VariableToBePlot -c Colors -s SamePlot   . ' )
parser.add_argument('-a',action='store',default='none',dest='args', nargs='+',help=' name and directory where the input SMILEI file is located')
parser.add_argument('-d',action='store',default='Fields.h5',dest='dir', nargs='+',help='directory where Fields.h5 is located, for ex. ../../Fields.h5 default=Fields.h5')
parser.add_argument('-D',action='store',default='plot_fields',dest='dir_plot', nargs='+',help='directory where plots will be located, for ex. ../../toto default=plot_fields') 
parser.add_argument('-f', action='store', default='.png',dest='format',help='plot format (default:True)')
what=parser.parse_args()

#-----------------------------------------------
# check of variables
#-----------------------------------------------
if(what.args=='none'):
    sys.exit('you must specify name and directory where the imput SMILEI file is located')
if(what.dir=='Fields.h5'):
    dir=str(what.dir)
else:
    dir=str(what.dir)[2:-2]
if(what.dir_plot=='plot_fields'):
    dir_plot=str(what.dir_plot)
else:
    dir_plot=str(what.dir_plot)[2:-2]
#-----------------------------------------------
# 