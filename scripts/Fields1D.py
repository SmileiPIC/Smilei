#! /usr/bin/env python
print '# Warning, this is a parallel python code to plot some information about the fields (Eright & Eleft). You need to install py27-mpi4py @1.3_3 or other version. To run the code: openmpirun -np ncpus ipython-2.7 ./par_diag_profiles.py #'
# -----------
# All imports
# -----------
import tables as tbl
import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import sys


from mpi4py import MPI
from pylab import *
# ---------------------------------------------
# Initialize MPI environment and folder
# ---------------------------------------------
comm=MPI.COMM_WORLD
my_rank=comm.Get_rank()
size=comm.Get_size()
# new directory where to save the plots
path='plot'
# format of the plots
file_format='.png'
# creation of the directory 
if(my_rank==0):
	if os.path.exists(path)==False:
		os.makedirs(path)
	else:
		shutil.rmtree(path)
		os.makedirs(path)
		
# -------------------------------------------------
# Open Fields.h5 in f 
# -------------------------------------------------
f = tbl.openFile("Fields.h5")	 
# -------------------------------------------------
# Density profile function
# -------------------------------------------------
# If you use a new density profile in SMILEI, you need to add it here
def define_k(plasma_geometry):
	if(plasma_geometry=='constant'):
		for k in range(0,n_points):
			if k<(n_v):
				k0[k]=k_void
			elif k<(n_v+n_pl):
				k0[k]=k_plasma
			elif k<n_points:
				k0[k]=k_void
	elif(plasma_geometry=='trap'):
	 	if 'slope_length' in globals():
 			for k in range(0,n_points):
 				if k<(n_v):
					k0[k]=k_void
				elif k<(n_v+n_sl):
					k0[k]=k_void+((k_plasma-k_void)/n_sl)*(k-n_v)
				elif k<(n_v-n_sl+n_pl):
					k0[k]=k_plasma
				elif k<(n_v+n_pl):
					k0[k]=k_plasma+((k_void-k_plasma)/n_sl)*(k-(n_v+n_pl-n_sl))
				elif k<n_points :
					k0[k]=k_void
		else:
			for k in range(0,n_points):
				if k<(n_v):
					k0[k]=k_void
				elif k<(n_v+n_sl_l):
					k0[k]=k_void+((k_plasma-k_void)/n_sl_l)*(k-n_v)
				elif k<(n_v+n_pl-n_sl_r):
					k0[k]=k_plasma
				elif k<(n_v+n_pl):
					k0[k]=k_plasma+((k_void-k_plasma)/n_sl_r)*(k-(n_v+n_pl-n_sl_r))
				elif k<n_points :
					k0[k]=k_void
	elif(plasma_geometry=='triangular'):
		right_slope_length=plasma_length-slope_length
		n_sl_r=right_slope_length*res_space
		for k in range(0,n_points):
			if k<(n_v):
				k0[k]=k_void
			elif k<(n_v+n_sl):
				k0[k]=k_void+((k_plasma-k_void)/n_sl)*(k-n_v)
			elif k<(n_v+n_pl):
				k0[k]=k_plasma+((k_void-k_plasma)/n_sl_r)*(k-(n_v+n_pl-n_sl_r))
			elif k<n_points :
				k0[k]=k_void
	else:
		sys.exit(plasma_geometry +' density profile not yet defined!!')
# -------------------------------------------------
# Simulation data
# -------------------------------------------------
plasma_geometry = 'constant'	
res_time = 10
sim_time = 300

res_space  = 8
sim_length = 130

vacuum_length	= 5
plasma_length	= 120

density=0.3
I_norm= 8.5e-2
fieldDump_every = 125


# -------------------------------------------------
# plotting
# -------------------------------------------------
n_time= res_time*sim_time
n_points = res_space*sim_length
n_v=vacuum_length*res_space
n_pl=plasma_length*res_space

if 'slope_length' in globals():
	n_sl=slope_length*res_space
if 'left_slope_length' in globals():
	n_sl_l=left_slope_length*res_space
if 'left_slope_length' in globals():
	n_sl_r=right_slope_length*res_space

k_void=1.0
k_plasma=sqrt(1.0/(1.0-density))


time = []	 # time axis
Ez   = []	 # Ez
By   = []	 # By
ni   = []	 # ni
ne   = []	 # ne
k0 = ndarray(n_points)
Ex  =[]

groups_name = []

# it defines the values of k0 on the domain
define_k(plasma_geometry)

for n in f.root:
	groups_name.append(n._v_name)

n_plot_tot= len(groups_name)
n_plot_pp=n_plot_tot/size
residual=n_plot_tot%size

for j in range(1,residual+1):
	if(my_rank==j-1):
		n_plot_pp+=1
	
my_iter=my_rank*n_plot_pp
my_groups=[]
for n in range(my_iter,my_iter+n_plot_pp):
	my_groups.append(groups_name[n])
	
for n in range(0,len(my_groups)):
	for g in f.walkGroups("/"+my_groups[n]):
		time.append(g._v_name)
		Ez.append(g.Ez)
		By.append(g.By_m)
		ni.append(g.rho_s0)
		ne.append(g.rho_s1)
		Ex.append(g.Ex)

# --------------------------------------------------------------------
# Calculate forward/backward moving waves & electron and ion densities
# --------------------------------------------------------------------
for t in range(0,len(time)):
	  Eleft		= []
	  Eright	= []
	  ni_ov_ni0 = []
	  ne_ov_ne0 = []
	  Ex_ = []
	  
	  for k in range(0,n_points):
		Eleft.append(	  ( 0.5*(Ez[t][k] - 0.5 *k0[k]* (By[t][k+1]+By[t][k])) )/I_norm )
		Eright.append(	  ( Ez[t][k] - I_norm*Eleft[k] )/I_norm )
		Ex_.append(Ex[t][k])
		ni_ov_ni0.append(  ni[t][k]/density )
		ne_ov_ne0.append( -ne[t][k]/density )

	  if(my_rank==0):
		print str(float(t)/float(len(time)-1))+'% done'
	  plot(ne_ov_ne0,color='y')
	  plot(ni_ov_ni0,color='g')
	  plot(Eleft,color='r')
	  plot(Eright,color='k')
	  name= str(my_groups[t])
	  savefig(path+'/diag_profiles_'+name+file_format)
	
f.close()