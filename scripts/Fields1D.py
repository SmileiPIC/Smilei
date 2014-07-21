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
parser.add_argument('-p',action='store',default='0',dest='l_norm', nargs='+',help='lesar which respect to make the normalization') 
parser.add_argument('-f', action='store', default='.png',dest='format',help='plot format (default:True)')
what=parser.parse_args()

#-----------------------------------------------
# check of variables
#-----------------------------------------------
if(what.args=='none'):
    sys.exit('you must specify name and directory where the imput SMILEI file is located')
else:
    input=str(what.args)[2:-2]
if(what.dir=='Fields.h5'):
    dir=str(what.dir)
else:
    dir=str(what.dir)[2:-2]
if(what.dir_plot=='plot_fields'):
    dir_plot=str(what.dir_plot)
else:
    dir_plot=str(what.dir_plot)[2:-2]
    
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
#-----------------------------------------------
# creation of the directory
#-----------------------------------------------
if(my_rank==0):
	if os.path.exists(dir_plot)==False:
		os.makedirs(dir_plot)
	else:
		shutil.rmtree(dir_plot)
		os.makedirs(dir_plot)
#--------------------------------------------
#read the input file
#---------------------------------------
with open(input) as file:
    input_data=file.readlines()
    

species=[]
I_lasers=[]
for line in input_data:
    words=line.split()
    if(len(words)!=0):
        if(words[0]!='#'):
            if(words[0]=='res_time'):
                res_time=float(words[2])
            elif(words[0]=='sim_time'):
                sim_time=float(words[2])
            elif(words[0]=='res_space'):
                res_space=float(words[2])
            elif(words[0]=='sim_length'):
                sim_length=float(words[2])
            elif(words[0]=='plasma_geometry'):
                plasma_geometry=str(words[2])    
            elif(words[0]=='vacuum_length'):
                vacuum_length=float(words[2])
            elif(words[0]=='plasma_length'):
                plasma_length=float(words[2])
            elif(words[0]=='slope_length'):
                slope_length=float(words[2])
            elif(words[0]=='left_slope_length'):
                left_slope_length=float(words[2])
            elif(words[0]=='right_slope_length'):
                right_slope_length=float(words[2])
            elif(words[0]=='species_type'):
                species.append(str(words[2]))
            elif(words[0]=='density'):
                n_norm=float(words[2])
            elif(words[0]=='a0'):
                I_lasers.append(float(words[2]))
            elif(words[0]=='fieldDump_every'):
                field_every=float(words[2])

n_species=len(species)
n_lasers=len(I_lasers)   

I_norm=I_lasers[int(str(what.l_norm))]

n_time= res_time*sim_time
n_points = int(res_space*sim_length)
n_v=int(vacuum_length*res_space)
n_pl=int(plasma_length*res_space)   

k_void=1.0
k_plasma=sqrt(1.0/(1.0-n_norm))
k0 = ndarray(n_points)
if 'slope_length' in globals():
	n_sl=slope_length*res_space
if 'left_slope_length' in globals():
	n_sl_l=left_slope_length*res_space
if 'right_slope_length' in globals():
	n_sl_r=right_slope_length*res_space

define_k(plasma_geometry)

f = tbl.openFile(dir)
n_plot_tot=int(n_time/field_every)+1

n_plot_pp=n_plot_tot/size
residual=n_plot_tot%size
for p in range(0,residual):
    if(my_rank==p):
        n_plot_pp+=1
n_plot_pp=comm.gather(n_plot_pp)
n_plot_pp=comm.bcast(n_plot_pp)
pattern=[0]*(size+1)
for r in range(1,size+1):
    pattern[r]=pattern[r-1]+n_plot_pp[r-1]
group_name=[str]*(pattern[my_rank+1]-pattern[my_rank])
for p in range(pattern[my_rank],pattern[my_rank+1]):
    group_name[p-pattern[my_rank]]=str(int(p*field_every)).zfill(10)

time=[]
Ez=[]
By=[]
ni=[]
ne=[]
Ex=[]

for n in range(0,len(group_name)):
	for g in f.walkGroups("/"+group_name[n]):
		time.append(g._v_name)
		Ez.append(g.Ez)
		By.append(g.By_m)
		ni.append(g.rho_s0)
		ne.append(g.rho_s1)
		# for nn in range(0,n_species):
# 		    if(species[nn]=='ion'):
# 		        ni.append(g.rho_s0)
# 		    elif(species[nn]=='eon'):
# 		        ne.append(g.rho_s1)
		    
for t in range(0,len(group_name)):
	Eleft		= []
	Eright	= []
	ni_ov_ni0 = []
	ne_ov_ne0 = []
	for k in range(0,n_points):
	    ni_ov_ni0.append(  ni[t][k]/n_norm )
	    ne_ov_ne0.append( -ne[t][k]/n_norm )
	    Eleft.append(( 0.5*(Ez[t][k] - 0.5 *k0[k]* (By[t][k+1]+By[t][k])) )/I_norm )
	    Eright.append(( Ez[t][k] - I_norm*Eleft[k] )/I_norm )
        
	
	plt.plot(Eleft,color='r',label='Eleft')
	plt.plot(Eright,color='b',label='Eright')
	plt.plot(ni_ov_ni0,color='y',label='ion density')
	plt.plot(ne_ov_ne0,color='g',label='electron density')
	plt.ylim([-1.5,1.5])
	ax=plt.subplot(111)
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.72, box.height])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plt.savefig(str(group_name[t])+what.format)
	
	clf()
	    
f.close()	