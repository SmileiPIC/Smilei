#Read the Fields.h5 file.
import os
import scipy
import matplotlib
from matplotlib import pyplot
import tables
import sys
import math
import time;
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
p = comm.Get_size()

time.sleep(rank)


######### INPUT THE FOLLOWING PARAMETERS ############################################
#directory = "/sps/beck/smilei/tst2d_accelec" #Directory of the simulation data
directory = "/sps/beck/smilei/benchmark" #Directory of the simulation data
#list_fields=['Ex','Ey','Ez','Bz_m','Rho_electron1','Rho_proton','Jx','Jy','Jz'] #List of the fields you want to extract (Ei,Bi,Ji_sn,Rho,rho_sn, where n is the number of the species and i the direction (x,y,z))
list_fields=['Ex','Ey','Ez','Bz_m','Bx_m','By_m','Rho_electron','Rho_proton','Jx_electron'] 
#list_fields=['Ex','Ey','Ez','Bz_m','Bx_m','By_m'] 
#list_fields=['Rho_electron','Rho_proton'] 
#list_fields=['Ey','Ex'] 
#list_fields=['Ex','Ey','Ez','Rho','rho_s1','rho_s0','Jx_s1'] 
#list_fields=['Bx_m','By_m','Bz_m','Jy_s1','Jz_s1'] 
first_cycle = 1
last_cycle = 4032
cycle_step = 1000 # step between two displayed cycles (minimum is given by outputcyle of the inputfile)
plot_on_axis = 0 # Also plot 1D graph of the quantities on axis if ==1.
suffix = "benchmark" #Suffix to be added in produced files name.
####################################################################################

filename = directory + "/Fields.h5" # Proc file name
h5file = tables.openFile(filename, mode = "r", title = "Fields_file")

existing_files = scipy.arange(0,last_cycle,cycle_step)
asked_files = scipy.arange(first_cycle, last_cycle)
list_files=scipy.intersect1d(existing_files,asked_files)
nb_files=len(list_files)

Filebyproc=nb_files/p #Minimum Nomber of files you have to read
reste=nb_files-Filebyproc*p
if (rank<reste):
    Filerange=scipy.arange(rank*(Filebyproc+1),(rank+1)*(Filebyproc+1))
else:
    Filerange=scipy.arange((reste*(Filebyproc+1)+(rank-reste)*Filebyproc),(reste*(Filebyproc+1)+(rank-reste+1)*Filebyproc))
#print nb_files, 'filerange =', Filerange, rank
list_cycle = list(list_files[Filerange])
print "cycle", list_cycle

#list_cycle=range(first_cycle,last_cycle+1,cycle_step) 

print h5file.root._f_listNodes
if (rank==0):
    print h5file.root._f_getChild(repr(list_cycle[0]).rjust(10,"0"))._f_listNodes

for field in list_fields:
    for cycle in list_cycle:    
            num = repr(cycle).rjust(10,"0")
            group = h5file.root._f_getChild(num)
            array = group._f_getChild(field).read()
            print field, cycle, scipy.shape(array)
            if cycle == list_cycle[0]:
                v_minl = array.min() #local
                v_maxl = array.max() #local
                v_minl = -0.008
                #v_min = comm.allreduce(v_minl,op=MPI.MIN)
                #v_max = comm.allreduce(v_maxl,op=MPI.MAX)
                v_min = comm.bcast(v_minl,root=0)
                v_max = comm.bcast(v_maxl,root=0)
            #    print v_minl,"  ",v_min,"  ",v_maxl,"  ",v_max
            pyplot.figure()
            #pyplot.imshow(array.T , origin=0,vmin=v_min,vmax=v_max,aspect="auto")
            pyplot.imshow(array.T , origin=0, aspect="auto")
            pyplot.title(field + " " + num)
            pyplot.colorbar()
            pyplot.savefig(directory + "/"+field+num+".png",format="png")
            pyplot.close()
    comm.Barrier()
    if (rank==0):
        os.system("convert -delay 10  -loop 0 "+directory+"/"+field+"*.png "+ directory+"/"+field+suffix+".gif")

if (rank==0):
    os.system("rm -f "+directory+"/*.png")

comm.Barrier()
if plot_on_axis == 1:
    for field in list_fields:
        for cycle in list_cycle:    
                print field, cycle
                num = repr(cycle).rjust(10,"0")
                group = h5file.root._f_getChild(num)
                array = group._f_getChild(field).read()
                pyplot.figure()
                pyplot.plot(array[:,scipy.shape(array)[1]/2],'r+')
                pyplot.title(field + " " + num)
                pyplot.savefig(directory + "/"+field+num+".png",format="png")
                pyplot.close()
        comm.Barrier()
        if (rank==0):
            os.system("convert -delay 10  -loop 0 "+directory+"/"+field+"*.png "+ directory+"/"+field+suffix+"onaxis.gif")
    
    if (rank==0):
        os.system("rm -f "+directory+"/*.png")

h5file.close()



