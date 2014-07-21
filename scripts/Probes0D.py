#! /usr/bin/env python

# -----------
# All imports
# -----------
import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import sys
from mpi4py import MPI
from pylab import *
import argparse
import tables as tbl
# # ---------------------------------------------
# # Initialize MPI environment 
# # ---------------------------------------------
 
comm=MPI.COMM_WORLD
my_rank=comm.Get_rank()
size=comm.Get_size()

# ---------------------------------------------
# Parsing of the arguments from the command line
# ---------------------------------------------
parser = argparse.ArgumentParser(description=' ' )

parser.add_argument('-a',action='store',default='all',dest='args',nargs='+',help='arg1 arg2..they will be plotted in function of time. default=it plots all the variables')
parser.add_argument('-p',action='store',default='all',dest='probes',nargs='+',help='it plots the data of a probe. default= it plots all the probes')
parser.add_argument('-d',action='store',default='Probes0D.h5',dest='dir', nargs='+',help='directory where Probes0D.h5 is located, for ex. ../../Probes0D.h5 default=Probes0D.h5')
parser.add_argument('-D',action='store',default='plot_probes0D',dest='dir_plot', nargs='+',help='directory where plots will be located, for ex. ../../toto default=plot_probes0D')
parser.add_argument('-c',action='store',default='rgbcmykw',dest='colors', nargs='+',help='color for each plot')
parser.add_argument('-o',action='store',default='',dest='lim',nargs='+',help='xmin xmax ymin ymax')
parser.add_argument('-s', action='store_false', default=True,dest='same_plot',help='plot on same plot (default:True)')
parser.add_argument('-n',action='store_false',default=True,dest='clean_folder',help='clean the folder (default:False)') 
parser.add_argument('-f', action='store', default='.png',dest='format',help='plot format (default:True)')
parser.add_argument('-fft',action='store',default='None',dest='fft_args',nargs='+',help='it is implementing the FFT on the args indicated')
what=parser.parse_args()

if(what.dir!='Probe0D.h5'):
    dir=str(what.dir)[2:-2]
else:
    dir=str(what.dir)
if(what.dir_plot!='plot_probes0D'):
    path=str(what.dir_plot)[2:-2]
    print path
else:
    path=str(what.dir_plot)

data_list=['time','Ex','Ey','Ez','Bx','By','Bz']
mapping=[0,1,2,3,4,5,6]
if(what.args=='all'):
    args=mapping[1:]
else:
    args=[]
    for w in range(0,len(what.args)):
        found=False
        for l in range(0,len(mapping)):
            if(what.args[w]==data_list[l]):
                args.append(mapping[l])
                found=True
        if(found==False):
            sys.exit(what.args[w]+' does not exist')

if(what.clean_folder==True):
     if(my_rank==0):
        if (os.path.exists(path)==False):
	        os.makedirs(path)
        else:
		    shutil.rmtree(path)
		    os.makedirs(path)




f=tbl.openFile(dir)
data_h=[]
for n in f.root:
    data_h.append(n.read())
data=array(data_h)        
id_probes=[]
if(what.probes=="all"):
    n_probes=len(data)
    for p in range(0,n_probes):
        id_probes.append(p)
else:
    n_probes=len(what.probes)
    for p in range(0,n_probes):
        id_probes.append(int(what.probes[p]))

if(what.colors!='rgbcmykw'):
    if(len(what.colors)!=len(args)):
        sys.exit('colors not well defined, there are '+ str(len(what.colors)) + ' colors specified for ' +str(len(args))+ ' arguments')
lim=what.lim
if(what.same_plot==True):
    n_pp=n_probes/size
    res=n_probes%size
    for k in range(0,res):
        if(my_rank==k):
            n_pp+=1
    mapping=[]
    n_pp=comm.gather(n_pp)
    n_pp=comm.bcast(n_pp)
    for r in range(0,size):
        for s in range(0,n_pp[r]):
            mapping.append(r)
    if(size!=1):
        for p in range(0,len(mapping)):
            if(my_rank==mapping[p]):
                title=path+'/p'+str(id_probes[p])+'$time$'
                for a in range(0,len(args)):
                    plot(data[id_probes[p]][:,args[a]],label=data_list[args[a]],color=what.colors[a])
                    title+='$'+data_list[args[a]]
                    legend(loc=1)
                xlabel('Time')
                if(len(lim)!=0):
                    xlim([float(lim[0]),float(lim[1])])
                    ylim([float(lim[2]),float(lim[3])])
                savefig(title+what.format)
                clf()
    else:  
        for p in range(0,n_probes):           
            title=path+'/p'+str(id_probes[p])+'$time$'
            for a in range(0,len(args)):
                plot(data[id_probes[p]][:,args[a]],label=data_list[args[a]],color=what.colors[a])
                title+='$'+data_list[args[a]]
                legend(loc=1)
            xlabel('Time')
            if(len(lim)!=0):
                xlim([float(lim[0]),float(lim[1])])
                ylim([float(lim[2]),float(lim[3])])
            savefig(title+what.format)
            clf()
else:
    if(size!=1):
        mapping_var=[]
        mapping_pr=[]
        mapping=[]
        n_pp=len(args)*n_probes/size
        res=(len(args)*n_probes)%size
        for k in range (0,res):
            if(my_rank==k):
                n_pp+=1
        for p in range(0,n_probes):
            for a in range(0,len(args)):
                mapping_pr.append(id_probes[p])
                mapping_var.append(args[a])
        n_pp=comm.gather(n_pp)
        n_pp=comm.bcast(n_pp)
        for r in range(0,size):
            for s in range(0,n_pp[r]):
                mapping.append(r)
        for m in range(0,len(mapping)):
            if(my_rank==mapping[m]):
                title=path+'/p'+str(mapping_pr[m])+'$time$'+data_list[mapping_var[m]]
                plot(data[mapping_pr[m]][:,mapping_var[m]],label=data_list[mapping_var[m]])
                legend(loc=1)
                xlabel('Time')
                if(len(lim)!=0):
                    xlim([float(lim[0]),float(lim[1])])
                    ylim([float(lim[2]),float(lim[3])])
                savefig(title+what.format)
                clf()
    else:
        for p in range(0, n_probes):
            for a in range(0,len(args)):
                title=path+'/p'+str(id_probes[p])+'$time$'+data_list[args[a]]
                plot(data[id_probes[p]][:,args[a]],label=data_list[args[a]])
                legend(loc=1)
                xlabel('Time')
                if(len(lim)!=0):
                    xlim([float(lim[0]),float(lim[1])])
                    ylim([float(lim[2]),float(lim[3])])
                savefig(title+what.format)
                clf()

               
if(what.fft_args!='None'):
    message='FFT evaluated on '
    for a in range(0,len(what.fft_args)):
        message+=what.fft_args[a]+' '
    if(my_rank==0):
        print message
        print 'WARNING: FFT is not implemented yet'

f.close()