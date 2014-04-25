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

# # ---------------------------------------------
# # Initialize MPI environment 
# # ---------------------------------------------
 
comm=MPI.COMM_WORLD
my_rank=comm.Get_rank()
size=comm.Get_size()


# ---------------------------------------------
# Parsing of the arguments from the command line
# ---------------------------------------------
parser = argparse.ArgumentParser(description=' this is a parallel python code to plot the information from the file scalars.txt. You need to install py27-mpi4py @1.3_3 or other version. To run the code :  openmpirun -np ncpus python scalars.py -a VariableToBePlot -c Colors -s SamePlot   . ' )
parser.add_argument('-a',action='store',default='none',dest='args', nargs='+',help='(arg1 arg2) (arg3 arg4) ..')
parser.add_argument('-d',action='store',default='scalars.txt',dest='dir', nargs='+',help='directory where scalars.txt is located, for ex. ../../scalars.txt default=scalars.txt')
parser.add_argument('-D',action='store',default='plot_scalars',dest='dir_plot', nargs='+',help='directory where plots will be located, for ex. ../../toto default=plot_scalars')
parser.add_argument('-c',action='store',default='rgbcmykw',dest='colors', nargs='+',help='color for each plot')
parser.add_argument('-o',action='store',default='',dest='lim',nargs='+',help='xmin xmax ymin ymax')
parser.add_argument('-s', action='store_false', default=True,dest='same_plot',help='plot on same plot (default:True)')
parser.add_argument('-n',action='store_false',default=True,dest='clean_folder',help='clean the folder (default:False)') 
parser.add_argument('-f', action='store', default='.png',dest='format',help='plot format (default:True)')
what=parser.parse_args()
args=what.args
colors=what.colors
lim=what.lim
same=what.same_plot
clean_folder=what.clean_folder
file_format=what.format
n_plots=len(args)/2

if(what.dir!='scalars.txt'):
    dir=str(what.dir)[2:-2]
else:
    dir=str(what.dir)
if(what.dir_plot!='plot_scalars'):
    path=str(what.dir_plot)[2:-2]
else:
    path=str(what.dir_plot)
    

if (len(args)%2!=0):
    sys.exit('arguments to be plotted not well defined')
elif(args=='none'):
    sys.exit('You must specify at least two arguments in the option -a')
if (my_rank==0):
    if (colors=='rgbcmykw'):
        print 'warning, default colors are used (rgbcmykw)'
    elif(len(color)!=n_plots):
        sys.exit('colors not well defined')
    if(same==True)&(len(lim)!=0):
        if(len(lim)!=4):
            sys.exit('axis limits missing in one direction, there must be 4 numbers (xmin, xmax,ymin,ymax)')
    



#creation of the directory 
if(clean_folder==True):
    if(my_rank==0):
        if (os.path.exists(path)==False):
	        os.makedirs(path)
        else:
		    shutil.rmtree(path)
		    os.makedirs(path)
# # ---------------------------------------------
# # Reading from scalars.txt
# # ---------------------------------------------
with open(dir) as file:
    data=file.readlines()
    n_lines=len(data)
# if(my_rank==0):
#     start_time = time.time()
 
found_x=[False]*n_plots
found_y=[False]*n_plots
found_l=False
first_line=0
n_data_lines=0
pattern=np.zeros(size+1) 
n_lines_pp=np.ones(size)
scalars=[]
col_x=[int]*n_plots
col_y=[int]*n_plots


for line in data:
    words=line.split()
    scalars.append(words)
    if(found_l==False):
        if(words[0]!='#'):
            found_l=True
        else:
            first_line+=1
            
    for p_iter in range(0,n_plots):
        if(found_x[p_iter]==False)|(found_y[p_iter]==False):
            for l in range(0,len(words)):
                if(found_x[p_iter]==False):
                    if(args[p_iter*2]==words[l]): 
                        col_x[p_iter]=int(words[l-1])-1
                        found_x[p_iter]=True
                if(found_y[p_iter]==False):
                    if(args[p_iter*2+1]==words[l]): 
                        col_y[p_iter]=int(words[l-1])-1                      
                        found_y[p_iter]=True
        
for p_iter in range(0,n_plots):       
    if(found_x[p_iter]==False):
        sys.exit('the arguments '+ args[p_iter*2]+ ' does not exist')
    elif(found_y[p_iter]==False):
        sys.exit('the arguments '+ args[p_iter*2+1]+ ' does not exist')
                    
n_data_lines=n_lines-first_line
n_lines_pp*=(n_data_lines/size)
residual=n_data_lines%size
for j in range(0,residual):
    n_lines_pp[j]+=1
pattern[0]=first_line
for r in range(1,size):
    pattern[r]=pattern[r-1]+n_lines_pp[r-1]
pattern[size]=n_lines
loc_x=np.zeros((n_plots,n_data_lines))
loc_y=np.zeros((n_plots,n_data_lines))
gl_x=np.zeros((n_plots,n_data_lines))
gl_y=np.zeros((n_plots,n_data_lines))
    
for p_iter in range(0,n_plots):   
    for l in range(int(pattern[my_rank]),int(pattern[my_rank+1])):
        loc_x[p_iter][l-first_line]=float(scalars[l][col_x[p_iter]])
        loc_y[p_iter][l-first_line]=float(scalars[l][col_y[p_iter]])

my_legend=[""]*n_plots
if(same==True):
    title="/"
else:
    title=["/"]*n_plots
    
for a in range(0,len(args)):
    if(same==True):
        title+='$'+args[a]
        if(a%2==0):
            my_legend[a/2]=args[a+1]
    else:
        title[a/2]+='$'+args[a]
        if(a%2==0):
            my_legend[a/2]=args[a+1]

        
        

gl_x=comm.reduce(loc_x,MPI.SUM) 
gl_y=comm.reduce(loc_y,MPI.SUM)                      
if(my_rank==0):
    if(same==True):
        for p_iter in range(0,n_plots):
            plt.plot(gl_x[p_iter],gl_y[p_iter],colors[p_iter],label=my_legend[p_iter])
            plt.legend(loc=2)
        if(len(lim)!=0):
            plt.xlim([float(lim[0]),float(lim[1])])
            plt.ylim([float(lim[2]),float(lim[3])])
        plt.xlabel(args[0])
        plt.savefig(path+title+file_format)
    else:
        for p_iter in range(0,n_plots):
            clf()
            figure()
            plot(gl_x[p_iter],gl_y[p_iter],colors[p_iter])
            legend([my_legend[p_iter]],loc=2)
            plt.xlabel(args[p_iter*2])
            plt.ylabel(args[p_iter*2+1])
            if(len(lim)!=0):
                xlim([float(lim[0]),float(lim[1])])
                plt.ylim([float(lim[2]),float(lim[3])])
            savefig(path+title[p_iter]+file_format)
            clf()
            