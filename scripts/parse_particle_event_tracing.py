import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import os,sys
import re
import math

# This script parses the Particle Event Tracing diagnostic, that can be 
# used with or without task parallelization. 
#
# This diagnostic stores the time when a particle event starts and ends, as well
# as the information on which MPI rank and OpenMP thread performed this event.
#
# For the purpose of this diagnostic, particle events are e.g. the PIC operators
# acting on a group of particles corresponding to given ibin-ispec-ipatch 
# combination.
#
# This diagnostic thus may help in visualizing the scheduling of macro-particle
# operations. Warning: with many bins, species or patches, the resulting plot 
# may become unreadable.
#
# Use of the script: select an iteration and if particle event data is available
# for that iteration a Gant chart plot will be created. Each rectangle represents
# a particle event for a given ibin-ispec-ipatch combination. The color-code 
# of the rectangles tells the operator treating those particles.
#
# Warning: the starting point of the axes is set reading the first event of
# OpenMP thread 0 of each rank. If this thread is slower than the others, negative
# times may appear in the plot of other threads' scheduling.

# ------- Inputs
iteration_to_analyze = 1500
Plot = True



# ------- Auxiliary functions
###########

# Hypothesis underlying these functions:
# the tracing file name has a structure particle_event_tracing_rank_X_thread_Y.txt

def get_rank_from_filename(filename):
	filename_components = filename.split("_")
	if len(filename_components)>3:
		return int(filename_components[3])
	else:	
		return -1

def get_thread_from_filename(filename):
	filename_components = filename.replace('.', '_').split("_")
	return int(filename_components[5])

def search_for(string_found,string_to_search,line):
        # this function tells if you should keep parsing for a certain parameter
    if string_found == True: # parameter already found, no need to parse it again
        return False
    else: # look for the parameter if present in the string
        return str(line).find(string_to_search)>0

###########
    
# Reads an event from two given lines about the same event

# Important hypothesis if OpenMP tasks are used:
# OpenMP tasks are tied (i.e. thread X starts an event, thread X ends this event)

# The event lines structures is RelativeTime Start/End EventName(1word)
    
dict_event_name = {0:'Interp', 1:'Push', 2:'BC', 3: 'Proj',4: 'DensityReduction',
                   5: 'Ionization',6: 'Radiation',7: 'MultiphotonBW',
                   8: 'IonizReduction',9: 'RadReduction',10: 'MBWReduction',11: 'VectoKeys'}    
    
def read_thread_event(line_start_event,line_end_event):
	words_line_start_event = line_start_event.split(" ")
	words_line_end_event   = line_end_event.split(" ")
	time_start_event = float(words_line_start_event[0])
	time_end_event = float(words_line_end_event[0])
	event_name = dict_event_name[int(words_line_start_event[2])]  # name of event
	# 0: Interp
	# 1: Push
	# 2: BC
	# 3: Proj
	# 4: Density Reduction
	# 5: Ionization
	# 6: Radiation
	# 7: Multiphoton Breit Wheeler
	# 8: Ionization Reduction
	# 9: Radiation Reduction
	#10: Multiphoton Breit Wheeler Reduction
	#11: Computation of keys and count for vectorization
	return time_start_event,time_end_event,event_name

# ------- Read data
start_path = os.getcwd()

# Find number of MPI
N_MPI = 1
MPI_found = False
with open("smilei.log") as f:
	lines = f.readlines()
	for iline,line in enumerate(lines):
		if search_for(MPI_found,"Number of MPI process",line):
			N_MPI = int(re.search(r'\d+', str(line)).group())
			MPI_found = True
			break
print("Number of MPI = ",N_MPI)

if Plot==True:
	plt.ion()
# For each rank, read and plot the tracing data
for rank_to_analyze in range(0,N_MPI):
	print("Reading tracing for rank ",rank_to_analyze)
	tracing_files = []
	# read only files for the desired MPI rank
	for file in os.listdir(start_path):
		if (file.startswith("particle_event_tracing_") and (get_rank_from_filename(file)==rank_to_analyze)):
			tracing_files.append(file)

	tracing_files = sorted(tracing_files)

	# for each file to read (one per thread, read the events in the desired iteration)
	threads = []

	Events_to_plot = []
	ref_time = 0.
	for tracing_file in tracing_files:
		event_lines_to_read = []
		thread_number = get_thread_from_filename(tracing_file)
		threads.append(thread_number)
		Iteration_Found = False
		with open(os.path.join(start_path,tracing_file)) as f: # read file corresponding to rank and thread
			lines = f.readlines()
			keepCurrentLine = False
			# Isolate the lines for the desired iteration
			for line in lines:
				if line == "Start Iteration "+str(iteration_to_analyze)+"\n":
					keepCurrentLine = True
					Iteration_Found = True
				if keepCurrentLine:
					event_lines_to_read.append(line)
				if line == "End Iteration "+str(iteration_to_analyze)+"\n":
					keepCurrentLine = False
			event_lines_to_read = [line[:-1] for line in event_lines_to_read]
			event_lines_to_read = event_lines_to_read[1:-1]
			if Iteration_Found == False:
				print("Error, iteration not found")
			if thread_number == 0:
				ref_time = float(event_lines_to_read[0].split(" ")[0])
		# Start extracting data from the desired iteration	
		Nevents = int(len(event_lines_to_read)/2)
		for i_event in range(0,Nevents):
			line_start_event = event_lines_to_read[2*i_event]
			line_end_event = event_lines_to_read[2*i_event+1]
			time_start_event,time_end_event,event_name = read_thread_event(line_start_event,line_end_event) 
			Events_to_plot.append([thread_number, time_start_event-ref_time,time_end_event-ref_time,event_name])	

	# for event in Events_to_plot:		 
	# 	print(event)	

	# ------- Plot the timeline bar graph

	if (Plot == False):
		quit()
	print("Plotting tracing for rank ",rank_to_analyze)


	colormapping = {'Interp' : "r", 'Push' : "yellow", 'BC' : "g", 'Proj' : "cyan", 'DensityReduction' : "b",'EnergyLost' : "purple",
                'Ionization' : "orange", 'Radiation' : "darkgreen", 'MultiphotonBW' : "blueviolet", 
                'IonizReduction' : "gold", 'RadReduction' : "forestgreen",'MBWReduction' : "fuchsia",'VectoKeys': "darkcyan"}
	vert_shift = {'Interp' : 0., 'Push' : 0.1, 'BC' : 0.2, 'Proj' : 0.3, 'DensityReduction' : 0.5, 
              'Ionization' : 0.6, 'Radiation' : 0.7, 'MultiphotonBW' : 0.8,
              'IonizReduction' : 0.6, 'RadReduction' : 0.7, 'MBWReduction' : 0.8, 'VectoKeys':0.5}

	### Plot all events

	fig, gnt = plt.subplots()
	fig.set_size_inches(7, 5.5)

	# Setting labels for x-axis and y-axis
	gnt.set_xlabel('Time (s)')
	gnt.set_ylabel('Thread')

	# Setting ticks on y-axis
	gnt.set_yticks(threads)

	gnt.set_yticklabels([str(thread) for thread in threads])

	lbl_pool = []
	for event in Events_to_plot:
		lbl = event[3]
		if lbl in lbl_pool:
			prefix = '_'
		else:
			lbl_pool.append(lbl)
			prefix=''
		gnt.broken_barh([(event[1], (event[2]-event[1]))], (event[0]-0.05-vert_shift[event[3]], 0.5), facecolors =(colormapping[event[3]]),label=prefix+lbl, edgecolors=('k'), linewidth=0.2)

	plt.legend( bbox_to_anchor=(0.2,-0.3,0.65,0.2), loc="lower center",
                mode="expand", borderaxespad=0, ncol=4, frameon=False)

	xlim = gnt.get_xlim()
	plt.title("Rank "+str(rank_to_analyze))
	plt.tight_layout()
	print("Rank ",rank_to_analyze," :")
	print("Xlim = ",xlim)



