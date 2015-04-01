# -----------------------------------------------------------------------
# HOW TO VIEW FIELDS                              -    F. Perez - 03/2015
# -----------------------------------------------------------------------
#
# >>>>>> Requirements
#   python2.7 with the following packages: numpy, matplotlib, pylab, h5py
#
# >>>>>> First step: invoke python and load this file
#      $ python -i Fields.py
#
# >>>>>> Second step: in the python shell, use the function ""
#
# >>>>>> Examples:



import h5py
import numpy as np
import os.path, glob, re
import matplotlib.pyplot as plt
import pylab
pylab.ion()


# Finds a parameter "param" in the input file
# Argument "after" is a string that must be found before "param"
def findParam(results_path, param, after=None):
	out = ""
	ok = True if after is None else False
	try:
		file = glob.glob(results_path+"/*.in")[0]
	except:
		print "Cannot open input file in "+results_path
		return False
	for line in open(file, 'r'):
		if "#" in line: line = line[:line.find("#")]
		if ok or (after in line and "=" in line):
			ok = True
		else:
			continue
		if param in line and "=" in line:
			out = line.split("=")[1]
			break
	return out.strip()
	

# get all available timesteps
def getAvailableTimesteps(results_path):
	try:
		file = results_path+'/Fields.h5'
		f = h5py.File(file, 'r')
	except:
		print "Cannot open file "+file
		return np.array([])
	times = np.double(f.keys())
	f.close()
	return times

# get all available fields
def getAvailableFields(results_path):
	try:
		file = results_path+'/Fields.h5'
		f = h5py.File(file, 'r')
	except:
		print "Cannot open file "+file
		return np.array([])
	try:
		fields = f.values()[0].keys() # list of fields
	except:
		fields = []
	f.close()
	return fields


# -------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------
def Fields(results_path, field=None, timesteps=None, slice=None,
           units="code", data_log=False, data_min=None, data_max=None,
           xmin=None, xmax=None, ymin=None, ymax=None,
           figure=None):
    
	# Prepare units
	try:
		sim_units  = findParam(results_path, "sim_units")
		if sim_units==False: return
	except:
		sim_units = "norm"
	try:
		cell_length = findParam(results_path, "cell_length")
		cell_length = np.double(cell_length.split())
	except:
		try:
			res_space = findParam(results_path, "res_space")
			res_space = np.double(res_space.split())
			sim_length = findParam(results_path, "sim_length")
			sim_length = np.double(sim_length.split())
			cell_length = sim_length/res_space
		except:
			print "Could not extract 'cell_length' from the input file"
			return
	if units == "nice":
		try:
			wavelength_SI = float( findParam(results_path, "wavelength_SI") )
		except:
			print "Could not extract 'wavelength_SI' from the input file"
			return
		if sim_units == "wavelength": cell_length *= wavelength_SI 
		else :                        cell_length *= wavelength_SI/(2.*np.pi)
		cell_length *= 1e2 # cm
		cell_volume = np.prod(cell_length)
		coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 * cell_volume
		coeff_current = coeff_density * 4.803e-9
	elif units == "code":
		if sim_units == "wavelength": cell_length *= 2.*np.pi
		coeff_density = 1.
		coeff_current = 1.
	if np.size(cell_length)==1:
		cell_length = np.hstack((cell_length,cell_length,cell_length))
	
	# Get available times
	times = getAvailableTimesteps(results_path)
	if np.size(times) == 0:
		print "No fields found in Fields.h5"
		return
	
	# Get available fields
	fields = getAvailableFields(results_path)
	
	# If no field requested, list available fields
	if field is None:
		print "Printing available fields:"
		print "-----------------------------------"
		if len(fields)%3==0:
			print '\n'.join(['\t\t'.join(list(i)) for i in np.reshape(fields,(-1,3))])
		else:
			print "%s\n" % '\n'.join()
		return
	
	# 1 - verifications, initialization
	# -------------------------------------------------------------------
	# Check value of field
	if field not in fields:
		fs = filter(lambda x:field in x, fields)
		if len(fs)==0:		
			print "No field `"+field+"` found in Fields.h5"
			return
		if len(fs)>1:
			print "Several fields match: "+' '.join(fs)
			print "Please be more specific and retry."
			return
		field = fs[0]
	fieldn = fields.index(field) # index of the requested field

	# Check slice is a dict
	if slice is not None  and  type(slice) is not dict:
		print "Argument `slice` must be a dictionary"
		return
	if slice is None: slice = {}
	
	if figure is not None:
		print "Using field "+field
	
	# Open hdf file
	file = results_path+'/Fields.h5'
	f = h5py.File(file, 'r')
	
	
	# 2 - Manage timesteps
	# -------------------------------------------------------------------
	# If timesteps is None, then keep all timesteps
	# otherwise, select timesteps
	if timesteps is not None:
		try:
			ts = np.array(np.double(timesteps),ndmin=1)
			if ts.size==2:
				times = times[ (times>=ts[0]) * (times<=ts[1]) ] # get all times in between bounds
			elif ts.size==1:
				times = np.array([times[(np.abs(times-ts)).argmin()]]) # get nearest time
			else:
				raise Exception()
		except:
			print "Argument `timesteps` must be one or two non-negative integers"
			f.close()
			return
	
	# Need at least one timestep
	if times.size < 1:
		print "Timesteps not found"
		f.close()
		return
	
		
	# 3 - Manage axes
	# -------------------------------------------------------------------
	# Fabricate all axes values
	sample = np.double(f.values()[0].values()[fieldn])
	shape = sample.shape
	naxes = len(shape)
	plot_label = []; plot_centers = []
	slices = [np.array([])]*3
	for iaxis in range(naxes):
		centers = np.linspace(0., shape[iaxis]*cell_length[iaxis], shape[iaxis])
		label = {0:"x", 1:"y", 2:"z"}[iaxis]
		axisunits = "[code units]"
		if units == "nice": axisunits = "[cm]"
		
		if label in slice:
			s = np.double(slice[label])
			if len(s)==1:
				indices = np.array([(np.abs(centers-s)).argmin()])
			elif len(s)==2:
				indices = np.nonzero( (centers>=s[0]) * (centers<=s[1]) )[0]
				if indices.size == 0:
					indices = np.array([(np.abs(centers-s[0])).argmin()])
			else:
				print "Slice along "+label+" must be one number or a list of two numbers"
				return
			if figure is not None:
				if indices.size == 0:
					print "Slice along "+label+" is out of the box range"
					return
				if indices.size == 1:
					print "Sliced at "+label+" = "+str(centers[indices])+" "+axisunits
				else:
					print "Sliced for "+label+" from "+str(centers[indices[ 0]])+" to "+str(centers[indices[-1]])+" "+axisunits
			slices[iaxis] = np.delete(np.arange(shape[iaxis]), indices)
		else:
			plot_centers.append(centers)
			plot_label  .append(label+" "+axisunits)
		
	if len(plot_centers) > 2:
		print "Cannot plot in "+str(len(plot_shape))+"d. You need to 'slice' some axes."
		return
	
	# Build units
	if units == "nice":
		fieldunits = {"B":"T"  ,"E":"V/m"  ,"J":"A"          ,"R":"1/cm$^3$"   }[field[0]]
		unitscoeff = {"B":10710,"E":3.21e12,"J":coeff_current,"R":coeff_density}[field[0]]
		title      = field + "("+fieldunits+")"
	else:
		fieldunits = {"B":"$m_e\omega/e$","E":"$m_ec\omega/e$","J":"$ecn_c$","R":"$n_c$"}[field[0]]
		unitscoeff = {"B":1              ,"E":1               ,"J":1        ,"R":1      }[field[0]]
		title      = field + " in units of "+fieldunits
	if data_log: title = "Log[ "+title+" ]"
	
	
	# 4 - Loop times
	# -------------------------------------------------------------------
	if figure is not None: fig = plt.figure(figure)
	for itime in range(times.size):
		
		time = times[itime]
		A = np.double(f.values()[itime].values()[fieldn])
		
		# apply the slicing
		for iaxis in range(naxes):
			if slices[iaxis].size == 0: continue
			A = np.delete(A, slices[iaxis], axis=iaxis) # remove parts outside of the slice
			A = np.sum(A, axis=iaxis, keepdims=True) # sum over the slice
		A = np.squeeze(A) # remove sliced axes
		A *= unitscoeff
		
		# log scale if requested
		if data_log: A = np.log10(A)
		
		if figure is None: break
		
		# plot
		if A.ndim == 1:
			fig.clf()
			ax = fig.add_subplot(1,1,1)
			ax.plot(plot_centers[0], A)
			ax.set_xlabel(plot_label[0])
			ax.set_xlim(xmin=xmin, xmax=xmax)
			ax.set_ylim(ymin=data_min, ymax=data_max)
			ax.set_title(title)
			fig.canvas.draw()
			plt.show()
		
		elif A.ndim == 2:
			fig.clf()
			ax = fig.add_subplot(1,1,1)
			extent = [plot_centers[0][0], plot_centers[0][-1], plot_centers[1][0], plot_centers[1][-1]]
			im = ax.imshow( np.flipud(A.transpose()),
				vmin = data_min, vmax = data_max, extent=extent,
				aspect="auto", interpolation="nearest")
			ax.set_xlabel(plot_label[0])
			ax.set_ylabel(plot_label[1])
			ax.set_xlim(xmin=xmin, xmax=xmax)
			ax.set_ylim(ymin=ymin, ymax=ymax)
			plt.colorbar(im)
			ax.set_title(title)
			fig.canvas.draw()
			plt.show()
		
		if figure is not None: print "timestep "+str(time)
	
	
	f.close()
	
	if figure is None:
		result = {"data":A}
		for i in range(naxes):
			label = {0:"x", 1:"y", 2:"z"}[i]
			result.update({ label:plot_centers[i] })
		return result


	
	
