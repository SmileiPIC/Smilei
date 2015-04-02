# -----------------------------------------------------------------------
# HOW TO VIEW PARTICLE DIAGNOSTICS                -    F. Perez - 03/2015
# -----------------------------------------------------------------------
#
# >>>>>> What can be done
#   During the simulation, each particle diagnostic collects the data from particles
#   into a N-dimensional histogram.
#   Each histogram axis can be: x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge.
#   In each bin of the histogram, several things may be summed: the weights (density), 
#     weight*charge (charge density) or weight*velocity (current density).
#   Examples:
#       +----------+------------+---------------------+
#       |   Rank   |   type     |         Axes        |
#       +----------+------------+---------------------+
#       |   1-D    |  density   |        'ekin'       | => energy distribution.
#       |   2-D    |  density   |     'x' and 'y'     | => density map.
#       |   2-D    |  x-current |     'x' and 'y'     | => x-current map.
#       |   2-D    |  density   |     'x' and 'px'    | => phase space.
#       |   3-D    |  density   | 'x', 'y' and 'ekin' | => density map for several energy ranges.
#       +----------+------------+---------------------+
#
# >>>>>> Requirements
#   python2.7 with the following packages: numpy, matplotlib, pylab, h5py
#
# >>>>>> First step: invoke python and load this file
#      $ python -i ParticleDiagnostic.py
#
# >>>>>> Second step: in the python shell, use the function "ParticleDiagnostic"
#
# ParticleDiagnostic(results_path, diagNumber=None, timesteps=None, slice=None,
#                    units="code", data_log=False, data_min=None, data_max=None,
#                    xmin=None, xmax=None, ymin=None, ymax=None,
#                    figure=None)
#
#      results_path = _string_
#                     Path to the directory where the outputs are stored
#                     (Also, this has to contain one and only one input file *.in)
#
#        diagNumber = _int_   (optional)
#                     Number of the diagnostic. `0` is the first diagnostic.
#                     If not given, then a list of available diagnostics is printed.
#
#         timesteps = _int_            (optional)
#         timesteps = [_int_, _int_]   (optional)
#                     If omitted, all timesteps are used.
#                     If one number  given, the nearest timestep available is used.
#                     If two numbers given, all the timesteps in between are used.
#
#             slice = { "axis" : "all", ... }                 (optional)
#             slice = { "axis" : _double_, ... }              (optional)
#             slice = { "axis" : [_double_, _double_], ... }  (optional)
#                     This parameter is used to reduce the number of dimensions of the array.
#                     If the `axis` key is present, then any axis of the same name will be removed.
#                     `axis` must be x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge.
#                      - If the value is "all", then a sum is performed over all the axis.
#                      - If the value is _double_, then only the bin closest to the value is kept.
#                      - If the value is [_double_,_double_], then a sum is performed between the two values.
#                     Example: {"x":[4,5]} will sum all the data for x in the range [4,5].
#
#             units = ("code") or "nice"    (optional)
#                     If "nice" is chosen, then units are converted into usual units.
#                     Distances in microns, density in 1/cm^3, energy in MeV.
#
#          data_log = True or (False)       (optional)
#                     If True, then log10 is applied to the output array before plotting.
#
#          data_min = _double_    (optional)
#          data_max = _double_    (optional)
#                     If present, output is rescaled before plotting.
#
#              xmin = _double_    (optional)
#              xmax = _double_    (optional)
#              ymin = _double_    (optional)
#              ymax = _double_    (optional)
#                     If present, axes are rescaled before plotting.
#
#            figure = None        (default)
#            figure = _int_       (optional)
#                     Choses the figure number that is passed to matplotlib.
#                     If absent or None, returns the first data without plotting.
#
# >>>>>> Examples:
#    ParticleDiagnostic('../test', diagNumber=1, slice={"y":"all"}, units="nice",data_min=0, data_max=3e14, figure=1)



import h5py
import numpy as np
import os.path, glob, re
import matplotlib.pyplot as plt
import pylab
pylab.ion()


# Prints the info obtained by the function "getInfo"
def printInfo(info):
	if info==False: return
	
	# 1 - diag number, type and list of species
	species = ""
	for i in range(len(info["species"])): species += str(info["species"][i])+" " # reconstitute species string
	print "Diag#"+str(info["#"])+" - "+info["output"]+" of species # "+species
	
	# 2 - period and time-averaging
	tavg = "no time-averaging"
	if (info["tavg"] > 1):
		tavg = "averaging over "+str(info["tavg"])+" timesteps"
	print "    Every "+ str(info["every"]) + " timesteps, "+tavg
	
	# 3 - axes
	for i in range(len(info["axes"])):
		axis = info["axes"][i];
		logscale = "" if not axis["log"] else " [ LOG SCALE ] "
		edges    = "" if not axis["edges_included"] else " [ INCLUDING EDGES ] "
		print ("    "+axis["type"]+" from "+str(axis["min"])+" to "+str(axis["max"])
		       +" in "+str(axis["size"])+" steps "+logscale+edges)
		       
	return True


# Gets info about diagnostic number "diagNumber" in the path "results_path"
def getInfo(results_path, diagNumber):	
	# path to the file
	file = results_path+'/ParticleDiagnostic'+str(diagNumber)+'.h5'
	# if no file, return
	if not os.path.isfile(file): return False
	# open file
	f = h5py.File(file, 'r')
	# get attributes from file
	attrs = f.attrs.items()
	axes = []
	# Parse each attribute
	for i in range(len(attrs)):
		name  = attrs[i][0]
		value = attrs[i][1]
		if (name == "output"): output = value
		if (name == "every" ): every  = int(value)
		if (name == "time_average"): time_average = int(value)
		if (name == "species"):
		    species = value.strip().split(" ") # get all species numbers
		    for i in range(len(species)):
		    	species[i] = int(species[i]) # convert string to int
		if (name[0:4] == "axis" ):
			n = int(name[4:]) # axis number
			sp = value.split(" ")
			axistype  = sp[0]
			axismin  = float(sp[1])
			axismax  = float(sp[2])
			axissize = int(sp[3])
			logscale = bool(int(sp[4]))
			edge_inclusive = bool(int(sp[5]))
			while len(axes)<n+1: axes.append({}) # extend the array to the adequate size
			axes[n] = {"type":axistype,"min":axismin,"max":axismax,"size":axissize,"log":logscale,"edges_included":edge_inclusive}
	f.close()
	
	info = {"#":diagNumber, "output":output, "every":every, "tavg":time_average, "species":species, "axes":axes}
	return info


# Finds a parameter "param" in the input file
# Argument "after" is a string that must be found before "param"
def findParam(results_path, param, after=None):
	out = ""
	ok = True if after is None else False
	file = glob.glob(results_path+"/*.in")[0]
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
	

# get all available timesteps for a given diagnostic
def getAvailableTimesteps(results_path, diagNumber):
	try:
		file = results_path+'/ParticleDiagnostic'+str(diagNumber)+'.h5'
		f = h5py.File(file, 'r')
	except:
		print "Cannot open file "+file
		return np.array([])
	items = f.items()
	ntimes = len(items)
	times = np.zeros(ntimes)
	for i in range(ntimes):
		times[i] = int(items[i][0].strip("timestep")) # fill the "times" array with the available timesteps
	f.close()
	return times
	

# -------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------
def ParticleDiagnostic(results_path, diagNumber=None, timesteps=None, slice=None,
                       units="code", data_log=False, data_min=None, data_max=None,
                       xmin=None, xmax=None, ymin=None, ymax=None,
                       figure=None):

	# If no diag requested, list available diags
	if diagNumber is None:
		print "Printing available particle diagnostics:"
		print "----------------------------------------"
		
		ok = True
		diagNumber = 0
		while True:
			ok = printInfo(getInfo(results_path, diagNumber))
			if not ok: break
			diagNumber += 1
		
		if diagNumber == 0:
			print "      No particle diagnostics found in "+results_path;
		
		return
	
	# Get info from the input file and prepare units
	try:
		dim = findParam(results_path, "dim")
		ndim = int(dim[0])
	except:
		print "Could not extract 'dim' from the input file"
		return
	try:
		sim_units  = findParam(results_path, "sim_units")
	except:
		print "Could not extract 'sim_units' from the input file"
		return
	try:
		sim_length = findParam(results_path, "sim_length")
		sim_length = np.double(sim_length.split())
	except:
		print "Could not extract 'sim_length' from the input file"
		return
	try:
		cell_length = findParam(results_path, "cell_length")
		cell_length = np.double(cell_length.split())
	except:
		try:
			res_space = findParam(results_path, "res_space")
			res_space = np.double(res_space.split())
			cell_length = sim_length/res_space
		except:
			print "Could not extract 'cell_length' or 'res_space' from the input file"
			return
	if   ndim == 1:
		sim_length  = sim_length [0]
		cell_length = cell_length[0]
	elif ndim == 2:
		if sim_length.size  == 1: sim_length  = np.array([sim_length,sim_length])
		else                    : sim_length  = sim_length[0:1]
		if cell_length.size == 1: cell_length = np.array([cell_length,cell_length])
		else                    : cell_length = cell_length[0:1]
	elif ndim == 3:
		if sim_length.size == 1: sim_length = np.array([sim_length,sim_length,sim_length])
		elif sim_length.size >2: sim_length = sim_length[0:2]
		else:
			print "In the input file, 'sim_length' should have 1 or 3 arguments for a 3d simulation"
			return
		if cell_length.size == 1: cell_length = np.array([cell_length,cell_length,cell_length])
		elif cell_length.size >2: cell_length = cell_length[0:2]
		else:
			print "In the input file, 'cell_length' or 'res_space' should have 1 or 3 arguments for a 3d simulation"
			return
	else:
		print "Could not understand simulation dimension 'dim="+dim+"' from the input file"
		return
	sim_length  = np.array(sim_length ,ndmin=1)
	cell_length = np.array(cell_length,ndmin=1)
	ncels = sim_length/cell_length
	if sim_units == "wavelength": cell_length *= 2.*np.pi
	cell_size = {"x":cell_length[0]}
	if ndim>1: cell_size.update({"y":cell_length[1]})
	if ndim>2: cell_size.update({"z":cell_length[2]})
	if units == "nice":
		try:
			wavelength_SI = float( findParam(results_path, "wavelength_SI") )
		except:
			print "Could not extract 'wavelength_SI' from the input file"
			return
		coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 # nc in cm^-3
		coeff_energy = 0.511
	elif units == "code":
		coeff_density = 1.
		coeff_energy = 1.
	else:
		print "Units type '"+units+"' not recognized. Use 'code' or 'nice'"
		return

	# 1 - verifications, initialization
	# -------------------------------------------------------------------
	# Check value of diagNumber
	if type(diagNumber) is not int or diagNumber<0:
		print "Argument 'diagNumber' must be a positive or zero integer"
		return
	
	# Check slice is a dict
	if slice is not None  and  type(slice) is not dict:
		print "Argument 'slice' must be a dictionary"
		return
	
	# Get the info on the requested diagNumber
	info = getInfo(results_path, diagNumber)
	
	# Leave if diag not found
	if info == False:
		print "Particle diagnostic #"+str(diagNumber)+" not found";
		return
		
	if figure is not None: printInfo(info)
	
	# Open hdf file
	file = results_path+'/ParticleDiagnostic'+str(diagNumber)+'.h5'
	f = h5py.File(file, 'r')
	
	# Make slice a dictionary
	if slice is None: slice = {}
	
	# 2 - Manage timesteps
	# -------------------------------------------------------------------
	# Get available timesteps
	items = f.items()
	ntimes = len(items)
	times = np.zeros(ntimes)
	data  = {}
	for i in range(ntimes):
		times[i] = int(items[i][0].strip("timestep")) # fill the "times" array with the available timesteps
		data.update({ times[i] : items[i][1] }) # fill the "data" dictionary with pointers to the data arrays

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
			print "Argument 'timesteps' must be one or two non-negative integers"
			f.close()
			return
	
	# Need at least one timestep
	if times.size < 1:
		print "Timesteps not found"
		f.close()
		return
	
		
	# 3 - Manage axes
	# -------------------------------------------------------------------
	# Fabricate all axes values and slices
	axes = info["axes"]
	naxes = len(axes)
	shape = []
	plot_shape = []; plot_type = []; plot_label = []; plot_centers = []; plot_log = []; plot_diff = []
	units_coeff = 1.
	unitsa = [0,0,0,0]
	spatialaxes = {"x":False, "y":False, "z":False}
	for iaxis in range(naxes):
		axis = axes[iaxis]
		shape.append(axis["size"])
		
		# Find the vector of values along the axis
		if axis["log"]:
			edges = np.linspace(np.log10(axis["min"]), np.log10(axis["max"]), axis["size"]+1)
			centers = edges + (edges[1]-edges[0])/2.
			centers = 10.**(centers[:-1])
		else:
			edges = np.linspace(axis["min"], axis["max"], axis["size"]+1)
			centers = edges + (edges[1]-edges[0])/2.
			centers = centers[:-1]
		axis.update({ "edges"   : edges   })
		axis.update({ "centers" : centers })
		
		# Find some quantities depending on the axis type
		overall_min = "-inf"; overall_max = "inf"
		axis_units = ""; axis_coeff = 1.
		if   axis["type"] in ["x","y","z"]:
			axis_units = " [ wavelength / 2Pi ]"
			if units == "nice":
				axis_units = " [ microns ]"
				axis_coeff = 1e6*wavelength_SI/(2.*np.pi)
			spatialaxes[axis["type"]] = True
		elif axis["type"] in ["px","py","pz","p"]:
			axis_units = " [ m c ]"
		elif axis["type"] in ["vx","vy","vz","v"]:
			axis_units = " [ c ]"
		elif axis["type"] == "gamma":
			overall_min = "1"
		elif axis["type"] == "ekin":
			axis_units = " [ m c^2 ]"
			if units == "nice":
				axis_units = " [ MeV ]"
				axis_coeff = 0.511
			overall_min = "0"
		elif axis["type"] == "charge":
			overall_min = "0"

		# if this axis has to be sliced, then select the slice
		if axis["type"] in slice:
					
			# if slice is "all", then all the axis has to be summed
			if slice[axis["type"]] == "all":
				indices = np.arange(axis["size"])
			
			# Otherwise, get the slice from the argument `slice`
			else:	
				try:
					s = np.double(slice[axis["type"]])
					if s.size>2 or s.size<1: raise Exception()
				except Exception(e):
					print "Slice along axis "+axis["type"]+" should be one or two floats"
					return
				# convert the slice into a range of indices
				if s.size == 1:
					indices = np.array([(np.abs(centers-s)).argmin()])
				else :
					indices = np.nonzero( (centers>=s[0]) * (centers<=s[1]) )[0]
			
			# calculate the size of the slice
			imin = indices.min()  ; emin = edges[imin]
			imax = indices.max()+1; emax = edges[imax]
			slice_size = emax - emin
			
			# print out the slicing
			if imin==0            and axis["edges_included"]: emin = overall_min
			if imin==axis["size"] and axis["edges_included"]: emax = overall_max
			if figure is not None:
				if indices.size == 1:
					print "   Slicing at "+axis["type"]+" = "+str(centers[indices])
				else:
					print "   Slicing at "+axis["type"]+" = "+str(edges[indices[0]])+" to "+str(edges[indices[-1]])
			
			# convert the range of indices into their "conjugate"
			indices = np.delete(np.arange(axis["size"]), indices)
			# put the slice in the dictionary
			axis.update({"slice":indices, "slice_size":slice_size})
			
			if axis["type"] in ["x","y","z"]:
				units_coeff *= cell_size[axis["type"]]/slice_size
		
		# if not sliced, then add this axis to the overall plot
		else:
			plot_shape.append(axis["size"])
			plot_type.append(axis["type"])
			plot_label.append(axis["type"]+axis_units)
			plot_centers.append(centers*axis_coeff)
			plot_log.append(axis["log"])
			plot_diff.append(np.diff(edges))
			if   axis["type"] in ["x","y","z"]:
				units_coeff *= cell_size[axis["type"]]
				unitsa[0] += 1
			elif axis["type"] in ["px","py","pz","p"]:
				unitsa[1] += 1
			elif axis["type"] in ["vx","vy","vz","v"]:
				unitsa[2] += 1
			elif axis["type"] == "ekin":
				units_coeff /= coeff_energy
				unitsa[3] += 1 
	

	if len(plot_shape) > 2:
		print "Cannot plot in "+str(len(plot_shape))+"d. You need to 'slice' some axes."
		return
	
	# Build units
	units_coeff *= coeff_density
	if   info["output"] == "density":               title = "Number density"
	elif info["output"] == "charge_density":        title = "Charge density"
	elif info["output"][:-1] == "current_density_": title = "J"+info["output"][-1]
	if units == "nice":
		if   info["output"] == "density":               unitss = "particles/cm$^3$"
		elif info["output"] == "charge_density":        unitss = "e/cm$^3$"
		elif info["output"][:-1] == "current_density_": unitss = "particles * c /cm$^3$"
		if unitsa[1]>0: unitss += "/(mc)"
		if unitsa[1]>1: unitss += "$^"+str(unitsa[1])+"$"
		if unitsa[2]>0: unitss += "/c"
		if unitsa[2]>1: unitss += "$^"+str(unitsa[2])+"$"
		if unitsa[3]>0: unitss += "/MeV"
		if unitsa[3]>1: unitss += "$^"+str(unitsa[3])+"$"
	elif units == "code":
		if   info["output"] == "density":               unitss = "$n_c$"
		elif info["output"] == "charge_density":        unitss = "e$\, n_c$"
		elif info["output"][:-1] == "current_density_": unitss = "particles * $c\, n_c$"
		if unitsa[1]>0: unitss += "/(mc)"
		if unitsa[1]>1: unitss += "$^"+str(unitsa[1])+"$"
		if unitsa[2]>0: unitss += "/c"
		if unitsa[2]>1: unitss += "$^"+str(unitsa[2])+"$"
		if unitsa[3]>0: unitss += "/(mc$^2$)"
		if unitsa[3]>1: unitss += "$^"+str(unitsa[3])+"$"
	title += " [ "+unitss+" ]"
	if data_log: title = "Log[ "+title+" ]"
	
	# If any spatial dimension did not appear, then count it for calculating the correct density
	if ndim>=1 and not spatialaxes["x"]: units_coeff /= ncels[0]
	if ndim>=2 and not spatialaxes["y"]: units_coeff /= ncels[1]
	if ndim==3 and not spatialaxes["z"]: units_coeff /= ncels[2]
	
	# Calculate the array that represents the bins sizes in order to get units right.
	# This array will be the same size as the plotted array
    if len(plot_diff)==1:
        bsize = np.prod( np.array( plot_diff ), axis=0)
    else:
        bsize = np.prod( np.array( np.meshgrid( *tuple(plot_diff) ) ), axis=0)    
	bsize /= units_coeff
	bsize = bsize.transpose()
	
	
	
	# 4 - Loop times
	# -------------------------------------------------------------------
	if figure is not None: fig = plt.figure(figure)
	for itime in range(times.size):
		
		time = times[itime]
		A = np.reshape(data[time],shape)
		
		# apply the slicing
		for iaxis in range(naxes):
			if "slice" in axes[iaxis]:
				A = np.delete(A, indices, axis=iaxis) # remove parts outside of the slice
				A = np.sum(A, axis=iaxis, keepdims=True) # sum over the slice
		A = np.squeeze(A) # remove sliced axes
		
		# Divide by the bins size
		A = A/bsize
		
		# log scale if requested
		if data_log: A = np.log10(A)
		
		if figure is None: break
		
		# plot
		if A.ndim == 1:
			fig.clf()
			ax = fig.add_subplot(1,1,1)
			ax.plot(plot_centers[0], A)
			if plot_log[0]: ax.set_xscale("log")
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
			if plot_log[0]: extent[0:2] = [np.log10(plot_centers[0][0]), np.log10(plot_centers[0][-1])]
			if plot_log[1]: extent[2:4] = [np.log10(plot_centers[1][0]), np.log10(plot_centers[1][-1])]
			im = ax.imshow( np.flipud(A.transpose()),
				vmin = data_min, vmax = data_max, extent=extent,
				aspect="auto", interpolation="nearest")
			if (plot_log[0]): ax.set_xlabel("Log[ "+plot_label[0]+" ]")
			else:             ax.set_xlabel(        plot_label[0]     )
			if (plot_log[1]): ax.set_ylabel("Log[ "+plot_label[1]+" ]")
			else:             ax.set_ylabel(        plot_label[1]     )
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
		for i in range(len(plot_type)):
			result.update({ plot_type[i]:plot_centers[i] })
		return result


	
	
