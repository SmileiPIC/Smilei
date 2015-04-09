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
# >>>>>> Second step: in the python shell, use the class "ParticleDiagnostic"
#                     to create a ParticleDiagnostic object
#
# ParticleDiagnostic(results_path, diagNumber=None, timesteps=None, slice=None,
#                    units="code", data_log=False):
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
#
# >>>>>> Third step: To plot the data, use the following method.
#
# ParticleDiagnostic(..., figure=1, data_min=None, data_max=None,
#                         xmin=None, xmax=None, ymin=None, ymax=None) .plot()
#
#            figure = _int_       (optional)
#                     The figure number that is passed to matplotlib.
#                     If absent, figure 1 is used.
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
#
# >>>>>> Instead of plotting, you can obtain the data as an array.
#        Note that only the first requested timestep is given.
#
# ParticleDiagnostic(...).getData()
#   This method returns only the data array.
#
# ParticleDiagnostic(...).get()
#   This method returns a dictionary containing the data, and the axes scales.
#
#
# >>>>>> To simultaneously plot multiple diagnostics in the same figure:
#
# multiPlot(diag1, diag2, ... , figure=1, shape=None)
#
#              diag1 = diagnostic prepared by ParticleDiagnostic(...)
#              diag2 = diagnostic prepared by ParticleDiagnostic(...)
#                      ...
#
#            figure = _int_             (optional)
#                     The figure number that is passed to matplotlib.
#                     If absent, figure 1 is used.
#
#             shape = [_int_ , _int_]   (optional)
#                     The arrangement of plots inside the figure.
#                     For instance, [2, 1] makes two plots stacked vertically,
#                      and [1, 2] makes two plots stacked horizontally.
#                     If absent, stacks plots vertically.
#
# >>>>>> Examples:
#    ParticleDiagnostic('../test', diagNumber=1, slice={"y":"all"}, units="nice").plot(figure=1, data_min=0, data_max=3e14)



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
# Main class
# -------------------------------------------------------------------
class ParticleDiagnostic(object):

	# We use __new__ to prevent object creation if no diagNumber requested
	def __new__(cls, *args, **kwargs):
		# Is there a "results_path" argument ?
		try:
			results_path = args[0]
		except:
			try:
				results_path = kwargs["results_path"]
			except:
				print "Argument 'results_path' required"
				return None
		if not os.path.isdir(results_path):
			print "Could not find directory "+results_path
			return None
		if len(glob.glob(results_path+"/*.in"))==0:
			print "Could not find an input file in directory "+results_path
			return None
		if len(glob.glob(results_path+"/*.in"))>1:
			print "Directory "+results_path+" contains more than one input file. There should be only one."
			return None
		# Is there a "diagNumber" argument ?
		try:
			diagNumber = args[1]
		except:
			try:
				diagNumber = kwargs["diagNumber"]
			except:
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
				return None
		# If everything ok, then we create the object
		return object.__new__(cls, *args, **kwargs)
	
	
	# This is the constructor, which creates the object
	def __init__(self, results_path, diagNumber=None, timesteps=None, slice=None,
	                    units="code", data_log=False, **kwargs):
		
		self.valid = False
		
		# Get info from the input file and prepare units
		try:
			dim = findParam(results_path, "dim")
			ndim = int(dim[0])
		except:
			print "Could not extract 'dim' from the input file"
			return None
		try:
			sim_units  = findParam(results_path, "sim_units")
		except:
			print "Could not extract 'sim_units' from the input file"
			return None
		try:
			sim_length = findParam(results_path, "sim_length")
			sim_length = np.double(sim_length.split())
		except:
			print "Could not extract 'sim_length' from the input file"
			return None
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
				return None
		if   ndim == 1:
			sim_length  = sim_length [0]
			cell_length = cell_length[0]
		elif ndim == 2:
			if sim_length.size  == 1: sim_length  = np.array([sim_length,sim_length])
			else                    : sim_length  = sim_length[0:2]
			if cell_length.size == 1: cell_length = np.array([cell_length,cell_length])
			else                    : cell_length = cell_length[0:2]
		elif ndim == 3:
			if sim_length.size == 1: sim_length = np.array([sim_length,sim_length,sim_length])
			elif sim_length.size >2: sim_length = sim_length[0:3]
			else:
				print "In the input file, 'sim_length' should have 1 or 3 arguments for a 3d simulation"
				return None
			if cell_length.size == 1: cell_length = np.array([cell_length,cell_length,cell_length])
			elif cell_length.size >2: cell_length = cell_length[0:3]
			else:
				print "In the input file, 'cell_length' or 'res_space' should have 1 or 3 arguments for a 3d simulation"
				return None
		else:
			print "Could not understand simulation dimension 'dim="+dim+"' from the input file"
			return None
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
				return None
			coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 # nc in cm^-3
			coeff_energy = 0.511
		elif units == "code":
			coeff_density = 1.
			coeff_energy = 1.
		else:
			print "Units type '"+units+"' not recognized. Use 'code' or 'nice'"
			return None
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check value of diagNumber
		if type(diagNumber) is not int or diagNumber<0:
			print "Argument 'diagNumber' must be a positive or zero integer"
			return None
		
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			print "Argument 'slice' must be a dictionary"
			return None
		
		# Get the info on the requested diagNumber
		self.info = getInfo(results_path, diagNumber)
		
		# Leave if diag not found
		if self.info == False:
			print "Particle diagnostic #"+str(diagNumber)+" not found";
			return None
		
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self.data_log = data_log
		
		# Default some other parameters
		self.figure   = 1
		self.data_min = None
		self.data_max = None
		self.xmin     = None
		self.xmax     = None
		self.ymin     = None
		self.ymax     = None
		self.figurekwargs = {}
		self.plotkwargs = {}
		self.imkwargs = {"interpolation":"nearest", "aspect":"auto"}
		self.colorbarkwargs = {}
		kwargs = self.setPlot(**kwargs)
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# Get available timesteps
		self.file = results_path+'/ParticleDiagnostic'+str(diagNumber)+'.h5'
		f = h5py.File(self.file, 'r')
		items = f.items()
		ntimes = len(items)
		self.times = np.zeros(ntimes)
		self.data  = {}
		for i in range(ntimes):
			self.times[i] = int(items[i][0].strip("timestep")) # fill the "times" array with the available timesteps
			self.data.update({ self.times[i] : i }) # fill the "data" dictionary with indices to the data arrays
		f.close()
	
		# If timesteps is None, then keep all timesteps
		# otherwise, select timesteps
		if timesteps is not None:
			try:
				ts = np.array(np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self.times = self.times[ (self.times>=ts[0]) * (self.times<=ts[1]) ]
				elif ts.size==1:
					# get nearest time
					self.times = np.array([self.times[(np.abs(self.times-ts)).argmin()]])
				else:
					raise Exception()
			except:
				print "Argument 'timesteps' must be one or two non-negative integers"
				f.close()
				return None
		
		# Need at least one timestep
		if self.times.size < 1:
			print "Timesteps not found"
			f.close()
			return None
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values and slices
		self.axes = self.info["axes"]
		self.naxes = len(self.axes)
		self.shape = []
		self.plot_shape = []; self.plot_type = []; plot_diff = []
		self.plot_label = []; self.plot_centers = []; self.plot_log = []
		units_coeff = 1.
		unitsa = [0,0,0,0]
		spatialaxes = {"x":False, "y":False, "z":False}
		for axis in self.axes:
			self.shape.append(axis["size"])
			
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
						return None
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
				if indices.size == 1:
					axis.update({ "sliceInfo" : "      Slicing at "+axis["type"]+" = "+str(centers[indices][0]) })
				else:
					axis.update({ "sliceInfo" : "      Slicing "+axis["type"]+" from "+str(edges[indices[0]])+" to "+str(edges[indices[-1]+1]) })
				
				# convert the range of indices into their "conjugate"
				indices = np.delete(np.arange(axis["size"]), indices)
				# put the slice in the dictionary
				axis.update({"slice":indices, "slice_size":slice_size})
				
				if axis["type"] in ["x","y","z"]:
					units_coeff *= cell_size[axis["type"]]/slice_size
			
			# if not sliced, then add this axis to the overall plot
			else:
				self.plot_type.append(axis["type"])
				self.plot_shape.append(axis["size"])
				self.plot_centers.append(centers*axis_coeff)
				self.plot_log.append(axis["log"])
				self.plot_label.append(axis["type"]+axis_units)
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
		
		
		if len(self.plot_shape) > 2:
			print "Cannot plot in "+str(len(self.plot_shape))+"d. You need to 'slice' some axes."
			return None
		
		# Build units
		units_coeff *= coeff_density
		self.title = "??"
		unitss = "??"
		if   self.info["output"] == "density":               self.title = "Number density"
		elif self.info["output"] == "charge_density":        self.title = "Charge density"
		elif self.info["output"][:-1] == "current_density_": self.title = "J"+self.info["output"][-1]
		if units == "nice":
			if   self.info["output"] == "density":               unitss = "particles/cm$^3$"
			elif self.info["output"] == "charge_density":        unitss = "e/cm$^3$"
			elif self.info["output"][:-1] == "current_density_": unitss = "particles * c /cm$^3$"
			if unitsa[1]>0: unitss += "/(mc)"
			if unitsa[1]>1: unitss += "$^"+str(unitsa[1])+"$"
			if unitsa[2]>0: unitss += "/c"
			if unitsa[2]>1: unitss += "$^"+str(unitsa[2])+"$"
			if unitsa[3]>0: unitss += "/MeV"
			if unitsa[3]>1: unitss += "$^"+str(unitsa[3])+"$"
		elif units == "code":
			if   self.info["output"] == "density":               unitss = "$n_c$"
			elif self.info["output"] == "charge_density":        unitss = "e$\, n_c$"
			elif self.info["output"][:-1] == "current_density_": unitss = "particles * $c\, n_c$"
			if unitsa[1]>0: unitss += "/(mc)"
			if unitsa[1]>1: unitss += "$^"+str(unitsa[1])+"$"
			if unitsa[2]>0: unitss += "/c"
			if unitsa[2]>1: unitss += "$^"+str(unitsa[2])+"$"
			if unitsa[3]>0: unitss += "/(mc$^2$)"
			if unitsa[3]>1: unitss += "$^"+str(unitsa[3])+"$"
		self.title += " [ "+unitss+" ]"
		if self.data_log: self.title = "Log[ "+self.title+" ]"
		
		# If any spatial dimension did not appear, then count it for calculating the correct density
		if ndim>=1 and not spatialaxes["x"]: units_coeff /= ncels[0]
		if ndim>=2 and not spatialaxes["y"]: units_coeff /= ncels[1]
		if ndim==3 and not spatialaxes["z"]: units_coeff /= ncels[2]
		
		# Calculate the array that represents the bins sizes in order to get units right.
		# This array will be the same size as the plotted array
		if len(plot_diff)==1:
			self.bsize = plot_diff[0]
		else:
			self.bsize = np.prod( np.array( np.meshgrid( *plot_diff ) ), axis=0)
			self.bsize = self.bsize.transpose()
		self.bsize /= units_coeff
		
		# Finish constructor
		self.valid = True
		
	
	# Method to verify everything was ok during initialization
	def validate(self):
		if not self.valid:
			print "Diagnostic is invalid"
			return False
		return True
	
	# Method to print info on the current diag
	def print_info(self):
		if not self.validate(): return
		printInfo(self.info)
		return
		
	# Method to set optional plotting arguments
	def setPlot(self, **kwargs):
		# First, we manage the main optional arguments
		self.figure   = kwargs.pop("figure"  ,self.figure  )
		self.data_min = kwargs.pop("data_min",self.data_min)
		self.data_max = kwargs.pop("data_max",self.data_max)
		self.xmin     = kwargs.pop("xmin"    ,self.xmin    )
		self.xmax     = kwargs.pop("xmax"    ,self.xmax    )
		self.ymin     = kwargs.pop("ymin"    ,self.ymin    )
		self.ymax     = kwargs.pop("ymax"    ,self.ymax    )
		# Second, we manage all the other arguments that are directly the ones
		#  of matplotlib
		# For each keyword argument provided, we save these arguments separately
		#  depending on the type: figure, plot, image, colorbar.
		for kwa in kwargs:
			val = kwargs[kwa]
			if kwa in ["figsize","dpi","facecolor","edgecolor"]:
				self.figurekwargs.update({kwa:val})
			if kwa in ["color","dashes","drawstyle","fillstyle","label","linestyle",
			           "linewidth","marker","markeredgecolor","markeredgewidth",
			           "markerfacecolor","markerfacecoloralt","markersize","markevery",
			           "visible","zorder"]:
				self.plotkwargs.update({kwa:val})
			if kwa in ["cmap","aspect","interpolation"]:
				self.imkwargs.update({kwa:val})
			if kwa in ["orientation","fraction","pad","shrink","anchor","panchor",
			           "extend","extendfrac","extendrect","spacing","ticks","format",
			           "drawedges"]:
				self.colorbarkwargs.update({kwa:val})
		# special case: "aspect" is ambiguous because it exists for both imshow and colorbar
		if "cbaspect" in kwargs: self.colorbarkwargs.update({"aspect":kwargs["cbaspect"]})
		return kwargs
		
	# Same Method, without returns
	def set(self, **kwargs):
		kwargs = self.setPlot(**kwargs)
	
	
	# Method to obtain the data only
	def getData(self, time=0):
		if not self.validate(): return
		if time not in self.data:
			print "Time "+time+" not found in this diagnostic"
			return []
		# Open file
		f = h5py.File(self.file, 'r')
		# get data
		A = np.reshape(f.items()[self.data[time]][1],self.shape)
		f.close()
		return A
	
	
	# Method to obtain the data and the axes
	def get(self, time=0):
		if not self.validate(): return
		print "timestep "+str(self.times[time])
		A = self.getData(time)		
		result = {"data":A}
		for i in range(len(self.plot_type)):
			result.update({ self.plot_type[i]:self.plot_centers[i] })
		return result
	
	# Method to obtain the plot limits
	def limits(self):
		l = []
		for i in range(len(self.plot_shape)):
			l.append([min(self.plot_centers[0]), max(self.plot_centers[0])])
		return l
	
	# Method to plot the data when axes are made
	def plotOnAxes(self, ax, time):
		if not self.validate(): return None
		# get data
		A = self.getData(time)
		# apply the slicing
		for iaxis in range(self.naxes):
			axis = self.axes[iaxis]
			if "slice" in axis:
				A = np.delete(A, axis["slice"], axis=iaxis) # remove parts outside of the slice
				A = np.sum(A, axis=iaxis, keepdims=True) # sum over the slice
		A = np.squeeze(A) # remove sliced axes
		# Divide by the bins size
		A /= self.bsize
		# log scale if requested
		if self.data_log: A = np.log10(A)
		# plot
		if A.ndim == 1:
			im, = ax.plot(self.plot_centers[0], A, **self.plotkwargs)
			if self.plot_log[0]: ax.set_xscale("log")
			ax.set_xlabel(self.plot_label[0])
			ax.set_xlim(xmin=self.xmin, xmax=self.xmax)
			ax.set_ylim(ymin=self.data_min, ymax=self.data_max)
			ax.set_title(self.title)
		elif A.ndim == 2:
			extent = [self.plot_centers[0][0], self.plot_centers[0][-1], self.plot_centers[1][0], self.plot_centers[1][-1]]
			if self.plot_log[0]: extent[0:2] = [np.log10(self.plot_centers[0][0]), np.log10(self.plot_centers[0][-1])]
			if self.plot_log[1]: extent[2:4] = [np.log10(self.plot_centers[1][0]), np.log10(self.plot_centers[1][-1])]
			im = ax.imshow( np.flipud(A.transpose()),
				vmin = self.data_min, vmax = self.data_max, extent=extent, **self.imkwargs)
			if (self.plot_log[0]): ax.set_xlabel("Log[ "+self.plot_label[0]+" ]")
			else:                  ax.set_xlabel(        self.plot_label[0]     )
			if (self.plot_log[1]): ax.set_ylabel("Log[ "+self.plot_label[1]+" ]")
			else:                  ax.set_ylabel(        self.plot_label[1]     )
			ax.set_xlim(xmin=self.xmin, xmax=self.xmax)
			ax.set_ylim(ymin=self.ymin, ymax=self.ymax)
			try: # if colorbar exists
				ax.cax.cla()
				plt.colorbar(mappable=im, cax=ax.cax, **self.colorbarkwargs)
			except AttributeError:
				ax.cax = plt.colorbar(mappable=im, ax=ax, **self.colorbarkwargs).ax
			ax.set_title(self.title)
		return im
	
	
	# Method to plot the current diagnostic
	def plot(self, **kwargs):
		if not self.validate(): return
		
		kwargs = self.setPlot(**kwargs)
		self.print_info()
		for ax in self.axes:
			if "sliceInfo" in ax: print ax["sliceInfo"]
		
		# Loop times
		fig = plt.figure(self.figure)
		fig.set(**self.figurekwargs)
		fig.clf()
		for timeindex in range(self.times.size):
			time = self.times[timeindex]
			print "timestep "+str(time)
			# plot
			ax = fig.add_subplot(1,1,1)
			ax.cla()
			artist = self.plotOnAxes(ax, time)
			fig.canvas.draw()
			plt.show()
			

# Function to plot multiple diags on the same figure
def multiPlot(*diags, **kwargs):
	ndiags = len(diags)
	# Verify diags are valid
	if ndiags == 0: return
	for diag in diags:
		if not diag.validate(): return
	# Gather all times
	alltimes = np.unique(np.concatenate([diag.times for diag in diags]))
	# Get keyword arguments
	figure = kwargs.pop("figure", 1)
	shape  = kwargs.pop("shape" , None)
	# Determine whether to plot all cases on the same axes
	sameAxes = False
	if shape is None or shape == [1,1]:
		sameAxes = True
		for diag in diags:
			if len(diag.plot_shape)!=1 or diag.plot_type!=diags[0].plot_type:
				sameAxes = False
				break
	if not sameAxes and shape == [1,1]:
		print "Cannot have shape=[1,1] with these diagnostics"
		return
	# Determine the shape
	if sameAxes: shape = [1,1]
	if shape is None: shape = [ndiags,1]
	nplots = np.array(shape).prod()
	if not sameAxes and nplots != ndiags:
		print "The 'shape' argument is incompatible with the number of diagnostics:"
		print "  "+str(ndiags)+" diagnostics do not fit "+str(nplots)+" plots"
		return
	# Make the figure
	fig = plt.figure(figure)
	fig.set(**diags[0].figurekwargs)
	fig.clf()
	fig.subplots_adjust(wspace=0.5, hspace=0.5, bottom=0.15)
	ax = []
	xmin =  float("inf")
	xmax = -float("inf")
	c = plt.matplotlib.rcParams['axes.color_cycle']
	for i in range(nplots):
		ax.append( fig.add_subplot(shape[0], shape[1], i) )
	for i, diag in enumerate(diags):
		if sameAxes: diag.ax = ax[0]
		else       : diag.ax = ax[i]
		diag.artist = None
		l = diag.limits()[0]
		xmin = min(xmin,l[0])
		xmax = max(xmax,l[1])
		if "color" not in diag.plotkwargs: diag.plotkwargs.update({ "color":c[i%len(c)] })
	# Loop all times
	for time in alltimes:
		for diag in diags:
			if time in diag.times:
				if sameAxes:
					if diag.artist is not None: diag.artist.remove()
				else:
					diag.ax.cla()
				diag.artist = diag.plotOnAxes(diag.ax, time)
				if sameAxes:
					diag.ax.set_xlim(xmin,xmax)
		fig.canvas.draw()
		plt.show()
	return
	

