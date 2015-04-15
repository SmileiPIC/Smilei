# -----------------------------------------------------------------------
# HOW TO VIEW DIAGNOSTICS                -    F. Perez - 03/2015
# -----------------------------------------------------------------------
#  +---------------------------+
#  | 1 . Particle diagnostics  |
#  +---------------------------+
# >>>>>> What can be done
#   During the simulation, each particle diagnostic collects the data from particles
#   into a N-dimensional histogram.
#   Each histogram axis can be: x, y, z, px, py, pz, p, gamma, ekin, vx, vy, vz, v or charge.
#   In each bin of the histogram, several things may be summed: the weights (density), 
#     weight*charge (charge density), weight*charge*velocity (current density),
#     or weight*momentum (momentum density)
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
#      $ python -i Diagnostics.py
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
#
#
#  +---------------------------+
#  | 2 . Fields                |
#  +---------------------------+
#
# You can view fields with a similar procedure:
#
# Field(results_path, field=None, timesteps=None, slice=None,
#       units="code", data_log=False)
#
# This works almost the same way as ParticleDiagnostic(...), with some exceptions:
#
# - `field` must be one of "Bx_m", "By_m", "Bz_m", "Ex", "Ey", "Ez", "Jx", "Jy", "Jz",
#   "Jx_[species]", "Jy_[species]", "Jz_[species]", "Rho" or "Rho_[species]"
#   where [species] is the name of one of the existing species.
#   If you omit the argument `field`, the list of available fields will be displayed.
#
# - `slice` can only accept three axes: "x", "y", "z".
#   For instance, slice={"x":"all"}.
#   Note that the slice does NOT calculate the sum of the axis, but the AVERAGE.
#
#
#
#  +---------------------------+
#  | 3 . Scalars               |
#  +---------------------------+
# TODO


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

# Method to manage Matplotlib-compatible kwargs
def matplotlibArgs(kwargs):
	figurekwargs   = {}
	plotkwargs     = {}
	imkwargs       = {}
	colorbarkwargs = {}
	for kwa in kwargs:
		val = kwargs[kwa]
		if kwa in ["figsize","dpi","facecolor","edgecolor"]:
			figurekwargs.update({kwa:val})
		if kwa in ["color","dashes","drawstyle","fillstyle","label","linestyle",
		           "linewidth","marker","markeredgecolor","markeredgewidth",
		           "markerfacecolor","markerfacecoloralt","markersize","markevery",
		           "visible","zorder"]:
			plotkwargs.update({kwa:val})
		if kwa in ["cmap","aspect","interpolation"]:
			imkwargs.update({kwa:val})
		if kwa in ["orientation","fraction","pad","shrink","anchor","panchor",
		           "extend","extendfrac","extendrect","spacing","ticks","format",
		           "drawedges"]:
			colorbarkwargs.update({kwa:val})
	return figurekwargs, plotkwargs, imkwargs, colorbarkwargs


# -------------------------------------------------------------------
# Mother class for all diagnostics
# -------------------------------------------------------------------
class Diagnostic(object):


	# When no action is performed on the object, this is what appears
	def __repr__(self):
		self.info()
		return ""
	
	# Various methods to extract stuff from the input file
	def read_sim_units(self):
		try:
			return findParam(self.results_path, "sim_units")
		except:
			print "Could not extract 'sim_units' from the input file"
			raise
	def read_ndim(self):
		try:
			dim = findParam(self.results_path, "dim")
			ndim = int(dim[0])
		except:
			print "Could not extract 'dim' from the input file"
			raise
		if ndim not in [1,2,3]:
			print "Could not understand simulation dimension 'dim="+dim+"' from the input file"
			raise
		return ndim
	def read_ncels_cell_length(self, ndim, sim_units):
		try:
			sim_length = findParam(self.results_path, "sim_length")
			sim_length = np.double(sim_length.split())
		except:
			print "Could not extract 'sim_length' from the input file"
			raise
		try:
			cell_length = findParam(self.results_path, "cell_length")
			cell_length = np.double(cell_length.split())
		except:
			try:
				res_space = findParam(self.results_path, "res_space")
				res_space = np.double(res_space.split())
				cell_length = sim_length/res_space
			except:
				print "Could not extract 'cell_length' or 'res_space' from the input file"
				raise
		if   ndim == 1:
			sim_length  = sim_length[0]
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
				raise
			if cell_length.size == 1: cell_length = np.array([cell_length,cell_length,cell_length])
			elif cell_length.size >2: cell_length = cell_length[0:3]
			else:
				print "In the input file, 'cell_length' or 'res_space' should have 1 or 3 arguments for a 3d simulation"
				raise
		sim_length  = np.array(sim_length ,ndmin=1)
		cell_length = np.array(cell_length,ndmin=1)
		ncels = sim_length/cell_length
		if sim_units == "wavelength": cell_length *= 2.*np.pi
		return ncels, cell_length
	def read_timestep(self,sim_units):
		try:
			timestep = np.double(findParam(self.results_path, "timestep"))
		except:
			try:
				res_time = np.double(findParam(self.results_path, "res_time"))
				sim_time = np.double(findParam(self.results_path, "sim_time"))
				timestep = sim_time/res_time
			except:
				print "Could not extract 'timestep' or 'res_time' from the input file"
				raise
		if sim_units == "wavelength": timestep *= 2.*np.pi
		return timestep
	def read_wavelength_SI(self):
		try:
			wavelength_SI = np.double( findParam(self.results_path, "wavelength_SI") )
		except:
			print "Could not extract 'wavelength_SI' from the input file"
			raise
		return wavelength_SI

	# Method to verify everything was ok during initialization
	def validate(self):
		if not self.valid:
			print "Diagnostic is invalid"
			return False
		return True

	# Method to verify that results_path is valid
	@staticmethod
	def validatePath(*args, **kwargs):
		try:
			results_path = args[0]
			if "results_path" in kwargs:
				print "Too many arguments 'results_path'"
				raise
		except:
			try:
				results_path = kwargs["results_path"]
			except:
				print "Argument 'results_path' required"
				raise
		if not os.path.isdir(results_path):
			print "Could not find directory "+results_path
			raise
		if len(glob.glob(results_path+"/*.in"))==0:
			print "Could not find an input file in directory "+results_path
			raise
		if len(glob.glob(results_path+"/*.in"))>1:
			print "Directory "+results_path+" contains more than one input file. There should be only one."
			raise
		return results_path
	
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
		# Second, we manage all the other arguments that are directly the ones of matplotlib
		# For each keyword argument provided, we save these arguments separately
		#  depending on the type: figure, plot, image, colorbar.
		args = matplotlibArgs(kwargs)
		self.figurekwargs  .update(args[0])
		self.plotkwargs    .update(args[1])
		self.imkwargs      .update(args[2])
		self.colorbarkwargs.update(args[3])
		# special case: "aspect" is ambiguous because it exists for both imshow and colorbar
		if "cbaspect" in kwargs: self.colorbarkwargs.update({"aspect":kwargs["cbaspect"]})
		return kwargs
	def set(self, **kwargs):
		kwargs = self.setPlot(**kwargs)
	
	# Method to set the default plot parameters
	def setDefaultPlot(self):
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
	
	# Method to obtain the plot limits
	def limits(self):
		l = []
		for i in range(len(self.plot_shape)):
			l.append([min(self.plot_centers[0]), max(self.plot_centers[0])])
		return l
	
	# Method to obtain the data and the axes
	# Note: this is overloaded in the case of scalars
	def get(self, **kwargs):
		if not self.validate(): return
		# if optional argument "time" not provided, find out which time to plot
		try:    time = kwargs["time"]
		except: time = self.times[0]
		# obtain the data array
		A = self.getData(time=time)		
		# print timestep
		print "timestep "+str(time)+ "   -   t = "+str(time*self.coeff_time)+self.time_units
		# format the results into a dictionary
		result = {"data":A}
		for i in range(len(self.plot_type)):
			result.update({ self.plot_type[i]:self.plot_centers[i] })
		return result
	
	# Method to plot the current diagnostic
	def plot(self, **kwargs):
		if not self.validate(): return
		kwargs = self.setPlot(**kwargs)
		self.info()
		
		# Make figure
		fig = plt.figure(self.figure)
		fig.set(**self.figurekwargs)
		fig.clf()
		ax = fig.add_subplot(1,1,1)
		# Animation if several dimensions
		if len(self.plot_shape) > 0:
			# Loop times
			for timeindex in range(self.times.size):
				time = self.times[timeindex]
				print "timestep "+str(time)+ "   -   t = "+str(time*self.coeff_time)+self.time_units
				# plot
				ax.cla()
				artist = self.animateOnAxes(ax, time)
				fig.canvas.draw()
				plt.show()
		# Static plot if 0 dimensions
		else:
			ax.cla()
			artist = self.plotVsTime(ax)
	
	# Method to plot the data when axes are made
	def animateOnAxes(self, ax, time):
		if not self.validate(): return None
		# get data
		A = self.getData(time=time)
		# plot
		if A.ndim == 0: # as a function of time
			times = self.times[self.times<=time]
			A = np.zeros(times.size)
			for i, time in enumerate(times):
				A[i] = self.getData(time=time)
			im, = ax.plot(times*self.coeff_time, A, **self.plotkwargs)
			ax.set_xlabel('Time ['+self.time_units+' ]')
			ax.set_xlim(xmax=self.times[-1]*self.coeff_time)
			if self.data_min is not None: ax.set_ylim(ymin=self.data_min)
			if self.data_max is not None: ax.set_ylim(ymax=self.data_max)
		elif A.ndim == 1:
			im, = ax.plot(self.plot_centers[0], A, **self.plotkwargs)
			if self.plot_log[0]: ax.set_xscale("log")
			ax.set_xlabel(self.plot_label[0])
			if self.data_min is not None: ax.set_ylim(ymin=self.data_min)
			if self.data_max is not None: ax.set_ylim(ymax=self.data_max)
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
			if self.ymin is not None: ax.set_ylim(ymin=self.ymin)
			if self.ymax is not None: ax.set_ylim(ymax=self.ymax)
			try: # if colorbar exists
				ax.cax.cla()
				plt.colorbar(mappable=im, cax=ax.cax, **self.colorbarkwargs)
			except AttributeError:
				ax.cax = plt.colorbar(mappable=im, ax=ax, **self.colorbarkwargs).ax
		if self.xmin is not None: ax.set_xlim(xmin=self.xmin)
		if self.xmax is not None: ax.set_xlim(xmax=self.xmax)
		if self.title is not None: ax.set_title(self.title)
		return im
	
	# If the sliced data has 0 dimension, this function can plot it 
	def plotVsTime(self, ax):
		if len(self.plot_shape) > 0:
			print "To plot vs. time, it is necessary to slice all axes in order to obtain a 0-D array"
			return None
		# Loop times to gather data
		A = np.zeros(self.times.size)
		for i, time in enumerate(self.times):
			A[i] = self.getData(time=time)
		im, = ax.plot(self.times*self.coeff_time, A, **self.plotkwargs)
		ax.set_xlabel('Time ['+self.time_units+' ]')
		if self.xmin is not None: ax.set_xlim(xmin=self.xmin)
		if self.xmax is not None: ax.set_xlim(xmax=self.xmax)
		if self.data_min is not None: ax.set_ylim(ymin=self.data_min)
		if self.data_max is not None: ax.set_ylim(ymax=self.data_max)
		if self.title is not None: ax.set_title(self.title)
		return im
	




# -------------------------------------------------------------------
# Class for particle diagnostics
# -------------------------------------------------------------------
class ParticleDiagnostic(Diagnostic):
	# We use __new__ to prevent object creation if no diagNumber requested
	def __new__(cls, *args, **kwargs):
		# Is there a "results_path" argument ?
		try   : results_path = cls.validatePath(*args, **kwargs)
		except: return None
		# Is there a "diagNumber" argument ?
		try:
			diagNumber = args[1]
		except:
			try:
				diagNumber = kwargs["diagNumber"]
			except:
				print "Printing available particle diagnostics:"
				print "----------------------------------------"
				diagNumber = 0
				while cls.printInfo(cls.getInfo(results_path,diagNumber)):
					diagNumber += 1
				if diagNumber == 0:
					print "      No particle diagnostics found in "+results_path;
				return None
		# If everything ok, then we create the object
		return Diagnostic.__new__(cls, *args, **kwargs)
	
	
	# This is the constructor, which creates the object
	def __init__(self, results_path, diagNumber=None, timesteps=None, slice=None,
	             units="code", data_log=False, **kwargs):
		
		self.valid = False
		self.results_path = results_path
		
		# Get info from the input file and prepare units
		try:
			ndim               = self.read_ndim()
			sim_units          = self.read_sim_units()
			ncels, cell_length = self.read_ncels_cell_length(ndim, sim_units)
			timestep           = self.read_timestep(sim_units)
			cell_size = {"x":cell_length[0]}
			if ndim>1: cell_size.update({"y":cell_length[1]})
			if ndim>2: cell_size.update({"z":cell_length[2]})
		except:
			return None

		if units == "nice":
			try   : wavelength_SI = self.read_wavelength_SI()
			except: return None
			coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 # nc in cm^-3
			coeff_energy = 0.511
			self.coeff_time = timestep * wavelength_SI/3.e8 # in seconds
			self.time_units = " s"
		elif units == "code":
			coeff_density = 1.
			coeff_energy = 1.
			self.coeff_time = timestep
			self.time_units = " 1/w"
		else:
			print "Units type '"+units+"' not recognized. Use 'code' or 'nice'"
			return None
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check the requested diags are ok
		if type(diagNumber) is int:
			if diagNumber<0:
				print "Argument 'diagNumber' cannot be a negative integer."
				return None
			self.operation = '#' + str(diagNumber)
		elif type(diagNumber) is str:
			self.operation = diagNumber
		else:
			print "Argument 'diagNumber' must be and integer or a string."
			return None
			
		# Get list of requested diags
		self.diags = sorted(set([ int(d[1:]) for d in re.findall('#\d+',self.operation) ]))
		try:
			exec(re.sub('#\d+','1.',self.operation)) in None
		except:
			print "Cannot understand operation '"+self.operation+"'"
			return None
		# Verify that all requested diags exist and they all have the same shape
		self.info_ = {}
		self.shape = {}
		self.axes = {}
		self.naxes = {}
		for d in self.diags:
			self.getMyInfo(d)
			try:
				self.info_.update({ d:self.getMyInfo(d) })
			except:
				print "Particle diagnostic #"+str(d)+" not found."
				return None
			self.axes .update ({ d:self.info_[d]["axes"] })
			self.naxes.update ({ d:len(self.axes[d]) })
			self.shape.update({ d:[ axis["size"] for axis in self.axes[d] ] })
			if self.naxes[d] != self.naxes[self.diags[0]]:
				print "All diagnostics in operation '"+self.operation+"' must have as many axes."
				print (" Diagnotic #"+str(d)+" has "+str(self.naxes[d])+" axes and #"+
					str(self.diags[0])+" has "+str(self.naxes[self.diags[0]])+" axes")
				return None
			for a in self.axes[d]:
				if self.axes[d] != self.axes[self.diags[0]]:
					print ("In operation '"+self.operation+"', diagnostics #"+str(d)+" and #"
						+str(self.diags[0])+" must have the same shape.")
					return None
		
		self.axes  = self.axes [self.diags[0]]
		self.naxes = self.naxes[self.diags[0]]
		self.shape = self.shape[self.diags[0]]
		
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			print "Argument 'slice' must be a dictionary"
			return None
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self.data_log = data_log
		
		# Default plot parameters
		self.setDefaultPlot()
		# Set user's plot parameters
		kwargs = self.setPlot(**kwargs)
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# Get available timesteps
		self.file = {}
		self.times = {}
		self.data  = {}
		for d in self.diags:
			# get all diagnostics files
			self.file.update({ d:results_path+'/ParticleDiagnostic'+str(d)+'.h5' })
			# get all diagnostics timesteps
			self.times.update({ d:self.getAvailableTimesteps(d) })
			# fill the "data" dictionary with indices to the data arrays
			self.data.update({ d:{} })
			for i,t in enumerate(self.times[d]):
				self.data[d].update({ t : i })
			# If timesteps is None, then keep all timesteps, otherwise, select timesteps
			if timesteps is not None:
				try:
					ts = np.array(np.double(timesteps),ndmin=1)
					if ts.size==2:
						# get all times in between bounds
						self.times[d] = self.times[d][ (self.times>=ts[0]) * (self.times<=ts[1]) ]
					elif ts.size==1:
						# get nearest time
						self.times[d] = np.array([self.times[d][(np.abs(self.times[d]-ts)).argmin()]])
					else:
						raise
				except:
					print "Argument 'timesteps' must be one or two non-negative integers"
					return None
			# Verify that timesteps are the same for all diagnostics
			if (self.times[d] != self.times[self.diags[0]]).any() :
				print "All diagnostics in operation '"+self.operation+"' must have the same timesteps."
				print (" Diagnotic #"+str(d)+" has "+str(len(self.times[d]))+ " timesteps and #"
					+str(self.diags[0])+" has "+str(len(self.times[self.diags[0]])))+ " timesteps"
				return None
		# Now we need to keep only one array of timesteps because they should be all the same
		self.times = self.times[self.diags[0]]
		
		# Need at least one timestep
		if self.times.size < 1:
			print "Timesteps not found"
			return None
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values for all diags
		self.plot_shape = []; self.plot_type = []; plot_diff = []
		self.plot_label = []; self.plot_centers = []; self.plot_log = []
		units_coeff = 1.
		unitsa = [0,0,0,0]
		spatialaxes = {"x":False, "y":False, "z":False}
		for axis in self.axes:
			
			# Find the vector of values along the axis
			if axis["log"]:
				edges = np.linspace(np.log10(axis["min"]), np.log10(axis["max"]), axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				edges = 10.**edges
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
						if s.size>2 or s.size<1: raise
					except:
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
				self.plot_type   .append(axis["type"])
				self.plot_shape  .append(axis["size"])
				self.plot_centers.append(centers*axis_coeff)
				self.plot_log    .append(axis["log"])
				self.plot_label  .append(axis["type"]+axis_units)
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
		self.titles = {}
		self.units = {}
		for d in self.diags:
			self.titles.update({ d:"??" })
			unitss = "??"
			output = self.info_[d]["output"]
			if   output == "density"                        : self.titles[d] = "Number density"
			elif output == "charge_density"                 : self.titles[d] = "Charge density"
			elif output[:-1] == "current_density_"          : self.titles[d] = "J"+output[-1]
			elif output == "p_density"                      : self.titles[d] = "P density"
			elif output[2:] == "_density" and output[0]=="p": self.titles[d] = "P"+output[1]+" density"
			if units == "nice":
				if   output == "density"                        : unitss = "particles/cm$^3$"
				elif output == "charge_density"                 : unitss = "$e$/cm$^3$"
				elif output[:-1] == "current_density_"          : unitss = "particles * $c$ /cm$^3$"
				elif output == "p_density"                      : unitss = "particles * $m\,c$ /cm$^3$"
				elif output[2:] == "_density" and output[0]=="p": unitss = "particles * $m\,c$ /cm$^3$"
				if unitsa[1]>0: unitss += "/(mc)"
				if unitsa[1]>1: unitss += "$^"+str(unitsa[1])+"$"
				if unitsa[2]>0: unitss += "/c"
				if unitsa[2]>1: unitss += "$^"+str(unitsa[2])+"$"
				if unitsa[3]>0: unitss += "/MeV"
				if unitsa[3]>1: unitss += "$^"+str(unitsa[3])+"$"
			elif units == "code":
				if   output == "density"                        : unitss = "$n_c$"
				elif output == "charge_density"                 : unitss = "$e\, n_c$"
				elif output[:-1] == "current_density_"          : unitss = "particles * $c\, n_c$"
				elif output == "p_density"                      : unitss = "particles * $m\,c\, n_c$"
				elif output[2:] == "_density" and output[0]=="p": unitss = "particles * $m\,c\, n_c$"
				if unitsa[1]>0: unitss += "/(mc)"
				if unitsa[1]>1: unitss += "$^"+str(unitsa[1])+"$"
				if unitsa[2]>0: unitss += "/c"
				if unitsa[2]>1: unitss += "$^"+str(unitsa[2])+"$"
				if unitsa[3]>0: unitss += "/(mc$^2$)"
				if unitsa[3]>1: unitss += "$^"+str(unitsa[3])+"$"
			self.units[d] = " [ "+unitss+" ]"
		# finish title creation
		if len(self.diags) == 1:
			self.title = self.titles[self.diags[0]] + self.units[self.diags[0]]
		else:
			self.title = self.operation
			for d in self.diags:
				self.title = self.title.replace("#"+str(d), self.titles[d])
		if self.data_log: self.title = "Log[ "+self.title+" ]"
		
		# If any spatial dimension did not appear, then count it for calculating the correct density
		if ndim>=1 and not spatialaxes["x"]: units_coeff /= ncels[0]
		if ndim>=2 and not spatialaxes["y"]: units_coeff /= ncels[1]
		if ndim==3 and not spatialaxes["z"]: units_coeff /= ncels[2]
		units_coeff *= coeff_density
			
		# Compute the total coefficient of units
		units_operation = re.sub("#\d+",str(units_coeff),self.operation)
		exec("units_coeff = " + units_operation) in None
		
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
	
	
	# Gets info about diagnostic number "diagNumber"
	@staticmethod
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
		return {"#":diagNumber, "output":output, "every":every, "tavg":time_average, "species":species, "axes":axes}
	def getMyInfo(self, diagNumber):
		return self.getInfo(self.results_path, diagNumber)
	
	
	# Prints the info obtained by the function "getInfo"
	@staticmethod
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
	
	# Method to print info on all included diags
	def info(self):
		if not self.validate(): return
		for d in self.diags:
			self.printInfo(self.info_[d])
		if len(self.operation)>2: print "Operation : "+self.operation
		for ax in self.axes:
			if "sliceInfo" in ax: print ax["sliceInfo"]
		return
	
	# get all available timesteps for a given diagnostic
	def getAvailableTimesteps(self, diagNumber=None):
		# if argument "diagNumber" not provided, return the times calculated in __init__
		if diagNumber is None:
			return self.times
		# Otherwise, get the timesteps specifically available for the single requested diagnostic
		else:
			try:
				file = self.results_path+'/ParticleDiagnostic'+str(diagNumber)+'.h5'
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
		
	# Method to obtain the data only
	def getData(self, **kwargs):
		if not self.validate(): return
		# if optional argument "time" not provided, find out which time to plot
		try:    time = kwargs["time"]
		except: time = self.times[0]
		# Verify that the timestep is valid
		if time not in self.times:
			print "Timestep "+time+" not found in this diagnostic"
			return []
		# Get arrays from all requested diagnostics
		A = {}
		for d in self.diags:
			# Open file
			f = h5py.File(self.file[d], 'r')
			# get data
			index = self.data[d][time]
			A.update({ d:np.reshape(f.items()[index][1],self.shape) })
			f.close()
			# Apply the slicing
			for iaxis in range(self.naxes):
				axis = self.axes[iaxis]
				if "slice" in axis:
					A[d] = np.delete(A[d], axis["slice"], axis=iaxis) # remove parts outside of the slice
					A[d][np.isnan(A[d])] = 0.
					A[d] = np.sum(A[d], axis=iaxis, keepdims=True) # sum over the slice
			A[d] = np.squeeze(A[d]) # remove sliced axes
			# Divide by the bins size
			A[d] /= self.bsize
		# Calculate operation
		data_operation = self.operation
		for d in self.diags:
			data_operation = data_operation.replace("#"+str(d),"A["+str(d)+"]")
		exec("A = "+data_operation) in None
		# log scale if requested
		if self.data_log: A = np.log10(A)
		return A




# -------------------------------------------------------------------
# Class for fields diagnostics
# -------------------------------------------------------------------
class Field(Diagnostic):
	# We use __new__ to prevent object creation if no field requested
	def __new__(cls, *args, **kwargs):
		# Is there a "results_path" argument ?
		try   : results_path = cls.validatePath(*args, **kwargs)
		except: return None
		# Is there a "field" argument ?
		try:
			field = args[1]
		except:
			try:
				field = kwargs["field"]
			except:
				fields = cls.getFieldsIn(results_path)
				if len(fields)>0:
					print "Printing available fields:"
					print "--------------------------"
					l = (len(fields)/3) * 3
					if l>0:
						print '\n'.join(['\t\t'.join(list(i)) for i in np.reshape(fields[:l],(-1,3))])
					print '\t\t'.join(fields[l:])
				else:
					print "No fields found in '"+results_path+"'"
				return None
		# If everything ok, then we create the object
		return Diagnostic.__new__(cls, *args, **kwargs)
	
	
	# This is the constructor, which creates the object
	def __init__(self,results_path, field=None, timesteps=None, slice=None,
	             units="code", data_log=False, **kwargs):

		self.valid = False
		self.results_path = results_path

		# Get info from the input file and prepare units
		try:
			ndim               = self.read_ndim()
			sim_units          = self.read_sim_units()
			ncels, cell_length = self.read_ncels_cell_length(ndim, sim_units)
			timestep           = self.read_timestep(sim_units)
		except:
			return None
		
		if units == "nice":
			try   : wavelength_SI = self.read_wavelength_SI()
			except: return None
			cell_length *= 1e2*wavelength_SI/(2.*np.pi) # in cm
			cell_volume = np.prod(cell_length)
			coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 * cell_volume
			coeff_current = coeff_density * 4.803e-9
			self.coeff_time = timestep * wavelength_SI/3.e8 # in seconds
			self.time_units = " s"
		elif units == "code":
			coeff_density = 1.
			coeff_current = 1.
			self.coeff_time = timestep
			self.time_units = " 1/w"
		
		# Get available times
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
			print "No fields found in Fields.h5"
			return
		
		# Get available fields
		fields = self.getFields()
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check value of field
		if field not in fields:
			fs = filter(lambda x:field in x, fields)
			if len(fs)==0:		
				print "No field `"+field+"` found in Fields.h5"
				return
			if len(fs)>1:
				print "Several fields match: "+(' '.join(fs))
				print "Please be more specific and retry."
				return
			field = fs[0]
		self.fieldn = fields.index(field) # index of the requested field
		self.fieldname = field
		
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			print "Argument `slice` must be a dictionary"
			return None
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self.data_log = data_log
		
		# Default plot parameters
		self.setDefaultPlot()
		# Set user's plot parameters
		kwargs = self.setPlot(**kwargs)
		
		# Get the shape of field
		self.file = results_path+'/Fields.h5'
		f = h5py.File(self.file, 'r')
		sample = np.double(f.values()[0].values()[self.fieldn])
		shape = sample.shape
		f.close()
		
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with indices to the data arrays
		self.data = {}
		for i,t in enumerate(self.times):
			self.data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
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
					raise
			except:
				print "Argument `timesteps` must be one or two non-negative integers"
				return
		
		# Need at least one timestep
		if self.times.size < 1:
			print "Timesteps not found"
			return
		
			
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values
		self.naxes = ndim
		self.plot_type = []
		self.plot_label = []
		self.plot_centers = []
		self.plot_shape = []
		self.plot_log   = []
		self.sliceinfo = {}
		self.slices = [None]*ndim
		for iaxis in range(self.naxes):
			centers = np.linspace(0., shape[iaxis]*cell_length[iaxis], shape[iaxis])
			label = {0:"x", 1:"y", 2:"z"}[iaxis]
			axisunits = "[code units]"
			if units == "nice": axisunits = "[cm]"
			
			if label in slice:
				# if slice is "all", then all the axis has to be summed
				if slice[label] == "all":
					indices = np.arange(shape[iaxis])
				# Otherwise, get the slice from the argument `slice`
				else:
					try:
						s = np.double(slice[label])
						if s.size>2 or s.size<1: raise
					except:
						print "Slice along axis "+label+" should be one or two floats"
						return None
					if s.size==1:
						indices = np.array([(np.abs(centers-s)).argmin()])
					elif s.size==2:
						indices = np.nonzero( (centers>=s[0]) * (centers<=s[1]) )[0]
					if indices.size == 0:
						print "Slice along "+label+" is out of the box range"
						return None
					if indices.size == 1:
						self.sliceinfo.update({ label:"Sliced at "+label+" = "+str(centers[indices])+" "+axisunits })
					else:
						self.sliceinfo.update({ label:"Sliced for "+label
							+" from "+str(centers[indices[ 0]])+" to "+str(centers[indices[-1]])+" "+axisunits })
				# convert the range of indices into their "conjugate"
				self.slices[iaxis] = np.delete(np.arange(shape[iaxis]), indices)
			else:
				self.plot_type   .append(label)
				self.plot_shape  .append(shape[iaxis])
				self.plot_centers.append(centers)
				self.plot_label  .append(label+" "+axisunits)
				self.plot_log    .append(False)
			
		if len(self.plot_centers) > 2:
			print "Cannot plot in "+str(len(self.plot_shape))+"d. You need to 'slice' some axes."
			return
		
		# Build units
		self.fieldunits = "??"
		self.unitscoeff = "??"
		self.title = "??"
		if units == "nice":
			self.fieldunits = {"B":"T"  ,"E":"V/m"  ,"J":"A"          ,"R":"1/cm$^3$"   }[field[0]]
			self.unitscoeff = {"B":10710,"E":3.21e12,"J":coeff_current,"R":coeff_density}[field[0]]
			self.title      = field + "("+self.fieldunits+")"
		else:
			self.fieldunits = {"B":"$m_e\omega/e$","E":"$m_ec\omega/e$","J":"$ecn_c$","R":"$n_c$"}[field[0]]
			self.unitscoeff = {"B":1              ,"E":1               ,"J":1        ,"R":1      }[field[0]]
			self.title      = field + " in units of "+self.fieldunits
		if data_log: self.title = "Log[ "+self.title+" ]"
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included fields
	def info(self):
		if not self.validate(): return
		print "Field "+self.fieldname,
		#todo
		return
	
	# get all available fields
	@staticmethod
	def getFieldsIn(results_path):
		try:
			file = results_path+'/Fields.h5'
			f = h5py.File(file, 'r')
		except:
			print "Cannot open file "+file
			return []
		try:
			fields = f.values()[0].keys() # list of fields
		except:
			fields = []
		f.close()
		return fields
	def getFields(self):
		return self.getFieldsIn(self.results_path)
	
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:
			file = self.results_path+'/Fields.h5'
			f = h5py.File(file, 'r')
		except:
			print "Cannot open file "+file
			return np.array([])
		times = np.double(f.keys())
		f.close()
		return times
	
	# Method to obtain the data only
	def getData(self, **kwargs):
		if not self.validate(): return
		# if optional argument "time" not provided, find out which time to plot
		try:    time = kwargs["time"]
		except: time = self.times[0]
		# Verify that the timestep is valid
		if time not in self.times:
			print "Timestep "+time+" not found in this diagnostic"
			return []
		# Get arrays from requested field
		# Open file
		f = h5py.File(self.file, 'r')
		# get data
		index = self.data[time]
		A = np.double(f.values()[index].values()[self.fieldn])
		f.close()
		# Apply the slicing
		for iaxis in range(self.naxes):
			if self.slices[iaxis] is None: continue
			A = np.delete(A, self.slices[iaxis], axis=iaxis) # remove parts outside of the slice
			A = np.mean(A, axis=iaxis, keepdims=True) # sum over the slice
		A = np.squeeze(A) # remove sliced axes
		A *= self.unitscoeff
		# log scale if requested
		if self.data_log: A = np.log10(A)
		return A




# -------------------------------------------------------------------
# Class for scalars
# -------------------------------------------------------------------
class Scalar(Diagnostic):
	# We use __new__ to prevent object creation if no field requested
	def __new__(cls, *args, **kwargs):
		# Is there a "results_path" argument ?
		try   : results_path = cls.validatePath(*args, **kwargs)
		except: return None
		# Is there a "scalar" argument ?
		try:
			scalar = args[1]
		except:
			try:
				scalar = kwargs["scalar"]
			except:
				scalars = cls.getScalarsIn(results_path)
				if len(scalars)>0:
					print "Printing available scalars:"
					print "---------------------------"
					l = [""]
					for s in scalars:
						if s[:2] != l[-1][:2] and s[-2:]!=l[-1][-2:]:
							if l!=[""]: print "\t".join(l)
							l = []
						l.append(s)
				else:
					print "No scalars found in '"+results_path+"'"
				return None
		# If everything ok, then we create the object
		return Diagnostic.__new__(cls, *args, **kwargs)
	
	
	# This is the constructor, which creates the object
	def __init__(self,results_path, scalar=None, timesteps=None,
	             units="code", data_log=False, **kwargs):
	
		self.valid = False
		self.results_path = results_path

		# Get info from the input file and prepare units
		try:
			sim_units          = self.read_sim_units()
			self.timestep      = self.read_timestep("")
		except:
			return None
		
		if units == "nice":
			try   : wavelength_SI = self.read_wavelength_SI()
			except: return None
			self.coeff_time = self.timestep * wavelength_SI/3.e8/(2.*np.pi) # in seconds
			self.time_units = " s"
		elif units == "code":
			self.coeff_time = self.timestep
			self.time_units = " 1/w"
		if sim_units == "wavelength": self.coeff_time *= 2. * np.pi
		
		# Get available scalars
		scalars = self.getScalars()
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check value of field
		if scalar not in scalars:
			fs = filter(lambda x:scalar in x, scalars)
			if len(fs)==0:		
				print "No scalar `"+scalar+"` found in scalars.txt"
				return
			if len(fs)>1:
				print "Several scalars match: "+(' '.join(fs))
				print "Please be more specific and retry."
				return
			scalar = fs[0]
		self.scalarn = scalars.index(scalar) # index of the requested scalar
		self.scalarname = scalar
		
		# Put data_log as object's variable
		self.data_log = data_log
		
		# Default plot parameters
		self.setDefaultPlot()
		# Set user's plot parameters
		kwargs = self.setPlot(**kwargs)
		
		# Already get the data from the file
		# Loop file line by line
		self.times = []
		self.values = []
		file = results_path+'/scalars.txt'
		f = open(file, 'r')
		for line in f:
			line = line.strip()
			if line[0]=="#": continue
			line = line.split()
			self.times .append( int( float(line[0]) / self.timestep ) )
			self.values.append( float(line[self.scalarn+1]) )
		self.times  = np.array(self.times )
		self.values = np.array(self.values)
		f.close()
		
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with the index to each time
		self.data = {}
		for i,t in enumerate(self.times):
			self.data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
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
					raise
			except:
				print "Argument `timesteps` must be one or two non-negative integers"
				return
		
		# Need at least one timestep
		if self.times.size < 1:
			print "Timesteps not found"
			return
		
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# There are no axes for scalars
		self.naxes = 0
		self.plot_type = []
		self.plot_label = []
		self.plot_centers = []
		self.plot_shape = []
		self.plot_log   = []
		self.slices = []
		# Build units
		self.scalarunits = "unknown units"
		self.unitscoeff = 1.
		self.title = "??"
		if units == "nice":
			self.title      = scalar +"( "+self.scalarunits+" )" # todo
		else:
			self.title      = scalar +"( "+self.scalarunits+" )" # todo
		if data_log: self.title = "Log[ "+self.title+" ]"
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included scalars
	def info(self):
		if not self.validate(): return
		print "Scalar "+self.scalarname,
		#todo
		return
	
	# get all available scalars
	@staticmethod
	def getScalarsIn(results_path):
		try:
			file = results_path+'/scalars.txt'
			f = open(file, 'r')
		except:
			print "Cannot open file "+file
			return []
		try:
			# Find last commented line 
			prevline = ""
			for line in f:
				line = line.strip()
				if line[0]!="#": break
				prevline = line[1:].strip()
			scalars = prevline.split() # list of scalars
			scalars = scalars[1:] # remove first, which is "time"
		except:
			scalars = []
		f.close()
		return scalars
	def getScalars(self):
		return self.getScalarsIn(self.results_path)
	
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		return self.times
	
	# Method to obtain the data only
	def getData(self, **kwargs):
		if not self.validate(): return
		# if optional argument "time" not provided, find out which time to plot
		try:    time = kwargs["time"]
		except: time = self.times[0]
		# Verify that the timestep is valid
		if time not in self.times:
			print "Timestep "+time+" not found in this diagnostic"
			return []
		# Get value at selected time
		A = self.values[ self.data[time] ]
		A *= self.unitscoeff
		# log scale if requested
		if self.data_log: A = np.log10(A)
		return A




# Function to plot multiple particle diags on the same figure
def multiPlot(*Pdiags, **kwargs):
	nPdiags = len(Pdiags)
	# Verify Pdiags are valid
	if nPdiags == 0: return
	for Pdiag in Pdiags:
		if not Pdiag.validate(): return
	# Gather all times
	alltimes = np.unique(np.concatenate([Pdiag.times for Pdiag in Pdiags]))
	# Get keyword arguments
	figure = kwargs.pop("figure", 1)
	shape  = kwargs.pop("shape" , None)
	# Determine whether to plot all cases on the same axes
	sameAxes = False
	if shape is None or shape == [1,1]:
		sameAxes = True
		for Pdiag in Pdiags:
			if len(Pdiag.plot_shape)==0 and len(Pdiags[0].plot_shape)==0:
				continue
			if len(Pdiag.plot_shape)!=1 or Pdiag.plot_type!=Pdiags[0].plot_type:
				sameAxes = False
				break
	if not sameAxes and shape == [1,1]:
		print "Cannot have shape=[1,1] with these diagnostics"
		return
	# Determine the shape
	if sameAxes: shape = [1,1]
	if shape is None: shape = [nPdiags,1]
	nplots = np.array(shape).prod()
	if not sameAxes and nplots != nPdiags:
		print "The 'shape' argument is incompatible with the number of diagnostics:"
		print "  "+str(nPdiags)+" diagnostics do not fit "+str(nplots)+" plots"
		return
	# Make the figure
	fig = plt.figure(figure)
	fig.set(**matplotlibArgs(kwargs)[0]) # Apply figure kwargs
	fig.clf()
	fig.subplots_adjust(wspace=0.5, hspace=0.5, bottom=0.15)
	ax = []
	xmin =  float("inf")
	xmax = -float("inf")
	c = plt.matplotlib.rcParams['axes.color_cycle']
	for i in range(nplots):
		ax.append( fig.add_subplot(shape[0], shape[1], i) )
	for i, Pdiag in enumerate(Pdiags):
		if sameAxes: Pdiag.ax = ax[0]
		else       : Pdiag.ax = ax[i]
		Pdiag.artist = None
		try:
			l = Pdiag.limits()[0]
			xmin = min(xmin,l[0])
			xmax = max(xmax,l[1])
		except:
			pass
		if "color" not in Pdiag.plotkwargs:
			Pdiag.plotkwargs.update({ "color":c[i%len(c)] })
	# Static plot
	if sameAxes and len(Pdiags[0].plot_shape)==0:
		for Pdiag in Pdiags:
			Pdiag.artist = Pdiag.plotVsTime(Pdiag.ax)
		fig.canvas.draw()
		plt.show()
	# Animated plot
	else:
		# Loop all times
		for time in alltimes:
			for Pdiag in Pdiags:
				if time in Pdiag.times:
					if sameAxes:
						if Pdiag.artist is not None: Pdiag.artist.remove()
					else:
						Pdiag.ax.cla()
					Pdiag.artist = Pdiag.animateOnAxes(Pdiag.ax, time)
					if sameAxes:
						Pdiag.ax.set_xlim(xmin,xmax)
			fig.canvas.draw()
			plt.show()
		return
	

