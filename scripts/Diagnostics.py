# -----------------------------------------------------------------------
# HOW TO VIEW DIAGNOSTICS                -    F. Perez - 03/2015
# -----------------------------------------------------------------------
# >>>>>> Requirements
#   python2.7 with the following packages: numpy, matplotlib, pylab, h5py
#
# >>>>>> First step: invoke python and load this file
#      $ python -i Diagnostics.py
#
#
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
# >>>>>> In the python shell, use the class "ParticleDiagnostic"
#        to create a ParticleDiagnostic object
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
# >>>>>> To plot the data, use the following method.
#
# ParticleDiagnostic(..., figure=1, vmin=None, vmax=None,
#                         xmin=None, xmax=None, ymin=None, ymax=None) .plot()
#
#            figure = _int_       (optional)
#                     The figure number that is passed to matplotlib.
#                     If absent, figure 1 is used.
#
#              vmin = _double_    (optional)
#              vmax = _double_    (optional)
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
#   This method returns only a list of the data arrays (for each timestep requested).
#
# ParticleDiagnostic(...).get()
#   This method returns a dictionary containing the data, the list of timesteps and the axes scales.
#
#
# >>>>>> Examples:
#    ParticleDiagnostic('../test', diagNumber=1, slice={"y":"all"}, units="nice").plot(figure=1, vmin=0, vmax=3e14)
#
#
#  +---------------------------+
#  | 2 . Fields                |
#  +---------------------------+
#
# You can view fields with a similar procedure:
#
# Field(results_path, field=None, timesteps=None, slice=None, units="code", data_log=False)
#
# This works almost the same way as ParticleDiagnostic(...), with some exceptions:
#
# - `field` must be one of "Bx_m", "By_m", "Bz_m", "Ex", "Ey", "Ez", "Jx", "Jy", "Jz",
#   "Jx_[species]", "Jy_[species]", "Jz_[species]", "Rho" or "Rho_[species]"
#   where [species] is the name of one of the existing species.
#   If you omit the argument `field`, the list of available fields will be displayed.
#   Additionally, you can write an operation instead of just one field. For instance,
#   you can have "Jx+Jy".
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
#
# For scalars, this is almost the same:
#
# Scalar(results_path, scalar=None, timesteps=None, units="code", data_log=False)
#
# - `scalar` must be an available scalar name. To get a list of available scalars, 
#    simply omit this argument.
#
# - there is no more `slice` argument of course!
#
#
#  +---------------------------+
#  | 4 . Multiple diagnostics  |
#  +---------------------------+
# To simultaneously plot multiple diagnostics in the same figure:
#
# multiPlot(diag1, diag2, ... , figure=1, shape=None, **kwargs)
#
#              diag1 = diagnostic prepared by ParticleDiagnostic(), Field() or Scalar()
#              diag2 = diagnostic prepared by ParticleDiagnostic(), Field() or Scalar()
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
#            kwargs = many other keyword-arguments can be used -> refer to the doc.


class Smilei(object):
	""" Smilei(results_path=".")
	
	Import Smilei simulation information.
	
	The argument `results_path` specifies where the simulation results are located.
	Omit this argument if you are already in the results path.
	"""
	
	valid = False
	
	def __init__(self, results_path="."):
		# Import packages
		import h5py
		import numpy as np
		import os.path, glob, re, sys
		import matplotlib.pyplot
		import pylab
		pylab.ion()
		# Transfer packages to local attributes
		self._results_path = results_path
		self._h5py = h5py
		self._np = np
		self._ospath = os.path
		self._glob = glob.glob
		self._re = re
		self._plt = matplotlib.pyplot
		# Verify that results_path is valid
		if not self._ospath.isdir(results_path):
			print "Could not find directory "+results_path
			return
		if len(self._glob(results_path+"/smilei.py"))==0:
			print "Could not find an input file in directory "+results_path
			return
		# Fetch the python namelist
		namespace={}
		execfile(results_path+'/smilei.py',namespace) # execute the namelist into an empty namespace
		class Namelist: pass # empty class to store the namelist variables
		self.namelist = Namelist() # create new empty object
		for key, value in namespace.iteritems(): # transfer all variables to this object
			if key[0]=="_": continue # skip builtins
			setattr(self.namelist, key, value)
		
		self.valid = True
	
	def __repr__(self):
		if not self.valid:
			return "Invalid Smilei simulation"
		file = self._glob(self._results_path+"/smilei.py")[0]
		return "Smilei simulation with input file located at `"+file+"`"
	
	# Methods to create instances of diagnostics
	def Scalar(self, *args, **kwargs):
		""" Scalar(scalar=None, timesteps=None, units="code", data_log=False)
		
		Import and analyze a scalar diagnostic from a Smilei simulation
		
		Parameters:
		-----------
		scalar : an available scalar name. (optional)
			To get a list of available scalars, simply omit this argument.
		timesteps : int or [int, int] (optional)
			If omitted, all timesteps are used.
			If one number  given, the nearest timestep available is used.
			If two numbers given, all the timesteps in between are used.
		units : "code" or "nice"    (optional)
			If "nice" is chosen, then units are converted into usual units:
			distances in microns, density in 1/cm^3, energy in MeV.
		data_log : True or False    (optional)
			If True, then log10 is applied to the output array before plotting.
		
		Usage:
		------
			scalar = S.Scalar(...)
		where S is a `Smilei` object.
		"""
		return Scalar(self, *args, **kwargs)
	def Field(self, *args, **kwargs):
		""" Field(field=None, timesteps=None, slice=None, units="code", data_log=False)
		
		Import and analyze a field diagnostic from a Smilei simulation
		
		Parameters:
		-----------
		field : an available field name. (optional)
			To get a list of available fields, simply omit this argument.
			You may write an operation instead of just one field, e.g. "Jx+Jy".
		timesteps : int or [int, int] (optional)
			If omitted, all timesteps are used.
			If one number  given, the nearest timestep available is used.
			If two numbers given, all the timesteps in between are used.
		slice : a python dictionary of the form { axis:range, ... } (optional)
			`axis` may be "x", "y" or "z".
			`range` may be "all", a float, or [float, float].
			For instance, slice={"x":"all", "y":[2,3]}.
			The average of all values within the 'slice' is computed.
		units : "code" or "nice"    (optional)
			If "nice" is chosen, then units are converted into usual units:
			distances in microns, density in 1/cm^3, energy in MeV.
		data_log : True or False    (optional)
			If True, then log10 is applied to the output array before plotting.
		
		Usage:
		------
			field = S.Field(...) # S is a Smilei object
			field.get()
			field.plot()
		"""
		return Field(self, *args, **kwargs)
	def Probe(self, *args, **kwargs):
		""" Probe(probeNumber=None, field=None, timesteps=None, slice=None, units="code", data_log=False)
		
		Import and analyze a probe diagnostic from a Smilei simulation
		
		Parameters:
		-----------
		probeNumber : index of an available probe. (optional)
			To get a list of available probes, simply omit this argument.
		field : name of an available field. (optional)
			To get a list of available fields, simply omit this argument.
		timesteps : int or [int, int] (optional)
			If omitted, all timesteps are used.
			If one number  given, the nearest timestep available is used.
			If two numbers given, all the timesteps in between are used.
		slice : a python dictionary of the form { axis:range, ... } (optional)
			`axis` may be "axis1" or "axis2" (the probe axes).
			`range` may be "all", a float, or [float, float].
			For instance, slice={"axis1":"all", "axis2":[2,3]}.
			The average of all values within the 'slice' is computed.
		units : "code" or "nice"    (optional)
			If "nice" is chosen, then units are converted into usual units:
			distances in microns, density in 1/cm^3, energy in MeV.
		data_log : True or False    (optional)
			If True, then log10 is applied to the output array before plotting.
		
		Usage:
		------
			probe = S.Probe(...) # S is a Smilei object
			probe.get()
			probe.plot()
		"""
		return Probe(self, *args, **kwargs)
	def ParticleDiagnostic(self, *args, **kwargs):
		""" ParticleDiagnostic(diagNumber=None, timesteps=None, slice=None, units="code", data_log=False)
		
		Import and analyze a particle diagnostic from a Smilei simulation
		
		Parameters:
		-----------
		diagNumber : index of an available particle diagnostic. (optional)
			To get a list of available diags, simply omit this argument.
		timesteps : int or [int, int] (optional)
			If omitted, all timesteps are used.
			If one number  given, the nearest timestep available is used.
			If two numbers given, all the timesteps in between are used.
		slice : a python dictionary of the form { axis:range, ... } (optional)
			`axis` may be "x", "y", "z", "px", "py", "pz", "p", "gamma", "ekin", "vx", "vy", "vz", "v" or "charge".
			`range` may be "all", a float, or [float, float].
			For instance, slice={"x":"all", "y":[2,3]}.
			The SUM of all values within the 'slice' is computed.
		units : "code" or "nice"    (optional)
			If "nice" is chosen, then units are converted into usual units:
			distances in microns, density in 1/cm^3, energy in MeV.
		data_log : True or False    (optional)
			If True, then log10 is applied to the output array before plotting.
		
		Usage:
		------
			diag = S.ParticleDiagnostic(...) # S is a Smilei object
			diag.get()
			diag.plot()
		"""
		return ParticleDiagnostic(self, *args, **kwargs)
	
	
class Options(object):
	""" Class to contain matplotlib plotting options """
	
	def __init__(self, **kwargs):
		self.figure = 1
		self.xmin   = None
		self.xmax   = None
		self.ymin   = None
		self.ymax   = None
		self.vmin   = None
		self.vmax   = None
		self.figure0 = {}
		self.figure1 = {"facecolor":"w"}
		self.axes = {}
		self.plot = {}
		self.image = {"interpolation":"nearest", "aspect":"auto"}
		self.colorbar = {}
		self.xtick = {"useOffset":False}
		self.ytick = {"useOffset":False}
		self.set(**kwargs)
	
	# Method to set optional plotting arguments
	def set(self, **kwargs):
		# First, we manage the main optional arguments
		self.figure0.update({ "num":kwargs.pop("figure",self.figure) })
		self.xmin     = kwargs.pop("xmin"   ,self.xmin  )
		self.xmax     = kwargs.pop("xmax"   ,self.xmax  )
		self.ymin     = kwargs.pop("ymin"   ,self.ymin  )
		self.ymax     = kwargs.pop("ymax"   ,self.ymax  )
		self.vmin     = kwargs.pop("vmin"   ,kwargs.pop("data_min",self.vmin))
		self.vmax     = kwargs.pop("vmax"   ,kwargs.pop("data_max",self.vmax))
		# Second, we manage all the other arguments that are directly the ones of matplotlib
		for kwa in kwargs:
			val = kwargs[kwa]
			if kwa in ["figsize"]:
				self.figure0.update({kwa:val})
			if kwa in ["dpi","facecolor","edgecolor"]:
				self.figure1.update({kwa:val})
			if kwa in ["aspect","axis_bgcolor",
					   "frame_on","position","title","visible","xlabel","xscale","xticklabels",
					   "xticks","ylabel","yscale","yticklabels","yticks","zorder"]:
				self.axes.update({kwa:val})
			if kwa in ["color","dashes","drawstyle","fillstyle","label","linestyle",
					   "linewidth","marker","markeredgecolor","markeredgewidth",
					   "markerfacecolor","markerfacecoloralt","markersize","markevery",
					   "visible","zorder"]:
				self.plot.update({kwa:val})
			if kwa in ["cmap","aspect","interpolation"]:
				self.image.update({kwa:val})
			if kwa in ["orientation","fraction","pad","shrink","anchor","panchor",
					   "extend","extendfrac","extendrect","spacing","ticks","format",
					   "drawedges"]:
				self.colorbar.update({kwa:val})
			if kwa in ["style_x","scilimits_x","useOffset_x"]:
				self.xtick.update({kwa[:-2]:val})
			if kwa in ["style_y","scilimits_y","useOffset_y"]:
				self.ytick.update({kwa[:-2]:val})
		# special case: "aspect" is ambiguous because it exists for both imshow and colorbar
		if "cbaspect" in kwargs:
			self.colorbar.update({"aspect":kwargs["cbaspect"]})


# -------------------------------------------------------------------
# Mother class for all diagnostics
# -------------------------------------------------------------------
class Diagnostic(object):
	
	# Initialize with "results_path" argument being either the `results_path` or
	# the parent `Smilei` object
	def __init__(self, results_path=None, *args, **kwargs):
		self.valid = False
		self._previousdata = None
		# if string, try to use it as a results_path
		if type(results_path) is str:
			self.Smilei = Smilei(results_path)
		# if given a Smilei object, use it directly
		elif type(results_path) is Smilei:
			self.Smilei = results_path
		# Otherwise, error
		else:
			print "Could not find information on the Smilei simulation"
			return
		if not self.Smilei.valid: return
		# pass packages to each diagnostic
		self._results_path = self.Smilei._results_path
		self._h5py = self.Smilei._h5py
		self._np = self.Smilei._np
		self._ospath = self.Smilei._ospath
		self._glob = self.Smilei._glob
		self._re = self.Smilei._re
		self._plt = self.Smilei._plt
		self.namelist = self.Smilei.namelist
		# Make the Options object
		self.options = Options(**kwargs)
		# Call the '_init' function of the child class
		self._init(*args, **kwargs)
	
	# When no action is performed on the object, this is what appears
	def __repr__(self):
		if not self.valid: return ""
		self.info()
		return ""
	
	# Various methods to extract stuff from the input file
	def _read_ndim(self):
		try:
			dim = self.namelist.dim
			ndim = int(dim[0])
		except:
			print "Could not extract 'dim' from the input file"
			raise
		if ndim not in [1,2,3]:
			print "Could not understand simulation dimension 'dim="+dim+"' from the input file"
			raise
		return ndim
	def _read_ncels_cell_length(self, ndim):
		try:
			sim_length = self._np.double( self.namelist.sim_length )
			if sim_length.size==0: raise
		except:
			print "Could not extract 'sim_length' from the input file"
			raise
		try:
			cell_length = self._np.double( self.namelist.cell_length )
			if cell_length.size==0: raise
		except:
			try:
				res_space = self._np.double( self.namelist.res_space )
				cell_length = 1./res_space
				if cell_length.size==0: raise
			except:
				print "Could not extract 'cell_length' or 'res_space' from the input file"
				raise
		if   ndim == 1:
			sim_length  = sim_length[0]
			cell_length = cell_length[0]
		elif ndim == 2:
			if sim_length.size  == 1: sim_length  = self._np.array([sim_length,sim_length])
			else                    : sim_length  = sim_length[0:2]
			if cell_length.size == 1: cell_length = self._np.array([cell_length,cell_length])
			else                    : cell_length = cell_length[0:2]
		elif ndim == 3:
			if sim_length.size == 1: sim_length = self._np.array([sim_length,sim_length,sim_length])
			elif sim_length.size >2: sim_length = sim_length[0:3]
			else:
				print "In the input file, 'sim_length' should have 1 or 3 arguments for a 3d simulation"
				raise
			if cell_length.size == 1: cell_length = self._np.array([cell_length,cell_length,cell_length])
			elif cell_length.size >2: cell_length = cell_length[0:3]
			else:
				print "In the input file, 'cell_length' or 'res_space' should have 1 or 3 arguments for a 3d simulation"
				raise
		sim_length  = self._np.array(sim_length ,ndmin=1)
		cell_length = self._np.array(cell_length,ndmin=1)
		ncels = sim_length/cell_length
		return ncels, cell_length
	def _read_timestep(self):
		try:
			timestep = self._np.double(self.namelist.timestep)
			if not self._np.isfinite(timestep): raise 
		except:
			try:
				res_time = self._np.double(self.namelist.res_time)
				timestep = 1./res_time
				if not self._np.isfinite(timestep): raise 
			except:
				print "Could not extract 'timestep' or 'res_time' from the input file"
				raise
		return timestep
	def _read_wavelength_SI(self):
		try:
			wavelength_SI = self._np.double( self.namelist.wavelength_SI )
		except:
			print "Could not extract 'wavelength_SI' from the input file"
			raise
		return wavelength_SI
	
	# Method to verify everything was ok during initialization
	def _validate(self):
		if not self.Smilei.valid or not self.valid:
			print "Diagnostic is invalid"
			return False
		return True
		
	# Method to set optional plotting arguments
	def set(self, **kwargs):
		self.options.set(**kwargs)
	
	# Method to obtain the plot limits
	def limits(self):
		l = []
		for i in range(len(self._plot_shape)):
			l.append([min(self._plot_centers[0]), max(self._plot_centers[0])])
		return l
	
	# Method to get only the arrays of data
	def getData(self):
		if not self._validate(): return
		data = []
		for t in self.times:
			data.append( self._getDataAtTime(t) )
		return data
	
	# Method to obtain the data and the axes
	def get(self):
		if not self._validate(): return
		# obtain the data arrays
		data = self.getData()
		# format the results into a dictionary
		result = {"data":data, "times":self.times}
		for i in range(len(self._plot_type)):
			result.update({ self._plot_type[i]:self._plot_centers[i] })
		return result
	
	# Method to plot the current diagnostic
	def plot(self, **kwargs):
		if not self._validate(): return
		self.set(**kwargs)
		self.info()
		
		# Make figure
		fig = self._plt.figure(**self.options.figure0)
		fig.set(**self.options.figure1)
		fig.clf()
		ax = fig.add_subplot(1,1,1)
		# Animation if several dimensions
		if len(self._plot_shape) > 0:
			# Loop times
			for timeindex in range(self.times.size):
				time = self.times[timeindex]
				print "timestep "+str(time)+ "   -   t = "+str(time*self._coeff_time)+self._time_units
				# plot
				ax.cla()
				artist = self._animateOnAxes(ax, time)
				fig.canvas.draw()
				self._plt.show()
			#return artist
		# Static plot if 0 dimensions
		else:
			ax.cla()
			artist = self._plotVsTime(ax)
	
	# Method to set limits to a plot
	def _setLimits(self, ax, xmin=None, xmax=None, ymin=None, ymax=None):
		if xmin is not None: ax.set_xlim(xmin=xmin)
		if xmax is not None: ax.set_xlim(xmax=xmax)
		if ymin is not None: ax.set_ylim(ymin=ymin)
		if ymax is not None: ax.set_ylim(ymax=ymax)
	
	# Method to plot the data when axes are made
	def _animateOnAxes(self, ax, t):
		if not self._validate(): return None
		# get data
		A = self._getDataAtTime(t)
		# plot
		if A.ndim == 0: # as a function of time
			if self._previousdata is None:
				self._previousdata = self._np.zeros(self.times.size)
				for i, t in enumerate(self.times):
					self._previousdata[i] = self._getDataAtTime(t)
			times = self.times[self.times<=t]
			A     = self._previousdata[self.times<=t]
			im, = ax.plot(times*self._coeff_time, A, **self.options.plot)
			ax.set_xlabel('Time ['+self._time_units+' ]')
			self._setLimits(ax, xmax=self.times[-1]*self._coeff_time, ymin=self.options.vmin, ymax=self.options.vmax)
		elif A.ndim == 1:
			im, = ax.plot(self._plot_centers[0], A, **self.options.plot)
			if self._plot_log[0]: ax.set_xscale("log")
			ax.set_xlabel(self._plot_label[0])
			self._setLimits(ax, ymin=self.options.vmin, ymax=self.options.vmax)
		elif A.ndim == 2:
			im = self._animateOnAxes_2D(ax, A)
			if (self._plot_log[0]): ax.set_xlabel("Log[ "+self._plot_label[0]+" ]")
			else:                   ax.set_xlabel(        self._plot_label[0]     )
			if (self._plot_log[1]): ax.set_ylabel("Log[ "+self._plot_label[1]+" ]")
			else:                   ax.set_ylabel(        self._plot_label[1]     )
			self._setLimits(ax, ymin=self.options.ymin, ymax=self.options.ymax)
			try: # if colorbar exists
				ax.cax.cla()
				self._plt.colorbar(mappable=im, cax=ax.cax, **self.options.colorbar)
			except AttributeError:
				ax.cax = self._plt.colorbar(mappable=im, ax=ax, **self.options.colorbar).ax
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax)
		if self._title is not None: ax.set_title(self._title)
		ax.set(**self.options.axes)
		try:
			if len(self.options.xtick)>0: ax.ticklabel_format(axis="x",**self.options.xtick)
		except:
			print "Cannot format x ticks (typically happens with log-scale)"
			self.xtickkwargs = []
		try:
			if len(self.options.ytick)>0: ax.ticklabel_format(axis="y",**self.options.ytick)
		except:
			print "Cannot format y ticks (typically happens with log-scale)"
			self.xtickkwargs = []
		return im
	
	# Special case: 2D plot
	# This is overloaded by class "Probe" because it requires to replace imshow
	def _animateOnAxes_2D(self, ax, A):
		extent = [self._plot_centers[0][0], self._plot_centers[0][-1], self._plot_centers[1][0], self._plot_centers[1][-1]]
		if self._plot_log[0]: extent[0:2] = [self._np.log10(self._plot_centers[0][0]), self._np.log10(self._plot_centers[0][-1])]
		if self._plot_log[1]: extent[2:4] = [self._np.log10(self._plot_centers[1][0]), self._np.log10(self._plot_centers[1][-1])]
		im = ax.imshow( self._np.flipud(A.transpose()),
			vmin = self.options.vmin, vmax = self.options.vmax, extent=extent, **self.options.image)
		return im
	
	# If the sliced data has 0 dimension, this function can plot it 
	def _plotVsTime(self, ax):
		if len(self._plot_shape) > 0:
			print "To plot vs. time, it is necessary to slice all axes in order to obtain a 0-D array"
			return None
		# Loop times to gather data
		A = self._np.squeeze(self.getData())
		im, = ax.plot(self.times*self._coeff_time, A, **self.options.plot)
		ax.set_xlabel('Time ['+self._time_units+' ]')
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.vmin, ymax=self.options.vmax)
		if self._title is not None: ax.set_title(self._title)
		ax.set(**self.options.axes)
		ax.ticklabel_format(axis="x",**self.options.xtick)
		ax.ticklabel_format(axis="y",**self.options.ytick)
		return im
	
	
	def dim(self):
		return len(self._plot_shape)
	


# -------------------------------------------------------------------
# Class for particle diagnostics
# -------------------------------------------------------------------
class ParticleDiagnostic(Diagnostic):

	# This is the constructor, which creates the object
	def _init(self, diagNumber=None, timesteps=None, slice=None,
				 units="code", data_log=False, **kwargs):
		
		if not self.Smilei.valid: return None
		if diagNumber is None:
			print "Printing available particle diagnostics:"
			print "----------------------------------------"
			diagNumber = 0
			while self._printInfo(self._getInfo(diagNumber)):
				diagNumber += 1
			if diagNumber == 0:
				print "      No particle diagnostics found in "+self._results_path;
			return None
		
		
		# Get info from the input file and prepare units
		try:
			ndim               = self._read_ndim()
			
			ncels, cell_length = self._read_ncels_cell_length(ndim)
			self.timestep      = self._read_timestep()
			cell_size = {"x":cell_length[0]}
			if ndim>1: cell_size.update({"y":cell_length[1]})
			if ndim>2: cell_size.update({"z":cell_length[2]})
		except:
			return None
		
		if units == "nice":
			try   : wavelength_SI = self._read_wavelength_SI()
			except: return None
			coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 # nc in cm^-3
			coeff_energy = 0.511
			self._coeff_time = self.timestep * wavelength_SI/3.e8 # in seconds
			self._time_units = " s"
		elif units == "code":
			coeff_density = 1.
			coeff_energy = 1.
			self._coeff_time = self.timestep
			self._time_units = " $1/\omega$"
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
			print"Argument 'diagNumber' must be and integer or a string."
			return None
		
		# Get list of requested diags
		self._diags = sorted(set([ int(d[1:]) for d in self._re.findall('#\d+',self.operation) ]))
		for diag in self._diags:
			if not self._getInfo(diag):
				print "No particle diagnostic #"+str(diag)
				return None
		try:
			exec(self._re.sub('#\d+','1.',self.operation)) in None
		except:
			print "Cannot understand operation '"+self.operation+"'"
			return None
		# Verify that all requested diags exist and they all have the same shape
		self._info = {}
		self._shape = {}
		self._axes = {}
		self._naxes = {}
		for d in self._diags:
			try:
				self._info.update({ d:self._getInfo(d) })
			except:
				print "Particle diagnostic #"+str(d)+" not found."
				return None
			self._axes .update ({ d:self._info[d]["axes"] })
			self._naxes.update ({ d:len(self._axes[d]) })
			self._shape.update({ d:[ axis["size"] for axis in self._axes[d] ] })
			if self._naxes[d] != self._naxes[self._diags[0]]:
				print ("All diagnostics in operation '"+self.operation+"' must have as many axes."
					+ " Diagnotic #"+str(d)+" has "+str(self._naxes[d])+" axes and #"+
					str(self._diags[0])+" has "+str(self._naxes[self._diags[0]])+" axes")
				return None
			for a in self._axes[d]:
				if self._axes[d] != self._axes[self._diags[0]]:
					print ("In operation '"+self.operation+"', diagnostics #"+str(d)+" and #"
						+str(self._diags[0])+" must have the same shape.")
					return None
		
		self._axes  = self._axes [self._diags[0]]
		self._naxes = self._naxes[self._diags[0]]
		self._shape = self._shape[self._diags[0]]
		
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			print "Argument 'slice' must be a dictionary"
			return None
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# Get available timesteps
		self._file = {}
		self.times = {}
		self._data  = {}
		self._h5items = {}
		for d in self._diags:
			# get all diagnostics files
			self._file.update({ d:self._h5py.File(self._results_path+'/ParticleDiagnostic'+str(d)+'.h5') })
			self._h5items.update({ d:self._file[d].items() })
			# get all diagnostics timesteps
			self.times.update({ d:self.getAvailableTimesteps(d) })
			# fill the "data" dictionary with indices to the data arrays
			self._data.update({ d:{} })
			for i,t in enumerate(self.times[d]):
				self._data[d].update({ t : i })
			# If timesteps is None, then keep all timesteps, otherwise, select timesteps
			if timesteps is not None:
				try:
					ts = self._np.array(self._np.double(timesteps),ndmin=1)
					if ts.size==2:
						# get all times in between bounds
						self.times[d] = self.times[d][ (self.times[d]>=ts[0]) * (self.times[d]<=ts[1]) ]
					elif ts.size==1:
						# get nearest time
						self.times[d] = self._np.array([self.times[d][(self._np.abs(self.times[d]-ts)).argmin()]])
					else:
						raise
				except:
					print "Argument 'timesteps' must be one or two non-negative integers"
					return None
			# Verify that timesteps are the same for all diagnostics
			if (self.times[d] != self.times[self._diags[0]]).any() :
				print ("All diagnostics in operation '"+self.operation+"' must have the same timesteps."
					+" Diagnotic #"+str(d)+" has "+str(len(self.times[d]))+ " timesteps and #"
					+str(self._diags[0])+" has "+str(len(self.times[self._diags[0]])))+ " timesteps"
				return None
		# Now we need to keep only one array of timesteps because they should be all the same
		self.times = self.times[self._diags[0]]
		
		# Need at least one timestep
		if self.times.size < 1:
			print "Timesteps not found"
			return None
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values for all diags
		self._plot_shape = []; self._plot_type = []; plot_diff = []
		self._plot_label = []; self._plot_centers = []; self._plot_log = []
		units_coeff = 1.
		unitsa = [0,0,0,0]
		spatialaxes = {"x":False, "y":False, "z":False}
		for axis in self._axes:
		
			# Find the vector of values along the axis
			if axis["log"]:
				edges = self._np.linspace(self._np.log10(axis["min"]), self._np.log10(axis["max"]), axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				edges = 10.**edges
				centers = 10.**(centers[:-1])
			else:
				edges = self._np.linspace(axis["min"], axis["max"], axis["size"]+1)
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
					axis_coeff = 1e6*wavelength_SI/(2.*self._np.pi)
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
					indices = self._np.arange(axis["size"])
				
				# Otherwise, get the slice from the argument `slice`
				else:
					try:
						s = self._np.double(slice[axis["type"]])
						if s.size>2 or s.size<1: raise
					except:
						print "Slice along axis "+axis["type"]+" should be one or two floats"
						return None
					# convert the slice into a range of indices
					if s.size == 1:
						indices = self._np.array([(self._np.abs(centers-s)).argmin()])
					else :
						indices = self._np.nonzero( (centers>=s[0]) * (centers<=s[1]) )[0]
				
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
				indices = self._np.delete(self._np.arange(axis["size"]), indices)
				# put the slice in the dictionary
				axis.update({"slice":indices, "slice_size":slice_size})
			
				if axis["type"] in ["x","y","z"]:
					units_coeff *= cell_size[axis["type"]]/slice_size
			
			# if not sliced, then add this axis to the overall plot
			else:
				self._plot_type   .append(axis["type"])
				self._plot_shape  .append(axis["size"])
				self._plot_centers.append(centers*axis_coeff)
				self._plot_log    .append(axis["log"])
				self._plot_label  .append(axis["type"]+axis_units)
				plot_diff.append(self._np.diff(edges))
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
		
		
		if len(self._plot_shape) > 2:
			print "Cannot plot in "+str(len(self._plot_shape))+"d. You need to 'slice' some axes."
			return None
		
		# Build units
		self._titles = {}
		self._units = {}
		for d in self._diags:
			self._titles.update({ d:"??" })
			unitss = "??"
			output = self._info[d]["output"]
			if   output == "density"                        : self._titles[d] = "Number density"
			elif output == "charge_density"                 : self._titles[d] = "Charge density"
			elif output[:-1] == "current_density_"          : self._titles[d] = "J"+output[-1]
			elif output == "p_density"                      : self._titles[d] = "P density"
			elif output[2:] == "_density" and output[0]=="p": self._titles[d] = "P"+output[1]+" density"
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
			self._units[d] = " [ "+unitss+" ]"
		# finish title creation
		if len(self._diags) == 1:
			self._title = self._titles[self._diags[0]] + self._units[self._diags[0]]
		else:
			self._title = self.operation
			for d in self._diags:
				self._title = self._title.replace("#"+str(d), self._titles[d])
		if self._data_log: self._title = "Log[ "+self._title+" ]"
	
		# If any spatial dimension did not appear, then count it for calculating the correct density
		if ndim>=1 and not spatialaxes["x"]: units_coeff /= ncels[0]
		if ndim>=2 and not spatialaxes["y"]: units_coeff /= ncels[1]
		if ndim==3 and not spatialaxes["z"]: units_coeff /= ncels[2]
		units_coeff *= coeff_density
	
		# Calculate the array that represents the bins sizes in order to get units right.
		# This array will be the same size as the plotted array
		if len(plot_diff)==0:
			self._bsize = 1.
		elif len(plot_diff)==1:
			self._bsize = plot_diff[0]
		else:
			self._bsize = self._np.prod( self._np.array( self._np.meshgrid( *plot_diff ) ), axis=0)
			self._bsize = self._bsize.transpose()
		self._bsize /= units_coeff
		
		# Finish constructor
		self.valid = True
	
	
	# Gets info about diagnostic number "diagNumber"
	def _getInfo(self,diagNumber):
		# path to the file
		file = self._results_path+'/ParticleDiagnostic'+str(diagNumber)+'.h5'
		# if no file, return
		if not self._ospath.isfile(file): return False
		# open file
		f = self._h5py.File(file, 'r')
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
	
	
	# Prints the info obtained by the function "getInfo"
	@staticmethod
	def _printInfo(info):
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
		if not self._validate(): return
		for d in self._diags:
			self._printInfo(self._info[d])
		if len(self.operation)>2: print "Operation : "+self.operation
		for ax in self._axes:
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
				file = self._results_path+'/ParticleDiagnostic'+str(diagNumber)+'.h5'
				f = self._h5py.File(file, 'r')
			except:
				print "Cannot open file "+file
				return self._np.array([])
			items = f.items()
			ntimes = len(items)
			times = self._np.zeros(ntimes)
			for i in range(ntimes):
				times[i] = int(items[i][0].strip("timestep")) # fill the "times" array with the available timesteps
			f.close()
			return times
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self.times:
			print "Timestep "+t+" not found in this diagnostic"
			return []
		# Get arrays from all requested diagnostics
		A = {}
		for d in self._diags:
			# get data
			index = self._data[d][t]
			A.update({ d:self._np.reshape(self._h5items[d][index][1],self._shape) })
			# Apply the slicing
			for iaxis in range(self._naxes):
				axis = self._axes[iaxis]
				if "slice" in axis:
					A[d] = self._np.delete(A[d], axis["slice"], axis=iaxis) # remove parts outside of the slice
					A[d][self._np.isnan(A[d])] = 0.
					A[d] = self._np.sum(A[d], axis=iaxis, keepdims=True) # sum over the slice
			A[d] = self._np.squeeze(A[d]) # remove sliced axes
			# Divide by the bins size
			A[d] /= self._np.squeeze(self._bsize)
		# Calculate operation
		data_operation = self.operation
		for d in reversed(self._diags):
			data_operation = data_operation.replace("#"+str(d),"A["+str(d)+"]")
		exec("A = "+data_operation) in None
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A




# -------------------------------------------------------------------
# Class for fields diagnostics
# -------------------------------------------------------------------
class Field(Diagnostic):
	
	# This is the constructor, which creates the object
	def _init(self, field=None, timesteps=None, slice=None,
				 units="code", data_log=False, **kwargs):
		
		if not self.Smilei.valid: return None
		if field is None:
			fields = self.getFields()
			if len(fields)>0:
				print "Printing available fields:"
				print "--------------------------"
				l = (len(fields)/3) * 3
				if l>0:
					print '\n'.join(['\t\t'.join(list(i)) for i in self._np.reshape(fields[:l],(-1,3))])
				print '\t\t'.join(fields[l:])
			else:
				print "No fields found in '"+self._results_path+"'"
			return None
		
		
		# Get info from the input file and prepare units
		try:
			ndim               = self._read_ndim()
			ncels, cell_length = self._read_ncels_cell_length(ndim)
			self.timestep      = self._read_timestep()
		except:
			return None
		
		if units == "nice":
			try   : wavelength_SI = self._read_wavelength_SI()
			except: return None
			cell_length *= 1e2*wavelength_SI/(2.*self._np.pi) # in cm
			cell_volume = self._np.prod(cell_length)
			coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 * cell_volume # in e/cm^3
			coeff_current = coeff_density * 4.803e-9 # in A/cm^2
			self._coeff_time = self.timestep * wavelength_SI/3.e8 # in seconds
			self._time_units = " s"
		elif units == "code":
			coeff_density = 1. # in nc
			coeff_current = 1. # in e*c*nc
			self._coeff_time = self.timestep # in 1/w
			self._time_units = " $1/\omega$"
		
		# Get available times
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
			print "No fields found in Fields.h5"
			return
		
		# Get available fields
		fields = self.getFields()
		sortedfields = reversed(sorted(fields, key = len))
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Parse the `field` argument
		self.operation = field
		for f in sortedfields:
			i = fields.index(f)
			self.operation = self.operation.replace(f,"#"+str(i))
		requested_fields = self._re.findall("#\d+",self.operation)
		if len(requested_fields) == 0:
			print "Could not find any existing field in `"+field+"`"
			return None
		self._fieldn = [ int(f[1:]) for f in requested_fields ] # indexes of the requested fields
		self._fieldn = list(set(self._fieldn))
		self._fieldname = [ fields[i] for i in self._fieldn ] # names of the requested fields
	
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			print "Argument `slice` must be a dictionary"
			return None
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self._data_log = data_log
	
		# Get the shape of fields
		self._file = self._results_path+'/Fields.h5'
		f = self._h5py.File(self._file, 'r')
		self._h5items = f.values()
		self._shape = self._np.double(self._h5items[0].values()[0]).shape
		for n in self._fieldn:
			s = self._np.double(self._h5items[0].values()[n]).shape
			self._shape = self._np.min((self._shape, s), axis=0)
		
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with indices to the data arrays
		self._data = {}
		for i,t in enumerate(self.times):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				ts = self._np.array(self._np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self.times = self.times[ (self.times>=ts[0]) * (self.times<=ts[1]) ]
				elif ts.size==1:
					# get nearest time
					self.times = self._np.array([self.times[(self._np.abs(self.times-ts)).argmin()]])
				else:
					raise
			except:
				print "Argument `timesteps` must be one or two non-negative integers"
				return None
	
		# Need at least one timestep
		if self.times.size < 1:
			print "Timesteps not found"
			return None
	
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values
		self._naxes = ndim
		self._plot_type = []
		self._plot_label = []
		self._plot_centers = []
		self._plot_shape = []
		self._plot_log   = []
		self._sliceinfo = {}
		self._slices = [None]*ndim
		for iaxis in range(self._naxes):
			centers = self._np.linspace(0., self._shape[iaxis]*cell_length[iaxis], self._shape[iaxis])
			label = {0:"x", 1:"y", 2:"z"}[iaxis]
			axisunits = "[code units]"
			if units == "nice": axisunits = "[cm]"
		
			if label in slice:
				# if slice is "all", then all the axis has to be summed
				if slice[label] == "all":
					indices = self._np.arange(self._shape[iaxis])
				# Otherwise, get the slice from the argument `slice`
				else:
					try:
						s = self._np.double(slice[label])
						if s.size>2 or s.size<1: raise
					except:
						print "Slice along axis "+label+" should be one or two floats"
						return None
					if s.size==1:
						indices = self._np.array([(self._np.abs(centers-s)).argmin()])
					elif s.size==2:
						indices = self._np.nonzero( (centers>=s[0]) * (centers<=s[1]) )[0]
					if indices.size == 0:
						print "Slice along "+label+" is out of the box range"
						return None
					if indices.size == 1:
						self._sliceinfo.update({ label:"Sliced at "+label+" = "+str(centers[indices])+" "+axisunits })
					else:
						self._sliceinfo.update({ label:"Sliced for "+label
							+" from "+str(centers[indices[ 0]])+" to "+str(centers[indices[-1]])+" "+axisunits })
				# convert the range of indices into their "conjugate"
				self._slices[iaxis] = self._np.delete(self._np.arange(self._shape[iaxis]), indices)
			else:
				self._plot_type   .append(label)
				self._plot_shape  .append(self._shape[iaxis])
				self._plot_centers.append(centers)
				self._plot_label  .append(label+" "+axisunits)
				self._plot_log    .append(False)
		
		if len(self._plot_centers) > 2:
			print "Cannot plot in "+str(len(self._plot_shape))+"d. You need to 'slice' some axes."
			return
	
		# Build units
		self._titles = {}
		self._fieldunits = {}
		self._unitscoeff = {}
		for f in self._fieldname:
			i = fields.index(f)
			self._fieldunits.update({ i:"??" })
			self._unitscoeff.update({ i:1 })
			self._titles    .update({ i:"??" })
			if units == "nice":
				self._fieldunits[i] = " ("+{"B":"T"  ,"E":"V/m"  ,"J":"A/cm$^2$"   ,"R":"e/cm$^3$"   }[f[0]]+")"
				self._unitscoeff[i] =      {"B":10710,"E":3.21e12,"J":coeff_current,"R":coeff_density}[f[0]]
				self._titles    [i] = f
			else:
				self._fieldunits[i] = " in units of "+{"B":"$m_e\omega/e$","E":"$m_ec\omega/e$","J":"$ecn_c$"    ,"R":"$n_c$"      }[f[0]]
				self._unitscoeff[i] =                 {"B":1              ,"E":1               ,"J":coeff_current,"R":coeff_density}[f[0]]
				self._titles    [i] = f
		# finish title creation
		if len(self._fieldname) == 1:
			self._title = self._titles[self._fieldn[0]] + self._fieldunits[self._fieldn[0]]
		else:
			self._title = self.operation
			for n in self._fieldn:
				self._title = self._title.replace("#"+str(n), self._titles[n])
	
		# Finish constructor
		self.valid = True

	# Method to print info on included fields
	def info(self):
		if not self._validate(): return
		print self._title
		#todo
		return

	# get all available fields, sorted by name length
	def getFields(self):
		try:
			file = self._results_path+'/Fields.h5'
			f = self._h5py.File(file, 'r')
		except:
			print "Cannot open file "+file
			return []
		try:
			fields = f.values()[0].keys() # list of fields
		except:
			fields = []
		f.close()
		return fields
	
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:
			file = self._results_path+'/Fields.h5'
			f = self._h5py.File(file, 'r')
		except:
			print "Cannot open file "+file
			return self._np.array([])
		times = self._np.double(f.keys())
		f.close()
		return times
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self.times:
			print "Timestep "+t+" not found in this diagnostic"
			return []
		# Get arrays from requested field
		# get data
		index = self._data[t]
		C = {}
		op = "A=" + self.operation
		for n in reversed(self._fieldn): # for each field in operation
			B = self._np.double(self._h5items[index].values()[n]) # get array
			B *= self._unitscoeff[n]
			for axis, size in enumerate(self._shape):
				l = self._np.arange(size, B.shape[axis])
				B = self._np.delete(B, l, axis=axis) # remove extra cells if necessary
			C.update({ n:B })
			op = op.replace("#"+str(n), "C["+str(n)+"]")
		# Calculate the operation
		exec op in None
		# Apply the slicing
		for iaxis in range(self._naxes):
			if self._slices[iaxis] is None: continue
			A = self._np.delete(A, self._slices[iaxis], axis=iaxis) # remove parts outside of the slice
			A = self._np.mean(A, axis=iaxis, keepdims=True) # sum over the slice
		A = self._np.squeeze(A) # remove sliced axes
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A




# -------------------------------------------------------------------
# Class for scalars
# -------------------------------------------------------------------
class Scalar(Diagnostic):
	
	# This is the constructor, which creates the object
	def _init(self, scalar=None, timesteps=None,
				 units="code", data_log=False, **kwargs):
		
		if not self.Smilei.valid: return None
		if scalar is None:
			scalars = self.getScalars()
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
				print "No scalars found in '"+self._results_path+"'"
			return None
		
		# Get info from the input file and prepare units
		try:
			self.timestep      = self._read_timestep()
		except:
			return None
		
		if units == "nice":
			try   : wavelength_SI = self._read_wavelength_SI()
			except: return None
			self._coeff_time = self.timestep * wavelength_SI/3.e8/(2.*self.np.pi) # in seconds
			self._time_units = " s"
		elif units == "code":
			self._coeff_time  = self.timestep
			self._time_units = " $1/\omega$"
	
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
		self._scalarn = scalars.index(scalar) # index of the requested scalar
		self._scalarname = scalar
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Already get the data from the file
		# Loop file line by line
		self.times = []
		self._values = []
		file = self._results_path+'/scalars.txt'
		f = open(file, 'r')
		for line in f:
			line = line.strip()
			if line[0]=="#": continue
			line = line.split()
			self.times .append( int( self._np.round(float(line[0]) / float(self.timestep)) ) )
			self._values.append( float(line[self._scalarn+1]) )
		self.times  = self._np.array(self.times )
		self._values = self._np.array(self._values)
		f.close()
		
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with the index to each time
		self._data = {}
		for i,t in enumerate(self.times):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				ts = self._np.array(self._np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self.times = self.times[ (self.times>=ts[0]) * (self.times<=ts[1]) ]
				elif ts.size==1:
					# get nearest time
					self.times = self._np.array([self.times[(self._np.abs(self.times-ts)).argmin()]])
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
		self._naxes = 0
		self._plot_type = []
		self._plot_label = []
		self._plot_centers = []
		self._plot_shape = []
		self._plot_log   = []
		self._slices = []
		# Build units
		self._scalarunits = "unknown units"
		self._unitscoeff = 1.
		self._title = "??"
		if units == "nice":
			self._title      = scalar +"( "+self._scalarunits+" )" # todo
		else:
			self._title      = scalar +"( "+self._scalarunits+" )" # todo
		if data_log: self._title = "Log[ "+self._title+" ]"
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included scalars
	def info(self):
		if not self._validate(): return
		print "Scalar "+self._scalarname,
		#todo
		return
	
	# get all available scalars
	def getScalars(self):
		try:
			file = self._results_path+'/scalars.txt'
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
	
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		return self.times
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self.times:
			print "Timestep "+t+" not found in this diagnostic"
			return []
		# Get value at selected time
		A = self._values[ self._data[t] ]
		A *= self._unitscoeff
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A




# -------------------------------------------------------------------
# Class for probe diagnostics
# -------------------------------------------------------------------
class Probe(Diagnostic):

	# This is the constructor, which creates the object
	def _init(self, probeNumber=None, field=None, timesteps=None, slice=None,
				 units="code", data_log=False, **kwargs):
		
		if not self.Smilei.valid: return None
		if probeNumber is None:
			probes = self.getProbes()
			if len(probes)>0:
				print "Printing available probes:"
				print "--------------------------"
				for p in probes:
					self._printInfo(self._getInfo(p))
			else:
				print "No probes found in '"+self._results_path+"'"
			return None
		if field is None:
			print "Printing available fields for probes:"
			print "-------------------------------------"
			print "Ex Ey Ez Bx By Bz Jx Jy Jz Rho"
			return None
		
		
		self.probeNumber  = probeNumber
		self._file = self._results_path+"/Probes.h5"
		f = self._h5py.File(self._file, 'r')
		self._h5items = f.values()
		
		# Get info from the input file and prepare units
		try:
			ndim               = self._read_ndim()
			ncels, cell_length = self._read_ncels_cell_length(ndim)
			self.timestep      = self._read_timestep()
		except:
			return None
		
		if units == "nice":
			try   : wavelength_SI = self._read_wavelength_SI()
			except: return None
			cell_length *= 1e2*wavelength_SI/(2.*self._np.pi) # in cm
			cell_volume = self._np.prod(cell_length)
			coeff_distance = 1e2*wavelength_SI/(2.*self._np.pi) # in cm
			coeff_density = 1.11e21 / (wavelength_SI/1e-6)**2 * cell_volume # in e/cm^3
			coeff_current = coeff_density * 4.803e-9 # in A/cm^2
			self._coeff_time = self.timestep * wavelength_SI/3.e8 # in seconds
			self._time_units = " s"
		elif units == "code":
			coeff_distance = 1 # in c/w
			coeff_density = 1. # in nc
			coeff_current = 1. # in e*c*nc
			self._coeff_time = self.timestep # in 1/w
			self._time_units = " $1/\omega$"
		
		# Get available times
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
			print "No probes found in Probes.h5"
			return
		
		# Get available fields
		fields = ["Ex","Ey","Ez","Bx","By","Bz","Jx","Jy","Jz","Rho"]
		sortedfields = reversed(sorted(fields, key = len))
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Parse the `field` argument
		self.operation = field
		for f in sortedfields:
			i = fields.index(f)
			self.operation = self.operation.replace(f,"#"+str(i))
		requested_fields = self._re.findall("#\d+",self.operation)
		if len(requested_fields) == 0:
			print "Could not find any existing field in `"+field+"`"
			return None
		self._fieldn = [ int(f[1:]) for f in requested_fields ] # indexes of the requested fields
		self._fieldn = list(set(self._fieldn))
		self._fieldname = [ fields[i] for i in self._fieldn ] # names of the requested fields
		
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			print "Argument `slice` must be a dictionary"
			return None
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Get the shape of the probe
		self._info = self._getMyInfo()
		self._shape = self._info["shape"]
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with indices to the data arrays
		self._data = {}
		for i,t in enumerate(self.times):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				ts = self._np.array(self._np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self.times = self.times[ (self.times>=ts[0]) * (self.times<=ts[1]) ]
				elif ts.size==1:
					# get nearest time
					self.times = self._np.array([self.times[(self._np.abs(self.times-ts)).argmin()]])
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
		self._naxes = self._shape.size
		self._plot_type = []
		self._plot_label = []
		self._plot_centers = []
		self._plot_shape = []
		self._plot_log   = []
		self._sliceinfo = {}
		self._slices = [None]*ndim
		for iaxis in range(self._naxes):
		
			# calculate grid points locations
			p0 = self._info["p0"            ] # reference point
			pi = self._info["p"+str(iaxis+1)] # end point of this axis
			centers = self._np.zeros((self._shape[iaxis],p0.size))
			for i in range(p0.size):
				centers[:,i] = self._np.linspace(p0[i],pi[i],self._shape[iaxis])
			centers *= coeff_distance
		
			label = {0:"axis1", 1:"axis2", 2:"axis3"}[iaxis]
			axisunits = "[code units]"
			if units == "nice": axisunits = "[cm]"
		
			if label in slice:
				# if slice is "all", then all the axis has to be summed
				if slice[label] == "all":
					indices = self._np.arange(self._shape[iaxis])
				# Otherwise, get the slice from the argument `slice`
				else:
					indices = self._np.arange(self._shape[iaxis])
					try:
						s = self._np.double(slice[label])
						if s.size>2 or s.size<1: raise
					except:
						print "Slice along axis "+label+" should be one or two floats"
						return None
					if s.size==1:
						indices = self._np.array([(self._np.abs(indices-s)).argmin()])
					elif s.size==2:
						indices = self._np.nonzero( (indices>=s[0]) * (indices<=s[1]) )[0]
					if indices.size == 0:
						print "Slice along "+label+" is out of the box range"
						return None
					if indices.size == 1:
						self._sliceinfo.update({ label:"Sliced at "+label+" = "+str(indices[0]) })
					else:
						self._sliceinfo.update({ label:"Sliced for "+label+" from "+str(indices[0])+" to "+str(indices[-1]) })
				# convert the range of indices into their "conjugate"
				self._slices[iaxis] = self._np.delete(self._np.arange(self._shape[iaxis]), indices)
			else:
				self._plot_type   .append(label)
				self._plot_shape  .append(self._shape[iaxis])
				self._plot_centers.append(centers)
				self._plot_label  .append(label+" "+axisunits)
				self._plot_log    .append(False)
			
		
		if len(self._plot_centers) > 2:
			print "Cannot plot in "+str(len(self._plot_shape))+"d. You need to 'slice' some axes."
			return
		
		# Special case in 1D: we convert the point locations to scalar distances
		if len(self._plot_centers) == 1:
			self._plot_centers[0] = self._np.sqrt(self._np.sum((self._plot_centers[0]-self._plot_centers[0][0])**2,axis=1))
		# Special case in 2D: we have to prepare for pcolormesh instead of imshow
		elif len(self._plot_centers) == 2:
			p1 = self._plot_centers[0] # locations of grid points along first dimension
			d = self._np.diff(p1, axis=0) # separation between the points
			p1 = self._np.vstack((p1, p1[-1,:])) # add last edges at the end of box
			p1[1:-1] -= d/2 # move points by one half
			p2 = self._plot_centers[1] # locations of grid points along second dimension
			d = self._np.diff(p2, axis=0) # separation between the points
			p2 = self._np.vstack((p2, p2[-1,:])) # add last edges at the end of box
			p2[1:-1] -= d/2 # move points by one half
			# Now p1 and p2 contain edges grid points along the 2 dimensions
			# We have to convert into X and Y 2D arrays (similar to meshgrid)
			X = self._np.zeros((p1.shape[0], p2.shape[0]))
			Y = self._np.zeros((p1.shape[0], p2.shape[0]))
			for i in range(p2.shape[0]):
				X[:,i] = p1[:,0] + p2[i,0]-p2[0,0]
				Y[:,i] = p1[:,1] + p2[i,1]-p2[0,1]
			self._plot_edges = [X, Y]
			self._plot_label = ["x "+axisunits, "y "+axisunits]
		
		
		# Build units
		self._titles = {}
		self._fieldunits = {}
		self._unitscoeff = {}
		for f in self._fieldname:
			i = fields.index(f)
			self._fieldunits.update({ i:"??" })
			self._unitscoeff.update({ i:1 })
			self._titles    .update({ i:"??" })
			if units == "nice":
				self._fieldunits[i] = " ("+{"B":"T"  ,"E":"V/m"  ,"J":"A/cm$^2$"   ,"R":"e/cm$^3$"   }[f[0]]+")"
				self._unitscoeff[i] =      {"B":10710,"E":3.21e12,"J":coeff_current,"R":coeff_density}[f[0]]
				self._titles    [i] = f
			else:
				self._fieldunits[i] = " in units of "+{"B":"$m_e\omega/e$","E":"$m_ec\omega/e$","J":"$ecn_c$"    ,"R":"$n_c$"      }[f[0]]
				self._unitscoeff[i] =                 {"B":1              ,"E":1               ,"J":coeff_current,"R":coeff_density}[f[0]]
				self._titles    [i] = f
		# finish title creation
		if len(self._fieldname) == 1:
			self._title = self._titles[self._fieldn[0]] + self._fieldunits[self._fieldn[0]]
		else:
			self._title = self.operation
			for n in self._fieldn:
				self._title = self._title.replace("#"+str(n), self._titles[n])
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included probe
	def info(self):
		if not self._validate(): return
		self._printInfo(self._getMyInfo())
	# Method to print info previously obtained with getInfo
	@staticmethod
	def _printInfo(info):
		print ("Probe #"+str(info["probeNumber"])+": "+str(info["dimension"])+"-dimensional,"
			+" every "+str(info["every"])+" timesteps")
		i = 0
		while "p"+str(i) in info:
			print "p"+str(i)+" = "+" ".join(info["p"+str(i)].astype(str).tolist())
			i += 1
		if info["shape"].size>0:
			print "number = "+" ".join(info["shape"].astype(str).tolist())
	
	# Method to get info on a given probe
	def _getInfo(self, probeNumber):
		try:
			file = self._results_path+'/Probes.h5'
			f = self._h5py.File(file, 'r')
		except:
			print "Cannot open file "+file
			return {}
		try:
			probes = [int(key.strip("p")) for key in f.keys()]
			k = probes.index(probeNumber)
			probe = f.values()[k]
		except:
			print "Cannot find probe "+str(probeNumber)+" in file "+file
			return {}
		out = {}
		out.update({"probeNumber":probeNumber, "dimension":probe.attrs["dimension"],
			"every":probe.attrs["every"]})
		i = 0
		while "p"+str(i) in probe.keys():
			k = probe.keys().index("p"+str(i))
			out.update({ "p"+str(i):self._np.array(probe.values()[k]) })
			i += 1
		k = probe.keys().index("number")
		out.update({ "shape":self._np.array(probe.values()[k]) })
		f.close()
		return out
	def _getMyInfo(self):
		return self._getInfo(self.probeNumber)
	
	# get all available fields, sorted by name length
	def getProbes(self):
		try:
			file = self._results_path+'/Probes.h5'
			f = self._h5py.File(file, 'r')
		except:
			print "Cannot open file "+file
			return []
		try:
			probes = [int(key.strip("p")) for key in f.keys()] # list of probes
		except:
			probes = []
		f.close()
		return probes
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:
			f = self._h5py.File(self._file, 'r')
		except:
			print "Cannot open file "+self._file
			return self._np.array([])
		try:
			probes = [int(key.strip("p")) for key in f.keys()]
			k = probes.index(self.probeNumber)
			probe = f.values()[k]
		except:
			print "Cannot find probe "+str(self.probeNumber)+" in file "+file
			return self._np.array([])
		times = []
		for key in probe.keys():
			try   : times.append( int(key) )
			except: pass
		f.close()
		return self._np.double(times)
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self.times:
			print "Timestep "+t+" not found in this diagnostic"
			return []
		# Get arrays from requested field
		# Open file
		f = self._h5py.File(self._file, 'r')
		# get data
		index = self._data[t]
		C = {}
		op = "A=" + self.operation
		for n in reversed(self._fieldn): # for each field in operation
			B = self._np.double(self._h5items[self.probeNumber].values()[index][:,n]) # get array
			B = self._np.reshape(B, self._shape) # reshape array because it is flattened in the file
			B *= self._unitscoeff[n]
			C.update({ n:B })
			op = op.replace("#"+str(n), "C["+str(n)+"]")
		f.close()
		# Calculate the operation
		exec op in None
		# Apply the slicing
		for iaxis in range(self._naxes):
			if self._slices[iaxis] is None: continue
			A = self._np.delete(A, self._slices[iaxis], axis=iaxis) # remove parts outside of the slice
			A = self._np.mean(A, axis=iaxis, keepdims=True) # average over the slice
		A = self._np.squeeze(A) # remove sliced axes
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
	
	# Overloading a plotting function in order to use pcolormesh instead of imshow
	def _animateOnAxes_2D(self, ax, A):
		# first, we remove kwargs that are not supported by pcolormesh
		kwargs = dict(self.options.image)
		for kwarg in self.options.image:
			if kwarg not in ["cmap"]: del kwargs[kwarg]
		im = ax.pcolormesh(self._plot_edges[0], self._plot_edges[1], self._np.flipud(A.transpose()),
			vmin = self.options.vmin, vmax = self.options.vmax, **kwargs)
		return im



def multiPlot(*Diags, **kwargs):
	""" multiplot(Diag1, Diag2, ..., figure=1, shape=None)
	
	Plots simultaneously several diagnostics.
	
	Parameters:
	-----------
	Diag1, Diag2, ... : Several objects of classes 'Scalar', 'Field', 'Probe' or 'ParticleDiagnostic' 
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	nDiags = len(Diags)
	# Verify Diags are valid
	if nDiags == 0: return
	for Diag in Diags:
		if not Diag.valid: return
	# Gather all times
	alltimes = np.unique(np.concatenate([Diag.times*Diag.timestep for Diag in Diags]))
	# Get keyword arguments
	shape  = kwargs.pop("shape" , None)
	# Determine whether to plot all cases on the same axes
	sameAxes = False
	if shape is None or shape == [1,1]:
		sameAxes = True
		for Diag in Diags:
			if Diag.dim()==0 and Diags[0].dim()==0:
				continue
			if Diag.dim()!=1 or Diag._plot_type!=Diags[0]._plot_type:
				sameAxes = False
				break
	if not sameAxes and shape == [1,1]:
		print "Cannot have shape=[1,1] with these diagnostics"
		return
	# Determine the shape
	if sameAxes: shape = [1,1]
	if shape is None: shape = [nDiags,1]
	nplots = np.array(shape).prod()
	if not sameAxes and nplots != nDiags:
		print "The 'shape' argument is incompatible with the number of diagnostics:"
		print "  "+str(nDiags)+" diagnostics do not fit "+str(nplots)+" plots"
		return
	# Make the figure
	if "facecolor" not in kwargs: kwargs.update({ "facecolor":"w" })
	options = Options(**kwargs)
	fig = plt.figure(**options.figure0)
	fig.set(**options.figure1) # Apply figure kwargs
	fig.clf()
	fig.subplots_adjust(wspace=0.5, hspace=0.5, bottom=0.15)
	ax = []
	xmin =  float("inf")
	xmax = -float("inf")
	c = plt.matplotlib.rcParams['axes.color_cycle']
	for i in range(nplots):
		ax.append( fig.add_subplot(shape[0], shape[1], i+1) )
	for i, Diag in enumerate(Diags):
		if sameAxes: Diag._ax = ax[0]
		else       : Diag._ax = ax[i]
		Diag._artist = None
		try:
			l = Diag.limits()[0]
			xmin = min(xmin,l[0])
			xmax = max(xmax,l[1])
		except:
			pass
		if "color" not in Diag.options.plot:
			Diag.options.plot.update({ "color":c[i%len(c)] })
	# Static plot
	if sameAxes and len(Diags[0]._plot_shape)==0:
		for Diag in Diags:
			Diag._artist = Diag._plotVsTime(Diag._ax)
		fig.canvas.draw()
		plt.show()
	# Animated plot
	else:
		# Loop all times
		for time in alltimes:
			for Diag in Diags:
				t = np.round(time/Diag.timestep) # convert time to timestep
				if t in Diag.times:
					if sameAxes:
						if Diag._artist is not None: Diag._artist.remove()
					else:
						Diag._ax.cla()
					Diag._artist = Diag._animateOnAxes(Diag._ax, t)
					if sameAxes:
						Diag._ax.set_xlim(xmin,xmax)
			fig.canvas.draw()
			plt.show()
		return
	

