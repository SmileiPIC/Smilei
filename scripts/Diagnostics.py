# -----------------------------------------------------------------------
# HOW TO VIEW DIAGNOSTICS                -    F. Perez - 2015
# -----------------------------------------------------------------------
# Check out the documentation in the doc/Sphinx directory.
# It must be compiled with the Sphinx software (sphinx-doc.org).


def setMatplotLibBackend(show=True):
	import matplotlib, sys
	usingAgg = (matplotlib.get_backend().lower() == "agg")
	if not show and not usingAgg:
		if "matplotlib.pyplot" in sys.modules:
			print "WARNING: 'show=False' requires you restart python."
		else:
			matplotlib.use("Agg")
	if show and usingAgg:
		if "matplotlib.pyplot" in sys.modules:
			print "WARNING: 'show=False' was set earlier. Restart python if you want figures to appear."
	#print matplotlib.get_backend()


class Smilei(object):
	""" Smilei(results_path=".", show=True)
	
	Import Smilei simulation information.
	
	* `results_path` specifies where the simulation results are located.
	Omit this argument if you are already in the results path.
	
	* `show` can be set to False to prevent figures to actually appear on screen.
	
	"""
	
	def __init__(self, results_path=".", show=True):
		self.valid = False
		# Import packages
		import h5py
		import numpy as np
		import os.path, glob, re, sys
		setMatplotLibBackend(show=show)
		import matplotlib.pyplot
		import matplotlib.pylab as pylab
		pylab.ion()
		# Transfer packages to local attributes
		self._results_path = results_path
		self._h5py = h5py
		self._np = np
		self._ospath = os.path
		self._glob = glob.glob
		self._re = re
		self._plt = matplotlib.pyplot
		self.reload()
	
	def reload(self):
		self.valid = False
		# Verify that results_path is valid
		if not self._ospath.isdir(self._results_path):
			print "Could not find directory "+self._results_path
			return
		if len(self._glob(self._results_path+"/smilei.py"))==0:
			print "Could not find an input file in directory "+self._results_path
			return
		# Fetch the python namelist
		namespace={}
		execfile(self._results_path+'/smilei.py',namespace) # execute the namelist into an empty namespace
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
		""" Scalar(scalar=None, timesteps=None, units=[""], data_log=False)
		
		Import and analyze a scalar diagnostic from a Smilei simulation
		
		Parameters:
		-----------
		scalar : an available scalar name. (optional)
			To get a list of available scalars, simply omit this argument.
		timesteps : int or [int, int] (optional)
			If omitted, all timesteps are used.
			If one number  given, the nearest timestep available is used.
			If two numbers given, all the timesteps in between are used.
		units : A units specification such as ["m","second"]
		data_log : True or False    (optional)
			If True, then log10 is applied to the output array before plotting.
		
		Usage:
		------
			scalar = S.Scalar(...)
		where S is a `Smilei` object.
		"""
		return Scalar(self, *args, **kwargs)
	def Field(self, *args, **kwargs):
		""" Field(field=None, timesteps=None, slice=None, units=[""], data_log=False)
		
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
		units : A units specification such as ["m","second"]
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
		""" Probe(probeNumber=None, field=None, timesteps=None, slice=None, units=[""], data_log=False)
		
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
		units : A units specification such as ["m","second"]
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
		""" ParticleDiagnostic(diagNumber=None, timesteps=None, slice=None, units=[""], data_log=False)
		
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
		units : A units specification such as ["m","second"]
		data_log : True or False    (optional)
			If True, then log10 is applied to the output array before plotting.
		
		Usage:
		------
			diag = S.ParticleDiagnostic(...) # S is a Smilei object
			diag.get()
			diag.plot()
		"""
		return ParticleDiagnostic(self, *args, **kwargs)
	def TrackParticles(self, *args, **kwargs):
		""" TrackParticles(species=None, select="", axes=[], timesteps=None, units=[""], skipAnimation=False)
		
		Import and analyze tracked particles from a Smilei simulation
		
		Parameters:
		-----------
		species : name of a tracked species. (optional)
			To get a list of available tracked species, simply omit this argument.
		select: Instructions for selecting particles among those available.
			Syntax 1: select="any(times, condition)"
			Syntax 2: select="all(times, condition)"
			`times` is a selection of timesteps t, for instance `t>50`.
			`condition` is a condition on particles properties (x, y, z, px, py, pz), for instance `px>0`.
			Syntax 1 selects particles satisfying `condition` for at least one of the `times`.
			Syntax 2 selects particles satisfying `condition` at all `times`.
			Example: select="all(t<40, px<0.1)" selects particles that kept px<0.1 until timestep 40.
			Example: select="any(t>0, px>1.)" selects particles that reached px>1 at some point.
			It is possible to make logical operations: + is OR; * is AND; - is NOT.
			Example: select="any((t>30)*(t<60), px>1) + all(t>0, (x>1)*(x<2))"
		timesteps : int or [int, int] (optional)
			If omitted, all timesteps are used.
			If one number  given, the nearest timestep available is used.
			If two numbers given, all the timesteps in between are used.
		units : A units specification such as ["m","second"]
		axes: A list of axes for plotting the trajectories.
			Each axis is "x", "y", "z", "px", "py" or "pz".
			Example: axes = ["x"] corresponds to x versus time.
			Example: axes = ["x","y"] correspond to 2-D trajectories.
			Example: axes = ["x","px"] correspond to phase-space trajectories.
		skipAnimation: when True, the plot() will directly show the full trajectory.
		
		Usage:
		------
			diag = S.ParticleDiagnostic(...) # S is a Smilei object
			diag.get()
			diag.plot()
		"""
		return TrackParticles(self, *args, **kwargs)


class Options(object):
	""" Class to contain matplotlib plotting options """
	
	def __init__(self, **kwargs):
		self.figure  = 1
		self.xfactor = None
		self.xmin    = None
		self.xmax    = None
		self.yfactor = None
		self.ymin    = None
		self.ymax    = None
		self.vfactor = None
		self.vmin    = None
		self.vmax    = None
		self.skipAnimation = False
		self.streakPlot = False
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
		self.xfactor  = kwargs.pop("xfactor",self.xfactor  )
		self.xmin     = kwargs.pop("xmin"   ,self.xmin  )
		self.xmax     = kwargs.pop("xmax"   ,self.xmax  )
		self.yfactor  = kwargs.pop("yfactor",self.yfactor  )
		self.ymin     = kwargs.pop("ymin"   ,self.ymin  )
		self.ymax     = kwargs.pop("ymax"   ,self.ymax  )
		self.vfactor  = kwargs.pop("vfactor",self.vfactor  )
		self.vmin     = kwargs.pop("vmin"   ,kwargs.pop("data_min",self.vmin))
		self.vmax     = kwargs.pop("vmax"   ,kwargs.pop("data_max",self.vmax))
		self.skipAnimation = kwargs.pop("skipAnimation", self.skipAnimation)
		self.streakPlot    = kwargs.pop("streakPlot"   , self.streakPlot   )
		# Second, we manage all the other arguments that are directly the ones of matplotlib
		for kwa, val in kwargs.iteritems():
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


class Units(object):
	""" Units()
	
	Class to handle units smartly. Based on the *pint* package.
	"""
	
	def __init__(self, *args, **kwargs):
		# All args are parsed
		self.requestedUnits = []
		self.requestedX = ""
		self.requestedY = ""
		self.requestedV = ""
		for a in args:
			if type(a) is str:
				self.requestedUnits.append( a )
			else:
				raise TypeError("Arguments of Units() should be strings")
		for kwa, val in kwargs.iteritems():
			if type(val) is not str:
				raise TypeError("Arguments of Units() should be strings")
			if   kwa == "x": self.requestedX = val
			elif kwa == "y": self.requestedY = val
			elif kwa == "v": self.requestedV = val
			else: raise TypeError("Units() got an unexpected keyword argument '"+kwa+"'")
		
		# We try to import the pint package
		self.UnitRegistry = None
		try:
			from pint import UnitRegistry
		except:
			print "WARNING: you do not have the *pint* package, so you cannot modify units."
			print "       : The results will stay in code units."
			return
		self.UnitRegistry = UnitRegistry
	
	def _divide(self,units1, units2):
		division = self.ureg("("+units1+") / ("+units2+")").to_base_units()
		if not division.dimensionless: raise
		return division.magnitude or 1., units2
	
	def _convert(self, knownUnits, requestedUnits):
		if knownUnits:
			if requestedUnits:
				try:
					return self._divide(knownUnits,requestedUnits)
				except:
					print "WARNING: cannot convert units to <"+requestedUnits+">"
					print "       : Conversion discarded."
			else:
				for units in self.requestedUnits:
					try   : return self._divide(knownUnits,units)
					except: pass
			val = self.ureg(knownUnits)
			return val.magnitude or 1., u"{0.units:P}".format(val)
		return 1., ""
	
	def prepare(self, wavelength_SI=None, xunits="", yunits="", vunits="", tunits=""):
		if self.UnitRegistry:
			if wavelength_SI:
				# Load pint's default unit registry
				self.ureg = self.UnitRegistry()
				# Define code units
				self.ureg.define("V_r = speed_of_light"                   ) # velocity
				self.ureg.define("W_r = 2*pi*V_r/("+str(wavelength_SI)+"*meter)") # frequency
				self.ureg.define("M_r = electron_mass"                    ) # mass
				self.ureg.define("Q_r = 1.602176565e-19 * coulomb"        ) # charge
			else:
				# Make blank unit registry
				self.ureg = self.UnitRegistry(None)
				self.ureg.define("V_r = [code_velocity]"                  ) # velocity
				self.ureg.define("W_r = [code_frequency]"                 ) # frequency
				self.ureg.define("M_r = [code_mass]"                      ) # mass
				self.ureg.define("Q_r = [code_charge]"                    ) # charge
				self.ureg.define("epsilon_0 = 1")
			self.ureg.define("L_r = V_r / W_r"                        ) # length
			self.ureg.define("T_r = 1   / W_r"                        ) # time
			self.ureg.define("P_r = M_r * V_r"                        ) # momentum
			self.ureg.define("K_r = M_r * V_r**2"                     ) # energy
			self.ureg.define("N_r = epsilon_0 * M_r * W_r**2 / Q_r**2") # density
			self.ureg.define("J_r = V_r * Q_r * N_r"                  ) # current
			self.ureg.define("B_r = M_r * W_r / Q_r "                 ) # magnetic field
			self.ureg.define("E_r = B_r * V_r"                        ) # electric field
			self.ureg.define("S_r = K_r * V_r * N_r"                  ) # poynting
			# Convert units if possible
			self.xcoeff, self.xname = self._convert(xunits, self.requestedX)
			self.ycoeff, self.yname = self._convert(yunits, self.requestedY)
			self.vcoeff, self.vname = self._convert(vunits, self.requestedV)
			self.tcoeff, self.tname = self._convert("T_r", None)
		else:
			self.xcoeff, self.xname = 1., "code units"
			self.ycoeff, self.yname = 1., "code units"
			self.vcoeff, self.vname = 1., "code units"
			self.tcoeff, self.tname = 1., "code units"



# -------------------------------------------------------------------
# Mother class for all diagnostics
# -------------------------------------------------------------------
class Diagnostic(object):
	
	# Initialize with "results_path" argument being either the `results_path` or
	# the parent `Smilei` object
	def __init__(self, results_path=None, *args, **kwargs):
		self.valid = False
		self._tmpdata = None
		self._animateOnAxes = None
		self._shape = []
		self._centers = []
		self._type = []
		self._label = []
		self._units = []
		self._log = []
		self._data_log = False
		
		# if results_path is string, try to use it as a results_path
		if type(results_path) is str:
			self.Smilei = Smilei(results_path)
		# if given a Smilei object, use it directly
		elif type(results_path) is Smilei:
			self.Smilei = results_path
		# if the Smilei object is outdated, reload
		elif hasattr(results_path,"_results_path") and hasattr(results_path,"reload") \
				and callable(getattr(results_path,"reload")):
			self.Smilei = results_path
			self.Smilei.reload()
		# Otherwise, error
		else:
			print "Could not find information on the Smilei simulation"
			return
		if not self.Smilei.valid:
			print "Invalid Smilei simulation"
			return
		
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
		
		# Make or retrieve the Units object
		self.units = kwargs.pop("units", [""])
		if type(self.units) in [list, tuple]: self.units = Units(*self.units)
		if type(self.units) is dict         : self.units = Units(**self.units)
		if type(self.units) is not Units:
			print "Could not understand the 'units' argument"
			return
		
		# Get some info on the simulation
		try:
			# get number of dimensions
			error = "Error extracting 'dim' from the input file"
			self._ndim = int(self.namelist.dim[0])
			if self._ndim not in [1,2,3]: raise
			# get box size
			error = "Error extracting 'sim_length' from the input file"
			sim_length = self._np.atleast_1d(self._np.double(self.namelist.sim_length))
			if sim_length.size != self._ndim: raise
			# get cell size
			error = "Error extracting 'cell_length' from the input file"
			self._cell_length = self._np.atleast_1d(self._np.double(self.namelist.cell_length))
			if self._cell_length.size != self._ndim: raise
			# calculate number of cells in each dimension
			self._ncels = sim_length/self._cell_length
			# extract time-step
			error = "Error extracting 'timestep' from the input file"
			self.timestep = self._np.double(self.namelist.timestep)
			if not self._np.isfinite(self.timestep): raise
		except:
			print error
			return None
		
		# Call the '_init' function of the child class
		self._init(*args, **kwargs)
		
		# Prepare units
		self._dim = len(self._shape)
		if self.valid:
			try:    wavelength_SI = self.namelist.wavelength_SI
			except: wavelength_SI = None
			xunits = None
			yunits = None
			if self._dim > 0: xunits = self._units[0]
			if self._dim > 1: yunits = self._units[1]
			self.units.prepare(wavelength_SI, xunits, yunits, self._vunits)
	
	# When no action is performed on the object, this is what appears
	def __repr__(self):
		if not self.valid: return ""
		self.info()
		return ""
	
	# Method to verify everything was ok during initialization
	def _validate(self):
		try:
			self.Smilei
		except:
			print "No valid Smilei simulation selected"
			return False
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
		for i in range(len(self._shape)):
			l.append([min(self._centers[i]), max(self._centers[i])])
		return l
	
	# Method to get only the arrays of data
	def getData(self):
		if not self._validate(): return
		self._prepare1() # prepare the vfactor
		data = []
		for t in self.times:
			data.append( self._vfactor*self._getDataAtTime(t) )
		return data
	
	# Method to obtain the data and the axes
	def get(self):
		if not self._validate(): return
		# obtain the data arrays
		data = self.getData()
		# format the results into a dictionary
		result = {"data":data, "times":self.times}
		for i in range(len(self._type)):
			result.update({ self._type[i]:self._centers[i] })
		return result
	
	# Method to plot the current diagnostic
	def plot(self, movie="", fps=15, dpi=200, saveAs=None, **kwargs):
		if not self._validate(): return
		self._prepare()
		self.set(**kwargs)
		self.info()
		
		# Make figure
		fig = self._plt.figure(**self.options.figure0)
		fig.set(**self.options.figure1)
		fig.clf()
		ax = fig.add_subplot(1,1,1)
		
		# Case of a streakPlot (no animation)
		if self.options.streakPlot:
			if len(self.times) < 2:
				print "ERROR: a streak plot requires at least 2 times"
				return
			if not hasattr(self,"_getDataAtTime"):
				print "ERROR: this diagnostic cannot do a streak plot"
				return
			if len(self._shape) != 1:
				print "ERROR: Diagnostic must be 1-D for a streak plot"
				return
			if not (self._np.diff(self.times)==self.times[1]-self.times[0]).all():
				print "WARNING: times are not evenly spaced. Time-scale not plotted"
				ylabel = "Unevenly-spaced times"
			else:
				ylabel = "Timesteps"
			# Loop times and accumulate data
			A = []
			for time in self.times: A.append(self._getDataAtTime(time))
			A = self._np.double(A)
			# Plot
			ax.cla()
			xmin = self._xfactor*self._centers[0][0]
			xmax = self._xfactor*self._centers[0][-1]
			extent = [xmin, xmax, self.times[0], self.times[-1]]
			if self._log[0]: extent[0:2] = [self._np.log10(xmin), self._np.log10(xmax)]
			im = ax.imshow(A, vmin = self.options.vmin, vmax = self.options.vmax, extent=extent, **self.options.image)
			ax.set_xlabel(self._xlabel)
			ax.set_ylabel(ylabel)
			self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.ymin, ymax=self.options.ymax)
			try: # if colorbar exists
				ax.cax.cla()
				self._plt.colorbar(mappable=im, cax=ax.cax, **self.options.colorbar)
			except AttributeError:
				ax.cax = self._plt.colorbar(mappable=im, ax=ax, **self.options.colorbar).ax
			self._setSomeOptions(ax)
			self._plt.draw()
			self._plt.pause(0.00001)
			return
		
		# Possible to skip animation
		if self.options.skipAnimation:
			self._animateOnAxes(ax, self.times[-1])
			self._plt.draw()
			self._plt.pause(0.00001)
			return
		
		# Otherwise, animation
		# Movie requested ?
		mov = Movie(fig, movie, fps, dpi)
		# Save to file requested ?
		save = SaveAs(saveAs, fig, self._plt)
		# Loop times for animation
		for time in self.times:
			print "timestep "+str(time)
			# plot
			ax.cla()
			if self._animateOnAxes(ax, time) is None: return
			self._plt.draw()
			self._plt.pause(0.00001)
			mov.grab_frame()
			save.frame(time)
		# Movie ?
		if mov.writer is not None: mov.finish()
	
	# Method to prepare some data before plotting
	def _prepare(self):
		self._prepare1()
		self._prepare2()
		self._prepare3()
		self._prepare4()
	
	# Methods to prepare stuff
	def _prepare1(self):
		# prepare the factors
		self._xfactor = (self.options.xfactor or 1.) * self.units.xcoeff
		self._yfactor = (self.options.yfactor or 1.) * self.units.ycoeff
		self._vfactor = self.units.vcoeff
		self._tfactor = (self.options.xfactor or 1.) * self.units.tcoeff * self.timestep
	def _prepare2(self):
		# prepare the animating function
		if not self._animateOnAxes:
			if   self._dim == 0: self._animateOnAxes = self._animateOnAxes_0D
			elif self._dim == 1: self._animateOnAxes = self._animateOnAxes_1D
			elif self._dim == 2: self._animateOnAxes = self._animateOnAxes_2D
			else:
				print "Cannot plot with more than 2 dimensions !"
				return
		# prepare t label
		self._tlabel = self.units.tname
		if self.options.xfactor: self._tlabel += "/"+str(self.options.xfactor)
		self._tlabel = 'Time ( '+self._tlabel+' )'
		# prepare x label
		if self._dim > 0:
			self._xlabel = self.units.xname
			if self.options.xfactor: self._xlabel += "/"+str(self.options.xfactor)
			self._xlabel = self._label[0] + " (" + self._xlabel + ")"
			if self._log[0]: self._xlabel = "Log[ "+self._xlabel+" ]"
		# prepare y label
		if self._dim > 1:
			self._ylabel = self.units.yname
			if self.options.yfactor: self._ylabel += "/"+str(self.options.yfactor)
			self._ylabel = self._label[1] + " (" + self._ylabel + ")"
			if self._log[1]: self._ylabel = "Log[ "+self._ylabel+" ]"
			# prepare extent for 2d plots
			self._extent = [self._xfactor*self._centers[0][0], self._xfactor*self._centers[0][-1], self._yfactor*self._centers[1][0], self._yfactor*self._centers[1][-1]]
			if self._log[0]:
				self._extent[0] = self._np.log10(self._extent[0])
				self._extent[1] = self._np.log10(self._extent[1])
			if self._log[1]:
				self._extent[2] = self._np.log10(self._extent[2])
				self._extent[3] = self._np.log10(self._extent[3])
		# prepare title
		self._vlabel = ""
		if self.units.vname: self._vlabel += " (" + self.units.vname + ")"
		if self._title     : self._vlabel = self._title + self._vlabel
		if self._data_log  : self._vlabel = "Log[ "+self._vlabel+" ]"
	def _prepare3(self):
		# prepare temporary data if zero-d plot
		if self._dim == 0 and self._tmpdata is None:
			self._tmpdata = self._np.zeros(self.times.size)
			for i, t in enumerate(self.times):
				self._tmpdata[i] = self._getDataAtTime(t)
	def _prepare4(self): pass
	
	# Method to set limits to a plot
	def _setLimits(self, ax, xmin=None, xmax=None, ymin=None, ymax=None):
		if xmin is not None: ax.set_xlim(xmin=xmin)
		if xmax is not None: ax.set_xlim(xmax=xmax)
		if ymin is not None: ax.set_ylim(ymin=ymin)
		if ymax is not None: ax.set_ylim(ymax=ymax)
	
	# Methods to plot the data when axes are made
	def _animateOnAxes_0D(self, ax, t):
		times = self.times[self.times<=t]
		A     = self._tmpdata[self.times<=t]
		im, = ax.plot(self._tfactor*times, self._vfactor*A, **self.options.plot)
		ax.set_xlabel(self._tlabel)
		self._setLimits(ax, xmax=self._tfactor*self.times[-1], ymin=self.options.vmin, ymax=self.options.vmax)
		self._setSomeOptions(ax)
		return im
	def _animateOnAxes_1D(self, ax, t):
		A = self._getDataAtTime(t)
		im, = ax.plot(self._xfactor*self._centers[0], self._vfactor*A, **self.options.plot)
		if self._log[0]: ax.set_xscale("log")
		ax.set_xlabel(self._xlabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.vmin, ymax=self.options.vmax)
		self._setSomeOptions(ax)
		return im
	def _animateOnAxes_2D(self, ax, t):
		A = self._getDataAtTime(t)
		im = self._animateOnAxes_2D_(ax, self._vfactor*A)
		ax.set_xlabel(self._xlabel)
		ax.set_ylabel(self._ylabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.ymin, ymax=self.options.ymax)
		try: # if colorbar exists
			ax.cax.cla()
			self._plt.colorbar(mappable=im, cax=ax.cax, **self.options.colorbar)
		except AttributeError:
			ax.cax = self._plt.colorbar(mappable=im, ax=ax, **self.options.colorbar).ax
		self._setSomeOptions(ax)
		return im
	
	# Special case: 2D plot
	# This is overloaded by class "Probe" because it requires to replace imshow
	def _animateOnAxes_2D_(self, ax, A):
		im = ax.imshow( self._np.flipud(A.transpose()),
			vmin = self.options.vmin, vmax = self.options.vmax, extent=self._extent, **self.options.image)
		return im
	
	# set options during animation
	def _setSomeOptions(self, ax):
		if self._vlabel: ax.set_title(self._vlabel)
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
	
	def dim(self):
		return len(self._shape)
	


# -------------------------------------------------------------------
# Class for particle diagnostics
# -------------------------------------------------------------------
class ParticleDiagnostic(Diagnostic):

	# This is the constructor, which creates the object
	def _init(self, diagNumber=None, timesteps=None, slice=None, data_log=False, **kwargs):
		
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
		
		cell_size = {"x":self._cell_length[0]}
		if self._ndim>1: cell_size.update({"y":self._cell_length[1]})
		if self._ndim>2: cell_size.update({"z":self._cell_length[2]})
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check the requested diags are ok
		if type(diagNumber) is int:
			if diagNumber<0:
				print "Argument 'diagNumber' cannot be a negative integer."
				return
			self.operation = '#' + str(diagNumber)
		elif type(diagNumber) is str:
			self.operation = diagNumber
		else:
			print"Argument 'diagNumber' must be and integer or a string."
			return
		
		# Get list of requested diags
		self._diags = sorted(set([ int(d[1:]) for d in self._re.findall('#\d+',self.operation) ]))
		for diag in self._diags:
			if not self._getInfo(diag):
				print "No particle diagnostic #"+str(diag)
				return
		try:
			exec(self._re.sub('#\d+','1.',self.operation)) in None
		except ZeroDivisionError: pass
		except:
			print "Cannot understand operation '"+self.operation+"'"
			return
		# Verify that all requested diags exist and they all have the same shape
		self._info = {}
		self._ishape = {}
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
			self._ishape.update({ d:[ axis["size"] for axis in self._axes[d] ] })
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
		self._ishape = self._ishape[self._diags[0]]
		
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
		plot_diff = []
		coeff = 1.
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
			axis_units = ""
			if   axis["type"] in ["x","y","z"]:
				axis_units = "L_r"
				spatialaxes[axis["type"]] = True
				coeff *= cell_size[axis["type"]]
			elif axis["type"] in ["px","py","pz","p"]:
				axis_units = "P_r"
			elif axis["type"] in ["vx","vy","vz","v"]:
				axis_units = "V_r"
			elif axis["type"] == "gamma":
				overall_min = "1"
			elif axis["type"] == "ekin":
				axis_units = "K_r"
				overall_min = "0"
			elif axis["type"] == "charge":
				axis_units = "Q_r"
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
				
				if axis["type"] in ["x","y","z"]: coeff /= slice_size
			
			# if not sliced, then add this axis to the overall plot
			else:
				self._type   .append(axis["type"])
				self._shape  .append(axis["size"])
				self._centers.append(centers)
				self._log    .append(axis["log"])
				self._label  .append(axis["type"])
				self._units  .append(axis_units)
				plot_diff.append(self._np.diff(edges))
		
		if len(self._shape) > 2:
			print "Cannot plot in "+str(len(self._shape))+"d. You need to 'slice' some axes."
			return None
		
		# Build units
		titles = {}
		units = {}
		for d in self._diags:
			titles.update({ d:"??" })
			units.update({ d:"??" })
			val_units = "??"
			output = self._info[d]["output"]
			if   output == "density":
				titles[d] = "Number density"
				val_units = "N_r"
			elif output == "charge_density":
				titles[d] = "Charge density"
				val_units = "N_r * Q_r"
			elif output[0] == "j":
				titles[d] = "J"+output[1]
				val_units = "J_r"
			elif output[0]=="p" and output[-8:] == "_density":
				titles[d] = "P"+output[1].strip("_")+" density"
				val_units = "N_r * P_r"
			elif output[:8]=="pressure":
				titles[d] = "Pressure "+output[-2]
				val_units = "N_r * K_r"
			axes_units = [unit for unit in self._units if unit!="L_r"]
			units[d] = val_units
			if len(axes_units)>0: units[d] += " / ( " + " * ".join(axes_units) + " )"
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for d in self._diags:
			self._vunits = self._vunits.replace("#"+str(d), "( "+units[d]+" )")
			self._title  = self._title .replace("#"+str(d), titles[d])
		
		# If any spatial dimension did not appear, then count it for calculating the correct density
		if self._ndim>=1 and not spatialaxes["x"]: coeff /= self._ncels[0]
		if self._ndim>=2 and not spatialaxes["y"]: coeff /= self._ncels[1]
		if self._ndim==3 and not spatialaxes["z"]: coeff /= self._ncels[2]
		
		# Calculate the array that represents the bins sizes in order to get units right.
		# This array will be the same size as the plotted array
		if len(plot_diff)==0:
			self._bsize = 1.
		elif len(plot_diff)==1:
			self._bsize = plot_diff[0]
		else:
			self._bsize = self._np.prod( self._np.array( self._np.meshgrid( *plot_diff ) ), axis=0)
			self._bsize = self._bsize.transpose()
		self._bsize /= coeff
		
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
			A.update({ d:self._np.reshape(self._h5items[d][index][1],self._ishape) })
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
	def _init(self, field=None, timesteps=None, slice=None, data_log=False, **kwargs):
		
		if not self.Smilei.valid: return None
		
		self._file = self._results_path+'/Fields.h5'
		try:
			self._f = self._h5py.File(self._file, 'r')
		except:
			print "Cannot open file "+self._file
			return
		self._h5items = self._f.values()

		if field is None:
			fields = self.getFields()
			if len(fields)>0:
				print "Printing available fields:"
				print "--------------------------"
				l = (len(fields)/3) * 3
				maxlength = str(self._np.max([len(f) for f in fields])+4)
				fields = [('%'+maxlength+'s')%f for f in fields]
				if l>0:
					print '\n'.join([''.join(list(i)) for i in self._np.reshape(fields[:l],(-1,3))])
				print ''.join(list(fields[l:]))
			else:
				print "No fields found in '"+self._results_path+"'"
			return None
		

		
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
		self._operation = "A="+self.operation
		self._fieldname = []
		for f in sortedfields:
			if f in self._operation:
				self._operation = self._re.sub(r"\b"+f+r"\b","C['"+f+"']",self._operation)
				self._fieldname.append(f)
		
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			print "Argument `slice` must be a dictionary"
			return None
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Get the shape of fields
		iterfields = self._h5items[0].itervalues();
		self._ishape = iterfields.next().shape;
		for fd in iterfields:
			self._ishape = self._np.min((self._ishape, fd.shape), axis=0)
		
		
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
		self._naxes = self._ndim
		self._sliceinfo = {}
		self._slices = [None]*self._ndim
		for iaxis in range(self._naxes):
			centers = self._np.linspace(0., self._ishape[iaxis]*self._cell_length[iaxis], self._ishape[iaxis])
			label = {0:"x", 1:"y", 2:"z"}[iaxis]
			axisunits = "L_r"
			
			if label in slice:
				# if slice is "all", then all the axis has to be summed
				if slice[label] == "all":
					indices = self._np.arange(self._ishape[iaxis])
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
				self._slices[iaxis] = self._np.delete(self._np.arange(self._ishape[iaxis]), indices)
			else:
				self._type   .append(label)
				self._shape  .append(self._ishape[iaxis])
				self._centers.append(centers)
				self._label  .append(label)
				self._units  .append(axisunits)
				self._log    .append(False)
		
		if len(self._centers) > 2:
			print "Cannot plot in "+str(len(self._shape))+"d. You need to 'slice' some axes."
			return
		
		# Build units
		units = {}
		for f in self._fieldname:
			units.update({ f:{"B":"B_r", "E":"E_r", "J":"J_r", "R":"N_r"}[f[0]] })
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for f in self._fieldname:
			self._vunits = self._vunits.replace(f, units[f])
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included fields
	def info(self):
		if not self._validate(): return
		print self._title
		return
	
	# get all available fields, sorted by name length
	def getFields(self):
		try:    fields = self._f.itervalues().next().keys() # list of fields
		except: fields = []
		return fields
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		times = self._np.double(self._f.keys())
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
		h5item = self._h5items[index]
		for field in self._fieldname: # for each field in operation
			B = self._np.double(h5item.get(field)) # get array
			for axis, size in enumerate(self._ishape):
				l = self._np.arange(size, B.shape[axis])
				B = self._np.delete(B, l, axis=axis) # remove extra cells if necessary
			C.update({ field:B })
		# Calculate the operation
		exec self._operation in None
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
	def _init(self, scalar=None, timesteps=None, data_log=False, **kwargs):
		
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
		self._slices = []
		# Build units
		self._vunits = "??"
		if   self._scalarname == "time":
			self._vunits = "T_r"
		elif self._scalarname == "Ubal_norm" or self._scalarname[0] in ["N","Z"]:
			self._vunits = ""
		else:
			self._vunits = {"U":"K_r", "E":"E_r", "B":"B_r", "J":"J_r",
									"R":"N_r", "P":"S_r"}[self._scalarname[0]]
		self._title =self._scalarname
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included scalars
	def info(self):
		if not self._validate(): return
		print "Scalar "+self._scalarname
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
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A




# -------------------------------------------------------------------
# Class for probe diagnostics
# -------------------------------------------------------------------
class Probe(Diagnostic):

	# This is the constructor, which creates the object
	def _init(self, probeNumber=None, field=None, timesteps=None, slice=None, data_log=False, **kwargs):
		
		if not self.Smilei.valid: return None
		
		# If no probeNumber, print available probes
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
		
		# Try to get the probe from the hdf5 file
		self.probeNumber  = probeNumber
		self._file = self._results_path+"/Probes.h5"
		f = self._h5py.File(self._file, 'r')
		self._h5probe = None
		for key in f.keys():
			if key[0] != "p": continue
			if int(key.strip("p"))==probeNumber:
				self._h5probe = f.get(key)
				break
		if self._h5probe is None:
			print "Cannot find probe "+str(probeNumber)+" in file "+file
			f.close()
			return None
		
		# Extract available fields
		fields = self._h5probe.attrs["fields"].split(",")
		if len(fields) == 0:
			print "Probe #"+probeNumber+" is empty"
			f.close()
			return None
		# If no field, print available fields
		if field is None:
			print "Printing available fields for probe #"+str(probeNumber)+":"
			print "----------------------------------------"
			print ", ".join(fields)
			f.close()
			return None
		
		# Get available times
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
			print "No probes found in Probes.h5"
			return
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Parse the `field` argument
		sortedfields = reversed(sorted(fields, key = len))
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
		self._ishape = self._info["shape"]
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with indices to the data arrays
		self._data = {}
		for t in self.times:
			self._data.update({ t : "%010i"%t })
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
		self._naxes = self._ishape.size
		self._sliceinfo = {}
		self._slices = [None]*self._ndim
		for iaxis in range(self._naxes):
		
			# calculate grid points locations
			p0 = self._info["p0"            ] # reference point
			pi = self._info["p"+str(iaxis+1)] # end point of this axis
			centers = self._np.zeros((self._ishape[iaxis],p0.size))
			for i in range(p0.size):
				centers[:,i] = self._np.linspace(p0[i],pi[i],self._ishape[iaxis])
			
			label = {0:"axis1", 1:"axis2", 2:"axis3"}[iaxis]
			axisunits = "L_r"
			
			if label in slice:
				# if slice is "all", then all the axis has to be summed
				if slice[label] == "all":
					indices = self._np.arange(self._ishape[iaxis])
				# Otherwise, get the slice from the argument `slice`
				else:
					indices = self._np.arange(self._ishape[iaxis])
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
				self._slices[iaxis] = self._np.delete(self._np.arange(self._ishape[iaxis]), indices)
			else:
				self._type   .append(label)
				self._shape  .append(self._ishape[iaxis])
				self._centers.append(centers)
				self._label  .append(label)
				self._units  .append(axisunits)
				self._log    .append(False)
			
		
		if len(self._shape) > 2:
			print "Cannot plot in "+str(len(self._shape))+"d. You need to 'slice' some axes."
			return
		
		# Special case in 1D: we convert the point locations to scalar distances
		if len(self._centers) == 1:
			self._centers[0] = self._np.sqrt(self._np.sum((self._centers[0]-self._centers[0][0])**2,axis=1))
		# Special case in 2D: we have to prepare for pcolormesh instead of imshow
		elif len(self._centers) == 2:
			p1 = self._centers[0] # locations of grid points along first dimension
			d = self._np.diff(p1, axis=0) # separation between the points
			p1 = self._np.vstack((p1, p1[-1,:])) # add last edges at the end of box
			p1[1:-1] -= d/2 # move points by one half
			p2 = self._centers[1] # locations of grid points along second dimension
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
			self._edges = [X, Y]
			self._label = ["x", "y"]
			self._units = [axisunits, axisunits]
		
		# Build units
		titles = {}
		fieldunits = {}
		for f in self._fieldname:
			i = fields.index(f)
			fieldunits.update({ i:{"B":"B_r","E":"E_r","J":"J_r","R":"N_r"}[f[0]] })
			titles    .update({ i:f })
		# Make total units and title
		self._title  = self.operation
		self._vunits = self.operation
		for n in self._fieldn:
			self._title  = self._title .replace("#"+str(n), titles    [n])
			self._vunits = self._vunits.replace("#"+str(n), fieldunits[n])
		
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
			+" every "+str(info["every"])+" timesteps, with fields "+info["fields"])
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
		for key in f.iterkeys():
			if key[0] != "p": continue
			if int(key.strip("p"))==probeNumber:
				probe = f[key]
		if probe is None:
			print "Cannot find probe "+str(probeNumber)+" in file "+file
			return {}
		out = {}
		out.update({"probeNumber":probeNumber, "dimension":probe.attrs["dimension"],
			"every":probe.attrs["every"], "shape":self._np.array(probe["number"]),
			"fields":probe.attrs["fields"] })
		i = 0
		while "p"+str(i) in probe.keys():
			k = probe.keys().index("p"+str(i))
			out.update({ "p"+str(i):self._np.array(probe.values()[k]) })
			i += 1
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
			probes = [int(key.strip("p")) for key in f.iterkeys()] # list of probes
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
		times = []
		for key in self._h5probe.iterkeys():
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
		# get data
		index = self._data[t]
		C = {}
		op = "A=" + self.operation
		for n in reversed(self._fieldn): # for each field in operation
			B = self._np.double(self._h5probe[index][n,:]) # get array
			B = self._np.reshape(B, self._ishape) # reshape array because it is flattened in the file
			C.update({ n:B })
			op = op.replace("#"+str(n), "C["+str(n)+"]")
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
	
	# We override _prepare4
	def _prepare4(self):
		# If 2D plot, we remove kwargs that are not supported by pcolormesh
		if self._dim == 2:
			authorizedKwargs = ["cmap"]
			for kwarg in self.options.image.keys():
				if kwarg not in authorizedKwargs: del self.options.image[kwarg]
	
	# Overloading a plotting function in order to use pcolormesh instead of imshow
	def _animateOnAxes_2D_(self, ax, A):
		im = ax.pcolormesh(self._xfactor*self._edges[0], self._yfactor*self._edges[1], self._np.flipud(A.transpose()),
			vmin = self.options.vmin, vmax = self.options.vmax, **self.options.image)
		return im




# -------------------------------------------------------------------
# Class for tracked particles diagnostics
# -------------------------------------------------------------------
class TrackParticles(Diagnostic):

	# This is the constructor, which creates the object
	def _init(self, species=None, select="", axes=[], timesteps=None, **kwargs):
		
		if not self.Smilei.valid: return None
		
		# If argument 'species' not provided, then print available species and leave
		if species is None:
			species = self.getTrackSpecies()
			if len(species)>0:
				print "Printing available tracked species:"
				print "-----------------------------------"
				for s in species: print s
			else:
				print "No tracked particles files found in '"+self._results_path+"'"
			return None
		
		# Get info from the hdf5 files + verifications
		# -------------------------------------------------------------------
		self.species  = species
		self._file = self._results_path+"/TrackParticles_"+species+".h5"
		f = self._h5py.File(self._file, 'r')
		self._h5items = f.values()
		self._every = f.attrs["every"]
		
		# Get available times in the hdf5 file
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
			print "No tracked particles found in "+self._file
			return
		alltimes = self.times
		# If specific timesteps requested, narrow the selection
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
		
		# Get available properties ("x", "y", etc.)
		self._properties = {}
		translateProperties = {"Id":"Id", "x":"Position-0", "y":"Position-1", "z":"Position-2",
							"px":"Momentum-0", "py":"Momentum-1", "pz":"Momentum-2"}
		availableProperties = f.keys()
		for k,v in translateProperties.iteritems():
			try:
				i = availableProperties.index(v)
				self._properties.update({ k:i })
				if k == "Id": self._Id = self._h5items[i]
			except:
				pass
		
		# Get number of particles
		self.nParticles = self._h5items[0].shape[1]
		
		# Select particles
		# -------------------------------------------------------------------
		if type(select) is not str:
			print "Error: the argument 'select' must be a string"
			return
		def findClosingCharacter(string, character, start=0):
			i = start
			stack = []
			associatedBracket = {")":"(", "]":"[", "}":"{"}
			while i < len(string):
				if string[i] == character and len(stack)==0: return i
				if string[i] in ["(", "[", "{"]:
					stack.append(string[i])
				if string[i] in [")", "]", "}"]:
					if len(stack)==0:
						raise Exception("Error in selector syntax: missing `"+character+"`")
					if stack[-1]!=associatedBracket[string[i]]:
						raise Exception("Error in selector syntax: missing closing parentheses or brackets")
					del stack[-1]
				i+=1
			raise Exception("Error in selector syntax: missing `"+character+"`")
		i = 0
		stack = []
		operation = ""
		while i < len(select):
			if i+4<len(select):
				if select[i:i+4] == "any(" or select[i:i+4] == "all(":
					if select[i:i+4] == "any(": function = self._np.logical_or
					if select[i:i+4] == "all(": function = self._np.logical_and
					comma = findClosingCharacter(select, ",", i+4)
					parenthesis = findClosingCharacter(select, ")", comma+1)
					timeSelector = select[i+4:comma]
					try:
						s = self._re.sub(r"\bt\b","alltimes",timeSelector)
						times = alltimes[eval(s)]
					except:
						raise Exception("Error in selector syntax: time selector not understood in "+select[i:i+3]+"()")
					try:
						particleSelector = select[comma+1:parenthesis]
						for prop in self._properties.keys():
							particleSelector = self._re.sub(r"\b"+prop+r"\b", "self._np.double(self._h5items["+str(self._properties[prop])+"][t,:])", particleSelector)
					except:
						raise Exception("Error in selector syntax: not understood: "+select[i:parenthesis+1])
					if select[i:i+4] == "any(": selection = self._np.array([False]*self.nParticles)
					if select[i:i+4] == "all(": selection = self._np.array([True]*self.nParticles)
					#try:
					ID = self._np.zeros((self.nParticles,), dtype=self._np.int16)
					for t in times:
						selectionAtTimeT = eval(particleSelector) # array of True or False
						self._Id.read_direct(ID, source_sel=self._np.s_[t,:], dest_sel=self._np.s_[:]) # read the particle Ids
						selectionAtTimeT = selectionAtTimeT[ID>0] # remove zeros, which are dead particles
						id = ID[ID>0]-1 # remove zeros, which are dead particles
						selection[id] = function( selection[id], selectionAtTimeT)
					#except:
					#	raise Exception("Error in selector syntax: not understood: "+select[i:parenthesis+1])
					stack.append(selection)
					operation += "stack["+str(len(stack)-1)+"]"
					i = parenthesis+1
					continue
			operation += select[i]
			i+=1
		if len(operation)==0.:
			self.selectedParticles = self._np.arange(self.nParticles)
		else:
			self.selectedParticles = eval(operation).nonzero()[0]
		self.selectedParticles.sort()
		self.selectedParticles += 1
		self.nselectedParticles = len(self.selectedParticles)
		
		# Manage axes
		# -------------------------------------------------------------------
		if type(axes) is not list:
			print "Error: Argument 'axes' must be a list"
			return
		if len(axes)==0:
			print "Error: must define at least one axis."
			return
		self.axes = axes
		self._axesIndex = []
		for axis in axes:
			if axis not in self._properties.keys():
				print "Error: Argument 'axes' has item '"+str(axis)+"' unknown."
				print "       Available axes are: "+(", ".join(sorted(self._properties.keys())))
				return
			self._axesIndex.append( self._properties[axis] ) # axesIndex contains the index in the hdf5 file
		self._type = axes
		for i, axis in enumerate(axes):
			axisi = self._axesIndex[i]
			vals = self._np.double(self._h5items[axisi])
			axisunits = ""
			if axis != "Id":
				if axis in ["x" , "y" , "z" ]: axisunits = "L_r"
				if axis in ["px", "py", "pz"]: axisunits = "P_r"
			self._centers.append( [vals.min(), vals.max()] )
			self._log.append( False )
			self._label.append( axis )
			self._units.append( axisunits )
		self._title = "Track particles '"+species+"'"
		self._shape = [0]*len(axes)
		# Hack to work with 1 axis
		if len(axes)==1: self._vunits = self._units[0]
		else: self._vunits = ""
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included probe
	def info(self):
		if not self._validate(): return
		print "Track particles: species '"+self.species+"' containing "+str(self.nParticles)+" particles"
		if len(self.selectedParticles) != self.nParticles:
			print "                with selection of "+str(len(self.selectedParticles))+" particles"
	
	# get all available tracked species
	def getTrackSpecies(self):
		files = self._glob(self._results_path+"/TrackParticles_*.h5")
		species = []
		for file in files:
			species.append(self._re.search("TrackParticles_(.*).h5",file).groups()[0])
		return species
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:
			ntimes = self._h5items[0].len()
		except:
			print "Unable to find tracked particle data in file "+self._file
			return self._np.array([])
		return self._np.arange(0,ntimes*self._every,self._every)
	
	# We override the get and getData methods
	def getData(self):
		if not self._validate(): return
		self._prepare1() # prepare the vfactor
		# create dictionary with info on the axes
		data = {}
		for axis in self.axes:
			data.update({ axis:self._np.zeros((len(self.times), self.nselectedParticles)) })
			data[axis].fill(self._np.nan)
		# loop times and fill up the data
		ID = self._np.zeros((self.nParticles,), dtype=self._np.int16)
		B = self._np.zeros((self.nParticles,))
		indices = self.selectedParticles - 1
		for ti in range(len(self.times)):
			self._Id.read_direct(ID, source_sel=self._np.s_[ti,:], dest_sel=self._np.s_[:]) # read the particle Ids
			deadParticles = (ID==0).nonzero()
			for i, axis in enumerate(self.axes):
				axisi = self._axesIndex[i]
				self._h5items[axisi].read_direct(B, source_sel=self._np.s_[ti,:], dest_sel=self._np.s_[:])
				B[deadParticles]=self._np.nan
				data[axis][ti, :] = B[indices].squeeze() * self._vfactor
		data.update({ "times":self.times })
		return data
	def get(self):
		return self.getData()
	
	# We override _prepare3
	def _prepare3(self):
		if self._tmpdata is None:
			A = self.getData()
			self._tmpdata = []
			for axis in self.axes: self._tmpdata.append( A[axis] )
	
	# We override the plotting methods
	def _animateOnAxes_0D(self, ax, t):
		pass
	def _animateOnAxes_1D(self, ax, t):
		times = self.times[self.times<=t]
		A     = self._tmpdata[0][self.times<=t,:]
		if times.size == 1:
			times = self._np.double([times, times]).squeeze()
			A = self._np.double([A, A]).squeeze()
		ax.plot(self._tfactor*times, self._vfactor*A, **self.options.plot)
		ax.set_xlabel(self._tlabel)
		ax.set_ylabel(self.axes[0]+" ("+self.units.vname+")")
		self._setLimits(ax, xmax=self._tfactor*self.times[-1], ymin=self.options.vmin, ymax=self.options.vmax)
		self._setSomeOptions(ax)
		ax.set_title(self._title) # override title
		return 1
	def _animateOnAxes_2D(self, ax, t):
		x = self._tmpdata[0][self.times<=t,:]
		y = self._tmpdata[1][self.times<=t,:]
		ax.plot(self._xfactor*x, self._yfactor*y, **self.options.plot)
		ax.set_xlabel(self._xlabel)
		ax.set_ylabel(self._ylabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.ymin, ymax=self.options.ymax)
		self._setSomeOptions(ax)
		return 1


class Movie:
	
	def __init__(self, fig, movie="", fps=15, dpi=200):
		import os.path as ospath
		self.writer = None
		if type(movie) is not str:
			print "ERROR: argument 'movie' must be a filename"
			return
		if len(movie)>0:
			if ospath.isfile(movie):
				print "ERROR: file '"+movie+"' already exists. Please choose another name"
				return
			if ospath.isdir(movie):
				print "ERROR: '"+movie+"' is a path, not a file"
				return
			try:
				import matplotlib.animation as anim
			except:
				print "ERROR: it looks like your version of matplotlib is too old for movies"
				return
			try:
				self.writer = anim.writers['ffmpeg'](fps=fps)
			except:
				print "ERROR: you need the 'ffmpeg' software to make movies"
				return
			self.writer.setup(fig, movie, dpi)
	
	def finish(self):
		if self.writer is not None:
			self.writer.finish()
	
	def grab_frame(self):
		if self.writer is not None:
			self.writer.grab_frame()



class SaveAs:
	
	def __init__(self, smartPath, fig, plt):
		import os.path as p
		from os import sep as sep
		default_extension = ".png"
		self.figure = fig
		
		self.prefix = False
		if type(smartPath) is str:
			path = smartPath
			if len(smartPath)==0 or smartPath[0].isalnum():
				path = "."+sep+path
			if p.isdir(path):
				self.prefix = p.normpath(path)+sep
				self.suffix = default_extension
			else:
				path, base = p.split(path)
				if p.isdir(path):
					basesplit = base.rsplit(".",1)
					self.prefix = path+sep+basesplit[0]
					self.suffix = ""
					if len(basesplit)>1: self.suffix = "."+basesplit[1]
				else:
					self.prefix = False
					self.suffix = "`"+path+"` is not a directory"
		else:
			return
		
		if not self.prefix:
			print "WARNING: "+self.suffix
			print "         Will not save figures to files"
		
		# For some reason, the following freezes the animation; removed temporarily
		#else:
		#	supportedTypes=plt.matplotlib.backend_bases.FigureCanvasBase(fig).get_supported_filetypes()
		#	if self.suffix.strip(".") not in supportedTypes.iterkeys():
		#		print "WARNING: file format not supported, will not save figures to files"
		#		print "         Supported formats: "+",".join(supportedTypes.iterkeys())
		#		self.prefix = False
	
	def frame(self, id):
		if self.prefix:
			self.figure.savefig(self.prefix + "%010d"%int(id) + self.suffix)


def multiPlot(*Diags, **kwargs):
	""" multiplot(Diag1, Diag2, ..., shape=None, movie="", fps=15, dpi=200, saveAs=None)
	
	Plots simultaneously several diagnostics.
	
	Parameters:
	-----------
	Diag1, Diag2, ... : Several objects of classes 'Scalar', 'Field', 'Probe' or 'ParticleDiagnostic'
	shape : 2D list giving the number of figures in x and y.
	movie : filename to create a movie.
	fps : frames per second for the movie.
	dpi : resolution of the movie.
	saveAs : path where to store individual frames as pictures.
	skipAnimation : toggle going directly to the last frame.
	"""
	# Verify Diags are valid
	nDiags = len(Diags)
	if nDiags == 0: return
	for Diag in Diags:
		if not Diag.valid: return
	np  = Diags[0]._np  # numpy
	plt = Diags[0]._plt # pyplot
	# Get keyword arguments
	shape  = kwargs.pop("shape" , None)
	movie  = kwargs.pop("movie" , ""  )
	fps    = kwargs.pop("fps"   , 15  )
	dpi    = kwargs.pop("dpi"   , 200 )
	saveAs = kwargs.pop("saveAs", None)
	skipAnimation = kwargs.pop("skipAnimation", False )
	# Gather all times
	if skipAnimation:
		alltimes = np.unique([Diag.times[-1]*Diag.timestep for Diag in Diags])
	else:
		alltimes = np.unique(np.concatenate([Diag.times*Diag.timestep for Diag in Diags]))
	# Determine whether to plot all cases on the same axes
	sameAxes = False
	if shape is None or shape == [1,1]:
		sameAxes = True
		for Diag in Diags:
			if type(Diag) is TrackParticles:
				sameAxes = False
				break
			if Diag.dim()==0 and Diags[0].dim()==0:
				continue
			if Diag.dim()!=1 or Diag._type != Diags[0]._type:
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
		Diag._prepare()
	# Static plot
	if sameAxes and len(Diags[0]._shape)==0:
		for Diag in Diags:
			Diag._artist = Diag._animateOnAxes(Diag._ax, Diag.times[-1])
			plt.draw()
			plt.pause(0.00001)
	# Animated plot
	else:
		# Loop all times
		mov = Movie(fig, movie, fps, dpi)
		save = SaveAs(saveAs, fig, plt)
		for i,time in enumerate(alltimes):
			t = None
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
			plt.draw()
			plt.pause(0.00001)
			mov.grab_frame()
			if t is not None: save.frame(int(t))
		mov.finish()
		return


