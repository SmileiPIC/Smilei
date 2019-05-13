from .._Utils import *

class Diagnostic(object):
	"""Mother class for all Diagnostics.
	To create a diagnostic, refer to the doc of the SmileiSimulation class.
	Once such object is created, you may get help on its diagnostics.

	Example:
		S = happi.Open("path/to/simulation") # Load a simulation
		help( S )                            # General help on the simulation's diagnostics
		help( S.Field )                      # Help on loading a Field diagnostic
	"""

	def __init__(self, simulation, *args, **kwargs):
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
		self._error = []
		self._xoffset = 0.
		
		# The 'simulation' is a SmileiSimulation object. It is passed as an instance attribute
		self.simulation = simulation

		# Transfer the simulation's packages to the diagnostic
		self._h5py    = self.simulation._h5py
		self._np      = self.simulation._np
		self._os      = self.simulation._os
		self._glob    = self.simulation._glob
		self._re      = self.simulation._re
		self._plt     = self.simulation._plt
		self._verbose = self.simulation._verbose
		
		# Reload the simulation, in case it has been updated
		self.simulation.reload()
		if not self.simulation.valid:
			self._error += ["Invalid Smilei simulation"]
			return
		
		# Copy some parameters from the simulation
		self._results_path = self.simulation._results_path
		self.namelist      = self.simulation.namelist
		self._ndim         = self.simulation._ndim
		self._cell_length  = self.simulation._cell_length
		self._ncels        = self.simulation._ncels
		self.timestep      = self.simulation._timestep
		
		# Make the Options object
		self.options = Options()
		kwargs = self.options.set(**kwargs)
		
		# Make or retrieve the Units object
		self.units = kwargs.pop("units", [""])
		if type(self.units) in [list, tuple]: self.units = Units(*self.units , verbose = self._verbose)
		if type(self.units) is dict         : self.units = Units(verbose = self._verbose, **self.units)
		if type(self.units) is not Units:
			self._error += ["Could not understand the 'units' argument"]
			return
		self.units.prepare(self.simulation._reference_angular_frequency_SI)
		
		# DEPRECATION ERRORS
		if "slice" in kwargs:
			self._error += ["Argument `slice` is deprecated: use `subset` instead."]
			return
		if "stride" in kwargs:
			self._error += ["Argument `stride` is deprecated: use `average` or `sum` instead."]
			return
		
		# Call the '_init' function of the child class
		remaining_kwargs = self._init(*args, **kwargs)
		if remaining_kwargs is not None and len(remaining_kwargs) > 0:
			self.valid = False
			self._error += ["The following keyword-arguments are unknown: "+", ".join(remaining_kwargs.keys())]
			return
		
		# Prepare units for axes
		self.dim = len(self._shape)
		if self.valid:
			xunits = None
			yunits = None
			if self.dim > 0: xunits = self._units[0]
			if self.dim > 1: yunits = self._units[1]
			self.units.convertAxes(xunits, yunits, self._vunits)
		
		# Prepare data_log output
		if self._data_log:
			self._dataAtTime = self._dataLogAtTime
		
	# When no action is performed on the object, this is what appears
	def __repr__(self):
		self.info()
		return ""

	# Method to verify everything was ok during initialization
	def _validate(self):
		try:
			self.simulation.valid
		except:
			print("No valid Smilei simulation selected")
			return False
		if not self.simulation.valid or not self.valid:
			print("Diagnostic is invalid")
			print("\n".join(self._error))
			return False
		return True

	# Method to set optional plotting arguments
	def set(self, **kwargs):
		self.options.set(**kwargs)
		return self

	# Method to set optional plotting arguments, but also checks all are known
	def _setAndCheck(self, **kwargs):
		kwargs = self.options.set(**kwargs)
		if len(kwargs)>0:
			unknown_kwargs = ", ".join(kwargs.keys())
			print("Error: The following arguments are unknown: "+unknown_kwargs)
			return False
		else:
			return True

	# Method to obtain the plot limits
	def limits(self):
		"""Gets the overall limits of the diagnostic along its axes

		Returns:
		--------
		A list of [min, max] for each axis.
		"""
		l = []
		factor = [self._xfactor, self._yfactor]
		for i in range(self.dim):
			l.append([min(self._centers[i])*factor[i], max(self._centers[i])*factor[i]])
		return l

	# Method to print info on this diag
	def info(self):
		if self._validate() and self._verbose:
			print(self._info())

	# Method to get only the arrays of data
	def getData(self, timestep=None):
		"""Obtains the data from the diagnostic.

		Parameters:
		-----------
		timestep: int (default: None, which means all available timesteps)

		Returns:
		--------
		A list of arrays: each array corresponding to the diagnostic data at a given
		timestep.
		"""
		if not self._validate(): return
		self._prepare1() # prepare the vfactor
		data = []

		if timestep is None:
			for t in self._timesteps:
				data.append( self._dataAtTime(t) )
		elif timestep not in self._timesteps:
			print("ERROR: timestep "+str(timestep)+" not available")
		else:
			data.append( self._dataAtTime(timestep) )

		return data

	def getTimesteps(self):
		"""Obtains the list of timesteps selected in this diagnostic"""
		if not self._validate(): return []
		return self._timesteps

	def getTimes(self):
		"""
		Obtains the list of times selected in this diagnostic.
		By default, times are in the code's units, but are converted to the diagnostic's
		units defined by the `units` argument, if provided.
		"""
		if not self._validate(): return []
		return self.units.tcoeff * self.timestep * self._np.array(self._timesteps)

	def getAxis(self, axis):
		"""
		Obtains the list of positions of the diagnostic data along the requested axis.
		By default, axis positions are in the code's units, but are converted to
		the diagnostic's units defined by the `units` argument, if provided.

		Parameters:
		-----------
		axis: str
			The name of the requested axis.

		Returns:
		--------
		A list of positions along the requested axis.
		(If the requested axis is not available, returns an empty list.)

		Example: if `x` is an available axis, `Diag.getAxis("x")` returns a list
		of the positions of the diagnostic data along x.
		"""
		try: axis_index = self._type.index(axis)
		except: return []
		factor = 1.
		if axis_index == 1: factor = (self.options.xfactor or 1.) * self.units.xcoeff
		if axis_index == 2: factor = (self.options.yfactor or 1.) * self.units.ycoeff
		return factor * self._np.array(self._centers[axis_index])

	# Method to obtain the data and the axes
	def get(self, timestep=None):
		"""Obtains the data from the diagnostic and some additional information.

		!!! Deprecated !!!
		Use functions `getData`, `getTimesteps`, `getTimes` and `getAxis` instead.
		"""
		if not self._validate(): return
		# obtain the data arrays
		data = self.getData(timestep=timestep)
		# format the results into a dictionary
		result = {"data":data, "times":self._timesteps}
		for i in range(len(self._type)):
			result.update({ self._type[i]:self._centers[i] })
		return result

	def _make_axes(self, axes):
		if axes is None:
			fig = self._plt.figure(**self.options.figure0)
			fig.set(**self.options.figure1)
			fig.clf()
			ax = fig.add_subplot(1,1,1)
			if self.options.side == "right":
				ax.yaxis.tick_right()
				ax.yaxis.set_label_position("right")
			return ax
		else:
			return axes

	def plot(self, timestep=None, saveAs=None, axes=None, **kwargs):
		""" Plots the diagnostic.

		Parameters:
		-----------
		timestep: int (default: None, which means the last timestep)
			The number of the timestep to plot
		figure: int (default: 1)
			The figure number that is passed to matplotlib.
		axes: (default: None)
			Matplotlib's axes handle on which to plot. If None, make new axes.
		vmin, vmax: floats (default to the current limits)
			Data value limits.
		xmin, xmax, ymin, ymax: floats (default to the current limits)
			Axes limits.
		xfactor, yfactor: floats (default: 1)
			Factors to rescale axes.
		saveAs: path string (default: None)
			Name of a directory where to save each frame as figures.
			You can even specify a filename such as mydir/prefix.png
			and it will automatically make successive files showing
			the timestep: mydir/prefix0.png, mydir/prefix1.png, etc.

		Example:
		--------
			S = happi.Open("path/to/my/results")
			S.ParticleBinning(1).plot(vmin=0, vmax=1e14)
		"""
		if not self._validate(): return
		if not self._prepare(): return
		if not self._setAndCheck(**kwargs): return
		self.info()
		ax = self._make_axes(axes)
		fig = ax.figure

		if timestep is None:
			timestep = self._timesteps[-1]
		elif timestep not in self._timesteps:
			print("ERROR: timestep "+str(timestep)+" not available")
			return

		save = SaveAs(saveAs, fig, self._plt)
		self._animateOnAxes(ax, timestep)
		self._plt.draw()
		self._plt.pause(0.00001)
		save.frame()
		return

	def streak(self, saveAs=None, axes=None, **kwargs):
		""" Plots the diagnostic with one axis being time.

		Parameters:
		-----------
		figure: int (default: 1)
			The figure number that is passed to matplotlib.
		axes: (default: None)
			Matplotlib's axes handle on which to plot. If None, make new axes.
		vmin, vmax: floats (default to the current limits)
			Data value limits.
		xmin, xmax, ymin, ymax: floats (default to the current limits)
			Axes limits.
		xfactor, yfactor: floats (default: 1)
			Factors to rescale axes.
		saveAs: path string (default: None)
			Name of a directory where to save each frame as figures.
			You can even specify a filename such as mydir/prefix.png
			and it will automatically make successive files showing
			the timestep: mydir/prefix0.png, mydir/prefix1.png, etc.

		Example:
		--------
			S = happi.Open("path/to/my/results")
			S.ParticleBinning(1).streak(vmin=0, vmax=1e14)
		"""
		if not self._validate(): return
		if not self._prepare(): return
		if not self._setAndCheck(**kwargs): return
		self.info()
		ax = self._make_axes(axes)
		fig = ax.figure

		if len(self._timesteps) < 2:
			print("ERROR: a streak plot requires at least 2 times")
			return
		if not hasattr(self,"_getDataAtTime"):
			print("ERROR: this diagnostic cannot do a streak plot")
			return
		if self.dim != 1:
			print("ERROR: Diagnostic must be 1-D for a streak plot")
			return
		if not (self._np.diff(self._timesteps)==self._timesteps[1]-self._timesteps[0]).all():
			print("WARNING: times are not evenly spaced. Time-scale not plotted")
			ylabel = "Unevenly-spaced times"
		else:
			self._yfactor = self.options.yfactor or 1.
			ylabel = "Timesteps"
			if self._yfactor != 1.:
				ylabel += " x "+str(self._yfactor)
		# Loop times and accumulate data
		A = self._np.double([self._dataAtTime(t) for t in self._timesteps])
		# Plot
		ax.cla()
		xmin = self._xfactor*self._centers[0][0]
		xmax = self._xfactor*self._centers[0][-1]
		extent = [xmin, xmax, self._yfactor*self._timesteps[0], self._yfactor*self._timesteps[-1]]
		if self._log[0]: extent[0:2] = [self._np.log10(xmin), self._np.log10(xmax)]
		im = ax.imshow(self._np.flipud(A), vmin = self.options.vmin, vmax = self.options.vmax, extent=extent, **self.options.image)
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

		# Save?
		save = SaveAs(saveAs, fig, self._plt)
		save.frame()

	def animate(self, movie="", fps=15, dpi=200, saveAs=None, axes=None, **kwargs):
		""" Animates the diagnostic over all its timesteps.
		If the data is 1D, it is plotted as a curve, and is animated for all requested timesteps.
		If the data is 2D, it is plotted as a map, and is animated for all requested timesteps.
		If the data is 0D, it is plotted as a curve as function of time.

		Parameters:
		-----------
		figure: int (default: 1)
			The figure number that is passed to matplotlib.
		axes: (default: None)
			Matplotlib's axes handle on which to plot. If None, make new axes.
		vmin, vmax: floats (default to the current limits)
			Data value limits.
		xmin, xmax, ymin, ymax: floats (default to the current limits)
			Axes limits.
		xfactor, yfactor: floats (default: 1)
			Factors to rescale axes.
		movie: path string (default: "")
			Name of a file to create a movie, such as "movie.avi"
			If movie="" no movie is created.
		fps: int (default: 15)
			Number of frames per second (only if movie requested).
		dpi: int (default: 200)
			Number of dots per inch (only if movie requested).
		saveAs: path string (default: None)
			Name of a directory where to save each frame as figures.
			You can even specify a filename such as mydir/prefix.png
			and it will automatically make successive files showing
			the timestep: mydir/prefix0.png, mydir/prefix1.png, etc.

		Example:
		--------
			S = happi.Open("path/to/my/results")
			S.ParticleBinning(1).animate(vmin=0, vmax=1e14)

			This takes the particle binning diagnostic #1 and plots the resulting array in figure 1 from 0 to 3e14.
		"""
		if not self._validate(): return
		if not self._prepare(): return
		if not self._setAndCheck(**kwargs): return
		self.info()
		ax = self._make_axes(axes)
		fig = ax.figure

		# Reset ctrl-C exception
		import sys
		if hasattr(sys,"last_type"): del sys.last_type

		# Movie requested ?
		mov = Movie(fig, movie, fps, dpi)
		# Save to file requested ?
		save = SaveAs(saveAs, fig, self._plt)
		# Loop times for animation
		for time in self._timesteps:
			if self._verbose: print("timestep "+str(time))
			# plot
			ax.cla()
			if self._animateOnAxes(ax, time) is None: return
			self._plt.draw()
			self._plt.pause(0.00001)
			# Catch ctrl-C
			if hasattr(sys,"last_type"):
				if sys.last_type is KeyboardInterrupt: break
			# Copy to movie or file
			mov.grab_frame()
			save.frame(time)
		# Movie ?
		if mov.writer is not None: mov.finish()

	# Method to select specific timesteps among those available in times
	def _selectTimesteps(self, timesteps, times):
		ts = self._np.array(self._np.double(timesteps),ndmin=1)
		if ts.size==2:
			# get all times in between bounds
			times = times[ (times>=ts[0]) * (times<=ts[1]) ]
		elif ts.size==1:
			# get nearest time
			times = self._np.array([times[(self._np.abs(times-ts)).argmin()]])
		else:
			raise
		return times

	# Method to select portion of a mesh based on a slice
	def _selectSubset(self, portion, meshpoints, axisname, axisunits, operation):
		try:
			s = self._np.double(portion)
			if s.size>3 or s.size<1: raise
		except:
			self._error += ["`"+operation+"` along axis "+axisname+" should be a list of 1 to 3 floats"]
			raise
		step = 1
		if s.size==1:
			indices = self._np.array([(self._np.abs(meshpoints-s)).argmin()])
		else:
			indices = self._np.nonzero( (meshpoints>=s[0]) * (meshpoints<=s[1]) )[0]
			if indices.size == 0:
				indices = self._np.array([(self._np.abs(meshpoints-s[:2].mean())).argmin()])
			if s.size==3:
				try:
					step = int(s[2])
					if step - s[2] != 0: raise
				except:
					self._error += ["`"+operation+"` along axis "+axisname+": third number must be an integer"]
					raise
				indices = indices[::step]
		if indices.size == 0:
			self._error += ["`"+operation+"` along "+axisname+" is out of range"]
			raise
		elif indices.size == 1:
			info = operation+" at "+axisname+" = "+str(meshpoints[indices])+" "+axisunits
			selection = self._np.s_[indices[0]]
			finalShape = 1
		else:
			info = operation+" for "+axisname+" from "+str(meshpoints[indices[0]])+" to "+str(meshpoints[indices[-1]])+" "+axisunits
			if step > 1: info += " every "+str(step)+" cells"
			selection = self._np.s_[indices[0]:indices[-1]+1:step]
			finalShape = len(indices)
		return info, selection, finalShape

	# Method to select portion of a mesh based on a range
	def _selectRange(self, portion, meshpoints, axisname, axisunits, operation, edgeInclusive=False):
		# if portion is "all", then select all the axis
		if portion == "all":
			info = operation+" for all "+axisname
			selection = self._np.s_[:]
			finalShape = meshpoints.size
		# Otherwise, parse the portion
		else:
			try:
				s = self._np.double(portion)
				if s.size>2 or s.size<1: raise
			except:
				self._error += ["`"+operation+"` along axis "+axisname+" should be one or two floats"]
				raise
			if s.size==1:
				indices = self._np.array([(self._np.abs(meshpoints-s)).argmin()])
			elif s.size==2:
				indices = self._np.nonzero( (meshpoints>=s[0]) * (meshpoints<=s[1]) )[0]
				if indices.size == 0:
					indices = self._np.array([(self._np.abs(meshpoints-s.mean())).argmin()])
			if indices.size == 0:
				self._error += ["`"+operation+"` along "+axisname+" is out of range"]
				raise
			elif indices.size == 1:
				info = operation+" at "+axisname+" = "+str(meshpoints[indices])+" "+axisunits
				selection = slice(indices[0],indices[0]+1)
				finalShape = 1
			else:
				if edgeInclusive:
					axismin = "-infinity" if indices[ 0]==0                 else str(meshpoints[indices[ 0]])+" "+axisunits
					axismax =  "infinity" if indices[-1]==len(meshpoints)-1 else str(meshpoints[indices[-1]])+" "+axisunits
					info = operation+" for "+axisname+" from "+axismin+" to "+axismax
				else:
					info = operation+" for "+axisname+" from "+str(meshpoints[indices[0]])+" to "+str(meshpoints[indices[-1]])+" "+axisunits
				selection = slice(indices[0],indices[-1])
				finalShape = indices[-1] - indices[0]
		return info, selection, finalShape

	# Method to prepare some data before plotting
	def _prepare(self):
		self._prepare1()
		if not self._prepare2(): return False
		if not self._prepare3(): return False
		self._prepare4()
		return True

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
			if   self.dim == 0: self._animateOnAxes = self._animateOnAxes_0D
			elif self.dim == 1: self._animateOnAxes = self._animateOnAxes_1D
			elif self.dim == 2: self._animateOnAxes = self._animateOnAxes_2D
			else:
				print("Cannot plot in "+str(self.dim)+" dimensions !")
				return False
		# prepare t label
		self._tlabel = self.units.tname
		if self.options.xfactor: self._tlabel += "/"+str(self.options.xfactor)
		self._tlabel = 'Time ( '+self._tlabel+' )'
		# prepare x label
		if self.dim > 0:
			self._xlabel = self.units.xname
			if self.options.xfactor: self._xlabel += "/"+str(self.options.xfactor)
			self._xlabel = self._label[0] + " (" + self._xlabel + ")"
			if self._log[0]: self._xlabel = "Log[ "+self._xlabel+" ]"
		# prepare y label
		if self.dim > 1:
			self._ylabel = self.units.yname
			if self.options.yfactor: self._ylabel += "/"+str(self.options.yfactor)
			self._ylabel = self._label[1] + " (" + self._ylabel + ")"
			if self._log[1]: self._ylabel = "Log[ "+self._ylabel+" ]"
			# prepare extent for 2d plots
			self._extent = [
				self._xfactor*self._centers[0][0],
				self._xfactor*self._centers[0][-1],
				self._yfactor*self._centers[1][0],
				self._yfactor*self._centers[1][-1]
			]
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
		if self.dim==1:
			self._ylabel = self._vlabel
			self._vlabel = ""
		return True
	def _prepare3(self):
		# prepare temporary data if zero-d plot
		if self.dim == 0 and self._tmpdata is None:
			self._tmpdata = self._np.zeros(self._timesteps.size)
			for i, t in enumerate(self._timesteps):
				self._tmpdata[i] = self._dataAtTime(t)
		# prepare the colormap if 2d plot
		if self.dim == 2 and self.options.transparent:
			cmap = self.options.image["cmap"]
			if type(cmap)==str: cmap = self._plt.matplotlib.cm.get_cmap(cmap)
			d = cmap._segmentdata
			new_cmap = self._plt.matplotlib.colors.LinearSegmentedColormap("tmp_cmap", cmap._segmentdata, N=256, gamma=1.0)
			if self.options.transparent in ["both", "under"]:
				new_cmap.set_under(color="white", alpha="0")
			if self.options.transparent in ["both", "over"]:
				new_cmap.set_over (color="white", alpha="0")
			self.options.image["cmap"] = new_cmap
		return True

	def _prepare4(self): pass

	# Method to set limits to a plot
	def _setLimits(self, ax, xmin=None, xmax=None, ymin=None, ymax=None):
		if xmin is not None: ax.set_xlim(xmin=xmin)
		if xmax is not None: ax.set_xlim(xmax=xmax)
		if ymin is not None: ax.set_ylim(ymin=ymin)
		if ymax is not None: ax.set_ylim(ymax=ymax)

	# Methods to plot the data when axes are made
	def _animateOnAxes_0D(self, ax, t, cax_id=0):
		times = self._timesteps[self._timesteps<=t]
		A     = self._tmpdata[self._timesteps<=t]
		im, = ax.plot(self._tfactor*times, A, **self.options.plot)
		ax.set_xlabel(self._tlabel)
		self._setLimits(ax, xmax=self._tfactor*self._timesteps[-1], ymin=self.options.vmin, ymax=self.options.vmax)
		self._setSomeOptions(ax, t)
		return im
	def _animateOnAxes_1D(self, ax, t, cax_id=0):
		A = self._dataAtTime(t)
		im, = ax.plot(self._xfactor*(self._xoffset+self._centers[0]), A, **self.options.plot)
		if self._log[0]: ax.set_xscale("log")
		ax.set_xlabel(self._xlabel)
		ax.set_ylabel(self._ylabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.vmin, ymax=self.options.vmax)
		self._setSomeOptions(ax, t)
		return im
	def _animateOnAxes_2D(self, ax, t, cax_id=0):
		A = self._dataAtTime(t)
		im = self._animateOnAxes_2D_(ax, A)
		ax.set_xlabel(self._xlabel)
		ax.set_ylabel(self._ylabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.ymin, ymax=self.options.ymax)
		try: ax.cax
		except: ax.cax = {}
		if "aspect" not in self.options.colorbar.keys() or self.options.colorbar["aspect"]>0:
			try: # if colorbar exists
				ax.cax[cax_id].cla()
				ax.figure.colorbar(mappable=im, cax=ax.cax[cax_id])
			except:
				ax.cax[cax_id] = ax.figure.colorbar(mappable=im, ax=ax, use_gridspec=False, **self.options.colorbar).ax
		self._setSomeOptions(ax, t)
		return im

	# Special case: 2D plot
	# This is overloaded by class "Probe" because it requires to replace imshow
	# Also overloaded by class "Performances" to add a line plot
	def _animateOnAxes_2D_(self, ax, A):
		im = ax.imshow( self._np.rot90(A),
			vmin = self.options.vmin, vmax = self.options.vmax, extent=self._extent, **self.options.image)
		return im

	# set options during animation
	def _setSomeOptions(self, ax, t=None):
		title = []
		if self._vlabel:
			title += [self._vlabel]
		if t is not None:
			title += ["t = %.2f "%(t*self.timestep*self.units.tcoeff)+self.units.tname]
		ax.set_title("  ".join(title))
		for option, value in self.options.axes.items():
			if type(value) is dict:
				getattr(ax, "set_"+option)( **value )
			else:
				getattr(ax, "set_"+option)( value )
		try:
			if len(self.options.xtick)>0: ax.ticklabel_format(axis="x",**self.options.xtick)
		except:
			if self._verbose: print("Cannot format x ticks (typically happens with log-scale)")
			self.xtickkwargs = []
		try:
			if len(self.options.ytick)>0: ax.ticklabel_format(axis="y",**self.options.ytick)
		except:
			if self._verbose: print("Cannot format y ticks (typically happens with log-scale)")
			self.xtickkwargs = []

	# Define and output directory in case of exporting
	def _setExportDir(self, diagName):
		if self.options.export_dir is not None:
			directory = self.options.export_dir
		elif len(self._results_path) == 1:
			directory = self._results_path[0]
		else:
			directory = self._results_path[0] +self._os.sep+ ".."
		directory += self._os.sep + diagName + self._os.sep
		return directory
	
	def _mkdir(self, dir):
		if not self._os.path.exists(dir): self._os.makedirs(dir)
	
	def _dataAtTime(self, t):
		return self._vfactor*self._getDataAtTime(t)
	def _dataLogAtTime(self, t):
		return self._np.log10( self._vfactor*self._getDataAtTime(t) )
	
	# Convert data to VTK format
	def toVTK(self, numberOfPieces=1):
		if not self._validate(): return
		# prepare vfactor
		self._prepare1()

		if self.dim<2 or self.dim>3:
			print ("Cannot export "+str(self.dim)+"D data to VTK")
			return

		self._mkdir(self._exportDir)
		fileprefix = self._exportDir + self._exportPrefix

		spacings = [c[1]-c[0] for c in self._centers]
		extent = []
		for i in range(self.dim): extent += [0,self._shape[i]-1]
		origin = [0.] * self.dim
		ntimes = len(self._timesteps)

		vtk = VTKfile()

		# If 2D data, then do a streak plot
		if self.dim == 2:
			dt = self._timesteps[1]-self._timesteps[0]

			# Get the data
			data = self._np.zeros(list(self._shape)+[ntimes])
			for itime in range(ntimes):
				data[:,:,itime] = self._dataAtTime(self._timesteps[itime])
			arr = vtk.Array(self._np.ascontiguousarray(data.flatten(order='F'), dtype='float32'), self._title)

			# If all timesteps are regularly spaced
			if (self._np.diff(self._timesteps)==dt).all():
				spacings += [dt]
				extent += [0, ntimes-1]
				origin += [self._timesteps[0]]
				vtk.WriteImage(arr, origin, extent, spacings, fileprefix+".pvti", numberOfPieces)
				if self._verbose: print("Successfully exported regular streak plot to VTK, folder='"+self._exportDir)

			# If timesteps are irregular, make an irregular grid
			else:
				vtk.WriteRectilinearGrid(
					(self._shape[0], self._shape[1], ntimes),
					vtk.Array(self._centers[0].astype('float32'), "x"),
					vtk.Array(self._centers[1].astype('float32'), "y"),
					vtk.Array(self._timesteps  .astype('float32'), "t"),
					arr,
					fileprefix+".vtk"
				)
				if self._verbose: print("Successfully exported irregular streak plot to VTK, folder='"+self._exportDir)

		# If 3D data, then do a 3D plot
		elif self.dim == 3:
			for itime in range(ntimes):
				data = self._np.ascontiguousarray(self._dataAtTime(self._timesteps[itime]).flatten(order='F'), dtype='float32')
				arr = vtk.Array(data, self._title)
				# Output using the diag number
				#vtk.WriteImage(arr, origin, extent, spacings, fileprefix+"_"+str(itime)+".pvti", numberOfPieces)
				# Output using the timestep number
				filename = fileprefix+"_{:08d}.pvti".format(int(self._timesteps[itime]))
				if self._verbose: print("* Processing {}".format(filename))
				vtk.WriteImage(arr, origin, extent, spacings, filename, numberOfPieces)
			if self._verbose: print("Successfully exported 3D plot to VTK, folder='"+self._exportDir)
