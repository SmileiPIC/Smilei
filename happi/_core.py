from ._Utils import *
from ._Diagnostics import Scalar, Field, Probe, ParticleBinning, RadiationSpectrum, Performances, Screen, TrackParticles


class ScalarFactory(object):
	"""Import and analyze a scalar diagnostic from a Smilei simulation
	
	Parameters:
	-----------
	scalar : string (optional)
		The name of the scalar ("Utot", "Ubal", etc.)
		To get a list of available scalars, simply omit this argument.
	timesteps : int or [int, int] (optional)
		If omitted, all timesteps are used.
		If one number  given, the nearest timestep available is used.
		If two numbers given, all the timesteps in between are used.
	units : A units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	
	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		scalar = S.Scalar(...)               # Load the scalar diagnostic
		scalar.get()                         # Obtain the data
	"""
	
	def __init__(self, simulation, scalar=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		if not simulation._scan: return
		
		# If not a specific scalar (root level), build a list of scalar shortcuts
		if scalar is None:
			if simulation._verbose: print("Scanning for Scalar diagnostics")
			# Create a temporary, empty scalar diagnostic
			tmpDiag = Scalar.Scalar(simulation)
			# Get a list of scalars
			scalars = tmpDiag.getScalars()
			# Create scalars shortcuts
			for scalar in scalars:
				setattr(self, scalar, ScalarFactory(simulation, scalar))
		
		else:
			# the scalar is saved for generating the object in __call__
			self._additionalArgs += (scalar, )

	def __call__(self, *args, **kwargs):
		return Scalar.Scalar(self._simulation, *(self._additionalArgs+args), **kwargs)


class FieldFactory(object):
	"""Import and analyze a Field diagnostic from a Smilei simulation

	Parameters:
	-----------
	field : string (optional)
		The name of the field ("Ex", "Ey", etc.)
		To get a list of available fields, simply omit this argument.
		You may write an operation instead of just one field, e.g. "Jx+Jy".
	timesteps : int or [int, int] (optional)
		If omitted, all timesteps are used.
		If one number  given, the nearest timestep available is used.
		If two numbers given, all the timesteps in between are used.
	subset: a python dictionary of the form { axis:range, ... } (optional)
		`axis` must be "x", "y" or "z".
		`range` must be a list of 1 to 3 floats such as [start, stop, step]
		WARNING: THE VALUE OF `step` IS A NUMBER OF CELLS.
		Only the data within the chosen axes' selections is extracted.
		Example: subset = {"y":[10, 80, 4]}
	average : a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y" or "z".
		`range` may be "all", a float, or [float, float].
		For instance, average={"x":"all", "y":[2,3]}.
		The average of all values within the bounds is computed.
	units : a units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	export_dir : the directory to export to VTK

	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		field = S.Field(...)                 # Load the field diagnostic
		field.get()                          # Obtain the data
	"""
	
	def __init__(self, simulation, diagNumber=None, field=None, timestep=None, availableTimesteps=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		self._children = []
		if not simulation._scan: return
		
		# If not a specific diag (root level), build a list of diag shortcuts
		if diagNumber is None:
			if simulation._verbose: print("Scanning for Field diagnostics")
			# Create diags shortcuts
			for diag in simulation._diag_numbers["Fields"]:
				child = FieldFactory(simulation, diag)
				setattr(self, "Field"+str(diag), child)
				self._children += [child]
		
		else:
			# the diag is saved for generating the object in __call__
			self._additionalArgs += (diagNumber, )
			
			# If not a specific field, build a list of field shortcuts
			if field is None:
				# Create a temporary, empty field diagnostic
				tmpDiag = Field.Field(simulation, diagNumber)
				# Get a list of fields
				fields = tmpDiag.getFields()
#				# Get a list of timesteps
#				timesteps = tmpDiag.getAvailableTimesteps()
#				# Create fields shortcuts
#				for field in fields:
#					child = FieldFactory(simulation, diagNumber, field, availableTimesteps=timesteps)
#					setattr(self, field, child)
				# Create fields shortcuts
				for field in fields:
					child = FieldFactory(simulation, diagNumber, field)
					setattr(self, field, child)

			else:
				# the field is saved for generating the object in __call__
				self._additionalArgs += (field, )

#				# If not a specific timestep, build a list of timesteps shortcuts
#				if timestep is None:
#					for timestep in availableTimesteps:
#						setattr(self, 't%0.10i'%timestep, FieldFactory(simulation, diagNumber, field, timestep))
#
#				else:
#					# the timestep is saved for generating the object in __call__
#					self._additionalArgs += (timestep, )

	def __call__(self, *args, **kwargs):
		return Field.Field(self._simulation, *(self._additionalArgs+args), **kwargs)


class ProbeFactory(object):
	"""Import and analyze a probe diagnostic from a Smilei simulation

	Parameters:
	-----------
	probeNumber : int (optional)
		Index of an available probe.
		To get a list of available probes, simply omit this argument.
	field : string (optional)
		Name of the field ("Ex", "Ey", etc)
		To get a list of available fields, simply omit this argument.
	timesteps : int or [int, int] (optional)
		If omitted, all timesteps are used.
		If one number  given, the nearest timestep available is used.
		If two numbers given, all the timesteps in between are used.
	subset: a python dictionary of the form { axis:range, ... } (optional)
		`axis` must be "axis1", "axis2" or "axis3".
		`range` must be a list of 1 to 3 floats such as [start, stop, step]
		WARNING: THE VALUE OF `step` IS A NUMBER OF PROBE POINTS.
		Only the data within the chosen axes' selections is extracted.
		Example: subset = {"axis2":[10, 80, 4]}
	average : a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "axis1" or "axis2" (the probe axes).
		`range` may be "all", a float, or [float, float].
		For instance, average={"axis1":"all", "axis2":[2,3]}.
		The average of all values within the bounds is computed.
	units : A units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	export_dir : the directory to export to VTK

	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		probe = S.Probe(...)                 # Load the probe diagnostic
		probe.get()                          # Obtain the data
	"""

	def __init__(self, simulation, probeNumber=None, field=None, timestep=None, availableTimesteps=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		if not simulation._scan: return

		# If not a specific probe, build a list of probe shortcuts
		if probeNumber is None:
			if simulation._verbose: print("Scanning for Probe diagnostics")
			# Create probe shortcuts
			for probe in simulation._diag_numbers["Probes"]:
				setattr(self, 'Probe'+str(probe), ProbeFactory(simulation, probe))

		else:
			# the probe is saved for generating the object in __call__
			self._additionalArgs += (probeNumber,)

			# If not a specific field, build a list of field shortcuts
			if field is None:
				# Create a temporary, empty probe diagnostic
				tmpDiag = Probe.Probe(simulation, probeNumber)
				# Get a list of fields
				fields = tmpDiag.getFields()
#				# Get a list of timesteps
#				timesteps = tmpDiag.getAvailableTimesteps()
#				# Create fields shortcuts
#				for field in fields:
#					setattr(self, field, ProbeFactory(simulation, probeNumber, field, availableTimesteps=timesteps))
				# Create fields shortcuts
				for field in fields:
					setattr(self, field, ProbeFactory(simulation, probeNumber, field))

			else:
				# the field is saved for generating the object in __call__
				self._additionalArgs += (field, )

#				# If not a specific timestep, build a list of timesteps shortcuts
#				if timestep is None:
#					for timestep in availableTimesteps:
#						setattr(self, 't%0.10i'%timestep, ProbeFactory(simulation, probeNumber, field, timestep))
#
#				else:
#					# the timestep is saved for generating the object in __call__
#					self._additionalArgs += (timestep, )

	def __call__(self, *args, **kwargs):
		return Probe.Probe(self._simulation, *(self._additionalArgs+args), **kwargs)



class ParticleBinningFactory(object):
	"""Import and analyze a ParticleBinning diagnostic from a Smilei simulation

	Parameters:
	-----------
	diagNumber : int (optional)
		Index of an available ParticleBinning diagnostic.
		To get a list of available diags, simply omit this argument.
	timesteps : int or [int, int] (optional)
		If omitted, all timesteps are used.
		If one number  given, the nearest timestep available is used.
		If two numbers given, all the timesteps in between are used.
	subset: a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y", "z", "px", "py", "pz", "p", "gamma", "ekin", "vx", "vy", "vz", "v" or "charge".
		`range` must be a list of 1 to 3 floats such as [start, stop, step]
		WARNING: THE VALUE OF `step` IS A NUMBER OF BINS.
		Only the data within the chosen axes' selections is extracted.
		Example: subset = {"y":[10, 80, 4]}
	sum : a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y", "z", "px", "py", "pz", "p", "gamma", "ekin", "vx", "vy", "vz", "v" or "charge".
		`range` may be "all", a float, or [float, float].
		For instance, sum={"x":"all", "y":[2,3]}.
		The sum of all values within the bounds is computed.
	units : A units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	export_dir : the directory to export to VTK

	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		part = S.ParticleBinning(...)        # Load the particle binning diagnostic
		part.get()                           # Obtain the data
	"""

	def __init__(self, simulation, diagNumber=None, timestep=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		if not simulation._scan: return

		# If not a specific diag (root level), build a list of diag shortcuts
		if diagNumber is None:
			if simulation._verbose: print("Scanning for ParticleBinning diagnostics")
			# Create diags shortcuts
			for diag in simulation._diag_numbers["ParticleBinning"]:
				setattr(self, 'Diag'+str(diag), ParticleBinningFactory(simulation, diag))

		else:
			# the diag is saved for generating the object in __call__
			self._additionalArgs += (diagNumber, )

			## If not a specific timestep, build a list of timesteps shortcuts
			#if timestep is None:
			#	# Create a temporary, empty particle binning diagnostic
			#	tmpDiag = ParticleBinning.ParticleBinning(simulation, diagNumber)
			#	# Get a list of timesteps
			#	timesteps = tmpDiag.getAvailableTimesteps()
			#	# Create timesteps shortcuts
			#	for timestep in timesteps:
			#		setattr(self, 't%0.10i'%timestep, ParticleBinningFactory(simulation, diagNumber, timestep))
			#
			#else:
			#	# the timestep is saved for generating the object in __call__
			#	self._additionalArgs += (timestep, )

	def __call__(self, *args, **kwargs):
		return ParticleBinning.ParticleBinning(self._simulation, *(self._additionalArgs+args), **kwargs)


class RadiationSpectrumFactory(object):
	"""Import and analyze a RadiationSpectrum diagnostic from a Smilei simulation

	Parameters:
	-----------
	diagNumber : int (optional)
		Index of an available RadiationSpectrum diagnostic.
		To get a list of available diags, simply omit this argument.
	timesteps : int or [int, int] (optional)
		If omitted, all timesteps are used.
		If one number  given, the nearest timestep available is used.
		If two numbers given, all the timesteps in between are used.
	subset: a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y", "z", "px", "py", "pz", "p", "gamma", "ekin", "vx", "vy", "vz", "v" or "charge".
		`range` must be a list of 1 to 3 floats such as [start, stop, step]
		WARNING: THE VALUE OF `step` IS A NUMBER OF BINS.
		Only the data within the chosen axes' selections is extracted.
		Example: subset = {"y":[10, 80, 4]}
	sum : a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y", "z", "px", "py", "pz", "p", "gamma", "ekin", "vx", "vy", "vz", "v" or "charge".
		`range` may be "all", a float, or [float, float].
		For instance, sum={"x":"all", "y":[2,3]}.
		The sum of all values within the bounds is computed.
	units : A units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.

	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		rad = S.RadiationSpectrum(...)      # Load the RadiationSpectrum diagnostic
		rad.get()                           # Obtain the data
	"""

	def __init__(self, simulation, diagNumber=None, timestep=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		if not simulation._scan: return
		
		# If not a specific diag (root level), build a list of diag shortcuts
		if diagNumber is None:
			if simulation._verbose: print("Scanning for RadiationSpectrum diagnostics")
			# Create diags shortcuts
			for diag in simulation._diag_numbers["RadiationSpectrum"]:
				setattr(self, 'Diag'+str(diag), RadiationSpectrumFactory(simulation, diag))

		else:
			# the diag is saved for generating the object in __call__
			self._additionalArgs += (diagNumber, )

			## If not a specific timestep, build a list of timesteps shortcuts
			#if timestep is None:
			#	# Create a temporary, empty radiation spectrum diagnostic
			#	tmpDiag = RadiationSpectrum.RadiationSpectrum(simulation, diagNumber)
			#	# Get a list of timesteps
			#	timesteps = tmpDiag.getAvailableTimesteps()
			#	# Create timesteps shortcuts
			#	for timestep in timesteps:
			#		setattr(self, 't%0.10i'%timestep, RadiationSpectrumFactory(simulation, diagNumber, timestep))
			#
			#else:
			#	# the timestep is saved for generating the object in __call__
			#	self._additionalArgs += (timestep, )

	def __call__(self, *args, **kwargs):
		return RadiationSpectrum.RadiationSpectrum(self._simulation, *(self._additionalArgs+args), **kwargs)



class PerformancesFactory(object):
	"""Import and analyze a Performances diagnostic from a Smilei simulation

	Parameters:
	-----------
	raw : a string
		(incompatible with `map` or `histogram`)
		The name of a quantity, or an operation between them (see quantities below).
		The requested quantity is obtained for each process.
	map : a string
		(incompatible with `raw` or `histogram`)
		The name of a quantity, or an operation between them (see quantities below).
		The requested quantity is obtained vs. space coordinates.
	histogram : ["quantity", min, max, nsteps]
		(incompatible with `raw` or `map`)
		Makes a histogram of the requested quantity between min an max, with nsteps bins.
		The quantity may be an operation between the quantities listed further below.
	timesteps : int or [int, int] (optional)
		If omitted, all timesteps are used.
		If one number  given, the nearest timestep available is used.
		If two numbers given, all the timesteps in between are used.
	units : A units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	export_dir : the directory to export to VTK

	Available "quantities":
	-----------------------
	hindex                     : the starting index of each proc in the hilbert curve
	number_of_cells            : the number of cells in each proc
	number_of_particles        : the number of particles in each proc (except frozen ones)
	number_of_frozen_particles : the number of frozen particles in each proc
	total_load                 : the `load` of each proc (number of particles and cells with cell_load coefficient)
	timer_global               : global simulation time (only available for proc 0)
	timer_particles            : time spent computing particles by each proc
	timer_maxwell              : time spent solving maxwell by each proc
	timer_densities            : time spent projecting densities by each proc
	timer_collisions           : time spent computing collisions by each proc
	timer_movWindow            : time spent handling the moving window by each proc
	timer_loadBal              : time spent balancing the load by each proc
	timer_syncPart             : time spent synchronzing particles by each proc
	timer_syncField            : time spent synchronzing fields by each proc
	timer_syncDens             : time spent synchronzing densities by each proc
	timer_total                : the sum of all timers above (except timer_global)

	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		part = S.Performances(...)           # Load the particle binning diagnostic
		part.get()                           # Obtain the data
	"""

	def __init__(self, simulation):
		self._simulation = simulation
		self._additionalArgs = tuple()
		if not simulation._scan: return
		if simulation._verbose: print("Scanning for Performance diagnostics")
	
	def __call__(self, *args, **kwargs):
		return Performances.Performances(self._simulation, *(self._additionalArgs+args), **kwargs)



class ScreenFactory(object):
	"""Import and analyze a screen diagnostic from a Smilei simulation

	Parameters:
	-----------
	diagNumber : int (optional)
		Index of an available screen diagnostic.
		To get a list of available diags, simply omit this argument.
	timesteps : int or [int, int] (optional)
		If omitted, all timesteps are used.
		If one number  given, the nearest timestep available is used.
		If two numbers given, all the timesteps in between are used.
	subset: a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y", "z", "px", "py", "pz", "p", "gamma", "ekin", "vx", "vy", "vz", "v" or "charge".
		`range` must be a list of 1 to 3 floats such as [start, stop, step]
		WARNING: THE VALUE OF `step` IS A NUMBER OF BINS.
		Only the data within the chosen axes' selections is extracted.
		Example: subset = {"y":[10, 80, 4]}
	sum : a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y", "z", "a", "b", "theta", "phi", "px", "py", "pz", "p", "gamma", "ekin", "vx", "vy", "vz", "v" or "charge".
		`range` may be "all", a float, or [float, float].
		For instance, sum={"x":"all", "y":[2,3]}.
		The sum of all values within the bounds is computed.
	units : A units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	export_dir : the directory to export to VTK

	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		screen = S.Screen(...)               # Load the Screen diagnostic
		screen.get()                         # Obtain the data
	"""

	def __init__(self, simulation, diagNumber=None, timestep=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		if not simulation._scan: return

		# If not a specific diag (root level), build a list of diag shortcuts
		if diagNumber is None:
			if simulation._verbose: print("Scanning for Screen diagnostics")
			# Create diags shortcuts
			for diag in simulation._diag_numbers["Screen"]:
				setattr(self, 'Screen'+str(diag), ScreenFactory(simulation, diag))

		else:
			# the diag is saved for generating the object in __call__
			self._additionalArgs += (diagNumber, )

			## If not a specific timestep, build a list of timesteps shortcuts
			#if timestep is None:
			#	# Create a temporary, empty Screen diagnostic
			#	tmpDiag = Screen.Screen(simulation, diagNumber)
			#	# Get a list of timesteps
			#	timesteps = tmpDiag.getAvailableTimesteps()
			#	# Create timesteps shortcuts
			#	for timestep in timesteps:
			#		setattr(self, 't%0.10i'%timestep, ScreenFactory(simulation, diagNumber, timestep))
			#
			#else:
			#	# the timestep is saved for generating the object in __call__
			#	self._additionalArgs += (timestep, )

	def __call__(self, *args, **kwargs):
		return Screen.Screen(self._simulation, *(self._additionalArgs+args), **kwargs)



class TrackParticlesFactory(object):
	"""Import and analyze tracked particles from a Smilei simulation

	Parameters:
	-----------
	species : string (optional)
		Name of a tracked species.
		To get a list of available tracked species, simply omit this argument.
	select: string (optional)
		Instructions for selecting particles among those available.
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
	length: int (optional)
		The length of each plotted trajectory, in number of timesteps.
	units : A units specification such as ["m","second"]
	axes: A list of axes for plotting the trajectories.
		Each axis is "x", "y", "z", "px", "py" or "pz", "chi".
		Example: axes = ["x"] corresponds to x versus time.
		Example: axes = ["x","y"] correspond to 2-D trajectories.
		Example: axes = ["x","px"] correspond to phase-space trajectories.
	skipAnimation: bool (default: False)
		When True, the plot() will directly show the full trajectory.
	export_dir : the directory to export to VTK

	Usage:
	------
		S = happi.Open("path/to/simulation") # Load the simulation
		track = S.TrackParticles(...)        # Load the tracked-particle diagnostic
		track.get()                          # Obtain the data
	"""

	def __init__(self, simulation, species=None, timestep=None):
		self._simulation = simulation
		self._additionalKwargs = dict()
		self._children = []
		if not simulation._scan: return

		# If not a specific species (root level), build a list of species shortcuts
		if species is None:
			if simulation._verbose: print("Scanning for Tracked particle diagnostics")
			# Create a temporary, empty tracked-particle diagnostic
			tmpDiag = TrackParticles.TrackParticles(simulation)
			# Get a list of species
			specs = tmpDiag.getTrackSpecies()
			# Create species shortcuts
			for spec in specs:
				child = TrackParticlesFactory(simulation, spec)
				setattr(self, spec, child)
				self._children += [child]

		else:
			# the species is saved for generating the object in __call__
			self._additionalKwargs.update( {"species":species} )

			# For now, the following block is de-activated
			# It is not possible to have pre-loaded timesteps because the file ordering
			# would take place, and it takes a long time

			## If not a specific timestep, build a list of timesteps shortcuts
			#if timestep is None:
			#	# Create a temporary, empty tracked-particle diagnostic
			#	tmpDiag = TrackParticles.TrackParticles(simulation, species)
			#	# Get a list of timesteps
			#	timesteps = tmpDiag.getAvailableTimesteps()
			#	# Create timesteps shortcuts
			#	for timestep in timesteps:
			#		setattr(self, 't%0.10i'%timestep, TrackParticlesFactory(simulation, species, timestep))
			#
			#else:
			#	# the timestep is saved for generating the object in __call__
			#	self._additionalKwargs.update( {"timesteps":timestep} )

	def __call__(self, *args, **kwargs):
		kwargs.update(self._additionalKwargs)
		return TrackParticles.TrackParticles(self._simulation, *args, **kwargs)




def Open(*args, **kwargs):
	""" Import a Smilei simulation

	Parameters:
	-----------
	results_path : string or list of strings (default '.').
		Directory containing simulation results, or list of directories.
		Omit this argument if you are already in the results directory.

	reference_angular_frequency_SI : float (default None)
		Sets or change the value of reference_angular_frequency_SI, which may
		be defined in the block Main() of any Smilei namelist.

	show : bool (default True)
		Can be set to False to prevent figures to actually appear on screen.

	verbose : bool (default True)
		If False, no warnings or information are printed.

	scan : bool (default True)
		If False, the HDF5 output files are not initially scanned.

	Returns:
	--------
	A SmileiSimulation object, i.e. a container that holds information about a simulation.

	Attributes of the returned object:
	----------------------------------
	namelist :
		An object that holds the information of the original user namelist.
	Scalar :
		A method to access the `DiagScalar` diagnostic.
	Field :
		A method to access the `DiagField` diagnostic.
	Probe :
		A method to access the `DiagProbe` diagnostic.
	ParticleBinning :
		A method to access the `DiagParticleBinning` diagnostic.
	TrackParticles :
		A method to access the tracked particles diagnostic.
	Performances :
		A method to access the `Performances` diagnostic.

	"""
	return SmileiSimulation(*args, **kwargs)


class SmileiSimulation(object):
	"""Object for handling the outputs of a Smilei simulation

	Attributes:
	-----------
	namelist :
		An object that holds the information of the original user namelist.
	Scalar :
		A method to access the `DiagScalar` diagnostic.
	Field :
		A method to access the `DiagField` diagnostic.
	Probe :
		A method to access the `DiagProbe` diagnostic.
	ParticleBinning :
		A method to access the `DiagParticleBinning` diagnostic.
	TrackParticles :
		A method to access the tracked particles diagnostic.
	Performances :
		A method to access the `Performances` diagnostic.

	"""

	def __init__(self, results_path=".", reference_angular_frequency_SI=None, show=True, verbose=True, scan=True):
		self.valid = False
		# Import packages
		import h5py
		import numpy as np
		import os, glob, re
		setMatplotLibBackend(show=show)
		updateMatplotLibColormaps()
		import matplotlib.pyplot
		import matplotlib.pylab as pylab
		pylab.ion()
		# Transfer packages to local attributes
		self._results_path = results_path
		self._h5py = h5py
		self._np = np
		self._os = os
		self._glob = glob.glob
		self._re = re
		self._plt = matplotlib.pyplot
		self._mtime = 0
		self._verbose = verbose
		self._reference_angular_frequency_SI = reference_angular_frequency_SI
		self._scan = scan

		# Load the simulation (verify the path, get the namelist)
		self.reload()

		# Load diagnostics factories
		if self.valid:
			self._diag_numbers = {}
			self._diag_names = {}
			for diagType in ["Fields", "Probes", "ParticleBinning", "Screen", "RadiationSpectrum"]:
				self._diag_numbers[diagType], self._diag_names[diagType] = None, None
				if self._scan:
					self.getDiags(diagType)
			
			self.Scalar = ScalarFactory(self)
			self.Field = FieldFactory(self)
			self.Probe = ProbeFactory(self)
			self.ParticleBinning = ParticleBinningFactory(self)
			self.RadiationSpectrum = RadiationSpectrumFactory(self)
			self.Performances = PerformancesFactory(self)
			self.Screen = ScreenFactory(self)
			self.TrackParticles = TrackParticlesFactory(self)


	def _openNamelist(self, path):
		# Fetch the python namelist
		namespace={}
		exec(open(path+self._os.sep+'smilei.py').read(), namespace) # execute the namelist into an empty namespace
		class Namelist: pass # empty class to store the namelist variables
		namelist = Namelist() # create new empty object
		for key, value in namespace.items(): # transfer all variables to this object
			if key[0]=="_": continue # skip builtins
			setattr(namelist, key, value)

		# Get some info on the simulation
		try:
			# get number of dimensions
			error = "Error extracting 'geometry' from the input file"
			ndim_fields, ndim_particles = {
				"1Dcartesian"   : (1,1),
				"2Dcartesian"   : (2,2),
				"3Dcartesian"   : (3,3),
				"AMcylindrical" : (2,3),
			} [namelist.Main.geometry]
			# get box size
			error = "Error extracting 'grid_length' from the input file"
			grid_length = self._np.atleast_1d(self._np.double(namelist.Main.grid_length))
			if grid_length.size != ndim_fields: raise
			# get cell size
			error = "Error extracting 'cell_length' from the input file"
			cell_length = self._np.atleast_1d(self._np.double(namelist.Main.cell_length))
			if cell_length.size != ndim_fields: raise
			# calculate number of cells in each dimension
			ncels = grid_length/cell_length
			# extract time-step
			error = "Error extracting 'timestep' from the input file"
			timestep = self._np.double(namelist.Main.timestep)
			if not self._np.isfinite(timestep): raise
		except:
			print(error)
			return
		try:
			reference_angular_frequency_SI = namelist.Main.reference_angular_frequency_SI
		except:
			reference_angular_frequency_SI = None
		return namelist, ndim_fields, ndim_particles, cell_length, ncels, timestep, reference_angular_frequency_SI

	def reload(self):
		"""Reloads the simulation, if it has been updated"""
		self.valid = False

		# Obtain the path(s) to the simulation(s) results
		if type(self._results_path) is not list:
			self._results_path = [self._results_path]
		allPaths = []
		for path in self._results_path:
			if type(path) is not str:
				print("The `results_path` parameter must be a string or a list of strings")
				return
			validPaths = []
			for match in self._glob(path):
				if self._os.path.isdir(match) and self._os.path.isfile(match+self._os.sep+"smilei.py"):
					validPaths.append(match)
			if len(validPaths)==0 and self._verbose:
				print("WARNING: `"+path+"` does not point to any valid Smilei simulation path")
			allPaths.extend( validPaths )
		self._results_path = allPaths

		if len(self._results_path)==0:
			print("No valid paths to Smilei simulation results have been provided")
			return

		# Check the last modification date and get paths which are newer
		lastmodif = 0
		newPaths = []
		for path in self._results_path:
			thismtime = self._os.path.getmtime(path+self._os.sep+"/smilei.py")
			if thismtime > self._mtime: newPaths.append(path)
			lastmodif = max(lastmodif, thismtime)

		# Reload if necessary
		if lastmodif > self._mtime:
			W_r = None
			# Loop paths and verify the namelist is compatible
			for path in newPaths:
				self.namelist, ndim_fields, ndim_particles, cell_length, ncels, timestep, reference_angular_frequency_SI = self._openNamelist(path)
				try:
					if (
						ndim_fields != self._ndim_fields or
						ndim_particles != self._ndim_particles or
						(cell_length != self._cell_length).any() or
						(ncels != self._ncels).any() or
						timestep != self._timestep or
						reference_angular_frequency_SI != W_r
					):
						print("The simulation in path '"+path+"' is not compatible with the other ones")
						return
				except:
					pass
				self._ndim_fields = ndim_fields
				self._ndim_particles = ndim_particles
				self._cell_length = cell_length
				self._ncels = ncels
				self._timestep = timestep
				W_r = reference_angular_frequency_SI
				if self._verbose: print("Loaded simulation '"+path+"'")
			if self._reference_angular_frequency_SI is None:
				self._reference_angular_frequency_SI = W_r
		
		self._mtime = lastmodif
		self.valid = True
	
	def getDiags(self, diagType):
		if self._diag_numbers[diagType] is None:
			self._diag_numbers[diagType], self._diag_names[diagType] = self.scanDiags(diagType)
		return self._diag_numbers[diagType], self._diag_names[diagType]
	
	def scanDiags(self, diagType):
		diags = []
		for path in self._results_path:
			files = self._glob(path+self._os.sep+diagType+'*.h5')
			these_diags = []
			for file in files:
				# get number
				number = int(self._re.findall(diagType+"([0-9]+).h5$",file)[0])
				# get name
				with self._h5py.File(file, 'r') as f:
					name = f.attrs["name"].decode() if "name" in f.attrs else ""
				these_diags += [(number, name)]
			# Update diags with those of previous paths
			if diags == []: diags = these_diags
			else          : diags = [ d for d in diags if d in these_diags ]
		if diags == []:
			return [], []
		else:
			return zip( *sorted( diags, key=lambda x:x[0] ) )
	
	def __repr__(self):
		if not self.valid:
			return "Invalid Smilei simulation"
		else:
			files = [self._glob(path+self._os.sep+"smilei.py")[0] for path in self._results_path]
			files = "\n\t".join(files)
			return "Smilei simulation with input file(s) located at:\n\t"+files
