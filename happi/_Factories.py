from ._Diagnostics.Scalar import Scalar
from ._Diagnostics.Field import Field
from ._Diagnostics.Probe import Probe
from ._Diagnostics.ParticleBinning import ParticleBinning
from ._Diagnostics.Screen import Screen
from ._Diagnostics.RadiationSpectrum import RadiationSpectrum
from ._Diagnostics.TrackParticles import TrackParticles
from ._Diagnostics.Performances import Performances

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
			tmpDiag = Scalar(simulation)
			# Get a list of scalars
			scalars = tmpDiag.getScalars()
			# Create scalars shortcuts
			for scalar in scalars:
				setattr(self, scalar, ScalarFactory(simulation, scalar))
		
		else:
			# the scalar is saved for generating the object in __call__
			self._additionalArgs += (scalar, )

	def __call__(self, *args, **kwargs):
		return Scalar(self._simulation, *(self._additionalArgs+args), **kwargs)


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
				# Get a list of fields
				fields = simulation.fieldInfo(diagNumber)["fields"]
				# Create fields shortcuts
				for field in fields:
					child = FieldFactory(simulation, diagNumber, field)
					setattr(self, field, child)
			
			else:
				# the field is saved for generating the object in __call__
				self._additionalArgs += (field, )
	
	def __call__(self, *args, **kwargs):
		return Field(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def __repr__(self):
		msg = object.__repr__(self)
		if len(self._additionalArgs) == 0 and self._simulation._scan:
			diag_numbers = self._simulation._diag_numbers["Fields"]
			diag_names = self._simulation._diag_names["Fields"]
			if diag_numbers:
				msg += "\nAvailable Field diagnostics:"
				for n, name in zip(diag_numbers, diag_names):
					msg += "\n\t#%d"%n + (" ('%s')"%name if name else "")
					fields = sorted(self._simulation.fieldInfo(n)["fields"])
					if fields:
						msg += ", with fields: " + ", ".join(fields)
					else:
						msg += ", no fields"
			else:
				msg += "\nNo Field diagnostics available"
		return msg

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
				# Get a list of fields
				fields = self._simulation.probeInfo(probeNumber)["fields"]
				# Create fields shortcuts
				for field in fields:
					setattr(self, field, ProbeFactory(simulation, probeNumber, field))

			else:
				# the field is saved for generating the object in __call__
				self._additionalArgs += (field, )
	
	def __call__(self, *args, **kwargs):
		return Probe(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def __repr__(self):
		msg = object.__repr__(self)
		if len(self._additionalArgs) == 0 and self._simulation._scan:
			diag_numbers = self._simulation._diag_numbers["Probes"]
			diag_names = self._simulation._diag_names["Probes"]
			if diag_numbers:
				msg += "\nAvailable Probe diagnostics:"
				for n, name in zip(diag_numbers, diag_names):
					msg += "\n\t#%d"%n + (" ('%s')"%name if name else "")
					fields = sorted(self._simulation.probeInfo(n)["fields"])
					if fields:
						msg += ", with fields: " + ", ".join(fields)
					else:
						msg += ", no fields"
			else:
				msg += "\nNo Probe diagnostics available"
		return msg


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
	
	def __call__(self, *args, **kwargs):
		return ParticleBinning(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def __repr__(self):
		msg = object.__repr__(self)
		if len(self._additionalArgs) == 0 and self._simulation._scan:
			diag_numbers = self._simulation._diag_numbers["ParticleBinning"]
			diag_names = self._simulation._diag_names["ParticleBinning"]
			if diag_numbers:
				msg += "\nAvailable ParticleBinning diagnostics:\n"
				for n, name in zip(diag_numbers, diag_names):
					msg += "\n" + ParticleBinning._printInfo(ParticleBinning._getInfo(self._simulation, "ParticleBinning", n))
			else:
				msg += "\nNo ParticleBinning diagnostics available"
		return msg

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
		
	def __call__(self, *args, **kwargs):
		return RadiationSpectrum(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def __repr__(self):
		msg = object.__repr__(self)
		if len(self._additionalArgs) == 0 and self._simulation._scan:
			diag_numbers = self._simulation._diag_numbers["RadiationSpectrum"]
			diag_names = self._simulation._diag_names["RadiationSpectrum"]
			if diag_numbers:
				msg += "\nAvailable RadiationSpectrum diagnostics:\n"
				for n, name in zip(diag_numbers, diag_names):
					msg += "\n" + RadiationSpectrum._printInfo(RadiationSpectrum._getInfo(self._simulation, "RadiationSpectrum", n))
			else:
				msg += "\nNo RadiationSpectrum diagnostics available"
		return msg


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
		return Performances(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def __repr__(self):
		msg = object.__repr__(self)
		info = self._simulation.performanceInfo()
		msg += "\nAvailable quantities:\n"+", ".join([str(q) for q in info["quantities_uint"]+info["quantities_double"]])
		return msg


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
	
	def __call__(self, *args, **kwargs):
		return Screen(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def __repr__(self):
		msg = object.__repr__(self)
		if len(self._additionalArgs) == 0 and self._simulation._scan:
			diag_numbers = self._simulation._diag_numbers["Screen"]
			diag_names = self._simulation._diag_names["Screen"]
			if diag_numbers:
				msg += "\nAvailable Screen diagnostics:\n"
				for n, name in zip(diag_numbers, diag_names):
					msg += "\n" + Screen._printInfo(Screen._getInfo(self._simulation, "Screen", n))
			else:
				msg += "\nNo Screen diagnostics available"
		return msg


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
			# Get a list of species
			specs = self._simulation.getTrackSpecies()
			# Create species shortcuts
			for spec in specs:
				child = TrackParticlesFactory(simulation, spec)
				setattr(self, spec, child)
				self._children += [child]

		else:
			# the species is saved for generating the object in __call__
			self._additionalKwargs.update( {"species":species} )
	
	def __call__(self, *args, **kwargs):
		kwargs.update(self._additionalKwargs)
		return TrackParticles(self._simulation, *args, **kwargs)
	
	def __repr__(self):
		msg = object.__repr__(self)
		if len(self._additionalKwargs) == 0 and self._simulation._scan:
			specs = self._simulation.getTrackSpecies()
			if specs:
				msg += "\nAvailable TrackParticles diagnostics:\n" + ", ".join(specs)
			else:
				msg += "\nNo TrackParticles diagnostics available"
		return msg
