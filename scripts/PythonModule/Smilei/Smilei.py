from ._Utils import *
from ._Diagnostics import Scalar, Field, Probe, ParticleDiagnostic, TrackParticles


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
		S = Smilei("path/to/simulation") # Load the simulation
		scalar = S.Scalar(...)           # Load the scalar diagnostic
		scalar.get()                     # Obtain the data
	"""
	
	def __init__(self, simulation, scalar=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		
		# If not a specific scalar (root level), build a list of scalar shortcuts
		if scalar is None:
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
	
	def toXDMF(self):
		pass
	
	
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
	slice : a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "x", "y" or "z".
		`range` may be "all", a float, or [float, float].
		For instance, slice={"x":"all", "y":[2,3]}.
		The average of all values within the 'slice' is computed.
	units : a units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	stride : int (default: 1)
		Step size for sampling the grid.
	
	Usage:
	------
		S = Smilei("path/to/simulation") # Load the simulation
		field = S.Field(...)             # Load the field diagnostic
		field.get()                      # Obtain the data
	"""
	
	def __init__(self, simulation, diagNumber=None, field=None, timestep=None, availableTimesteps=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		self._children = []
		
		# If not a specific diag (root level), build a list of diag shortcuts
		if diagNumber is None:
			# Create a temporary, empty field diagnostic
			tmpDiag = Field.Field(simulation)
			# Get a list of diag
			diags = tmpDiag.getDiags()
			# Create diags shortcuts
			for diag in diags:
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
				# Get a list of timesteps
				timesteps = tmpDiag.getAvailableTimesteps()
				# Create fields shortcuts
				for field in fields:
					child = FieldFactory(simulation, diagNumber, field, availableTimesteps=timesteps)
					setattr(self, field, child)
			
			else:
				# the field is saved for generating the object in __call__
				self._additionalArgs += (field, )
				
				# If not a specific timestep, build a list of timesteps shortcuts
				if timestep is None:
					for timestep in availableTimesteps:
						setattr(self, 't%0.10i'%timestep, FieldFactory(simulation, diagNumber, field, timestep))
				
				else:
					# the timestep is saved for generating the object in __call__
					self._additionalArgs += (timestep, )
	
	def __call__(self, *args, **kwargs):
		return Field.Field(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def toXDMF(self):
		if len(self._children) > 0:
			for child in self._children:
				child.toXDMF()
		else:
			self().toXDMF()


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
	slice : a python dictionary of the form { axis:range, ... } (optional)
		`axis` may be "axis1" or "axis2" (the probe axes).
		`range` may be "all", a float, or [float, float].
		For instance, slice={"axis1":"all", "axis2":[2,3]}.
		The average of all values within the 'slice' is computed.
	units : A units specification such as ["m","second"]
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	
	Usage:
	------
		S = Smilei("path/to/simulation") # Load the simulation
		probe = S.Probe(...)             # Load the probe diagnostic
		probe.get()                      # Obtain the data
	"""
	
	def __init__(self, simulation, probeNumber=None, field=None, timestep=None, availableTimesteps=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		
		# If not a specific probe, build a list of probe shortcuts
		if probeNumber is None:
			# Create a temporary, empty probe diagnostic
			tmpDiag = Probe.Probe(simulation)
			# Get a list of probes
			probes = tmpDiag.getProbes()
			# Create probe shortcuts
			for probe in probes:
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
				# Get a list of timesteps
				timesteps = tmpDiag.getAvailableTimesteps()
				# Create fields shortcuts
				for field in fields:
					setattr(self, field, ProbeFactory(simulation, probeNumber, field, availableTimesteps=timesteps))
			
			else:
				# the field is saved for generating the object in __call__
				self._additionalArgs += (field, )
				
				# If not a specific timestep, build a list of timesteps shortcuts
				if timestep is None:
					for timestep in availableTimesteps:
						setattr(self, 't%0.10i'%timestep, ProbeFactory(simulation, probeNumber, field, timestep))
				
				else:
					# the timestep is saved for generating the object in __call__
					self._additionalArgs += (timestep, )
	
	def __call__(self, *args, **kwargs):
		return Probe.Probe(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def toXDMF(self):
		pass


class ParticleDiagnosticFactory(object):
	"""Import and analyze a particle diagnostic from a Smilei simulation
	
	Parameters:
	-----------
	diagNumber : int (optional)
		Index of an available particle diagnostic.
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
	data_log : bool (default: False)
		If True, then log10 is applied to the output array before plotting.
	stride : int (default: 1)
		Step size for sampling the grid.
	
	Usage:
	------
		S = Smilei("path/to/simulation") # Load the simulation
		part = S.ParticleDiagnostic(...) # Load the particle diagnostic
		part.get()                       # Obtain the data
	"""
	
	def __init__(self, simulation, diagNumber=None, timestep=None):
		self._simulation = simulation
		self._additionalArgs = tuple()
		
		# If not a specific diag (root level), build a list of diag shortcuts
		if diagNumber is None:
			# Create a temporary, empty particle diagnostic
			tmpDiag = ParticleDiagnostic.ParticleDiagnostic(simulation)
			# Get a list of diags
			diags = tmpDiag.getDiags()
			# Create diags shortcuts
			for diag in diags:
				setattr(self, 'Diag'+str(diag), ParticleDiagnosticFactory(simulation, diag))
		
		else:
			# the diag is saved for generating the object in __call__
			self._additionalArgs += (diagNumber, )
			
			# If not a specific timestep, build a list of timesteps shortcuts
			if timestep is None:
				# Create a temporary, empty particle diagnostic
				tmpDiag = ParticleDiagnostic.ParticleDiagnostic(simulation, diagNumber)
				# Get a list of timesteps
				timesteps = tmpDiag.getAvailableTimesteps()
				# Create timesteps shortcuts
				for timestep in timesteps:
					setattr(self, 't%0.10i'%timestep, ParticleDiagnosticFactory(simulation, diagNumber, timestep))
			
			else:
				# the timestep is saved for generating the object in __call__
				self._additionalArgs += (timestep, )
	
	def __call__(self, *args, **kwargs):
		return ParticleDiagnostic.ParticleDiagnostic(self._simulation, *(self._additionalArgs+args), **kwargs)
	
	def toXDMF(self):
		pass


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
		Each axis is "x", "y", "z", "px", "py" or "pz".
		Example: axes = ["x"] corresponds to x versus time.
		Example: axes = ["x","y"] correspond to 2-D trajectories.
		Example: axes = ["x","px"] correspond to phase-space trajectories.
	skipAnimation: bool (default: False)
		When True, the plot() will directly show the full trajectory.
	
	Usage:
	------
		S = Smilei("path/to/simulation") # Load the simulation
		track = S.TrackParticles(...)    # Load the tracked-particle diagnostic
		track.get()                      # Obtain the data
	"""
	
	def __init__(self, simulation, species=None, timestep=None):
		self._simulation = simulation
		self._additionalKwargs = dict()
		self._children = []
		
		# If not a specific species (root level), build a list of species shortcuts
		if species is None:
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
	
	def toXDMF(self):
		if len(self._children) > 0:
			for child in self._children:
				child.toXDMF()
		else:
			self().toXDMF()





class Smilei(object):
	""" Import a Smilei simulation
	
	Parameters:
	-----------
	results_path : string or list of strings (default '.').
		Directory containing simulation results, or list of directories.
		Omit this argument if you are already in the results directory.
	
	show : bool (default True)
		Can be set to False to prevent figures to actually appear on screen.
	
	Returns:
	--------
	A Smilei object, i.e. a container that holds information about a simulation.
	
	Attributes:
	-----------
	namelist :
		An object that holds the information of the original user namelist.
	Scalar :
		A callable object to access the `DiagScalar` diagnostic.
	Field :
		A callable object to access the `DiagField` diagnostic.
	Probe :
		A callable object to access the `DiagProbe` diagnostic.
	ParticleDiagnostic :
		A callable object to access the `DiagParticle` diagnostic.
	TrackParticles :
		A callable object to access the tracked particles diagnostic.
		
	"""
	
	def __init__(self, results_path=".", show=True, verbose=True):
		self.valid = False
		# Import packages
		import h5py
		import numpy as np
		import os, glob, re, sys
		setMatplotLibBackend(show=show)
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
		
		# Load the simulation (verify the path, get the namelist)
		self.reload()
		
		# Load diagnostics factories
		if self.valid:
			self.Scalar = ScalarFactory(self)
			self.Field = FieldFactory(self)
			self.Probe = ProbeFactory(self)
			self.ParticleDiagnostic = ParticleDiagnosticFactory(self)
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
			error = "Error extracting 'dim' from the input file"
			ndim = int(namelist.Main.geometry[0])
			if ndim not in [1,2,3]: raise
			# get box size
			error = "Error extracting 'sim_length' from the input file"
			sim_length = self._np.atleast_1d(self._np.double(namelist.Main.sim_length))
			if sim_length.size != ndim: raise
			# get cell size
			error = "Error extracting 'cell_length' from the input file"
			cell_length = self._np.atleast_1d(self._np.double(namelist.Main.cell_length))
			if cell_length.size != ndim: raise
			# calculate number of cells in each dimension
			ncels = sim_length/cell_length
			# extract time-step
			error = "Error extracting 'timestep' from the input file"
			timestep = self._np.double(namelist.Main.timestep)
			if not self._np.isfinite(timestep): raise
		except:
			print(error)
			return
		try:
			referenceAngularFrequency_SI = namelist.Main.referenceAngularFrequency_SI
		except:
			referenceAngularFrequency_SI = None
		return namelist, ndim, cell_length, ncels, timestep, referenceAngularFrequency_SI
	
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
			if len(validPaths)==0:
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
			# Get the previous simulation parameters
			try:    prevArgs = (self._ndim, self._cell_length, self._ncels, self._timestep, self._referenceAngularFrequency_SI)
			except: prevArgs = ()
			# Loop paths and verify the namelist is compatible
			for path in newPaths:
				args = self._openNamelist(path)
				if len(prevArgs)==0:
					prevArgs = args[1:]
				elif args[1]!=prevArgs[0] or (args[2]!=prevArgs[1]).any() or (args[3]!=prevArgs[2]).any() or args[4:]!=prevArgs[3:]:
					print("The simulation in path '"+path+"' is not compatible with the other ones")
					return
				if self._verbose: print("Loaded simulation '"+path+"'")
			# Update the simulation parameters
			self._ndim, self._cell_length, self._ncels, self._timestep, self._referenceAngularFrequency_SI = args[1:]
			self.namelist = args[0]
		
		self._mtime = lastmodif
		self.valid = True
	
	def __repr__(self):
		if not self.valid:
			return "Invalid Smilei simulation"
		else:
			files = [self._glob(path+self._os.sep+"smilei.py")[0] for path in self._results_path]
			files = "\n\t".join(files)
			return "Smilei simulation with input file(s) located at:\n\t"+files
	
	def toXDMF(self):
		if not self.valid: return
		
		self.Scalar            .toXDMF()
		self.Field             .toXDMF()
		self.Probe             .toXDMF()
		self.ParticleDiagnostic.toXDMF()
		self.TrackParticles    .toXDMF()


