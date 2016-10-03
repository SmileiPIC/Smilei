from ._Utils import *
from ._Diagnostics.Scalar import Scalar
import _Diagnostics.Field
from ._Diagnostics.Probe import Probe
from ._Diagnostics.ParticleDiagnostic import ParticleDiagnostic
from ._Diagnostics.TrackParticles import TrackParticles


class FieldFactory(object):
	"""Import and analyze a Field diagnostic from a Smilei simulation
	
	Parameters:
	-----------
	field : string (optional)
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
	data_log : bool (default False)
		If True, then log10 is applied to the output array before plotting.
	stride : int (default 1)
		Step size for sampling the grid.
	
	Usage:
	------
		S = Smilei("path/to/simulation") # Load the simulation
		field = S.Field(...)             # Load the field diagnostic
		field.get()                      # Obtain the data
	"""
	
	def __init__(self, simulation):
		self._simulation = simulation
		# Create a list of shortcuts for fields
		def createField(field):
			def F(*args, **kwargs):
				return _Diagnostics.Field.Field(self._simulation, field=field, *args, **kwargs)
			F.__doc__ = "Alias of Smilei().Field('"+field+"')"
			F.__name__ = str(field)
			return F
		for field in _Diagnostics.Field.Field(self._simulation).getFields():
			setattr(self, field, createField(field))
	
	def __call__(self, *args, **kwargs):
		return _Diagnostics.Field.Field(self._simulation, *args, **kwargs)


class Smilei(object):
	""" Import a Smilei simulation
	
	Parameters:
	-----------
	results_path : string (default '.').
		Directory containing simulation results.
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
	
	def __init__(self, results_path=".", show=True):
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
		
		# Load the simulation (verify the path, get the namelist)
		self.reload()
		
		self.Field = FieldFactory(self)
	
	def reload(self):
		self.valid = False
		# Verify that results_path is valid
		if not self._os.path.isdir(self._results_path):
			print("Could not find directory "+self._results_path)
			return
		if len(self._glob(self._results_path+"/smilei.py"))==0:
			print("Could not find an input file in directory "+self._results_path)
			return
		# Check the last modification date
		lastmodif = self._os.path.getmtime(self._results_path+"/smilei.py")
		if self._mtime < lastmodif:
			# Fetch the python namelist
			namespace={}
			exec(open(self._results_path+'/smilei.py').read(), namespace) # execute the namelist into an empty namespace
			class Namelist: pass # empty class to store the namelist variables
			self.namelist = Namelist() # create new empty object
			for key, value in namespace.items(): # transfer all variables to this object
				if key[0]=="_": continue # skip builtins
				setattr(self.namelist, key, value)
		
		self._mtime = lastmodif
		self.valid = True
	
	def __repr__(self):
		if not self.valid:
			return "Invalid Smilei simulation"
		file = self._glob(self._results_path+"/smilei.py")[0]
		return "Smilei simulation with input file located at `"+file+"`"
	
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
		""" ParticleDiagnostic(diagNumber=None, timesteps=None, slice=None, units=[""], data_log=False, stride=1)
		
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
		stride : step size for reading the grid.
			If the grid is too large, use a stride > 1 to reduce the amount of data.
		
		Usage:
		------
			diag = S.ParticleDiagnostic(...) # S is a Smilei object
			diag.get()
			diag.plot()
		"""
		return ParticleDiagnostic(self, *args, **kwargs)
	
	
	def TrackParticles(self, *args, **kwargs):
		""" TrackParticles(species=None, select="", axes=[], timesteps=None, length=None, units=[""], skipAnimation=False)
		
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
		length: The length of each plotted trajectory, in number of timesteps.
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
	
