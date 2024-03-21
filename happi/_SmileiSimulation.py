from ._Factories import ScalarFactory, FieldFactory, ProbeFactory, ParticleBinningFactory, RadiationSpectrumFactory, PerformancesFactory, ScreenFactory, TrackParticlesFactory, NewParticlesFactory
from ._Utils import *


PintWarningIssued = False

class SmileiSimulation(object):
	"""Object for handling the outputs of a Smilei simulation

	Attributes:
	-----------
	namelist :
		An object that holds the information of the original user namelist.
	Scalar :
		A method to access the `DiagScalar` diagnostic.
	Field :
		A method to access the `DiagField` diagnostics.
	Probe :
		A method to access the `DiagProbe` diagnostics.
	ParticleBinning :
		A method to access the `DiagParticleBinning` diagnostics.
	Screen :
		A method to access the `Screen` diagnostics.
	RadiationSpectrum :
		A method to access the `RadiationSpectrum` diagnostics.
	TrackParticles :
		A method to access the tracked particles diagnostics.
	Performances :
		A method to access the `Performances` diagnostic.

	"""

	def __init__(self, results_path=".", reference_angular_frequency_SI=None, show=True, verbose=True, scan=True, pint=True):
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
		self._ureg = None
		
		# Load the simulation (verify the path, get the namelist)
		self.reload()
		
		# Load diagnostic factories
		if self.valid:
			
			# Manage units with the pint package
			global PintWarningIssued
			if pint:
				try:
					from pint import UnitRegistry
				except Exception as e:
					if self._verbose and not PintWarningIssued:
						print("WARNING: you do not have the *Pint* package, so you cannot modify units.")
						print("       : The results will stay in code units.")
						PintWarningIssued = True
				else:
					self._ureg = UnitRegistry()
					if self._reference_angular_frequency_SI:
						self._ureg.define("W_r = "+str(self._reference_angular_frequency_SI)+"*hertz") # frequency
					else:
						self._ureg.define("W_r = [reference_frequency]"                 ) # frequency
					self._ureg.define("V_r = speed_of_light"                   ) # velocity
					self._ureg.define("M_r = electron_mass"                    ) # mass
					self._ureg.define("Q_r = 1.602176565e-19 * coulomb"        ) # charge
					self._ureg.define("L_r = V_r / W_r"                        ) # length
					self._ureg.define("T_r = 1   / W_r"                        ) # time
					self._ureg.define("P_r = M_r * V_r"                        ) # momentum
					self._ureg.define("K_r = M_r * V_r**2"                     ) # energy
					self._ureg.define("N_r = epsilon_0 * M_r * W_r**2 / Q_r**2") # density
					self._ureg.define("J_r = V_r * Q_r * N_r"                  ) # current
					self._ureg.define("B_r = M_r * W_r / Q_r"                  ) # magnetic field
					self._ureg.define("E_r = B_r * V_r"                        ) # electric field
					self._ureg.define("S_r = K_r * V_r * N_r"                  ) # poynting
			elif self._verbose and not PintWarningIssued:
				print("WARNING: *Pint* package disabled. The results will stay in code units.")
				PintWarningIssued = True
			
			self.cylindrical = self.namelist.Main.geometry == "AMcylindrical"
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
			self.NewParticles = NewParticlesFactory(self)
			
	def _openNamelist(self, path):
		# empty class to store the namelist variables
		class Namelist: pass
		namelist = Namelist()
		
		# Fetch the python namelist
		if self._scan:
			namespace={}
			with open(path+self._os.sep+'smilei.py') as f:
				exec(f.read(), namespace) # execute the namelist into an empty namespace
			for key, value in namespace.items(): # transfer all variables to this object
				if key[0]=="_": continue # skip builtins
				setattr(namelist, key, value)
		else:
			import shelve
			class Block(object):
				def __init__(self, **kwargs):
					for k,v in kwargs.items():
						setattr(self, k, v)
			with shelve.open(path+self._os.sep+'info.shelf') as f:
				for k in f:
					if k == "_singletons":
						for singletonName, singletonDict in f[k].items():
							setattr(Namelist, singletonName, Block(**singletonDict))
					elif k == "_components":
						for componentName, componentList in f[k].items():
							setattr(Namelist, componentName, [Block(**component) for component in componentList])
					else:
						setattr(Namelist, k, f[k])
		
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
		except Exception as e:
			print(error)
			return
		try:
			reference_angular_frequency_SI = namelist.Main.reference_angular_frequency_SI
		except Exception as e:
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
		self._results_path = sorted(allPaths)

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
				except Exception as e:
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
			diags = list(set(diags+these_diags)) # unique diags
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
	
	def _getParticleListSpecies(self, filePrefix):
		""" List the available species in diagnostics of type ParticleList """
		species = []
		for path in self._results_path:
			files = self._glob(path+self._os.sep+filePrefix+"*.h5")
			for file in files:
				s = self._re.search("^"+filePrefix+"_(.+).h5",self._os.path.basename(file))
				if s: species += [ s.groups()[0] ]
		return list(set(species)) # unique species
	
	# get all available scalars
	def getScalars(self):
		allScalars = None
		for path in self._results_path:
			try:
				file = path+'/scalars.txt'
				f = open(file, 'r')
			except Exception as e:
				continue
			try:
				# Find last commented line
				prevline = ""
				for line in f:
					line = line.strip()
					if line[0]!="#": break
					prevline = line[1:].strip()
				scalars = str(prevline).split() # list of scalars
				scalars = scalars[1:] # remove first, which is "time"
			except Exception as e:
				scalars = []
			f.close()
			if allScalars is None:
				allScalars = scalars
			else:
				allScalars = self._np.intersect1d(allScalars, scalars)
		if allScalars is None:
			return []
		return allScalars
	
	def getTrackSpecies(self):
		""" List the available tracked species """
		return self._getParticleListSpecies("TrackParticlesDisordered")
	
	def getNewParticlesSpecies(self):
		""" List the available NewParticles species """
		return self._getParticleListSpecies("NewParticles")
	
	def fieldInfo(self, diag):
		""" Information on a specific Field diagnostic
		
		Parameters:
		-----------
		diag: the number or name of a Field diagnostic
		
		Returns:
		--------
		A dictionnary containing:
		* "diagNumber": the diagnostic number
		* "diagName": the diagnostic name
		* "fields": list of the available fields in this diagnostic. In the case of
		  `AMcylindrical` geometry, this is a dictionnary with a list of modes for each field.
		"""
		diag_numbers, diag_names = self.getDiags("Fields")
		
		if type(diag) is str:
			if diag not in diag_names:
				raise Exception("No Field diagnostic `"+diag+"` found")
			i = diag_names.index( diag )
		else:
			if diag not in diag_numbers:
				raise Exception("No Field diagnostic #"+str(diag)+" found")
			i = diag_numbers.index( diag )
		diagNumber = diag_numbers[i]
		diagName = diag_names[i]
		
		raw_fields = set()
		for path in self._results_path:
			file = path+self._os.sep+'Fields'+str(diagNumber)+'.h5'
			try:
				f = self._h5py.File(file, 'r')
			except Exception as e:
				continue
			values = f["data"].values()
			if len(values)==0:
				continue
			these_fields =  set(next(iter(values)).keys())
			raw_fields = (raw_fields & these_fields) or these_fields
			f.close()
		
		# Case of a cylindrical geometry
		if self.cylindrical:
			from ._Diagnostics.Field import Field
			fields = {}
			for f in raw_fields:
				fname, imode, envelope = Field._cylindricalMode(f)
				if fname not in fields:
					fields[fname] = []
				fields[fname] += [int(imode)]
		else:
			fields = sorted(list(raw_fields))
			
		return dict( diagNumber=diagNumber, diagName=diagName, fields=fields )
	
	def probeInfo(self, diag):
		""" Information on a specific Probe diagnostic
		
		Parameters:
		-----------
		diag: the number or name of a Probe diagnostic
		
		Returns:
		--------
		A dictionnary containing:
		* "probeNumber": the diagnostic number
		* "probeName": the diagnostic name
		* "fields": list of the available fields in this diagnostic
		"""
		diag_numbers, diag_names = self.getDiags("Probes")
		
		if type(diag) is str:
			if diag not in diag_names:
				raise Exception("No probe diagnostic `"+diag+"` found")
			i = diag_names.index( diag )
		else:
			if diag not in diag_numbers:
				raise Exception("No probe diagnostic #"+str(diag)+" found")
			i = diag_numbers.index( diag )
		probeNumber = diag_numbers[i]
		probeName = diag_names[i]
		
		fields = []
		for path in self._results_path:
			file = path+self._os.sep+"Probes"+str(probeNumber)+".h5"
			try:
				f = self._h5py.File(file, 'r')
			except Exception as e:
				continue
			# Verify that this file is compatible with the previous ones
			try:
				for key, val in verifications.items():
					if f[key][()] != val:
						raise Exception("Probe #"+str(probeNumber)+" in path '"+path+"' is incompatible with the other ones")
			except Exception as e:
				verifications = {"number":f["number"][()]}
				npoints = f["number"].size
				if f["number"][()].prod() > 1:
					npoints += 1
				for i in range(npoints):
					verifications["p"+str(i)] = f["p"+str(i)][()]
			# Get list of fields
			fields_here = _decode(f.attrs["fields"]).split(",")
			if fields and fields != fields_here:
				raise Exception("Probe #"+str(probeNumber)+" in path '"+path+"' is incompatible with the other ones")
			fields = fields_here
			
			f.close()
		
		if not verifications:
			raise Exception("Error opening probe #"+str(probeNumber))
		if not fields:
			raise Exception("No fields found for probe #"+str(probeNumber))
		
		return dict( probeNumber=probeNumber, probeName=probeName, fields=list(fields) )
	
	def performanceInfo(self):
		""" Information on the available quantities in the performance diagnostic
		
		Returns:
		--------
		A dictionnary containing:
		* "quantities_uint": a list of the available integer quantities
		* "quantities_double": a list of the available float quantities
		* "patch_arrangement": the type of patch arrangement
		"""
		
		available_uint   = []
		available_double = []
		patch_arrangement = "?"
		timesteps = set()
		for path in self._results_path:
			file = path+self._os.sep+'Performances.h5'
			try:
				f = self._h5py.File(file, 'r')
			except Exception as e:
				continue
			# Verify all simulations have all quantities
			try:
				quantities_uint   = [_decode(a) for a in f.attrs["quantities_uint"  ]]
				quantities_double = [_decode(a) for a in f.attrs["quantities_double"]]
				if available_uint   and available_uint   != quantities_uint  : raise
				if available_double and available_double != quantities_double: raise
				available_uint   = quantities_uint
				available_double = quantities_double
				if "patch_arrangement" in f.attrs:
					patch_arrangement = _decode(f.attrs["patch_arrangement"])
				timesteps = timesteps.union([int(k) for k in f])
				f.close()
			except Exception as e:
				raise Exception("File '"+file+"' does not seem to contain correct data")
			
		return dict(
			quantities_uint = available_uint,
			quantities_double = available_double,
			patch_arrangement = patch_arrangement,
			timesteps = sorted(timesteps)
			)