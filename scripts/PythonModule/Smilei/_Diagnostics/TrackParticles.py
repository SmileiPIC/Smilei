from .Diagnostic import Diagnostic
from .._Utils import *

# -------------------------------------------------------------------
# Class for tracked particles diagnostics
# -------------------------------------------------------------------
class TrackParticles(Diagnostic):
	# This is the constructor, which creates the object
	def _init(self, species=None, select="", axes=[], timesteps=None, length=None, **kwargs):
		
		# If argument 'species' not provided, then print available species and leave
		if species is None:
			species = self.getTrackSpecies()
			if len(species)>0:
				self._error += "Printing available tracked species:\n"
				self._error += "-----------------------------------\n"
				self._error += "\n".join(species)
			else:
				self._error = "No tracked particles files found in '"+self._results_path+"'"
			return None
		
		# Get info from the hdf5 files + verifications
		# -------------------------------------------------------------------
		self.species  = species
		self._file = self._results_path+"/TrackParticles_"+species+".h5"
		try:
			f = self._h5py.File(self._file, 'r')
		except:
			self._orderFile( self._results_path+"/TrackParticlesDisordered_"+species+".h5", self._file )
			f = self._h5py.File(self._file, 'r')
		self._h5items = list(f.values())
		
		# Get available times in the hdf5 file
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
			self._error = "No tracked particles found in "+self._file
			return
		alltimes = self.times
		# If specific timesteps requested, narrow the selection
		if timesteps is not None:
			try:
				ts = self._np.array(self._np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self._itimes = self._np.nonzero((self.times>=ts[0]) * (self.times<=ts[1]))[0]
					self.times = self.times[ self._itimes ]
				elif ts.size==1:
					# get nearest time
					self._itimes = self._np.array([(self._np.abs(self.times-ts)).argmin()])
					self.times = self._np.array(self.times[ self._itimes ])
				else:
					raise
			except:
				self._error = "Argument `timesteps` must be one or two non-negative integers"
				return
		else:
			self._itimes = self._np.arange(len(self.times))
		# Need at least one timestep
		if self.times.size < 1:
			self._error = "Timesteps not found"
			return
		
		# Get available properties ("x", "y", etc.)
		self._properties = {}
		translateProperties = {"Id":"Id", "x":"Position-0", "y":"Position-1", "z":"Position-2",
			"px":"Momentum-0", "py":"Momentum-1", "pz":"Momentum-2"}
		availableProperties = list(f.keys())
		for k,v in translateProperties.items():
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
		if type(select) is str:
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
							time_indices = self._np.nonzero(eval(s))[0]
						except:
							raise Exception("Error in selector syntax: time selector not understood in "+select[i:i+3]+"()")
						try:
							particleSelector = select[comma+1:parenthesis]
							for prop in self._properties.keys():
								particleSelector = self._re.sub(r"\b"+prop+r"\b", "self._np.double(self._h5items["+str(self._properties[prop])+"][ti,:])", particleSelector)
						except:
							raise Exception("Error in selector syntax: not understood: "+select[i:parenthesis+1])
						if select[i:i+4] == "any(": selection = self._np.array([False]*self.nParticles)
						if select[i:i+4] == "all(": selection = self._np.array([True]*self.nParticles)
						#try:
						ID = self._np.zeros((self.nParticles,), dtype=self._np.int32)
						for ti in time_indices:
							selectionAtTimeT = eval(particleSelector) # array of True or False
							self._Id.read_direct(ID, source_sel=self._np.s_[ti,:], dest_sel=self._np.s_[:]) # read the particle Ids
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
		
		else:
			try:
				self.selectedParticles = self._np.array(select,dtype=int)
			except:
				self._error = "Error: argument 'select' must be a string or a list of particle IDs"
				return
		
		self.nselectedParticles = len(self.selectedParticles)
		
		# Manage axes
		# -------------------------------------------------------------------
		if type(axes) is not list:
			self._error = "Error: Argument 'axes' must be a list"
			return
		if len(axes)==0:
			self._error = "Error: must define at least one axis."
			return
		self.axes = axes
		self._axesIndex = []
		for axis in axes:
			if axis not in self._properties.keys():
				self._error += "Error: Argument 'axes' has item '"+str(axis)+"' unknown.\n"
				self._error += "       Available axes are: "+(", ".join(sorted(self._properties.keys())))
				return
			self._axesIndex.append( self._properties[axis] ) # axesIndex contains the index in the hdf5 file
		self._type = axes
		for i, axis in enumerate(axes):
			axisi = self._axesIndex[i]
			axisunits = ""
			if axis == "Id":
				self._centers.append( [0, self._h5items[axisi][0,-1]] )
			if axis in ["x" , "y" , "z" ]:
				axisunits = "L_r"
				self._centers.append( [0., self.namelist.Main.sim_length[{"x":0,"y":1,"z":2}[axis]]] )
			if axis in ["px", "py", "pz"]:
				axisunits = "P_r"
				self._centers.append( [-1., 1.] )
			if axis == "Charge":
				axisunits = "Q_r"
				self._centers.append( [-10., 10.] )
			self._log.append( False )
			self._label.append( axis )
			self._units.append( axisunits )
		self._title = "Track particles '"+species+"'"
		self._shape = [0]*len(axes)
		# Hack to work with 1 axis
		if len(axes)==1: self._vunits = self._units[0]
		else: self._vunits = ""
		
		self._rawData = None
		
		# Finish constructor
		self.length = length or self.times[-1]
		self.valid = True
	
	# Method to print info on included probe
	def info(self):
		if not self._validate():
			print(self._error)
		else:
			print("Track particles: species '"+self.species+"' containing "+str(self.nParticles)+" particles")
			if len(self.selectedParticles) != self.nParticles:
				print("                with selection of "+str(len(self.selectedParticles))+" particles")
	
	# get all available tracked species
	def getTrackSpecies(self):
		files = self._glob(self._results_path+"/TrackParticles*.h5")
		species = []
		for file in files:
			species_ = self._re.search("_(.+).h5",self._os.path.basename(file)).groups()[0]
			if species_ not in species: species.append(species_)
		return species
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:
			ntimes = self._h5items[0].len()
		except:
			print("Unable to find tracked particle data in file "+self._file)
			return self._np.array([])
		for item in self._h5items:
			if item.name == "/Times":
				return item.value
		print("Unable to find the list of timesteps in file "+self._file)
		return self._np.array([])
	
	# Make the particles ordered by Id in the file, in case they are not
	def _orderFile( self, fileDisordered, fileOrdered ):
		print("Ordering particles ... (this could take a while)")
		# Copy the disordered file
		from platform import system
		s = system()
		if s in ['Windows']:
			status = self._os.system('xcopy "%s" "%s"' % (fileDisordered, fileOrdered))
		elif s in ['Linux','Darwin']:
			status = self._os.system('cp -fr %s %s' % (fileDisordered, fileOrdered) )
		else:
			status = 0
			try:
				from shutil import copyfile
				copyfile(fileDisordered, fileOrdered)
			except:
				status = 1
		if status != 0:
			raise Exception("New file could not be created: "+str(fileOrdered))
		# Open the file which will become ordered
		print("    Created new file "+fileOrdered)
		f = self._h5py.File(fileOrdered)
		# Get list of properties
		properties = [p.name[1:] for p in f.values() if len(p.shape)==2]
		# For each time
		ntimes, npart = f["Id"].shape
		times = f["Times"]
		A = self._np.zeros((npart,))
		for i in range(ntimes):
			print("    Ordering @ timestep = "+str(times[i]))
			# Get the indices for sorting arrays
			ids = f["Id"][i,:]
			remaining_particles = ids>0
			read_indices  = self._np.nonzero(remaining_particles)
			write_indices = ids[read_indices]-1
			B = self._np.zeros((npart,))
			# Sort arrays
			for property in properties:
				f[property].read_direct (A, source_sel=self._np.s_[i,:])
				B[write_indices] = A[read_indices]
				f[property].write_direct(B, dest_sel  =self._np.s_[i,:])
		# Close files
		f.close()
		print("Ordering succeeded")
	
	# We override the get and getData methods
	def getData(self):
		if not self._validate(): return
		self._prepare1() # prepare the vfactor
		
		if self._rawData is None:
			print("Preparing data ...")
			# create dictionary with info on the axes
			self._rawData = {}
			ntimes = len(self.times)
			for axis in self.axes:
				self._rawData.update({ axis:self._np.zeros((ntimes, self.nselectedParticles)) })
				self._rawData[axis].fill(self._np.nan)
			print("Loading data ...")
			# loop times and fill up the data
			ID = self._np.zeros((self.nParticles,), dtype=self._np.int16)
			B = self._np.zeros((self.nParticles,))
			indices = self.selectedParticles - 1
			for it, ti in enumerate(self._itimes):
				print("     iteration "+str(it+1)+"/"+str(ntimes))
				self._Id.read_direct(ID, source_sel=self._np.s_[ti,:], dest_sel=self._np.s_[:]) # read the particle Ids
				deadParticles = (ID==0).nonzero()
				for i, axis in enumerate(self.axes):
					axisi = self._axesIndex[i]
					self._h5items[axisi].read_direct(B, source_sel=self._np.s_[ti,:], dest_sel=self._np.s_[:])
					B[deadParticles]=self._np.nan
					self._rawData[axis][it, :] = B[indices].squeeze()
			self._rawData.update({ "times":self.times })
			print("... done")
		
		# Multiply by the vfactor
		data = {}
		data.update({ "times":self.times })
		for axis in self.axes:
			data.update({ axis:self._rawData[axis]*self._vfactor })
		
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
		timeSelection = (self.times<=t)*(self.times>=t-self.length)
		times = self.times[timeSelection]
		A     = self._tmpdata[0][timeSelection,:]
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
		timeSelection = (self.times<=t)*(self.times>=t-self.length)
		x = self._tmpdata[0][timeSelection,:]
		y = self._tmpdata[1][timeSelection,:]
		ax.plot(self._xfactor*x, self._yfactor*y, **self.options.plot)
		ax.set_xlabel(self._xlabel)
		ax.set_ylabel(self._ylabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.ymin, ymax=self.options.ymax)
		self._setSomeOptions(ax)
		return 1
