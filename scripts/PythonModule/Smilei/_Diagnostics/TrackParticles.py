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
				self._error = "No tracked particles files found"
			return None
		
		# Get info from the hdf5 files + verifications
		# -------------------------------------------------------------------
		self.species  = species
		translateProperties = {"Id":"Id", "x":"Position-0", "y":"Position-1", "z":"Position-2",
			"px":"Momentum-0", "py":"Momentum-1", "pz":"Momentum-2"}
		self._h5items = None
		for pathNumber, path in enumerate(self._results_path):
			# Order the particles by ID
			file = path+self._os.sep+"TrackParticles_"+species+".h5"
			try:
				f = self._h5py.File(file, 'r')
			except:
				self._orderFile( path+self._os.sep+"TrackParticlesDisordered_"+species+".h5", file )
				f = self._h5py.File(file, 'r')
			# Create arrays to store h5 items
			if self._h5items is None:
				self._h5items = {}
				self._locationForTime = {}
				for prop, name in translateProperties.items():
					if name in f.keys():
						self._h5items[prop] = []
			# Memorize the locations of each property (x, y, etc) in the files
			for prop, val in self._h5items.items():
				val.append( f[translateProperties[prop]] )
			# Memorize the locations of timesteps in the files
			for it, t in enumerate(f["Times"]):
				self._locationForTime[t] = [pathNumber, it]
		self.times = self._np.array(sorted(self._locationForTime.keys()))
		self._times = self.times[:]
		
		# Get available times in the hdf5 file
		if self.times.size == 0:
			self._error = "No tracked particles found"
			return
		# If specific timesteps requested, narrow the selection
		if timesteps is not None:
			try:
				ts = self._np.array(self._np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self.times = self.times[ self._np.nonzero((self.times>=ts[0]) * (self.times<=ts[1]))[0] ]
				elif ts.size==1:
					# get nearest time
					self.times = self._np.array(self.times[ self._np.array([(self._np.abs(self.times-ts)).argmin()]) ])
				else:
					raise
			except:
				self._error = "Argument `timesteps` must be one or two non-negative integers"
				return
		# Need at least one timestep
		if self.times.size < 1:
			self._error = "Timesteps not found"
			return
		
		# Get number of particles
		self.nParticles = self._h5items["Id"][0].shape[1]
		
		# Select particles
		# -------------------------------------------------------------------
		if type(select) is str:
			# Define a function that finds the next closing character in a string
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
			# Define a function that gets some requested data
			def getData(property, time, buffer=None):
				pathNumber, it = self._locationForTime[time]
				dataset = self._h5items[property][pathNumber]
				if buffer is None: buffer = self._np.zeros((self.nParticles,))
				dataset.read_direct(buffer, source_sel=self._np.s_[it,:])
				return buffer
			# Start reading the selector
			i = 0
			stack = []
			operation = ""
			while i < len(select):
				if i+4<len(select) and select[i:i+4] in ["any(","all("]:
					if select[i:i+4] == "any(": function = self._np.logical_or
					if select[i:i+4] == "all(": function = self._np.logical_and
					comma = findClosingCharacter(select, ",", i+4)
					parenthesis = findClosingCharacter(select, ")", comma+1)
					timeSelector = select[i+4:comma]
					try:
						s = self._re.sub(r"\bt\b","self._times",timeSelector)
						times = self._times[eval(s)]
					except:
						raise Exception("Error in selector syntax: time selector not understood in "+select[i:i+3]+"()")
					try:
						particleSelector = select[comma+1:parenthesis]
						for prop in self._h5items.keys():
							particleSelector = self._re.sub(r"\b"+prop+r"\b", "getData('"+prop+"',time)", particleSelector)
					except:
						raise Exception("Error in selector syntax: not understood: "+select[i:parenthesis+1])
					if select[i:i+4] == "any(": selection = self._np.array([False]*self.nParticles)
					if select[i:i+4] == "all(": selection = self._np.array([True]*self.nParticles)
					#try:
					ID = self._np.zeros((self.nParticles,), dtype=self._np.int32)
					for time in times:
						selectionAtTimeT = eval(particleSelector) # array of True or False
						getData("Id", time, ID)
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
		for axis in axes:
			if axis not in self._h5items.keys():
				self._error += "Error: Argument 'axes' has item '"+str(axis)+"' unknown.\n"
				self._error += "       Available axes are: "+(", ".join(sorted(self._h5items.keys())))
				return
		self._type = axes
		for axis in axes:
			axisunits = ""
			if axis == "Id":
				self._centers.append( [0, self._h5items[axis][0][0,-1]] )
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
		
		# Set the directory in case of exporting
		self._exportPrefix = "TrackParticles_"+self.species+"_"+"".join(self.axes)
		self._exportDir = self._setExportDir(self._exportPrefix)
		
		# Finish constructor
		self.length = length or self.times[-1]
		self.valid = True
	
	# Method to print info on included probe
	def _info(self):
		info = "Track particles: species '"+self.species+"' containing "+str(self.nParticles)+" particles"
		if len(self.selectedParticles) != self.nParticles:
			info += "\n                with selection of "+str(len(self.selectedParticles))+" particles"
		return info
	
	# get all available tracked species
	def getTrackSpecies(self):
		for path in self._results_path:
			files = self._glob(path+self._os.sep+"TrackParticles*.h5")
			species_here = [self._re.search("_(.+).h5",self._os.path.basename(file)).groups()[0] for file in files]
			try   : species = [ s for s in species if s in species_here ]
			except: species = species_here
		return species
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		return self._times
	
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
			for it, time in enumerate(self.times):
				print("     iteration "+str(it+1)+"/"+str(ntimes))
				pathNumber, timeIndexInPath = self._locationForTime[time]
				self._h5items["Id"][pathNumber].read_direct(ID, source_sel=self._np.s_[timeIndexInPath,:]) # read the particle Ids
				deadParticles = (ID==0).nonzero()
				for axis in self.axes:
					self._h5items[axis][pathNumber].read_direct(B, source_sel=self._np.s_[timeIndexInPath,:])
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
	
	# Convert to XDMF format for ParaView
	def toXDMF(self):
		
		self._mkdir(self._exportDir)
		
		# Make the XDMF for usual time collections
		with open(self._exportDir+sep+"TrackParticles_"+str(self.species)+".xmf",'w') as f:
			f.write('<?xml version="1.0" ?>\n')
			f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
			f.write('<Xdmf Version="3.0">\n')
			f.write('	<Domain>\n')
			npoints = self._h5items['Id'][0].shape[1]
			f.write('		<DataItem Name="Zeroes" ItemType="Uniform" NumberType="Float" Dimensions="'+str(npoints)+'" Format="XML">'+"0. "*npoints+'</DataItem>\n')
			f.write('		<Grid GridType="Collection" CollectionType="Temporal">\n')
			nfiles = len(self._h5items['Id'])
			for ifile in range(nfiles):
				file = self._h5items['Id'][ifile].file
				filename = self._os.path.abspath(file.filename)
				ntimes = len(file['Times'])
				for itime in range(ntimes):
					selection = "%d,%d:%d,%d:%d,%d:%d,%d" % (itime,0, 1,1, 1,npoints, 1,npoints)
					f.write('			<Grid Name="Timestep_'+str(itime)+'" GridType="Uniform">\n')
					f.write('				<Time Value="'+str(file['Times'][itime])+'" />\n')
					f.write('				<Topology TopologyType="Polyvertex" NumberOfElements="'+str(npoints)+'"/>\n')
					f.write('				<Geometry Name="geometry" GeometryType="VXVYVZ">\n')
					f.write('					<DataItem ItemType="Uniform" NumberType="Float" Dimensions="'+str(npoints)+'" Precision="8" Format="HDF">'+filename+':/Position-0|'+selection+'</DataItem>\n')
					if self._ndim < 2:
						f.write('					<DataItem ItemType="Uniform" NumberType="Float" Dimensions="'+str(npoints)+'" Precision="8" Format="XML" Reference="XML">/Xdmf/Domain/DataItem[@Name="Zeroes"]</DataItem>\n')
					else:
						f.write('					<DataItem ItemType="Uniform" NumberType="Float" Dimensions="'+str(npoints)+'" Precision="8" Format="HDF">'+filename+':/Position-1|'+selection+'</DataItem>\n')
					if self._ndim < 3:
						f.write('					<DataItem ItemType="Uniform" NumberType="Float" Dimensions="'+str(npoints)+'" Precision="8" Format="XML" Reference="XML">/Xdmf/Domain/DataItem[@Name="Zeroes"]</DataItem>\n')
					else:
						f.write('					<DataItem ItemType="Uniform" NumberType="Float" Dimensions="'+str(npoints)+'" Precision="8" Format="HDF">'+filename+':/Position-2|'+selection+'</DataItem>\n')
					f.write('				</Geometry>\n')
					f.write('				<Attribute Name="Px" Center="Node" AttributeType="Scalar">\n')
					f.write('					<DataItem ItemType="Uniform" NumberType="Float" Precision="8" Dimensions="'+str(npoints)+'" Format="HDF">'+filename+':/Momentum-0|'+selection+'</DataItem>\n')
					f.write('				</Attribute>\n')
					f.write('				<Attribute Name="Py" Center="Node" AttributeType="Scalar">\n')
					f.write('					<DataItem ItemType="Uniform" NumberType="Float" Precision="8" Dimensions="'+str(npoints)+'" Format="HDF">'+filename+':/Momentum-1|'+selection+'</DataItem>\n')
					f.write('				</Attribute>\n')
					f.write('				<Attribute Name="Pz" Center="Node" AttributeType="Scalar">\n')
					f.write('					<DataItem ItemType="Uniform" NumberType="Float" Precision="8" Dimensions="'+str(npoints)+'" Format="HDF">'+filename+':/Momentum-2|'+selection+'</DataItem>\n')
					f.write('				</Attribute>\n')
					f.write('			</Grid>\n')
			f.write('		</Grid>\n')
			f.write('	</Domain>\n')
			f.write('</Xdmf>\n')