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
			return
		
		# Get info from the hdf5 files + verifications
		# -------------------------------------------------------------------
		self.species  = species
		# If the first path does not contain the ordered file, we must create it
		orderedfile = self._results_path[0]+self._os.sep+"TrackParticles_"+species+".h5"
		if not self._os.path.isfile(orderedfile):
			disorderedfiles = []
			for path in self._results_path:
				file = path+self._os.sep+"TrackParticlesDisordered_"+species+".h5"
				if not self._os.path.isfile(file):
					self._error = "Missing TrackParticles file in directory "+path
					return
				disorderedfiles += [file]
			self._orderFiles(disorderedfiles, orderedfile)
		# Create arrays to store h5 items
		f = self._h5py.File(orderedfile)
		self._h5items = {}
		self._locationForTime = {}
		for prop in ["Id", "x", "y", "z", "px", "py", "pz"]:
			if prop in f:
				self._h5items[prop] = f[prop]
		# Memorize the locations of timesteps in the files
		for it, t in enumerate(f["Times"]):
			self._locationForTime[t] = it
		self.times = self._np.array(sorted(f["Times"]))
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
		self.nParticles = self._h5items["Id"].shape[1]
		
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
				it = self._locationForTime[time]
				dataset = self._h5items[property]
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
						loc = self._np.flatnonzero(ID>0) # indices of existing particles
						selection[loc] = function( selection[loc], selectionAtTimeT[loc])
					#except:
					#	raise Exception("Error in selector syntax: not understood: "+select[i:parenthesis+1])
					stack.append(selection)
					operation += "stack["+str(len(stack)-1)+"]"
					i = parenthesis+1
				else:
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
				self._centers.append( [0, self._h5items[axis][0,-1]] )
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
	def _orderFiles( self, filesDisordered, fileOrdered ):
		print("Ordering particles ... (this could take a while)")
		# Obtain the list of all times in all disordered files
		time_locations = {}
		for fileIndex, fileD in enumerate(filesDisordered):
			f = self._h5py.File(fileD, "r")
			for t in f["data"].keys():
				try   : time_locations[int(t)] = (fileIndex, t)
				except: pass
			f.close()
		times = sorted(time_locations.keys())
		# Open the last file and get the number of particles from each MPI
		last_file_index, tname = time_locations[times[-1]]
		f = self._h5py.File(filesDisordered[last_file_index], "r")
		number_of_particles = (f["data"][tname]["latest_IDs"].value % (2**32)).astype('uint32')
		# Calculate the offset that each MPI needs
		offset = self._np.cumsum(number_of_particles)
		total_number_of_particles = offset[-1]
		offset = self._np.roll(offset, 1)
		offset[0] = 0
		# Make new (ordered) file
		f0 = self._h5py.File(fileOrdered, "w")
		# Make new datasets
		properties = {"id":"Id", "position/x":"x", "position/y":"y", "position/z":"z",
		              "momentum/x":"px", "momentum/y":"py", "momentum/z":"pz"}
		for k, name in properties.items():
			try   : f0.create_dataset(name, (len(times), total_number_of_particles), f["data"][tname]["particles"][self.species][k].dtype, fillvalue=self._np.nan)
			except: pass
		f.close()
		# Loop times and fill arrays
		for it, t in enumerate(times):
			print("    Ordering @ timestep = "+str(t))
			file_index, tname = time_locations[t]
			f = self._h5py.File(filesDisordered[file_index], "r")
			group = f["data"][tname]["particles"][self.species]
			if group["id"].size == 0: continue
			# Get the Ids and find where they should be stored in the final file
			locs = group["id"].value % 2**32 + offset[ group["id"].value>>32 ] -1
			# Loop datasets and order them
			for k, name in properties.items():
				if k not in group: continue
				disordered = group[k].value
				ordered = self._np.zeros((total_number_of_particles, ), dtype=disordered.dtype)
				ordered[locs] = disordered
				f0[name].write_direct(ordered, dest_sel=self._np.s_[it,:])
			f.close()
		# Create the "Times" dataset
		f0.create_dataset("Times", data=times)
		# Close file
		f0.close()
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
				self._rawData[axis] = self._np.zeros((ntimes, self.nselectedParticles))
				self._rawData[axis].fill(self._np.nan)
			print("Loading data ...")
			# loop times and fill up the data
			ID = self._np.zeros((self.nParticles,), dtype=self._np.int16)
			B = self._np.zeros((self.nParticles,))
			indices = self.selectedParticles - 1
			for it, time in enumerate(self.times):
				print("     iteration "+str(it+1)+"/"+str(ntimes))
				timeIndex = self._locationForTime[time]
				self._h5items["Id"].read_direct(ID, source_sel=self._np.s_[timeIndex,:]) # read the particle Ids
				deadParticles = (ID==0).nonzero()
				for axis in self.axes:
					self._h5items[axis].read_direct(B, source_sel=self._np.s_[timeIndex,:])
					B[deadParticles]=self._np.nan
					self._rawData[axis][it, :] = B[indices].squeeze()
			# Add the lineBreaks array which indicates where lines are broken (e.g. loop around the box)
			self._rawData['brokenLine'] = self._np.zeros((self.nselectedParticles,), dtype=bool)
			self._rawData['lineBreaks'] = {}
			dt = self._np.diff(self.times)*self.timestep
			for axis in ["x","y","z"]:
				if axis in self.axes:
					dudt = self._np.diff(self._rawData[axis],axis=0)
					for i in range(dudt.shape[1]): dudt[:,i] /= dt
					self._rawData['brokenLine'] += self._np.abs(dudt).max(axis=0) > 1.
					broken_particles = self._np.flatnonzero(self._rawData['brokenLine'])
					for broken_particle in broken_particles:
						broken_times = list(self._np.flatnonzero(self._np.abs(dudt[:,broken_particle]) > 1.)+1)
						if broken_particle in self._rawData['lineBreaks'].keys():
							self._rawData['lineBreaks'][broken_particle] += broken_times
						else:
							self._rawData['lineBreaks'][broken_particle] = broken_times
			# Add the times array
			self._rawData["times"] = self.times
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
		tmin = t-self.length
		tmax = t
		timeSelection = (self.times<=tmax)*(self.times>=tmin)
		selected_times = self._np.flatnonzero(timeSelection)
		itmin = selected_times[0]
		itmax = selected_times[-1]
		# Plot first the non-broken lines
		x = self._tmpdata[0][timeSelection,:][:,~self._rawData["brokenLine"]]
		y = self._tmpdata[1][timeSelection,:][:,~self._rawData["brokenLine"]]
		ax.plot(self._xfactor*x, self._yfactor*y, **self.options.plot)
		# Then plot the broken lines
		ax.hold("on")
		for line, breaks in self._rawData['lineBreaks'].items():
			x = self._tmpdata[0][:, line]
			y = self._tmpdata[1][:, line]
			prevline = None
			for ibrk in range(len(breaks)):
				if breaks[ibrk] <= itmin: continue
				iti = itmin
				if ibrk>0: iti = max(itmin, breaks[ibrk-1])
				itf = min( itmax, breaks[ibrk] )
				if prevline:
					ax.plot(self._xfactor*x[iti:itf], self._yfactor*y[iti:itf], color=prevline.get_color(), **self.options.plot)
				else:
					prevline, = ax.plot(self._xfactor*x[iti:itf], self._yfactor*y[iti:itf], **self.options.plot)
				if breaks[ibrk] > itmax: break
		ax.hold("off")
		# Add labels and options
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
			npoints = self._h5items['Id'].shape[1]
			f.write('		<DataItem Name="Zeroes" ItemType="Uniform" NumberType="Float" Dimensions="'+str(npoints)+'" Format="XML">'+"0. "*npoints+'</DataItem>\n')
			f.write('		<Grid GridType="Collection" CollectionType="Temporal">\n')
			file = self._h5items['Id'].file
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
	
	
	
	# Convert data to VTK format
	def toVTK(self, numberOfPieces=1):
		if not self._validate(): return
		
		if self._ndim!=3:
			print "Cannot export tracked particles of a "+str(self._ndim)+"D simulation to VTK"
			return
		
		self._mkdir(self._exportDir)
		fileprefix = self._exportDir + self._exportPrefix
		
		ntimes = len(self.times)
		
		vtk = VTKfile()
		
		# If 3D simulation, then do a 3D plot
		if self._ndim == 3:
			if "x" not in self.axes or "y" not in self.axes or "z" not in self.axes:
				print("Error exporting tracked particles to VTK: axes 'x', 'y' and 'z' are required")
				return
			data = self.getData()
			pcoords = self._np.stack((data["x"],data["y"],data["z"])).transpose()
			npoints, nt, nd = pcoords.shape
			pcoords = self._np.reshape(pcoords, (npoints*nt, nd))
			pcoords = self._np.ascontiguousarray(pcoords, dtype='float32')
			pcoords = vtk.Array(pcoords, "")
			connectivity = self._np.ascontiguousarray([[nt]+[nt*i+j for j in range(nt)] for i in range(npoints)])
			
			attributes = []
			for ax in self.axes:
				if ax!="x" and  ax!="y" and  ax!="z":
					attributes += [vtk.Array(self._np.ascontiguousarray(data[ax].flatten(),'float32'),ax)]
			
			vtk.WriteLines(pcoords, connectivity, attributes, fileprefix+".vtk")
			print("Successfully exported tracked particles to VTK, folder='"+self._exportDir)
		
