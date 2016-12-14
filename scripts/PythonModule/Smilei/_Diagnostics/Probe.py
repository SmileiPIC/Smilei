from .Diagnostic import Diagnostic
from .._Utils import *

# -------------------------------------------------------------------
# Class for probe diagnostics
# -------------------------------------------------------------------
class Probe(Diagnostic):
	# This is the constructor, which creates the object
	def _init(self, probeNumber=None, field=None, timesteps=None, slice=None, data_log=False, **kwargs):
		
		self._h5probe = []
		self._times = []
		
		# Get the available probes
		for path in self._results_path:
			files = self._glob(path+self._os.sep+"Probes*.h5")
			probes = [self._re.findall(r"Probes([0-9]+)[.]h5$",file)[0] for file in files]
			try   : self._probes = [p for p in probes if p in self._probes]
			except: self._probes = probes
		
		# If no probeNumber, print available probes
		if probeNumber is None:
			if len(self._probes)>0:
				self._error += "Printing available probes:\n"
				self._error += "--------------------------\n"
				for p in self._probes:
					self._error += self._info(self._getInfo(p))
			else:
				self._error += "No probes found"
			return
		
		# Try to get the probe from the hdf5 file
		self.probeNumber = probeNumber
		for path in self._results_path:
			# Open file
			file = path+self._os.sep+"Probes"+str(self.probeNumber)+".h5"
			try:
				self._h5probe.append( self._h5py.File(file, 'r') )
			except:
				self._error = "Error opening probe #"+str(probeNumber)+" in path '"+path+"'"
				return
			# Verify that this file is compatible with the previous ones
			try:
				for key, val in verifications.items():
					if self._h5probe[-1][key].value != val:
						self._error = "Probe #"+str(probeNumber)+" in path '"+path+"' is incompatible with the other ones"
						return
			except:
				verifications = {"number":self._h5probe[-1]["number"].value}
				npoints = self._h5probe[-1]["number"].size
				if self._h5probe[-1]["number"].value.prod() > 1:
					npoints += 1
				for i in range(npoints):
					verifications["p"+str(i)] = self._h5probe[-1]["p"+str(i)].value
		
		# Get available times
		self._dataForTime = {}
		for file in self._h5probe:
			for key, val in file.items():
				try   : self._dataForTime[int(key)] = val
				except: break
		self._times = self._np.double(self._dataForTime.keys())
		if self._times.size == 0:
			self._error = "No timesteps found"
			return
		
		# Extract available fields
		fields = self.getFields()
		if len(fields) == 0:
			self._error = "No fields found for probe #"+probeNumber
			return
		# If no field, print available fields
		if field is None:
			self._error += "Printing available fields for probe #"+str(probeNumber)+":\n"
			self._error += "----------------------------------------\n"
			self._error += ", ".join(fields)+"\n"
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
			self._error += "Could not find any existing field in `"+field+"`"
			return
		self._fieldn = [ int(f[1:]) for f in requested_fields ] # indexes of the requested fields
		self._fieldn = list(set(self._fieldn))
		self._fieldname = [ fields[i] for i in self._fieldn ] # names of the requested fields
		
		# Check slice is a dict
		if slice is not None  and  type(slice) is not dict:
			self._error += "Argument `slice` must be a dictionary"
			return
		# Make slice a dictionary
		if slice is None: slice = {}
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Get the shape of the probe
		self._myinfo = self._getMyInfo()
		self._initialShape = self._myinfo["shape"]
		if self._initialShape.prod()==1: self._initialShape=self._np.array([])
		self.numpoints = self._h5probe[0]["positions"].shape[0]
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		self.times = self._times
		if timesteps is not None:
			try:
				self.times = self._selectTimesteps(timesteps, self.times)
			except:
				self._error = "Argument `timesteps` must be one or two non-negative integers"
				return
		
		# Need at least one timestep
		if self.times.size < 1:
			self._error = "Timesteps not found"
			return
		
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values
		self._naxes = self._initialShape.size
		self._finalShape = []
		self._sliceinfo = {}
		self._slices = []
		self._selection = ()
		p = []
		for iaxis in range(self._naxes):
		
			# calculate grid points locations
			p0 = self._myinfo["p0"            ] # reference point
			pi = self._myinfo["p"+str(iaxis+1)] # end point of this axis
			p.append( pi-p0 )
			centers = self._np.zeros((self._initialShape[iaxis],p0.size))
			for i in range(p0.size):
				centers[:,i] = self._np.linspace(p0[i],pi[i],self._initialShape[iaxis])
			
			label = {0:"axis1", 1:"axis2", 2:"axis3"}[iaxis]
			axisunits = "L_r"
			
			if label in slice:
				
				self._slices.append(True)
				
				# if slice is "all", then all the axis has to be summed
				if slice[label] == "all":
					self._selection += ( self._np.s_[:], )
					self._finalShape .append( self._initialShape[iaxis] )
				
				# Otherwise, get the slice from the argument `slice`
				else:
					try:
						s = self._np.double(slice[label])
						if s.size>2 or s.size<1: raise
					except:
						self._error += "Slice along axis "+label+" should be one or two floats"
						return
					indices = self._np.arange(self._initialShape[iaxis])
					if s.size==1:
						indices = self._np.array([(self._np.abs(indices-s)).argmin()])
					elif s.size==2:
						indices = self._np.nonzero( (indices>=s[0]) * (indices<=s[1]) )[0]
					if indices.size == 0:
						self._error += "Slice along "+label+" is out of range"
						return
					if indices.size == 1:
						self._sliceinfo.update({ label:"Sliced at "+label+" = "+str(indices[0]) })
						self._selection += ( self._np.s_[indices[0]], )
						self._finalShape .append( 1 )
					else:
						self._sliceinfo.update({ label:"Sliced for "+label+" from "+str(indices[0])+" to "+str(indices[-1]) })
						self._selection += ( self._np.s_[indices[0]:indices[-1]], )
						self._finalShape.append( indices[-1] - indices[0] )
			else:
				self._selection += ( self._np.s_[:], )
				self._finalShape.append( self._initialShape[iaxis] )
				self._slices .append(False)
				self._type   .append(label)
				self._shape  .append(self._initialShape[iaxis])
				self._centers.append(centers)
				self._label  .append(label)
				self._units  .append(axisunits)
				self._log    .append(False)
		
		# Special case in 1D: we convert the point locations to scalar distances
		if len(self._centers) == 1:
			self._centers[0] = self._np.sqrt(self._np.sum((self._centers[0]-self._centers[0][0])**2,axis=1))
			self._centers[0] = self._np.maximum( self._centers[0], 0.)
			self._centers[0] = self._np.minimum( self._centers[0], self._ncels[0]*self._cell_length[0])
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
			X = self._np.maximum( X, 0.)
			X = self._np.minimum( X, self._ncels[0]*self._cell_length[0])
			Y = self._np.maximum( Y, 0.)
			Y = self._np.minimum( Y, self._ncels[1]*self._cell_length[1])
			self._edges = [X, Y]
			self._label = ["x", "y"]
			self._units = [axisunits, axisunits]
		
		# Prepare the reordering of the points for patches disorder
		positions = self._h5probe[0]["positions"].value # actual probe points positions
		self._ordering = None
		if self._naxes>0:
			p = self._np.array(p) # matrix of the probe generating vectors
			# Subtract by p0
			p0 = self._myinfo["p0"]
			for i in range(p0.size):
				positions[:,i] -= p0[i]
			# If 1D probe, convert positions to distances
			if self._naxes==1:
				p  = self._np.sqrt(self._np.sum(p**2))
				invp = self._np.array(1./p, ndmin=1)
				positions = self._np.sqrt(self._np.sum(positions**2,1))
			# If 2D or 3D probe, must calculate matrix inverse
			else:
				invp = self._np.linalg.inv(p.transpose())
			# Make the ordering vector
			self._ordering = self._np.zeros((self._initialShape.prod(),), dtype=int)-1
			for indexInFile in range(positions.shape[0]):
				posInFile = positions[indexInFile]
				ijk = self._np.dot(invp, posInFile)*(self._initialShape-1) # find the indices of the point
				indexInArray = ijk[0]
				for l in range(1,len(ijk)): indexInArray = indexInArray*self._initialShape[l]+ijk[l] # linearized index
				indexInArray = int(round(indexInArray))
				#print indexInFile, indexInArray, posInFile, ijk
				try   : self._ordering[indexInArray] = indexInFile
				except: pass
		
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
		
		self._buffer = self._np.zeros((self.numpoints,), dtype="double")
		
		# Finish constructor
		self.valid = True
	
	# destructor
	def __del__(self):
		for file in self._h5probe:
			file.close()
	
	# Method to print info previously obtained with getInfo
	def _info(self, info=None):
		if info is None: info = self._getMyInfo()
		printedInfo = "Probe #"+str(info["probeNumber"])+": "+str(info["dimension"])+"-dimensional,"+" with fields "+str(info["fields"])+"\n"
		i = 0
		while "p"+str(i) in info:
			printedInfo += "p"+str(i)+" = "+" ".join(info["p"+str(i)].astype(str).tolist())+"\n"
			i += 1
		if info["shape"].size>0:
			printedInfo += "number = "+" ".join(info["shape"].astype(str).tolist())+"\n"
		return printedInfo
	
	# Method to get info on a given probe
	def _getInfo(self, probeNumber):
		try:
			file = self._results_path[0]+"/Probes"+str(probeNumber)+".h5"
			probe = self._h5py.File(file, 'r')
		except:
			print("Cannot open file "+file)
			return {}
		out = {}
		out.update({"probeNumber":probeNumber, "dimension":probe.attrs["dimension"],
			"shape":self._np.array(probe["number"]),"fields":probe.attrs["fields"] })
		i = 0
		while "p"+str(i) in probe.keys():
			out.update({ "p"+str(i):self._np.array(probe["p"+str(i)]) })
			i += 1
		probe.close()
		return out
	def _getMyInfo(self):
		return self._getInfo(self.probeNumber)
	
	# get all available probes
	def getProbes(self):
		return self._probes
	
	# get all available fields
	def getFields(self):
		for file in self._h5probe:
			fields_here = str(file.attrs["fields"]).split(",")
			try   : fields = [f for f in fields_here if f in fields]
			except: fields = fields_here
		return fields
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		return self._times
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self.times:
			print("Timestep "+t+" not found in this diagnostic")
			return []
		# Get arrays from requested field
		# get data
		C = {}
		op = self.operation
		for n in reversed(self._fieldn): # for each field in operation
			self._dataForTime[t].read_direct(self._buffer, source_sel=self._np.s_[n,:])
			C.update({ n:self._buffer })
			op = op.replace("#"+str(n), "C["+str(n)+"]")
		# Calculate the operation
		A = eval(op)
		# Reorder probes for patch disorder
		if self._ordering is not None: A = A[self._ordering]
		# Reshape array because it is flattened in the file
		A = self._np.reshape(A, self._initialShape)
		# Extract the selection
		A = self._np.reshape(A[self._selection], self._finalShape)
		# Apply the slicing
		for iaxis in range(self._naxes):
			if self._slices[iaxis]:
				A = self._np.mean(A, axis=iaxis, keepdims=True) # average over the slice
		A = self._np.squeeze(A) # remove sliced axes
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
	
	# We override _prepare4
	def _prepare4(self):
		# If 2D plot, we remove kwargs that are not supported by pcolormesh
		if self.dim == 2:
			authorizedKwargs = ["cmap"]
			newoptionsimage = {}
			for kwarg in self.options.image.keys():
				if kwarg in authorizedKwargs: newoptionsimage[kwarg]=self.options.image[kwarg]
			self.options.image = newoptionsimage
	
	# Overloading a plotting function in order to use pcolormesh instead of imshow
	def _animateOnAxes_2D_(self, ax, A):
		im = ax.pcolormesh(self._xfactor*self._edges[0], self._yfactor*self._edges[1], (A),
			vmin = self.options.vmin, vmax = self.options.vmax, **self.options.image)
		return im
