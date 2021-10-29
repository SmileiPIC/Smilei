from .Diagnostic import Diagnostic
from .._Utils import *

class Probe(Diagnostic):
	"""Class for loading a Probe diagnostic"""

	def _init(self, probeNumber=None, field=None, timesteps=None, subset=None, average=None, data_log=False, chunksize=10000000, data_transform=None, **kwargs):

		self._h5probe = []
		self._alltimesteps = []
		self._chunksize = chunksize
		self._subsetinfo = {}
		
		# Search available diags
		diag_numbers, diag_names = self.simulation.getDiags("Probes")
		
		# If no probeNumber, print available probes
		if probeNumber is None:
			if len(diag_numbers)>0:
				error = ["Argument `probeNumber` not provided"]
				error += ["Printing available probes:"]
				error += ["--------------------------"]
				for p in diag_numbers:
					error += [self._info(self._getInfo(p))]
			else:
				error += ["No probes found"]
			raise Exception("\n".join(error))
		
		info = self.simulation.probeInfo(probeNumber)
		self.probeNumber = info["probeNumber"]
		self.probeName   = info["probeName"]
		self._fields = info["fields"]
		
		# Try to get the probe from the hdf5 file
		for path in self._results_path:
			# Open file
			file = path+self._os.sep+"Probes"+str(self.probeNumber)+".h5"
			try:
				self._h5probe.append( self._h5py.File(file, 'r') )
			except Exception as e:
				continue
		
		# If no field, print available fields
		if field is None:
			error = ["Argument `field` not provided"]
			error += ["Printing available fields for probe #"+str(probeNumber)+":"]
			error += ["----------------------------------------"]
			error += [str(", ".join(self._fields))]
			raise Exception("\n".join(error))
		
		# Get available times
		self._dataForTime = {}
		for file in self._h5probe:
			for key, val in file.items():
				try   : self._dataForTime[int(key)] = val
				except Exception as e: break
		self._alltimesteps = self._np.double(sorted(self._dataForTime.keys()))
		if self._alltimesteps.size == 0:
			raise Exception("No timesteps found")
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Parse the `field` argument
		sortedfields = reversed(sorted(self._fields, key = len))
		self.operation = field
		for f in sortedfields:
			i = self._fields.index(f)
			self.operation = self.operation.replace(f,"#"+str(i))
		requested_fields = self._re.findall("#\d+",self.operation)
		if len(requested_fields) == 0:
			raise Exception("Could not find any existing field in `"+field+"`")
		self._fieldn = [ int(f[1:]) for f in requested_fields ] # indexes of the requested fields
		self._fieldn = list(set(self._fieldn))
		self._fieldname = [ self._fields[i] for i in self._fieldn ] # names of the requested fields

		# Check subset
		if subset is None:
			subset = {}
		elif type(subset) is not dict:
			raise Exception("Argument `subset` must be a dictionary")
		
		# Check average
		if average is None:
			average = {}
		elif type(average) is not dict:
			raise Exception("Argument `average` must be a dictionary")

		# Put data_log as object's variable
		self._data_log = data_log
		self._data_transform = data_transform

		# Get the shape of the probe
		self._myinfo = self._getMyInfo()
		self._initialShape = self._myinfo["shape"]
		if self._initialShape.prod()==1:
			self._initialShape = self._np.array([], dtype=int)
		self.numpoints = self._h5probe[0]["positions"].shape[0]
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		self._timesteps = self._alltimesteps
		if timesteps is not None:
			try:
				self._timesteps = self._selectTimesteps(timesteps, self._timesteps)
			except Exception as e:
				raise Exception("Argument `timesteps` must be one or two non-negative integers")

		# Need at least one timestep
		if self._timesteps.size < 1:
			raise Exception("Timesteps not found")

		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values
		self._naxes = self._initialShape.size
		self._finalShape = self._np.copy(self._initialShape)
		self._averages = [False]*self._naxes
		self._selection = [self._np.s_[:]]*self._naxes
		p = []
		self.p_plot = []
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

			# If averaging over this axis
			if label in average:
				if label in subset:
					raise Exception("`subset` not possible on the same axes as `average`")

				self._averages[iaxis] = True

				distances = self._np.sqrt(self._np.sum((centers-centers[0])**2,axis=1))
				self._subsetinfo[label], self._selection[iaxis], self._finalShape[iaxis] \
					= self._selectRange(average[label], distances, label, axisunits, "average")
			
			# Otherwise
			else:
				# If taking a subset of this axis
				if label in subset:
					distances = self._np.sqrt(self._np.sum((centers-centers[0])**2,axis=1))
					self._subsetinfo[label], self._selection[iaxis], self._finalShape[iaxis] \
						= self._selectSubset(subset[label], distances, label, axisunits, "subset")
				
				# If subset has more than 1 point (or no subset), use this axis in the plot
				if type(self._selection[iaxis]) is slice:
					self._type   .append(label)
					self._shape  .append(self._initialShape[iaxis])
					self._centers.append(centers[self._selection[iaxis],:])
					self._label  .append(label)
					self._units  .append(axisunits)
					self._log    .append(False)
					self.p_plot  .append(p[-1])
		
		self._selection = tuple(s if type(s) is slice else slice(s,s+1) for s in self._selection)
		
		# Special case in 1D: we convert the point locations to scalar distances
		if len(self._centers) == 1:
			self._centers[0] = self._np.sqrt(self._np.sum((self._centers[0]-self._centers[0][0])**2,axis=1))
		# Special case in 2D: we have to prepare for pcolormesh instead of imshow
		elif len(self._centers) == 2:
			p1 = self._centers[0] # locations of grid points along first dimension
			d = (p1[1,:] - p1[0,:])/2. # half separation between the points
			p1[0,:] += d
			p1 = self._np.vstack((p1, p1[-1,:]+d)) # add last edges at the end of box
			offset = p1[0,:]
			p1 = self._np.apply_along_axis(lambda x: x-offset, 1, p1) # move points
			p2 = self._centers[1] # locations of grid points along second dimension
			d = (p2[1,:] - p2[0,:])/2. # half separation between the points
			p2[0,:] += d
			p2 = self._np.vstack((p2, p2[-1,:]+d)) # add last edges at the end of box
			offset = p2[0,:]
			p2 = self._np.apply_along_axis(lambda x: x-offset, 1, p2) # move points
			# Trick in a 3D simulation (the probe has to be projected)
			if self._ndim_particles==3:
				# unit vectors in the two dimensions + perpendicular
				u1 = self.p_plot[0] / self._np.linalg.norm(self.p_plot[0])
				u2 = self.p_plot[1] / self._np.linalg.norm(self.p_plot[1])
				# Distances along first direction
				p1[:,0] = self._np.dot(p1, u1)
				p1[:,1:] = 0.
				# Distances along second direction
				p2x = self._np.dot(p2, u1)
				p2[:,1] = self._np.dot(p2, u2)
				p2[:,0] = p2x
				p2[:,2:] = 0.
			# Now p1 and p2 contain edges grid points along the 2 dimensions
			# We have to convert into X and Y 2D arrays (similar to meshgrid)
			X = self._np.zeros((p1.shape[0], p2.shape[0]))
			Y = self._np.zeros((p1.shape[0], p2.shape[0]))
			for i in range(p2.shape[0]):
				X[:,i] = p1[:,0] + (p2[i,0]-p2[0,0])
				Y[:,i] = p1[:,1] + (p2[i,1]-p2[0,1])
			#FOLLOWING LINES CREATE PB IN THE SELECTION OF min/max OF AXES
			#X = self._np.maximum( X, 0.)
			#X = self._np.minimum( X, self._ncels[0]*self._cell_length[0])
			#Y = self._np.maximum( Y, 0.)
			#Y = self._np.minimum( Y, self._ncels[1]*self._cell_length[1])
			self._edges = [X, Y]

		# Prepare the reordering of the points for patches disorder
		tmpShape = self._initialShape
		if self._naxes == 0:
			self._ordering = self._np.array([0.], dtype=int)
		
		else:
			self._ordering = self._np.zeros((self._finalShape.prod(),), dtype=int)-1
			
			# calculate matrix inverse
			if self._naxes==1:
				p  = self._np.sqrt(self._np.sum(self._np.array(p)**2))
				invp = self._np.array(1./p, ndmin=2)
			else:
				if (self._naxes==2 and self._ndim_particles==3):
					pp = self._np.cross(p[0],p[1])
					p.append(pp/self._np.linalg.norm(pp))
					tmpShape = self._np.hstack((tmpShape, 1))
				invp = self._np.linalg.inv(self._np.array(p).transpose())
			
			# calculate ordering
			p0 = self._myinfo["p0"]
			for first, last, npart in ChunkedRange(self.numpoints, chunksize):
				positions = self._h5probe[0]["positions"][first:last,:].T # actual probe points positions
				# Subtract by p0
				for i in range(p0.size):
					positions[i,:] -= p0[i]
				# In 1D convert positions to distances
				if self._naxes==1:
					positions = self._np.sqrt(self._np.sum(positions**2,0))[self._np.newaxis,:]
				# Find the indices of the points
				ijk = (self._np.dot(invp, positions)*(tmpShape-1)[:,self._np.newaxis]).round().astype(int)
				keep = self._np.ones( (npart,), dtype=bool )
				for d,sel in enumerate(self._selection): # keep only points in selection
					start = sel.start or 0
					stop = sel.stop or tmpShape[d]
					step = sel.step or 1
					keep *= (ijk[d] >= start) * (ijk[d] < stop) * ((ijk[d]-start) % step == 0)
				ijk = ijk[:,keep]
				indexInFile = self._np.arange(first, last, dtype=int)[keep]
				# Convert to indices in the selection
				for d,sel in enumerate(self._selection):
					ijk[d] = (ijk[d] - (sel.start or 0)) // (sel.step or 1)
				# Linearize index
				indexInArray = ijk[0]
				for d in range(1,len(self._finalShape)):
					indexInArray = indexInArray*self._finalShape[d] + ijk[d]
				# Store ordering
				self._ordering[indexInArray] = indexInFile
		
		# Build units
		titles = {}
		fieldunits = {}
		unitsForField = {"B":"B_r","E":"E_r","J":"J_r","R":"N_r","P":"V_r*K_r*N_r"}
		self.time_integral = self._myinfo["time_integral"]
		
		for f in self._fieldname:
			i = self._fields.index(f)
			if self.time_integral:
				fieldunits.update({ i:unitsForField[f[0]] + "*T_r" })
				titles    .update({ i:"Time-integrated "+f })
			else:
				fieldunits.update({ i:unitsForField[f[0]] })
				titles    .update({ i:f })
		# Make total units and title
		self._title  = self.operation
		self._vunits = self.operation
		for n in self._fieldn:
			self._title  = self._title .replace("#"+str(n), titles    [n])
			self._vunits = self._vunits.replace("#"+str(n), fieldunits[n])
		self._vunits = self.units._getUnits(self._vunits)

		# Set the directory in case of exporting
		self._exportPrefix = "Probe"+str(probeNumber)+"_"+"".join(self._fieldname)
		self._exportDir = self._setExportDir(self._exportPrefix)

		# Finish constructor
		self.valid = True
		return kwargs

	# destructor
	def __del__(self):
		if hasattr(self, "_h5probe"):
			for file in self._h5probe:
				try:
					file.close()
				except Exception as e:
					pass

	# Method to print info previously obtained with getInfo
	def _info(self, info=None):
		if info is None: info = self._getMyInfo()
		printedInfo = "Probe #%s: "%info["probeNumber"]
		if "dimension" in info:
			printedInfo += str(info["dimension"])+"-dimensional,"+" with fields "+_decode(info["fields"])
			i = 0
			while "p"+str(i) in info:
				printedInfo += "\n\tp"+str(i)+" = "+" ".join(info["p"+str(i)].astype(str).tolist())
				i += 1
			if info["shape"].size>0:
				printedInfo += "\n\tnumber = "+" ".join(info["shape"].astype(str).tolist())
		else:
			printedInfo += "\n\tFile not found or not readable"
		for l in self._subsetinfo:
			printedInfo += "\n\t"+self._subsetinfo[l]
		return printedInfo

	# Method to get info on a given probe
	def _getInfo(self, probeNumber):
		out = {}
		out["probeNumber"] = probeNumber
		for path in self._results_path:
			try:
				file = path+"/Probes"+str(probeNumber)+".h5"
				probe = self._h5py.File(file, 'r')
			except Exception as e:
				continue
			out["dimension"] = probe.attrs["dimension"]
			out["shape"] = self._np.array(probe["number"], dtype=int)
			out["fields"] = probe.attrs["fields"]
			if "time_integral" in probe.attrs:
				out["time_integral"] = probe.attrs["time_integral"]
			else:
				out["time_integral"] = False
			i = 0
			while "p"+str(i) in probe.keys():
				out["p"+str(i)] = self._np.array(probe["p"+str(i)])
				i += 1
			probe.close()
			return out
		raise("Cannot open any file Probes"+str(probeNumber)+".h5")
	
	def _getMyInfo(self):
		return self._getInfo(self.probeNumber)
	
	# get all available fields
	def getFields(self):
		return self._fields
	
	# get the value of x_moved for a requested timestep
	def getXmoved(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# get h5 iteration group
		h5item = self._dataForTime[t]
		# Change units
		factor, _ = self.units._convert("L_r", None)
		return h5item.attrs["x_moved"]*factor if "x_moved" in h5item.attrs else 0.

	# get all available timesteps
	def getAvailableTimesteps(self):
		return self._alltimesteps

	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+t+" not found in this diagnostic")
			return []
		# Get arrays from requested field
		# get data
		C = {}
		op = self.operation
		for n in reversed(self._fieldn): # for each field in operation
			buffer = self._np.zeros((self._ordering.size,), dtype="double")
			for first, last, npart in ChunkedRange(self.numpoints, self._chunksize):
				o = self._ordering - first
				keep = self._np.flatnonzero((o>=0) * (o<npart))
				data = self._dataForTime[t][n,first:last]
				buffer[keep] = data[o[keep]]
			C.update({ n:buffer })
			op = op.replace("#"+str(n), "C["+str(n)+"]")
		# Calculate the operation
		A = eval(op)
		# Reshape array because it is flattened in the file
		A = self._np.reshape(A, self._finalShape)
		# Apply the averaging
		for iaxis in range(self._naxes):
			if self._averages[iaxis]:
				A = self._np.mean(A, axis=iaxis, keepdims=True)
		A = self._np.squeeze(A) # remove averaged axes

		if callable(self._data_transform): A = self._data_transform(A)
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
	def _plotOnAxes_2D_(self, ax, A):
		vmin = self.options.vmin
		vmax = self.options.vmax
		if self.options.vsym:
			if vmin or vmax:
				print("WARNING: vsym set on the same Diagnostic as vmin and/or vmax. Ignoring vmin/vmax.")
		        
			if self.options.vsym is True:
				vmax = self._np.abs(A).max()
			else:
				vmax = self._np.abs(self.options.vsym)

			vmin = -vmax
		self._plot = ax.pcolormesh(self._xfactor*self._edges[0], self._yfactor*self._edges[1], (A),
			vmin = vmin, vmax = vmax, **self.options.image)
		return self._plot
	def _animateOnAxes_2D_(self, ax, A):
		self._plot.set_array( A.flatten() )
		return self._plot
