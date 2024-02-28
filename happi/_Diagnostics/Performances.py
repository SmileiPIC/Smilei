from .Diagnostic import Diagnostic
from .._Utils import *


# Method to create a matrix containing the hindex of a 2D Hilbert curve
def HilbertCurveMatrix2D(m1, m2=None, oversize=0):
	import numpy as np
	if m2 is None: m2=m1
	o = oversize
	m = min(m1, m2)
	M = abs(m1-m2)
	A = np.empty((2**m2+2*o, 2**m1+2*o), dtype="uint32")
	A[o,o] = 0
	# For each Hilbert iteration
	for i in range(0,m):
		I = 2**i
		J = 2*I
		Npoints = I**2
		K1 = I + o
		K2 = J + o
		# Second quadrant is a transpose of the first
		A[o:K1 , K1:K2] = np.transpose(A[o:K1,o:K1])
		A[o:K1 , K1:K2] += Npoints
		Npoints *= 4
		# Second half is the same as the first half, but flipped
		A[K1:K2, o:K2] = (Npoints-1) - np.flipud(A[o:K1, o:K2])
		# For next iteration, the whole matrix is transposed
		A[o:K2, o:K2] = np.transpose(A[o:K2, o:K2])
	# Duplicate the Hilbert curve several times if non-square
	if m2>m1:
		A[o:K2, o:K2] = np.transpose(A[o:K2, o:K2])
		for j in range(1, 2**M):
			A[(j*J+o):((j+1)*J+o), o:K2] = A[o:K2, o:K2] + (j*Npoints)
	elif m1>m2:
		for j in range(1, 2**M):
			A[o:K2, (j*J+o):((j+1)*J+o)] = A[o:K2, o:K2] + (j*Npoints)
	return A

# Method to create a matrix containing the hindex of a 2D linYX curve
def LinYXCurveMatrix2D(n, oversize=0):
	import numpy as np
	o = oversize
	A = np.zeros((n[1]+2*o, n[0]+2*o), dtype="uint32")
	A[o:-o, o:-o] = np.arange(n[0]*n[1]).reshape(n[1],n[0])
	return A

# Method to create a matrix containing the hindex of a 2D linXY curve
def LinXYCurveMatrix2D(n, oversize=0):
	import numpy as np
	o = oversize
	A = np.zeros((n[1]+2*o, n[0]+2*o), dtype="uint32")
	A[o:-o, o:-o] = np.arange(n[0]*n[1]).reshape(n[0],n[1]).T
	return A

# Method to partition a matrix depending on a list of values
def PartitionMatrix( matrix, listOfValues, oversize=0 ):
	import numpy as np
	partitioned = np.empty(matrix.shape, dtype="uint32")
	max_index = 0
	for v in listOfValues:
		partitioned[matrix >= v] = max_index
		max_index += 1
	if oversize>0:
		partitioned[ :oversize,:] = np.uint32(-1)
		partitioned[:, :oversize] = np.uint32(-1)
		partitioned[-oversize:,:] = np.uint32(-1)
		partitioned[:,-oversize:] = np.uint32(-1)
	return partitioned


class Performances(Diagnostic):
	"""Class for loading a Performances diagnostic"""

	def _init(self, raw=None, map=None, histogram=None, timesteps=None, data_log=False, data_transform=None, species=None, cumulative=True, **kwargs):

		info = self.simulation.performanceInfo()
		self._availableQuantities_uint   = info["quantities_uint"]
		self._availableQuantities_double = info["quantities_double"]
		self.patch_arrangement = info["patch_arrangement"]
		
		# Open the file(s) and load the data
		self._h5items = {}
		for path in self._results_path:
			file = path+self._os.sep+'Performances.h5'
			try:
				f = self._h5py.File(file, 'r')
			except Exception as e:
				continue
			self._h5items.update( dict(f) )
		if not self._h5items:
			raise Exception("Could not open any file Performances.h5")
		# Converted to ordered list
		self._h5items = sorted(self._h5items.values(), key=lambda x:int(x.name[1:]))

		nargs = (raw is not None) + (map is not None) + (histogram is not None)

		if nargs>1:
			raise Exception("Choose only one of `raw`, `map` or `histogram`")

		if nargs == 0:
			raise Exception(
				"Must define raw='quantity', map='quantity' or histogram=['quantity',min,max,nsteps]\n"\
				+"Available quantities: "+", ".join([str(q) for q in self.getAvailableQuantities()])
			)
		
		# Get available times
		self._timesteps = self.getAvailableTimesteps()
		if self._timesteps.size == 0:
			raise Exception("No timesteps found")
		
		# Get the number of procs of the data
		self._nprocs = self._h5items[0]["quantities_uint"].shape[1]
		for item in self._h5items:
			if item["quantities_uint"].shape[1] != self._nprocs:
				raise Exception("Incompatible simulations")
		
		# Get the shape of patches
		self._number_of_patches = self.simulation.namelist.Main.number_of_patches
		self._tot_number_of_patches = self._np.prod( self._number_of_patches )
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Parse the `map` or `histogram` arguments
		if raw is not None:
			if type(raw) is not str:
				raise Exception("Argument `raw` must be a string")
			self.operation = raw
			self._mode = "raw"

		elif map is not None:
			if type(map) is not str:
				raise Exception("Argument `map` must be a string")
			if self._ndim_fields > 2:
				raise Exception("Argument `map` not available in "+str(self._ndim_fields)+"D")
			self.operation = map
			self._mode = "map"
			self._m = [int(self._np.log2(n)) for n in self._number_of_patches]

		elif histogram is not None:
			if type(histogram) is not list or len(histogram) != 4:
				raise Exception("Argument `histogram` must be a list with 4 elements")
			if type(histogram[0]) is not str:
				raise Exception("Argument `histogram` must be a list with first element being a string")
			self.operation = histogram[0]
			try:
				histogram_min    = float(histogram[1])
				histogram_max    = float(histogram[2])
				histogram_nsteps = int  (histogram[3])
			except Exception as e:
				raise Exception("Diagnostic not loaded: argument `histogram` must be a list like ['quantity',min,max,nsteps]")
			self._mode = "hist"

		# Parse the operation
		self._operation = self.operation
		self._operationunits = self.operation
		index_in_output = 0
		self._quantities_uint = []
		used_quantities = []
		for index_in_file, q in enumerate(self._availableQuantities_uint):
			if self._re.search(r"\b%s\b"%q,self._operation):
				self._operation = self._re.sub(r"\b%s\b"%q,"C["+str(index_in_output)+"]",self._operation)
				units = "seconds" if q.startswith("timer") else "1"
				self._operationunits = self._operationunits.replace(q, units)
				self._quantities_uint.append([index_in_file, q])
				used_quantities.append( q )
				index_in_output += 1
		self._quantities_double = []
		for index_in_file, q in enumerate(self._availableQuantities_double):
			if self._re.search(r"\b%s\b"%q,self._operation):
				self._operation = self._re.sub(r"\b%s\b"%q,"C["+str(index_in_output)+"]",self._operation)
				units = "seconds" if q.startswith("timer") else "1"
				self._operationunits = self._operationunits.replace(q, units)
				self._quantities_double.append([index_in_file, q])
				used_quantities.append( q )
				index_in_output += 1
		
		# Put data_log as object's variable
		self._data_log = data_log
		self._data_transform = data_transform
		self._cumulative = cumulative
		
		# In case of "vecto" quantity, get the species
		if species is not None:
			if self.operation != "vecto":
				raise Exception("Argument `species` only valid with quantity 'vecto'")
			self._species = str(species)
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with indices to the data arrays
		self._data = {}
		for i,t in enumerate(self._timesteps):
			self._data.update({ t : i })
		# Select timesteps if requested
		timestep_indices = kwargs.pop("timestep_indices", None)
		self._timesteps = self._selectTimesteps(timesteps, timestep_indices, self._timesteps)
		assert self._timesteps.size > 0, "Timesteps not found"
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		if raw is not None:
			self._type   .append("index of each process")
			self._shape  .append(self._nprocs)
			self._centers.append(self._np.arange(self._nprocs))
			self._label  .append("index of each process")
			self._units  .append("")
			self._log    .append(False)
			self._vunits = self._operationunits
			self._title  = self.operation

		elif map is not None:
			self._patch_length = [0]*self._ndim_fields
			for iaxis in range(self._ndim_fields):
				sim_length = self._ncels[iaxis]*self._cell_length[iaxis]
				self._patch_length[iaxis] = sim_length / float(self._number_of_patches[iaxis])
				centers = self._np.linspace(self._patch_length[iaxis]*0.5, sim_length-self._patch_length[iaxis]*0.5, self._number_of_patches[iaxis])
				self._type   .append("xyz"[iaxis])
				self._shape  .append(self._number_of_patches[iaxis])
				self._centers.append(centers)
				self._label  .append("xyz"[iaxis])
				self._units  .append("L_r")
				self._log    .append(False)
			self._vunits = self._operationunits
			self._title  = self.operation

		elif histogram is not None:
			bin_length = (histogram_max - histogram_min) / histogram_nsteps
			self._edges = self._np.linspace(histogram_min, histogram_max, histogram_nsteps+1)
			centers = self._edges[:-1] + bin_length * 0.5
			self._type   .append("quantity")
			self._shape  .append(histogram_nsteps)
			self._centers.append(centers)
			self._label  .append("quantity")
			self._units  .append("")
			self._log    .append(False)
			self._vunits = "1"
			self._title  = "number of processes"

		# Set the directory in case of exporting
		self._exportPrefix = "Performances"
		if len(used_quantities):
			self._exportPrefix += "_".join(used_quantities)
		if species is not None:
			self._exportPrefix += "_{}".format(species)
		self._exportDir = self._setExportDir(self._exportPrefix)

		# Finish constructor
		self.valid = True
		return kwargs

	# Method to print info
	def _info(self):
		return "Performances diagnostic "+self._title

	# get all available timesteps
	def getAvailableTimesteps(self):
		try:    times = [float(a.name[1:]) for a in self._h5items]
		except Exception as e: times = []
		return self._np.double(times)

	# get all available quantities
	def getAvailableQuantities(self):
		return self._availableQuantities_uint + self._availableQuantities_double

	# Method to obtain the data only
	def _getDataAtTime(self, t):
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		
		# Get arrays from requested field
		# get data
		index = self._data[t]
		C = []
		for dtype, quantities in (("uint",self._quantities_uint), ("double",self._quantities_double)):
			h5item = self._h5items[index]["quantities_"+dtype]
			for index_in_file, quantity in quantities:
				B = self._np.empty((self._nprocs,), dtype=dtype)
				h5item.read_direct( B, source_sel=self._np.s_[index_in_file,:] )
				# If not cumulative, make the difference with the previous time
				if not self._cumulative and quantity.startswith("timer") and index > 0:
					prevh5item = self._h5items[index-1]["quantities_"+dtype]
					prevB = self._np.empty((self._nprocs,), dtype=dtype)
					prevh5item.read_direct( prevB, source_sel=self._np.s_[index_in_file,:] )
					B -= prevB
				C.append( B )
		
		# Calculate the operation
		# First patch performance information
		if  self.operation in ["vecto", "mpi_rank"]:
			if self._mode != "raw":
				print("With quantities `vecto` or `mpi_rank`, only mode `raw` is supported")
				return []
			
			if "patches" not in self._h5items[index].keys():
				print("No patches group in timestep {}".format(str(t)))
				return []

			if self.operation=="vecto":

				if self._species not in self._h5items[index]["patches"].keys():
					print("Requested species {} does not have a group".format(self._species))
					return []
				patches_buffer = self._np.array(self._h5items[index]["patches"][self._species]["vecto"])

			elif self.operation=="mpi_rank":

				if "mpi_rank" not in self._h5items[index]["patches"].keys():
					print("Requested mpi_rank does not have a dataset")
					return []

				patches_buffer = self._np.array(self._h5items[index]["patches"]["mpi_rank"])

			# Get the position of the patches
			x_patches = self._np.array(self._h5items[index]["patches"]["x"][:])
			y_patches = self._np.array(self._h5items[index]["patches"]["y"][:])
			if "z" in self._h5items[index]["patches"]:
				z_patches = self._np.array(self._h5items[index]["patches"]["z"][:])
			else:
				z_patches = self._np.zeros_like(x_patches)	
			i_patch = self._np.ravel_multi_index(
				(x_patches, y_patches, z_patches),
				(x_patches.max()+1, y_patches.max()+1, z_patches.max()+1)
			)
			
			# Matrix of patches reconstituted
			A = self._np.empty_like(i_patch)
			A[i_patch] = patches_buffer
			A = self._np.squeeze(A.reshape([x_patches.max()+1, y_patches.max()+1, z_patches.max()+1]))

		# Or global performance information
		else:
			A = eval(self._operation)
		
		# If raw requested
		if self._mode == "raw":
			return A
		
		# If map requested
		elif self._mode == "map":
			
			# Extract the array "hindex"
			hindices = self._np.empty((self._nprocs,), dtype="uint")
			h5item = self._h5items[index]["quantities_uint"]
			index_in_file = self._availableQuantities_uint.index("hindex")
			h5item.read_direct( hindices, source_sel=self._np.s_[index_in_file,:] )

			if self._ndim_fields == 1:
				# Make a matrix with MPI ranks at each patch location
				ranks = self._np.empty((self._number_of_patches[0],), dtype=self._np.uint32)
				previous_h = 0
				rank = 0
				for h in hindices[1:]:
					ranks[previous_h:h] = rank
					rank += 1
					previous_h = h
				ranks[previous_h:] = rank

				# For each patch, associate the data of corresponding MPI rank
				return A[ranks]

			elif self._ndim_fields == 2:
				# Make a matrix with patch indices on the Hilbert curve
				if not hasattr(self, "_curvematrix"):
					if self.patch_arrangement == 'hilbertian':
						self._curvematrix = HilbertCurveMatrix2D(self._m[0], self._m[1], oversize=1)
					elif self.patch_arrangement == 'linearized_XY':
						self._curvematrix = LinXYCurveMatrix2D(self._number_of_patches, oversize=1)
					elif self.patch_arrangement == 'linearized_YX':
						self._curvematrix = LinYXCurveMatrix2D(self._number_of_patches, oversize=1)
					else:
						print("Error: patch arrangement "+str(self.patch_arrangement)+" not implemented")
						return []
					
				# Make a matrix with MPI ranks at each patch location
				self._ranks = PartitionMatrix( self._curvematrix, hindices, oversize=1 )

				# For each patch, associate the data of corresponding MPI rank
				return A[self._ranks[1:-1,1:-1]]

		# If histogram requested
		elif self._mode == "hist":
			histogram, _ = self._np.histogram( A, self._edges )
			return histogram

	# Convert data to VTK format
	def toVTK(self,numberOfPieces=1,axis_quantity="patch"):
		"""
		Export the performance data to Vtk
		"""

		# Checking of the arguments
		if axis_quantity not in ["patch","grid"]:
			print("axis_quantity must be `patch` or `grid`")
			return []

		# Creation of the directory and base name
		self._mkdir(self._exportDir)
		fileprefix = self._exportDir + self._exportPrefix + self.operation

		vtk = VTKfile()

		# If 2D data, then do a streak plot
		if self._ndim_fields == 2:

			if self._verbose: print("2D non implemented for VTK conversion")

		# Else 3D data
		elif self._ndim_fields == 3:

			# Loop over the requested time steps
			for istep,step in enumerate(self._timesteps):

				if self._verbose: print("Step: {}, file {}_{}.pvti".format(istep, fileprefix, int(step)))

				raw = self._dataAtTime(self._timesteps[istep])
				shape = list(raw.shape)
				origin = [0,0,0]
				extent = []
				for i in range(self._ndim_fields):
					extent += [0,shape[i]-1]
				# The axis of the gird are used
				if (axis_quantity == "grid"):
					spacings = [0,0,0]
					for i in range(self._ndim_fields):
						spacings[i] = self.simulation.namelist.Main.grid_length[i] / shape[i]
				# Axis are the patch x, y, z indexes
				else:
					spacings = [1,1,1]

				data = self._np.ascontiguousarray(raw.flatten(order='F'), dtype='float32')
				arr = vtk.Array(data, self._title)
				vtk.WriteImage(arr, origin, extent, spacings, fileprefix+"_{:08d}.pvti".format(int(step)), numberOfPieces)

			if self._verbose: print("Successfully exported to VTK, folder='"+self._exportDir)

	def _prepareExtent(self):
		self._extent = [
			self._xfactor*self._xoffset, self._xfactor*( self._xoffset + self._number_of_patches[0]*self._patch_length[0] ),
			self._yfactor*self._yoffset, self._yfactor*( self._yoffset + self._number_of_patches[1]*self._patch_length[1] )
		]
	
	def _calculateMPIcontours_2D(self):
		# Add lines to visualize MPI contours
		# Vertical lines
		vlines_i    = []
		vlines_jmin = []
		vlines_jmax = []
		vdiff = (self._np.diff(self._ranks, axis=1) != 0)
		vindices = self._np.flatnonzero(self._np.any(vdiff, axis=0)) # i-indices where vertical lines occur
		for i in vindices:
			j = self._np.flatnonzero(self._np.diff(vdiff[:,i])) # starts and ends of vertical lines
			vlines_i += [ self._np.full((len(j)//2), i, dtype=self._np.uint32) ]
			vlines_jmin += [ j[ 0::2 ] ]
			vlines_jmax += [ j[ 1::2 ] ]
		vlines_i    = (self._xoffset + self._np.concatenate( vlines_i    )*self._patch_length[0])*self._xfactor
		vlines_jmin = (self._yoffset + self._np.concatenate( vlines_jmin )*self._patch_length[1])*self._yfactor
		vlines_jmax = (self._yoffset + self._np.concatenate( vlines_jmax )*self._patch_length[1])*self._yfactor

		# Horizontal lines
		hlines_j    = []
		hlines_imin = []
		hlines_imax = []
		hdiff = (self._np.diff(self._ranks, axis=0) != 0)
		hindices = self._np.flatnonzero(self._np.any(hdiff, axis=1)) # j-indices where horizontal lines occur
		for j in hindices:
			i = self._np.flatnonzero(self._np.diff(hdiff[j,:])) # starts and ends of horizontal lines
			hlines_j += [ self._np.full((len(i)//2), j, dtype=self._np.uint32) ]
			hlines_imin += [ i[ 0::2 ] ]
			hlines_imax += [ i[ 1::2 ] ]
		hlines_j    = (self._yoffset + self._np.concatenate( hlines_j    )*self._patch_length[1])*self._yfactor
		hlines_imin = (self._xoffset + self._np.concatenate( hlines_imin )*self._patch_length[0])*self._xfactor
		hlines_imax = (self._xoffset + self._np.concatenate( hlines_imax )*self._patch_length[0])*self._xfactor
		
		return vlines_i, vlines_jmin, vlines_jmax, hlines_j, hlines_imin, hlines_imax
	
	def _plotOnAxes_2D_(self, ax, A):
		# Display the data
		self._plot = ax.imshow( self._np.flipud(A), extent=self._extent, **self.options.image)
		vlines_i, vlines_jmin, vlines_jmax, hlines_j, hlines_imin, hlines_imax = self._calculateMPIcontours_2D()
		self._vlines = ax.vlines( vlines_i, vlines_jmin, vlines_jmax, **self.options.plot)
		self._hlines = ax.hlines( hlines_j, hlines_imin, hlines_imax, **self.options.plot)
		return self._plot

	def _animateOnAxes_2D_(self, ax, A):
		# Display the data
		self._plot.set_data( self._np.flipud(A))
		vlines_i, vlines_jmin, vlines_jmax, hlines_j, hlines_imin, hlines_imax = self._calculateMPIcontours_2D()
		try:
			self._vlines.remove()
			self._hlines.remove()
		except:
			ax.collections.remove(self._vlines)
			ax.collections.remove(self._hlines)
		self._vlines = ax.vlines( vlines_i, vlines_jmin, vlines_jmax, **self.options.plot)
		self._hlines = ax.hlines( hlines_j, hlines_imin, hlines_imax, **self.options.plot)
		return self._plot
