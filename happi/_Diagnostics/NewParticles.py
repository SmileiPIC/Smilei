from .Diagnostic import Diagnostic
from .ParticleList import ParticleList
from .._Utils import *

class NewParticles(ParticleList):
	"""Class for loading a NewParticles diagnostic"""
	
	_diagType = "NewParticles"
	
	def _init(self, species=None, select="", axes=[], chunksize=20000000, **kwargs):
		
		# If argument 'species' not provided, then print available species and leave
		if species is None:
			species = self.simulation.getNewParticlesSpecies()
			error = ["Argument `species` not provided"]
			if len(species)>0:
				error += ["Printing available species for NewParticles:"]
				error += ["--------------------------------------------"]
				error += ["\n".join(species)]
			else:
				error += ["No files found for NewParticles"]
			raise Exception("\n".join(error))
		
		# Get info from the hdf5 files + verifications
		# -------------------------------------------------------------------
		self.species  = species
		files = self._findFiles("NewParticles")
		self._h5files = [self._h5py.File(file, "r") for file in files]
		
		# List available properties
		self._short_properties_from_raw["birth_time"] = "t"
		self._raw_properties_from_short = {v:k for k,v in self._short_properties_from_raw.items()}
		T0 = self._h5files[0]["data/0/particles/"+self.species]
		self.available_properties = {v for k,v in self._short_properties_from_raw.items() if k in T0}
		
		# Find the overlap between files
		if len(self._h5files) == 1:
			# If only 1 file, there is 1 range containing the whole file
			self.nParticles = self._h5files[0]["iteration_npart"][-1,1]
			self.ranges_in_file = {0:[self._np.s_[0:self.nParticles]]}
		else:
			# If several files, this is more complicated
			iterations_in_file = []
			npart_in_file = []
			all_iterations = self._np.array([], dtype=int)
			# for each file, get the list of iterations dumped and the corresponding number of particles for each iteration
			for i,f in enumerate(self._h5files):
				iterations_in_file += [f["iteration_npart"][:,0]]
				npart_in_file += [f["iteration_npart"][:,1]]
				# all_iteration is the union of all iterations in all files
				all_iterations = self._np.union1d( all_iterations, iterations_in_file[-1] )
			# Now, assign one file to each existing iteration. If 2 files contain the same iteration, the last file prevails
			file_vs_iteration = self._np.empty_like( all_iterations, dtype = int )
			for i,f in enumerate(self._h5files):
				locations = self._np.isin( all_iterations, iterations_in_file[i], assume_unique=True )
				file_vs_iteration[locations] = i
			# find for which iterations there is a change in the file number
			change_locs = 1 + self._np.flatnonzero(self._np.diff(file_vs_iteration))
			change_locs = self._np.concatenate(([0], change_locs, [len(all_iterations)]))
			# Now the goal is to calculate several ranges of particles that each file should handle
			# Instead of looping every iteration, we loop on the iterations where a file change occurs
			self.nParticles = 0
			self.ranges_in_file = {i:[] for i in range(len(self._h5files))}
			for i in range(len(change_locs)-1):
				# Get the locations of relevant iterations
				istart, istop = change_locs[i], change_locs[i+1]
				# Get the corresponding file
				ifile = file_vs_iteration[istart]
				# Get the relevant iterations
				iterstart = all_iterations[istart]
				iterstop = all_iterations[-1]+1 if istop >= len(all_iterations) else all_iterations[istop]
				# Get the iteration locations in the file
				iterindexstart, iterindexstop = self._np.searchsorted(iterations_in_file[ifile], [iterstart, iterstop])
				# Get the corresponding particle indices in the file
				ipartstart = 0 if iterindexstart == 0 else npart_in_file[ifile][iterindexstart-1]
				ipartstop  = npart_in_file[ifile][iterindexstop -1]
				# Add this range of particles
				self.nParticles += ipartstop - ipartstart
				self.ranges_in_file[ifile] += [self._np.s_[ipartstart:ipartstop]]
		
		# Manage axes
		self._initAxes(axes)
		
		# Select particles
		self.selectedParticles = None
		self.nselectedParticles = self.nParticles
		if type(select) is str:
			# Transform the operation
			try:
				particleSelector = select
				doubleProps = []
				int16Props = []
				for prop in self.available_properties:
					particleSelector, nsubs = self._re.subn(r"\b"+prop+r"\b", "particle_data['"+prop+"']", particleSelector)
					if nsubs > 0:
						if   prop == "q" : int16Props  += [prop]
						else             : doubleProps += [prop]
			except:
				raise Exception("Error: `select` not understood")
			# Nothing to select if empty `select` does not contain any particle property
			if len(int16Props) + len(doubleProps) > 0.:
				# Get particle selection chunk by chunk
				self.selectedParticles = []
				for particle_data in self._iterAllParticles(chunksize):
					self.selectedParticles += [self._np.flatnonzero( eval(particleSelector) )]
				# Concatenate all chunks
				self.selectedParticles = self._np.concatenate( self.selectedParticles )
				self.nselectedParticles = len(self.selectedParticles)
		
		# Otherwise, the selection can be a list of particle IDs
		elif select:
			# Loop only IDs
			self.selectedParticles = []
			for particle_data in self._iterAllParticles(chunksize, axes=["Id"]):
				self.selectedParticles += [self._np.flatnonzero(self._np.in1d(particle_data["Id"], select))]
			self.selectedParticles = self._np.sort( self._np.concatenate( self.selectedParticles ) )
			self.nselectedParticles = len(self.selectedParticles)
		
		# Now split the selected particles per the different files where they belong
		if self.nselectedParticles < self.nParticles:
			sel = {}
			offset = 0
			for ifile, ranges in self.ranges_in_file.items():
				sel[ifile] = []
				for r in ranges:
					istart, istop = self._np.searchsorted( self.selectedParticles, [offset + r.start, offset + r.stop] )
					sel[ifile] += [self.selectedParticles[istart:istop]-offset]
					offset += r.stop - r.start
				sel[ifile] = self._np.concatenate( sel[ifile] )
			self.selectedParticles = sel
		
		if select and self._verbose:
			print("Kept "+str(self.nselectedParticles)+" particles")
		
		# Finish constructor
		assert "timesteps" not in kwargs and "timestep_indices" not in kwargs,\
			"The NewParticles diagnostic does not accept 'timesteps' or 'timestep_indices' as arguments."
		self._timesteps = [None]
		self.valid = True
		return kwargs
	
	# We override the get and getData methods
	def getData(self):
		for data in self.iterParticles(chunksize=self.nParticles):
			return data
	
	# Iterator on particles
	def iterParticles(self, chunksize=1):
		if self.nselectedParticles < self.nParticles:
			return self._iterSelectedParticles(chunksize, self.axes)
		else:
			return self._iterAllParticles(chunksize, self.axes)
	
	# Iterator for all particles
	def _iterAllParticles(self, chunksize=1, axes=None):
		if axes is None:
			axes = self.axes
		# Make buffers
		properties = self._raw_properties_from_short
		data = {}
		for axis in axes:
			if axis == "Id":
				data[axis] = self._np.empty((chunksize,), dtype=self._np.uint64)
			elif axis == "q":
				data[axis] = self._np.empty((chunksize,), dtype=self._np.int16 )
			else:
				data[axis] = self._np.empty((chunksize,), dtype=self._np.double)
		
		# Loop files
		npart_in_chunk = 0
		for ifile,f in enumerate(self._h5files):
			# Get the main group containing all datasets
			try:
				group = f["data/0/particles/"+self.species]
			except:
				continue
			# Loop on the ranges of particles
			for r in self.ranges_in_file[ifile]:
				nremaining_in_chunk = chunksize - npart_in_chunk
				start = r.start
				# Yield chunks until the current range of particles is complete
				while r.stop - start > nremaining_in_chunk:
					for axis in axes:
						group[properties[axis]].read_direct(data[axis], source_sel=self._np.s_[start:start+nremaining_in_chunk], dest_sel=self._np.s_[npart_in_chunk:])
					yield data
					# Start new chunk
					start += nremaining_in_chunk
					npart_in_chunk = 0
					nremaining_in_chunk = chunksize
				# If the range still has a few particles, start filling the next chunk
				if start < r.stop:
					for axis in axes:
						group[properties[axis]].read_direct(data[axis], source_sel=self._np.s_[start:], dest_sel=self._np.s_[npart_in_chunk:npart_in_chunk+r.stop-start])
					npart_in_chunk += r.stop-start
		# If the last chunk was not complete, yield that
		if npart_in_chunk > 0:
			# Resize buffers
			for axis in axes:
				data[axis].resize((npart_in_chunk,))
			yield data
	
	# Iterator for selected particles
	def _iterSelectedParticles(self, chunksize=1, axes=None):
		# Make buffers
		properties = self._raw_properties_from_short
		data = {}
		for axis in axes:
			if axis == "Id":
				data[axis] = self._np.empty((chunksize,), dtype=self._np.uint64)
			elif axis == "q":
				data[axis] = self._np.empty((chunksize,), dtype=self._np.int16 )
			else:
				data[axis] = self._np.empty((chunksize,), dtype=self._np.double)
		
		# Loop files
		npart_in_chunk = 0
		for ifile,f in enumerate(self._h5files):
			# Get the main group containing all datasets
			try:
				group = f["data/0/particles/"+self.species]
			except:
				continue
			#selection for this file
			sel = self.selectedParticles[ifile]
			# Yield chunks until the current selection of particles is complete
			nremaining_in_chunk = chunksize - npart_in_chunk
			start = 0; stop = sel.size
			while stop - start > nremaining_in_chunk:
				print("chunk %d"%chunkstart)
				for axis in axes:
					group[properties[axis]].read_direct(data[axis], source_sel=sel[start:start+nremaining_in_chunk], dest_sel=self._np.s_[npart_in_chunk:])
				yield data
				# Start new chunk
				start += nremaining_in_chunk
				npart_in_chunk = 0
				nremaining_in_chunk = chunksize
			# If the selection still has a few particles, start filling the next chunk
			if start < stop:
				for axis in axes:
					group[properties[axis]].read_direct(data[axis], source_sel=sel[start:], dest_sel=self._np.s_[npart_in_chunk:npart_in_chunk+stop-start])
				npart_in_chunk += stop-start
		# If the last chunk was not complete, yield that
		if npart_in_chunk > 0:
			# Resize buffers
			for axis in axes:
				data[axis].resize((npart_in_chunk,))
			yield data
	
	# Method to get info
	def _info(self):
		info = "New particles diagnostic: species '"+self.species+"' containing "+str(self.nParticles)+" particles"
		if self.selectedParticles:
			info += "\n\twith selection of "+str(self.nselectedParticles)+" particles"
		info += "\n\tAxes: " + ", ".join(self.axes)
		return info
	
	# Prevent using some plotting functions
	def streak(self, *args, **kwargs):
		raise("Error: Cannot use streak in NewParticles")
	def animate(self, *args, **kwargs):
		raise("Error: Cannot use animate in NewParticles")
	def slide(self, *args, **kwargs):
		raise("Error: Cannot use slide in NewParticles")
	def toVTK(self, *args, **kwargs):
		raise("Error: Cannot use toVTK in NewParticles")
	
	# Override plotting functions
	def _prepare3(self):
		if self._tmpdata is None:
			A = self.getData()
			self._tmpdata = []
			for axis in self.axes: self._tmpdata.append( A[axis] )
		return True
	
	def _setTitle(ax, t=None):
		ax.set_title(self.options.title.format(self._vlabel), self.options.labels_font["title"])
	
	def _plotOnAxes_0D(self, ax, t, cax_id=0):
		pass
	
	def _plotOnAxes_1D(self, ax, t, cax_id=0):
		pass
	
	def _plotOnAxes_2D(self, ax, t, cax_id=0):
		xfactor = self.options.xfactor or 1.
		yfactor = self.options.yfactor or 1.
		try   : ax.set_prop_cycle (None)
		except:	ax.set_color_cycle(None)
		ax.plot(xfactor*self._tmpdata[0], yfactor*self._tmpdata[1], '.', **self.options.plot)
		# Add labels and options
		ax.set_xlabel(self._xlabel)
		ax.set_ylabel(self._ylabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.ymin, ymax=self.options.ymax)
		self._setTitle(ax)
		self._setAxesOptions(ax)
		return 1
	
