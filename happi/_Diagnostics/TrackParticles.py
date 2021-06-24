from .Diagnostic import Diagnostic
from .._Utils import *

# Define a function that finds the next closing character in a string
def findClosingCharacter(string, character, start=0):
	stack = []
	associatedBracket = {")":"(", "]":"[", "}":"{"}
	for i in range(start, len(string)):
		if string[i] == character and len(stack)==0:
			return i
		elif string[i] in ("(", "[", "{"):
			stack.append(string[i])
		elif string[i] in (")", "]", "}"):
			if len(stack)==0:
				raise Exception("Error in selector syntax: missing `"+character+"`")
			if stack[-1]!=associatedBracket[string[i]]:
				raise Exception("Error in selector syntax: missing closing parentheses or brackets")
			del stack[-1]
	raise Exception("Error in selector syntax: missing `"+character+"`")


class TrackParticles(Diagnostic):
	"""Class for loading a TrackParticles diagnostic"""

	def _init(self, species=None, select="", axes=[], timesteps=None, sort=True, sorted_as="", length=None, chunksize=20000000, **kwargs):

		# If argument 'species' not provided, then print available species and leave
		if species is None:
			species = self.simulation.getTrackSpecies()
			error = ["Argument `species` not provided"]
			if len(species)>0:
				error += ["Printing available tracked species:"]
				error += ["-----------------------------------"]
				error += ["\n".join(species)]
			else:
				error += ["No tracked particles files found"]
			raise Exception("\n".join(error))
		
		if type(sort) not in [bool, str]:
			raise Exception("Argument `sort` must be `True` or `False` or a string")
		if not sort and select!="":
			raise Exception("Cannot select particles if not sorted")
		self._sort = sort
		
		
		# Get info from the hdf5 files + verifications
		# -------------------------------------------------------------------
		self.species  = species
		self._h5items = {}
		disorderedfiles = self._findDisorderedFiles()
		self._short_properties_from_raw = {
			"id":"Id", "position/x":"x", "position/y":"y", "position/z":"z",
			"momentum/x":"px", "momentum/y":"py", "momentum/z":"pz",
			"charge":"q", "weight":"w", "chi":"chi",
			"E/x":"Ex", "E/y":"Ey", "E/z":"Ez", "B/x":"Bx", "B/y":"By", "B/z":"Bz"
		}
		
		# If sorting allowed, find out if ordering needed
		needsOrdering = False
		if sort:
			if type(sort) is str:
				# The sorted file gets a name from `sorted_as`
				if type(sorted_as) is not str or self._re.search(r"[^a-zA-Z0-9_]","_"+sorted_as):
					raise Exception("Argument `sorted_as` must be a keyword composed of letters and numbers")
				if not sorted_as:
					raise Exception("Argument `sorted_as` is required when `sort` is a selection")
			if sorted_as:
				sorted_as = "_"+sorted_as
			orderedfile = self._results_path[0]+self._os.sep+"TrackParticles_"+species+sorted_as+".h5"
			needsOrdering = self._needsOrdering(orderedfile)
			if sorted_as and not needsOrdering and type(sort) is str:
				print("WARNING: ordered file `"+"TrackParticles_"+species+sorted_as+".h5"+"` already exists.")
				print("         Skipping sorting operation.")
		
		# Find times in disordered files
		if not sort or needsOrdering:
			self._locationForTime = {}
			for file in disorderedfiles:
				f = self._h5py.File(file, "r")
				self._locationForTime.update( {int(t):[f,it] for it, t in enumerate(f["data"].keys())} )
			self._lastfile = f
			self._timesteps = self._np.array(sorted(self._locationForTime))
			self._alltimesteps = self._np.copy(self._timesteps)
			
			# List available properties
			try: # python 2
				self._raw_properties_from_short = {v:k for k,v in self._short_properties_from_raw.iteritems()}
				T0 = next(self._lastfile["data"].itervalues())["particles/"+self.species]
			except: # python 3
				self._raw_properties_from_short = {v:k for k,v in self._short_properties_from_raw.items()}
				T0 = next(iter(self._lastfile["data"].values()))["particles/"+self.species]
			self.available_properties = [v for k,v in self._short_properties_from_raw.items() if k in T0]
		
		# If sorting allowed, then do the sorting
		if sort:
			# If the first path does not contain the ordered file (or it is incomplete), we must create it
			if needsOrdering:
				self._orderFiles(orderedfile, chunksize, sort)
				if self._needsOrdering(orderedfile):
					raise Exception("Ordering not succesful")
			# Create arrays to store h5 items
			self._lastfile = self._h5py.File(orderedfile, "r")
			for prop in ["Id", "x", "y", "z", "px", "py", "pz", "q", "w", "chi",
			             "Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
				if prop in self._lastfile:
					self._h5items[prop] = self._lastfile[prop]
			self.available_properties = list(self._h5items.keys())
			# Memorize the locations of timesteps in the files
			self._locationForTime = {t:it for it, t in enumerate(self._lastfile["Times"])}
			self._timesteps = self._np.array(sorted(self._lastfile["Times"]))
			self._alltimesteps = self._np.copy(self._timesteps)
			self.nParticles = self._h5items["Id"].shape[1]
		
		# Add moving_x in the list of properties
		if "x" in self.available_properties:
			file = disorderedfiles[0]
			with self._h5py.File(file, "r") as f:
				try: # python 2
					D = next(f["data"].itervalues())
				except: # python 3
					D = next(iter(f["data"].values()))
				if "x_moved" in D.attrs:
					self.available_properties += ["moving_x"]
		
		# Get available times in the hdf5 file
		if self._timesteps.size == 0:
			raise Exception("No tracked particles found")
		# If specific timesteps requested, narrow the selection
		if timesteps is not None:
			try:
				ts = self._np.array(self._np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self._timesteps = self._timesteps[ self._np.nonzero((self._timesteps>=ts[0]) * (self._timesteps<=ts[1]))[0] ]
				elif ts.size==1:
					# get nearest time
					self._timesteps = self._np.array(self._timesteps[ self._np.array([(self._np.abs(self._timesteps-ts)).argmin()]) ])
				else:
					raise
			except:
				raise Exception("Argument `timesteps` must be one or two non-negative integers")
		# Need at least one timestep
		if self._timesteps.size < 1:
			raise Exception("Timesteps not found")
		
		# Select particles
		# -------------------------------------------------------------------
		if sort:
			self.selectedParticles = self._selectParticles( select, True, chunksize )
			if self.selectedParticles is None:
				raise Exception("Error: argument 'select' must be a string or a list of particle IDs")
			
			# Remove particles that are not actually tracked during the requested timesteps
			if self._verbose: print("Removing dead particles ...")
			if type(self.selectedParticles) is not slice and len(self.selectedParticles) > 0:
				first_time = self._locationForTime[self._timesteps[ 0]]
				last_time  = self._locationForTime[self._timesteps[-1]]+1
				IDs = self._readUnstructuredH5(self._h5items["Id"], self.selectedParticles, first_time, last_time)
				dead_particles = self._np.flatnonzero(self._np.all( self._np.isnan(IDs) + (IDs==0), axis=0 ))
				self.selectedParticles = self._np.delete( self.selectedParticles, dead_particles )

			# Calculate the number of selected particles
			if type(self.selectedParticles) is slice:
				self.nselectedParticles = self.nParticles
			else:
				self.nselectedParticles = len(self.selectedParticles)
			if self.nselectedParticles == 0:
				raise Exception("No particles found")
			if self._verbose: print("Kept "+str(self.nselectedParticles)+" particles")

		# Manage axes
		# -------------------------------------------------------------------
		if type(axes) is not list:
			raise Exception("Error: Argument 'axes' must be a list")
		# if axes provided, verify them
		if len(axes)>0:
			self.axes = axes
			for axis in axes:
				if axis not in self.available_properties:
					raise Exception(
						"Error: Argument 'axes' has item '"+str(axis)+"' unknown.\n"
						+ "       Available axes are: "+(", ".join(sorted(self.available_properties)))
					)
		# otherwise use default
		else:
			self.axes = self.available_properties
		
		# Get x_moved if necessary
		if "moving_x" in self.axes:
			self._XmovedForTime = {}
			for file in disorderedfiles:
				with self._h5py.File(file, "r") as f:
					for t in f["data"].keys():
						self._XmovedForTime[int(t)] = f["data"][t].attrs["x_moved"]
		
		# Then figure out axis units
		self._type = self.axes
		self._factors = []
		for axis in self.axes:
			axisunits = ""
			if axis == "Id":
				self._centers.append( [0, 281474976710655] )
			elif axis in ["x" , "y" , "z", "moving_x"]:
				axisunits = "L_r"
				self._centers.append( [0., self.namelist.Main.grid_length[{"x":0,"y":1,"z":-1}[axis[-1]]]] )
			elif axis in ["px", "py", "pz"]:
				axisunits = "P_r"
				self._centers.append( [-1., 1.] )
			elif axis == "w":
				axisunits = "N_r * L_r^%i" % self._ndim_particles
				self._centers.append( [0., 1.] )
			elif axis == "q":
				axisunits = "Q_r"
				self._centers.append( [-10., 10.] )
			elif axis == "chi":
				axisunits = "1"
				self._centers.append( [0., 2.] )
			elif axis[0] == "E":
				axisunits = "E_r"
				self._centers.append( [-1., 1.] )
			elif axis[0] == "B":
				axisunits = "B_r"
				self._centers.append( [-1., 1.] )
			self._log += [False]
			self._label += [axis]
			self._units += [axisunits]
			if axis == "Id":
				self._factors += [1]
			else:
				factor, _ = self.units._convert(axisunits, None)
				self._factors += [factor]
		self._title = "Track particles '"+species+"'"
		self._shape = [0]*len(self.axes)
		self._centers = [self._np.array(c) for c in self._centers]

		# Hack to work with 1 axis
		if len(axes)==1: self._vunits = self._units[0]
		else: self._vunits = ""

		# Set the directory in case of exporting
		self._exportPrefix = "TrackParticles_"+self.species+"_"+"".join(self.axes)
		self._exportDir = self._setExportDir(self._exportPrefix)

		self._rawData = None

		# Finish constructor
		self.length = length or self._timesteps[-1]
		self.valid = True
		return kwargs

	def _needsOrdering(self, orderedfile):
		if not self._os.path.isfile(orderedfile):
			return True
		else:
			try:
				f = self._h5py.File(orderedfile, "r")
				if "finished_ordering" not in f.attrs.keys():
					return True
			except:
				self._os.remove(orderedfile)
				return True
			finally:
				f.close()
		return False

	def _selectParticles( self, select, already_sorted, chunksize ):
		if type(select) is str:
			# Parse the selector
			i = 0
			operation = ""
			seltype = []
			selstr = []
			timeSelector = []
			particleSelector = []
			doubleProps = []
			int16Props = []
			while i < len(select):
				if i+4<len(select) and select[i:i+4] in ["any(","all("]:
					seltype += [select[i:i+4]]
					if seltype[-1] not in ["any(","all("]:
						raise Exception("Error in selector syntax: unknown argument "+seltype[-1][:-1])
					comma = findClosingCharacter(select, ",", i+4)
					parenthesis = findClosingCharacter(select, ")", comma+1)
					timeSelector += [select[i+4:comma]]
					selstr += [select[i:parenthesis]]
					try:
						timeSelector[-1] = "self._alltimesteps["+self._re.sub(r"\bt\b","self._alltimesteps",timeSelector[-1])+"]"
						eval(timeSelector[-1])
					except:
						raise Exception("Error in selector syntax: time selector not understood in "+select[i:i+3]+"()")
					try:
						particleSelector += [select[comma+1:parenthesis]]
						doubleProps += [[]]
						int16Props += [[]]
						for prop in self.available_properties:
							(particleSelector[-1], nsubs) = self._re.subn(r"\b"+prop+r"\b", "properties['"+prop+"'][:actual_chunksize]", particleSelector[-1])
							if nsubs > 0:
								if   prop == "q" : int16Props [-1] += [prop]
								else             : doubleProps[-1] += [prop]
					except:
						raise Exception("Error in selector syntax: not understood: "+select[i:parenthesis+1])
					operation += "stack["+str(len(seltype)-1)+"]"
					i = parenthesis+1
				elif not already_sorted and not select[i].isspace():
					raise Exception("Complex selection operations not allowed for unsorted files (bad character %s)"%select[i])
				else:
					operation += select[i]
					i+=1
			nOperations = len(seltype)
			
			# Nothing to select if empty operation
			if len(operation)==0.:
				return self._np.s_[:]
			
			# Execute the selector
			if self._verbose: print("Selecting particles ... (this may take a while)")
			
			def makeBuffers(size):
				properties = {}
				for k in range(nOperations):
					for prop in int16Props[k]:
						if prop not in properties:
							properties[prop] = self._np.empty((size,), dtype=self._np.int16)
					for prop in doubleProps[k]:
						if prop not in properties:
							properties[prop] = self._np.empty((size,), dtype=self._np.double)
				properties["Id"] = self._np.empty((size,), dtype=self._np.uint64)
				return properties
			
			if already_sorted:
				# Setup the chunks of particles (if too many particles)
				chunks = ChunkedRange(self.nParticles, chunksize)
				# Allocate buffers
				selectedParticles = self._np.array([], dtype=self._np.uint64)
				properties = makeBuffers(chunks.adjustedchunksize)
				# Loop on chunks
				for chunkstart, chunkstop, actual_chunksize in chunks:
					# Execute each of the selector items
					stack = []
					for k in range(nOperations):
						selection = self._np.empty((chunks.adjustedchunksize,), dtype=bool)
						if   seltype[k] == "any(": selection.fill(False)
						elif seltype[k] == "all(": selection.fill(True )
						requiredProps = doubleProps[k] + int16Props[k] + ["Id"]
						# Loop times
						for time in eval(timeSelector[k]):
							if self._verbose: print("   Selecting block `"+selstr[k]+")`, at time "+str(time))
							# Extract required properties from h5 files
							it = self._locationForTime[time]
							for prop in requiredProps:
								self._h5items[prop].read_direct(properties[prop], source_sel=self._np.s_[it,chunkstart:chunkstop], dest_sel=self._np.s_[:actual_chunksize])
							# Calculate the selector
							selectionAtTimeT = eval(particleSelector[k]) # array of True or False
							# Combine with selection of previous times
							selectionAtTimeT[self._np.isnan(selectionAtTimeT)] = False
							existing = properties["Id"][:actual_chunksize]>0 # existing particles at that timestep
							if   seltype[k] == "any(": selection[existing] += selectionAtTimeT[existing]
							elif seltype[k] == "all(": selection *= selectionAtTimeT * existing
						stack.append(selection)
					# Merge all stack items according to the operations
					selectedParticles = self._np.union1d( selectedParticles, eval(operation).nonzero()[0] )
			else:
				# Execute the selector item
				selectedParticles = self._np.array([], dtype="uint64")
				k = 0
				requiredProps = doubleProps[k] + int16Props[k] + ["Id"]
				# Loop times
				for time in eval(timeSelector[k]):
					if self._verbose: print("   Selecting block `"+selstr[k]+")`, at time "+str(time))
					# Get group in file
					[f, it] = self._locationForTime[time]
					group = f["data/"+"%010i"%time+"/particles/"+self.species]
					npart = group["id"].shape[0]
					# Loop on chunks
					selectionAtTimeT = []
					for chunkstart, chunkstop, actual_chunksize in ChunkedRange(npart, chunksize):
						# Allocate buffers
						properties = makeBuffers(actual_chunksize)
						# Extract required properties from h5 files
						for prop in requiredProps:
							group[self._raw_properties_from_short[prop]].read_direct(properties[prop], source_sel=self._np.s_[chunkstart:chunkstop], dest_sel=self._np.s_[:actual_chunksize])
						# Calculate the selector
						sel = eval(particleSelector[k]) # array of True or False
						selectionAtTimeT.append(properties["Id"][sel])
					selectionAtTimeT = self._np.concatenate(selectionAtTimeT)
					# Combine with selection of previous times
					if   seltype[k] == "any(": selectedParticles = self._np.union1d(selectedParticles, selectionAtTimeT)
					elif seltype[k] == "all(": selectedParticles = self._np.intersect1d(selectedParticles, selectionAtTimeT)
			selectedParticles.sort()
			return selectedParticles

		# Otherwise, the selection can be a list of particle IDs
		else:
			try:
				IDs = self._lastfile["unique_Ids"] # get all available IDs
				return self._np.flatnonzero(self._np.in1d(IDs, select)) # find the requested IDs
			except:
				return

	# Method to get info
	def _info(self):
		info = "Track particles: species '"+self.species+"'"
		if self._sort:
			info += " containing "+str(self.nParticles)+" particles"
			if self.nselectedParticles != self.nParticles:
				info += "\n\twith selection of "+str(self.nselectedParticles)+" particles"
		info += "\n\tAxes: " + ", ".join(self.axes)
		return info

	# Read hdf5 dataset faster with unstrusctured list of indices
	def _readUnstructuredH5(self, dataset, indices, first_time, last_time=None):
		if last_time is None:
			last_time = first_time + 1
		cs = 1000
		if type(indices) is slice or len(indices) < cs:
			return dataset[first_time:last_time, indices]
		else:
			n = len(indices)
			result = self._np.empty(( last_time - first_time, n ), dtype=dataset.dtype)
			chunksize = min(cs,n)
			nchunks = int(n/chunksize)
			chunksize = int(n / nchunks)
			chunkstop = 0
			for ichunk in range(nchunks):
				chunkstart = chunkstop
				chunkstop  = min(chunkstart + chunksize, n)
				result[:,chunkstart:chunkstop] = dataset[first_time:last_time, indices[chunkstart:chunkstop]]
			return result

	# get all available timesteps
	def getAvailableTimesteps(self):
		return self._alltimesteps

	# Get a list of disordered files
	def _findDisorderedFiles(self):
		disorderedfiles = []
		for path in self._results_path:
			file = path+self._os.sep+"TrackParticlesDisordered_"+self.species+".h5"
			if self._os.path.isfile(file):
				disorderedfiles += [file]
		if not disorderedfiles:
			raise Exception("No TrackParticles files")
		return disorderedfiles

	# Make the particles ordered by Id in the file, in case they are not
	def _orderFiles( self, fileOrdered, chunksize, sort ):
		if self._verbose:
			print("Ordering particles ... (this could take a while)")
			if type(sort) is str:
				print("    Selecting particles according to "+sort)
		try:
			# If ordered file already exists, find out which timestep was done last
			latestOrdered = -1
			if self._os.path.isfile(fileOrdered):
				f0 = self._h5py.File(fileOrdered, "r+")
				try:    latestOrdered = f0.attrs["latestOrdered"]
				except: pass
			# otherwise, make new (ordered) file
			else:
				f0 = self._h5py.File(fileOrdered, "w")
			# Open the last file and get the number of particles from each MPI
			last_time = self._timesteps[-1]
			last_file, _ = self._locationForTime[last_time]
			number_of_particles = (last_file["data/"+"%010i/"%last_time+"latest_IDs"][()] % (2**32)).astype('uint32')
			if self._verbose: print("Number of particles: "+str(number_of_particles.sum()))
			# Calculate the offset that each MPI needs
			offset = self._np.cumsum(number_of_particles, dtype='uint64')
			total_number_of_particles = offset[-1]
			offset = self._np.roll(offset, 1)
			offset[0] = 0
			# Do the particle selection if requested
			selectedIds = None
			selectedIndices = self._np.s_[:]
			nparticles_to_write = total_number_of_particles
			if type(sort) is str:
				selectedIds = self._selectParticles( sort, False, chunksize )
				nparticles_to_write = len(selectedIds)
			# Make datasets if not existing already
			size = (len(self._timesteps), nparticles_to_write)
			group = last_file["data/"+"%010i/"%last_time+"particles/"+self.species]
			for k, name in self._short_properties_from_raw.items():
				try   : f0.create_dataset(name, size, group[k].dtype, fillvalue=(0 if name=="Id" else self._np.nan))
				except: pass
			# Loop times and fill arrays
			for it, t in enumerate(self._timesteps):
				
				# Skip previously-ordered times
				if it<=latestOrdered: continue
				
				if self._verbose: print("    Ordering @ timestep = "+str(t))
				f, _ = self._locationForTime[t]
				group = f["data/"+"%010i/"%t+"particles/"+self.species]
				nparticles = group["id"].size
				if nparticles == 0: continue
				
				# If not too many particles, sort all at once
				if nparticles_to_write < chunksize and nparticles < chunksize:
					# Get the Ids and find where they should be stored in the final file
					if selectedIds is None:
						locs = (
							group["id"][()].astype("uint32") # takes the second hald of id (meaning particle number)
							+ offset[ (group["id"][()]>>32).astype("uint32") & 0b111111111111111111111111 ]
							-1
						)
					else:
						_,selectedIndices,locs = self._np.intersect1d( group["id"][()], selectedIds, return_indices=True )
					# Loop datasets and order them
					if len(locs) > 0:
						for k, name in self._short_properties_from_raw.items():
							if k not in group: continue
							ordered = self._np.empty((nparticles_to_write, ), dtype=group[k].dtype)
							if k == "id": ordered.fill(0)
							else        : ordered.fill(self._np.nan)
							ordered[locs] = group[k][()][selectedIndices]
							f0[name].write_direct(ordered, dest_sel=self._np.s_[it,:])
				
				# If too many particles, sort by chunks
				else:
					data = {}
					for k, name in self._short_properties_from_raw.items():
						data[k] = self._np.empty((chunksize,), dtype=self._np.int16 if k == "charge" else self._np.double)
					# Loop chunks of the output
					for first_o, last_o, npart_o in ChunkedRange(nparticles_to_write, chunksize):
						for k, name in self._short_properties_from_raw.items():
							if k not in group: continue
							if k == "id": data[k].fill(0)
							else        : data[k].fill(self._np.nan)
						# Loop chunks of the input
						for first_i, last_i, npart_i in ChunkedRange(nparticles, chunksize):
							# Obtain IDs
							ID = group["id"][first_i:last_i]
							# Extract useful IDs for this output chunk
							if selectedIds is None:
								loc_in_output = ID.astype("uint32") + offset[ (ID>>32).astype("uint32") & 0b111111111111111111111111 ] - 1
								keep = self._np.flatnonzero((loc_in_output >= first_o) * (loc_in_output < last_o))
								loc_in_output = loc_in_output[keep] - first_o
							else:
								_,keep,loc_in_output = self._np.intersect1d( ID, selectedIds[first_o:last_o], return_indices=True )
							# Fill datasets with this chunk
							for k, name in self._short_properties_from_raw.items():
								if k not in group: continue
								data[k][loc_in_output] = group[k][first_i:last_i][keep]
						# Accumulated data is written out
						for k, name in self._short_properties_from_raw.items():
							if k not in group: continue
							f0[name][it, first_o:last_o] = data[k][:npart_o]
						
				# Indicate that this iteration was succesfully ordered
				f0.attrs["latestOrdered"] = it
				f0.flush()
			if self._verbose: print("    Finalizing the ordering process")
			# Create the "Times" dataset
			f0.create_dataset("Times", data=self._timesteps)
			# Create the "unique_Ids" dataset
			if selectedIds is None:
				unique_Ids = self._np.empty((nparticles_to_write,), dtype=f0["Id"].dtype)
				for iMPI in range(number_of_particles.size):
					for first, last, npart in ChunkedRange(number_of_particles[iMPI], chunksize):
						o = int(offset[iMPI])
						unique_Ids[o+first : o+last] = ((iMPI<<32) + 1) + self._np.arange(first,last,dtype='uint64')
				f0.create_dataset("unique_Ids", data=unique_Ids)
			else:
				f0.create_dataset("unique_Ids", data=selectedIds)
			# Indicate that the ordering is finished
			f0.attrs["finished_ordering"] = True
			# Close file
			f0.close()
		except Exception as e:
			print("Error in the ordering of the tracked particles")
			if self._verbose:
				print(e)
			raise
		finally:
			# Close disordered files
			for t in self._locationForTime:
				self._locationForTime[t][0].close()
		if self._verbose: print("Ordering succeeded")

	# Method to generate the raw data (only done once)
	def _generateRawData(self, times=None):
		if not self._validate(): return
		self._prepare1() # prepare the vfactor

		if self._sort:
			if self._rawData is None:
				self._rawData = {}
				first_time = self._locationForTime[self._timesteps[0]]
				last_time  = self._locationForTime[self._timesteps[-1]] + 1
				if self._verbose: print("Loading data ...")
				# fill up the data
				ID = self._readUnstructuredH5(self._h5items["Id"], self.selectedParticles, first_time, last_time)
				deadParticles = (ID==0).nonzero()
				for axis in self.axes:
					if self._verbose: print("   axis: "+axis)
					if axis == "Id":
						self._rawData[axis] = ID
					else:
						if axis=="moving_x":
							data = self._readUnstructuredH5(self._h5items["x"], self.selectedParticles, first_time, last_time)
							for it, time in enumerate(self._timesteps):
								data[it,:] -= self._XmovedForTime[time]
						else:
							data = self._readUnstructuredH5(self._h5items[axis], self.selectedParticles, first_time, last_time)
						data[deadParticles] = self._np.nan
						self._rawData[axis] = data

				if self._verbose: print("Process broken lines ...")
				# Add the lineBreaks array which indicates where lines are broken (e.g. loop around the box)
				self._rawData['brokenLine'] = self._np.zeros((self.nselectedParticles,), dtype=bool)
				self._rawData['lineBreaks'] = {}
				if self._timesteps.size > 1:
					dt = self._np.diff(self._timesteps)*self.timestep
					for axis in ["x","y","z"]:
						if axis in self.axes:
							dudt = self._np.diff(self._rawData[axis],axis=0)
							for i in range(dudt.shape[1]): dudt[:,i] /= dt
							dudt[~self._np.isnan(dudt)] = 0. # NaNs already break lines
							# Line is broken if velocity > c
							self._rawData['brokenLine'] += self._np.abs(dudt).max(axis=0) > 1.
							broken_particles = self._np.flatnonzero(self._rawData['brokenLine'])
							for broken_particle in broken_particles:
								broken_times = list(self._np.flatnonzero(self._np.abs(dudt[:,broken_particle]) > 1.)+1)
								if broken_particle in self._rawData['lineBreaks'].keys():
									self._rawData['lineBreaks'][broken_particle] += broken_times
								else:
									self._rawData['lineBreaks'][broken_particle] = broken_times
				# Add the times array
				self._rawData["times"] = self._timesteps
				if self._verbose: print("... done")

		# If not sorted, get different kind of data
		else:
			if self._rawData is None:
				self._rawData = {}

			if self._verbose: print("Loading data ...")
			properties = dict(self._raw_properties_from_short, moving_x="position/x")
			if times is None: times = self._timesteps
			for time in times:
				if time in self._rawData: continue
				[f, timeIndex] = self._locationForTime[time]
				group = f["data/"+"%010i"%time+"/particles/"+self.species]
				self._rawData[time] = {}
				for axis in self.axes:
					self._rawData[time][axis] = group[properties[axis]][()]
				if "moving_x" in self.axes:
					self._rawData[time]["moving_x"] -= self._XmovedForTime[time]

			if self._verbose: print("... done")

	# We override the get and getData methods
	def getData(self, timestep=None):
		if not self._validate(): return
		self._prepare1() # prepare the vfactor

		if timestep is None:
			ts = self._timesteps
		elif timestep not in self._timesteps:
			print("ERROR: timestep "+str(timestep)+" not available")
			return {}
		else:
			ts = [timestep]
			indexOfRequestedTime = self._np.where(self._timesteps==timestep)

		if len(ts)==1 and not self._sort:
			self._generateRawData(ts)
		else:
			self._generateRawData()

		data = {}
		data.update({ "times":ts })

		if self._sort:
			for axis, factor in zip(self.axes, self._factors):
				if timestep is None:
					data[axis] = self._rawData[axis]
				else:
					data[axis] = self._rawData[axis][indexOfRequestedTime]
				data[axis] *= factor
		else:
			for t in ts:
				data[t] = {}
				for axis, factor in zip(self.axes, self._factors):
					data[t][axis] = self._rawData[t][axis] * factor
		return data

	def get(self):
		return self.getData()

	# Iterator on UNSORTED particles for a given timestep
	def iterParticles(self, timestep, chunksize=1):
		if not self._validate(): return
		self._prepare1() # prepare the vfactor

		if timestep not in self._timesteps:
			print("ERROR: timestep "+str(timestep)+" not available")
			return

		properties = {"moving_x":"x"}
		properties.update( self._raw_properties_from_short )

		disorderedfiles = self._findDisorderedFiles()
		for file in disorderedfiles:
			f = self._h5py.File(file, "r")
			# This is the timestep for which we want to produce an iterator
			try:
				group = f["data/"+("%010d"%timestep)+"/particles/"+self.species]
			except:
				f.close()
				continue
			npart = group["id"].size
			ID          = self._np.empty((chunksize,), dtype=self._np.uint64)
			data_double = self._np.empty((chunksize,), dtype=self._np.double)
			data_int16  = self._np.empty((chunksize,), dtype=self._np.int16 )
			for chunkstart in range(0, npart, chunksize):
				chunkend = chunkstart + chunksize
				if chunkend > npart:
					chunkend = npart
					ID          = self._np.empty((chunkend-chunkstart,), dtype=self._np.uint64)
					data_double = self._np.empty((chunkend-chunkstart,), dtype=self._np.double)
					data_int16  = self._np.empty((chunkend-chunkstart,), dtype=self._np.int16 )
				data = {}
				for axis in self.axes:
					if axis == "Id":
						group[properties[axis]].read_direct(ID, source_sel=self._np.s_[chunkstart:chunkend])
						data[axis] = ID.copy()
					elif axis == "q":
						group[properties[axis]].read_direct(data_int16, source_sel=self._np.s_[chunkstart:chunkend])
						data[axis] = data_int16.copy()
					elif axis == "moving_x":
						group[properties["x"]].read_direct(data_double, source_sel=self._np.s_[chunkstart:chunkend])
						data[axis] = data_double.copy()
					else:
						group[properties[axis]].read_direct(data_double, source_sel=self._np.s_[chunkstart:chunkend])
						data[axis] = data_double.copy()
				yield data
			f.close()

	# We override _prepare3
	def _prepare3(self):
		if not self._sort:
			print("Cannot plot non-sorted data")
			return False
		if self._tmpdata is None:
			A = self.getData()
			self._tmpdata = []
			for axis in self.axes: self._tmpdata.append( A[axis] )
		return True


	# We override the plotting methods
	def _animateOnAxes_0D(self, ax, t, cax_id=0):
		pass
	
	def _animateOnAxes_1D(self, ax, t, cax_id=0):
		timeSelection = (self._timesteps<=t)*(self._timesteps>=t-self.length)
		times = self._timesteps[timeSelection]
		A     = self._tmpdata[0][timeSelection,:]
		if times.size == 1:
			times = self._np.double([times, times]).squeeze()
			A = self._np.double([A, A]).squeeze()
		try   : ax.set_prop_cycle (None)
		except:	ax.set_color_cycle(None)
		self._plot = ax.plot(self._tfactor*times, self._vfactor*A, **self.options.plot)
		ax.set_xlabel(self._tlabel)
		ax.set_ylabel(self.axes[0]+" ("+self.units.vname+")")
		self._setLimits(ax, xmax=self._tfactor*self._timesteps[-1], ymin=self.options.vmin, ymax=self.options.vmax)
		ax.set_title(self._title) # override title
		self._setAxesOptions(ax)
		return self._plot
	
	def _animateOnAxes_2D(self, ax, t, cax_id=0):
		if hasattr(ax, "_lines"):
			if self in ax._lines:
				for line in ax._lines[self]:
					line.remove()
				del ax._lines[self]
		else:
			ax._lines = {}
		tmin = t-self.length
		tmax = t
		timeSelection = (self._timesteps<=tmax)*(self._timesteps>=tmin)
		selected_times = self._np.flatnonzero(timeSelection)
		itmin = selected_times[0]
		itmax = selected_times[-1]
		# Plot first the non-broken lines
		x = self._tmpdata[0][timeSelection,:][:,~self._rawData["brokenLine"]]
		y = self._tmpdata[1][timeSelection,:][:,~self._rawData["brokenLine"]]
		try   : ax.set_prop_cycle (None)
		except:	ax.set_color_cycle(None)
		ax._lines[self] = ax.plot(self._xfactor*x, self._yfactor*y, **self.options.plot)
		# Then plot the broken lines
		try   : ax.hold("on")
		except: pass
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
					ax._lines[self] += ax.plot(self._xfactor*x[iti:itf], self._yfactor*y[iti:itf], color=prevline.get_color(), **self.options.plot)
				else:
					prevline, = ax.plot(self._xfactor*x[iti:itf], self._yfactor*y[iti:itf], **self.options.plot)
				ax._lines[self] += [prevline]
				if breaks[ibrk] > itmax: break
		try   : ax.hold("off")
		except: pass
		# Add labels and options
		ax.set_xlabel(self._xlabel)
		ax.set_ylabel(self._ylabel)
		self._setLimits(ax, xmin=self.options.xmin, xmax=self.options.xmax, ymin=self.options.ymin, ymax=self.options.ymax)
		self._setTitle(ax, t)
		self._setAxesOptions(ax)
		return 1
	
	_plotOnAxes_0D = _animateOnAxes_0D
	_plotOnAxes_1D = _animateOnAxes_1D
	_plotOnAxes_2D = _animateOnAxes_2D
	
	# Convert data to VTK format
	def toVTK(self, rendering="trajectory", data_format="xml"):
		"""
		Export the data to Vtk
		"""
		if not self._validate(): return

		if not self._sort:
			print("Cannot export non-sorted data")
			return

		if self._ndim_particles != 3:
			print ("Cannot export tracked particles of a "+str(self._ndim_particles)+"D simulation to VTK")
			return

		# The specified rendering option is checked
		if rendering not in ["trajectory","cloud"]:
			print ("Rendering of type {} is not valid. It should be `trajectory` or `cloud`.".format(rendering))
			return

		# The specified data format is checked
		if data_format not in ["xml","vtk"]:
			print ("Format of type {} is not valid. Should be `xml` or `vtk` ".format(data_format))
			return

		self._mkdir(self._exportDir)
		fileprefix = self._exportDir + self._exportPrefix + "_" + rendering 

		ntimes = len(self._timesteps)

		# Determine the correct file extension according to the given data format
		if data_format == "xml":
			extension = "vtp"
		else:
			extension = "vtk"

		# Creation of a customed vtk object
		vtk = VTKfile()
		
		# Require x, y and z
		xaxis = "x"
		if "x" not in self.axes:
			xaxis = "moving_x"
		if xaxis not in self.axes or "y" not in self.axes or "z" not in self.axes:
			print("Error exporting tracked particles to VTK: axes 'x', 'y' and 'z' are required")
			return
		
		# Cloud mode: each time step is a separated cloud of particles
		# If there is only one timestep, the trajectory mode becomes a cloud
		if (ntimes == 1)or(rendering == "cloud"):

			data = self.getData()

			for istep,step in enumerate(self._timesteps):
				
				data_clean_step = {}
				
				# Clean data at istep: remove NaN
				mask = self._np.ones(len(data[self.axes[0]][istep]), dtype=bool)
				for ax in self.axes:
					mask = self._np.logical_and(mask,self._np.logical_not(self._np.isnan(self._np.asarray(data[ax][istep]))))
				for ax in self.axes:
					#print(ax,data[ax][istep])
					data_clean_step[ax] = self._np.asarray(data[ax][istep])[mask]
				
				pcoords_step = self._np.stack((data_clean_step[xaxis],data_clean_step["y"],data_clean_step["z"])).transpose()
				pcoords_step = self._np.ascontiguousarray(pcoords_step, dtype='float32')

				# Convert pcoords that is a numpy array into vtkFloatArray
				pcoords_step = vtk.Array(pcoords_step, "")

				# List of scalar arrays
				attributes = []
				for ax in self.axes:
					if ax not in ["x", "y", "z", "moving_x", "Id"]:
						attributes += [vtk.Array(self._np.ascontiguousarray(data_clean_step[ax].flatten(),'float32'),ax)]
					# Integer arrays
					elif ax == "Id":
						attributes += [vtk.Array(self._np.ascontiguousarray(data_clean_step[ax].flatten(),'int32'),ax)]

				vtk.WriteCloud(pcoords_step, attributes, data_format, fileprefix+"_{:06d}.{}".format(step,extension))
				print("Exportation of {}_{:06d}.{}".format(fileprefix,step,extension))

			print("Successfully exported tracked particles to VTK, folder='"+self._exportDir)

		# Trajectory mode
		elif (rendering == "trajectory"):

			data = self.getData()
			pcoords = self._np.stack((data[xaxis],data["y"],data["z"])).transpose()
			npoints, nt, nd = pcoords.shape

			pcoords = self._np.reshape(pcoords, (npoints*nt, nd))
			pcoords = self._np.ascontiguousarray(pcoords, dtype='float32')

			# Convert pcoords that is a numpy array into vtkFloatArray
			pcoords = vtk.Array(pcoords, "")

			# Segments between points to describe the trajectories
			connectivity = self._np.ascontiguousarray([[nt]+[nt*i+j for j in range(nt)] for i in range(npoints)])

			# List of scalar arrays
			attributes = []
			for ax in self.axes:
				if ax not in ["x", "y", "z", "moving_x", "Id"]:
					attributes += [vtk.Array(self._np.ascontiguousarray(data[ax].flatten(),'float32'),ax)]
				# Integer arrays
				elif ax == "Id":
					attributes += [vtk.Array(self._np.ascontiguousarray(data[ax].flatten(),'int32'),ax)]

			vtk.WriteLines(pcoords, connectivity, attributes, data_format, fileprefix+".{}".format(extension))
			print("Successfully exported tracked particles to VTK, folder='"+self._exportDir)
