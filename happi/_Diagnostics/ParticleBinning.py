from .Diagnostic import Diagnostic
from .._Utils import *

class ParticleBinning(Diagnostic):
	"""Class for loading a particle binning diagnostic"""
	
	_diagName = "ParticleBinning"
	hasComposite = False
	
	def _init(self, diagNumber=None, timesteps=None, subset=None, sum=None, data_log=False, include={}, **kwargs):
		
		# Search available diags
		diag_numbers, diag_names = self.simulation.getDiags(self._diagName)
		
		if diagNumber is None:
			self._error += ["Printing available %s:" % self._diagName]
			self._error += ["------------------------------------------------"]
			for diagNumber in diag_numbers:
				self._error += [self._printInfo(self._getInfo(diagNumber))]
			if len(diag_numbers)==0:
				self._error += ["      No %s found" % self._diagName]
			return
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check the requested diags are ok
		if type(diagNumber) is int:
			if diagNumber<0:
				self._error += ["Argument 'diagNumber' cannot be a negative integer."]
				return
			self.operation = '#' + str(diagNumber)
		elif type(diagNumber) is str:
			if diagNumber in diag_names:
				i = diag_names.index(diagNumber)
				self.operation = '#' + str(diag_numbers[i])
			else:
				self.operation = diagNumber
		else:
			self._error += ["Argument 'diagNumber' must be and integer or a string."]
			return
		
		# Get list of requested diags
		self._myinfo = {}
		self._diags = sorted(set([ int(d[1:]) for d in self._re.findall('#\d+',self.operation) ]))
		for d in self._diags:
			try:
				info = self._getInfo(d)
				if info is False: raise
				self._myinfo.update({ d:info })
			except:
				self._error += ["%s #%d invalid or non-existent" % (self._diagName,d)]
				return
		# Test the operation
		self._include = include
		try:
			exec(self._re.sub('#\d+','1.',self.operation), self._include, {"t":0})
		except ZeroDivisionError: pass
		except:
			self._error += ["Cannot understand operation '"+self.operation+"'"]
			return
		# Verify that all requested diags all have the same shape
		self._axes = {}
		self._naxes = {}
		for d in self._diags:
			self._axes .update ({ d:self._myinfo[d]["axes"] })
			self._naxes.update ({ d:len(self._axes[d]) })
			if self._naxes[d] != self._naxes[self._diags[0]]:
				self._error += [
					"All diagnostics in operation '%s' must have as many axes. %s #%d has %d axes and #%d has %d axes"
					% (self.operation, self._diagName, d, self._naxes[d], self._diags[0], self._naxes[self._diags[0]])
				]
				return
			for a in self._axes[d]:
				if self._axes[d] != self._axes[self._diags[0]]:
					self._error += [
						"In operation '%s', %s #%d and #%d must have the same shape."
						% (self.operation, self._diagName, d, self._diags[0])
					]
					return
		
		self._axes  = self._axes [self._diags[0]]
		self._naxes = self._naxes[self._diags[0]]
		
		# Check subset
		if subset is None: subset = {}
		elif type(subset) is not dict:
			self._error += ["Argument `subset` must be a dictionary"]
			return
		
		# Check sum
		if sum is None: sum = {}
		elif type(sum) is not dict:
			self._error += ["Argument 'sum' must be a dictionary"]
			return
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# Get available timesteps
		self._timesteps = {}
		self._alltimesteps = {}
		self._indexOfTime  = {}
		self._h5items = {}
		for d in self._diags:
			# Gather data from all timesteps, and the list of timesteps
			items = {}
			for path in self._results_path:
				f = self._h5py.File(path+self._os.sep+self._diagName+str(d)+'.h5')
				items.update( dict(f) )
			items = sorted(items.items())
			self._h5items[d] = [it[1] for it in items]
			self._timesteps[d] = self._np.array([ int(it[0].strip("timestep")) for it in items ])
			self._alltimesteps[d] = self._np.copy(self._timesteps[d])
			# fill the "_indexOfTime" dictionary with indices to the data arrays
			self._indexOfTime.update({ d:{} })
			for i,t in enumerate(self._timesteps[d]):
				self._indexOfTime[d].update({ t : i })
			# If timesteps is None, then keep all timesteps, otherwise, select timesteps
			if timesteps is not None:
				try:
					self._timesteps[d] = self._selectTimesteps(timesteps, self._timesteps[d])
				except:
					self._error += ["Argument 'timesteps' must be one or two non-negative integers"]
					return
			# Verify that timesteps are the same for all diagnostics
			if (self._timesteps[d] != self._timesteps[self._diags[0]]).any() :
				self._error += [
					"All diagnostics in operation '%s' must have the same timesteps. Diagnotic #%d has %d timesteps and #%d has %d timesteps"
					% (self.operation, d, len(self._timesteps[d]), self._diags[0], len(self._timesteps[self._diags[0]]))
				]
				return
		# Now we need to keep only one array of timesteps because they should be all the same
		self._timesteps  = self._timesteps [self._diags[0]]
		self._alltimesteps = self._alltimesteps[self._diags[0]]
		
		# Need at least one timestep
		if self._timesteps.size < 1:
			self._error += ["Timesteps not found"]
			return
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values for all diags
		plot_diff = []
		coeff = 1.
		spatialaxes = {"x":False, "y":False, "z":False}
		self._finalShape = [[]]*self._naxes
		self._sums = [False]*self._naxes
		self._selection = [self._np.s_[:]]*self._naxes
		uniform = True
		
		for iaxis in range(self._naxes):
			axis = self._axes[iaxis]
			axistype = axis["type"]
			
			# Find the vector of values along the axis
			if axis["log"]:
				edges = self._np.linspace(self._np.log10(axis["min"]), self._np.log10(axis["max"]), axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				edges = 10.**edges
				centers = 10.**(centers[:-1])
			else:
				edges = self._np.linspace(axis["min"], axis["max"], axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				centers = centers[:-1]
			axis.update({ "edges"   : edges   })
			axis.update({ "centers" : centers })
			
			# Find some quantities depending on the axis type
			overall_min = "-inf"; overall_max = "inf"
			axis_units = ""
			if   axistype in ["x","y","z","moving_x"]:
				axis_units = "L_r"
				spatialaxes[axistype[-1]] = True
			elif axistype in ["a","b"]:
				axis_units = "L_r"
				self.hasComposite = True
			elif axistype == "theta" and self._ndim_particles==2:
				axis_units = "rad"
				overall_min = "-3.141592653589793"
				overall_max = "3.141592653589793"
			elif axistype == "theta" and self._ndim_particles==3:
				axis_units = "rad"
				overall_min = "0"
				overall_max = "3.141592653589793"
			elif axistype == "phi":
				axis_units = "rad"
				overall_min = "-3.141592653589793"
				overall_max = " 3.141592653589793"
			elif axistype in ["px","py","pz","p"]:
				axis_units = "P_r"
			elif axistype in ["vx","vy","vz","v"]:
				axis_units = "V_r"
			elif axistype in ["vperp2"]:
				axis_units = "V_r**2"
				overall_min = "0"
			elif axistype == "gamma":
				overall_min = "1"
			elif axistype == "ekin":
				axis_units = "K_r"
				overall_min = "0"
			elif axistype == "charge":
				axis_units = "Q_r"
				overall_min = "0"
			elif axistype == "chi":
				overall_min = "0"
			
			# if this axis has to be summed, then select the bounds
			if axistype in sum:
				if axistype in subset:
					self._error += ["`subset` not possible on the same axes as `sum`"]
					return
				
				self._sums[iaxis] = True
				
				try:
					axis["sumInfo"], self._selection[iaxis], self._finalShape[iaxis] \
						= self._selectRange(sum[axistype], centers, axistype, axis_units, "sum", axis["edges_included"])
				except:
					return
				
				if axistype in ["x","y","z","moving_x"]:
					first_edge = edges[self._selection[iaxis].start or 0]
					last_edge  = edges[(self._selection[iaxis].stop or len(centers))]
					coeff /= last_edge - first_edge
				
				plot_diff.append( self._np.ones((self._finalShape[iaxis],)) )
			
			# if not summed
			else:
				# If taking a subset of this axis
				if axistype in subset:
					try:
						axis["subsetInfo"], self._selection[iaxis], self._finalShape[iaxis] \
							= self._selectSubset(subset[axistype], centers, axistype, axis_units, "subset")
					except:
						return
					# If selection is not a slice (meaning only one element) then axis removed from plot
					if type(self._selection[iaxis]) is not slice and axistype in ["x","y","z","moving_x"]:
						first_edge = edges[self._selection[iaxis]]
						last_edge  = edges[self._selection[iaxis]+1]
						coeff /= last_edge - first_edge
				
				# If no subset, or subset has more than 1 point, use this axis in the plot
				if type(self._selection[iaxis]) is slice:
					self._type   .append(axistype)
					self._shape  .append(axis["size"])
					self._centers.append(centers[self._selection[iaxis]])
					self._log    .append(axis["log"])
					self._label  .append(axistype)
					self._units  .append(axis_units)
					if axistype == "theta" and self._ndim_particles==3:
						uniform = False
						plot_diff.append(self._np.diff(self._np.cos(edges))[self._selection[iaxis]])
					else:
						plot_diff.append(self._np.diff(edges)[self._selection[iaxis]])
					self._finalShape[iaxis] = len(self._centers[-1])
					if axis["log"]:
						uniform = False
				else:
					plot_diff.append( self._np.ones((self._finalShape[iaxis],)) )
		
		self._selection = tuple(self._selection)
		
		# Build units
		titles = {}
		units = {}
		axes_units = [unit or "1" for unit in self._units if (self.hasComposite or unit!="L_r")]
		axes_units = (" / ( " + " * ".join(axes_units) + " )") if axes_units else ""
		for d in self._diags:
			titles.update({ d:"??" })
			val_units = "1"
			deposited_quantity = self._myinfo[d]["deposited_quantity"]
			if   deposited_quantity == "weight":
				titles[d] = "Number" + ("" if self.hasComposite else " density")
				val_units = "1" if self.hasComposite else "N_r"
			elif deposited_quantity == "weight_charge":
				titles[d] = "Charge" + ("" if self.hasComposite else " density")
				val_units = "Q_r" if self.hasComposite else "N_r * Q_r"
			elif deposited_quantity == "weight_ekin":
				titles[d] = "Energy" + ("" if self.hasComposite else " density")
				val_units = "K_r" if self.hasComposite else "N_r * K_r"
			elif deposited_quantity[:15] == "weight_charge_v":
				titles[d] = "J"+deposited_quantity[-1] + (" x Volume" if self.hasComposite else "")
				val_units = "J_r/N_r" if self.hasComposite else "J_r"
			elif deposited_quantity[:8] == "weight_p":
				titles[d] = "P"+deposited_quantity[8:] + ("" if self.hasComposite else " density")
				val_units = "P_r" if self.hasComposite else "N_r * P_r"
			elif deposited_quantity[:8] == "weight_v":
				titles[d] = "Pressure "+deposited_quantity[8]+deposited_quantity[11] + (" x Volume" if self.hasComposite else "")
				val_units = "K_r" if self.hasComposite else "N_r * K_r"
			elif deposited_quantity[:13] == "weight_ekin_v":
				titles[d] = "Energy ("+deposited_quantity[-1]+") flux" + (" x Volume" if self.hasComposite else " density")
				val_units = "K_r" if self.hasComposite else "N_r * K_r"
			elif deposited_quantity[:13] == "weight_power": # for radiation spectrum
				titles[d] = "Power" + ("" if self.hasComposite else " density")
				val_units = "K_r / T_r" if self.hasComposite else "N_r * K_r / T_r"
			units[d] = val_units + axes_units
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for d in self._diags:
			self._vunits = self._vunits.replace("#"+str(d), "( "+units[d]+" )")
			self._title  = self._title .replace("#"+str(d), titles[d])
		self._vunits = self.units._getUnits(self._vunits)
		
		# If any spatial dimension did not appear, then count it for calculating the correct density
		if self._ndim_particles>=1 and not spatialaxes["x"]: coeff /= self._ncels[ 0]*self._cell_length[ 0]
		if self._ndim_particles>=2 and not spatialaxes["y"]: coeff /= self._ncels[ 1]*self._cell_length[ 1]
		if self._ndim_particles==3 and not spatialaxes["z"]: coeff /= self._ncels[-1]*self._cell_length[-1]
		
		# Calculate the array that represents the bins sizes in order to get units right.
		# This array will be the same size as the plotted array
		if uniform:
			self._bsize = 1.
			for d in plot_diff:
				self._bsize *= d[0]
		else:
			if len(plot_diff)==0:
				self._bsize = 1.
			elif len(plot_diff)==1:
				self._bsize = plot_diff[0]
			else:
				self._bsize = self._np.prod( self._np.array( self._np.meshgrid( *plot_diff ) ), axis=0)
				self._bsize = self._bsize.transpose([1,0]+list(range(2,len(plot_diff))))
		self._bsize = 1. / self._bsize
		if not self.hasComposite: self._bsize *= coeff
		
		# Set the directory in case of exporting
		self._exportPrefix = self._diagName+"_"+"-".join([str(d) for d in self._diags])
		self._exportDir = self._setExportDir(self._exportPrefix)
		
		# Finish constructor
		self.valid = True
		return kwargs
	
	# Gets info about diagnostic number "diagNumber"
	def _getInfo(self,diagNumber):
		info = {}
		for path in self._results_path:
			# Open file
			try:
				file = path+self._os.sep+self._diagName+str(diagNumber)+'.h5'
				f = self._h5py.File(file, 'r')
			except:
				return False
			# get attributes from file
			axes = []
			deposited_quantity = "weight_power" # necessary for radiation spectrum
			time_average = None
			# Parse each attribute
			for name, value in f.attrs.items():
				if name == "deposited_quantity":
					try:
						deposited_quantity = bytes.decode(value)
					except:
						deposited_quantity = "user_function"
				elif name == "time_average":
					time_average = int(value)
				elif name == "species":
					species = bytes.decode(value.strip()).split() # get all species numbers
					species = [int(s) for s in species]
				elif name[0:4] == "axis":
					n = int(name[4:]) # axis number
					sp = bytes.decode(value).split()
					while len(axes)<n+1: axes.append({}) # extend the array to the adequate size
					axes[n] = dict(
						type = sp[0], min = float(sp[1]), max = float(sp[2]), size = int(sp[3]),
						log = bool(int(sp[4])), edges_included = bool(int(sp[5])), coefficients = 0 if sp[6]=="[]" else eval(sp[6])
					)
				elif name == "photon_energy_axis":
					sp = bytes.decode(value).split()
					axes.append( dict(
						type = "gamma", min = float(sp[0]), max = float(sp[1]), size = int(sp[2]),
						log = bool(int(sp[3])), edges_included = bool(int(sp[4]))
					))
			f.close()
			# Verify that the info corresponds to the diag in the other paths
			if info == {}:
				info = {"#":diagNumber, "deposited_quantity":deposited_quantity, "tavg":time_average, "species":species, "axes":axes}
			else:
				if deposited_quantity!=info["deposited_quantity"] or axes!=info["axes"]:
					print(self._diagName+" #"+str(diagNumber)+" in path '"+path+"' is incompatible with the other ones")
					return False
		return info

	
	# Prints the info obtained by the function "getInfo"
	@staticmethod
	def _printInfo(info):
		if not info:
			return "Error while reading file(s)"
		
		# 1 - diag number, type and list of species
		species = ""
		for i in range(len(info["species"])): species += str(info["species"][i])+" " # reconstitute species string
		printedInfo = "#"+str(info["#"])+" - "+info["deposited_quantity"]+" of species # "+species+"\n"
		
		# 2 - period and time-averaging
		if info["tavg"] and info["tavg"] > 1:
			printedInfo += "    Averaging over "+str(info["tavg"])+" timesteps\n"
		
		# 3 - axes
		for i in range(len(info["axes"])):
			axis = info["axes"][i]
			logscale = "" if not axis["log"] else " [ LOG SCALE ] "
			edges    = "" if not axis["edges_included"] else " [ INCLUDING EDGES ] "
			printedInfo += "    "+axis["type"]+" from "+str(axis["min"])+" to "+str(axis["max"]) \
				   +" in "+str(axis["size"])+" steps "+logscale+edges+"\n"
		return printedInfo
	
	# Method to print info on all included diags
	def _info(self):
		info = ""
		for d in self._diags:
			info += self._printInfo(self._myinfo[d])+"\n"
		if len(self.operation)>2: info += "Operation : "+self.operation+"\n"
		for ax in self._axes:
			if "sumInfo" in ax: info += ax["sumInfo"]+"\n"
			if "subsetInfo" in ax: info += ax["subsetInfo"]+"\n"
		return info
	
	# get all available timesteps for a given diagnostic
	def getAvailableTimesteps(self, diagNumber=None):
		# if argument "diagNumber" not provided, return the times calculated in __init__
		if diagNumber is None:
			return self._alltimesteps
		# Otherwise, get the timesteps specifically available for the single requested diagnostic
		else:
			times = set()
			for path in self._results_path:
				try:
					file = path+self._os.sep+self._diagName+str(diagNumber)+'.h5'
					f = self._h5py.File(file, 'r')
				except:
					print("Cannot open file "+file)
					return self._np.array([])
				times.update( set(f.keys()) )
				f.close()
			times = [int(t.strip("timestep")) for t in times]
			return self._np.array(times)
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Get arrays from all requested diagnostics
		A = {}
		for d in self._diags:
			# find the index of the array corresponding to the requested timestep
			try:
				index = self._indexOfTime[d][t]
			except:
				print("Timestep "+str(t)+" not found in this diagnostic")
				return []
			# get data
			B = self._np.empty(self._finalShape)
			try:
				self._h5items[d][index].read_direct(B, source_sel=self._selection) # get array
			except:
				B = self._np.squeeze(B)
				self._h5items[d][index].read_direct(B, source_sel=self._selection) # get array
				B = self._np.reshape(B, self._finalShape)
			B[self._np.isnan(B)] = 0.
			# Divide by the bins size
			B *= self._bsize
			# Append this diag's array for the operation
			A.update({ d:B })
		# Calculate operation
		data_operation = self.operation
		for d in reversed(self._diags):
			data_operation = data_operation.replace("#"+str(d),"A["+str(d)+"]")
		A = eval(data_operation, self._include, locals())
		# Apply the summing
		for iaxis in range(self._naxes):
			if self._sums[iaxis]:
				A = self._np.sum(A, axis=iaxis, keepdims=True)
		# remove summed axes
		A = self._np.squeeze(A)
		return A
