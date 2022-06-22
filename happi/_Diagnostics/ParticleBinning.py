from .Diagnostic import Diagnostic
from .._Utils import *

class ParticleBinning(Diagnostic):
	"""Class for loading a particle binning diagnostic"""
	
	_diagType = "ParticleBinning"
	hasComposite = False
	
	def _init(self, diagNumber=None, timesteps=None, subset=None, average=None, sum=None, data_log=False, data_transform=None, include={}, **kwargs):
		
		# Search available diags
		diag_numbers, diag_names = self.simulation.getDiags(self._diagType)
		
		if diagNumber is None:
			error = ["`diagNumber` is not defined"]
			error += ["Printing available %s:" % self._diagType]
			error += ["------------------------------------------------"]
			for diagNumber in diag_numbers:
				error += [self._printInfo(self._getInfo(self.simulation, self._diagType, diagNumber))]
			if len(diag_numbers)==0:
				error += ["      No %s found" % self._diagType]
			raise Exception("\n".join(error))
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check the requested diags are ok
		if type(diagNumber) is int:
			if diagNumber<0:
				raise Exception("Argument 'diagNumber' cannot be a negative integer.")
			self.operation = '#' + str(diagNumber)
		elif type(diagNumber) is str:
			if diagNumber in diag_names:
				i = diag_names.index(diagNumber)
				self.operation = '#' + str(diag_numbers[i])
			else:
				self.operation = diagNumber
		else:
			raise Exception("Argument 'diagNumber' must be and integer or a string.")
		
		# Get list of requested diags
		self._myinfo = {}
		self._diags = sorted(set([ int(d[1:]) for d in self._re.findall('#\d+',self.operation) ]))
		for d in self._diags:
			try:
				info = self._getInfo(self.simulation, self._diagType, d)
				if info is False: raise
				self._myinfo.update({ d:info })
			except Exception as e:
				raise Exception("%s #%d invalid or non-existent" % (self._diagType,d))
		# Test the operation
		self._include = include
		try:
			exec(self._re.sub('#\d+','1.',self.operation), self._include, {"t":0})
		except ZeroDivisionError:
			pass
		except Exception as e:
			raise Exception("Cannot understand operation '"+self.operation+"'")
		# Verify that all requested diags all have the same shape
		self._axes = {}
		self._naxes = {}
		for d in self._diags:
			self._axes .update ({ d:self._myinfo[d]["axes"] })
			self._naxes.update ({ d:len(self._axes[d]) })
			if self._naxes[d] != self._naxes[self._diags[0]]:
				raise Exception(
					"All diagnostics in operation '%s' must have as many axes. %s #%d has %d axes and #%d has %d axes"
					% (self.operation, self._diagType, d, self._naxes[d], self._diags[0], self._naxes[self._diags[0]])
				)
			if self._axes[d] != self._axes[self._diags[0]]:
				raise Exception(
					"In operation '%s', %s #%d and #%d must have the same shape."
					% (self.operation, self._diagType, d, self._diags[0])
				)
		
		self._axes  = self._axes [self._diags[0]]
		self._naxes = self._naxes[self._diags[0]]
		
		# Check subset
		if subset is None:
			subset = {}
		elif type(subset) is not dict:
			raise Exception("Argument `subset` must be a dictionary")
		
		# Check average
		if average is None:
			if type(sum) is dict:
				average = sum
			else:
				average = {}
		elif type(average) is not dict:
			raise Exception("Argument `average` must be a dictionary")
		
		# Put data_log as object's variable
		self._data_log = data_log
		self._data_transform = data_transform

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
				try:
					f = self._h5py.File(path+self._os.sep+self._diagType+str(d)+'.h5', 'r')
				except:
					continue
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
				except Exception as e:
					raise Exception("Argument 'timesteps' must be one or two non-negative integers")
			# Verify that timesteps are the same for all diagnostics
			if (self._timesteps[d] != self._timesteps[self._diags[0]]).any() :
				raise Exception(
					"All diagnostics in operation '%s' must have the same timesteps. Diagnotic #%d has %d timesteps and #%d has %d timesteps"
					% (self.operation, d, len(self._timesteps[d]), self._diags[0], len(self._timesteps[self._diags[0]]))
				)
		# Now we need to keep only one array of timesteps because they should be all the same
		self._timesteps  = self._timesteps [self._diags[0]]
		self._alltimesteps = self._alltimesteps[self._diags[0]]
		
		# Need at least one timestep
		if self._timesteps.size < 1:
			raise Exception("Timesteps not found")
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values for all diags
		self._spatialaxes = {"x":False, "y":False, "z":False}
		self.auto_axes = False
		user_axes = []
		for iaxis, axis in enumerate(self._axes):
			# Find some quantities depending on the axis type
			overall_min = "-inf"
			overall_max = "inf"
			axis["units"] = ""
			if   axis["type"] in ["x","y","z","moving_x"]:
				axis["units"] = "L_r"
				self._spatialaxes[axis["type"][-1]] = True
			elif axis["type"] in ["a","b"]:
				axis["units"] = "L_r"
				self.hasComposite = True
			elif axis["type"] == "theta" and self._ndim_particles==2:
				axis["units"] = "rad"
				overall_min = "-3.141592653589793"
				overall_max = "3.141592653589793"
			elif axis["type"] == "theta" and self._ndim_particles==3:
				axis["units"] = "rad"
				overall_min = "0"
				overall_max = "3.141592653589793"
			elif axis["type"] == "phi":
				axis["units"] = "rad"
				overall_min = "-3.141592653589793"
				overall_max = " 3.141592653589793"
			elif axis["type"] in ["px","py","pz","p"]:
				axis["units"] = "P_r"
			elif axis["type"] in ["vx","vy","vz","v"]:
				axis["units"] = "V_r"
			elif axis["type"] in ["vperp2"]:
				axis["units"] = "V_r**2"
				overall_min = "0"
			elif axis["type"] == "gamma":
				overall_min = "1"
			elif axis["type"] == "ekin":
				axis["units"] = "K_r"
				overall_min = "0"
			elif axis["type"] == "charge":
				axis["units"] = "Q_r"
				overall_min = "0"
			elif axis["type"] == "chi":
				overall_min = "0"
			elif axis["type"][:4] == "user":
				axis["units"] = ""
				user_axes += [axis["type"]]
			
			# Store average/subset info
			if axis["type"] in average:
				if axis["type"] in subset:
					raise Exception("`subset` not possible on the same axes as `average`")
				axis["average"] = average[axis["type"]]
			elif axis["type"] in subset:
				axis["subset"] = subset[axis["type"]]
			
			# Get limits when auto limits
			if axis["min"] == "auto":
				axis["auto_min"] = [it.attrs["min%d"%iaxis] for it in self._h5items[self._diags[0]]]
				self.auto_axes = True
			if axis["max"] == "auto":
				axis["auto_max"] = [it.attrs["max%d"%iaxis] for it in self._h5items[self._diags[0]]]
				self.auto_axes =  True
		
		# Set all axes average/subset
		self._updateAxes(self._timesteps[0])
		
		# Build units
		titles = {}
		units = {}
		axes_units = [unit or "1" for unit in self._units if (self.hasComposite or unit!="L_r")]
		axes_units = (" / ( " + " * ".join(axes_units) + " )") if axes_units else ""
		for d in self._diags:
			deposited_quantity = self._myinfo[d]["deposited_quantity"]
			titles[d], val_units = self._make_units(deposited_quantity, self.hasComposite)
			units[d] = val_units + axes_units
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for d in self._diags:
			self._vunits = self._vunits.replace("#"+str(d), "( "+units[d]+" )")
			self._title  = self._title .replace("#"+str(d), titles[d])
		self._vunits = self.units._getUnits(self._vunits)
		if user_axes:
			self._title = "(%s)/(%s)"%(self._title, " x ".join(user_axes))
		if deposited_quantity == "user_function":
			self._title += "/volume"+("  x" if self._vunits else "")
		
		# Set the directory in case of exporting
		self._exportPrefix = self._diagType+"_"+"-".join([str(d) for d in self._diags])
		self._exportDir = self._setExportDir(self._exportPrefix)
		
		# Finish constructor
		self.valid = True
		return kwargs
	
	@staticmethod
	def _make_units(deposited_quantity, hasComposite):
		if   deposited_quantity == "weight":
			title = "Number" + ("" if hasComposite else " density")
			units = "1" if hasComposite else "N_r"
		elif deposited_quantity == "weight_charge":
			title = "Charge" + ("" if hasComposite else " density")
			units = "Q_r" if hasComposite else "N_r * Q_r"
		elif deposited_quantity == "weight_ekin":
			title = "Energy" + ("" if hasComposite else " density")
			units = "K_r" if hasComposite else "N_r * K_r"
		elif deposited_quantity[:15] == "weight_charge_v":
			title = "J"+deposited_quantity[-1] + (" x Volume" if hasComposite else "")
			units = "J_r/N_r" if hasComposite else "J_r"
		elif deposited_quantity[:8] == "weight_p":
			title = "P"+deposited_quantity[8:] + ("" if hasComposite else " density")
			units = "P_r" if hasComposite else "N_r * P_r"
		elif deposited_quantity[:8] == "weight_v":
			title = "Pressure "+deposited_quantity[8]+deposited_quantity[11] + (" x Volume" if hasComposite else "")
			units = "K_r" if hasComposite else "N_r * K_r"
		elif deposited_quantity[:13] == "weight_ekin_v":
			title = "Energy ("+deposited_quantity[-1]+") flux" + (" x Volume" if hasComposite else " density")
			units = "K_r" if hasComposite else "N_r * K_r"
		elif deposited_quantity[:13] == "weight_power": # for radiation spectrum
			title = "Power" + ("" if hasComposite else " density")
			units = "K_r / T_r" if hasComposite else "N_r * K_r / T_r"
		elif deposited_quantity == "user_function":
			title = "user_function"
			units = "1"
		else:
			title = ""
			units = "1"
		return title, units
	
	# Gets info about diagnostic number "diagNumber"
	@staticmethod
	def _getInfo(simulation, diagType, diagNumber):
		info = {}
		for path in simulation._results_path:
			# Open file
			try:
				file = path+simulation._os.sep+diagType+str(diagNumber)+'.h5'
				f = simulation._h5py.File(file, 'r')
			except Exception as e:
				continue
			# get attributes from file
			axes = []
			deposited_quantity = "weight_power" # necessary for radiation spectrum
			time_average = None
			# Parse each attribute
			for name, value in f.attrs.items():
				if name == "deposited_quantity":
					try:
						deposited_quantity = bytes.decode(value)
					except Exception as e:
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
						type = sp[0], min = sp[1], max = sp[2], size = int(sp[3]),
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
					print(diagType+" #"+str(diagNumber)+" in path '"+path+"' is incompatible with the other ones")
					return False
		if not info:
			return False
		return info
	
	# Prints the info obtained by the function "getInfo"
	@staticmethod
	def _printInfo(info):
		if not info:
			return "Error while reading file(s)"
		
		# 1 - diag number, type and list of species
		species = ",".join([str(s) for s in info["species"]])
		title, _ = ParticleBinning._make_units(info["deposited_quantity"], False)
		printedInfo = "#"+str(info["#"])+" - "+title+" of species # "+species+"\n"
		
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
		info = "\n"
		for d in self._diags:
			info += self._printInfo(self._myinfo[d]) + "\n"
		if len(self.operation)>2:
			info += "Operation : "+self.operation+"\n"
		for ax in self._axes:
			if "averageInfo" in ax: info += ax["averageInfo"]+"\n"
			if "subsetInfo" in ax: info += ax["subsetInfo"]+"\n"
		info += self._units_explanation+"\n"
		return info
	
	def _updateAxes(self, timestep):
		uniform = True
		plot_diff = []
		coeff = 1.
		self._finalShape = [[]]*self._naxes
		self._selection = [self._np.s_[:]]*self._naxes
		self._type    = []
		self._shape   = []
		self._centers = []
		self._log     = []
		self._label   = []
		self._units   = []
		i = self._indexOfTime[self._diags[0]][timestep]
		
		for iaxis, axis in enumerate(self._axes):
			axismin = axis["auto_min"][i] if axis["min"]=="auto" else float(axis["min"])
			axismax = axis["auto_max"][i] if axis["max"]=="auto" else float(axis["max"])
			
			# Find the vector of values along the axis
			if axis["log"]:
				edges = self._np.linspace(self._np.log10(axismin), self._np.log10(axismax), axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				edges = 10.**edges
				centers = 10.**(centers[:-1])
			else:
				edges = self._np.linspace(axismin, axismax, axis["size"]+1)
				centers = edges + (edges[1]-edges[0])/2.
				centers = centers[:-1]
			axis["edges"  ] = edges
			axis["centers"] = centers
			
			# if this axis has to be averaged, then select the bounds
			if "average" in axis:
				axis["averageInfo"], self._selection[iaxis], self._finalShape[iaxis] \
					= self._selectRange(axis["average"], centers, axis["type"], axis["units"], "average", axis["edges_included"])
				
				if axis["type"] in ["x","y","z","moving_x"]:
					first_edge = edges[self._selection[iaxis].start or 0]
					last_edge  = edges[(self._selection[iaxis].stop or len(centers))]
					coeff /= last_edge - first_edge
				
				plot_diff.append( self._np.ones((self._finalShape[iaxis],)) )
			
			# if not averaged
			else:
				# If taking a subset of this axis
				if "subset" in axis:
					axis["subsetInfo"], self._selection[iaxis], self._finalShape[iaxis] \
						= self._selectSubset(axis["subset"], centers, axis["type"], axis["units"], "subset")
					
					# If selection is not a slice (meaning only one element) then axis removed from plot
					if type(self._selection[iaxis]) is not slice and axis["type"] in ["x","y","z","moving_x"]:
						first_edge = edges[self._selection[iaxis]]
						last_edge  = edges[self._selection[iaxis]+1]
						coeff /= last_edge - first_edge
				
				# If no subset, or subset has more than 1 point, use this axis in the plot
				if type(self._selection[iaxis]) is slice:
					self._type   .append(axis["type"])
					self._shape  .append(axis["size"])
					self._centers.append(centers[self._selection[iaxis]])
					self._log    .append(axis["log"])
					self._label  .append(axis["type"])
					self._units  .append(axis["units"])
					if axis["type"] == "theta" and self._ndim_particles==3:
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
		
		# If any spatial dimension did not appear, then count it for calculating the correct density
		self._units_explanation = []
		if self._ndim_particles>=1 and not self._spatialaxes["x"]: 
			coeff /= self._ncels[ 0]*self._cell_length[ 0]
			self._units_explanation += ["grid_length[0]"]
		if self._ndim_particles>=2 and not self._spatialaxes["y"]:
			coeff /= self._ncels[ 1]*self._cell_length[ 1]
			self._units_explanation += ["grid_length[1]"]
		if self._ndim_particles==3 and not self._spatialaxes["z"]:
			coeff /= self._ncels[-1]*self._cell_length[-1]
			self._units_explanation += ["grid_length[-1]"]
		if self._units_explanation:
			self._units_explanation = " and by " + " * ".join(self._units_explanation)
		self._units_explanation = "The value in each bin is the sum of the `deposited_quantity` divided by the bin size"\
			+ self._units_explanation
		
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
		if not self.hasComposite:
			self._bsize *= coeff
		
	def _getCenters(self, axis_index, timestep):
		if self.auto_axes:
			self._updateAxes(timestep)
		return self._np.array(self._centers[axis_index])
	
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
					file = path+self._os.sep+self._diagType+str(diagNumber)+'.h5'
					f = self._h5py.File(file, 'r')
				except Exception as e:
					print("Cannot open file "+file)
					return self._np.array([])
				times.update( set(f.keys()) )
				f.close()
			times = [int(t.strip("timestep")) for t in times]
			return self._np.array(times)
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Auto axes require recalculation of bin size and centers
		if self.auto_axes:
			self._updateAxes(t)
			if len(self._shape) > 1:
				# prepare extent for 2d plots
				self._extent = [
					self._xfactor*self._centers[0][0],
					self._xfactor*self._centers[0][-1],
					self._yfactor*self._centers[1][0],
					self._yfactor*self._centers[1][-1]
				]
				if self._log[0]:
					self._extent[0] = self._np.log10(self._extent[0])
					self._extent[1] = self._np.log10(self._extent[1])
				if self._log[1]:
					self._extent[2] = self._np.log10(self._extent[2])
					self._extent[3] = self._np.log10(self._extent[3])
		# Get arrays from all requested diagnostics
		A = {}
		for d in self._diags:
			# find the index of the array corresponding to the requested timestep
			try:
				index = self._indexOfTime[d][t]
			except Exception as e:
				print("Timestep "+str(t)+" not found in this diagnostic")
				return []
			# get data
			B = self._np.empty(self._finalShape)
			try:
				self._h5items[d][index].read_direct(B, source_sel=self._selection) # get array
			except Exception as e:
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
		# Apply the averaging
		for iaxis in range(self._naxes):
			if "average" in self._axes[iaxis]:
				A = self._np.sum(A, axis=iaxis, keepdims=True)
		# remove averaged axes
		A = self._np.squeeze(A)
		# transform if requested
		if callable(self._data_transform): A = self._data_transform(A)
		return A
