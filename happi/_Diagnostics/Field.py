from .Diagnostic import Diagnostic
from .._Utils import *

class Field(Diagnostic):
	"""Class for loading a Field diagnostic"""
	
	def _init(self, diagNumber=None, field=None, timesteps=None, subset=None, average=None, data_log=False, **kwargs):
		
		self.cylindrical = self.namelist.Main.geometry == "3drz"
		
		# Search available diags
		diags = self.getDiags()
		
		# Return directly if no diag number provided
		if diagNumber is None:
			self._error += "Diagnostic not loaded: diagNumber is not defined\n"
			if len(diags)>0:
				self._error += "Please choose among: "+", ".join([str(d) for d in diags])
			else:
				self._error += "(No Field diagnostics existing anyways)"
			return
		else:
			self.diagNumber = diagNumber
			if diagNumber not in diags:
				self._error = "Diagnostic not loaded: no field diagnostic #"+str(diagNumber)+" found"
				return
		
		# Open the file(s) and load the data
		self._h5items = {}
		self._fields = []
		for path in self._results_path:
			file = path+self._os.sep+'Fields'+str(diagNumber)+'.h5'
			try:
				f = self._h5py.File(file, 'r')
			except:
				self._error = "Diagnostic not loaded: Could not open '"+file+"'"
				return
			self._h5items.update( dict(f["data"]) )
			# Select only the fields that are common to all simulations
			values = f["data"].values()
			if len(values)==0:
				self._fields = []
			elif len(self._fields)==0:
				self._fields = list(next(iter(values)).keys())
			else:
				self._fields = [f for f in next(iter(values)).keys() if f in self._fields]
		# Remove "tmp" dataset
		if "tmp" in self._h5items: del self._h5items["tmp"]
		# Converted to ordered list
		self._h5items = sorted(self._h5items.values(), key=lambda x:int(x.name[6:]))
		
		# Case of a cylindrical geometry
		# Build the list of fields that can be reconstructed
		if self.cylindrical:
			all_fields = list(self._fields)
			self._fields = {}
			for f in ["Ex","Er","Et","Bx","Br","Bt","Jx","Jr","Jt","Rho"]:
				imode = 0
				while f+"_mode_"+str(imode) in all_fields:
					self._fields[f] = imode
					imode += 1
		
		# If no field selected, print available fields and leave
		if field is None:
			if len(self._fields)>0:
				self._error += "Printing available fields:\n"
				self._error += "--------------------------\n"
				l = int(len(self._fields)/3) * 3
				maxlength = str(self._np.max([len(f) for f in self._fields])+4)
				fields = [('%'+maxlength+'s')%f for f in self._fields]
				if l>0:
					self._error += '\n'.join([''.join(list(i)) for i in self._np.reshape(fields[:l],(-1,3))])
				self._error += '\n'+''.join(list(fields[l:]))
			else:
				self._error += "No fields found"
			return
		
		# Get available times
		self._timesteps = self.getAvailableTimesteps()
		if self._timesteps.size == 0:
			self._error = "Diagnostic not loaded: No fields found"
			return
		
		# Get available fields
		sortedfields = reversed(sorted(self._fields, key = len))
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Parse the `field` argument
		self.operation = field
		self._operation = self.operation
		self._fieldname = []
		for f in sortedfields:
			if self._re.search(r"\b"+f+r"\b",self._operation):
				self._operation = self._re.sub(r"\b"+f+r"\b","C['"+f+"']",self._operation)
				self._fieldname.append(f)
		
		# Check subset
		if subset is None: subset = {}
		elif type(subset) is not dict:
			self._error = "Argument `subset` must be a dictionary"
			return
		
		# Check average
		if average is None: average = {}
		elif type(average) is not dict:
			self._error = "Argument `average` must be a dictionary"
			return
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Case of a cylindrical geometry
		# Check whether "theta" or "build3d" option is chosen
		if self.cylindrical:
			self._theta   = kwargs.pop("theta"  , None)
			self._build3d = kwargs.pop("build3d", None)
			if (self._theta is None) == (self._build3d is None):
				self._error = "In cylindrical geometry, one (and only one) option `theta` or `build3d` is required"
				return
			self._mode_theta = self._theta != None
			if self._mode_theta:
				try:
					self._theta = float(self._theta)
				except:
					self._error = "Option `theta` must be a number"
					return
				self._getDataAtTime = self._theta_getDataAtTime
			else:
				try:
					self._build3d = self._np.array(self._build3d)
					if self._build3d.shape != (3,3): raise
				except:
					self._error = "Option `build3d` must be a list of three lists"
					return
				self._getDataAtTime = self._build3d_getDataAtTime
		
		# Get the shape of fields
		fields = [f for f in self._h5items[0].values() if f]
		self._initialShape = fields[0].shape
		for fd in fields:
			self._initialShape = self._np.min((self._initialShape, fd.shape), axis=0)
		axis_start = [0.]*3
		axis_stop  = (self._initialShape-0.5)*self._cell_length
		axis_step  = self._cell_length
		if self.cylindrical:
			if self._mode_theta:
				self._initialShape[1] /= 2
				axis_stop = (self._initialShape-0.5)*self._cell_length
			else:
				self._initialShape = [int((s[1]-s[0])/s[2])+1 for s in self._build3d]
				axis_start = self._build3d[:,0]
				axis_stop  = self._build3d[:,1]
				axis_step  = self._build3d[:,2]
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with indices to the data arrays
		self._data = {}
		for i,t in enumerate(self._timesteps):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				self._timesteps = self._selectTimesteps(timesteps, self._timesteps)
			except:
				self._error = "Argument `timesteps` must be one or two non-negative integers"
				return
		
		# Need at least one timestep
		if self._timesteps.size < 1:
			self._error = "Timesteps not found"
			return
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		self._naxes = len(self._initialShape)
		self._subsetinfo = {}
		self._finalShape = self._np.copy(self._initialShape)
		self._averages = [False]*self._ndim
		self._selection = [self._np.s_[:]]*self._ndim
		axis_name = "xr" if self.cylindrical and self._mode_theta else "xyz"
		for iaxis in range(self._naxes):
			centers = self._np.arange(axis_start[iaxis], axis_stop[iaxis], axis_step[iaxis])
			label = axis_name[iaxis]
			axisunits = "L_r"
			
			# If averaging over this axis
			if label in average:
				if label in subset:
					self._error = "`subset` not possible on the same axes as `average`"
					return
				
				self._averages[iaxis] = True
				
				try:
					self._subsetinfo[label], self._selection[iaxis], self._finalShape[iaxis] \
						= self._selectRange(average[label], centers, label, axisunits, "average")
				except:
					return
			# Otherwise
			else:
				# If taking a subset of this axis
				if label in subset:
					try:
						self._subsetinfo[label], self._selection[iaxis], self._finalShape[iaxis] \
							= self._selectSubset(subset[label], centers, label, axisunits, "subset")
					except:
						return
				# If subset has more than 1 point (or no subset), use this axis in the plot
				if type(self._selection[iaxis]) is slice:
					self._type     .append(label)
					self._shape    .append(self._finalShape[iaxis])
					self._centers  .append(centers[self._selection[iaxis]])
					self._label    .append(label)
					self._units    .append(axisunits)
					self._log      .append(False)
		
		self._selection = tuple(self._selection)
		
		# Special selection in cylindrical geometry due to complex numbers
		if self.cylindrical:
			self._complex_selection_real = list(self._selection)
			self._complex_selection_imag = list(self._selection)
			if type(self._selection[1]) is slice:
				self._complex_selection_real[1] = slice(
					None if self._selection[1].start is None else self._selection[1].start*2,
					None if self._selection[1].stop  is None else self._selection[1].stop *2,
					(self._selection[1].step  or 1)*2
				)
				self._complex_selection_imag[1] = slice(
					(self._selection[1].start or 0)*2 + 1,
					None if self._selection[1].stop is None else self._selection[1].stop*2 + 1,
					(self._selection[1].step  or 1)*2
				)
			else:
				self._complex_selection_real[1] = self._selection[1]*2
				self._complex_selection_imag[1] = self._selection[1]*2 + 1
			self._complex_selection_real = tuple(self._complex_selection_real)
			self._complex_selection_imag = tuple(self._complex_selection_imag)
		
		# Build units
		units = {}
		for f in self._fieldname:
			units.update({ f:{"B":"B_r", "E":"E_r", "J":"J_r", "R":"N_r"}[f[0]] })
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for f in self._fieldname:
			self._vunits = self._vunits.replace(f, units[f])
		
		# Set the directory in case of exporting
		self._exportPrefix = "Field"+str(diagNumber)+"_"+"".join(self._fieldname)
		self._exportDir = self._setExportDir(self._exportPrefix)
		
		# Finish constructor
		self.valid = True
		return kwargs
	
	# Method to print info on included fields
	def _info(self):
		return "Field diagnostic "+self._title
	
	# get all available field diagnostics
	def getDiags(self):
		diags = []
		for path in self._results_path:
			files = self._glob(path+self._os.sep+'Fields*.h5')
			if len(files)==0:
				self._error = "No fields found in '"+path+"'"
				return []
			diagNumbers = [ int(self._re.findall("Fields([0-9]+).h5$",file)[0]) for file in files ]
			if diags == []: diags = diagNumbers
			else          : diags = [ d for d in diags if d in diagNumbers ]
		return diags
	
	# get all available fields, sorted by name length
	def getFields(self):
		return self._fields
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:    times = [float(a.name[6:]) for a in self._h5items]
		except: times = []
		return self._np.double(times)
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Get arrays from requested field
		# get data
		index = self._data[t]
		C = {}
		h5item = self._h5items[index]
		for field in self._fieldname: # for each field in operation
			B = self._np.empty(self._finalShape)
			try:
				h5item[field].read_direct(B, source_sel=self._selection) # get array
			except:
				B = self._np.squeeze(B)
				h5item[field].read_direct(B, source_sel=self._selection) # get array
				B = self._np.reshape(B, self._finalShape)
			C.update({ field:B })
		# Calculate the operation
		A = eval(self._operation)
		# Apply the averaging
		A = self._np.reshape(A,self._finalShape)
		for iaxis in range(self._naxes):
			if self._averages[iaxis]:
				A = self._np.mean(A, axis=iaxis, keepdims=True)
		# remove averaged axes
		A = self._np.squeeze(A)
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
	
	# Method to obtain the data only
	# Specific to cylindrical geometry, to reconstruct the plane at some theta
	def _theta_getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Get arrays from requested field
		# get data
		index = self._data[t]
		C = {}
		h5item = self._h5items[index]
		self._finalShape
		for field in self._fieldname: # for each field in operation
			nmode = self._fields[field]
			F = self._np.zeros(self._finalShape)
			f = field + "_mode_"
			for imode in range(nmode):
				B_real = self._np.empty(self._finalShape)
				B_imag = self._np.empty(self._finalShape)
				try:
					h5item[f+str(imode)].read_direct(B_real, source_sel=self._complex_selection_real)
					h5item[f+str(imode)].read_direct(B_imag, source_sel=self._complex_selection_imag)
				except:
					B_real = self._np.squeeze(B_real)
					B_imag = self._np.squeeze(B_imag)
					h5item[f+str(imode)].read_direct(B_real, source_sel=self._complex_selection_real)
					h5item[f+str(imode)].read_direct(B_imag, source_sel=self._complex_selection_imag)
					B_real = self._np.reshape(B_real, self._finalShape)
					B_imag = self._np.reshape(B_imag, self._finalShape)
				F += (2.*self._np.cos(imode*self._theta)) * B_real + (2.*self._np.sin(imode*self._theta)) * B_imag
				
			C.update({ field:F })
		
		# Calculate the operation
		A = eval(self._operation)
		# Apply the averaging
		A = self._np.reshape(A,self._finalShape)
		for iaxis in range(self._naxes):
			if self._averages[iaxis]:
				A = self._np.mean(A, axis=iaxis, keepdims=True)
		# remove averaged axes
		A = self._np.squeeze(A)
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
		
	# Method to obtain the data only
	# Specific to cylindrical geometry, to reconstruct a 3d box
	def _build3d_getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Calculate the positions of points in the box
		x = self._np.arange(*self._build3d[0])
		y = self._np.arange(*self._build3d[1])
		z = self._np.arange(*self._build3d[2])
		yy, zz = self._np.meshgrid(x,y)
		r = self._np.sqrt(yy**2 + zz**2)
		theta = self._np.atan2(zz, yy)
		
		# Get arrays from requested field
		index = self._data[t]
		C = {}
		h5item = self._h5items[index]
		self._finalShape
		for field in self._fieldname: # for each field in operation
			nmode = self._fields[field]
			F = self._np.zeros(self._finalShape)
			f = field + "_mode_"
			for imode in range(nmode):
				B_real = self._np.empty(self._finalShape)
				B_imag = self._np.empty(self._finalShape)
				try:
					h5item[f+str(imode)].read_direct(B_real, source_sel=self._complex_selection_real)
					h5item[f+str(imode)].read_direct(B_imag, source_sel=self._complex_selection_imag)
				except:
					B_real = self._np.squeeze(B_real)
					B_imag = self._np.squeeze(B_imag)
					h5item[f+str(imode)].read_direct(B_real, source_sel=self._complex_selection_real)
					h5item[f+str(imode)].read_direct(B_imag, source_sel=self._complex_selection_imag)
					B_real = self._np.reshape(B_real, self._finalShape)
					B_imag = self._np.reshape(B_imag, self._finalShape)
				F += (2.*self._np.cos(imode*self._theta)) * B_real + (2.*self._np.sin(imode*self._theta)) * B_imag
				
			C.update({ field:F })
		
		# Calculate the operation
		A = eval(self._operation)
		# Apply the averaging
		A = self._np.reshape(A,self._finalShape)
		for iaxis in range(self._naxes):
			if self._averages[iaxis]:
				A = self._np.mean(A, axis=iaxis, keepdims=True)
		# remove averaged axes
		A = self._np.squeeze(A)
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
