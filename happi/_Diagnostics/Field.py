from .Diagnostic import Diagnostic
from .._Utils import *

class Field(Diagnostic):
	"""Class for loading a Field diagnostic"""
	
	def _init(self, diagNumber=None, field=None, timesteps=None, subset=None, average=None, data_log=False, moving=False, **kwargs):
		
		self.moving = moving
		self._subsetinfo = {}
		self.cylindrical = self.namelist.Main.geometry == "AMcylindrical"
		
		# Search available diags
		diag_numbers, diag_names = self.simulation.getDiags("Fields")
		
		# Return directly if no diag number provided
		if diagNumber is None:
			self._error += ["Diagnostic not loaded: diagNumber is not defined"]
			if len(diag_numbers)>0:
				self._error += ["Please choose among: "+", ".join([str(d) for d in diag_numbers])]
				self._error += ["(names: "+", ".join([d for d in diag_names if d])+" )"]
			else:
				self._error += ["(No Field diagnostics existing anyways)"]
			return
		elif type(diagNumber) is str:
			if diagNumber not in diag_names:
				self._error += ["Diagnostic not loaded: no field diagnostic #"+str(diagNumber)+" found"]
				return
			i = diag_names.index( diagNumber )
		else:
			if diagNumber not in diag_numbers:
				self._error += ["Diagnostic not loaded: no field diagnostic #"+str(diagNumber)+" found"]
				return
			i = diag_numbers.index( diagNumber )
		self.diagNumber = diag_numbers[i]
		self.diagName = diag_names[i]
		
		# Open the file(s) and load the data
		self._h5items = {}
		self._fields = []
		for path in self._results_path:
			file = path+self._os.sep+'Fields'+str(self.diagNumber)+'.h5'
			try:
				f = self._h5py.File(file, 'r')
			except:
				self._error += ["Diagnostic not loaded: Could not open '"+file+"'"]
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
			self._is_complex = True
			has_complex = False
			for f in all_fields:
				try:
					if f[:5] in ["Bl_m_","Br_m_","Bt_m_"]:
						fname = f[:4]
						f = f[5:]
						has_complex = True
					elif f[:4] in ["Rho_"]:
						fname = f[:3]
						f = f[4:]
						has_complex = True
					elif f[:3] in ["El_","Er_","Et_","Bl_","Br_","Bt_","Jl_","Jr_","Jt_"]:
						fname = f[:2]
						f = f[3:]
						has_complex = True
					elif f[:10] in ["Env_A_abs_","Env_E_abs_"]:
						fname = f[:9]
						f = f[10:]
						self._is_complex = False
						build3d = None
					elif f[:11] in ["Env_Ex_abs_"]:
						fname = f[:10]
						f = f[11:]
						self._is_complex = False
						build3d = None
					elif f[:8] in ["Env_Chi_"]:
						fname = f[:7]
						f = f[8:]
						self._is_complex = False
						build3d = None
					else:
						raise
					try:
						wordmode, imode = f.split('_')
						species_name = ""
					except:
						ff = f.split('_')
						species_name = "_" + "_".join(ff[:-2])
						wordmode = ff[-2]
						imode = ff[-1]
					if wordmode != "mode":
						raise
					fname += species_name
					if fname not in self._fields:
						self._fields[fname] = []
					self._fields[fname] += [int(imode)]
				except:
					print('WARNING: found unknown field '+f)
					continue
			if has_complex and not self._is_complex:
				self._error += ["Cannot use both envelope and normal fields"]
				return
		
		# If no field selected, print available fields and leave
		if field is None:
			if len(self._fields)>0:
				self._error += ["Error: no field chosen"]
				self._error += ["Printing available fields:"]
				self._error += ["--------------------------"]
				l = int(len(self._fields)/3) * 3
				maxlength = str(self._np.max([len(f) for f in self._fields])+4)
				fields = [('%'+maxlength+'s')%f for f in self._fields]
				if l>0:
					self._error += ['\n'.join([''.join(list(i)) for i in self._np.reshape(fields[:l],(-1,3))])]
				self._error += ['\n'+''.join(list(fields[l:]))]
			else:
				self._error += ["No fields found"]
			return
		
		# Get available times
		self._timesteps = self.getAvailableTimesteps()
		if self._timesteps.size == 0:
			self._error += ["Diagnostic not loaded: No fields found"]
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
		if not self._fieldname:
			self._error += ["String "+self.operation+" does not seem to include any field"]
			return
		
		# Check subset
		if subset is None: subset = {}
		elif type(subset) is not dict:
			self._error += ["Argument `subset` must be a dictionary"]
			return
		
		# Check average
		if average is None: average = {}
		elif type(average) is not dict:
			self._error += ["Argument `average` must be a dictionary"]
			return
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Case of a cylindrical geometry
		# Check whether "theta" or "build3d" option is chosen
		if self.cylindrical :
			theta   = kwargs.pop("theta"  , None)
			build3d = kwargs.pop("build3d", None)
			modes   = kwargs.pop("modes"  , None)
			if (theta is None) == (build3d is None):
				self._error += ["In cylindrical geometry, one (and only one) option `theta` or `build3d` is required"]
				return
			if theta is not None:
				try:
					self._theta = float(theta)
				except:
					self._error += ["Option `theta` must be a number"]
					return
				self._getDataAtTime = self._theta_getDataAtTime
			else:
				try:
					build3d = self._np.array(build3d)
					if build3d.shape != (3,3): raise
				except:
					self._error += ["Option `build3d` must be a list of three lists"]
					return
				self._getDataAtTime = self._build3d_getDataAtTime
			# Test whether "modes" is an int or an iterable on ints
			max_nmode = max([max(self._fields[f]) for f in self._fieldname])+1
			if modes is None:
				self._modes = range(max_nmode)
			else:
				try:
					self._modes = [int(modes)]
				except:
					try:
						self._modes = [int(imode) for imode in modes]
					except:
						self._error += ["Option `modes` must be a number or an iterable on numbers"]
						return
				for imode in self._modes:
					if imode >= max_nmode:
						self._error += ["Option `modes` cannot contain modes larger than "+str(max_nmode-1)]
						return
			# Warning if some modes are not available
			for f in self._fieldname:
				unavailable_modes = [str(i) for i in set(self._modes) - set(self._fields[f])]
				if unavailable_modes:
					print("WARNING: field "+f+" does not have mode(s): "+','.join(unavailable_modes))
		
		# Get the shape of fields
		fields = [f for f in self._h5items[0].values() if f]
		self._raw_shape = fields[0].shape
		self._initialShape = fields[0].shape
		for fd in fields:
			self._initialShape = self._np.min((self._initialShape, fd.shape), axis=0)
		self._offset    = fields[0].attrs['gridGlobalOffset']
		self._spacing   = fields[0].attrs['gridSpacing']
		axis_start = self._offset
		axis_stop  = self._offset + (self._initialShape-0.5)*self._spacing
		axis_step  = self._spacing
		if self.cylindrical :
			if build3d is not None:
				self._initialShape = [int(self._np.ceil( (s[1]-s[0])/float(s[2]) )) for s in build3d]
				axis_start = build3d[:,0]
				axis_stop  = build3d[:,1]
				axis_step  = build3d[:,2]
			else:
				if self._is_complex:
					self._initialShape[1] /= 2
				axis_stop = self._offset + (self._initialShape-0.5)*self._spacing
		
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
				self._error += ["Argument `timesteps` must be one or two non-negative integers"]
				return
		
		# Need at least one timestep
		if self._timesteps.size < 1:
			self._error += ["Timesteps not found"]
			return
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		self._naxes = len(self._initialShape)
		self._finalShape = self._np.copy(self._initialShape)
		self._averages = [False]*self._naxes
		self._selection = [self._np.s_[:]]*self._naxes
		self._offset  = fields[0].attrs['gridGlobalOffset']
		self._spacing = fields[0].attrs['gridSpacing']
		axis_name = "xyz" if not self.cylindrical or not self._is_complex or build3d is not None else "xr"
		for iaxis in range(self._naxes):
			centers = self._np.arange(axis_start[iaxis], axis_stop[iaxis], axis_step[iaxis])
			label = axis_name[iaxis]
			axisunits = "L_r"
			
			# If averaging over this axis
			if label in average:
				if label in subset:
					self._error += ["`subset` not possible on the same axes as `average`"]
					return
				
				self._averages[iaxis] = True
				
				try:
					self._subsetinfo[label], self._selection[iaxis], self._finalShape[iaxis] \
						= self._selectRange(average[label], centers, label, axisunits, "average")
				except Exception as e:
					if not self._error:
						self._error += ["Error handling average:"]
						self._error += [str(e)]
					return
			# Otherwise
			else:
				# If taking a subset of this axis
				if label in subset:
					try:
						self._subsetinfo[label], self._selection[iaxis], self._finalShape[iaxis] \
							= self._selectSubset(subset[label], centers, label, axisunits, "subset")
					except Exception as e:
						if not self._error:
							self._error += ["Error handling subset:"]
							self._error += [str(e)]
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
			complex_factor = 2 #cylindrical default is complex
			if not self._is_complex :
				complex_factor = 1
			if type(self._selection[1]) is slice:
				self._complex_selection_real[1] = slice(
					None if self._selection[1].start is None else self._selection[1].start*complex_factor,
					None if self._selection[1].stop  is None else self._selection[1].stop *complex_factor,
					(self._selection[1].step  or 1)*complex_factor
				)
				self._complex_selection_imag[1] = slice(
					(self._selection[1].start or 0)*complex_factor + 1,
					None if self._selection[1].stop is None else self._selection[1].stop*complex_factor + 1,
					(self._selection[1].step  or 1)*complex_factor
				)
			else:
				self._complex_selection_real[1] = self._selection[1]*complex_factor
				self._complex_selection_imag[1] = self._selection[1]*complex_factor + 1
			self._complex_selection_real = tuple(self._complex_selection_real)
			self._complex_selection_imag = tuple(self._complex_selection_imag)
			
			# In the case of "build3d", prepare some data for the 3D construction
			if build3d is not None:
				# Calculate the raw data positions
				self._raw_positions = (
					self._np.arange(self._offset[0], self._offset[0] + (self._raw_shape[0]  -0.5)*self._spacing[0], self._spacing[0] ),
					self._np.arange(self._offset[1], self._offset[1] + (self._raw_shape[1]/complex_factor-0.5)*self._spacing[1], self._spacing[1] ),
				)
				# Calculate the positions of points in the final box
				x = self._np.arange(*build3d[0])
				y = self._np.arange(*build3d[1])
				z = self._np.arange(*build3d[2])
				if len(x)==0 or len(y)==0 or len(z)==0:
					self._error += ["Error: The array shape to be constructed seems to be empty"]
					return
				y2, z2 = self._np.meshgrid(y,z)
				r2 = self._np.sqrt(y2**2 + z2**2)
				self._theta = self._np.arctan2(z2, y2)
				del y2, z2
				r3 = self._np.tile(r2, (len(x), 1, 1))
				self._theta = self._np.tile(self._theta, (len(x), 1, 1))[self._selection]
				x3 = self._np.rollaxis(self._np.tile(x,(len(y),len(z),1)),2)
				self._xr = self._np.stack((x3, r3), axis=-1)
				del x3, r3, r2
		
		# Build units
		units = {}
		for f in self._fieldname:
			units.update({ f:{"B":"B_r", "E":"E_r", "J":"J_r", "R":"N_r"}[f[0]] })
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for f in self._fieldname:
			self._vunits = self._vunits.replace(f, units[f])
		self._vunits = self.units._getUnits(self._vunits)
		
		# Set the directory in case of exporting
		self._exportPrefix = "Field"+str(self.diagNumber)+"_"+"".join(self._fieldname)
		self._exportDir = self._setExportDir(self._exportPrefix)
		
		# Finish constructor
		self.valid = True
		return kwargs
	
	# Method to print info on included fields
	def _info(self):
		s = "Field diagnostic #"+str(self.diagNumber)+": "+self._title
		tavg = self.namelist.DiagFields[self.diagNumber].time_average
		if tavg > 1:
			s += "\n\tTime_average: " + str(tavg) + " timesteps"
		if any(self._offset > 0.):
			s += "\n\tGrid offset: " + ", ".join([str(a) for a in self._offset])
		if any(self._spacing != self._cell_length):
			s += "\n\tGrid spacing: " + ", ".join([str(a) for a in self._spacing])
		for l in self._subsetinfo:
			s += "\n\t"+self._subsetinfo[l]
		return s
	
	# get all available fields, sorted by name length
	def getFields(self):
		return self._fields
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:    times = [float(a.name[6:]) for a in self._h5items]
		except: times = []
		return self._np.double(times)
	
	# get the value of x_moved for a requested timestep
	def getXmoved(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# get h5 iteration group
		index = self._data[t]
		h5item = self._h5items[index]
		return h5item.attrs["x_moved"] if "x_moved" in h5item.attrs else 0.
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Get arrays from requested field
		index = self._data[t]
		C = {}
		h5item = self._h5items[index]
		
		# Handle moving window
		if self.moving and "x_moved" in h5item.attrs and 'x' in self._type:
			self._xoffset = h5item.attrs["x_moved"]
			if self.dim>1 and hasattr(self,"_extent"):
				self._extent[0] = self._xfactor*(self._xoffset + self._centers[0][ 0])
				self._extent[1] = self._xfactor*(self._xoffset + self._centers[0][-1])
		
		 # for each field in operation, obtain the data
		for field in self._fieldname:
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
		for field in self._fieldname: # for each field in operation
			available_modes = self._fields[field]
			F = self._np.zeros(self._finalShape)
			f = field + "_mode_"
			for imode in self._modes:
				if imode not in available_modes: continue
				B_real = self._np.empty(self._finalShape)
				B_imag = self._np.empty(self._finalShape)
				try:
					h5item[f+str(imode)].read_direct(B_real, source_sel=self._complex_selection_real)
					if imode > 0:
						h5item[f+str(imode)].read_direct(B_imag, source_sel=self._complex_selection_imag)
				except:
					B_real = self._np.squeeze(B_real)
					h5item[f+str(imode)].read_direct(B_real, source_sel=self._complex_selection_real)
					B_real = self._np.reshape(B_real, self._finalShape)
					if imode > 0:
						B_imag = self._np.squeeze(B_imag)
						h5item[f+str(imode)].read_direct(B_imag, source_sel=self._complex_selection_imag)
						B_imag = self._np.reshape(B_imag, self._finalShape)
				F += (self._np.cos(imode*self._theta)) * B_real
				if imode > 0:
					F += (self._np.sin(imode*self._theta)) * B_imag
				
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
		return A
	
	# Method to obtain the data only
	# Specific to cylindrical geometry, to reconstruct a 3d box
	def _build3d_getDataAtTime(self, t):
		from scipy.interpolate import RegularGridInterpolator
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		
		# Get arrays from requested field
		index = self._data[t]
		C = {}
		h5item = self._h5items[index]
		step = 2 if self._is_complex else 1
		for field in self._fieldname: # for each field in operation
			available_modes = self._fields[field]
			F = self._np.zeros(self._finalShape)
			f = field + "_mode_"
			for imode in self._modes:
				if imode not in available_modes: continue
				B = self._np.empty(self._raw_shape)
				try:
					h5item[f+str(imode)].read_direct(B)
				except:
					B = self._np.squeeze(B)
					h5item[f+str(imode)].read_direct(B)
					B = self._np.reshape(B, self._raw_shape)
				B_real = RegularGridInterpolator(self._raw_positions, B[:,0::step], bounds_error=False, fill_value=0.)(self._xr)
				F += self._np.cos(imode*self._theta) * B_real[self._selection]
				if imode > 0.:
					B_imag = RegularGridInterpolator(self._raw_positions, B[:,1::step], bounds_error=False, fill_value=0.)(self._xr)
					F += self._np.sin(imode*self._theta) * B_imag[self._selection]
				
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
		return A
