from .Diagnostic import Diagnostic
from .._Utils import *


class Field(Diagnostic):
	""" The Field diagnostic of a Smilei simulation"""
	
	def _init(self, field=None, timesteps=None, slice=None, data_log=False, stride=1, **kwargs):
		
		# Open the file
		self._file = self._results_path+'/Fields.h5'
		try:
			f = self._h5py.File(self._file, 'r')
		except:
			self._error = "Diagnostic not loaded: No fields found"
			return
		self._h5items = list(f.values())
		
		# If no field selected, print available fields and leave
		if field is None:
			fields = self.getFields()
			if len(fields)>0:
				self._error += "Printing available fields:\n"
				self._error += "--------------------------\n"
				l = int(len(fields)/3) * 3
				maxlength = str(self._np.max([len(f) for f in fields])+4)
				fields = [('%'+maxlength+'s')%f for f in fields]
				if l>0:
					self._error += '\n'.join([''.join(list(i)) for i in self._np.reshape(fields[:l],(-1,3))])
				self._error += '\n'+''.join(list(fields[l:]))
			else:
				self._error += "No fields found in '"+self._results_path
			return
		
		# Get available times
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
			self._error = "Diagnostic not loaded: No fields found in Fields.h5"
			return
		
		# Get available fields
		fields = self.getFields()
		sortedfields = reversed(sorted(fields, key = len))
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Parse the `field` argument
		self.operation = field
		self._operation = self.operation
		self._fieldname = []
		for f in sortedfields:
			if f in self._operation:
				self._operation = self._re.sub(r"\b"+f+r"\b","C['"+f+"']",self._operation)
				self._fieldname.append(f)
		
		# Check slice is a dict
		if slice is None: slice = {}
		if type(slice) is not dict:
			self._error = "Diagnostic not loaded: Argument `slice` must be a dictionary"
			return
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Get the shape of fields
		fields = list(self._h5items[0].values());
		self._initialShape = fields[0].shape;
		for fd in fields:
			self._initialShape = self._np.min((self._initialShape, fd.shape), axis=0)
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with indices to the data arrays
		self._data = {}
		for i,t in enumerate(self.times):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				self.times = self._selectTimesteps(timesteps, self.times)
			except:
				self._error = "Diagnostic not loaded: Argument `timesteps` must be one or two non-negative integers"
				return
		
		# Need at least one timestep
		if self.times.size < 1:
			self._error = "Diagnostic not loaded: Timesteps not found"
			return
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# Fabricate all axes values
		self._naxes = self._ndim
		self._sliceinfo = {}
		self._finalShape = self._initialShape
		self._slices = [False]*self._ndim
		self._selection = ()
		for iaxis in range(self._naxes):
			centers = self._np.linspace(0., self._initialShape[iaxis]*self._cell_length[iaxis], self._initialShape[iaxis])
			label = {0:"x", 1:"y", 2:"z"}[iaxis]
			axisunits = "L_r"
			
			if label in slice:
				self._slices[iaxis] = True
				# if slice is "all", then all the axis has to be summed
				if slice[label] == "all":
					self._sliceinfo.update({ label:"Sliced for all "+label })
					self._selection += ( self._np.s_[:], )
				# Otherwise, get the slice from the argument `slice`
				else:
					try:
						s = self._np.double(slice[label])
						if s.size>2 or s.size<1: raise
					except:
						self._error = "Diagnostic not loaded: Slice along axis "+label+" should be one or two floats"
						return
					if s.size==1:
						indices = self._np.array([(self._np.abs(centers-s)).argmin()])
					elif s.size==2:
						indices = self._np.nonzero( (centers>=s[0]) * (centers<=s[1]) )[0]
					if indices.size == 0:
						self._error = "Diagnostic not loaded: Slice along "+label+" is out of the box range"
						return None
					if indices.size == 1:
						self._sliceinfo.update({ label:"Sliced at "+label+" = "+str(centers[indices])+" "+axisunits })
						self._selection += ( self._np.s_[indices[0]], )
						self._finalShape[iaxis] = 1
					else:
						self._sliceinfo.update({ label:"Sliced for "+label
							+" from "+str(centers[indices[0]])+" to "+str(centers[indices[-1]])+" "+axisunits })
						self._selection += ( self._np.s_[indices[0]:indices[-1]], )
						self._finalShape[iaxis] = indices[-1] - indices[0]
			else:
				centers = centers[:self._initialShape[iaxis]:stride]
				self._selection += ( self._np.s_[:self._initialShape[iaxis]:stride], )
				self._finalShape[iaxis] = len(centers)
				self._type     .append(label)
				self._shape    .append(self._finalShape[iaxis])
				self._centers  .append(centers)
				self._label    .append(label)
				self._units    .append(axisunits)
				self._log      .append(False)
		
		if len(self._centers) > 2:
			self._error = "Diagnostic not loaded: Cannot plot in "+str(len(self._shape))+"d. You need to 'slice' some axes."
			return
		
		# Build units
		units = {}
		for f in self._fieldname:
			units.update({ f:{"B":"B_r", "E":"E_r", "J":"J_r", "R":"N_r"}[f[0]] })
		# Make total units and title
		self._vunits = self.operation
		self._title  = self.operation
		for f in self._fieldname:
			self._vunits = self._vunits.replace(f, units[f])
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included fields
	def _info(self):
		return "Field diagnostic "+self._title
	
	# get all available fields, sorted by name length
	def getFields(self):
		try:    fields = self._h5items[0].keys() # list of fields
		except: fields = []
		return list(fields)
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		try:    times = [float(a.name[1:]) for a in self._h5items[:-1]]
		except: times = []
		return self._np.double(times)
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self.times:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Get arrays from requested field
		# get data
		index = self._data[t]
		C = {}
		h5item = self._h5items[index]
		for field in self._fieldname: # for each field in operation
			B = self._np.zeros(self._finalShape)
			h5item[field].read_direct(B, source_sel=self._selection) # get array
			C.update({ field:B })
		# Calculate the operation
		A = eval(self._operation)
		# Apply the slicing
		for iaxis in range(self._naxes):
			if self._slices[iaxis]:
				A = self._np.mean(A, axis=iaxis, keepdims=True) # mean over the slice
		# remove sliced axes
		A = self._np.squeeze(A)
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
