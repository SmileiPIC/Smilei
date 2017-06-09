from .Diagnostic import Diagnostic
from .._Utils import *


class Field(Diagnostic):
	""" The Field diagnostic of a Smilei simulation"""
	
	def _init(self, diagNumber=None, field=None, timesteps=None, average=None, data_log=False, stride=1, **kwargs):
		
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
				self._fields = values[0].keys()
			else:
				self._fields = [f for f in values[0].keys() if f in self._fields]
		# Remove "tmp" dataset
		if "tmp" in self._h5items: del self._h5items["tmp"]
		# Converted to ordered list
		self._h5items = sorted(self._h5items.values(), key=lambda x:int(x.name[6:]))
		
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
		self.times = self.getAvailableTimesteps()
		if self.times.size == 0:
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
		
		# Check average is a dict
		if average is None: average = {}
		if type(average) is not dict:
			self._error = "Diagnostic not loaded: Argument `average` must be a dictionary"
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
		self._naxes = self._ndim
		self._averageinfo = {}
		self._finalShape = self._np.copy(self._initialShape)
		self._averages = [False]*self._ndim
		self._selection = []
		for iaxis in range(self._naxes):
			centers = self._np.linspace(0., (self._initialShape[iaxis]-1)*self._cell_length[iaxis], self._initialShape[iaxis])
			label = {0:"x", 1:"y", 2:"z"}[iaxis]
			axisunits = "L_r"
			self._selection += [ self._np.s_[:self._initialShape[iaxis]:stride] ]
			
			if label in average:
				self._averages[iaxis] = True
				
				self._averageinfo[label], self._selection[iaxis], self._finalShape[iaxis] \
					= self._selectRange(average[label], centers, label, axisunits, "average")
				
			else:
				centers = centers[:self._initialShape[iaxis]:stride]
				self._finalShape[iaxis] = len(centers)
				self._type     .append(label)
				self._shape    .append(self._finalShape[iaxis])
				self._centers  .append(centers)
				self._label    .append(label)
				self._units    .append(axisunits)
				self._log      .append(False)
		self._selection = tuple(self._selection)
		
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
				self._error = "Diagnostic not loaded: No fields found in '"+path+"'"
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
		if t not in self.times:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Get arrays from requested field
		# get data
		index = self._data[t]
		C = {}
		h5item = self._h5items[index]
		for field in self._fieldname: # for each field in operation
			B = self._np.squeeze(self._np.zeros(self._finalShape))
			h5item[field].read_direct(B, source_sel=self._selection) # get array
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
	
	# Convert to XDMF format for ParaView
	def toXDMF(self):
		
		# Calculate a few things
		ndim = self._ndim
		shape = self._h5items[0].values()[0].shape
		cell_length = list(self._cell_length)
		if ndim == 1:
			ndim = 2
			shape += (1,)
			cell_length += [1.]
		shapestr = " ".join([str(a) for a in shape])
		#field_axis = "xyz".index(field[1])
		#magnetic_field = (field[0]=="B")
		#origin = [ str(((field_axis==i)^magnetic_field)*(-0.5)) for i in range(ndim) ]
		origin = [ 0. for i in range(ndim) ]
		try:    requestedfields = self._fieldname
		except: requestedfields = False
		
		self._mkdir(self._exportDir)
		fileprefix = self._exportDir+sep+"Fields"+str(self.diagNumber)
		if requestedfields: fileprefix += "".join(requestedfields)
		
		# Make the XDMF for usual time collections
		with open(fileprefix+".xmf",'w') as f:
			f.write('<?xml version="1.0" ?>\n')
			f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
			f.write('<Xdmf Version="3.0">\n')
			f.write('	<Domain>\n')
			topology='				<Topology Name="Fields topology" TopologyType="'+str(ndim)+'DCoRectMesh" Dimensions="'+shapestr+'"/>\n'
			geometry=('				<Geometry Name="Fields geometry" GeometryType="ORIGIN_'+"".join(["DX","DY","DZ"][0:ndim])+'">\n'
				+'					<DataItem Format="XML" NumberType="Float" Dimensions="'+str(ndim)+'">'+" ".join([str(o) for o in origin])+'</DataItem>\n'
				+'					<DataItem Format="XML" NumberType="Float" Dimensions="'+str(ndim)+'">'+" ".join([str(o) for o in cell_length])+'</DataItem>\n'
				+'				</Geometry>\n')
			XYZ = self._np.meshgrid(*[[origin[dim]+i*self._cell_length[dim] for i in range(shape[dim])] for dim in range(self._ndim)])
			for dim in range(self._ndim):
				f.write('		<DataItem Name="Space'+"XYZ"[dim]+'" ItemType="Uniform" NumberType="Float" Dimensions="'+shapestr+'" Format="XML">\n')
				f.write('			'+" ".join([str(i) for i in XYZ[dim].flatten()])+'\n')
				f.write('		</DataItem>\n')
			f.write('		<Grid GridType="Collection" CollectionType="Temporal">\n')
			for item in self._h5items:
				f.write('			<Grid Name="Timestep_'+str(item.name[6:])+'" GridType="Uniform">\n')
				f.write('				<Time Value="'+str(float(item.name[6:])*self.timestep)+'"/>\n')
				f.write(topology)
				f.write(geometry)
				for dim in range(self._ndim):
					f.write('				<Attribute Name="'+"XYZ"[dim]+'" Center="Node" AttributeType="Scalar">\n')
					f.write('					<DataItem ItemType="Uniform" NumberType="Float" Dimensions="'+shapestr+'" Format="XML" Reference="XML">/Xdmf/Domain/DataItem[@Name="Space'+"XYZ"[dim]+'"]</DataItem>\n')
					f.write('				</Attribute>\n')
				for field in item.values():
					if requestedfields and field.name not in requestedfields: continue
					location = self._os.path.abspath(item.file.filename)+':'+field.name
					f.write('				<Attribute Name="'+self._os.path.basename(field.name)+'" Center="Node" AttributeType="Scalar">\n')
					f.write('					<DataItem ItemType="Uniform" NumberType="Float" Precision="8" Dimensions="'+shapestr+'" Format="HDF">'+location+'</DataItem>\n')
					f.write('				</Attribute>\n')
				f.write('			</Grid>\n')
			f.write('		</Grid>\n')
			f.write('	</Domain>\n')
			f.write('</Xdmf>\n')
		
		# Make the XDMF for time streak
		if self._ndim < 3:
			ndim = self._ndim + 1
			shape = self._h5items[0].values()[0].shape + (len(self._h5items),)
			shapestr = " ".join([str(a) for a in shape])
			origin = [ 0. for i in range(ndim) ]
			cell_length = list(self._cell_length) + [self.timestep]
			axes = "XYZ"[0:self._ndim] + "T"
			fields = self._h5items[0].keys()
			if requestedfields:
				fields = [F for F in fields if F in requestedfields]
			with open(fileprefix+"_streak.xmf",'w') as f:
				f.write('<?xml version="1.0" ?>\n')
				f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
				f.write('<Xdmf Version="3.0">\n')
				f.write('	<Domain>\n')
				f.write('		<Grid GridType="Uniform">\n')
				f.write('			<Topology Name="Fields topology" TopologyType="'+str(ndim)+'DCoRectMesh" Dimensions="'+shapestr+'"/>\n')
				f.write('			<Geometry Name="Fields geometry" GeometryType="ORIGIN_'+"".join(["DX","DY","DZ"][0:ndim])+'">\n')
				f.write('				<DataItem Format="XML" NumberType="Float" Dimensions="'+str(ndim)+'">'+" ".join([str(o) for o in origin])+'</DataItem>\n')
				f.write('				<DataItem Format="XML" NumberType="Float" Dimensions="'+str(ndim)+'">'+" ".join([str(o) for o in cell_length])+'</DataItem>\n')
				f.write('			</Geometry>\n')
				XYZ = self._np.meshgrid(*[[origin[dim]+i*cell_length[dim] for i in range(shape[dim])] for dim in reversed(range(ndim))])
				for dim in range(ndim):
					f.write('			<Attribute Name="'+axes[ndim-dim-1]+'" Center="Node" AttributeType="Scalar">\n')
					f.write('				<DataItem ItemType="Uniform" NumberType="Float" Precision="8" Dimensions="'+shapestr+'" Format="XML">\n')
					f.write('					'+" ".join([str(i) for i in XYZ[ndim-dim-1].flatten()])+'\n')
					f.write('				</DataItem>\n')
					f.write('			</Attribute>\n')
				for field in fields:
					f.write('			<Attribute Name="'+self._os.path.basename(field)+'" Center="Node" AttributeType="Scalar">\n')
					f.write('				<DataItem ItemType="Function" Function="'+"|".join(["$"+str(i) for i in range(len(self._h5items))])+'" Dimensions="'+shapestr+'">\n')
					for item in self._h5items:
						location = self._os.path.abspath(item.file.filename)+':'+item.name+"/"+field
						f.write('					<DataItem ItemType="Uniform" NumberType="Float" Precision="8" Dimensions="'+ " ".join([str(a) for a in shape[:-1]])+'" Format="HDF">'+location+'</DataItem>\n')
					f.write('				</DataItem>\n')
					f.write('			</Attribute>\n')
				f.write('		</Grid>\n')
				f.write('	</Domain>\n')
				f.write('</Xdmf>\n')
