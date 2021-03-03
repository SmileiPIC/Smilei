from .Diagnostic import Diagnostic
from .._Utils import *

class Scalar(Diagnostic):
	"""Class for loading a Scalar diagnostic"""
	
	def _init(self, scalar=None, timesteps=None, data_log=False, data_transform=None, **kwargs):
		# Get available scalars
		scalars = self.getScalars()
		
		# If no scalar chosen, only print the available scalars
		if scalar is None:
			if len(scalars)>0:
				self._error += ["Error: no scalar chosen"]
				self._error += ["Printing available scalars:"]
				self._error += ["---------------------------"]
				l = [""]
				for s in scalars:
					if len(s)>4 and s[:2]!=l[-1][:2] and s[-2:]!=l[-1][-2:]:
						if l!=[""]: self._error += ["\t".join(l)]
						l = []
					l.append(s)
				if l!=[""]: self._error += ["\t".join(l)]
			else:
				self._error += ["No scalars found"]
			return
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Find which scalar is requested
		if scalar not in scalars:
			fs = list(filter(lambda x:scalar in x, scalars))
			if len(fs)==0:
				self._error += ["No scalar `"+scalar+"` found"]
				return
			if len(fs)>1:
				self._error += ["Several scalars match: "+(' '.join(fs))]
				self._error += ["Please be more specific and retry."]
				return
			scalar = fs[0]
		self._scalarname = scalar
		
		# Put data_log as object's variable
		self._data_log = data_log
		self._data_transform = data_transform
		
		# Already get the data from the file
		# Loop file line by line
		self._alltimesteps = []
		self._values = []
		times_values = {}
		for path in self._results_path:
			with open(path+'/scalars.txt') as f:
				for line in f:
					line = line.strip()
					if line[0]!="#": break
					prevline = line
				scalars = prevline[1:].strip().split() # list of scalars
				scalarindex = scalars.index(scalar) # index of the requested scalar
				line = str(line.strip()).split()
				times_values[ int( self._np.round(float(line[0]) / float(self.timestep)) ) ] = float(line[scalarindex])
				for line in f:
					line = str(line.strip()).split()
					times_values[ int( self._np.round(float(line[0]) / float(self.timestep)) ) ] = float(line[scalarindex])
		self._alltimesteps  = self._np.array(sorted(times_values.keys()))
		self._values = self._np.array([times_values[k] for k in self._alltimesteps])
		self._timesteps = self._np.copy(self._alltimesteps)
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "_data" dictionary with the index to each time
		self._data = {}
		for i,t in enumerate(self._timesteps):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				self._timesteps = self._selectTimesteps(timesteps, self._timesteps)
			except Exception as e:
				self._error += ["Argument `timesteps` must be one or two non-negative integers"]
				return
		
		# Need at least one timestep
		if self._timesteps.size < 1:
			self._error += ["Timesteps not found"]
			return
		
		
		# 3 - Build units
		# -------------------------------------------------------------------
		self._vunits = "??"
		if   self._scalarname == "time":
			self._vunits = "T_r"
		elif self._scalarname == "Ubal_norm":
			self._vunits = ""
		else:
			self._vunits = {
				"U":"K_r * N_r * L_r^%i" % self._ndim_particles,
				"P":"K_r * N_r * L_r^%i" % self._ndim_particles,
				"D":"N_r * L_r^%i" % self._ndim_particles,
				"E":"E_r",
				"B":"B_r",
				"J":"J_r",
				"R":"N_r",
				"Z":"Q_r",
				"N":"",
				}[self._scalarname[0]]
		self._title =self._scalarname
		
		# Finish constructor
		self.valid = True
		return kwargs
	
	# Method to print info on included scalars
	def _info(self):
		return "Scalar "+self._scalarname
	
	# get all available scalars
	def getScalars(self):
		for path in self._results_path:
			try:
				file = path+'/scalars.txt'
				f = open(file, 'r')
			except Exception as e:
				self._error += ["Cannot open 'scalars.txt' in directory '"+path+"'"]
				return []
			try:
				# Find last commented line
				prevline = ""
				for line in f:
					line = line.strip()
					if line[0]!="#": break
					prevline = line[1:].strip()
				scalars = str(prevline).split() # list of scalars
				scalars = scalars[1:] # remove first, which is "time"
			except Exception as e:
				scalars = []
			f.close()
			try:    allScalars = self._np.intersect1d(allScalars, scalars)
			except Exception as e: allScalars = scalars
		return allScalars
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		return self._alltimesteps
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Get value at selected time
		A = self._values[ self._data[t] ]

		if callable(self._data_transform): A = self._data_transform(A)
		return A
