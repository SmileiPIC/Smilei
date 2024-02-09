from .Diagnostic import Diagnostic
from .._Utils import *

class Scalar(Diagnostic):
	"""Class for loading a Scalar diagnostic"""
	
	def _init(self, scalar=None, timesteps=None, data_log=False, data_transform=None, **kwargs):
		# Get available scalars
		scalars = self.simulation.getScalars()
		
		# If no scalar chosen, only print the available scalars
		if scalar is None:
			error = []
			if len(scalars)>0:
				error += ["Error: no scalar chosen"]
				error += ["Printing available scalars:"]
				error += ["---------------------------"]
				l = [""]
				for s in scalars:
					if len(s)>4 and s[:2]!=l[-1][:2] and s[-2:]!=l[-1][-2:]:
						if l!=[""]: error += ["\t".join(l)]
						l = []
					l.append(s)
				if l!=[""]: error += ["\t".join(l)]
			else:
				error += ["No scalars found"]
			raise Exception("\n".join(error))
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Find which scalar is requested & associated units
		self.operation = scalar
		def scalarTranslator(s):
			units = "??"
			if   s == "time":
				units = "T_r"
			elif s == "Ubal_norm":
				units = ""
			else:
				units = {
					"U":"K_r * N_r * L_r^%i" % self._ndim_particles,
					"P":"K_r * N_r * L_r^%i" % self._ndim_particles,
					"D":"N_r * L_r^%i" % self._ndim_particles,
					"E":"E_r",
					"B":"B_r",
					"J":"J_r",
					"R":"Q_r * N_r",
					"Z":"Q_r",
					"N":"",
					}[s[0]]
			return units, "S['%s']"%s, s
		self._operation = Operation(scalar, scalarTranslator, self._ureg)
		self._scalarname = self._operation.variables
		self._vunits = self._operation.translated_units
		self._title  = self._operation.title
		if not self._scalarname:
			raise Exception("String "+self.operation+" does not seem to include any scalar")
		
		# Put data_log as object's variable
		self._data_log = data_log
		self._data_transform = data_transform
		
		# Already get the data from the file
		# Loop file line by line
		self._alltimesteps = []
		S = { s:[] for s in self._scalarname }
		for path in self._results_path:
			try:
				f = open(path+'/scalars.txt')
			except:
				continue
			for line in f:
				if line[0]!="#": break
				prevline = line
			scalars = prevline[1:].strip().split() # list of scalars
			scalarindexes = {s:scalars.index(s) for s in self._scalarname} # indexes of the requested scalars
			line = str(line.strip()).split()
			self._alltimesteps += [ int( self._np.round(float(line[0]) / float(self.timestep)) ) ]
			for s,i in scalarindexes.items():
				S[s] += [float(line[i])]
			for line in f:
				line = str(line.strip()).split()
				self._alltimesteps += [ int( self._np.round(float(line[0]) / float(self.timestep)) ) ]
				for s,i in scalarindexes.items():
					S[s] += [float(line[i])]
			f.close()
		self._alltimesteps  = self._np.array(self._alltimesteps)
		for s in S:
			S[s] = self._np.array(S[s])
		self._values = self._operation.eval(locals())
		self._timesteps = self._np.copy(self._alltimesteps)
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "_data" dictionary with the index to each time
		self._data = {}
		for i,t in enumerate(self._timesteps):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		timestep_indices = kwargs.pop("timestep_indices", None)
		self._timesteps = self._selectTimesteps(timesteps, timestep_indices, self._timesteps)
		
		# Need at least one timestep
		assert self._timesteps.size > 0, "Timesteps not found"
		
		# Finish constructor
		self.valid = True
		return kwargs
	
	# Method to print info on included scalars
	def _info(self):
		return "Scalar "+self.operation
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		return self._alltimesteps
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		# Verify that the timestep is valid
		if t not in self._timesteps:
			print("Timestep "+str(t)+" not found in this diagnostic")
			return []
		# Get value at selected time
		A = self._values[ self._data[t] ]
		return A
