from .Diagnostic import Diagnostic
from .._Utils import *

class Scalar(Diagnostic):
	"""Class for loading a scalar diagnostic"""
	
	def _init(self, scalar=None, timesteps=None, data_log=False, **kwargs):
		# Get available scalars
		scalars = self.getScalars()
		
		# If no scalar chosen, only print the available scalars
		if scalar is None:
			if len(scalars)>0:
				self._error += "Printing available scalars:\n"
				self._error += "---------------------------\n"
				l = [""]
				for s in scalars:
					if s[:2] != l[-1][:2] and s[-2:]!=l[-1][-2:]:
						if l!=[""]: self._error += "\t".join(l)+"\n"
						l = []
					l.append(s)
			else:
				self._error += "No scalars found in '"+self._results_path+"'"
			return
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Find which scalar is requested
		if scalar not in scalars:
			fs = list(filter(lambda x:scalar in x, scalars))
			if len(fs)==0:
				self._error += "No scalar `"+scalar+"` found in scalars.txt"
				return
			if len(fs)>1:
				self._error += "Several scalars match: "+(' '.join(fs))+"\n"
				self._error += "Please be more specific and retry.\n"
				return
			scalar = fs[0]
		scalarindex = scalars.index(scalar) # index of the requested scalar
		self._scalarname = scalar
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Already get the data from the file
		# Loop file line by line
		self._times = []
		self._values = []
		with open(self._results_path+'/scalars.txt') as f:
			for line in f:
				line = line.strip()
				if line[0]=="#": continue
				line = str(line).split()
				self._times.append( int( self._np.round(float(line[0]) / float(self.timestep)) ) )
				self._values.append( float(line[scalarindex+1]) )
			self._times  = self._np.array(self._times )
			self._values = self._np.array(self._values)
			self.times = self._times
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "_data" dictionary with the index to each time
		self._data = {}
		for i,t in enumerate(self.times):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				self.times = self._selectTimesteps(timesteps, self.times)
			except:
				self._error += "Argument `timesteps` must be one or two non-negative integers"
				return
		
		# Need at least one timestep
		if self.times.size < 1:
			self._error += "Timesteps not found"
			return
		
		
		# 3 - Build units
		# -------------------------------------------------------------------
		self._vunits = "??"
		if   self._scalarname == "time":
			self._vunits = "T_r"
		elif self._scalarname == "Ubal_norm" or self._scalarname[0] in ["N","Z"]:
			self._vunits = ""
		else:
			self._vunits = {"U":"K_r", "E":"E_r", "B":"B_r", "J":"J_r",
									"R":"N_r", "P":"S_r"}[self._scalarname[0]]
		self._title =self._scalarname
		
		# Finish constructor
		self.valid = True
	
	# Method to print info on included scalars
	def _info(self):
		return "Scalar "+self._scalarname
	
	# get all available scalars
	def getScalars(self):
		try:
			file = self._results_path+'/scalars.txt'
			f = open(file, 'r')
		except:
			print("Cannot open file "+file)
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
		except:
			scalars = []
		f.close()
		return scalars
	
	# get all available timesteps
	def getAvailableTimesteps(self):
		return self._times
	
	# Method to obtain the data only
	def _getDataAtTime(self, t):
		if not self._validate(): return
		# Verify that the timestep is valid
		if t not in self.times:
			print("Timestep "+t+" not found in this diagnostic")
			return []
		# Get value at selected time
		A = self._values[ self._data[t] ]
		# log scale if requested
		if self._data_log: A = self._np.log10(A)
		return A
