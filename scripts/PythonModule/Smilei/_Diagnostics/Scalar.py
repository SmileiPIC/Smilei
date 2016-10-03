from .Diagnostic import Diagnostic
from .._Utils import *

# -------------------------------------------------------------------
# Class for scalars
# -------------------------------------------------------------------
class Scalar(Diagnostic):
	# This is the constructor, which creates the object
	def _init(self, scalar=None, timesteps=None, data_log=False, **kwargs):
		
		if not self.Smilei.valid: return None
		if scalar is None:
			scalars = self.getScalars()
			if len(scalars)>0:
				print("Printing available scalars:")
				print("---------------------------")
				l = [""]
				for s in scalars:
					if s[:2] != l[-1][:2] and s[-2:]!=l[-1][-2:]:
						if l!=[""]: print("\t".join(l))
						l = []
					l.append(s)
			else:
				print("No scalars found in '"+self._results_path+"'")
			return None
		
		# Get available scalars
		scalars = self.getScalars()
		
		# 1 - verifications, initialization
		# -------------------------------------------------------------------
		# Check value of field
		if scalar not in scalars:
			fs = list(filter(lambda x:scalar in x, scalars))
			if len(fs)==0:
				print("No scalar `"+scalar+"` found in scalars.txt")
				return
			if len(fs)>1:
				print("Several scalars match: "+(' '.join(fs)))
				print("Please be more specific and retry.")
				return
			scalar = fs[0]
		self._scalarn = scalars.index(scalar) # index of the requested scalar
		self._scalarname = scalar
		
		# Put data_log as object's variable
		self._data_log = data_log
		
		# Already get the data from the file
		# Loop file line by line
		self.times = []
		self._values = []
		file = self._results_path+'/scalars.txt'
		f = open(file, 'r')
		for line in f:
			line = line.strip()
			if line[0]=="#": continue
			line = str(line).split()
			self.times .append( int( self._np.round(float(line[0]) / float(self.timestep)) ) )
			self._values.append( float(line[self._scalarn+1]) )
		self.times  = self._np.array(self.times )
		self._values = self._np.array(self._values)
		f.close()
		
		# 2 - Manage timesteps
		# -------------------------------------------------------------------
		# fill the "data" dictionary with the index to each time
		self._data = {}
		for i,t in enumerate(self.times):
			self._data.update({ t : i })
		# If timesteps is None, then keep all timesteps otherwise, select timesteps
		if timesteps is not None:
			try:
				ts = self._np.array(self._np.double(timesteps),ndmin=1)
				if ts.size==2:
					# get all times in between bounds
					self.times = self.times[ (self.times>=ts[0]) * (self.times<=ts[1]) ]
				elif ts.size==1:
					# get nearest time
					self.times = self._np.array([self.times[(self._np.abs(self.times-ts)).argmin()]])
				else:
					raise
			except:
				print("Argument `timesteps` must be one or two non-negative integers")
				return
		
		# Need at least one timestep
		if self.times.size < 1:
			print("Timesteps not found")
			return
		
		
		# 3 - Manage axes
		# -------------------------------------------------------------------
		# There are no axes for scalars
		self._naxes = 0
		self._slices = []
		# Build units
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
	def info(self):
		if not self._validate(): return
		print("Scalar "+self._scalarname)
		return
	
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
		return self.times
	
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
