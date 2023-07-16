from .Diagnostic import Diagnostic
from .._Utils import *

class ParticleList(Diagnostic):
	"""Generic class for a diagnostic with listed particles (trackParticles or newParticles)"""
	
	_short_properties_from_raw = {
		"id":"Id", "position/x":"x", "position/y":"y", "position/z":"z",
		"momentum/x":"px", "momentum/y":"py", "momentum/z":"pz",
		"charge":"q", "weight":"w", "chi":"chi",
		"E/x":"Ex", "E/y":"Ey", "E/z":"Ez", "B/x":"Bx", "B/y":"By", "B/z":"Bz",
		"W/x":"Wx", "W/y":"Wy", "W/z":"Wz"
	}
	
	def _initAxes(self, axes):
		
		if type(axes) is not list:
			raise Exception("Error: Argument 'axes' must be a list")
			
		# if axes provided, verify them
		if len(axes)>0:
			self.axes = axes
			unknown_axes = set(axes) - self.available_properties
			if unknown_axes:
				raise Exception(
					"Error: Argument 'axes' has unknown item(s): "+", ".join(unknown_axes)+".\n"
					+ "       Available axes are: "+", ".join(sorted(self.available_properties))
				)
		# otherwise use default
		else:
			self.axes = list(self.available_properties)
		
		# Then figure out axis units
		self._factors = []
		for axis in self.axes:
			axisunits = ""
			if axis == "Id":
				self._centers.append( [0, 281474976710655] )
			elif axis in ["x" , "y" , "z", "moving_x"]:
				axisunits = "L_r"
				self._centers.append( [0., self.namelist.Main.grid_length[{"x":0,"y":1,"z":-1}[axis[-1]]]] )
			elif axis in ["px", "py", "pz"]:
				axisunits = "P_r"
				self._centers.append( [-1., 1.] )
			elif axis == "w":
				axisunits = "N_r * L_r^%i" % self._ndim_particles
				self._centers.append( [0., 1.] )
			elif axis == "q":
				axisunits = "Q_r"
				self._centers.append( [-10., 10.] )
			elif axis == "chi":
				axisunits = "1"
				self._centers.append( [0., 2.] )
			elif axis[0] == "E":
				axisunits = "E_r"
				self._centers.append( [-1., 1.] )
			elif axis[0] == "B":
				axisunits = "B_r"
				self._centers.append( [-1., 1.] )
			elif axis[0] == "W":
				axisunits = "K_r"
				self._centers.append( [0., 1.] )
			self._log += [False]
			self._label += [axis]
			self._units += [axisunits]
			if axis == "Id":
				self._factors += [1]
			else:
				factor, _ = self.units._convert(axisunits, None)
				self._factors += [factor]
		self._title = self._diagType+" '"+self.species+"'"
		self._shape = [0]*len(self.axes)
		self._centers = [self._np.array(c) for c in self._centers]
		self._type = self.axes
		
		# Hack to work with 1 axis
		self._vunits = self._units[0] if len(axes)==1 else ""
	
	def get(self):
		return self.getData()

	# Get a list of files
	def _findFiles(self, filePrefix):
		files = []
		for path in self._results_path:
			file = path+self._os.sep+filePrefix+"_"+self.species+".h5"
			if self._os.path.isfile(file):
				files += [file]
		if not files:
			raise Exception("No files "+filePrefix+" found")
		return files
