
def setMatplotLibBackend(show=True):
	import matplotlib, sys
	usingAgg = (matplotlib.get_backend().lower() == "agg")
	if not show and not usingAgg:
		if "matplotlib.pyplot" in sys.modules:
			print("WARNING: 'show=False' requires you restart python.")
		else:
			matplotlib.use("Agg")
	if show and usingAgg:
		if "matplotlib.pyplot" in sys.modules:
			print("WARNING: 'show=False' was set earlier. Restart python if you want figures to appear.")
	matplotlib.rcParams['font.family'] = 'serif'
	matplotlib.rcParams['font.serif'] = 'Times New Roman'
	#print(matplotlib.get_backend())


def updateMatplotLibColormaps():
	# Add smilei's colormap to the available matplotlib colormaps
	import matplotlib.colors, matplotlib.pyplot
	if "smilei" in matplotlib.pyplot.colormaps(): return
	def register(name, d):
		cmap = matplotlib.colors.LinearSegmentedColormap(name, d, N=256, gamma=1.0)
		matplotlib.pyplot.register_cmap(cmap=cmap)
	register(u"smilei", {
			'red'  :((0., 0., 0.), (0.0625 , 0.091, 0.091), (0.09375, 0.118, 0.118), (0.125 , 0.127, 0.127), (0.1875 , 0.135, 0.135), (0.21875, 0.125, 0.125), (0.28125, 0.034, 0.034), (0.3125 , 0.010, 0.010), (0.34375, 0.009, 0.009), (0.4375 , 0.049, 0.049), (0.46875, 0.057, 0.057), (0.5 , 0.058, 0.058), (0.59375, 0.031, 0.031), (0.625 , 0.028, 0.028), (0.65625, 0.047, 0.047), (0.71875, 0.143, 0.143), (0.78125, 0.294, 0.294), (0.84375, 0.519, 0.519), (0.90625, 0.664, 0.664), (0.9375 , 0.760, 0.760), (0.96875, 0.880, 0.880), (1., 1., 1. )),
			'green':((0., 0., 0.), (0.21875, 0.228, 0.228), (0.78125, 0.827, 0.827), (0.8125 , 0.852, 0.852), (0.84375, 0.869, 0.869), (0.9375 , 0.937, 0.937), (0.96875, 0.967, 0.967), (1. , 1. , 1. )),
			'blue' :((0., 0., 0.), (0.0625 , 0.174, 0.184), (0.09375, 0.207, 0.207), (0.21875, 0.447, 0.447), (0.25 , 0.500, 0.500), (0.5 , 0.507, 0.507), (0.5625 , 0.502, 0.502), (0.625 , 0.485, 0.485), (0.6875 , 0.452, 0.452), (0.75 , 0.398, 0.398), (0.78125, 0.363, 0.363), (0.8125 , 0.345, 0.345), (0.84375, 0.377, 0.377), (0.90625, 0.534, 0.534), (0.9375 , 0.660, 0.660), (0.96875, 0.790, 0.790), (1. , 1. , 1. ))
		})
	register(u"smilei_r", {
			'red'  :((0.0, 1.0, 1.0), (0.03125, 0.88, 0.88), (0.0625, 0.76, 0.76), (0.09375, 0.664, 0.664), (0.15625, 0.519, 0.519), (0.21875, 0.294, 0.294), (0.28125, 0.143, 0.143), (0.34375, 0.047, 0.047), (0.375, 0.028, 0.028), (0.40625, 0.031, 0.031), (0.5, 0.058, 0.058), (0.53125, 0.057, 0.057), (0.5625, 0.049, 0.049), (0.65625, 0.009, 0.009), (0.6875, 0.01, 0.01), (0.71875, 0.034, 0.034), (0.78125, 0.125, 0.125), (0.8125, 0.135, 0.135), (0.875, 0.127, 0.127), (0.90625, 0.118, 0.118), (0.9375, 0.091, 0.091), (1.0, 0.0, 0.0)),
			'green':((0.0, 1.0, 1.0), (0.03125, 0.967, 0.967), (0.0625, 0.937, 0.937), (0.15625, 0.869, 0.869), (0.1875, 0.852, 0.852), (0.21875, 0.827, 0.827), (0.78125, 0.228, 0.228), (1.0, 0.0, 0.0)),
			'blue' :((0.0, 1.0, 1.0), (0.03125, 0.79, 0.79), (0.0625, 0.66, 0.66), (0.09375, 0.534, 0.534), (0.15625, 0.377, 0.377), (0.1875, 0.345, 0.345), (0.21875, 0.363, 0.363), (0.25, 0.398, 0.398), (0.3125, 0.452, 0.452), (0.375, 0.485, 0.485), (0.4375, 0.502, 0.502), (0.5, 0.507, 0.507), (0.75, 0.5, 0.5), (0.78125, 0.447, 0.447), (0.90625, 0.207, 0.207), (0.9375, 0.174, 0.184), (1.0, 0.0, 0.0))
		})
	register(u"smileiD", {
			'red'  :((0.000,0.786,0.786),(0.033,0.807,0.807),(0.067,0.868,0.868),(0.100,0.932,0.932),(0.133,0.927,0.927),(0.167,0.919,0.919),(0.200,0.913,0.913),(0.233,0.909,0.909),(0.267,0.908,0.908),(0.300,0.909,0.909),(0.333,0.913,0.913),(0.367,0.920,0.920),(0.400,0.930,0.930),(0.433,0.943,0.943),(0.467,0.961,0.961),(0.500,1.000,1.000),(0.533,0.972,0.972),(0.567,0.914,0.914),(0.600,0.853,0.853),(0.633,0.789,0.789),(0.667,0.721,0.721),(0.700,0.648,0.648),(0.733,0.570,0.570),(0.767,0.483,0.483),(0.800,0.383,0.383),(0.833,0.261,0.261),(0.867,0.116,0.116),(0.900,0.215,0.215),(0.933,0.280,0.280),(0.967,0.288,0.288),(1.000,0.237,0.237),),
			'green':((0.000,0.025,0.025),(0.033,0.167,0.167),(0.067,0.144,0.144),(0.100,0.004,0.004),(0.133,0.248,0.248),(0.167,0.364,0.364),(0.200,0.452,0.452),(0.233,0.527,0.527),(0.267,0.594,0.594),(0.300,0.656,0.656),(0.333,0.714,0.714),(0.367,0.770,0.770),(0.400,0.824,0.824),(0.433,0.877,0.877),(0.467,0.928,0.928),(0.500,1.000,1.000),(0.533,0.987,0.987),(0.567,0.965,0.965),(0.600,0.945,0.945),(0.633,0.926,0.926),(0.667,0.909,0.909),(0.700,0.893,0.893),(0.733,0.877,0.877),(0.767,0.862,0.862),(0.800,0.846,0.846),(0.833,0.830,0.830),(0.867,0.810,0.810),(0.900,0.774,0.774),(0.933,0.733,0.733),(0.967,0.697,0.697),(1.000,0.671,0.671),),
			'blue' :((0.000,0.681,0.681),(0.033,0.701,0.701),(0.067,0.736,0.736),(0.100,0.834,0.834),(0.133,0.854,0.854),(0.167,0.861,0.861),(0.200,0.867,0.867),(0.233,0.872,0.872),(0.267,0.877,0.877),(0.300,0.885,0.885),(0.333,0.894,0.894),(0.367,0.906,0.906),(0.400,0.921,0.921),(0.433,0.939,0.939),(0.467,0.960,0.960),(0.500,1.000,1.000),(0.533,0.971,0.971),(0.567,0.914,0.914),(0.600,0.853,0.853),(0.633,0.789,0.789),(0.667,0.723,0.723),(0.700,0.654,0.654),(0.733,0.581,0.581),(0.767,0.503,0.503),(0.800,0.419,0.419),(0.833,0.323,0.323),(0.867,0.199,0.199),(0.900,0.119,0.119),(0.933,0.206,0.206),(0.967,0.272,0.272),(1.000,0.284,0.284),),
		})
	register(u"smileiD_r", {
			'red'  :((0.000,0.237,0.237),(0.033,0.288,0.288),(0.067,0.280,0.280),(0.100,0.215,0.215),(0.133,0.116,0.116),(0.167,0.261,0.261),(0.200,0.383,0.383),(0.233,0.483,0.483),(0.267,0.570,0.570),(0.300,0.648,0.648),(0.333,0.721,0.721),(0.367,0.789,0.789),(0.400,0.853,0.853),(0.433,0.914,0.914),(0.467,0.972,0.972),(0.500,1.000,1.000),(0.533,0.961,0.961),(0.567,0.943,0.943),(0.600,0.930,0.930),(0.633,0.920,0.920),(0.667,0.913,0.913),(0.700,0.909,0.909),(0.733,0.908,0.908),(0.767,0.909,0.909),(0.800,0.913,0.913),(0.833,0.919,0.919),(0.867,0.927,0.927),(0.900,0.932,0.932),(0.933,0.868,0.868),(0.967,0.807,0.807),(1.000,0.786,0.786),),
			'green':((0.000,0.671,0.671),(0.033,0.697,0.697),(0.067,0.733,0.733),(0.100,0.774,0.774),(0.133,0.810,0.810),(0.167,0.830,0.830),(0.200,0.846,0.846),(0.233,0.862,0.862),(0.267,0.877,0.877),(0.300,0.893,0.893),(0.333,0.909,0.909),(0.367,0.926,0.926),(0.400,0.945,0.945),(0.433,0.965,0.965),(0.467,0.987,0.987),(0.500,1.000,1.000),(0.533,0.928,0.928),(0.567,0.877,0.877),(0.600,0.824,0.824),(0.633,0.770,0.770),(0.667,0.714,0.714),(0.700,0.656,0.656),(0.733,0.594,0.594),(0.767,0.527,0.527),(0.800,0.452,0.452),(0.833,0.364,0.364),(0.867,0.248,0.248),(0.900,0.004,0.004),(0.933,0.144,0.144),(0.967,0.167,0.167),(1.000,0.025,0.025),),
			'blue' :((0.000,0.284,0.284),(0.033,0.272,0.272),(0.067,0.206,0.206),(0.100,0.119,0.119),(0.133,0.199,0.199),(0.167,0.323,0.323),(0.200,0.419,0.419),(0.233,0.503,0.503),(0.267,0.581,0.581),(0.300,0.654,0.654),(0.333,0.723,0.723),(0.367,0.789,0.789),(0.400,0.853,0.853),(0.433,0.914,0.914),(0.467,0.971,0.971),(0.500,1.000,1.000),(0.533,0.960,0.960),(0.567,0.939,0.939),(0.600,0.921,0.921),(0.633,0.906,0.906),(0.667,0.894,0.894),(0.700,0.885,0.885),(0.733,0.877,0.877),(0.767,0.872,0.872),(0.800,0.867,0.867),(0.833,0.861,0.861),(0.867,0.854,0.854),(0.900,0.834,0.834),(0.933,0.736,0.736),(0.967,0.701,0.701),(1.000,0.681,0.681),),
		})

class ChunkedRange:
	def __init__(self, size, chunksize):
		self.size = int(size)
		self.nchunks = int( (size-1) // int(chunksize) + 1 )
		self.adjustedchunksize = int( (size-1) // self.nchunks + 1 )
		self.ichunk = 0
	def __iter__(self):
		return self
	def __next__(self):
		if self.ichunk < self.nchunks:
			first = self.ichunk * self.adjustedchunksize
			last  = min( first + self.adjustedchunksize, self.size )
			self.ichunk += 1
			return (first, last, last-first)
		else:
			raise StopIteration
	def next(self): # for python 2
		return self.__next__()


def openNamelist(namelist):
	"""
	Function to execute a namelist and store all its content in the returned object.

	Example:
		namelist = happi.openNamelist("path/no/my/namelist.py")
		print( namelist.Main.timestep)
	"""

	from . import happi_directory
	from os import sep
	from os.path import isdir, exists
	smilei_python_directory = happi_directory + sep + ".." + sep + "src" + sep + "Python"
	if not isdir(smilei_python_directory):
		raise Exception("Cannot find the Smilei/src directory")
	for file in ["pyinit.py", "pycontrol.py", "pyprofiles.py"]:
		if not exists(smilei_python_directory + sep + file):
			raise Exception("Cannot find the Smilei/src/Python/"+file+" file")
	namespace={}
	exec(open(smilei_python_directory+sep+"pyinit.py").read(), namespace)
	exec(open(smilei_python_directory+sep+"pyprofiles.py").read(), namespace)
	exec(open(namelist).read(), namespace) # execute the namelist
	exec(open(smilei_python_directory+sep+"pycontrol.py").read(), namespace)
	class Namelist: pass # empty class to store the namelist variables
	namelist = Namelist() # create new empty object
	for key, value in namespace.items(): # transfer all variables to this object
		if key[0]=="_": continue # skip builtins
		setattr(namelist, key, value)
	return namelist


class Options(object):
	""" Class to contain matplotlib plotting options """

	def __init__(self, **kwargs):
		self.figure  = 1
		self.xfactor = None
		self.xmin    = None
		self.xmax    = None
		self.yfactor = None
		self.ymin    = None
		self.ymax    = None
		self.vfactor = None
		self.vmin    = None
		self.vmax    = None
		self.figure0 = {}
		self.figure1 = {"facecolor":"w"}
		self.axes = {}
		self.plot = {}
		self.image = {"cmap":"smilei", "interpolation":"nearest", "aspect":"auto"}
		self.colorbar = {}
		self.xtick = {"useOffset":False}
		self.ytick = {"useOffset":False}
		self.side = "left"
		self.transparent = None
		self.export_dir = None

	# Method to set optional plotting arguments
	def set(self, **kwargs):
		# First, we manage the main optional arguments
		self.figure0["num"] = kwargs.pop("figure", self.figure)
		self.xfactor     = kwargs.pop("xfactor"    , self.xfactor  )
		self.xmin        = kwargs.pop("xmin"       , self.xmin  )
		self.xmax        = kwargs.pop("xmax"       , self.xmax  )
		self.yfactor     = kwargs.pop("yfactor"    , self.yfactor  )
		self.ymin        = kwargs.pop("ymin"       , self.ymin  )
		self.ymax        = kwargs.pop("ymax"       , self.ymax  )
		self.vfactor     = kwargs.pop("vfactor"    , self.vfactor  )
		self.vmin        = kwargs.pop("vmin"       , self.vmin )
		self.vmax        = kwargs.pop("vmax"       , self.vmax )
		self.side        = kwargs.pop("side"       , self.side )
		self.transparent = kwargs.pop("transparent", self.transparent )
		self.export_dir  = kwargs.pop("export_dir", self.export_dir )
		# Second, we manage all the other arguments that are directly the ones of matplotlib
		for kwa, val in kwargs.copy().items():
			if kwa in ["figsize"]:
				self.figure0[kwa] = val
			elif kwa in ["facecolor","edgecolor"]:
				self.figure1[kwa] = val
			elif kwa in ["aspect","axis_bgcolor",
					   "frame_on","position","title","visible","xlabel","xscale","xticklabels",
					   "xticks","ylabel","yscale","yticklabels","yticks","zorder"]:
				self.axes[kwa] = val
			elif kwa in ["color","dashes","drawstyle","fillstyle","label","linestyle",
					   "linewidth","marker","markeredgecolor","markeredgewidth",
					   "markerfacecolor","markerfacecoloralt","markersize","markevery",
					   "visible","zorder"]:
				self.plot[kwa] = val
			elif kwa in ["cmap","aspect","interpolation"]:
				self.image[kwa] = val
			elif kwa in ["orientation","fraction","pad","shrink","anchor","panchor",
					   "extend","extendfrac","extendrect","spacing","ticks","format",
					   "drawedges"]:
				self.colorbar[kwa] = val
			elif kwa in ["style_x","scilimits_x","useOffset_x"]:
				self.xtick[kwa] = val
			elif kwa in ["style_y","scilimits_y","useOffset_y"]:
				self.ytick[kwa] = val
			else:
				continue
			kwargs.pop(kwa)
		# special case: "aspect" is ambiguous because it exists for both imshow and colorbar
		if "cbaspect" in kwargs:
			self.colorbar["aspect"] = kwargs.pop("cbaspect")
		if self.side=="right" and "pad" not in self.colorbar:
			self.colorbar["pad"] = 0.15
		return kwargs

PintWarningIssued = False

class Units(object):
	""" Units()

	Class to handle units smartly. Based on the *pint* package.
	"""

	def __init__(self, *args, **kwargs):
		# All args are parsed
		self.requestedUnits = []
		self.requestedX = ""
		self.requestedY = ""
		self.requestedV = ""
		self.verbose = True
		for a in args:
			if type(a) is str:
				self.requestedUnits.append( a )
			else:
				raise TypeError("Arguments of Units() should be strings")
		for kwa, val in kwargs.items():
			if kwa == "verbose": self.verbose = val
			else:
				if type(val) is not str:
					raise TypeError("Arguments of Units() should be strings")
				if   kwa == "x": self.requestedX = val
				elif kwa == "y": self.requestedY = val
				elif kwa == "v": self.requestedV = val
				else: raise TypeError("Units() got an unexpected keyword argument '"+kwa+"'")

		# We try to import the pint package
		self.UnitRegistry = None
		try:
			from pint import UnitRegistry
			self.UnitRegistry = UnitRegistry
		except:
			global PintWarningIssued
			if self.verbose and not PintWarningIssued:
				print("WARNING: you do not have the *pint* package, so you cannot modify units.")
				print("       : The results will stay in code units.")
				PintWarningIssued = True
			return
	
	def _getUnits(self, units):
		if self.UnitRegistry:
			u = self.ureg(units)
			try: u = u.units.format_babel()
			except: u = ""
			return u
		else:
			return "1"	
	
	def _divide(self,units1, units2):
		division = self.ureg("("+units1+") / ("+units2+")").to_base_units()
		if not division.dimensionless: raise
		return division.magnitude or 1., units2

	def _convert(self, knownUnits, requestedUnits):
		if knownUnits:
			if requestedUnits:
				try:
					return self._divide(knownUnits,requestedUnits)
				except:
					if self.verbose:
						print("WARNING: cannot convert units to <"+requestedUnits+">")
						print("       : Conversion discarded.")
			else:
				for units in self.requestedUnits:
					try   : return self._divide(knownUnits,units)
					except: pass
			try:
				if knownUnits=="1": knownUnits=""
				val = self.ureg(knownUnits)
				return 1., u"{0.units:P}".format(val)
			except:
				if self.verbose:
					print("WARNING: units unknown: "+str(knownUnits))
				return 1., ""
		return 1., ""

	def prepare(self, reference_angular_frequency_SI=None):
		if self.UnitRegistry:
			if reference_angular_frequency_SI:
				# Load pint's default unit registry
				self.ureg = self.UnitRegistry()
				# Define code units
				self.ureg.define("V_r = speed_of_light"                   ) # velocity
				self.ureg.define("W_r = "+str(reference_angular_frequency_SI)+"*hertz") # frequency
				self.ureg.define("M_r = electron_mass"                    ) # mass
				self.ureg.define("Q_r = 1.602176565e-19 * coulomb"        ) # charge
			else:
				# Make blank unit registry
				self.ureg = self.UnitRegistry(None)
				self.ureg.define("V_r = [code_velocity]"                  ) # velocity
				self.ureg.define("W_r = [code_frequency]"                 ) # frequency
				self.ureg.define("M_r = [code_mass]"                      ) # mass
				self.ureg.define("Q_r = [code_charge]"                    ) # charge
				self.ureg.define("epsilon_0 = 1")
				# Add radians and degrees
				self.ureg.define("radian    = [] = rad"               )
				self.ureg.define("degree    = pi/180*radian = deg"    )
				self.ureg.define("steradian = radian ** 2 = sr"       )
			self.ureg.define("L_r = V_r / W_r"                        ) # length
			self.ureg.define("T_r = 1   / W_r"                        ) # time
			self.ureg.define("P_r = M_r * V_r"                        ) # momentum
			self.ureg.define("K_r = M_r * V_r**2"                     ) # energy
			self.ureg.define("N_r = epsilon_0 * M_r * W_r**2 / Q_r**2") # density
			self.ureg.define("J_r = V_r * Q_r * N_r"                  ) # current
			self.ureg.define("B_r = M_r * W_r / Q_r"                  ) # magnetic field
			self.ureg.define("E_r = B_r * V_r"                        ) # electric field
			self.ureg.define("S_r = K_r * V_r * N_r"                  ) # poynting
	
	def convertAxes(self, xunits="", yunits="", vunits="", tunits=""):
		if self.UnitRegistry:
			self.xcoeff, self.xname = self._convert(xunits, self.requestedX)
			self.ycoeff, self.yname = self._convert(yunits, self.requestedY)
			self.vcoeff, self.vname = self._convert(vunits, self.requestedV)
			self.tcoeff, self.tname = self._convert("T_r", None)
		else:
			self.xcoeff, self.xname = 1., "code units"
			self.ycoeff, self.yname = 1., "code units"
			self.vcoeff, self.vname = 1., "code units"
			self.tcoeff, self.tname = 1., "code units"





class Movie:

	def __init__(self, fig, movie="", fps=15, dpi=200):
		import os.path as ospath
		self.writer = None
		if type(movie) is not str:
			print("ERROR: argument 'movie' must be a filename")
			return
		if len(movie)>0:
			if ospath.isfile(movie):
				print("ERROR: file '"+movie+"' already exists. Please choose another name")
				return
			if ospath.isdir(movie):
				print("ERROR: '"+movie+"' is a path, not a file")
				return
			try:
				import matplotlib.animation as anim
			except:
				print("ERROR: it looks like your version of matplotlib is too old for movies")
				return
			filename, file_extension = ospath.splitext(movie)
			if file_extension == ".gif":
				writer = 'imagemagick_file'
			else:
				writer = 'ffmpeg'
			try:
				self.writer = anim.writers[writer](fps=fps)
			except:
				print("ERROR: you need the '"+writer+"' software to make movies")
				return

			self.writer.setup(fig, movie, dpi)

	def finish(self):
		if self.writer is not None:
			self.writer.finish()

	def grab_frame(self):
		if self.writer is not None:
			self.writer.grab_frame()



class SaveAs:

	def __init__(self, smartPath, fig, plt):
		import os.path as p
		from os import sep as sep
		default_extension = ".png"
		self.figure = fig

		self.prefix = False
		if type(smartPath) is str:
			path = smartPath
			if len(smartPath)==0 or smartPath[0].isalnum():
				path = "."+sep+path
			if p.isdir(path):
				self.prefix = p.normpath(path)+sep
				self.suffix = default_extension
			else:
				path, base = p.split(path)
				if p.isdir(path):
					basesplit = base.rsplit(".",1)
					self.prefix = path+sep+basesplit[0]
					self.suffix = ""
					if len(basesplit)>1: self.suffix = "."+basesplit[1]
				else:
					self.prefix = False
					self.suffix = "`"+path+"` is not a directory"
		else:
			return

		if not self.prefix:
			print("WARNING: "+self.suffix)
			print("         Will not save figures to files")

		# For some reason, the following freezes the animation; removed temporarily
		#else:
		#	supportedTypes=plt.matplotlib.backend_bases.FigureCanvasBase(fig).get_supported_filetypes()
		#	if self.suffix.strip(".") not in supportedTypes.keys():
		#		print("WARNING: file format not supported, will not save figures to files")
		#		print("         Supported formats: "+",".join(supportedTypes.keys()))
		#		self.prefix = False

	def frame(self, id=None):
		if self.prefix:
			file = self.prefix + ("" if id is None else "%010d"%int(id)) + self.suffix
			self.figure.savefig(file)



def multiPlot(*Diags, **kwargs):
	""" multiplot(Diag1, Diag2, ...,
	              shape=None,
	              movie="", fps=15, dpi=200, saveAs=None,
	              skipAnimation=False
	              )

	Plots simultaneously several diagnostics.

	Parameters:
	-----------
	Diag1, Diag2, ... : Several objects of classes 'Scalar', 'Field', 'Probe' or 'ParticleBinning'
	shape : 2-element list giving the number of figures in x and y.
	movie : filename to create a movie, e.g. "my/path/mov.avi" or "my/path/mov.gif"
	fps : frames per second for the movie.
	dpi : resolution of the movie.
	saveAs : path where to store individual frames as pictures, e.g. "my/path/fig.png"
	skipAnimation : if True, plots only the last frame.
	"""

	from ._Diagnostics import TrackParticles

	# Verify Diags are valid
	nDiags = len(Diags)
	if nDiags == 0: return
	for Diag in Diags:
		if not Diag.valid: return
	np  = Diags[0]._np  # numpy
	plt = Diags[0]._plt # pyplot
	# Get keyword arguments
	shape  = kwargs.pop("shape" , None)
	movie  = kwargs.pop("movie" , ""  )
	fps    = kwargs.pop("fps"   , 15  )
	dpi    = kwargs.pop("dpi"   , 200 )
	saveAs = kwargs.pop("saveAs", None)
	skipAnimation = kwargs.pop("skipAnimation", False )
	timesteps = kwargs.pop("timesteps", None )
	# Gather all times
	alltimes = []
	for Diag in Diags:
		diagtimes = Diag.getTimesteps()
		if timesteps is not None:
			diagtimes = Diag._selectTimesteps(timesteps, diagtimes)
		diagtimes = list( diagtimes*Diag.timestep )
		if skipAnimation: alltimes += [diagtimes[-1]]
		else            : alltimes += diagtimes
	alltimes = np.sort(np.unique(alltimes))
	# Determine whether to plot all cases on the same axes
	sameAxes = False
	if shape is None or shape == [1,1]:
		sameAxes = True
		for d in Diags:
			if type(d) is TrackParticles or d._type!=Diags[0]._type or d.dim!=Diags[0].dim:
				sameAxes = False
				break
	if not sameAxes and shape == [1,1] and nDiags>1:
		print("Cannot have shape=[1,1] with these diagnostics")
		return
	# Determine the shape
	if sameAxes: shape = [1,1]
	if shape is None: shape = [nDiags,1]
	nplots = np.array(shape).prod()
	if not sameAxes and nplots != nDiags:
		print("The 'shape' argument is incompatible with the number of diagnostics:")
		print("  "+str(nDiags)+" diagnostics do not fit "+str(nplots)+" plots")
		return
	# Make the figure
	if "facecolor" not in kwargs: kwargs.update({ "facecolor":"w" })
	options = Options()
	options.set(**kwargs)
	fig = plt.figure(**options.figure0)
	fig.set(**options.figure1) # Apply figure kwargs
	fig.clf()
	fig.subplots_adjust(wspace=0.5, hspace=0.5, bottom=0.15)
	ax = []
	xmin =  float("inf")
	xmax = -float("inf")
	ymin =  float("inf")
	ymax = -float("inf")
	option_xmin = []
	option_xmax = []
	option_ymin = []
	option_ymax = []
	try:
		c = plt.matplotlib.rcParams['axes.color_cycle']
	except:
		c = plt.matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
	for i in range(nplots):
		ax.append( fig.add_subplot(shape[0], shape[1], i+1) )
	rightside = [d.options.side=="right" for d in Diags]
	allright  = all(rightside)
	bothsides = any(rightside) and not allright
	for i, Diag in enumerate(Diags):
		Diag._cax_id = 0
		if sameAxes:
			Diag._ax = ax[0]
			if Diag.dim==2: Diag._cax_id = i
		else:
			Diag._ax = ax[i]
		if Diag.options.side == "right":
			if sameAxes and not allright:
				try   : Diag._ax.twin # check if twin exists
				except: Diag._ax.twin = Diag._ax.twinx()
				Diag._ax = Diag._ax.twin
			else:
				Diag._ax.yaxis.tick_right()
				Diag._ax.yaxis.set_label_position("right")
		Diag._plot = None
		if Diag.options.xmin is not None: option_xmin += [Diag.options.xmin]
		if Diag.options.xmax is not None: option_xmax += [Diag.options.xmax]
		if Diag.options.ymin is not None: option_ymin += [Diag.options.ymin]
		if Diag.options.ymax is not None: option_ymax += [Diag.options.ymax]
		if "color" not in Diag.options.plot:
			Diag.options.plot.update({ "color":c[i%len(c)] })
		Diag._prepare()
		l = Diag.limits()
		if len(l) > 0:
			if Diag.options.xmin is None: xmin = min(xmin,l[0][0])
			if Diag.options.xmax is None: xmax = max(xmax,l[0][1])
			if len(l) > 1:
				if Diag.options.ymin is None: ymin = min(ymin,l[1][0])
				if Diag.options.ymax is None: ymax = max(ymax,l[1][1])
	# Find min max
	if option_xmin: xmin = min([xmin]+option_xmin)
	if option_xmax: xmax = max([xmax]+option_xmax)
	if option_ymin: ymin = min([ymin]+option_ymin)
	if option_ymax: ymax = max([ymax]+option_ymax)
	
	# Static plot
	if sameAxes and Diags[0].dim==0:
		for Diag in Diags:
			Diag._plotOnAxes(Diag._ax, Diag.getTimesteps()[-1])
			plt.draw()
			plt.pause(0.00001)
	# Animated plot
	else:
		# Loop all times
		mov = Movie(fig, movie, fps, dpi)
		save = SaveAs(saveAs, fig, plt)
		for i,time in enumerate(alltimes):
			t = None
			for Diag in Diags:
				t = np.round(time/Diag.timestep) # convert time to timestep
				if t in Diag.getTimesteps():
					if Diag._plot is None:
						Diag._plotOnAxes(Diag._ax, t, cax_id = Diag._cax_id)
					else:
						Diag._animateOnAxes(Diag._ax, t, cax_id = Diag._cax_id)
					if sameAxes:
						Diag._ax.set_xlim(xmin,xmax)
						if Diag.dim<2 and bothsides:
							color = Diag._plot.get_color()
							Diag._ax.yaxis.label.set_color(color)
							Diag._ax.tick_params(axis='y', colors=color)
							if Diag.options.side == "right":
								Diag._ax.spines['right'].set_color(color)
								Diag._ax.spines['left'].set_color((1.,1.,1.,0.))
							else:
								Diag._ax.spines['left'].set_color(color)
					try: Diag._ax.set_position(Diag._ax.twin.get_position())
					except: pass
			plt.legend()
			plt.draw()
			plt.pause(0.00001)
			mov.grab_frame()
			if t is not None: save.frame(int(t))
		mov.finish()
		return


class VTKfile:

	def __init__(self):
		try:
			import vtk
		except:
			print("Python module 'vtk' not found. Could not export to VTK format")
			return

		self.vtk = vtk

	def Array(self, data, name):
		"""
		Convert a numpy array in a vtkFloatArray
		"""

		from numpy import float32
		from numpy import int32

		shape = data.shape
		if len(shape)==1:   npoints, nComponents = shape[0], 1
		elif len(shape)==2: npoints, nComponents = shape
		else: raise Exception("impossible")

		if (data.dtype == int32):

			arr = self.vtk.vtkIntArray()
			arr.SetNumberOfTuples(npoints)
			arr.SetNumberOfComponents(nComponents)
			# Replace the pointer in arr by the pointer to the data
			arr.SetVoidArray(data, npoints*nComponents, 1)
			arr.SetName(name)
			# keep reference to "data" to not have a null reference
			# else pcoords would be deleted
			# (see the trick: http://vtk.1045678.n5.nabble.com/More-zero-copy-array-support-for-Python-td5743662.html)
			arr.array = data

			return arr

		elif (data.dtype == float32):

			arr = self.vtk.vtkFloatArray()
			arr.SetNumberOfTuples(npoints)
			arr.SetNumberOfComponents(nComponents)
			# Replace the pointer in arr by the pointer to the data
			arr.SetVoidArray(data, npoints*nComponents, 1)
			arr.SetName(name)
			# keep reference to "data" to not have a null reference
			# else pcoords would be deleted
			# (see the trick: http://vtk.1045678.n5.nabble.com/More-zero-copy-array-support-for-Python-td5743662.html)
			arr.array = data

			return arr

		else:
			raise Exception("In Array: Unknown data type for data {data.dtype}".format())


	def WriteImage(self, array, origin, extent, spacings, file, numberOfPieces):
		img = self.vtk.vtkImageData()
		img.SetOrigin(origin)
		img.SetExtent(extent)
		img.SetSpacing(spacings)
		img.GetPointData().SetScalars(array)
		writer = self.vtk.vtkXMLPImageDataWriter()
		writer.SetFileName(file)
		writer.SetNumberOfPieces(numberOfPieces)
		writer.SetEndPiece(numberOfPieces-1)
		writer.SetStartPiece(0);
		if float(self.vtk.VTK_VERSION[:3]) < 6:
		    writer.SetInput(img)
		else:
		    writer.SetInputData(img)
		writer.Write()

	def WriteRectilinearGrid(self, dimensions, xcoords, ycoords, zcoords, array, file):
		grid = self.vtk.vtkRectilinearGrid()
		grid.SetDimensions(dimensions)
		grid.SetXCoordinates(xcoords)
		grid.SetYCoordinates(ycoords)
		grid.SetZCoordinates(zcoords)
		grid.GetPointData().SetScalars(array)
		writer = self.vtk.vtkRectilinearGridWriter()
		writer.SetFileName(file)
		if float(self.vtk.VTK_VERSION[:3]) < 6:
		    writer.SetInput(grid)
		else:
		    writer.SetInputData(grid)
		writer.Write()

	def WriteCloud(self, pcoords, attributes, data_format, file):
		"""
		Create a vtk file that describes a cloud of points (using vtkPolyData)

		Inputs:

			pcoodrs: vtk array that describes the point coordinates
			attributes: Vtk arrays containing additional values for each point
			data_format: the output data format
			file: output file path

		"""

		points = self.vtk.vtkPoints()
		#points.SetDataTypeToFloat()
		points.SetData(pcoords)

		pdata = self.vtk.vtkPolyData()
		pdata.SetPoints(points)

		# Add scalars for xml
		if data_format == "xml":

			for attribute in attributes:
				# SetScalar allows only one scalar
				# pdata.GetPointData().SetScalars(attribute)

				# AddArray creates scalar and then fields
				pdata.GetPointData().AddArray(attribute)

			# The first attribute (first scalar) is the main one
			if len(attributes) > 0:
				pdata.GetPointData().SetActiveScalars(attributes[0].GetName())

		# Add scalars for vtk
		else:

			if len(attributes) > 0:
				pdata.GetPointData().SetScalars(attributes[0])
				pdata.GetPointData().SetActiveScalars(attributes[0].GetName())



		writer = self.vtk.vtkPolyDataWriter()

		# For XML output
		if data_format == "xml":
			writer = self.vtk.vtkXMLDataSetWriter()

		writer.SetFileName(file)
		if float(self.vtk.VTK_VERSION[:3]) < 6:
		    writer.SetInput(pdata)
		else:
		    writer.SetInputData(pdata)
		writer.Write()

		# Add the following attributes by hand because the API
		# limits vtk to 1 scalar
		if data_format == "vtk":
			file_object = open(file, 'a')
			for attribute in attributes[1:]:
				if (attribute.GetDataType() == 6):
					data_type = "int"
				elif (attribute.GetDataType() == 10):
					data_type = "float"
				size = attribute.GetSize()
				file_object.write("POINT_DATA {} \n".format(pdata.GetNumberOfPoints()))
				file_object.write("SCALARS {} {} \n".format(attribute.GetName(),data_type))
				file_object.write("LOOKUP_TABLE default \n")
				for i in range(0,size,8):
					remaining = min(size - i,8)
					for j in range(remaining):
						file_object.write("{} ".format(attribute.GetValue(i + j)))
					file_object.write("\n")


	def WritePoints(self, pcoords, file):
		points = self.vtk.vtkPoints()
		points.SetData(pcoords)
		grid = self.vtk.vtkUnstructuredGrid()
		grid.SetPoints(points)
		writer = self.vtk.vtkUnstructuredGridWriter()
		writer.SetFileName(file)
		if float(self.vtk.VTK_VERSION[:3]) < 6:
		    writer.SetInput(grid)
		else:
		    writer.SetInputData(grid)
		writer.Write()

	def WriteLines(self, pcoords, connectivity, attributes, data_format, file):
		"""
		Create a vtk file that describes lines such as trajectories

		Inputs:

			pcoodrs: vtk array that describes the point coordinates
			connectivity: connection betwwen coordiantes in pcoords to form trajectories
			attributes: Vtk arrays containing additional values for each point
			data_format: the output data format
			file: output file path
		"""
		ncel = len(connectivity)
		connectivity = connectivity.flatten()

		id = self.vtk.vtkIdTypeArray()
		id.SetNumberOfTuples(connectivity.size)
		id.SetNumberOfComponents(1)
		id.SetVoidArray(connectivity, connectivity.size, 1)
		connec = self.vtk.vtkCellArray()
		connec.SetCells(ncel, id)

		points = self.vtk.vtkPoints()
		#points.SetDataTypeToFloat()
		points.SetData(pcoords)

		pdata = self.vtk.vtkPolyData()
		pdata.SetPoints(points)
		pdata.SetLines(connec)

		# Add scalars
		for attribute in attributes:
			# SetScalar allows only one scalar
			#pdata.GetPointData().SetScalars(attribute)
			pdata.GetPointData().AddArray(attribute)

		# The first attribute (first scalar) is the main one
		if len(attributes) > 0:
			pdata.GetPointData().SetActiveScalars(attributes[0].GetName())

		writer = self.vtk.vtkPolyDataWriter()

		# For XML output
		if data_format == "xml":
			writer = self.vtk.vtkXMLDataSetWriter()

		writer.SetFileName(file)
		if float(self.vtk.VTK_VERSION[:3]) < 6:
		    writer.SetInput(pdata)
		else:
		    writer.SetInputData(pdata)
		writer.Write()
