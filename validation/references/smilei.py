# coding: utf-8

"""
    GENERAL DEFINITIONS FOR SMILEI
"""

import math, os, gc, operator

def _add_metaclass(metaclass):
    """Class decorator for creating a class with a metaclass."""
    # Taken from module "six" for compatibility with python 2 and 3
    def wrapper(cls):
        orig_vars = cls.__dict__.copy()
        slots = orig_vars.get('__slots__')
        if slots is not None:
            if isinstance(slots, str): slots = [slots]
            for slots_var in slots: orig_vars.pop(slots_var)
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper


class SmileiComponentType(type):
    """Metaclass to all Smilei components"""
    
    # Constructor of classes
    def __init__(self, name, bases, attrs):
        self._list = []
        self._verify = True
        # Run standard metaclass init
        super(SmileiComponentType, self).__init__(name, bases, attrs)
    
    # Functions to define the iterator
    def __iter__(self):
        self.current = 0
        return self
    def next(self):
        if self.current >= len(self._list):
            raise StopIteration
        self.current += 1
        return self._list[self.current - 1]
    __next__ = next #python3
    
    # Function to return one given instance, for example DiagParticles[0]
    # Special case: species can also be indexed by their name: Species["ion1"]
    def __getitem__(self, key):
        if self.__name__ == "Species" and type(key) is str:
            for obj in self._list:
                if obj.species_type == key:
                    return obj
        else:
            return self._list[key]
    
    # Function to return the number of instances, for example len(Species)
    def __len__(self):
        return len(self._list)
    
    # Function to display the list of instances
    def __repr__(self):
        if len(self._list)==0:
            return "<Empty list of "+self.__name__+">"
        else:
            l = []
            for obj in self._list: l.append(str(obj))
            return "["+", ".join(l)+"]"

@_add_metaclass(SmileiComponentType)
class SmileiComponent(object):
    """Smilei component generic class"""
    
    # Function to initialize new components
    def _fillObjectAndAddToList(self, cls, obj, **kwargs):
        # add all kwargs as internal class variables
        if kwargs is not None:
            for key, value in kwargs.items():
                if key=="_list":
                    print("Python warning: in "+cls.__name__+": cannot have argument named '_list'. Discarding.")
                elif not hasattr(cls, key):
                    raise Exception("ERROR in the namelist: cannot define `"+key+"` in block "+cls.__name__+"()");
                else:
                    setattr(obj, key, value)
        # add the new component to the "_list"
        cls._list.append(obj)
    
    # Constructor for all SmileiComponents
    def __init__(self, **kwargs):
        self._fillObjectAndAddToList(type(self), self, **kwargs)
    
    def __repr__(self):
        return "<Smilei "+type(self).__name__+">"


class SmileiSingletonType(SmileiComponentType):
    """Metaclass to all Smilei singletons"""
    
    def __repr__(self):
        return "<Smilei "+str(self.__name__)+">"

@_add_metaclass(SmileiSingletonType)
class SmileiSingleton(SmileiComponent):
    """Smilei singleton generic class"""
    
    # Prevent having two instances
    def __new__(cls, **kwargs):
        if len(cls._list) >= 1:
            raise Exception("ERROR in the namelist: cannot define block "+cls.__name__+"() twice")
        return super(SmileiSingleton, cls).__new__(cls)
    
    # Constructor for all SmileiSingletons
    def __init__(self, **kwargs):
        self._fillObjectAndAddToList(type(self), type(self), **kwargs)


class Main(SmileiSingleton):
    """Main parameters"""
    
    # Default geometry info
    geometry = None
    cell_length = []
    sim_length = []
    number_of_cells = []
    timestep = None
    sim_time = None
    number_of_timesteps = None
    interpolation_order = 2
    number_of_patches = None
    clrw = 1
    timestep = None
    timestep_over_CFL = None
    
    # Poisson tuning
    solve_poisson = True
    poisson_iter_max = 50000
    poisson_error_max = 1.e-14
    
    # Default fields
    maxwell_sol = 'Yee'
    bc_em_type_x = []
    bc_em_type_y = []
    bc_em_type_z = []
    time_fields_frozen = 0.
    
    # Default Misc
    referenceAngularFrequency_SI = 0.
    print_every = None
    output_dir = None
    random_seed = None
    
    def __init__(self, **kwargs):
        super(Main, self).__init__(**kwargs)
        #initialize timestep if not defined based on timestep_over_CFL
        if Main.timestep is None:
            if Main.timestep_over_CFL is None:
                raise Exception("timestep and timestep_over_CFL not defined")
            else:
                if Main.cell_length is None:
                    raise Exception("Need cell_length to calculate timestep")
                if Main.maxwell_sol == 'Yee':
                    if Main.geometry == '1d3v':
                        Main.timestep = Main.timestep_over_CFL*Main.cell_length[0]
                    elif Main.geometry == '2d3v':         
                        Main.timestep = Main.timestep_over_CFL/math.sqrt(1.0/(Main.cell_length[0]**2)+1.0/(Main.cell_length[1]**2))
                    elif Main.geometry == '3d3v':         
                        Main.timestep = Main.timestep_over_CFL/math.sqrt(1.0/(Main.cell_length[0]**2)+1.0/(Main.cell_length[1]**2)+1.0/(Main.cell_length[2]**2))
                    else: 
                        raise Exception("timestep: geometry not implemented "+Main.geometry)
                else:
                    raise Exception("timestep: maxwell_sol not implemented "+Main.maxwell_sol)

        #initialize sim_length if not defined based on number_of_cells and cell_length
        if len(Main.sim_length) is 0:
            if len(Main.number_of_cells) is 0:
                raise Exception("sim_length and number_of_cells not defined")
            elif len(Main.number_of_cells) != len(Main.cell_length):
                raise Exception("sim_length and number_of_cells not defined")
            else :
                Main.sim_length = [a*b for a,b in zip(Main.number_of_cells, Main.cell_length)]

        #initialize sim_time if not defined based on number_of_timesteps and timestep
        if Main.sim_time is None:
            if Main.number_of_timesteps is None:
                raise Exception("sim_time and number_of_timesteps not defined")
            else:
                Main.sim_time = Main.number_of_timesteps * Main.timestep

class LoadBalancing(SmileiSingleton):
    """Load balancing parameters"""
    
    every = 150
    initial_balance = True
    coef_cell = 1.0
    coef_frozen = 0.1


class MovingWindow(SmileiSingleton):
    """Moving window parameters"""
    
    time_start = 0.
    velocity_x = 1.


class DumpRestart(SmileiSingleton):
    """Dump and restart parameters"""
    
    restart_dir = None
    restart_number = None
    dump_step = 0
    dump_minutes = 0.
    dump_file_sequence = 2
    dump_deflate = 0
    exit_after_dump = True
    file_grouping = None


class Species(SmileiComponent):
    """Species parameters"""
    species_type = None
    initPosition_type = None
    initMomentum_type = ""
    n_part_per_cell = None
    c_part_max = 1.0
    mass = None
    charge = None
    charge_density = None
    nb_density = None
    mean_velocity = [0.]
    temperature = [1e-10]
    thermT = None
    thermVelocity = None
    dynamics_type = "norm"
    time_frozen = 0.0
    radiating = False
    bc_part_type_xmin = None
    bc_part_type_xmax = None
    bc_part_type_ymax = None
    bc_part_type_ymin = None
    bc_part_type_zmin = None
    bc_part_type_zmax = None
    ionization_model = "none"
    ionization_electrons = None
    atomic_number = None
    isTest = False
    track_every = 0
    track_ordered = False
    track_flush_every = 1

class Laser(SmileiComponent):
    """Laser parameters"""
    boxSide = "xmin"
    omega = 1.
    chirp_profile = 1.
    time_envelope = 1.
    space_envelope = [1., 0.]
    phase = [0., 0.]
    space_time_profile = None

class Collisions(SmileiComponent):
    """Collisions parameters"""
    species1 = None
    species2 = None
    coulomb_log = 0.
    debug_every = 0
    ionizing = False


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    every = None
    number = []
    pos = []
    pos_first = []
    pos_second = []
    pos_third = []
    fields = []
    flush_every = 1

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    output = None
    every = None
    time_average = 1
    species = None
    axes = []
    flush_every = 1

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    every = None
    precision = 10
    vars = []

class DiagFields(SmileiComponent):
    """Diagnostic Fields"""
    every = None
    fields = []
    time_average = 1
    flush_every = 1

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    field = None
    profile = None

# external current (antenna)
class Antenna(SmileiComponent):
    """Antenna"""
    field = None
    time_profile  = None
    space_profile = None

# Particle wall
class PartWall(SmileiComponent):
    """Particle Wall"""
    kind = "none"
    x = None
    y = None
    z = None

# Smilei-defined
smilei_mpi_rank = 0
smilei_mpi_size = 1
smilei_rand_max = 2**31-1
smilei_version='v3.0-84-g92ff763-CI'
# Some predefined profiles (see doc)

def constant(value, xvacuum=-float("inf"), yvacuum=-float("inf"), zvacuum=-float("inf")):
    global Main
    if len(Main)==0:
        raise Exception("constant profile has been defined before `Main()`")
    if Main.geometry == "1d3v":
        f = lambda x  : value if x>=xvacuum else 0.
    if Main.geometry == "2d3v":
        f = lambda x,y: value if (x>=xvacuum and y>=yvacuum) else 0.
        f.yvacuum = yvacuum
    if Main.geometry == "3d3v":
        f = lambda x,y,z: value if (x>=xvacuum and y>=yvacuum and z>=zvacuum) else 0.
        f.yvacuum = yvacuum
        f.zvacuum = zvacuum
    f.profileName = "constant"
    f.value   = value
    f.xvacuum = xvacuum
    return f
constant._reserved = True

def trapezoidal(max,
                xvacuum=0., xplateau=None, xslope1=0., xslope2=0.,
                yvacuum=0., yplateau=None, yslope1=0., yslope2=0.,
                zvacuum=0., zplateau=None, zslope1=0., zslope2=0. ):
    global Main
    if len(Main)==0:
        raise Exception("trapezoidal profile has been defined before `Main()`")
    if len(Main.sim_length)>0 and xplateau is None: xplateau = Main.sim_length[0]-xvacuum
    if len(Main.sim_length)>1 and yplateau is None: yplateau = Main.sim_length[1]-yvacuum
    if len(Main.sim_length)>2 and zplateau is None: zplateau = Main.sim_length[2]-zvacuum
    def trapeze(max, vacuum, plateau, slope1, slope2):
        def f(position):
            # vacuum region
            if position < vacuum: return 0.
            # linearly increasing density
            elif position < vacuum+slope1: return max*(position-vacuum) / slope1
            # density plateau
            elif position < vacuum+slope1+plateau: return max
            # linearly decreasing density
            elif position < vacuum+slope1+plateau+slope2:
                return max*(1. - ( position - (vacuum+slope1+plateau) ) / slope2)
            # beyond the plasma
            else: return 0.0
        return f
    if   Main.geometry == "1d3v": dim = 1
    elif Main.geometry == "2d3v": dim = 2
    elif Main.geometry == "3d3v": dim = 3
    fx = trapeze(max, xvacuum, xplateau, xslope1, xslope2)
    f = fx
    if dim > 1:
        fy = trapeze(1. , yvacuum, yplateau, yslope1, yslope2)
        f = lambda x,y: fx(x)*fy(y)
    if dim > 2:
        fz = trapeze(1. , zvacuum, zplateau, zslope1, zslope2)
        f = lambda x,y,z: fx(x)*fy(y)*fz(z)
    f.profileName = "trapezoidal"
    f.value    = max
    f.xvacuum  = xvacuum
    f.xplateau = xplateau
    f.xslope1  = xslope1
    f.xslope2  = xslope2
    if dim > 1:
        f.yvacuum  = yvacuum
        f.yplateau = yplateau
        f.yslope1  = yslope1
        f.yslope2  = yslope2
    if dim > 2:
        f.zvacuum  = zvacuum
        f.zplateau = zplateau
        f.zslope1  = zslope1
        f.zslope2  = zslope2
    return f

def gaussian(max,
             xvacuum=0., xlength=float("inf"), xfwhm=None, xcenter=None, xorder=2,
             yvacuum=0., ylength=float("inf"), yfwhm=None, ycenter=None, yorder=2,
             zvacuum=0., zlength=float("inf"), zfwhm=None, zcenter=None, zorder=2 ):
    import math
    global Main
    if len(Main)==0:
        raise Exception("gaussian profile has been defined before `Main()`")
    if len(Main.sim_length)>0:
        if xlength is None: xlength = Main.sim_length[0]-xvacuum
        if xfwhm   is None: xfwhm   = (Main.sim_length[0]-xvacuum)/3.
        if xcenter is None: xcenter = xvacuum + (Main.sim_length[0]-xvacuum)/2.
    if len(Main.sim_length)>1: 
        if ylength is None: ylength = Main.sim_length[1]-yvacuum
        if yfwhm   is None: yfwhm   = (Main.sim_length[1]-yvacuum)/3.
        if ycenter is None: ycenter = yvacuum + (Main.sim_length[1]-yvacuum)/2.
    if len(Main.sim_length)>2: 
        if zlength is None: zlength = Main.sim_length[2]-zvacuum
        if zfwhm   is None: zfwhm   = (Main.sim_length[2]-zvacuum)/3.
        if zcenter is None: zcenter = zvacuum + (Main.sim_length[2]-zvacuum)/2.
    def gauss(max, vacuum, length, sigma, center, order):
        def f(position):
            if order == 0: return max
            # vacuum region
            if position < vacuum: return 0.
            # gaussian
            elif position < vacuum+length: return max*math.exp( -(position-center)**order / sigma )
            # beyond
            else: return 0.0
        return f
    if Main.geometry == "1d3v": dim = 1
    if Main.geometry == "2d3v": dim = 2
    if Main.geometry == "3d3v": dim = 3
    xsigma = (0.5*xfwhm)**xorder/math.log(2.0)
    fx = gauss(max, xvacuum, xlength, xsigma, xcenter, xorder)
    f = fx
    if dim > 1:
        ysigma = (0.5*yfwhm)**yorder/math.log(2.0)
        fy = gauss(1., yvacuum, ylength, ysigma, ycenter, yorder)
        f = lambda x,y: fx(x)*fy(y)
    if dim > 2:
        zsigma = (0.5*zfwhm)**zorder/math.log(2.0)
        fz = gauss(1., zvacuum, zlength, zsigma, zcenter, zorder)
        f = lambda x,y,z: fx(x)*fy(y)*fz(z)
    f.profileName = "gaussian"
    f.value   = max
    f.xvacuum = xvacuum
    f.xlength = xlength
    f.xsigma  = xsigma
    f.xcenter = xcenter
    f.xorder  = xorder
    if dim > 1:
        f.yvacuum = yvacuum
        f.ylength = ylength
        f.ysigma  = ysigma
        f.ycenter = ycenter
        f.yorder  = yorder
    if dim > 2:
        f.zvacuum = zvacuum
        f.zlength = zlength
        f.zsigma  = zsigma
        f.zcenter = zcenter
        f.zorder  = zorder
    return f


def polygonal(xpoints=[], xvalues=[]):
    global Main
    if len(Main)==0:
        raise Exception("polygonal profile has been defined before `Main()`")
    if len(xpoints)!=len(xvalues):
        raise Exception("polygonal profile requires as many points as values")
    if len(Main.sim_length)>0 and len(xpoints)==0:
        xpoints = [0., Main.sim_length[0]]
        xvalues = [1., 1.]
    N = len(xpoints)
    xpoints = [float(x) for x in xpoints]
    xvalues = [float(x) for x in xvalues]
    xslopes = [0. for i in range(1,N)]
    for i in range(1,N):
        if xpoints[i] == xpoints[i-1]: continue
        xslopes[i-1] = (xvalues[i]-xvalues[i-1])/(xpoints[i]-xpoints[i-1])
    def f(x,y=0.):
        if x < xpoints[0]: return 0.0;
        for i in range(1,N):
            if x<xpoints[i]: return xvalues[i-1] + xslopes[i-1] * ( x-xpoints[i-1] )
        return 0.
    f.profileName = "polygonal"
    f.xpoints = xpoints
    f.xvalues = xvalues
    f.xslopes = xslopes
    return f

def cosine(base,
           xamplitude=1., xvacuum=0., xlength=None, xphi=0., xnumber=2,
           yamplitude=1., yvacuum=0., ylength=None, yphi=0., ynumber=2,
           zamplitude=1., zvacuum=0., zlength=None, zphi=0., znumber=2):
    import math
    global Main
    if len(Main)==0:
        raise Exception("cosine profile has been defined before `Main()`")
    
    if len(Main.sim_length)>0 and xlength is None: xlength = Main.sim_length[0]-xvacuum
    if len(Main.sim_length)>1 and ylength is None: ylength = Main.sim_length[1]-yvacuum
    if len(Main.sim_length)>2 and zlength is None: zlength = Main.sim_length[2]-zvacuum
    
    def cos(base, amplitude, vacuum, length, phi, number):
        def f(position):
            #vacuum region
            if position < vacuum: return 0.
            # profile region
            elif position < vacuum+length:
                return base + amplitude * math.cos(phi + 2.*math.pi * number * (position-vacuum)/length)
            # beyond
            else: return 0.
        return f
    if Main.geometry == "1d3v": dim = 1
    if Main.geometry == "2d3v": dim = 2
    if Main.geometry == "3d3v": dim = 3
    fx = cos(base, xamplitude, xvacuum, xlength, xphi, xnumber)
    f = fx
    if dim > 1:
        fy = cos(base, yamplitude, yvacuum, ylength, yphi, ynumber)
        f = lambda x,y: fx(x)*fy(y)
    if dim > 2:
        fz = cos(base, zamplitude, zvacuum, zlength, zphi, znumber)
        f = lambda x,y,z: fx(x)*fy(y)*fz(z)
    f.profileName = "cosine"
    f.base        = base
    f.xamplitude  = xamplitude
    f.xvacuum     = xvacuum
    f.xlength     = xlength
    f.xphi        = xphi
    f.xnumber     = float(xnumber)
    if dim > 1:
        f.yamplitude  = yamplitude
        f.yvacuum     = yvacuum
        f.ylength     = ylength
        f.yphi        = yphi
        f.ynumber     = float(ynumber)
    if dim > 2:
        f.zamplitude  = zamplitude
        f.zvacuum     = zvacuum
        f.zlength     = zlength
        f.zphi        = zphi
        f.znumber     = float(znumber)
    return f

def polynomial(**kwargs):
    global Main
    if len(Main)==0:
        raise Exception("polynomial profile has been defined before `Main()`")
    x0 = 0.
    y0 = 0.
    z0 = 0.
    coeffs = dict()
    for k, a in kwargs.items():
        if   k=="x0":
            x0 = a
        elif k=="y0":
            y0 = a
        elif k[:5]=="order":
            if type(a) is not list: a = [a]
            order = int(k[5:])
            coeffs[ order ] = a
            if Main.geometry=="1d3v":
                if len(a)!=1:
                    raise Exception("1D polynomial profile must have one coefficient at order "+str(order))
            elif Main.geometry=="2d3v":
                if len(a)!=order+1:
                    raise Exception("2D polynomial profile must have "+str(order+1)+" coefficients at order "+str(order))
            elif Main.geometry=="3d3v":
                if len(a)!=(order+1)*(order+2)/2:
                    raise Exception("3D polynomial profile must have "+str((order+1)*(order+2)/2)+" coefficients at order "+str(order))
    if Main.geometry=="1d3v":
        def f(x):
            r = 0.
            xx0 = x-x0
            xx = 1.
            currentOrder = 0
            for order, c in sorted(coeffs.items()):
                while currentOrder<order:
                    currentOrder += 1
                    xx *= xx0
                r += c[0] * xx
            return r
    elif Main.geometry=="2d3v":
        def f(x,y):
            r = 0.
            xx0 = x-x0
            yy0 = y-y0
            xx = [1.]
            currentOrder = 0
            for order, c in sorted(coeffs.items()):
                while currentOrder<order:
                    currentOrder += 1
                    yy = xx[-1]*yy0
                    xx = [ xxx * xx0 for xxx in xx ] . append(yy)
                for i in range(order+1): r += c[i]*xx[i]
            return r
    elif Main.geometry=="3d3v":
        def f(x,y,z):
            r = 0.
            xx0 = x-x0
            yy0 = y-y0
            zz0 = z-z0
            xx = [1.]
            currentOrder = 0
            for order, c in sorted(coeffs.items()):
                while currentOrder<order:
                    currentOrder += 1
                    zz = xx[-1]*zz0
                    yy = [ xxx * yy0 for xxx in xx[-currentOrder-1:] ] . append(zz)
                    xx = [ xxx * xx0 for xxx in xx ] . extend(yy)
                for i in range(order+1): r += c[i]*xx[i]
            return r
    else:
        raise Exception("polynomial profiles are not available in this geometry yet")
    f.profileName = "polynomial"
    f.x0 = x0
    f.y0 = y0
    f.z0 = z0
    f.orders = []
    f.coeffs = []
    for order, c in sorted(coeffs.items()):
        f.orders.append( order )
        f.coeffs.append( c     )
    return f



def tconstant(start=0.):
    def f(t):
        return 1. if t>=start else 0.
    f.profileName = "tconstant"
    f.start       = start
    return f
tconstant._reserved = True

def ttrapezoidal(start=0., plateau=None, slope1=0., slope2=0.):
    global Main
    if len(Main)==0:
        raise Exception("ttrapezoidal profile has been defined before `Main()`")
    if plateau is None: plateau = Main.sim_time - start
    def f(t):
        if t < start: return 0.
        elif t < start+slope1: return (t-start) / slope1
        elif t < start+slope1+plateau: return 1.
        elif t < start+slope1+plateau+slope2:
            return 1. - ( t - (start+slope1+slope2) ) / slope2
        else: return 0.0
    f.profileName = "ttrapezoidal"
    f.start       = start
    f.plateau     = plateau
    f.slope1      = slope1
    f.slope2      = slope2
    return f

def tgaussian(start=0., duration=None, fwhm=None, center=None, order=2):
    import math
    global Main
    if len(Main)==0:
        raise Exception("tgaussian profile has been defined before `Main()`")
    if duration is None: duration = Main.sim_time-start
    if fwhm     is None: fwhm     = (Main.sim_time-start)/3.
    if center   is None: center   = start + (Main.sim_time-start)/2.
    sigma = (0.5*fwhm)**order/math.log(2.0)
    def f(t):
        if t < start: return 0.
        elif t < start+duration: return math.exp( -(t-center)**order / sigma )
        else: return 0.0
    f.profileName = "tgaussian"
    f.start       = start
    f.duration    = duration
    f.sigma       = sigma
    f.center      = center
    f.order       = order
    return f

def tpolygonal(points=[], values=[]):
    global Main
    if len(Main)==0:
        raise Exception("tpolygonal profile has been defined before `Main()`")
    if len(points)==0:
        points = [0., Main.sim_time]
        values = [1., 1.]
    N = len(points)
    points = [float(x) for x in points]
    values = [float(x) for x in values]
    slopes = [0. for i in range(1,N)]
    for i in range(1,N):
        if points[i] == points[i-1]: continue
        slopes[i-1] = (values[i]-values[i-1])/(points[i]-points[i-1])
    def f(t):
        if t < points[0]: return 0.0;
        for i in range(1,N):
            if t<points[i]: return values[i-1] + slopes[i-1] * ( t-points[i-1] )
        return 0.
    f.profileName = "tpolygonal"
    f.points      = points
    f.values      = values
    f.slopes      = slopes
    return f

def tcosine(base=0., amplitude=1., start=0., duration=None, phi=0., freq=1.):
    import math
    global Main
    if len(Main)==0:
        raise Exception("tcosine profile has been defined before `Main()`")
    if duration is None: duration = Main.sim_time-start
    def f(t):
        if t < start: return 0.
        elif t < start+duration:
            return base + amplitude * math.cos(phi + freq*(t-start))
        else: return 0.
    f.profileName = "tcosine"
    f.base        = base
    f.amplitude   = amplitude
    f.start       = start
    f.duration    = duration
    f.phi         = phi
    f.freq        = freq
    return f

def tpolynomial(**kwargs):
    t0 = 0.
    coeffs = dict()
    for k, a in kwargs.items():
        if   k=="t0":
            t0 = a
        elif k[:5]=="order":
            order = int(k[5:])
            try: coeffs[ order ] = a*1.
            except: raise Exception("tpolynomial profile must have one coefficient per order")
    def f(t):
        r = 0.
        tt0 = t-t0
        tt = 1.
        currentOrder = 0
        for order, c in sorted(coeffs.items()):
            while currentOrder<order:
                currentOrder += 1
                tt *= tt0
            r += c * tt
        return r
    f.profileName = "tpolynomial"
    f.t0 = t0
    f.orders = []
    f.coeffs = []
    for order, c in sorted(coeffs.items()):
        f.orders.append( order )
        f.coeffs.append( c     )
    return f


def transformPolarization(polarizationPhi, ellipticity):
    import math
    p = (1.-ellipticity**2)*math.sin(2.*polarizationPhi)/2.
    if abs(p) < 1e-10:
        if abs(ellipticity**2-1.)<1e-10: polarizationPhi=0.
        dephasing = math.pi/2.
        amplitude = math.sqrt(1./(1.+ellipticity**2))
        amplitudeY = amplitude * (math.cos(polarizationPhi)+math.sin(polarizationPhi)*ellipticity)
        amplitudeZ = amplitude * (math.sin(polarizationPhi)+math.cos(polarizationPhi)*ellipticity)
    else:
        dephasing = math.atan(ellipticity/p)
        theta = 0.5 * math.atan( math.tan(2.*polarizationPhi) / math.cos(dephasing) )
        while theta<0.: theta += math.pi/2.
        amplitudeY = math.sqrt(2.) * math.cos(theta)
        amplitudeZ = math.sqrt(2.) * math.sin(theta)
    return [dephasing, amplitudeY, amplitudeZ]


def LaserPlanar1D( boxSide="xmin", a0=1., omega=1.,
        polarizationPhi=0., ellipticity=0., time_envelope=tconstant()):
    import math
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarizationPhi, ellipticity)
    amplitudeY *= a0
    amplitudeZ *= a0
    # Create Laser
    Laser(
        boxSide        = boxSide,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ amplitudeZ, amplitudeY ],
        phase          = [ dephasing, 0. ],
    )



def LaserGaussian2D( boxSide="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=0.,
        polarizationPhi=0., ellipticity=0., time_envelope=tconstant()):
    import math
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarizationPhi, ellipticity)
    amplitudeY *= a0
    amplitudeZ *= a0
    # Space and phase envelopes
    Zr = omega * waist**2/2.
    phaseZero = 0.
    if incidence_angle == 0.:
        Y1 = focus[1]
        w  = math.sqrt(1./(1.+(focus[0]/Zr)**2))
        invWaist2 = (w/waist)**2
        coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
        def spatial(y):
            return w * math.exp( -invWaist2*(y-focus[1])**2 )
        def phase(y):
            return coeff * (y-focus[1])**2
    else:
        invZr  = math.sin(incidence_angle) / Zr
        invZr2 = invZr**2
        invZr3 = (math.cos(incidence_angle) / Zr)**2 / 2.
        invWaist2 = (math.cos(incidence_angle) / waist)**2
        omega_ = omega * math.sin(incidence_angle)
        Y1 = focus[1] + focus[0]/math.tan(incidence_angle)
        Y2 = focus[1] - focus[0]*math.tan(incidence_angle)
        def spatial(y):
            w2 = 1./(1. + invZr2*(y-Y1)**2)
            return math.sqrt(w2) * math.exp( -invWaist2*w2*(y-Y2)**2 )
        def phase(y):
            dy = y-Y1
            return omega_*dy*(1.+ invZr3*(y-Y2)**2/(1.+invZr2*dy**2)) + math.atan(invZr*dy)
        phaseZero = phase(Y2)
    # Create Laser
    Laser(
        boxSide        = boxSide,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y:amplitudeZ*spatial(y), lambda y:amplitudeY*spatial(y) ],
        phase          = [ lambda y:phase(y)-phaseZero+dephasing, lambda y:phase(y)-phaseZero ],
    )

def LaserGaussian3D( boxSide="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=[0.,0.],
        polarizationPhi=0., ellipticity=0., time_envelope=tconstant()):
    import math
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarizationPhi, ellipticity)
    amplitudeY *= a0
    amplitudeZ *= a0
    # Space and phase envelopes
    Zr = omega * waist**2/2.
    phaseZero = 0.
    if incidence_angle == [0.,0.]:
        w  = math.sqrt(1./(1.+(focus[0]/Zr)**2))
        invWaist2 = (w/waist)**2
        coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
        def spatial(y,z):
            return w * math.exp( -invWaist2*((y-focus[1])**2 + (z-focus[2])**2 )  )
        def phase(y,z):
            return coeff * ( (y-focus[1])**2 + (z-focus[2])**2 )
    else:
        invZr = 1./Zr
        invW  = 1./waist
        alpha = omega * Zr
        cy = math.cos(incidence_angle[0]); sy = math.sin(incidence_angle[0])
        cz = math.cos(incidence_angle[1]); sz = math.sin(incidence_angle[1])
        cycz = cy*cz; cysz = cy*sz; sycz = sy*cz; sysz = sy*sz
        def spatial(y,z):
            X = invZr * (-focus[0]*cycz + (y-focus[1])*cysz - (z-focus[2])*sy )
            Y = invW  * ( focus[0]*sz   + (y-focus[1])*cz                     )
            Z = invW  * (-focus[0]*sycz + (y-focus[1])*sysz + (z-focus[2])*cy )
            invW2 = 1./(1.+X**2)
            return math.sqrt(invW2) * math.exp(-(Y**2+Z**2)*invW2)
        def phase(y,z):
            X = invZr * (-focus[0]*cycz + (y-focus[1])*cysz - (z-focus[2])*sy )
            Y = invZr * ( focus[0]*sz   + (y-focus[1])*cz                     )
            Z = invZr * (-focus[0]*sycz + (y-focus[1])*sysz + (z-focus[2])*cy )
            return alpha * X*(1.+0.5*(Y**2+Z**2)/(1.+X**2)) - math.atan(X)
        phaseZero = phase(focus[1]-sz/cz*focus[0], focus[2]+sy/cy/cz*focus[0])
        
    # Create Laser
    Laser(
        boxSide        = boxSide,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y,z:amplitudeZ*spatial(y,z), lambda y,z:amplitudeY*spatial(y,z) ],
        phase          = [ lambda y,z:phase(y,z)-phaseZero+dephasing, lambda y,z:phase(y,z)-phaseZero ],
    )

"""
-----------------------------------------------------------------------
    BEGINNING OF THE USER NAMELIST
"""

dx = 0.125
dtrans = 3.
dt = 0.124
nx = 896
ntrans = 40
Lx = nx * dx
Ltrans = ntrans*dtrans
npatch_x = 128
laser_fwhm = 19.80

Main(
    geometry = "3d3v",
    
    interpolation_order = 2,
    
    timestep = dt,
    sim_time = int(2*Lx/dt)*dt,
    
    cell_length  = [dx, dtrans, dtrans],
    sim_length = [ Lx,  Ltrans, Ltrans],
    
    number_of_patches = [npatch_x, 4, 4],
    
    clrw = nx/npatch_x,
    
    bc_em_type_x = ["silver-muller","silver-muller"],
    bc_em_type_y = ["silver-muller","silver-muller"],
    bc_em_type_z = ["silver-muller","silver-muller"],
    
    random_seed = 0,
    solve_poisson = False,
    print_every = 100
)

MovingWindow(
    time_start = Main.sim_length[0],
    velocity_x = 0.9997
)

LoadBalancing(
    initial_balance = False,
    every = 20,
    coef_cell = 1.,
    coef_frozen = 0.1
)

Species( 
    species_type = "electron",
    initPosition_type = "regular",
    initMomentum_type = "cold",
    n_part_per_cell = 8,
    c_part_max = 1.0,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.000494,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.0],
    dynamics_type = "norm",    
    time_frozen = 0.0,
    radiating = False,
    bc_part_type_xmin = "supp",
    bc_part_type_xmax = "supp",
    bc_part_type_ymin ="supp",
    bc_part_type_ymax ="supp",
    bc_part_type_zmin ="supp",
    bc_part_type_zmax ="supp"
)

LaserGaussian3D(
    boxSide         = "xmin",
    a0              = 2.,
    focus           = [0., Main.sim_length[1]/2., Main.sim_length[2]/2.],
    waist           = 26.16,
    time_envelope   = tgaussian(center=2**0.5*laser_fwhm, fwhm=laser_fwhm)
)

DumpRestart(
    dump_step = 0,
    dump_minutes = 0.0,
    exit_after_dump = False,
)

list_fields = ['Ex','Ey','Rho','Jx']

DiagFields(
    every = 100,
    fields = list_fields
)

DiagProbe(
	every = 10,
	pos = [0., Main.sim_length[1]/2., Main.sim_length[2]/2.],
	pos_first = [Main.sim_length[0], Main.sim_length[1]/2., Main.sim_length[2]/2.],
	number = [nx],
	fields = ['Ex','Ey','Rho','Jx']
)

DiagProbe(
	every = 10,
	pos = [0., Main.sim_length[1]/4., Main.sim_length[2]/2.],
	pos_first = [0., 3*Main.sim_length[1]/4., Main.sim_length[2]/2.],
	pos_second = [Main.sim_length[0], Main.sim_length[1]/4., Main.sim_length[2]/2.],
	number = [nx, ntrans],
	fields = ['Ex','Ey','Rho','Jx']
)

DiagScalar(every = 10, vars=['Uelm','Ukin_electron','ExMax','ExMaxCell','EyMax','EyMaxCell', 'RhoMin', 'RhoMinCell'])



"""
    END OF THE USER NAMELIST
-----------------------------------------------------------------------
"""

gc.collect()
import math

def _mkdir(role, path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except:
            raise Exception("ERROR in the namelist: "+role+" "+path+" cannot be created")
    elif not os.path.isdir(path):
        raise Exception("ERROR in the namelist: "+role+" "+path+" exists but is not a directory")

def _smilei_check():
    """Do checks over the script"""
    # Verify classes were not overriden
    for CheckClassName in ["SmileiComponent","Species", "Laser","Collisions",
            "DiagProbe","DiagParticles", "DiagScalar","DiagFields","ExtField",
            "SmileiSingleton","Main","DumpRestart","LoadBalancing","MovingWindow"]:
        CheckClass = globals()[CheckClassName]
        try:
            if not CheckClass._verify: raise Exception("")
        except:
            raise Exception("ERROR in the namelist: it seems that the name `"+CheckClassName+"` has been overriden")
    # Verify the output_dir
    if smilei_mpi_rank == 0 and Main.output_dir:
        _mkdir("output_dir", Main.output_dir)
    # Checkpoint: prepare dir tree
    if smilei_mpi_rank == 0 and (DumpRestart.dump_step>0 or DumpRestart.dump_minutes>0.):
        checkpoint_dir = (Main.output_dir or ".") + os.sep + "checkpoints" + os.sep
        if DumpRestart.file_grouping :
            ngroups = smilei_mpi_size/DumpRestart.file_grouping+1
            ngroups_chars = int(math.log10(ngroups))+1
            for group in range(ngroups):
                group_dir = checkpoint_dir + '%*s'%(ngroups_chars,group)
                _mkdir("checkpoint", group_dir)
        else:
            _mkdir("checkpoint", checkpoint_dir)
    # Checkpoint: Verify the restart_dir
    if len(DumpRestart)==1 and DumpRestart.restart_dir:
        if not os.path.isdir(DumpRestart.restart_dir):
            raise Exception("ERROR in the namelist: restart_dir = `"+DumpRestart.restart_dir+"` is not a directory")
    # Verify that constant() and tconstant() were not redefined
    if not hasattr(constant, "_reserved") or not hasattr(tconstant, "_reserved"):
        raise Exception("Names `constant` and `tconstant` cannot be overriden")
    # Convert float profiles to constant() or tconstant()
    def toSpaceProfile(input):
        try   : return constant(input*1.)
        except: return input
    def toTimeProfile(input):
        try:
            input*1.
            return tconstant()
        except: return input
    for s in Species:
        s.nb_density      = toSpaceProfile(s.nb_density      )
        s.charge_density  = toSpaceProfile(s.charge_density  )
        s.n_part_per_cell = toSpaceProfile(s.n_part_per_cell )
        s.charge          = toSpaceProfile(s.charge          )
        s.mean_velocity   = [ toSpaceProfile(p) for p in s.mean_velocity ]
        s.temperature     = [ toSpaceProfile(p) for p in s.temperature   ]
    for e in ExtField:
        e.profile         = toSpaceProfile(e.profile         )
    for a in Antenna:
        a.space_profile   = toSpaceProfile(a.space_profile   )
        a.time_profile    = toTimeProfile (a.time_profile    )
    for l in Laser:
        l.chirp_profile   = toTimeProfile( l.chirp_profile )
        l.time_envelope   = toTimeProfile( l.time_envelope )
        l.space_envelope  = [ toSpaceProfile(p) for p in l.space_envelope ]
        l.phase           = [ toSpaceProfile(p) for p in l.phase          ]

# this function will be called after initialising the simulation, just before entering the time loop
# if it returns false, the code will call a Py_Finalize();
def _keep_python_running():
    ps = [[las.time_envelope, las.chirp_profile] for las in Laser]
    ps += [[ant.time_profile] for ant in Antenna]
    if len(MovingWindow)>0 or len(LoadBalancing)>0:
        ps += [[s.nb_density, s.charge_density, s.n_part_per_cell, s.charge] + s.mean_velocity + s.temperature for s in Species]
    profiles = []
    for p in ps: profiles += p
    for prof in profiles:
        if callable(prof) and not hasattr(prof,"profileName"):
            return True
    return False

# Prevent creating new components (by mistake)
def _noNewComponents(cls, *args, **kwargs):
    print("Please do not create a new "+cls.__name__)
    return None
SmileiComponent.__new__ = staticmethod(_noNewComponents)



