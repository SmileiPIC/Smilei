# Some predefined profiles (see doc)

def constant(value, xvacuum=-float("inf"), yvacuum=-float("inf"), zvacuum=-float("inf")):
    global Main
    if len(Main)==0:
        raise Exception("constant profile has been defined before `Main()`")
    if Main.geometry == "1Dcartesian":
        f = lambda x  : value if x>=xvacuum else 0.
    if Main.geometry == "2Dcartesian":
        f = lambda x,y: value if (x>=xvacuum and y>=yvacuum) else 0.
        f.yvacuum = yvacuum
    if Main.geometry == "3Dcartesian":
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
    if len(Main.grid_length)>0 and xplateau is None: xplateau = Main.grid_length[0]-xvacuum
    if len(Main.grid_length)>1 and yplateau is None: yplateau = Main.grid_length[1]-yvacuum
    if len(Main.grid_length)>2 and zplateau is None: zplateau = Main.grid_length[2]-zvacuum
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
    if   Main.geometry == "1Dcartesian": dim = 1
    elif Main.geometry == "2Dcartesian": dim = 2
    elif Main.geometry == "3Dcartesian": dim = 3
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
    if len(Main.grid_length)>0:
        if xlength is None: xlength = Main.grid_length[0]-xvacuum
        if xfwhm   is None: xfwhm   = (Main.grid_length[0]-xvacuum)/3.
        if xcenter is None: xcenter = xvacuum + (Main.grid_length[0]-xvacuum)/2.
    if len(Main.grid_length)>1: 
        if ylength is None: ylength = Main.grid_length[1]-yvacuum
        if yfwhm   is None: yfwhm   = (Main.grid_length[1]-yvacuum)/3.
        if ycenter is None: ycenter = yvacuum + (Main.grid_length[1]-yvacuum)/2.
    if len(Main.grid_length)>2: 
        if zlength is None: zlength = Main.grid_length[2]-zvacuum
        if zfwhm   is None: zfwhm   = (Main.grid_length[2]-zvacuum)/3.
        if zcenter is None: zcenter = zvacuum + (Main.grid_length[2]-zvacuum)/2.
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
    if Main.geometry == "1Dcartesian": dim = 1
    if Main.geometry == "2Dcartesian": dim = 2
    if Main.geometry == "3Dcartesian": dim = 3
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
    if len(Main.grid_length)>0 and len(xpoints)==0:
        xpoints = [0., Main.grid_length[0]]
        xvalues = [1., 1.]
    N = len(xpoints)
    xpoints = [float(x) for x in xpoints]
    xvalues = [float(x) for x in xvalues]
    xslopes = [0. for i in range(1,N)]
    for i in range(1,N):
        if xpoints[i] == xpoints[i-1]: continue
        xslopes[i-1] = (xvalues[i]-xvalues[i-1])/(xpoints[i]-xpoints[i-1])
    def f(x,y=0.,z=0.):
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
    
    if len(Main.grid_length)>0 and xlength is None: xlength = Main.grid_length[0]-xvacuum
    if len(Main.grid_length)>1 and ylength is None: ylength = Main.grid_length[1]-yvacuum
    if len(Main.grid_length)>2 and zlength is None: zlength = Main.grid_length[2]-zvacuum
    
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
    if Main.geometry == "1Dcartesian": dim = 1
    if Main.geometry == "2Dcartesian": dim = 2
    if Main.geometry == "3Dcartesian": dim = 3
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
        elif k=="z0":
            z0 = a
        elif k[:5]=="order":
            if type(a) is not list: a = [a]
            order = int(k[5:])
            coeffs[ order ] = a
            if Main.geometry=="1Dcartesian":
                if len(a)!=1:
                    raise Exception("1D polynomial profile must have one coefficient at order "+str(order))
            elif Main.geometry=="2Dcartesian":
                if len(a)!=order+1:
                    raise Exception("2D polynomial profile must have "+str(order+1)+" coefficients at order "+str(order))
            elif Main.geometry=="3Dcartesian":
                if len(a)!=(order+1)*(order+2)/2:
                    raise Exception("3D polynomial profile must have "+str((order+1)*(order+2)/2)+" coefficients at order "+str(order))
    if Main.geometry=="1Dcartesian":
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
    elif Main.geometry=="2Dcartesian":
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
                    xx = [ xxx * xx0 for xxx in xx ] + [yy]
                for i in range(order+1): r += c[i]*xx[i]
            return r
    elif Main.geometry=="3Dcartesian":
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
                    yy = [ xxx * yy0 for xxx in xx[-currentOrder-1:] ] + [zz]
                    xx = [ xxx * xx0 for xxx in xx ] + yy
                for i in range(len(c)): r += c[i]*xx[i]
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
    if plateau is None: plateau = Main.simulation_time - start
    def f(t):
        if t < start: return 0.
        elif t < start+slope1: return (t-start) / slope1
        elif t < start+slope1+plateau: return 1.
        elif t < start+slope1+plateau+slope2:
            return 1. - ( t - (start+slope1+plateau) ) / slope2
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
    if duration is None: duration = Main.simulation_time-start
    if fwhm     is None: fwhm     = (Main.simulation_time-start)/3.
    if center   is None: center   = start + (Main.simulation_time-start)/2.
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
        points = [0., Main.simulation_time]
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
    if duration is None: duration = Main.simulation_time-start
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

def tsin2plateau(start=0., fwhm=0., plateau=None, slope1=None, slope2=None):
    import math
    global Main
    if len(Main)==0:
        raise Exception("tsin2plateau profile has been defined before `Main()`")
    if plateau is None: plateau = 0 # default is a simple sin2 profile (could be used for a 2D or 3D laserPulse too)
    if slope1 is None: slope1 = fwhm
    if slope2 is None: slope2 = slope1
    def f(t):
        if t < start:
            return 0.
        elif (t < start+slope1) and (slope1!=0.):
            return math.pow( math.sin(0.5*math.pi*(t-start)/slope1) , 2 )
        elif t < start+slope1+plateau:
            return 1.
        elif t < start+slope1+plateau+slope2 and (slope2!=0.):
            return math.pow(  math.cos(0.5*math.pi*(t-start-slope1-plateau)/slope2) , 2 )
        else:
            return 0.
    f.profileName = "tsin2plateau"
    f.start       = start
    #f.fwhm        = fwhm
    f.plateau     = plateau
    f.slope1      = slope1
    f.slope2      = slope2
    return f


def transformPolarization(polarization_phi, ellipticity):
    import math
    p = (1.-ellipticity**2)*math.sin(2.*polarization_phi)/2.
    if abs(p) < 1e-10:
        if abs(ellipticity**2-1.)<1e-10: polarization_phi=0.
        dephasing = math.pi/2.
        amplitude = math.sqrt(1./(1.+ellipticity**2))
        amplitudeY = amplitude * (math.cos(polarization_phi)+math.sin(polarization_phi)*ellipticity)
        amplitudeZ = amplitude * (math.sin(polarization_phi)+math.cos(polarization_phi)*ellipticity)
    else:
        dephasing = math.atan(ellipticity/p)
        theta = 0.5 * math.atan( math.tan(2.*polarization_phi) / math.cos(dephasing) )
        while theta<0.: theta += math.pi/2.
        amplitudeY = math.sqrt(2.) * math.cos(theta)
        amplitudeZ = math.sqrt(2.) * math.sin(theta)
    return [dephasing, amplitudeY, amplitudeZ]

def LaserPlanar1D( box_side="xmin", a0=1., omega=1.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(),phaseZero=0.):
    import math
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0
    amplitudeZ *= a0
    # Create Laser
    Laser(
        box_side        = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ amplitudeZ, amplitudeY ],
        phase          = [ dephasing-phaseZero, -phaseZero ],
    )


def LaserGaussian2D( box_side="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=0.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phaseZero=0.):
    import math
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0
    amplitudeZ *= a0
    # Space and phase envelopes
    Zr = omega * waist**2/2.
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
        phaseZero += phase(Y2)
    # Create Laser
    Laser(
        box_side        = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y:amplitudeZ*spatial(y), lambda y:amplitudeY*spatial(y) ],
        phase          = [ lambda y:phase(y)-phaseZero+dephasing, lambda y:phase(y)-phaseZero ],
    )

def LaserGaussian3D( box_side="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=[0.,0.],
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phaseZero=0.):
    import math
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0
    amplitudeZ *= a0
    # Space and phase envelopes
    Zr = omega * waist**2/2.
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
        phaseZero += phase(focus[1]-sz/cz*focus[0], focus[2]+sy/cy/cz*focus[0])
        
    # Create Laser
    Laser(
        box_side        = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y,z:amplitudeZ*spatial(y,z), lambda y,z:amplitudeY*spatial(y,z) ],
        phase          = [ lambda y,z:phase(y,z)-phaseZero+dephasing, lambda y,z:phase(y,z)-phaseZero ],
    )

# Define the tools for the MPI fft and the propagation of a laser profile
_missing_packages = []
try   : import numpy as np
except: _missing_packages += ["numpy"]
try   : from mpi4py import MPI
except: _missing_packages += ["mpi4py"]
try   :
    import h5py
    if not h5py.h5.get_config().mpi:
        _missing_packages += ["h5py built with MPI"]
except:
    _missing_packages += ["h5py"]

if _missing_packages == []:
    _N_LaserOffset = 0
    
    try   : from inspect import getfullargspec as getargspec
    except: from inspect import getargspec
    
    # Gets the number of arguments of a function
    def _nargs(f):
        return len(getargspec(f)[0])
    
    # Gets the smallest multiple of M greater or equal to N
    def _optimalFFTsize(M, N):
        return int(np.ceil(N/float(M))*M)
    
    class MPI_FFT(object):
        """
        2D & 3D FFT using MPI domain decomposition
        The first two axes of the global domain must be divisible by the MPI size
        """
        
        def __init__(self, N, L, comm):
            from numpy.fft import fftfreq
            self.comm = comm
            self.MPI_rank = comm.Get_rank()
            self.MPI_size = comm.Get_size()
            self.N = np.array(N, dtype=np.int   ) # Arrays shape
            self.L = np.array(L, dtype=np.double) # Physical dimensions of the arrays
            assert self.N.ndim == 1
            assert self.L.ndim == 1
            assert self.N.size == self.L.size
            assert self.N.size in [2,3]
            assert self.N[0] % self.MPI_size == 0
            self.Nlocal = self.N // self.MPI_size # Size in the current proc
            
            # Array shape in the 'real space', for the current processor
            self.X_shape = (self.Nlocal[0], self.N[1])
            # Array shape in the 'frequency space', for the current processor
            self.K_shape = (self.N[0], self.Nlocal[1])
            # Slice of the global array in 'real space', for the current processor
            self.X_slice = (
                slice(self.MPI_rank*self.Nlocal[0], (self.MPI_rank+1)*self.Nlocal[0], 1),
                slice(0, self.N[1], 1),
            )
            # Slice of the global array in 'frequency space', for the current processor
            self.K_slice = (
                slice(0, self.N[0], 1),
                slice(self.MPI_rank*self.Nlocal[1], (self.MPI_rank+1)*self.Nlocal[1], 1),
            )
            # Arrays of wavenumbers (kx, ky, ky) owning to the current processor
            self.local_k = (
                fftfreq(self.N[0], self.L[0]/self.N[0]) *2.*np.pi ,
                fftfreq(self.N[1], self.L[1]/self.N[1])[self.K_slice[1]] *2.*np.pi,
            )
            
            if self.N.size == 3:
                assert self.N[1] % self.MPI_size == 0
                self.X_shape += (self.N[2],)
                self.K_shape += (self.N[2],)
                self.X_slice += (slice(0, self.N[2], 1),)
                self.K_slice += (slice(0, self.N[2], 1),)
                self.local_k += (fftfreq(self.N[2], self.L[2]/self.N[2]) *2.*np.pi,)
                self.fft  = self._fft3
                self.ifft = self._ifft3
        
        def fft(self, A):
            from numpy.fft import fft
            # Last direction
            Hslab = fft(A, axis=1)
            # Transpose MPI decomposition
            blocks = np.ascontiguousarray(Hslab.reshape((self.Nlocal[0], self.MPI_size, self.Nlocal[1])).transpose((1,0,2)))
            Vslab = np.ascontiguousarray(np.empty((self.N[0], self.Nlocal[1]), dtype=np.complex))
            self.comm.Alltoall([blocks, MPI.C_DOUBLE_COMPLEX], [Vslab, MPI.C_DOUBLE_COMPLEX])
            # First direction
            return fft(Vslab, axis=0)
        
        def ifft(self, A, last_axis=True):
            from numpy.fft import ifft
            # First direction
            Vslab = np.ascontiguousarray(ifft(A, axis=0))
            if last_axis:
                # Transpose MPI decomposition
                blocks = np.ascontiguousarray(np.empty( (self.MPI_size, self.Nlocal[0], self.Nlocal[1]), dtype=np.complex))
                self.comm.Alltoall([Vslab, MPI.C_DOUBLE_COMPLEX], [blocks, MPI.C_DOUBLE_COMPLEX])
                Hslab = blocks.transpose((1,0,2)).reshape((self.Nlocal[0], self.N[1]))
                # Other direction
                return ifft(Hslab, axis=1)
            else:
                return Vslab
        
        def _fft3(self, A):
            from numpy.fft import fft, fft2
            # Last two directions
            Hslab = fft2(A, axes=(1,2))
            # Transpose MPI decomposition
            blocks = np.ascontiguousarray(Hslab.reshape((self.Nlocal[0], self.MPI_size, self.Nlocal[1], self.N[2])).transpose((1,0,2,3)))
            Vslab = np.ascontiguousarray(np.empty((self.N[0], self.Nlocal[1], self.N[2]), dtype=np.complex))
            self.comm.Alltoall([blocks, MPI.C_DOUBLE_COMPLEX], [Vslab, MPI.C_DOUBLE_COMPLEX])
            # First direction
            return fft(Vslab, axis=0)
        
        def _ifft3(self, A, last_axis=True):
            from numpy.fft import ifft, ifft2
            N2 = A.shape[2]
            # First direction
            Vslab = np.ascontiguousarray(ifft(A, axis=0))
            # Transpose MPI decomposition
            blocks = np.ascontiguousarray(np.empty( (self.MPI_size, self.Nlocal[0], self.Nlocal[1], N2), dtype=np.complex))
            self.comm.Alltoall([Vslab, MPI.C_DOUBLE_COMPLEX], [blocks, MPI.C_DOUBLE_COMPLEX])
            Hslab = blocks.transpose((1,0,2,3)).reshape((self.Nlocal[0], self.N[1], N2))
            # Other directions
            if last_axis:
                return ifft2(Hslab, axes=(1,2))
            else:
                return ifft2(Hslab, axes=(1,))
    
    def LaserOffset(box_side="xmin", B_profiles=[], offset=0.):
        global _N_LaserOffset
        
        profiles = []
        profiles_n = []
        for i,p in enumerate(B_profiles):
            if p:
                profiles.append( p )
                profiles_n.append( str(i+1) )
        if len(profiles) == 0:
            raise Exception("LaserOffset requires at least one profile defined")
        
        # Obtain the size of the array from Main()
        global Main
        if len(Main)==0:
            raise Exception("LaserOffset has been defined before `Main()`")
        N = np.array(Main.number_of_cells) + 8 # includes oversize of 4 cells
        L = np.array(Main.grid_length) + np.array(Main.cell_length)*8
        Nt = int( Main.simulation_time / Main.timestep )
        Lt = Main.simulation_time
        side = {"x":0, "y":1, "z":2}[box_side[0]]
        ndim = N.size
        if ndim not in [2,3]:
            raise Exception("LaserOffset cannot be defined in "+str(ndim)+"D")
        _2D = ndim == 2
        _3D = ndim == 3
        if _2D:
            N = np.array([N[1-side], Nt])
            L = np.array([L[1-side], Lt])
        else:
            N = np.array([N[(side+1)%3], N[(side+2)%3], Nt])
            L = np.array([L[(side+1)%3], L[(side+2)%3], Lt])
        
        # Define the file name
        file = 'LaserOffset'+str(_N_LaserOffset)+'.h5'
        
        # Test the profiles
        for p in profiles:
            if _nargs(p) != ndim:
                raise Exception("LaserOffset requires profiles with "+str(ndim)+" arguments ("+("y,t" if ndim==2 else "y,z,t")+"), found "+str(n))
        # Convert profiles for numpy
        profiles = [np.vectorize(p) for p in profiles]
        
        if not _test_mode:
            
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            
            # Make the array bigger in order to accomodate for the FFT
            N = np.vectorize(_optimalFFTsize)(size, N)
            
            # Define the FFT procedure
            FFT = MPI_FFT(N, L, comm)
            
            # Calculate grid coordinates
            coordinates = [np.linspace(0., L[i], N[i])[FFT.X_slice[i]] for i in range(len(N))]
            mesh = np.meshgrid(*coordinates, indexing="ij")
            
            # Calculate the profiles as arrays
            B = [p(*mesh) for p in profiles]
            
            # Fourier transform of the fields at destination (keep only positive time frequencies)
            B_FT = [FFT.fft(b) for b in B]
            
            # Get lists of ky, kz and omega
            k_coordinates = list(FFT.local_k)
            
            # Select only interesting omegas
            local_spectrum = sum([np.abs(b)**2 for b in B_FT])
            while local_spectrum.ndim>1: local_spectrum = np.sum( local_spectrum, axis=0 )
            if _2D:
                # In 2D, the spectrum is scattered across processors, so we gather
                spectrum = np.empty((size, len(local_spectrum)))
                comm.Allgather( local_spectrum, spectrum )
                spectrum = spectrum.reshape((-1,))
            else:
                # In 3D, each processor has the full spectrum, so we sum all contributions
                spectrum = np.empty_like(local_spectrum)
                comm.Allreduce( local_spectrum, spectrum )
            indices = np.sort(np.argsort( spectrum[:Nt//2] )[-100:]) # select 100 most intense omegas
            N_FT = N.copy()
            N_FT[-1] = len(indices)
            if _2D:
                # Keep only indices in this proc
                indices -= rank*FFT.Nlocal[1]
                this_proc = (0<indices) * (indices<FFT.Nlocal[1])
                try: MPI_omega_offset = np.flatnonzero(this_proc)[0]
                except: MPI_omega_offset = 0
                indices = indices[this_proc]
            k_coordinates[-1] = k_coordinates[-1][indices] # select omegas
            B_FT = [b[..., indices] for b in B_FT]
            
            # Calculate kx (direction of propagation)
            k_mesh = np.meshgrid(*k_coordinates, indexing="ij")
            kx = k_mesh[-1]**2 # omega^2
            for i in range(0,len(k_mesh)-1): kx -= k_mesh[i]**2;
            kx[kx<0.] = 0.
            kx = np.sqrt(kx)
            
            # Backwards propagation
            P = np.exp( 1j * kx * offset )
            B_FT = [b*P for b in B_FT]
            
            # Fourier transform back to real space
            P = np.exp( -1j * k_mesh[-1] * offset ) / Lt # add time offset and divide by omega increment
            B = [ P * FFT.ifft(b, last_axis=False) for b in B_FT ]
            
            # Find the file region where each proc will write
            if _2D:
                region = (slice(None), slice(MPI_omega_offset, MPI_omega_offset+len(indices), 1))
            else:
                region = (FFT.X_slice[0], FFT.X_slice[1], slice(None))
            
            # Store omega in hdf5 file
            if _2D:
                with h5py.File(file, 'w', driver='mpio', comm=comm) as f:
                    dataset = f.create_dataset('omega', (N_FT[-1],), dtype=k_coordinates[-1].dtype)
                    if len(indices)>0:
                        dataset[region[1]] = k_coordinates[-1]
            elif rank == 0: # only rank 0 needed in 3D
                with h5py.File(file, 'w') as f:
                    f.create_dataset('omega', data=k_coordinates[-1])
            comm.barrier()
            
            # Now store the absolute value and the phase, via MPI
            with h5py.File(file, 'r+', driver='mpio', comm=comm) as f:
                i = 0
                for b in B:
                    A = np.abs(b)
                    dataset = f.create_dataset('magnitude'+profiles_n[i], N_FT, dtype=A.dtype)
                    if len(indices)>0:
                        dataset[region] = A
                    A = np.angle(b)
                    dataset = f.create_dataset('phase'+profiles_n[i], N_FT, dtype=A.dtype)
                    if len(indices)>0:
                         dataset[region] = A
                    i += 1
        
        # Create the Laser object
        Laser(
            box_side = "xmin",
            file = file
        )
        
        _N_LaserOffset += 1

else:
    
    def LaserOffset(box_side="xmin", B_profiles=[], offset=0.):
        print("WARNING: LaserOffset unavailable due to missing packages: "+", ".join(_missing_packages))

"""
-----------------------------------------------------------------------
    BEGINNING OF THE USER NAMELIST
"""
