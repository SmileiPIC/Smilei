# Some predefined profiles (see doc)

def constant(value, xvacuum=-float("inf"), yvacuum=-float("inf"), zvacuum=-float("inf")):
    global Main
    if len(Main)==0:
        raise Exception("constant profile has been defined before `Main()`")
    if Main.geometry == "1Dcartesian":
        f = lambda x,t=0: value if x>=xvacuum else 0.
    if (Main.geometry == "2Dcartesian" or Main.geometry == "AMcylindrical"):
        f = lambda x,y,t=0: value if (x>=xvacuum and y>=yvacuum) else 0.
    if Main.geometry == "3Dcartesian":
        f = lambda x,y,z,t=0: value if (x>=xvacuum and y>=yvacuum and z>=zvacuum) else 0.
    f.profileName = "constant"
    f.value   = value
    f.xvacuum = xvacuum
    f.yvacuum = yvacuum
    f.zvacuum = zvacuum
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
    elif (Main.geometry == "2Dcartesian" or Main.geometry == "AMcylindrical"): dim = 2
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
    if (Main.geometry == "2Dcartesian" or Main.geometry == "AMcylindrical"): dim = 2
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
    if (Main.geometry == "2Dcartesian" or Main.geometry == "AMcylindrical"): dim = 2
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
            elif (Main.geometry=="2Dcartesian" or Main.geometry == "AMcylindrical"):
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
    elif (Main.geometry=="2Dcartesian" or Main.geometry == "AMcylindrical"):
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
    if fwhm     is None: fwhm     = duration/3.
    if center   is None: center   = start + duration/2.
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
    from math import pi, sqrt, sin, cos, tan, atan2
    e2 = ellipticity**2
    p = (1.-e2)*sin(2.*polarization_phi)/2.
    dephasing = atan2(ellipticity, p)
    amplitude = sqrt(1./(1.+e2))
    c2 = cos(polarization_phi)**2
    s2 = 1. - c2
    amplitudeY = amplitude * sqrt(c2 + e2*s2)
    amplitudeZ = amplitude * sqrt(s2 + e2*c2)
    return [dephasing, amplitudeY, amplitudeZ]

def LaserPlanar1D( box_side="xmin", a0=1., omega=1.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(),phase_offset=0.):
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
    # Create Laser
    Laser(
        box_side        = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ amplitudeZ, amplitudeY ],
        phase          = [ dephasing-phase_offset, -phase_offset ],
        delay_phase    = [ 0., dephasing ]
    )

def LaserEnvelopePlanar1D( a0=1., omega=1., time_envelope=tconstant(),
        envelope_solver = "explicit",box_side = "inside",Envelope_boundary_conditions = [["reflective"]],
        polarization_phi = 0.,ellipticity = 0.):
    from numpy import vectorize, sqrt

    def spatial_envelope1D(x):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        return (a0*polarization_amplitude_factor)
    
    if (box_side=="inside"):
        def envelope_profile(x,t):
            return spatial_envelope1D(x)*complex( vectorize(time_envelope)(t) )
    elif (box_side=="xmin"):
        def envelope_profile(t):
            return spatial_envelope1D(0)*complex( vectorize(time_envelope)(t) )
    else:
        print("LaserEnvelope error: box_side must be either 'inside' or 'xmin'. ")

    # Create Laser Envelope
    LaserEnvelope(
        omega                        = omega,
        envelope_profile             = envelope_profile,
        envelope_solver              = envelope_solver,
        box_side                     = box_side,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi             = polarization_phi,
        ellipticity                  = ellipticity
    )


def LaserGaussian2D( box_side="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=0.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phase_offset=0.):
    from math import pi, cos, sin, tan, atan, sqrt, exp
    assert len(focus)==2, "LaserGaussian2D: focus must be a list of length 2."
    global Main
    assert len(Main)==1, "LaserGaussian2D profile has been defined before `Main()`"
    grid_length = Main.grid_length
    # Polarization and amplitude
    dephasing, amplitudeZ, amplitudeY = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
    delay_phase = [0., dephasing]
    # Injection on ymin/ymax
    if box_side[0] == "y":
        focus = focus[::-1]
        grid_length = grid_length[::-1]
        amplitudeY = -amplitudeY
    # Injection on max boundary
    if box_side.endswith("max"):
        focus[0] = grid_length[0] - focus[0]
    # Space and phase envelopes
    Zr = omega * waist**2/2.
    if incidence_angle == 0.:
        w  = sqrt(1./(1.+(focus[0]/Zr)**2))
        invWaist2 = (w/waist)**2
        coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
        def spatial(y):
            return sqrt(w) * exp( -invWaist2*(y-focus[1])**2 )
        def phase(y):
            return coeff * (y-focus[1])**2
    else:
        invZr  = sin(incidence_angle) / Zr
        invZr2 = invZr**2
        invZr3 = (cos(incidence_angle) / Zr)**2 / 2.
        invWaist2 = (cos(incidence_angle) / waist)**2
        omega_ = omega * sin(incidence_angle)
        Y1 = focus[1] + focus[0]/tan(incidence_angle)
        Y2 = focus[1] - focus[0]*tan(incidence_angle)
        amplitudeY *= cos(incidence_angle)
        def spatial(y):
            w2 = 1./(1. + invZr2*(y-Y1)**2)
            return sqrt(sqrt(w2)) * exp( -invWaist2*w2*(y-Y2)**2 )
        def phase(y):
            dy = y-Y1
            return omega_*dy*(1.+ invZr3*(y-Y2)**2/(1.+invZr2*dy**2)) - 0.5*atan(invZr*dy)
        # Adjust the phase to match that of a laser that could come from another face
        if Y2 < 0:
            distance_to_boundary = focus[1] / sin(incidence_angle)
        elif Y2 < grid_length[1]:
            distance_to_boundary = focus[0] / cos(incidence_angle)
        else:
            distance_to_boundary = (focus[1] - grid_length[1]) / sin(incidence_angle)
        phase_offset -= omega * distance_to_boundary - 0.5*atan(distance_to_boundary/Zr)
    # Create Laser
    Laser(
        box_side       = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y:amplitudeY*spatial(y), lambda y:amplitudeZ*spatial(y) ],
        phase          = [ lambda y:phase(y)-phase_offset+delay_phase[1], lambda y:phase(y)-phase_offset+delay_phase[0] ],
        delay_phase    = delay_phase
    )

def LaserSquareFlattenedGaussian2D( box_side="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=0.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phase_offset=0.,N=10,flattened_intensity_position="far_from_focus"):
    import numpy as np
    import scipy.special as sp
    assert len(focus)==2, "LaserSquareFlattenedGaussian2D: focus must be a list of length 2."
    assert incidence_angle == 0, "LaserSquareFlattenedGaussian2D: currently only incidence_angle=0 is supported."
    assert box_side == "xmin", "LaserSquareFlattenedGaussian2D: currently only box_side=xmin is supported."
    assert isinstance(N, int), "LaserSquareFlattenedGaussian2D: N must be an integer."
    assert (
            flattened_intensity_position == "far_from_focus"
            or flattened_intensity_position == "at_focus"
           ), (
           "LaserSquareFlattenedGaussian2D: flattened_intensity_position must be either 'at_focus' or 'far_from_focus'."
           )

    # Polarization and amplitude
    dephasing, amplitudeZ, amplitudeY = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
    delay_phase = [0., dephasing]

    # Waist and Rayleigh length
    waist_corrected = (
                       waist * np.sqrt(N + 1) if flattened_intensity_position == "far_from_focus"
                       else waist / np.sqrt(N + 1)
                       )
    # Effective Rayleigh length, normalized
    x_R  = omega * waist_corrected**2/2.

    # Store the Hermite-Gauss (HG) mode coefficients
    cn = np.zeros(N+1)
    for n in range(N+1):
        m_values = np.arange(n, N+1)
        # computing this in log scale and then using the exponential
        # avoids overflow for high N
        # remember that gamma(n+1)=n!
        log_terms = (
                    -3*m_values*np.log(2.)
                    + sp.gammaln(2*m_values+1)
                    - sp.gammaln(m_values+1)
                    - sp.gammaln(m_values-n+1))
        cn[n] = np.sum(np.exp(log_terms - np.max(log_terms))) \
                    * np.exp(np.max(log_terms) - sp.gammaln(2*n+1))

    # The HG mode coefficients have alternating signs
    # if the flattened profile is far from focus
    if flattened_intensity_position == "far_from_focus":
        cn = cn*(-1.)**np.arange(N+1)

    # Normalizing quantity to have a normalized peak field before the multiplication by a0
    # This recursive way of computing it avoids the overflow given by a brute force calculation
    # of the factorial
    def S_N(N):
            S = 0.0
            k = 1.0  # k_0 = 1
            for m in range(0, N+1):
                    S += k
                    # update a -> a_{m+1}
                    k *= (2*m + 1) / (2*(m + 1))
            return S

    normalization_constant = 1. if (flattened_intensity_position=="at_focus") else S_N(N)

    # Store Hermite polynomials
    def even_hermite_polynomials(x, N):
        # Returns an array of Hermite polynomials with even index using recursion relations
        H_even = np.empty((N+1,) + np.shape(x), dtype=float)
        H0 = np.ones_like(x)
        H_even[0] = H0
        if N == 0:
            return H_even
        H1   = 2*x
        Hnm2 = H0
        Hnm1 = H1
        even = 1

        for k in range(1, 2*N):
            Hn = 2*x*Hnm1 - 2*k*Hnm2
            Hnm2 = Hnm1
            Hnm1 = Hn
            if (k+1) % 2 == 0:
                H_even[even] = Hn
                even += 1

        return H_even

    # Compute constant terms at x=0
    x              = 0.
    # HG mode waist at x
    w              = waist_corrected * np.sqrt(1 + ((x-focus[0]) /x_R)**2)
    # HG mode Gouy phase argument at x, to multiply by the the mode factor
    Gouy_phase_arg = np.arctan2((x-focus[0]),x_R)
    # HG mode curved wavefront at x
    one_ov_R       = (x - focus[0]) / ((x - focus[0])**2 + x_R**2)

    # Square flattened Gauss definition in 2D Cartesian geometry
    # Compared to a 3D definition, the amplitude decreases as 1/sqrt(w(x)/w0) instead of 1/(w(x)/w0)
    # and the Gouy phase is different
    def rectangular_flattened_Gaussian_beam2D(y):
        y              = np.asarray(y)
        curved_phase_y = np.exp(1j * omega * (y-focus[1])**2 * one_ov_R / 2 )
        # HG mode exponential decay
        exp_along_y    = np.exp(-(y-focus[1])**2 / w**2)
        # Precompute the Hermite polynomials
        mask           = exp_along_y > 0.
        y_scaled       = np.sqrt(2) * (y-focus[1]) / w
        H = np.zeros((N+1,) + y.shape, dtype=float)
        if np.any(mask):
            # Only evaluate Hermite polynomials where the Gaussian has not underflown to zero
            H[:, mask] = even_hermite_polynomials(y_scaled[mask], N)
            # Convert the nan to zero, that can happen only when the exponential is ~0
            H          = np.nan_to_num(H,nan=0.0,posinf=0.0,neginf=0.0)
        # Sum the HG modes parts that change for each mode
        HG_field_along_y = np.zeros_like(y,dtype=complex)
        for n in range(0, N+1):
            HG_field_along_y += cn[n] * H[n] * np.exp(-1j*(2*n+1/2.)*Gouy_phase_arg)
        # Multiply the result by the part in common for all modes
        HG_field_along_y = HG_field_along_y * exp_along_y * curved_phase_y * np.sqrt(waist_corrected/w)

        return HG_field_along_y/normalization_constant

    # define the Laser at x=0 through the space_time profile of By and Bz
    def complex_envelope(y,t):
        return rectangular_flattened_Gaussian_beam2D(y)*time_envelope(t)*np.exp(-1j*phase_offset)
    def By_profile(y,t):
        return amplitudeY*np.real(complex_envelope(y,t)*np.exp(1j*delay_phase[1])*np.exp(-1j*omega*t))
    def Bz_profile(y,t):
        return amplitudeZ*np.real(complex_envelope(y,t)*np.exp(1j*delay_phase[0])*np.exp(-1j*omega*t))

    # Create Laser
    Laser(
        box_side = "xmin",
        space_time_profile = [ By_profile, Bz_profile ]
    )

def rotation(x,y,ang) :
    '''
    Lineare tranformation: Rotation matrix
    (x,y)->(x',y')
    '''
    from math import cos, sin
    xrot = +cos(ang)*x + sin(ang)*y
    yrot = -sin(ang)*x + cos(ang)*y
    return xrot,yrot

def transform(x,y,xf,yf,L,ang) :
    '''
    Function to transform coordinate of laser-RPP formula
    x,y : Lab/Simulation box coordinate
    X,Y : Coordinate where X is the propagation axis of the laser, Y is transvers axis
    X,Y are rotated and translated coordinate in order the user define the 'focal spot' xf,yf in box coordinate and the angle of incidence with respect of x-axis of the simulated box
    '''
    from math import cos,sin,tan
    X,Y = rotation(x,y,ang)
    X = X+(L-xf/cos(ang))-(yf-tan(ang)*xf)*sin(ang)
    Y = Y-(yf-tan(ang)*xf)*cos(ang)
    return X,Y

def LaserSmoothing2D(box_side="xmin", a0=1., omega=1., focus=None, incidence_angle=0.,polarization_phi=0.,ellipticity=0.,phase_zero=0.,
               Lf=3.00e6,fnumber=8.00,
               N=6,rpp_random_seed=42,
               temporal_smoothing=None,
               omega_m=0.,modulation_depth=0,spectral_profile=lambda w: 1.,frequency_comb=False,mode_locking=False,
               temporal_freq_random_seed=1789,temporal_phi_random_seed=1793,
               rpp_per_mode=False,rpp_seed_per_mode=[42],
               omega_m_trans=0.,modulation_depth_trans=0,mode2generate_trans=None,chirp_profile=tconstant(),
               omega_m_longi=0.,modulation_depth_longi=0,mode2generate_longi=None,
               space_envelope=lambda y:1.,time_envelope=tconstant()):
    '''
    Default values are in code units
    incidence_angle in radian
    a0                     : Maximum of the envelope at focal spot for 1 speckle (i.e. N=40 and no random phase between element, or N=1). Otherwise, for N=40, a = a0/sqrt(N=40) in the simulation box.
    Lf                     : Longueur focale without SSD
    fnumber                : F-number
    N                      : Number of phase plate element
    rpp_random_seed        : Seed in order to have a Random Phase Plate (None is = no random, all element have zero phase-shift),
    temporal_smoothing     : None/'Broadband'/'TSSD'/'LSSD'
    omega_m                : modulation frequency for Broadband Laser
    modulation_depth       : depth 'm' of modulation and frequency bandwith = 2m for Broadband Laser
    frequency_comb         : False/True mode are equally spaced (modulation time) or not (no modulation time)
    mode_locking           : False/True mode have randomly chosen phase or not
    rpp_per_mode           : False/'Stardriver'/'ISI' : Same RPP for each mode / Change the RPP for each mode / Phase delay per mode
    rpp_seed_per_mode      : Seed for RRP or [parameter] for ISI echelon
    omega_m_trans          : modulation frequency for transverse TSSD
    omega_m_longi          : modulation frequency for longitudinal LSSD
    modulation_depth_trans : depth 'm' of modulation and frequency bandwith = 2m for transverse SSD
    modulation_depth_longi : depth 'm' of modulation and frequency bandwith = 2m for longitudinal SSD
    '''
    import numpy as np
    from math import pi, sqrt, cos, sin, tan, fabs
    from cmath import exp,rect,polar
    from scipy.special import erf,jv
    
    global Main

    if temporal_smoothing==None:
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='Broadband':
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='TSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='LSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0

    if len(Main)==0:
        raise Exception("LaserSmoothing2D profile has been defined before `Main()`") 
        
    k0 = omega
    waist = fnumber*N*(2.00*pi/omega)
    D = Lf/fnumber #taille lame de phase
    d = D/N #taille element lame de phase
    R = sqrt(waist*d/pi)

    x_focus,y_focus = focus[0],focus[1]

    El = (a0*omega/N)/np.sqrt(k0*d**2/(2*Lf*np.pi**2))
    # Polarization and amplitude
    dephasing, amplitudeZ, amplitudeY = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= El * cos(incidence_angle)
    amplitudeZ *= El
    delay_phase = [0., dephasing]

    krpp = np.linspace(-D/2,D/2,N+1)
    phik = np.zeros(N)
    if rpp_random_seed != None :
        np.random.seed(rpp_random_seed) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
        phikinit = np.random.rand(N)
        for i in range(0,N):
            phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi

    m_broadband = modulation_depth
    m_trans = modulation_depth_trans
    m_longi = modulation_depth_longi
    alpha_t = 2*pi/D
    alpha_x = 1/omega

    modes_trans = range(-m_trans,m_trans+1,1)
    modes_longi = range(-m_longi,m_longi+1,1)
    modes_broadband = range(-m_broadband,m_broadband+1,1)
    
    if temporal_smoothing=='Broadband':
        if frequency_comb :
            mode_offsets = np.zeros(len(modes_broadband))
        else :
            np.random.seed(temporal_freq_random_seed)
            mode_offsets = np.random.uniform(-0.5,0.5,len(modes_broadband))
        if mode_locking :
            phase_w = np.zeros(len(modes_broadband)) 
        else :
            np.random.seed(temporal_phi_random_seed)
            phase_w = 2*pi*np.random.rand(len(modes_broadband))

        # Legacy
        # Ebb = np.ones(len(modes_broadband))/np.sqrt(len(modes_broadband))

        if not callable(spectral_profile):
            raise Exception("spectral_profile must be a callable function of omega, e.g. lambda w: 1.")
        
        # Example 1 for intensity spectrum : Gaussian
        # sigma = 10 * omega_m
        # spectral_profile = lambda w: np.exp(-(w - omega)**2 / (2*sigma**2))
        # Example 2 for intensity spectrum : Lorentzian
        # gamma = 5 * omega_m
        # spectral_profile = lambda w: gamma / ((w - omega)**2 + gamma**2)
        # Default : Flat
        # spectral_profile = lambda w: 1.
        
        Ebb = np.array([spectral_profile(omega * (1. + (mode + mode_offsets[int(mode+m_broadband)]) * omega_m / omega)) for mode in modes_broadband])
        # Normalization : sum of all spectral intensity line = 1
        Ebb = Ebb / np.sqrt(np.sum(Ebb**2))
            

    def ERPP(y,imode,imode_t,imode_l,phik) :
        '''
        Formula (24) of "Cross-beam energy transfer between spatially smoothed laser beams" [A. Oudin, A. Debayle, C. Ruyer]
        Spatial envelope definition at x=0 of the simulation domain.
        For a fixed SSD mode. SSD is treated as a Laser() superposition at different frequency.
        '''

        X,Y = transform(0,y,x_focus,y_focus,Lf,incidence_angle)

        Lfw = Lf*(1+alpha_x*imode_l*omega_m_longi)
        # K   = sqrt( fabs( k0/(2*X) - 1/R**2 ) )
        K = sqrt( fabs( k0*(Lfw-X)/(2*X*Lfw) ) )
        sum_erfm = 0+0*1j
        factm    = 0+0*1j

        if (temporal_smoothing == None) & ((imode_t != 0) | (imode_l != 0)):
            raise Exception("Input inconsistency : No temporal smoothing selected but non-zero 'modulation depth'")

        # Valid for X<L
        if X<Lfw :
            factm = (1/np.sqrt(pi))*0.5*sqrt(Lfw/(Lfw-X))*exp( -1j*k0*Y**2/(2*(Lfw-X)) + 1j*(imode_t*alpha_t*Y-(imode_t*alpha_t)**2*X/(2*k0))*Lfw/(Lfw-X)) # Phase RPP+TSSD
            for n in range (0,N):
                  # sum_erfm += (erf(exp(-1j*pi/4)*K*(krpp[n+1]-(Y-imode_t*alpha_t*X/k0)*k0/2/X/K**2)) - erf(exp(-1j*pi/4)*K*(krpp[n]-(Y-imode_t*alpha_t*X/k0)*k0/2/X/K**2)))*exp(1j*phik[n])
                  sum_erfm += (erf(exp(-1j*pi/4)*K*(krpp[n+1]-(Y-imode_t*alpha_t*X/k0)*Lfw/(Lfw-X))) - erf(exp(-1j*pi/4)*K*(krpp[n]-(Y-imode_t*alpha_t*X/k0)*Lfw/(Lfw-X))))*exp(1j*phik[n])
        # At focus
        elif X==Lfw :
           factm = exp(-1j*pi*0.25)*sqrt(k0*d*d/(2*pi*pi*Lfw))*exp(1j*k0*Y**2/(2*Lfw))*sin(k0*d/(2*Lfw)*(Y-imode_t*alpha_t*X/k0))/(k0*d/(2*Lfw)*(Y-imode_t*alpha_t*X/k0))
           for n in range (0,N):
               sum_erfm += exp(1j*phik[n]-1j*k0/Lfw*(Y-imode_t*alpha_t*X/k0)*(krpp[n+1]+krpp[n])/2)
        # Beyond focus
        else :
            factm = (1/np.sqrt(pi))*0.5*exp(-1j*pi*0.5)*sqrt(Lfw/fabs(Lfw-X))*exp(-1j*k0*Y**2/(2*(Lfw-X)) + 1j*(imode_t*alpha_t*Y-(imode_t*alpha_t)**2*X/(2*k0))*Lfw/(Lfw-X)) # Phase RPP+TSSD
            for n in range (0,N):
                # sum_erfm += (erf(exp(+1j*pi/4)*fabs(K)*(krpp[n+1]-k0*(Y-imode_t*alpha_t*X/k0)/2/X/K**2)) - erf(exp(+1j*pi/4)*fabs(K)*(krpp[n]-k0*(Y-imode_t*alpha_t*X/k0)/2/X/K**2)))*exp(1j*phik[n]) # Somme des erf
                sum_erfm += (erf(exp(+1j*pi/4)*fabs(K)*(krpp[n+1]-(Y-imode_t*alpha_t*X/k0)*Lfw/(Lfw-X))) - erf(exp(+1j*pi/4)*fabs(K)*(krpp[n]-(Y-imode_t*alpha_t*X/k0)*Lfw/(Lfw-X))))*exp(1j*phik[n])

        if temporal_smoothing==None:
            Einit = sum_erfm*factm*exp(1j*k0*X)#*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_t*omega_m_trans/k0)
        elif (temporal_smoothing=='TSSD') | (temporal_smoothing=='LSSD'):
            Einit = jv(imode_t,m_trans)*jv(imode_l,m_longi)*sum_erfm*factm*exp(1j*k0*X)*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_t*omega_m_trans/k0)
        elif temporal_smoothing=='Broadband':
            Einit = Ebb[imode]*sum_erfm*factm*exp(1j*k0*X)*exp(1j*phase_w[imode])*exp(1j*k0*X*imode*omega_m/k0)
        else :
            raise Exception("Temporal_smoothing method not implemented yet")

        Amp,Phase = polar(Einit)
        return Amp,Phase

    def ERPP_ampBz(y,imode_,imode_t_,imode_l_,phik_) :
        return amplitudeZ*space_envelope(y)*ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_ampBy(y,imode_,imode_t_,imode_l_,phik_) :
        return amplitudeY*space_envelope(y)*ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_phaseBz(y,imode_,imode_t_,imode_l_,phik_) :
        return ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[0]
    def ERPP_phaseBy(y,imode_,imode_t_,imode_l_,phik_) :
        return ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[1]

    # Encapsulated function to avoid warning
    def _mk_amp_By(im, it, il, pk):
        return lambda y: ERPP_ampBy(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)
    def _mk_amp_Bz(im, it, il, pk):
        return lambda y: ERPP_ampBz(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)
    def _mk_phase_By(im, it, il, pk):
        return lambda y: ERPP_phaseBy(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)
    def _mk_phase_Bz(im, it, il, pk):
        return lambda y: ERPP_phaseBz(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)

    fct_amp_By = []
    fct_amp_Bz = []
    fct_phase_By = []
    fct_phase_Bz = []

    if temporal_smoothing=='Broadband':
        phik_base = phik.copy()
        for mode in modes_broadband :
            if rpp_per_mode=='Stardriver':
                if len(rpp_seed_per_mode)!=len(modes_broadband):
                    raise Exception("len(rpp_seed_per_mode): "+str(len(rpp_seed_per_mode))+". len(modes_broadband): "+str(len(modes_broadband))+". Length of rpp_seed_per_mode have to be equal to 2 x modulation_depth + 1 ")
                phik = np.zeros(N)
                np.random.seed(rpp_seed_per_mode[mode+m_broadband]) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
                phikinit = np.random.rand(N)
                for i in range(0,N):
                    phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi
            elif rpp_per_mode=='ISI':
                phik_mode = phik_base.copy()
                delta_t_echelon = rpp_seed_per_mode[0] * 2*pi / (2 * m_broadband * omega_m)
                omega_k = omega * (1. + (mode+mode_offsets[mode+m_broadband]) * omega_m / omega)
                for i in range(0,N):
                    phik_mode[i] += omega_k * i * delta_t_echelon
                phik_mode = phik_mode % (2*pi)
                phik = phik_mode
            else :
                phik = 1*phik
            #fct_amp_By.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_amp_Bz.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_By.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_Bz.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            fct_amp_By.append(_mk_amp_By(mode, 0, 0, phik))
            fct_amp_Bz.append(_mk_amp_Bz(mode, 0, 0, phik))
            fct_phase_By.append(_mk_phase_By(mode, 0, 0, phik))
            fct_phase_Bz.append(_mk_phase_Bz(mode, 0, 0, phik))
        fct_amp_By = np.array(fct_amp_By)
        fct_amp_Bz = np.array(fct_amp_Bz)
        fct_phase_By = np.array(fct_phase_By)
        fct_phase_Bz = np.array(fct_phase_Bz)
    else :
        for mode_t in modes_trans :
            for mode_l in modes_longi :
                #fct_amp_By.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                #fct_amp_Bz.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                #fct_phase_By.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                #fct_phase_Bz.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                fct_amp_By.append(_mk_amp_By(0, mode_t, mode_l, phik))
                fct_amp_Bz.append(_mk_amp_Bz(0, mode_t, mode_l, phik))
                fct_phase_By.append(_mk_phase_By(0, mode_t, mode_l, phik))
                fct_phase_Bz.append(_mk_phase_Bz(0, mode_t, mode_l, phik))
        fct_amp_By = np.reshape(np.array(fct_amp_By),(2*m_trans+1,2*m_longi+1))
        fct_amp_Bz = np.reshape(np.array(fct_amp_Bz),(2*m_trans+1,2*m_longi+1))
        fct_phase_By = np.reshape(np.array(fct_phase_By),(2*m_trans+1,2*m_longi+1))
        fct_phase_Bz = np.reshape(np.array(fct_phase_Bz),(2*m_trans+1,2*m_longi+1))

    if temporal_smoothing=='Broadband':
        for mode in modes_broadband :
            im = int(mode+m_broadband)
            Laser(
                box_side       = box_side,
                omega          = omega*(1.+(mode+mode_offsets[im])*omega_m/omega),
                # omega          = omega,
                # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                time_envelope  = time_envelope,
                space_envelope = [fct_amp_By[im],fct_amp_Bz[im]],
                phase          = [fct_phase_By[im],fct_phase_Bz[im]],
                delay_phase    = delay_phase
            )
    else :
        if mode2generate_trans != None :
            mode_t = 1.*mode2generate_trans
            if mode2generate_longi != None :
                mode_l = 1.*mode2generate_longi
                im_t = int(mode_t+m_trans)
                im_l = int(mode_l+m_longi)
                Laser(
                    box_side       = box_side,
                    omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                    # omega          = omega,
                    # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                    time_envelope  = time_envelope,
                    space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                    phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                    delay_phase    = delay_phase
                )
            else :
                for mode_l in modes_longi :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
        else :
            if mode2generate_longi != None :
                mode_l = mode2generate_longi
                for mode_t in modes_trans :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
            else :
                for mode_t in modes_trans :
                    for mode_l in modes_longi :
                        im_t = int(mode_t+m_trans)
                        im_l = int(mode_l+m_longi)
                        Laser(
                            box_side       = box_side,
                            omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                            # omega          = omega,
                            # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                            time_envelope  = time_envelope,
                            space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                            phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                            delay_phase    = delay_phase
                        )

def LaserSmoothingPeriodic2D(box_side="xmin", a0=1., omega=1., focus=None, incidence_angle=0.,polarization_phi=0.,ellipticity=0.,phase_zero=0.,
               Lf=3.00e6,fnumber=8.00,
               N=6,rpp_random_seed=42,
               temporal_smoothing=None,
               omega_m=0.,modulation_depth=0,spectral_profile=lambda w: 1.,frequency_comb=False,mode_locking=False,
               temporal_freq_random_seed=1789,temporal_phi_random_seed=1793,
               rpp_per_mode=False,rpp_seed_per_mode=[42],
               omega_m_trans=0.,modulation_depth_trans=0,mode2generate_trans=None,chirp_profile=tconstant(),
               omega_m_longi=0.,modulation_depth_longi=0,mode2generate_longi=None,
               space_envelope=lambda y:1.,time_envelope=tconstant()):
    '''
    Default values are in code units
    incidence_angle in radian
    a0                     : Maximum of the envelope at focal spot for 1 speckle (i.e. N=40 and no random phase between element, or N=1). Otherwise, for N=40, a = a0/sqrt(N=40) in the simulation box.
    Lf                     : Longueur focale without SSD
    fnumber                : F-number
    N                      : Number of phase plate element
    rpp_random_seed        : Seed in order to have a Random Phase Plate (None is = no random, all element have zero phase-shift),
    temporal_smoothing     : None/'Broadband'/'TSSD'/'LSSD'
    omega_m                : modulation frequency for Broadband Laser
    modulation_depth       : depth 'm' of modulation and frequency bandwith = 2m for Broadband Laser
    frequency_comb         : False/True mode are equally spaced (modulation time) or not (no modulation time)
    mode_locking           : False/True mode have randomly chosen phase or not
    rpp_per_mode           : False/'Stardriver'/'ISI' : Same RPP for each mode / Change the RPP for each mode / Phase delay per mode
    rpp_seed_per_mode      : Seed for RRP or [parameter] for ISI echelon
    omega_m_trans          : modulation frequency for transverse TSSD
    omega_m_longi          : modulation frequency for longitudinal LSSD
    modulation_depth_trans : depth 'm' of modulation and frequency bandwith = 2m for transverse SSD
    modulation_depth_longi : depth 'm' of modulation and frequency bandwith = 2m for longitudinal SSD
    '''
    import numpy as np
    from math import pi, sqrt, cos, sin, tan, fabs
    from cmath import exp,rect,polar
    from scipy.special import erf,jv
    
    global Main

    if temporal_smoothing==None:
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='Broadband':
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='TSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='LSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0

    if len(Main)==0:
        raise Exception("LaserSmoothingPeriodic2D profile has been defined before `Main()`") 
        
    k0 = omega
    waist = fnumber*N*(2.00*pi/omega)
    D = Lf/fnumber #taille lame de phase
    d = D/N #taille element lame de phase
    R = sqrt(waist*d/pi)

    x_focus,y_focus = focus[0],focus[1]

    El = (a0*omega/N)/np.sqrt(k0*d**2/(2*Lf*np.pi**2))
    # Polarization and amplitude
    dephasing, amplitudeZ, amplitudeY = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= El * cos(incidence_angle)
    amplitudeZ *= El
    delay_phase = [0., dephasing]

    krpp = np.linspace(-D/2,D/2,N+1)
    phik = np.zeros(N)
    if rpp_random_seed != None :
        np.random.seed(rpp_random_seed) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
        phikinit = np.random.rand(N)
        for i in range(0,N):
            phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi

    m_broadband = modulation_depth
    m_trans = modulation_depth_trans
    m_longi = modulation_depth_longi
    alpha_t = 2*pi/D
    alpha_x = 1/omega

    modes_trans = range(-m_trans,m_trans+1,1)
    modes_longi = range(-m_longi,m_longi+1,1)
    modes_broadband = range(-m_broadband,m_broadband+1,1)
    
    if temporal_smoothing=='Broadband':
        if frequency_comb :
            mode_offsets = np.zeros(len(modes_broadband))
        else :
            np.random.seed(temporal_freq_random_seed)
            mode_offsets = np.random.uniform(-0.5,0.5,len(modes_broadband))
        if mode_locking :
            phase_w = np.zeros(len(modes_broadband)) 
        else :
            np.random.seed(temporal_phi_random_seed)
            phase_w = 2*pi*np.random.rand(len(modes_broadband))
        
        # Legacy
        # Ebb = np.ones(len(modes_broadband))/np.sqrt(len(modes_broadband))

        if not callable(spectral_profile):
            raise Exception("spectral_profile must be a callable function of omega, e.g. lambda w: 1.")

        # Example 1 for intensity spectrum : Gaussian
        # sigma = 10 * omega_m
        # spectral_profile = lambda w: np.exp(-(w - omega)**2 / (2*sigma**2))
        # Example 2 for intensity spectrum : Lorentzian
        # gamma = 5 * omega_m
        # spectral_profile = lambda w: gamma / ((w - omega)**2 + gamma**2)
        # Default : Flat
        # spectral_profile = lambda w: 1.

        Ebb = np.array([spectral_profile(omega * (1. + (mode + mode_offsets[int(mode+m_broadband)]) * omega_m / omega)) for mode in modes_broadband])
        # Normalization : sum of all spectral intensity line = 1
        Ebb = Ebb / np.sqrt(np.sum(Ebb**2))


    def ERPP(y,imode,imode_t,imode_l,phik) :
        '''
        Formula (24) of "Cross-beam energy transfer between spatially smoothed laser beams" [A. Oudin, A. Debayle, C. Ruyer]
        Spatial envelope definition at x=0 of the simulation domain.
        For a fixed SSD mode. SSD is treated as a Laser() superposition at different frequency.
        '''

        X,Y = transform(0,y,x_focus,y_focus,Lf,incidence_angle)

        Lfw = Lf*(1+alpha_x*imode_l*omega_m_longi)
        # K   = sqrt( fabs( k0/(2*X) - 1/R**2 ) )
        K = sqrt( fabs( k0*(Lfw-X)/(2*X*Lfw) ) )
        sum_erfm = 0+0*1j
        factm    = 0+0*1j

        if (temporal_smoothing == None) & ((imode_t != 0) | (imode_l != 0)):
            raise Exception("Input inconsistency : No temporal smoothing selected but non-zero 'modulation depth'")

        factm = exp(-1j*pi*0.25)*sqrt(k0*d*d/(2*pi*pi*Lfw))
        for n in range (0,N):
            sum_erfm += exp(1j*phik[n]-1j*k0/Lfw*(Y-imode_t*alpha_t*X/k0)*(krpp[n+1]+krpp[n])/2)

        if temporal_smoothing==None:
            Einit = sum_erfm*factm*exp(1j*k0*X)#*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_t*omega_m_trans/k0)
        elif (temporal_smoothing=='TSSD') | (temporal_smoothing=='LSSD'):
            Einit = jv(imode_t,m_trans)*jv(imode_l,m_longi)*sum_erfm*factm*exp(1j*k0*X)*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_t*omega_m_trans/k0)
        elif temporal_smoothing=='Broadband':
            Einit = Ebb[imode]*sum_erfm*factm*exp(1j*k0*X)*exp(1j*phase_w[imode])*exp(1j*k0*X*imode*omega_m/k0)
        else :
            raise Exception("Temporal_smoothing method not implemented yet")

        Amp,Phase = polar(Einit)
        return Amp,Phase

    def ERPP_ampBz(y,imode_,imode_t_,imode_l_,phik_) :
        return amplitudeZ*space_envelope(y)*ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_ampBy(y,imode_,imode_t_,imode_l_,phik_) :
        return amplitudeY*space_envelope(y)*ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_phaseBz(y,imode_,imode_t_,imode_l_,phik_) :
        return ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[0]
    def ERPP_phaseBy(y,imode_,imode_t_,imode_l_,phik_) :
        return ERPP(y,imode=imode_,imode_t=imode_t_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[1]

    # Encapsulated function to avoid warning
    def _mk_amp_By(im, it, il, pk):
        return lambda y: ERPP_ampBy(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)
    def _mk_amp_Bz(im, it, il, pk):
        return lambda y: ERPP_ampBz(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)
    def _mk_phase_By(im, it, il, pk):
        return lambda y: ERPP_phaseBy(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)
    def _mk_phase_Bz(im, it, il, pk):
        return lambda y: ERPP_phaseBz(y, imode_=im, imode_t_=it, imode_l_=il, phik_=pk)

    fct_amp_By = []
    fct_amp_Bz = []
    fct_phase_By = []
    fct_phase_Bz = []

    if temporal_smoothing=='Broadband':
        phik_base = phik.copy()
        for mode in modes_broadband :
            if rpp_per_mode=='Stardriver':
                if len(rpp_seed_per_mode)!=len(modes_broadband):
                    raise Exception("len(rpp_seed_per_mode): "+str(len(rpp_seed_per_mode))+". len(modes_broadband): "+str(len(modes_broadband))+". Length of rpp_seed_per_mode have to be equal to 2 x modulation_depth + 1 ")
                phik = np.zeros(N)
                np.random.seed(rpp_seed_per_mode[mode+m_broadband]) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
                phikinit = np.random.rand(N)
                for i in range(0,N):
                    phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi
            elif rpp_per_mode=='ISI':
                phik_mode = phik_base.copy()
                delta_t_echelon = rpp_seed_per_mode[0] * 2*pi / (2 * m_broadband * omega_m)
                omega_k = omega * (1. + (mode+mode_offsets[mode+m_broadband]) * omega_m / omega)
                for i in range(0,N):
                    phik_mode[i] += omega_k * i * delta_t_echelon
                phik_mode = phik_mode % (2*pi)
                phik = phik_mode
            else :
                phik = 1*phik
            #fct_amp_By.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_amp_Bz.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_By.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_Bz.append(lambda y,imode_tmp=mode,imode_t_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            fct_amp_By.append(_mk_amp_By(mode, 0, 0, phik))
            fct_amp_Bz.append(_mk_amp_Bz(mode, 0, 0, phik))
            fct_phase_By.append(_mk_phase_By(mode, 0, 0, phik))
            fct_phase_Bz.append(_mk_phase_Bz(mode, 0, 0, phik))
        fct_amp_By = np.array(fct_amp_By)
        fct_amp_Bz = np.array(fct_amp_Bz)
        fct_phase_By = np.array(fct_phase_By)
        fct_phase_Bz = np.array(fct_phase_Bz)
    else :
        for mode_t in modes_trans :
            for mode_l in modes_longi :
                #fct_amp_By.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                #fct_amp_Bz.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                #fct_phase_By.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBy(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                #fct_phase_Bz.append(lambda y,imode_tmp=0,imode_t_tmp=mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBz(y,imode_=imode_tmp,imode_t_=imode_t_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                fct_amp_By.append(_mk_amp_By(0, mode_t, mode_l, phik))
                fct_amp_Bz.append(_mk_amp_Bz(0, mode_t, mode_l, phik))
                fct_phase_By.append(_mk_phase_By(0, mode_t, mode_l, phik))
                fct_phase_Bz.append(_mk_phase_Bz(0, mode_t, mode_l, phik))
        fct_amp_By = np.reshape(np.array(fct_amp_By),(2*m_trans+1,2*m_longi+1))
        fct_amp_Bz = np.reshape(np.array(fct_amp_Bz),(2*m_trans+1,2*m_longi+1))
        fct_phase_By = np.reshape(np.array(fct_phase_By),(2*m_trans+1,2*m_longi+1))
        fct_phase_Bz = np.reshape(np.array(fct_phase_Bz),(2*m_trans+1,2*m_longi+1))

    if temporal_smoothing=='Broadband':
        for mode in modes_broadband :
            im = int(mode+m_broadband)
            Laser(
                box_side       = box_side,
                omega          = omega*(1.+(mode+mode_offsets[im])*omega_m/omega),
                # omega          = omega,
                # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                time_envelope  = time_envelope,
                space_envelope = [fct_amp_By[im],fct_amp_Bz[im]],
                phase          = [fct_phase_By[im],fct_phase_Bz[im]],
                delay_phase    = delay_phase
            )
    else :
        if mode2generate_trans != None :
            mode_t = 1.*mode2generate_trans
            if mode2generate_longi != None :
                mode_l = 1.*mode2generate_longi
                im_t = int(mode_t+m_trans)
                im_l = int(mode_l+m_longi)
                Laser(
                    box_side       = box_side,
                    omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                    # omega          = omega,
                    # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                    time_envelope  = time_envelope,
                    space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                    phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                    delay_phase    = delay_phase
                )
            else :
                for mode_l in modes_longi :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
        else :
            if mode2generate_longi != None :
                mode_l = mode2generate_longi
                for mode_t in modes_trans :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
            else :
                for mode_t in modes_trans :
                    for mode_l in modes_longi :
                        im_t = int(mode_t+m_trans)
                        im_l = int(mode_l+m_longi)
                        Laser(
                            box_side       = box_side,
                            omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                            # omega          = omega,
                            # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                            time_envelope  = time_envelope,
                            space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                            phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                            delay_phase    = delay_phase
                        )

def LaserEnvelopeGaussian2D( a0=1., omega=1., focus=None, waist=3., time_envelope=tconstant(),
        envelope_solver = "explicit",box_side = "inside",Envelope_boundary_conditions = [["reflective"]],
        polarization_phi = 0.,ellipticity = 0.):
    import cmath
    from numpy import exp, sqrt, arctan, vectorize
    assert len(focus)==2, "LaserEnvelopeGaussian2D: focus must be a list of length 2."

    def gaussian_beam2D(x,y):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        Zr = omega * waist**2/2.
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        coeff = omega * (x-focus[0]) * w**2 / (2.*Zr**2)
        phase = coeff * ( (y-focus[1])**2 )
        exponential_with_total_phase = exp(1j*(phase-arctan( (x-focus[0])/Zr )))
        invWaist2 = (w/waist)**2
        spatial_amplitude = a0 *polarization_amplitude_factor * sqrt(w) * exp( -invWaist2*(y-focus[1])**2)
        return spatial_amplitude * exponential_with_total_phase
        
    if (box_side=="inside"):
        def envelope_profile(x,y,t):
            return gaussian_beam2D(x,y)*vectorize(time_envelope)(t)
    elif (box_side=="xmin"):
        def envelope_profile(y,t):
            return gaussian_beam2D(0,y)*vectorize(time_envelope)(t)
    else:
        print("LaserEnvelope error: box_side must be either 'inside' or 'xmin'. ")

    # Create Laser Envelope
    LaserEnvelope(
        omega                        = omega,
        envelope_profile             = envelope_profile,
        envelope_solver              = envelope_solver,
        box_side      = box_side,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi             = polarization_phi,
        ellipticity                  = ellipticity
    )

def LaserEnvelopeSquareFlattenedGaussian2D( a0=1., omega=1., focus=None, waist=3., N=10, time_envelope=tconstant(),
        envelope_solver = "explicit",box_side = "inside",Envelope_boundary_conditions = [["reflective"]],
        polarization_phi = 0.,ellipticity = 0.,flattened_intensity_position="far_from_focus"):
    import numpy as np
    import scipy.special as sp
    assert len(focus)==2, "LaserEnvelopeSquareFlattenedGaussian2D: focus must be a list of length 2."
    assert isinstance(N, int), "LaserEnvelopeSquareFlattenedGaussian2D: N must be an integer."
    assert (
            flattened_intensity_position == "far_from_focus"
            or flattened_intensity_position == "at_focus"
           ), (
           "LaserEnvelopeSquareFlattenedGaussian2D: flattened_intensity_position must be either 'at_focus' or 'far_from_focus'."
           )

    # Effective waist
    waist_corrected = (
                       waist * np.sqrt(N + 1) if flattened_intensity_position == "far_from_focus"
                       else waist / np.sqrt(N + 1)
                       )

    # Effective Rayleigh length, normalized
    x_R  = omega * waist_corrected**2/2.

    # Polarization amplitude factor
    polarization_amplitude_factor = 1/np.sqrt(1.+ellipticity**2)

    # Store the Hermite-Gauss (HG) mode coefficients
    cn = np.zeros(N+1)
    for n in range(N+1):
        m_values = np.arange(n, N+1)
        # computing this in log scale and then using the exponential
        # avoids overflow for high N
        # remember that gamma(n+1)=n!
        log_terms = (
                    -3*m_values*np.log(2.)
                    + sp.gammaln(2*m_values+1)
                    - sp.gammaln(m_values+1)
                    - sp.gammaln(m_values-n+1))
        cn[n] = np.sum(np.exp(log_terms - np.max(log_terms))) \
                    * np.exp(np.max(log_terms) - sp.gammaln(2*n+1))

    # The HG mode coefficients have alternating signs
    # if the flattened profile is far from focus
    if flattened_intensity_position == "far_from_focus":
        cn = cn*(-1.)**np.arange(N+1)

    # Normalizing quantity to have a normalized peak before the multiplication by a0
    # This recursive way of computing it avoids the overflow given by a brute force calculation
    # of the factorial
    def S_N(N):
        S = 0.0
        k = 1.0  # k_0 = 1
        for m in range(0, N+1):
            S += k
            # update a -> a_{m+1}
            k *= (2*m + 1) / (2*(m + 1))
        return S

    normalization_constant = 1. if (flattened_intensity_position=="at_focus") else S_N(N)

    # Recursive definition of the HG modes with even index
    def even_hermite_polynomials(x, N):
        # Returns an array of Hermite polynomials with even index using recursion relations
        H_even = np.empty((N+1,) + np.shape(x), dtype=float)
        H0 = np.ones_like(x)
        H_even[0] = H0
        if N == 0:
            return H_even
        H1   = 2*x
        Hnm2 = H0
        Hnm1 = H1
        even = 1

        for k in range(1, 2*N):
            Hn = 2*x*Hnm1 - 2*k*Hnm2
            Hnm2 = Hnm1
            Hnm1 = Hn
            if (k+1) % 2 == 0:
                H_even[even] = Hn
                even += 1

        return H_even

    # Square flattened Gauss definition in 2D Cartesian geometry
    def rectangular_flattened_Gaussian_beam2D(x,y):
        x, y             = np.broadcast_arrays(x, y)
        # HG mode waist at x
        w                = waist_corrected * np.sqrt(1 + ((x-focus[0]) /x_R)**2)
        # HG mode Gouy phase argument at x, to multiply by the the mode factor
        Gouy_phase_arg   = np.arctan2((x-focus[0]),x_R)
        # HG mode curved wavefront at x
        one_ov_R         = (x - focus[0]) / ((x - focus[0])**2 + x_R**2)
        curved_phase_y   = np.exp(1j * omega * (y-focus[1])**2 * one_ov_R / 2 )
        # HG mode exponential decay
        exp_along_y      = np.exp(-(y-focus[1])**2 / w**2)
        # Precompute the Hermite polynomials
        mask             = exp_along_y > 0.
        y_scaled         = np.sqrt(2) * (y-focus[1]) / w
        H = np.zeros((N+1,) + y.shape, dtype=float)
        if np.any(mask):
            # Only evaluate Hermite polynomials where the Gaussian has not underflown to zero
            H[:, mask]   = even_hermite_polynomials(y_scaled[mask], N)
            # Convert the nan to zero, that can happen only when the exponential is ~0
            H            = np.nan_to_num(H,nan=0.0,posinf=0.0,neginf=0.0)
        # Sum the HG modes parts that change for each mode
        HG_field_along_y = np.zeros_like(y,dtype=complex)
        for n in range(0, N+1):
            HG_field_along_y += cn[n] * H[n] * np.exp(-1j*(2*n+1/2.)*Gouy_phase_arg)
        # Multiply by the part in common for all modes
        HG_field_along_y = HG_field_along_y * exp_along_y * curved_phase_y * np.sqrt(waist_corrected/w)

        return a0*omega*polarization_amplitude_factor*HG_field_along_y/normalization_constant

    # Return the complex envelope, including the temporal envelope
    if (box_side=="inside"):
        def envelope_profile(x,y,t):
            return rectangular_flattened_Gaussian_beam2D(x,y)*(time_envelope)(t))
    elif (box_side=="xmin"):
        def envelope_profile(y,t):
            return rectangular_flattened_Gaussian_beam2D(0,y)*(time_envelope)(t))
    else:
        print("LaserEnvelope error: box_side must be either 'inside' or 'xmin'. ")

    # Create Laser Envelope
    LaserEnvelope(
        omega                        = omega,
        envelope_profile             = envelope_profile,
        envelope_solver              = envelope_solver,
        box_side                     = box_side,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi             = polarization_phi,
        ellipticity                  = ellipticity
    )


def LaserGaussian3D( box_side="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=[0.,0.],
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phase_offset=0.):
    from math import pi, cos, sin, tan, atan, sqrt, exp
    assert len(focus)==3, "LaserGaussian3D: focus must be a list of length 3."
    global Main
    assert len(Main)==1, "LaserGaussian3D profile has been defined before `Main()`"
    grid_length = Main.grid_length
    # Polarization and amplitude
    [dephasing, amplitudeZ, amplitudeY] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
    # Injection on ymin/ymax or zmin/zmax
    if box_side[0] == "y":
        focus = [focus[1],focus[0],focus[2]]
        grid_length = [grid_length[1],grid_length[0],grid_length[2]]
        amplitudeY = -amplitudeY
    elif box_side[0] == "z":
        focus = [focus[2],focus[0],focus[1]]
        grid_length = [grid_length[2],grid_length[0],grid_length[1]]
    # Injection on max boundary
    if box_side.endswith("max"):
        focus[0] = grid_length[0] - focus[0]
    # Space and phase envelopes
    Zr = omega * waist**2/2.
    if incidence_angle == [0.,0.]:
        w  = sqrt(1./(1.+(focus[0]/Zr)**2))
        invWaist2 = (w/waist)**2
        coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
        def spatial(y,z):
            return w * exp( -invWaist2*((y-focus[1])**2 + (z-focus[2])**2 )  )
        def phase(y,z):
            return coeff * ( (y-focus[1])**2 + (z-focus[2])**2 )
    else:
        invZr = 1./Zr
        invW  = 1./waist
        alpha = omega * Zr
        cy = cos(incidence_angle[0]); sy = sin(incidence_angle[0])
        cz = cos(incidence_angle[1]); sz = sin(incidence_angle[1])
        cycz = cy*cz; cysz = cy*sz; sycz = sy*cz; sysz = sy*sz
        amplitudeZ = sysz * amplitudeY + cy * amplitudeZ
        amplitudeY *= cz
        def spatial(y,z):
            X = invZr * (-focus[0]*cycz + (y-focus[1])*cysz - (z-focus[2])*sy )
            Y = invW  * ( focus[0]*sz   + (y-focus[1])*cz                     )
            Z = invW  * (-focus[0]*sycz + (y-focus[1])*sysz + (z-focus[2])*cy )
            invW2 = 1./(1.+X**2)
            return sqrt(invW2) * exp(-(Y**2+Z**2)*invW2)
        def phase(y,z):
            X = invZr * (-focus[0]*cycz + (y-focus[1])*cysz - (z-focus[2])*sy )
            Y = invZr * ( focus[0]*sz   + (y-focus[1])*cz                     )
            Z = invZr * (-focus[0]*sycz + (y-focus[1])*sysz + (z-focus[2])*cy )
            return alpha * X*(1.+0.5*(Y**2+Z**2)/(1.+X**2)) - atan(X)
        # Adjust the phase to match that of a laser that could come from another face
        faces = (focus[0],focus[1],focus[2],focus[1]-grid_length[1],focus[2]-grid_length[2])
        denominators = (cycz, cysz, -sy, cysz, -sy)
        distance_to_boundary = min([N/D for N,D in zip(faces,denominators) if D != 0 and N/D > 0])
        phase_offset -= omega * distance_to_boundary - atan(distance_to_boundary/Zr)
    # Create Laser
    Laser(
        box_side       = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y,z:amplitudeY*spatial(y,z), lambda y,z:amplitudeZ*spatial(y,z) ],
        phase          = [ lambda y,z:phase(y,z)-phase_offset+dephasing, lambda y,z:phase(y,z)-phase_offset ],
        delay_phase    = [ 0., dephasing ]
    )

# We will assume in 3D that angle is only in the (x,y) plane
# Transverse SSD along Y or Z in order to have a transverse SSD in the angle plane or perpendicular to the plane

def rotation_3d(x,y,z,ang) :
    '''
    Lineare tranformation: Rotation matrix
    (x,y)->(x',y')
    '''
    from math import cos, sin
    xrot = +cos(ang)*x + sin(ang)*y
    yrot = -sin(ang)*x + cos(ang)*y
    return xrot,yrot,z

def transform_3d(x,y,z,xf,yf,zf,L,ang) :
    '''
    Function to transform coordinate of laser-RPP formula
    x,y : Lab/Simulation box coordinate
    X,Y : Coordinate where X is the propagation axis of the laser, Y is transvers axis
    X,Y are rotated and translated coordinate in order the user define the 'focal spot' xf,yf in box coordinate and the angle of incidence with respect of x-axis of the simulated box
    '''
    from math import cos,sin,tan
    X,Y,Z = rotation_3d(x,y,z,ang)
    X = X+(L-xf/cos(ang))-(yf-tan(ang)*xf)*sin(ang)
    Y = Y-(yf-tan(ang)*xf)*cos(ang)
    Z = Z-zf
    return X,Y,Z

def LaserSmoothing3D(box_side="xmin", a0=1., omega=1., focus=None, incidence_angle=0.,polarization_phi=0.,ellipticity=0.,phase_zero=0.,
               Lf=3.00e6,fnumber=8.00,
               N=[6,6],rpp_random_seed=42,
               temporal_smoothing=None,
               omega_m=0.,modulation_depth=0,spectral_profile=lambda w: 1.,frequency_comb=False,mode_locking=False,
               temporal_freq_random_seed=1789,temporal_phi_random_seed=1793,
               rpp_per_mode=False,rpp_seed_per_mode=[42],
               omega_m_trans=0.,modulation_depth_trans=0,mode2generate_trans=None,direction='y',chirp_profile=tconstant(),
               omega_m_longi=0.,modulation_depth_longi=0,mode2generate_longi=None,
               space_envelope=lambda y,z:1.,time_envelope=tconstant()):
    '''
    Default values are in code units
    incidence_angle in radian ONLY IN (X,Y) PLANE
    a0                     : Maximum of the envelope at focal spot for 1 speckle (i.e. N=40 and no random phase between element, or N=1). Otherwise, for N=40, a = a0/sqrt(N=40) in the simulation box.
    Lf                     : Longueur focale without SSD
    fnumber                : F-number
    N                      : List of number of phase plate element per direction (for Ntot=36, then N=[6,6])
    rpp_random_seed        : Seed in order to have a Random Phase Plate (None is = no random, all element have zero phase-shift),
    temporal_smoothing     : None/'Broadband'/'TSSD'/'LSSD'
    omega_m                : modulation frequency for Broadband Laser
    modulation_depth       : depth 'm' of modulation and frequency bandwith = 2m for Broadband Laser
    frequency_comb         : False/True mode are equally spaced (modulation time) or not (no modulation time)
    mode_locking           : False/True mode have randomly chosen phase or not
    rpp_per_mode           : False/'Stardriver'/'ISI' : Same RPP for each mode / Change the RPP for each mode / Phase delay per mode
    rpp_seed_per_mode      : Seed for RRP or [parameter] for ISI echelon
    omega_m_trans          : modulation frequency for transverse TSSD
    omega_m_longi          : modulation frequency for longitudinal LSSD
    modulation_depth_trans : depth 'm' of modulation and frequency bandwith = 2m for transverse SSD
    modulation_depth_longi : depth 'm' of modulation and frequency bandwith = 2m for longitudinal SSD
    direction              : direction of transverse TSSD : 'y' or 'z'
    '''

    import numpy as np
    from math import pi, sqrt, cos, sin, tan, fabs
    from cmath import exp,rect,polar
    from scipy.special import erf,jv
    
    global Main

    if temporal_smoothing==None:
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='Broadband':
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='TSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='LSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0

    if len(Main)==0:
        raise Exception("LaserSmoothing3D profile has been defined before `Main()`") 
        
    k0 = omega
    Ny,Nz = N[0],N[1]
    Ntot = int(Ny*Nz)
    waist_y = fnumber*Ny*(2.00*pi/omega)
    waist_z = fnumber*Nz*(2.00*pi/omega)
    D = Lf/fnumber #taille lame de phase : On fait l'hypothese/le choix ici que l'ouverture du faisceau est identique avant la lame de phase (D et fnumber pareil pour y et z)
    dy = D/Ny #taille element lame de phase y
    dz = D/Nz #taille element lame de phase z
    Ry = sqrt(waist_y*dy/pi)
    Rz = sqrt(waist_z*dz/pi)

    x_focus,y_focus,z_focus = focus[0],focus[1],focus[2]

    El = (a0*omega/Ntot)/np.sqrt(k0*dy**2/(2*Lf*np.pi**2))/np.sqrt(k0*dz**2/(2*Lf*np.pi**2))
    # Polarization and amplitude
    dephasing, amplitudeZ, amplitudeY = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= El * cos(incidence_angle)
    amplitudeZ *= El
    delay_phase = [0., dephasing]

    krpp_y = np.linspace(-D/2,D/2,Ny+1)
    krpp_z = np.linspace(-D/2,D/2,Nz+1)

    #phik_y = np.zeros(Ny)
    #phik_z = np.zeros(Nz)
    #if rpp_random_seed != None :
    #    np.random.seed(rpp_random_seed) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
    #    phikinit = np.random.rand(int(Ny+Nz))
    #    for i in range(0,Ny):
    #        phik_y[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi
    #    for i in range(Ny,Ny+Nz):
    #        phik_z[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi

    phik = np.zeros(Ntot)
    if rpp_random_seed != None :
        np.random.seed(rpp_random_seed) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
        phikinit = np.random.rand(Ntot)
        for i in range(0,Ntot):
            phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi

    m_broadband = modulation_depth
    m_trans = modulation_depth_trans
    m_longi = modulation_depth_longi
    alpha_t = 2*pi/D
    alpha_x = 1/omega

    modes_trans = range(-m_trans,m_trans+1,1)
    modes_longi = range(-m_longi,m_longi+1,1)
    modes_broadband = range(-m_broadband,m_broadband+1,1)
    
    if temporal_smoothing=='Broadband':
        if frequency_comb :
            mode_offsets = np.zeros(len(modes_broadband))
        else :
            np.random.seed(temporal_freq_random_seed)
            mode_offsets = np.random.uniform(-0.5,0.5,len(modes_broadband))
        if mode_locking :
            phase_w = np.zeros(len(modes_broadband)) 
        else :
            np.random.seed(temporal_phi_random_seed)
            phase_w = 2*pi*np.random.rand(len(modes_broadband))
        
        # Legacy
        # Ebb = np.ones(len(modes_broadband))/np.sqrt(len(modes_broadband))

        if not callable(spectral_profile):
            raise Exception("spectral_profile must be a callable function of omega, e.g. lambda w: 1.")

        # Example 1 for intensity spectrum : Gaussian
        # sigma = 10 * omega_m
        # spectral_profile = lambda w: np.exp(-(w - omega)**2 / (2*sigma**2))
        # Example 2 for intensity spectrum : Lorentzian
        # gamma = 5 * omega_m
        # spectral_profile = lambda w: gamma / ((w - omega)**2 + gamma**2)
        # Default : Flat
        # spectral_profile = lambda w: 1.

        Ebb = np.array([spectral_profile(omega * (1. + (mode + mode_offsets[int(mode+m_broadband)]) * omega_m / omega)) for mode in modes_broadband])
        # Normalization : sum of all spectral intensity line = 1
        Ebb = Ebb / np.sqrt(np.sum(Ebb**2))


    def ERPP(y,z,imode,imode_tY,imode_tZ,imode_l,phik) :
        '''
        Formula (24) of "Cross-beam energy transfer between spatially smoothed laser beams" [A. Oudin, A. Debayle, C. Ruyer]
        Spatial envelope definition at x=0 of the simulation domain.
        For a fixed SSD mode. SSD is treated as a Laser() superposition at different frequency.
        '''

        X,Y,Z = transform_3d(0,y,z,x_focus,y_focus,z_focus,Lf,incidence_angle)

        Lfw = Lf*(1+alpha_x*imode_l*omega_m_longi)
        # K   = sqrt( fabs( k0/(2*X) - 1/R**2 ) )
        K = sqrt( fabs( k0*(Lfw-X)/(2*X*Lfw) ) )
        sum_erfm = 0+0*1j
        sum_erfm_y = 0+0*1j
        sum_erfm_z = 0+0*1j
        factm    = 0+0*1j

        if (temporal_smoothing == None) & ((imode_tY != 0) | (imode_tZ != 0)  | (imode_l != 0)):
            raise Exception("Input inconsistency : No temporal smoothing selected but non-zero 'modulation depth'")

        # Valid for X<L
        if X<Lfw :
            factm = (1/np.sqrt(pi))*0.5*sqrt(Lfw/(Lfw-X))*exp(-1j*k0*Y**2/(2*(Lfw-X)) + 1j*(imode_tY*alpha_t*Y-(imode_tY*alpha_t)**2*X/(2*k0))*Lfw/(Lfw-X)) # Phase RPP+TSSD
            factm *= (1/np.sqrt(pi))*0.5*sqrt(Lfw/(Lfw-X))*exp(-1j*k0*Z**2/(2*(Lfw-X)) + 1j*(imode_tZ*alpha_t*Z-(imode_tZ*alpha_t)**2*X/(2*k0))*Lfw/(Lfw-X))
            for ny in range (0,Ny):
                for nz in range (0,Nz):
                    sum_erfm_y = (erf(exp(-1j*pi/4)*K*(krpp_y[ny+1]-(Y-imode_tY*alpha_t*X/k0)*Lfw/(Lfw-X))) - erf(exp(-1j*pi/4)*K*(krpp_y[ny]-(Y-imode_tY*alpha_t*X/k0)*Lfw/(Lfw-X))))
                    sum_erfm_z = (erf(exp(-1j*pi/4)*K*(krpp_z[nz+1]-(Z-imode_tZ*alpha_t*X/k0)*Lfw/(Lfw-X))) - erf(exp(-1j*pi/4)*K*(krpp_z[nz]-(Z-imode_tZ*alpha_t*X/k0)*Lfw/(Lfw-X))))
                    sum_erfm += sum_erfm_y*sum_erfm_z*exp(1j*phik[ny*Nz+nz])
        # At focus
        elif X==Lfw :
           factm = exp(-1j*pi*0.25)*sqrt(k0*dy*dy/(2*pi*pi*Lfw))*exp(1j*k0*Y**2/(2*Lfw))*sin(k0*d/(2*Lfw)*(Y-imode_tY*alpha_t*X/k0))/(k0*d/(2*Lfw)*(Y-imode_tY*alpha_t*X/k0))
           factm *= exp(-1j*pi*0.25)*sqrt(k0*dz*dz/(2*pi*pi*Lfw))*exp(1j*k0*Z**2/(2*Lfw))*sin(k0*d/(2*Lfw)*(Z-imode_tZ*alpha_t*X/k0))/(k0*d/(2*Lfw)*(Z-imode_tZ*alpha_t*X/k0))
           for ny in range (0,Ny):
               for nz in range (0,Nz):
                   sum_erfm_y = exp(-1j*k0/Lfw*(Y-imode_tY*alpha_t*X/k0)*(krpp_y[ny+1]+krpp_y[ny])/2)
                   sum_erfm_z = exp(-1j*k0/Lfw*(Z-imode_tZ*alpha_t*X/k0)*(krpp_z[nz+1]+krpp_z[nz])/2)
                   sum_erfm += sum_erfm_y*sum_erfm_z*exp(1j*phik[ny*Nz+nz])
        # Beyond focus
        else :
            factm = (1/np.sqrt(pi))*0.5*sqrt(Lfw/fabs(Lfw-X))*exp(-1j*k0*Y**2/(2*(Lfw-X)) + 1j*(imode_tY*alpha_t*Y-(imode_tY*alpha_t)**2*X/(2*k0))*Lfw/(Lfw-X))*exp(-1j*pi*0.5)# Phase RPP+TSSD
            factm *= (1/np.sqrt(pi))*0.5*sqrt(Lfw/fabs(Lfw-X))*exp(-1j*k0*Z**2/(2*(Lfw-X)) + 1j*(imode_tZ*alpha_t*Z-(imode_tZ*alpha_t)**2*X/(2*k0))*Lfw/(Lfw-X))*exp(-1j*pi*0.5)
            for ny in range (0,Ny):
                for nz in range (0,Nz):
                    sum_erfm_y = (erf(exp(+1j*pi/4)*fabs(K)*(krpp_y[ny+1]-(Y-imode_tY*alpha_t*X/k0)*Lfw/(Lfw-X))) - erf(exp(+1j*pi/4)*fabs(K)*(krpp_y[ny]-(Y-imode_tY*alpha_t*X/k0)*Lfw/(Lfw-X))))
                    sum_erfm_z = (erf(exp(+1j*pi/4)*fabs(K)*(krpp_z[nz+1]-(Z-imode_tZ*alpha_t*X/k0)*Lfw/(Lfw-X))) - erf(exp(+1j*pi/4)*fabs(K)*(krpp_z[nz]-(Z-imode_tZ*alpha_t*X/k0)*Lfw/(Lfw-X))))
                    sum_erfm += sum_erfm_y*sum_erfm_z*exp(1j*phik[ny*Nz+nz])

        if temporal_smoothing==None:
            Einit = sum_erfm*factm*exp(1j*k0*X)#*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_t*omega_m_trans/k0)
        elif (temporal_smoothing=='TSSD') | (temporal_smoothing=='LSSD'):
            if direction=='y':
                Einit = jv(imode_tY,m_trans)*jv(imode_l,m_longi)*sum_erfm*factm*exp(1j*k0*X)*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_tY*omega_m_trans/k0)
            elif direction=='z':
                Einit = jv(imode_tZ,m_trans)*jv(imode_l,m_longi)*sum_erfm*factm*exp(1j*k0*X)*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_tZ*omega_m_trans/k0)
            else :
                raise Exception("'direction' parameter for TSSD unknown")
        elif temporal_smoothing=='Broadband':
            Einit = Ebb[imode]*sum_erfm*factm*exp(1j*k0*X)*exp(1j*phase_w[imode])*exp(1j*k0*X*imode*omega_m/k0)
        else :
            raise Exception("Temporal_smoothing method not implemented yet")

        Amp,Phase = polar(Einit)
        return Amp,Phase

    def ERPP_ampBz(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return amplitudeZ*space_envelope(y,z)*ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_ampBy(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return amplitudeY*space_envelope(y,z)*ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_phaseBz(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[0]
    def ERPP_phaseBy(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[1]

    # Encapsulated function to avoid warning
    def _mk_amp_By(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_ampBy(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)
    def _mk_amp_Bz(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_ampBz(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)
    def _mk_phase_By(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_phaseBy(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)
    def _mk_phase_Bz(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_phaseBz(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)

    fct_amp_By = []
    fct_amp_Bz = []
    fct_phase_By = []
    fct_phase_Bz = []

    if temporal_smoothing=='Broadband':
        phik_base = phik.copy()
        for mode in modes_broadband :
            if rpp_per_mode=='Stardriver':
                if len(rpp_seed_per_mode)!=len(modes_broadband):
                    raise Exception("len(rpp_seed_per_mode): "+str(len(rpp_seed_per_mode))+". len(modes_broadband): "+str(len(modes_broadband))+". Length of rpp_seed_per_mode have to be equal to 2 x modulation_depth + 1 ")
                phik = np.zeros(Ntot)
                np.random.seed(rpp_seed_per_mode[mode+m_broadband]) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
                phikinit = np.random.rand(Ntot)
                for i in range(0,Ntot):
                    phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi
            elif rpp_per_mode=='ISI':
                phik_mode = phik_base.copy()
                delta_t_y = rpp_seed_per_mode[0] * 2*pi / (2 * m_broadband * omega_m)
                delta_t_z = rpp_seed_per_mode[1] * 2*pi / (2 * m_broadband * omega_m) # Separable 1D echelon. Ny plus Nz independant beamlet and not Ny x Nz
                omega_k = omega * (1. + (mode+mode_offsets[mode+m_broadband]) * omega_m / omega)
                for ny in range(0, Ny):
                    for nz in range(0, Nz):
                        phik_mode[ny*Nz+nz] += omega_k * (ny * delta_t_y + nz * delta_t_z)
                phik_mode = phik_mode % (2*pi)
                phik = phik_mode
            else :
                phik = 1*phik
            #fct_amp_By.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_amp_Bz.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_By.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_Bz.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            fct_amp_By.append(_mk_amp_By(mode, 0, 0, 0, phik))
            fct_amp_Bz.append(_mk_amp_Bz(mode, 0, 0, 0, phik))
            fct_phase_By.append(_mk_phase_By(mode, 0, 0, 0, phik))
            fct_phase_Bz.append(_mk_phase_Bz(mode, 0, 0, 0, phik))
        fct_amp_By = np.array(fct_amp_By)
        fct_amp_Bz = np.array(fct_amp_Bz)
        fct_phase_By = np.array(fct_phase_By)
        fct_phase_Bz = np.array(fct_phase_Bz)
    else :
        for mode_t in modes_trans :
            for mode_l in modes_longi :
                if (m_trans==0) | (direction=='y') :
                    #fct_amp_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_amp_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    fct_amp_By.append(_mk_amp_By(0, mode_t, 0, mode_l, phik))
                    fct_amp_Bz.append(_mk_amp_Bz(0, mode_t, 0, mode_l, phik))
                    fct_phase_By.append(_mk_phase_By(0, mode_t, 0, mode_l, phik))
                    fct_phase_Bz.append(_mk_phase_Bz(0, mode_t, 0, mode_l, phik))
                elif direction=='z':
                    #fct_amp_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_amp_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    fct_amp_By.append(_mk_amp_By(0, 0, mode_t, mode_l, phik))
                    fct_amp_Bz.append(_mk_amp_Bz(0, 0,mode_t, 0, mode_l, phik))
                    fct_phase_By.append(_mk_phase_By(0, 0, mode_t, mode_l, phik))
                    fct_phase_Bz.append(_mk_phase_Bz(0, 0, mode_t, mode_l, phik))
                else :
                     raise Exception(" m_trans != 0 and direction for transverse SSD unknown : Have to be 'y' or 'z' ")
        fct_amp_By = np.reshape(np.array(fct_amp_By),(2*m_trans+1,2*m_longi+1))
        fct_amp_Bz = np.reshape(np.array(fct_amp_Bz),(2*m_trans+1,2*m_longi+1))
        fct_phase_By = np.reshape(np.array(fct_phase_By),(2*m_trans+1,2*m_longi+1))
        fct_phase_Bz = np.reshape(np.array(fct_phase_Bz),(2*m_trans+1,2*m_longi+1))

    if temporal_smoothing=='Broadband':
        for mode in modes_broadband :
            im = int(mode+m_broadband)
            Laser(
                box_side       = box_side,
                omega          = omega*(1.+(mode+mode_offsets[im])*omega_m/omega),
                # omega          = omega,
                # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                time_envelope  = time_envelope,
                space_envelope = [fct_amp_By[im],fct_amp_Bz[im]],
                phase          = [fct_phase_By[im],fct_phase_Bz[im]],
                delay_phase    = delay_phase
            )
    else :
        if mode2generate_trans != None :
            mode_t = 1.*mode2generate_trans
            if mode2generate_longi != None :
                mode_l = 1.*mode2generate_longi
                im_t = int(mode_t+m_trans)
                im_l = int(mode_l+m_longi)
                Laser(
                    box_side       = box_side,
                    omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                    # omega          = omega,
                    # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                    time_envelope  = time_envelope,
                    space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                    phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                    delay_phase    = delay_phase
                )
            else :
                for mode_l in modes_longi :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
        else :
            if mode2generate_longi != None :
                mode_l = mode2generate_longi
                for mode_t in modes_trans :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
            else :
                for mode_t in modes_trans :
                    for mode_l in modes_longi :
                        im_t = int(mode_t+m_trans)
                        im_l = int(mode_l+m_longi)
                        Laser(
                            box_side       = box_side,
                            omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                            # omega          = omega,
                            # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                            time_envelope  = time_envelope,
                            space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                            phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                            delay_phase    = delay_phase
                        )

def LaserSmoothingPeriodic3D(box_side="xmin", a0=1., omega=1., focus=None, incidence_angle=0.,polarization_phi=0.,ellipticity=0.,phase_zero=0.,
               Lf=3.00e6,fnumber=8.00,
               N=[6,6],rpp_random_seed=42,
               temporal_smoothing=None,
               omega_m=0.,modulation_depth=0,spectral_profile=lambda w: 1.,frequency_comb=False,mode_locking=False,
               temporal_freq_random_seed=1789,temporal_phi_random_seed=1793,
               rpp_per_mode=False,rpp_seed_per_mode=[42],
               omega_m_trans=0.,modulation_depth_trans=0,mode2generate_trans=None,direction='y',chirp_profile=tconstant(),
               omega_m_longi=0.,modulation_depth_longi=0,mode2generate_longi=None,
               space_envelope=lambda y,z:1.,time_envelope=tconstant()):
    '''
    Default values are in code units
    incidence_angle in radian ONLY IN (X,Y) PLANE
    a0                     : Maximum of the envelope at focal spot for 1 speckle (i.e. N=40 and no random phase between element, or N=1). Otherwise, for N=40, a = a0/sqrt(N=40) in the simulation box.
    Lf                     : Longueur focale without SSD
    fnumber                : F-number
    N                      : List of number of phase plate element per direction (for Ntot=36, then N=[6,6])
    rpp_random_seed        : Seed in order to have a Random Phase Plate (None is = no random, all element have zero phase-shift),
    temporal_smoothing     : None/'Broadband'/'TSSD'/'LSSD'
    omega_m                : modulation frequency for Broadband Laser
    modulation_depth       : depth 'm' of modulation and frequency bandwith = 2m for Broadband Laser
    frequency_comb         : False/True mode are equally spaced (modulation time) or not (no modulation time)
    mode_locking           : False/True mode have randomly chosen phase or not
    rpp_per_mode           : False/'Stardriver'/'ISI' : Same RPP for each mode / Change the RPP for each mode / Phase delay per mode
    rpp_seed_per_mode      : Seed for RRP or [parameter] for ISI echelon
    omega_m_trans          : modulation frequency for transverse TSSD
    omega_m_longi          : modulation frequency for longitudinal LSSD
    modulation_depth_trans : depth 'm' of modulation and frequency bandwith = 2m for transverse SSD
    modulation_depth_longi : depth 'm' of modulation and frequency bandwith = 2m for longitudinal SSD
    direction              : direction of transverse TSSD : 'y' or 'z'
    '''
    import numpy as np
    from math import pi, sqrt, cos, sin, tan, fabs
    from cmath import exp,rect,polar
    from scipy.special import erf,jv
    
    global Main

    if temporal_smoothing==None:
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='Broadband':
        omega_m_trans=0.
        modulation_depth_trans=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='TSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_longi=0.
        modulation_depth_longi=0
    elif temporal_smoothing=='LSSD':
        omega_m=0.
        modulation_depth=0
        omega_m_trans=0.
        modulation_depth_trans=0

    if len(Main)==0:
        raise Exception("LaserSmoothingPeriodic3D profile has been defined before `Main()`") 
        
    k0 = omega
    Ny,Nz = N[0],N[1]
    Ntot = int(Ny*Nz)
    waist_y = fnumber*Ny*(2.00*pi/omega)
    waist_z = fnumber*Nz*(2.00*pi/omega)
    D = Lf/fnumber #taille lame de phase : On fait l'hypothese/le choix ici que l'ouverture du faisceau est identique avant la lame de phase (D et fnumber pareil pour y et z)
    dy = D/Ny #taille element lame de phase y
    dz = D/Nz #taille element lame de phase z
    Ry = sqrt(waist_y*dy/pi)
    Rz = sqrt(waist_z*dz/pi)

    x_focus,y_focus,z_focus = focus[0],focus[1],focus[2]

    El = (a0*omega/Ntot)/np.sqrt(k0*dy**2/(2*Lf*np.pi**2))/np.sqrt(k0*dz**2/(2*Lf*np.pi**2))
    # Polarization and amplitude
    dephasing, amplitudeZ, amplitudeY = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= El * cos(incidence_angle)
    amplitudeZ *= El
    delay_phase = [0., dephasing]

    krpp_y = np.linspace(-D/2,D/2,Ny+1)
    krpp_z = np.linspace(-D/2,D/2,Nz+1)

    #phik_y = np.zeros(Ny)
    #phik_z = np.zeros(Nz)
    #if rpp_random_seed != None :
    #    np.random.seed(rpp_random_seed) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
    #    phikinit = np.random.rand(int(Ny+Nz))
    #    for i in range(0,Ny):
    #        phik_y[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi
    #    for i in range(Ny,Ny+Nz):
    #        phik_z[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi

    phik = np.zeros(Ntot)
    if rpp_random_seed != None :
        np.random.seed(rpp_random_seed) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
        phikinit = np.random.rand(Ntot)
        for i in range(0,Ntot):
            phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi

    m_broadband = modulation_depth
    m_trans = modulation_depth_trans
    m_longi = modulation_depth_longi
    alpha_t = 2*pi/D
    alpha_x = 1/omega

    modes_trans = range(-m_trans,m_trans+1,1)
    modes_longi = range(-m_longi,m_longi+1,1)
    modes_broadband = range(-m_broadband,m_broadband+1,1)
    
    if temporal_smoothing=='Broadband':
        if frequency_comb :
            mode_offsets = np.zeros(len(modes_broadband))
        else :
            np.random.seed(temporal_freq_random_seed)
            mode_offsets = np.random.uniform(-0.5,0.5,len(modes_broadband))
        if mode_locking :
            phase_w = np.zeros(len(modes_broadband)) 
        else :
            np.random.seed(temporal_phi_random_seed)
            phase_w = 2*pi*np.random.rand(len(modes_broadband))
        
        # Legacy
        # Ebb = np.ones(len(modes_broadband))/np.sqrt(len(modes_broadband))

        if not callable(spectral_profile):
            raise Exception("spectral_profile must be a callable function of omega, e.g. lambda w: 1.")

        # Example 1 for intensity spectrum : Gaussian
        # sigma = 10 * omega_m
        # spectral_profile = lambda w: np.exp(-(w - omega)**2 / (2*sigma**2))
        # Example 2 for intensity spectrum : Lorentzian
        # gamma = 5 * omega_m
        # spectral_profile = lambda w: gamma / ((w - omega)**2 + gamma**2)
        # Default : Flat
        # spectral_profile = lambda w: 1.

        Ebb = np.array([spectral_profile(omega * (1. + (mode + mode_offsets[int(mode+m_broadband)]) * omega_m / omega)) for mode in modes_broadband])
        # Normalization : sum of all spectral intensity line = 1
        Ebb = Ebb / np.sqrt(np.sum(Ebb**2))

    def ERPP(y,z,imode,imode_tY,imode_tZ,imode_l,phik) :
        '''
        Formula (24) of "Cross-beam energy transfer between spatially smoothed laser beams" [A. Oudin, A. Debayle, C. Ruyer]
        Spatial envelope definition at x=0 of the simulation domain.
        For a fixed SSD mode. SSD is treated as a Laser() superposition at different frequency.
        '''

        X,Y,Z = transform_3d(0,y,z,x_focus,y_focus,z_focus,Lf,incidence_angle)

        Lfw = Lf*(1+alpha_x*imode_l*omega_m_longi)
        # K   = sqrt( fabs( k0/(2*X) - 1/R**2 ) )
        K = sqrt( fabs( k0*(Lfw-X)/(2*X*Lfw) ) )
        sum_erfm = 0+0*1j
        sum_erfm_y = 0+0*1j
        sum_erfm_z = 0+0*1j
        factm    = 0+0*1j

        if (temporal_smoothing == None) & ((imode_tY != 0) | (imode_tZ != 0)  | (imode_l != 0)):
            raise Exception("Input inconsistency : No temporal smoothing selected but non-zero 'modulation depth'")

        factm = exp(-1j*pi*0.25)*sqrt(k0*dy*dy/(2*pi*pi*Lfw))
        factm *= exp(-1j*pi*0.25)*sqrt(k0*dz*dz/(2*pi*pi*Lfw))
        for ny in range (0,Ny):
            for nz in range (0,Nz):
                sum_erfm_y = exp(-1j*k0/Lfw*(Y-imode_tY*alpha_t*X/k0)*(krpp_y[ny+1]+krpp_y[ny])/2)
                sum_erfm_z = exp(-1j*k0/Lfw*(Z-imode_tZ*alpha_t*X/k0)*(krpp_z[nz+1]+krpp_z[nz])/2)
                sum_erfm += sum_erfm_y*sum_erfm_z*exp(1j*phik[ny*Nz+nz])

        if temporal_smoothing==None:
            Einit = sum_erfm*factm*exp(1j*k0*X)#*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_t*omega_m_trans/k0)
        elif (temporal_smoothing=='TSSD') | (temporal_smoothing=='LSSD'):
            if direction=='y':
                Einit = jv(imode_tY,m_trans)*jv(imode_l,m_longi)*sum_erfm*factm*exp(1j*k0*X)*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_tY*omega_m_trans/k0)
            elif direction=='z':
                Einit = jv(imode_tZ,m_trans)*jv(imode_l,m_longi)*sum_erfm*factm*exp(1j*k0*X)*exp(1j*k0*X*imode_l*omega_m_longi/k0)*exp(1j*k0*X*imode_tZ*omega_m_trans/k0)
            else :
                raise Exception("'direction' parameter for TSSD unknown")
        elif temporal_smoothing=='Broadband':
            Einit = Ebb[imode]*sum_erfm*factm*exp(1j*k0*X)*exp(1j*phase_w[imode])*exp(1j*k0*X*imode*omega_m/k0)
        else :
            raise Exception("Temporal_smoothing method not implemented yet")

        Amp,Phase = polar(Einit)
        return Amp,Phase

    def ERPP_ampBz(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return amplitudeZ*space_envelope(y,z)*ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_ampBy(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return amplitudeY*space_envelope(y,z)*ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[0]
    def ERPP_phaseBz(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[0]
    def ERPP_phaseBy(y,z,imode_,imode_tY_,imode_tZ_,imode_l_,phik_) :
        return ERPP(y,z,imode=imode_,imode_tY=imode_tY_,imode_tZ=imode_tZ_,imode_l=imode_l_,phik=phik_)[1]-phase_zero+delay_phase[1]

    # Encapsulated function to avoid warning
    def _mk_amp_By(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_ampBy(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)
    def _mk_amp_Bz(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_ampBz(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)
    def _mk_phase_By(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_phaseBy(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)
    def _mk_phase_Bz(im, itY, itZ, il, pk):
        return lambda y,z: ERPP_phaseBz(y, z, imode_=im, imode_tY_=itY, imode_tZ_=itZ, imode_l_=il, phik_=pk)

    fct_amp_By = []
    fct_amp_Bz = []
    fct_phase_By = []
    fct_phase_Bz = []

    if temporal_smoothing=='Broadband':
        phik_base = phik.copy()
        for mode in modes_broadband :
            if rpp_per_mode=='Stardriver':
                if len(rpp_seed_per_mode)!=len(modes_broadband):
                    raise Exception("len(rpp_seed_per_mode): "+str(len(rpp_seed_per_mode))+". len(modes_broadband): "+str(len(modes_broadband))+". Length of rpp_seed_per_mode have to be equal to 2 x modulation_depth + 1 ")
                phik = np.zeros(Ntot)
                np.random.seed(rpp_seed_per_mode[mode+m_broadband]) # Meme suite de nombre aleatoire utilisee pour comparer des cas avec meme lame de phase
                phikinit = np.random.rand(Ntot)
                for i in range(0,Ntot):
                    phik[i] = 2*pi*phikinit[i] # remplissage avec phi entre 0 et 2pi
            elif rpp_per_mode=='ISI':
                phik_mode = phik_base.copy()
                delta_t_y = rpp_seed_per_mode[0] * 2*pi / (2 * m_broadband * omega_m)
                delta_t_z = rpp_seed_per_mode[1] * 2*pi / (2 * m_broadband * omega_m) # Separable 1D echelon. Ny plus Nz independant beamlet and not Ny x Nz
                omega_k = omega * (1. + (mode+mode_offsets[mode+m_broadband]) * omega_m / omega)
                for ny in range(0, Ny):
                    for nz in range(0, Nz):
                        phik_mode[ny*Nz+nz] += omega_k * (ny * delta_t_y + nz * delta_t_z)
                phik_mode = phik_mode % (2*pi)
                phik = phik_mode
            else :
                phik = 1*phik
            #fct_amp_By.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_amp_Bz.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_ampBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_By.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            #fct_phase_Bz.append(lambda y,z,imode_tmp=mode,imode_tY_tmp=0,imode_tZ_tmp=0,imode_l_tmp=0,phik_tmp=phik: ERPP_phaseBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
            fct_amp_By.append(_mk_amp_By(mode, 0, 0, 0, phik))
            fct_amp_Bz.append(_mk_amp_Bz(mode, 0, 0, 0, phik))
            fct_phase_By.append(_mk_phase_By(mode, 0, 0, 0, phik))
            fct_phase_Bz.append(_mk_phase_Bz(mode, 0, 0, 0, phik))
        fct_amp_By = np.array(fct_amp_By)
        fct_amp_Bz = np.array(fct_amp_Bz)
        fct_phase_By = np.array(fct_phase_By)
        fct_phase_Bz = np.array(fct_phase_Bz)
    else :
        for mode_t in modes_trans :
            for mode_l in modes_longi :
                if (m_trans==0) | (direction=='y') :
                    #fct_amp_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_amp_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=1*mode_t,imode_tZ_tmp=0*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    fct_amp_By.append(_mk_amp_By(0, mode_t, 0, mode_l, phik))
                    fct_amp_Bz.append(_mk_amp_Bz(0, mode_t, 0, mode_l, phik))
                    fct_phase_By.append(_mk_phase_By(0, mode_t, 0, mode_l, phik))
                    fct_phase_Bz.append(_mk_phase_Bz(0, mode_t, 0, mode_l, phik))
                elif direction=='z':
                    #fct_amp_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_amp_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_ampBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_By.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBy(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    #fct_phase_Bz.append(lambda y,z,imode_tmp=0,imode_tY_tmp=0*mode_t,imode_tZ_tmp=1*mode_t,imode_l_tmp=mode_l,phik_tmp=phik: ERPP_phaseBz(y,z,imode_=imode_tmp,imode_tY_=imode_tY_tmp,imode_tZ_=imode_tZ_tmp,imode_l_=imode_l_tmp,phik_=phik_tmp))
                    fct_amp_By.append(_mk_amp_By(0, 0, mode_t, mode_l, phik))
                    fct_amp_Bz.append(_mk_amp_Bz(0, 0,mode_t, 0, mode_l, phik))
                    fct_phase_By.append(_mk_phase_By(0, 0, mode_t, mode_l, phik))
                    fct_phase_Bz.append(_mk_phase_Bz(0, 0, mode_t, mode_l, phik))
                else :
                     raise Exception(" m_trans != 0 and direction for transverse SSD unknown : Have to be 'y' or 'z' ")
        fct_amp_By = np.reshape(np.array(fct_amp_By),(2*m_trans+1,2*m_longi+1))
        fct_amp_Bz = np.reshape(np.array(fct_amp_Bz),(2*m_trans+1,2*m_longi+1))
        fct_phase_By = np.reshape(np.array(fct_phase_By),(2*m_trans+1,2*m_longi+1))
        fct_phase_Bz = np.reshape(np.array(fct_phase_Bz),(2*m_trans+1,2*m_longi+1))

    if temporal_smoothing=='Broadband':
        for mode in modes_broadband :
            im = int(mode+m_broadband)
            Laser(
                box_side       = box_side,
                omega          = omega*(1.+(mode+mode_offsets[im])*omega_m/omega),
                # omega          = omega,
                # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                time_envelope  = time_envelope,
                space_envelope = [fct_amp_By[im],fct_amp_Bz[im]],
                phase          = [fct_phase_By[im],fct_phase_Bz[im]],
                delay_phase    = delay_phase
            )
    else :
        if mode2generate_trans != None :
            mode_t = 1.*mode2generate_trans
            if mode2generate_longi != None :
                mode_l = 1.*mode2generate_longi
                im_t = int(mode_t+m_trans)
                im_l = int(mode_l+m_longi)
                Laser(
                    box_side       = box_side,
                    omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                    # omega          = omega,
                    # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                    time_envelope  = time_envelope,
                    space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                    phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                    delay_phase    = delay_phase
                )
            else :
                for mode_l in modes_longi :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega))
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
        else :
            if mode2generate_longi != None :
                mode_l = mode2generate_longi
                for mode_t in modes_trans :
                    im_t = int(mode_t+m_trans)
                    im_l = int(mode_l+m_longi)
                    Laser(
                        box_side       = box_side,
                        omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                        # omega          = omega,
                        # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                        time_envelope  = time_envelope,
                        space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                        phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                        delay_phase    = delay_phase
                    )
            else :
                for mode_t in modes_trans :
                    for mode_l in modes_longi :
                        im_t = int(mode_t+m_trans)
                        im_l = int(mode_l+m_longi)
                        Laser(
                            box_side       = box_side,
                            omega          = omega*(1.+mode_t*omega_m_trans/omega+mode_l*omega_m_longi/omega),
                            # omega          = omega,
                            # chirp_profile  = tpolynomial(t0=0.0, order0=(1.+mode*omega_m/omega)),
                            time_envelope  = time_envelope,
                            space_envelope = [fct_amp_By[im_t,im_l],fct_amp_Bz[im_t,im_l]],
                            phase          = [fct_phase_By[im_t,im_l],fct_phase_Bz[im_t,im_l]],
                            delay_phase    = delay_phase
                        )

def LaserEnvelopeGaussian3D( a0=1., omega=1., focus=None, waist=3., time_envelope=tconstant(),
        envelope_solver = "explicit",Envelope_boundary_conditions = [["reflective"]], box_side = "inside",
        polarization_phi = 0.,ellipticity = 0.):
    import cmath
    from numpy import exp, sqrt, arctan, vectorize
    assert len(focus)==3, "LaserEnvelopeGaussian3D: focus must be a list of length 3."
    
    def gaussian_beam3D(x,y,z):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        Zr = omega * waist**2/2.
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        coeff = omega * (x-focus[0]) * w**2 / (2.*Zr**2)
        phase = coeff * ( (y-focus[1])**2 + (z-focus[2])**2 )
        exponential_with_total_phase = exp(1j*(phase-arctan( (x-focus[0])/Zr )))
        invWaist2 = (w/waist)**2
        spatial_amplitude = a0*polarization_amplitude_factor* w * exp( -invWaist2*(  (y-focus[1])**2 + (z-focus[2])**2 )  )
        return spatial_amplitude * exponential_with_total_phase
    
    if (box_side=="inside"):
        def envelope_profile(x,y,z,t):
            return gaussian_beam3D(x,y,z)*vectorize(time_envelope)(t)
    elif (box_side=="xmin"):
        def envelope_profile(y,z,t):
            return gaussian_beam3D(0,y,z)*vectorize(time_envelope)(t)
    else:
        print("LaserEnvelope error: box_side must be either 'inside' or 'xmin'. ")

    # Create Laser Envelope
    LaserEnvelope(
        omega                        = omega,
        envelope_profile             = envelope_profile,
        envelope_solver              = envelope_solver,
        box_side                     = box_side,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi             = polarization_phi,
        ellipticity                  = ellipticity
    )


def LaserGaussianAM( box_side="xmin", a0=1., omega=1., focus=None, waist=3.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phase_offset=0.):
    from math import cos, sin, tan, atan, sqrt, exp
    if (len(focus)<1) or (len(focus)>2): 
        print("ERROR: focus should be a list of length 1")
        exit(1)
    elif (len(focus)==2):
        print("WARNING: deprecated focus in LaserGaussianAM should be a list of length 1")
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
    # Space and phase envelopes
    Zr = omega * waist**2/2.
    w  = sqrt(1./(1.+(focus[0]/Zr)**2))
    invWaist2 = (w/waist)**2
    coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
    def spatial(r):
        return w * exp( -invWaist2*(r)**2 )
    def phase(r):
        return coeff * (r)**2
    # Create Laser
    Laser(
        box_side        = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda r:amplitudeZ*spatial(r), lambda r:amplitudeY*spatial(r) ],
        phase          = [ lambda r:phase(r)-phase_offset+dephasing, lambda r:phase(r)-phase_offset ],
        delay_phase    = [ 0., dephasing ]
    )


def LaserEnvelopeGaussianAM( a0=1., omega=1., focus=None, waist=3., time_envelope=tconstant(),
        envelope_solver = "explicit",box_side = "inside",Envelope_boundary_conditions = [["reflective"]],
        Env_pml_sigma_parameters = [[0.90,2],[10.0,2],[10.0,2]],
        Env_pml_kappa_parameters = [[1.00,1.00,2],[1.00,1.00,2],[1.00,1.00,2]],
        Env_pml_alpha_parameters = [[0.90,0.90,1],[0.75,0.75,1],[0.75,0.75,1]],
        polarization_phi = 0.,ellipticity = 0.):
    import cmath
    from numpy import exp, sqrt, arctan, vectorize
    if (len(focus)<1) or (len(focus)>2): 
        print("ERROR: focus should be a list of length 1")
        exit(1)
    elif (len(focus)==2):
        print("WARNING: deprecated focus in LaserEnvelopeGaussianAM should be a list of length 1")

    def gaussian_beamAM(x,r):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        Zr = omega * waist**2/2.
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        coeff = omega * (x-focus[0]) * w**2 / (2.*Zr**2)
        phase = coeff * ( r**2 )
        exponential_with_total_phase = exp(1j*(phase-arctan( (x-focus[0])/Zr )))
        invWaist2 = (w/waist)**2
        spatial_amplitude = a0 * polarization_amplitude_factor * w * exp( -invWaist2*(  r**2  ) )
        return spatial_amplitude  * exponential_with_total_phase
        
    if (box_side=="inside"):
        def envelope_profile(x,r,t):
            return gaussian_beamAM(x,r)*vectorize(time_envelope)(t)
    elif (box_side=="xmin"):
        def envelope_profile(r,t):
            return gaussian_beamAM(0,r)*vectorize(time_envelope)(t)
    else:
        print("LaserEnvelope error: box_side must be either 'inside' or 'xmin'. ")
            
    # Create Laser Envelope
    LaserEnvelope(
        omega                        = omega,
        envelope_profile             = envelope_profile,
        envelope_solver              = envelope_solver,
        box_side                     = box_side,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        Env_pml_sigma_parameters     = Env_pml_sigma_parameters,
        Env_pml_kappa_parameters     = Env_pml_kappa_parameters,
        Env_pml_alpha_parameters     = Env_pml_alpha_parameters,
        polarization_phi             = polarization_phi,
        ellipticity                  = ellipticity
    )

# Define the tools for the propagation of a laser profile
try:
    import numpy as np
    
    _N_LaserOffset = 0
    
    def LaserOffset(box_side="xmin", space_time_profile=[], offset=0., angle=0., extra_envelope=lambda *a:1.,
            fft_time_window=None, fft_time_step=None, keep_n_strongest_modes=100,
            number_of_processes=None, file=None):
        global _N_LaserOffset
        
        file_ = file or ('LaserOffset'+str(_N_LaserOffset)+'.h5')
        
        L = Laser(
            box_side = box_side,
            file = file_,
        )
        
        L._offset = offset
        L._extra_envelope = extra_envelope
        L._profiles = space_time_profile
        L._fft_time_window = fft_time_window or Main.simulation_time
        L._fft_time_step = fft_time_step or Main.timestep
        L._keep_n_strongest_modes = keep_n_strongest_modes
        L._angle = angle
        L._number_of_processes = number_of_processes
        if file:
            if not os.path.exists(file):
                raise Exception("File not found or not accessible: "+file)
            L._propagate = False
        else:
            L._propagate = True
        
        _N_LaserOffset += 1

except:
    
    def LaserOffset(box_side="xmin", space_time_profile=[], offset=0., fft_time_window=None, extra_envelope=lambda *a:1., keep_n_strongest_modes=100, angle=0., number_of_processes=None, file=None):
        L = Laser(
            box_side = box_side,
            file = "none",
            time_envelope = extra_envelope
        )
        print("WARNING: LaserOffset unavailable because numpy was not found")

