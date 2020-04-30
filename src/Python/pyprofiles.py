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
    from math import pi, sqrt, sin, cos, tan, atan
    e2 = ellipticity**2
    p = (1.-e2)*sin(2.*polarization_phi)/2.
    if abs(p) < 1e-10:
        dephasing = pi/2.
    else:
        dephasing = atan(ellipticity/p)
    amplitude = sqrt(1./(1.+e2))
    c2 = cos(polarization_phi)**2
    s2 = 1. - c2
    amplitudeY = amplitude * sqrt(c2 + e2*s2)
    amplitudeZ = amplitude * sqrt(s2 + e2*c2)
    return [dephasing, amplitudeY, amplitudeZ]

def LaserPlanar1D( box_side="xmin", a0=1., omega=1.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(),phaseZero=0.):
    import math
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
        phase          = [ dephasing-phaseZero, -phaseZero ],
        delay_phase    = [ 0., dephasing ]
    )

def LaserEnvelopePlanar1D( a0=1., omega=1., focus=None, time_envelope=tconstant(),
        envelope_solver = "explicit",Envelope_boundary_conditions = [["reflective"]],
        polarization_phi = 0.,ellipticity = 0.):
    import cmath
    from numpy import vectorize, sqrt

    def space_time_envelope(x,t):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        return (a0*omega*polarization_amplitude_factor) * complex( vectorize(time_envelope)(t) )

    # Create Laser Envelope
    LaserEnvelope(
        omega               = omega,
        envelope_profile    = space_time_envelope,
        envelope_solver     = envelope_solver,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi    = polarization_phi,
        ellipticity         = ellipticity
    )


def LaserGaussian2D( box_side="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=0.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phaseZero=0.):
    from math import cos, sin, tan, atan, sqrt, exp
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
    # Space and phase envelopes
    Zr = omega * waist**2/2.
    if incidence_angle == 0.:
        Y1 = focus[1]
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
        amplitudeZ *= cos(incidence_angle)
        def spatial(y):
            w2 = 1./(1. + invZr2*(y-Y1)**2)
            return sqrt(sqrt(w2)) * exp( -invWaist2*w2*(y-Y2)**2 )
        def phase(y):
            dy = y-Y1
            return omega_*dy*(1.+ invZr3*(y-Y2)**2/(1.+invZr2*dy**2)) + atan(invZr*dy)
        phaseZero += phase(Y2)
    # Create Laser
    Laser(
        box_side        = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y:amplitudeZ*spatial(y), lambda y:amplitudeY*spatial(y) ],
        phase          = [ lambda y:phase(y)-phaseZero+dephasing, lambda y:phase(y)-phaseZero ],
        delay_phase    = [ 0., dephasing ]
    )

def LaserEnvelopeGaussian2D( a0=1., omega=1., focus=None, waist=3., time_envelope=tconstant(),
        envelope_solver = "explicit",Envelope_boundary_conditions = [["reflective"]],
        polarization_phi = 0.,ellipticity = 0.):
    import cmath
    from numpy import exp, sqrt, arctan, vectorize

    def gaussian_beam_with_temporal_profile(x,y,t):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        Zr = omega * waist**2/2.
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        coeff = omega * (x-focus[0]) * w**2 / (2.*Zr**2)
        phase = coeff * ( (y-focus[1])**2 )
        exponential_with_total_phase = exp(1j*(phase-arctan( (x-focus[0])/Zr )))
        invWaist2 = (w/waist)**2
        spatial_amplitude = a0*polarization_amplitude_factor*omega * sqrt(w) * exp( -invWaist2*(y-focus[1])**2)
        space_time_envelope = spatial_amplitude * vectorize(time_envelope)(t)
        return space_time_envelope * exponential_with_total_phase

    # Create Laser Envelope
    LaserEnvelope(
        omega               = omega,
        envelope_profile    = gaussian_beam_with_temporal_profile,
        envelope_solver     = envelope_solver,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi    = polarization_phi,
        ellipticity         = ellipticity
    )

def LaserGaussian3D( box_side="xmin", a0=1., omega=1., focus=None, waist=3., incidence_angle=[0.,0.],
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phaseZero=0.):
    import math
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
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
        amplitudeZ = sysz * amplitudeY + cy * amplitudeZ
        amplitudeY *= cz
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
        delay_phase    = [ 0., dephasing ]
    )

def LaserEnvelopeGaussian3D( a0=1., omega=1., focus=None, waist=3., time_envelope=tconstant(),
        envelope_solver = "explicit",Envelope_boundary_conditions = [["reflective"]],
        polarization_phi = 0.,ellipticity = 0.):
    import cmath
    from numpy import exp, sqrt, arctan, vectorize

    def gaussian_beam_with_temporal_profile(x,y,z,t):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        Zr = omega * waist**2/2.
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        coeff = omega * (x-focus[0]) * w**2 / (2.*Zr**2)
        phase = coeff * ( (y-focus[1])**2 + (z-focus[2])**2 )
        exponential_with_total_phase = exp(1j*(phase-arctan( (x-focus[0])/Zr )))
        invWaist2 = (w/waist)**2
        spatial_amplitude = a0*omega*polarization_amplitude_factor* w * exp( -invWaist2*(  (y-focus[1])**2 + (z-focus[2])**2 )  )
        space_time_envelope = spatial_amplitude * vectorize(time_envelope)(t)
        return space_time_envelope * exponential_with_total_phase

    # Create Laser Envelope
    LaserEnvelope(
        omega               = omega,
        envelope_profile    = gaussian_beam_with_temporal_profile,
        envelope_solver     = envelope_solver,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi    = polarization_phi,
        ellipticity         = ellipticity
    )


def LaserGaussianAM( box_side="xmin", a0=1., omega=1., focus=None, waist=3.,
        polarization_phi=0., ellipticity=0., time_envelope=tconstant(), phaseZero=0.):
    from math import cos, sin, tan, atan, sqrt, exp
    assert len(focus)==2, "LaserGaussianAM: focus must be a list of length 2."
    # Polarization and amplitude
    [dephasing, amplitudeY, amplitudeZ] = transformPolarization(polarization_phi, ellipticity)
    amplitudeY *= a0 * omega
    amplitudeZ *= a0 * omega
    # Space and phase envelopes
    Zr = omega * waist**2/2.
    w  = sqrt(1./(1.+(focus[0]/Zr)**2))
    invWaist2 = (w/waist)**2
    coeff = -omega * focus[0] * w**2 / (2.*Zr**2)
    def spatial(y):
        return w * exp( -invWaist2*(y-focus[1])**2 )
    def phase(y):
        return coeff * (y-focus[1])**2
    # Create Laser
    Laser(
        box_side        = box_side,
        omega          = omega,
        chirp_profile  = tconstant(),
        time_envelope  = time_envelope,
        space_envelope = [ lambda y:amplitudeZ*spatial(y), lambda y:amplitudeY*spatial(y) ],
        phase          = [ lambda y:phase(y)-phaseZero+dephasing, lambda y:phase(y)-phaseZero ],
        delay_phase    = [ 0., dephasing ]
    )


def LaserEnvelopeGaussianAM( a0=1., omega=1., focus=None, waist=3., time_envelope=tconstant(),
        envelope_solver = "explicit",Envelope_boundary_conditions = [["reflective"]],
        polarization_phi = 0.,ellipticity = 0.):
    import cmath
    from numpy import exp, sqrt, arctan, vectorize

    def gaussian_beam_with_temporal_profile(x,r,t):
        polarization_amplitude_factor = 1/sqrt(1.+ellipticity**2)
        Zr = omega * waist**2/2.
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        coeff = omega * (x-focus[0]) * w**2 / (2.*Zr**2)
        phase = coeff * ( r**2 )
        exponential_with_total_phase = exp(1j*(phase-arctan( (x-focus[0])/Zr )))
        invWaist2 = (w/waist)**2
        spatial_amplitude = a0*omega *polarization_amplitude_factor* w * exp( -invWaist2*(  r**2  ) )
        space_time_envelope = spatial_amplitude * vectorize(time_envelope)(t)
        return space_time_envelope * exponential_with_total_phase

    # Create Laser Envelope
    LaserEnvelope(
        omega               = omega,
        envelope_profile    = gaussian_beam_with_temporal_profile,
        envelope_solver     = envelope_solver,
        Envelope_boundary_conditions = Envelope_boundary_conditions,
        polarization_phi    = polarization_phi,
        ellipticity         = ellipticity
    )

# Define the tools for the propagation of a laser profile
try:
    import numpy as np
    
    _N_LaserOffset = 0
    
    def LaserOffset(box_side="xmin", space_time_profile=[], offset=0., extra_envelope=lambda *a:1., keep_n_strongest_modes=100, angle=0.):
        global _N_LaserOffset
        
        file = 'LaserOffset'+str(_N_LaserOffset)+'.h5'
        
        L = Laser(
            box_side = "xmin",
            file = file,
        )
        
        L._offset = offset
        L._extra_envelope = extra_envelope
        L._profiles = space_time_profile
        L._keep_n_strongest_modes = keep_n_strongest_modes
        L._angle = angle
        L._propagate = True
        
        _N_LaserOffset += 1

except:
    
    def LaserOffset(box_side="xmin", space_time_profile=[], offset=0., time_envelope=1.):
        L = Laser(
            box_side = "xmin",
            file = "none",
            time_envelope = time_envelope
        )
        print("WARNING: LaserOffset unavailable because numpy was not found")

