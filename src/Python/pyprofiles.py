# Some predefined profiles (see doc)

def constant(value, xvacuum=-float("inf"), yvacuum=-float("inf")):
    global dim
    if not dim:
        raise Exception("constant profile has been defined before `dim`")
    if dim == "1d3v":
        f = lambda x  : value if x>=xvacuum else 0.
    if dim == "2d3v":
        f = lambda x,y: value if (x>=xvacuum and y>=yvacuum) else 0.
        f.yvacuum = yvacuum
    f.profileName = "constant"
    f.value   = value
    f.xvacuum = xvacuum
    return f
constant._reserved = True

def trapezoidal(max,
                xvacuum=0., xplateau=None, xslope1=0., xslope2=0.,
                yvacuum=0., yplateau=None, yslope1=0., yslope2=0. ):
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("trapezoidal profile has been defined before `dim` or `sim_length`")
    if len(sim_length)>0 and xplateau is None: xplateau = sim_length[0]-xvacuum
    if len(sim_length)>1 and yplateau is None: yplateau = sim_length[1]-yvacuum
    def fx(x):
        # vacuum region
        if x < xvacuum: return 0.
        # linearly increasing density
        elif x < xvacuum+xslope1: return max*(x-xvacuum) / xslope1
        # density plateau
        elif x < xvacuum+xslope1+xplateau: return max
        # linearly decreasing density
        elif x < xvacuum+xslope1+xplateau+xslope2:
            return max*(1. - ( x - (xvacuum+xslope1+xslope2) ) / xslope2)
        # beyond the plasma
        else: return 0.0
    if dim == "1d3v":
        f = fx
    if dim == "2d3v":
        def fy(y):
            # vacuum region
            if y < yvacuum: return 0.
            # linearly increasing density
            elif y < yvacuum+yslope1: return (y-yvacuum) / yslope1
            # density plateau
            elif y < yvacuum+yslope1+yplateau: return 1.
            # linearly decreasing density
            elif y < yvacuum+yslope1+yplateau+yslope2:
                return 1. - ( y - (yvacuum+yslope1+yslope2) ) / yslope2
            # beyond
            else: return 0.0
        f = lambda x,y: fx(x)*fy(y)
        f.yvacuum  = yvacuum
        f.yplateau = yplateau
        f.yslope1  = yslope1
        f.yslope2  = yslope2
    f.profileName = "trapezoidal"
    f.value    = max
    f.xvacuum  = xvacuum
    f.xplateau = xplateau
    f.xslope1  = xslope1
    f.xslope2  = xslope2
    return f

def gaussian(max,
             xvacuum=0., xlength=float("inf"), xfwhm=None, xcenter=None, xorder=2,
             yvacuum=0., ylength=float("inf"), yfwhm=None, ycenter=None, yorder=2 ):
    import math
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("gaussian profile has been defined before `dim` or `sim_length`")
    if len(sim_length)>0:
        if xlength is None: xlength = sim_length[0]-xvacuum
        if xfwhm   is None: xfwhm   = xlength/3.
        if xcenter is None: xcenter = xvacuum + xlength/2.
    if len(sim_length)>1: 
        if ylength is None: ylength = sim_length[1]-yvacuum
        if yfwhm   is None: yfwhm   = ylength/3.
        if ycenter is None: ycenter = yvacuum + ylength/2.
    xsigma = (0.5*xfwhm)**xorder/math.log(2.0)
    def fx(x):
        # vacuum region
        if x < xvacuum: return 0.
        # gaussian
        elif x < xvacuum+xlength: return max*math.exp( -(x-xcenter)**xorder / xsigma )
        # beyond
        else: return 0.0
    if dim == "1d3v":
        f = fx
    if dim == "2d3v":
        ysigma = (0.5*yfwhm)**yorder/math.log(2.0)
        def fy(y):
            if yorder == 0: return 1.
            # vacuum region
            if y < yvacuum: return 0.
            # gaussian
            elif y < yvacuum+ylength: return math.exp( -(y-ycenter)**yorder / ysigma )
            # beyond
            else: return 0.0
        f = lambda x,y: fx(x)*fy(y)
        f.yvacuum = yvacuum
        f.ylength = ylength
        f.ysigma  = ysigma
        f.ycenter = ycenter
        f.yorder  = yorder
    f.profileName = "gaussian"
    f.value   = max
    f.xvacuum = xvacuum
    f.xlength = xlength
    f.xsigma  = xsigma
    f.xcenter = xcenter
    f.xorder  = xorder
    return f


def polygonal(xpoints=[], xvalues=[]):
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("polygonal profile has been defined before `dim` or `sim_length`")
    if len(xpoints)!=len(xvalues):
        raise Exception("polygonal profile requires as many points as values")
    if len(sim_length)>0 and len(xpoints)==0:
        xpoints = [0., sim_length[0]]
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
           yamplitude=1., yvacuum=0., ylength=None, yphi=0., ynumber=2):
    import math
    global dim, sim_length
    if not dim or not sim_length:
        raise Exception("cosine profile has been defined before `dim` or `sim_length`")
    
    if len(sim_length)>0 and xlength is None: xlength = sim_length[0]-xvacuum
    if len(sim_length)>1 and ylength is None: ylength = sim_length[1]-yvacuum
    
    def fx(x):
        #vacuum region
        if x < xvacuum: return 0.
        # profile region
        elif x < xvacuum+xlength:
            return base + xamplitude * math.cos(xphi + 2.*math.pi * xnumber * (x-xvacuum)/xlength)
        # beyond
        else: return 0.
    if dim == "1d3v":
        f = fx
    if dim == "2d3v":
        def fy(y):
            #vacuum region
            if y < yvacuum: return 0.
            # profile region
            elif y < yvacuum+ylength:
                return base + yamplitude * math.cos(yphi + 2.*math.pi * ynumber * (y-yvacuum)/ylength)
            # beyond
            else: return 0.
        f = lambda x,y: fx(x)*fy(y)
        f.yamplitude  = yamplitude
        f.yvacuum     = yvacuum
        f.ylength     = ylength
        f.yphi        = yphi
        f.ynumber     = float(ynumber)
    f.profileName = "cosine"
    f.base        = base
    f.xamplitude  = xamplitude
    f.xvacuum     = xvacuum
    f.xlength     = xlength
    f.xphi        = xphi
    f.xnumber     = float(xnumber)
    return f

def polynomial(**kwargs):
    global dim
    if not dim:
        raise Exception("polynomial profile has been defined before `dim`")
    x0 = 0.
    y0 = 0.
    coeffs = dict()
    for k, a in kwargs.iteritems():
        if   k=="x0":
            x0 = a
        elif k=="y0":
            y0 = a
        elif k[:5]=="order":
            if type(a) is not list: a = [a]
            order = int(k[5:])
            coeffs[ order ] = a
            if dim=="1d3v" and len(a)!=1:
                raise Exception("1D polynomial profile must have one coefficient per order")
    if dim=="1d3v":
        def f(x):
            r = 0.
            xx0 = x-x0
            xx = 1.
            currentOrder = 0
            for order, c in sorted(coeffs.iteritems()):
                while currentOrder<order:
                    currentOrder += 1
                    xx *= xx0
                r += c[0] * xx
            return r
    if dim=="2d3v":
        def f(x,y):
            r = 0.
            xx0 = x-x0
            yy0 = y-y0
            xx = [1.]
            currentOrder = 0
            for order, c in sorted(coeffs.iteritems()):
                while currentOrder<order:
                    currentOrder += 1
                    last = xx[-1]*yy0
                    xx = [ xxx * xx0 for xxx in xx ]
                    xx.append(last)
                for i in range(order+1): r += c[i]*xx[i]
            return r
    f.profileName = "polynomial"
    f.x0 = x0
    f.y0 = y0
    f.orders = []
    f.coeffs = []
    for order, c in sorted(coeffs.iteritems()):
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
    global sim_time
    if not sim_time:
        raise Exception("ttrapezoidal profile has been defined before `sim_time`")
    if plateau is None: plateau = sim_time - start
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
    global sim_time
    if not sim_time:
        raise Exception("tgaussian profile has been defined before `sim_time`")
    if duration is None: duration = sim_time-start
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
    global sim_time
    if not sim_time:
        raise Exception("tpolygonal profile has been defined before `sim_time`")
    if len(points)==0:
        points = [0., sim_time]
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
    global sim_time
    if not sim_time:
        raise Exception("tcosine profile has been defined before `sim_time`")
    if duration is None: duration = sim_time-start
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
    for k, a in kwargs.iteritems():
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
        for order, c in sorted(coeffs.iteritems()):
            while currentOrder<order:
                currentOrder += 1
                tt *= tt0
            r += c * tt
        return r
    f.profileName = "tpolynomial"
    f.t0 = t0
    f.orders = []
    f.coeffs = []
    for order, c in sorted(coeffs.iteritems()):
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


def LaserPlanar1D( boxSide="west", a0=1., omega=1.,
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



def LaserGaussian2D( boxSide="west", a0=1., omega=1., focus=None, waist=3., incidence_angle=0.,
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





