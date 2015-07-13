# Some predefined profiles (see doc)

def constant(value, xvacuum=0., yvacuum=0.):
    global dim, sim_length
    if dim == "1d3v": return lambda x: value if x>=xvacuum else 0.
    if dim == "2d3v": return lambda x,y: value if (x>=xvacuum and y>=yvacuum) else 0.

def trapezoidal(max,
                xvacuum=0., xplateau=None, xslope1=0., xslope2=0.,
                yvacuum=0., yplateau=None, yslope1=0., yslope2=0. ):
    global dim, sim_length
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
    if dim == "1d3v": return fx
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def gaussian(max,
             xvacuum=0., xlength=None, xfwhm=None, xcenter=None, xorder=2,
             yvacuum=0., ylength=None, yfwhm=None, ycenter=None, yorder=2 ):
    import math
    global dim, sim_length
    if len(sim_length)>0:
        if xlength is None: xlength = sim_length[0]-xvacuum
        if xfwhm   is None: xfwhm   = xlength/3.
        if xcenter is None: xcenter = xvacuum + xlength/2.
    if len(sim_length)>1: 
        if ylength is None:ylength = sim_length[1]-yvacuum
        if yfwhm   is None:yfwhm   = ylength/3.
        if ycenter is None:ycenter = yvacuum + ylength/2.
    def fx(x):
        sigmaN = (0.5*xfwhm)**xorder/math.log(2.0)
        # vacuum region
        if x < xvacuum: return 0.
        # gaussian
        elif x < xvacuum+xlength: return max*math.exp( -(x-xcenter)**xorder / sigmaN )
        # beyond
        else: return 0.0
    def fy(y):
        if yorder == 0: return 1.
        sigmaN = (0.5*yfwhm)**yorder/math.log(2.0)
        # vacuum region
        if y < yvacuum: return 0.
        # gaussian
        elif y < yvacuum+ylength: return math.exp( -(y-ycenter)**yorder / sigmaN )
        # beyond
        else: return 0.0
    if dim == "1d3v": return fx
    if dim == "2d3v": return lambda x,y: fx(x)*fy(y)

def polygonal(xpoints=[], xvalues=[]):
    global dim, sim_length
    if len(xpoints)!=len(xvalues):
        raise Exception("polygonal profile requires as many points as values")
    if len(sim_length)>0 and len(xpoints)==0:
        xpoints = [0., sim_length[0]]
        xvalues = [1., 1.]
    N = len(xpoints)
    def f(x,y=0.):
        # vacuum region
        if x < xpoints[0]: return 0.0;
        # polygon region (defined over N segments)
        elif x < xpoints[-1]:
            for i in range(1,len(xpoints)):
                if xpoints[i-1]==xpoints[i]: return xvalues[i-1]
                if x>=xpoints[i-1] and x<xpoints[i]:
                    m = (xvalues[i]-xvalues[i-1])/(xpoints[i]-xpoints[i-1])
                    return xvalues[i-1] + m * ( x-xpoints[i-1] )
        # beyond
        else: return 0.
    return f

def cosine(base, amplitude=1.,
           xvacuum=0., xlength=None, phi=0., xnumber=1):
    import math
    global sim_length
    if len(sim_length)>0 and xlength is None: xlength = sim_length[0]-xvacuum
    def f(x,y=0.):
        #vacuum region
        if x < xvacuum: return 0.
        # profile region
        elif x < xvacuum+xlength:
            return base + amplitude * math.cos(phi + 2.*math.pi * xnumber * (x-xvacuum)/xlength)
        # beyond
        else: return 0.
    return f
