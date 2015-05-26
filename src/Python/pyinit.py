# this is useful to kill the running simulation
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


class Smilei():
    """Smilei Main class"""
    def __init__(self, **kwargs):
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                setattr(self, key, value)
    species=[]
    laser=[]
    diag_probe=[]
    diag_particles=[]
    diag_phase=[]
    diag_scalar=[]
    
    # We define methods that create components for the current Smilei object
    def SmileiComponent(self, **kwargs): return SmileiComponent(self, **kwargs)
    def Species        (self, **kwargs): return Species        (self, **kwargs)
    def Laser          (self, **kwargs): return Laser          (self, **kwargs)
    def DiagProbe      (self, **kwargs): return DiagProbe      (self, **kwargs)
    def DiagParticles  (self, **kwargs): return DiagParticles  (self, **kwargs)
    def DiagPhase      (self, **kwargs): return DiagPhase      (self, **kwargs)
    def DiagScalar     (self, **kwargs): return DiagScalar     (self, **kwargs)

def get_smilei():
    """ get first Smilei object declared"""
    for key,value in globals().items() :
        if isinstance(value,Smilei) :
            return value
    return None


class SmileiComponent():
    """Species generic class"""
    def __init__(self, *args, **kwargs):
        # Try to find the Smilei object as an argument of the constructor
        self.mysim = None
        for arg in args:
            if isinstance(arg,Smilei):
                self.mysim=arg
        # If not found, get the first available Smilei object
        if self.mysim is None:
            self.mysim = get_smilei()
        # If still not found, raise error
        if self.mysim is None:
            raise Exception("Error in component `"+type(self).__name__+"`: cannot find associated simulation")
        # Collect all kwargs and adds them as internal class variables
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                setattr(self, key, value)

class Species(SmileiComponent):
    """Species parameters"""
    species_type='None'
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        self.mysim.species.append(self)

class Laser(SmileiComponent):
    """Laser parameters"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        self.mysim.laser.append(self)


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        self.mysim.diag_probe.append(self)

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        self.mysim.diag_particles.append(self)

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        self.mysim.diag_phase.append(self)

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        self.mysim.diag_scalar.append(self)




