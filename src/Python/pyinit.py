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
    
    def __init__(self):
        """"default constructor"""

def get_smilei():
    """ get first Smilei object declared"""
    for key,value in globals().items() :
        if isinstance(value,Smilei) :
            return value
    return None

        
class SmileiComponent():
    """Species generic class"""
    def __init__(self, *args, **kwargs):
        self.mysim = get_smilei() # set mysim to the first Smilei object available
        for key in args:
            if isinstance(key,Smilei) :
                self.mysim=key
        if kwargs is not None: # add all kwargs as internal class variables
            for key, value in kwargs.iteritems():
                setattr(self, key, value)

class Species(SmileiComponent):
    """Species parameters"""
    species_type='None'
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        if isinstance(self.mysim,Smilei):
            self.mysim.species.append(self)
    
                    
class Laser(SmileiComponent):
    """Laser parameters"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        if isinstance(self.mysim,Smilei):
            self.mysim.laser.append(self)


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        if isinstance(self.mysim,Smilei):
            self.mysim.diag_probe.append(self)

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        if isinstance(self.mysim,Smilei):
            self.mysim.diag_particles.append(self)

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        if isinstance(self.mysim,Smilei):
            self.mysim.diag_phase.append(self)

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    def __init__(self, *args, **kwargs):
        SmileiComponent.__init__(self, *args, **kwargs)
        if isinstance(self.mysim,Smilei):
            self.mysim.diag_scalar.append(self)

    



