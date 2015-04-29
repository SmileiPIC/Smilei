# this is useful to kill the running simulation
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


print "hello from pyinit.py"

class smilei():
    """smilei parameters"""
    length_conversion=1.0
    time_conversion=1.0
    species=[]
    laser=[]
    diag_probe=[]
    diag_particles=[]
    diag_phase=[]
    diag_scalar=[]
    
    def __init__(self):
        """"default constructor"""
        print "smilei constructor"
    def add_species(self,spec) :
        if isinstance(spec,species):
            self.species.append(spec)

class species():
    """species parameters"""
    species_type='None'
    initPosition_type=''
    def __init__(self, _parent_smilei=None):
        if isinstance(_parent_smilei,smilei):
            _parent_smilei.species.append(self)

class laser():
    """laser parameters"""
    def __init__(self, _parent_smilei=None):
        if isinstance(_parent_smilei,smilei):
            _parent_smilei.laser.append(self)


#diagnostics
class diag_probe():
    """diagnostic probe"""
    def __init__(self, _parent_smilei=None):
        if isinstance(_parent_smilei,smilei):
            _parent_smilei.diag_probe.append(self)

class diag_particles():
    """diagnostic particles"""
    def __init__(self, _parent_smilei=None):
        if isinstance(_parent_smilei,smilei):
            _parent_smilei.diag_particles.append(self)

class diag_phase():
    """diagnostic phase"""
    def __init__(self, _parent_smilei=None):
        if isinstance(_parent_smilei,smilei):
            _parent_smilei.diag_phase.append(self)

class diag_scalar():
    """diagnostic scalar"""
    def __init__(self, _parent_smilei=None):
        if isinstance(_parent_smilei,smilei):
            _parent_smilei.diag_scalar.append(self)

    

#useful functions
def get_smilei() :
    for i,j in globals().items():
        if isinstance(j,smilei) :
            return j
    return None


