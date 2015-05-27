"""@package pyinit
    Definition of Smilei components
"""

# \fixme Please explain what this is for
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)


class SmileiComponent(object):
    """Species generic class"""
    verify = True
    # This constructor is used always for all child classes
    def __init__(self, *args, **kwargs):
        if kwargs is not None: # add all kwargs as internal class variables
            for key, value in kwargs.iteritems():
                setattr(self, key, value)
        type(self).list.append(self) # add the current object to the static list "list"

class Species(SmileiComponent):
    """Species parameters"""
    species_type='None'
    list = []

class Laser(SmileiComponent):
    """Laser parameters"""
    list = []

class Collisions(SmileiComponent):
    """Collisions parameters"""
    list = []


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    list = []

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    list = []

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    list = []

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    list = []

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    list = []



