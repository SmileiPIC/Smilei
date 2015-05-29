"""@package pyinit
    Definition of Smilei components
"""

# This import allows for using ctrl+C during the simulation
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class Meta(type):
    """Metaclass to all Smilei components"""
    
    # This function occurs before the constructor of any Smilei component.
    # Sometimes, we don't want the class constructor to actually create a new instance.
    # Here we use __call__ to access to some Smilei components.
    # For instance, Species(2) returns the third existing Species instead of creating a new one.
    def __call__(cls, *args, **kwargs):
        # if no standalone arguments provided, then usual constructor is called
        if len(args)==0:
            self = cls.__new__(cls, **kwargs) # create object
            cls.__init__(self, **kwargs) # initialize object
            return self
        # if one standalone argument is provided, we return an existing instance
        elif len(args)==1:
            if len(kwargs)>0:
                print "ERROR: "+cls.__name__+"() cannot have both standalone and keyword arguments"
                return None
            if type(args[0]) is not int:
                print "ERROR: Argument to "+cls.__name__+"() should be an integer"
                return None
            if len(cls.list)<=args[0]:
                print "ERROR: There are only "+str(len(cls.list))+" "+cls.__name__+" objects, but requested #"+str(args[0])
                return None
            return cls.list[args[0]]
        # else, error
        else:
            print "ERROR: "+cls.__name__+"() accepts only one standalone argument"
            return None

class SmileiComponent(object):
    """Smilei component generic class"""
    __metaclass__ = Meta
    verify = True
    
    # This constructor is used always for all child classes
    def __init__(self, **kwargs):
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



