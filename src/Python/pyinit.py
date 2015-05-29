"""@package pyinit
    Definition of Smilei components
"""

# Since the pytohn interpreter grabs key keyboards,
# we have to filter the ctrl-c kill command:
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

class SmileiComponentType(type):
    """Metaclass to all Smilei components"""
    
    # Constructor of classes
    def __init__(self, name, bases, attrs):
        self.list = []
        self.verify = True
        self.current = 0
        # Run standard metaclass init
        super(SmileiComponentType, self).__init__(name, bases, attrs)
    
    # Functions to define the iterator
    def __iter__(self):
        return self
    def next(self):
        if self.current >= len(self.list):
            raise StopIteration
        self.current += 1
        return self.list[self.current - 1]
    
    # Function to return one given instance, for example DiagParticles[0]
    # Special case: species can also be indexed by their name: Species["ion1"]
    def __getitem__(self, key):
        if self.__name__ == "Species" and type(key) is str:
            for obj in self.list:
                if obj.species_type == key:
                    return obj
        else:
            return self.list[key]
    
    # Function to return the number of instances, for example len(Species)
    def __len__(self):
        return len(self.list)
    
    # Function to display the content of the component
    def __repr__(self):
        if len(self.list)==0:
            return "<Empty list of "+self.__name__+">"
        else:
            l = []
            for obj in self.list: l.append(str(obj))
            return "["+", ".join(l)+"]"


class SmileiComponent(object):
    """Smilei component generic class"""
    __metaclass__ = SmileiComponentType
    
    # This constructor is used always for all child classes
    def __init__(self, **kwargs):
        if kwargs is not None: # add all kwargs as internal class variables
            for key, value in kwargs.iteritems():
                setattr(self, key, value)
        type(self).list.append(self) # add the current object to the static list "list"


class Species(SmileiComponent):
    """Species parameters"""
    species_type='None'

class Laser(SmileiComponent):
    """Laser parameters"""
    pass

class Collisions(SmileiComponent):
    """Collisions parameters"""
    pass


#diagnostics
class DiagProbe(SmileiComponent):
    """Diagnostic probe"""
    pass

class DiagParticles(SmileiComponent):
    """Diagnostic particles"""
    pass

class DiagPhase(SmileiComponent):
    """Diagnostic phase"""
    pass

class DiagScalar(SmileiComponent):
    """Diagnostic scalar"""
    pass

# external fields
class ExtField(SmileiComponent):
    """External Field"""
    pass



