print "hello Python"

class smilei:
    """smilei parameters"""
    length_conversion=1.0
    time_conversion=1.0
    species=[]
    def __init__(self):
        """"default constructor"""
        print "constructor"
        

class species:
    """species parameters"""
    name='unknown'
    def __init__(self, _parent_smilei=None):
        if _parent_smilei:
            _parent_smilei.species.append(self)


def get_smilei() :
    for i,j in globals().items():
        if isinstance(j,smilei) :
            return j
    return None


def check_namelist() :
    if get_smilei() is None:
        raise NameError('No smilei defined') 