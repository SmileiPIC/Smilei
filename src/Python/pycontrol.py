"""@package pycontrol
    here we check if the namelist is clean and without errors
"""

import gc 
gc.collect()

def smilei_check():
    """Do checks over the script"""
    
    # Verify classes were not overriden
    for CheckClassName,CheckClass in {"SmileiComponent":SmileiComponent,"Species":Species,
            "Laser":Laser,"Collisions":Collisions,"DiagProbe":DiagProbe,"DiagParticles":DiagParticles,
            "DiagScalar":DiagScalar,"DiagPhase":DiagPhase}.iteritems():
        try:
            if not CheckClass.verify: raise
        except:
            raise Exception("ERROR in the namelist: it seems that the name `"+CheckClassName+"` has been overriden")
    
    
smilei_check()


def keep_python_running():
    retval=True
    for las in Laser:
        for prof in (las.time_profile, las.transv_profile):
            if callable(prof):
                print "Don't stop me now! I'm having such a good time... "
                retval=False
    return retval
        