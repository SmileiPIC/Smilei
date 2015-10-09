"""@package pycontrol
    here we check if the namelist is clean and without errors
"""

import gc 
gc.collect()

import os

def _smilei_check():
    """Do checks over the script"""
    
    # Verify classes were not overriden
    for CheckClassName,CheckClass in {"SmileiComponent":SmileiComponent,"Species":Species,
            "Laser":Laser,"Collisions":Collisions,"DiagProbe":DiagProbe,"DiagParticles":DiagParticles,
            "DiagScalar":DiagScalar,"DiagPhase":DiagPhase,"ExtField":ExtField}.iteritems():
        try:
            if not CheckClass.verify: raise
        except:
            raise Exception("ERROR in the namelist: it seems that the name `"+CheckClassName+"` has been overriden")
    
    # Check species for undefined/duplicate type
    all_species=[]
    for spec in Species:
        if spec.species_type in all_species:
            raise Exception("ERROR in the namelist: there is duplicate species_type")
        else:
            all_species.append(spec.species_type)

    if dump_step and dump_dir:
        if not os.path.exists(dump_dir):
            try:
                os.makedirs(dump_dir)
            except OSError as exception:
                raise Exception("ERROR in the namelist: dump_dir "+dump_dir+" does not exists and cannot be created")
        elif not os.path.isdir(dump_dir):
                raise Exception("ERROR in the namelist: dump_dir "+dump_dir+" exists and is not a dir")

    if restart and restart_dir:
        if not os.path.isdir(restart_dir):
            raise Exception("ERROR in the namelist: restart_dir "+restart_dir+" is not a dir")



# this function will be called after initialising the simulation, just before entering the time loop
# if it returns false, the code will call a Py_Finalize();
def _keep_python_running():
    for las in Laser:
        for prof in (las.time_profile, las.transv_profile):
            if callable(prof): return True
    for ant in Antenna:
        if callable(ant.time_profile): return True

    if not nspace_win_x == 0:
        return True

    return False

# Prevent creating new components (by mistake)
def _noNewComponents(cls, *args, **kwargs):
    print "Please do not create a new "+cls.__name__
    return None
SmileiComponent.__new__ = staticmethod(_noNewComponents)



