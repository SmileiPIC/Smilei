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
    
    # check species for undefined/duplicate species_type
    all_species=[]
    for spec in Species:
        this_spec=spec.species_type
        try:
            if this_spec==None: raise
        except:
            raise Exception("ERROR in the namelist: there is a species without species_type")

        all_species.append(spec.species_type)
    try:
        if len(all_species)!=len(set(all_species)): raise
    except:
        raise Exception("ERROR in the namelist: there is duplicate species_type")
    
    
# this function will be called after initialising the simulation, just before entering the time loop
# if it returns false, the code will call a Py_Finalize();
def keep_python_running():
    retval=False
    for las in Laser:
        for prof in (las.time_profile, las.transv_profile):
            if callable(prof):
                retval=True
    return retval
        