
"""
    END OF THE USER NAMELIST
-----------------------------------------------------------------------
"""

gc.collect()
import math

def _mkdir(role, path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except:
            raise Exception("ERROR in the namelist: "+role+" "+path+" cannot be created")
    elif not os.path.isdir(path):
        raise Exception("ERROR in the namelist: "+role+" "+path+" exists but is not a directory")

def _smilei_check():
    """Do checks over the script"""
    # Verify classes were not overriden
    for CheckClassName in ["SmileiComponent","Species", "Laser","Collisions",
            "DiagProbe","DiagParticles", "DiagScalar","DiagFields","ExtField",
            "SmileiSingleton","Main","DumpRestart","LoadBalancing","MovingWindow"]:
        CheckClass = globals()[CheckClassName]
        try:
            if not CheckClass._verify: raise Exception("")
        except:
            raise Exception("ERROR in the namelist: it seems that the name `"+CheckClassName+"` has been overriden")
    # Verify the output_dir
    if smilei_mpi_rank == 0 and Main.output_dir:
        _mkdir("output_dir", Main.output_dir)
    # Checkpoint: prepare dir tree
    if smilei_mpi_rank == 0 and (DumpRestart.dump_step>0 or DumpRestart.dump_minutes>0.):
        checkpoint_dir = (Main.output_dir or ".") + os.sep + "checkpoints" + os.sep
        if DumpRestart.file_grouping :
            ngroups = smilei_mpi_size/DumpRestart.file_grouping+1
            ngroups_chars = int(math.log10(ngroups))+1
            for group in range(ngroups):
                group_dir = checkpoint_dir + '%*s'%(ngroups_chars,group)
                _mkdir("checkpoint", group_dir)
        else:
            _mkdir("checkpoint", checkpoint_dir)
    # Checkpoint: Verify the restart_dir
    if len(DumpRestart)==1 and DumpRestart.restart_dir:
        if not os.path.isdir(DumpRestart.restart_dir):
            raise Exception("ERROR in the namelist: restart_dir = `"+DumpRestart.restart_dir+"` is not a directory")
    # Verify that constant() and tconstant() were not redefined
    if not hasattr(constant, "_reserved") or not hasattr(tconstant, "_reserved"):
        raise Exception("Names `constant` and `tconstant` cannot be overriden")
    # Convert float profiles to constant() or tconstant()
    def toSpaceProfile(input):
        try   : return constant(input*1.)
        except: return input
    def toTimeProfile(input):
        try:
            input*1.
            return tconstant()
        except: return input
    for s in Species:
        s.nb_density      = toSpaceProfile(s.nb_density      )
        s.charge_density  = toSpaceProfile(s.charge_density  )
        s.n_part_per_cell = toSpaceProfile(s.n_part_per_cell )
        s.charge          = toSpaceProfile(s.charge          )
        s.mean_velocity   = [ toSpaceProfile(p) for p in s.mean_velocity ]
        s.temperature     = [ toSpaceProfile(p) for p in s.temperature   ]
    for e in ExtField:
        e.profile         = toSpaceProfile(e.profile         )
    for a in Antenna:
        a.space_profile   = toSpaceProfile(a.space_profile   )
        a.time_profile    = toTimeProfile (a.time_profile    )
    for l in Laser:
        l.chirp_profile   = toTimeProfile( l.chirp_profile )
        l.time_envelope   = toTimeProfile( l.time_envelope )
        l.space_envelope  = [ toSpaceProfile(p) for p in l.space_envelope ]
        l.phase           = [ toSpaceProfile(p) for p in l.phase          ]

# this function will be called after initialising the simulation, just before entering the time loop
# if it returns false, the code will call a Py_Finalize();
def _keep_python_running():
    ps = [[las.time_envelope, las.chirp_profile] for las in Laser]
    ps += [[ant.time_profile] for ant in Antenna]
    if len(MovingWindow)>0 or len(LoadBalancing)>0:
        ps += [[s.nb_density, s.charge_density, s.n_part_per_cell, s.charge] + s.mean_velocity + s.temperature for s in Species]
    profiles = []
    for p in ps: profiles += p
    for prof in profiles:
        if callable(prof) and not hasattr(prof,"profileName"):
            return True
    return False

# Prevent creating new components (by mistake)
def _noNewComponents(cls, *args, **kwargs):
    print("Please do not create a new "+cls.__name__)
    return None
SmileiComponent.__new__ = staticmethod(_noNewComponents)



