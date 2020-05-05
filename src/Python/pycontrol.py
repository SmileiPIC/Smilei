
"""
    END OF THE USER NAMELIST
-----------------------------------------------------------------------
"""

gc.collect()
import math
import glob, re

def _mkdir(role, path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except:
            raise Exception("ERROR in the namelist: "+role+" "+path+" cannot be created")
    elif not os.path.isdir(path):
        raise Exception("ERROR in the namelist: "+role+" "+path+" exists but is not a directory")

def _prepare_checkpoint_dir():
    # Checkpoint: prepare dir tree
    if smilei_mpi_rank == 0 and (Checkpoints.dump_step>0 or Checkpoints.dump_minutes>0.):
        checkpoint_dir = "." + os.sep + "checkpoints" + os.sep
        if Checkpoints.file_grouping:
            ngroups = int((smilei_mpi_size-1)/Checkpoints.file_grouping + 1)
            ngroups_chars = int(math.log10(ngroups))+1
            for group in range(ngroups):
                group_dir = checkpoint_dir + '%0*d'%(ngroups_chars,group)
                _mkdir("checkpoint", group_dir)
        else:
            _mkdir("checkpoint", checkpoint_dir)

def _smilei_check():
    """Do checks over the script"""
    # Verify classes were not overriden
    for CheckClassName in ["SmileiComponent","Species", "Laser","Collisions",
            "DiagProbe","DiagParticleBinning", "DiagScalar","DiagFields",
            "DiagTrackParticles","DiagPerformances","ExternalField","PrescribedField",
            "SmileiSingleton","Main","Checkpoints","LoadBalancing","MovingWindow",
            "RadiationReaction", "ParticleData", "MultiphotonBreitWheeler",
            "Vectorization"]:
        CheckClass = globals()[CheckClassName]
        try:
            if not CheckClass._verify: raise Exception("")
        except:
            raise Exception("ERROR in the namelist: it seems that the name `"+CheckClassName+"` has been overriden")

    # Checkpoint: Verify the restart_dir and find possible restart file for each rank
    if len(Checkpoints)==1 and Checkpoints.restart_dir:
        if len(Checkpoints.restart_files) == 0 :
            Checkpoints.restart = True
            pattern = Checkpoints.restart_dir + os.sep + "checkpoints" + os.sep
            if Checkpoints.file_grouping:
                pattern += "*"+ os.sep
            pattern += "dump-*-*.h5";
            # pick those file that match the mpi rank
            files = filter(lambda a: smilei_mpi_rank==int(re.search(r'dump-[0-9]*-([0-9]*).h5$',a).groups()[-1]), glob.glob(pattern))
            
            if Checkpoints.restart_number is not None:
                # pick those file that match the restart_number
                files = filter(lambda a: Checkpoints.restart_number==int(re.search(r'dump-([0-9]*)-[0-9]*.h5$',a).groups()[-1]), files)
            
            Checkpoints.restart_files = list(files)
            
            if len(Checkpoints.restart_files) == 0:
                raise Exception(
                "ERROR in the namelist: cannot find valid restart files for processor "+str(smilei_mpi_rank) +
                "\n\t\trestart_dir = '" + Checkpoints.restart_dir +
                "'\n\t\trestart_number = " + str(Checkpoints.restart_number) +
                "\n\t\tmatching pattern: '" + pattern + "'" )
            
        else :
            raise Exception("restart_dir and restart_files are both not empty")

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
        s.number_density      = toSpaceProfile(s.number_density)
        s.charge_density  = toSpaceProfile(s.charge_density)
        s.particles_per_cell = toSpaceProfile(s.particles_per_cell)
        s.charge          = toSpaceProfile(s.charge)
        s.mean_velocity   = [ toSpaceProfile(p) for p in s.mean_velocity ]
        s.temperature     = [ toSpaceProfile(p) for p in s.temperature   ]
    for e in ExternalField:
        e.profile         = toSpaceProfile(e.profile)
    for e in PrescribedField:
        e.profile         = toSpaceProfile(e.profile)
    for a in Antenna:
        a.space_profile   = toSpaceProfile(a.space_profile   )
        a.time_profile    = toTimeProfile (a.time_profile    )
    for l in Laser:
        l.chirp_profile   = toTimeProfile( l.chirp_profile )
        l.time_envelope   = toTimeProfile( l.time_envelope )
        l.space_envelope  = [ toSpaceProfile(p) for p in l.space_envelope ]
        l.phase           = [ toSpaceProfile(p) for p in l.phase          ]
    for s in ParticleInjector:
        s.number_density = toSpaceProfile(s.number_density)
        s.charge_density = toSpaceProfile(s.charge_density)
        s.time_envelope  = toTimeProfile(s.time_envelope)
        s.particles_per_cell = toSpaceProfile(s.particles_per_cell )
        s.mean_velocity   = [ toSpaceProfile(p) for p in s.mean_velocity ]
        s.temperature     = [ toSpaceProfile(p) for p in s.temperature   ]
# this function will be called after initialising the simulation, just before entering the time loop
# if it returns false, the code will call a Py_Finalize();
def _keep_python_running():
    # Verify all temporal profiles, and all profiles that depend on the moving window or on the load balancing
    profiles = []
    for las in Laser:
        profiles += [las.time_envelope]
        profiles += [las.chirp_profile]
        if type(las.space_time_profile) is list:
            profiles += las.space_time_profile
        if type(las.space_time_profile_AM) is list:
            profiles += las.space_time_profile_AM
        if hasattr(las, "_extra_envelope"):
            profiles += [las._extra_envelope]
    profiles += [ant.time_profile for ant in Antenna]
    if len(MovingWindow)>0 or len(LoadBalancing)>0:
        for s in Species:
            profiles += [s.number_density, s.charge_density, s.particles_per_cell, s.charge] + s.mean_velocity + s.temperature
    if len(MovingWindow)>0:
        for e in ExternalField:
            profiles += [e.profile]
    if len(LoadBalancing)>0:
        if (Main.uncoupled_grids):
            return True
    for e in PrescribedField:
        profiles += [e.profile]
    for s in ParticleInjector:
        profiles += [s.time_envelope, s.number_density, s.charge_density, s.particles_per_cell] + s.mean_velocity + s.temperature
    for prof in profiles:
        if callable(prof) and not hasattr(prof,"profileName"):
            return True
    # Verify the tracked species that require a particle selection
    for d in DiagTrackParticles:
        if d.filter is not None:
            return True
    # Verify the particle binning having a function for deposited_quantity or axis type
    for d in DiagParticleBinning._list + DiagScreen._list:
        if type(d.deposited_quantity) is not str:
            return True
        for ax in d.axes:
            if type(ax[0]) is not str:
                return True
    # Verify the radiation spectrum having a function for deposited_quantity or axis type
    for d in DiagRadiationSpectrum._list:
        if type(d.photon_energy_axis[0]) is not str:
            return True
        for ax in d.axes:
            if type(ax[0]) is not str:
                return True
    # Verify if ionization from rate is used
    for s in Species:
        if s.ionization_rate is not None:
            return True
    # else False
    return False

# Prevent creating new components (by mistake)
def _noNewComponents(cls, *args, **kwargs):
    print("Please do not create a new "+cls.__name__)
    return None
SmileiComponent.__new__ = staticmethod(_noNewComponents)
