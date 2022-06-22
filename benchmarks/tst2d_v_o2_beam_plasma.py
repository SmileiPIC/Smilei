import math as m
import numpy as np

def arrondiChiffre(chiffre, nombreChiffreApresLaVirgule):
    LeChiffreArrondi = chiffre*10**nombreChiffreApresLaVirgule
    LeChiffreArrondi = int(LeChiffreArrondi)
    LeChiffreArrondi = float(LeChiffreArrondi)/float(10**nombreChiffreApresLaVirgule)
    return LeChiffreArrondi
##############################################
#######
#######  personnal values
#######
##############################################
fluct = 'OFF'

TeV = 200                      # electron & ion temperature in eV
Te = TeV/511.e3                # electron & ion temperature in me c^2

Tex=Te
Tey=Te
Tez=Te

Tix=0.1*Te
Tiy=0.1*Te
Tiz=0.1*Te

n0  = 1                         #in unit of n_c
nb = 0.0005*n0                      #in unit of n_c


Lde = m.sqrt(Te/n0)                #-- Debye length in units of c/\omega_{pe}
VT =  m.sqrt(2*Te/n0)               #-- VT ain unit of 1/c
vb = 9*VT                           # -- in unit of 1/c
v0 = -(nb/(n0-nb))*vb
dx  =m.sqrt(2)*Lde #0.5*Lde                  # cell length (same in x & y)
dy  = dx
dt  = 0.5 * dx/m.sqrt(2.)

number_part_per_cell = 180 # 1800
rapport_de_masse=1836.

NumberOfCell_inX=1024#
NumberOfCell_inY=1024
OutputEveryNumberOfCellForFieldData=1 # on sauvegarde les point grille des fichiers grilles

SizePatch=16 #8*16 #2*8

#################################
NumberOfTimeStep=60000

Scalar_save = True
NumberOfTimeStepForSnapshotSCALAR = 10 #-> 6000 sorties # sortie valeurs scalaire tout les NumberOfTimeStepForSnapshotSCALAR pas de temps

Fields_save = True
NumberOfTimeStepForSnapshotFIELD = 20 # sortie des fichiers champs tout les NumberOfTimeStepForSnapshotFIELD pas de temps

Particles_save = True
NumberOfTimeStepForSnapshotPART = 3000 # -> # sortie des fichiers particules tout les NumberOfTimeStepForSnapshotPART pas de temps

Distribution_save = True # -> # sortie des focntions de distribution tout les NumberOfTimeStepForSnapshotDISTRI pas de temps
NumberOfTimeStepForSnapshotDISTRI = 600
# Velocity limit for the distribution histogram
Distri_VXmin = -6.*VT
Distri_VXmax = vb+15.*VT

Distri_VYmin = -15.*VT
Distri_VYmax = 15.*VT

Distri_VZmin = -15.*VT
Distri_VZmax = 15.*VT
# number of bins for the histogram
Nbins = 400

NumberOfPrintOUTInformation=10# signifie que le code ecrit dans le .out 10 fois durant le run
#################################

dt=arrondiChiffre(dt,3)# pour avoir des snapshots espaces correctement
                       # et toujours verifier la condition CFL

number_of_patchX=int(m.pow(2, int(m.log(NumberOfCell_inX/SizePatch)/m.log(2))))
number_of_patchY=int(m.pow(2, int(m.log(NumberOfCell_inY/SizePatch)/m.log(2))))

NumberOfCell_inX=int(number_of_patchX*SizePatch)
NumberOfCell_inY=int(number_of_patchY*SizePatch)

Lx = NumberOfCell_inX*dx
Ly = NumberOfCell_inY*dy


DNlevel=0.#0.05   # average level of fluctuations
lambdaX=9      # average wavelength of fluctuations along x
lambdaY=9      # average wavelength of fluctuations along y

if (fluct=='ON'):

    #Limites du plot en (kx,ky)
#    kxmin=-0.15
#    kymin=-0.15
#    kxmax=0.15
#    kymax=0.15
    ##################################
    x=np.linspace(0,Lx,NumberOfCell_inX)
    y=np.linspace(0,Ly,NumberOfCell_inY)
    kx=2*np.pi*np.linspace(-NumberOfCell_inX/2,NumberOfCell_inX/2-1,NumberOfCell_inX)/Lx
    ky=2*np.pi*np.linspace(-NumberOfCell_inY/2,NumberOfCell_inY/2-1,NumberOfCell_inY)/Ly
 
    densikk=np.zeros((NumberOfCell_inX, NumberOfCell_inY))
    densi=np.zeros((NumberOfCell_inX, NumberOfCell_inY))
    nnx=Lx*kx/2/np.pi
    nny=Ly*ky/2/np.pi
    kappaX=2*np.pi/lambdaX
    kappaY=2*np.pi/lambdaY
    
    np.random.seed(seed=1)
    rm=np.random.rand(NumberOfCell_inX, NumberOfCell_inY) # Nx*Ny matrix with random numbers

    for ii in range(0,NumberOfCell_inX):
       for jj in range(0,NumberOfCell_inY):
           densikk[ii,jj]=kx[ii]**2/kappaX**2+ky[jj]**2/kappaY**2

    #Gaussian spectrum
    
    densikk=np.exp(-densikk)*np.exp(1j*2*np.pi*rm) # product term by term of matrices
    densi=np.real(np.fft.ifft2(np.fft.ifftshift(densikk)))
    DN=densi*densi
    DNl=np.sqrt(DN.mean())
    densi=densi*DNlevel/DNl
    densikk=np.fft.ifftshift(np.fft.fft2(densi))

    dn=np.fft.ifftshift(densi)


def n_ion(x,y):
    if (fluct=='ON'):
        #with open('/scratch/cnt0026/lpp0106/ggauthier/SmileiKtest_Occi/toto'+str(x)+'.txt', 'w') as f:
            #sys.stdout = f # Change the standard output to the file we created.
            #print(x,y, x/dx, y/dy, int(x/dx), int(y/dy))
            #print(dn[int(x/dx), int(y/dy)])
            #sys.stdout = original_stdout
        return n0+dn[int(x/dx),int(y/dy)]
    else:
        if (x<=Lx):
            return n0

def n_elec(x,y):
    if (fluct=='ON'):
        return n0+dn[int(x/dx),int(y/dy)]-nb
    else:
        if (x<=Lx):
            return n0-nb
#
###############################
######
######  DATA FOR SMILEI
######
###############################
# Tsim : enorme chiffre que l'on atteint jamais

Main(
    geometry = "2Dcartesian",
    interpolation_order=2,
    cell_length = [dx,dy],
    number_of_cells= [NumberOfCell_inX,NumberOfCell_inY],
    number_of_patches = [number_of_patchX,number_of_patchY],
    timestep = dt,
    #print_every=NumberOfTimeStep/NumberOfPrintOUTInformation,
    simulation_time = 100*dt,
    EM_boundary_conditions = [["periodic"],["periodic"]],
    random_seed = 0,
    print_every = 2
)

LoadBalancing(
    initial_balance = True,
    every = 150,
    cell_load = 1.,
    frozen_particle_load = 0.1,
)

Vectorization(
   mode = "on"
)

Species(
    name = "ion",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    mean_velocity = [0.,0.,0.],
    temperature = [Tix, Tiy, Tiz],
    particles_per_cell = number_part_per_cell,
    mass = rapport_de_masse,
    charge = 1.,
    number_density = n_ion,
    boundary_conditions = [["periodic", "periodic"],["periodic", "periodic"]]
)
Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    mean_velocity = [v0,0.,0.],
    temperature = [Tex, Tey, Tez],
    particles_per_cell = number_part_per_cell,
    mass = 1,
    charge = -1.,
    number_density = n_elec,
    boundary_conditions = [["periodic", "periodic"],["periodic", "periodic"]]
)
Species(
    name = "electron-beam",
    position_initialization = "random",
    momentum_initialization = "maxwell-juettner",
    mean_velocity = [vb,0.,0.],
    temperature = [Tex, Tey, Tez],
    particles_per_cell = number_part_per_cell,
    mass = 1,
    charge = -1.,
    number_density = nb,
    boundary_conditions = [["periodic", "periodic"],["periodic", "periodic"]]
)


ApplyWhenFiltreIsON='ON' # toujours ON ici sert pour le decoupage (voir lecture_namelist.py)

restart_run='initial'
# initial pour un premier run
# restart pour lire un fichier restart pour continuer
# all job have to be in $SCRATCHDIR directory (see below for the eact path)
chemin_restart='SmileiGG12v2_Occi-016'


if (restart_run == 'initial'):
    Checkpoints(
            dump_step = NumberOfTimeStep,
#            dump_minutes = 240.,
#            dump_deflate = 0,
            exit_after_dump = True,
            keep_n_dumps = 2,
)
if (restart_run == 'restart'):
    Checkpoints(
            restart_dir = chemin_restart,
            dump_step = NumberOfTimeStep,
            exit_after_dump = True,
            keep_n_dumps = 2,
)
####### DIAGNOSTICS######
if (Scalar_save):
    DiagScalar(
        every = NumberOfTimeStepForSnapshotSCALAR ,
        vars = ["Utot","Uelm","Ukin","Uelm_Ex","Uelm_Ey","Uelm_Ez","Uelm_Bx_m","Uelm_By_m","Uelm_Bz_m", "Ukin_electron-beam","Ukin_electron", "Ukin_ion"],
        precision = 10
    )

if (Fields_save):
    DiagFields(
        every = NumberOfTimeStepForSnapshotFIELD,
        flush_every = NumberOfTimeStepForSnapshotFIELD,
        fields = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Rho_electron', 'Rho_electron-beam', 'Rho_ion']
    )


if (Distribution_save):
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins]
                ]
        )
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron"],
        axes = [
                ["x", 0., Ly, Nbins],
                ["vx", Distri_VXmin, Distri_VXmax , Nbins]
                ]
        )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["x", 0., Lx, Nbins],
                ["vx", Distri_VXmin, Distri_VXmax , Nbins]
                ]
        )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins]
                ]
        )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron", "electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins]
                ]
        )
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron", "electron-beam"],
        axes = [
                ["x", 0., Lx, Nbins],
                ["vx", Distri_VXmin, Distri_VXmax , Nbins]
                ]
        )

#******

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins],
                ["vy", Distri_VYmin, Distri_VYmax , Nbins]
                ]
        )
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins],
                ["vz", Distri_VZmin, Distri_VZmax , Nbins]
                ]
        )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["vy", Distri_VYmin, Distri_VYmax , Nbins],
                ["vz", Distri_VZmin, Distri_VZmax , Nbins]
                ]
        )
	

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins],
                ["vy", Distri_VYmin, Distri_VYmax , Nbins]
                ]
        )
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins],
                ["vz", Distri_VZmin, Distri_VZmax , Nbins]
                ]
        )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron"],
        axes = [
                ["vy", Distri_VYmin, Distri_VYmax , Nbins],
                ["vz", Distri_VZmin, Distri_VZmax , Nbins]
                ]
        )
	

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron", "electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins],
                ["vy", Distri_VYmin, Distri_VYmax , Nbins]
                ]
        )
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron","electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , Nbins],
                ["vz", Distri_VZmin, Distri_VZmax , Nbins]
                ]
        )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron","electron-beam"],
        axes = [
                ["vy", Distri_VYmin, Distri_VYmax , Nbins],
                ["vz", Distri_VZmin, Distri_VZmax , Nbins]
                ]
        )


if (Particles_save):
    DiagTrackParticles(
        species = "electron-beam",
        every = NumberOfTimeStepForSnapshotPART,
        attributes = ["x", "y", "px", "py", "pz"]
        )
