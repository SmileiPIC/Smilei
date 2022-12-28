import math as m
import numpy as np

def arrondiChiffre(chiffre, nombreChiffreApresLaVirgule):
    LeChiffreArrondi = chiffre*10**nombreChiffreApresLaVirgule
    LeChiffreArrondi = int(LeChiffreArrondi)
    LeChiffreArrondi = float(LeChiffreArrondi)/float(10**nombreChiffreApresLaVirgule)
    return LeChiffreArrondi
#def my_filter(particles):
#    global vb, VT
#    randomNumber = np.random(len(particules.px), size=NtrajectPart)
#    return(particles[randomNumber].px>vb+VT)*(particles[randomNumber].px<vb-VT)

##############################################
#######
#######  personnal values
#######
##############################################
fluct = 'OFF'
mag = 'OFF'
PROBE = 'ON' # utilise les sorties paralleles avec PROBES

TeV = 200                      # electron & ion temperature in eV
Te = TeV/511.e3                # electron & ion temperature in me c^2

Tex = Te
Tey = Te
Tez = Te

TexBeam = 1*Tex 
TeyBeam = 1*Tey
TezBeam = 1*Tez

Tix = 0.1*Te
Tiy = 0.1*Te
Tiz = 0.1*Te

if (mag == 'ON'):
    bxo = 0.33 #0.35
    byo = 0
    bzo = 0

DNlevel=0.05   # average level of density fluctuations
n0  = 1                        #in unit of n_c
nb = 0.0005*n0                      #in unit of n_c


Lde = m.sqrt(Te/n0)                #-- Debye length in units of c/\omega_{pe}
VT =  m.sqrt(2*Te/n0)              #-- VT ain unit of 1/c
vb = 9*VT  

                        # -- in unit of 1/c

dx  =m.sqrt(2.)*Lde #0.5*Lde                  # cell length (same in x & y)
dy  = dx
dt  = 0.5 * dx/m.sqrt(2.)
#dt = dt/3.

number_part_per_cell = 1800
rapport_de_masse=1836.

NumberOfCell_inX=128*32 # = 4096
NumberOfCell_inY=128*16 # = 2048

OutputEveryNumberOfCellForFieldData=1 # on sauvegarde les point grille des fichiers grilles

SizePatch=128

#################################
NumberOfTimeStep=1000

Scalar_save = 'ON'
NumberOfTimeStepForSnapshotSCALAR = 100

Fields_save = 'OFF'
NumberOfTimeStepForSnapshotFIELD = 100 # sortie des fichiers champs tout les NumberOfTimeStepForSnapshotFIELD pas de temps
FLUSH_EVERY_SnapshotFIELD = 10*NumberOfTimeStepForSnapshotFIELD # flush les donnees tout FLUSH_EVERY

Particles_save = 'OFF'
NumberOfTimeStepForSnapshotPART = 3000 # -> # sortie des fichiers particules tout les NumberOfTimeStepForSnapshotPART pas de temps
NumberOfTimeStepForTrajectoriesPART=10 # sortie des trajectoires de NtrajectPart particules tout les NumberOfTimeStepForTrajectoriesPART pas de temps
NtrajectPart = 20000

Distribution_save = 'OFF'
NumberOfTimeStepForSnapshotDISTRI = 6000
# Velocity limit for the distribution histogram
Distri_VXmin = -20.*VT
Distri_VXmax = vb+20.*VT

Distri_VYmin = -20.*VT
Distri_VYmax = 20.*VT

Distri_VZmin = Distri_VYmin
Distri_VZmax = Distri_VYmax

Distri_VXminION = -0.06*VT
Distri_VXmaxION = 0.06*VT

Distri_VYminION = -0.08*VT
Distri_VYmaxION = 0.08*VT

Distri_VZminION = Distri_VYminION
Distri_VZmaxION = Distri_VYmaxION

# number of bins for the histogram
NbinsX = 800
NbinsY = 400
NbinsV = 600
PrintOUTInformationEvery=2

#################################

dt=arrondiChiffre(dt,3)# pour avoir des snapshots espaces correctement
                       # et toujours verifier la condition CFL

number_of_patchX=int(m.pow(2, int(m.log(NumberOfCell_inX/SizePatch)/m.log(2))))
number_of_patchY=int(m.pow(2, int(m.log(NumberOfCell_inY/SizePatch)/m.log(2))))

NumberOfCell_inX=int(number_of_patchX*SizePatch)
NumberOfCell_inY=int(number_of_patchY*SizePatch)

Lx = NumberOfCell_inX*dx
Ly = NumberOfCell_inY*dy


#DNlevel=0.05   # average level of fluctuations
SeedRandom=1
lambdaX=12      # average wavelength of fluctuations along x
lambdaY=12      # average wavelength of fluctuations along y

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
    
    np.random.seed(seed=SeedRandom)
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
    # v0 = -(nb/(n0+DNlevel-nb))*vb
    v0 = -(nb/(n0-nb))*vb
else:
    v0 = -(nb/(n0-nb))*vb

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
    print_every=PrintOUTInformationEvery,
    simulation_time = NumberOfTimeStep*dt,
    EM_boundary_conditions = [["periodic"],["periodic"]],
    random_seed = 0,
    gpu_computing=True,
)

# # Disabled for GPUs as of 2022/06
# LoadBalancing(
#     initial_balance = True,
#     every = 150,
#     cell_load = 1.,
#     frozen_particle_load = 0.1,
# )

# Always disabled for GPU (counter intuitive)
Vectorization(
   mode = "off"
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
    position_initialization = "ion",
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
    position_initialization = "ion",
    momentum_initialization = "maxwell-juettner",
    mean_velocity = [vb,0.,0.],
    temperature = [TexBeam, TeyBeam, TezBeam],
    particles_per_cell = number_part_per_cell,
    mass = 1,
    charge = -1.,
    number_density = nb,
    boundary_conditions = [["periodic", "periodic"],["periodic", "periodic"]]
)

# EXTERNAL MAGNETIC FIELDS
if (mag == 'ON'):
    ExternalField(
        field = "Bx",
        profile = bxo
    )
    ExternalField(
        field = "By",
        profile = byo
    )
    ExternalField(
        field = "Bz",
        profile = bzo
    )

ApplyWhenFiltreIsON='ON' # toujours ON ici sert pour le decoupage (voir lecture_namelist.py)

# 'disable' to ... disable
# 'initial' pour un premier run
# 'restart' pour lire un fichier restart pour continuer

restart_run='disable'
# restart_run='initial'
# restart_run='restart'

chemin_restart='.'

if restart_run == 'initial':
    Checkpoints(
            dump_step = NumberOfTimeStep,
            exit_after_dump = True,
            keep_n_dumps = 2)
elif restart_run == 'restart':
    Checkpoints(
            restart_dir = chemin_restart,
            dump_step = NumberOfTimeStep,
            exit_after_dump = True,
            keep_n_dumps = 2)
else:
    print("Restart disabled !")

####### DIAGNOSTICS######
if (Scalar_save == 'ON'):
    DiagScalar(
        every = NumberOfTimeStepForSnapshotSCALAR ,
        vars = ["Utot","Uelm","Ukin","Uelm_Ex","Uelm_Ey","Uelm_Ez","Uelm_Bx_m","Uelm_By_m","Uelm_Bz_m", "Ukin_electron-beam","Ukin_electron", "Ukin_ion"],
        precision = 15
    )

if (Fields_save == 'ON'):
    if (PROBE == 'ON'):
        DiagProbe(
            every = NumberOfTimeStepForSnapshotFIELD,
            origin = [0., 0.],
	    flush_every = FLUSH_EVERY_SnapshotFIELD,
            corners = [ [(NumberOfCell_inX+1)*dx, 0 ], [0, (NumberOfCell_inY+1)*dy]],
            number = [(NumberOfCell_inX+1)/OutputEveryNumberOfCellForFieldData, (NumberOfCell_inY+1)/OutputEveryNumberOfCellForFieldData],
#            fields = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez',
#                      'Rho_electron', 'Rho_electron-beam', 'Rho_ion']
            fields = ['Bz', 'Ex', 'Ey',
                      'Rho_electron', 'Rho_electron-beam', 'Rho_ion']
        )
    
    else:
        DiagFields(
            every = NumberOfTimeStepForSnapshotFIELD,
            flush_every = NumberOfTimeStepForSnapshotFIELD,
            fields = ['Bz', 'Ex', 'Ey','Jx_electron-beam','Jy_electron-beam', 'Jz_electron-beam', 'Rho_electron', 'Rho_electron-beam', 'Rho_ion']
        )

if (Distribution_save == 'ON'):

# # ******* PLASMA
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["vx", Distri_VXminION, Distri_VXmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["vy", Distri_VYminION, Distri_VYmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["vz", Distri_VZminION, Distri_VZmaxION , NbinsV]
#                 ]
#         )

#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["y", 0., Ly, NbinsY]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["vx", Distri_VXminION, Distri_VXmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["vy", Distri_VYminION, Distri_VYmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["vz", Distri_VZminION, Distri_VZmaxION , NbinsV]
#                 ]
#         )

#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["y", 0., Ly, NbinsY],
#                 ["vx", Distri_VXminION, Distri_VXmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["y", 0., Ly, NbinsY],
#                 ["vy", Distri_VYminION, Distri_VYmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["y", 0., Ly, NbinsY],
#                 ["vz", Distri_VZminION, Distri_VZmaxION , NbinsV]
#                 ]
#         )

#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["vx", Distri_VXminION, Distri_VXmaxION , NbinsV],
#                 ["vy", Distri_VYminION, Distri_VYmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"],
#         axes = [
#                 ["vx", Distri_VXminION, Distri_VXmaxION , NbinsV],
#                 ["vz", Distri_VZminION, Distri_VZmaxION , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["ion"], 
#         axes = [
#                 ["vy", Distri_VYminION, Distri_VYmaxION , NbinsV],
#                 ["vz", Distri_VZminION, Distri_VZmaxION , NbinsV]
#                 ]
#         )


# # ******* PLASMA
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["y", 0., Ly, NbinsY]
#                 ]
#         )

#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["x", 0., Lx, NbinsX],
#                 ["vz", Distri_VZmin, Distri_VZmax , NbinsV]

#                 ]
#         )

#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["y", 0., Ly, NbinsY],
#                 ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["y", 0., Ly, NbinsY],
#                 ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["y", 0., Ly, NbinsY],
#                 ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
#                 ]
#         )

#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["vx", Distri_VXmin, Distri_VXmax , NbinsV],
#                 ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["vx", Distri_VXmin, Distri_VXmax , NbinsV],
#                 ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
#                 ]
#         )
#     DiagParticleBinning(
#         deposited_quantity = "weight",
#         every = NumberOfTimeStepForSnapshotDISTRI,
#         time_average = 1,
#         species = ["electron"],
#         axes = [
#                 ["vy", Distri_VYmin, Distri_VYmax , NbinsV],
#                 ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
#                 ] 
#         )



#    ******  PLASMA-BEAM

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
                ]
        )
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
                ]
        )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron-beam"],
    #     axes = [
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["x", 0., Lx, NbinsX],
                ["y", 0., Ly, NbinsY]
                ]
        )
	
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["x", 0., Lx, NbinsX],
                ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
                ]
        )
    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["x", 0., Lx, NbinsX],
                ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
                ]
        )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron-beam"],
    #     axes = [
    #             ["x", 0., Lx, NbinsX],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )

    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron-beam"],
    #     axes = [
    #             ["y", 0., Ly, NbinsY],
    #             ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron-beam"],
    #     axes = [
    #             ["y", 0., Ly, NbinsY],
    #             ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron-beam"],
    #     axes = [
    #             ["y", 0., Ly, NbinsY],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )

    DiagParticleBinning(
        deposited_quantity = "weight",
        every = NumberOfTimeStepForSnapshotDISTRI,
        time_average = 1,
        species = ["electron-beam"],
        axes = [
                ["vx", Distri_VXmin, Distri_VXmax , NbinsV],
                ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
                ]
        )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron-beam"],
    #     axes = [
    #             ["vx", Distri_VXmin, Distri_VXmax , NbinsV],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron-beam"],
    #     axes = [
    #             ["vy", Distri_VYmin, Distri_VYmax , NbinsV],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )


#    ******  PLASMA + PLASMA-BEAM

    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )

    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["x", 0., Lx, NbinsX],
    #             ["y", 0., Ly, NbinsY]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["x", 0., Lx, NbinsX],
    #             ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["x", 0., Lx, NbinsX],
    #             ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["x", 0., Lx, NbinsX],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )

    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["y", 0., Ly, NbinsY],
    #             ["vx", Distri_VXmin, Distri_VXmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["y", 0., Ly, NbinsY],
    #             ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["y", 0., Ly, NbinsY],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )

    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["vx", Distri_VXmin, Distri_VXmax , NbinsV],
    #             ["vy", Distri_VYmin, Distri_VYmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["vx", Distri_VXmin, Distri_VXmax , NbinsV],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )
    # DiagParticleBinning(
    #     deposited_quantity = "weight",
    #     every = NumberOfTimeStepForSnapshotDISTRI,
    #     time_average = 1,
    #     species = ["electron", "electron-beam"],
    #     axes = [
    #             ["vy", Distri_VYmin, Distri_VYmax , NbinsV],
    #             ["vz", Distri_VZmin, Distri_VZmax , NbinsV]
    #             ]
    #     )

    

if (Particles_save=='ON'):
    DiagTrackParticles(
        species = "electron-beam",
        every = NumberOfTimeStepForSnapshotPART,
        attributes = ["x", "y", "px", "py", "pz"]
        )
#    DiagTrackParticles(
#        species = "electron-beam",
#	filter = my_filter,
#	every = NumberOfTimeStepForTrajectoriesPART,
#        flush_every = NumberOfTimeStepForSnapshotPART,
#        attributes = ["x", "y", "px", "py", "pz"]
#	)	
