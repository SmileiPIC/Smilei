# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

import numpy as np

l0 = 2.*np.pi               # laser wavelength
t0 = l0                     # optical cycle
Lsim = [20.*l0,20.*l0]      # length of the simulation
Tsim = 30.*t0               # duration of the simulation
resx = 8.                  # nb of cells in on laser wavelength
rest = 1.1*np.sqrt(2)*resx  # time of timestep in one optical cycle 

Main(
    geometry = "2Dcartesian",
    maxwell_solver = "Yee",
    interpolation_order = 2 ,
    cell_length = [l0/resx,l0/resx],
    grid_length  = Lsim,
    number_of_patches = [ 16, 16 ],
    timestep = t0/rest,
    simulation_time = Tsim,
    EM_boundary_conditions = [
                              ['PML','PML'],
                              ['PML','PML'],
                             ],
    number_of_pml_cells             = [[10,10],[10,10]],
)

Antenna(
    field='Jz',
    time_profile = lambda t: np.sin(2.*np.pi*t/t0)*np.sin(2.*np.pi*t/(4*t0))*(1.-np.heaviside(t-2*t0, 1)),
    space_profile=gaussian(1., xfwhm=l0, yfwhm=l0, xcenter=Main.grid_length[0]*0.5, ycenter=Main.grid_length[1]*0.5)
)

globalEvery=1.

DiagScalar(
    every=globalEvery,
    vars=['Uelm']
)

DiagFields(
    every = 100,
    fields = ['Ez'])
