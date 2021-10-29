from math import pi, cos, sin, tan, asin, atan

l0 = 2.0*pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [10.*l0,8.*l0,5.*l0]  # length of the simulation
Tsim = 18.*t0                 # duration of the simulation
resx = 12.                    # nb of cells in one laser wavelength
rest = 22.                    # nb of timesteps in one optical cycle 

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4,4,4 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    
)

ang = [-pi/7., pi/6.]

LaserGaussian3D(
    box_side        = "xmin",
    a0              = 1.,
    omega           = 1.,
    focus           = [0.7*Lsim[0], 0.5*Lsim[1], 0.6*Lsim[2]],
    waist           = l0,
    incidence_angle = ang,
    time_envelope   = tgaussian(fwhm=t0*6)
)

LaserGaussian3D(
    box_side        = "ymin",
    a0              = 1.,
    omega           = 1.,
    focus           = [0.7*Lsim[0], 0.5*Lsim[1], 0.6*Lsim[2]],
    waist           = l0,
    incidence_angle = [ang[0], pi/2.-ang[1]],
    time_envelope   = tgaussian(fwhm=t0*6)
)

LaserGaussian3D(
    box_side        = "zmin",
    a0              = 1.,
    omega           = 1.,
    focus           = [0.7*Lsim[0], 0.5*Lsim[1], 0.6*Lsim[2]],
    waist           = l0,
    incidence_angle = [-asin(cos(ang[0])*sin(ang[1])), -atan(cos(ang[1])/tan(ang[0]))],
    time_envelope   = tgaussian(fwhm=t0*6),
    polarization_phi = pi/2.,
)


globalEvery = int(rest*4)

DiagScalar(
    every=globalEvery
)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez','Bz']
)
