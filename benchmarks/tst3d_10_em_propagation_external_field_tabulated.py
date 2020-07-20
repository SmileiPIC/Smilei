import math
from numpy import exp, sqrt, arctan, vectorize, real, sin, cos, arctan, zeros_like, arange, meshgrid, linspace
from scipy import integrate
from math import log
import h5py

#Simulation parameters
dx = 0.393*0.5
dtrans = 4.712
dt = 0.96 * dx
nx =  1000*2 
ntrans = 100*2 
Lx = nx * dx
Ltrans = ntrans * dtrans
npatch_x = 8
npatch_trans = 8
Nit =  50

#Laser parameters and functions
a0    = 0.01
laser_FWHM_E = 58.81
focus = [ Lx/2., Ltrans/2., Ltrans/2.]
laser_initial_position = Lx/2.
c_vacuum = 1.
waist       = 157.
omega       = 1.
Zr          = omega * waist**2/2.  # Rayleigh length

# time gaussian function
def time_gaussian(fwhm, center, order=2):
    import math
    sigma = (0.5*fwhm)**order/log(2.0)
    def f(t):
        return exp( -(t-center)**order / sigma )
    return f

time_envelope_t              = time_gaussian(center=laser_initial_position                  , fwhm=laser_FWHM_E)
time_envelope_t_plus_half_dt = time_gaussian(center=(laser_initial_position+c_vacuum*0.5*dt), fwhm=laser_FWHM_E)

# laser waist function
def w(x):
        w  = sqrt(1./(1.+   ( (x-focus[0])/Zr  )**2 ) )
        return w

# Define Ex as solution of div(E)=0 for linear polarization along Y.
def fullEx(x,y,z):
    integration_constant = 0.
    nxl, nyl, nzl = x.shape
    A = zeros_like(x)
    y2d = y[0,:,:]-focus[1]
    z2d = z[0,:,:]-focus[2]
    xpatch = x[:,0,0]
    rsquare2d = y2d**2+z2d**2
    
    def r2overRC(xp):
        if xp != 0. :
            return -0.5*rsquare2d/(xp + Zr**2/xp)
        else:
            return 0.
    def yoverRC(xp):
        if xp != 0. :
            return y2d/(xp + Zr**2/xp)
        else:
            return 0.
    def invWaist2(xp):
        return  (w(xp)/waist)**2
    def spatial_amplitude(xp):
        return   w(xp) * exp( -invWaist2(xp)*rsquare2d)
    def modified_phase(xp):
        return xp + r2overRC(xp) + arctan(xp/Zr)
    def integrand_hyperslab(xp):
        return a0*time_envelope_t(xp)*spatial_amplitude(xp)*( 2*y2d*invWaist2(xp)*sin(modified_phase(xp-focus[0])) + yoverRC(xp-focus[0]) * cos(modified_phase(xp-focus[0]) ))
    
    for ix in range(nxl-1):
        slab = integrand_hyperslab(xpatch[ix])*dx/2.
        A[ix,:,:] += slab
        A[ix+1,:,:] = A[ix,:,:] + slab
    A[ 0,:,:]  = 0 # integrate from 0 to 0 equals 0.
    A[-1,:,:] += integrand_hyperslab(xpatch[nxl-1])*dx/2. #Last term

    for iy in range(nyl):
        for iz in range(nzl):
            mean = A[:,iy,iz].mean()
            A[:,iy,iz] -= mean

    return A

# Define Ey with paraxial approximation for linear polarization along Y.
def fullEy(x,y,z):
    nxl, nyl, nzl = x.shape
    A = zeros_like(x)
    y2d = y[0,:,:]-focus[1]
    z2d = z[0,:,:]-focus[2]
    xpatch = x[:,0,0]
    rsquare2d = y2d**2+z2d**2
    
    def r2overRC(xp):
        if xp != 0. :
            return -0.5*rsquare2d/(xp + Zr**2/xp)
        else:
            return 0.
    def invWaist2(xp):
        return  (w(xp)/waist)**2
    def spatial_amplitude(xp):
        return   w(xp) * exp( -invWaist2(xp)*rsquare2d)
    def modified_phase(xp):
        return xp + r2overRC(xp) + arctan(xp/Zr)
    def integrand_hyperslab(xp):
        return a0*time_envelope_t(xp)*spatial_amplitude(xp)*( sin(modified_phase(xp-focus[0])) )
    
    for ix in range(nxl):
        slab = integrand_hyperslab(xpatch[ix])
        A[ix,:,:] += slab
    return A

# Define Bx as solution of div(B)=0. Identical to Ex but with a dt/2 shift in time and dependency along z.
def fullBx(x,y,z):
    integration_constant = 0.
    nxl, nyl, nzl = x.shape
    A = zeros_like(x)
    y2d = y[0,:,:]-focus[1]
    z2d = z[0,:,:]-focus[2]
    xpatch = x[:,0,0]
    rsquare2d = y2d**2+z2d**2
    
    def r2overRC(xp):
        if xp != 0. :
            return -0.5*rsquare2d/(xp + Zr**2/xp)
        else:
            return 0.
    def zoverRC(xp):
        if xp != 0. :
            return z2d/(xp + Zr**2/xp)
        else:
            return 0.
    def invWaist2(xp):
        return  (w(xp)/waist)**2
    def spatial_amplitude(xp):
        return   w(xp) * exp( -invWaist2(xp)*rsquare2d)
    def modified_phase(xp):
        return xp -dt/2. + r2overRC(xp) + arctan(xp/Zr)
    def integrand_hyperslab(xp):
        return a0*time_envelope_t_plus_half_dt(xp)*spatial_amplitude(xp)*( 2*z2d*invWaist2(xp)*sin(modified_phase(xp-focus[0])) + zoverRC(xp-focus[0]) * cos(modified_phase(xp-focus[0]) ))
    
    for ix in range(nxl-1):
        slab = integrand_hyperslab(xpatch[ix])*dx/2.
        A[ix,:,:] += slab
        A[ix+1,:,:] = A[ix,:,:] + slab
    A[ 0,:,:]  = 0 # integrate from 0 to 0 equals 0.
    A[-1,:,:] += integrand_hyperslab(xpatch[nxl-1])*dx/2. #Last term

    for iy in range(nyl):
        for iz in range(nzl):
            mean = A[:,iy,iz].mean()
            A[:,iy,iz] -= mean
    return A

# Define Bz with paraxial approximation for linear polarization along Y. (identical to Ey with time shift of dt/2.)
def fullBz(x,y,z):
    nxl, nyl, nzl = x.shape
    A = zeros_like(x)
    y2d = y[0,:,:]-focus[1]
    z2d = z[0,:,:]-focus[2]
    xpatch = x[:,0,0]
    rsquare2d = y2d**2+z2d**2
    
    def r2overRC(xp):
        if xp != 0. :
            return -0.5*rsquare2d/(xp + Zr**2/xp)
        else:
            return 0.
    def invWaist2(xp):
        return  (w(xp)/waist)**2
    def spatial_amplitude(xp):
        return   w(xp) * exp( -invWaist2(xp)*rsquare2d)
    def modified_phase(xp):
        return xp - dt/2. + r2overRC(xp) + arctan(xp/Zr)
    def integrand_hyperslab(xp):
        return a0*time_envelope_t_plus_half_dt(xp)*spatial_amplitude(xp)*( sin(modified_phase(xp-focus[0])) )
    
    for ix in range(nxl):
        slab = integrand_hyperslab(xpatch[ix])
        A[ix,:,:] += slab
    return A

def tabulateEx():
    #dual, primal, primal
    xarray = linspace(-5*dx/2.,Lx+5*dx/2.,nx+2*2+2)
    yarray = linspace(-2.*dtrans,Ltrans+2.*dtrans,ntrans+2*2+1)
    zarray = linspace(-2.*dtrans,Ltrans+2.*dtrans,ntrans+2*2+1)
    X,Y,Z = meshgrid(xarray, yarray, zarray, indexing="ij")
    A = fullEx(X,Y,Z)
    h5f = h5py.File('tabulated_ex.h5', 'w')
    h5f.create_dataset('dataset_1', data=A)
    h5f.close()

def tabulateEy():
    #primal, dual, primal
    xarray = linspace(-2.*dx,Lx+2.*dx,nx+2*2+1)
    yarray = linspace(-5.*dtrans/2.,Ltrans+5.*dtrans/2.,ntrans+2*2+2)
    zarray = linspace(-2.*dtrans,Ltrans+2.*dtrans,ntrans+2*2+1)
    X,Y,Z = meshgrid(xarray, yarray, zarray, indexing="ij")
    A = fullEy(X,Y,Z)
    h5f = h5py.File('tabulated_ey.h5', 'w')
    h5f.create_dataset('dataset_1', data=A)
    h5f.close()

def tabulateBx():
#primal, dual, dual
    xarray = linspace(-2.*dx,Lx+2.*dx,nx+2*2+1)
    yarray = linspace(-5.*dtrans/2.,Ltrans+5.*dtrans/2.,ntrans+2*2+2)
    zarray = linspace(-5.*dtrans/2.,Ltrans+5.*dtrans/2.,ntrans+2*2+2)
    X,Y,Z = meshgrid(xarray, yarray, zarray, indexing="ij")
    A = fullBx(X,Y,Z)
    h5f = h5py.File('tabulated_bx.h5', 'w')
    h5f.create_dataset('dataset_1', data=A)
    h5f.close()

def tabulateBz():
#dual, dual, primal
    xarray = linspace(-5*dx/2.,Lx+5*dx/2.,nx+2*2+2)
    yarray = linspace(-5.*dtrans/2.,Ltrans+5.*dtrans/2.,ntrans+2*2+2)
    zarray = linspace(-2.*dtrans,Ltrans+2.*dtrans,ntrans+2*2+1)
    X,Y,Z = meshgrid(xarray, yarray, zarray, indexing="ij")
    A = fullBz(X,Y,Z)
    h5f = h5py.File('tabulated_bz.h5', 'w')
    h5f.create_dataset('dataset_1', data=A)
    h5f.close()

def preprocess():
    if smilei_mpi_size >= 4:
        if smilei_mpi_rank == 0:
            tabulateEx()
        if smilei_mpi_rank == 1:
            tabulateEy()
        if smilei_mpi_rank == 2:
            tabulateBx()
        if smilei_mpi_rank == 3:
            tabulateBz()
    else:
        if smilei_mpi_rank == 0:
            tabulateEx()
            tabulateEy()
            tabulateBx()
            tabulateBz()

def cleanup():
    if smilei_mpi_rank == 0:
        import os
        os.remove("tabulated_ex.h5")
        os.remove("tabulated_ey.h5")
        os.remove("tabulated_bx.h5")
        os.remove("tabulated_bz.h5")

Main(
    geometry = "3Dcartesian",
    interpolation_order = 2 ,
    cell_length  = [dx, dtrans, dtrans],
    grid_length = [ Lx,  Ltrans, Ltrans],
    number_of_patches = [ npatch_x, npatch_trans, npatch_trans ],
    timestep = dt,
    simulation_time = dt*Nit,
    clrw = 5,
    EM_boundary_conditions = [ ['silver-muller'] ],
    print_every = 1,
    random_seed = smilei_mpi_rank,
    save_magnectic_fields_for_SM = False,
)

MovingWindow(
    time_start = 0.,
    velocity_x = 1.0
)



################ Laser gaussian pulse, defined through external fields ###################

# Read tabulated Ex
def Ex(x,y,z):
    A = zeros_like(x)
    if x[0,0,0] < Lx-2.5*dx:
        nxl, nyl, nzl = x.shape
        ix = (x[0,0,0]+2.5*dx)/dx 
        iy = (y[0,0,0]+2.*dtrans)/dtrans 
        iz = (z[0,0,0]+2.*dtrans)/dtrans 
        ix = int(round(ix))
        iy = int(round(iy))
        iz = int(round(iz))
        h5f = h5py.File('tabulated_ex.h5','r') 
        A = h5f['dataset_1'][ix:ix+nxl,iy:iy+nyl,iz:iz+nzl].copy()
        h5f.close()
    return A

# Read tabulated Ey
def Ey(x,y,z):
    A = zeros_like(x)
    if x[0,0,0] < Lx-2*dx:
        ix = (x[0,0,0]+2*dx)/dx 
        iy = (y[0,0,0]+2.5*dtrans)/dtrans 
        iz = (z[0,0,0]+2*dtrans)/dtrans 
        nxl, nyl, nzl = x.shape
        ix = int(round(ix))
        iy = int(round(iy))
        iz = int(round(iz))
        h5f = h5py.File('tabulated_ey.h5','r') 
        A = h5f['dataset_1'][ix:ix+nxl,iy:iy+nyl,iz:iz+nzl].copy()
        h5f.close()
    return A

def Ez(x,y,z):
        return 0.*x

# Magnetic field
# Read tabulated Bx
def Bx(x,y,z):
    A = zeros_like(x)
    if x[0,0,0] < Lx-2.*dx:
        nxl, nyl, nzl = x.shape
        ix = (x[0,0,0]+2.*dx)/dx 
        iy = (y[0,0,0]+2.5*dtrans)/dtrans 
        iz = (z[0,0,0]+2.5*dtrans)/dtrans 
        ix = int(round(ix))
        iy = int(round(iy))
        iz = int(round(iz))
        h5f = h5py.File('tabulated_bx.h5','r') 
        A = h5f['dataset_1'][ix:ix+nxl,iy:iy+nyl,iz:iz+nzl].copy()
        h5f.close()
    return A

def By(x,y,z):
        return 0.*x

# Read tabulated Bz
def Bz(x,y,z):
    A = zeros_like(x)
    if x[0,0,0] < Lx-2.5*dx:
        nxl, nyl, nzl = x.shape
        ix = (x[0,0,0]+2.5*dx)/dx 
        iy = (y[0,0,0]+2.5*dtrans)/dtrans 
        iz = (z[0,0,0]+2.*dtrans   )/dtrans 
        ix = int(round(ix))
        iy = int(round(iy))
        iz = int(round(iz))
        h5f = h5py.File('tabulated_bz.h5','r') 
        A = h5f['dataset_1'][ix:ix+nxl,iy:iy+nyl,iz:iz+nzl].copy()
        h5f.close()
    return A

field_profile = {'Ex': Ex, 'Ey': Ey, 'Ez': Ez, 'Bx': Bx, 'By': By, 'Bz': Bz}

for field in ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']:
        ExternalField(
                field = field,
                profile = field_profile[field],
        )


##########################################################################################

globalEvery = int(50)

DiagScalar(
    every=globalEvery
)

#DiagFields(
#    every = globalEvery,
#    fields = ['Ex','Ey','Ez']
#)
#from numpy import s_
#DiagFields(
#    every = globalEvery,
#    fields = ['Ex','Ey','Ez'],
#    subgrid = s_[4:100:3, 5:400:10, 6:300:80]
#)

DiagProbe(
    every = 50,
    origin = [0.1*Main.grid_length[0], 0.5*Main.grid_length[1], 0.5*Main.grid_length[2]],
    fields = []
)

DiagProbe(
    every = 50,
    number = [nx],
    origin = [0.1*Main.grid_length[0], 0.5*Main.grid_length[1], 0.5*Main.grid_length[2]],
    corners = [[0.9*Main.grid_length[0], 0.5*Main.grid_length[1], 0.5*Main.grid_length[2]]],
    fields = []
)

DiagProbe(
    every = 50,
    number = [nx/2, ntrans/2],
    origin = [0.*Main.grid_length[0], 0.*Main.grid_length[1], 0.5*Main.grid_length[2]],
    corners = [
        [1.*Main.grid_length[0], 0. *Main.grid_length[1], 0.5*Main.grid_length[2]],
        [0.*Main.grid_length[0], 1.*Main.grid_length[1], 0.5*Main.grid_length[2]],
    ],
    fields = []
)

