import os, re, numpy as np, math
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter, maximum_filter1d
from h5py import File
from glob import glob
import happi

S = happi.Open(["./restart*"], verbose=False)



# COMPARE THE Ey FIELD
Ey = S.Field.Field0.Ey(timesteps=1000).getData()[0]
Validate("Ey field at iteration 1000", Ey, 0.001)

# COMPARE THE Ez FIELD
Ez = S.Field.Field0.Ez(timesteps=1000).getData()[0]
Validate("Ez field at iteration 1000", Ez, 1e-9)

# VERIFY THE IONIZATION RATE Vs THEORY
w_r = S.namelist.Main.reference_angular_frequency_SI
au_to_w0 = 4.134137172e+16 / w_r
Ec_to_au = 3.314742578e-15 * w_r
a0 = S.namelist.Laser[0].space_envelope[1]

def calculate_ionization(Ip, l):
	Zat = len(Ip)
	Z = np.arange(Zat)
	nstar = (Z+1.) * (2.*Ip)**-0.5
	alpha = 2.*nstar - 1.
	gc = np.array([math.gamma(c) for c in 2.*nstar])
	beta = 2.**(2.*nstar-1.)/(2.*nstar)/gc * (8.*l+4.) * Ip * au_to_w0
	gamma = 2.*(2.*Ip)**1.5
	
	tmax  = 0.3*2.*np.pi
	dt    = 2.*np.pi/200.
	prev_n = np.zeros((Zat+1,)); prev_n[0]=1.
	times = []
	Zstar = []
	for t in np.arange(0, tmax, dt):
		E = abs(a0 * (math.sin(t-0.05) if t>0.05 else 0.)) *Ec_to_au
		if E > 1e-5:
			n = np.zeros((Zat+1,))
			for z in range(0,Zat+1):
				Wp = 0.
				if z < Zat:
					deltap = gamma[z]/E
					if deltap>1e-18:
						Wp = beta[z] * deltap**alpha[z] * math.exp(-deltap/3.)
				n[z] = (1.-Wp*dt/2.)/(1.+Wp*dt/2.)*prev_n[z]
				if z > 0:
					deltam = gamma[z-1]/E
					if deltam>1e-18:
						Wm = beta[z-1] * deltam**alpha[z-1] * math.exp(-deltam/3.)
						n[z] += Wm*dt/(1.+Wm*dt/2.)*prev_n[z-1]
			prev_n = n
		times += [t]
		Zstar += [ np.sum(prev_n*np.arange(Zat+1))/np.sum(prev_n) ]
	return times, Zstar

for i, model in enumerate(("tunnel", "tunnel_full_PPT", "tunnel_TL", "tunnel_BSI")):
    # hydrogen
    charge_distribution = S.ParticleBinning("hydrogen_"+model).getData()
    charge_distribution /= charge_distribution[0].sum()
    n = charge_distribution[0].size
    mean_charge = [np.sum(d*np.arange(n)) for d in charge_distribution]
    # # theory
    # Ip = np.array([13.5984])/27.2114
    # l  = np.array([0])
    # t, Zs = calculate_ionization(Ip, l)
    # times = charge["times"]*S.namelist.Main.timestep
    # Zs_theory = interp1d(t, Zs) (times)
    Validate(model+": Hydrogen mean charge vs time", mean_charge, 0.000001)

    # carbon
    charge_distribution = S.ParticleBinning("carbon_"+model).getData()
    charge_distribution /= charge_distribution[0].sum()
    n = charge_distribution[0].size
    mean_charge = [np.sum(d*np.arange(n)) for d in charge_distribution]
    # # theory
    # Ip  = np.array([11.2602,24.3845,47.8877,64.4935,392.0905,489.9931]) /27.2114
    # l   = np.array([1,1,0,0,0,0])
    # t, Zs = calculate_ionization(Ip, l)
    # times = charge["times"]*S.namelist.Main.timestep
    # Zs_theory = interp1d(t, Zs) (times)
    Validate(model+": Carbon mean charge vs time", mean_charge, 0.0001)

    # SCALARS RELATED TO SPECIES
    Validate(model+": Scalar Dens_electron", S.Scalar("Dens_electron_"+model).getData(), 0.0000001)
    Validate(model+": Scalar Ntot_electron", S.Scalar("Ntot_electron_"+model).getData(), 1.)
    Validate(model+": Scalar Zavg_carbon"  , S.Scalar("Zavg_carbon_"+model).getData(), 0.0000001)

    # TRACKING DIAGNOSTIC
    d = S.TrackParticles("electron_"+model, axes=["Id","x","Wy"], timesteps=1000).getData()
    keep = d["Id"] > 0
    order = np.argsort(d["x"][keep])
    Validate(model+": Track electron x", d["x"][keep][order][::200], 1e-4)
    Validate(model+": Track electron Wy", gaussian_filter(maximum_filter1d(d["Wy"][keep][order],20),200)[::200], 1e-5)

    # NEW PARTICLES DIAGNOSTIC
    d = S.NewParticles("electron_"+model).get()
    t = d["t"]
    q = d["q"]
    Validate(model+": DiagNewParticles: number of particles", t.size, 5. )
    tavg = [np.mean(t[q==i]) for i in [0,1,2,3]]
    Validate(model+": DiagNewParticles: time vs ionization state", tavg, 0.01 )
