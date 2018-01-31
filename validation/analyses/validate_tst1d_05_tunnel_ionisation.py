import os, re, numpy as np, math
from scipy.interpolate import interp1d as interp
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
au_to_w0 = 4.134137172e+16 / w_r;
Ec_to_au = 3.314742578e-15 * w_r;
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

# hydrogen
charge = S.ParticleBinning.Diag0().get()
charge_distribution = np.array( charge["data"] )
charge_distribution /= charge_distribution[0,:].sum()
n1, n2 = charge_distribution.shape
mean_charge = (charge_distribution * np.outer(np.ones((n1,)), np.arange(n2))).sum(axis=1)
# # theory
# Ip = np.array([13.5984])/27.2114
# l  = np.array([0])
# t, Zs = calculate_ionization(Ip, l)
# times = charge["times"]*S.namelist.Main.timestep
# Zs_theory = interp(t, Zs) (times)
Validate("Hydrogen mean charge vs time", mean_charge, 0.1)

# carbon (does not work yet)
charge = S.ParticleBinning.Diag1().get()
charge_distribution = np.array( charge["data"] )
charge_distribution /= charge_distribution[0,:].sum()
n1, n2 = charge_distribution.shape
mean_charge = (charge_distribution * np.outer(np.ones((n1,)), np.arange(n2))).sum(axis=1)
# # theory
# Ip  = np.array([11.2602,24.3845,47.8877,64.4935,392.0905,489.9931]) /27.2114
# l   = np.array([1,1,0,0,0,0])
# t, Zs = calculate_ionization(Ip, l)
# times = charge["times"]*S.namelist.Main.timestep
# Zs_theory = interp(t, Zs) (times)
Validate("Carbon mean charge vs time", mean_charge, 0.1)

# SCALARS RELATED TO SPECIES
Validate("Scalar Dens_electron", S.Scalar.Dens_electron().getData(), 0.003)
Validate("Scalar Ntot_electron", S.Scalar.Ntot_electron().getData(), 100.)
Validate("Scalar Zavg_carbon"  , S.Scalar.Zavg_carbon  ().getData(), 0.2)




