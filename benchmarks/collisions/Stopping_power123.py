import happi
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import kv
from scipy.integrate import quad
plt.ion()

def tot(p,v1,a):
	g2 = (p**2+1.)**0.5
	v2 = p/g2
	g1 = (1.-v1**2)**-0.5
	Ap = g1*g2*(1.+v1*v2)
	Am = g1*g2*(1.-v1*v2)
	Bp = (Ap**2-1.)**0.5
	Bm = (Am**2-1.)**0.5
	F  = (g1-g2)*( -Ap/Bp + np.log(np.abs((Ap+Bp)/(Am+Bm))) + (Ap**2-2.)/Bp )
	return F/(v1*v2*g1**2*g2**2) * p**2*np.exp(-a*np.sqrt(1.+p**2))

def tot1(p,v1,a):
	g2 = (p**2+1.)**0.5
	v2 = p/g2
	g1 = (1.-v1**2)**-0.5
	Am = g1*g2*(1.-v1*v2)
	Bm = (Am**2-1.)**0.5
	F = (g1-g2)*(Am/Bm-(Am**2-2.)/Bm)
	return F/(v1*v2*g1**2*g2**2) * p**2*np.exp(-a*np.sqrt(1.+p**2))

def FrankelStoppingPower(E0,T):
	# Frankel (PRA 1979) - formula 5.10.
	# Only ee collisions.
	# E0 : initial electron energy (MeV)
	# T  : background plasma temperature (MeV)
	E0_ = np.double(E0)
	n = E0_.size
	v1 = np.double(np.sqrt(1.-(E0_/.511+1.)**-2))
	g1 = np.double((1.-v1**2)**(-0.5))
	a = .511/T
	y = np.zeros(n)
	for k in range(n):
		if E0_[k] < 0.1:
			pmin=g1[k]*v1[k]*(1.-g1[k]**2/4.)
			pmax=g1[k]*v1[k]*(1.+g1[k]**2/4.)
			dp = (pmax-pmin)/1000000.
			y[k] = quad(lambda p: tot1(p,v1[k],a), 0.   ,pmin   ,epsrel=3.e-14)[0]
			y[k]+= quad(lambda p: tot1(p,v1[k],a), pmin ,pmax-dp,epsrel=3.e-14)[0]
			y[k]+= quad(lambda p: tot1(p,v1[k],a), pmax ,np.inf ,epsrel=3.e-14)[0]
			y[k]+= quad(lambda p: tot (p,v1[k],a), 0.   ,np.inf ,epsrel=3.e-14)[0]
		else:
			y[k] = quad(lambda p: tot (p,v1[k],a), 0.   ,np.inf ,epsrel=3.e-14)[0]
			y[k]+= quad(lambda p: tot1(p,v1[k],a), 0.   ,np.inf ,epsrel=3.e-14)[0]
	y *= (a/(4.*np.pi*kv(2,a)))/v1
	y *= 3.204e-24 # 8*pi^2*me*c^2*re^2 in MeV*cm^2
	return y


for path in ["Stopping_power1","Stopping_power2","Stopping_power3"]:

	sim = happi.Open(path)
	temperature_electron = np.double(sim.namelist.Species["backgroundelectron"].temperature)
	density_electron     = np.double(sim.namelist.Species["backgroundelectron"].charge_density)
	coulomb_log          = np.double(sim.namelist.Collisions[0].coulomb_log)
	dt                   = np.double(sim.namelist.Main.timestep)/(2*np.pi)
	
	re = 2.8179403267e-15 # meters
	wavelength = 1e-6 # meters
	c = 3e8
	
	times = np.double(sim.ParticleBinning(diagNumber=0).getAvailableTimesteps())
	nx = sim.ParticleBinning(diagNumber=0,timesteps=0).get()["x"].size
	
	Ekin = np.zeros((nx,len(times)))
	electrons = sim.ParticleBinning(0).get()
	ekin = electrons["ekin"]*0.511
	dekin = np.diff(ekin)
	dekin = np.hstack((dekin, dekin[-1]))
	
	fig = None
	#fig = plt.figure(1)
	if fig: fig.clf()
	if fig: ax = fig.add_subplot(1,1,1)
	for i,t in enumerate(times):
		A = electrons["data"][i]
		A = A*dekin
		for k in range(nx): Ekin[k][i] = (A[k,:]*ekin).sum()/A[k,:].sum()
		
		if fig:
			ax.cla()
			ax.plot(ekin,A[0,:],'b')
			plt.draw()
	
	
	times *= 3.33*dt # fs
	
	fig = plt.figure(2)
	#fig.clf()
	fig.set_facecolor('w')
	ax = fig.add_subplot(1,1,1)
	for k in range(nx): ax.plot(times, Ekin[k]-Ekin[k][0])
	
	
	fig = plt.figure(3)
	fig.set_facecolor('w')
	ax = fig.add_subplot(1,1,1)
	for k in range(nx):
		Q = (Ekin[k][-1]-Ekin[k][1]) / (times[-1]-times[1])
		E0 = Ekin[k][0]
		v0 = np.sqrt(1.-(E0/0.511+1.)**-2)*c
		Q /= -1e-15*v0*coulomb_log # MeV/m
		Q *= 0.01 # MeV/cm
		Q /= density_electron*1.11e21 # MeV*cm^2
		ax.loglog(E0,Q, 'ok')
		
E = np.logspace(-2,1,1000)
ax.loglog(E,FrankelStoppingPower(E,temperature_electron*0.511), 'k', label='Frankel theory' )
if not ax.legend_: ax.legend()
ax.set_xlim(0.01,10)
ax.set_xlabel('Electron energy (MeV)')
ax.set_title('Stopping power $Q/(n_e \ln\Lambda)$ (MeV cm$^2$)')

plt.show()
