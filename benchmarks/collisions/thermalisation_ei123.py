
from happi import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf as erf

for path in ["thermalisation_ei1","thermalisation_ei2","thermalisation_ei3"]:

	sim = Smilei(path)
	mass_ion             = np.double(sim.namelist.Species["ion1"].mass)
	charge_ion           = np.double(sim.namelist.Species["ion1"].charge)
	density_ion          = np.double(sim.namelist.Species["ion1"].charge_density)/charge_ion
	temperature_ion      = np.double(sim.namelist.Species["ion1"].temperature)
	velocity_electron    = np.double(sim.namelist.Species["electron1"].mean_velocity)[0]
	temperature_electron = np.double(sim.namelist.Species["electron1"].temperature)
	coulomb_log          = np.double(sim.namelist.Collisions[0].coulomb_log)
	dt                   = np.double(sim.namelist.Main.timestep)/(2*np.pi)
	
	re_ = 2.8179403267e-15 # meters
	wavelength = 1e-6 # meters
	c = 3e8
	coeff = (2.*np.pi/wavelength)**2*re_*c / 8.
	
	times = sim.ParticleBinning(diagNumber=0).getAvailableTimesteps()
	
	electrons0 = sim.ParticleBinning(0, sum={"x":"all"}).get()
	evx = electrons0["vx"]
	electrons1 = sim.ParticleBinning(1, sum={"x":"all"}).get()
	evy = electrons1["vy"]
	electrons2 = sim.ParticleBinning(2, sum={"x":"all"}).get()
	evz = electrons2["vz"]
	ions0 = sim.ParticleBinning(3, sum={"x":"all"}).get()
	ivx = ions0["vx"]
	ions1 = sim.ParticleBinning(4, sum={"x":"all"}).get()
	ivy = ions1["vy"]
	ions2 = sim.ParticleBinning(5, sum={"x":"all"}).get()
	ivz = ions2["vz"]
	
	e_T_mean = np.zeros(len(times))
	i_T_mean = np.zeros(len(times))
	
	fig = None
	fig = plt.figure(1)
	if fig: fig.clf()
	if fig: ax = fig.add_subplot(1,1,1)
	for i,t in enumerate(times):

		A = electrons0["data"][i]
		vx0 = (A*evx).sum() / A.sum()
		e_T_mean[i] = (A*(evx-vx0)**2).sum() / A.sum()
		
		A = electrons1["data"][i]
		vy0 = (A*evy).sum() / A.sum()
		e_T_mean[i] += (A*(evy-vy0)**2).sum() / A.sum()
		
		A = electrons2["data"][i]
		vz0 = (A*evz).sum() / A.sum()
		e_T_mean[i] += (A*(evz-vz0)**2).sum() / A.sum()
		if fig:
			ax.cla()
			ax.plot(evx,A,'b')
		
		A = ions0["data"][i]
		vx0 = (A*ivx).sum() / A.sum()
		i_T_mean[i] = (A*(ivx-vx0)**2).sum() / A.sum()
		A = ions1["data"][i]
		vy0 = (A*ivy).sum() / A.sum()
		i_T_mean[i] += (A*(ivy-vy0)**2).sum() / A.sum()
		A = ions2["data"][i]
		vz0 = (A*ivz).sum() / A.sum()
		i_T_mean[i] += (A*(ivz-vz0)**2).sum() / A.sum()
		if fig:
			ax.cla()
			ax.plot(ivx,A,'b')
	
	e_T_mean *= 1./3.
	i_T_mean *= 1./3.*mass_ion
	
	times *= 3.33*dt # fs
	
	# NRL relaxation	
	Te = temperature_electron*1
	Ti = temperature_ion*1
	t = np.linspace(0,300,3000)
	dt = np.diff(t).mean()
	Te_theory = t*0.
	Ti_theory = t*0.
	for i in range(len(t)):
		Te_theory[i] = Te
		Ti_theory[i] = Ti
		nu0 = coeff * charge_ion**2 * density_ion * coulomb_log  *np.sqrt(mass_ion)/((2./3.)**1.5*(Te+mass_ion*Ti)**1.5)
		Te -= nu0*(Te-Ti)* dt*1e-15
		Ti += nu0*(Te-Ti)* dt*1e-15
	
	
	fig = plt.figure(2)
	#fig.clf()
	fig.set_facecolor('w')
	
	ax = fig.add_subplot(1,1,1)
	ax.plot(times, e_T_mean*511, 'b', label='electrons')
	ax.plot(times, i_T_mean*511, 'r', label='ions')
	ax.plot(t, Te_theory*511,'k')
	ax.plot(t, Ti_theory*511,'k')
	ax.set_xlim(0.,300.)
	ax.set_ylim(temperature_ion*511,temperature_electron*511)
	ax.set_xlabel('time in fs')
	ax.set_ylabel('Temperature [keV]')
	if not ax.legend_: ax.legend()
	
	plt.show()


