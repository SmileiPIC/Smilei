
import happi
from math import sqrt, pi
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf as erf

for path in ["thermalisation_ei1","thermalisation_ei2","thermalisation_ei3"]:

	sim = happi.Open(path)
	mass_ion             = np.double(sim.namelist.Species["ion1"].mass)
	charge_ion           = np.double(sim.namelist.Species["ion1"].charge)
	density_ion          = np.double(sim.namelist.Species["ion1"].charge_density)/charge_ion
	temperature_ion      = np.double(sim.namelist.Species["ion1"].temperature)
	velocity_electron    = np.double(sim.namelist.Species["electron1"].mean_velocity)[0]
	temperature_electron = np.double(sim.namelist.Species["electron1"].temperature)
	coulomb_log          = np.double(sim.namelist.Collisions[0].coulomb_log)
	dt                   = np.double(sim.namelist.Main.timestep)
	Wr = sim.namelist.Main.reference_angular_frequency_SI
	
	re_ = 2.8179403267e-15 # meters
	c = 3e8
	coeff = (2./3.) * sqrt(2./pi) * Wr**2 * re_/c
	
	times = np.double(sim.ParticleBinning(diagNumber=0).getAvailableTimesteps())
	
	electrons0 = sim.ParticleBinning(0).get()
	evx = electrons0["vx"]
	electrons1 = sim.ParticleBinning(1).get()
	evy = electrons1["vy"]
	electrons2 = sim.ParticleBinning(2).get()
	evz = electrons2["vz"]
	ions0 = sim.ParticleBinning(3).get()
	ivx = ions0["vx"]
	ions1 = sim.ParticleBinning(4).get()
	ivy = ions1["vy"]
	ions2 = sim.ParticleBinning(5).get()
	ivz = ions2["vz"]
	
	e_T_mean = np.zeros(len(times))
	i_T_mean = np.zeros(len(times))
	
	fig = None
	#fig = plt.figure(1)
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
	
	times *= (1e15/Wr)*dt # fs
	
	# NRL relaxation	
	Te = temperature_electron*1
	Ti = temperature_ion*1
	t = np.linspace(0,times.max(),3000)
	dt = np.diff(t).mean()
	Te_theory = t*0.
	Ti_theory = t*0.
	for i in range(len(t)):
		Te_theory[i] = Te
		Ti_theory[i] = Ti
		nu0 = coeff * charge_ion**2 * density_ion * coulomb_log  *np.sqrt(mass_ion)/((Te*mass_ion+Ti)**1.5)
		Te -= nu0*(Te-Ti)* dt*1e-15
		Ti += nu0*(Te-Ti)* dt*1e-15
	
	
	fig = plt.figure(2)
	#fig.clf()
	fig.set_facecolor('w')
	
	ax = fig.add_subplot(1,1,1)
	
	#ax.plot(times[:-1], np.diff(e_T_mean)/np.diff(times)*511, 'b', label='electrons')
	#ax.plot(times[:-1], np.diff(i_T_mean)/np.diff(times)*511, 'r', label='ions')
	#ax.plot(t[:-1], np.diff(Te_theory)/np.diff(t)*511,'k')
	#ax.plot(t[:-1], np.diff(Ti_theory)/np.diff(t)*511,'k')
	#ax.set_xlabel('time in fs')
	
	ax.plot(times, e_T_mean*511, 'b', label='electrons')
	ax.plot(times, i_T_mean*511, 'r', label='ions')
	ax.plot(t, Te_theory*511,'k')
	ax.plot(t, Ti_theory*511,'k')
	#ax.set_ylim(temperature_ion*511,temperature_electron*511)
	ax.set_xlabel('time in fs')
	ax.set_ylabel('Temperature [keV]')

	
	if not ax.legend_: ax.legend()
	
	plt.show()


