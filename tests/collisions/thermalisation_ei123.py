
execfile("../../scripts/Diagnostics.py")
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
	dt                   = np.double(sim.namelist.timestep)/(2*np.pi)
	
	re_ = 2.8179403267e-15 # meters
	wavelength = 1e-6 # meters
	c = 3e8
	coeff = (2.*np.pi/wavelength)**2*re_*c / 8.
	
	times = sim.ParticleDiagnostic(diagNumber=0).getAvailableTimesteps()
	
	e_T_mean = np.zeros(len(times))
	i_T_mean = np.zeros(len(times))
	
	fig = None
	fig = plt.figure(1)
	if fig: fig.clf()
	if fig: ax = fig.add_subplot(1,1,1)
	for i,t in enumerate(times):
		electrons = sim.ParticleDiagnostic(0, units="nice", slice={"x":"all"}, timesteps=t).get()
		vx = electrons["vx"]
		A = electrons["data"][0]
		vx0 = (A*vx).sum() / A.sum()
		e_T_mean[i] = (A*(vx-vx0)**2).sum() / A.sum()
		electrons = sim.ParticleDiagnostic(1, units="nice", slice={"x":"all"}, timesteps=t).get()
		vy = electrons["vy"]
		A = electrons["data"][0]
		vy0 = (A*vy).sum() / A.sum()
		e_T_mean[i] += (A*(vy-vy0)**2).sum() / A.sum()
		electrons = sim.ParticleDiagnostic(2, units="nice", slice={"x":"all"}, timesteps=t).get()
		vz = electrons["vz"]
		A = electrons["data"][0]
		vz0 = (A*vz).sum() / A.sum()
		e_T_mean[i] += (A*(vz-vz0)**2).sum() / A.sum()
		if fig:
			ax.cla()
			ax.plot(vx,A,'b')
		
		ions = sim.ParticleDiagnostic(3, units="nice", slice={"x":"all"}, timesteps=t).get()
		vx = ions["vx"]
		A = ions["data"][0]
		vx0 = (A*vx).sum() / A.sum()
		i_T_mean[i] = (A*(vx-vx0)**2).sum() / A.sum()
		ions = sim.ParticleDiagnostic(4, units="nice", slice={"x":"all"}, timesteps=t).get()
		vy = ions["vy"]
		A = ions["data"][0]
		vy0 = (A*vy).sum() / A.sum()
		i_T_mean[i] += (A*(vy-vy0)**2).sum() / A.sum()
		ions = sim.ParticleDiagnostic(5, units="nice", slice={"x":"all"}, timesteps=t).get()
		vz = ions["vz"]
		A = ions["data"][0]
		vz0 = (A*vz).sum() / A.sum()
		i_T_mean[i] += (A*(vz-vz0)**2).sum() / A.sum()
		if fig:
			ax.cla()
			ax.plot(vx,A,'b')
	
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


