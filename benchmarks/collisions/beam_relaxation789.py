import happi
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf as erf

for path in ["beam_relaxation7","beam_relaxation8","beam_relaxation9"]:

	sim = happi.Open(path)
	mass_ion             = np.double(sim.namelist.Species["ion1"].mass)
	charge_ion           = np.double(sim.namelist.Species["ion1"].charge)
	density_ion          = np.double(sim.namelist.Species["ion1"].charge_density)
	temperature_ion      = np.double(sim.namelist.Species["ion1"].temperature)
	velocity_electron    = np.double(sim.namelist.Species["electron1"].mean_velocity)[0]
	temperature_electron = np.double(sim.namelist.Species["electron1"].temperature)
	coulomb_log          = np.double(sim.namelist.Collisions[0].coulomb_log)
	dt                   = np.double(sim.namelist.Main.timestep)/(2*np.pi)
	
	re = 2.8179403267e-15 # meters
	wavelength = 1e-6 # meters
	c = 3e8
	coeff = (2.*np.pi/wavelength)**2*re*c
	
	times = sim.ParticleBinning(diagNumber=0).getAvailableTimesteps()
	
	e_vx_mean = np.zeros(len(times))
	e_vperp2  = np.zeros(len(times))
	i_vx_mean = np.zeros(len(times))
	Ti        = np.zeros(len(times))
	
	electrons0 = sim.ParticleBinning(0, sum={"x":"all"}).get()
	electrons1 = sim.ParticleBinning(1, sum={"x":"all"}).get()
	ions       = sim.ParticleBinning(2, sum={"x":"all"}).get()
	
	fig = None
	#fig = plt.figure(1)
	if fig: fig.clf()
	if fig: ax = fig.add_subplot(1,1,1)
	for i,t in enumerate(times):
		vx = electrons0["vx"]
		A = electrons0["data"][i]
		e_vx_mean[i] = (A*vx).sum() / A.sum()
	
		if fig:
			ax.cla()
			ax.plot(vx,A,'b')
	
		vperp2 = electrons1["vperp2"]
		A = electrons1["data"][i]
		e_vperp2[i] = (A*vperp2).sum() / A.sum()
	
		vx = ions["vx"]
		A = ions["data"][i]
		i_vx_mean[i] = (A*vx).sum() / A.sum()
		Ti[i] = (A*(vx-i_vx_mean[i])**2).sum() / A.sum() * mass_ion
		
		if fig: 
			ax.plot(vx,A,'r')
			ax.set_ylim(ymax=1e24)
			fig.canvas.draw()
	
	
	times = times * 3.33*dt # fs
	
	# NRL relaxation	
	v = velocity_electron
	t = np.linspace(0,1,1000)
	dt = np.diff(t).mean()
	vi = np.interp(t, times, i_vx_mean)
	ti = np.interp(t, times, Ti)
	ve_theory = t*0.
	vperp2_theory = t*0.
	vp2 = 0.0006**2
	for i in range(len(t)):
		ve_theory[i] = v
		vperp2_theory[i] = vp2
		nu0 = coeff * charge_ion * density_ion * coulomb_log / (v-vi[i])**3
		x = mass_ion*(v-vi[i])**2/(2.*ti[i])
		phi = erf(np.sqrt(x))-2./np.sqrt(np.pi)*np.exp(-x)*np.sqrt(x)
		v -= nu0*(1.+1./mass_ion)*phi* dt*1e-15*v
		if v<=vi[i]: break
		vp2 += nu0*2.*(1-1./(2.*x))*phi* dt*1e-15*v**2
	fig = plt.figure(2)
	#fig.clf()
	fig.set_facecolor('w')
	
	ax = fig.add_subplot(1,1,1)
	ax.plot(times, e_vx_mean, 'b', label='electrons $v_x$')
	ax.plot(times, np.sqrt(e_vperp2),'--b', label='electrons $v_\perp$')
	ax.plot(times, i_vx_mean, 'r', label='ions')
	ax.plot(t, ve_theory, '-k')
	ax.plot(t, np.sqrt(vperp2_theory), '--k')
	#ax.plot(times, ve_theory, '--k')
	#ax.plot(times, vi_theory, '--k')
	ax.set_xlim(0.,0.06)
	ax.set_ylim(0.,velocity_electron)
	ax.set_xlabel('time in fs')
	ax.set_ylabel('$v_x / c$')
	if not ax.legend_: ax.legend(loc='lower right')
	
	plt.show()
	
	
	