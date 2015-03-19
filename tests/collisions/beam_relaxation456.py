
execfile("../../scripts/ParticleDiagnostic.py")
from scipy.special import erf as erf

for path in ["beam_relaxation4","beam_relaxation5","beam_relaxation6"]:
	
	mass_ion             = np.double(findParam(path, "mass"         ,"ion1"))
	charge_ion           = np.double(findParam(path, "charge"       ,"ion1"))
	density_ion          = np.double(findParam(path, "density"      ,"ion1"))
	temperature_ion      = np.double(findParam(path, "temperature"  ,"ion1"))
	velocity_electron    =           findParam(path, "mean_velocity","electron1")
	velocity_electron    = np.double(velocity_electron.split()[0])
	temperature_electron = np.double(findParam(path, "temperature"  ,"electron1"))
	coulomb_log          = np.double(findParam(path, "coulomb_log"  ,"ion1"))
	dt                   = np.double(findParam(path, "timestep"))
	
	re = 2.8179403267e-15 # meters
	wavelength = 1e-6 # meters
	c = 3e8
	coeff = (2.*np.pi/wavelength)**2*re*c
	
	times = getAvailableTimesteps(path, diagNumber=0)
	
	e_vx_mean = np.zeros(len(times))
	e_vperp2  = np.zeros(len(times))
	i_vx_mean = np.zeros(len(times))
	Ti        = np.zeros(len(times))
	
	fig = None
	#fig = plt.figure(1)
	if fig: fig.clf()
	if fig: ax = fig.add_subplot(1,1,1)
	for i,t in enumerate(times):
		electrons = ParticleDiagnostic(path,0, units="nice", slice={"x":"all"}, timesteps=t)
		vx = electrons["vx"]
		A = electrons["data"]
		e_vx_mean[i] = (A*vx).sum() / A.sum()
	
		if fig:
			ax.cla()
			ax.plot(vx,A,'b')
	
		electrons = ParticleDiagnostic(path,1, units="nice", slice={"x":"all"}, timesteps=t)
		vperp2 = electrons["vperp2"]
		A = electrons["data"]
		e_vperp2[i] = (A*vperp2).sum() / A.sum()
	
		ions = ParticleDiagnostic(path,2, units="nice", slice={"x":"all"}, timesteps=t)
		vx = ions["vx"]
		A = ions["data"]
		i_vx_mean[i] = (A*vx).sum() / A.sum()
		Ti[i] = (A*(vx-i_vx_mean[i])**2).sum() / A.sum() * mass_ion
		
		if fig: 
			ax.plot(vx,A,'r')
			ax.set_ylim(ymax=1e24)
			fig.canvas.draw()
	
	
	times *= 3.33*dt # fs
	
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
		nu0 = coeff * charge_ion**2 * density_ion * coulomb_log / (v-vi[i])**3
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
	ax.set_xlim(0.,0.3)
	ax.set_ylim(0.,velocity_electron)
	ax.set_xlabel('time in fs')
	ax.set_ylabel('$v_x / c$')
	if not ax.legend_: ax.legend()
	
	plt.show()
	
	
	