
execfile("../../scripts/ParticleDiagnostic.py")
from scipy.special import erf as erf

for path in ["thermalisation_ei1","thermalisation_ei2","thermalisation_ei3"]:

	mass_ion             = np.double(findParam(path, "mass"         ,"ion1"))
	charge_ion           = np.double(findParam(path, "charge"       ,"ion1"))
	density_ion          = np.double(findParam(path, "density"      ,"ion1"))/charge_ion
	temperature_ion      = np.double(findParam(path, "temperature"  ,"ion1"))
	velocity_electron    =           findParam(path, "mean_velocity","electron1")
	velocity_electron    = np.double(velocity_electron.split()[0])
	temperature_electron = np.double(findParam(path, "temperature"  ,"electron1"))
	coulomb_log          = np.double(findParam(path, "coulomb_log"  ,"ion1"))
	dt                   = np.double(findParam(path, "timestep"))
	
	re = 2.8179403267e-15 # meters
	wavelength = 1e-6 # meters
	c = 3e8
	coeff = (2.*np.pi/wavelength)**2*re*c / 8.
	
	times = getAvailableTimesteps(path, diagNumber=0)
	
	e_T_mean = np.zeros(len(times))
	i_T_mean = np.zeros(len(times))
	
	fig = None
	fig = plt.figure(1)
	if fig: fig.clf()
	if fig: ax = fig.add_subplot(1,1,1)
	for i,t in enumerate(times):
		electrons = ParticleDiagnostic(path,0, units="nice", slice={"x":"all"}, timesteps=t)
		vx = electrons["vx"]
		A = electrons["data"]
		vx0 = (A*vx).sum() / A.sum()
		e_T_mean[i] = (A*(vx-vx0)**2).sum() / A.sum()
		electrons = ParticleDiagnostic(path,1, units="nice", slice={"x":"all"}, timesteps=t)
		vy = electrons["vy"]
		A = electrons["data"]
		vy0 = (A*vy).sum() / A.sum()
		e_T_mean[i] += (A*(vy-vy0)**2).sum() / A.sum()
		electrons = ParticleDiagnostic(path,2, units="nice", slice={"x":"all"}, timesteps=t)
		vz = electrons["vz"]
		A = electrons["data"]
		vz0 = (A*vz).sum() / A.sum()
		e_T_mean[i] += (A*(vz-vz0)**2).sum() / A.sum()
		if fig:
			ax.cla()
			ax.plot(vx,A,'b')
		
		ions = ParticleDiagnostic(path,3, units="nice", slice={"x":"all"}, timesteps=t)
		vx = ions["vx"]
		A = ions["data"]
		vx0 = (A*vx).sum() / A.sum()
		i_T_mean[i] = (A*(vx-vx0)**2).sum() / A.sum()
		ions = ParticleDiagnostic(path,4, units="nice", slice={"x":"all"}, timesteps=t)
		vy = ions["vy"]
		A = ions["data"]
		vy0 = (A*vy).sum() / A.sum()
		i_T_mean[i] += (A*(vy-vy0)**2).sum() / A.sum()
		ions = ParticleDiagnostic(path,5, units="nice", slice={"x":"all"}, timesteps=t)
		vz = ions["vz"]
		A = ions["data"]
		vz0 = (A*vz).sum() / A.sum()
		i_T_mean[i] += (A*(vz-vz0)**2).sum() / A.sum()
		if fig:
			ax.cla()
			ax.plot(vx,A,'b')
	
	e_T_mean *= 1./3.
	i_T_mean *= 1./3.*mass_ion
	
	times *= 3.33*dt # fs
	
	# NRL relaxation	
	Te = temperature_electron
	Ti = temperature_ion
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


