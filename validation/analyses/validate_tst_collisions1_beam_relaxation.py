import os, re, numpy as np
from scipy.special import erf
import happi

S = happi.Open(["./restart*"], verbose=False)



for i in range(3):
	ion = "ion"+str(i)
	eon = "eon"+str(i)
	ion_mass = S.namelist.Species[ion].mass
	
	times = np.double(S.ParticleBinning(0).getAvailableTimesteps())
	ones = np.ones_like(times)
	
	eon_vx = S.ParticleBinning(i*3+0, sum={"x":"all"}).get()
	eon_mean_vx = np.array(eon_vx["data"])
	eon_mean_vx = (np.outer(ones, eon_vx["vx"])*eon_mean_vx).sum(axis=1) / eon_mean_vx.sum(axis=1)
	
	eon_vp2 = S.ParticleBinning(i*3+1, sum={"x":"all"}).get()
	eon_mean_vp2 = np.array(eon_vp2["data"])
	eon_mean_vp2 = (np.outer(ones, eon_vp2["vperp2"])*eon_mean_vp2).sum(axis=1) / eon_mean_vp2.sum(axis=1)
	
	ion_vx = S.ParticleBinning(i*3+2, sum={"x":"all"}).get()
	ion_vxd = np.array(ion_vx["data"])
	ion_mean_vx = (np.outer(ones, ion_vx["vx"])*ion_vxd).sum(axis=1) / ion_vxd.sum(axis=1)
	ion_T = ( np.outer(ones, ion_vx["vx"]) - np.outer(ion_mean_vx, np.ones_like(ion_vx["vx"])) )**2 * ion_vxd
	ion_T = ion_T.sum(axis=1) / ion_vxd.sum(axis=1) * ion_mass
	
#	# NRL relaxation
#	ion_charge      = S.namelist.Species[ion].charge
#	ion_density     = S.namelist.Species[ion].number_density
#	eon_velocity    = S.namelist.Species[eon].mean_velocity[0]
#	eon_temperature = S.namelist.Species[eon].temperature[0]
#	coulomb_log     = S.namelist.Collisions[i].coulomb_log
#	W_r             = S.namelist.Main.reference_angular_frequency_SI
#	r_e = 2.8179403267e-15 # meters
#	c = 3e8
#	coeff = W_r**2*r_e/c * ion_charge * ion_density * coulomb_log
#	t = np.linspace(0., S.namelist.Main.simulation_time, 1000)
#	dt = np.diff(t).mean() / W_r
#	S_times = times * S.namelist.Main.timestep
#	vi = np.interp(t, S_times, ion_mean_vx)
#	ti = np.interp(t, S_times, ion_T)
#	ve_theory = t*np.nan
#	vperp2_theory = t*np.nan
#	vp2 = 0.
#	v = eon_velocity
#	for i in range(len(t)):
#		ve_theory[i] = v
#		vperp2_theory[i] = vp2
#		nu0 = coeff * (v-vi[i])**-3
#		x = ion_mass*(v-vi[i])**2/(2.*ti[i])
#		phi = erf(np.sqrt(x))-2.*np.exp(-x)*np.sqrt(x/np.pi)
#		v -= nu0*(1.+1./ion_mass)*phi* dt*v
#		if v**2<=vperp2_theory[i]*0.8: break
#		vp2 += nu0*2.*(1.-1./(2.*x))*phi* dt*v**2
#	
#	plt.plot(S_times/W_r*1e15, eon_mean_vx, 'b', label='electrons $v_x$')
#	plt.plot(S_times/W_r*1e15, np.sqrt(eon_mean_vp2),'--b', label='electrons $v_\perp$')
#	plt.plot(S_times/W_r*1e15, ion_mean_vx, 'r', label='ions')
#	plt.plot(t      /W_r*1e15, ve_theory, '-k')
#	plt.plot(t      /W_r*1e15, np.sqrt(vperp2_theory), '--k')
#	plt.show()
	
	Validate(eon+" mean vx", eon_mean_vx, 0.001)
	Validate(eon+" mean vperp", np.sqrt(eon_mean_vp2), 0.001)
	Validate(ion+" mean vx", ion_mean_vx, 0.001)

