execfile("../../scripts/Diagnostics.py")
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf as erf

path = "temperature_isotropization1"

sim = Smilei(path)
density_electron     = np.double(sim.namelist.Species["electron1"].density)
coulomb_log          = np.double(sim.namelist.Collisions[0].coulomb_log)
dt                   = np.double(sim.namelist.timestep)

re_ = 2.8179403267e-15 # meters
wavelength = 1e-6 # meters
c = 3e8
coeff = (2.*np.pi/wavelength)**2*re_*c / (2.*np.sqrt(np.pi))

times = sim.ParticleDiagnostic(diagNumber=0).getAvailableTimesteps()

e_Tpar  = np.zeros(len(times))
e_Tperp = np.zeros(len(times))

fig = None
#fig = plt.figure(1)
if fig: fig.clf()
if fig: ax = fig.add_subplot(1,1,1)
for i,t in enumerate(times):
	electrons = sim.ParticleDiagnostic(0, units="nice", slice={"x":"all"}, timesteps=t).get()
	vx = electrons["vx"]
	A = electrons["data"][0]
	vx0 = (A*vx).sum() / A.sum()
	e_Tpar[i] = (A*(vx-vx0)**2).sum() / A.sum()
	if fig:
		ax.cla()
		ax.plot(vx,A,'b')

	electrons = sim.ParticleDiagnostic(1, units="nice", slice={"x":"all"}, timesteps=t).get()
	vy = electrons["vy"]
	A = electrons["data"][0]
	vy0 = (A*vy).sum() / A.sum()
	e_Tperp[i] = (A*(vy-vy0)**2).sum() / A.sum()
	electrons = sim.ParticleDiagnostic(2, units="nice", slice={"x":"all"}, timesteps=t).get()
	vz = electrons["vz"]
	A = electrons["data"][0]
	vz0 = (A*vz).sum() / A.sum()
	e_Tperp[i] += (A*(vz-vz0)**2).sum() / A.sum()
	if fig:
		ax.plot(vz,A,'r')
		
		fig.canvas.draw()

e_Tperp *= 1./2.

times *= 3.33*dt # fs

# NRL relaxation	
Tpar  = 0.0003
Tperp = 0.0002
t = np.linspace(0,60,3000)
dt = np.diff(t).mean()
Tpar_theory  = t*0.
Tperp_theory = t*0.
for i in range(len(t)):
	Tpar_theory[i] = Tpar
	Tperp_theory[i] = Tperp
	A = Tperp/Tpar - 1.
	if A>0: break
	nu0 = coeff * density_electron * coulomb_log /(2./3.*Tpar)**1.5 * A**-2 *(
		-3. + (A+3.) * np.arctanh((-A)**0.5)/(-A)**0.5 )
	#print A, Tpar, Tperp, nu0
	Tpar  -= 2.*nu0*(Tpar-Tperp)* dt*1e-15
	Tperp +=    nu0*(Tpar-Tperp)* dt*1e-15


fig = plt.figure(2)
#fig.clf()
fig.set_facecolor('w')

ax = fig.add_subplot(1,1,1)
ax.plot(times, e_Tpar *511, 'b', label='longitudinal temperature')
ax.plot(times, e_Tperp*511, 'r', label='perpendicular temperature')
ax.plot(t, Tpar_theory*511,'k')
ax.plot(t, Tperp_theory*511,'k')
ax.set_xlim(0.,60.)
ax.set_xlabel('time in fs')
ax.set_ylabel('Temperature [keV]')
if not ax.legend_: ax.legend()

plt.show()


