
execfile("../../scripts/ParticleDiagnostic.py")
from scipy.special import erf as erf

path = "temperature_isotropization1"

density_electron    = np.double(findParam(path, "density","electron1"))
coulomb_log         = np.double(findParam(path, "coulomb_log"  ,"electron1"))
dt                  = np.double(findParam(path, "timestep"))

re = 2.8179403267e-15 # meters
wavelength = 1e-6 # meters
c = 3e8
coeff = (2.*np.pi/wavelength)**2*re*c / (2.*np.sqrt(np.pi))

times = getAvailableTimesteps(path, diagNumber=0)

e_Tpar  = np.zeros(len(times))
e_Tperp = np.zeros(len(times))

fig = None
#fig = plt.figure(1)
if fig: fig.clf()
if fig: ax = fig.add_subplot(1,1,1)
for i,t in enumerate(times):
	electrons = ParticleDiagnostic(path,0, units="nice", slice={"x":"all"}, timesteps=t)
	vx = electrons["vx"]
	A = electrons["data"]
	vx0 = (A*vx).sum() / A.sum()
	e_Tpar[i] = (A*(vx-vx0)**2).sum() / A.sum()
	if fig:
		ax.cla()
		ax.plot(vx,A,'b')

	electrons = ParticleDiagnostic(path,1, units="nice", slice={"x":"all"}, timesteps=t)
	vy = electrons["vy"]
	A = electrons["data"]
	vy0 = (A*vy).sum() / A.sum()
	e_Tperp[i] = (A*(vy-vy0)**2).sum() / A.sum()
	electrons = ParticleDiagnostic(path,2, units="nice", slice={"x":"all"}, timesteps=t)
	vz = electrons["vz"]
	A = electrons["data"]
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


