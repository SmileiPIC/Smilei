
execfile("../../scripts/ParticleDiagnostic.py")
from scipy.special import erf as erf

path = "Maxwellianization1"

density_electron    = np.double(findParam(path, "density","electron1"))
coulomb_log         = np.double(findParam(path, "coulomb_log"  ,"electron1"))
dt                  = np.double(findParam(path, "timestep"))

re = 2.8179403267e-15 # meters
wavelength = 1e-6 # meters
c = 3e8
coeff = (2.*np.pi/wavelength)**2*re*c / (2.*np.sqrt(np.pi))

times = getAvailableTimesteps(path, diagNumber=0)

fig = None
fig = plt.figure(1)
if fig: fig.clf()
if fig: ax = fig.add_subplot(1,1,1)
for i,t in enumerate(times):
	if i%3>0: continue
	electrons = ParticleDiagnostic(path,0, units="nice", slice={"x":"all"}, timesteps=t)
	vx = electrons["vx"]
	A = electrons["data"]
	if fig:
		#ax.cla()
		ax.plot(vx,A,'b')

	if fig:
		fig.canvas.draw()
		#plt.pause(0.1)
	

fig.set_facecolor('w')

ax.plot(vx, 0.92e24*np.exp(-(vx/0.0068)**2), 'r', linewidth=2)

ax.set_xlabel('$v_x / c$')
ax.set_title('Velocity distribution at different times')

fig.canvas.draw()
plt.show()


