
execfile("../../scripts/Diagnostics.py")
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf as erf

path = "Maxwellianization1"

sim = Smilei(path)
coulomb_log          = np.double(sim.namelist.Collisions[0].coulomb_log)
dt                   = np.double(sim.namelist.timestep)/(2*np.pi)


re_ = 2.8179403267e-15 # meters
wavelength = 1e-6 # meters
c = 3e8
coeff = (2.*np.pi/wavelength)**2*re_*c / (2.*np.sqrt(np.pi))

times = ParticleDiagnostic(path, diagNumber=0).getAvailableTimesteps()

fig = None
fig = plt.figure(1)
if fig: fig.clf()
if fig: ax = fig.add_subplot(1,1,1)
for i,t in enumerate(times):
	if i%3>0: continue
	electrons = ParticleDiagnostic(path,0, units="nice", slice={"x":"all"}, timesteps=t).get()
	vx = electrons["vx"]
	A = electrons["data"][0]
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


