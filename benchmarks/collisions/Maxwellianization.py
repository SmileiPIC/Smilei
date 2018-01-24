
import happi
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf as erf
plt.ion()

path = "Maxwellianization1"

sim = happi.Open(path)
coulomb_log = np.double(sim.namelist.Collisions[0].coulomb_log)
dt          = np.double(sim.namelist.Main.timestep)/(2*np.pi)

re_ = 2.8179403267e-15 # meters
wavelength = 1e-6 # meters
c = 3e8

times = sim.ParticleBinning(0).getAvailableTimesteps()
electrons = sim.ParticleBinning(0, sum={"x":"all"}).get()
vx = electrons["vx"]

fig = None
fig = plt.figure(1)
if fig: fig.clf()
if fig: ax = fig.add_subplot(1,1,1)
for i,t in enumerate(times):
	if i%3>0: continue
	A = electrons["data"][i]
	if fig:
		#ax.cla()
		ax.plot(vx,A,'b')
		plt.draw()
		plt.pause(0.0001)

fig.set_facecolor('w')

ax.plot(vx, 820.*np.exp(-(vx/0.0068)**2), 'r', linewidth=2)

ax.set_xlabel('$v_x / c$')
ax.set_title('Velocity distribution at different times')

fig.canvas.draw()
plt.show()


