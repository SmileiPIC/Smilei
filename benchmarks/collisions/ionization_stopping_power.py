import happi

import numpy as np
import matplotlib.pyplot as plt
ln = np.log

def RohrlichStoppingPower(Z, E): # MeV*cm^2
	# Formula of Rohrlich&Carlson (1954)
	#  taken from thesis of F Perez
	g = E/511.+1.
	v2 = 1.-g**-2
	I0 = ( 9.76*Z + 58.8*Z**-0.19 )/1000.
	return 2.55e-25*Z/(1.-g**-2) * (
		ln((g+1.)/2.*(E/I0)**2)
		+(1.-(2.*g-1.)*ln(2.) + (g-1.)**2/8.)/g**2  )



plt.figure(1,figsize=(8,3.5))

# First, slowing down for a few examples
########################################
D = []
sims = ["2","3","4"]
for sim in sims:
	S=happi.Open("ionization_stopping_power"+sim)
	
	D.append(S.ParticleBinning("#0/#1",sum={"x":"all"},units=["fs"],
		linestyle="None", marker='o', markersize=4, markeredgewidth=0.))
happi.multiPlot(*D, skipAnimation=True)
fig = plt.gcf()
ax = plt.gca()
ax.xaxis.labelpad = 0
ax.yaxis.labelpad = 0
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(10.,1e6)
ax.set_ylim(1e-2,3.)
ax.set_position([0.1, 0.19, 0.35,  0.77])

# Theory
for sim in sims:
	S=happi.Open("ionization_stopping_power"+sim)
	
	Zi = np.double(S.namelist.Species["ion1"].atomic_number)
	ni = np.double(S.namelist.Species["ion1"].number_density) * 1.11e21 # cm^-3
	Ee = S.namelist.E * 511. # electron energy in keV
	vel = S.namelist.vel * 3e10 # electron velocity in cm/s
	dt = 1e-15 * S.namelist.Main.timestep*4. # in s
	t = 0.
	energy = []
	time = []
	for k in range(30000):
		time.append(t)
		energy.append(Ee)
		Ee -= ni*vel*dt*RohrlichStoppingPower(Zi, Ee)*1e3
		if Ee<0.: break
		vel = np.sqrt(1.-1./(1.+Ee/511.)**2)*3e10
		t += dt
	momentum = np.sqrt((1.+np.array(energy)/511.)**2-1.)
	ax.plot(np.array(time)*1e15, momentum)

# Make nicer plot
plt.plot([0],"sk", markersize=4, markeredgewidth=0., label="Smilei")
plt.plot([0],"-k",label="Rohrlich & Carlson")
plt.legend(loc="best", prop={'size':10})
ax.set_ylabel("Momentum / ($m_e c$)")
ax.set_title("")



## Second, stopping power as a function of initial energy
#########################################################
S=happi.Open("ionization_stopping_power1")

reference_angular_frequency_SI = np.double(S.namelist.Main.reference_angular_frequency_SI)
timestep = np.double(S.namelist.Main.timestep)

# get ion stuff
ni = np.double(S.namelist.Species["ion1"].number_density)
mass = np.double(S.namelist.Species["ion1"].mass)*9.11e-28 # g
rho = mass*ni*1.1e21 # g/cm^3
Zi = np.double(S.namelist.Species["ion1"].atomic_number)


# Loop electron species
npoints = S.namelist.npoints
energy = []
Qsmilei = []
for i in range(npoints):
	print("electron"+str(i))
	# get electron velocity
	ve = np.double(S.namelist.Species["electron"+str(i)].mean_velocity[0])
	Ee = (1./np.sqrt(1.-ve**2)-1.)*511. # energy in keV
	energy.append(Ee)
	
	# Get data
	D = S.ParticleBinning("#"+str(2*i)+"/#"+str(2*i+1)+"",sum={"x":"all"})
	P = np.array(D.getData()) # momentum
	# stopping power in MeV/cm
	times = D.getTimesteps()
	diffP = -(P[-1]-P[0])/(times[-1]-times[0]) * (0.511*1e-2*reference_angular_frequency_SI/3e8/timestep)
	# stopping power in MeV*cm2/g
	Qsmilei.append( diffP/rho )
	
# Plot simulation result
ax = fig.add_axes([0.57, 0.19, 0.37,  0.77])
ax.loglog(energy, Qsmilei, 'ok', label="Smilei", markersize=4)
# Plot theory
E = np.logspace(-1,10,1000)
ax.loglog(E, RohrlichStoppingPower(Zi,E)/mass,'-k', label="Rohrlich & Carlson")
# Make nicer plot
ax.set_ylim(1.,3e2)
fig.set_facecolor("w")
plt.legend(loc="best", prop={'size':10})
ax.xaxis.labelpad = 0
ax.yaxis.labelpad = 0
ax.set_xlabel("Incident electron energy (keV)")
ax.set_ylabel("Stopping power (Mev cm$^2$/g)")
ax.set_xticks([1e-2, 1e0, 1e2, 1e4, 1e6, 1e8, 1e10])
