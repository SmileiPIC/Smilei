
import happi
try:
	execfile("resparis.py")
except:
	exec(open("resparis.py").read(), globals(), locals())
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf as erf
plt.ion()

v0  = { "conductivity1":[-0.00064 ,-0.00041,-0.00025 ],
        "conductivity2":[-0.00088,-0.00108,-0.00145],
        "conductivity3":[-0.0023  ,-0.005            ]}

dv0 = { "conductivity1":[0.00008, 0.00002, 0        ],
        "conductivity2":[-0.0001,-0.00015,-0.0002  ],
        "conductivity3":[-0.0005,-0.0012           ]}

style = { "conductivity1": 'b', "conductivity2":'g', "conductivity3":'r' }

velocity = []
temperature = []
density = []

for path in ["conductivity1","conductivity2","conductivity3"]:

	S = happi.Open(path)

	ncases = 0
	for d in S.namelist.DiagParticleBinning:
		if d.deposited_quantity == "weight_charge_vx":
			ncases += 1
	if ncases == 0: continue
	
	print("simulation "+path+" has %d cases"%ncases)
	coulomb_log          = np.double(S.namelist.Collisions[0].coulomb_log)
	dt                   = np.double(S.namelist.Main.timestep)/(2*np.pi)
	
	times = np.double(S.ParticleBinning(0).getAvailableTimesteps())
	
	vx_mean = np.zeros((ncases,len(times)))
	
	fig = None
	#fig = plt.figure(1)
	if fig: fig.clf()
	if fig: ax = fig.add_subplot(1,1,1)
	for k in range(ncases):
		evx_density = -np.array(S.ParticleBinning(k*3).getData())
		edensity = np.array(S.ParticleBinning(k*3+1).getData())
		vx_mean[k,:] = evx_density/edensity
	
	
	times = times * 3.33*dt # fs
	
	fig = plt.figure(2)
	#fig.clf()
	fig.set_facecolor('w')
	
	ax = fig.add_subplot(1,1,1)
	
	
	for k in range(ncases):
		ax.plot(times, vx_mean[k,:], style[path], label='electrons $v_x$  #' + str(k))
		ax.plot(times, v0[path][k]+times*dv0[path][k], "--"+style[path])
		
		velocity.append(v0[path][k])
		temperature.append( np.double(S.namelist.Species["electron"+str(k+1)].temperature))
		density    .append( np.double(S.namelist.Species["electron"+str(k+1)].charge_density(20*(2*np.pi))))
	
	ax.set_xlabel('time in fs')
	ax.set_ylabel('$v_x / c$')
	#if not ax.legend_: ax.legend()
	
	plt.show()
	


velocity    = np.double(velocity)
temperature = np.double(temperature)
density     = np.double(density)
Ex          = np.array([0.002, 0.002, 0.002, 0.01, 0.01, 0.01, 0.01, 0.01])

# Simulation results
coeff = 1.66782e4 # e0*(2*pi*c/(1um)) in S/m
conductivity = (-density*velocity/Ex) * coeff

# Lee-More theory
t = np.logspace(-3,0,1000)
eta,teff,ze,xne,efev,aa,zlam=resparis(t,t,29,63.5,8.9,1,1,0.11,7.73,1,1,1,0,'-k',0)

# Perez formula
Te = 2.*t/511./np.pi
l  = (76.*1.11e27*(2.818e-15)**3 * np.sqrt(3./4./np.pi))**(-1./3.)
a  = 2.*ze*np.sqrt(coulomb_log/np.pi)/(l*Te)
perez = 2.367432e13 * ( # 8*pi*e0*c/re in S/m
	2.*ze/(3.*np.pi*l**2*np.sqrt(Te))*(1.-(1.+a)*np.exp(-a))
	+ Te**1.5/(ze*coulomb_log)*(1.+a+a**2/2.+a**3/6.)*np.exp(-a)
	)

fig = None
fig = plt.figure(3)
if fig:
	fig.clf()
	fig.set_facecolor('w')
	ax = fig.add_subplot(1,1,1)
	ax.loglog(temperature*511000. , conductivity, 'ok', label="SMILEI")
	ax.loglog(t*1000., 1./eta, '-k', label="Lee & More theory")
	ax.loglog(t*1000., perez , '--k', label="Formula from Ref [1]")
	
	ax.legend(loc='upper left')
	ax.set_xlim(1,1050)
	ax.set_ylim(4e5,4e7)
	
	ax.set_xlabel('Temperature [eV]')
	ax.set_ylabel('Conductivity [S/m]')
	
