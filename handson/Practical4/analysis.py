#
#  ANALYSIS OF TUNNEL IONISATION SIMULATION (Practical 3)
#

simulation_to_analyse = 'to be defined'

# IMPORT SMILEI's POST-PROCESSING TOOL
# ------------------------------------

import happi


# IMPORT OTHER PYTHON PACKAGES
# ----------------------------

import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# DEFINE MATPLOTLIB PREFERENCES
# -----------------------------

### mpl.rcParams['text.usetex'] = True
plt.matplotlib.rcParams.update({
       'font.family':'serif',
       'font.serif':'Times',
        'font.size':16
})
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['ytick.minor.size'] = 5


# LOADING SIMULATION & IMPORTANT VARIABLES FROM NAMELIST
# ------------------------------------------------------

S  = happi.Open(simulation_to_analyse ,verbose=False)

l0   = S.namelist.l0
t0   = S.namelist.t0

tmax = S.namelist.tmax
Lv = S.namelist.Lv
Lp = S.namelist.Lp

rest = S.namelist.rest 
dt   = S.namelist.Main.timestep; nt = int(tmax/dt)
resx = S.namelist.resx

Lmu  = S.namelist.Lmu
I18  = S.namelist.I18
aL   = S.namelist.aL
Zat  = S.namelist.Species[0].atomic_number
reference_angular_frequency_SI = S.namelist.Main.reference_angular_frequency_SI

laser_time_envelope = S.namelist.Laser[0].time_envelope

print '- vector potential    aL = ',aL
print '- ref. ang. frequency w0 = ',reference_angular_frequency_SI



# SOLVE THE RATE EQUATION NUMERICALLY & PLOT THE RESULTS
# ------------------------------------------------------

execfile('solve_rate_eqs.py')

fig=plt.figure(1); 
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
plt.hold('on')
for Z in range(0,5):
    ax.plot(t/t0,n[Z,:],color='0.60',linewidth=1)
plt.fill_between(t/t0,Env, 0, interpolate=True, color='0.90')
plt.xlim(4,10)
plt.ylim(0,1.)
plt.show()



# SIMULATION ANALYSIS & COMPARISON TO RATE EQUATIONS
# --------------------------------------------------

# read n(Z,t): get the density of each charge state from the ParticleBinning diagnostics
n    = np.array( S.ParticleBinning(0).getData() )  
n00  = n[0,0]
n    = n/n00

# get corresponding time-steps
t    = dt*np.array(S.ParticleBinning(0).get()['times'])
t    = t - Lv - Lp # centering time axis 

# check conservation
nsum = 0
for Z in range(0,Zat+1):
    nsum += n[n.shape[0]-1,Z]
print '- checking conservation of the particle number: should give 1, returns:', nsum


# assign a color to each charge state
lcolor=['k--','r--','b--','g--','c--','m--','k--']

# plot the density of each charge state as a function of time
for Z in range(5):
    plt.plot(t/t0,n[:,Z],lcolor[Z],linewidth=2,label='Z = %d'%Z)
    ### plt.plot(t/t0,n[:,Z],lcolor[Z],linewidth=2,label=r'$Z^{\star}=$ %d'%Z)
    plt.legend(loc='upper left', frameon=False, borderpad=0.1, handletextpad=0.1, labelspacing=0.1)

plt.xlim(4,10)
plt.ylim(0,1.)
plt.xlabel('t [in optical cycles]') ### plt.xlabel(r'$c t/\lambda_0$')
plt.ylabel('N[Z]/N0') ### plt.ylabel(r'$N_{Z^{\star}}(t)/N_{0}(t=0) $')
plt.xticks( (4.,6.,8.,10.), color='k')
plt.yticks( (0.,0.25,0.5,0.75,1.), color='k')
ax.tick_params('both', length=4)
plt.savefig('Figure_Carbon_Ionization.eps')
plt.show()



