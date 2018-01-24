import happi
import math as m
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True
plt.matplotlib.rcParams.update({
       'font.family':'serif',
       'font.serif':'Times',
        'font.size':20
})
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['ytick.minor.size'] = 5


S  = happi.Open('/Users/mica/RESULTS/SMILEI/thermalPlasmaPxDrift/')
T  = S.namelist.Te
mu = 1./T
print mu
v0 = S.namelist.v0
g0 = 1./m.sqrt(1.-v0**2)

# read p-distribution fct
fp  = np.array(S.ParticleBinning(1).getData())[0]
p   = np.array(S.ParticleBinning(1).get()['px'])
print 'int over all px:', (p[1]-p[0])*np.sum(fp)


# compute theoretical distribution fct
fth  = np.zeros(p.shape)
fth = (g0*np.sqrt(1.+p**2)+T) * np.exp( -g0*mu* (np.sqrt(1.+p**2) - v0*p) )
itg = (p[1]-p[0])*np.sum(fth)
fth = fth/itg

# plot
plt.figure(1)
plt.semilogy(p,fp,'b')
plt.hold('on')
plt.semilogy(p,fth,'r')
plt.xlabel(r'$p_x$')
plt.ylabel(r'$f(p_x)$')
plt.xlim([p.min(),p.max()])
plt.ylim([fp.max()/1.0e4,fp.max()*10.])
plt.tight_layout()