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

S  = happi.Open('/Users/mica/RESULTS/SMILEI/thermalPlasmaNoDrift/')
T  = S.namelist.Te
mu = 1./T
v0 = S.namelist.v0
g0 = 1./m.sqrt(1.-v0**2)

# read gamma-distribution fct
fg  = np.array(S.ParticleBinning(0).getData())[0]
g   = np.array(S.ParticleBinning(0).get()['gamma'])
print 'int over all g:', (g[1]-g[0])*np.sum(fg)


# compute theoretical distribution fct
fth = np.zeros(g.shape)
fth = g*np.sqrt(g**2-1.)*np.exp(-mu*g)
itg = (g[1]-g[0])*np.sum(fth)
fth = fth/itg

# read p-distribution fct
fp  = np.array(S.ParticleBinning(1).getData())[0]
p   = np.array(S.ParticleBinning(1).get()['px'])
print 'int over all px:', (p[1]-p[0])*np.sum(fp)


# compute theoretical distribution fct
fpth  = np.zeros(p.shape)
fpth = (g0*np.sqrt(1.+p**2)+T) * np.exp( -g0*mu* (np.sqrt(1.+p**2) - v0*p) )
itg = (p[1]-p[0])*np.sum(fpth)
fpth = fpth/itg

# plot
plt.figure(1)

plt.subplot(211)
plt.semilogy(g,fg,'b')
plt.hold('on')
plt.semilogy(g,fth,'r')
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$f(\gamma)$')
plt.xlim([g.min(),g.max()])
plt.ylim([fg.max()/1.0e3,fg.max()*10.])

plt.subplot(212)
plt.semilogy(p,fp,'b')
plt.hold('on')
plt.semilogy(p,fpth,'r')
plt.xlabel(r'$p_x$')
plt.ylabel(r'$f(p_x)$')
plt.xlim([p.min(),p.max()])
plt.ylim([fp.max()/1.0e3,fp.max()*10.])

plt.tight_layout()