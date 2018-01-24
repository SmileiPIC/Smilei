# _____________________________________________________
#
# Particle tracking analysis of trajectories in 
# a circular plane wave
#
# Non-relativistic (classical) regime
# _____________________________________________________

import happi
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
import math as m


mpl.rcParams['font.size'] = 20
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['figure.facecolor'] = 'white'

# ________________________________________________
# Functions



# ________________________________________________
# Open results

res = happi.Open("./tst1d_cir_plane_wave_clas",verbose=False)

# ________________________________________________
# Parameters

#a0 = res_Boris.namelist.LaserGaussian3D[0].a0
a0 = 0.1
step = 2  
 
# ____________________________________________
# Fields

if False:

  Ey = res_Boris.Field(0, "Ey", timesteps=1300, average = {"z":[5.]}).get()
 
  Ey["data"] = np.array(Ey["data"][0].T)
 
  fig = plt.figure(figsize=(12,6))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :]) 

  ax.plot(Ey["x"],Ey["data"])
 

# ____________________________________________
# Track

Track_Borisnr = res.TrackParticles("electron_borisnr", axes=["x","px","py","pz"]).get()
Track_Boris = res.TrackParticles("electron_norm", axes=["x","px","py","pz"]).get()
Track_Vay = res.TrackParticles("electron_vay", axes=["x","px","py","pz"]).get()
Track_HC = res.TrackParticles("electron_higueracary", axes=["x","px","py","pz"]).get()

#gamma_Boris = m.sqrt(1 + np.array(Track['px'][:,0])**2 + np.array(Track['py'][:,0])**2 + np.array(Track['px'][:,0])**2)

# Momentum

x_Borisnr = np.array(Track_Borisnr['x'][::step,0])
px_Borisnr = np.array(Track_Borisnr['px'][::step,0])
py_Borisnr = np.array(Track_Borisnr['py'][::step,0])
pz_Borisnr = np.array(Track_Borisnr['pz'][::step,0])
portho_Borisnr = np.sqrt(py_Borisnr**2 + pz_Borisnr**2)
gf_Borisnr = np.sqrt(1. + px_Borisnr**2 +  py_Borisnr**2 + pz_Borisnr**2)

x_Boris = np.array(Track_Boris['x'][::step,0])
px_Boris = np.array(Track_Boris['px'][::step,0])
py_Boris = np.array(Track_Boris['py'][::step,0])
pz_Boris = np.array(Track_Boris['pz'][::step,0])
portho_Boris = np.sqrt(py_Boris**2 + pz_Boris**2)
gf_Boris = np.sqrt(1. + px_Boris**2 +  py_Boris**2 + pz_Boris**2)

x_Vay = np.array(Track_Vay['x'][::step,0])
px_Vay = np.array(Track_Vay['px'][::step,0])
py_Vay = np.array(Track_Vay['py'][::step,0])
pz_Vay = np.array(Track_Vay['pz'][::step,0])
portho_Vay = np.sqrt(py_Vay**2 + pz_Vay**2)
gf_Vay = np.sqrt(1. + px_Vay**2 +  py_Vay**2 + pz_Vay**2)

x_HC = np.array(Track_HC['x'][::step,0])
px_HC = np.array(Track_HC['px'][::step,0])
py_HC = np.array(Track_HC['py'][::step,0])
pz_HC = np.array(Track_HC['pz'][::step,0])
portho_HC = np.sqrt(py_HC**2 + pz_HC**2)
gf_HC = np.sqrt(1. + px_HC**2 +  py_HC**2 + pz_HC**2)

# x - P_ortho
if True:

  fig = plt.figure(figsize=(12,6))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :]) 

  ax.plot(x_Borisnr,portho_Borisnr,
        color='r',lw=2,label='Boris nr')
  ax.plot(x_Boris,portho_Boris,
        color='b',lw=2,label='Boris')
  ax.plot(x_Vay,portho_Vay,
        color='green',ls='--',lw=2,label='Vay')
  ax.plot(x_HC,portho_HC,
        color='orange',ls=':',lw=2,label='HC')

  ax.set_xlabel(r"$x\omega_0/c$")
  ax.set_ylabel(r"$p_{\perp}$")

  ax.legend(loc="best")

  plt.tight_layout()

# x - Px
if True:

  fig = plt.figure(figsize=(12,6))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :]) 

  ax.plot(x_Borisnr,px_Borisnr,
        color='r',lw=2,label='Boris nr')
  ax.plot(x_Boris,px_Boris,
        color='b',lw=2,label='Boris')
  ax.plot(x_Vay,px_Vay,
        color='green',ls='--',lw=2,label='Vay')
  ax.plot(x_HC,px_HC,
        color='orange',ls=':',lw=2,label='HC')

  ax.set_xlabel(r"$x\omega_0/c$")
  ax.set_ylabel(r"$p_{\parallel}$")

  ax.legend(loc="best")

  plt.tight_layout()

# x - Gamma
if True:

  fig = plt.figure(figsize=(12,6))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :]) 

  ax.plot(x_Borisnr,gf_Borisnr,
        color='r',lw=2,label='Boris nr')
  ax.plot(x_Boris,gf_Boris,
        color='b',lw=2,label='Boris')
  ax.plot(x_Vay,gf_Vay,
        color='green',ls='--',lw=2,label='Vay')
  ax.plot(x_HC,gf_HC,
        color='orange',ls=':',lw=2,label='HC')

  ax.set_xlabel(r"$x\omega_0/c$")
  ax.set_ylabel(r"$\gamma$")

  ax.legend(loc="lower center")

  plt.tight_layout()

print
print('P_{ortho, Boris}:' + str(max(portho_Boris)) + ' ' + str(a0/m.sqrt(2.)))
print('P_{ortho, Vay}: ' + str(max(portho_Vay)) + ' ' + str(a0/m.sqrt(2.)))
print('P_{ortho, HC}: ' + str(max(portho_HC)) + ' ' + str(a0/m.sqrt(2.)))

print
print('P_{parallel, Boris}: ' + str(max(px_Boris)) + ' ' + str(a0**2/2))
print('P_{parallel, Vay}: ' + str(max(px_Vay)) + ' ' + str(a0**2/2.))
print('P_{parallel, HC}: ' + str(max(px_HC)) + str(a0**2/2.))

print
print('Gamma_{Boris}: ' + str(max(gf_Boris)) + ' ' + str((a0**2+2.)/2.))
print('Gamma_{Vay}: ' + str(max(gf_Vay)) + ' ' + str((a0**2+2.)/2.))
print('Gamma_{HC}: ' + str(max(gf_HC)) + ' ' + str((a0**2+2.)/2.))

plt.tight_layout()

plt.show()
