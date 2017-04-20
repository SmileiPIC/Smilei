from Smilei import *
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

res_Boris = Smilei("./tst1d_cir_plane_wave_Boris")
res_Vay = Smilei("./tst1d_cir_plane_wave_Vay")
res_HC = Smilei("./tst1d_cir_plane_wave_HC")

# ________________________________________________
# Parameters

#a0 = res_Boris.namelist.LaserGaussian3D[0].a0
a0 = 2.

# _________________________________________
# Scalar

Ukin = res_Boris.Scalar("Ukin").get()
Ukin_Vay = res_Vay.Scalar("Ukin").get()
Ukin_HC = res_HC.Scalar("Ukin").get()
  
fig = plt.figure(figsize=(12,6))
gs = gridspec.GridSpec(2, 2)
ax = plt.subplot(gs[:, :]) 
  
#print scalar_HC["times"],scalar_HC["data"]

#Ukin["data"] = np.array(Ukin["data"])/(192.*1e-8)
#Ukin_Vay["data"] = np.array(Ukin_Vay["data"])/(192.*1e-8)
#Ukin_HC["data"] = np.array(Ukin_HC["data"])/(192.*1e-8)
 
ax.plot(Ukin["times"],Ukin["data"],
        color='b',
        ls='-',
        lw=2,
        marker='v',
        mec='b',
        markevery=100)
  
ax.plot(Ukin_Vay["times"],Ukin_Vay["data"],
        color='green',
        ls='--',
        lw=2)

ax.plot(Ukin_HC["times"],Ukin_HC["data"],
        color='orange',
        ls=':',
        lw=2)

ax.set_xlabel('time steps')
ax.set_ylabel('Normalized Kinetic Energy')
  
plt.tight_layout()
   
# ____________________________________________
# Fields

if False:

  Ey = res_Boris.Field(0, "Ey", timesteps=1300,slice = {"z":[5.]}).get()
 
  Ey["data"] = np.array(Ey["data"][0].T)
 
  fig = plt.figure(figsize=(12,6))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :]) 
 

# ____________________________________________
# Track

Track = res_Boris.TrackParticles("electron", axes=["x","px","py","pz"]).get()
Track_Vay = res_Vay.TrackParticles("electron", axes=["x","px","py","pz"]).get()
Track_HC = res_HC.TrackParticles("electron", axes=["x","px","py","pz"]).get()

#gamma_Boris = m.sqrt(1 + np.array(Track['px'][:,0])**2 + np.array(Track['py'][:,0])**2 + np.array(Track['px'][:,0])**2)

# Momentum

px_Boris = np.array(Track['px'][:,0])
py_Boris = np.array(Track['py'][:,0])
pz_Boris = np.array(Track['pz'][:,0])
portho_Boris = np.sqrt(py_Boris**2 + pz_Boris**2)

px_Vay = np.array(Track_Vay['px'][:,0])
py_Vay = np.array(Track_Vay['py'][:,0])
pz_Vay = np.array(Track_Vay['pz'][:,0])
portho_Vay = np.sqrt(py_Vay**2 + pz_Vay**2)

px_HC = np.array(Track_HC['px'][:,0])
py_HC = np.array(Track_HC['py'][:,0])
pz_HC = np.array(Track_HC['pz'][:,0])
portho_HC = np.sqrt(py_HC**2 + pz_HC**2)

if True:

  fig = plt.figure(figsize=(12,6))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :]) 

  ax.plot(Track['x'][:,0],portho_Boris,
        color='b',lw=2)
  ax.plot(Track_Vay['x'][:,0],portho_Vay,
        color='green',ls='--',lw=2)
  ax.plot(Track_HC['x'][:,0],portho_HC,
        color='orange',ls=':',lw=2)

if True:

  fig = plt.figure(figsize=(12,6))
  gs = gridspec.GridSpec(2, 2)
  ax = plt.subplot(gs[:, :]) 

  ax.plot(Track['x'][:,0],px_Boris,
        color='b',lw=2)
  ax.plot(Track_Vay['x'][:,0],px_Vay,
        color='green',ls='--',lw=2)
  ax.plot(Track_HC['x'][:,0],px_HC,
        color='orange',ls=':',lw=2)

print('P_{ortho, Boris}:' + str(max(portho_Boris)) + ' ' + str(a0/m.sqrt(2.)))
print('P_{ortho, Vay}: ' + str(max(portho_Vay)) + ' ' + str(a0/m.sqrt(2.)))
print('P_{ortho, HC}: ' + str(max(portho_HC)) + ' ' + str(a0/m.sqrt(2.)))

print('P_{parallel, Boris}: ' + str(max(px_Boris)) + ' ' + str(a0**2/2))
print('P_{parallel, Vay}: ' + str(max(px_Vay)) + ' ' + str(a0**2/2.))
print('P_{parallel, HC}: ' + str(max(px_HC)) + str(a0**2/2.))

print('Gamma_{Boris}: ' + str(max(np.sqrt(1.+px_Boris**2+py_Boris**2+pz_Boris**2))) + ' ' + str((a0**2+2.)/2.))
print('Gamma_{Vay}: ' + str(max(np.sqrt(1.+px_Vay**2+py_Vay**2+pz_Vay**2))) + ' ' + str((a0**2+2.)/2.))
print('Gamma_{HC}: ' + str(max(np.sqrt(1.+px_HC**2+py_HC**2+pz_HC**2))) + ' ' + str((a0**2+2.)/2.))

plt.tight_layout()

plt.show()
