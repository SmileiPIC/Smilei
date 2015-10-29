#This scriptts simulate the behaviour of init_patch_count
#We assume the load per patch is constant (no vacuum)
#We assume all MPI processes run on device with the same capability

import scipy
import matplotlib.pyplot as plt

#   USer defined parameters
nmpi = 12
npatchx = 256
npatchy = 32
###########################
Ndpi=1000

mat = scipy.fromfile("data.txt",sep=" ",dtype=int)
mat=mat.reshape(-1,3)
x=mat[:,1]
y=mat[:,2]

mat_plot = scipy.zeros((npatchy,npatchx),dtype=int)
patch_count = scipy.zeros((nmpi),dtype=int)

total_patches = npatchx*npatchy
target_load = total_patches/float(nmpi)

Tcur = target_load #Current target
Ncur = 0           #Current number of patches
Lcur = 0           #current Load
r = 0              #MPI rank

for patch in range(total_patches):
    Lcur += 1  #Load of 1 patch = 1.
    Ncur += 1  #Give one more patch to current MPI rank
    if (Lcur > Tcur or nmpi-r >= total_patches-patch):
        above_target = Lcur -  Tcur
        below_target = Tcur - (Lcur-1)
        if (above_target > below_target):
            patch_count[r] = Ncur - 1
            Ncur = 1
        else:
            patch_count[r] = Ncur
            Ncur = 0
        r = r + 1
        Tcur = Tcur + target_load

print patch_count

for j in range(nmpi):
    for i in scipy.arange(patch_count[j])+patch_count[:j].sum():
        mat_plot[y[i],x[i]] = j


fig, ax1 = plt.subplots(figsize=(3.48,3.543),dpi=200)
#plt.matshow(mat_plot,aspect="auto")
#ax1.axis('off')
ax1.plot(x,y,color='black',lw=0.1,marker='o',markersize=1.)
image=plt.imshow(mat_plot, origin=0, aspect='auto',interpolation='nearest',cmap="gray", vmin=-2, vmax=mat_plot.max()+2 )
#ax1.set_position([0.05,0.05,0.85,0.85])
#ax1.axes.set_visible(False)
ax1.xaxis.set_tick_params(direction='out',labelsize=7)
ax1.yaxis.set_tick_params(direction='out',labelsize=7)
#ax1.yaxis.set_ticks([0.,8.,16.,24.,31.])
#ax1.xaxis.set_ticks([0.,8.,16.,24.,31.])
ax1.set_ylabel(r'$\rm Patches\ Y\ coordinate$',fontsize=8)
ax1.set_xlabel(r'$\rm Patches\ X\ coordinate$',fontsize=8)

#box = image.axes.get_window_extent()

plt.savefig("Hilbert_curve.eps",format='eps',dpi=Ndpi, frameon=False)
#plt.show()
