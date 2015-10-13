#This scriptts simulate the behaviour of init_patch_count
#We assume the load per patch is constant (no vacuum)
#We assume all MPI processes run on device with the same capability

import scipy
import matplotlib.pyplot as plt

#   USer defined parameters
nmpi = 24
npatchx = 32
npatchy = 32
###########################

mat = scipy.fromfile("data.txt",sep=" ",dtype=int)
mat=mat.reshape(-1,3)
x=mat[:,1]
y=mat[:,2]

mat_plot = scipy.zeros((npatchx,npatchy))
patch_count = scipy.zeros((nmpi))

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
        mat_plot[mat[i,1],mat[i,2]] = j

print "rank 7 has patches from ", patch_count[:7].sum(), " to ", patch_count[:8].sum()-1


plt.matshow(mat_plot,aspect="auto")
plt.plot(x,y,color='black',lw=2)
plt.show()
