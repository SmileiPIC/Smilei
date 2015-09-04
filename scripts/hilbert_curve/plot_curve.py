import scipy
import matplotlib.pyplot as plt

mat = scipy.fromfile("data.txt",sep=" ",dtype=int)
mat=mat.reshape(-1,3)
nmpi=15.
step = int(mat[:,0].size/nmpi)
print "step = ", step

mat_plot=scipy.zeros((32,8))

for i in mat[:,0]:
    mat_plot[mat[i,1],mat[i,2]] = min(i/step,nmpi-1)

plt.matshow(mat_plot,aspect="auto")
plt.show()
