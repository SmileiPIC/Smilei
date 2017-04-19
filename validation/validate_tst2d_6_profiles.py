import os, re, numpy as np
from Smilei import *

S = Smilei(".", verbose=False)

for name, profile in S.namelist.profiles.items():
	A = S.Field.Field0("Rho_"+name)
	data = A.get()
	values = data["data"][0]
#	y,x = np.meshgrid( A.get()["x"], A.get()["y"] )
#	v = x[:,:]
#	for i in range(x.shape[0]):
#		for j in range(x.shape[1]):
#			v[i,j] = profile(x[i,j],y[i,j])
	Validate("Profile "+name, values[::2,::2], 0.01 )

#
#fig=plt.figure(1)
#for name, profile in S.namelist.profiles.items():
#	fig.clf()
#	ax1=fig.add_subplot(3,1,1)
#	print "Rho_"+name
#	A=S.Field.Field0("Rho_"+name, slice={"z":np.pi})
#	data = A.get()
#	v0 = data["data"][0]
#	plt.colorbar( ax1.imshow(v0) )
#	y,x = np.meshgrid( data["x"], data["y"] )
#	v = x[:,:]
#	for i in range(x.shape[0]):
#		for j in range(x.shape[1]):
#			v[i,j] = profile(x[i,j],y[i,j])
#	ax2=fig.add_subplot(3,1,2)
#	plt.colorbar( ax2.imshow(v) )
#	ax3=fig.add_subplot(3,1,3)
#	plt.colorbar( ax3.imshow(np.log10(np.abs(v-v0))) )
#	plt.draw()
#	plt.waitforbuttonpress()
#	print "--------"
#