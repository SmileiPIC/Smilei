import os, re, numpy as np
import happi

S = happi.Open(["./restart*"], verbose=False)



for name, profile in S.namelist.profiles.items():
	A = S.Field.Field0("Rho_"+name)
	data = A.get()
	values = data["data"][0]
#	z0 = np.pi
#	y,x = np.meshgrid( A.get()["x"], A.get()["y"] )
#	v = x[:,:]
#	for i in range(x.shape[0]):
#		for j in range(x.shape[1]):
#			v[i,j] = profile(x[i,j],y[i,j], z0)
	Validate("Profile "+name, values[::10,::10,::10], 0.01 )


#fig=plt.figure(1)
#for name, profile in S.namelist.profiles.items():
#	fig.clf()
#	ax1=fig.add_subplot(3,1,1)
#	print "Rho_"+name
#	A=S.Field.Field0("Rho_"+name, average={"z":3.})
#	z0 = A.getData()[0]
#	plt.colorbar( ax1.imshow(z0) )
#	y,x = np.meshgrid( A.get()["x"], A.get()["y"] )
#	z = x[:,:]
#	for i in range(x.shape[0]):
#		for j in range(x.shape[1]):
#			z[i,j] = profile(x[i,j],y[i,j], 3.)
#	ax2=fig.add_subplot(3,1,2)
#	plt.colorbar( ax2.imshow(z) )
#	ax3=fig.add_subplot(3,1,3)
#	plt.colorbar( ax3.imshow(np.log10(np.abs(z-z0))) )
#	plt.draw()
#	plt.waitforbuttonpress()
#	print "--------"
#


# Verify external fields
Lz = S.namelist.Main.grid_length[2]
for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
	F = S.Field.Field0(field, timesteps=0, subset={"z":Lz/2.}).getData()[0][::10,::10]
	Validate(field+" field", F, 0.01)

# Verify prescribed fields
for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
	F = S.Field.Field0(field, timesteps=4, subset={"z":Lz/2.}).getData()[0][::10,::10]
	Validate("prescribed "+field+" field", F, 0.01)

# Verify initialization from numpy
A = ["x","y","z","px","py","pz"]
T = S.TrackParticles("from_numpy",axes=A)
Validate("Number of particles from numpy", T.nParticles)
D = T.getData()
Validate("Particles from numpy", np.concatenate( [D[a] for a in A] ), 0.00001)
