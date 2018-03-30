import os, re, numpy as np
from scipy.signal import butter, filtfilt
import happi

S = happi.Open(["./restart*"], verbose=False)



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
#	A=S.Field.Field0("Rho_"+name, average={"z":np.pi})
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




# Maxwell-Juttner initialization
p = S.ParticleBinning(0).get()
p_distr = p["data"][0]
p = p["px"]
b, a = butter(8, 0.15, btype='low', analog=False)
p_filt = filtfilt(b, a, p_distr)
# # theory
# Te = S.namelist.Te
# fth = (np.sqrt(1.+p**2)+Te) * np.exp( -1./Te* np.sqrt(1.+p**2) )
# itg = (p[1]-p[0])*np.sum(fth)
# fth = fth/itg
# plt.plot(p, p_filt, '.k')
# plt.plot(p, fth, '-k')
Validate("Maxwell-Juttner Momentum distribution", p_filt, p_filt.max()*1e-2)


# Verify external fields
#for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
for field in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]:
	F = S.Field.Field0(field, timesteps=0).getData()[0][::4,::4]
	Validate(field+" field", F, 0.01)


