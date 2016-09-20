import numpy as np
import h5py
import sys
import mpi4py

f_in  = h5py.File("Fields.h5")
f_out = h5py.File("Fields_folded.h5", "w")

# Size of a patch (in each dimension)
patch_size = f_in.attrs["patch_size"]
patch_size = patch_size.astype(int)
# Total number of cells in a patch
total_patch_size = patch_size.prod()
# Global number of cells
total_size = f_in.itervalues().next().itervalues().next().shape[0]
# Dimension of the simulation
dim = len(patch_size)
# Global simulation shape
simulation_shape = f_in.attrs["sim_length"]/f_in.attrs["cell_length"]
simulation_shape = simulation_shape.astype(int) + 1
# Global log2 of the number of patches in each direction
m = np.round(np.log2(simulation_shape / patch_size)).astype(int)


# Hilbert curve utilities
def bit(i, k):
	return (i>>k)&1
def setbit(i, k, value):
	return (i & ~(1<<k)) | (value<<k)
def gc(i):
	return i^(i>>1)
def rotl(value, shift, dim):
	sd = shift%dim
	dd = 2**dim-1
	return (value << sd) & dd | ((value & dd) >> (dim-sd))
def tedinv(e, d, b, dim):
	return rotl( b , d+1, dim ) ^ e
def tsb(i):
	k = 0
	while i & 1:
		i = i>>1
		k+=1
	return k
def direction(i, dim):
	if  i == 0: return 0
	elif i & 1: return tsb(i)%dim
	else      : return tsb(i-1)%dim
def entry(i):
	if i == 0: return 0
	else     : return gc(2*((i-1)/2))

# 2D Hilbert curve
def generalHilbert(m0, m1, h):
	shift = 0
	mmin = min(m0,m1)
	e = 0
	d = 1 if m0 < m1 else 0
	# First define in which sub-hypercube of side 2^mmin the point is.
	for i in range(m0+m1-1, mmin+mmin-1, -1):
		l = bit(h,i) 
		shift += l*(1<<(i-mmin))
		h -= l*(1<<i)
	# Run the cubic inversion algorithm in the sub hypercube.
	x = 0
	y = 0
	for i in range(mmin-1, -1 ,-1):
		w = ((bit(h,2*i+1))<<1) + bit(h,2*i)
		l = gc(w)
		l = tedinv(e,d,l,2)
		x = setbit(x,i, bit(l,0))
		y = setbit(y,i, bit(l,1))
		e = e ^ (rotl(entry(w), d+1, 2))
		d = (d + direction(w, 2) +1 )%2
	# Shift the appropriate coordinate by the necessary value.
	if m0 >= m1: x += shift
	else       : y += shift
	return (x,y)


# Function that tells the patch coordinates
def patchSlices1D( hindex ):
	return (
		# slice in output
		(slice(hindex*(patch_size[0]-1)+int(hindex>0),(hindex+1)*(patch_size[0]-1)+1),),
		# slice in patch
		(slice(int(hindex>0), patch_size[0]),),
	)
def patchSlices2D( hindex ):
	px, py = generalHilbert(m[0], m[1], hindex)
	return (
		# slice in output
		(slice(px*(patch_size[0]-1)+int(px>0),(px+1)*(patch_size[0]-1)+1),
		 slice(py*(patch_size[1]-1)+int(py>0),(py+1)*(patch_size[1]-1)+1),),
		# slice in patch
		(slice(int(px>0), patch_size[0]),
		 slice(int(py>0), patch_size[1]),),
	)
if   dim == 1: patchSlices = patchSlices1D
elif dim == 2: patchSlices = patchSlices2D
else         : raise Exception(str(dim)+"D unavailable")


# Loop timesteps
for t in f_in.values():
	t_out = f_out.create_group(t.name)
        print "t = ", t
	
	# Loop fields
	for field in t.values():

                print "field = ", field
		
		# Loop all patches and fill the array
		hindex = 0
		data = np.zeros(simulation_shape)
                rawdata = field.value
		for istart in range(0, total_size, total_patch_size):
			# Obtain the field in current patch
			#A = np.reshape( field.value[istart:istart+total_patch_size], patch_size )
			A = np.reshape( rawdata[istart:istart+total_patch_size], patch_size )
			# Put the patch in the data
			slice_in, slice_out = patchSlices( hindex )
			data[slice_in] = A[slice_out]
			hindex += 1
		
		# Write the data out
		t_out.create_dataset(field.name, data=data)


