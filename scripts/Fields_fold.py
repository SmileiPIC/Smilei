import numpy as np
import h5py

f_in  = h5py.File("Fields.h5")
f_out = h5py.File("Fields_unfolded.h5", "w")

# Size of a patch (in each dimension)
patch_size = f_in.attrs["patch_size"]
patch_size = patch_size.astype(int)
# Total number of cells in a patch
total_patch_size = patch_size.prod()
# Global number of cells
total_size = f_in.itervalues().next().itervalues().next().shape[0]
# Dimension of the simulation
dim = len(patch_size)
# Global number of patches
npatches = total_size / total_patch_size
# Global simulation shape
simulation_shape = f_in.attrs["sim_length"]/f_in.attrs["cell_length"]
simulation_shape = simulation_shape.astype(int) + 1

# Function that tells the patch coordinates
def patchSlices1D( hindex ):
	bonus = int(hindex>0)
	return (
		# slice in output
		(slice(hindex*(patch_size[0]-1)+bonus,(hindex+1)*(patch_size[0]-1)+1),),
		# slice in patch
		(slice(bonus, patch_size[0]),),
	)
if   dim == 1: patchSlices = patchSlices1D
else         : raise Exception(str(dim)+"D unavailable")

# Loop timesteps
for time in f_in.values():
	time_out = f_out.create_group(time.name)
	
	# Loop fields
	for field in time.values():
		
		# Loop all patches and fill the array
		hindex = 0
		data = np.zeros(simulation_shape)
		for istart in range(0, total_size, total_patch_size):
			# Obtain the field in current patch
			A = np.reshape( field.value[istart:istart+total_patch_size], patch_size )
			# Put the patch in the data
			slice_in, slice_out = patchSlices( hindex );
			data[slice_in] = A[slice_out]
			hindex += 1
		
		# Write the data out
		time_out.create_dataset(field.name, data=data)
		



