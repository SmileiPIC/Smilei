from numpy import log
import h5py

with h5py.File("./restart000/BinaryProcesses0.h5") as f:
	debye = f["t00000001/debyelength"][()].reshape((16,16))

logdebye = log( debye )

Validate("log( debye length )", logdebye, 0.1)

