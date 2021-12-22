from ._Utils import multiPlot, multiSlide, Units, openNamelist

import os as _os
happi_directory = _os.path.dirname(_os.path.abspath(__file__))


def Open(*args, **kwargs):
	""" Import a Smilei simulation

	Parameters:
	-----------
	results_path : string or list of strings (default '.').
		Directory containing simulation results, or list of directories.
		Omit this argument if you are already in the results directory.

	reference_angular_frequency_SI : float (default None)
		Sets or change the value of reference_angular_frequency_SI, which may
		be defined in the block Main() of any Smilei namelist.

	show : bool (default True)
		Can be set to False to prevent figures to actually appear on screen.

	verbose : bool (default True)
		If False, no warnings or information are printed.

	scan : bool (default True)
		If False, the HDF5 output files are not initially scanned.

	Returns:
	--------
	A SmileiSimulation object, i.e. a container that holds information about a simulation.

	"""
	
	from ._SmileiSimulation import SmileiSimulation
	return SmileiSimulation(*args, **kwargs)


