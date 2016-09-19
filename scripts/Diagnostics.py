
# For backwards compatibility, this scripts imports Smilei's python module
#  located in scripts/PythonModule

import os
import inspect
this_script_dir = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))
s = os.sep
module_file = this_script_dir+s+'PythonModule'+s+'Smilei'+s+'__init__.py'
if not os.path.isfile(module_file):
	print("The Diagnostics.py file cannot be open this way.")
	print("Use the function `execfile`, or the `%run` command in ipython.")
else:
	exec(open(module_file,'r').read())