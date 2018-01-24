
# For backwards compatibility, this scripts imports Smilei's python module "happi"

import os
import inspect

this_script_dir = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))
s = os.sep
module_dir = this_script_dir+s+'..'+s
if not os.path.isfile(module_dir+'happi'+s+'__init__.py'):
	print("The Diagnostics.py file cannot be open this way.")
	print("Use the function `execfile`, or the `%run` command in ipython.")
	print("In python 3, replace the command `execfile(file)` by `exec(compile(open(file).read(), file, 'exec'))`")
else:
	import sys
	sys.path.insert(0, module_dir)
	import happi
