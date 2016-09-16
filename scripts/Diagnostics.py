
# For backwards compatibility, this scripts imports the python module
import os
import inspect
exec(open(os.path.dirname(os.path.abspath(inspect.stack()[0][1]))+os.sep+'PythonModule'+os.sep+'Smilei'+os.sep+'__init__.py','r').read())