
# --------------------------------------------
# Different tools for the validation process
# --------------------------------------------

import os, sys
from time import ctime, strftime

try:
    execfile = execfile
except: # python3
    def execfile(file, globals=globals()):
        exec(compile(open(file).read(), file, 'exec'), globals)

try:
    raw_input
except: # python3
    raw_input = input

def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def date(BIN_NAME):
    return os.stat(BIN_NAME).st_ctime

def date_string(BIN_NAME):
    date_integer = date(BIN_NAME)
    date_time = ctime(date_integer)
    return date_time.replace(" ","-")
