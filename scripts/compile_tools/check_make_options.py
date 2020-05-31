#!/usr/bin/env python 

import sys
import subprocess

def option_is_in_the_list( keyword, option, options_registered ):
    for listed in options_registered.split() :
       if (option == listed):
           sys.exit(0)

    print(option, "do not exist for", keyword)
    print("Use one of : ", options_registered)
    sys.exit(-1)


# if no options are passed, return
if (len(sys.argv)==2):
    sys.exit(0)

options_registered = ""

# look for known options for config
if ( (sys.argv)[1] == "config" ):
    cmd = subprocess.Popen('grep config makefile|grep findstring', shell=True, stdout=subprocess.PIPE)
    for line in cmd.stdout :
        options_registered +=  (line.split(",")[1]).split()[1] + " "

    for i in range(2,len(sys.argv)):
        option_is_in_the_list( "config", (sys.argv)[i], options_registered)

# look for known options for machine
if ( (sys.argv)[1] == "machine" ):
    # a single machine can be passed
    if (len(sys.argv)!=3):
        sys.exit(-1)

    cmd = subprocess.Popen('ls -1 scripts/compile_tools/machine', shell=True, stdout=subprocess.PIPE)
    for line in cmd.stdout :
        options_registered += (line.split())[0].decode() + " "

    option_is_in_the_list( "machine", (sys.argv)[2], options_registered)

