from os.path import join, dirname, realpath, isfile
from subprocess import check_output

try:
	out1 = check_output('git describe --tags --always 2>/dev/null', shell=True).decode().strip('v\n')
	out2 = check_output('git rev-parse --abbrev-ref HEAD  2>/dev/null', shell=True).decode().strip('v\n')
	print(str(out1.strip())+"-"+str(out2.strip()))
except Exception as e:
	f_name = join( dirname( realpath(__file__) ), "version" )
	if isfile(f_name):
		with open(f_name, 'r') as f:
			print("v-"+f.read())
	else:
		print("??-??")


