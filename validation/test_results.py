from easi import Validation
V = Validation()
V.compare("benchname.py","path_to_bench")
#Directory "path_to_bench" must contain a directory named restart containing a valid Smilei simulation which will be compared to the reference results of "benchname.py".
