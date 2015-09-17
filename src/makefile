SMILEICXX     ?=  mpixlcxx_r
#SMILEICXX     ?= skin mpixlcxx_r
#SMILEICXX     ?= mpixlcxx_r -pg -cxx=powerpc64-bgq-linux-g++ 
#SMILEICXX     ?= mpixlcxx_r -cxx=powerpc64-bgq-linux-g++ 
#SMILEICXX     ?= mpixlcxx_r 
HDF5_ROOT_DIR ?=

PYTHONCONFIG ?= python-config
##################################################

VERSION:="$(shell git describe 2>/dev/null || echo '??')"

#CFLAGS += -D__VERSION=\"$(VERSION)\" -I${HDF5_ROOT_DIR}/include
#LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5 /bgsys/linux/RHEL6.5_V1R2M2-7/opt/ibmcmp/vac/bg/12.1/bglib64/libxl.a /bgsys/linux/RHEL6.5_V1R2M2-7/opt/ibmcmp/vac/bg/12.1/bglib64/libxlopt.a
CFLAGS += -D__VERSION=\"$(VERSION)\" -I${HDF5_ROOT_DIR}/include -I/bglocal/cn/pub/HDF5/1.8.14/par/include 
LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5 /bgsys/linux/RHEL6.5_V1R2M2-7/opt/ibmcmp/vac/bg/12.1/bglib64/libxl.a /bgsys/linux/RHEL6.5_V1R2M2-7/opt/ibmcmp/vac/bg/12.1/bglib64/libxlopt.a -L/bglocal/cn/pub/HDF5/1.8.14/par/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lm -L/bglocal/cn/pub/zlib/1.2.5/lib -lz

ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz
endif

#add subdirs
DIRS = $(shell find * -type d)
#add include directives for subdirs
CFLAGS += $(DIRS:%=-I%)

EXEC = smilei

default: $(EXEC)

help: 
	@echo 'Usage: make [clean] [config=OPTIONS]'
	@echo '		OPTIONS is a string composed of one or more of:'
	@echo '	        debug    : to compile in debug mode (code runs really slow)'
	@echo '	        openmp   : to compile with openmp enabled'
	@echo '         scalasca : to compile using scalasca'
	@echo ' e.g. make config="debug openmp"'
	@echo ''
	@echo 'Environment variables : SMILEICXX (c++ compiler) HDF5_ROOT_DIR (HDF5 dir with /include and /lib)'


####################################################
DESCRIBE:=$(shell git describe 2>/dev/null || echo '??')
BRANCH:=$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '??')
COMMITDATE:=$(shell git show -s --pretty="%ci" 2>/dev/null || echo '??')

VERSION:="$(DESCRIBE)-$(BRANCH)"

CFLAGS += -D__VERSION=\"$(VERSION)\" -D__COMMITDATE=\""$(COMMITDATE)"\" -I${HDF5_ROOT_DIR}/include -std=c++0x 
LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5 -lz


ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz -L/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/lib
endif

#add subdirs
DIRS = $(shell find * -type d)
#add include directives for subdirs
CFLAGS += $(DIRS:%=-I%)

PYSCRIPTS := $(shell find Python -name \*.py)
PYHEADERS := $(PYSCRIPTS:.py=.pyh)

#collect all cpp files
SRCS := $(shell find . -name \*.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

CFLAGS+=$(shell $(PYTHONCONFIG) --includes)
LDFLAGS+=$(shell $(PYTHONCONFIG) --ldflags)

# check for variable config
ifneq (,$(findstring debug,$(config)))
	CFLAGS += -g -pg -Wall -D__DEBUG -O0# -shared-intel 
else
	CFLAGS += -O3              
endif

ifneq (,$(findstring scalasca,$(config)))
    SMILEICXX = scalasca -instrument mpic++
endif

ifneq (,$(findstring openmp,$(config)))
    ifneq (,$(findstring icpc,$(shell $(SMILEICXX) --showme:command)))
        CFLAGS += -qsmp=omp -D_OMP
    else
        CFLAGS += -fopenmp -D_OMP
    endif
endif


clean:
	rm -f $(OBJS) $(DEPS) $(PYHEADERS)

distclean: clean
	rm -f $(EXEC)

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

# this generates a .h file containing a char[] with the python script in binary then
#you can just include this file to get the contents
%.pyh: %.py
	@ echo "Creating binary char for $< : $@"
	@ xxd -i $< > $@

#%.d: %.cpp $(PYHEADERS) 
#	@ echo "Checking dependencies for $<"
# create and modify dependecy file .d to take into account the location subdir
#	@ $(SMILEICXX) $(CFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(shell  dirname $<)/\1.d $(shell  dirname $<)/\1.o:@" > $@  

%.o : %.cpp %.d
	$(SMILEICXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(PYHEADERS) $(OBJS)
	$(SMILEICXX) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

# these are kept for backward compatibility and might be removed (see make help)

obsolete:
	@echo "[WARNING] Please consider using make config=\"$(MAKECMDGOALS)\""

debug: obsolete
	make config=debug

openmp: obsolete
	make config=openmp

scalasca: obsolete
	make config=scalasca


# DO NOT DELETE
