SMILEICXX     ?=  mpixlcxx_r
#SMILEICXX     ?= skin mpixlcxx_r
#SMILEICXX     ?= mpixlcxx_r -pg -cxx=powerpc64-bgq-linux-g++ 
#SMILEICXX     ?= mpixlcxx_r -cxx=powerpc64-bgq-linux-g++ 
#SMILEICXX     ?= mpixlcxx_r 
HDF5_ROOT_DIR ?=

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

default: smilei

help: 
	@echo 'Usage                 : make [clean] [openmp={gnu|intel}] [config=debug]'
	@echo 'Environment variables : SMILEICXX (c++ compiler) HDF5_ROOT_DIR (HDF5 dir with /include and /lib)'


# these are kept for backward compatibility and might be removed (see make help)
release:
	make

openmpgnu:
	make openmp=gnu

openmpintel:
	make openmp=intel

debug:
	make config=debug

scalasca: 
	make config=scalasca



#collect all cpp files
SRCS := $(shell find . -name \*.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

ifeq ($(config),debug)
	CFLAGS += -g -pg -Wall -D__DEBUG
else
	CFLAGS += -O3              
endif

ifeq ($(openmp),gnu)
	CFLAGS += -fopenmp -D_OMP
else ifeq ($(openmp),intel)
	CFLAGS += -qsmp=omp -D_OMP
endif

ifeq ($(config),scalasca)
	CFLAGS += -qsmp=omp -D_OMP
    SMILEICXX = scalasca -instrument mpic++
endif

clean:
	rm -f $(OBJS)
	
distclean: clean
	rm -f $(EXEC)
	
NODEPS:=clean	
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
    #Chances are, these files don't exist.  GMake will create them and
    #clean up automatically afterwards
	-include $(DEPS)
endif


#%.d: %.cpp
	#@ echo "Checking dependencies for $<"
	#$(eval DIR := $(shell  dirname $<))
# create and modify dependecy file .d to take into account the location subdir	
	#@ $(SMILEICXX) -MM $(CFLAGS) $*.cpp 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(DIR)/\1.d $(DIR)/\1.o:@" > $@  

%.o : %.cpp
	$(SMILEICXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(SMILEICXX) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

# DO NOT DELETE
