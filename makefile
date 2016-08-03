
#-----------------------------------------------------
# Default variables

MPIVERSION = $(shell mpirun --version 2>&1| head -n 1)
ifneq (,$(findstring Open MPI,$(MPIVERSION)))
    SMILEICXX ?= mpicxx
else
    SMILEICXX ?= mpiicpc
endif

HDF5_ROOT_DIR ?=

BUILD_DIR ?= build

PYTHONCONFIG ?= python-config

EXEC = smilei

#-----------------------------------------------------
# Get some git information
DESCRIBE:=$(shell git describe 2>/dev/null || echo '??')
BRANCH:=$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '??')
VERSION="$(DESCRIBE)-$(BRANCH)"

#-----------------------------------------------------
# Collect directories and files
DIRS := $(shell find src -type d)
SRCS := $(shell find src/* -name \*.cpp)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.d))
PYSCRIPTS := $(shell find src/Python -name \*.py)
SITEDIR=$(shell python -c 'import site; site._script()' --user-site)

#-----------------------------------------------------
# Build flags 

# Smilei version
CXXFLAGS += -D__VERSION=\"$(VERSION)\"
# C++ version
CXXFLAGS += -std=c++0x 
# HDF5 library
ifneq ($(strip $(HDF5_ROOT_DIR)),)
CXXFLAGS += -I${HDF5_ROOT_DIR}/include 
LDFLAGS += -L${HDF5_ROOT_DIR}/lib 
endif
LDFLAGS += -lhdf5 
# Include subdirs
CXXFLAGS += $(DIRS:%=-I%)
# Python-related flags
CXXFLAGS += -I$(BUILD_DIR)/src/Python
PYHEADERS := $(addprefix $(BUILD_DIR)/, $(PYSCRIPTS:.py=.pyh))
PY_CXXFLAGS:=$(shell $(PYTHONCONFIG) --includes)
CXXFLAGS+=$(PY_CXXFLAGS)
PY_LDFLAGS:=$(shell $(PYTHONCONFIG) --ldflags)
LDFLAGS+=$(PY_LDFLAGS)
ifneq ($(strip $(PYTHONHOME)),)
LDFLAGS+=-L$(PYTHONHOME)/lib
endif 

#-----------------------------------------------------
# Machine-specific configuration

ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz -L/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/lib
endif

ifneq (,$(findstring turing,$(config)))
	CXXFLAGS += -I$(BG_PYTHONHOME)/include/python2.7 -qlanglvl=extended0x
	LDFLAGS  += -qnostaticlink -L(BG_PYTHONHOME)/lib64 -lpython2.7 -lutil
endif

ifneq (,$(findstring icpc,$(SMILEI_COMPILER)))
    CXXFLAGS += -xHost -no-vec
endif

#-----------------------------------------------------
# Manage options in the "config" parameter

ifneq (,$(findstring debug,$(config)))
	CXXFLAGS += -g -pg -Wall -D__DEBUG -O0 # -shared-intel 
else
	CXXFLAGS += -O3 #-ipo
endif

ifneq (,$(findstring scalasca,$(config)))
    SMILEICXX = scalasca -instrument $(SMILEICXX)
endif

ifeq (,$(findstring noopenmp,$(config)))
    OPENMP_FLAG ?= -fopenmp 
	LDFLAGS += -lm
    OPENMP_FLAG += -D_OMP
    LDFLAGS += $(OPENMP_FLAG)
    #LDFLAGS += -mt_mpi
    CXXFLAGS += $(OPENMP_FLAG)
endif


#-----------------------------------------------------
# Set verbosity prefix
ifeq (,$(findstring verbose,$(config)))
Q := @
else
Q := 
endif


#-----------------------------------------------------
# Rules

default: $(EXEC)

clean:
	@echo "Cleaning $(BUILD_DIR)"
	$(Q) rm -rf $(BUILD_DIR) 
	$(Q) rm -rf smilei-$(VERSION).tgz
	$(Q) make -C doc clean

distclean: clean
	$(Q) rm -f $(EXEC)

env:
	@echo "$(MPIVERSION)"

# Deprecated rules
obsolete:
	@echo "[WARNING] Please consider using make config=\"$(MAKECMDGOALS)\""

debug: obsolete
	make config=debug

scalasca: obsolete
	make config=scalasca

# Create python header files
$(BUILD_DIR)/%.pyh: %.py
	@echo "Creating binary char for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) python scripts/CompileTools/hexdump.py "$<" "$@"

# Calculate dependencies
$(BUILD_DIR)/%.d: %.cpp
	@echo "Checking dependencies for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
# create and modify dependecy file .d to take into account the location subdir
	$(Q) $(SMILEICXX) $(CXXFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(BUILD_DIR)/$(shell  dirname $<)/\1.d $(BUILD_DIR)/$(shell  dirname $<)/\1.o:@" > $@  

$(BUILD_DIR)/%.o : %.cpp
	@echo "Compiling $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	@echo "Linking $@"
	$(Q) $(SMILEICXX) $(OBJS) -o $(BUILD_DIR)/$@ $(LDFLAGS) 
	$(Q) cp $(BUILD_DIR)/$@ $@

# ????
ifeq ($(filter clean help doc tar,$(MAKECMDGOALS)),) 
# Let's try to make the next lines clear: we include $(DEPS) and pygenerator
-include $(DEPS) pygenerator
# and pygenerator will create all the $(PYHEADERS) (which are files)
pygenerator : $(PYHEADERS)
endif


.PHONY: pygenerator doc help clean default tar env

doc:
	make -C doc all

sphinx:
	make -C doc/Sphinx html
tar:
	git archive -o smilei-$(VERSION).tgz --prefix smilei-$(VERSION)/ HEAD

# Install the python module in the python path
install_python:
	@echo "Installing $(SITEDIR)/smilei.pth"
	$(Q) mkdir -p "$(SITEDIR)"
	$(Q) echo "$(PWD)/scripts/PythonModule" > "$(SITEDIR)/smilei.pth"

uninstall_python:
	@echo "Uninstalling $(SITEDIR)/smilei.pth"
	$(Q) rm "$(SITEDIR)/smilei.pth"

help: 
	@echo 'Usage: make config=OPTIONS'
	@echo '	    OPTIONS is a string composed of one or more of:'
	@echo '	        verbose    : to print compile command lines'
	@echo '	        debug      : to compile in debug mode (code runs really slow)'
	@echo '         scalasca   : to compile using scalasca'
	@echo '         noopenmp   : to compile without openmp'
	@echo ' examples:'
	@echo '     make config=verbose'
	@echo '     make config=debug'
	@echo '     make config=noopenmp'
	@echo '     make config="debug noopenmp"'
	@echo 
	@echo ' other usage:' 
	@echo '      make doc     : builds the documentation'
	@echo '      make tar     : creates an archive of the sources'
	@echo '      make clean   : remove build'
	@echo ''
	@echo 'Environment variables :'
	@echo '     SMILEICXX     : mpi c++ compiler'
	@echo '     HDF5_ROOT_DIR : HDF5 dir'
	@echo '     BUILD_DIR     : directory used to store build files (default: "build")'
	@echo '     OPENMP_FLAG   : flag to use openmp (default: "-fopenmp")'

