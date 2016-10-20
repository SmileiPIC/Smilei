#-----------------------------------------------------
# Variables that can be defined by the user:
# 
# SMILEICXX     : the MPI C++ executable (for instance mpicxx, mpiicpc, etc.)
# HDF5_ROOT_DIR : the local path to the HDF5 library
# BUILD_DIR     : the path to the build directory (default: ./build)
# PYTHON_CONFIG : the executable `python-config` usually shipped with python installation
# EXEC          : the name of the created executable for Smilei


MPIVERSION = $(shell mpirun --version 2>&1| head -n 1)
ifneq (,$(findstring Open MPI,$(MPIVERSION)))
    SMILEICXX ?= mpicxx
else
    SMILEICXX ?= mpiicpc
endif
HDF5_ROOT_DIR ?=
BUILD_DIR ?= build


#-----------------------------------------------------
# check if python-config exists
ifeq (,$(shell command -v python-config 2> /dev/null))
	PYTHONCONFIG := python-config
else
	PYTHONCONFIG := python scripts/CompileTools/python-config.py
endif
    

EXEC = smilei


#-----------------------------------------------------
# Set the verbosity prefix
ifeq (,$(findstring verbose,$(config)))
    Q := @
else
    Q := 
endif

#-----------------------------------------------------
# Git information
DESCRIBE:=$(shell git describe 2>/dev/null || echo '??')
BRANCH:=$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '??')
VERSION="$(DESCRIBE)-$(BRANCH)"

#-----------------------------------------------------
# Directories and files
DIRS := $(shell find src -type d)
SRCS := $(shell find src/* -name \*.cpp)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.d))
SITEDIR = $(shell python -c 'import site; site._script()' --user-site)

#-----------------------------------------------------
# Flags 

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
PYSCRIPTS = $(shell find src/Python -name \*.py)
PYHEADERS := $(addprefix $(BUILD_DIR)/, $(PYSCRIPTS:.py=.pyh))
PY_CXXFLAGS := $(shell $(PYTHONCONFIG) --includes)
CXXFLAGS += $(PY_CXXFLAGS)
PY_LDFLAGS := $(shell $(PYTHONCONFIG) --ldflags)
LDFLAGS += $(PY_LDFLAGS)
ifneq ($(strip $(PYTHONHOME)),)
    LDFLAGS += -L$(PYTHONHOME)/lib
endif 

# Machine-specific configuration
ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz -L/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/lib
endif

ifneq (,$(findstring turing,$(config)))
	CXXFLAGS += -I$(BG_PYTHONHOME)/include/python2.7 -qlanglvl=extended0x
	LDFLAGS  += -qnostaticlink -L(BG_PYTHONHOME)/lib64 -lpython2.7 -lutil
endif

# Compiler-specific configuration
ifneq (,$(findstring icpc,$(SMILEI_COMPILER)))
    CXXFLAGS += -xHost -no-vec
endif

# Manage options in the "config" parameter
ifneq (,$(findstring debug,$(config)))
	CXXFLAGS += -g -pg -Wall -D__DEBUG -O0 # -shared-intel 
else
	CXXFLAGS += -O3 # -g #-ipo
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
#LDFLAGS += -mt_mpi

#-----------------------------------------------------
# Rules for building Smilei

default: $(EXEC)

clean:
	@echo "Cleaning $(BUILD_DIR)"
	$(Q) rm -rf $(BUILD_DIR) 
	$(Q) rm -rf $(EXEC)-$(VERSION).tgz
	@echo "Cleaning doc/html"
	$(Q) rm -rf doc/html

distclean: clean uninstall_python
	$(Q) rm -f $(EXEC)
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
#	$(Q) $(SMILEICXX) $(CXXFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(BUILD_DIR)/$(shell  dirname $<)/\1.d $(BUILD_DIR)/$(shell  dirname $<)/\1.o:@" > $@  
	$(Q) $(SMILEICXX) $(CXXFLAGS) -MF"$@" -MG -MM -MP -MT"$@ $(@:.d=.o)" $<

$(BUILD_DIR)/%.o : %.cpp
	@echo "Compiling $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	@echo "Linking $@"
	$(Q) $(SMILEICXX) $(OBJS) -o $(BUILD_DIR)/$@ $(LDFLAGS) 
	$(Q) cp $(BUILD_DIR)/$@ $@

# Avoid to check dependencies and to create .pyh if not necessary
ifeq ($(filter-out $(wildcard print-*),$(MAKECMDGOALS)),) 
    FILTER_RULES=clean distclean help env obsolete debug scalasca doc doxygen sphinx tar install_python uninstall_python
    ifeq ($(filter $(FILTER_RULES),$(MAKECMDGOALS)),) 
        # Let's try to make the next lines clear: we include $(DEPS) and pygenerator
        -include $(DEPS) pygenerator
        # and pygenerator will create all the $(PYHEADERS) (which are files)
        pygenerator : $(PYHEADERS)
    endif
endif

# these are not file-related rules
.PHONY: pygenerator $(FILTER_RULES)

#-----------------------------------------------------
# Doc rules

doc: sphinx doxygen

sphinx:
	@echo "Compiling sphinx documentation in doc/html/Sphinx/html"
	$(Q) if type "sphinx-build" >/dev/null 2>&1; then\
		make -C doc/Sphinx BUILDDIR=../html/Sphinx html;\
	else \
		echo "Cannot build Sphinx doc because Sphinx is not installed";\
	fi
	
doxygen:
	@echo "Compiling doxygen documentation in doc/html/Doxygen/html"
	$(Q) if type "doxygen" >/dev/null 2>&1; then\
		mkdir -p doc/html/Doxygen; (echo "PROJECT_NUMBER=${VERSION}\nOUTPUT_DIRECTORY=doc/html/Doxygen"; cat doc/Doxygen/smilei.dox) | doxygen - ;\
	else \
		echo "Cannot build doxygen doc because doxygen is not installed";\
	fi	

tar:
	@echo "Creating archive $(EXEC)-$(VERSION).tgz"
	$(Q) git archive -o $(EXEC)-$(VERSION).tgz --prefix $(EXEC)-$(VERSION)/ HEAD


#-----------------------------------------------------
# Python module rules

# Install the python module in the user python path
install_python:
	@echo "Installing $(SITEDIR)/smilei.pth"
	$(Q) mkdir -p "$(SITEDIR)"
	$(Q) echo "$(PWD)/scripts/PythonModule" > "$(SITEDIR)/smilei.pth"

uninstall_python:
	@echo "Uninstalling $(SITEDIR)/smilei.pth"
	$(Q) rm -f "$(SITEDIR)/smilei.pth"


#-----------------------------------------------------
# Info rules

# print internal makefile variable (used to debug compilation problems)
print-% :
	$(info $* : $($*)) @true

# print a set of important internal variables
env: print-SMILEICXX print-MPIVERSION print-VERSION print-OPENMP_FLAG print-HDF5_ROOT_DIR print-SITEDIR print-PY_CXXFLAGS print-PY_LDFLAGS print-CXXFLAGS


help: 
	@echo 'TO BUILD SMILEI:'
	@echo '----------------'
	@echo 'Usage:'
	@echo '  make'
	@echo 'or, to compile with 4 cpus (for instance):'
	@echo '  make -j 4'
	@echo
	@echo 'More options:'
	@echo '  make config="[ verbose ] [ debug ] [ scalasca ] [ noopenmp ]"'
	@echo '    verbose    : to print compile command lines'
	@echo '    debug      : to compile in debug mode (code runs really slow)'
	@echo '    scalasca   : to compile using scalasca'
	@echo '    noopenmp   : to compile without openmp'
	@echo
	@echo 'Examples:'
	@echo '  make config=verbose'
	@echo '  make config=debug'
	@echo '  make config="debug noopenmp"'
	@echo
	@echo 'OTHER PURPOSES:'
	@echo '---------------'
	@echo '  make doc              : builds all the documentation'
	@echo '  make sphinx           : builds the `sphinx` documentation only (for users)'
	@echo '  make doxygen          : builds the `doxygen` documentation only (for developers)'
	@echo '  make tar              : creates an archive of the sources'
	@echo '  make clean            : cleans the build directory'
	@echo "  make install_python   : install Smilei's python module"
	@echo "  make uninstall_python : remove Smilei's python module"
	@echo '  make print-XXX        : print internal makefile variable XXX'
	@echo '  make env              : print important internal makefile variables'
	@echo ''
	@echo 'Environment variables :'
	@echo '  SMILEICXX     : mpi c++ compiler'
	@echo '  HDF5_ROOT_DIR : HDF5 dir'
	@echo '  BUILD_DIR     : directory used to store build files [$(BUILD_DIR)]'
	@echo '  OPENMP_FLAG   : flag to use openmp [$(OPENMP_FLAG)]'
	@echo 
	@echo 'http://www.maisondelasimulation.fr/smilei'

