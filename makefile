#-----------------------------------------------------
# Variables that can be defined by the user:
# 
# SMILEICXX     : the MPI C++ executable (for instance mpicxx, mpiicpc, etc.)
# HDF5_ROOT_DIR : the local path to the HDF5 library
# BUILD_DIR     : the path to the build directory (default: ./build)
# PYTHON_CONFIG : the executable `python-config` usually shipped with python installation

SMILEICXX ?= mpicxx
HDF5_ROOT_DIR ?= 
BUILD_DIR ?= build
PYTHONEXE ?= python

PYTHONCONFIG := $(PYTHONEXE) scripts/CompileTools/python-config.py

#-----------------------------------------------------
# Git information
VERSION:=$(shell $(PYTHONEXE) scripts/CompileTools/get-version.py )

#-----------------------------------------------------
# Directories and files
DIRS := $(shell find src -type d)
SRCS := $(shell find src/* -name \*.cpp)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.d))
SITEDIR = $(shell $(PYTHONEXE) -c 'import site; site._script()' --user-site)

#-----------------------------------------------------
# Flags 

# Smilei version
CXXFLAGS += -D__VERSION=\"$(VERSION)\"
# C++ version
CXXFLAGS += -std=c++11 -Wall 
# HDF5 library
ifneq ($(strip $(HDF5_ROOT_DIR)),)
CXXFLAGS += -I${HDF5_ROOT_DIR}/include 
LDFLAGS := -L${HDF5_ROOT_DIR}/lib $(LDFLAGS)
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


PICSAR=FALSE
ifeq ($(PICSAR),TRUE)
        # New environment variable
	FFTW3_LIB ?= $(FFTW_LIB_DIR)
	LIBPXR ?= picsar/lib
	# Set Picsar link environment
	CXXFLAGS += -D_PICSAR
	LDFLAGS += -L$(LIBPXR) -lpxr
	LDFLAGS += -L$(FFTW3_LIB) -lfftw3_mpi
	LDFLAGS += -L$(FFTW3_LIB) -lfftw3_threads
	LDFLAGS += -L$(FFTW3_LIB) -lfftw3
	LDFLAGS += -lgfortran
endif

# Manage options in the "config" parameter
ifneq (,$(findstring debug,$(config)))
    CXXFLAGS += -g -pg -D__DEBUG -O0
# With gdb
else ifneq (,$(findstring gdb,$(config)))
    CXXFLAGS += -v -da -Q

# With valgrind
else ifneq (,$(findstring valgrind,$(config)))
    CXXFLAGS += -g  -O3

# Scalasca
else ifneq (,$(findstring scalasca,$(config)))
    CXXFLAGS += -g  -O3
    SMILEICXX = scalasca -instrument $(SMILEICXX)

# With Intel Advisor / Vtune
else ifneq (,$(findstring advisor,$(config)))
    CXXFLAGS += -g -O3 -debug inline-debug-info -shared-intel -parallel-source-info=2

# Optimization report
else ifneq (,$(findstring opt-report,$(config)))
    CXXFLAGS += -O3 -qopt-report5

# Default configuration
else
    CXXFLAGS += -O3 #-xHost -no-prec-div -ipo
endif


ifeq (,$(findstring noopenmp,$(config)))
    OPENMP_FLAG ?= -fopenmp 
    LDFLAGS += -lm
    OPENMP_FLAG += -D_OMP
    LDFLAGS += $(OPENMP_FLAG)
    CXXFLAGS += $(OPENMP_FLAG)
#else 
#    LDFLAGS += -mt_mpi # intelmpi only
endif


#-----------------------------------------------------
# check whether to use a machine specific definitions
ifneq ($(machine),)
	ifneq ($(wildcard scripts/CompileTools/machine/$(machine)),)
	-include scripts/CompileTools/machine/$(machine)
	endif
endif
#-----------------------------------------------------
# Set the verbosity prefix
ifeq (,$(findstring verbose,$(config)))
    Q := @
else
    Q := 
endif

#-----------------------------------------------------
# Rules for building the excutable smilei

EXEC = smilei

default: $(EXEC) $(EXEC)_test

clean:
	@echo "Cleaning $(BUILD_DIR)"
	$(Q) rm -rf $(BUILD_DIR) 
	$(Q) rm -rf $(EXEC)-$(VERSION).tgz

distclean: clean uninstall_happi
	$(Q) rm -f $(EXEC) $(EXEC)_test
	

# Create python header files
$(BUILD_DIR)/%.pyh: %.py
	@echo "Creating binary char for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(PYTHONEXE) scripts/CompileTools/hexdump.py "$<" "$@"

# Calculate dependencies
$(BUILD_DIR)/%.d: %.cpp
	@echo "Checking dependencies for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(SMILEICXX) $(CXXFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

$(BUILD_DIR)/src/Diagnostic/DiagnosticScalar.o : src/Diagnostic/DiagnosticScalar.cpp
	@echo "SPECIAL COMPILATION FOR $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -O2 -c $< -o $@

# Compile cpps
$(BUILD_DIR)/%.o : %.cpp
	@echo "Compiling $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -c $< -o $@

# Link the main program
$(EXEC): $(OBJS)
	@echo "Linking $@"
	$(Q) $(SMILEICXX) $(OBJS) -o $(BUILD_DIR)/$@ $(LDFLAGS) 
	$(Q) cp $(BUILD_DIR)/$@ $@

# Compile the the main program again for test mode
$(BUILD_DIR)/src/Smilei_test.o: src/Smilei.cpp $(EXEC)
	@echo "Compiling src/Smilei.cpp for test mode"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -DSMILEI_TESTMODE -c src/Smilei.cpp -o $@

# Link the main program for test mode
$(EXEC)_test : $(OBJS:Smilei.o=Smilei_test.o)
	@echo "Linking $@ for test mode"
	$(Q) $(SMILEICXX) $(OBJS:Smilei.o=Smilei_test.o) -o $(BUILD_DIR)/$@ $(LDFLAGS)
	$(Q) cp $(BUILD_DIR)/$@ $@

# Avoid to check dependencies and to create .pyh if not necessary
FILTER_RULES=clean distclean help env debug doc tar happi uninstall_happi
ifeq ($(filter-out $(wildcard print-*),$(MAKECMDGOALS)),) 
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

doc:
	$(Q) if type "sphinx-build" >/dev/null 2>&1; then\
		make -C doc/Sphinx BUILDDIR=../../$(BUILD_DIR) html;\
		echo "Sphinx documentation in $(BUILD_DIR)/html/index.html";\
	else \
		echo "Cannot build Sphinx doc because Sphinx is not installed";\
	fi
	
#-----------------------------------------------------
# Archive in tgz file

tar:
	@echo "Creating archive $(EXEC)-$(VERSION).tgz"
	$(Q) git archive -o $(EXEC)-$(VERSION).tgz --prefix $(EXEC)-$(VERSION)/ HEAD
	$(Q) tar -zxf $(EXEC)-$(VERSION).tgz
	$(Q) echo $(VERSION) > $(EXEC)-$(VERSION)/.version
	$(Q) tar -czf $(EXEC)-$(VERSION).tgz $(EXEC)-$(VERSION) && rm -R $(EXEC)-$(VERSION)


#-----------------------------------------------------
# Python module rules

# Install the python module in the user python path

happi:
	@echo "Installing $(SITEDIR)/smilei.pth"
	$(Q) mkdir -p "$(SITEDIR)"
	$(Q) echo "$(CURDIR)" > "$(SITEDIR)/smilei.pth"

uninstall_happi:
	@echo "Uninstalling $(SITEDIR)/smilei.pth"
	$(Q) rm -f "$(SITEDIR)/smilei.pth"


#-----------------------------------------------------
# Info rules
print-% :
	$(info $* : $($*)) @true

env: print-SMILEICXX print-PYTHONEXE print-MPIVERSION print-VERSION print-OPENMP_FLAG print-HDF5_ROOT_DIR print-SITEDIR print-PY_CXXFLAGS print-PY_LDFLAGS print-CXXFLAGS print-LDFLAGS	


#-----------------------------------------------------
# help

help:
	@echo 'TO BUILD SMILEI:'
	@echo '----------------'
	@echo 'Usage:'
	@echo '  make'
	@echo 'or, to compile with 4 cpus (for instance):'
	@echo '  make -j 4'
	@echo
	@echo 'Config options:'
	@echo '  make config="[ verbose ] [ debug ] [ scalasca ] [ noopenmp ]"'
	@echo '    verbose              : to print compile command lines'
	@echo '    debug                : to compile in debug mode (code runs really slow)'
	@echo '    scalasca             : to compile using scalasca'
	@echo '    noopenmp             : to compile without openmp'
	@echo
	@echo 'Examples:'
	@echo '  make config=verbose'
	@echo '  make config=debug'
	@echo '  make config="debug noopenmp"'
	@echo
	@echo 'Machine options:'
	@echo '  make machine=XXX      : include machine file in scripts/CompileTools/machine/XXX'
	@echo '  make machine=XXX help : print help for machine'
	@echo
	@echo 'OTHER PURPOSES:'
	@echo '---------------'
	@echo '  make doc              : builds the documentation'
	@echo '  make tar              : creates an archive of the sources'
	@echo '  make clean            : cleans the build directory'
	@echo "  make happi            : install Smilei's python module"
	@echo "  make uninstall_happi  : remove Smilei's python module"
	@echo '  make env              : print important internal makefile variables'
	@echo '  make print-XXX        : print internal makefile variable XXX'
	@echo ''
	@echo 'Environment variables :'
	@echo '  SMILEICXX             : mpi c++ compiler [$(SMILEICXX)]'
	@echo '  HDF5_ROOT_DIR         : HDF5 dir [$(HDF5_ROOT_DIR)]'
	@echo '  BUILD_DIR             : directory used to store build files [$(BUILD_DIR)]'
	@echo '  OPENMP_FLAG           : openmp flag [$(OPENMP_FLAG)]'
	@echo '  PYTHONEXE             : python executable [$(PYTHONEXE)]'	
	@echo '  FFTW3_LIB             : FFTW3 libraries directory [$(FFTW3_LIB)]'
	@echo '  LIB PXR               : Picsar library directory [$(LIBPXR)]'
	@echo 
	@echo 'http://www.maisondelasimulation.fr/smilei'
	@echo 'https://github.com/SmileiPIC/Smilei'
	@echo
	@if [ -f  scripts/CompileTools/machine/$(machine) ]; then echo "Machine comments for $(machine):"; grep '^#' scripts/CompileTools/machine/$(machine); fi

