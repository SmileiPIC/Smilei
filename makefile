#-----------------------------------------------------
# Variables that can be defined by the user:
#
# SMILEICXX     : the MPI C++ executable (for instance mpicxx, mpiicpc, etc.)
# HDF5_ROOT_DIR : the local path to the HDF5 library
# BUILD_DIR     : the path to the build directory (default: ./build)
# PYTHON_CONFIG : the executable `python-config` usually shipped with python installation

SMILEICXX ?= mpicxx
HDF5_ROOT_DIR ?= $(HDF5_ROOT)
BOOST_ROOT_DIR ?= $(BOOST_ROOT)
BUILD_DIR ?= build
PYTHONEXE ?= python
TABLES_BUILD_DIR ?= tools/tables/build

#-----------------------------------------------------
# check whether to use a machine specific definitions
ifneq ($(machine),)
    ifneq ($(wildcard scripts/compile_tools/machine/$(machine)),)
    -include scripts/compile_tools/machine/$(machine)
    else
define errormsg
ERROR: Cannot find machine file for "$(machine)"
Available machines are:
$(shell ls -1 scripts/compile_tools/machine)
endef
    $(error $(errormsg))
    endif
endif

PYTHONCONFIG := $(PYTHONEXE) scripts/compile_tools/python-config.py

#-----------------------------------------------------
# Git information
VERSION:=$(shell $(PYTHONEXE) scripts/compile_tools/get-version.py )

#-----------------------------------------------------
# Directories and files

# Smilei
DIRS := $(shell find src -type d)
SRCS := $(shell find src/* -name \*.cpp)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.d))
SITEDIR = $(shell $(PYTHONEXE) -c 'import site; site._script()' --user-site)

# Smilei tools
TABLES_DIR := tools/tables
TABLES_SRCS := $(shell find tools/tables/* -name \*.cpp | rev | cut -d '/' -f1 | rev)
TABLES_DEPS := $(addprefix $(TABLES_BUILD_DIR)/, $(SRCS:.cpp=.d))
TABLES_OBJS := $(addprefix $(TABLES_BUILD_DIR)/, $(TABLES_SRCS:.cpp=.o))
TABLES_SRCS := $(shell find tools/tables/* -name \*.cpp)


#-----------------------------------------------------
# Flags

# Smilei version
CXXFLAGS += -D__VERSION=\"$(VERSION)\" -D_VECTO
# C++ version
CXXFLAGS += -std=c++11 -Wall #-Wshadow
# HDF5 library
ifneq ($(strip $(HDF5_ROOT_DIR)),)
CXXFLAGS += -I$(HDF5_ROOT_DIR)/include
LDFLAGS := -L$(HDF5_ROOT_DIR)/lib  $(LDFLAGS)
endif
# Boost library
ifneq ($(strip $(BOOST_ROOT_DIR)),)
CXXFLAGS += -I$(BOOST_ROOT_DIR)/include
LDFLAGS := -L$(BOOST_ROOT_DIR)/lib $(LDFLAGS)
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

my_config:=$(config)
define parse_config
$(findstring $(1),$(config))$(eval my_config:=$(filter-out $(1),$(my_config)))
endef

# Manage options in the "config" parameter
ifneq (,$(call parse_config,debug))
    CXXFLAGS += -g -pg -D__DEBUG -O0
# With gdb
else ifneq (,$(call parse_config,gdb))
    CXXFLAGS += -g -D__DEBUG -O0

# With valgrind
else ifneq (,$(call parse_config,valgrind))
    CXXFLAGS += -g  -O3

# Scalasca
else ifneq (,$(call parse_config,scalasca))
    CXXFLAGS += -g  -O3
    SMILEICXX = scalasca -instrument $(SMILEICXX)

# With Intel Advisor / Vtune
else ifneq (,$(call parse_config,advisor))
    CXXFLAGS += -g -O3 -shared-intel -debug inline-debug-info -qopenmp-link dynamic -parallel-source-info=2

# With Intel Inspector
else ifneq (,$(call parse_config,inspector))
    CXXFLAGS += -g -O0 -I$(INSPECTOR_ROOT_DIR)/include/
    LDFLAGS += $(INSPECTOR_ROOT_DIR)/lib64/libittnotify.a

# Optimization report
else ifneq (,$(call parse_config,opt-report))
    CXXFLAGS += -O3 -qopt-report5

# Default configuration
else
    CXXFLAGS += -O3 -g #-xHost -no-prec-div -ipo
endif

# Manage options in the "config" parameter
ifneq (,$(call parse_config,detailed_timers))
    CXXFLAGS += -D__DETAILED_TIMERS
endif

#activate openmp unless noopenmp flag
ifeq (,$(call parse_config,noopenmp))
    OPENMP_FLAG ?= -fopenmp
    LDFLAGS += -lm
    OPENMP_FLAG += -D_OMP
    LDFLAGS += $(OPENMP_FLAG)
    CXXFLAGS += $(OPENMP_FLAG)
endif

ifneq (,$(call parse_config,picsar))
    # New environment variable
	FFTW3_LIB ?= $(FFTW_LIB_DIR)
	LIBPXR ?= picsar/lib
	# Set Picsar link environment
	CXXFLAGS += -D_PICSAR
	LDFLAGS += -L$(LIBPXR) -lpxr
	LDFLAGS += -L$(FFTW3_LIB) -lfftw3_mpi  -lopenblas

	LDFLAGS += -L$(FFTW3_LIB) -lfftw3_threads
	LDFLAGS += -L$(FFTW3_LIB) -lfftw3
	LDFLAGS += -lgfortran
endif


# Manage MPI communications by a single thread (master in MW)
ifneq (,$(call parse_config,no_mpi_tm))
    CXXFLAGS += -D_NO_MPI_TM
endif

#last: check remaining arguments and raise error
ifneq ($(strip $(my_config)),)
$(error "Unused parameters in config : $(my_config)")
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

check:
	$(Q) $(PYTHONEXE) scripts/compile_tools/check_make_options.py config $(config)
	$(Q) $(PYTHONEXE) scripts/compile_tools/check_make_options.py machine $(machine)

# Create python header files
$(BUILD_DIR)/%.pyh: %.py
	@echo "Creating binary char for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(PYTHONEXE) scripts/compile_tools/hexdump.py "$<" "$@"

# Calculate dependencies
$(BUILD_DIR)/%.d: %.cpp
	@echo "Checking dependencies for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(SMILEICXX) $(CXXFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

$(BUILD_DIR)/src/Diagnostic/DiagnosticScalar.o : src/Diagnostic/DiagnosticScalar.cpp
	@echo "SPECIAL COMPILATION FOR $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -O1 -c $< -o $@

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
# astyle
style:
	@echo "Astyle is applied on all files"
	$(Q) astyle --style=1tbs --fill-empty-lines --pad-comma --unpad-paren --pad-paren-in --align-pointer=name --align-reference=name -n -r src/*.cpp,*.h

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
# Smilei tables

TABLES_EXEC = smilei_tables

tables: tables_folder $(TABLES_EXEC)
	
tables_folder:
	@echo "Installing smilei_tables tool"
	@mkdir -p $(TABLES_BUILD_DIR)

tables_clean:
	@echo "Cleaning $(TABLES_BUILD_DIR)"
	@rm -r $(TABLES_BUILD_DIR)

# Calculate dependencies
$(TABLES_BUILD_DIR)/%.d: %.cpp
	@echo "Checking dependencies for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(SMILEICXX) $(CXXFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

# Compile cpps
$(TABLES_BUILD_DIR)/%.o : $(TABLES_DIR)/%.cpp
	@echo "Compiling $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -c $< -o $@

# Link the main program
$(TABLES_EXEC): $(TABLES_OBJS)
	@echo "Linking $@"
	$(Q) $(SMILEICXX) $(TABLES_OBJS) -o $(TABLES_BUILD_DIR)/$@ $(LDFLAGS)
	$(Q) cp $(TABLES_BUILD_DIR)/$@ $@

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
	@echo '    detailed_timers      : to compile the code with more refined timers (refined time report)'
	@echo '    noopenmp             : to compile without openmp'
	@echo '    no_mpi_tm            : to compile with a MPI library without MPI_THREAD_MULTIPLE support'
	@echo '    opt-report           : to generate a report about optimization, vectorization and inlining (Intel compiler)'
	@echo '    scalasca             : to compile using scalasca'
	@echo '    advisor              : to compile for Intel Advisor analysis'
	@echo '    vtune                : to compile for Intel Vtune analysis'
	@echo '    inspector            : to compile for Intel Inspector analysis'
	@echo
	@echo 'Examples:'
	@echo '  make config=verbose'
	@echo '  make config=debug'
	@echo '  make config="debug noopenmp"'
	@echo
	@echo 'Machine options:'
	@echo '  make machine=XXX      : include machine file in scripts/compile_tools/machine/XXX'
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
	@echo 'SMILEI TABLES:'
	@echo '---------------'
	@echo '  make tables           : compilation of the tool smilei_tables'
	@echo ''
	@echo 'Environment variables:'
	@echo '  SMILEICXX             : mpi c++ compiler [$(SMILEICXX)]'
	@echo '  HDF5_ROOT_DIR         : HDF5 dir. Defaults to the value of HDF5_ROOT [$(HDF5_ROOT_DIR)]'
	@echo '  BUILD_DIR             : directory used to store build files [$(BUILD_DIR)]'
	@echo '  OPENMP_FLAG           : openmp flag [$(OPENMP_FLAG)]'
	@echo '  PYTHONEXE             : python executable [$(PYTHONEXE)]'
	@echo '  FFTW3_LIB             : FFTW3 libraries directory [$(FFTW3_LIB)]'
	@echo '  LIB PXR               : Picsar library directory [$(LIBPXR)]'
	@echo
	@echo 'Intel Inspector environment:'
	@echo '  INSPECTOR_ROOT_DIR    : only needed to use the inspector API (__itt functions) [$(INSPECTOR_ROOT_DIR)]'
	@echo
	@echo 'http://www.maisondelasimulation.fr/smilei'
	@echo 'https://github.com/SmileiPIC/Smilei'
	@echo
	@if [ -f scripts/compile_tools/machine/$(machine) ]; then echo "Machine comments for $(machine):"; grep '^#' scripts/compile_tools/machine/$(machine) || echo "None"; else echo "Available machines:"; ls -1 scripts/compile_tools/machine; fi
