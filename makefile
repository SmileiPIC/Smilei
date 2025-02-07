#
# ----------------------------
# >>>>>>>> READ THIS <<<<<<<<<
# ----------------------------
#
# Instead of modifying this file, you should set the following
# environment variables on your system
#
# BUILD_DIR        : the path to the build directory (default: ./build)
# SMILEICXX        : the MPI C++ executable (for instance mpicxx, mpiicpc, etc.)
# PYTHONEXE        : the python executable to be used in smilei
# HDF5_ROOT_DIR    : the local path to the HDF5 library
# BOOST_ROOT_DIR   : the local path to the boost library
# TABLES_BUILD_DIR : build directory for databases (default ./tools/tables/build)

BUILD_DIR ?= build
SMILEICXX ?= mpicxx
PYTHONEXE ?= python
HDF5_ROOT_DIR ?= $(HDF5_ROOT)
BOOST_ROOT_DIR ?= $(BOOST_ROOT)
TABLES_BUILD_DIR ?= tools/tables/build

#-----------------------------------------------------
# Machines scripts may need that
my_config:=$(config)
define parse_config
$(findstring $(1),$(config))$(eval my_config:=$(filter-out $(1),$(my_config)))
endef

PYTHONCONFIG := $(PYTHONEXE) scripts/compile_tools/python-config.py

#-----------------------------------------------------
# Git information
VERSION:=$(shell $(PYTHONEXE) scripts/compile_tools/get-version.py )

#-----------------------------------------------------
# Compiler specific flags

COMPILER_INFO := $(shell $(SMILEICXX) -show | cut -d' ' -f1)

ifeq ($(findstring g++, $(COMPILER_INFO)), g++)
	CXXFLAGS += -Wno-reorder
else ifeq ($(findstring clang++, $(COMPILER_INFO)), clang++)
	CXXFLAGS += -Wdeprecated-register
endif

#-----------------------------------------------------
# Directories and files

# Smilei
DIRS := $(shell find src -type d)
SRCS := $(shell find src/* -name \*.cpp)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.d))
SITEDIR = $(shell d=`$(PYTHONEXE) -m site --user-site` && echo $$d || $(PYTHONEXE) -c "import sysconfig; print(sysconfig.get_path('purelib'))")

# Smilei tools
TABLES_DIR := tools/tables
TABLES_SRCS := $(shell find tools/tables/* -name \*.cpp | rev | cut -d '/' -f1 | rev)
TABLES_DEPS := $(addprefix $(TABLES_BUILD_DIR)/, $(SRCS:.cpp=.d))
TABLES_OBJS := $(addprefix $(TABLES_BUILD_DIR)/, $(TABLES_SRCS:.cpp=.o))
TABLES_SRCS := $(shell find tools/tables/* -name \*.cpp)

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

#-----------------------------------------------------
# Flags

# Smilei version
CXXFLAGS += -D__VERSION=\"$(VERSION)\"
# Remove OpenMPI warnings
CXXFLAGS += -DOMPI_SKIP_MPICXX
# C++ version
ifeq ($(findstring armclang++, $(COMPILER_INFO)), armclang++)
	CXXFLAGS += -std=c++11 -Wall
else ifeq ($(findstring clang++, $(COMPILER_INFO)), clang++)
	CXXFLAGS += -std=c++11 -Wall -Wno-unused-command-line-argument
else ifeq ($(findstring g++, $(COMPILER_INFO)), g++)
	CXXFLAGS += -std=c++11 -Wall -Wextra
else ifeq ($(findstring FCC, $(COMPILER_INFO)), FCC)
	CXXFLAGS += -std=c++11
else ifeq ($(findstring FCC, $(COMPILER_INFO)), FCCpx)
	CXXFLAGS += -std=c++11
else
	CXXFLAGS += -std=c++14 #-Wall #not recognized by nvcc, make an exception
endif

# HDF5 library
ifneq ($(strip $(HDF5_ROOT_DIR)),)
CXXFLAGS += -I$(HDF5_ROOT_DIR)/include
LDFLAGS := -L$(HDF5_ROOT_DIR)/lib  $(LDFLAGS)
DEPSFLAGS += -I$(HDF5_ROOT_DIR)/include
endif
# Boost library
ifneq ($(strip $(BOOST_ROOT_DIR)),)
CXXFLAGS += -I$(BOOST_ROOT_DIR)/include
LDFLAGS := -L$(BOOST_ROOT_DIR)/lib $(LDFLAGS)
endif
LDFLAGS += -lhdf5
# Include subdirs
CXXFLAGS += $(DIRS:%=-I%)
DEPSFLAGS += $(DIRS:%=-I%)
# Python-related flags
CXXFLAGS += -I$(BUILD_DIR)/src/Python
DEPSFLAGS += -I$(BUILD_DIR)/src/Python

PYSCRIPTS = $(shell find src/Python -name \*.py)
PYHEADERS := $(addprefix $(BUILD_DIR)/, $(PYSCRIPTS:.py=.pyh))
PY_CXXFLAGS := $(shell $(PYTHONCONFIG) --includes)
CXXFLAGS += $(PY_CXXFLAGS)
DEPSFLAGS += $(PY_CXXFLAGS)

PY_LDFLAGS := $(shell $(PYTHONCONFIG) --ldflags)
LDFLAGS += $(PY_LDFLAGS)
ifneq ($(strip $(PYTHONHOME)),)
	LDFLAGS += -L$(PYTHONHOME)/lib
endif

# Manage options in the "config" parameter
ifneq (,$(call parse_config,debug))
	CXXFLAGS += -g -pg -D__DEBUG -O0
# With gdb
else ifneq (,$(call parse_config,gdb))
	CXXFLAGS += -g -D__DEBUG -O0
# With gdb
else ifneq (,$(call parse_config,ddt))
	# -g
	CXXFLAGS += -O0 -g
# With valgrind
else ifneq (,$(call parse_config,valgrind))
	CXXFLAGS += -g -O3

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

# Default configuration
else
	ifeq ($(findstring clang++, $(COMPILER_INFO)), clang++)
		CXXFLAGS += -O3 -g -fno-math-errno
	else ifeq ($(findstring armclang++, $(COMPILER_INFO)), armclang++)
		CXXFLAGS += -Ofast -g
	else ifeq ($(findstring FCC, $(COMPILER_INFO)), FCC)
		CXXFLAGS += -O3 -Kfast -g
	else ifeq ($(findstring FCCpx, $(COMPILER_INFO)), FCCpx)
		CXXFLAGS += -O3 -Kfast -g
	else ifeq ($(findstring pgi, $(COMPILER_INFO)), pgi)
		CXXFLAGS += -O3
	else
		CXXFLAGS += -O3 -g
	endif
endif

# Optimization report
ifneq (,$(call parse_config,opt-report))
	# Clang compiler
	ifeq ($(findstring clang++, $(COMPILER_INFO)), clang++)
		CXXFLAGS += -fsave-optimization-record -Rpass-analysis=loop-vectorize
	else ifeq ($(findstring armclang++, $(COMPILER_INFO)), armclang++)
		CXXFLAGS += -fsave-optimization-record -Rpass-analysis=loop-vectorize
	else ifeq ($(findstring FCC, $(COMPILER_INFO)), FCC)
		CXXFLAGS += -Koptmsg=2 -Nlst=t
	else ifeq ($(findstring FCCpx, $(COMPILER_INFO)), FCCpx)
		CXXFLAGS += -Koptmsg=2 -Nlst=t
	else ifeq ($(findstring g++, $(COMPILER_INFO)), g++)
		CXXFLAGS += -fopt-info
	# Intel compiler
	else ifeq ($(findstring icpc, $(COMPILER_INFO)), icpc)
		CXXFLAGS += -qopt-report5
	endif
endif

# Detailed timers
ifneq (,$(call parse_config,detailed_timers))
	CXXFLAGS += -D__DETAILED_TIMERS
endif

# NVIDIA GPUs
ifneq (,$(call parse_config,gpu_nvidia))
	override config += noopenmp # Prevent openmp for nvidia
	
	CXXFLAGS += -DSMILEI_ACCELERATOR_GPU -DSMILEI_ACCELERATOR_GPU_OACC
	GPU_COMPILER ?= nvcc
	GPU_COMPILER_FLAGS += -x cu -DSMILEI_ACCELERATOR_GPU -DSMILEI_ACCELERATOR_GPU_OACC $(DIRS:%=-I%)
	GPU_COMPILER_FLAGS += -I$(BUILD_DIR)/src/Python $(PY_CXXFLAGS)
	GPU_KERNEL_SRCS := $(shell find src/* -name \*.cu)
	GPU_KERNEL_OBJS := $(addprefix $(BUILD_DIR)/, $(GPU_KERNEL_SRCS:.cu=.o))
	
	OBJS += $(GPU_KERNEL_OBJS)
endif

# AMD GPUs
ifneq (,$(call parse_config,gpu_amd))
	CXXFLAGS += -DSMILEI_ACCELERATOR_GPU -DSMILEI_ACCELERATOR_GPU_OMP
	GPU_COMPILER ?= $(CC)
	GPU_COMPILER_FLAGS += -x hip -DSMILEI_ACCELERATOR_GPU -DSMILEI_ACCELERATOR_GPU_OMP -std=c++14 $(DIRS:%=-I%)
	GPU_COMPILER_FLAGS += -I$(BUILD_DIR)/src/Python $(PY_CXXFLAGS)
	GPU_KERNEL_SRCS := $(shell find src/* -name \*.cu)
	GPU_KERNEL_OBJS := $(addprefix $(BUILD_DIR)/, $(GPU_KERNEL_SRCS:.cu=.o))
	
	OBJS += $(GPU_KERNEL_OBJS)
endif

#activate openmp unless noopenmp flag
# For Fujitsu compiler: -Kopenmp
ifeq (,$(call parse_config,noopenmp))
	ifeq ($(findstring FCC, $(COMPILER_INFO)), FCC)
		OPENMP_FLAG ?= -Kopenmp -Kopenmp_simd
	else ifeq ($(findstring FCCpx, $(COMPILER_INFO)), FCCpx)
		OPENMP_FLAG ?= -Kopenmp -Kopenmp_simd
	else
		OPENMP_FLAG ?= -fopenmp
	endif
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
	#LDFLAGS += -lgfortran
endif


# Manage MPI communications by a single thread (master in MW)
ifneq (,$(call parse_config,no_mpi_tm))
	CXXFLAGS += -D_NO_MPI_TM
endif

ifneq (,$(call parse_config,part_event_tracing))
	CXXFLAGS += -D_PARTEVENTTRACING
endif

CXXFLAGS0 = $(shell echo $(CXXFLAGS)| sed "s/O3/O0/g")

#-----------------------------------------------------
# Set the verbosity prefix
ifeq (,$(call parse_config,verbose))
	Q := @
else
	Q :=
endif

#last: check remaining arguments and raise error
ifneq ($(strip $(my_config)),)
$(error "Unused parameters in config : $(my_config)")
endif

SMILEICXX_DEPS ?= $(SMILEICXX)

#-----------------------------------------------------
# Rules for building the excutable smilei

EXEC = smilei

default: $(PYHEADERS) $(EXEC) $(EXEC)_test

#-----------------------------------------------------
# Header
header:
	@echo " _____________________________________"
	@echo ""
	@echo " SMILEI compilation"
	@echo ""
	@if [ $(call parse_config,debug) ]; then echo "- Debug option requested"; fi;
	@if [ $(call parse_config,gdb) ]; then echo "- Compilation for GDB requested"; fi;
	@if [ $(call parse_config,picsar) ]; then echo "- SMILEI linked to PICSAR requested"; fi;
	@if [ $(call parse_config,opt-report) ]; then echo "- Optimization report requested"; fi;
	@if [ $(call parse_config,detailed_timers) ]; then echo "- Detailed timers option requested"; fi;
	@if [ $(call parse_config,no_mpi_tm) ]; then echo "- Compiled without MPI_THREAD_MULTIPLE"; fi;
	@if [ $(call parse_config,part_event_tracing) ]; then echo "- Compiled with particle events tracing"; fi;
	@echo " _____________________________________"
	@echo ""

clean:
	@echo "Cleaning $(BUILD_DIR)"
	$(Q) rm -rf $(EXEC)
	$(Q) rm -rf $(EXEC)_test
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
	$(Q) $(SMILEICXX_DEPS) $(DEPSFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

# Calculate dependencies: special for Params.cpp which needs pyh files
$(BUILD_DIR)/src/Params/Params.d: src/Params/Params.cpp $(PYHEADERS)
	@echo "Checking dependencies for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(SMILEICXX_DEPS) $(DEPSFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

ifeq ($(findstring icpc, $(COMPILER_INFO)), icpc)
$(BUILD_DIR)/src/Diagnostic/DiagnosticScalar.o : src/Diagnostic/DiagnosticScalar.cpp
	@echo "SPECIAL COMPILATION FOR $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -O1 -c $< -o $@
endif

$(BUILD_DIR)/src/MultiphotonBreitWheeler/MultiphotonBreitWheelerTablesDefault.o : src/MultiphotonBreitWheeler/MultiphotonBreitWheelerTablesDefault.cpp
	@echo "SPECIAL COMPILATION FOR $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS0) -c $< -o $@

$(BUILD_DIR)/src/Radiation/RadiationTablesDefault.o : src/Radiation/RadiationTablesDefault.cpp
	@echo "SPECIAL COMPILATION FOR $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS0) -c $< -o $@

# Compile cpps
$(BUILD_DIR)/%.o : %.cpp
	@echo "Compiling $<"
	$(Q) $(SMILEICXX) $(CXXFLAGS) -c $< -o $@

# Compile cus
$(BUILD_DIR)/%.o : %.cu
	@echo "Compiling $<"
	$(Q) $(GPU_COMPILER) $(GPU_COMPILER_FLAGS) -c $< -o $@

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

# these are not file-related rules
PHONY_RULES=clean distclean help env debug doc tar happi uninstall_happi
.PHONY: $(PHONY_RULES)

# Check dependencies only when necessary
GOALS = $(if $(MAKECMDGOALS), $(MAKECMDGOALS), default)
ifneq ($(filter-out $(PHONY_RULES) print-%, $(GOALS)),)
	-include $(DEPS)
endif

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

env:  print-VERSION print-SMILEICXX print-OPENMP_FLAG print-HDF5_ROOT_DIR print-FFTW3_LIB_DIR print-SITEDIR print-PYTHONEXE print-PY_CXXFLAGS print-PY_LDFLAGS print-CXXFLAGS print-LDFLAGS print-GPU_COMPILER print-GPU_COMPILER_FLAGS print-COMPILER_INFO

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

# TODO(Etienne M): This does not get the dependencies of the .cu files !
# Calculate dependencies
$(TABLES_BUILD_DIR)/%.d: %.cpp
	@echo "Checking dependencies for $<"
	$(Q) if [ ! -d "$(@D)" ]; then mkdir -p "$(@D)"; fi;
	$(Q) $(SMILEICXX) $(DPSFLAGS) -MF"$@" -MM -MP -MT"$@ $(@:.d=.o)" $<

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
	@echo 'or, to compile with 4 threads (for instance):'
	@echo '  make -j 4'
	@echo
	@echo 'Config options:'
	@echo '  make config="[ verbose ] [ debug ] [ scalasca ] [ noopenmp ]"'
	@echo '    verbose                      : to print compile command lines'
	@echo '    noopenmp                     : to compile without openmp'
	@echo '    no_mpi_tm                    : to compile with a MPI library without MPI_THREAD_MULTIPLE support'
	@echo '    gpu_nvidia                   : to compile for NVIDIA GPU (uses OpenACC)'
	@echo '    gpu_amd                      : to compile for AMP GPU (uses OpenMP)'
	@echo '    detailed_timers              : to compile the code with more refined timers (refined time report)'
	@echo '    debug                        : to compile in debug mode (code runs really slow)'
	@echo '    opt-report                   : to generate a report about optimization, vectorization and inlining (Intel compiler)'
	@echo '    scalasca                     : to compile using scalasca'
	@echo '    advisor                      : to compile for Intel Advisor analysis'
	@echo '    vtune                        : to compile for Intel Vtune analysis'
	@echo '    inspector                    : to compile for Intel Inspector analysis'
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
	@echo 'Environment variables needed for compilation:'
	@echo '  SMILEICXX         : mpi c++ compiler (possibly GPU-aware) [mpicxx]'
	@echo '  SMILEICXX_DEPS    : c++ compiler for calculating dependencies [$$SMILEICXX]'
	@echo '  CXXFLAGS          : FLAGS for $$SMILEICXX []'
	@echo '  LDFLAGS           : FLAGS for the linker []'
	@echo '  HDF5_ROOT_DIR     : folder where the HDF5 library was installed [$$HDF5_ROOT_DIR]'
	@echo '  BUILD_DIR         : custom folder for building Smilei [build]'
	@echo '  PYTHONEXE         : python executable [python]'
	@echo '  GPU_COMPILER      : compiler for cuda-like files [$$CC]'
	@echo '  GPU_COMPILER_FLAGS: flags for the $$GPU_COMPILER []'
#    @echo '  FFTW3_LIB_DIR  : FFTW3 libraries directory [$(FFTW3_LIB_DIR)]'
#    @echo '  LIBPXR         : Picsar library directory [$(LIBPXR)]'
	@echo
	@echo 'Intel Inspector environment:'
	@echo '  INSPECTOR_ROOT_DIR    : only needed to use the inspector API (__itt functions) [$(INSPECTOR_ROOT_DIR)]'
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
	@echo 
	@echo 'https://smileipic.github.io/Smilei/'
	@echo 'https://github.com/SmileiPIC/Smilei'
	@echo
	@if [ -f scripts/compile_tools/machine/$(machine) ]; then echo "Machine comments for $(machine):"; grep '^#' scripts/compile_tools/machine/$(machine) || echo "None"; else echo "Available machines:"; ls -1 scripts/compile_tools/machine; fi
