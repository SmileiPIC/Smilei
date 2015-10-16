SMILEICXX     ?= mpic++
HDF5_ROOT_DIR ?=

BUILD_DIR ?= build

PYTHONCONFIG ?= python-config

EXEC = smilei

default: $(EXEC)

####################################################
DESCRIBE:=$(shell git describe 2>/dev/null || echo '??')
BRANCH:=$(shell git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '??')
VERSION="$(DESCRIBE)-$(BRANCH)"
COMMITDATE:="$(shell git show -s --pretty="%ci" 2>/dev/null || echo '??')"

CXXFLAGS += -D__VERSION=\"$(VERSION)\" -D__COMMITDATE=\"$(COMMITDATE)\" -D__CONFIG=\""$(config)"\"

CXXFLAGS += -I${HDF5_ROOT_DIR}/include -std=c++0x 
LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5 -lz


ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz -L/gpfslocal/pub/python/anaconda/Anaconda-2.1.0/lib
endif

#add subdirs
DIRS := $(shell find src -type d)
#add include directives for subdirs
CXXFLAGS += $(DIRS:%=-I%)

#collect all cpp files
SRCS := $(shell find src/* -name \*.cpp)
OBJS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.o))
DEPS := $(addprefix $(BUILD_DIR)/, $(SRCS:.cpp=.d))
PYSCRIPTS := $(shell find src/Python -name \*.py)
CXXFLAGS += -I$(BUILD_DIR)/src/Python
PYHEADERS := $(addprefix $(BUILD_DIR)/, $(PYSCRIPTS:.py=.pyh))

PY_CXXFLAGS:=$(shell $(PYTHONCONFIG) --includes)
CXXFLAGS+=$(PY_CXXFLAGS)
PY_LDFLAGS:=$(shell $(PYTHONCONFIG) --ldflags)
LDFLAGS+=$(PY_LDFLAGS)


# check for variable config
ifneq (,$(findstring debug,$(config)))
	CXXFLAGS += -g -pg -Wall -D__DEBUG -O0 # -shared-intel 
else
	CXXFLAGS += -O3 # -xHost -ipo
endif

ifneq (,$(findstring scalasca,$(config)))
    SMILEICXX = scalasca -instrument mpic++
endif

ifneq (,$(findstring turing,$(config)))
	CXXFLAGS += -I$(BG_PYTHONHOME)/include/python2.7 -qlanglvl=extended0x
	LDFLAGS  += -qnostaticlink -L(BG_PYTHONHOME)/lib64 -lpython2.7 -lutil
endif

ifeq (,$(findstring noopenmp,$(config)))
	SMILEI_COMPILER:=$(shell $(SMILEICXX) --showme:command)
    ifneq (,$(findstring icpc,$(SMILEI_COMPILER)))
        OPENMPFLAGS = -openmp
    else
        OPENMPFLAGS = -fopenmp 
    endif
    OPENMPFLAGS += -D_OMP
    LDFLAGS += $(OPENMPFLAGS)
    CXXFLAGS += $(OPENMPFLAGS)
endif

clean:
	rm -f $(OBJS) $(DEPS) $(PYHEADERS)
	rm -rf $(BUILD_DIR) 
	rm -rf smilei-$(VERSION).tgz
	make -C doc clean
	
distclean: clean
	rm -f $(EXEC)
	
	
# this generates a .h file containing a char[] with the python script in binary then
# you can just include this file to get the contents (in Params/Params.cpp)
$(BUILD_DIR)/%.pyh: %.py
	@ echo "Creating binary char for $< : $@"
	@ mkdir -p "$(@D)" 
	@ cd "$(<D)" && xxd -i "$(<F)" > "$(@F)"
	@ mv "$(<D)/$(@F)" "$@"

$(BUILD_DIR)/%.d: %.cpp
	@ echo "Checking dependencies for $<"
# create and modify dependecy file .d to take into account the location subdir
	@ $(SMILEICXX) $(CXXFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(BUILD_DIR)/$(shell  dirname $<)/\1.d $(BUILD_DIR)/$(shell  dirname $<)/\1.o:@" > $@  

$(BUILD_DIR)/%.o : %.cpp
	$(SMILEICXX) $(CXXFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(SMILEICXX) $(OBJS) -o $(BUILD_DIR)/$@  $(LDFLAGS)
	cp $(BUILD_DIR)/$@ $@

# these are kept for backward compatibility and might be removed (see make help)
obsolete:
	@echo "[WARNING] Please consider using make config=\"$(MAKECMDGOALS)\""

debug: obsolete
	make config=debug

scalasca: obsolete
	make config=scalasca


ifeq ($(filter clean help doc,$(MAKECMDGOALS)),) 
# Let's try to make the next lines clear: we include $(DEPS) and pygenerator
-include $(DEPS) pygenerator
# we specify that pygenerator is not a file
.PHONY : pygenerator buildtree
# create the tree to store .d .o .pyh files
buildtree:
	@echo "Creating build dirtree in $(BUILD_DIR)"
	@mkdir -p $(addprefix $(BUILD_DIR)/, $(DIRS))
# and pygenerator will create all the $(PYHEADERS) (which are files)
pygenerator : buildtree $(PYHEADERS)
endif


.PHONY: doc src help clean default tar

doc:
	make -C doc all

tar:
	git archive -o smilei-$(VERSION).tgz --prefix smilei-$(VERSION)/ HEAD
	
help: 
	@echo 'Usage: make config=OPTIONS'
	@echo '	    OPTIONS is a string composed of one or more of:'
	@echo '	        debug      : to compile in debug mode (code runs really slow)'
	@echo '         scalasca   : to compile using scalasca'
	@echo '         noopenmp   : to compile without openmp'
	@echo ' examples:'
	@echo '     make config=debug'
	@echo '     make config=noopenmp'
	@echo '     make config="debug noopenmp"'
	@echo 
	@echo 'Environment variables :'
	@echo '     SMILEICXX     : mpi c++ compiler'
	@echo '     HDF5_ROOT_DIR : HDF5 dir'
	@echo '     BUILD_DIR         :directory used to store build files (default: "build")'
	@echo 
	@echo '      make doc     : builds the documentation'
	@echo '      make tar     : creates an archive of the sources'
	@echo '      make clean   : remove build'

