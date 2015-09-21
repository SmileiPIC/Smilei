SMILEICXX     ?= mpic++
HDF5_ROOT_DIR ?=

PYTHONCONFIG ?= python-config

EXEC = ../smilei

default: $(EXEC)

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

BUILD ?= ../build
	

#collect all cpp files
SRCS := $(shell find * -name \*.cpp)
OBJS := $(addprefix $(BUILD)/, $(SRCS:.cpp=.o))
DEPS := $(addprefix $(BUILD)/, $(SRCS:.cpp=.d))
PYSCRIPTS := $(shell find Python -name \*.py)
CFLAGS += -I$(BUILD)/Python
PYHEADERS := $(addprefix $(BUILD)/, $(PYSCRIPTS:.py=.pyh))

CFLAGS+=$(shell $(PYTHONCONFIG) --includes)
LDFLAGS+=$(shell $(PYTHONCONFIG) --ldflags)


# check for variable config
ifneq (,$(findstring debug,$(config)))
	CFLAGS += -g -pg -Wall -D__DEBUG -O0# -shared-intel 
else
	CFLAGS += -O3#  -xHost -ipo
endif

ifneq (,$(findstring scalasca,$(config)))
    SMILEICXX = scalasca -instrument mpic++
endif

ifneq (,$(findstring openmp,$(config)))
    ifneq (,$(findstring icpc,$(shell $(SMILEICXX) --showme:command)))
        CFLAGS += -openmp
    else
        CFLAGS += -fopenmp 
    endif
    CFLAGS += -D_OMP
endif

clean:
	rm -f $(OBJS) $(DEPS) $(PYHEADERS)
	rm -rf $(BUILD) 
	
distclean: clean
	rm -f $(EXEC)
	
	
# this generates a .h file containing a char[] with the python script in binary then
#you can just include this file to get the contents
$(BUILD)/%.pyh: %.py
	@ echo "Creating binary char for $< : $@"
	@ mkdir -p $(@D)
	@ cd $(<D) && xxd -i $(<F) > ../$@

$(BUILD)/%.d: %.cpp
	@ echo "Checking dependencies for $<"
# create and modify dependecy file .d to take into account the location subdir
	@ $(SMILEICXX) $(CFLAGS) -MM $< 2>/dev/null | sed -e "s@\(^.*\)\.o:@$(BUILD)/$(shell  dirname $<)/\1.d $(BUILD)/$(shell  dirname $<)/\1.o:@" > $@  

$(BUILD)/%.o : %.cpp
	$(SMILEICXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(SMILEICXX) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# these are kept for backward compatibility and might be removed (see make help)
obsolete:
	@echo "[WARNING] Please consider using make config=\"$(MAKECMDGOALS)\""

debug: obsolete
	make config=debug

openmp: obsolete
	make config=openmp

scalasca: obsolete
	make config=scalasca


ifeq ($(filter clean,$(MAKECMDGOALS)),) 
# Let's try to make the next lines clear: we include $(DEPS) and pygenerator
-include $(DEPS) pygenerator
# we specify that pygenerator is not a file
.PHONY : pygenerator
# create the tree to store .d .o .pyh files
buildtree:
	@echo "Creating build dirtree in $(BUILD)"
	@mkdir -p $(addprefix $(BUILD)/, $(DIRS))
# and pygenerator will create all the $(PYHEADERS) (which are files)
pygenerator : buildtree $(PYHEADERS)
endif



# DO NOT DELETE
