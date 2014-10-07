VERSION:="$(shell git describe 2>/dev/null || echo '??')"
CFLAGS += -D__VERSION=\"$(VERSION)\" -I${HDF5_ROOT_DIR}/include
LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5

SMILEICXX     ?= mpic++
HDF5_ROOT_DIR ?=

ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz
endif

DIRS = ./ $(wildcard */)
CFLAGS +=  $(DIRS:%=-I%) # add includes flags for each dir

##################################################

all: release

EXEC = smilei

SRCS := $(foreach dir,$(DIRS),$(wildcard $(dir)*.cpp))
OBJS := $(patsubst %.cpp,%.o,$(SRCS))
DEPS := $(patsubst %.cpp,%.d,$(SRCS))

release: CFLAGS += -O3
release: $(EXEC)

openmpgnu: CFLAGS += -O3 -fopenmp -D_OMP
openmpgnu: $(EXEC)

openmpintel: CFLAGS += -O3 -openmp -D_OMP
openmpintel: $(EXEC)

debug: CFLAGS += -g -pg -Wall -D__DEBUG
debug: $(EXEC)

-include $(OBJS:.o=.d)

%.d: %.cpp
	@ echo "Checking dependencies for $<"
	@../scripts/depend.sh $(shell dirname $<) $(SMILEICXX) -MM -MG $(CFLAGS) $*.cpp > $@ 2>/dev/null

%.o : %.cpp
	$(SMILEICXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(SMILEICXX) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

clean:
	rm -f $(EXEC) $(OBJS) $(DEPS) .depend
	
# DO NOT DELETE
