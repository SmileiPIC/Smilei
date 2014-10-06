VERSION:="$(shell git describe 2>/dev/null || echo '??')"

SMILEICXX     ?= mpic++
HDF5_ROOT_DIR ?=


ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz
endif

SUBDIRS = ./ $(wildcard */)
INC_SUB = $(foreach d, $(SUBDIRS), -I./$d)
override CFLAGS  += -D__VERSION=\"$(VERSION)\" $(INC_SUB) -I${HDF5_ROOT_DIR}/include
LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5

##################################################

all: release

EXEC = smilei

SRCS := $(foreach dir,$(SUBDIRS),$(wildcard $(dir)*.cpp))
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

release: CFLAGS += -O3
release: $(EXEC)

openmpgnu: CFLAGS += -O3 -fopenmp -D_OMP
openmpgnu: $(EXEC)

openmpintel: CFLAGS += -O3 -openmp -D_OMP
openmpintel: $(EXEC)

debug: CFLAGS += -g -pg -Wall -D__DEBUG
debug: $(EXEC)

.depend: $(SRCS)
	$(SMILEICXX) $(CFLAGS) -MM $^ > .depend

include .depend 

%.o : %.cpp
	$(SMILEICXX) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(SMILEICXX) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

clean:
	rm -f $(EXEC) $(OBJS) .depend
