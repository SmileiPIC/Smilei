VERSION:="$(shell git describe 2>/dev/null || echo '??')"

SMILEICXX     ?= mpic++
HDF5_ROOT_DIR ?=

CFLAGS  += -D__VERSION=\"$(VERSION)\" -I${HDF5_ROOT_DIR}/include #-Wstrict-aliasing
LDFLAGS += -lm -L${HDF5_ROOT_DIR}/lib -lhdf5

ifneq (,$(findstring poincare,$(HOSTNAME)))
    LDFLAGS += -lgpfs -lz
endif



all: release

EXEC = smilei

SRCS = $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

release: CFLAGS += -O3
release: $(EXEC)

openmpgnu: CFLAGS += -O3 -fopenmp
openmpgnu: $(EXEC)

openmpintel: CFLAGS += -O3 -openmp
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
