VERSION:="v$(shell git rev-list HEAD --count):$(shell git log -1 --pretty=format:%h)"

CC      = mpic++
#CC      = g++
CFLAGS  = -D__VERSION=\"$(VERSION)\" #-Wstrict-aliasing
LDFLAGS = -lm 

all: release install

INSTALLDIR = ..
EXEC = smilei

SRCS = $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

release: CFLAGS += -O3
release: $(EXEC)

debug: CFLAGS += -g -pg -Wall -D__DEBUG
debug: $(EXEC)

.depend: $(SRCS)
	rm -rf .depend
	$(CC) $(CFLAGS) -MM $^ >> .depend;

include .depend 

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJS) $(LDFLAGS)

install: $(EXEC)
	install $(EXEC) $(INSTALLDIR)
	
clean:
	rm -f $(EXEC) $(OBJS) .depend


