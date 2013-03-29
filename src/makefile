VERSION:="v$(shell git rev-list HEAD --count):$(shell git log -1 --pretty=format:%h)"

CC      = g++
CFLAGS  = -D__VERSION=\"$(VERSION)\"
LDFLAGS = -lm 

all: release
EXEC=../spice

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

clean:
	rm -f $(EXEC) $(OBJS) .depend


