INPUT = $(wildcard *.in)

PROCS ?= ''

all: $(INPUT)

FORCE :

%.in : FORCE
	../scripts/smilei.sh $(PROCS) $@
	
clean:
	rm -f $(OUTPUT)