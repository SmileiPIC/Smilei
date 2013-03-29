INPUT = $(wildcard *.in)
OUTPUT = $(patsubst %.in,%.out,$(INPUT))

.PHONY : all

all: $(OUTPUT)

%.out : %.in
	../spice $< > $@ 2>&1
	
clean:
	rm -f $(OUTPUT)