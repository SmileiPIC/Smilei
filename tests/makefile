INPUT = $(wildcard *.in)
OUTPUT = $(patsubst %.in,%.out,$(INPUT))

.PHONY : all

all: $(OUTPUT)

%.out : %.in
	../smilei $< > $@ 2>&1
	
clean:
	rm -f $(OUTPUT)