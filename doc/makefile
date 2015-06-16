VERSION:="$(shell git describe 2>/dev/null || echo '??')-$(shell git rev-parse --abbrev-ref HEAD || echo '??')"

TEXFILES=Collisions/Collisions.tex Diagnostics/Diagnostics.tex Profiles/Profiles.tex

PDFS:=$(TEXFILES:.tex=.pdf)

.PHONY: html

default: $(PDFS) 

all: $(PDFS) html

%.pdf: %.tex
	latexmk -cd -pdf $< && latexmk -cd -c $<

html:
#	convert -resize 100x100 doc/smileiLogo/smileiLogo.png doc/logo.png
	cd doc; (echo "PROJECT_NUMBER=${VERSION}"; cat smilei.dox) | doxygen -
#   pdf generation is buggy... dropped
#	latexmk -cd -pdf latex/refman.tex && latexmk -cd -c latex/refman.tex
	
clean:
	for fname in $(TEXFILES); do latexmk -cd -C $$fname; done
	rm -rf html smilei_warn.log

