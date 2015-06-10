VERSION:="$(shell git describe 2>/dev/null || echo '??')"

TEXFILES=Collisions/Collisions.tex Diagnostics/Diagnostics.tex Profiles/Profiles.tex latex/refman.tex

PDFS:=$(TEXFILES:.tex=.pdf)



.PHONY: doc


default: doc

clean:
	for fname in $(TEXFILES); do cd `dirname $$fname` && echo `pwd` && latexmk -C && cd ..; done
	rm -rf latex html

doc: html $(PDFS)

PDFLATEX_FLAGS=-interaction=batchmode -halt-on-error -file-line-error

%.pdf: %.tex
	@ echo "Creating doc for $< : $@"
	cd $(shell  dirname $<) && latexmk -pdf $(shell basename $<)


html: smilei.dox
#	convert -resize 100x100 doc/smileiLogo/smileiLogo.png doc/logo.png
	cd doc; (cat smilei.dox; echo "PROJECT_NUMBER=${VERSION}") | doxygen -

