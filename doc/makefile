VERSION:="$(shell git describe 2>/dev/null || echo '??')-$(shell git rev-parse --abbrev-ref HEAD || echo '??')"

TEXFILES=

PDFS:=$(TEXFILES:.tex=.pdf)

.PHONY: html

default: all

all: $(PDFS) html

%.pdf: %.tex
	latexmk -cd -pdf $< && latexmk -cd -c $<

html:
#	convert -resize 100x100 doc/smileiLogo/smileiLogo.png doc/logo.png
	@echo "-------------";
	@if type "doxygen" >/dev/null 2>&1; then\
		cd Doxygen; (echo "PROJECT_NUMBER=${VERSION}"; cat smilei.dox) | doxygen - ;\
	else \
		echo "Cannot build doxygen doc because doxygen is not installed";\
	fi	
	@echo "-------------";
#   pdf generation is buggy... dropped
#	latexmk -cd -pdf latex/refman.tex && latexmk -cd -c latex/refman.tex
	@if type "sphinx-build" >/dev/null 2>&1; then\
		cd Sphinx; make html;\
	else \
		echo "Cannot build Sphinx doc because Sphinx is not installed";\
	fi
	@echo "-------------";
	
	
clean:
	for fname in $(TEXFILES); do latexmk -cd -C $$fname; done;
	rm -rf Doxygen/html Doxygen/smilei_warn.log
	rm -rf Sphinx/_build
