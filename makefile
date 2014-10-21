VERSION:="$(shell git describe 2>/dev/null || echo '??')"

.PHONY: doc src

default: release

release:
	make -C src

openmpgnu:
	make -C src openmp=gnu

openmpintel:
	make -C src openmp=intel

debug:
	make -C src config=debug

scalasca:
	make -C src config=scalasca

clean:
	make -C src clean
	rm -rf doc/latex doc/html
	rm -rf build

doc: htmldoc pdfdoc

htmldoc:
#	convert -resize 100x100 doc/smileiLogo/smileiLogo.png doc/logo.png
	cd doc; (cat smilei.dox; echo "PROJECT_NUMBER=${VERSION}") | doxygen -

pdfdoc:
	cd doc/latex; pdflatex refman.tex; bibtex refman; pdflatex refman.tex
