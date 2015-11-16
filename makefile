VERSION:="$(shell git describe 2>/dev/null || echo '??')"

.PHONY: doc src

default: release

release:
	make -C src

openmp:
	make -C src config=openmp

openmpintelmpi:
	make -C src openmp=intelmpi

debug:
	make -C src config=debug

scalasca:
	make -C src config=scalasca

clean:
	make -C src clean
	make -C doc clean

doc:
	make -C doc all

tar:
	git archive -o smilei-$(VERSION).tgz --prefix smilei-$(VERSION)/ HEAD
	
