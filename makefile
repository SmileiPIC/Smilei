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
	make -C doc clean

doc:
	make -C doc all

