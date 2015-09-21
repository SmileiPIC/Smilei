VERSION:="$(shell git describe 2>/dev/null || echo '??')"

.PHONY: doc src

default: 
	@echo "Making in src"
	make -C src $(MAKECMDGOALS) $(MAKEFLAGS)
	@echo "Going out of src"

clean:
	make -C src clean
	make -C doc clean
	rm -rf smilei-$(VERSION).tgz

doc:
	make -C doc all

tar:
	git archive -o smilei-$(VERSION).tgz --prefix smilei-$(VERSION)/ HEAD
	
help: 
	@echo 'Usage: make [clean] [config=OPTIONS]'
	@echo '		OPTIONS is a string composed of one or more of:'
	@echo '	        debug    : to compile in debug mode (code runs really slow)'
	@echo '	        openmp   : to compile with openmp enabled'
	@echo '         scalasca : to compile using scalasca'
	@echo ' e.g.  make config="debug openmp"'
	@echo ''
	@echo 'Environment variables :'
	@echo '      SMILEICXX (mpi c++ compiler)'
	@echo '      HDF5_ROOT_DIR (HDF5 dir)'
	@echo ''
	@echo '       make doc : builds the documentation'
	@echo '       make tar : creates an archive of the sources'


%::
	@echo "Making in src"
	make -C src $(MAKECMDGOALS) $(MAKEFLAGS)
	@echo "Going out of src"
