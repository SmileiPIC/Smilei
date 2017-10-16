INSTRUCTIONS FOR COMPILING SMILEI WITH PICSAR 

1)clone Picsar
git clone https://kallalahaythem@bitbucket.org/berkeleylab/picsar.git

2)Move to gpstd_group branch
git fetch && git checkout gpstd_group

3)set prod=library into the Makefile

4) Launch:

make lib
Now the libpxr.so and libpxr.a are in the lib/ file

5)Link Smilei with libpxr.so
LIBPXR = ~/picsar/lib
LDFLAGS += -L$(LIBPXR) -lpxr

6) Link fftw to smilei
LDFLAGS += -I$(FFTW3_INCLUDE)
LDFLAGS += -L$(FFTW3_LIB) -lfftw3_mpi
LDFLAGS += -L$(FFTW3_LIB) -lfftw3_threads
LDFLAGS += -L$(FFTW3_LIB) -lfftw3

7) Add lgfortran flag to compliling options 
LDFLAGS += -lgfortran


