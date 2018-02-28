INSTRUCTIONS FOR COMPILING SMILEI WITH PICSAR 

1)clone Picsar
git clone git@bitbucket.org:berkeleylab/picsar.git

2)Move to master branch of picsar 

3)set MODE=library into the Makefile

  set COMP=intel
  set FARGS= -g -check bound -O3 -fopenmp -module ./
  set FFTW3_LIB=$($FFTW_LIB_DIR)
  set FFTW3_INCLUDE=$(FFTW_INC_DIR)

4) Launch:

make lib
Now the libpxr.so and libpxr.a are in the lib/ file

5)Link Smilei with libpxr.so
Set PICSAR=TRUE in SMILEI makefile
LIBPXR = ~/picsar/lib
LDFLAGS += -L$(LIBPXR) -lpxr

6) Link fftw to smilei
LDFLAGS += -I$(FFTW3_INCLUDE)
LDFLAGS += -L$(FFTW3_LIB) -lfftw3_mpi
LDFLAGS += -L$(FFTW3_LIB) -lfftw3_threads
LDFLAGS += -L$(FFTW3_LIB) -lfftw3

7) Add lgfortran flag to compliling options 
LDFLAGS += -lgfortran


