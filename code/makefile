# programming environment
COMPILER     := mpic++
INCLUDE      := -I/home/vilasini/Thesis/progs/LATfield2 -I/usr/local/include -I/usr/include/hdf5/serial # add the path to LATfield2 and other libraries (if necessary)
LIB          := -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lm -lgsl -lgslcblas -lfftw3 #-lm  -lgsl -lgslcblas


# target and source
EXEC         := global
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   :=  -DHDF5 -DH5_HAVE_PIXIE -DFFT3D

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2)  $(INCLUDE) $(LIB)



clean:
	-rm -f $(EXEC) lccat lcmap
