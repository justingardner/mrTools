MLAB = /Applications/MATLAB7.7.app####Example of a Matlab path on a Mac
# MLAB = /usr/local/MATLAB/R2013b####Example of a Matlab path on a linux machine
MXSFX = mexmaci64####MEX extension for a MAC Intel 64bit
# MXSFX = mexa64####MEX extension for a Linux 64bit platform
MEX = ${MLAB}/bin/mex
ARCH = x86_64
MXARCH = -maci64
# MXARCH = -glnxa64

MFLAGS = -largeArrayDims -cxx ${MXARCH}
INC = -I ${MLAB}/extern/include
LIB = -L ${MLAB}/extern/lib

CC = cc -Wall -pedantic -no-cpp-precomp -arch ${ARCH} # some of these options might be obsolete
# CC = cc -Wall -pedantic -fPIC # This worked on linux 64 bit
C_OPTIMIZE_SWITCH = -O2    ## For GCC
CFLAGS = ${C_OPTIMIZE_SWITCH} ${INC} ${LIB}

all: corrDn.${MXSFX} upConv.${MXSFX} 

clean:
	/bin/rm *.o

corrDn.${MXSFX}: corrDn.o wrap.o convolve.o edges.o
	${MEX} ${MFLAGS} corrDn.o wrap.o convolve.o edges.o

upConv.${MXSFX}: upConv.o wrap.o convolve.o edges.o
	${MEX} ${MFLAGS} upConv.o wrap.o convolve.o edges.o

convolve.o wrap.o edges.o: convolve.h 

%.o : %.c
	${CC} -c ${CFLAGS} $<		

