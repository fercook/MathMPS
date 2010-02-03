
include ARmake.inc

#
# issue make all
#

# To link, this works from the command line
#
#   ifort -O3 -m32 -arch i386 FortranArpackDriver.o arpackformps.o mathematicatemplatetm.o -L/Users/fercook/ARPACK/ -larpack_MAC -I/Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/CompilerAdditions -L/Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/CompilerAdditions -lMLi3 -o TEST -nofor-main  -lstdc++
#
#

ARCH:=$(shell uname)

ifeq ($(ARCH),Darwin)
VERSION=6.0
MLINKDIR = /Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/CompilerAdditions
EXTRA_CFLAGS = -m32 #-arch i386  #-Wno-long-double 
INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}
MPREP = ${CADDSDIR}/mprep
ARPACKDIR=$(HOME)/ARPACK
ARPACKlibb=arpack_MAC
RM=rm
MLINKLIB=MLi3
CC=g++
LINKER=ifort
LINKFLAGS=-nofor-main -lstdc++
endif

ifeq ($(ARCH),Linux)
VERSION=6.0
MLINKDIR =  #/usr/local/Wolfram/Mathematica/${VERSION}/SystemFiles/Links/MathLink/DeveloperKit
FFLAGS = -O3 -m64
SYS = Linux-x86-64
CADDSDIR =.  #${MLINKDIR}/${SYS}/CompilerAdditions
EXTRA_CFLAGS = -m64 #-arch i386  #-Wno-long-double
INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}
MPREP = ./mprep #${CADDSDIR}/mprep
ARPACKDIR=$(HOME)/ARPACK
ARPACKlibb=arpack_X86
RM=rm
MLINKLIB=ML64i3
CC=g++
LINKER=ifort
LINKFLAGS=-nofor-main -lm -lpthread -lrt -lstdc++
endif

all: arpackformps
obj: object
mlink: mathlink
exec: executable
test: secondtest
debug: secondtestdebug


#--------  HERE START THE USEFUL BITS

object: ArpackForMPS.f90
	 $(FC) $(FFLAGS) -c ArpackForMPS.f90

mathlink: mathematicatemplate.tm
	${MPREP} $? -o $@

arpackformps: mathematicatemplatetm.o arpackformps.o FortranArpackDriver.o
	${LINKER} ${EXTRA_CFLAGS} -I${INCDIR} mathematicatemplatetm.o arpackformps.o FortranArpackDriver.o ${LINKFLAGS} -L${LIBDIR}  -L${ARPACKDIR} -l${ARPACKlibb} -l${MLINKLIB} -o $@
	cp arpackformps $(SYS)/.

mathematicatemplatetm.c: mathematicatemplate.tm 
	${MPREP} $? -o $@

FortranArpackDriver.o: FortranArpackDriver.f90
	 $(FC) $(FFLAGS) -c FortranArpackDriver.f90

arpackformps.o: arpackformps.cpp
	${CC} -c ${EXTRA_CFLAGS} -I${INCDIR} arpackformps.cpp

mathematicatemplatetm.o : mathematicatemplatetm.c
	${CC} -c ${EXTRA_CFLAGS} -I${INCDIR} mathematicatemplatetm.c


clean :
	@ ${RM} -rf *.o *tm.c $(BINARIES)

#-----------------------------------------------------------------------
#
firsttest: complexdriver.o
	$(FC) $(FFLAGS) complexdriver.f90 $(ALIBS) -o firsttest
#
secondtest: spinchaindriverTEST.o drivermodule.f90
	$(FC) $(FFLAGS) drivermodule.f90 spinchaindriverTEST.f90 $(ALIBS) -o secondtest

secondtestdebug:
	$(FC) $(FFLAGS) drivermodule.f90 spinchaindriverTEST.f90 $(ALIBS) -g -o secondtest

