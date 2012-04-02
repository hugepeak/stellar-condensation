#///////////////////////////////////////////////////////////////////////////////
# The following lines must be edited according to where you have
# placed your libnucnet and wn_matrix modules.  You may also
# edit to choose a different compiler (e.g., g++), to use Valgrind or
# not, or to profile:
#///////////////////////////////////////////////////////////////////////////////

GC=gcc -g
FC=gfortran -g

MODULES = /Users/bradleymeyer/modules

SPARSKITDIR = $(MODULES)/SPARSKIT2/

LIBNUCNETDIR = $(MODULES)/libnucnet/0.13/src/
MATRIXSRCDIR = $(MODULES)/wn_matrix/0.13/src/
SOLVEDIR = $(MODULES)/wn_sparse_solve/0.4/src/
VALGRIND= yes
PROFILE= no

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

CINCLUDE= `xml2-config --cflags` `gsl-config --cflags` -I$(LIBNUCNETDIR) -I$(MATRIXSRCDIR) -I$(SOLVEDIR)
CLIBS= `xml2-config --libs` `gsl-config --libs` -lm -L$(SPARSKITDIR) -lskit

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS1= -ansi -Werror -Wall -W \
         -Wconversion -Wshadow \
         -Wpointer-arith -Wcast-qual -Wcast-align \
         -Wwrite-strings \
         -fshort-enums -fno-common -Dinline= -g

ifeq ($(GC), gcc) 
	CFLAGS2= $(CFLAGS1) -Wmissing-prototypes -Wstrict-prototypes -Wnested-externs
else
	CFLAGS2= $(CFLAGS1)
endif

ifeq ($(VALGRIND), yes)
	CFLAGS3= $(CFLAGS2) -O0
else
	CFLAGS3= $(CFLAGS2) -O2
endif

ifeq ($(PROFILE), yes)
	CFLAGS= $(CFLAGS3) -pg
else
	CFLAGS= $(CFLAGS3)
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)


#===============================================================================
# TMPDIR is the temporary directory for codes compilation, this is where
# object files are created. 
#===============================================================================

ifndef TMPDIR
TMPDIR = ./tmp/
TMP_DIR := $(shell mkdir -p $(TMPDIR))
endif

#===============================================================================
# Compile matrix codes.
#===============================================================================

$(TMPDIR)WnMatrix.o: $(MATRIXSRCDIR)WnMatrix.c
	$(CC) -c -o $(TMPDIR)WnMatrix.o $(MATRIXSRCDIR)WnMatrix.c

$(TMPDIR)WnSparseSolve.o: $(SOLVEDIR)WnSparseSolve.c
	$(CC) -c -o $(TMPDIR)WnSparseSolve.o $(SOLVEDIR)WnSparseSolve.c

#===============================================================================
# Compile unsupported sparse kit codes.
#===============================================================================
$(TMPDIR)itaux.o: $(SOLVEDIR)itaux.f
	$(FC) -c -o $(TMPDIR)itaux.o $(SOLVEDIR)itaux.f

$(TMPDIR)blas1.o: $(SOLVEDIR)blas1.f
	$(FC) -c -o $(TMPDIR)blas1.o $(SOLVEDIR)blas1.f

$(TMPDIR)exppro.o: $(SOLVEDIR)exppro.f
	$(FC) -c -o $(TMPDIR)exppro.o $(SOLVEDIR)exppro.f

$(TMPDIR)phipro.o: $(SOLVEDIR)phipro.f
	$(FC) -c -o $(TMPDIR)phipro.o $(SOLVEDIR)phipro.f

#===============================================================================
# Compile Libnucnet codes.
#===============================================================================
$(TMPDIR)Libnucnet__Nuc.o: $(LIBNUCNETDIR)Libnucnet__Nuc.c
	$(CC) -c -o $(TMPDIR)Libnucnet__Nuc.o $(LIBNUCNETDIR)Libnucnet__Nuc.c

$(TMPDIR)Libnucnet__Nse.o: $(LIBNUCNETDIR)Libnucnet__Nse.c
	$(CC) -c -o $(TMPDIR)Libnucnet__Nse.o $(LIBNUCNETDIR)Libnucnet__Nse.c

$(TMPDIR)Libnucnet__Reac.o: $(LIBNUCNETDIR)Libnucnet__Reac.c
	$(CC) -c -o $(TMPDIR)Libnucnet__Reac.o $(LIBNUCNETDIR)Libnucnet__Reac.c

$(TMPDIR)Libnucnet.o: $(LIBNUCNETDIR)Libnucnet.c
	$(CC) -c -o $(TMPDIR)Libnucnet.o $(LIBNUCNETDIR)Libnucnet.c

#===============================================================================
# Compile main codes.
#===============================================================================
$(TMPDIR)testC.o: testC.c
	$(CC) -c -o $(TMPDIR)testC.o testC.c

$(TMPDIR)testF.o: testF.f 
	$(FC) -c -o $(TMPDIR)testF.o testF.f

$(TMPDIR)con1123a.o: con1123a.f 
	$(FC) -c -o $(TMPDIR)con1123a.o con1123a.f

$(TMPDIR)con1115b.o: con1115b.f
	$(FC) -c -o $(TMPDIR)con1115b.o con1115b.f 

$(TMPDIR)con1123c.o: con1123c.f 
	$(FC) -c -o $(TMPDIR)con1123c.o con1123c.f

$(TMPDIR)con1115d.o: con1115d.f 
	$(FC) -c -o $(TMPDIR)con1115d.o con1115d.f

$(TMPDIR)con0330e.o: con0330e.f 
	$(FC) -c -o $(TMPDIR)con0330e.o con0330e.f

$(TMPDIR)con1115f.o: con1115f.f 
	$(FC) -c -o $(TMPDIR)con1115f.o con1115f.f

#--------------------------------------------------------------------------

COBJ = $(TMPDIR)con1115b.o              \
       $(TMPDIR)con1123c.o       \
       $(TMPDIR)con1115d.o      \
       $(TMPDIR)con0330e.o      \
       $(TMPDIR)con1115f.o

NOBJ =  $(TMPDIR)Libnucnet__Nuc.o   \
        $(TMPDIR)testC.o           

ROBJ =  $(TMPDIR)Libnucnet__Reac.o    

OBJ =   $(NOBJ)                     \
        $(ROBJ)                     \
        $(TMPDIR)WnMatrix.o         \
        $(TMPDIR)WnSparseSolve.o    \
        $(TMPDIR)exppro.o           \
        $(TMPDIR)phipro.o           \
        $(TMPDIR)blas1.o            \
        $(TMPDIR)itaux.o            \
        $(TMPDIR)Libnucnet.o        \
        $(TMPDIR)con1115b.o       \
        $(TMPDIR)con1123c.o       \
        $(TMPDIR)con1115d.o      \
        $(TMPDIR)con0330e.o      \
        $(TMPDIR)con1115f.o

#===============================================================================
# Compile codes depending only on Libnucnet__Nuc.
#===============================================================================

testF : $(TMPDIR)testF.o $(OBJ)
	$(FC) $(TMPDIR)testF.o $(OBJ) $(CLIBS) -o testF

cwin: $(TMPDIR)con1123a.o $(COBJ)
	$(FC) $(TMPDIR)con1123a.o $(COBJ) $(CLIBS) -o cwin


#===============================================================================
# Clean up.
#===============================================================================

clean: 
	rm -f $(TMPDIR)*.c
	rm -f $(TMPDIR)*.h
	rm -f $(TMPDIR)*.o
	rm -f *.exe

cleanall: clean
	rm testF


