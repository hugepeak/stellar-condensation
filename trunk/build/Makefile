ifndef BUILD_DEF

GC=gcc -g
FC=gfortran -g

#===============================================================================
# Versions and urls.
#===============================================================================

LIBNUCNET_VERSION = 0.13
WNMATRIX_VERSION = 0.13
WNSPARSE_VERSION = 0.4

SOURCEFORGE_URL = http://sourceforge.net/projects/
WEBNUCLEO_URL = http://www.webnucleo.org/home/modules
SKIT_URL = http://www-users.cs.umn.edu/~saad/software/SPARSKIT/

LIBNUCNET_URL = $(SOURCEFORGE_URL)/libnucnet/files/libnucnet
WNMATRIX_URL = $(SOURCEFORGE_URL)/wnmatrix/files/wn_matrix
WNSPARSE_URL = $(WEBNUCLEO_URL)/wn_sparse_solve/downloads

#===============================================================================
# VENDORDIR is the directory for storing webnucleo codes.
#===============================================================================

ifndef VENDORDIR
VENDORDIR = ../vendor
VENDOR_DIR := $(shell mkdir -p $(VENDORDIR))
endif

#===============================================================================
# Directories.
#===============================================================================

LIBNUCNETDIR = $(VENDORDIR)/libnucnet/$(LIBNUCNET_VERSION)/src
WNMATRIXDIR = $(VENDORDIR)/wn_matrix/$(WNMATRIX_VERSION)/src
WNSPARSEDIR = $(VENDORDIR)/wn_sparse_solve/$(WNSPARSE_VERSION)/src
SKITDIR = $(VENDORDIR)/SPARSKIT2

VALGRIND= yes
PROFILE= no

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

#===============================================================================
# Sparsekit.
#===============================================================================

include ../build/Makefile.sparsekit

#===============================================================================
# Include directories and libraries.
#===============================================================================

VPATH= $(WNMATRIXDIR):$(LIBNUCNETDIR):$(WNSPARSEDIR):$(SKITDIR):$(SKITDIR)/ITSOL:$(SKITDIR)/MATGEN:$(SKITDIR)/UNSUPP
CINCLUDE= `xml2-config --cflags` `gsl-config --cflags` -I$(LIBNUCNETDIR) -I$(WNMATRIXDIR) -I$(WNSPARSEDIR) -I$(SKITDIR)
CLIBS= `xml2-config --libs` `gsl-config --cflags --libs`
FLIBS= -L$(SKITDIR) -lskit

#===============================================================================
# Compiler flags.
#===============================================================================

CFLAGS1= -ansi -Werror -Wall -W -pedantic\
         -Wpointer-arith -Wcast-qual -Wcast-align \
         -Wwrite-strings -Wconversion \
         -fshort-enums -fno-common -g \

ifeq ($(VALGRIND), yes)
	CFLAGS2= $(CFLAGS1) -O0
else
	CFLAGS2= $(CFLAGS1) -O2
endif

ifeq ($(PROFILE), yes)
	CFLAGS= $(CFLAGS2) -pg
else
	CFLAGS= $(CFLAGS2)
endif

CC=$(GC) $(CFLAGS) $(CINCLUDE)

#===============================================================================
# OBJDIR is the temporary directory for codes compilation, this is where
# object files are created. 
#===============================================================================

ifndef OBJDIR
OBJDIR = ../obj
OBJ_DIR := $(shell mkdir -p $(OBJDIR))
endif

#===============================================================================
# BINDIR is the temporary directory for code executables.
#===============================================================================

ifndef BINDIR
BINDIR = .
BIN_DIR := $(shell mkdir -p $(BINDIR))
endif

#===============================================================================
# Build dependencies.
#===============================================================================

$(VENDORDIR)/SPARSKIT2.tar.gz:
	wget $(SKIT_URL)/SPARSKIT2.tar.gz -O $@
	tar xz -C $(VENDORDIR) -f $@

lib: $(VENDORDIR)/SPARSKIT2.tar.gz libskit

$(VENDORDIR)/wn_matrix.tar.gz:
	wget $(WNMATRIX_URL)/wn_matrix_$(WNMATRIX_VERSION).tar.gz -O $@

$(OBJDIR)/WnMatrix.o: $(VENDORDIR)/wn_matrix.tar.gz 
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(WNMATRIXDIR)/WnMatrix.c -o $@

$(VENDORDIR)/wn_sparse_solve.tar.gz:
	wget $(WNSPARSE_URL)/wn_sparse_solve_$(WNSPARSE_VERSION).tar.gz -O $@

$(OBJDIR)/WnSparseSolve.o: $(VENDORDIR)/wn_sparse_solve.tar.gz
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(WNSPARSEDIR)/WnSparseSolve.c -o $@

$(VENDORDIR)/libnucnet.tar.gz:
	wget $(LIBNUCNET_URL)/libnucnet_$(LIBNUCNET_VERSION).tar.gz -O $@

$(OBJDIR)/Libnucnet__Nuc.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet__Nuc.c -o $@

$(OBJDIR)/Libnucnet__Reac.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet__Reac.c -o $@

$(OBJDIR)/Libnucnet.o: $(VENDORDIR)/libnucnet.tar.gz $(OBJDIR)/WnMatrix.o
	tar xz -C $(VENDORDIR) -f $<
	$(CC) -c $(LIBNUCNETDIR)/Libnucnet.c -o $@

#--------------------------------------------------------------------------

WN_OBJ =$(OBJDIR)/WnMatrix.o		\
	$(OBJDIR)/Libnucnet__Nuc.o	\
	$(OBJDIR)/Libnucnet__Reac.o	\
	$(OBJDIR)/Libnucnet.o		\
	$(OBJDIR)/WnSparseSolve.o	\

#
# non-library and unsupported objects
#
SP_OBJ = $(OBJDIR)/exppro.o	\
	 $(OBJDIR)/phipro.o	\
	 $(OBJDIR)/blas1.o	\
	 $(OBJDIR)/itaux.o	\

$(SP_OBJ): $(OBJDIR)/%.o : $(WNSPARSEDIR)/%.f
	$(FC) -c -o $@ $<

#===============================================================================
# Get data.
#===============================================================================

ifdef NUCNET_TARGET

DATA_URL=http://www.webnucleo.org/experimental/mbradle/download/network_data
DATA=astr820.tar
DATA_DIR=$(NUCNET_TARGET)/data_pub

data:
	mkdir -p $(DATA_DIR)
	wget $(DATA_URL)/$(DATA).gz
	mv ${DATA}.gz $(DATA_DIR)
	tar xz -C ${DATA_DIR} -f ${DATA_DIR}/${DATA}.gz
	rm ${DATA_DIR}/${DATA}.gz

endif

#===============================================================================
# Get data.
#===============================================================================

DATA_URL=http://www.webnucleo.org/experimental/mbradle/download/special
DATA=stellar_cond_data.tar
DATA_DIR=../data_pub

data:
	mkdir -p $(DATA_DIR)
	wget $(DATA_URL)/$(DATA).gz
	mv ${DATA}.gz $(DATA_DIR)
	tar xz -C ${DATA_DIR} -f ${DATA_DIR}/${DATA}.gz
	rm ${DATA_DIR}/${DATA}.gz

#===============================================================================
# Clean up. 
#===============================================================================

.PHONY: clean clean_all clean_dist

clean:
	rm -fr $(OBJDIR)

clean_dist: clean
	rm -fr $(VENDORDIR)

clean_data:
	rm -fr $(DATA_DIR)

#===============================================================================
# End include.
#===============================================================================

BUILD_DEF = yes

endif

