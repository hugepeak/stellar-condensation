BUILD_DIR = ../build

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

VPATH = $(BUILD_DIR)

#===============================================================================
# Compile main codes.
#===============================================================================

$(OBJDIR)/testF.o: testF.f 
	$(FC) -c -o $(OBJDIR)/testF.o testF.f

$(OBJDIR)/testC.o: testC.c 
	$(CC) -c -o $(OBJDIR)/testC.o testC.c

#--------------------------------------------------------------------------

CON_OBJ = $(OBJDIR)/con1115b.o         		\
          $(OBJDIR)/con1123c.o       		\
          $(OBJDIR)/con1115d.o      		\
          $(OBJDIR)/con0330e.o      		\
          $(OBJDIR)/con1115f.o			\

$(CON_OBJ): $(OBJDIR)/%.o : %.f
	$(FC) -c -o $@ $<

OBJ =   $(WN_OBJ)				\
	$(CON_OBJ)				\
	$(SP_OBJ)				\
        $(OBJDIR)/testC.o			\

#===============================================================================
# Compile codes depending only on Libnucnet__Nuc.
#===============================================================================

testF : lib $(OBJDIR)/testF.o $(OBJ)
	$(FC) $(OBJDIR)/testF.o $(OBJ) $(CLIBS) $(FLIBS) -o testF

cwin: $(OBJDIR)/con1123a.o $(COBJ)
	$(FC) $(OBJDIR)/con1123a.o $(COBJ) $(CLIBS) -o cwin

#===============================================================================
# Clean up.
#===============================================================================

cleanall: clean
	rm testF


