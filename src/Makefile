BUILD_DIR = ../build

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

VPATH = $(BUILD_DIR)

#--------------------------------------------------------------------------

$(OBJDIR)/testC.o: testC.c 
	$(CC) -c -o $(OBJDIR)/testC.o testC.c

CON_OBJ = $(OBJDIR)/con1115b.o         		\
          $(OBJDIR)/con1123c.o       		\
          $(OBJDIR)/con1115d.o      		\
          $(OBJDIR)/con0330e.o      		\
          $(OBJDIR)/con1115f.o			\

$(CON_OBJ): $(OBJDIR)/%.o : %.f
	$(FC) -c -o $@ $<

COBJ =  $(WN_OBJ)				\
	$(SP_OBJ)				\
        $(OBJDIR)/testC.o			\

FOBJ =  $(CON_OBJ)				\



#===============================================================================
# Compile codes.
#===============================================================================

CEXEC = average_zones

FEXEC = testF

$(CEXEC) : lib $(COBJ)
	cp $(LIBNUCNETDIR)/../examples/average_zones.c .
	$(CC) -c -o $(OBJDIR)/$@.o $@.c
	$(FC) -o $@ $(CLIBS) $(FLIBS) $(OBJDIR)/$@.o $(COBJ)

$(FEXEC) : lib $(COBJ) $(FOBJ)
	$(FC) -c -o $(OBJDIR)/$@.o $@.f
	$(FC) -o $@ $(CLIBS) $(FLIBS) $(OBJDIR)/$@.o $(COBJ) $(FOBJ)

cwin: $(OBJDIR)/con1123a.o $(COBJ)
	$(FC) $(OBJDIR)/con1123a.o $(COBJ) $(CLIBS) -o cwin

#===============================================================================
# Clean up.
#===============================================================================

cleanall: clean
	rm -fr $(CEXEC) $(FEXEC)


