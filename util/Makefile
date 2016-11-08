#*******************************************************************
# FILENAME  : util/Makefile
# DATE      : 10/03
# DESIGNER  : CCDAS group
# PURPOSE   : generate code and build CCDAS objects for the UTILity programs
#*******************************************************************

# This is the Makefile of the UTIL source code. 

# include Make.inc (specific Makefile with macros)
#include $(TOP)/Make.inc

# options to include local settings and overwrite the above (HEW: better in $(TOP)/Make.inc)
#include $(TOP)/Makeicl_local
#include $(TOP)/Makeicl_eiche_lf95
# include $(TOP)/Makeicl_eiche_ifc

#VPATH=$(BETHY)

#---------------------------------------------------------
# suffix rules
#---------------------------------------------------------
.SUFFIXES : .o .f90 .mod .f .d 

.f.o :
	$(F90) $(F90FLAGS) $(DBL) -c $(@D)/$(<F) -o $(OBJ)/$(@F)

.f90.o :
	$(F90) $(F90FLAGS) $(DBL) -c $<

.f90.mod :
	$(F90) $(F90FLAGS) $(DBL) -c $<
#	mv $@ $(OBJDIR)

.f90.d :
	$(F90) $(F90FLAGS) $(DBL) -c $<

#---------------------------------------------------------
# rules for utility programms
#---------------------------------------------------------
pfwdunc.f90: pfwdunc.g90 
	$(CPP) $(CPPFLAGS) $< >$@

../pfwdunc: pfwdunc.o
	$(F90) $(DBL)  -o $@ pfwdunc.o  

../func: func.o
	$(F90) $(DBL) $(LIBDIRS) -L$(LAPACK) -L$(BLAS1) -o $@ \
	func.o $(LIBS) -l$(LAPACKLIB) -l$(BLAS1LIB) 

../dpdb: dpdb.o
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	dpdb.o $(LIBS) 

../printp: printp.o
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	printp.o $(LIBS) 

#---------------------------------------------------------
# remove .o and .mod files
#---------------------------------------------------------
# remove .o, .mod rsp .d files
clearo:
	@$(RM) pfwdunc.f90
	@$(RM) *.o *.d *.mod 

rmall:
	@$(RM) *~ \#*\# *.out *.o *.mod *.d 


