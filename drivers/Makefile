#*******************************************************************
# FILENAME  : Makefile
# DATE      : 10/03
# DESIGNER  : CCDAS group
# PURPOSE   : Makefile for the DRIVERS
#*******************************************************************

# This is the Makefile of the DRIVERS 

# include Make.inc (specific makefile with macros)
# include $(TOP)/Make.inc

#---------------------------------------------------------
# rules for removing
#---------------------------------------------------------
# remove .o, .mod rsp .d files
clearo:
	@$(RM) *.o *.d *.mod

# remove intermediate files, warmstart still possible
clean:	clearo
	@$(RM) *~ \#*\# 

scratch: clean

rmall: scratch
