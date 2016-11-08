#***********************************************************************
# FILENAME  : simple/Makefile
# DATE      : 10/03
# DESIGNER  : CCDAS group
# PURPOSE   : generate code
#		 and build CCDAS objects of  SIMPLE test code
#***********************************************************************

# include Make.inc (specific makefile with macros)
#include $(TOP)/Make.inc
# set top level directory
ifeq ($(TOP),)
   TOP = $(PWD)
endif

#include $(TOP)/Makeicl_local

#---------------------------------------------------------
# simple test function, code, adcode and objects
#	all starts with letter T
#---------------------------------------------------------

SIMPLE = $(TOP)/simple# simple test program
EXT    = .f90#		extension of code files

# source code (and generated code)
TSRC = simple.f90
TADSRC  = simple_ad.f90
TTLSRC  = simple_tl.f90
THESSRC	= taf_hes.f90 simple_ad_hes.f90 simple_func_hes.f90
TPTLSRC = taf_ptl.f90 simple_ptl.f90 simple_pfunc_ptl.f90 
TPADSRC = taf_pad.f90 simple_pad.f90 simple_pfunc_pad.f90 
TSPADSRC = ttaf_spad.f90 simple_spad.f90 simple_pfunc_spad.f90 
TSPTLSRC = ttaf_sptl.f90 simple_sptl.f90 simple_pfunc_sptl.f90 
TITLSRC = taf_itl.f90 simple_ad_itl.f90 simple_ifunc_itl.f90 

# object code
TOBJ    = $(TSRC:.f90=.o)
TTLOBJ  = $(TTLSRC:.f90=.o)
TADOBJ  = $(TADSRC:.f90=.o)
THESOBJ = $(THESSRC:.f90=.o)
TPTLOBJ = $(TPTLSRC:.f90=.o)
TPADOBJ = $(TPADSRC:.f90=.o)
TSPTLOBJ = $(TSPTLSRC:.f90=.o)
TSPADOBJ = $(TSPADSRC:.f90=.o)
TITLOBJ = $(TITLSRC:.f90=.o)

# test code
TEXE 	= $(TOP)/tcost $(TOP)/tcost2 $(TOP)/ttsttlm \
	$(TOP)/ttstadm $(TOP)/ttsthesscol  $(TOP)/ttstptlm $(TOP)/ttstpadm \
	$(TOP)/tnextx $(TOP)/tfgrad $(TOP)/trev $(TOP)/ttstspadm

#---------------------------
# pattern rules for drivers
#---------------------------
$(DRIVERS)/%.o: $(DRIVERS)/%.f90
	$(F90) $(F90FLAGS) $(DBL) -c -o $@ $<

#---------------------------------------------------------
# suffix rules
#---------------------------------------------------------
.SUFFIXES : .o .f90 .mod .f .d .F90

.f.o :
	$(F90) $(F90FLAGS) $(DBL) -c $<

.f90.o :
	$(F90) $(F90FLAGS) $(DBL) -c $<

.f90.mod :
	$(F90) $(F90FLAGS) $(DBL) -c $<

.f90.d :
	$(F90) $(F90FLAGS) $(DBL) -c $<

.F90.f90 :
	$(CPP) $(CPPFLAGS) $< > $@

# disable default suffix rule for .mod files
.mod.o:	
	true

#---------------------------------------------------------
# generate derivative code:
#---------------------------------------------------------
adm $(TADSRC): $(TSRC)
	$(TAF) -l - $(TTAFFLAGS) $(MADFLAGS) -reverse $(INCDIRS) $(TSRC)

tlm $(TTLSRC): $(TSRC)
	$(TAF) -l - $(TTAFFLAGS) $(MADFLAGS) -forward $(INCDIRS) $(TSRC)

hess $(THESSRC): $(TADSRC) simple_func.f90
	$(TAF) -l - $(TTAFFLAGS) $(VADFLAGS) -forward -ftlmark _hes -ldg -vecfor 2	\
	   $(TADSRC) simple_func.f90

padm $(TPADSRC): $(TSRC) simple_pfunc.f90
	$(TAF) -l - $(TTAFFLAGS) $(VADFLAGS) -reverse -admark _pad $(INCDIRS)	\
	-vecrev 2 -ldg $(INCDIRS) $(TSRC) simple_pfunc.f90

ptlm $(TPTLSRC): $(TSRC) simple_pfunc.f90
	$(TAF) -l - $(TTAFFLAGS) $(VADFLAGS) -forward -ftlmark _ptl $(INCDIRS)	\
	-vecfor 2 -ldg $(INCDIRS) $(TSRC) simple_pfunc.f90

spadm $(TSPADSRC): $(TSRC) simple_pfunc.f90
	$(TAF) -l - $(TTAFFLAGS) $(TAF_SPAD) $(VADFLAGS) -reverse $(INCDIRS)	\
	$(TSRC) simple_pfunc.f90
	- mv taf_spad.f90 ttaf_spad.f90

sptlm $(TSPTLSRC): $(TSRC) simple_pfunc.f90
	$(TAF) -l - $(TTAFFLAGS) $(TAF_SPTL) $(VADFLAGS) -forward $(INCDIRS)	\
	$(TSRC) simple_pfunc.f90
	- mv taf_sptl.f90 ttaf_sptl.f90

itlm $(TITLSRC): $(TADSRC) simple_ifunc.f90
	$(TAF) -l - -v2 -toplevel func -input xbgr -output y -forward -ftlmark _itl $(INCDIRS)	\
	-vecfor 1 -ldg $(INCDIRS) $(TADSRC) simple_ifunc.f90

#---------------------------------------------------------
# targets for executables
#---------------------------------------------------------
#HEW-D../tcost: $(DRIVERS)/prgcost.o $(DRIVERS)/storef.o $(TOBJ)
#HEW-D	$(F90) $(DBL) $(F90FLAGS) $(LIBDIRS) -o $@	\
#HEW-D	$(DRIVERS)/prgcost.o $(DRIVERS)/storef.o $(TOBJ) $(LIBS)

../cost: $(DRIVERS)/prgcost.o $(DRIVERS)/storef.o $(TOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgcost.o $(DRIVERS)/storef.o $(TOBJ) $(LIBS)

# calculate costfunction twice
../cost2: $(DRIVERS)/prgcost2.o $(DRIVERS)/storef.o $(TOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgcost2.o $(DRIVERS)/storef.o $(TOBJ) $(LIBS)

# test tangent linear * vector against finite differences
../tsttlm: $(DRIVERS)/prgtsttlm.o $(DRIVERS)/storef.o $(TTLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtsttlm.o $(DRIVERS)/storef.o $(TTLOBJ) $(LIBS)

# test gradient from adjoint against finite differences
../tstadm: $(DRIVERS)/prgtstadm.o $(DRIVERS)/storef.o $(TADOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstadm.o $(DRIVERS)/storef.o $(TADOBJ) $(LIBS)

# compute columns of Jacobian in forward mode and check against finite differences
../tsthesscol: $(DRIVERS)/prgtstjaccol.o $(DRIVERS)/storef.o $(THESOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccol.o $(DRIVERS)/storef.o $(THESOBJ) $(LIBS)

# determine next control vector x in current optimisation iteration step
../nextx: $(DRIVERS)/prgopti.o $(DRIVERS)/storef.o $(TADOBJ) # $(LSOPT)/liblsopt.a
	$(F90) $(OPTLIBDIRS)	\
		-o $@ $(DRIVERS)/prgopti.o $(DRIVERS)/storef.o $(TADOBJ) $(OPTLIBS)

# compute the gradient of the costfunction
../fgrad: $(DRIVERS)/prggrad.o $(DRIVERS)/storef.o $(TADOBJ) 
	$(F90) $(OPTLIBDIRS)	\
		-o $@ $(DRIVERS)/prggrad.o $(DRIVERS)/storef.o $(TADOBJ) $(OPTLIBS)

# columns of Jacobian in forward mode, check against finite differences
../tstptlm: $(DRIVERS)/prgtstjaccolptl.o $(DRIVERS)/storef.o $(TPTLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccolptl.o $(DRIVERS)/storef.o $(TPTLOBJ) $(LIBS) 

# columns of Jacobian in reverse mode and check against finite differences
../tstpadm: $(DRIVERS)/prgtstjaccolpad.o $(DRIVERS)/storef.o $(TPADOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccolpad.o $(DRIVERS)/storef.o $(TPADOBJ) $(LIBS) 

# 
../rev: $(DRIVERS)/prgrev.o $(DRIVERS)/storef.o $(TADOBJ)
	$(F90) $(LIBDIRS)	\
		-o $@ $(DRIVERS)/prgrev.o $(DRIVERS)/storef.o $(TADOBJ) $(LIBS)

#
../hesscol: $(DRIVERS)/prgjaccol.o $(DRIVERS)/storef.o $(TADOBJ) $(THESOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgjaccol.o $(DRIVERS)/storef.o $(TADOBJ) $(THESOBJ) $(LIBS)

# 
../tstspadm: $(DRIVERS)/prgasd.o $(DRIVERS)/storef.o $(TSPTLOBJ) $(TSPADOBJ) $(TPADOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgasd.o $(DRIVERS)/storef.o $(TSPTLOBJ) $(TSPADOBJ) $(TPADOBJ) $(LIBS) 

# columns of Jacobian in forward mode, check against finite differences
../tstitlm: $(DRIVERS)/prgtstjaccolitl.o $(DRIVERS)/storef.o $(TITLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccolitl.o $(DRIVERS)/storef.o $(TITLOBJ) $(LIBS) 

#---------------------------------------------------------
# run an optimisation
#---------------------------------------------------------

lsopt: ../nextx ../fgrad $(TOP)/trunopti
	cd $(TOP) ; $(TOP)/trunopti

../dfp: $(DRIVERS)/prgdfpmin_no_netcdf.o $(DRIVERS)/storef.o $(TADOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgdfpmin_no_netcdf.o $(DRIVERS)/storef.o $(TADOBJ) $(LIBS) ../libnumrec.a

../line: $(DRIVERS)/line.o $(TOBJ) 
	$(F90) $(LIBDIRS)	\
		-o $@ $(DRIVERS)/line.o $(TOBJ) $(LIBS)

../dline: $(DRIVERS)/dline.o $(TTLOBJ) 
	$(F90) $(LIBDIRS)	\
		-o $@ $(DRIVERS)/dline.o $(TTLOBJ) $(LIBS)

#---------------------------------------------------------
# for simple not implemented targets (e.g. prognostic targets) :
#---------------------------------------------------------

../pfunc2 ../prev ../pfwd ../printp ../parmap :



#---------------------------------------------------------
# rules for removing
#---------------------------------------------------------
# remove .o, .mod rsp .d files
clearo:
	@$(RM) *.o *.d *.mod

# remove intermediate files, warmstart still possible
clean:	clearo
	@$(RM) taf_command taf_output

# remove generated files, only coldstart possible
scratch: clean
	@$(RM) *.out
	@$(RM) -r adcode  hesscode  padm  ptlcode spadcode  tlcode
	@$(RM) *_ad.f90 *tl.f90 *_hes.f90 *_ptl.f90 *_pad.f90 *_spad.f90 *_sptl.f90 *_itl.f90 

rmall: scratch
	@$(RM) *~ \#*\# *.out 
	@$(RM) -r tstoutput* 
	@$(RM) OPWARM* lsopt.tmp

# dependencies test environment
simple_sptl.o : ttaf_sptl.o
simple_spad.o : ttaf_spad.o










