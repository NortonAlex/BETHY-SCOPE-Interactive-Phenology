#*******************************************************************
# FILENAME  : $(TOP)/Makefile
# DATE      : 10/02
# DESIGNER  : CCDAS group
# PURPOSE   : build the targets of CCDAS
#*******************************************************************
# This is the CENTRAL Makefile of CCDAS. Changes for your specific shell environment
# or specific target architecture should NOT be made here, but in "Makeicl_local",
# which is included in ALL Makefiles of the program package.

# subpress echoing of command executing
#.SILENT :

export
# set top level directory
ifeq ($(TOP),)
   TOP = $(PWD)
endif
#
###### MACROS ######
#
#*** Macros, which depend e.g. on your platform and 
#    on where needed libraries and TAF are located:
#
# netCDF stuff
NETCDFINCDIR = /client/include
NETCDFLIB    = netcdf
NETCDFLIBDIR = /client/lib
INCDIRS = -I$(NETCDFINCDIR)     # include directories
LIBDIRS = -L$(NETCDFLIBDIR)     # libraries (only netCDF stuff)
LIBS    = $(NETCDFLIB)
# libraries for optimisation
LSOPT 	= /scratch/local1/adbethy/
LSOPT_tst	= /home/thomas/CCDAS/bethy/lsopt/head/src
BLAS1 	= /scratch/local1/adbethy/blas1
LAPACK 	= /home/users/rayner/lib/ifort64/lib
LSOPTLIB = lsopt_v2
LAPACKLIB = lapack
BLASLIB = blas
BLAS1LIB        = blas
OPTLIBDIRS = -L$(NETCDFLIBDIR) -L$(LSOPT) -L$(LAPACK) -L$(BLAS) -L$(TAFLIB)# libraries (for optimisation targets)
OPTLIBS    = -l$(NETCDFLIB) -l$(LSOPTLIB) -l$(LAPACKLIB) -l$(BLASLIB) # libraries (for optimisation targets)
#OPTLIBS    = -l$(NETCDFLIB) -llsopt_v2 -lblas1# libraries (for optimisation targets) 
# TAF stuff
TAFDIR	= /scratch/local1/m211024/taf
TAF	= staf -include netcdf.inc
#TAF	= taf
TAFLINK	= $(TAFDIR)/bin/taflink
TAFLIB	= $(TAFDIR)/../taf-1.4.7_ifc/lib
TAFFLAGS	= -r8 -f95 -f90decl -active-file -closetapes -v2 -newv# -showrecomp 
TTAFFLAGS	= -r8 -closetapes -active-file -v2
MADFLAGS	= -toplevel model -input x -output fc
GTBFLAGS	= -toplevel func -input x -output y 
VADFLAGS	= -toplevel func -input x -output y 
SVADFLAGS	= -toplevel func -input x -output y -pure -modmark mdasd

#TAFTAPE = FILE_TAPE
TAFTAPE = STATIC_TAPE

INPUT = FUTURE
INPUT = OBS


#*** Macros, set and changed by the users of CCDAS
# according to their purposes and environment
#
#HEW-ADD-050310 : Switch between different models rsp. source code directories
#/////////////////////////////
#* Model source code directory: 
#          bethy (default)  CARBON-BETHY coupled to TM2
#	   simple 	    Simple test	
CASE		= bethy
#CASE		= simple
#CASE		= structure
MODDIR		= $(TOP)/$(CASE)

#RESOLUTION      = hires
#RESOLUTION      = lores
RESOLUTION      = local

#* PCASE determines target quantities for prognostic mode
#PCASE		= GRID# yet to be implemented
#PCASE		= REGIONS# yet to be implemented
#PCASE		= BANDS
#PCASE		= NPP#
PCASE		= LOREGIONS
#PCASE		= NULL#
NGRID		= 2#
NYEAR		= 1#

#* Dimensions
# first for CASE = simple or structure
NN     =  6
NPADM   = 10
NPTLM	= $(NN)
ifeq ($(CASE),bethy)# extended BETHY version has 72 parameters
NN	= 6
NPTLM	= $(NN)
NPADM   = 1
ifeq ($(PCASE),NPP)#
NPADM  =  $(shell expr 12 \* $(NYEAR) \* $(NGRID))
endif
ifeq ($(PCASE),REGIONS)#
NPADM  =  9
endif
ifeq ($(PCASE),LOREGIONS)#
NPADM  =  12
endif
ifeq ($(PCASE),BANDS)#
NPADM  =  10
endif
# settings for background sensitivities
ILATB	= 1# 				initial lat band
ILATE	= 24# 				final lat band
NILAT	= 2# 				# of lat bands per itlm run
NITLM	= $(shell expr 36 \* $(NILAT))# # of independents
endif
# settings for future prog uncertainties
NFTLM = 11

# ASD options, dimensions can be computed from Jacobian dimensions
KIND		= 4
BITSIZE		= $(shell expr $(KIND) \* 8)
NSPAD		= $(shell expr \( 1 \+ \( $(NPADM) \- 1 \) / $(BITSIZE) \) )
NSPTL		= $(shell expr \( 1 \+ \( $(NPTLM) \- 1 \) / $(BITSIZE) \) )
TAF_SPAD	= -asd -orgmark _soad -admark _spad -asdkind $(KIND) -jacobian $(NSPAD) 
TAF_SPTL	= -asd -orgmark _sotl -ftlmark _sptl -asdkind $(KIND) -jacobian $(NSPTL) 

ifeq ($(FCASE),)
FCASE	= NULL
endif

# macro FCASE not set: normal diagnostic runs
# macro FCASE=PREPROG: normal run plus saving cs on disk for future prognostic run
# macro FCASE=FUTPROG: future prognostic run reading cs from disk plus saving 
#			npp, resf, interress on disk for net flux uncertainties  

DIAGOUT	= POST
# macro DIAGOUT not set: no diagnostic output (all targets except post)
# macro DIAGOUT=POST   : diagnostic output written (only for post)

# number of optimisations performed by mopti
BOPTI   = -1
EOPTI   = 1
EXP     = standard

#////////////////////////////////////////
#* Compiler :
#
# Preprocessing
CPP		= cpp
CPPFLAGS	= -P -traditional -DNFTLMY=$(NFTLMY) -DNPADM=$(NPADM) -DNPROJ=$(NGRID) -DNPTLM=$(NPTLM) -DNITLM=$(NITLM) -D$(PCASE) -D$(FCASE) -D$(RESOLUTION) -D$(DIAGOUT) -D$(TAFTAPE) -D$(INPUT)
#
# Choose the compiler: either "NAG", "FUJITSU" or "INTEL"
#COMPILER        = FUJITSU
COMPILER        = INTEL
ifeq ($(COMPILER),FUJITSU)# Fujitsu Fortran90/95 Compiler
	F90		= lf95 
	FC              = $(F90)
	DBL		= --dbl
	DEBUGFLAG	= -g --trace --ap --chk #a,e,s,u,x#--info
	RTECHECK	= --chkglobal
	OPTFLAG	        = -O3 #--tpp 
	F90FLAGS	=   $(INCDIRS) $(RTECHECK) $(DEBUGFLG) # --prefetch 
	F90FLAGS	=   $(INCDIRS) $(OPTFLAG)
	FFLAGS        	=  $(F90FLAGS)  
#       ADDFLAGS  = --wide -Cpp
else
ifeq ($(COMPILER),INTEL)# INTEL Compiler
#		F90             = ifc #$(DEBUG)
	F90             = ifort #$(DEBUG)
	FC              = $(F90)
	DBL             = -autodouble# turn off double precision
	DEBUGFLAG       = -g # -mp -C -inline_debug_info
	RTECHECK        = -C
#	OPTFLAG         = -O3#-xW -tpp7 
	F90FLAGS        =  $(INCDIRS) $(OPTFLAG)
	F90FLAGS        =  $(INCDIRS) $(DEBUGFLAG) 
	FFLAGS        	=  $(F90FLAGS)  
	LSOPT 	= /scratch/local1/adbethy/
	TAFDIR	= /scratch/local1/m211024/taf
else
	ERROR_MESSAGE = "Invalid Compiler $COMPILER" 
endif
endif

VERSION=$(shell svnversion)# current ccdas revison in svn
DIST=ccdas$(VERSION)
# Sources
SRCDIRS = $(DRIVERS) $(BETHY) $(UTIL) $(SIMPLE) $(DOCDIR) structure# fgtb source directories
#
BETHY = $(TOP)/bethy# bethy source code
SIMPLE = $(TOP)/simple# simple test program
UTIL = $(TOP)/util# utility programs
DRIVERS = $(TOP)/drivers# driver program code
DOCDIR = $(TOP)/doc# documents

IDLROUTINES=../idlroutines

#*** You can include an additioanal include Makefile
#    and override there Macros set above
#    MAKEICL should *MUST* be predefined and
#    allows easy switching when invoking make
ifdef MAKEICL
  $(shell ln -sf config/$(MAKEICL) Makeicl_local)
endif 
include Makeicl_local

# remove !
RM	= rm -f
MKDIR   = mkdir -p


# CCDAS excecutables 
EXE	= cost cost2 tstadm tsttlm tstvtlm tstgtb rev tsthesscol hesscol \
	func-gtb2 pfunc2 tstptlm tstpadm prev pfwd tstspadm \
	nextx fgrad opti-gtb post line dline ifunc2 tstitlm ifunc2 ifwd dfp dfptl \
	ffwd tstftlm ffunc2

# excecutables for utils
UEXE	= printp pfwdunc invhess parmap

#########################################################
#####  TARGETS ##########################################
#########################################################

#---------------------------------------------------------
# default rule print targets
#---------------------------------------------------------
default:
	@echo ''
	@echo '  Targets	Description'
	@echo '--------------------------------------------------------------'
	@echo ''
	@echo ' documentation'
	@echo '  [default] 	show all targets'
	@echo '  doc 		generates documentation of current version'
	@echo '  showInfo 	display current settings in the make description files'
	@echo ''
	@echo ' prepare tests and runs'
	@echo '  init		initialise run (create ouput directories)'
	@echo
	@echo ' sensitivities:'
	@echo '  rev		make executable to compute parameter sensitivities by adjoint'
	@echo ''
	@echo ' do optimisation:'
	@echo '  coldstart	prepare for new optimization'
	@echo '  nextx		make executable to compute next control vector'
	@echo '  fgrad		make executable to compute function and gradient'
	@echo '  post		make executable for prognostic model run after optimisation'
	@echo '  line		compute cost function between two control vectors'
	@echo '  dline		compute directional derivative between two control vectors'
	@echo '  dfp		run optimisation using DFP'
	@echo
	@echo ' evaluate Hessian matrix:'
	@echo '  hesscol	make executable to compute columns of Hessian'
	@echo
	@echo ' do prognostic step:'
	@echo '  prev		make executable to compute Jacobian for prognostic step (reverse)'
	@echo '  pfwd		make executable to compute Jacobian for prognostic step (forward)'
	@echo
	@echo ' sensitivities of optimal parameters wrt background:'
	@echo '  ifwd		make executable for forward mode calculation'
	@echo '  runifwd	run forward mode calculation'
	@echo
	@echo ' future prognostic runs:'
	@echo '  preprog	run post plus save cs (slow carbon pool) on disk'
	@echo '  postprog	run post for future period, reads cs from disk and'
	@echo ' 			saves npp, resfv,ressv, lfac on disk'
	@echo '  ffwd		make executable for forward mode calculation, ffwd reads'
	@echo '  			cs, npp, resfv, ressv, lfac from disk'
	@echo
	@echo
	@echo ' utility programs:'
	@echo '  printp	make executable to print parameter vector from OPWARMD'
	@echo '  parmap	make executable to check mapping of parameters read from OPWARMD'
	@echo '  invhess	make executable to invert and calculate eigen spectrum of the Hessian'
	@echo '  invhess_pseudo	make executable to add pseudo uncertainties to Hessian'
	@echo '  pfwdunc	make executable to calculate prognostic uncertainties from Jacobian'
	@echo '                 	computed in forward mode and Hessian'
	@echo '  func		make executable to calculate future prognostic uncertainties from'
	@echo '                 	Jacobians and Hessian'
	@echo
	@echo ' generate derivative code:'
	@echo '  adm		make adjoint model'
	@echo '  tlm		make tangent linear model'
	@echo '  hess		make Hessian times matrix code'
	@echo '  padm		make adjoint model for prognostic step'
	@echo '  ptlm		make tangent linear model for prognostic step'
	@echo '  spadm		make ASD reverse code for prognostic step'
	@echo '  sptlm		make ASD forward code for prognostic step'
	@echo '  ftlm		make tangent linear code for future prog uncertainties'
	@echo ''
	@echo ' testing/debugging for individual codes'
	@echo '  cost		make executable to compute cost function'
	@echo '  cost2		make executable to compute cost function twice'
	@echo '  pfunc2		make executable to evaluate function for prognostic step'
	@echo '  ifunc2		make executable to evaluate function for background sens'
	@echo '  tstadm  	combine timing and derivative check for adjoint'
	@echo '  tsttlm  	combine timing and derivative check for tangent'
	@echo '  tstspadm	combine timing and ASD check for prognostic adjoint'
	@echo '  tstpadm	combine timing and derivative check for prognostic adjoint'
	@echo '  tstptlm	combine timing and derivative check for prognostic tangent'
	@echo '  tstitlm	combine timing and derivative check for background sens tangent'
	@echo '  tstftlm	combine timing and derivative check for future prog tangent'
	@echo '  tsthesscol	combine timing and derivative check for columns of hessian'
	@echo
	@echo ' run series of test'
	@echo '  testall	generate and verify all derivative code'
	@echo
	@echo ' cleaning up:'
	@echo '  clearo	remove object files'
	@echo '  clean		remove intermediate files'
	@echo '  scratch	remove all generated files'
	@echo '  rmall		remove all except source code'
	@echo
	@echo ' pack a distribution:'
	@echo '  dist		pack tarballs ccdas$(VERSION).tgz and input.tgz'
#
#---------------------------------------------------------
# targets to create and show the CCDAS documention
#---------------------------------------------------------
#
# create documentation
doc: init
	cd $(DOCDIR) ; $(MAKE) pdf

showInfo:
	@echo '== Current settings: ========'
	@echo '- Version         $(VERSION) '
	@echo '- Model dir       $(MODDIR) '
	@echo '- Compiler        $(COMPILER) '
	@echo '- F90 =      $(F90)'
	@echo '- F90FLAGS = $(F90FLAGS)'
	@echo '- File control_bethy/control:'
	@cat control_bethy/control
	@echo '- !!!! NOTE: File adbethy.options NOT LONGER USED, HEW050307 !!!!:'
####	@cat control_bethy/adbethy.options


#---------------------------------------------------------
# generate derivative code:
#---------------------------------------------------------

# adjoint model
adm:
	cd $(MODDIR) ; $(MAKE) $@

# tangent linear model
tlm:
	cd $(MODDIR) ; $(MAKE) $@

# adjoint model for GTB
gtb:
	cd $(MODDIR) ; $(MAKE) $@

# Hessian times matrix code
hess: 
	cd $(MODDIR) ; $(MAKE) $@ 

# adjoint model for prognostic step
padm: 
	cd $(MODDIR) ; $(MAKE) $@

# tangent linear model for prognostics 
ptlm:
	cd $(MODDIR) ; $(MAKE) $@

# ASD reverse code for prognostic step
spadm :
	  cd $(MODDIR) ; $(MAKE) $@

# ASD forward code for prognostic step
sptlm:
	  cd $(MODDIR) ; $(MAKE) $@

# implicit tangent code
itlm: 
	cd $(MODDIR) ; $(MAKE) $@

# future prognostic uncertaities tangent code
ftlm: 
	cd $(MODDIR) ; $(MAKE) $@

# force rebuilt for targets with equal named executables
force_rebuild:

# initialise run (ouput directories)
init:
	@ ( if [ ! -d output ]; then mkdir output; fi; )
	@ ( if [ ! -d optiout ]; then mkdir optiout; fi; )

#---------------------------------------------------------
# testing/debugging for individual codes
#---------------------------------------------------------

# compute cost function
cost: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# compute cost function twice
cost2: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# evaluate future prog uncertainties function forward
ffunc2: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# evaluate implicit function forward
ifunc2: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# evaluate scalar valued function/prognostic?? twice
pfunc2: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# evaluate scalar valued function/prognostic?? twice
func-gtb2: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

#   compute function and gradient by adjoint
rev: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative check for adjoint
tstadm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative check for tangent linear
tsttlm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative check for tangent linear in vector mode
tstvtlm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative check for adjoint
tstgtb: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative check for progn adjoint
tstpadm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and ASD check for progn adjoint
tstspadm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative check for progn tangent
tstptlm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative check for columns of hessian
tsthesscol: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative for implicit tangent
tstitlm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# combine timing and derivative for future prog  tangent
tstftlm: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

#---------------------------------------------------------
# testing/debugging for individual codes
#---------------------------------------------------------
include Make-tst.inc
alltests: 
	@echo local test results  > alltests.log
	$(MAKE) testlocal
	@echo lores test results  >> alltests.log
	$(MAKE) testlores
	@cat alltests.log

testlocal: 
	ln -fs control_bethy/control_hainich_short control
	$(MAKE) xall "RESOLUTION=local" "VERSION=hainich_short" "PCASE=NPP"
#	$(MAKE) xonly "XTARGET=xcost2 xadm xhess" "RESOLUTION=local" "PCASE=NPP" "VERSION=hainich_short" "TAFTAPE=FILE_TAPE"

testlores: 
	ln -fs control_bethy/control.lores.1year control
#	ln -fs control_bethy/control.lores.2years control
#	$(MAKE) xonly "XTARGET=xcost2 xtlm xadm xhess xvectlm xvecadm" "RESOLUTION=lores" "VERSION=lores_short" "PCASE=BANDS"
#	$(MAKE) xonly "XTARGET=xtlm xadm" "RESOLUTION=lores" "VERSION=lores_short" "PCASE=BANDS"
	$(MAKE) xall "RESOLUTION=lores" "VERSION=lores_short" "PCASE=BANDS"

testhires: 
	ln -fs control_bethy/control.hires.1year control
	$(MAKE) xonly "XTARGET=xcost2 xtlm xadm" "RESOLUTION=hires" "VERSION=hires_short" "PCASE=BANDS"
#	$(MAKE) xonly "XTARGET=xcost2" "RESOLUTION=hires" "VERSION=hires_short" "PCASE=BANDS"
#	$(MAKE) xall "RESOLUTION=hires" "VERSION=hires_short" "PCASE=BANDS"

#---------------------------------------------------------
# do optimisation
#---------------------------------------------------------

# prepare for new optimization
coldstart:
	$(RM) OPWARMI OPWARMD lsopt.tmp result 

# compute next control vector
nextx: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

#  compute function and gradient
fgrad: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

#  compute function and gradient
opti-gtb: init libfgtb.a
	cd $(MODDIR) ; $(MAKE) ../$@ 
# build gtb lib
libfgtb.a: 
	cd fgtb; $(MAKE) FC=$(FC) "FFLAGS=$(FFLAGS) $(DBL)"

# build dfp lib
libnumrec.a: 
	@ ( if [ ! -d libnumrec.a ]; then cd numrec; $(MAKE) FC=$(FC) "FFLAGS=$(FFLAGS) $(DBL)"; fi; )

# print current control vector and gradient
post: init
	$(RM) post $(MODDIR)/bethy.f90 $(MODDIR)/diagnostics.f90 $(MODDIR)/mo_diagnostics.f90  
	cd $(MODDIR) ; $(MAKE) CPPFLAGS+=" -DPOST" ../post  
	#$(RM) $(MODDIR)/bethy.f90 $(MODDIR)/diagnostics.f90 $(MODDIR)/mo_diagnostics.f90  

# do a linesearch between two control vectors
line: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

dline: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# run an optimiasation
lsopt: init
	cd $(MODDIR) ; $(MAKE) $@ 

# run an optimiasation with DFP
dfp: init libnumrec.a
	cd $(MODDIR) ; $(MAKE) ../$@ 

# run an ensemble of optimisations
loop2: dfp
	$(MKDIR) outmopti
	iopti=$(BOPTI); while test $$iopti -le $(EOPTI) ; do \
		outdir=outmopti/opti-$(EXP)-$$iopti ; \
		$(MKDIR) $$outdir ; \
		echo $$iopti > pert.dat ; \
		echo running opti no: $$iopti ; \
		iloop=0; while test ! -e xf.b ; do \
			echo running loop pass: $$iloop ; \
			dfp 2>&1 >> dfp.out ; \
			cp px.b x.b; \
			$(RM) pert.dat ; \
			let "iloop = iloop + 1"; \
		done ; \
		mv dfp.out h.dat x.b xf.b lx.b OPWARM* $$outdir ; \
		let "iopti = iopti + 1"; \
	done

mopti: dfp
	$(MKDIR) outmopti
	iopti=$(BOPTI); while test $$iopti -le $(EOPTI) ; do \
		outdir=outmopti/opti-$(EXP)-$$iopti ; \
		$(MKDIR) $$outdir ; \
		echo $$iopti > pert.dat ; \
		echo running opti no: $$iopti ; \
		dfp >& dfp.out ; \
		mv dfp.out h.dat x.b xf.b lx.b OPWARM* $$outdir ; \
		let "iopti = iopti + 1"; \
	done

popti: dfp
	- $(MKDIR) outmopti
	iopti=$(BOPTI); while test $$iopti -le $(EOPTI) ; do \
		outdir=outmopti/opti-$(EXP)-$$iopti ; \
		$(MKDIR) $$outdir/output ; \
		cp control $$outdir ; \
		cd $$outdir ; cp ../../dfp ./; \
		ln -s ../../input ; \
		ln -s ../../control_bethy ; \
		ln -s ../../dfpmin.par ; \
		echo $$iopti > pert.dat ; \
		echo now starting opti no: $$iopti ; \
		(dfp >& dfp.out &); \
		let "iopti = iopti + 1"; \
		cd ../.. ; \
	done

# run an opti with pseudo fluxes
pseudo: init libnumrec.a
	make scratch
	cd $(MODDIR) ; $(MAKE) CPPFLAGS+=" -DDIAGNOSTIC" ../cost  
	echo 'pseudo: generating pseudo data'
	./cost
	make scratch
	cd $(MODDIR) ; $(MAKE) CPPFLAGS+=" -DPSEUDOFLUX" ../cost
	echo 'pseudo: checking pseudo data'
	./cost >& cost.log
	grep cost_flux cost.log
	cd $(MODDIR) ; $(MAKE) CPPFLAGS+=" -DPSEUDOFLUX" ../tstadm
	echo 'pseudo: checking ADM'
	./tstadm >& tstadm.log
	-grep DIFFERENT tstadm.log
	cd $(MODDIR) ; $(MAKE) CPPFLAGS+=" -DPSEUDOFLUX" ../dfp
	echo 'pseudo: running optimisation'
	./dfp 2>&1 | tee $@.log

# run an optimiasation with DFP
dfptl: init libnumrec.a
	cd $(MODDIR) ; $(MAKE) ../$@ 

#-------------------------------------------------------
# evaluate Hessian matrix:
#---------------------------------------------------------

# compute columns of the Hessian
hesscol: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

#-------------------------------------------------------
# do prognostic step:
#---------------------------------------------------------

# compute Jacobian for prognostic step (reverse)
prev: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# compute Jacobian for prognostic step (forward)
pfwd: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

#------------------------------------------------------
# future prognostics:
#--------------------------------------------------------

# do post for optimisation period and save cs
preprog-do: post

preprog: init
	$(RM) $(MODDIR)/bethy.f90
	$(MAKE) FCASE=PREPROG preprog-do

# do post for whole period (incl. prognostic period), need to read in cs
# and saves npp,ress,interress for needed for calculating future prognostic uncertainties
postprog-do: post

postprog: init
	$(RM) $(MODDIR)/bethy.f90
	$(MAKE) FCASE=FUTPROG postprog-do

# executable for computing dNEP/dcs for future prognostic uncertainties 
ffwd: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

#-------------------------------------------------------
# implicit function:
#---------------------------------------------------------

# executable for computing d2Jdxdb 
ifwd: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# executable for computing d2Jdxdb 
runifwd: ifwd init
	ilat=$(ILATB); while test $$ilat -le $(ILATE) ; do \
		echo $$ilat > latb.par ; \
		cat latb.par ; \
		./ifwd ; \
		cp output/d2Jdxdb.bin output/d2Jdxdb-$$ilat.bin ; \
		let "ilat = ilat + $(NILAT)"; \
	done

# compute d2Jdxdb 
tifwd: init
	cd $(MODDIR) ; $(MAKE) ../ifwd

#---------------------------------------------------------
# rules for utility programms
#---------------------------------------------------------

#  print current parameter vector (old version, without params transformation)
printparams: force_rebuild
	cd $(MODDIR) ; $(MAKE) ../$@ 

#  plot a priori parameter (PDF, Trafo, Gaussian and Costfunction)
plotparams: force_rebuild
	cd $(IDLROUTINES) ; ./batch_boundmap

#  print current parameter vector
printp: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# test parameter mapping
parmap: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# inverte and calculate eigen spectrum of the Hessian
invhess: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# inverte and calculate eigen spectrum of the Hessian
invhess_pseudo: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

jachess: init
	cd $(MODDIR) ; $(MAKE) ../$@ 

# calculates prognostic uncertainties from Jacobian and Hessian
pfwdunc: init
	cd $(UTIL) ; $(MAKE) ../$@ 

# calculates future prognostic uncertainties from Jacobians and Hessian
func: init
	cd $(UTIL) ; $(MAKE) ../$@ 

# calculates background sensitivities 
dpdb: init
	cd $(UTIL) ; $(MAKE) ../$@ 

#---------------------------------------------------------
# targets for removing and cleaning
#---------------------------------------------------------

# remove .o, .mod rsp .d files
clearo:
	@for DIR in $(BETHY) $(SIMPLE) $(DRIVERS) util structure;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) clearo ; cd $$back ;) ;\
	  done
#	@$(RM) *.o *.d *.mod

# remove intermediate files, warmstart still possible
clean:	clearo
	@$(RM) monthfile_* g_monthfile_* _hesmonthfile* flux_temp.bin 
	@$(RM) outfile_* _hesoutfile_?_bethy_z* fort.* 
	@$(RM) *carbon_* *scale_* *day_* *diurnal_* *cbalance_* *state.b
	@$(RM) pert.dat

# remove generated files, only coldstart possible
scratch: clean
	@for DIR in $(BETHY) $(SIMPLE) $(DRIVERS) structure;\
	  do \
	  (cd $$DIR ;\
	  $(MAKE) scratch ; cd $$back ;) ;\
	  done
	@$(RM) output/* optiout/*
	@$(RM) OPWARM* *.tmp result params_out.txt
	@$(RM) error.dat debugging.dat
	@$(RM) x.b xf.b lx.b h.dat line.set dline.set parlist.dat parset.dat
	@$(RM) $(EXE) $(UEXE) $(TEXE)

rmall:	scratch
	@for DIR in $(SRCDIRS) ;\
	do \
	( cd $$DIR ; $(MAKE) rmall ; cd $$back ;) ;\
	done
	@$(RM) *~ \#*\# *.out outopti* doc/*~  doc/\#*\# 
#	@$(RM) -r optiout output obj ccdas$(VERSION) $(DIST)
	@$(RM) input.tgz ccdas$(VERSION).tgz
	@$(RM) *.bin
	@$(RM) numrec/*.o
	@$(RM) libnumrec.a


#---------------------------------------------------------
# target for packing a distribution (with script pack_dist)
#---------------------------------------------------------
dist: adm tlm hess padm ptlm spadm sptlm
	@$(RM) -r $(DIST)
	@$(MKDIR) $(DIST)/input/gv
	@$(MKDIR) $(DIST)/input/jactd
	@$(MKDIR) $(DIST)/input/background
	@$(MKDIR) $(DIST)/input/climate

	@cp CHANGES PROBLEMS QUESTIONS CHECKS README_dist Makefile Makeicl* \
		all_targets \
		opti_save linesearch hesse pjac trunopti \
		jaccol.par pjaccol.par lsopt.par tst.par $(DIST)
	@cp -a control_bethy $(DIST)
	@cp -a bethy doc drivers simple util $(DIST) 
	@cp input/background/co2_flx.atm.cdf input/background/fossil*.nc input/background/fossil_amp.txt \
		input/background/ocean_tak.nc input/background/luc_flux.nc $(DIST)/input/background 
	@cp input/climate/b_climate3veg_79-00_bethy13.5.nc input/climate/b_lores_climate3veg_79-00.nc \
		input/climate/b_lores_spinup_climate3veg_79-00.nc input/climate/b_spinup_climate3veg_79-00.nc \
		$(DIST)/input/climate
	@(for stat in $$(cat input/gv41_data); do \
		cp input/gv/$${stat}*.co2 $(DIST)/input/gv ; \
		cp input/jactd/mat_$${stat}*.nc $(DIST)/input/jactd; \
	done)

	@tar -zcvf $(DIST).tgz $(DIST)

#---------------------------------------------------------
# target for packing a distribution (with script pack_dist)
#---------------------------------------------------------
dist-gtb: gtb
	@$(RM) -r $(DIST)
	@$(MKDIR) $(DIST)/input/gv
	@$(MKDIR) $(DIST)/input/jactd
	@$(MKDIR) $(DIST)/input/background
	@$(MKDIR) $(DIST)/input/climate
	@cp README_dist Makefile Make-tst.inc Makeicl_local \
		tstgtb.par $(DIST)
	@cp -a control_bethy $(DIST)
	@cp -a bethy doc drivers util $(DIST) 
	@cp input/background/co2_flx.atm.cdf input/background/fossil*.nc input/background/fossil_amp_multi.txt \
		input/background/ocean_tak.nc input/background/luc_flux.nc $(DIST)/input/background 
	@cp input/climate/b_lores_spinup* input/climate/b_lores_climate* \
		$(DIST)/input/climate
	@(for stat in $$(cat control_bethy/gv41_data); do \
		cp input/gv/$${stat}*.co2 $(DIST)/input/gv ; \
		cp input/jactd/mat_$${stat}*.nc $(DIST)/input/jactd; \
	done)
	@tar -zcvf $(DIST).tgz $(DIST)
