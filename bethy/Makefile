#*******************************************************************
# FILENAME  : bethy/Makefile
# DATE      : 10/03
# DESIGNER  : CCDAS group
# PURPOSE   : generate code and build CCDAS objects from BETHY sourcecode
#*******************************************************************

# This is the Makefile of the BETHY source code. 

#---------------------------------------------------------
# code and objects
#---------------------------------------------------------
SRC = 	cbalance.f90 	 	budget.f90	\
	diagnostics.f90 	\
	surface.f90 		numbmod.f90		initmod.f90		\
	bethy.f90		model.f90		postmod.f90		\
	postadm.f90		\
	veg_allocate.f90 		veg_deallocate.f90

MSRC	=			\
	constants_pjr.f90	mo_taf.f90		mo_netcdf.f90		\
	mo_constants.f90	mo_grid.f90		mo_namelist.f90		\
	mo_calendar.f90		\
	mo_carparams.f90	mo_carvar.f90				\
	mo_ctrl.f90 		mo_beta.f90		mo_vegetation.f90	\
	tm.f90 			mo_climate.f90					\
	mo_diagnostics.f90	mo_helper.f90		mo_hydro.f90		\
	mo_pheno.f90		mo_mapping.f90 		mo_trafo.f90 		\
	fluo_func.f90		fluo_param.f90		mo_rtmo.f90 		\
	chemical.f90		mo_rtmf.f90		fluo.f90		\
	integration.f90 	costf.f90 			\
	mo_prog.f90		mo_carbon.f90		\
	mo_io.f90		mo_surface.f90	\
	bethy2tm.f90		bgr.f90			mo_config.f90 \

OBJ	= $(SRC:.f90=.o)
MOBJ	= $(MSRC:.f90=.o)

# DEPENDENCIES of Module-FILES
include ../Makeicl_local

#---------------------------
# pattern rules for bethy sources
#---------------------------
%.o: %.f
	$(F90) $(F90FLAGS) $(DBL) -c $(@D)/$(<F) -o $(OBJ)/$(@F)

%.o: %.f90
	$(F90) $(F90FLAGS) $(DBL) -c -o $@ $<

%.mod: %.f90
	$(F90) $(F90FLAGS) $(DBL) -c -o $@ $<

%.d: %.f90
	$(F90) $(F90FLAGS) $(DBL) -c -o $@ $<

%.f90: %.g90
	$(F90) $(F90FLAGS) $(DBL) -c -o $@ $<

# disable default suffix rule for .mod files
%.o: %.mod
	true

#---------------------------
# pattern rules for drivers
#---------------------------
$(DRIVERS)/%.o: $(DRIVERS)/%.f90
	$(F90) $(F90FLAGS) $(DBL) -c -o $@ $<

#---------------------------------------------------------
# AD code and objects
#---------------------------------------------------------
TLSRC = mo_constants_tl.f90	mo_grid_tl.f90		 mo_namelist_tl.f90	\
	mo_vegetation_tl.f90					\
	mo_surface_tl.f90	mo_beta_tl.f90		\
	mo_carvar_tl.f90	mo_carparams_tl.f90	mo_netcdf_tl.f90	\
	mo_calendar_tl.f90	mo_ctrl_tl.f90 		constants_pjr_tl.f90	\
        tm_tl.f90		\
	mo_prog_tl.f90		mo_climate_tl.f90	mo_helper_tl.f90	\
	mo_hydro_tl.f90		mo_diagnostics_tl.f90	mo_pheno_tl.f90		\
	mo_mapping_tl.f90	mo_trafo_tl.f90		\
	costf_tl.f90		veg_deallocate_tl.f90	\
	veg_allocate_tl.f90 	surface_tl.f90		budget_tl.f90		\
	mo_io_tl.f90		mo_config_tl.f90 	postmod_tl.f90		\
	postadm_tl.f90		bgr_tl.f90		bethy2tm_tl.f90		\
	mo_taf_tl.f90		numbmod_tl.f90		mo_carbon_tl.f90	\
	model_tl.f90		bethy_tl.f90		initmod_tl.f90		\
	diagnostics_tl.f90	cbalance_tl.f90		

VTLSRC = taf_vtl.f90 $(TLSRC:_tl.f90=_vtl.f90)

ADSRC = mo_constants_ad.f90	mo_grid_ad.f90		 mo_namelist_ad.f90	\
	mo_vegetation_ad.f90					\
	mo_surface_ad.f90	mo_beta_ad.f90		\
	mo_carvar_ad.f90	mo_carparams_ad.f90	mo_netcdf_ad.f90	\
	mo_calendar_ad.f90	mo_ctrl_ad.f90 		constants_pjr_ad.f90	\
        tm_ad.f90		\
	mo_prog_ad.f90		mo_climate_ad.f90	mo_helper_ad.f90	\
	mo_hydro_ad.f90		mo_diagnostics_ad.f90	mo_pheno_ad.f90		\
	mo_mapping_ad.f90	mo_trafo_ad.f90		costf_ad.f90		\
	veg_deallocate_ad.f90	\
	veg_allocate_ad.f90 	surface_ad.f90		budget_ad.f90		\
	mo_io_ad.f90		mo_config_ad.f90 	postmod_ad.f90		\
	postadm_ad.f90		bgr_ad.f90		bethy2tm_ad.f90		\
	mo_taf_ad.f90		numbmod_ad.f90		mo_carbon_ad.f90	\
	cbalance_ad.f90		bethy_ad.f90		model_ad.f90		initmod_ad.f90		\
	diagnostics_ad.f90						

GTBSRC	= \
	mo_constants_gtb.f90	mo_grid_gtb.f90		 mo_namelist_gtb.f90	\
	mo_vegetation_gtb.f90	\
	mo_surface_gtb.f90	mo_beta_gtb.f90		\
	mo_carvar_gtb.f90	mo_carparams_gtb.f90	mo_netcdf_gtb.f90	\
	mo_calendar_gtb.f90	mo_ctrl_gtb.f90 		constants_pjr_gtb.f90	\
        tm_gtb.f90		mo_climate_gtb.f90	mo_helper_gtb.f90	mo_hydro_gtb.f90		\
	mo_diagnostics_gtb.f90	mo_pheno_gtb.f90	mo_mapping_gtb.f90	\
	mo_trafo_gtb.f90	costf_gtb.f90		veg_deallocate_gtb.f90	\
	veg_allocate_gtb.f90 	surface_gtb.f90		mo_prog_gtb.f90		\
	mo_io_gtb.f90		mo_config_gtb.f90 	postmod_gtb.f90		\
	postadm_gtb.f90		bgr_gtb.f90		bethy2tm_gtb.f90		\
	mo_taf_gtb.f90		numbmod_gtb.f90		mo_carbon_gtb.f90	\
	bethy_gtb.f90		initmod_gtb.f90		\
	diagnostics_gtb.f90				\
	cbalance_gtb.f90	model_gtb.f90		func-gtb_gtb.f90

PTLSRC	= taf_ptl.f90 \
	mo_constants_ptl.f90	mo_grid_ptl.f90		 mo_namelist_ptl.f90	\
	mo_vegetation_ptl.f90	\
	mo_surface_ptl.f90	mo_beta_ptl.f90		\
	mo_carvar_ptl.f90	mo_carparams_ptl.f90	mo_netcdf_ptl.f90	\
	mo_calendar_ptl.f90	mo_ctrl_ptl.f90 	constants_pjr_ptl.f90	\
        tm_ptl.f90		mo_climate_ptl.f90	mo_helper_ptl.f90	mo_hydro_ptl.f90	\
	mo_diagnostics_ptl.f90	mo_pheno_ptl.f90	mo_mapping_ptl.f90	\
	mo_trafo_ptl.f90	costf_ptl.f90		veg_deallocate_ptl.f90	\
	veg_allocate_ptl.f90 	surface_ptl.f90		mo_prog_ptl.f90		\
	mo_io_ptl.f90		mo_config_ptl.f90 	postmod_ptl.f90		\
	postadm_ptl.f90		bgr_ptl.f90		bethy2tm_ptl.f90	\
	mo_taf_ptl.f90		numbmod_ptl.f90		mo_carbon_ptl.f90	\
	model_ptl.f90		bethy_ptl.f90		initmod_ptl.f90		\
	diagnostics_ptl.f90				budget_ptl.f90		\
	cbalance_ptl.f90	pfunc_ptl.f90

PADSRC	= taf_pad.f90 \
	mo_constants_pad.f90	mo_grid_pad.f90		 mo_namelist_pad.f90	\
	mo_vegetation_pad.f90	\
	mo_surface_pad.f90	mo_beta_pad.f90		\
	mo_carvar_pad.f90	mo_carparams_pad.f90	mo_netcdf_pad.f90	\
	mo_calendar_pad.f90	mo_ctrl_pad.f90 		constants_pjr_pad.f90	\
        tm_pad.f90		mo_climate_pad.f90	mo_helper_pad.f90	mo_hydro_pad.f90		\
	mo_diagnostics_pad.f90	mo_pheno_pad.f90	mo_mapping_pad.f90	\
	mo_trafo_pad.f90	costf_pad.f90		veg_deallocate_pad.f90	\
	veg_allocate_pad.f90 	surface_pad.f90		mo_prog_pad.f90		\
	mo_io_pad.f90		mo_config_pad.f90 	postmod_pad.f90		\
	postadm_pad.f90		bgr_pad.f90		bethy2tm_pad.f90		\
	mo_taf_pad.f90		numbmod_pad.f90		mo_carbon_pad.f90	\
	bethy_pad.f90		initmod_pad.f90		\
	diagnostics_pad.f90				budget_pad.f90		\
	cbalance_pad.f90	pfunc_pad.f90		model_pad.f90		

HESSRC=   taf_hes.f90 $(ADSRC:_ad.f90=_ad_hes.f90) func_hes.f90	


SPADSRC	= taf_spad.f90 mo_vegetation_spad.f90    mo_diagnostics_spad.f90		\
	  mo_carvar_spad.f90        mo_carbon_spad.f90	\
	  mo_beta_spad.f90          mo_surface_spad.f90       \
	  cbalance_spad.f90         mo_climate_spad.f90	climinday_spad.f90	\
          tm_spad.f90               mo_pheno_spad.f90		\
	  mo_hydro_spad.f90	    mo_mapping_spad.f90		mo_trafo_spad.f90	\
	  diagnostics_spad.f90      costf_spad.f90            surface_spad.f90	\
	  bethy2tm_spad.f90         bgr_spad.f90		\
	  veg_allocate_spad.f90	  	veg_deallocate_spad.f90			\
	  model_spad.f90 		  climate_allocate_spad.f90	\
	  climate_deallocate_spad.f90 pfunc_spad.f90

SPTLSRC	= taf_sptl.f90 mo_vegetation_sptl.f90    mo_diagnostics_sptl.f90		\
	  mo_carvar_sptl.f90        mo_carbon_sptl.f90	\
	  mo_beta_sptl.f90          mo_surface_sptl.f90       \
	  cbalance_sptl.f90         mo_climate_sptl.f90	climinday_sptl.f90	\
          tm_sptl.f90               mo_pheno_sptl.f90		\
	  mo_hydro_sptl.f90	   mo_mapping_sptl.f90		mo_trafo_sptl.f90	\
	  diagnostics_sptl.f90      costf_sptl.f90            surface_sptl.f90	\
	  bethy2tm_sptl.f90         bgr_sptl.f90		\
	  veg_allocate_sptl.f90	  	veg_deallocate_sptl.f90			\
	  model_sptl.f90 		  climate_allocate_sptl.f90	\
	  climate_deallocate_sptl.f90 pfunc_sptl.f90

ITLSRC=   taf_itl.f90 $(ADSRC:_ad.f90=_ad_itl.f90) ifunc_itl.f90	

FSRC	=			\
	mo_constants.f90	mo_diagnostics.f90 	mo_grid.f90		\
	mo_namelist.f90		mo_ctrl.f90		mo_vegetation.f90	\
	mo_beta.f90		mo_mapping.f90 		mo_trafo.f90		\
	mo_netcdf.f90 		\
	mo_calendar.f90		constants_pjr.f90	tm.f90 			\
	costf.f90 		mo_carparams.f90	mo_climate.f90		\
	mo_hydro.f90		mo_pheno.f90		mo_config.f90		\
	budget.f90		

FTLSRC=   taf_ftl.f90 $(FSRC:.f90=_ftl.f90) ffunc_ftl.f90	


TLOBJ	= $(TLSRC:.f90=.o)
VTLOBJ	= $(VTLSRC:.f90=.o)
ADOBJ	= $(ADSRC:.f90=.o)
GTBOBJ	= $(GTBSRC:.f90=.o)

PTLOBJ	= $(PTLSRC:.f90=.o)
PADOBJ	= $(PADSRC:.f90=.o)

HESOBJ= $(HESSRC:.f90=.o)

SPADOBJ	= $(SPADSRC:.f90=.o)
SPTLOBJ	= $(SPTLSRC:.f90=.o)

ITLOBJ= $(ITLSRC:.f90=.o)

FOBJ	= $(FSRC:.f90=.o)
FTLOBJ= $(FTLSRC:.f90=.o)

#---------------------------------------------------------
# make netcf code
# (mo_netcdf must be compiled without the double precision flag)
#---------------------------------------------------------
mo_netcdf.o : mo_netcdf.f90
	$(F90) $(F90FLAGS) -c mo_netcdf.f90

mo_netcdf_tl.o : mo_netcdf_tl.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_tl.f90

mo_netcdf_vtl.o : mo_netcdf_vtl.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_vtl.f90

mo_netcdf_ad.o : mo_netcdf_ad.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_ad.f90

mo_netcdf_ad_hes.o : mo_netcdf_ad_hes.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_ad_hes.f90

mo_netcdf_ad_itl.o : mo_netcdf_ad_itl.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_ad_itl.f90

mo_netcdf_ptl.o : mo_netcdf_ptl.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_ptl.f90

mo_netcdf_ftl.o : mo_netcdf_ftl.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_ftl.f90

mo_netcdf_gtb.o : mo_netcdf_gtb.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_gtb.f90

mo_netcdf_pad.o : mo_netcdf_pad.f90
	$(F90) $(F90FLAGS) -c mo_netcdf_pad.f90

mo_netcdf.mod : mo_netcdf.f90
	$(F90) $(F90FLAGS) -c mo_netcdf.f90

mo_netcdf.d : mo_netcdf.f90
	$(F90) $(F90FLAGS) -c mo_netcdf.f90

# disable default suffix rule for .mod files
.mod.o:
	true

#---------------------------------------------------------
# generate derivative code:
#---------------------------------------------------------

# tangent and adjoint models
tlm $(TLSRC):   $(MSRC) $(SRC)
	$(TAF) -l taf-tlm.log $(TAFFLAGS) $(MADFLAGS) -forward $(INCDIRS)	\
	 $(MSRC) $(SRC) 

vtlm $(VTLSRC):   $(MSRC) $(SRC)
	$(TAF) -l  taf-vtlm.log $(TAFFLAGS) $(MADFLAGS) -forward $(INCDIRS)	\
	 -ldg -vecfor $(NN) -ftlmark _vtl $(MSRC) $(SRC) 

adm $(ADSRC): $(MSRC) $(SRC)
	$(TAF) -l taf-adm.log $(TAFFLAGS) $(MADFLAGS) -reverse $(INCDIRS) \
	$(MSRC) $(SRC) 

# adm for GTB
gtb $(GTBSRC): $(MSRC) $(SRC) func-gtb.f90
	$(TAF) -l taf-gtb.log $(TAFFLAGS) $(GTBFLAGS) -reverse  -admark _gtb $(INCDIRS) \
	$(MSRC) $(SRC) func-gtb.f90 

# Hessian times matrix code
hess $(HESSRC):  $(ADSRC) func.f90 
	$(TAF) -l taf-hess.log $(TAFFLAGS) -toplevel func -input x -output y 	\
	-forward  -ftlmark _hes $(INCDIRS) -ldg -vecfor $(NN)	\
	 taf_ad.f90 $(ADSRC) func.f90  

# adjoint model for prognostic step
padm $(PADSRC):   $(MSRC) $(SRC) pfunc.f90
	$(TAF) -l taf-padm.log $(TAFFLAGS) $(VADFLAGS) -reverse  -admark _pad $(INCDIRS)	\
	-vecrev $(NPADM) -ldg  $(MSRC) $(SRC) pfunc.f90 

# tangent linear model for prognostics 
ptlm $(PTLSRC):   $(MSRC) $(SRC) pfunc.f90
	$(TAF) -l taf-ptlm.log $(TAFFLAGS) $(VADFLAGS) -forward  -ftlmark _ptl $(INCDIRS)	\
	-vecfor $(NPTLM) -ldg  $(MSRC) $(SRC) pfunc.f90 

# ASD reverse code for prognostic step
spadm $(SPADSRC):   $(MSRC) $(SRC) pfunc.f90
	$(TAF) -l taf-spadm.log $(TAFFLAGS) $(TAF_SPAD) $(SVADFLAGS) -reverse $(INCDIRS)	\
	 $(MSRC) $(SRC) pfunc.f90 

# ASD forward code for prognostic step
sptlm $(SPTLSRC):   $(MSRC) $(SRC) pfunc.f90
	$(TAF) -l taf-sptlm.log $(TAFFLAGS) $(TAF_SPTL) $(SVADFLAGS) -forward $(INCDIRS)	\
	 $(MSRC) $(SRC) pfunc.f90 

# implicit tangent code
itlm $(ITLSRC):  ifunc.f90 $(ADSRC)
	$(TAF) -l taf-itlm.log $(TAFFLAGS) -toplevel func -input xbgr -output y 	\
	-forward  -ftlmark _itl $(INCDIRS) -ldg -vecfor $(NITLM)	\
	taf_ad.f90  $(ADSRC) ifunc.f90   

# tangent linear model for future prognostic uncertainties
ftlm $(FTLSRC):   $(FSRC) ffunc.f90
	$(TAF) -l taf-ftlm.log $(TAFFLAGS) $(VADFLAGS) -forward  -ftlmark _ftl $(INCDIRS)	\
	-vecfor $(NFTLM) -ldg  $(FSRC) ffunc.f90 

#---------------------------------------------------------
# testing/debugging for individual codes
#---------------------------------------------------------

# compute cost function
../cost: $(DRIVERS)/prgcost.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ)
	$(F90) $(LIBDIRS) -o $@	\
	$(DRIVERS)/storef.o $(DRIVERS)/prgcost.o \
	$(MOBJ) $(OBJ) $(LIBS)

# compute cost function twice
../cost2: $(DRIVERS)/prgcost2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgcost2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) $(LIBS)

# compute function and gradient by adjoint
../rev: $(DRIVERS)/prgrev.o $(DRIVERS)/storef.o $(ADOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@      \
	$(DRIVERS)/prgrev.o $(DRIVERS)/storef.o $(ADOBJ) $(LIBS)

# combine timing and derivative check for adjoint
../tstadm: $(DRIVERS)/prgtstadm.o $(DRIVERS)/storef.o $(ADOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@      \
	$(DRIVERS)/prgtstadm.o $(DRIVERS)/storef.o $(ADOBJ) $(LIBS)

# combine timing and derivative check for tangent
../tsttlm: $(DRIVERS)/prgtsttlm.o $(DRIVERS)/storef.o $(TLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtsttlm.o $(DRIVERS)/storef.o $(TLOBJ) $(LIBS)

# combine timing and derivative check for tangent
../tstvtlm: $(DRIVERS)/prgtstvtlm.o $(DRIVERS)/storef.o $(VTLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstvtlm.o $(DRIVERS)/storef.o $(VTLOBJ) $(LIBS)

# combine timing and derivative check for gtb adjoint
../tstgtb: $(DRIVERS)/prgtstgtb.o $(DRIVERS)/storef.o $(GTBOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@      \
	$(DRIVERS)/prgtstgtb.o $(DRIVERS)/storef.o $(GTBOBJ) $(LIBS)

# combine timing and ASD check for progn adjoint
../tstspadm: $(DRIVERS)/prgasd.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) $(SPTLOBJ) $(SPADOBJ) $(PADOBJ) pfunc.o
	$(F90) $(DBL) $(LIBDIRS) -o $@ \
	$(DRIVERS)/prgasd.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) \
	$(SPTLOBJ) $(SPADOBJ)  $(PADOBJ) pfunc.o $(LIBS) 

# combine timing and derivative check for progn adjoint
../tstpadm: $(DRIVERS)/prgtstjaccolpad.o $(DRIVERS)/storef.o $(PADOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@      \
	$(DRIVERS)/prgtstjaccolpad.o $(DRIVERS)/storef.o $(PADOBJ) $(LIBS)

# combine timing and derivative check for progn tangent
../tstptlm: $(DRIVERS)/prgtstjaccolptl.o $(DRIVERS)/storef.o $(PTLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccolptl.o $(DRIVERS)/storef.o $(PTLOBJ) $(LIBS) 

# combine timing and derivative check for columns of hessian
../tsthesscol: $(DRIVERS)/prgtstjaccol.o $(DRIVERS)/storef.o $(HESOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccol.o $(DRIVERS)/storef.o $(HESOBJ) $(LIBS) 

# combine timing and derivative check for ift tangent
../tstitlm: $(DRIVERS)/prgtstjaccolitl.o $(DRIVERS)/storef.o $(ITLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccolitl.o  $(DRIVERS)/storef.o $(ITLOBJ) $(LIBS) 

# combine timing and derivative check for columns of hessian
../tstftlm: $(DRIVERS)/prgtstjaccolftl.o $(DRIVERS)/storef.o $(FTLOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgtstjaccolftl.o $(DRIVERS)/storef.o $(FTLOBJ) $(LIBS) 

#---------------------------------------------------------
# do optimisation
#---------------------------------------------------------

#  compute function and gradient
../fgrad: $(DRIVERS)/prggrad.o $(DRIVERS)/storef.o $(ADOBJ)
	$(F90) $(LIBDIRS) $(OPTLIBDIRS)	\
		-o $@ $(DRIVERS)/prggrad.o $(DRIVERS)/storef.o $(ADOBJ) $(LIBS)

# compute next control vector
../nextx: $(DRIVERS)/prgopti.o $(DRIVERS)/storef.o $(OBJ) $(MOBJ)
	$(F90) $(OPTLIBDIRS)	\
		-o $@ $(DRIVERS)/prgopti.o $(DRIVERS)/storef.o $(OBJ) $(MOBJ) $(OPTLIBS)

#  optimisation with gtb
../opti-gtb: $(DRIVERS)/prggtb.o $(OBJ) $(MOBJ) $(GTBOBJ) func-gtb.o
	$(F90) $(LIBDIRS) \
		-o $@ $(DRIVERS)/prggtb.o $(OBJ) $(MOBJ) $(GTBOBJ) func-gtb.o $(LIBS) -lfgtb

../dfp: $(DRIVERS)/storef.o $(ADOBJ) $(DRIVERS)/prgdfpmin.o 
	$(F90) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgdfpmin.o $(DRIVERS)/storef.o $(ADOBJ) $(LIBS) ../libnumrec.a

../dfptl: $(DRIVERS)/prgdfpmintl.o $(DRIVERS)/storef.o $(TLOBJ) 
	$(F90) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgdfpmintl.o $(TLOBJ) $(LIBS) ../libnumrec.a

# print current control vector and gradient
../post: $(DRIVERS)/prgpost.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) 
	$(F90) $(LIBDIRS)	\
		-o $@ $(DRIVERS)/prgpost.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) $(LIBS)

# do a linesearch between two control vectors
../line: $(DRIVERS)/line.o $(MOBJ) $(OBJ) 
	$(F90) $(LIBDIRS)	\
		-o $@ $(DRIVERS)/line.o $(MOBJ) $(OBJ) $(LIBS)

../dline: $(DRIVERS)/dline.o $(TLOBJ) 
	$(F90) $(LIBDIRS)	\
		-o $@ $(DRIVERS)/dline.o $(TLOBJ) $(LIBS)

#---------------------------------------------------------
# evaluate Hessian matrix:
#---------------------------------------------------------

# compute columns of the Hessian
../hesscol: $(DRIVERS)/prgjaccol.o $(DRIVERS)/storef.o $(HESOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgjaccol.o $(DRIVERS)/storef.o $(HESOBJ) $(LIBS) 

#---------------------------------------------------------
# rules for utility programms with bethy modules
#---------------------------------------------------------

# print current parameter vector
#../printp: mo_namelist.o mo_mapping.o $(UTIL)/printp.o 
#	$(F90) $(DBL) $(LIBDIRS) -o $@ \
#	$(UTIL)/printp.o mo_namelist.o mo_mapping.o 
../printp: $(MOBJ) $(OBJ)  $(UTIL)/printp.o 
	$(F90) $(DBL) $(LIBDIRS) -o $@ \
	$(MOBJ) $(OBJ) $(UTIL)/printp.o $(LIBS) 

# print current parameter vector old version (without params transformation)
../printparams: mo_namelist.o mo_constants.o $(UTIL)/printparams.o 
	$(F90) $(DBL) $(LIBDIRS) -o $@ \
	$(UTIL)/printparams.o mo_namelist.o mo_constants.o

# test program for parameter transformation 
../parmap: $(MOBJ) $(OBJ) $(UTIL)/parmap.o  
	$(F90) $(DBL) $(LIBDIRS) -o $@ \
	$(UTIL)/parmap.o $(MOBJ) $(OBJ) $(LIBS)

../invhess: mo_constants.o mo_grid.o mo_namelist.o mo_trafo.o $(UTIL)/invhess.o  
	$(F90) $(DBL) $(LIBDIRS) -L$(LAPACK) -L$(BLAS) -o $@	\
	$(UTIL)/invhess.o mo_constants.o mo_grid.o mo_namelist.o mo_trafo.o $(LIBS) -l$(LAPACKLIB) -l$(BLAS1LIB) 

../invhess_pseudo: mo_constants.o mo_namelist.o $(UTIL)/invhess_pseudo.o  
	$(F90) $(DBL) $(LIBDIRS) -L$(LAPACK) -L$(BLAS) -o $@	\
	$(UTIL)/invhess_pseudo.o mo_constants.o mo_namelist.o $(LIBS) -l$(LAPACKLIB) -l$(BLAS1LIB) 

../jachess: mo_netcdf.o mo_constants.o mo_grid.o mo_namelist.o $(UTIL)/jac_hess.o  
	$(F90) $(DBL) $(LIBDIRS) -L$(LAPACK) -L$(BLAS) -o $@	\
	$(UTIL)/jac_hess.o mo_netcdf.o mo_constants.o mo_grid.o mo_namelist.o $(LIBS) -l$(LAPACKLIB) -l$(BLASLIB) 


#---------------------------------------------------------
# do Prognostic Step
#---------------------------------------------------------

# adapt pfunc.f90 and ifunc
mo_grid.f90: mo_grid.g90 
	$(CPP) $(CPPFLAGS) $< >$@

pfunc.f90: pfunc.g90 
	$(CPP) $(CPPFLAGS) $< >$@

ifunc.f90: ifunc.g90 
	$(CPP) $(CPPFLAGS) $< >$@

ffunc.f90: ffunc.g90 
	$(CPP) $(CPPFLAGS) $< >$@

model.f90: model.g90 
	$(CPP) $(CPPFLAGS) $< >$@

bethy.f90: bethy.g90 
	$(CPP) $(CPPFLAGS) $< >$@

costf.f90: costf.g90 
	$(CPP) $(CPPFLAGS) $< >$@

diagnostics.f90: diagnostics.g90 
	$(CPP) $(CPPFLAGS) $< >$@

mo_diagnostics.f90: mo_diagnostics.g90 
	$(CPP) $(CPPFLAGS) $< >$@

mo_climate.f90: mo_climate.g90 
	$(CPP) $(CPPFLAGS) $< >$@

cbalance.f90: cbalance.g90 
	$(CPP) $(CPPFLAGS) $< >$@

# evaluate scalar valued function/prognostic?? twice
../ffunc2: $(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(FOBJ) ffunc.o 
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(FOBJ) ffunc.o $(LIBS)

# evaluate scalar valued function/prognostic?? twice
../pfunc2: $(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) pfunc.o
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) pfunc.o $(LIBS)

# evaluate scalar valued function/prognostic?? twice
../func-gtb2: $(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) func-gtb.o
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) func-gtb.o $(LIBS)

# evaluate scalar valued function/prognostic?? twice
../ifunc2: $(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) $(ADOBJ) ifunc.o
	$(F90) $(DBL) $(LIBDIRS) -o $@	\
	$(DRIVERS)/prgfunc2.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) $(ADOBJ) ifunc.o $(LIBS)

# compute Jacobian for prognostic step (reverse)
../prev: $(DRIVERS)/prgjacrevpad.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) $(PADOBJ) pfunc.o
	$(F90) $(DBL) $(LIBDIRS) -o $@ 	\
	$(DRIVERS)/prgjacrevpad.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) $(PADOBJ)  pfunc.o $(LIBS) 

# compute Jacobian for prognostic step (forward)
../pfwd: $(DRIVERS)/prgjacfwdptl.o $(DRIVERS)/storef.o $(PTLOBJ)
	$(F90) $(DBL) $(LIBDIRS) -o $@ 	\
	$(DRIVERS)/prgjacfwdptl.o $(DRIVERS)/storef.o $(PTLOBJ) $(LIBS) 

# compute additional second derivative term for background sensitivity
../ifwd: $(DRIVERS)/prgjaccolitl.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) ifunc.o $(ADOBJ) $(ITLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@ 	\
	$(DRIVERS)/prgjaccolitl.o $(DRIVERS)/storef.o $(MOBJ) $(OBJ) ifunc.o $(ADOBJ) $(ITLOBJ) $(LIBS) 

# compute additional second derivative term for background sensitivity
../ffwd: $(DRIVERS)/prgjaccolftl.o $(DRIVERS)/storef.o $(FTLOBJ) 
	$(F90) $(DBL) $(LIBDIRS) -o $@ 	\
	$(DRIVERS)/prgjaccolftl.o $(DRIVERS)/storef.o $(FTLOBJ) $(LIBS) 

#---------------------------------------------------------
# rules for removing
#---------------------------------------------------------
# remove .o, .mod rsp .d files
clearo:
	@$(RM) *.o *.d *.mod

# remove intermediate files, warmstart still possible
clean:	clearo
# remove generated files, only coldstart possible
scratch: clean
	@$(RM) *_ad.flow *_ftl.flow  taf_command taf_output 
	@$(RM) mo_grid.f90 pfunc.f90 ifunc.f90 ffunc.f90 model.f90 bethy.f90 costf.f90 diagnostics.f90 mo_diagnostics.f90 mo_climate.f90 cbalance.f90 
	@$(RM) *~ \#*\# work.pc work.pcl

rmall: scratch
	@$(RM) *_tl.f90 *_vtl.f90 *_ad.f90 *_gtb.f90 *_pad.f90 *_ptl.f90 *_hes.f90
	@$(RM) *_spad.f90 *_sptl.f90 *_itl.f90 *_iad.f90 tafad.f90
	@$(RM) *_ad.f *_ftl.f *_pad.f *_ptl.f *_spad.f *_sptl.f *_hes.f 
	@$(RM) taf*.log

#
#-------------------------------------------------------
# DEPENDENCIES SRC-FILES
#---------------------------------------------------------
imbethy.o: mo_imbethy.o
adbethy.o : mo_constants.o mo_carparams.o mo_carbon.o mo_grid.o mo_io.o mo_climate.o mo_diagnostics.o mo_vegetation.o costf.o model.o mo_surface.o mo_carvar.o
cbalance.o : mo_constants.o mo_namelist.o  mo_grid.o mo_vegetation.o mo_beta.o mo_diagnostics.o mo_calendar.o mo_climate.o
budget.o : mo_grid.o
climate_allocate.o : mo_climate.o mo_constants.o
climate_deallocate.o : mo_climate.o mo_constants.o
climatereadday.o : mo_constants.o mo_io.o mo_climate.o mo_grid.o
#climinday.o : mo_constants.o  mo_carparams.o mo_grid.o mo_io.o mo_climate.o
climinday.o : mo_constants.o  mo_carparams.o mo_grid.o mo_namelist.o mo_netcdf.o mo_climate.o  mo_io.o 
climinsubday.o : mo_constants.o mo_io.o mo_climate.o
daytemp.o : mo_constants.o 
diagnostics.o : mo_constants.o mo_diagnostics.o mo_grid.o mo_hydro.o mo_climate.o mo_pheno.o
func.o : mo_namelist.o
init_site.o : mo_constants.o mo_grid.o mo_netcdf.o mo_carparams.o mo_vegetation.o costf.o
initialize.o : mo_constants.o mo_grid.o mo_netcdf.o mo_carparams.o mo_vegetation.o
initmod.o : costf.o mo_trafo.o
#model.o : mo_constants.o mo_vegetation.o mo_climate.o mo_io.o costf.o mo_surface.o  surface.o mo_carbon.o climinday.o climinsubday.o bgr.o tm.o bethy2tm.o mo_taf.o mo_mapping.o mo_namelist.o mo_beta.o mo_ctrl.o mo_carparams.o mo_carvar.o mo_diagnostics.o mo_grid.o
bethy.o : mo_constants.o mo_vegetation.o mo_climate.o mo_io.o costf.o mo_surface.o  surface.o mo_carvar.o mo_taf.o mo_carbon.o bgr.o tm.o bethy2tm.o mo_taf.o mo_mapping.o mo_namelist.o mo_beta.o mo_ctrl.o mo_carparams.o mo_carvar.o mo_diagnostics.o mo_grid.o mo_pheno.o mo_hydro.o mo_calendar.o mo_helper.o
model.o : mo_config.o mo_trafo.o
#numbmod.o : mo_constants.o mo_vegetation.o mo_climate.o mo_io.o costf.o mo_surface.o  surface.o climinday.o climinsubday.o bgr.o tm.o bethy2tm.o mo_mapping.o mo_namelist.o mo_beta.o
numbmod.o : mo_constants.o mo_vegetation.o mo_climate.o mo_io.o mo_ctrl.o mo_beta.o constants_pjr.o tm.o mo_mapping.o costf.o mo_surface.o  surface.o bgr.o tm.o bethy2tm.o mo_mapping.o mo_namelist.o mo_beta.o mo_prog.o mo_config.o
pfunc.o : mo_constants.o constants_pjr.o tm.o mo_namelist.o mo_grid.o tm.o
postadm.o : mo_diagnostics.o mo_io.o costf.o tm.o mo_namelist.o constants_pjr.o mo_trafo.o
postmod.o : mo_diagnostics.o mo_io.o costf.o tm.o mo_namelist.o constants_pjr.o mo_prog.o mo_trafo.o
rad.o : mo_constants.o
surface.o : mo_constants.o mo_carparams.o mo_surface.o mo_grid.o
veg_allocate.o : mo_constants.o mo_vegetation.o  mo_surface.o mo_carvar.o mo_carparams.o mo_beta.o mo_pheno.o mo_hydro.o
veg_deallocate.o : mo_constants.o mo_vegetation.o  mo_surface.o mo_carvar.o mo_carparams.o mo_beta.o mo_pheno.o mo_hydro.o
wgen.o : mo_constants.o

#---------------------------------------------------------
# DEPENDENCIES MODULES
#---------------------------------------------------------
#HEW-ADD:
mo_namelist.o :constants_pjr.o  mo_constants.o mo_grid.o
#mo_namelist.o: constants_pjr.o 
bethy2tm.o : constants_pjr.o mo_constants.o mo_netcdf.o
bgr.o : constants_pjr.o mo_netcdf.o
costf.o : mo_constants.o mo_vegetation.o mo_carparams.o constants_pjr.o mo_mapping.o mo_climate.o mo_trafo.o
mo_beta.o : mo_constants.o
mo_carbon.o : mo_constants.o mo_carparams.o mo_carvar.o mo_climate.o mo_taf.o mo_grid.o
mo_carparams.o : mo_constants.o
mo_climate.o: mo_constants.o mo_constants.o mo_netcdf.o mo_calendar.o mo_carparams.o
mo_ctrl.o : mo_constants.o
mo_grid.o : mo_constants.o
mo_io.o : mo_constants.o mo_netcdf.o mo_grid.o mo_climate.o mo_namelist.o
mo_mapping.o : mo_constants.o mo_vegetation.o mo_beta.o mo_ctrl.o tm.o mo_grid.o mo_pheno.o
tm.o : constants_pjr.o mo_netcdf.o
mo_pheno.o : mo_constants.o mo_climate.o mo_namelist.o mo_grid.o mo_netcdf.o mo_diagnostics.o mo_helper.o
mo_hydro.o : mo_constants.o mo_climate.o mo_namelist.o mo_grid.o mo_netcdf.o mo_helper.o
mo_diagnostics.o : mo_constants.o

# DEPENDENCIES SRC-FILES
#pfunc.f90 : pfunc.F90

