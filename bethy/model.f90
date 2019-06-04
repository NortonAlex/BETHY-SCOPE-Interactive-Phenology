
SUBROUTINE model( nvar, x, fc )

  USE mo_constants
  USE mo_carparams
  USE mo_namelist, ONLY : nrun, year0, year1, yearin0, yearin1, nspin, dayint, &
       & year0_site, year1_site, site_file, grid_file, optpftl, optpftg, optbsmg, optbsml
  USE tm ! transport moduleOD
  USE bgr ! background fluxes
  USE costf
  USE bethy2tm
  USE mo_mapping
  USE mo_calendar 
  USE mo_climate
  USE mo_grid, ONLY : ng, vp, sp, gridp, nrs, frac, vg, sumfrac
  USE mo_hydro
  USE mo_pheno
  USE mo_config
  USE mo_prog
  USE mo_trafo
  USE mo_taf
  USE mo_diagnostics, ONLY : rlai
  USE mo_vegetation, ONLY : pft

  IMPLICIT NONE

! Arguments ..
  INTEGER              , INTENT(in)  :: nvar
  REAL, DIMENSION(nvar), INTENT(inout)  :: x
  REAL                 , INTENT(out) :: fc

! Local ..
  real, parameter :: mult = 1.e8
  REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: faparl
  REAL, ALLOCATABLE, DIMENSION(:, :, :)    :: faparg
  REAL :: fc_p, fc_c, fc_faparl, fc_faparg, fc_flux, fc_wmax, fc_pft!, fc_straf

  REAL, ALLOCATABLE, DIMENSION(:) :: flux_magnitudes 
  REAL, DIMENSION(nvar) :: par
  INTEGER :: i_x, i_flux ! index variables for direct flux part
  INTEGER :: nchk, nxp
  INTEGER :: k, j, np, i, y, m
  INTEGER :: scale
  REAL :: h1, pnlty

!$taf init top_tape = static, 1
! maxscale = 2
!$taf init scale_tape = static, 2
!$taf init day_tape = static,  maxkeyday
!$taf init dayint_tape = static,  maxkeydayint
!$taf init diurnal_tape = static,  maxkeydiurnal
!$taf init carbon_tape = static,  maxkeydiurnal
!$taf init carbon_tape_nl = static,  nl*maxkeydiurnal
!FastOpt !$taf init cbalance_tape = static,  maxdpy*(aspin+maxyears)*2 

  fc = 0.
  nchk = 1
  nxp = size(xpf)

! MAS-RM-100407 removed xsf as scaling flag
!  x(1:nxp) = x(1:nxp) / ABS(xsf)

! parameter transformation from x (unitless and scaled) to p
  par = x2p(x, x0, xpf, xpfa, xpfb, px0, p0su, a, nvar)

  IF (n_stats>0) THEN
     PRINT*
     PRINT*,'BETHY run on global grid for flask stations.'
     scale = 1
! .. initialize model
     CALL config_global (grid_file) 
     
     CALL init_config (ng)

! allocate the model variable fields     
     CALL model_allocate(outt, nrun, ng, vp)
     CALL pheno_allocate(vp)
     CALL hydro_allocate(ng,vp,scale)

! distributes physical parameter values on grid
     CALL mapping( par, nvar, nxp, scale) 

! generates calender and time management
     CALL calendar(year0,year1,yearin0,yearin1,nspin,dayint)

     call init_faparg(ng, nrun, outt)

! allocates variables for climate input data (dprecip, dtmin, dtmax, dswdown)
     CALL climate_allocate(ng, vp, sdays)

! reads in climate input data
     CALL get_global_climate
     
     ALLOCATE(faparg(nrun,outt,ng))

!------------------------------------------------------------------
! first bethy call, global grid
!------------------------------------------------------------------
     CALL bethy (nchk, dayint, ng, vp, nrun, outt, scale, netflux, faparg)

!------------------------------------------------------------------
! do interpolation from bethy to TM grid
!------------------------------------------------------------------
!$taf store netflux = top_tape, rec = 1 
! netflux only needed for its 'size', should be improved
     CALL run_bethy2tm(f_tm) 

! pjr 04/04/28 add parts for direct fluxes
     ALLOCATE( flux_magnitudes( COUNT( xpf == direct_flux)))
! collect all the direct_flux parameters from the param vector and scale them by their unc
     i_flux = 1
     DO i_x = 1, size(xpf)
        IF( xpf( i_x) == direct_flux) THEN
           flux_magnitudes( i_flux) = x( i_x) * p0su(  i_x) 
           i_flux = i_flux + 1
        ENDIF
     END DO

!------------------------------------------------------------------
! add background fluxes
!------------------------------------------------------------------
!$taf store f_tm = top_tape, rec = 1 
! f_tm only needed for its 'size', should be improved
     CALL run_bgr ( f_tm , flux_magnitudes)
     DEALLOCATE( flux_magnitudes)

!------------------------------------------------------------------
! now transport it
!------------------------------------------------------------------
     CALL run_tm (glob_offset, f_tm, conc) 

!------------------------------------------------------------------
! deallocate climate variables
!------------------------------------------------------------------
     CALL climate_deallocate (ng, vp, sdays)

!------------------------------------------------------------------
! deallocate calendar variables
!------------------------------------------------------------------
     CALL deallocate_cal  

!      if (optpftg) then
!!         fc = fc + cost_pft(x(nxp+1:nxp+vp),frac_u)
!         call cost_pft(fc_pft,(x(nxp+1:nxp+vp),frac_u)
!         fc = fc + fc_pft
!         if (optbsmg) then
!            call cost_wmax(fc_wmax,x(nxp+vp+1:nxp+vp+vp),pasmmax_u)
!!            fc = fc + cost_wmax(x(nxp+vp+1:nxp+vp+vp),pasmmax_u)
!            fc = fc + fc_wmax
!!         else
!!            fc = fc
!         endif
!      else
!         if (optbsmg) then
!            call cost_wmax(cf_wmax,x(nxp+vp+1:nxp+vp+vp),pasmmax_u)
!            fc = fc + fc_wmax
!!            fc = fc + cost_wmax(x(nxp+1:nxp+vp),pasmmax_u)           
!!         else
!!            fc = fc
!         endif
!      endif

!     fc_straf= 0.
!     do i=1,vp
!        j = gridp (i)
!        if (pft(i)==1 .and. frac(i)>0.5) then
!           do y=1,nrun
!              do m=1,12
!                 if (rlai(y,m,j) < 3. ) then
!                    fc_straf = fc_straf + (rlai(y,m,j) - 3.)**2/(0.5)**2
!                 endif
!              enddo
!           enddo
!        endif
!     enddo
!     print*,'cost_straf: ',fc_straf
!     fc = fc + fc_straf
!------------------------------------------------------------------
! deallocate init variables
!------------------------------------------------------------------
!     CALL config_deallocate
     CALL config_deallocate(ng,vp) ! note: use these if passing dimension to TAF
     CALL model_deallocate(outt, nrun, ng, vp)
     CALL hydro_deallocate(ng, vp)
     CALL pheno_deallocate(vp)

!------------------------------------------------------------------
! .. concentration and FAPAR costfunction global
!------------------------------------------------------------------
!$taf store conc = top_tape, rec = 1
     call cost_c(fc_c,conc)
     call cost_faparg(fc_faparg,faparg,nrun,outt,ng)
     fc = fc + fc_c + fc_faparg

     DEALLOCATE(faparg,faparg_obs,faparg_unc)

  ENDIF

!==================================================================
!==================================================================
! until here global run
!==================================================================
!==================================================================

  IF (n_sites>0) THEN

     PRINT*
     PRINT*,' BETHY run at flux sites.'
     scale = 2

! initialise site specific values
     CALL config_sites (site_file, n_sites) 
     
     CALL init_config (n_sites) 

     CALL model_allocate(366, nrs, n_sites, sp)
     CALL pheno_allocate(sp)
     CALL hydro_allocate(n_sites,sp,scale)

! distributes site specific physical parameter values
     CALL mapping( par, nvar, nxp, scale) 

! generates calender and time management
     CALL calendar(year0_site,year1_site,yearin0,yearin1,nspin,1)

     CALL init_faparl(n_sites, nrs, 366)

! allocates variables for climate input data (dprecip, dtmin, dtmax, dswdown)
     CALL climate_allocate (n_sites, sp, sdays)

! reads in climate input data
     CALL get_local_climate

     ALLOCATE(faparl(nrs,366,n_sites))

!------------------------------------------------------------------
! add term to costfunction if PFT frac <0 or sum(PFT frac) >1, disallow negative fractions
!------------------------------------------------------------------
!$taf loop = parallel
!      np = 1
!      DO k=1,n_sites
!         DO j=1,vg(k)
!            IF (frac(np) < 0.) THEN
!               fc = fc + mult*frac(np)**2
!               print*,'Model: negative fraction of ',frac(np),' at sub-gridcell ',np, ' with penalty ',mult*frac(np)**2
!            ENDIF
!            sumfrac(k)=sumfrac(k)+frac(np)
!            np = np + 1
!         ENDDO
!      ENDDO

!      DO k=1,n_sites
!         IF (sumfrac(k) > 1.) THEN
!            print*,'Model: PFT sum of ',sumfrac(k),' with penalty ', 1.e+3*(1-sumfrac(k))**4
!            fc = fc + 1.e+3*(1-sumfrac(k))**4
!         ENDIF
!      ENDDO

! some diagnostic output for optimsing fractional cover frac
!------------------------------------------------------------------
! second bethy call, 1 sites
!------------------------------------------------------------------
     CALL bethy (nchk, dayint, n_sites, sp, nrs, 366, scale, netflux, faparl)

!------------------------------------------------------------------
! handling of flux observations: add cost of fluxes
!------------------------------------------------------------------
     call cost_flux(fc_flux, netflux, nrs, outt, n_sites)
     fc = fc + fc_flux

!------------------------------------------------------------------
! deallocate climate variables
!------------------------------------------------------------------
     CALL climate_deallocate(n_sites, sp, sdays)

!------------------------------------------------------------------
! deallocate calendar variables
!------------------------------------------------------------------
     CALL deallocate_cal  

!------------------------------------------------------------------
! deallocate init variables
!------------------------------------------------------------------
!     CALL config_deallocate
     CALL config_deallocate(n_sites,sp) ! note: use this if passing dimension to TAF
     CALL model_deallocate(366, nrs, n_sites, sp)
     CALL hydro_deallocate(n_sites,sp)
     CALL pheno_deallocate(sp)

!------------------------------------------------------------------
! .. fapar costfunction at site level: add cost of FAPAR
!------------------------------------------------------------------

!     fc = fc + cost_faparl(faparl,nrs,366,n_sites) 
     call cost_faparl(fc_faparl,faparl,nrs,366,n_sites) 
     fc = fc + fc_faparl 
     DEALLOCATE(faparl,faparl_obs,faparl_unc)

  ENDIF

!------------------------------------------------------------------
! .. parameter contribution to costfunction
!------------------------------------------------------------------
  call cost_p(fc_p,x)
  fc = fc + fc_p

END SUBROUTINE model

!------------------------------------------------------------------
! do all the model allocations (except hydrology & phenology)
!------------------------------------------------------------------
SUBROUTINE model_allocate(outt, nrs, ng, vp)
  USE bethy2tm
  USE mo_diagnostics
  IMPLICIT NONE
  INTEGER, INTENT(in) :: outt,nrs, ng, vp

  ALLOCATE(netflux(nrs,outt,ng))
  CALL veg_allocate(vp)
  CALL diagnostics_allocate(nrs,outt, ng, vp)
END SUBROUTINE model_allocate


!------------------------------------------------------------------
! do all the model deallocations (except hydrology & phenology)
!------------------------------------------------------------------
SUBROUTINE model_deallocate(outt, nrs, ng, vp)
  USE bethy2tm
  USE mo_diagnostics
  IMPLICIT NONE
  INTEGER, INTENT(in) :: outt,nrs, ng, vp

!$taf next required = nrs,outt,ng
  DEALLOCATE(netflux)
  CALL veg_deallocate(vp)
  CALL diagnostics_deallocate(nrs, outt, ng, vp)
END SUBROUTINE model_deallocate


