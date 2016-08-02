!$taf module mo_io SUBROUTINE savefnc INPUT =1,2,3,4,5,6,7
!$TAF module mo_io subroutine savefnc output=
!$TAF module mo_io subroutine savefnc adname=ignore
!*********************************************************
!*********************************************************
MODULE mo_io

IMPLICIT NONE
!  mtmp monthly mean temperature (degree C)


  REAL, ALLOCATABLE, DIMENSION (:,:) :: mtran, mpp, mrt, mtmp
  REAL, ALLOCATABLE, DIMENSION (:,:) :: mvp, mwind
! REAL, ALLOCATABLE, DIMENSION (:,:,:) :: mlai, mpasm

CONTAINS

!HEW : initopt removed, now only read in from file control !!!
!*********************************************************
!*  SUBROUTINE initopt
!*  reads option file, initializes and allocates memory
!*********************************************************

!HEW-DEL SUBROUTINE initopt

!HEW-DEL ...


!HEW-DEL END SUBROUTINE initopt


!***********************************************************
!*   SUBROUTINE setio
!*   opens input/output files 
!***********************************************************

SUBROUTINE setio(outdir, scale)

! .. Use Statements ..
  USE mo_constants

  IMPLICIT NONE

! .. Arguments
  CHARACTER(len=80),INTENT(in) :: outdir
  INTEGER, INTENT(in) :: scale

! .. Local Scalars ..
  CHARACTER (len=4) :: cyin
  CHARACTER (len=80) :: ofname
  INTEGER :: uu

!cccc OPEN OUTPUT FILES:

!     error messages 
  uu=8
  ofname='error.dat'
  OPEN (uu,file=ofname,status='unknown',form='formatted')

!     debugging messages 
  uu=9
  ofname='debugging.dat'
  OPEN (uu,file=ofname,status='unknown',form='formatted')

  IF (scale==1) THEN

!     GPP [gC / m**2 month]
     uu=140
     ofname=TRIM(outdir) // 'gpp_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     NPP [gC / m**2 month]
     uu=141
     ofname=TRIM(outdir)//'npp_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     NEP [gC / m**2 month]
     uu=142
     ofname=TRIM(outdir)//'rnep_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     slow soil pool respiration [gC / m**2 month]
     uu=143
     ofname=TRIM(outdir)//'ress_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     fast soil pool respiration [gC / m**2 month]
     uu=144
     ofname=TRIM(outdir)//'resf_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     growth respiration [gC / m**2 month]
     uu=145
     ofname=TRIM(outdir)//'rgrw_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     maintenance respiration [gC / m**2 month]
     uu=146
     ofname=TRIM(outdir)//'rmnt_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!MAS-ADD-BEG-040220 output of npp per PFT 
!     npp per PFT  [gC / m**2 month]
     uu=147
     ofname=TRIM(outdir)//'npp1_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=148
     ofname=TRIM(outdir)//'npp2_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=149
     ofname=TRIM(outdir)//'npp3_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
!MAS-ADD-END-040220

!WOK-ADD-070725  various diagnostics from hydrology and phenology
!     soil reflectance
     uu=150
     ofname=TRIM(outdir)//'rhos_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     plant-available soil moisture [mm]
     uu=151
     ofname=TRIM(outdir)//'pasm_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     leaf-area index
     uu=152
     ofname=TRIM(outdir)//'lai_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     absorbed fraction of photosynthetically active radation
     uu=153
     ofname=TRIM(outdir)//'fpar_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     runoff [mm/month]
     uu=154
     ofname=TRIM(outdir)//'runoff_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     soil evaporation [mm/month]
     uu=155
     ofname=TRIM(outdir)//'sevp_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     snow melt [mm/month]
     uu=156
     ofname=TRIM(outdir)//'snmelt_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     actual evapotranspiration [mm/month]
     uu=157
     ofname=TRIM(outdir)//'aet_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     potential evapotranspiration [mm/month]
     uu=158
     ofname=TRIM(outdir)//'pet_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

     uu=163
     ofname=TRIM(outdir)//'fpar1_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=164
     ofname=TRIM(outdir)//'fpar2_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=165
     ofname=TRIM(outdir)//'fpar3_global.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

  ELSEIF (scale==2) THEN

!     GPP [gC / m**2 month]
     uu=240
     ofname=TRIM(outdir) // 'gpp_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     NPP [gC / m**2 month]
     uu=241
     ofname=TRIM(outdir)//'npp_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     NEP [gC / m**2 month]
     uu=242
     ofname=TRIM(outdir)//'rnep_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     slow soil pool respiration [gC / m**2 month]
     uu=243
     ofname=TRIM(outdir)//'ress_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     fast soil pool respiration [gC / m**2 month]
     uu=244
     ofname=TRIM(outdir)//'resf_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     growth respiration [gC / m**2 month]
     uu=245
     ofname=TRIM(outdir)//'rgrw_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     maintenance respiration [gC / m**2 month]
     uu=246
     ofname=TRIM(outdir)//'rmnt_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!MAS-ADD-BEG-040220 output of npp per PFT 
!     npp per PFT  [gC / m**2 month]
     uu=247
     ofname=TRIM(outdir)//'npp1_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=248
     ofname=TRIM(outdir)//'npp2_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=249
     ofname=TRIM(outdir)//'npp3_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=260
     ofname=TRIM(outdir)//'gpp1_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=261
     ofname=TRIM(outdir)//'gpp2_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=262
     ofname=TRIM(outdir)//'gpp3_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=263
     ofname=TRIM(outdir)//'fpar1_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=264
     ofname=TRIM(outdir)//'fpar2_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=265
     ofname=TRIM(outdir)//'fpar3_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=266
     ofname=TRIM(outdir)//'nep1_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=267
     ofname=TRIM(outdir)//'nep2_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
     uu=268
     ofname=TRIM(outdir)//'nep3_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
!MAS-ADD-END-040220

!WOK-ADD-070725  various diagnostics from hydrology and phenology
!     soil reflectance
     uu=250
     ofname=TRIM(outdir)//'rhos_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     plant-available soil moisture [mm]
     uu=251
     ofname=TRIM(outdir)//'pasm_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     leaf-area index
     uu=252
     ofname=TRIM(outdir)//'lai_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     absorbed fraction of photosynthetically active radation
     uu=253
     ofname=TRIM(outdir)//'fpar_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     runoff [mm/month]
     uu=254
     ofname=TRIM(outdir)//'runoff_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     soil evaporation [mm/month]
     uu=255
     ofname=TRIM(outdir)//'sevp_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     snow melt [mm/month]
     uu=256
     ofname=TRIM(outdir)//'snmelt_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     actual evapotranspiration [mm/month]
     uu=257
     ofname=TRIM(outdir)//'aet_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')

!     potential evapotranspiration [mm/month]
     uu=258
     ofname=TRIM(outdir)//'pet_site.dat'
     OPEN (uu,file=ofname,status='unknown',form='formatted')
  ENDIF
END SUBROUTINE setio

SUBROUTINE closeio
  close (8)
  close (9)
END SUBROUTINE closeio

!***********************************************************
!*   SUBROUTINE diagout
!*   writes diagnostic output files
!***********************************************************

!MAS-ADD-040220 nppp
!WOK-ADD-070725  various diagnostics from hydrology and phenology
SUBROUTINE diagout (ng,vp,scale,outint)

! .. Use Statements ..
  USE mo_constants, ONLY: outt, nv, rdays
  USE mo_grid, ONLY: lon, lat, nrs
!HEW-ADD: added nrun, ...
  USE mo_namelist, ONLY: outdir, nrun, year0, year0_site
  USE mo_diagnostics, ONLY: rgpp,rnpp,rnep,raet,rpet,ress,resf,rgrw,rmnt,nppp, &
                       &   rrhos,rpasm,rlai,rfpar,rrunoff,rsevp,rsnmelt,rnppp, &
                       &   rgppp, rfparp, rnepp, rfluo, rgppfluo, rrad, rpar, &
                       &   PAR_scope, PAR_scope_cab, rfluo_diurnal, rgppfluo_diurnal

  IMPLICIT NONE


! .. Arguments ..
  INTEGER, INTENT(in) :: ng, vp, scale, outint

!!MS$  REAL, DIMENSION(0:nrun,outt,ng), INTENT(inout) :: rgpp,rnpp,rnep,raet
!!MS$  REAL, DIMENSION(0:nrun,outt,ng), INTENT(inout) :: rpet,ress,resf,rgrw,rmnt
!!MS$!MAS-ADD-040220 nppp
!!MS$  REAL, DIMENSION(0:nrun,outt,ng,nv), INTENT(inout) :: nppp
!!MS$  REAL, DIMENSION(0:nrun,outt,ng), INTENT(inout) :: rrhos, rpasm, rlai, rfpar
!!MS$  REAL, DIMENSION(0:nrun,outt,ng), INTENT(inout) :: rrunoff, rsevp, rsnmelt

! .. Local Scalars ..
  INTEGER :: sp,i,fm, y, v, p
  REAL, DIMENSION(0:nrs,outt,ng,nv) :: rneppm
  sp=nrun*outt

  IF (scale==1) THEN
!     CALL savefield(140,'f8.2',rgpp,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(141,'f8.2',rnpp,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(142,'f8.2',rnep,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(143,'f8.2',ress,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(144,'f8.2',resf,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(145,'f8.2',rgrw,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(146,'f8.2',rmnt,lat,lon,nrun,outt,ng,year0)
!MAS-ADD-BEG-040220 nppp
!     CALL savefield(147,'f8.2',nppp(:,:,:,1),lat,lon,nrun,outt,ng,year0)
!     CALL savefield(148,'f8.2',nppp(:,:,:,2),lat,lon,nrun,outt,ng,year0)
!     CALL savefield(149,'f8.2',nppp(:,:,:,3),lat,lon,nrun,outt,ng,year0)
!MAS-ADD-END-040220 nppp
!     CALL savefield(150,'f6.3',rrhos,lat,lon,nrun,outt,ng,year0)
     CALL savefield(151,'f8.1',rpasm,lat,lon,nrun,outt,ng,year0)
     CALL savefield(152,'f6.2',rlai,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(153,'f6.3',rfpar,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(154,'f8.2',rrunoff,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(155,'f8.2',rsevp,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(156,'f8.2',rsnmelt,lat,lon,nrun,outt,ng,year0)
!WOK-CHG-071031 additional output for actual and potential evapotranspiration
!     CALL savefield(157,'f8.2',raet,lat,lon,nrun,outt,ng,year0)
!     CALL savefield(158,'f8.2',rpet,lat,lon,nrun,outt,ng,year0)

!     CALL savefield(163,'f6.3',rfparp(:,:,:,1),lat,lon,nrun,outt,ng,year0)
!     CALL savefield(164,'f6.3',rfparp(:,:,:,2),lat,lon,nrun,outt,ng,year0)
!     CALL savefield(165,'f6.3',rfparp(:,:,:,3),lat,lon,nrun,outt,ng,year0)
! and now netcdf gridded versions
     CALL savefnc(TRIM(outdir)//'rgpp.nc', rgpp,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rnpp.nc', rnpp,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rnep.nc', rnep,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'ress.nc', ress,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'resf.nc', resf,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rgrw.nc', rgrw,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rmnt.nc', rmnt,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rpasm.nc', rpasm,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rlai.nc', rlai,lat,lon,nrun,outt,ng,sp)

!ANorton. Fluorescence output fields
     CALL savefnc(TRIM(outdir)//'rfluo.nc', rfluo,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rgppfluo.nc', rgppfluo,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rrad.nc', rrad,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'rpar.nc', rpar,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'par_scope.nc', PAR_scope,lat,lon,nrun,outt,ng,sp)
     CALL savefnc(TRIM(outdir)//'par_scope_cab.nc', PAR_scope_cab,lat,lon,nrun,outt,ng,sp)
     CALL savefnc_diurnal(TRIM(outdir)//'rfluo_diurnal.nc', rfluo_diurnal,nrun,365,24,vp,sp)
     CALL savefnc_diurnal(TRIM(outdir)//'rgppfluo_diurnal.nc', rgppfluo_diurnal,nrun,365,24,vp,sp)


  ELSEIF (scale==2) THEN ! SITE SCALE !

     rneppm = 0.

     ! sum up daily nep values to monthly
!     print*, 'diagout: evaluating rnepp with these dimensions:'
!  do i = 1,4
!     print*, 'bounds of rnepp component ',i,' = ',lbound(rnepp,i),ubound(rnepp,i)
!  enddo
     do y=1,nrs
        do v=1,ng
           DO p=1,nv
              fm = 1
              do i=1,12
!                 print*, 'rnepp(',y,',',fm,':',fm+rdays(i)-1,',',v,',',p,')'
                 rneppm(y,i,v,p) = sum(rnepp(y,fm:fm+rdays(i)-1,v,p))
                 fm = fm + rdays(i)
              enddo
           enddo
        enddo
     enddo

     CALL savefield(240,'f8.3',rgpp,lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(241,'f8.3',rnpp,lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(242,'f8.3',rnep,lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(243,'f8.3',ress,lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(244,'f8.3',resf,lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(245,'f8.2',rgrw,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(246,'f8.2',rmnt,lat,lon,nrs,outt,ng,year0_site)
!MAS-ADD-BEG-040220 nppp
     CALL savefield(247,'f8.2',rnppp(:,:,:,1),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(248,'f8.2',rnppp(:,:,:,2),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(249,'f8.2',rnppp(:,:,:,3),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(260,'f8.2',rgppp(:,:,:,1),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(261,'f8.2',rgppp(:,:,:,2),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(262,'f8.2',rgppp(:,:,:,3),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(263,'f6.3',rfparp(:,:,:,1),lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(264,'f6.3',rfparp(:,:,:,2),lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(265,'f6.3',rfparp(:,:,:,3),lat,lon,nrs,outint,ng,year0_site)
     CALL savefield(266,'f8.2',rneppm(:,:,:,1),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(267,'f8.2',rneppm(:,:,:,2),lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(268,'f8.2',rneppm(:,:,:,3),lat,lon,nrs,outt,ng,year0_site)
!MAS-ADD-END-040220 nppp
     CALL savefield(250,'f6.3',rrhos,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(251,'f8.1',rpasm,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(252,'f6.2',rlai,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(253,'f6.3',rfpar,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(254,'f8.2',rrunoff,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(255,'f8.2',rsevp,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(256,'f8.2',rsnmelt,lat,lon,nrs,outt,ng,year0_site)
!WOK-CHG-071031 additional output for actual and potential evapotranspiration
     CALL savefield(257,'f8.2',raet,lat,lon,nrs,outt,ng,year0_site)
     CALL savefield(258,'f8.2',rpet,lat,lon,nrs,outt,ng,year0_site)
  ENDIF


  RETURN

END SUBROUTINE diagout


SUBROUTINE savefnc( filename, field,lat,lon,nm,ot,ng,binnm)
  ! pjr april 2002 for saving as netcdf output gridded
  USE mo_netcdf
  USE mo_namelist, ONLY : year0  ! avoiding name conflict
  USE mo_constants, ONLY: rdays
  USE mo_grid, ONLY: nlon, nlat

  IMPLICIT NONE

  CHARACTER(len=*) :: filename
! .. Scalar Arguments ..
  INTEGER, INTENT (IN) :: ot, ng, nm,  binnm

! .. Array Arguments ..
  REAL, DIMENSION (0:nm,ot,ng), INTENT (IN) :: field
  REAL, DIMENSION (ng), INTENT (IN) :: lon, lat

! local variables
  TYPE(ncfile) :: outfile
  TYPE(ncvar) :: nc_lats, nc_lons, nc_time, nc_flux
  TYPE(ncdim), TARGET, DIMENSION(3) :: nc_dims
  REAL, DIMENSION( nlon, nlat, nm* ot) :: outfield
  REAL, DIMENSION(:), ALLOCATABLE :: lats, lons, time
  INTEGER i, j, ig, i_m ! index variables

  outfield = 0.
  DO ig =1, ng
     i = NINT( MODULO( lon( ig) , 360.) *float( nlon) /360. +0.5)
     j = NINT((lat( ig) +90.)*float( nlat) /180. +0.5)
     DO i_m = 1, nm
        outfield( i, j,(i_m -1)* ot +1 : i_m * ot ) =  field( i_m, :, ig)
     END DO

  END DO

  ! if the output is monthly then normalize to per year for fluxes
!  IF( ot == 12) THEN
!     DO i = 1, nm * ot
!        j = MODULO(i-1,12)+ 1
!
!        outfield(:,:,i) = outfield(:,:,i) * 365./float(rdays(j))
!     END DO
!  END IF

  ! make lat, lon and time arrays
  ALLOCATE( lats( nlat))
  ALLOCATE( lons( nlon))
  ALLOCATE( time( ot * nm))
  DO i= 1, nm * ot
     time(i) = year0 + (i-0.5) / ot ! equal month lenghts, slight bug
  END DO
  DO i= 1,nlat
     lats( i) = -90. + (i-0.5) *180. /float (nlat)
  END DO
  DO i= 1, nlon
     lons( i) = (i -0.5) * 360./float( nlon)
  END DO
  outfile%name =TRIM( filename)
!  PRINT*,filename
  CALL nccreate(outfile, incmode = nf_clobber)
  nc_lats%name = 'latitude';  nc_lats%xtype  = nf_float
  nc_lats%ndims = 1; nc_lats%dims(1)%p => nc_dims(2)
  nc_lons%name = 'longitude';  nc_lons%xtype = nf_float
  nc_lons%ndims = 1; nc_lons%dims(1)%p =>  nc_dims(1)
  nc_time%name = 'time'; nc_time%xtype= nf_float
  nc_time%ndims = 1; nc_time%dims(1)%p => nc_dims(3)
  ! construct field name from filename
  i = INDEX(filename, "/", back=.TRUE.) +1
  j = INDEX( filename, ".", back=.TRUE.) -1
  nc_flux%name = filename( i: j); nc_flux%xtype = nf_float
  nc_flux%ndims = 3;
  DO i=1,3
     nc_flux%dims(i)%p => nc_dims(i)
  END DO

  nc_dims(1)%name = 'longitude'; nc_dims(1)%len = nlon
  nc_dims(3)%name = 'time'; nc_dims(3)%len = nm*ot
  nc_dims(2)%name ='latitude';  nc_dims(2)%len = nlat
  DO i=1,3
     CALL ncdimdef(outfile, nc_dims(i))
  END DO
  CALL ncvardef(outfile, nc_lons)

  CALL ncvardef(outfile, nc_lats)
  CALL ncvardef(outfile, nc_time)
  CALL ncvardef(outfile, nc_flux)
  CALL ncputatt(outfile, nc_lats, "units", "degrees_N")
  CALL ncputatt(outfile, nc_lons, "units", "degrees_E")
  CALL ncputatt(outfile, nc_time, 'units', 'years')
  CALL ncenddef(outfile)
  CALL ncprint(outfile, nc_lats, lats)
  CALL ncprint(outfile, nc_lons, lons)
  CALL ncprint(outfile, nc_time, time)
  CALL ncprint(outfile, nc_flux, outfield)
  CALL ncclose(outfile)
  DEALLOCATE( lats, lons, time)

  RETURN

END SUBROUTINE savefnc

SUBROUTINE savefnc_diurnal( filename, field,nm,ot,nhrs,nvp,binnm)
  ! pjr april 2002 for saving as netcdf output gridded
  USE mo_netcdf
  USE mo_namelist, ONLY : year0,day0  ! avoiding name conflict
  USE mo_constants, ONLY: rdays
  USE mo_grid, ONLY: nlon, nlat

  IMPLICIT NONE

  CHARACTER(len=*) :: filename
  CHARACTER(len=80) :: attval
! .. Scalar Arguments ..
  INTEGER, INTENT (IN) :: ot, nhrs, nvp, nm,  binnm

! .. Array Arguments ..
  REAL, DIMENSION (0:nm,ot,nhrs,nvp), INTENT (IN) :: field
!  REAL, DIMENSION (ng), INTENT (IN) :: lon, lat

! local variables
  TYPE(ncfile) :: outfile
  TYPE(ncvar) :: nc_vegp, nc_time, nc_flux
  TYPE(ncdim), TARGET, DIMENSION(2) :: nc_dims
  REAL, DIMENSION( nm* ot* nhrs, nvp ) :: outfield
  REAL, DIMENSION(:), ALLOCATABLE :: vegp, time
  INTEGER i, j, i_m, i_day, idoy, ihr  ! index variables
!  REAL, DIMENSION(:), ALLOCATABLE :: hrs

  outfield = 0.

  DO i_m= 1, nm
     DO i_day= 1,ot
        DO ihr= 1,nhrs
           i = (i_day-1)*24 + ihr
           outfield(i,:) = field(i_m,i_day,ihr,:)
        END DO
     END DO
  END DO

  ! make veg-point and time arrays
  ALLOCATE( vegp( nvp))
  ALLOCATE( time( ot * nm * nhrs))
!  ALLOCATE( doy( ot))
!  ALLOCATE( hrs(nhrs))

  ! determine approx mid-day of the month (Day Of Year)
!  DO i= 1,ot
!     doy(i) = SUM( rdays(1:i))       ! Note: not the same as the actual mid-period of forcings for which SCOPE is simulated
!  END DO

!  DO i= 1,nhrs
!     hrs(i) = i*1.
!  END DO

  ! time array as fractional year
!  DO i= 1, nm * ot * nhrs
!     doy = (i-1)/24 + 1
!     ihr = MOD(i-1, nhrs)+1 
!     time(i) = year0 + doy/365. + hrs(ihr)/24/365.     !/365. + hrs(ihr)/24/365. !
!  END DO
  DO i= 1, nm * ot * nhrs
     time(i) = i*1.
  END DO

  DO i= 1, nvp
     vegp(i) = i*1.  !field(1,1,1,i) 
  END DO
!  DO i= 1,nlat
!     lats( i) = -90. + (i-0.5) *180. /float (nlat)
!  END DO
!  DO i= 1, nlon
!     lons( i) = (i -0.5) * 360./float( nlon)
!  END DO
  outfile%name =TRIM( filename)
!  PRINT*,filename
  CALL nccreate(outfile, incmode = nf_clobber)
  nc_time%name = 'time'; nc_time%xtype= nf_float
  nc_time%ndims = 1; nc_time%dims(1)%p => nc_dims(1)
  nc_vegp%name = 'veg-point'; nc_vegp%xtype= nf_float
  nc_vegp%ndims = 1; nc_vegp%dims(1)%p => nc_dims(2)
  ! construct field name from filename
  i = INDEX(filename, "/", back=.TRUE.) +1
  j = INDEX( filename, ".", back=.TRUE.) -1
  nc_flux%name = filename( i: j); nc_flux%xtype = nf_float
  nc_flux%ndims = 2;
  nc_flux%dims(1)%p => nc_dims(1)
  nc_flux%dims(2)%p => nc_dims(2)

  nc_dims(1)%name = 'time'; nc_dims(1)%len = nm*ot*nhrs
  nc_dims(2)%name = 'veg-point'; nc_dims(2)%len = nvp
  DO i=1,2
     CALL ncdimdef(outfile, nc_dims(i))
  END DO

  write(attval,'(a40,2i4)') 'hours since midnight on (year,day) ',year0,day0

  CALL ncvardef(outfile, nc_time)
  CALL ncvardef(outfile, nc_vegp)
  CALL ncvardef(outfile, nc_flux)
  CALL ncputatt(outfile, nc_vegp, "units", "vegetation-point number")
  CALL ncputatt(outfile, nc_time, "comment", "first day of the averaging period")
  CALL ncputatt(outfile, nc_time, "units", attval)
  CALL ncenddef(outfile)
  CALL ncprint(outfile, nc_time, time)
  CALL ncprint(outfile, nc_vegp, vegp)
  CALL ncprint(outfile, nc_flux, outfield)
  CALL ncclose(outfile)
  DEALLOCATE( vegp, time)
!  DEALLOCATE( doy, hrs )

  RETURN

END SUBROUTINE savefnc_diurnal


SUBROUTINE savefbin(unit,field,lat,lon,nm,ot,ng,binnm)

  IMPLICIT NONE

! .. Scalar Arguments ..
  INTEGER, INTENT (IN) :: ot, ng, nm, unit, binnm


! .. Array Arguments ..
  REAL, DIMENSION (0:nm,ot,ng), INTENT (IN) :: field
  REAL, DIMENSION (ng), INTENT (IN) :: lon, lat
  REAL, DIMENSION(ng,binnm) :: binfield

! .. Local Scalars ..
  INTEGER :: i, ye, ts, mo
  CHARACTER (len=80) :: ft, nf  

!cccc unit        output channel
!cccc form        output format for one value (i.e. 'f5.1')
!cccc             must be 4 characters long
!cccc field       field of variables to be saved
!cccc lat         list of latitudes [deg N]
!cccc lon         list of longitudes [deg E]
!cccc nm          first dimension of 'field', number of run years
!cccc ot          output intervall, number of outvalues per year
!cccc ng          third dimension of 'field', number of grid points

  
  ts=0
  DO ye=1,nm
     DO mo=1,ot
        ts = ts + 1
        DO i=1,ng
           binfield(i,ts) = field(ye,mo,i)
        ENDDO
     ENDDO
  ENDDO

  IF (ts/=binnm) THEN
     PRINT*,'warning: something went wrong while writing the bin file'
  ENDIF

  WRITE (unit) ng,binnm
  WRITE (unit) binfield


  RETURN

END SUBROUTINE savefbin



SUBROUTINE savefield(unit,form,field,lat,lon,nm,ot,ng,year0)

  IMPLICIT NONE

! .. Scalar Arguments ..
  INTEGER, INTENT (IN) :: ot, ng, nm, unit, year0
  CHARACTER (len=4), INTENT (IN) :: form

! .. Array Arguments ..
  REAL, DIMENSION (0:nm,ot,ng), INTENT (INout) :: field
  REAL, DIMENSION (ng), INTENT (IN) :: lon, lat

! .. Local Scalars ..
  INTEGER :: i, ye, ts
  CHARACTER (len=80) :: ft, nf  

!cccc unit        output channel
!cccc form        output format for one value (i.e. 'f5.1')
!cccc             must be 4 characters long
!cccc field       field of variables to be saved
!cccc lat         list of latitudes [deg N]
!cccc lon         list of longitudes [deg E]
!cccc nm          first dimension of 'field', number of run years
!cccc ot          output intervall, number of outvalues per year
!cccc ng          third dimension of 'field', number of grid points

  
  IF (ot == 12) THEN 
     ft = '(12' // TRIM(form) // ')'
  ELSEIF (ot == 366) THEN 
     ft = '(366' // TRIM(form) // ')'
  ENDIF

  DO i = 1,ng
     WRITE(unit,*) lat(i),lon(i), nm, year0
     DO ye=1,nm
        WRITE (unit,TRIM(ft)) field(ye,:,i)
     ENDDO
  END DO

  RETURN

END SUBROUTINE savefield

END MODULE mo_io

SUBROUTINE  ignore
  RETURN
END SUBROUTINE ignore
