SUBROUTINE postadm( nvar, x, fc, adx )
  USE mo_diagnostics
  USE mo_io
  USE costf
  USE tm, ONLY : conc
  USE mo_namelist, ONLY : year0, year1, outdir, year1_site, year0_site
  USE constants_pjr
  USE mo_trafo

  IMPLICIT NONE

  INTEGER              , INTENT(in) :: nvar
  REAL, DIMENSION(nvar), INTENT(inout) :: x
! HEW_s
  REAL, DIMENSION(nvar) :: par
  REAL,                  INTENT(inout) :: fc
  REAL, DIMENSION(nvar), INTENT(inout) :: adx
  ! local variables
  INTEGER i_stat, i_year, i_month , i
  CHARACTER (len=80) :: filnam

  WRITE(6,*)
  WRITE(6,'(a,1x,e17.10)') ' function value ', fc
  WRITE(6,*)
  WRITE(6,'(a,1x,e17.10)') ' gradient       ', SQRT( SUM( adx*adx ) )
  WRITE(6,*)
  WRITE(6,*) 'optimized parameters in physical units' 
  WRITE(6,'(a)') 'number        initial      initial_unc up     predicted     gradient'

  par= x2p(x,x0,xpf,xpfa,xpfb,px0,p0su,a,nvar)   
  DO i = 1,nvar
     if (i>size(px0)) then
        WRITE(6,'(i3,5(e17.10))') i, 0., 0., 0., par(i),adx(i)        
     else
        WRITE(6,'(i3,5(e17.10))') i, px0(i), p0su(i), par(i),adx(i)
     endif
  END DO

  CALL model( nvar, x, fc)
  WRITE(6,*) 'value of the cf after call model ',fc


  IF (n_stats>0) THEN
     WRITE(6,*) ' CO_2 concentrations for different stations from global run'
     WRITE(6,*) ' # of months simulated:', (year1 -year0 + 1) * months_per_year
     WRITE(6,'(a15)') 'Name of station'
     WRITE(6,*) '  year.frac   observed  simulated  uncertainties'

     DO i_stat  = 1,n_stats
        !MAS output only for mlo
        IF (stat_names(i_stat) == 'mlo') THEN
           WRITE(6,*) 'station: ',stat_names( i_stat)
           ! FastOpt: index i replaced by i_stat
           ! compiler check had complained about this
           DO i_year = year0, year1
              DO i_month = 1, months_per_year
                 i = ( i_year - year0)* n_stats * months_per_year + (i_stat -1) &
                      & * months_per_year  + i_month
                 WRITE(6,'(f10.2,3f12.4)') REAL(i_year) +(i_month -0.5)/ &
                      & months_per_year, c_obs(i), conc(i), c_unc(i)
              END DO
           END DO
        ENDIF
        filnam=TRIM(outdir)//TRIM(stat_names(i_stat))//'_modconc.dat'
        PRINT*,filnam
        OPEN (99,file=filnam,status='unknown',form='formatted')
        WRITE(99,'(a20)') stat_names( i_stat)
        WRITE(99,*) (year1 -year0 + 1) * months_per_year, 3
        DO i_year = year0, year1
           DO i_month = 1, months_per_year
              i = ( i_year - year0)* n_stats * months_per_year + (i_stat -1) &
                   & * months_per_year  + i_month
              WRITE(99,'(f10.2,3f10.4)') REAL(i_year) +(i_month -0.5)/ &
                   & months_per_year,c_obs(i), conc(i), c_unc(i)
           END DO
        END DO
        CLOSE(99)        
     END DO
  ENDIF

  IF (n_sites>0) THEN
     WRITE(6,*)
     WRITE(6,*) ' Local run at eddy flux sites'
     WRITE(6,*) ' # of months simulated:', (year1_site -year0_site + 1) * months_per_year
     WRITE(6,*) ' sites:'
     
     DO i_stat  = 1,n_sites
        WRITE(6,*) site_names( i_stat)
     ENDDO     
  ENDIF

  print *, 'reached end of postadm subroutine'
  
END SUBROUTINE postadm
