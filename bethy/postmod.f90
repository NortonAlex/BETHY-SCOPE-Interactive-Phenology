SUBROUTINE postmod( nx, x, fc)
  USE mo_diagnostics
  USE mo_io
  USE costf
  USE tm, ONLY : conc, f_tm, tm_deallocate
  USE mo_namelist, ONLY : year0, year1, year0_site,year1_site
  USE constants_pjr
  USE mo_pheno
  USE mo_hydro
  USE mo_climate
  USE mo_trafo
  USE mo_prog
  USE mo_grid, ONLY: sp, vp
  
  IMPLICIT NONE
  
  INTEGER              , INTENT(in) :: nx
  REAL, DIMENSION(nx), INTENT(inout) :: x
  REAL,                  INTENT(inout) :: fc
  REAL, DIMENSION(nx) :: par
  ! local variables
  INTEGER i_stat, i_year, i_month , i
  
  WRITE(6,*)
  WRITE(6,'(a,1x,e17.10)') ' function value ', fc
  WRITE(6,*)

  IF (n_stats>0) THEN
     WRITE(6,*) ' CO_2 concentrations for different stations from global run'
     WRITE(6,*) ' # of months simulated:', (year1 -year0 + 1) * months_per_year
!     WRITE(6,'(a15)') 'Name of station'
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
     END DO
     CALL tm_deallocate
     DEALLOCATE(paramap_global,prog_global)
  ELSE
     DEALLOCATE(conc,f_tm)
  ENDIF
  IF (n_sites>0) THEN
     WRITE(6,*)
     WRITE(6,*) ' Local run at eddy flux sites'
     WRITE(6,*) ' # of months simulated:', (year1_site -year0_site + 1) * months_per_year
     WRITE(6,*) ' sites:'
     
     DO i_stat  = 1,n_sites
        WRITE(6,*) site_names( i_stat)
     ENDDO
     DEALLOCATE(paramap_site,prog_sites)
     
  ENDIF

  call closeio
END SUBROUTINE postmod
