MODULE mo_calendar
  IMPLICIT NONE

  INTEGER, ALLOCATABLE  :: amonth(:), ayear(:)
  INTEGER, ALLOCATABLE  :: fmonth(:), lmonth(:)
  INTEGER, ALLOCATABLE  :: fyear(:), lyear(:)
  INTEGER, ALLOCATABLE  :: idayint(:), daycount(:)
  INTEGER, ALLOCATABLE  :: spin(:), doy(:)
  INTEGER :: tdays, sdays, nds0, ndays_in, dspin


CONTAINS

  SUBROUTINE calendar(year0,year1,yearin0,yearin1,nspin,phint)    
    USE mo_constants
    IMPLICIT NONE
    INTEGER, INTENT(in) :: year0, year1
    INTEGER, INTENT(in) :: yearin0, yearin1
    INTEGER, INTENT(in) :: nspin, phint
    INTEGER :: i, j, m, yr, d, im, iy, ip, icount
    INTEGER :: nleap, nleap_before, nleap_in, rdays0


! ... calculate numbers of days for various input / model periods
    ! number of leap years over model  run
    nleap = 0
    DO yr = year0, year1
      IF (MODULO (yr,4)==0) nleap = nleap + 1
    END DO
    sdays = (year1-year0+1) * 365 + nleap
    ! number of leap years in input data before model start
    nleap_before = 0
    DO yr = yearin0, year0-1
      IF (MODULO (yr,4)==0) nleap_before = nleap_before + 1
    END DO
    nds0 = (year0-yearin0) * 365 + nleap_before
    ! number of leap years in input data
    nleap_in = 0
    DO yr = yearin0, yearin1
      IF (MODULO (yr,4)==0) nleap_in = nleap_in + 1
    END DO
    ndays_in = (yearin1-yearin0+1) * 365 + nleap_in

    ! add number of days for spin-up period, no leap years for spin-up
    dspin = nspin * 365
    tdays = sdays + dspin

! ... allocate calendar variables
    ALLOCATE (amonth(tdays),ayear(tdays))
    ALLOCATE (fmonth(tdays),lmonth(tdays))
    ALLOCATE (fyear(tdays),lyear(tdays))
    ALLOCATE (idayint(tdays), daycount(tdays))
    ALLOCATE (spin(tdays))
    ALLOCATE (doy(tdays))

    i = 1
    DO yr = 1, nspin
       DO m = 1, 12
          rdays0 = rdays(m)
          DO d = 1, rdays0
             spin(i) = i
             ayear(i) = 0
             IF (d==1) fmonth(i:i+rdays0-1)=i
             IF (d==rdays0) lmonth(i-rdays0+1:i)=i
             amonth(i)=m
             IF (m==1.and.d==1) fyear(i:i+365-1)=i
             IF (m==12.and.d==rdays0) lyear(i-365+1:i)=i
             i = i + 1
          ENDDO
       END DO
    END DO

    DO yr = year0, year1
       DO m = 1, 12
          rdays0 = rdays(m)
          IF (m==2 .AND. MODULO (yr,4)==0) rdays0=rdays(m)+1
          DO d = 1, rdays0
             spin(i) = 0
             ayear(i) = yr-year0+1
             IF (d==1) fmonth(i:i+rdays0-1)=i
             IF (d==rdays0) lmonth(i-rdays0+1:i)=i
             amonth(i)=m
             IF (MODULO (yr,4)==0) THEN
                IF (m==1.and.d==1) fyear(i:i+366-1)=i
                IF (m==12.and.d==rdays0) lyear(i-366+1:i)=i
             ELSE
                IF (m==1.and.d==1) fyear(i:i+365-1)=i
                IF (m==12.and.d==rdays0) lyear(i-365+1:i)=i
             ENDIF
             i = i + 1
          ENDDO
       END DO
    END DO

    icount=0
    DO i=1,dspin
       IF (MOD(i-1,phint)==0) THEN
          ip=i+(phint-1)
          if (ip>dspin) ip=dspin
          idayint(i:ip)=i         
          icount=icount+1
          daycount(i:ip)=icount         
       ENDIF       
    ENDDO    
    DO i=1,sdays
       IF (MOD(i-1,phint)==0) THEN
          ip=i+(phint-1)
          if (ip>sdays) ip=sdays
          idayint(i+dspin:ip+dspin)=i+dspin         
          icount=icount+1
          daycount(i+dspin:ip+dspin)=icount         
       ENDIF       
    ENDDO

    DO i=1,tdays
       IF (fyear(i)==i) j = 1
       doy(i) = j
       j = j + 1
    ENDDO



  END SUBROUTINE calendar

  SUBROUTINE deallocate_cal    
    DEALLOCATE (amonth, ayear)
    DEALLOCATE (fmonth, lmonth)
    DEALLOCATE (fyear, lyear)
    DEALLOCATE (idayint, daycount, spin)    
    DEALLOCATE (doy)
  END SUBROUTINE deallocate_cal

END MODULE mo_calendar
