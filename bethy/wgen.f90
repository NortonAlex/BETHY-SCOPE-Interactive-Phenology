    SUBROUTINE wgen(pp,fwet,gday,gpp,wet)

!HEW-CHG: rflg:      USE mo_constants, ONLY : ng, jdpyear, rflg, jdlast
      USE mo_constants, ONLY : ng, jdpyear, jdlast
      USE mo_namelist, ONLY : rflg
      IMPLICIT NONE

! .. Parameters ..
      REAL, PARAMETER :: e = 2.7182818
      INTEGER, PARAMETER :: ia = 16807, im = 2147483647, iq = 127773
      INTEGER, PARAMETER :: ir = 2836
      REAL, PARAMETER :: am = 1./im
!      REAL, PARAMETER :: mtrand = 2., rtmin = 0.3

! .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: gday

! .. Array Arguments ..
      REAL, DIMENSION (ng), INTENT (IN) :: fwet, pp
      REAL, DIMENSION (ng), INTENT (OUT) :: gpp
      INTEGER, DIMENSION (ng), INTENT (INOUT) :: wet

! .. Local Scalars ..
      REAL :: alpha, beta, c, fr, ppw, pwd, pwet, pww, rr, x
!      REAL :: rtd, rtw, fw
      INTEGER :: i, k


! .. Intrinsic Functions ..
      INTRINSIC alog, exp, float, max, mod

! .. Common blocks ..
      INTEGER :: idum

      COMMON /vran/idum

! ..
!CCCC Rainfall is modeled with a 1st order Markov chain model combined
!CCCC with gamma-distribution following Geng et al. (1986), Agr. For. Met.
!CCCC 36, 363-376.
!CCCC INPUT:
!cccc pp        mean precipitation for current day [mm / day]
!cccc fwet      fraction of wet days per month for current day [frac]
!cccc gday      Julian day
!cccc rmode     0: stochastic rainfall
!cccc           1 or more: repeated rainfall every rmode/fwet days
!cccc jdpyear   number of days per year (e.g. 365)
!cccc wet       array to store wet/dry flag per grid point
!CCCC OUTPUT:
!cccc gpp       daily preicipation for day 'gday' for grid points 1...ng [mm]
!CCCC INTERNAL:
!cccc pwd       probality of wet day following dry day |
!cccc pww       probality of wet day following wet day | -> Markov chain
!cccc pwet      actual probability of wet day
!cccc ppw       precipitation per wet day [mm]
!cccc wet(1..n) array storing last days wet state (1:wet, 0:dry)
!cccc alpha     [unitless] | parameters of gamma-distribution to model wet-day
!cccc beta      [mm/day]   | precipitation (<x>=alpha*beta, s^2=alpha*beta^2)
!      if (jdlast.eq.0) idum = 10248098
      IF (jdlast==0) idum = 13458010
      IF (rflg==0 .AND. gday/=MOD(jdlast+1,jdpyear)) THEN
        DO i = 1, ng
          wet(i) = 0
          k = idum/iq
          idum = ia*(idum-k*iq) - ir*k
          IF (idum<0) idum = idum + im
          IF (am*idum<=fwet(i)) wet(i) = 1
        END DO
      END IF
      jdlast = gday

      IF (rflg>0) THEN
        DO i = 1, ng
          fr = rflg/MAX(fwet(i),.033)
          IF (MOD(float(gday),fr)<1.) THEN
            gpp(i) = pp(i)*fr
          ELSE
            gpp(i) = 0.
          END IF
        END DO
      ELSE IF (rflg==0) THEN
        DO i = 1, ng
!           Markov chain model for precipitation
          pwd = 0.75*fwet(i)
          pww = 0.25 + pwd
          pwet = pww*wet(i) + pwd*(1-wet(i))
          k = idum/iq
          idum = ia*(idum-k*iq) - ir*k
          IF (idum<0) idum = idum + im
          IF (am*idum<=pwet .AND. pp(i)>0.) THEN
            ppw = pp(i)/MAX(fwet(i),0.033)
            beta = -2.16 + 1.83*ppw
            IF (ppw<3.) beta = 1.11*ppw
            alpha = ppw/beta
!              Gamma-distribution for wet days:
            c = e/(alpha+e)
1           k = idum/iq
            idum = ia*(idum-k*iq) - ir*k
            IF (idum<0) idum = idum + im
            x = am*idum
            IF (x<=c) THEN
              x = (x/c)**(1./alpha)
              rr = EXP(-x)
            ELSE
              x = -alog((x-c)/(1.-c)/e)
              rr = x**(alpha-1.)
            END IF
            k = idum/iq
            idum = ia*(idum-k*iq) - ir*k
            IF (idum<0) idum = idum + im
            IF (am*idum>=rr) GO TO 1
            gpp(i) = beta*x
            wet(i) = 1
          ELSE
            gpp(i) = 0.
            wet(i) = 0
          END IF
        END DO
      END IF

      RETURN
    END SUBROUTINE wgen



    SUBROUTINE wgen2(fwet,gpp,dtran,drt,mth)

! .. Use Statements ..
      USE mo_constants, ONLY : ng

      IMPLICIT NONE

! .. Parameters ..
      REAL, PARAMETER :: mtrand = 2., rtmin = 0.3

! .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: mth

! .. Array Arguments ..
      REAL, DIMENSION (ng), INTENT (IN) :: gpp
      REAL, DIMENSION (ng), INTENT (INOUT) :: drt, dtran
      REAL, DIMENSION (12,ng), INTENT (IN) :: fwet

! .. Local Scalars ..
      REAL :: fw, rtd, rtw
      INTEGER :: i

! .. Intrinsic Functions ..
      INTRINSIC max, min

      DO i = 1, ng
        fw = MAX(fwet(mth,i),0.033)
        rtd = (drt(i)-rtmin*fw)/MAX(1.-fw,0.001)
        rtd = MIN(rtd,1.)
        rtw = (drt(i)-rtd*(1.-fw))/fw
        IF (gpp(i)>0.) THEN
          dtran(i) = dtran(i)/(mtrand+(1.-mtrand)*fw)
          drt(i) = rtw
        ELSE
          dtran(i) = dtran(i)/(1.+(1./mtrand-1.)*fw)
          drt(i) = rtd
        END IF
      END DO

      RETURN
    END SUBROUTINE wgen2


