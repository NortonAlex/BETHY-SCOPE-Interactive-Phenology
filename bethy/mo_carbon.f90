MODULE mo_carbon
!------------------------------------------------------------------------
! WOLFGANG KNORR/GEORG HOFFMANN 050398
! MATTHIAS O'CUNTZ since 11/98
! Kalle Schnitzler MPI Hamburg, May 2000
!------------------------------------------------------------------------
  IMPLICIT NONE

    REAL, PARAMETER :: OX       = 0.21     ! OX      OXYGEN CONCENTRATION [MOL(O2)/MOL(AIR)]  
    REAL, PARAMETER :: FRDC3    = 0.011      ! FRDC3   RATIO OF DARK RESPIRATION TO "PVM" 
                               !         AT 25C for C3
    REAL, PARAMETER :: FRDC4    = 0.042      ! FRDC4   RATIO OF DARK RESPIRATION TO "PVM" 
                               !         AT 25C for C4
    REAL, PARAMETER :: THETA    = 0.83       ! THETA   CURVATURE PARAMETER
    REAL, PARAMETER :: eta      = 0.99       ! eta     colimitation shape parameter
    REAL, PARAMETER :: pvmmin   = 1.e-20     ! pvmmin  security parameter for pvm
    REAL, PARAMETER :: jmt      = 0.04       ! jmt     temperature scaling for Jmax

  
CONTAINS
!--------------------------------------------------------------------------
!
! photo
!
!--------------------------------------------------------------------------


  SUBROUTINE photo1 (ng,vp,psradsd,temp,paps, &
                  & zrhos,ppar,zpdir,zlai,zgc,pgs,zfpar,zfc,zassc,zraut,ts, &
                  & c4flg,ph,class,vm,jmf,zrphc,fautleaf,ccost, &
                  & EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0) 

!-----------------------------------------------------------------------
! used Module variables
!-----------------------------------------------------------------------
    USE mo_carparams
    USE mo_constants
!    USE mo_namelist
    USE mo_carvar
    USE mo_climate, ONLY: mu, spds, cpds
    USE mo_taf
    USE mo_grid, ONLY : gridp

! .. Arguments
    INTEGER, INTENT(in) :: ng, vp, ts
    REAL, DIMENSION(ng), INTENT(in) :: psradsd, ppar, zpdir, paps
    REAL, DIMENSION(ng), INTENT(in) :: temp
    REAL, DIMENSION(vp), INTENT(in) :: EC, EO, EV, ER, EK, tgam
    REAL, DIMENSION(vp), INTENT(in) :: alpha, alc4, kc0, ko0
    REAL, DIMENSION(vp), INTENT(in) :: zlai,fautleaf,ccost, zrhos,zfc
    REAL, DIMENSION(vp), INTENT(inout) :: zgc
    REAL, DIMENSION(vp,nl), INTENT(inout) :: pgs
    REAL, DIMENSION(vp), INTENT(out) :: zfpar,zassc,zraut,zrphc
    INTEGER, DIMENSION (vp), INTENT(in) :: c4flg, ph, class
    REAL, DIMENSION(vp), INTENT(in) :: jmf
    REAL, DIMENSION(vp), INTENT(in) :: vm

! .. Local variables
    INTEGER :: klon, nveglist, jl, il, jj
    INTEGER, DIMENSION(vp) :: niveglst
    REAL, DIMENSION(vp) :: zfcmax, zrdc
    REAL, DIMENSION(vp) :: zrhospar
    REAL, DIMENSION(vp,nl) :: PLAIL, ZNSCL, ZAPAR
    REAL, DIMENSION(vp,nl) :: pass, prd, zrph
    REAL, DIMENSION(ng) :: omega, asp
    REAL, DIMENSION(ng) :: zirrin, zpar, coszen
    REAL, DIMENSION(0:nl) :: dl
    REAL :: zrmnt, zrcon

    coszen(:)=mu(ts,:)
    klon=vp
    nveglist=vp
    ! single scattering albedo of leaves
    omega=0.12
    ! aspect ratio of vegetation
    asp=0.
    
    DO JL=1,klon
       ZFCMAX(JL)     = 0.
       ZFPAR(JL)      = 0.
       ZASSC(JL)      = 0.
       ZRDC(JL)       = 0.
       ZRPHC(JL)      = 0.
       ZRAUT(JL)      = 0.
    END DO

    DO  IL = 1, NL
       DO JL = 1, KLON
          PLAIL(JL,IL) = 0.
          ZNSCL(JL,IL) = 0.
          ZAPAR(JL,IL) = 0.
          PASS(JL,IL)  = 0.
          PRD(JL,IL)   = 0.
          ZRPH(JL,IL)  = 0.
       END DO       
    END DO

!----------------------------------------------------------------
! borders of the NL canopy layers in terms of LAI (for integration)
!----------------------------------------------------------------

    DO jl = 0, nl
       DL(jl) = FLOAT(jl) / FLOAT(nl)
    END DO

!----------------------------------------------------------------
! 1.2  SET UP LIST OF VEGETATED GRID POINTS 
!----------------------------------------------------------------
    NVEGLIST=klon
    DO JL = 1, klon
       niveglst(JL) = JL
    ENDDO

    CALL NSCALE (ng,gridp,znscl, zlai, dl, klon, NVEGLIST, &
         & niVEGLST, nl, class, ph, spds, cpds, lailim0,laimax0)

!-----------------------------------------------------------------
! 2.2     SET SOIL TYPE TO 'MEDIUM' AND DETERMINE SOIL REFLECTANCE
!------------------------------------------------------------------     
!     COMPUTE SOIL REFLECTANCE IN PAR FROM SOIL LINE               
!     BY PRICE AND BAUSCH, REMOTE SENS. ENVIRON. 52, 55-63, 1995   
!     UNTIL NOW, ONLY MEDIUM SOIL TYPE AVAILABLE
!-----------------------------------------------------------------

    DO Jl = 1, NVEGLIST
       Jj = gridp(Jl)
       
!-----------------------------------------------------------------
!      SOIL SEEMS TO BE DARKER IN THE PAR BAND
!-----------------------------------------------------------------
       ZRHOSPAR(JL) = 0.92 * ZRHOS(JL) - 0.015
       ZRHOSPAR(JL) = MAX( ZRHOSPAR(JL), ZRHOSPARMIN )

!----------------------------------------------------------------
!     CONVERT ZPAR FROM W/M^2 TO MOL(PHOTONS)/(M^2 S)
!----------------------------------------------------------------
       ZPAR(Jj) = PPAR(Jj) / EPAR 
       ZIRRIN(Jj) = PSRADSD(Jj) / EPAR

!---------------------------------------------------------------------
!     SET FRACTIONAL COVER "ZFC"                                    
!     COMPUTE ABSORBED PAR PER LEAF AREA "ZAPAR" IN CANOPY LAYERS   
!      AND FPAR "ZFPAR"                                             
!---------------------------------------------------------------------
       ZFCMAX(JL) = FCMAX0

 END DO
    

!---------------------------------------------------------------------- 
!    CALCULATES THE ABSORBED PAR PER LEAF AERA PER CANOPY LAYER, ZAPAR    
!     WITH TWO FLUX APPROXIMATION, NORMALIZED TO INCOMING RADIATION = 1   
!    CALCULATES LAI PER CANOPY LAYER ZLAIL                                
!----------------------------------------------------------------------

    CALL FAPARL (ng,gridp,klon,nl,ZLAI,ZFC,ZFCMAX,ZRHOSPAR,coszen,ZPDIR, &
         & PLAIL,ZAPAR,OMEGA,ASP,DL, &
         & ZZENITMIN,ZFCMIN,ZLAIMIN,NIVEGLST,NVEGLIST)

!$TAF store  PLAIL,ZAPAR  = carbon_tape, key = keydiurnal

!----------------------------------------------------------------------
!    CALCULATES FPAR OUT OF ABSORBED PAR PER CANOPY LAYER ZAPAr
!---------------------------------------------------------------------- 
!$taf loop = parallel
    DO IL = 1, NL
!$taf loop = parallel
       DO JJ = 1, NVEGLIST
          JL = NIVEGLST(JJ)
          ZFPAR(JL) = ZFPAR(JL) + ZAPAR(JL,IL)
       ENDDO
    ENDDO     
 
!RG split loop for TAF
!----------------------------------------------------------------------
!    CONVERT 'ZAPAR' FROM RELATIVE PAR ABSORPTION (I.E. NORMALIZED TO
!    INCOMING RADIATION = 1) TO MOL(ABSORBED PHOTONS)/(M^2(LEAF AREA) S)
!----------------------------------------------------------------------
!$taf loop = parallel
    DO IL = 1, NL
!$taf loop = parallel
       DO Jl = 1, NVEGLIST
          Jj = gridp(Jl)
          ZAPAR(JL,IL) = ZAPAR(JL,IL) * ZPAR(Jj) &
               & / PLAIL(JL,IL) 
          
       ENDDO
    ENDDO     
            
!----------------------------------------------------------------------
!                                                                            
!     CALCULATE UNSTRESSED PHOTOSYNTHESIS, I.E. WITHOUT WATER LIMITATION   
!     SET CANOPY TEMPERATURE "PTV" TO FIRST-LAYER SOIL TEMPERATURE         
!     SET LEAF-INTERNAL CO2 CONCENTRATION "PCI" TO STANDARD VALUE "CI0"    
!     CALCULATE UNSTRESSED STOMATAL CONDUCTANCE "ZGS" BY CALCULATING       
!     FARQUHAR/COLLATZ MODEL AT FIXED "PCI"                               
!                                                                          
!----------------------------------------------------------------------
! SETS CI TO PREDEFINED VALUES FOR FIRST CALL OF 'PHYSN'
            
    zgc=0.
    zci=0.
    
    DO IL = 1, NL
       DO JJ = 1, NVEGLIST
          JL = NIVEGLST(JJ)
             
!--------------------------------------------------------------------------
!C CI IS CI0 WITH 0.87CA FOR C3 (BEERLING AND QUICK, 1995) 
!C  AND 0.67CA FOR C4 (SCHULZE ET AL., 1994)
!C CIO VALUES SET IN 'INITVEGDATA',    KNORR, 1997, S.44 (41E)
!-------------------------------------------------------------------------
          PCI(JL,IL) = (FCI1C3 * (1.-C4FLG(JL)) & 
               & + FCI1C4 * C4FLG(JL)) * CA
       END DO
    END DO
      
!$taf loop = parallel
    DO IL = 1, nl

!------------------------------------------------------------------------    
       CALL  PHSYN1 (ng,gridp,ZAPAR(1:klon,IL), PCI(1:klon,IL), PGS(1:klon,IL), &
            & temp, PAPS, ZNSCL(1:klon,IL), PASS(1:klon,IL), &
            & PRD(1:klon,IL), ZRPH(1:klon,IL), VM(1:klon), &
            & JMF(1:klon), C4FLG(1:klon), JMAXMIN, &
            & ZGSMIN, CA, ar, KLON, NIVEGLST, NVEGLIST, &
            & ZIRRIN,ts,EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0,coszen)

!$TAF STORE pgs (:,IL)    = carbon_tape_nl, key = IL+(keydiurnal-1)*nl

!------------------------------------------------------------------------
! GC0, CANOPY CONDUCTANCE AS THE INTEGRAL OF THE STOMATAL CONDUCTANCE OVER
!      THE LEAF AREa
!------------------------------------------------------------------------
!$taf loop = parallel
       DO JJ = 1, NVEGLIST
          JL = NIVEGLST(JJ)
          ZGC(JL) = ZGC(JL) + PGS(JL,IL) * PLAIL(JL,IL) 
       ENDDO

    ENDDO


!FastOpt !$TAF STORE ZGC    = carbon_tape, key = keydiurnal

!!$taf loop = parallel
! WOK-CHG-090401 should be taken care of through ZGSMIN
!    DO JJ = 1, NVEGLIST
!       JL = NIVEGLST(JJ)
!       ZGC(JL) = MAX( ZGC(JL), ZGCMIN )
!    ENDDO

  END SUBROUTINE photo1



  SUBROUTINE photo2 (ng,vp, psradsd,ptv,paps, &
                  & zrhos,ppar,zpdir,zlai,zgc,pgs,zfpar,zfc,zassc,zraut,ts, &
                  & c4flg,ph,class,vm,jmf,zrphc,fautleaf,ccost, &
                  & EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0,zgrowth,zmaint,dapar) 

!-----------------------------------------------------------------------
! used Module variables
!-----------------------------------------------------------------------
    USE mo_carparams
    USE mo_constants
    USE mo_carvar
    USE mo_climate, ONLY: mu, spds, cpds
    USE mo_taf
    USE mo_grid, ONLY : gridp

! .. Arguments
    INTEGER, INTENT(in) :: ng, vp, ts
    REAL, DIMENSION(ng), INTENT(in) :: psradsd, ppar, zpdir, paps
    REAL, DIMENSION(vp), INTENT(in) :: EC, EO, EV, ER, EK, tgam
    REAL, DIMENSION(vp), INTENT(in) :: alpha, alc4, kc0, ko0
    REAL, DIMENSION(vp), INTENT(in) :: zlai,fautleaf,ccost, zrhos, zfc
    REAL, DIMENSION(vp), INTENT(inout) :: zgc
    REAL, DIMENSION(vp,nl), INTENT(inout) :: pgs
    REAL, DIMENSION(vp), INTENT(out) :: zfpar,zassc,zraut,zrphc,zgrowth,zmaint
    REAL, DIMENSION(vp), INTENT(out) :: dapar
    INTEGER, DIMENSION (vp), INTENT(in) :: c4flg, ph, class
    REAL, DIMENSION(vp), INTENT(in) :: jmf, ptv
    REAL, DIMENSION(vp), INTENT(in) :: vm

! .. Local variables
    INTEGER :: klon, nveglist, jl, il, jj
    INTEGER, DIMENSION(vp) :: niveglst
    REAL, DIMENSION(vp) :: zfcmax, zrdc
    REAL, DIMENSION(vp) :: zrhospar
    REAL, DIMENSION(vp,nl) :: PLAIL, ZNSCL, ZAPAR
    REAL, DIMENSION(vp,nl) :: pass, prd, zrph
    REAL, DIMENSION(ng) :: omega, asp
    REAL, DIMENSION(ng) :: zirrin, zpar,coszen
    REAL, DIMENSION(0:nl) :: dl
    REAL :: zrmnt, zrcon, fcinh

    coszen(:)=mu(ts,:)
    klon=vp
    nveglist=vp
    ! single scattering albedo of leaves
    omega=0.12
    ! aspect ratio of vegetation
    asp=0.
    
    DO JL=1,klon
       ZFCMAX(JL)     = 0.
       ZFPAR(JL)      = 0.
       ZASSC(JL)      = 0.
       ZRDC(JL)       = 0.
       ZRPHC(JL)      = 0.
       ZRAUT(JL)      = 0.
       DAPAR(JL)      = 0.
    END DO

    DO  IL = 1, NL
       DO JL = 1, KLON
          PLAIL(JL,IL) = 0.
          ZNSCL(JL,IL) = 0.
          ZAPAR(JL,IL) = 0.
          PASS(JL,IL)  = 0.
          PRD(JL,IL)   = 0.
          ZRPH(JL,IL)  = 0.
       END DO       
    END DO

!----------------------------------------------------------------
! borders of the NL canopy layers in terms of LAI (for integration)
!----------------------------------------------------------------

    DO jl = 0, nl
       DL(jl) = FLOAT(jl) / FLOAT(nl)
    END DO

!----------------------------------------------------------------
! 1.2  SET UP LIST OF VEGETATED GRIP POINTS 
!----------------------------------------------------------------
    NVEGLIST=klon
    DO JL = 1, klon
       niveglst(JL) = JL
    ENDDO


    CALL NSCALE (ng,gridp,znscl, zlai, dl, klon, NVEGLIST, &
         & niVEGLST, nl, class, ph, spds, cpds, lailim0,laimax0)

!-----------------------------------------------------------------
! 2.2     SET SOIL TYPE TO 'MEDIUM' AND DETERMINE SOIL REFLECTANCE
!------------------------------------------------------------------     
!     COMPUTE SOIL REFLECTANCE IN PAR FROM SOIL LINE               
!     BY PRICE AND BAUSCH, REMOTE SENS. ENVIRON. 52, 55-63, 1995   
!     UNTIL NOW, ONLY MEDIUM SOIL TYPE AVAILABLE
!-----------------------------------------------------------------

    DO Jl = 1, NVEGLIST
       Jj = gridp(Jl)

!-----------------------------------------------------------------
!      SOIL SEEMS TO BE DARKER IN THE PAR BAND
!-----------------------------------------------------------------
       ZRHOSPAR(JL) = 0.92 * ZRHOS(JL) - 0.015
       ZRHOSPAR(JL) = MAX( ZRHOSPAR(JL), ZRHOSPARMIN )

!----------------------------------------------------------------
!     CONVERT ZPAR FROM W/M^2 TO MOL(PHOTONS)/(M^2 S)
!----------------------------------------------------------------
       ZPAR(Jj) = PPAR(Jj) / EPAR 
       ZIRRIN(Jj) = PSRADSD(Jj) / EPAR

!---------------------------------------------------------------------
!     SET FRACTIONAL COVER "ZFC"                                    
!     COMPUTE ABSORBED PAR PER LEAF AREA "ZAPAR" IN CANOPY LAYERS   
!      AND FPAR "ZFPAR"                                             
!---------------------------------------------------------------------
       ZFCMAX(JL) = FCMAX0

!---------------------------------------------------------------------
! TRIES TO AVOID WORST CASE, IF FC<< => LAI/FC >> 
! Fractional vegetation cover,fc
!  fc = fcMax (0.9) * Lai (max of the year, not here) / Lai0 (3) IF Lai<Lai0
!  fc = fcMax                                                    IF Lai>Lai0
!  PVGRAT (YN(198)='VGRAT' ! VEGETATION RATIO  US151093.6)
!         ZFC(JL)=AMAX1(PVGRAT(JL),RLAI(JL)/LAILIM0*FCMAX0)
!         ZFC(JL)=AMIN1(ZFC(JL) , FCMAX0)
!  (0.1 for numerical reasons: often (3+3)/2=2.99...)
!--------------------------------------------------------------------
!!$       IF ( ZLAI(JL).LT.( LAILIM0 - 0.1 ) ) THEN
!!$          ZFC(JL) = ZLAI(JL) / LAILIM0 * FCMAX0
!!$       ELSE 
!!$          ZFC(JL) = FCMAX0
!!$       ENDIF
!!$       ZFC(JL) = MAX( ZFC(JL), ZFCMIN )
!WOK-070705 handling of frational cover moved to mo_pheno
!       zfc(jl)=zfcmax(jl) 
    END DO
    

!---------------------------------------------------------------------- 
!    CALCULATES THE ABSORBED PAR PER LEAF AERA PER CANOPY LAYER, ZAPAR    
!     WITH TWO FLUX APPROXIMATION, NORMALIZED TO INCOMING RADIATION = 1   
!    CALCULATES LAI PER CANOPY LAYER ZLAIL                                
!----------------------------------------------------------------------

    CALL FAPARL (ng,gridp,klon,nl,ZLAI,ZFC,ZFCMAX,ZRHOSPAR,coszen,ZPDIR, &
         & PLAIL,ZAPAR,OMEGA,ASP,DL, &
         & ZZENITMIN,ZFCMIN,ZLAIMIN,NIVEGLST,NVEGLIST)

!$TAF store  PLAIL,ZAPAR  = carbon_tape, key = keydiurnal

!----------------------------------------------------------------------
!    CALCULATES FPAR OUT OF ABSORBED PAR PER CANOPY LAYER ZAPAr
!---------------------------------------------------------------------- 
!$taf loop = parallel
    DO IL = 1, NL
!$taf loop = parallel
       DO JJ = 1, NVEGLIST
          JL = NIVEGLST(JJ)
          ZFPAR(JL) = ZFPAR(JL) + ZAPAR(JL,IL)
       ENDDO
    ENDDO     

    ! Calculate APAR per vegetation point (umol m-2 s-1)
    DO IL = 1,NL
       DO JL = 1,NVEGLIST
          JJ = gridp(JL)
          DAPAR(JL) = DAPAR(JL) + ZAPAR(JL,IL) * ZPAR (JJ) *1e6  ! x1e6 to convert mol photons to umol photons
       ENDDO
    ENDDO

 
!RG split loop for TAF
!----------------------------------------------------------------------
!    CONVERT 'ZAPAR' FROM RELATIVE PAR ABSORPTION (I.E. NORMALIZED TO
!    INCOMING RADIATION = 1) TO MOL(ABSORBED PHOTONS)/(M^2(LEAF AREA) S)
!----------------------------------------------------------------------
!$taf loop = parallel
    DO IL = 1, NL
!$taf loop = parallel
       DO Jl = 1, NVEGLIST
          Jj = gridp(jl)
          ZAPAR(JL,IL) = ZAPAR(JL,IL) * ZPAR(Jj) &
               & / PLAIL(JL,IL) 
          
       ENDDO
    ENDDO     
    
!------------------------------------------------------------------------
!                                                                         
! COMPUTE ACTUAL PHOTOSYNTHESIS "ZASS", DARK RESPIRATION "ZRD"                
! AND PHOTORESPIRATION "ZRPH"                                              
! FOR EACH CANOPY LAYER AT CANOPY TEMP. "PTV" AND STOMATAL CONDUCTANCE "ZGS". 
! "PASSC" AND "PRDC" AND "PRPHC" ARE THE RESPECTIVE INTEGRALS OVER THE LEAF 
! AREA INDEX.  
!
!------------------------------------------------------------------------
!$TAF STORE pgs    = carbon_tape, key = keydiurnal

!$taf loop = parallel
    DO IL = 1, NL

!------------------------------------------------------------------------
! Calculate ASS, RD, RPH, CI with input GS, TC for each canopy layer IL
!        ACTUAL PHOTOSYNTHESIS:
!------------------------------------------------------------------------
       CALL  PHSYN2 (ng,gridp,ZAPAR(1:klon,IL),PCI(1:klon,IL),PGS(1:klon,IL), &
            & PTV,PAPS,ZNSCL(1:klon,IL),PASS(1:klon,IL), &
            & PRD(1:klon,IL),ZRPH(1:klon,IL),VM(1:klon), &
            & JMF(1:klon),C4FLG(1:klon),JMAXMIN, &
            & ZGSMIN,CA,ar,KLON,NIVEGLST,NVEGLIST, &
            & ZIRRIN,ts,EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0,coszen)



!$TAF STORE PASS(:,IL)   = carbon_tape_nl, key = il+(keydiurnal-1)*nl
!$TAF STORE PRD (:,IL)   = carbon_tape_nl, key = il+(keydiurnal-1)*nl

!------------------------------------------------------------------------
       DO JJ = 1, NVEGLIST
          JL = NIVEGLST(JJ)

!------------------------------------------------------------------------
! Integrate over Canopy
!------------------------------------------------------------------------
          ZASSC(JL) = ZASSC(JL) + PASS(JL,IL) * PLAIL(JL,IL)
          ZRDC(JL)  = ZRDC(JL)  + PRD(JL,IL)  * PLAIL(JL,IL)
          ZRPHC(JL) = ZRPHC(JL) + ZRPH(JL,IL) * PLAIL(JL,IL)
                  
!------------------------------------------------------------------------
          ZCI(JL)   = ZCI(JL)   + PCI(JL,IL)  * PLAIL(JL,IL)

       END DO  ! Longitudes for Integrartion
    END DO     ! Canopy layers for calling PHYSN


!$TAF STORE ZASSC, ZRDC  = carbon_tape, key = keydiurnal

!------------------------------------------------------------------------
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!------------------------------------------------------------------------
!*
!*     ACCUMULATE IMPORTANT OUTPUT FIELDS:             
!*       RGPP AS "GGP"                                 
!*       AUTOTROPHIC RESPIRATION AS "RAUTO"            
!*       LEAF RESPIRATION AS "RLEAF"                   
!*       PHOTORESPIRATION AS "RPHOTO"                  
!*       SOIL OR HETEROTROPHIC RESPIRATION AS "RHET"   
!------------------------------------------------------------------------
!$taf loop = parallel
    DO JJ = 1, NVEGLIST
       JL = NIVEGLST(JJ)

!------------------------------------------------------------------------
!* PRDC, LEAF OR DARK RESPIRATION, CALCULATED WITH Arrhenius formula
!*      40% (FAUTLEAF) OF THE MAINTENANCE RESPIRATION (ZRMNT)
!*       TAKES PLACE IN THE LEAVES (RYAN, 1991; KNORR, 1997, S.53FF)
!*        (ALL STILL IN MOL CO2/M^2/S), i.e. the dark respiration is
!*       part of the maintenance respiration
!------------------------------------------------------------------------ 
       ZRMNT = ZRDC(JL) / FAUTLEAF(JL)

!------------------------------------------------------------------------
!* ZRCON, GROWTH RESPIRATION, 25gC PER gC BIOMASS PRODUCED (CCOST - 1)(RYAN, 1991)
!*      I.E., ZRCON=(CCOST-1)*NPP=(CCOST-1)*(GPP-ZRMNT-ZRCON)
!*            ZRCON=(CCOST-1)*(GPP-ZRMNT) / CCOST
!------------------------------------------------------------------------
      ! WOK 2008-11-08 problems with minumum function here
      ! NOTE: removed it makes 'ccost' redundant.
!       ZRCON = ( ZASSC(JL) - ZRMNT ) * ( CCOST(jl) - 1. ) / CCOST(jl)
!------------------------------------------------------------------------
!* REDEFINITION OF PARAMETER 'CCOST' AS CCOST_NEW = (CCOST_OLD-1) / CCOST_OLD
!*      REASON: PROBLEMS WITH CCOST_OLD -> 0 OR < 1
!* DEFAULT VALUE OF CCOST_NEW = (1.25 - 1) / 1.25 = (1/4) / (5/4) = 0.2
!------------------------------------------------------------------------
       ZRCON = ( ZASSC(JL) - ZRMNT ) * CCOST(jl)
      ! ZRCON = MAX( ZRCON, ZGROWTHRESMIN )
      ! WOK 2008-11-08 all these did not work:
      ! ZRCON = ZRCON / (1. + exp(-ZRCON/1e-7))
      ! ZRCON = (ZRCON + SQRT(ZRCON**2+(2*ZGROWTHRESMIN*ZRMNT)**2)) / 2.
      ! ZRCON = (ZRCON + SQRT(ZRCON**2+(2*0.1*ZRMNT)**2)) / 2.
      ! if (jl==353) write (15,*) ZRCON, ZRMNT, ZGROWTHRESMIN

!------------------------------------------------------------------------
!* RAUT, AUTOTROPHIC RESPIRATION, SUM OF MAINTENANCE, GROWTH RESP.
!------------------------------------------------------------------------
       ZRAUT(JL) = ZRMNT + ZRCON

!------------------------------------------------------------------------
! convert from mol CO2 / m2 / s to g C / m2 / s
!------------------------------------------------------------------------
       fcinh = coldinhib(ptv(jl))
       zassc(jl)=(zassc(jl))*mmc*1000.  *fcinh
       zraut(jl)=zraut(jl)*mmc*1000.  *fcinh
       zgrowth(jl)=zrcon*mmc *1000. *fcinh
       zmaint(jl)=zrmnt*mmc *1000. *fcinh

!MAS-ADD-BEG-040220  reset resp, if lai is 0.            
       if (zlai(jl)==0) then
          zassc(jl)=0.
          zraut(jl)=0.
          zgrowth(jl)=0.
          zmaint(jl)=0.
       endif
!MAS-ADD-END-040220  reset resp, if lai is 0.
               

    END DO
     
  END SUBROUTINE photo2




!----------------------------------------------------------------------------
  SUBROUTINE PHSYN1 (ng,gridl,PAR, CI, GS, TC, P, NSCL, A, RD, RPH, PVM, PJMF, &
       &  PC4FLG, PJMAXMIN, PGSMIN, PCA, R, N, ILIST, &
       &  NLIST, PIRRIN,ts,EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0,coszen)

! .. Use Statements ..
  USE mo_helper, ONLY: maxx

!---------------------------------------------------------------------------- 
!    INPUT:                                                                      
!    PAR     ABSORBED PAR [MOL(PHOTONS) / M^2 S]                                 
!!    CI IF KFLG=0     CO2 CONCENTRATION INSIDE THE LEAF [MOL(CO2) / MOL(AIR)]    
!!    GS IF KFLG=1    STOMATAL CONDUCTANCE FOR WATER VAPOUR [M / S]               
!    TC      LEAF AND CANOPY TEMPERATURE [DEG C]                                 
!    P       AIR PRESSURE [PA]                                                   
!    NSCL    SCALING FACTOR FOR "PVM" AND "PJM"                                  
!    KFLG    FLAG FOR DIFFERENT COMPUTATION MODES                                
!            0: COMPUTE PHOTOS. & STOMATAL COND. AT GIVEN 'CI'                   
!            1: COMPUTE PHOTOS. & 'CI'  AT GIVEN  STOMATAL COND.                 
!    PVM     MAXIMUM RATE OF CARBOXYLATION [MOL(CO2) / M^2 S]                    
!    PJMF     Jmax/Vmax
!    PC4FLG  PHOTOSYNTHETIC PATHWAY                                              
!            0: C3; 1: C4                                                        
!    PJMAXMIN MINIMUM OF MAXIMUM CARBOXYLATION RATE               
!    PGSMIN   MINIMUM STOMATAL CONDUCTANCE                                       
!    PCA      ATMOSPHERIC CO2 MIXING RATIO [MOL(CO2)/MOL(AIR)]                  
!    R        GAS CONSTANT [J / K MOL]                                          
!    OUTPUT:                                                                    
!    A       GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2 S]                             
!    RD      DARK RESPIRATION OF LEAF [MOL(CO2) / M^2 S]                         
!    RPH     PHOTORESPIRATION [MOL(CO2) / M^2 S]                                 
!    GS IF KFLG=0    STOMATAL CONDUCTANCE FOR WATER VAPOUR [M / S]               
!    CI IF KFLG=1    CO2 CONCENTRATION INSIDE THE LEAF [MOL(CO2) / MOL(AIR)]     
!                                                                                
!    INTERNAL:                                                                   
!    PCA      CO2 CONCENTRATION OF AMBIENT AIR [MOL(CO2) / MOL(AIR)]             
!                                                                                
!    PHYSIOLOGY:                                                                 
!    C3 PLANTS: FARQHUAR, G.D., S. VON CAEMMERER AND J.A. BERRY, 1980.           
!               A BIOCHEMICAL MODEL OF PHOTOYNTHESIS IN LEAVES OF C3 SPECIES.    
!               PLANTA 149, 78-90.                                               
!    PVM      MAXIMUM CARBOXYLATION RATE AT 25C [MOL(CO2)/M^2 S]                 
!    PJM      MAXIMUM RATE OF ELECTRON TRANSPORT AT 25C [MOL(CO2)/M^2 S]         
!    VCMAX   MAXIMUM CARBOXYLATION RATE [MOL(CO2)/M^2 S]                         
!    JMAX    MAXIMUM RATE OF ELECTRON TRANSPORT [MOL(CO2)/M^2 S]                 
!    ALPHA   EFFICIENCY OF OF PHOTON CAPTURE                                     
!    OX      OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)]                           
!    KC0     MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)]      
!    KO0     MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)]        
!    GAM     COMPENSATION POINT WITHOUT DARK RESPIRATION [MOL(CO2) / MOL(AIR)]   
!    EC      ACTIVATION ENERGY FOR KC [J / MOL]                                  
!    EO      ACTIVATION ENERGY FOR KO [J / MOL]                                  
!    EV      ACTIVATION ENERGY FOR VCMAX [J / MOL]                               
!    ER      ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL]                    
!    AV      CONSTANT FOR HIGH-TEMPERATURE INHIBITION OF VCMAX [J / MOL]         
!    BV      CONSTANT FOR HIGH-TEMPERATURE INHIBITION OF VCMAX [J / MOL K]       
!    FRDC3   RATIO OF DARK RESPIRATION TO "PVM" AT 25C                           
!    C4 PLANTS: COLLATZ, G.J., M. RIBAS-CARBO AND J.A. BERRY, 1992.              
!               COUPLED PHOTOSYNTHESIS-STOMATAL CONDUCTANCE MODEL FOR LEAVES     
!               OF C4 PLANTS. AUST. J. PLANT PHYSIOL. 19, 519-538.               
!    ALC4    EFFECTIVE QUANTUM EFFICIENCY                                        
!    K       CO2 SPECIFICITY OF PEPCASE                                          
!    THETA   CURVATURE PARAMETER                                                 
!                                                                                
!----------------------------------------------------------------------------      


    INTEGER :: n, ts,ng
    REAL    PAR(N), CI(N), GS(N), TC(ng)
    REAL    P(ng), NSCL(N)
    REAL    A(N), RD(N), RPH(N)
    REAL    PVM(N), PJMF(N)
    REAL    PJMAXMIN, PGSMIN, PCA, R
    INTEGER  pc4flg(n)
    INTEGER ILIST(N), NLIST,gridl(n)
    REAL    PIRRIN(ng), coszen(ng)
    REAL, DIMENSION(n) :: EC, EO, EV, ER, EK, tgam
    REAL, DIMENSION(n) :: alpha, alc4, kc0, ko0

    REAL    dummy
    REAL    K
    REAL    VCMAX, KC, KO, GAM, JE, JC, JMAX, J0, J1
    REAL    T, T0, T1
    REAL    G0, K1, W1, K2, W2, B, C
    INTEGER JL, JJ


!-----------------------------------------------------------------------
! dummy assignments for TAF
!-----------------------------------------------------------------------
    A(1) = A(1)
    RD(1) = RD(1)

!---------------------------------------------------------------------------------
!                  FIRST ENTRY, TC=TA, CI0 => GC0, AC0                 
!---------------------------------------------------------------------------------  


!$taf loop = parallel
!!    IF (KFLG.EQ.0) THEN          ! First ENTRY, calculate gs0 from TC = T
!    print *,'            DO LOOP over vp: ', NLIST
       
    DO Jl = 1, NLIST
       Jj = gridl(Jl)
       T = TC(Jj) + 273.16    !  Canopy or Vegetation Temperature in Kelvin
       T0 = T - 298.16        !  relative T to 25 degree Celcius, means T - 25
       T1 = 298.16            !  25 degree Celsius in Kelvin

!RG       IF (pvm(jl)==0.) pvm(jl)=pvmmin
       
       IF (PC4FLG(JL).EQ.0) THEN   ! C3

!---------------------------------------------------------------------------------
!   C3     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
!
! Rate (with Aktivationenergy) vegetation temperature dependance is
!  k = k(25C) * EXP((Veg Temp - 25) * aktivation energy
!                    / 298.16 * R * (Veg Temp + 273.16))
! => k = k0 * EXP( T0 * E / T1 / R / T ),  WHERE R is the gas constant (8.314)
! This holds for OX-Oxygen partial pressure, KC-Michaelis-Menten constant for CO2,
!    KO-Michaelis-Menten constant for O2, PVM-carboxylation capacity,
!    RD-Dark respiration, K-PEPcase CO2 specivity
!    Knorr (106)
!---------------------------------------------------------------------------------
          KC = KC0(JL) * EXP(EC(JL) * T0 / T1 / R / T)
          KO = KO0(JL) * EXP(EO(JL) * T0 / T1 / R / T)

!---------------------------------------------------------------------------------!
! CO2 compensation point without leaf respiration, Gamma* is assumed to be linearly 
! dependant on vegetation temperature, Gamma* = 1.7 * TC (IF Gamma* in microMol/Mol)
! Here, Gam in Mol/Mol,       Knorr (105)
!---------------------------------------------------------------------------------
          GAM = MAX (tgam(JL) * TC(Jj), 0.)

!---------------------------------------------------------------------------------
!
! PVM and PJM are not only temperature dependant but also differ inside the canopy. 
! This is due to the fact that the plant distributes its Nitrogen content and 
! therefore Rubisco content so that, the place with the most incoming light got the 
! most Rubisco. Therefore, it is assumed that the Rubisco content falls 
! exponentially inside the canopy. This is reflected directly in the values of PVM 
! and PJM at 25 Celsius (PVM * nscl),  Knorr (107/108)
!---------------------------------------------------------------------------------
          VCMAX = PVM(JL) * NSCL(JL) * EXP(EV(JL) * T0 / T1 / R / T)

!---------------------------------------------------------------------------------
!
! The temperature dependance of the electron transport capacity follows
! Farqhuar(1988) with a linear temperature dependance according to the vegetation 
! temperature
!  J = J(25C) * TC / 25 WHERE J(25) = J0 * NSCL
! ?????????????????????/ PJMAXMIN=1E-12 vielleicht etwas gering
!---------------------------------------------------------------------------------
          JMAX = PVM(JL)*PJMF(JL)*NSCL(JL)*TC(Jj)*jmt + pjmaxmin
!Mas09072010          JMAX = MAX(JMAX,PJMAXMIN)

!---------------------------------------------------------------------------------
!                                                                                
!  C3       GROSS PHOTOSYNTHESIS AT GIVEN CI                                     
!                                                                                
!---------------------------------------------------------------------------------
!c  The assimilation follows the Farqhuar (1980) formulation for C3 plants
!  A = min{JC, JE} - RD
!  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
!  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)      with
!   J = alpha * I * PJM / sqrt(PJM^2 + alpha^2 * I^2) with I=PAR in Mol(Photons)
!        Knorr (102a-c, 103)
!  Here J = J1 and A is the gross photosynthesis, i.e. still including the
!          respiratory part RD
!---------------------------------------------------------------------------------
!Mas09072010          IF ( JMAX .GT. PJMAXMIN) THEN
             J1 = ALPHA(JL) * PAR(JL) * JMAX &
                  &  / SQRT(JMAX**2 + (ALPHA(JL) * PAR(JL))**2)
!Mas09072010          ELSE
!Mas09072010             J1 = 0.
!Mas09072010          ENDIF

          JE = J1 * (CI(JL) - GAM) / 4. / (CI(JL) + 2. * GAM)
          JC = VCMAX * (CI(JL) - GAM) / ( CI(JL) + KC * (1. + OX / KO) )

!RG
          IF ((je+jc)**2.0-4.0*eta*je*jc .LT. 1.e-24) THEN
             A(JL) =  (JE + JC)/(2.0*eta)
          ELSE
             A(JL) =  (JE + JC - SQRT((je+jc)**2.0-4.0*eta*je*jc))/(2.0*eta)
          ENDIF
!RG             A(JL) =  (JE + JC - SQRT((je+jc)**2.0-4.0*eta*je*jc))/(2.0*eta)


          A(JL) = A(JL) * HITINHIB(TC(Jj))
             
!---------------------------------------------------------------------------------
!                                                                                
!  C3      COMPUTE 'DARK' RESPIRATION, PHOTORESPIRATION, STOMATAL CONDUCTANCE    
!                                                                                
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
! Following Farqhuar et al. (1980), the dark respiration at 25C is proportional to
! PVM at 25C, therefore RD = const * PVM, but the temperature dependance goes with
! ER (for respiration) and not with EV (for PVM)
!---------------------------------------------------------------------------------
          RD(JL) = FRDC3 * PVM(JL) * NSCL(JL) * EXP(ER(JL)*T0/T1/R/T) * &
               & HITINHIB(TC(Jj)) * DARKINHIB(PIRRIN(Jj))

!---------------------------------------------------------------------------------
! Photorespiration is anothter thing !!!!!!
!---------------------------------------------------------------------------------
          RPH(JL) = VCMAX * GAM / ( CI(JL) + KC * ( 1. + OX / KO ) ) * &
               & HITINHIB(TC(Jj)) 

!---------------------------------------------------------------------------------
! Diffusion equation Flux = (PCA - CI) / resistence, rs
!   conductance gs = 1 / rs  =>  Flux = (PCA-CI) * gs
!   Flux of CO2 is A * amount, Assimilation rate * amount
!   A is here, Gross Assimilation, though A-RD = (net) Assimilation rate
!   the amount comes from the ideal gas equation pV=nRT => n/V = p / RT
!   the stomatal conductance for CO2 is less THEN the conductance of H2O by
!   the factor of 1.6: gs(CO2) = gs(H2O) / 1.6, due to its lower mobiblity due
!   to its higher mass
!   => A (net) = gs/1.6 * (PCA-CI) * p/RT
!   => gs = A(net)*1.6*RT/p/(PCA-CI)
!---------------------------------------------------------------------------------
!          GS(JL) = MAX (1.6 * (A(JL) - RD(JL)) / &
!               & (PCA - CI(JL)) * R * T / P(Jj), PGSMIN)
          GS(JL) = 1.6 * (A(JL) - RD(JL)) /  (PCA - CI(JL)) * R * T / P(Jj)
          GS(JL) = MAXX (GS(JL), PGSMIN, PGSMIN*10)

!---------------------------------------------------------------------------------
       ELSE IF (pc4flg(jl) .EQ. 1) THEN           ! C3 -> C4 Plants
                 
!---------------------------------------------------------------------------------
!
!                               C 4                                              
!
!---------------------------------------------------------------------------------
!                                                                                
!   C4     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
!           AND 'DARK' RESPIRATION                                               
!                                                                                
!---------------------------------------------------------------------------------
! For C4 plants the Farquhar equations are replaced by the set of equations of
!  Collatz et al. 1992:
!  A = min{JC, JE} - RD
!  JC = k * CI
!  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]      with
!  Ji = alphai * Ipar / Epar with Ji=PAR in Mol(Photons)
!        Knorr (114a-d)
!  alphai is the integrated quantum efficiency for C4 plants (ALC4 = 0.04, 
!    compared to the efficiency of C3 plants, ALPHA = 0.28)
!  Theta is the curve PARAMETER (0.83) which makes the change between
!   PVM and K limitation smooth
!  K is the PEPcase CO2 specifity instead of the electron transport capacity
!   within C3 plants
!  Ci is the stomatal CO2 concentration = Cimin + (Ci0 - Cimin)* GC/GC0 with
!    Cimin = 0.3 and 0.15 CA respectivly
!
! The factor 1E3 comes that PJM for C3 is in microMol and K is in milliMol,
!   which is not considered in INITVEGDATA
! K scales of course with EK
!---------------------------------------------------------------------------------
          K = PJMF(JL) * 1.E3 * NSCL(JL) * EXP(EK(JL)*T0/T1/R/T)

!---------------------------------------------------------------------------------
! same as C3
!---------------------------------------------------------------------------------
          VCMAX = PVM(JL) * NSCL(JL) * EXP(EV(JL)*T0/T1/R/T)

!---------------------------------------------------------------------------------
!        same as C4, just the 25 degree Celsius proportional factor is different
!    0.011 for C3,  0.0042 for C4 
!---------------------------------------------------------------------------------
          RD(JL) = FRDC4 * PVM(JL) * NSCL(JL) * EXP(ER(jl)*T0/T1/R/T) * &
               & HITINHIB(TC(Jj)) * DARKINHIB(PIRRIN(Jj))

!---------------------------------------------------------------------------------
!
!  C4       GROSS PHOTOSYNTHESIS AT GIVEN CI                                     C
!
!---------------------------------------------------------------------------------
!  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]
!    Ji = ALC4 * PAR
!  J0 is the sum of the first two terms in JE
!---------------------------------------------------------------------------------
          J0 = (ALC4(JL) * PAR(JL) + VCMAX) /  2. / THETA

!---------------------------------------------------------------------------------
!
!  last 2 terms:  with J0^2 = 1/4/Theta^2*(PVM+Ji)^2
!       sqrt(1/4/Theta^2)*sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji))
!   = sqrt (J0^2 - PVM*Ji/Theta)
!---------------------------------------------------------------------------------
          JE = J0 - SQRT (J0**2 - VCMAX * ALC4(JL) * PAR(JL) / THETA)

!---------------------------------------------------------------------------------
!         see above
!---------------------------------------------------------------------------------
          JC = K * CI(JL)

!---------------------------------------------------------------------------------
! same as C3, Farquhar Assimilation
!---------------------------------------------------------------------------------
!RG
          IF ((je+jc)**2.0-4.0*eta*je*jc .LT. 1.e-24) THEN
             A(JL) = (JE + JC)/(2.0*eta)
          ELSE
             A(JL) = (JE + JC - SQRT((je+jc)**2.0-4.0*eta*je*jc))/(2.0*eta)
          ENDIF
!RG             A(JL) = (JE + JC - SQRT((je+jc)**2.0-4.0*eta*je*jc))/(2.0*eta) 

          A(JL) = A(JL) * HITINHIB(TC(Jj))

!---------------------------------------------------------------------------------
!                                                                                
!   C4     COMPUTE PHOTORESPIRATION, STOMATAL CONDUCTANCE                        
!
!---------------------------------------------------------------------------------
! same as C3 (diffusion equation)
!---------------------------------------------------------------------------------
          GS(JL) = 1.6 * (A(JL) - RD(JL)) / (PCA - CI(JL)) * R * T / P(Jj)
!          GS(JL) = MAX(GS(JL), PGSMIN)
          GS(JL) = MAXX (GS(JL), PGSMIN, PGSMIN*10)

!---------------------------------------------------------------------------------
! Photorespiration is 0 for C4 plants
!---------------------------------------------------------------------------------
          RPH(JL) = 0.
          
       ENDIF                                   !  C3 / C4  Plants
    END DO                                     !  Veglist
!---------------------------------------------------------------------------------
    
  END SUBROUTINE PHSYN1




!----------------------------------------------------------------------------
  SUBROUTINE PHSYN2 (ng,gridl,PAR, CI, GS, TC, P, NSCL, A, RD, RPH, PVM, PJMF, &
       &  PC4FLG, PJMAXMIN, PGSMIN, PCA, R, N, ILIST, &
       &  NLIST, PIRRIN,ts,EC,EO,EV,ER,EK,tgam,alpha,alc4,kc0,ko0,coszen)


!---------------------------------------------------------------------------- 
!    INPUT:                                                                      
!    PAR     ABSORBED PAR [MOL(PHOTONS) / M^2 S]                                 
!!    CI IF KFLG=0     CO2 CONCENTRATION INSIDE THE LEAF [MOL(CO2) / MOL(AIR)]    
!!    GS IF KFLG=1    STOMATAL CONDUCTANCE FOR WATER VAPOUR [M / S]               
!    TC      LEAF AND CANOPY TEMPERATURE [DEG C]                                 
!    P       AIR PRESSURE [PA]                                                   
!    NSCL    SCALING FACTOR FOR "PVM" AND "PJM"                                  
!    KFLG    FLAG FOR DIFFERENT COMPUTATION MODES                                
!            0: COMPUTE PHOTOS. & STOMATAL COND. AT GIVEN 'CI'                   
!            1: COMPUTE PHOTOS. & 'CI'  AT GIVEN  STOMATAL COND.                 
!    PVM     MAXIMUM RATE OF CARBOXYLATION [MOL(CO2) / M^2 S]                    
!    PJMF     Jmax/Vmax
!    PC4FLG  PHOTOSYNTHETIC PATHWAY                                              
!            0: C3; 1: C4                                                        
!    PJMAXMIN MINIMUM OF MAXIMUM CARBOXYLATION RATE               
!    PGSMIN   MINIMUM STOMATAL CONDUCTANCE                                       
!    PCA      ATMOSPHERIC CO2 MIXING RATIO [MOL(CO2)/MOL(AIR)]                  
!    R        GAS CONSTANT [J / K MOL]                                          
!    OUTPUT:                                                                    
!    A       GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2 S]                             
!    RD      DARK RESPIRATION OF LEAF [MOL(CO2) / M^2 S]                         
!    RPH     PHOTORESPIRATION [MOL(CO2) / M^2 S]                                 
!    GS IF KFLG=0    STOMATAL CONDUCTANCE FOR WATER VAPOUR [M / S]               
!    CI IF KFLG=1    CO2 CONCENTRATION INSIDE THE LEAF [MOL(CO2) / MOL(AIR)]     
!                                                                                
!    INTERNAL:                                                                   
!    PCA      CO2 CONCENTRATION OF AMBIENT AIR [MOL(CO2) / MOL(AIR)]             
!                                                                                
!    PHYSIOLOGY:                                                                 
!    C3 PLANTS: FARQHUAR, G.D., S. VON CAEMMERER AND J.A. BERRY, 1980.           
!               A BIOCHEMICAL MODEL OF PHOTOYNTHESIS IN LEAVES OF C3 SPECIES.    
!               PLANTA 149, 78-90.                                               
!    PVM      MAXIMUM CARBOXYLATION RATE AT 25C [MOL(CO2)/M^2 S]                 
!    PJM      MAXIMUM RATE OF ELECTRON TRANSPORT AT 25C [MOL(CO2)/M^2 S]         
!    VCMAX   MAXIMUM CARBOXYLATION RATE [MOL(CO2)/M^2 S]                         
!    JMAX    MAXIMUM RATE OF ELECTRON TRANSPORT [MOL(CO2)/M^2 S]                 
!    ALPHA   EFFICIENCY OF OF PHOTON CAPTURE                                     
!    OX      OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)]                           
!    KC0     MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)]      
!    KO0     MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)]        
!    GAM     COMPENSATION POINT WITHOUT DARK RESPIRATION [MOL(CO2) / MOL(AIR)]   
!    EC      ACTIVATION ENERGY FOR KC [J / MOL]                                  
!    EO      ACTIVATION ENERGY FOR KO [J / MOL]                                  
!    EV      ACTIVATION ENERGY FOR VCMAX [J / MOL]                               
!    ER      ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL]                    
!    AV      CONSTANT FOR HIGH-TEMPERATURE INHIBITION OF VCMAX [J / MOL]         
!    BV      CONSTANT FOR HIGH-TEMPERATURE INHIBITION OF VCMAX [J / MOL K]       
!    FRDC3   RATIO OF DARK RESPIRATION TO "PVM" AT 25C                           
!    C4 PLANTS: COLLATZ, G.J., M. RIBAS-CARBO AND J.A. BERRY, 1992.              
!               COUPLED PHOTOSYNTHESIS-STOMATAL CONDUCTANCE MODEL FOR LEAVES     
!               OF C4 PLANTS. AUST. J. PLANT PHYSIOL. 19, 519-538.               
!    ALC4    EFFECTIVE QUANTUM EFFICIENCY                                        
!    K       CO2 SPECIFICITY OF PEPCASE                                          
!    THETA   CURVATURE PARAMETER                                                 
!                                                                                
!----------------------------------------------------------------------------      


    INTEGER :: n, ts,ng
    REAL    PAR(N), CI(N), GS(N), TC(N)
    REAL    P(ng), NSCL(N)
    REAL    A(N), RD(N), RPH(N)
    REAL    PVM(N), PJMF(N)
    REAL    PJMAXMIN, PGSMIN, PCA, R
    INTEGER  pc4flg(n)
    INTEGER ILIST(N), NLIST,gridl(n)
    REAL    PIRRIN(ng), coszen(ng)
    REAL, DIMENSION(n) :: EC, EO, EV, ER, EK, tgam
    REAL, DIMENSION(n) :: alpha, alc4, kc0, ko0

    REAL    K
    REAL    VCMAX, KC, KO, GAM, JE, JC, JMAX, J0, J1
    REAL    T, T0, T1
    REAL    G0, K1, W1, K2, W2, B, C
    INTEGER JL, JJ
    


!-----------------------------------------------------------------------
! dummy assignments for TAF
!-----------------------------------------------------------------------
    A(1) = A(1)
    RD(1) = RD(1)

!---------------------------------------------------------------------------------
!                                                                                
!                  SECOND ENTRY, TC,GC => AC, CI                                 
!
!---------------------------------------------------------------------------------
!$taf loop = parallel
!    print *,'            DO LOOP over vp: ', NLIST
    DO Jl = 1, NLIST
       Jj = gridl(Jl)
       T = TC(JL) + 273.16             !  Canopy or Vegetation Temperature K
       T0 = T - 298.16                 !  relative T to 25 degree C, T - 25
       T1 = 298.16                     !  25 degree Celsius in Kelvin

!RG       IF (pvm(jl)==0.) pvm(jl)=pvmmin


       IF (PC4FLG(JL) .EQ. 0.) THEN    !  C3

!---------------------------------------------------------------------------------
!                                C 3                                             
!                                                                                
!   C3     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
!           And 'DARK' RESPIRATION                                               
!
!---------------------------------------------------------------------------------
          KC = KC0(JL) * EXP(EC(JL) * T0 / T1 / R / T)
          KO = KO0(JL) * EXP(EO(JL) * T0 / T1 / R / T)

!             if (ts==30) print*,'phy: ',jl,tc(jl)
!---------------------------------------------------------------------------------
!  Tc is calculated in main PROGRAM before
!---------------------------------------------------------------------------------
          GAM = MAX (tgam(JL) * TC(JL), 0.)        ! Gam=1.7*TC
          VCMAX = PVM(JL) * NSCL(JL) &
               * EXP(EV(JL) * T0 / T1 / R / T)     ! PVM=PVM0*NSCL*EXP
                                                   ! (Energie scaled relative 25C)
          RD(JL) = FRDC3 * PVM(JL) * NSCL(JL) &    ! RD=RD0*NSCL*EXP
                                                   ! (Energie scaled relative 25C)
               * EXP(ER(JL) * T0 / T1 / R / T)  &  ! with RD0 proportional PVM0 
                                                   ! at 25C
               * HITINHIB(TC(JL)) &
               * DARKINHIB(pirrin(jj))
          JMAX = PVM(JL)*PJMF(JL) * NSCL(JL) &     ! PJM=PJM0*NSCL*EXP
                                                   ! (Energie scaled relative 25C)
               * TC(JL) * jmt + pjmaxmin
!Mas09072010          JMAX = MAX(JMAX, PJMAXMIN)

!---------------------------------------------------------------------------------
!
!  C3       GROSS PHOTOSYNTHESIS AT GIVEN TC                                     
!
!---------------------------------------------------------------------------------
!  Remember:
!  A = min{JC, JE} - RD
!  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
!  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)      with
!   J = alpha * I * PJM / sqrt(PJM^2 + alpha^2 * I^2) with I=PAR in Mol(Photons)
!        Knorr (102a-c, 103)
! J = J1
!---------------------------------------------------------------------------------
!Mas09072010          IF ( JMAX .GT. PJMAXMIN ) THEN
             J1 = ALPHA(JL) * PAR(JL) * JMAX &
                  / SQRT(JMAX**2 + (ALPHA(JL) * PAR(JL))**2)
!Mas09072010          ELSE
!Mas09072010             J1 = 0.
!Mas09072010          ENDIF

!---------------------------------------------------------------------------------
!         Helping friends K1, W1, W2, K2
!---------------------------------------------------------------------------------
          K1 = 2. * GAM
          W1 = J1 / 4.
          W2 = VCMAX
          K2 = KC * (1. + OX / KO)
          
!---------------------------------------------------------------------------------
! A = gs / 1.6 * (PCA - Ci) * p / R / T
! <=> Ci = PCA - 1.6 * R * T / p / gs * A = PCA - A / G0
!---------------------------------------------------------------------------------
          G0 = GS(JL) / 1.6 / R / T * P(jj)

!---------------------------------------------------------------------------------
! A = min{JC, JE} - RD
! => A = JC - RD ^ A = JE - RD
! Set this (A =) in Ci formula above
! Set Ci in
!  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)
! => quadratic formula in JE
! 0 = JE^2 -(RD+J/4+G0*(PCA+2*GAMMA))*JE +J/4*G0*(PCA-GAMMA)+J/4*RD
!---------------------------------------------------------------------------------

          B = RD(JL) + W1 + G0 * (PCA + K1)
          C = W1 * G0 * (PCA - GAM) + W1 * RD(JL)

!---------------------------------------------------------------------------------
! with 0 = x^2 + bx + c
!      x1/2 = -b/2 +/- sqrt(b^2/4 - c)
!  take '-' as minimum value of formula
!---------------------------------------------------------------------------------
!RG
          IF (B**2 / 4. - C .LT. 1.e-24) THEN
             JE = B / 2.
          ELSE
             JE = B / 2. - SQRT ( MAX (B**2 / 4. - C, 0.))
          ENDIF
!RG             JE = B / 2. - SQRT ( MAX (B**2 / 4. - C, 0.))

!---------------------------------------------------------------------------------
! Set Ci in
!  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
! WRITE JC = PVM * (Ci - Gam) / (Ci + K2)
! => quadratic formula in JC
! 0 = JC^2 -(RD+PVM+G0*(PCA+K2))*JC +PVM*G0*(PCA-GAMMA)+RD*PVM
!---------------------------------------------------------------------------------
          B = RD(JL) + W2 + G0 * (PCA + K2)
          C = W2 * G0 * (PCA - GAM) + W2 * RD(JL)
!RG
          IF (B**2 / 4. - C .LT. 1.e-24) THEN
             JC = B / 2.
          ELSE
             JC = B / 2. - SQRT ( MAX(B**2 / 4. - C, 0.))
          ENDIF
!RG             JC = B / 2. - SQRT ( MAX(B**2 / 4. - C, 0.))

!---------------------------------------------------------------------------------
! A = min{JC, JE} - RD
!  but here A = Gross Photosynthesis = GPP
!  => A = min{JC, JE}
!---------------------------------------------------------------------------------
!RG
          IF ((je+jc)**2.0-4.0*eta*je*jc .LT. 1.e-24) THEN
             A(JL) = (JE + JC)/(2.0*eta)
          ELSE
             A(JL) = (JE + JC - SQRT((je+jc)**2.0-4.0*eta*je*jc))/(2.0*eta)
          ENDIF
!RG             A(JL) = (JE + JC - SQRT((je+jc)**2.0-4.0*eta*je*jc))/(2.0*eta)

          A(JL) = A(JL) * HITINHIB(TC(JL))
   
!---------------------------------------------------------------------------------
!
!   C3     COMPUTE PHOTORESPIRATION AND INTER STOMATE CO2 CONCENTRATION          
!
!---------------------------------------------------------------------------------
! A = gs / 1.6 * (PCA - Ci) * p / R / T
! <=> Ci = PCA - 1.6 * R * T / p / gs * A = PCA - A / G0
!   with A = Net assimilation = NPP
!---------------------------------------------------------------------------------
          CI(JL) = PCA - MAX((A(JL)-RD(JL)) / MAX(G0, 1.E-6), 0.)

!---------------------------------------------------------------------------------
! Photorespiration 
! Carboxylilation controlled assimilation
! JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO)) =  PVM * (Ci - Gam) / (Ci + K2)
! JC = Assimilation - Photorespiration
! JC = PVM * (Ci - Gam) / (Ci + K2)
! Photorespiration = PVM * Gam / (Ci + K2)
          RPH(JL) = VCMAX * GAM / ( CI(JL) + KC * ( 1. + OX / KO )) * &
               & HITINHIB(TC(JL)) 

!---------------------------------------------------------------------------------

       ELSE                          ! C3 -> C4 Plants

!---------------------------------------------------------------------------------
!
!                               C 4                                              
!
!                                                                                
!   C4     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT          
!           AND 'DARK' RESPIRATION                                               
!
!---------------------------------------------------------------------------------
! Remind:
!  Collatz et al. 1992:
!  A = min{JC, JE} - RD
!  JC = k * Ci
!  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]      with
!   Ji = alphai * Ipar / Epar with Ji=PAR in Mol(Photons)           and
!   PAR = Ipar / Epar;  ALC4=alphai; J0=1/2/Theta *(PVM + Ji);
!   Ci = PCA - 1.6 * R * T / p / gs * A = PCA - A / G0                and
!   => A = JC - RD ^ A = JE - RD
!---------------------------------------------------------------------------------
          K = PJMF(JL) * 1.E3 * NSCL(JL) &          ! K=K0*NSCL*EXP
                                                    ! (Energie scaled relative 25C)
               * EXP(EK(JL) * T0 / T1 / R / T)      ! in mmol CO2 PepCase for C4
          VCMAX = PVM(JL) * NSCL(JL) &              ! PVM=PVM0*NSCL*EXP
                                                    ! (Energie scaled relative 25C)
               * EXP(EV(JL) * T0 / T1 / R / T)
          RD(JL) = FRDC4 * PVM(JL) * NSCL(JL) &     ! RD=RD0*NSCL*EXP
                                                    ! (Energie scaled relative 25C)
               * EXP(ER(JL) * T0 / T1 / R / T)   &  ! with RD0 proportional PVM0 
                                                    ! at 25C
               * HITINHIB(TC(JL)) &
               * DARKINHIB(pirrin(jj))

!---------------------------------------------------------------------------------
!
!  C4       GROSS PHOTOSYNTHESIS AT GIVEN CI                                     
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!  Ci = PCA - 1.6 * R * T / p / gs * A = PCA - A / G0    
!--------------------------------------------------------------------------------- 
          G0 = GS(JL) / 1.6 / R / T * P(jj)

!---------------------------------------------------------------------------------
!  J0=1/2/Theta *(PVM + Ji) = (alphai * PAR + PVM) / 2 / Theta
!---------------------------------------------------------------------------------
          J0 = (ALC4(JL) * PAR(JL) + VCMAX) /  2. / THETA

!---------------------------------------------------------------------------------
!  JE = J0 - sqrt( J0^2 - PVM*alphai*PAR/Theta)
!---------------------------------------------------------------------------------
          JE = J0 - SQRT (J0**2 - VCMAX * ALC4(JL) * PAR(JL) / THETA)

!---------------------------------------------------------------------------------
!  JC = (G0*PCA + RD) / (1 + G0/K)
!---------------------------------------------------------------------------------
          JC = (G0 * PCA + RD(JL)) / (1. + G0 / K) 

!---------------------------------------------------------------------------------
!  A = min{JC, JE} - RD, but here: A = GPP
!---------------------------------------------------------------------------------
          A(JL) = (JE + JC - SQRT((je+jc)**2.0-4.0*eta*je*jc))/(2.0*eta)
          A(JL) = A(JL) * HITINHIB(TC(JL))

!---------------------------------------------------------------------------------
!
!   C4     COMPUTE PHOTORESPIRATION AND INTER STOMATE CO2 CONCENTRATION          
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!  Ci = PCA - 1.6 * R * T / p / gs * A = PCA - A / G0, with A = NPP
!---------------------------------------------------------------------------------
          CI(JL) = PCA - MAX((A(JL)-RD(JL)) / MAX(G0, 1.E-6), 0.)

!---------------------------------------------------------------------------------
! No Photorespiration with C4 plants
!---------------------------------------------------------------------------------
          RPH(JL) = 0.

       ENDIF                                      !  C3 Plants = 1
    END DO                                        !  IVEGLIST
!!    ENDIF                                           !  IPHLG=0

    RETURN
  END SUBROUTINE PHSYN2



!-----------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------
  SUBROUTINE FAPARL (ng,gridl,N,nnl,PLAI,FC,FCMAX,RHOS,COSZEN,FDIR,LAIL, &
       & APAR,OMEGA,ASPE,DDL,ZENMIN,ZFMIN,ZLMIn,iLIST,nLIST)

!-----------------------------------------------------------------------
! ON ENTRY:
!   PLAI   : Total LAI
!   FC     : Canopy fraction
!   FCMAX  : Maximum canopy fraction (0.9)
!   RHOS   : soil albedo
!   COSZEN : cosine of solar zenith angle, mue
!   FDIR   : direct part of total Irradiance (e.g. 0.8)
!   OMEGA  : single leaf scattering albedo = 0.12
!   ASPE   : aspect ration (diameter/height) = 0
!   DDL    : leaf layer borders in units of LAI (0, 1/3, 2/3, 1)
!   ZENMIN : minimum of zenit winkel for which there is an PAR calculated (1e-3)
!   ZFMIN  : minimum of vegetation coverage
!   ZLMIN  : minimum LAI or LAI per layer
!   NNL    : # of canopy layers
!   N      : # of Longitude, DIMENSION OF THE IN AND OUTPUT FIELDS
!   ILIST  : Vegetation grid point list
!   NLIST  : # of grid points including vegetation on this latitude 
! ON EXIT:
!   LAIL   : LAI per canopy layer
!   APAR   : THE ABSORBED PAR PER LEAF AERA PER CANOPY LAYER
!            NORMALIZED TO INCOMING RADIATION = 1
!------------------------------------------------------------------------

    use mo_taf
    USE mo_constants, ONLY: pi


    INTEGER :: n, nnl,ng
    REAL      PLAI(N),  FC(N), FCMAX(N), RHOS(N) 
    REAL      COSZEN(ng), FDIR(ng)
    REAL      LAIL(N,NNL), APAR(N,NNL)
    REAL      OMEGA(ng) , ASPE(ng) , DDL(0:NNL)
    REAL      ZENMIN , ZFMIN , ZLMIN
    INTEGER   ILIST(N), NLIST, gridl(n)
    
!---------------------------------------------------------------------
! Local Variables
!---------------------------------------------------------------------
    REAL      ZH(N), K0(N), ZP0, ZP1(N)
    REAL      Q0(N), Q1(N), B0(N), B1(N)
    REAL      EKL(N), EHL(N), EKL0, EHL0, B4(N)
    REAL      X0, X1, X2, F, FC0, FC1, SHD, FSHD(N)
    INTEGER   JL, JJ, K
    
!-----------------------------------------------------------------------
! dummy assignments for TAF
!-----------------------------------------------------------------------
    LAIL(1,1) = LAIL(1,1)
    APAR(1,1) = APAR(1,1)

!$TAF STORE fdir  = carbon_tape, key = keydiurnal
!----------------------------------------------------------------------------
!                                                                    
!     COMPUTE ABSORBED PAR PER LEAF AREA (!) AND PLAI PARTITIONING   
!                                                                    
!----------------------------------------------------------------------------
! Initialize ABSORBED PAR PER LEAF AREA,  IF COSZEN(JL)<1e-3
! Spread LAI equally around the three layers
!----------------------------------------------------------------------------
    DO K = 1, NNL
       DO JJ = 1, NLIST
          JL = ILIST(JJ)
          APAR(JL,K) = 0.

!----------------------------------------------------------------------------
! LAIL=1/3*PLAI
!----------------------------------------------------------------------------
          LAIL(JL,K) = PLAI(JL) * (DDL(K) - DDL(K-1))
          LAIL(JL,K) = MAX(LAIL(JL,K), ZLMIN)
       ENDDO
    ENDDO

    b0 = 0.
    b1 = 0.
    b4 = 0.
    q0 = 0.
    q1 = 0.
    ekl = 0.
    ehl = 0.
    zh = 0.
    k0 = 0.
    zp1 = 0.
    fshd = 0.

!--------------------------------------------------------------------
! The Apsorbed Par per laef Area which is used later for the Net Assimilation,
! is calculated via the two stream approximation of Sellers (1985):
!  muq * d(Rdown)/dl + (1-(1-b)*omega)*Rdown - omega*b*Rup   = omega*muq*k*(1-b0)*R
! -muq * d(Rup)/dl   + (1-(1-b)*omega)*Rup   - omega*b*Rdown = omega*muq*k*b0*R
!  with
!   Rdown - downwards diffusive flux
!   Rup   - upwards diffusive Flux
!   R     - direct flux, with R(0) = dPAR * RPAR with RPAR incoming irridiance in PAR
!           and R = R0 * EXP(-kl) with R0=R(0) - exponential extinction law of 
!           Monsi and Racki (1953)
!           or similar Lambert-Beer's law
!   b     - forward scatter fraction of diffusive flux
!   b0    - forward scatter fraction of direct flux
!   k     - extinction coefficient for the direct flux
!   muq = int(1/K(mu))dmu|_0^1 - the integral of 1/K over the downward hemisphere
!   omega - the single leaf albedo
!
!  The general solutions are (kl,hl=k*l,h*l):
!  Rup   = q2*R0*EXP(-kl) + p1*B1*EXP(hl) + p2*B2*EXP(-hl)
!  Rdown =-q1*R0*EXP(-kl) +    B1*EXP(hl) +    B2*EXP(-hl)
!   with
!    h  = sqrt( (1-(1-b)*omega)^2/muq^2 - omega^2*b^2/muq^2 )
!    p1 = ( (1-(1-b)*omega) + muq * h )/omega/b
!    p2 = ( (1-(1-b)*omega) - muq * h )/omega/b
!   -q1 = (omega^2*b*muq*k* (1-b0) + omega*    b0  *muq*k*(1-(1-b)*omega - muq*k))/
!         ((1-(1-b)*omega)^2-muq^2*k^2-omega^2*b^2)
!    q2 = (omega^2*b*muq*k*    b0  + omega* (1-b0) *muq*k*(1-(1-b)*omega + muq*k))/
!         ((1-(1-b)*omega)^2-muq^2*k^2-omega^2*b^2)
!    B1/B2 from boundary conditions
!-------------------------------------------------------------------
!  Make two assumptions:
!  1) the distribution of leaf angles is isotropic
!  2) the leaf reflectivity and transmissivity are equal (the sum = omega)
!  => b=0.5, b0=0.5, k=0.5/mu with mu=cos(theta) solar zenith angle => muq=1
!
!  => k  = 1/2/mu
!     h  = sqrt( 1 - omega )
!     p1 = ( 1-omega/2 + h )/omega/2
!     p2 = ( 1-omega/2 - h )/omega/2
!   ! p2 = 1 / p1 !
!     q1 = ( (1 + 2*mu)*omega/2 )/(1-4*mu^2*(1-omega)) = ( k*(k + 1)*omega/2 )/
!                                                        (k^2-1-omega)
!     q2 = ( (1 - 2*mu)*omega/2 )/(1-4*mu^2*(1-omega)) = ( k*(k - 1)*omega/2 )/
!                                                        (k^2-1-omega)
!    
! Determine B1 and B2 from the boundary conditions:
!  1) Rdown(0) equals the incoming diffuse radiation
!     Rdown(0) = (1-dPAR)*RPAR
!  => Rdown(0) + R(0) = (1-dPAR)*RPAR + dPAR*RPAR = RPAR as total incoming PAR
!  2) the reflection at the lower boundary of the canopy is given by the soil 
!     reflectance
!     Rup(LAI) = RhosPAR * (R(LAI) + Rdown(LAI))
!  Here: FAPARL gets RhosPAR as Variable RHOS, LAI is the total canopy LAI, PLAI
! 
!  => B1 = + ( eta*R0 - (Rd+q1*R0) * gamma1 )/(gamma1 - gamma2)
!     B2 = - ( eta*R0 - (Rd+q1*R0) * gamma2 )/(gamma1 - gamma2)
!   with
!     eta    = RhosPAR * (1-q1)-q2) * EXP(-k*LAI)
!     gamma1 = ( p1 - RhosPAR) * EXP( + h*LAI)
!     gamma2 = ( p2 - RhosPAR) * EXP( - h*LAI)
!     Rd     = Rdown(0) = (1-dPAR)*RPAR
!------------------------------------------------------------------
! THAT IS THE COMPLETE SOLUTION OF THE TWO STREAM APPROXIMATION UNDER THE BOUNDARY
! CONDITIONS AND ASSUMPTIONS MENTIONED ABOVE !!!!!!!!!!!!!!!!!!
!
! Therefore, the absorbed Radiation inside the canopy is:
!   APAR = -d/dl ( R(l) + Rdown - Rup)
!        = (1-q1-q2)*k*R0*EXP(-kl) - (1-p1)*h*B1*EXP(hl) + (1-p2)*h*B2*EXP(-hl)
! But the absorbed PAR per canopy layer in discrete steps is:
!   APAR(JL) = 1/(D(JL)-D(JL-1)) * int(-d/dl(R(l) + Rdown - Rup))dl|_D(JL)^D(JL-1)
!            = (R(D(JL-1)) + Rdown(D(JL-1)) - Rup(D(JL-1)) - R(D(JL)) + Rdown(D(JL))
!              - Rup(D(JL)) / ((D(JL)-D(JL-1))
!  and  R(l)+Rdown-Rup = (1-q1-q2)*  R0*EXP(-kl) + (1-p1)*  B1*EXP(hl) + (1-p2) * 
!                        B2*EXP(-hl)
!------------------------------------------------------------------
! The clumping of the vegetation is taken into account, defining LAIc = LAI / fc
! as an effective LAI.
! Taken this into account, l = l/fc but the solutions stay as they are because 
! of the differentiations are take d/dl according to the NEW l=l/fc
! Only APAR has to be multiplied with fc at the END because D(JL)-D(JL-1) is still
! the old l
!------------------------------------------------------------------
!$taf loop = parallel
    DO Jl = 1, NLIST
       Jj = gridl(Jl)
       
       IF (COSZEN(Jj).GE.ZENMIN) THEN

!----------------------------------
!c h = sqrt( 1 - omega )
!----------------------------------
          ZH(JL) = SQRT (1. - OMEGA(Jj))

!----------------------------------
!c p1 = ( 1-omega/2 + h )/omega/
!----------------------------------
          ZP1(JL) = (1. - OMEGA(Jj) / 2. + ZH(JL)) / OMEGA(Jj) * 2. 

!----------------------------------
! k = 0.5/mu
!----------------------------------
          K0(JL) = 0.5 / COSZEN(Jj)
          IF (K0(JL).EQ.ZH(JL)) K0(JL) = K0(JL) + 1E-12
          IF (K0(JL).EQ.-ZH(JL)) K0(JL) = K0(JL) + 1E-12
!--------------------------------------------------------------
! denominator of q1 and q2
!--------------------------------------------------------------
          X0 = (1. - 4. * COSZEN(Jj)**2 * ZH(JL)**2)
!--------------------------------------------------------------
! q1 = ( (1 + 2*mu)*omega/2 )/(1-4*mu^2*(1-omega))
!    = ( k*(k + 1)*omega/2 )/(k^2-1-omega)
!--------------------------------------------------------------
          Q1(JL) = ((1. + 2. * COSZEN(Jj)) * OMEGA(Jj) / 2.) / X0

!--------------------------------------------------------------
! q2 = ( (1 - 2*mu)*omega/2 )/(1-4*mu^2*(1-omega))
!    = ( k*(k - 1)*omega/2 )/(k^2-1-omega)
!--------------------------------------------------------------
          Q0(JL) = ((1. - 2. * COSZEN(Jj)) * OMEGA(Jj) / 2.) / X0
          FC0 = MAX (FC(JL), ZFMIN)
          FC1 = MAX (FCMAX(JL), ZFMIN)

!---------------------------------------------------------------
! IF the vegetation clumps, the lower layers of the canopy get more direct,
! i.e. overall sunlight
! Not used yet, i.e. ASPE=0 => SHD=0 => FSHD = fc, just the simple effective LAI
!---------------------------------------------------------------
          SHD = MIN( 4. / PI * ASPE(Jj) &
               * SQRT( 1. / COSZEN(Jj)**2 - 1. ), 10. )
          FSHD(JL) = FC1 - (FC1 - FC0) * EXP(-FC0 * SHD / FC1)
       else
          fshd(JL) = 0.
          q0(JL) = 0.
          q1(JL) = 0.
          k0(JL) = 0.
          zh(JL) = 0.
          zp1(JL) = 0.
       ENDIF  ! mue>0.001
    END DO ! VEGLIST

!FastOpt !$TAF STORE fshd,k0,zh  = carbon_tape, key = keydiurnal
!$taf loop = parallel
    DO Jl = 1, NLIST
       Jj = gridl(Jl)
       
       IF (COSZEN(Jj).GE.ZENMIN) THEN
!-----------------------------------------------------------
! EXP(-k*LAI/fc)
!-----------------------------------------------------------
          EKL(JL) = EXP(-K0(JL) / FSHD(JL) * PLAI(JL))
              
!------------------------------------------------------------
! EXP(-h*LAI/fc)
!-----------------------------------------------------------
          EHL(JL) = EXP(-ZH(JL) / FSHD(JL) * PLAI(JL))

       else
          ehl(JL) = 0.
          ekl(JL) = 0.
       ENDIF  ! mue>0.001
    END DO ! VEGLIST

!FastOpt !$TAF STORE ehl,q0,q1  = carbon_tape, key = keydiurnal
!$taf loop = parallel
    DO Jl = 1, NLIST
       Jj = gridl(Jl)
       IF (COSZEN(Jj).GE.ZENMIN) THEN
!----------------------------------
! p2 = ( 1-omega/2 - h )/omega/2 = 1 / p1
!----------------------------------
          ZP0 = 1. / ZP1(JL)

!-----------------------------------------------------------
! gamma1 = ( p1 - RhosPAR) * EXP( + h*LAI)
!-----------------------------------------------------------

          X1 = (ZP1(JL) - RHOS(JL)) / EHL(JL)

!-----------------------------------------------------------
! gamma2 = ( p2 - RhosPAR) * EXP( - h*LAI)
!-----------------------------------------------------------          
          X0 = (ZP0 - RHOS(JL)) * EHL(JL)

!-----------------------------------------------------------
! eta = RhosPAR * (1-q1)-q2) * EXP(-k*LAI)
!-----------------------------------------------------------
          X2 = (RHOS(JL) * (1. - Q1(JL)) - Q0(JL)) * EKL(JL)

!------------------------------------------------------------
! F = 1 - dPAR + dPAR * q1
! => F*RPAR = Rd + q1*R0
! i.e. calculation takes RPAR=1
!-----------------------------------------------------------           
          F = 1. - FDIR(Jj) + Q1(JL) * FDIR(Jj)

!-------------------------------------------------------------
! B1(JL)*RPAR = B1, B2(JL)*RPAR = B2, B4(JL)*RPAR = R(0) + Rdown(0) - Rup(0)
!  B1 = + ( eta*R0 - (Rd+q1*R0) * gamma1 )/(gamma1 - gamma2)
!------------------------------------------------------------
          B1(JL) = (X2 * FDIR(Jj) - F * X0) / (X1 - X0)

!-----------------------------------------------------------
!  B2 = - ( eta*R0 - (Rd+q1*R0) * gamma2 )/(gamma1 - gamma2)
!-----------------------------------------------------------
          B0(JL) = (X2 * FDIR(Jj) - F * X1) / (X0 - X1)

       else
          b1(jl) = 0.
          b0(jl) = 0.
       ENDIF  ! mue>0.001
    END DO ! VEGLIST

!$taf loop = parallel
    DO Jl = 1, NLIST
       Jj = gridl(Jl)
       IF (COSZEN(Jj).GE.ZENMIN) THEN
!----------------------------------
! p2 = ( 1-omega/2 - h )/omega/2 = 1 / p1
!----------------------------------
          ZP0 = 1. / ZP1(JL)
!-----------------------------------------------------------
!  R(l)+Rdown-Rup = (1-q1-q2)* R0*EXP(-kl) + (1-p1)* B1*EXP(hl) + (1-p2)* B2*EXP(-hl)
!----------------------------------------------------------------------------
          B4(JL) = (1. - Q0(JL) - Q1(JL)) * FDIR(Jj) &
               + (1. - ZP1(JL)) * B1(JL) + (1. - ZP0) * B0(JL)
       else
          B4(JL) = 0.
       ENDIF  ! mue>0.001
    END DO ! VEGLIST

!FastOpt !$TAF STORE q0,q1,k0,b0,b1,fshd  = carbon_tape, key = keydiurnal
!------------------------------------------------------------
! Only Layers-1, i.e. except the lower boundary calculation
!----------------------------------------------------------------------------
!FastOpt !$taf loop = parallel
    DO K = 1, NNL - 1
!$taf loop = parallel
       DO Jl = 1, NLIST
          Jj = gridl(Jl)
          
          IF (COSZEN(Jj).GE.ZENMIN) THEN

!----------------------------------------------------------------------------
! p2=1/p1
!----------------------------------------------------------------------------
             ZP0 = 1. / ZP1(JL)
                 
!----------------------------------------------------------------------------
! EXP(-k*l/fc)
! with l=LAI*DDL(K), i.e. l is element of [0,LAI], i.e. 0,LAI/3,2LAI/3,LAI
!----------------------------------------------------------------------------
             EKL0 =  EXP(-K0(JL) / FSHD(JL) * DDL(K) * PLAI(JL))

!!----------------------------------------------------------------------------
! EXP(-h*l/fc)
!----------------------------------------------------------------------------
             EHL0 = EXP(-ZH(JL) / FSHD(JL) * DDL(K) * PLAI(JL))

!!----------------------------------------------------------------------------
! R(D(JL))+ Rdown(D(JL))- Rup(D(JL))=
!      (1-q1-q2)*R0*EXP(-kl)+(1-p1)*B1*EXP(hl)+(1-p2)*B2*EXP(-hl)
!  i.e. X0*RPAR = above
!----------------------------------------------------------------------------
             X0 = (1. - Q0(JL) - Q1(JL)) * EKL0 * FDIR(Jj) &
                  + (1. - ZP1(JL)) * B1(JL) / EHL0 &
                  + (1. - ZP0) * B0(JL) * EHL0

!----------------------------------------------------------------
! APAR(JL) = (R(D(JL-1))+Rdown(D(JL-1))-Rup(D(JL-1)) - R(D(JL))+
! Rdown(D(JL))-Rup(D(JL)) / ((D(JL)-D(JL-1))
! Here APAR only the nominator; the division is made outside FAPARL
!---------------------------------------------------------------
             APAR(JL,K) = B4(JL) - X0

!!---------------------------------------------------------------
! Partition evenly LAI, LAIL=LAI/3
!             LAIL(JL,K) = (DDL(K) - DDL(K-1)) * PLAI(JL)
! R(D(JL-1))+Rdown(D(JL-1))-Rup(D(JL-1) in next step
!---------------------------------------------------------------
             B4(JL) = X0
          else
             apar(jl,k) = 0.
             B4(JL) = 0.
          ENDIF ! mue>0.001
       ENDDO   ! VEGLIST
    ENDDO     ! # canopy layers - 1

!---------------------------------------------------------------
! Now the same for the last layer
!---------------------------------------------------------------
!FastOpt !$TAF STORE q0,q1,ekl,ehl,b0,b1  = carbon_tape, key = keydiurnal
    DO Jl = 1, NLIST
       Jj = gridl(Jl)
       IF (COSZEN(Jj).GE.ZENMIN) THEN

!----------------------------------------------------------------------------
! p2=1/p1
!----------------------------------------------------------------------------
             ZP0 = 1. / ZP1(JL)
!-------------------------------------------------------------
! R(D(NL))+ Rdown(D(NL))- Rup(D(NL))=
!      (1-q1-q2)*R0*EXP(-kLAI)+(1-p1)*B1*EXP(hLAI)+(1-p2)*B2*EXP(-hLAI)
!  i.e. X0*RPAR = above for the lower boundary
!------------------------------------------------------------
          X0 = (1. - Q0(JL) - Q1(JL)) * EKL(JL) * FDIR(Jj) &
               + (1. - ZP1(JL)) * B1(JL) / EHL(JL) &
               + (1. - ZP0) * B0(JL) * EHL(JL)

!------------------------------------------------------------!
! APAR(NL) = (R(D(NL-1))+Rdown(D(NL-1))-Rup(D(NL-1)) - R(D(NL))+
! Rdown(D(NL))-Rup(D(NL)) / ((D(NL)-D(NL-1))
! Here APAR only the nominator; the division is made outside FAPARL
!------------------------------------------------------------
          APAR(JL,NNL) = B4(JL) - X0

!------------------------------------------------------------
!             LAIL(JL,NNL) = (DDL(NNL) - DDL(NNL-1)) * PLAI(JL)
!------------------------------------------------------------
       else
          apar(jl,nnl) = 0.
       ENDIF ! mue>0.001
    END DO   ! VEGLIST

!FastOpt !$TAF STORE apar  = carbon_tape, key = keydiurnal
!$taf loop = parallel
    DO K = 1, NNL
!$taf loop = parallel
       DO Jl = 1, NLIST
          Jj = gridl(Jl)

!!------------------------------------------------------------
! APAR multiplied with fc because D(JL)-D(JL-1) 
! division is made outside FAPARL and is still the old l
!------------------------------------------------------------
          IF (COSZEN(Jj).GE.ZENMIN) APAR(JL,K)=APAR(JL,K)*FSHD(JL)
       ENDDO
    ENDDO
    
  END SUBROUTINE FAPARL


!----------------------------------------------------------------------------







!------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                           @
!                  Nitrogen scaling within canopy layers (nl)           @
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                           @
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!------------------------------------------------------------------------

  SUBROUTINE NSCALE (ng, gridl, NSCL, PLAI, DDL, KLON, nLIST, &
       & iLIST, nl, hclass, pheno, spds, cpds, lailim,laimax)


    USE mo_helper

! 1. --------------------------------------------------------------------
! NSCL, Nitrogen scaling factor per canopy layer per longitude
! PLAI, Leaf area INDEX per longitude (!!!not ZLAIL=lai per canopy layer)
! LAT, #Latitude
!------------------------------------------------------------------------        

    INTEGER :: klon, nl, ng
    REAL ::   NSCL(KLON,nl), PLAI(KLON), spds(ng), cpds(ng)
    REAL, dimension(0:nl) ::   DDL
    REAL ::   lailim, laimax, xlai
    INTEGER :: nlist, pheno(klon), hclass(klon)
    INTEGER :: ILIST(KLON), gridl(klon)

!------------------------------------------------------------------------
! Internal Variables 
!------------------------------------------------------------------------
    REAL    K12(ng)
    INTEGER JL, IL, JJ

!-----------------------------------------------------------------------
! dummy assignment for TAF
!-----------------------------------------------------------------------
    NSCL(1,1) = NSCL(1,1)

!-----------------------------------------------------------------------
!     LEAF-NITROGEN SCALING
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  K = 1 / 2mue     Knorr (119c)
!-----------------------------------------------------------------------
    k12 = 0.5 / MAX(spds+cpds,1e-3)
    
    DO IL = 1, NL
       DO Jl = 1, NLIST
          Jj = gridl(Jl)
!------------------------------------------------------------------------
!  IF LAI>LAI0  and (Vegetation either Tree or Shrub  OR  Irrigated or arable crop 
!    (0.1 for numerical reasons: often (3=3)/2=2.99...)
!-------------------------------------------------------------------------
!MS06082010          IF (PLAI(JL).GE.(LAILIM).AND.(HCLASS(JL).LE.2. &
!MS06082010               .OR.PHENO(JL).GE.6.)) THEN
          IF ((HCLASS(JL).LE.2. &
               .OR.PHENO(JL).GE.6.)) THEN

!-------------------------------------------------------------------------
!             IF (HCLASS(JL).LE.2) THEN
!  nscl = EXP(- K12 * l) WHERE l runs from 0 to LAI, DDL=0,1/3,2/3,1
!   (DDL(IL-1)+DDL(IL))/2. is the middle of each canopy layer (1/6, 1/2, 5/6)
!               NSCL(JL,IL) =
!      =           EXP (-K12*(DDL(IL-1)+DDL(IL))/2.*PLAI(JL))
!  Take Top of Canopy layer instead of middle
!-------------------------------------------------------------------------
!!MS$ 30.11.02    NSCL(JL,IL) =  EXP( -K12(jj) * DDL(IL-1) * PLAI(JL) )
!MS06082010             nscl(jl,il)=exp( -k12(jj)*ddl(il-1)*plai(jl) )* &
!MS06082010                  & (plai(jl)-lailim) / (laimax-lailim)+ &
!MS06082010                  & (laimax-plai(jl)) / (laimax-lailim)

             xlai = maxx(plai(jl),lailim,0.1)
             nscl(jl,il)=exp( -k12(jj)*ddl(il-1)*plai(jl) )* &
                  & (xlai-lailim) / (laimax-lailim)+ &
                  & (laimax-xlai) / (laimax-lailim)

          ELSE
             NSCL(JL,IL) = 1.
          ENDIF
       ENDDO
    ENDDO
    
    RETURN
        
  END SUBROUTINE NSCALE


     

!-----------------------------------------
  REAL FUNCTION HITINHIB(T)
    REAL T
!-----------------------------------------
! T : Canopy temperature in Degree CELSIUS
! FUNCTION which inhibits Assimilation and Respiration at temperatures above 
! 55 Celsius from Collatz et al., Physiological and environmental regulation
! of stomatal conductance, photosynthesis and transpiration: a model that 
! includes a laminar boundary layer,
! Agricultural and Forest Meteorology, 54, pp. 107-136, 1991
!------------------------------------------
    HITINHIB = 1. / ( 1. + EXP( 1.3 * ( T - 55. ) ) )

  END FUNCTION HITINHIB
!-----------------------------------------

!-----------------------------------------
  REAL FUNCTION COLDINHIB(T)
    REAL T
!------------------------------------------

!WOK-2009-11-10 Replaced step function by logistic function with range -0.5...+0.5C   
!    COLDINHIB =  SIGN(0.5,t) + 0.5
    COLDINHIB =  1./(1.+exp(-t*10))

  END FUNCTION COLDINHIB
!-----------------------------------------


!-----------------------------------------      
  REAL FUNCTION DARKINHIB(IRR)
    REAL IRR
!----------------------------------------------
! IRR : Total irridiance at the surface in mol/m^2s
!
! FUNCTION which inhibits Dark-Respiration in light
!  after Brooks and Farquhar, Effect of temperature on the CO2/O2 specifity on RBisCO
!  and the rate of respiration in the light, Planta 165, 397-406, 1985
!
! It is fitted to inhibit the dark-respiration to 50% of it's uninhibited value 
!   up from 50 umol/m^2s.
! The FUNCTION is an own creation.
!---------------------------------------------
    IF ( IRR .LT. 0 ) THEN
       DARKINHIB = 0.
    ELSE
       DARKINHIB = 0.5 + EXP( -1. * ALOG( 2. ) - IRR * 1.e6 / 10. )
    ENDIF

  END FUNCTION DARKINHIB
!-----------------------------------------------------

END MODULE mo_carbon
