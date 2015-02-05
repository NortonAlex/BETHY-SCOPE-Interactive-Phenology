!***********************************************************
!*   SUBROUTINE surfcae
!*   does prelimanary energy and water balance
!***********************************************************

SUBROUTINE surface (ng,vp,temp,nc,u,p,psoilst,lintw,zfpar,zfc,zlai,psradsd,zgc, &
                  & pgs,ptv,ea0,fe,ztrans,zptrans,zpcevp,zpsevp,ts,hv,day1p,zrhos,keydiurnal)


! .. Use Statements ..
  USE mo_constants
!HEW-ADD-050307:
  USE mo_namelist, ONLY : dtime
  USE mo_carparams
  USE mo_surface
  USE mo_grid, ONLY : gridp
  USE mo_helper, ONLY: mins, maxx

  IMPLICIT NONE
  
  INTEGER, INTENT(in) :: ng, vp
  REAL, DIMENSION(ng), INTENT(in) :: u, p, temp, nc, psradsd, ea0
  REAL, DIMENSION(vp), INTENT(in) :: psoilst, lintw, hv, zrhos,fe
!WOK-ADD-070626 zpevp, i.e. potential evaporation
  REAL, DIMENSION(vp), INTENT(out) :: ptv,ztrans,zptrans,zpcevp,zpsevp
  REAL, DIMENSION(vp), INTENT(inout) :: zgc, zfpar, zfc, zlai
  REAL, DIMENSION(vp,nl), INTENT(inout) :: pgs
  INTEGER, INTENT(in) :: day1p, keydiurnal

! .. Local Scalars ..

  INTEGER :: nveglist,ts
  INTEGER :: jl, jj
  REAL :: zfacbe, zpsoilst, eta, zgcold
  REAL :: ztk, zta, zrhoa, zea, zs, zrln
  REAL :: zlt, ztwodt, zxrn, zlcontent, zrnc, zrns
  REAL :: zunfrozen, roots, rootd, rootc
  REAL :: zlambda, zgamma, zav, zas, zsfl, ztlv
  REAL :: zfx, zfthw, a0, b0, esta, afx
  REAL, DIMENSION(vp) :: td3m1m, td4m1m, td5m1m, zga
  REAL, DIMENSION(vp) :: zdemand, rbeold
  REAL, DIMENSION(vp) :: zsupp, zltrans, zde


  nveglist = vp
  ztwodt = 2*dtime*60 ! in seconds
  eta=0.99

!$TAF store zgc, rbe  = diurnal_tape, rec = keydiurnal

  DO jl =1,nveglist
     jj=gridp(jl)
     rbeold(jl)=rbe(jl)
! .. do some additional calculations for bethy
     ! aerodynamic conductance
     zga(jl) = karmann**2*u(jj)/alog(href/(rz0*MAX(hv(jl),1E-3))+az0)**2    
  ENDDO
  

     
!------------------------------------------------------------------------
!   ADJUSTMENT OF STOMATAL CONTROL PARAMETER "RBE"                              
!    1. Set some used constants                                                 
!    2. Get Air temperature ZTK and density ZRHOA,                              
!       actual water vapour pressure ZEA                                        
!    3. CALCULATE AERODYNAMIC CONDUCTANCE "RGA" (KARMANN FORMULA)               
!    4. Soil water supply rate = "CW0"(MAX.ROOT SUPPLY RATE)                    
!                                *"PSOILST"(REL.PLANT AVAILABLE SOIL WATER)     
!    5. Calculate vegetation temperature and evapotranspiration                 
!      5.1 CALCULATE VAPOUR PRESSURE DEFICIT "ZDE" with water vapour pressure   
!           at AIR TEMPERATURE                                                  
!      5.2 Radiation balance at the surface                                     
!      5.3 ADJUST "RBE"                                                         
!         5.3.1 Trade off between demand of plants = non-water limited          
!               Photosynthesis and supply rate                                  
!         5.3.2 Adjust BE to the water limitation                               
!      5.4 CALCULATE CANOPY TEMPERATURE "PTV" with energy budget                
!      5.4 Adjust Canopy Conductance "ZGC" with control PARAMETER "RBE"         
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!   factor for opening (widening) the stomata per timestep
!------------------------------------------------------------------------
  ZFACBE = 0.98

!------------------------------------------------------------------------
!   Actual calculation over latitudes  C
!------------------------------------------------------------------------
!$taf loop = parallel
  DO Jl = 1, NVEGLIST              ! 1) to 2.4)
     jj=gridp(jl)
!------------------------------------------------------------------------
!  2.)   Air Temp, Air Density,               C
!        Air water vapour pressure            C
!
!------------------------------------------------------------------------
!    Air Kelvin Temperature full level, Air Celsius
!------------------------------------------------------------------------
     ZTK = temp(Jj) + tmelt
     ZTA = temp(jj)

!------------------------------------------------------------------------
!
!    air density = Mass*Pressure / (R*T)  (unit: kg/m3)
!    from pV=nRT
!         pM/rho=nRT
!------------------------------------------------------------------------
     ZRHOA = amd /1000. / ar / ZTK * P(Jj) 
     
     zfthw = SIGN(0.5,zta) + 0.5
     a0 = 17.269 * zfthw + 22.33 * (1.-zfthw)
     b0 = 237.3 * zfthw + 271.15 * (1.-zfthw)
     esta = 610.78*EXP(a0*zta/(b0+zta))
     zs = a0*b0*esta/(b0+zta)**2
     zea = ea0(jj) + reld0*fe(jl)*(esta-ea0(jj))
     zea = MIN(zea,esta)
     zde(jL) = esta - zea     

!------------------------------------------------------------------------
!
!   4.)  Soil supply rate                          
!
!------------------------------------------------------------------------
     ZPSOILST = PSOILST(JL)
     
!------------------------------------------------------------------------
! the supply is the rate cw0 (0.5 mm/hour) * an effective Soil water content / maximal soil water
!   cw0 is a free PARAMETER
!------------------------------------------------------------------------
     ZSUPP(JL) = CW0 / 3600. * ZPSOILST
!     ZSUPP(JL) = MAX( ZSUPP(JL), ZWATERMIN)
!     ZSUPP(JL) = MAXX( ZSUPP(JL), ZWATERMIN, ZWATERMIN*1E3)
             
!------------------------------------------------------------------------
!
!   5.) Calculating Evaporation and Leaf temperature                     
!
!------------------------------------------------------------------------
!
!   5.1) VAPOUR PRESSURE DEFICIT  C
!
!------------------------------------------------------------------------
!    latent evaporation heat, degree C
!     i.e. (2.501e6 - 2.38e3 * TV) IF TV>0 Celcius, 2.834e6 IF TV<0 Celcius
!      lambda = (2.501e6 - 2.38e3 * ta) * fthw + 2.834e6 * (1.-fthw)
!     here RLAM0 = 3.1512 because TV in Kelvin, i.e.
!      (3.1512e6 - 2.38e3 * 273.16)=2.501e6
!------------------------------------------------------------------------
     ZFTHW = SIGN( 0.5, ZTA ) + 0.5
     ZLAMBDA = ( RLAM0 + RLAM1 * ZTK ) * ZFTHW + RLAM2 &
              &  * ( 1. - ZFTHW )

!------------------------------------------------------------------------
!   psychometric constant for Penman-Monteith formula
!------------------------------------------------------------------------
     ZGAMMA = P(Jj) * cpd / ZEPS / ZLAMBDA


     zde(jL) = MAX(zde(jL),0.)

!------------------------------------------------------------------------
!   5.2) Radiation balance at surface              C
!      not depending on vegetation temperature     C
!      because thermal equilibrium is expected     C
!
!------------------------------------------------------------------------
!  Net thermal Radiatian = rad downwards - rad upwards
!    Net = sigma * Tkelvin^4 * (emissivity of the cloudless atmosphere, epsA
!           * cloudiness correction, rEps - surface emissivity, eps0)
!          Knorr, 1997, (40)-(39), S.34 (32e)
!         epsA = epsA0 (0.64) * (vapour pressure/Tkelvin)^(1/7)   Knorr (41)
!         rEps = 1 + 0.22 * cloud fraction^2                     Knorr (42)
!         eps0 = 0.97                                            Knorr S.35 (33e)
!         PACLCVM, cloud fraction t-dt (Total Cloud Cover:   US151189.228)
!------------------------------------------------------------------------
!WOK-CHG-070712 Formulation straightened in accordance with full BETHY 'dayphtrans'
     ZRLN = stbo * ZTK**4 &
          & * (( 1. + 0.22 * nc(Jj)**2 ) * ALBS0 &
          & * ( MAX(ZEA , 0.) / ZTK )**0.143 &
          & - 0.97)
     
!------------------------------------------------------------------------
!     -         - 0.97 * SIG * ( ZFC(JL) * (PTV(JL)+TMELT)**4
!     +                         + ( 1 - ZFC(JL) ) * PTSM1M(JL)**4 )
!------------------------------------------------------------------------
! Absorptivity of vegetation, aV
!    aV = (1 - albedo of dense vegetation, albV0 (0.15) - fraction absorbed by soil
!          under dense vegetation, aSoil0 (0.05) ) * fraction of absorbed PAR
!       Knorr (46)
!------------------------------------------------------------------------
     ZAV = ( 1. - ALBV0 - ASOIL0 ) * ZFPAR(JL)

!------------------------------------------------------------------------
! Absorptivity of soil, aS
!   aS = 1 - soil albedo, rS - ( 1 - rS - aSoil0) FPAR
!         Knorr (45)
!------------------------------------------------------------------------
     ZAS = 1. - ZRHOS(JL) - ( 1. - ZRHOS(JL) - ASOIL0 ) * ZFPAR(JL)

!------------------------------------------------------------------------
! Soil heat flux, G
!   G = 0.036 * total net Radiation, Rn                            Knorr (38)
!     Rn = thermal rad down - thermal rad up + 
!           (1 - surface albedo, albS) * total solar irridiance, Rs  Knorr (36)
!    albS = 1 - aV - aS => 1-albS=aV+aS                            Knorr (47)
!------------------------------------------------------------------------
     ZSFL = 0.036 * ( ZRLN + ( ZAV + ZAS ) * PSRADSD(Jj) )

!------------------------------------------------------------------------
! Transmission of radiation through canopy according to Beer's Law of radiation
!        absorption, tLV
! tLV = (1 - fractional vegetation cover, fc) +
!       fc * EXP( - LAI / fc)
!------------------------------------------------------------------------

!WOK-CHG-070712 corrected for uniform leaf-angle distribution
!     ZTLV =  (1. - zfc(jl)) + ZFC(JL) * EXP( -0.5 * ZLAI(JL)/ MAX(ZFC(JL),1e-3) )
!     ZTLV =  (1. - zfc(jl)) + ZFC(JL) * EXP( -1.0 * ZLAI(JL)/ MAX(ZFC(JL),1e-3) )
!WOK-CHG-081217 inconsistency in FC removed
!     ZTLV =  (1. - zfc(jl)) + MAX (ZFC(JL), ZFCMIN) * EXP( -0.5 * ZLAI(JL) / MAX (ZFC(JL), ZFCMIN) )
     ZFX = MAXX (ZFC(JL), ZFCMIN, ZFCMIN*10)
     ZTLV =  (1. - zfc(jl)) +  ZFX * EXP( -0.5 * ZLAI(JL) / ZFX )

!------------------------------------------------------------------------
! Vegetation part of total net radiation, Rnv
!  Rnv = (1 - tLV) * (rad down - rad up - G) + aV * Rs    Knorr (37a)
!  This is the Radiation which goes into the calculation of the canopy temperature
!------------------------------------------------------------------------

     ZRNC = ( 1.-ZTLV ) * ( ZRLN - ZSFL ) + ZAV * PSRADSD(Jj)

!------------------------------------------------------------------------
!   5.3) ADJUSTMENT OF PARAMETER "RBE"             C
!
!------------------------------------------------------------------------
!
!   5.3.1) Vegetation Demand of water              C
!        and trade off between Demand and Supply   C
!
!------------------------------------------------------------------------
! Denominator in Penman-Monteith formula, XRN
!  XRN = slope of water vapour pressure curve, zs * Rnv +
!         air density, rhoA * specific heat of air, cp * water vapour deficit, de
!          * aerodynamic conductance, GA
!------------------------------------------------------------------------
     ZXRN = ZS * ZRNC + ZRHOA * cpd * zde(jL) * ZGA(jl)
 
!------------------------------------------------------------------------
!         ZXRN = ZRHOA * CP * ZDE(JL) * ZGA(JL)
!------------------------------------------------------------------------
! PENMAN MONTEITH GIVES LAMBDA * EVAPOTRANSPIRATIONRATE
!  ZPTRANS = XRN / (zs + psychometric constant, gamma * 
!            (1 + GA / Not-water limited Canopy conductance GC0)) / evaporation heat, lambda
!                        Knorr (81)
! Because GC is still GC0 with no water limitation, ZPTRANS is the potential
!   transpiration rate, called the 'demand', i.e., IF there wouldn't be water limitation
!   the stomata conductance would go with the not-water limited photosynthesis rate, which is
!   limited by other factors (FARQUHAR-MODEL); 
!   demand means that this is the transpiration that the plant
!   would like to have, but there is only a special supply, ZSUPP of water out of the roots,
!   i.e. the plant has to CLOSE its stomate to make the demand smaller and demand=supply
!------------------------------------------------------------------------
     ZPTRANS(JL) = ZXRN / ( ZS + ZGAMMA &
                 & * ( 1. + ZGA(JL) / ZGC(JL) ) ) &
                 & / ZLAMBDA
        
!------------------------------------------------------------------------
!         ZPTRANS(JL) = ZXRN / ( ZGAMMA &
!             * ( 1. + ZGA(JL) / ZGC(JL) ) ) &
!             / ZLAMBDA
!------------------------------------------------------------------------
!     ZPTRANS(JL) = MAX( ZPTRANS(JL), ZWATERMIN)
!     ZPTRANS(JL) = MAXX( ZPTRANS(JL), ZWATERMIN, ZWATERMIN*1E3)

!WOK-ADD-070626 ZRNS & ZPSEVP from 'dayphtrans.f': 'rns' and 'pevap'
     ZRNS = ZTLV * ( ZRLN - ZSFL ) + ZAS * PSRADSD(JJ)
     ZPSEVP(JL) = ZS * ZRNS / (ZS + ZGAMMA) / ZLAMBDA
!     ZPSEVP(JL) = MAX( ZPSEVP(JL), ZWATERMIN)
     ZPSEVP(JL) = MAXX( ZPSEVP(JL), ZWATERMIN, ZWATERMIN*1E3)

!------------------------------------------------------------------------
! Potential canopy evaporation at the maximum evporation rate
!   which means in Penman-Monteith GC -> infinity
!------------------------------------------------------------------------
     ZPCEVP(JL) = ZXRN / (ZS + ZGAMMA ) / ZLAMBDA
!     ZPCEVP(JL) = MAX( ZPCEVP(JL), ZWATERMIN)

!------------------------------------------------------------------------
! IF de<0 THEN ea > es and there is a lot of water in the air, though that there is no
!  evatranspiration in the air out of the plant (not accounted for is interceptional water)
!
! The content of water in the interception reservoir in Echam is given in m equivalent water depth
!  which gives the water in kg/m^2 WHEN multiplied by the density of water = 1000 kg/m^3, AE -> kg: AE * Rho
!  => ZLCONTENT is the content of interception water in kg/m^2
! Therefore, ZLT is the time in which the content of water in the skin reservoir is evaporated
!  it is restricted to the length of the time step (ZTWODT)
! The actual evaporation is for the time ZTL with the rate for the skin reservoir and for the rest of
!  time with the Stomata evaporation rate
!
! Content of Interception water (1000 = water density)
!------------------------------------------------------------------------
!     ZLCONTENT = 0.     ! MS set to zero as we don't have interception here, 10.1.2002

!------------------------------------------------------------------------
! Time to evaporate the interception water at the transpiration rate 
! calculated with Penman-Monteith with Gc->infinity
!------------------------------------------------------------------------
! WOK-CHG-070612 replaced by new input field 'lintw'
!     ZLT = LINTW(JL) / ZLTRANS(JL)
!     ZLT = ZLCONTENT / ZLTRANS(JL)
!     ZLT = MAX( MIN( ZLT, ZTWODT ), 0. )
          
!------------------------------------------------------------------------
!  Demand through the stomata
!------------------------------------------------------------------------
     ZDEMAND(JL) = ZPTRANS(JL)

!------------------------------------------------------------------------
!  Actual Evaporation through the stomata is limited by the supply from the soil water
!------------------------------------------------------------------------
! WOK-CHG-07712 revert to minimum function, otherwise photosynthesis persists at 0 soil moisture
!     ZTRANS(JL) = (ZDEMAND(JL)+ZSUPP(JL)-SQRT((ZDEMAND(JL)+ZSUPP(JL))**2- &
!                &  4.0*eta*ZDEMAND(JL)*ZSUPP(JL)))/(2.*eta)
!     ZTRANS(JL) = Min( Zdemand(JL)  , Zsupp(jL) )
!if (jl==1) write (9,'(9E16.8)') Zdemand(JL), Zsupp(jL), ZPSOILST
! WOK-CHG-090402 another try, not believing the above
     ZTRANS(JL) = Mins( Zdemand(JL)  , Zsupp(jL), eta )

     ZTRANS(JL) = MAXX( ZTRANS(JL)  , ZWATERMIN, ZWATERMIN*1E3 )
!     ZTRANS(JL) = MAX( ZTRANS(JL)  , ZWATERMIN)

!------------------------------------------------------------------------
!  Overall Evaporation, i.e. the evapotranspiration through the stomata and the evaporation of the 
!   interception reservoir
!------------------------------------------------------------------------
!     ZCTRANS(JL) = ZTRANS(JL) + ZLTRANS(JL) * ZLT / ZTWODT
   
!------------------------------------------------------------------------
!
!   5.3.2) GC now from inverting penman/monteith with actual transpiration
!
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! 1.  fx = (ZXRN / (lambda * ZTRANS) - s) 
! IF demand < supply, THEN ZTRANS=ZPTRANS and ZTRANS
!  is the potential evatranspirationRATE and fx = gamma * (1 + GA / GC0)
! IF demand > supply, ZTRANS=ZSUPP, that's all water the soil can give
!  and ZTRANS is the actual evatranspirationRATE, fx = gamma * (1 + GA / GC)
!  (see GC instead of GC0)
!         ZFX =  ZXRN / ZLAMBDA / ( ZTRANS(JL) 
!     +                            + ZLTRANS(JL) * ZLT / ZTWODT ) - ZS
!------------------------------------------------------------------------
!     ZFX =  ZXRN / ZLAMBDA / ZTRANS(JL) - ZS
     ZFX =  ZXRN / ZLAMBDA / (ZTRANS(JL) + 1e-12) - ZS
 
!------------------------------------------------------------------------
! ZFX =  ZXRN / ZLAMBDA / ZCTRANS(JL)
!------------------------------------------------------------------------
! 2.  zgc = zga / (zfx/zgamm -1)
!------------------------------------------------------------------------
     zfx=(zfx/zgamma -1.) / ZGA(JL) * ZGC(JL)

!------------------------------------------------------------------------
! 3.  fx = (% - 1) / de = be
!------------------------------------------------------------------------
     ZFX = ( ZFX - 1. ) / (Zde(JL) + 1.0) ! NOTE: ZDE cannot be negative
!     ZFX = MAX( ZFX, 0. )
     ZFX = MAXX( ZFX, 0., 1E-6 )
!     IF ( Zde(JL) .GT. 1. ) THEN
!        ZFX = ( ZFX - 1. ) / Zde(JL)
!        ZFX = MAX( ZFX, 0. )
!     ELSE
!        ZFX = 0.
!     ENDIF
         

!------------------------------------------------------------------------
!
! 4.  Adjustment of be over the day with factor facbe (~0.85)
!     IF there is a water limitation, the stomata closes. It THEN starts
!     openeing again with the factor facbe. This interrupts sudden openings
!       and closings of the stomata.
!------------------------------------------------------------------------
!!     RBE(JL) = MAX( RBEold(JL) *zfacbe, ZFX )
     IF (day1p==0) rbe(jl)=zfx
!     IF (day1p==0) write (9,*) zfx


!------------------------------------------------------------------------
!
!     5.5)  Adjustment of "ZGC", Canopy Cond. C
!
!------------------------------------------------------------------------
! Actual canopy conductance GC = GC0 / ( 1 + be * de)    Knorr (82)
!------------------------------------------------------------------------
     ZGC(JL)  = ZGC(JL) / ( 1. + RBE(JL) * Zde(JL) )


!------------------------------------------------------------------------
!
!     5.4)  Adjustment of "PTV", Canopy Temp. C
!
!------------------------------------------------------------------------ 
! TV, CANOPY TEMPERATURE, ENERGY BALANCE, KNORR, 1997, (110), S.53 (5OE)
! TV = T + (NET RADIATION - REAL (actual) EVAPORATION)/
!              AIR DENSITY/AIR SPECIFIC HEAT/CANOPY CONDUCTANCE
!------------------------------------------------------------------------

!WOK-ADD-070712
!------------------------------------------------------------------------ 
! RECALCULATION OF TRANSPIRATION RATE AT GIVEN STOMATAL PARAMETER
!------------------------------------------------------------------------ 
     ZTRANS(JL) = ZXRN / ( ZS + ZGAMMA &
                 & * ( 1. + ZGA(JL) / ZGC(JL) ) ) &
                 & / ZLAMBDA

!WOK-CHG-070712 canopy temperature balance needs to take into account the effective area fraction
!     PTV(JL) = ZTA +( zrnc - ZLAMBDA * ZTRANS(JL))  &
!             & / ZRHOA / cpd / ZGA(JL)
!     ZFX = MAX (ZAV, 1E-3) / MAX (1. - ZTLV, 1E-3)
     ZFX = MAXX (ZAV, 1E-6, 1E-3) / MAXX (1. - ZTLV, 1E-6, 1E-3)
     PTV(JL) = ZTA +( ZRLN - ZSFL + ZFX * PSRADSD(Jj) - ZLAMBDA * ZTRANS(JL))  &
              & / ZRHOA / cpd / ZGA(JL)


  ENDDO

 
!$TAF store zde, rbe  = diurnal_tape, rec = keydiurnal
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!  ADJUSTMENT OF STOMATAL CONDUCTANCE FOR EACH CANOPY LAYER "ZGS",  C
!   "ZGC" IS "ZGS" INTEGRATED OVER THE LEAF AREA INDEX.             C
!------------------------------------------------------------------------
  DO jj = 1, NL
     DO Jl = 1, NVEGLIST
        
!------------------------------------------------------------------------
! stomatal conductance, gs = gs0 / (1+ be * de)
!------------------------------------------------------------------------
  
      PGS(JL,jj) = PGS(JL,jj) / ( 1. + RBE(JL) * Zde(JL) )
     END DO
  END DO
  

END SUBROUTINE surface
