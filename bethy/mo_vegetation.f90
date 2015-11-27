MODULE mo_vegetation

!CCCC vegetation ini parameters 



  IMPLICIT NONE

! .. Variables

  INTEGER, ALLOCATABLE, DIMENSION (:) :: c4flg, class, ph, pft
! WOK-ADD-070629 sla
  REAL, ALLOCATABLE, DIMENSION(:) ::  vm, jmf, fautleaf, ccost, hv, sla
  REAL, ALLOCATABLE,DIMENSION(:) :: aw, q10f, q10s, fracs, tauf,EC, EO, EV, ER, EK, tgam
  REAL, ALLOCATABLE,DIMENSION(:) :: alpha,alc4,kc0,ko0
  REAL, ALLOCATABLE, DIMENSION(:) :: Chl

!WOK-ADD-070612 ztrans
!WOK-CHG-070626 removed zctrans
!WOK-CHG-070626 removed zctrans
  REAL, ALLOCATABLE, DIMENSION(:) :: zfpar, zgc, ptv, zfc, zlai, lai
  REAL, ALLOCATABLE, DIMENSION(:) :: zassc, zraut, zrphc,zmaint,zgrowth
  REAL, ALLOCATABLE, DIMENSION(:) :: ztrans, zptrans, zpcevp, zpsevp
  REAL, ALLOCATABLE, DIMENSION(:,:) :: pgs
  
  REAL, PARAMETER :: frtdh = 0.3 ! herbaceous rooting depth as fraction of tree rooting depth

! hv         vegetation height [m]
! c4flg      photosynthetic pathway flag
! jm         Max electron transport in [Mol(CO2)/m^2/s]      
! vm         Max Carboxylation rate in [Mol(CO2)/m^2/s]
! class      height class


END MODULE mo_vegetation
