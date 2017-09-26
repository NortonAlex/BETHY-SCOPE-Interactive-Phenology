!------------------------------------------------------------------
!------------------------------------------------------------------
! module of array indices for control vector mapping
! also sets number of physical model variables used in mapping 'nvar'
!------------------------------------------------------------------
MODULE mo_ctrl

  IMPLICIT NONE

  INTEGER, PARAMETER :: ivm       =  1
  INTEGER, PARAMETER :: ijmf      =  2
  INTEGER, PARAMETER :: ifautleaf =  3
  INTEGER, PARAMETER :: iccost    =  4
  INTEGER, PARAMETER :: iq10f     =  5
  INTEGER, PARAMETER :: iq10s     =  6
  INTEGER, PARAMETER :: itauf     =  7
  INTEGER, PARAMETER :: iaw       =  8
  INTEGER, PARAMETER :: ifracs    =  9
  INTEGER, PARAMETER :: ier       = 10
  INTEGER, PARAMETER :: iev       = 11
  INTEGER, PARAMETER :: ieo       = 12
  INTEGER, PARAMETER :: iec       = 13
  INTEGER, PARAMETER :: iek       = 14
  INTEGER, PARAMETER :: ialpha    = 15
  INTEGER, PARAMETER :: ialc4     = 16
  INTEGER, PARAMETER :: ikc0      = 17
  INTEGER, PARAMETER :: iko0      = 18
  INTEGER, PARAMETER :: itgam     = 19
  INTEGER, PARAMETER :: ibeta     = 20
  INTEGER, PARAMETER :: iplaimax  = 21
  INTEGER, PARAMETER :: iptphen   = 22
  INTEGER, PARAMETER :: iptphenr  = 23
  INTEGER, PARAMETER :: ipdphen   = 24
  INTEGER, PARAMETER :: ipdphenr  = 25
  INTEGER, PARAMETER :: iplgr     = 26
  INTEGER, PARAMETER :: ipkl      = 27
  INTEGER, PARAMETER :: iptauw    = 28
  INTEGER, PARAMETER :: iChl      = 29
  INTEGER, PARAMETER :: iCdm      = 30
  INTEGER, PARAMETER :: iCsm      = 31
  INTEGER, PARAMETER :: iLIDFa      = 32
  INTEGER, PARAMETER :: iLIDFb      = 33
  INTEGER, PARAMETER :: ihc       = 34
  INTEGER, PARAMETER :: ileafwidth       = 35
  INTEGER, PARAMETER :: ivms      = 36
  INTEGER, PARAMETER :: ikc0s     = 37
  INTEGER, PARAMETER :: iko0s     = 38
  INTEGER, PARAMETER :: ivomf     = 39
  INTEGER, PARAMETER :: irdf      = 40

  INTEGER, PARAMETER :: npar      = 40

  REAL  :: dummy_in
  REAL  :: dummy_out


END MODULE mo_ctrl

