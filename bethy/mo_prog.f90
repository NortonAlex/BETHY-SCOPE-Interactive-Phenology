module mo_prog
  ! passes prognostic info for Jacobian runs
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE     :: prog_global
  REAL, DIMENSION(:,:,:), ALLOCATABLE     :: prog_sites
end module mo_prog
