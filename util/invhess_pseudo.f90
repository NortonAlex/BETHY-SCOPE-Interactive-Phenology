PROGRAM inv_hess_pseudo
  ! inverts a symmetric matrix and prints out roots of diagonals as well as saving matrix
  ! really a wrapper for symmetric inverse routines from lapack
  ! author Peter Rayner
  ! date:  Sep 2001
  ! modified, Aug 2004, Marko Scholze
  ! modified, Oct 2003, HEW
  
USE mo_namelist

  IMPLICIT NONE 
  
  INTEGER, PARAMETER :: n_jac = 288
  
  REAL, ALLOCATABLE, DIMENSION( :, :) :: hess, hessin, ht, jt, a, b, e ,jac, c_pseudo
  REAL, ALLOCATABLE, DIMENSION( :, :) :: jac1, jac2, c_pseudo1, c_pseudo2
  REAL, DIMENSION(:), ALLOCATABLE :: evals
  REAL, DIMENSION(1000000) ::  work
  REAL, ALLOCATABLE, DIMENSION(:) :: x, adx, p0, p, p0sl, p0su, xsf, lb, ub
  REAL :: fc
  INTEGER, ALLOCATABLE, DIMENSION(:) :: pflag
  INTEGER, DIMENSION(1000000) :: lwork
  INTEGER n_vars, i, j, info, k ,l, cutoff
  INTEGER, EXTERNAL :: IARGC
  CHARACTER :: header, w
  CHARACTER (LEN= 20) :: filename, progName

! Determine name of OPWARMD file
   j=IARGC()
   IF(j.GT.1) THEN
      CALL getarg(0,progName)       ! determine absolute program name, i.e. including path)
      PRINT*,"Usage: ",TRIM(progName)," [ <opwarmd file> ] "
      STOP
   END IF
   IF(j.EQ.0) THEN
      filename="OPWARMD"       ! Default name of opwarmd file
   ELSE
      CALL getarg(1,filename)  ! Name of opwarmd file
   END IF

! .. read namelist variables
  CALL get_namelist

! .. read parameter file
  OPEN(1,file=param_file,status='old')
  READ(1,*)
  READ(1,*) n_vars
!HEW-CHG-050105: now p0sl(n_vars) and p0su(n_vars) instead p0s readin, but only p0s=p0sl=p0sl for pflag=0 needed
  ALLOCATE(x(n_vars), adx(n_vars), p0(n_vars), p(n_vars), p0sl(n_vars), p0su(n_vars), xsf(n_vars), &
       & pflag(n_vars), lb(n_vars), ub(n_vars))
  DO i=1,n_vars
     READ(1,*) p0(i),p0sl(i),p0su(i),xsf(i), pflag(i), lb(i), ub(i)
  END DO
  CLOSE(1)


  WRITE(*,*) filename
  OPEN(2,file=filename,form='unformatted',access='direct',recl=8*n_vars)
  READ(2,rec=1) x
  READ(2,rec=2) adx
  CLOSE(2)

  ALLOCATE( hess(n_vars, n_vars), a(n_vars, n_vars), b(n_vars, n_vars), e(n_vars, n_vars))
  ALLOCATE( evals(n_vars), jt(n_vars,n_vars), ht(n_vars,n_vars), hessin(n_vars,n_vars))
  ALLOCATE( jac1(n_vars, n_jac),jac2(n_vars, n_jac), c_pseudo1(n_jac,n_jac),c_pseudo2(n_jac,n_jac))

  jt=0.
  ht=0.

  OPEN(unit=1,file='output/hesscol.bin',form='unformatted')
  READ (1) hessin
  CLOSE(1)

  OPEN(unit=1,file='output/nep_jac_57x288.bin',form='unformatted')
  READ (1) jac1
  CLOSE(1)

  OPEN(unit=1,file='output/aet_jac_57x288.bin',form='unformatted')
  READ (1) jac2
  CLOSE(1)
  
  c_pseudo1 = 0.
  c_pseudo2 = 0.
  do i=1,n_jac
     c_pseudo1(i,i) = .04
     c_pseudo2(i,i) = 24.0
  enddo


!!MS$  jac=0.
!!MS$  do i=1,1
!!MS$     jac(i,i)= 1
!!MS$  enddo
!!MS$
!!MS$  hessin=0.
  
  hessin = hessin +  MATMUL( TRANSPOSE(jac1), MATMUL( c_pseudo1, jac1))+  MATMUL( TRANSPOSE(jac2), MATMUL( c_pseudo2, jac2))

 PRINT*,'invert Hessian with respect to x (unitless parameter) or p (physical parameter): '
  READ(5,*) w
  

  IF (w == 'p') THEN
! For calculating posterior parameter uncertainties the Hessian has to be 
! remapped to be a function of p (hesscol calculates the Hessian as a function
! of x). t denotes the mapping between p and x (coded in mo_mapping.f90); 
! t: p->x
! The transformation of the Hessian is then given by:
! H_p = transpose(JT) * H_x * JT + sum ((adx(i) * HT_i),i=1,n_vars)
! with JT = del(t(p))/del(p) -> Jacobian of t
! and HT_i = del^2(t_i(p))/del(p_j)del(p_k) -> Hessian of i-th component of t
!
! Here HT is taken as a 2-D field only because the components of the 
! transformation are linear indepandent, e.g. each HT_i has only one element 
! on the diagonal (i,i) and for the transformation we only need the sum of all
! HT_i Hessians.
!HEW-CHG-050105: p0s -> p0sl
     DO i = 1,n_vars
        IF (pflag(i) == 0) THEN
           IF (xsf(i) > 0.) THEN
              p(i) = x(i)*p0sl(i)/ABS(xsf(i))           
              jt(i,i) = 1./p0sl(i)
           ELSE
              p(i) = EXP((x(i)/ABS(xsf(i)))*LOG(p0sl(i)))
              jt(i,i) = 1./(p(i)*LOG(p0sl(i)))
              ht(i,i) = adx(i)*(-1.)/(p(i)*p(i)*LOG(p0sl(i)))
           ENDIF
           print*,i,p(i)
        ELSE
           STOP 'so far, this works only if all parameters are not bounded'
        END IF        
     ENDDO  
    
     hess = MATMUL( TRANSPOSE(jt), MATMUL( hessin, jt)) + ht

  ELSEIF (w == 'x') THEN
! If the Hessian should be evaluated/inverted with respect to x (optimisation 
! parameters, useful) nothing has to be changed. This is neccessary for
! calculating prognostic uncertainties as the Jacobians of a prognostic 
! quantity is also evaluated with respect to x.
     hess = hessin 
  ELSE
     STOP 'you either have to type in p or x.'
  ENDIF

  a = hess

  PRINT*,'symmetry check ',MINVAL(a -TRANSPOSE(a)),MAXVAL(a -TRANSPOSE(a))

  CALL syminv( hess, info)
  b = MATMUL( a, hess) -RESHAPE ( (/ (1.0, (0.0, k= 1,n_vars), j= 1,n_vars-1), 1.0 /), (/ n_vars, n_vars /) )

  PRINT*,'inversion check',MINVAL(b),MAXVAL(b)

  IF(info /= 0) STOP 'inversion failure'

  PRINT*,'sqrt of diagonal of inverse Hessian (parameter uncertainties)'
  DO i=1,n_vars
     WRITE(6,*) i, SQRT(ABS(hess( i,i)))
  END DO
  ! calculate eigen spectrum 
  PRINT*,'eigenvalues of the Hessian in x-space'
  IF (w == 'p') THEN
     a=hessin
  ENDIF
  CALL dsyev('v', 'u', n_vars, a, n_vars,&
       evals, work, 100000, info)
  ! write out e-vals
  DO i=1,n_vars
     WRITE(6,*) i,evals(i)
  END DO

  ! make hessian positive definite, assume that largest neg eigenvalue represent noise level
!  cutoff = SUM( MINLOC( evals, (evals > ABS(evals(1))))) ! smallest eval  > size of largest negative eval
!  if (minval(evals) .gt. 0.) cutoff=1
  cutoff = SUM( MINLOC( evals, (evals > 1))) !  all eval  < 1
  PRINT*,'noise floor at eval ',cutoff
  ! now reconstruct hessian as observed part (> cutoff) + prior part *< cutoff)
  e = 0.
  ! insert eigen-values into diagonal matrix
  DO i = cutoff, n_vars
     e( i, i) = evals( i)
  END DO
!  b = MATMUL( a, MATMUL( e, TRANSPOSE( a)))
  ! and now add prior part, evals here are always 1
!  e = 0.
  DO i = 1, cutoff -1
     e( i, i) = 1.
  END DO
  b =  MATMUL( a, MATMUL( e, TRANSPOSE( a)))
 
  IF (w == 'p') THEN
    
     hess = MATMUL( TRANSPOSE(jt), MATMUL( b, jt)) + ht

     OPEN(unit=1, file='output/denoised_hess_p.bin', form='unformatted')
  ELSEIF (w == 'x') THEN
     OPEN(unit=1, file='output/denoised_hess_x_pseudo2.bin', form='unformatted')
     hess=b
  ENDIF
  REWIND 1
  WRITE (1) hess
  CLOSE(1)

  CALL syminv( hess, info)


  PRINT*,'sqrt diagonal elements of denoised Hessian'
  DO i=1,n_vars
     WRITE(6,*) i,SQRT(ABS(hess( i,i)))
  END DO


  STOP
  
CONTAINS 

  ! syminv calls inverting routines from LAPACK ??
  SUBROUTINE syminv(a, info)
    IMPLICIT NONE
    REAL :: a(:,:)
    INTEGER :: info ! status return
    ! local variables
    INTEGER :: n ! dimension
    INTEGER :: i,j ! index variables
    REAL, DIMENSION(:), ALLOCATABLE :: work    
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwork,ipiv
    
    REAL anorm ! 1-norm of A from slansy
    REAL rcond ! reciprocal condition number
    REAL slansy ! anorm function

    n = SIZE(a, 1)
    ALLOCATE(work(n*n),ipiv(n*n),iwork(n*n))
    
    ! routine for inverting a real symmetric positive definite matrix
    ! e.g. a covariance matrix (I hope)
    
    CALL dsytrf ('u', n, A, n, IPIV, WORK, n*n, info)
    CALL dsytri ('u', n, a, n, IPIV, WORK, info)
    DO j =2,n
       DO i=1,j
          a(j,i) = a(i,j)
       END DO
    END DO
    RETURN
    
  END SUBROUTINE syminv
END PROGRAM inv_hess_pseudo
