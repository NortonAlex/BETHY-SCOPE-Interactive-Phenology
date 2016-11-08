program line
  implicit none
  integer :: n

! get the number of control variables
  call numbmod( n )

! action
  call linerun( n )

end program line

subroutine linerun( n )
  implicit none

  integer :: n
  integer, parameter  :: uopt = 7        ! unit for options
  integer, parameter  :: usect = 22      ! unit for section
  integer, parameter  :: upar1 = 23      ! unit for parameter values 
  integer, parameter  :: upar2 = 24      ! unit for parameter values 
  integer, parameter  :: udx = 25        ! unit for search direction
  integer, parameter  :: uh = 26         ! unit for section with hessian
  logical :: hessian = .false.
  logical :: direction = .false.

  real    :: fc, fc0, fh
  real    :: x0(n), xend(n), x(n), t, dx(n), ndx, h(n,n)

  integer :: nt = 10                     ! number of steps in search direction
  integer :: its=-2,ite=10               ! start and end point 
  integer :: ixs=1,ixe                   ! use search direction projected onto range ixs to ixe
  real :: a = 1.                         ! scaling factor default: 1.

  integer :: i, j, it, ix
  namelist/nlist/a,nt,its,ite,ixs,ixe

! check for too large n
  if (n.gt.99) print *, 'Line: n = ',n,' too large, change print format'

! default values for options
  ixe = n

! read namelist for options
      open( uopt, file   = 'line.par', access = 'sequential', form = 'formatted'  )
      read( uopt, nlist, end=9000, err=8000 )
      print *, ' Line : options have been read'
 8000 continue
 9000 continue
      close( uopt )

  ixe = min(n,ixe)
  ixs = max(1,ixs)

! initialisize the model
  call initmod( n, x0 )

! get test interval
  
  print*, 'read final point '
  open (39,file='xf.b',status='unknown',form='unformatted')
  read(39) x0
  close(39)

  if (direction) then
     print*, 'read search direction '
     open (39,file='dx.b',status='unknown',form='unformatted')
     read(39) dx
     close(39)     
  else
     print*, 'read previous point from line search '
     open (39,file='lx.b',status='unknown',form='unformatted')
     read(39) xend
     close(39)
     ! compute delta x 
     dx = 0.
     do ix=ixs,ixe
        dx(ix) = xend(ix)-x0(ix)
     enddo
  endif

! compute norm of delta x 
  ndx = 0.
  do ix=1,n
     ndx=ndx+(dx(ix))**2
  enddo
  ndx = sqrt(ndx)

  print*, '|dx| = ',ndx
  print*, 'a |dx| = ',a*ndx

! get hessian
  if (hessian) then
     INQUIRE (file='output/hess.bin', exist=hessian)
     open (39,file='output/hess.bin',status='old',form='unformatted')
     read(39) h
     close(39)
     fh = 0.
     do i = 1, n
        do j = 1, n
           fh = fh + dx(i)*h(i,j)*dx(j)
        enddo
     enddo
     print*, 'Curvature in dx direction = ', fh/ndx/ndx
     fh = fh/2.
     open (uh,file='hess.set',form='formatted')
     write (uh,'(A)') '@WITH G00'
     write (uh,'(A)') '@G00 ON'
     write (uh,'(A)') '@TYPE xy'
     write (uh,'(A,E16.8)') '@WORLD XMAX', ndx
     write (uh,'(A)') '@s0 symbol 1'
     write (uh,'(A)') '@s0 line type 0'
     write (uh,'(A)') '@xaxis  label "norm of difference of control vector to base point"'
     write (uh,'(A)') '@TITLE "cost function"'
  endif

!  open units
  open (usect,file='line.set',form='formatted')
  open (upar1,file='parset.dat',form='formatted')
  open (upar2,file='parlist.dat',form='formatted')
  open (udx,file='dx.dat',form='formatted')

!  dx
  do ix=1,n
     write (udx,*) ix, dx(ix)
  enddo

!  xmgrace header
  write (usect,'(A)') '@WITH G00'
  write (usect,'(A)') '@G00 ON'
  write (usect,'(A)') '@TYPE xy'
  write (usect,'(A,E16.8)') '@WORLD XMAX', ndx
  write (usect,'(A)') '@s0 symbol 1'
  write (usect,'(A)') '@s0 line type 0'
  write (usect,'(A)') '@xaxis  label "norm of difference of control vector to base point"'
  write (usect,'(A)') '@TITLE "cost function"'

  write (upar1,'(A,A19,99I20)') '#', 't', (i,i=1,n)

! evaluate at end point of optimisation
  fc0 = 0.
  IF (ite-its.GT.1) call model( n, x0, fc0)

! this is the loop that does the section
  do it=its,ite
     t=a*float(it)/nt
     x=x0+t*dx
     fc = 0.
     call model( n, x, fc)
     write (usect,*) t*ndx,fc-fc0
     write (upar1,'(100e20.6)') t*ndx,x(:)
     write (upar2,*) '# parameter vector for it = ',it
     if (hessian) write (uh,*) t*ndx,fh * t**2
     do i=1,n
        write (upar2,*) i, x(i)
     enddo
  enddo

!  close units
  close (usect)
  close (upar1)
  close (upar2)
  close (udx)
  if (hessian) close (uh)
end subroutine linerun
