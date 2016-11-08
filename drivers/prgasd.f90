!--------------------
!     set bit sizes
!---------------------
      module machine
      integer    ibit_kind
      parameter( ibit_kind = 4 )
      integer    bitsize
      parameter( bitsize = ibit_kind*8 )
      end module machine

!----------------------------------
!     ASD driver
!---------------------------------------
      program asd
!     adapted to new TAF generated parameterlist
!     adapted for different n and m
!        01-03 TK
      use machine
      implicit none
      integer    n, m, nsad, nsgd

!------------------------------------------------
!
!------------------------------------------------
      call setfunc( n, m )

      nsad = 1 + (m-1) / bitsize
      nsgd = 1 + (n-1) / bitsize





      print *, 'n, m   = ', n, m
      print *, 'nsad   = ', nsad
      print *, 'nsgd   = ', nsgd


      call doit( n, m, nsad, nsgd )

      end


      subroutine doit( n, m, nsad, nsgd )
      use machine
      implicit none

      integer    n, m, nsad, nsgd

      real       eps
      parameter( eps = .00001 )
      integer    maxiter
      parameter( maxiter = 1 )
! HEW-TEST      parameter( maxiter = 100000 )

      real x(n)
      real f(m)

      real adx(m,n)
      real adf(m,m)

      real tfunc, tad, tsad, tsgd

      integer(kind=ibit_kind)   :: sadx(nsad,n)
      integer(kind=ibit_kind)   :: sadf(nsad,m)
      integer(kind=ibit_kind)   :: sgdx(nsgd,n)
      integer(kind=ibit_kind)   :: sgdf(nsgd,m)
      integer(kind=ibit_kind)   :: zero = 0

      integer i, j, k, it
      integer ia, ib, ja, jb
      integer ijac(n,m)
      integer cnt0, cnt1, cr, cm

      real t0, t1

      real second
      external second

!------------------------------------------------
!------------------------------------------------

      call initfunc( n, x )
      call initfunc_sotl( n, x )
      call initfunc_soad( n, x )

!------------------------------------------------
! run and time function
!------------------------------------------------

      t0 = second()

      do it = 1, maxiter
         call func( n, x, m, f )
      end do

      t1 = second()
      tfunc = t1 - t0



!------------------------------------------------
! run and time AD
!------------------------------------------------
         do i = 1, n
            do k = 1, m
               adx(k,i) = 0.0
            end do
         end do
         do j = 1, m
            do k = 1, m
               adf(k,j) = 0.0
            end do
            adf(j,j) = 1.0
         end do


      t0 = second()

      do it = 1, maxiter
         call func_pad( m, n, x, adx, m, m, f, adf, m )
!         call adfunc( n, x, adx, m, f, adf )
      end do

      t1 = second()
      tad = t1 - t0


      print*,'Jacobian from reverse'
      do j=1,m
         print'(5(x,e8.1))',(adx(j,i),i=1,n)
      enddo
      print*,'Jacobian from reverse converted to integer'
      call matprint( m, n, adx )

!------------------------------------------------
! run and time ASD
!------------------------------------------------
         do i = 1, n
            do ja = 1, nsad
               sadx(ja,i) = 0
            end do
         end do
         do j = 1, m
            do ja = 1, nsad
               sadf(ja,j) = 0
            end do
         end do
         do j = 1, m
            ja = 1 + (j -1)/bitsize
            jb = j - (ja-1)*bitsize
            sadf(ja,j) = ibset( zero, jb-1 )
         end do
!      print*,'initialisation of transposed Jacobian from reverse'
!      call bitprint( nsad, n, sadx, m )
!      print*,'initialisation of sadf'
!      call bitprint( nsad, m, sadf, m )


      t0 = second()

      do it = 1, maxiter
!         call sadfunc( n, x, sadx, m, f, sadf )
         call func_spad( n, x, sadx, m, f, sadf )
      end do

      t1 = second()
      tsad = t1 - t0


      print*,'sparsity structure of Jacobian from reverse'
      call bitprint( nsad, n, sadx, m )

!------------------------------------------------
! run and time SGD
!------------------------------------------------
         do j = 1, n
            do ja = 1, nsgd
               sgdx(ja,j) = 0
            end do
         end do
         do j = 1, n
            ja = 1 + (j -1)/bitsize
            jb = j - (ja-1)*bitsize
            sgdx(ja,j) = ibset( zero, jb-1 )
         end do
         sgdf=0.


      t0 = second()

      do it = 1, maxiter
!         call sgdfunc( n, x, sgdx, m, f, sgdf )
         call func_sptl( n, x, sgdx, m, f, sgdf )
      end do

      t1 = second()
      tsgd = t1 - t0


      print*,'sparsity structure of Jacobian from forward'
      call bitprintt( nsgd, m, sgdf, n )
!------------------------------------------------
!------------------------------------------------
      tfunc= tfunc / float(maxiter)
      tad  = tad   / float(maxiter)
      tsgd = tsgd  / float(maxiter)
      tsad = tsad  / float(maxiter)

!------------------------------------------------
!------------------------------------------------
      write(6,*) 'absolute func  ', tfunc
      write(6,*) 'absolute AD    ', tad
      write(6,*) 'absolute SAD   ', tsad
      write(6,*) 'absolute SGD   ', tsgd

      if (tfunc .lt. eps ) then
            print*,'tfunc too small, no output of timing for AD, SAD and SGD relatve to the function  '
      else
         write(6,*) 'relative AD/func  ', tad/tfunc
         write(6,*) 'relative SAD/func ', tsad/tfunc
         write(6,*) 'relative SGD/func ', tsgd/tfunc
      end if
      if (tsad .lt. eps ) then
            print*,'tsad too small, no output for relative timing of AD realtive to ASD '
      else
         write(6,*) 'relative AD/SAD ', tad/tsad
         write(6,*) 'relative AD/SGD ', tad/tsgd
      end if

!------------------------------------------------
! output for xvgr, xmgrace, actually
!------------------------------------------------

      if (tfunc .lt. eps ) then
            print*,'tfunc too small, no output of timing for AD, SAD and SGD relatve to the function  '
      else
         write(6,'(i8,3f8.2,a)') m,tad/tfunc,tsad/tfunc,tsgd/tfunc,'    xvgr'
      end if

!------------------------------------------------
! compare sparsity structure
!------------------------------------------------

      do i = 1, m
         ia = 1 + (i -1)/bitsize
         ib = i - (ia-1)*bitsize
         do j = 1, n
            ja = 1 + (j -1)/bitsize
            jb = j - (ja-1)*bitsize
            if ((adx(i,j) .ne. 0. ) .neqv. btest(sadx(ia,j),ib-1)) then
      print '(a,2i3,xe12.4,l)','DIFFERENT Rev ',i,j,adx(i,j),btest(sadx(ia,j),ib-1)
            end if
            if ((adx(i,j) .ne. 0. ) .neqv. btest(sgdf(ja,i),jb-1)) then
      print '(a,2i3,xe12.4,l)','DIFFERENT Fwd ',i,j,adx(i,j),btest(sgdf(ja,i),jb-1)
      print*,'--------'
            end if
         end do
      end do


      end



      subroutine bitprintt( nsgd, m, ibmat, n )
!     output of transposed of integer matrix ibmat,
!     where bits in first argument are transformed to integers
      use machine
      implicit none
      integer n, m, nsgd
      integer(kind=ibit_kind)   :: ibmat(nsgd,m)
      integer ijac(n)
      integer i,j,ia,ib

      write(6,'(a)') '---------------------------------------'
      do j = 1, m
         do i = 1, n
            ia = 1 + (i -1)/bitsize
            ib = i - (ia-1)*bitsize
            if ( btest(ibmat(ia,j),ib-1) ) then
!            print*,'bitprint: ia, ib, m, n =',ia, ib, m, n
               ijac(i) = 1
            else
               ijac(i) = 0
            end if
         end do

         write(*,'(99i1)') (ijac(i),i=1,n)
      end do

      end
      subroutine bitprint( nsad, n, ibmat, m )
!     output of integer matrix ibmat,
!     where bits in first argument are transformed to integers
      use machine
      implicit none
      integer n, m, nsad
      integer(kind=ibit_kind)   :: ibmat(nsad,n)
      integer ijac(m,n)
      integer i,j,ja,jb

      write(6,'(a)') '---------------------------------------'
      do j = 1, m
            ja = 1 + (j -1)/bitsize
            jb = j - (ja-1)*bitsize
            do i = 1, n
               if ( btest(ibmat(ja,i),jb-1) ) then
!                  print*,'bitprint: ja, jb, j, i, m, n =',ja, jb, j, i,m, n
                  ijac(j,i) = 1
               else
                  ijac(j,i) = 0
               end if
            end do
      end do
      do j = 1,m
         write(*,'(99i1)') (ijac(j,i),i=1,n)
      enddo

      end
      subroutine matprint( n, m, xb )
!     output of integer matrix xb,
      implicit none
      integer n, m
      real xb(n,m)
      integer ijac(m)
      integer i,j

      write(6,'(a)') '---------------------------------------'
      do i = 1, n
         do j = 1, m
            if ( abs(xb(i,j)) .gt. 1.e-32 ) then
               ijac(j) = 1
            else
               ijac(j) = 0
            end if
         end do

         write(*,'(99i1)') (ijac(j),j=1,m)
      end do

      end
