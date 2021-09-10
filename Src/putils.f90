module putils

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains

      subroutine pwrite(ptcls,n1in,nrayspin,nunit,idmin,idmax,ncols,nprecision)
! routine to perform parallel write of ptcls array
      use mpi
      implicit none
!     real*8, dimension(n1,*) :: ptcls
      real*8, dimension(:,:) :: ptcls
      real*8, allocatable, dimension(:,:) :: tblock
      integer :: n1in,nrayspin,nunit,idmin,idmax,ncols,nprecision
      integer :: n1,nraysp
      character*17 :: myformat='(99(1x,1pe16.9 ))'
      character*2 :: a2
      integer :: nrayspmax,nlength,mycount,nraysw,i,l,mprocs,myrank,ierr
      integer ::  mpistat(MPI_STATUS_SIZE)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      n1=ubound(ptcls,1)
      nraysp=ubound(ptcls,2)
!     if(myrank.eq.0)write(6,*)'(pwrite) n1,nraysp=',n1,nraysp
!
!---format:
      nlength=nprecision+7
      call num2stringit(nlength,a2,2)
      myformat(11:12)=a2
      call num2stringit(nprecision,a2,2)
      myformat(14:15)=a2
!---
!
      if(idmax.le.nraysp)then
        do i=1,idmax
          if(myrank.eq.0)write(nunit,myformat)ptcls(1:ncols,i)
        enddo
        return
      endif
!
!------------------------------------------------------------------
      call MPI_REDUCE(nraysp,nrayspmax,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      if(myrank.eq.0)then
        allocate(tblock(n1,nrayspmax))
        mycount=0
        do i=1,nraysp
         mycount=mycount+1
         if(mycount.lt.idmin)cycle
         if(mycount.gt.idmax)exit
         write(nunit,myformat)ptcls(1:ncols,i)
        enddo
        do l=1,mprocs-1
          call MPI_RECV(tblock,n1*nrayspmax,MPI_DOUBLE,l,98,MPI_COMM_WORLD,mpistat,ierr)
          call MPI_GET_COUNT(mpistat,MPI_DOUBLE,nraysw,ierr)
          nraysw=nraysw/n1
          do i=1,nraysw
           mycount=mycount+1
           if(mycount.lt.idmin)cycle
           if(mycount.gt.idmax)exit
           write(nunit,myformat)tblock(1:ncols,i)
          enddo
        enddo
        deallocate(tblock)
      else
        call MPI_SEND(ptcls,n1*nraysp,MPI_DOUBLE,0,98,MPI_COMM_WORLD,ierr)
      endif
      return
      end
!
      subroutine num2stringit(num,a,ndigits)
! converts an integer "num" with at most "ndigits" digits into a character string "a"
! if ndigits exceeds the actual number of digits in num, then leading zeroes are inserted in the string
      implicit none
      integer num,ndigits
      character a*(*)
      integer n,m,k
      m=num
      do n=1,ndigits
      k=m/10**(ndigits-n)
      if(k.eq.0)a(n:n)='0'
      if(k.eq.1)a(n:n)='1'
      if(k.eq.2)a(n:n)='2'
      if(k.eq.3)a(n:n)='3'
      if(k.eq.4)a(n:n)='4'
      if(k.eq.5)a(n:n)='5'
      if(k.eq.6)a(n:n)='6'
      if(k.eq.7)a(n:n)='7'
      if(k.eq.8)a(n:n)='8'
      if(k.eq.9)a(n:n)='9'
      m=m-k*10**(ndigits-n)
      enddo
      return
      end

end module putils
