module numerical_distributions

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains


      subroutine initdist(ptcls,iseed,indevar,npart_gbl,mimicmprocs)
      use mpi
      use lwprobcons, only : getmatrixandcentroid
      implicit none
      real*8, dimension(:,:) :: ptcls
      integer*8 :: npart_gbl
      integer :: iseed,indevar,mimicmprocs
      integer, allocatable, dimension(:) :: seedarray
      integer :: npart_lcl,nbatchsize,nstart,seedarraysize
      real*8 :: ekin
      real*8, dimension(6) :: cent
      real*8, dimension(6,6) :: sigmat
      integer :: mprocs,myrank,mpierr,m,n
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      if(npart_gbl.eq.1)then
        if(mprocs.ne.1.and.myrank.eq.0)write(6,*)'error: 1-partcle run must use 1 MPI process'
        if(mprocs.ne.1)stop
        call getmatrixandcentroid(ekin,cent,sigmat)
        ptcls(1:6,1)=cent(1:6)
        if(myrank.eq.0)write(6,*)'1-particle run with initial condition='
        if(myrank.eq.0)write(6,'(6(1pe14.7,1x))')ptcls(1:6,1)
        return
      endif
      if(mprocs.eq.1)seedarraysize=mimicmprocs
      if(mprocs.ne.1)seedarraysize=mprocs
      allocate(seedarray(0:seedarraysize-1)) !every proc owns entire seedarray (not necessary)
      npart_lcl=size(ptcls,2)
      n=iseed
      do m=0,seedarraysize-1
        n=mod(4096*n+150889,714025)  !for use if mprocs < 714025
        seedarray(m)=n
      enddo
      if(mprocs.gt.1)then
        call initrandom(seedarray(myrank))
        nstart=1
        call init_particles(ptcls,npart_lcl,nstart,indevar)
      else
        write(6,*)'serial job with initial conditions that mimic ',mimicmprocs,' procs'
        if(mod(npart_gbl,mimicmprocs).ne.0)then
          if(myrank.eq.0)write(6,*)'error: npart_gbl must be divisible by mimicmprocs'
          call myexit
        endif
        nbatchsize=npart_gbl/mimicmprocs
        nstart=1
        do m=1,mimicmprocs
          call initrandom(seedarray(m-1))
          call init_particles(ptcls,nbatchsize,nstart,indevar)
          nstart=nstart+nbatchsize
        enddo
      endif
! < if desired this is a good place to set particle#1 on PE0 to the centroid value >
      return
      end


      subroutine init_particles(ptcls,ntot,nstart,indevar)
! initialize the particle array
      use mpi
      use lwprobcons, only : getmatrixandcentroid !!!!!,getmatrixandcentroid2,fracinbunch2
      implicit none
      real*8, dimension(:,:) :: ptcls
      integer :: ntot,nstart,indevar,nnn,ntotarg,nstartarg
      integer :: icholesky

      real*8 :: ekin,gam0,gb0
      integer :: disttype                ! =0 for uniform, =1 for Gaussian
      real*8 :: gaussiancutoff           !cutoff if a Gaussian initial condition is used
      real*8, dimension(6) :: cent       !initial beam centroid
      real*8, dimension(6,6) :: sigmat   !initial beam sigma (i.e. 2nd moment) matrix
      real*8 :: theta

      integer :: n

      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

!this is a mess. fix later.

!commented out old code that allowed for 2 bunches
!!!!!!do nnn=1,2
!!!!!!if(nnn.eq.1)call getmatrixandcentroid(ekin,cent,sigmat)
!!!!!!if(nnn.eq.2 .and. fracinbunch2.eq.0.d0)return
!!!!!!if(nnn.eq.2)call getmatrixandcentroid2(ekin,cent,sigmat)

      call getmatrixandcentroid(ekin,cent,sigmat)

! generate the initial distribution:
      disttype=1
      gaussiancutoff=5.d0
      icholesky=1 ! =0 for 2x2 block diagonal case; =1 for a general sigmat but only use if emittance .ne. 0
      ! use the cholesky method when the matrix is not 2x2 block diagonal:
      if(sigmat(1,3).ne.0.d0 .or. sigmat(1,4).ne.0.d0 .or. sigmat(1,5).ne.0.d0 .or. sigmat(1,6).ne.0.d0)icholesky=1
!     if(sigmat(2,3)... should check the other rows of sigmat to see if 2x2 block diagonal or not
      if(myrank.eq.0)write(6,*)'disttype=',disttype
      if(myrank.eq.0)write(6,*)'gaussiancutoff=',gaussiancutoff
      if(myrank.eq.0)write(6,*)'icholesky=',icholesky
!!!!!!if(nnn.eq.1)then
        nstartarg=1
        ntotarg=ntot !!!!!nint((1.d0-fracinbunch2)*ntot)
!!!!!!endif
!!!!!!if(nnn.eq.2)then
!!!!!!  nstartarg=nint((1.d0-fracinbunch2)*ntot)+1
!!!!!!  ntotarg=ntot-nint((1.d0-fracinbunch2)*ntot)
!!!!!!endif
      if(myrank.eq.0)write(6,*)'calling gendist with nstartarg,ntotarg=',nstartarg,ntotarg
      call gendist(ptcls,cent,sigmat,gaussiancutoff,ntotarg,nstartarg,disttype,icholesky) !generate the initial distribution
      if(myrank.eq.0)write(6,*)'back from gendist'

! columns 7 and 8 (not used in this version since nprop is hardwired to 6 in probcons.f90):
      if(size(ptcls,1).gt.6)then
        write(6,*)'error: code should not get here in this version'
        do n=nstart,nstart+ntot-1
          ptcls(7,n)=1.d0*(myrank*size(ptcls,2)+n) !particle id#
          ptcls(8,n)=0.d0 !lost particle flag
        enddo
      endif

!!!!!!enddo
      return
      end

      subroutine gendist(ptcls,cent,sigmat,gaussiancutoff,ntot,nstart,disttype,icholesky)
      use mpi
      implicit none
      real*8, dimension(:,:) :: ptcls !particle array; lower indices are 1
      real*8, dimension(6) :: cent,av
      real*8, dimension(6,6) :: sigmat
      real*8 :: gaussiancutoff
      integer :: ntot,nstart,disttype,icholesky
      real*8, dimension(6,6) :: sq
      real*8, dimension(6) :: vec,newvec
      integer :: n
      real*8 :: r1,r2,r3,r4,r5,r6,arg,u5,u6
      real*8 :: epsx2,epsy2,epsz2
      integer :: mprocs,myrank,mpierr,ierr

      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr) !mprocs and myrank would be needed for I/O (e.g. error messages) and to
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr) !compute global quantities (e.g. moments) in this routine if desired
      if(myrank.eq.0)write(6,*)'inside gendist with ntot,nstart=',ntot,nstart

      epsx2=sigmat(1,1)*sigmat(2,2)-sigmat(1,2)**2
      epsy2=sigmat(3,3)*sigmat(4,4)-sigmat(3,4)**2
      epsz2=sigmat(5,5)*sigmat(6,6)-sigmat(5,6)**2
      if(epsx2.lt.0.d0.and.myrank.eq.0)write(6,*)'error: x-rms emittance squared <0'
      if(epsy2.lt.0.d0.and.myrank.eq.0)write(6,*)'error: y-rms emittance squared <0'
      if(epsz2.lt.0.d0.and.myrank.eq.0)write(6,*)'error: z-rms emittance squared <0'
      if(epsx2.lt.0.d0.or.epsy2.lt.0.d0.or.epsz2.lt.0.d0)call myexit


      if(icholesky.eq.0)then
      do n=nstart,nstart+ntot-1
    9   continue
        call normdv(r1,r2)
        call normdv(r3,r4)
        call normdv(r5,r6)
        if(r1*r1+r3*r3+r5*r5.gt.gaussiancutoff**2)goto 9
! this version assumes that sigmat is 2x2 block diagonal
!x-gbx:
        ptcls(1,n)=r1*sqrt(sigmat(1,1))
        ptcls(2,n)=r1*sigmat(1,2)/sqrt(sigmat(1,1))+r2*sqrt(sigmat(2,2)-sigmat(1,2)**2/sigmat(1,1))
!y-gby:
        ptcls(3,n)=r3*sqrt(sigmat(3,3))
        ptcls(4,n)=r3*sigmat(3,4)/sqrt(sigmat(3,3))+r4*sqrt(sigmat(4,4)-sigmat(3,4)**2/sigmat(3,3))
!z-gbz:
        ptcls(5,n)=r5*sqrt(sigmat(5,5))
        ptcls(6,n)=r5*sigmat(5,6)/sqrt(sigmat(5,5))+r6*sqrt(sigmat(6,6)-sigmat(5,6)**2/sigmat(5,5))
      enddo
      endif

      if(icholesky.ne.0)then
      call cholesky(sigmat,sq,6,6,ierr)
      if(ierr.ne.0)then
        if(myrank.eq.0)write(6,*)'error return from cholesky; ierr=',ierr
        stop
      endif
      do n=nstart,nstart+ntot-1
        if(disttype.eq.0)then
   10     continue
          call random_number(r1); r1=2.d0*r1-1.d0
          call random_number(r2); r2=2.d0*r2-1.d0
          call random_number(r3); r3=2.d0*r3-1.d0
          call random_number(r4); r4=2.d0*r4-1.d0
          call random_number(r5); r5=2.d0*r5-1.d0
          call random_number(r6); r6=2.d0*r6-1.d0
          if(r1*r1+r3*r3+r5*r5.gt.1.d0)goto 10
        endif
        if(disttype.eq.1)then
   11     continue
          call normdv(r1,r2); call normdv(r3,r4); call normdv(r5,r6)
          call random_number(u5); u5=2.d0*u5-1.d0
          if(r1*r1+r3*r3.gt.gaussiancutoff**2)goto 11
        endif
        vec(1)=r1;vec(2)=r2;vec(3)=r3;vec(4)=r4;vec(5)=r5;vec(6)=r6
        newvec=matmul(sq,vec)
        ptcls(1:6,n)=newvec(1:6)
      enddo
      endif
! impose the desired centroid:
      do n=nstart,nstart+ntot-1
        ptcls(1:6,n)=ptcls(1:6,n)+cent(1:6)
      enddo
      return
      end subroutine gendist
!
      subroutine normdv(d1,d2)
! routine to generate 2 normal deviates
      implicit none
      real*8 :: u1,u2
      real*8 twopi,d1,d2
      twopi=4.0d0*asin(1.0d0)
      call random_number(u1)
      call random_number(u2)
      d1=sqrt(-2.0d0*log(u1))*cos(twopi*u2)
      d2=sqrt(-2.0d0*log(u1))*sin(twopi*u2)
      return
      end subroutine normdv
!
!not sure where I got cholesky from
      subroutine cholesky(a,b,n,np,ierr)
      implicit none
      integer :: n,np,ierr
      real*8, dimension(np,np) :: a,b
      real*8, dimension(n) :: p
      integer :: i,j,k
      real*8 :: sumx
      ierr=-1
      if(n.gt.np)return
      b(:,:)=a(:,:)
      do i=1,n
        do j=i,n
          sumx=b(i,j)
          do k=i-1,1,-1
            sumx=sumx-b(i,k)*b(j,k)
          enddo
          if(i.eq.j)then
            if(sumx.le.0.d0)ierr=-i
            if(sumx.le.0.d0)return
            p(i)=sqrt(sumx)
          else
            b(j,i)=sumx/p(i)
          endif
        enddo
      enddo
      do i=1,n
        do j=i,n
          b(i,j)=0.d0
          if(i.eq.j)b(i,i)=p(i)
        enddo
      enddo
      ierr=0
      return
      end

end module
