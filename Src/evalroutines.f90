module evalroutines

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains

      subroutine evalt(t,y,f) ! t is the independent variable
      use mpi
      use lwprobcons, only : selfconsistent,nstartselfcon,nstepcounter
      use lienwiech, only : getselfconsistentfield
      use externalfield, only : getexternalfield
      implicit none
      real*8 :: t
      real*8, dimension(:,:) :: y,f
      integer :: n1,np,n,ierr
      real*8, dimension(3) :: e0,b0 !external field, calculated on-the-fly for each particle
      real*8, dimension(3) :: eslf,bslf
      real*8, allocatable, dimension(:,:) :: ebself
      logical, allocatable, dimension(:) :: msk
      real*8, dimension(6) :: av
      integer :: ngoodlcl,ngood
      real*8 :: gambet_ave_lcl,gambet_ave,gamma_ave,beta_ave !average values involving gamma,beta
      real*8 :: clite,cliteinv
      real*8 :: xmc2,qbymcc,qbymc,qbym
      real*8 :: gam,gam2,gbz2,gbz,den1
      real*8 :: gbx0,gby0,gbz0,gam0
      integer, save :: ncount=1
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      n1=ubound(y,1)
      np=ubound(y,2)
!     write(6,*)'n1,np=',n1,np
      clite=299792458.d0
      cliteinv=1.d0/clite
      xmc2=0.510998910d6
      qbymcc=1.d0/xmc2
      qbymc=1.d0/xmc2*clite
      qbym=1.d0/xmc2*clite**2
!
!space-charge calc:
      allocate(ebself(6,np))
      ebself(:,:)=0.d0
      if(selfconsistent.and.nstepcounter.ge.nstartselfcon)call getselfconsistentfield(t,np,y,ebself)
 
!loop over particles:
      f(1:n1,1:np)=0.d0 !ensures f=0 for lost particles and for columns beyond 6
      do n=1,np
!!!!!!!!if(y(n1,n).ne.0.d0)cycle !no lost particle flag in this version
        gam2=1.d0+y(2,n)**2+y(4,n)**2+y(6,n)**2
        gam=sqrt(gam2)
        eslf(1:3)=ebself(1:3,n)
        bslf(1:3)=ebself(4:6,n)
        call getexternalfield(y(1,n),y(3,n),y(5,n),t,e0,b0,ierr)
!       if(ierr.ne.0)then
!         y(n1,n)=-1.d0
!         cycle
!       endif
        e0(1:3)=e0(1:3)+eslf(1:3)
        b0(1:3)=b0(1:3)+bslf(1:3)
        f(1,n)=y(2,n)*clite/gam
        f(2,n)=qbymc*e0(1)+qbym*(y(4,n)*b0(3)-y(6,n)*b0(2))/gam
        f(3,n)=y(4,n)*clite/gam
        f(4,n)=qbymc*e0(2)+qbym*(y(6,n)*b0(1)-y(2,n)*b0(3))/gam
        f(5,n)=y(6,n)*clite/gam
        f(6,n)=qbymc*e0(3)+qbym*(y(2,n)*b0(2)-y(4,n)*b0(1))/gam
!skip!  f(7,n)=0.d0
      enddo
      return
      end

end module
