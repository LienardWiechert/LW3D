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

      subroutine evalt1(t,y,f) ! 1particle, t is the independent variable
      use mpi
      implicit none
      real*8 :: t
      real*8, dimension(:,:) :: y,f
      integer :: n1,np,n,ierr
      real*8, dimension(3) :: e0,b0 !external field, calculated on-the-fly for each particle
      real*8 :: clite,cliteinv
      real*8 :: xmc2,qbymcc,qbymc,qbym
      real*8 :: gam,gam2,gbz2,gbz,den1
      real*8 :: gbx0,gby0,gbz0,gam0
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      n1=ubound(y,1)
!!!   np=ubound(y,2)
      clite=299792458.d0
      cliteinv=1.d0/clite
      xmc2=0.510998910d6
      qbymcc=1.d0/xmc2
      qbymc=1.d0/xmc2*clite
      qbym=1.d0/xmc2*clite**2
!
      f(1:n1,1)=0.d0 !ensures f=0 for lost particles and for columns beyond 6
        gam2=1.d0+y(2,1)**2+y(4,1)**2+y(6,1)**2
        gam=sqrt(gam2)
        call getexternalfield(y(1,1),y(3,1),y(5,1),t,e0,b0,ierr)
        f(1,1)=y(2,1)*clite/gam
        f(2,1)=qbymc*e0(1)+qbym*(y(4,1)*b0(3)-y(6,1)*b0(2))/gam
        f(3,1)=y(4,1)*clite/gam
        f(4,1)=qbymc*e0(2)+qbym*(y(6,1)*b0(1)-y(2,1)*b0(3))/gam
        f(5,1)=y(6,1)*clite/gam
        f(6,1)=qbymc*e0(3)+qbym*(y(2,1)*b0(2)-y(4,1)*b0(1))/gam
!skip!  f(7,1)=0.d0
      return
      end

      subroutine getexternalfield(x,y,z,t,e0,b0,ierr)
      implicit none
      real*8 :: x,y,z,t
      real*8, dimension(3) :: e0,b0
      integer :: ierr
!     real*8, parameter :: bfielddip=0.3352642711907601d0 !B field for 100 MeV ss dipole case with rho=1m
      real*8, parameter :: pi=2.d0*asin(1.d0)
      real*8 :: bval,dz1bval,dz2bval,dz3bval,dz4bval
      real*8 :: eps,byplus,byminus,zshft,z1,z2,z3,z4
      integer :: n
      real*8 :: frntanh1,dz1frntanh1,dz2frntanh1,dz3frntanh1,dz4frntanh1
      real*8 :: frntanh2,dz1frntanh2,dz2frntanh2,dz3frntanh2,dz4frntanh2
      real*8 :: frntanh3,dz1frntanh3,dz2frntanh3,dz3frntanh3,dz4frntanh3
      real*8 :: frntanh4,dz1frntanh4,dz2frntanh4,dz3frntanh4,dz4frntanh4
      real*8, parameter :: twopi=4.d0*asin(1.d0)
!     real*8, parameter :: bcon= 0.0669164882226981d0    !corresponds to K=0.1
      real*8, parameter :: bcon=0.00669164882226981d0    !corresponds to K=0.01
      real*8, parameter :: lambda_u=16.d-3
      real*8, parameter :: n_u=50.d0
!     real*8, parameter :: k_u=twopi/lambda_u
      e0(1:3)=0.d0
      b0(1:3)=0.d0
!dipole hack HACK:
      if(t.ne.-12345678.)then
        b0(2)=0.85225274079636459d0 ! rho=1m for an electron with gamma=500.
!       b0(2)=0.3352642711907601d0 ! 0.33526427T gives rho=1m for 100MeV electron
        ierr=0
        return
      endif
!
!     b0(2)=bfielddip
!     if(z.ge.zmin.and.z.le.zmax)then
!       b0(2)=bcon*cos(k_u*z)
!     endif
      byplus=bcon*1.00d0
      byminus=-bcon*1.00d0
      bval=0.d0; dz1bval=0.d0; dz2bval=0.d0; dz3bval=0.d0; dz4bval=0.d0
      n=1
      eps=1.d-2/pi
      zshft=(n-1)*lambda_u
      z1=zshft+0.125*lambda_u; z2=zshft+0.375*lambda_u; z3=zshft+0.625*lambda_u; z4=zshft+0.875*lambda_u
      z1=z1+0.125*lambda_u
      call frntanh01234(z-z1,eps,frntanh1,dz1frntanh1,dz2frntanh1,dz3frntanh1,dz4frntanh1)
      call frntanh01234(z-z2,eps,frntanh2,dz1frntanh2,dz2frntanh2,dz3frntanh2,dz4frntanh2)
      call frntanh01234(z-z3,eps,frntanh3,dz1frntanh3,dz2frntanh3,dz3frntanh3,dz4frntanh3)
      call frntanh01234(z-z4,eps,frntanh4,dz1frntanh4,dz2frntanh4,dz3frntanh4,dz4frntanh4)
      bval=bval+byplus*frntanh1-byplus*frntanh2+byminus*frntanh3-byminus*frntanh4
      dz1bval=dz1bval+byplus*dz1frntanh1-byplus*dz1frntanh2+byminus*dz1frntanh3-byminus*dz1frntanh4
      dz2bval=dz2bval+byplus*dz2frntanh1-byplus*dz2frntanh2+byminus*dz2frntanh3-byminus*dz2frntanh4
      dz3bval=dz3bval+byplus*dz3frntanh1-byplus*dz3frntanh2+byminus*dz3frntanh3-byminus*dz3frntanh4
      dz4bval=dz4bval+byplus*dz4frntanh1-byplus*dz4frntanh2+byminus*dz4frntanh3-byminus*dz4frntanh4

      do n=2,nint(n_u)-1
        zshft=(n-1)*lambda_u
        z1=zshft+0.125*lambda_u; z2=zshft+0.375*lambda_u; z3=zshft+0.625*lambda_u; z4=zshft+0.875*lambda_u
        call frntanh01234(z-z1,eps,frntanh1,dz1frntanh1,dz2frntanh1,dz3frntanh1,dz4frntanh1)
        call frntanh01234(z-z2,eps,frntanh2,dz1frntanh2,dz2frntanh2,dz3frntanh2,dz4frntanh2)
        call frntanh01234(z-z3,eps,frntanh3,dz1frntanh3,dz2frntanh3,dz3frntanh3,dz4frntanh3)
        call frntanh01234(z-z4,eps,frntanh4,dz1frntanh4,dz2frntanh4,dz3frntanh4,dz4frntanh4)
        bval=bval+byplus*frntanh1-byplus*frntanh2+byminus*frntanh3-byminus*frntanh4
        dz1bval=dz1bval+byplus*dz1frntanh1-byplus*dz1frntanh2+byminus*dz1frntanh3-byminus*dz1frntanh4
        dz2bval=dz2bval+byplus*dz2frntanh1-byplus*dz2frntanh2+byminus*dz2frntanh3-byminus*dz2frntanh4
        dz3bval=dz3bval+byplus*dz3frntanh1-byplus*dz3frntanh2+byminus*dz3frntanh3-byminus*dz3frntanh4
        dz4bval=dz4bval+byplus*dz4frntanh1-byplus*dz4frntanh2+byminus*dz4frntanh3-byminus*dz4frntanh4
      enddo

      n=n_u
      zshft=(n-1)*lambda_u
      z1=zshft+0.125*lambda_u; z2=zshft+0.375*lambda_u; z3=zshft+0.625*lambda_u; z4=zshft+0.875*lambda_u
      z2=z2-0.125*lambda_u
      call frntanh01234(z-z1,eps,frntanh1,dz1frntanh1,dz2frntanh1,dz3frntanh1,dz4frntanh1)
      call frntanh01234(z-z2,eps,frntanh2,dz1frntanh2,dz2frntanh2,dz3frntanh2,dz4frntanh2)
      bval=bval+byplus*frntanh1-byplus*frntanh2!+byminus*frntanh(z-z3,eps)-byminus*frntanh(z-z4,eps)
      dz1bval=dz1bval+byplus*dz1frntanh1-byplus*dz1frntanh2!+byminus*dz1frntanh(z-z3,eps)-byminus*dz1frntanh(z-z4,eps)
      dz2bval=dz2bval+byplus*dz2frntanh1-byplus*dz2frntanh2!+byminus*dz2frntanh(z-z3,eps)-byminus*dz2frntanh(z-z4,eps)
      dz3bval=dz3bval+byplus*dz3frntanh1-byplus*dz3frntanh2!+byminus*dz3frntanh(z-z3,eps)-byminus*dz3frntanh(z-z4,eps)
      dz4bval=dz4bval+byplus*dz4frntanh1-byplus*dz4frntanh2!+byminus*dz4frntanh(z-z3,eps)-byminus*dz4frntanh(z-z4,eps)
!hardwired for no x-dependence:
      b0(2)=bval-0.5d0*y**2*dz2bval+y**4/24.d0*dz4bval
      b0(3)=y*dz1bval-y**3/6.d0*dz3bval
      ierr=0
      return
      end

      real*8 function frntanh(y,eps)
      implicit none
      real*8 :: y,eps
      frntanh = tanh(y/eps)/2.d0 +.5d0
      return
      end

      real*8 function dz1frntanh(y,eps)
      implicit none
      real*8 :: y,eps
      if (abs(y/eps).lt.10.d0) then
        dz1frntanh = 2*Exp(2*y/eps)/(eps*(1 + Exp(2*y/eps))**2)
      else
        dz1frntanh = 0.d0
      endif
      return
      end
!
      real*8 function dz2frntanh(y,eps)
      implicit none
      real*8 :: y,eps
      if (abs(y/eps).lt.10.d0) then
        dz2frntanh = (-4*Exp((2*y)/eps)*                       &
     &  (-1 + Exp((2*y)/eps)))/((1 + Exp((2*y)/eps))**3*eps**2)
      else
        dz2frntanh = 0.d0
      endif
      return
      end
!
      real*8 function dz3frntanh(y,eps)
      implicit none
      real*8 :: y,eps
      if (abs(y/eps).lt.10.d0) then
        dz3frntanh = 8*Exp(2*y/eps)*                           &
     &      (1 - 4*Exp(2*y/eps) + Exp(4*y/eps))/               &
     &   (eps**3*(1 + Exp(2*y/eps))**4)
      else
        dz3frntanh = 0.0d0
      endif
      return
      end
!
      real*8 function dz4frntanh(y,eps)
      implicit none
      real*8 :: y,eps
      if (abs(y/eps).lt.10.d0) then
        dz4frntanh =                                                    &
     &(-16*Exp((2*y)/eps)*(-1 + 11*Exp((2*y)/eps) - 11*Exp((4*y)/eps) + &
     &      Exp((6*y)/eps)))/((1 + Exp((2*y)/eps))**5*eps**4)
      else
        dz4frntanh = 0.d0
      endif
      return
      end

      subroutine frntanh01234(y,eps,frntanh,dz1frntanh,dz2frntanh,dz3frntanh,dz4frntanh)
      implicit none
      real*8 :: y,eps,frntanh,dz1frntanh,dz2frntanh,dz3frntanh,dz4frntanh
      real*8 :: exp2yeps,exp4yeps,exp6yeps

      frntanh = tanh(y/eps)/2.d0 +.5d0

      if(abs(y/eps).ge.10.d0)then
        dz1frntanh = 0.d0
        dz2frntanh = 0.d0
        dz3frntanh = 0.d0
        dz4frntanh = 0.d0
        return
      endif

      exp2yeps=Exp(2.d0*y/eps)
!     exp4yeps=Exp(4.d0*y/eps)
!     exp6yeps=Exp(6.d0*y/eps)
      exp4yeps=exp2yeps**2
      exp6yeps=exp2yeps*exp4yeps

        dz1frntanh = 2*exp2yeps/(eps*(1 + exp2yeps)**2)
 
        dz2frntanh = (-4*exp2yeps*                                      &
     &  (-1 + exp2yeps))/((1 + exp2yeps)**3*eps**2)
 
        dz3frntanh = 8*exp2yeps*                                        &
     &      (1 - 4*exp2yeps + exp4yeps)/                                &
     &   (eps**3*(1 + exp2yeps)**4)
 
        dz4frntanh =                                                    &
     &(-16*exp2yeps*(-1 + 11*exp2yeps - 11*exp4yeps +                   &
     &      exp6yeps))/((1 + exp2yeps)**5*eps**4)

      return
      end

end module
