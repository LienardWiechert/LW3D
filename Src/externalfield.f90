module externalfield

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains

      subroutine getexternalfield(x,y,z,t,e0,b0,ierr)
      implicit none
      real*8 :: x,y,z,t
      real*8, dimension(3) :: e0,b0
      integer :: ierr
      real*8 :: bval,dz1bval,dz2bval,dz3bval,dz4bval
      integer :: n
      integer, parameter :: ndipoles=4
      real*8, parameter :: thet0=0.349065850398866  !20 degrees
      real*8, parameter :: drouter=1.4d0*cos(thet0)
      real*8, parameter :: drmiddle=2.563d0
!     real*8, parameter :: effarclen=0.20d0*0.979815536051016d0
!     real*8, parameter :: dd=effarclen*sin(thet0)/thet0 ! dipole z-length
      real*8, parameter :: dd=0.3d0 ! dipole z-length
      real*8, parameter :: effarclen=dd*thet0/sin(thet0)
      real*8, parameter :: z0=0.0d0   !drift space before first dipole
      real*8, parameter :: gam0=120    !1.d0+60.d6/0.510998910d6   !NOTE WELL: fix later to get energy from probcons.f90
      real*8, parameter :: gb0=sqrt((gam0+1.d0)*(gam0-1.d0))
      real*8, parameter :: brho=gb0/299792458.d0*0.510998910d6
      real*8, parameter :: bchic=thet0*abs(brho)/effarclen
      real*8, parameter, dimension(ndipoles) :: z1=(/z0,      dd+drouter+z0, 2*dd+drmiddle+drouter+z0, 3*dd+drmiddle+2*drouter+z0/) !magnet leading edge
      real*8, parameter, dimension(ndipoles) :: z2=(/dd+z0, 2*dd+drouter+z0, 3*dd+drmiddle+drouter+z0, 4*dd+drmiddle+2*drouter+z0/) !magnet trailing edge
!     real*8, parameter, dimension(ndipoles) :: b00=(/-bchic,bchic,bchic,-bchic/)  !chicane magnet B field
      real*8, parameter, dimension(ndipoles) :: b00=(/+bchic,-bchic,-bchic,+bchic/)  !chicane magnet B field
!     real*8, parameter :: eps2=0.100d0/(2.d0*asin(1.d0)) !chicane magnet fringe field parameter
      real*8, parameter :: eps2=0.025  !0.060d0/(2.d0*asin(1.d0)) !chicane magnet fringe field parameter

      e0(1:3)=0.d0
      b0(1:3)=0.d0

!     call getquadbfield(x,y,z,b0)
!     call getsextbfield(x,y,z,b0)

      bval=0.d0
      dz1bval=0.d0
      dz2bval=0.d0
      dz3bval=0.d0
      dz4bval=0.d0
      do n=1,ndipoles
        bval=bval+b00(n)*frntanh(z-z1(n),eps2)-b00(n)*frntanh(z-z2(n),eps2)
        dz1bval=dz1bval+b00(n)*dz1frntanh(z-z1(n),eps2)-b00(n)*dz1frntanh(z-z2(n),eps2)
        dz2bval=dz2bval+b00(n)*dz2frntanh(z-z1(n),eps2)-b00(n)*dz2frntanh(z-z2(n),eps2)
        dz3bval=dz3bval+b00(n)*dz3frntanh(z-z1(n),eps2)-b00(n)*dz3frntanh(z-z2(n),eps2)
        dz4bval=dz4bval+b00(n)*dz4frntanh(z-z1(n),eps2)-b00(n)*dz4frntanh(z-z2(n),eps2)
      enddo
!hardwired for no x-dependence:
      b0(2)=b0(2) + bval   -0.5d0*y**2*dz2bval+y**4/24.d0*dz4bval
      b0(3)=b0(3) + y*dz1bval   -y**3/6.d0*dz3bval
      ierr=0
      return
      end

      subroutine undulatorgetexternalfield(x,y,z,t,e0,b0,ierr)
      implicit none
      real*8 :: x,y,z,t
      real*8, dimension(3) :: e0,b0
      integer :: ierr
      real*8 :: sig,ccon,g0,hinv
      real*8, parameter :: pi=2.d0*asin(1.d0)
      e0(1:3)=0.d0
      b0(1:3)=0.d0

!     sig= 1.0d-2
!     ccon=1.d0/(sqrt(2.d0*pi)*sig)
      hinv=1.d0/2.5d-2 !inverse of half-width of nearly square pulse
! 4.0 and 1.d-2 are pretty good
! 7.0 and 9.d-2 are better
        g0= 4.60d-0*exp(-0.5d0*(hinv*(z+0.20d0))**20) &
     &    + 4.60d-0*exp(-0.5d0*(hinv*(z+0.20d0-1.d0*2.56d0))**20) &
     &     +4.60d-0*exp(-0.5d0*(hinv*(z+0.20d0-2.d0*2.56d0))**20) &
     &    + 4.60d-0*exp(-0.5d0*(hinv*(z+0.20d0-3.d0*2.56d0))**20) &
     &     +4.60d-0*exp(-0.5d0*(hinv*(z+0.20d0-4.d0*2.56d0))**20) &
     &    + 4.60d-0*exp(-0.5d0*(hinv*(z+0.20d0-5.d0*2.56d0))**20) &
     &     +4.60d-0*exp(-0.5d0*(hinv*(z+0.20d0-6.d0*2.56d0))**20)
        if(abs(g0).lt.1.d-18)g0=0.d0 !not necessary but prevents a gnuplot plotting artifact
        b0(1)=g0*y
        b0(2)=g0*x


      call getundulatorfield(x,y,z,t,e0,b0,0.d0,ierr)
      call getundulatorfield(x,y,z,t,e0,b0,2.56d0,ierr)
      call getundulatorfield(x,y,z,t,e0,b0,5.12d0,ierr)
      call getundulatorfield(x,y,z,t,e0,b0,7.68d0,ierr)
      call getundulatorfield(x,y,z,t,e0,b0,10.24d0,ierr)
      call getundulatorfield(x,y,z,t,e0,b0,12.80d0,ierr)
      return
      end

      subroutine getundulatorfield(x,y,z,t,e0,b0,zstart,ierr)
      implicit none
      real*8 :: x,y,z,t,zstart
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
      real*8, parameter :: bcon=0.7918d0
      real*8, parameter :: lambda_u=2.8d-2
      real*8, parameter :: n_u=77.d0
!     real*8, parameter :: k_u=twopi/lambda_u
!     real*8, parameter :: zmax= (n_u/2.d0-0.25d0)*lambda_u
!     real*8, parameter :: zmax= (n_u/2.d0-0.00d0)*lambda_u
!     real*8, parameter :: zmax= (n_u/2.d0-0.001d0)*lambda_u
!     real*8, parameter :: zmin=-zmax
! IN THIS VERSION E0 AND B0 ARE ADDED TO
!     e0(1:3)=0.d0
!     b0(1:3)=0.d0
!dipole hack HACK:
!x    if(t.ne.-12345678.)then
!x      b0(2)=0.3352642711907601d0 ! 0.33526427T gives rho=1m for 100MeV electron
!x      ierr=0
!x      return
!x    endif
!
!     b0(2)=bfielddip
!     if(z.ge.zmin.and.z.le.zmax)then
!       b0(2)=bcon*cos(k_u*z)
!     endif
      byplus=bcon*1.00d0
      byminus=-bcon*1.00d0
      bval=0.d0; dz1bval=0.d0; dz2bval=0.d0; dz3bval=0.d0; dz4bval=0.d0
      n=1
      eps=1.5d-2/pi ! 4.d-3/pi
      zshft=(n-1)*lambda_u + zstart
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

!!!   eps=4.d-3/pi
!!!!!!do n=2,n_u-1
      do n=2,nint(n_u)-1
        zshft=(n-1)*lambda_u + zstart
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

!!!   eps=4.d-3/pi
      n=n_u
      zshft=(n-1)*lambda_u + zstart
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
      b0(2)=b0(2) + bval  -0.5d0*y**2*dz2bval+y**4/24.d0*dz4bval ! add to the value that b0 had on entry to routine
      b0(3)=b0(3)         + y*dz1bval-y**3/6.d0*dz3bval
      ierr=0
      return
      end

      subroutine oldgetexternalfield(x,y,z,t,e0,b0,ierr)
      implicit none
      real*8 :: x,y,z,t
      real*8, dimension(3) :: e0,b0
      integer :: ierr
!     real*8, parameter :: bfielddip=0.3352642711907601d0 !B field for 100 MeV ss dipole case with rho=1m
      real*8, parameter :: pi=2.d0*asin(1.d0)
      real*8 :: bval,dz1bval,dz2bval,dz3bval,dz4bval
      real*8 :: eps,byplus,byminus,zshft,z1,z2,z3,z4
      integer :: n
      real*8, parameter :: twopi=4.d0*asin(1.d0)
      real*8, parameter :: bcon=0.1d0          ! 0.1 Tesla
      real*8, parameter :: lambda_u=0.01d0     ! 1 cm
      real*8, parameter :: k_u=twopi/lambda_u
      real*8, parameter :: n_u=16.d0
!     real*8, parameter :: zmax= (n_u/2.d0-0.25d0)*lambda_u
!     real*8, parameter :: zmax= (n_u/2.d0-0.00d0)*lambda_u
      real*8, parameter :: zmax= (n_u/2.d0-0.001d0)*lambda_u
      real*8, parameter :: zmin=-zmax
      e0(1:3)=0.d0
      b0(1:3)=0.d0
!     b0(2)=bfielddip
!     if(z.ge.zmin.and.z.le.zmax)then
!       b0(2)=bcon*cos(k_u*z)
!     endif
      byplus=0.1d0*1.35d0
      byminus=-0.1d0*1.35d0
      bval=0.d0; dz1bval=0.d0; dz2bval=0.d0; dz3bval=0.d0; dz4bval=0.d0
      n=1
      eps=4.d-3/pi
      zshft=(n-1)*1.d-2
      z1=zshft+0.125d-2; z2=zshft+0.375d-2; z3=zshft+0.625d-2; z4=zshft+0.875d-2
      z1=z1+0.125d-2
      bval=bval+byplus*frntanh(z-z1,eps)-byplus*frntanh(z-z2,eps)+byminus*frntanh(z-z3,eps)-byminus*frntanh(z-z4,eps)
 dz1bval=dz1bval+byplus*dz1frntanh(z-z1,eps)-byplus*dz1frntanh(z-z2,eps)+byminus*dz1frntanh(z-z3,eps)-byminus*dz1frntanh(z-z4,eps)
 dz2bval=dz2bval+byplus*dz2frntanh(z-z1,eps)-byplus*dz2frntanh(z-z2,eps)+byminus*dz2frntanh(z-z3,eps)-byminus*dz2frntanh(z-z4,eps)
 dz3bval=dz3bval+byplus*dz3frntanh(z-z1,eps)-byplus*dz3frntanh(z-z2,eps)+byminus*dz3frntanh(z-z3,eps)-byminus*dz3frntanh(z-z4,eps)
 dz4bval=dz4bval+byplus*dz4frntanh(z-z1,eps)-byplus*dz4frntanh(z-z2,eps)+byminus*dz4frntanh(z-z3,eps)-byminus*dz4frntanh(z-z4,eps)

      eps=4.d-3/pi
!!!!!!do n=2,n_u-1
      do n=2,nint(n_u)-1
        zshft=(n-1)*1.d-2
        z1=zshft+0.125d-2; z2=zshft+0.375d-2; z3=zshft+0.625d-2; z4=zshft+0.875d-2
        bval=bval+byplus*frntanh(z-z1,eps)-byplus*frntanh(z-z2,eps)+byminus*frntanh(z-z3,eps)-byminus*frntanh(z-z4,eps)
   dz1bval=dz1bval+byplus*dz1frntanh(z-z1,eps)-byplus*dz1frntanh(z-z2,eps)+byminus*dz1frntanh(z-z3,eps)-byminus*dz1frntanh(z-z4,eps)
   dz2bval=dz2bval+byplus*dz2frntanh(z-z1,eps)-byplus*dz2frntanh(z-z2,eps)+byminus*dz2frntanh(z-z3,eps)-byminus*dz2frntanh(z-z4,eps)
   dz3bval=dz3bval+byplus*dz3frntanh(z-z1,eps)-byplus*dz3frntanh(z-z2,eps)+byminus*dz3frntanh(z-z3,eps)-byminus*dz3frntanh(z-z4,eps)
   dz4bval=dz4bval+byplus*dz4frntanh(z-z1,eps)-byplus*dz4frntanh(z-z2,eps)+byminus*dz4frntanh(z-z3,eps)-byminus*dz4frntanh(z-z4,eps)
      enddo

      eps=4.d-3/pi
      n=n_u
      zshft=(n-1)*1.d-2
      z1=zshft+0.125d-2; z2=zshft+0.375d-2; z3=zshft+0.625d-2; z4=zshft+0.875d-2
      z2=z2-0.125d-2
      bval=bval+byplus*frntanh(z-z1,eps)-byplus*frntanh(z-z2,eps)!+byminus*frntanh(z-z3,eps)-byminus*frntanh(z-z4,eps)
 dz1bval=dz1bval+byplus*dz1frntanh(z-z1,eps)-byplus*dz1frntanh(z-z2,eps)!+byminus*dz1frntanh(z-z3,eps)-byminus*dz1frntanh(z-z4,eps)
 dz2bval=dz2bval+byplus*dz2frntanh(z-z1,eps)-byplus*dz2frntanh(z-z2,eps)!+byminus*dz2frntanh(z-z3,eps)-byminus*dz2frntanh(z-z4,eps)
 dz3bval=dz3bval+byplus*dz3frntanh(z-z1,eps)-byplus*dz3frntanh(z-z2,eps)!+byminus*dz3frntanh(z-z3,eps)-byminus*dz3frntanh(z-z4,eps)
 dz4bval=dz4bval+byplus*dz4frntanh(z-z1,eps)-byplus*dz4frntanh(z-z2,eps)!+byminus*dz4frntanh(z-z3,eps)-byminus*dz4frntanh(z-z4,eps)
!hardwired for no x-dependence:
      b0(2)=bval-0.5d0*y**2*dz2bval+y**4/24.d0*dz4bval
      b0(3)=y*dz1bval-y**3/6.d0*dz3bval
      ierr=0
      return
      end

      subroutine stupakovgetexternalfield(x,y,z,t,e0,b0,ierr)
      implicit none
      real*8 :: x,y,z,t
      real*8, dimension(3) :: e0,b0
      integer :: ierr
      real*8 :: bval,dz1bval,dz2bval,dz3bval,dz4bval
      integer :: n
      integer, parameter :: ndipoles=4
      real*8, parameter :: thet0=2.14d0*2.d0*asin(1.d0)/180.d0
      real*8, parameter :: drouter=1.0d0*cos(thet0)
      real*8, parameter :: drmiddle=0.20d0
!!!!! real*8, parameter :: dd=0.20d0 ! dipole z-length
!!!!! real*8, parameter :: z0=0.10d0  !drift space before first dipole
!!!!! real*8, parameter :: bchic=0.74766244446414987d0
      real*8, parameter :: dd=10.d0  !   0.50d0  !0.50 is for Stupakov/Emma test    !!!!dd=0.20d0 ! dipole z-length
      real*8, parameter :: z0=0.d0 ! 0.10d0   !drift space before first dipole
!     real*8, parameter :: bchic=0.33469850546399954d0  !334... for Stupakov/Emma  !!!   0.74766244446414987d0
      real*8, parameter :: bchic=0.33526427119076008d0
      real*8, parameter, dimension(ndipoles) :: z1=(/z0,      dd+drouter+z0, 2*dd+drmiddle+drouter+z0, 3*dd+drmiddle+2*drouter+z0/) !magnet leading edge
      real*8, parameter, dimension(ndipoles) :: z2=(/dd+z0, 2*dd+drouter+z0, 3*dd+drmiddle+drouter+z0, 4*dd+drmiddle+2*drouter+z0/) !magnet trailing edge
      real*8, parameter, dimension(ndipoles) :: b00=(/+bchic,-bchic,-bchic,+bchic/)  !chicane magnet B field
!!!!!!real*8, parameter :: eps2=0.075d0/(2.d0*asin(1.d0)) !chicane magnet fringe field parameter
      real*8, parameter :: eps2=10.d0*0.0075d0/(2.d0*asin(1.d0)) !chicane magnet fringe field parameter
!hack-------
!     if(x.ne.-1234567.d0)then
!     e0(1:3)=0.d0
!     b0(1:3)=0.d0
!     b0(2)=0.33526427119076008d0   !0.335264T corresponds to rho=1m for 100MeV electron
!     b0(2)=0.33469850546399954d0   !corresponds to rho=1.5m for 150MeV electron
!     ierr=0
!     return
!     endif
!-----------
      e0(1:3)=0.d0
      b0(1:3)=0.d0
      bval=0.d0
      dz1bval=0.d0
      dz2bval=0.d0
      dz3bval=0.d0
      dz4bval=0.d0
      do n=1,ndipoles
        bval=bval+b00(n)*frntanh(z-z1(n),eps2)-b00(n)*frntanh(z-z2(n),eps2)
        dz1bval=dz1bval+b00(n)*dz1frntanh(z-z1(n),eps2)-b00(n)*dz1frntanh(z-z2(n),eps2)
        dz2bval=dz2bval+b00(n)*dz2frntanh(z-z1(n),eps2)-b00(n)*dz2frntanh(z-z2(n),eps2)
        dz3bval=dz3bval+b00(n)*dz3frntanh(z-z1(n),eps2)-b00(n)*dz3frntanh(z-z2(n),eps2)
        dz4bval=dz4bval+b00(n)*dz4frntanh(z-z1(n),eps2)-b00(n)*dz4frntanh(z-z2(n),eps2)
      enddo
!hardwired for no x-dependence:
      b0(2)=bval
!     b0(2)=bval-0.5d0*y**2*dz2bval+y**4/24.d0*dz4bval
!     b0(3)=y*dz1bval-y**3/6.d0*dz3bval
      ierr=0
      return
      end


      subroutine dipolegetexternalfield(x,y,z,t,e0,b0,ierr)
      implicit none
      real*8 :: x,y,z,t
      real*8, dimension(3) :: e0,b0
      integer :: ierr
      real*8, parameter :: bfielddip=0.74766244446414987d0
      e0(1:3)=0.d0
      b0(1:3)=0.d0
      b0(2)=bfielddip
      ierr=0
      return
      end

      real*8 function frntanh(y,eps)
      implicit none
      real*8 :: y,eps
      frntanh = tanh(y/eps)/2.d0 +.5d0
      return
      end

!
      real*8 function dz1frntanh(y,eps)
      implicit none
      real*8 :: y,eps
      if (abs(y/eps).lt.32.d0) then
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
      if (abs(y/eps).lt.32.d0) then
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
      if (abs(y/eps).lt.32.d0) then
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
      if (abs(y/eps).lt.32.d0) then
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

      if(abs(y/eps).ge.32.d0)then
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

      subroutine getquadbfield(xin,yin,zin,barg)
      implicit none
      real*8 :: xin,yin,zin
      real*8, dimension(3) :: barg,b,bp
      integer, parameter :: nquads=9
      real*8, parameter :: pi180=2.d0*asin(1.d0)/180.d0
      real*8, parameter, dimension(nquads) :: xcent=(/-1.1125774E-01,-1.9599797E-01,-3.2267838E-01,-4.3716955E-01,-4.3716955E-01,-4.3716955E-01,-3.2267838E-01,-1.9599797E-01,-1.1125774E-01/)
      real*8, parameter, dimension(nquads) :: xlen=(/0.49,0.49,0.49,0.49,0.49,0.49,0.49,0.49,0.49/)
      real*8, parameter, dimension(nquads) :: epsx=(/0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02/)
      real*8, parameter, dimension(nquads) :: zcent=(/5.4562402E+00,9.2071771E+00,1.4814551E+01,2.1615190E+01,2.2788973E+01,2.3962756E+01,3.0763395E+01,3.6370769E+01,4.0121706E+01/)
      real*8, parameter, dimension(nquads) :: zlen=(/0.714238,0.714238,1.0,0.714238,0.714238,0.714238,1.0,0.714238,0.714238/)
      real*8, parameter, dimension(nquads) :: epsz=(/0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02/)
      real*8, parameter, dimension(nquads) :: tilt=(/+0.022587909d0,+0.022587909d0,+0.022587909d0,0.d0,0.d0,0.d0,-0.022587909d0,-0.022587909d0,-0.022587909d0/)
      real*8, parameter, dimension(nquads) :: bcon=(/0.733338004,-0.485812703,0.284195941,0.396952876,-0.968972588,0.396952876,0.284195941,-0.485812703,0.733338004/) !these are k1 values

      real*8, parameter :: gam0=1.d0+10.d9/0.510998910d6   !NOTE WELL: fix later to get energy from probcons.f90
      real*8, parameter :: gb0=sqrt((gam0+1.d0)*(gam0-1.d0))
      real*8, parameter :: brho=gb0/299792458.d0*0.510998910d6

      real*8 :: x1,x2,z1,z2,ang,x,y,z,xp,yp,zp,b0,bval,dz1bval,dz2bval
!!!   real*8 :: dz3bval,dz4bval
      integer :: n
!     real*8, external :: bump,dz1bump,dz2bump,dz3bump,dz4bump
      b(1:3)=0.d0
      do n=1,nquads
!shift:
        x=xin-xcent(n)
        y=yin
        z=zin-zcent(n)
!rotate:
        ang=tilt(n) !-tilt(n) to try with minus sign
        xp=x*cos(ang)+z*sin(ang)
        yp=y
        zp=-x*sin(ang)+z*cos(ang)
!compute:
        b0=bcon(n)*brho ! multiply by brho to convert k1 values to B field gradient in T/m
        x1=-0.5d0*xlen(n); x2=0.5d0*xlen(n)
        z1=-0.5d0*zlen(n); z2=0.5d0*zlen(n)
        bval=b0*bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
        dz1bval=b0*dz1bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
        dz2bval=b0*dz2bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
!!!     dz3bval=b0*dz3bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
!!!     dz4bval=b0*dz4bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
        bp(1)=bval*yp  -dz2bval*(yp**3+3.d0*xp**2*yp)/12.d0
        bp(2)=bval*xp  -dz2bval*(xp**3+3.d0*yp**2*xp)/12.d0
        bp(3)=0.d0     +dz1bval*xp*yp
!rotate the field back to one aligned with the original frame:
        b(1)=b(1)+(bp(1)*cos(-ang)+bp(3)*sin(-ang))
        b(2)=b(2)+bp(2)
        b(3)=b(3)+(-bp(1)*sin(-ang)+bp(3)*cos(-ang))
      enddo
      barg(1:3)=barg(1:3)+b(1:3)
      return
      end

      subroutine getsextbfield(xin,yin,zin,barg)
      implicit none
      real*8 :: xin,yin,zin
      real*8, dimension(3) :: barg,b,bp
      integer, parameter :: nsexts=4
      real*8, parameter :: pi180=2.d0*asin(1.d0)/180.d0
      real*8, parameter :: gam0=1.d0+10.d9/0.510998910d6   !NOTE WELL: fix later to get energy from probcons.f90
      real*8, parameter :: gb0=sqrt((gam0+1.d0)*(gam0-1.d0))
      real*8, parameter :: brho=gb0/299792458.d0*0.510998910d6
      real*8, parameter, dimension(nsexts) :: xcent=(/-8.9019667E-02,-1.7483542E-01,-1.7483542E-01,-8.9019667E-02/)
      real*8, parameter, dimension(nsexts) :: xlen=(/0.49,0.49,0.49,0.49/)
      real*8, parameter, dimension(nsexts) :: epsx=(/0.02,0.02,0.02,0.02/)
      real*8, parameter, dimension(nsexts) :: zcent=(/4.4718954E+00,8.2704391E+00,3.7307507E+01,4.1106051E+01/)
      real*8, parameter, dimension(nsexts) :: zlen=(/0.762d0,0.250d0,0.250d0,0.762d0/)
      real*8, parameter, dimension(nsexts) :: epsz=(/0.02,0.02,0.02,0.02/)
      real*8, parameter, dimension(nsexts) :: tilt=(/+0.022587909d0,+0.022587909d0,-0.022587909d0,-0.022587909d0/)
!!!   real*8, parameter, dimension(nsexts) :: bcon=(/28.3284704,-0.636171592,-0.636171592,28.3284704/) !these are K2 values
!     real*8, parameter, dimension(nsexts) :: bcon=(/1416.1737483450d0/brho*2.d0,-0.636171592d0,-0.636171592d0,1416.1737483450d0/brho*2.d0/) !these are K2 values
      real*8, parameter, dimension(nsexts) :: bcon=(/42.6987867d0,-8.03335029d0,-8.03335029d0,42.6987867d0/) !these are K2 values


      real*8 :: x1,x2,z1,z2,ang,x,y,z,xp,yp,zp,b0,bval,dz1bval,dz2bval
!!!   real*8 :: dz3bval,dz4bval
      integer :: n
!     real*8, external :: bump,dz1bump,dz2bump,dz3bump,dz4bump
      b(1:3)=0.d0
      do n=1,nsexts
!shift:
        x=xin-xcent(n)
        y=yin
        z=zin-zcent(n)
!rotate:
        ang=tilt(n) !-tilt(n) to try with minus sign
        xp=x*cos(ang)+z*sin(ang)
        yp=y
        zp=-x*sin(ang)+z*cos(ang)
!compute:
        b0=bcon(n)*brho/2.d0
        x1=-0.5d0*xlen(n); x2=0.5d0*xlen(n)
        z1=-0.5d0*zlen(n); z2=0.5d0*zlen(n)
        bval=b0*bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
        dz1bval=b0*dz1bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
        dz2bval=b0*dz2bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
!!!     dz3bval=b0*dz3bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
!!!     dz4bval=b0*dz4bump(zp,z1,z2,epsz(n)) !*bump(xp,x1,x2,epsx(n))
        bp(1)=2.d0*bval*xp*yp
        bp(2)=bval*(xp**2-yp**2)
        bp(3)=0.d0
!rotate the field back to one aligned with the original frame:
        b(1)=b(1)+(bp(1)*cos(-ang)+bp(3)*sin(-ang))
        b(2)=b(2)+bp(2)
        b(3)=b(3)+(-bp(1)*sin(-ang)+bp(3)*cos(-ang))
      enddo
      barg(1:3)=barg(1:3)+b(1:3)
      return
      end
!
      real*8 function bump(s,s1,s2,eps)
      implicit none
      real*8 :: s,s1,s2,eps
!     real*8, external :: frntanh
      bump=(frntanh(s-s1,eps)-frntanh(s-s2,eps))
      return
      end
!
      real*8 function dz1bump(s,s1,s2,eps)
      implicit none
      real*8 :: s,s1,s2,eps
!     real*8, external :: dz1frntanh
      dz1bump=(dz1frntanh(s-s1,eps)-dz1frntanh(s-s2,eps))
      return
      end
!
      real*8 function dz2bump(s,s1,s2,eps)
      implicit none
      real*8 :: s,s1,s2,eps
!     real*8, external :: dz2frntanh
      dz2bump=(dz2frntanh(s-s1,eps)-dz2frntanh(s-s2,eps))
      return
      end
!
      real*8 function dz3bump(s,s1,s2,eps)
      implicit none
      real*8 :: s,s1,s2,eps
!     real*8, external :: dz3frntanh
      dz3bump=(dz3frntanh(s-s1,eps)-dz3frntanh(s-s2,eps))
      return
      end
!
      real*8 function dz4bump(s,s1,s2,eps)
      implicit none
      real*8 :: s,s1,s2,eps
!     real*8, external :: dz4frntanh
      dz4bump=(dz4frntanh(s-s1,eps)-dz4frntanh(s-s2,eps))
      return
      end

end module
