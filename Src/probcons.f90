module lwprobcons
      implicit none
      real*8, parameter :: tauini=0.d0     !integrate from time tauini to taufin; the distribution is created at tauini
      real*8, parameter :: taufin=2.0e-08
      integer, parameter :: nsteps=2000 ! 100000   !# of integration steps from tauini to taufin; if too few, LW search might fail
      real*8, parameter :: deltatau=(taufin-tauini)/(nsteps-1)
      real*8, parameter :: tauhistmin=tauini !store the history only if the time is between tauhistmin and tauhistmax
      real*8, parameter :: tauhistmax=taufin
      integer, parameter :: nhistpoints=nint((tauhistmax-tauhistmin)/deltatau)+5  !# of points in the stored history
      real*8, parameter :: tauminus=-2.d-9 !particles assumed to travel in free space (no field) from tauminus to tauini

      integer, parameter :: imax=1  !imax=3  !size of the x output grid
      integer, parameter :: jmax=1  !size of the y output grid
      integer, parameter :: kmax=129 ! 60000 !size of the z output grid
      real*8, parameter :: xoffmin=0.d0 !-1.d-5
      real*8, parameter :: xoffmax=0.d0 ! 1.d-5
      real*8, parameter :: yoffmin=-0.d0
      real*8, parameter :: yoffmax= 0.d0 
      real*8, parameter :: zoffmin=-5.d-5 ! 50.30d0 !50.400403d0 !5.d-5 !z grid range relative to the z centroid
      real*8, parameter :: zoffmax= 5.d-5 !50.45d0 !50.400404d0 !1.65d-4 
      integer :: irotategrid=1 !must be 0 or 1, this controls the x-z grid rotation
      real*8 :: relativetocentroid=1.d0 !must be 0.d0 or 1.d0

      real*8, allocatable, dimension(:,:,:) :: refptclhist
      real*8, allocatable, dimension(:) :: refptcltime
      integer :: laststoredrefpart=-1  !gets incremented immediately, so the history arrays run from 0 to nhistpoints

      integer, parameter :: indevar=0  !independent variable (0=t, 1=z)
      integer, parameter :: nprop=6    !6-vector of coords and momenta (later version could include loss flag, particle id,etc)
      integer*8 :: npart_gbl=5000000  !changed to integer*8 to enable npart_gbl > 2147483647
!!!!!!real*8 :: fracinbunch2=0.d0 ! the fraction of npart_gbl particles that is in bunch#2. Set to 0.d0 for a single bunch.
      real*8 :: chrgperbunch=1.0d-9 ! results are scaled to this value of charge/bunch (in Coulomb)
      integer :: lwseed=3147228 !random number seed
      integer :: iwritezeroes=0 ! =1 to write LW field even if the value is zero
      integer :: istopafter1fieldcalc=0 !=1 to stop after 1 field calc, otherwise the code keeps running
      logical :: rowmajor=.true.
      logical :: blankafter2dblock=.true.

save
contains

      subroutine getmatrixandcentroid(ekin,cent,sigmat)
      use mpi
      implicit none
      real*8 :: ekin !design kinetic energy
      real*8, dimension(6) :: cent       !initial beam centroid
      real*8, dimension(6,6) :: sigmat   !initial beam 2nd moment matrix
      real*8 :: gam0,gb0
      real*8 :: alfx,betx,gamx,epsx,sigx
      real*8 :: alfy,bety,gamy,epsy,sigy
      real*8 :: alfz,betz,gamz,epsz,sigz,h,r56,zpz,sigpz
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      ekin=70.d6
      gam0=1.d0+ekin/0.510998910d6       !gamma
      gb0=sqrt((gam0+1.d0)*(gam0-1.d0))  !gamma*beta
      if(myrank.eq.0)write(6,"('gam0,gb0,beta0=',3(1pe19.12,1x))")gam0,gb0,gb0/gam0
      if(myrank.eq.0)write(6,*)'brho=',gb0/299792458.d0*0.510998910d6

!beam centroid:
      cent(1:6)=0.d0
      cent(5)=-0.2d0 !simulation begins at -20cm, outside the dipole fringe field
      cent(6)=gb0    !the 6th variable is gamma*beta_z; set the centroid to the design value

!2nd moment matrix:
!initialize 2nd moment matrix to zero:
      sigmat(1:6,1:6)=0.d0
!fill in the other nonzero matrix elements:
!for example, upright in phase space with zero emittance in all planes:
      sigmat(1,1)=(1.d-5)**2
      sigmat(3,3)=(1.d-5)**2
      sigmat(5,5)=(1.d-5)**2

!here is old code in an attempt to initialize sigmat to achieve compression:
!     alfx= 0.75784343
!     betx= 3.17323883
!     gamx= (1.d0+alfx**2)/betx
!     alfy= 0.771049709
!     bety= 4.001700150
!     gamy=  (1.d0+alfy**2)/bety
!     sigx=21.75477d-6
!     sigy=24.4306d-6
!     epsx=sigx**2/betx
!     epsy=sigy**2/bety
!static to dynamic (i.e. unnormalized to normalized):
!     epsx=epsx*gam0
!     epsy=epsy*gam0
!     betx=betx/gam0
!     bety=bety/gam0
!     gamx=gamx*gam0
!     gamy=gamy*gam0
!
!     sigmat(1,1)= betx*epsx
!     sigmat(1,2)=-alfx*epsx*0.d0
!     sigmat(2,2)= gamx*epsx*1.d-18*0.d0
!     sigmat(2,1)=sigmat(1,2)
!     if(myrank.eq.0)then
!       write(6,*)'sigmat11=',sigmat(1,1)
!       write(6,*)'sigmat12=',sigmat(1,2)
!       write(6,*)'sigmat22=',sigmat(2,2)
!     endif
!     sigmat(3,3)= bety*epsy
!     sigmat(3,4)=-alfy*epsy*0.d0
!     sigmat(4,4)= gamy*epsy*1.d-18*0.d0
!     sigmat(4,3)=sigmat(3,4)
!     if(myrank.eq.0)then
!       write(6,*)'sigmat33=',sigmat(3,3)
!       write(6,*)'sigmat34=',sigmat(3,4)
!       write(6,*)'sigmat44=',sigmat(4,4)
!     endif
!     sigz=40.d-6  ! zrms=40 micron
!     r56=-0.43d0
!     h=-1.d0/r56 !chirp set for maximum compression
!     epsz=1.d-6
!     sigmat(5,5)=sigz**2
!     sigmat(5,6)=h*sigz**2*9.d3
!     sigmat(6,5)=sigmat(5,6)
!     sigmat(6,6)=(epsz**2+sigmat(5,6)**2)/sigmat(5,5)
!
!!!!  sigz=40.d-6
!!!!  sigpz=0.1d0
!!!!  zpz=3.87298334620742e-06
!!!!  sigmat(5,5)=sigz**2
!!!!  sigmat(5,6)=-zpz
!!!!  sigmat(6,6)=sigpz**2
!!!!  sigmat(6,5)=sigmat(5,6)
!
!     if(myrank.eq.0)then
!       write(6,*)'sigmat55=',sigmat(5,5)
!       write(6,*)'sigmat56=',sigmat(5,6)
!       write(6,*)'sigmat66=',sigmat(6,6)
!     endif
      return
      end

      subroutine writeprobcons
!this writes the problem constants so there is a record of them in the output file
      use mpi
      implicit none
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      if(myrank.eq.0)write(6,*)' '
      if(myrank.eq.0)write(6,*)'PROBLEM CONSTANTS:'
      if(myrank.eq.0)write(6,*)'tauini,taufin=',tauini,taufin
      if(myrank.eq.0)write(6,*)'nsteps=',nsteps
      if(myrank.eq.0)write(6,*)'deltatau=',deltatau
      if(myrank.eq.0)write(6,*)'tauhistmin,tauhistmax=',tauhistmin,tauhistmax
      if(myrank.eq.0)write(6,*)'nhistpoints=',nhistpoints
      if(myrank.eq.0)write(6,*)'tauminus=',tauminus
      if(myrank.eq.0)write(6,*)'imax,kmax=',imax,kmax
      if(myrank.eq.0)write(6,*)'xoffmin,xoffmax=',xoffmin,xoffmax
      if(myrank.eq.0)write(6,*)'zoffmin,zoffmax=',zoffmin,zoffmax
      if(myrank.eq.0)write(6,*)'npart_gbl=',npart_gbl
      if(myrank.eq.0)write(6,*)'chrgperbunch=',chrgperbunch
      if(myrank.eq.0)write(6,*)'lwseed=',lwseed
      if(myrank.eq.0)write(6,*)'iwritezeroes=',iwritezeroes
      if(myrank.eq.0)write(6,*)'istopafter1fieldcalc=',istopafter1fieldcalc
      if(myrank.eq.0)write(6,*)'irotategrid=',irotategrid
      if(myrank.eq.0)write(6,*)'relativetocentroid=',relativetocentroid
      if(myrank.eq.0)write(6,*)' '
      if(irotategrid.ne.0.and.irotategrid.ne.1)then
        if(myrank.eq.0)write(6,*)'error: irotategrid must be 0 or 1'
        stop
      endif
      if(relativetocentroid.ne.0.d0.and.relativetocentroid.ne.1.d0)then
        if(myrank.eq.0)write(6,*)'error: relativetocentroid must be 0.d0 or 1.d0'
        stop
      endif
      if(irotategrid.eq.1.and.relativetocentroid.eq.0.d0)then
       if(myrank.eq.0)write(6,*)'error: irotategrid=1, so relativetocentroid 0.d0 not allowed; change 1 to 0 or change 0.d0 to 1.d0'
       stop
      endif
      return
      end

!this is where the user specifies the test to decide whether or not to perform a LW field calc at the end of a step
      subroutine testtoperformfieldcalc(cent,nstep,idofieldcalc)
      implicit none
      real*8, dimension(*) :: cent !centroid
      integer :: nstep,idofieldcalc
      idofieldcalc=0
!     if(cent(5).gt.5.d0)idofieldcalc=1
!     if(nstep.ge.1000)idofieldcalc=1
      if(nstep.le.100 .and. mod(nstep,25).eq.0)idofieldcalc=1
      return
      end
end module lwprobcons
