module lwprobcons
      implicit none
      real*8, parameter :: tauini=0.d0     !integrate from time tauini to taufin; the distribution is created at tauini
      real*8, parameter :: taufin=2.0958492136543e-08 !  5.d-9 ! 172.d-9
      integer, parameter :: nsteps=1000 ! 100000   !# of integration steps from tauini to taufin; if too few, LW search might fail
      real*8, parameter :: deltatau=(taufin-tauini)/(nsteps-1)
      real*8, parameter :: tauhistmin=tauini !store the history only if the time is between tauhistmin and tauhistmax
      real*8, parameter :: tauhistmax=taufin
      integer, parameter :: nhistpoints=nint((tauhistmax-tauhistmin)/deltatau)+5  !# of points in the stored history
      real*8, parameter :: tauminus=-2.d-9 !particles assumed to travel in free space (no field) from tauminus to tauini

      integer, parameter :: imax=11 !size of the x output grid
      integer, parameter :: jmax=1  !size of the y output grid
      integer, parameter :: kmax=129 ! 60000 !size of the z output grid
      real*8, parameter :: xoffmin=-5.d-5
      real*8, parameter :: xoffmax= 5.d-5
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
      integer :: npart_gbl=624200000.d0 ! 624220000.d0 !note: code needs to be modified for npart_gbl > 2147483647
      real*8 :: fracinbunch2=0.d0 ! the fraction of npart_gbl particles that is in bunch#2. Set to 0.d0 for a single bunch.
      real*8 :: chrgperbunch=1.0d-10 !  1.d6 * 1.602176634d-19     !results are scaled to this value of charge/bunch (in nC)
      integer :: lwseed=3147228 !random number seed
      integer :: iwritezeroes=0 ! =1 to write LW field even if the value is zero
      integer :: istopafter1fieldcalc=1 !=1 to stop after 1 field calc, otherwise the code keeps running
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
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      ekin=254988456.09d0
      gam0=1.d0+ekin/0.510998910d6       !gamma
      gb0=sqrt((gam0+1.d0)*(gam0-1.d0))  !gamma*beta
      if(myrank.eq.0)write(6,"('gam0,gb0,beta0=',3(1pe19.12,1x))")gam0,gb0,gb0/gam0
      if(myrank.eq.0)write(6,*)'brho=',gb0/299792458.d0*0.510998910d6
      sigmat(1:6,1:6)=0.d0    !beam 2nd moments
      sigmat(1,1)=(1.d-5)**2
      sigmat(3,3)=(1.d-5)**2
      sigmat(5,5)=(1.d-5)**2

      cent(1:6)=0.d0          !beam centroid
      cent(5)=-1.e-15 ! -0.03d0
!     if(indevar.eq.0)then      !if integrating in time; check sign; need to fix if using deviation variables
        cent(6)=gb0
!     endif
!     if(indevar.eq.1)then      !if integrating in z; check sign; need to fix if using deviation variables
!       cent(6)=gam0
!     endif
      return
      end

      subroutine getmatrixandcentroid2(ekin,cent,sigmat)
      use mpi
      implicit none
      real*8 :: ekin !design kinetic energy
      real*8, dimension(6) :: cent       !initial beam centroid
      real*8, dimension(6,6) :: sigmat   !initial beam 2nd moment matrix
      real*8 :: gam0,gb0
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      ekin=254988456.09d0
      gam0=1.d0+ekin/0.510998910d6       !gamma
      gb0=sqrt((gam0+1.d0)*(gam0-1.d0))  !gamma*beta
      if(myrank.eq.0)write(6,"('gam0,gb0,beta0=',3(1pe19.12,1x))")gam0,gb0,gb0/gam0
      if(myrank.eq.0)write(6,*)'brho=',gb0/299792458.d0*0.510998910d6
      sigmat(1:6,1:6)=0.d0    !beam 2nd moments
      sigmat(1,1)=(1.d-7)**2
      sigmat(3,3)=(1.d-7)**2
      sigmat(5,5)=(1.d-7)**2

      cent(1:6)=0.d0          !beam centroid
!     cent(5)=-0.03d0
      cent(5)=-0.03d0 - 1.d-6 !test
!     if(indevar.eq.0)then      !if integrating in time; check sign; need to fix if using deviation variables
        cent(6)=gb0
!     endif
!     if(indevar.eq.1)then      !if integrating in z; check sign; need to fix if using deviation variables
!       cent(6)=gam0
!     endif
      if(cent(6).ne.-12345)then
              write(6,*)'error: code should not get here'
              stop
      endif
      return
      end

!this writes the problem constants so there is a record of them in the output file
      subroutine writeprobcons
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
      if(nstep.ge.1000)idofieldcalc=1
      return
      end
end module lwprobcons
