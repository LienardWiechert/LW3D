      program lw3d
! Parallel 3D Lienard-Wiechert code
! Robert D Ryne -- version 10 September, 2021
      use mpi
      use lwprobcons, only : nsteps,nhistpoints,indevar,refptclhist,refptcltime,  &
      &                      nprop,npart_gbl,chrgperbunch,tauini,taufin,tauminus,lwseed,writeprobcons
      use numerical_distributions, only : initdist
      use integrators, only : propagate,tdriftexact,afterstep,rk4cpy
      use evalroutines, only : getexternalfield,evalt
      use putils, only : pwrite
      use diagnostics, only : moments,profile1d
      implicit none

      real*8, allocatable, dimension(:,:) :: ptcls
      integer :: iseed

      integer :: mimicmprocs !used by 1-proc run to generate identical gendist random #s to mimic a mimicmprocs-proc run

      integer :: mprocs,myrank,mpierr ! # of MPI processes, rank of the proc, error flag
      integer :: npart_lcl,nremain !local variables used to compute the # of particles per MPI proc
      integer :: n,m,ierr,ihertz,iticks0,iticks1,iticks2,irtnerr

      real*8 :: zmin,zmax,hz
      real*8, dimension(3) :: e0,b0

      call MPI_INIT(mpierr)      !initialize MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      if(myrank.eq.0)then
        write(6,*)'******************************************************************'
        write(6,*)'3D Lienard-Wiechert Code version 10 Sept 2021 by Robert Ryne, LBNL'
        write(6,*)'******************************************************************'
      endif

      call writeprobcons
      if(myrank.eq.0)write(6,*)'STARTING SIMULATION. USING ',mprocs,' MPI PROCESSES'

      if(myrank.eq.0)then
        open(unit=91,file='refhist.out',status='unknown')
      endif

      call system_clock(count_rate=ihertz)
      call system_clock(count=iticks0)

      iseed=lwseed

      mimicmprocs=1  !  24 ! 1 ! if mprocs=1, this simulation is serial but with the same random numbers as "mimicmprocs" procs
      if(myrank.eq.0.and.mprocs.eq.1.and.mimicmprocs.ne.1)then
        write(6,*)'mimicmprocs=',mimicmprocs
      endif

!determine the local array size:
      npart_lcl=npart_gbl/mprocs
      nremain=mod(npart_gbl,mprocs)
      if(nremain.ne.0.and.myrank.le.nremain-1)npart_lcl=1+npart_gbl/mprocs
      if(myrank.eq.0)write(6,*)'(PE0) npart_lcl=',npart_lcl
!allocate arrays:
      allocate(ptcls(1:nprop,1:npart_lcl))
      allocate(refptclhist(6,0:nhistpoints,npart_lcl))
      allocate(refptcltime(0:nhistpoints))

! create the initial distribution of particles:
      call initdist(ptcls,iseed,indevar,npart_gbl,mimicmprocs)
      call tdriftexact(tauini,tauminus,ptcls)
      call afterstep(tauminus,ptcls)
      call tdriftexact(tauminus,tauini,ptcls)
      call afterstep(tauini,ptcls)

      if(myrank.eq.0)open(unit=37,file='initp.out',status='unknown')
      if(myrank.eq.0)open(unit=47,file='initzprof.out',status='unknown')
      if(myrank.eq.0)open(unit=20,file='amom.out',status='unknown')
      if(myrank.eq.0)open(unit=24,file='xrms.out',status='unknown')
      if(myrank.eq.0)open(unit=25,file='yrms.out',status='unknown')
      if(myrank.eq.0)open(unit=26,file='zrms.out',status='unknown')
      call pwrite(ptcls,nprop,npart_lcl,37,1,9600,6,15)
      if(myrank.eq.0)close(37)
      call profile1d(ptcls,-1.d0,5,1024,47)
      if(myrank.eq.0)close(47)
      call moments(ptcls,nprop,npart_lcl,0.d0,20,24,25,26,8,0,0,0,0,irtnerr)

      if(myrank.eq.0)then
        open(unit=93,file='extern.out',status='unknown')
        zmin=ptcls(5,1)
        zmax=zmin+299792458.*(taufin-tauini) ! 0.85d0
        hz=(zmax-zmin)/(nsteps-1)
        do n=1,nsteps
          call getexternalfield(0.d0,0.d0,zmin+(n-1)*hz,0.d0,e0,b0,ierr)
          if(abs(b0(2)).lt.1.d-20)cycle
          write(93,'(12(1pe17.10,1x))')0.d0,0.d0,zmin+(n-1)*hz,b0(1),b0(2),b0(3)
        enddo
        close(93)
      endif

! propagate the particles:
      call system_clock(count=iticks1)
      call propagate(ptcls,tauini,taufin,nsteps,indevar)
      call system_clock(count=iticks2)

      if(myrank.eq.0)open(unit=38,file='finalp.out',status='unknown')
      call pwrite(ptcls,nprop,npart_lcl,38,1,9600,6,15)
      if(myrank.eq.0)close(38)
      if(myrank.eq.0)open(unit=48,file='finalzprof.out',status='unknown')
      call profile1d(ptcls,-1.d0,5,1024,48)
      if(myrank.eq.0)close(48)

      if(myrank.eq.0)write(6,*)'initialization time=',dble(iticks1-iticks0)/dble(ihertz),' sec'
      if(myrank.eq.0)write(6,*)'propagation time=',dble(iticks2-iticks1)/dble(ihertz),' sec'
      if(myrank.eq.0)write(6,*)'total execution time=',dble(iticks2-iticks0)/dble(ihertz),' sec'
      call myexit
      end

      subroutine myexit
      use mpi
      implicit none
      integer :: myrank,mpierr
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      if(myrank.eq.0)write(6,*)'job finished'
      call MPI_FINALIZE(mpierr)
      call exit()
      end
