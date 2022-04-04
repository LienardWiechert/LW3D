module integrators

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains

      subroutine propagate(ptcls,tini,tfin,nsteps,indevar)
      use mpi
      use evalroutines, only : evalt
      implicit none
      real*8, dimension(:,:) :: ptcls !particle array
      real*8 :: tini,tfin
      integer :: nsteps,indevar
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      if(indevar.eq.0)call adams11(ptcls,tini,tfin,nsteps,evalt)
      return
      end

      subroutine rk4(y,tini,tfin,nsteps,eval)
      use mpi
      use lienwiech, only : getfield_atallpts
      use putils, only : pwrite
      use diagnostics, only : getcentroid,openfile
      use lwprobcons, only : tauhistmin,tauhistmax,refptcltime,laststoredrefpart,testtoperformfieldcalc,istopafter1fieldcalc,&
     &                       testtowriteptcls,npart_gbl
      implicit none
      real*8, dimension(:,:) :: y
      real*8 :: tini,tfin
      integer :: nsteps,mlo,mhi,mlomin,mhimax,gblmlomin,gblmhimax
      real*8, allocatable, dimension(:,:) :: yt,f,a,b,c,d
      integer :: n1,np,n,irtn,k,i,m,lth,idofieldcalc,iwriteptcls
      real*8 :: h,t,tt,xdum,zdum
      integer :: mprocs,myrank,mpierr,ierr
      character*32 :: pname
      integer, parameter :: nfiledigits=4
      integer, parameter :: nptclswritten=10000

      real*8, dimension(6) :: av  ! fix (6) later !!!!!!

      interface
          subroutine eval(t,y,f)
            real*8 :: t
            real*8, dimension(:,:) :: y,f
          end subroutine
      end interface

      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      n1=ubound(y,1)
      np=ubound(y,2)
      if(myrank.eq.0)write(6,*)'(rk4) n1,np=',n1,np
      allocate(yt(n1,np),f(n1,np),a(n1,np),b(n1,np),c(n1,np),d(n1,np))
      h=(tfin-tini)/nsteps
      t=tini
!     call moments(y,n1,np,t,20,24,25,26,9,0,0,0,0,irtn)
      do n=1,nsteps
        if(laststoredrefpart.eq.-1)then
          if(t+h.ge.tauhistmin)then
            call afterstep(t,y) !store the initial values in the history arrays
          endif
        endif
        if(myrank.eq.0.and.mod(n,1000).eq.0)write(6,*)'STARTING STEP ',n
        call eval(t,y,f)
        a(1:n1,1:np)=h*f(1:n1,1:np)
        yt(1:n1,1:np)=y(1:n1,1:np)+0.5d0*a(1:n1,1:np)
        tt=t+0.5d0*h
        call eval(tt,yt,f)
        b(1:n1,1:np)=h*f(1:n1,1:np)
        yt(1:n1,1:np)=y(1:n1,1:np)+0.5d0*b(1:n1,1:np)
        tt=t+0.5d0*h
        call eval(tt,yt,f)
        c(1:n1,1:np)=h*f(1:n1,1:np)
        yt(1:n1,1:np)=y(1:n1,1:np)+c(1:n1,1:np)
        tt=t+h
        call eval(tt,yt,f)
        d(1:n1,1:np)=h*f(1:n1,1:np)
        y(1:n1,1:np)=y(1:n1,1:np)+(a(1:n1,1:np)+2.d0*b(1:n1,1:np)+2.d0*c(1:n1,1:np)+d(1:n1,1:np))/6.d0
        t=tini+n*h
!       call moments(y,n1,np,t,20,24,25,26,9,0,0,0,0,irtn)
        if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y) !store latest values in the history arrays

!       at this point the particle history has been stored up to time t.
!       getfieldsatpoint() can be called to compute the LW fields.

        call getcentroid(y,av,n1,np) !do this to compute the LW fields when
        call testtoperformfieldcalc(av,n,idofieldcalc)
        if(idofieldcalc.eq.1)then
          if(myrank.eq.0)write(6,*)'Step ',n,' complete; computing fields using history values'
          if(myrank.eq.0)write(6,*)'refptcltime(0)=',refptcltime(0)
          if(myrank.eq.0)write(6,*)'refptcltime(',laststoredrefpart,')=',refptcltime(laststoredrefpart)
          call getfield_atallpts(t,av,np,n)
          if(myrank.eq.0)write(6,*)'Done computing fields.'
          if(istopafter1fieldcalc.eq.1)exit !test to exit from loop
        endif
        call testtowriteptcls(av,n,iwriteptcls)
      if(iwriteptcls.eq.1)then
        pname='ptcls'
        call openfile(pname,490,i,nfiledigits)
        call pwrite(y,n1,np,490,1,nptclswritten,6,15)
        close(490)
!       pname='zprof'
!       call openfile(pname,490,i,nfiledigits)
!       call profile1d(y,-1.d0,5,1024,490)
!       close(490)
      endif
      enddo
      return
      end

      subroutine adams11(y,tini,tfin,nsteps,eval)
      use mpi
      use lienwiech, only : getfield_atallpts
      use putils, only : pwrite
      use diagnostics, only : getcentroid,moments,openfile
      use lwprobcons, only : tauhistmin,tauhistmax,refptcltime,laststoredrefpart,testtoperformfieldcalc,istopafter1fieldcalc,&
     &                       testtowriteptcls,npart_gbl,nstepcounter
      implicit none
      real*8 :: tini,tfin
      integer :: nsteps
      real*8, dimension(:,:) :: y
      real*8, allocatable, dimension(:,:) :: yp,yc,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,yminusav
      real*8, dimension(10) :: am,bm
      real*8, parameter, dimension(10) :: a(1:10)= &
     &  (/57281.d0,-583435.d0,2687864.d0,-7394032.d0,13510082.d0,-17283646.d0,16002320.d0,-11271304.d0,9449717.d0,2082753.d0/)
      real*8, parameter, dimension(10) :: b(1:10)= &
     &  (/-2082753.d0,20884811.d0,-94307320.d0,252618224.d0,-444772162.d0,538363838.d0,-454661776.d0,265932680.d0,-104995189.d0,&
     &     30277247.d0/)
      integer :: n1,np,nsa,i,j,idofieldcalc,iwriteptcls,k
      real*8 :: h,t,tint,hdiv
      integer :: mprocs,myrank,mpierr,irtnerr
      character*32 :: pname
      integer, parameter :: nfiledigits=4
      integer, parameter :: nptclswritten=10000

      real*8, dimension(6) :: av  ! fix (6) later !!!!!!

      interface
          subroutine eval(t,y,f)
            real*8 :: t
            real*8, dimension(:,:) :: y,f
          end subroutine
      end interface

      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      n1=ubound(y,1)
      np=ubound(y,2)
      allocate(yp(n1,np),yc(n1,np))
      allocate(f1(n1,np),f2(n1,np),f3(n1,np),f4(n1,np),f5(n1,np),f6(n1,np),f7(n1,np),f8(n1,np),f9(n1,np),f10(n1,np),f11(n1,np))
      allocate(yminusav(n1,np))

      call getcentroid(y,av,n1,np)
      if(myrank.eq.0)write(6,*)'Entered Adams routine...'

      h=(tfin-tini)/nsteps
      t=tini
      if(laststoredrefpart.eq.-1)then
        if(t+h.ge.tauhistmin)then
          call afterstep(t,y) !store the initial values in the history arrays
        endif
      endif
      nstepcounter=0
      call eval(t,y,f1)
      call rk4cpy(y,tini,tini+h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=1
      call eval(t,y,f2)
      call rk4cpy(y,tini+h,tini+2.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=2
      call eval(t,y,f3)
      call rk4cpy(y,tini+2.d0*h,tini+3.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=3
      call eval(t,y,f4)
      call rk4cpy(y,tini+3.d0*h,tini+4.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=4
      call eval(t,y,f5)
      call rk4cpy(y,tini+4.d0*h,tini+5.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=6
      call eval(t,y,f6)
      call rk4cpy(y,tini+5.d0*h,tini+6.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=6
      call eval(t,y,f7)
      call rk4cpy(y,tini+6.d0*h,tini+7.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=7
      call eval(t,y,f8)
      call rk4cpy(y,tini+7.d0*h,tini+8.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=8
      call eval(t,y,f9)
      call rk4cpy(y,tini+8.d0*h,tini+9.d0*h,10,eval)
      t=t+h; if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y)
      nstepcounter=9
      call eval(t,y,f10)
      nsa=nsteps-9
      hdiv=h/7257600.0d+00

      am(1:10)=hdiv*a(1:10)
      bm(1:10)=hdiv*b(1:10)
      tint=t
      call moments(y,n1,np,tint,20,24,25,26,8,0,0,0,0,irtnerr)

      do 100 i=1,nsa
      nstepcounter=i+9
      yp(1:n1,1:np)=y(1:n1,1:np)+bm(1)*f1(1:n1,1:np)+bm(2)*f2(1:n1,1:np)+bm(3)*f3(1:n1,1:np) &
     & +bm(4)*f4(1:n1,1:np)+bm(5)*f5(1:n1,1:np)+bm(6)*f6(1:n1,1:np)+bm(7)*f7(1:n1,1:np) &
     & +bm(8)*f8(1:n1,1:np)+bm(9)*f9(1:n1,1:np)+bm(10)*f10(1:n1,1:np)
      call eval(t+h,yp,f11)
      yp(1:n1,1:np)=y(1:n1,1:np)+am(1)*f2(1:n1,1:np)+am(2)*f3(1:n1,1:np)+am(3)*f4(1:n1,1:np)+am(4)*f5(1:n1,1:np)+ &
     &                   am(5)*f6(1:n1,1:np)+am(6)*f7(1:n1,1:np)+am(7)*f8(1:n1,1:np)+am(8)*f9(1:n1,1:np)+am(9)*f10(1:n1,1:np)
      yc(1:n1,1:np)=yp(1:n1,1:np)+am(10)*f11(1:n1,1:np)
      call eval(t+h,yc,f11)
      y(1:n1,1:np)=yp(1:n1,1:np)+am(10)*f11(1:n1,1:np)
      f1(1:n1,1:np)=f2(1:n1,1:np)
      f2(1:n1,1:np)=f3(1:n1,1:np)
      f3(1:n1,1:np)=f4(1:n1,1:np)
      f4(1:n1,1:np)=f5(1:n1,1:np)
      f5(1:n1,1:np)=f6(1:n1,1:np)
      f6(1:n1,1:np)=f7(1:n1,1:np)
      f7(1:n1,1:np)=f8(1:n1,1:np)
      f8(1:n1,1:np)=f9(1:n1,1:np)
      f9(1:n1,1:np)=f10(1:n1,1:np)
      f10(1:n1,1:np)=f11(1:n1,1:np)
      t=tint+i*h
      call moments(y,n1,np,t,20,24,25,26,8,0,0,0,0,irtnerr)
      if(t.ge.tauhistmin .and. t.le.tauhistmax)call afterstep(t,y) !store latest values in the history arrays
      call getcentroid(y,av,n1,np)
      call testtoperformfieldcalc(av,i+9,idofieldcalc)
      if(idofieldcalc.eq.1)then
        if(myrank.eq.0)write(6,*)'Step ',i+9,' complete; computing fields at t=',t
        call getfield_atallpts(t,av,np,i+9)
        if(myrank.eq.0)write(6,*)'Done computing fields.'
        if(istopafter1fieldcalc.eq.1)exit !test to exit from loop
      endif
      call testtowriteptcls(av,i+9,iwriteptcls)
      if(iwriteptcls.eq.1)then
        pname='ptcls'
        call openfile(pname,490,i+9,nfiledigits)
        if(myrank.eq.0)write(6,*)'particle output to file ',i+9,' with <z>=',av(5)
        do k=1,np
          yminusav(1,k)=y(1,k)-av(1)
          yminusav(2,k)=y(2,k)-av(2)
          yminusav(3,k)=y(3,k)-av(3)
          yminusav(4,k)=y(4,k)-av(4)
          yminusav(5,k)=y(5,k)-av(5)
          yminusav(6,k)=y(6,k)-av(6)
        enddo
        call pwrite(yminusav,n1,np,490,1,nptclswritten,6,15)
        close(490)
!       pname='zprof'
!       call openfile(pname,490,i+9,nfiledigits)
!       call profile1d(y,-1.d0,5,1024,490)
!       close(490)
      endif
  100 continue
      return
      end

      subroutine afterstep(t,y)
      use mpi
      use lwprobcons, only : refptclhist,refptcltime,laststoredrefpart,nhistpoints
      implicit none
      real*8 :: t
      real*8, dimension(:,:) :: y
      real*8, dimension(6) :: av
      integer :: n1,np,n,nstep
      integer :: myrank,mpierr

      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      n1=ubound(y,1)
      np=ubound(y,2)
      nstep=laststoredrefpart+1 !assumes initilized to -1 in module  !!!!!!

      if(nstep.eq.nhistpoints+1)then
        if(myrank.eq.0)write(6,*)'shifting the arrays of stored quantities'
        do n=1,nhistpoints
          refptcltime(n-1)=refptcltime(n)
          refptclhist(1:6,n-1,1:np)=refptclhist(1:6,n,1:np)
        enddo
        nstep=nhistpoints
      endif

      refptcltime(nstep)=t

      do n=1,np
      refptclhist(1:6,nstep,n)=y(1:6,n)
      enddo
      if(myrank.eq.0)write(91,'(7(1pe16.9,1x))')refptcltime(nstep),refptclhist(1:6,nstep,1)
      laststoredrefpart=nstep
      return
      end


      subroutine rk4cpy(y,tini,tfin,nsteps,eval)
      use mpi
      implicit none
      real*8, dimension(:,:) :: y
      real*8 :: tini,tfin
      integer :: nsteps,mlo,mhi,mlomin,mhimax,gblmlomin,gblmhimax
      real*8, allocatable, dimension(:,:) :: yt,f,a,b,c,d
      integer :: n1,np,n,irtn,k,i,m,lth
      real*8 :: h,t,tt,xdum,zdum
      integer :: mprocs,myrank,mpierr,ierr

      real*8, dimension(6) :: av  ! fix (6) later !!!!!!

      interface
          subroutine eval(t,y,f)
            real*8 :: t
            real*8, dimension(:,:) :: y,f
          end subroutine
      end interface

      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      n1=ubound(y,1)
      np=ubound(y,2)
      allocate(yt(n1,np),f(n1,np),a(n1,np),b(n1,np),c(n1,np),d(n1,np))
      h=(tfin-tini)/nsteps
      t=tini
      do n=1,nsteps
        call eval(t,y,f)
        a(1:n1,1:np)=h*f(1:n1,1:np)
        yt(1:n1,1:np)=y(1:n1,1:np)+0.5d0*a(1:n1,1:np)
        tt=t+0.5d0*h
        call eval(tt,yt,f)
        b(1:n1,1:np)=h*f(1:n1,1:np)
        yt(1:n1,1:np)=y(1:n1,1:np)+0.5d0*b(1:n1,1:np)
        tt=t+0.5d0*h
        call eval(tt,yt,f)
        c(1:n1,1:np)=h*f(1:n1,1:np)
        yt(1:n1,1:np)=y(1:n1,1:np)+c(1:n1,1:np)
        tt=t+h
        call eval(tt,yt,f)
        d(1:n1,1:np)=h*f(1:n1,1:np)
        y(1:n1,1:np)=y(1:n1,1:np)+(a(1:n1,1:np)+2.d0*b(1:n1,1:np)+2.d0*c(1:n1,1:np)+d(1:n1,1:np))/6.d0
        t=tini+n*h
      enddo
      return
      end

      subroutine tdriftexact(tauini,taufin,y)
      implicit none
      real*8 :: tauini,taufin
      real*8, dimension(:,:) :: y
      real*8 :: gam2,gam
      integer :: n1,np,n
      real*8, parameter :: clite=299792458.d0
      n1=ubound(y,1)
      np=ubound(y,2)
      do n=1,np
        gam2=1.d0+y(2,n)**2+y(4,n)**2+y(6,n)**2
        gam=sqrt(gam2)
        y(1,n)=y(1,n)+(taufin-tauini)*y(2,n)*clite/gam
        y(3,n)=y(3,n)+(taufin-tauini)*y(4,n)*clite/gam
        y(5,n)=y(5,n)+(taufin-tauini)*y(6,n)*clite/gam
      enddo
      return
      end

end module
