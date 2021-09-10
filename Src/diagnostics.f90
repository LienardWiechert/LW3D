module diagnostics

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains
      subroutine moments(coordsin,n1,nraysp,tors,nfile0,nfile1,nfile2,nfile3,nprecision,ncross,includepi,ncent,icorr,irtrnerr)
! routine that does DOES rotate into a frame aligned with the direction of motion, then computes and writes 2nd moments
! compute 2nd moment matrix and some auxilliary quantites (e.g. emittances)
!   tors = the independent variable (usually "t" or "s")
!   nprecision = precision at which to write data
!   ncross = 0/1 to report cross-terms (<xpx>,...) as is/as ratios
!   includepi = 0/1 to no/yes divide emittances by pi
!   ncent = 0/1 to remove/keep beam centroid in 2nd moments computation
!   icorr = 0/1 to include/omit correction term to deal with roundoff
      use mpi
      implicit none
      integer :: n1,nraysp
      real*8, dimension(n1,*) :: coordsin
      real*8, dimension(n1,nraysp) :: coords
      real*8 :: tors
      integer :: nfile0,nfile1,nfile2,nfile3,nprecision,ncross,includepi,ncent,icorr,irtrnerr
      real*8 rgood
      real*8 :: pi
!     real*8 :: wex,wey,wet
!
      logical msk(nraysp)
      real*8 avloc(6),av(6),rv(6)     !local,global first moments (i.e. average values); rv are the "raw" values
      real*8 vecm(21),gvecm(21) !local,global vector of 2nd moments
      real*8 correc(6)
      real*8 gcorr(6)
      integer ierr,n,nlength,nparm
      integer i,j,k
      integer ngood,ngoodloc
      real*8 den1,den2
      real*8 xbar,ybar,tbar,pxbar,pybar,ptbar
      real*8 xrms,yrms,trms,pxrms,pyrms,ptrms
      real*8 sq1,sq2,sq3,sq4,sq5,sq6
      real*8 xpx,ypy,tpt,xpxfac,ypyfac,tptfac
      real*8 epsx,epsy,epst,epsx2,epsy2,epst2
      character*17 :: myformat='(35(1x,1pe18.11))'
      character*2 :: a2
      real*8 :: theta,x00,z00
      integer :: mprocs,myrank,mpierr
!
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      irtrnerr=0
      if(ncent.ne.0)then
        write(6,*)'error: this version must have ncent=0'
        call myexit
      endif

      coords(1:n1,1:nraysp)=coordsin(1:n1,1:nraysp) ! 17 June 2018

!
! set up particle mask
      msk(1:nraysp)=.true. !in the future this can be used to omit lost particles
!     call flaglostparticles(coords) !assumes phase cut is for array in radians
!     call setlostparticlemask(msk)
      ngoodloc=count(msk(1:nraysp))
! count # of "good" particles:
      call MPI_ALLREDUCE(ngoodloc,ngood,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!     ngood=ngoodloc
      if(ngood.eq.0)then
        write(6,*)'moments error: ngood=0'
        irtrnerr=-1
        return
      endif
      rgood=ngood
      den1=1.d0/(1.d0*ngood)  ! 17 June 2018
      den2=den1*den1
!     if(myrank.eq.0)write(6,*)'rgood=',rgood
!
!   rotate so the direction of motion is aligned with the z-direction.
!   this version ignores the y-plane.
!
!   compute the beam centroid:
!     avloc(1)=sum(coords(1,1:nraysp),msk(1:nraysp))*den1
!     avloc(2)=sum(coords(2,1:nraysp),msk(1:nraysp))*den1
!     avloc(3)=sum(coords(3,1:nraysp),msk(1:nraysp))*den1
!     avloc(4)=sum(coords(4,1:nraysp),msk(1:nraysp))*den1
!     avloc(5)=sum(coords(5,1:nraysp),msk(1:nraysp))*den1
!     avloc(6)=sum(coords(6,1:nraysp),msk(1:nraysp))*den1
      avloc(1:6)=0.d0
      do n=1,nraysp
        avloc(1:6)=avloc(1:6)+coords(1:6,n)
      enddo
      avloc(1:6)=avloc(1:6)*den1

!     if(myrank.eq.0)write(6,*)'avloc:'
!     if(myrank.eq.0)write(6,*)avloc(1)
!     if(myrank.eq.0)write(6,*)avloc(2)
!     if(myrank.eq.0)write(6,*)avloc(3)
!     if(myrank.eq.0)write(6,*)avloc(4)
!     if(myrank.eq.0)write(6,*)avloc(5)
!     if(myrank.eq.0)write(6,*)avloc(6)
      call MPI_ALLREDUCE(avloc,rv,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!     rv(1:6)=avloc(1:6) !keep the "raw" centroid values to print later
!     if(myrank.eq.0)write(6,*)'rv:'
!     if(myrank.eq.0)write(6,*)rv(1)
!     if(myrank.eq.0)write(6,*)rv(2)
!     if(myrank.eq.0)write(6,*)rv(3)
!     if(myrank.eq.0)write(6,*)rv(4)
!     if(myrank.eq.0)write(6,*)rv(5)
!     if(myrank.eq.0)write(6,*)rv(6)
!
!   remove the x- and z- centroid:
      x00=rv(1)
      z00=rv(5)
      coords(1,1:nraysp)=coords(1,1:nraysp)-x00
      coords(5,1:nraysp)=coords(5,1:nraysp)-z00
!
!   since this version uses gbx and gbz, the angle is atan2(<gbx>,<gbz>)
      theta=-atan2(rv(2),rv(6)) !why is there a minus sign???
!     write(6,*)'theta=',theta
!   rotate:
      call xzrot(theta,coords,msk,n1,nraysp)
!
!   re-compute beam centroid
      avloc(1)=sum(coords(1,1:nraysp),msk(1:nraysp))*den1
      avloc(2)=sum(coords(2,1:nraysp),msk(1:nraysp))*den1
      avloc(3)=sum(coords(3,1:nraysp),msk(1:nraysp))*den1
      avloc(4)=sum(coords(4,1:nraysp),msk(1:nraysp))*den1
      avloc(5)=sum(coords(5,1:nraysp),msk(1:nraysp))*den1
      avloc(6)=sum(coords(6,1:nraysp),msk(1:nraysp))*den1
      call MPI_ALLREDUCE(avloc,av,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!     av(1:6)=avloc(1:6)
!
!correction term to deal with roundoff issues:
      if(icorr.eq.0)then
        correc(1)=sum(coords(1,1:nraysp)-av(1),msk(1:nraysp))*den1
        correc(2)=sum(coords(2,1:nraysp)-av(2),msk(1:nraysp))*den1
        correc(3)=sum(coords(3,1:nraysp)-av(3),msk(1:nraysp))*den1
        correc(4)=sum(coords(4,1:nraysp)-av(4),msk(1:nraysp))*den1
        correc(5)=sum(coords(5,1:nraysp)-av(5),msk(1:nraysp))*den1
        correc(6)=sum(coords(6,1:nraysp)-av(6),msk(1:nraysp))*den1
        call MPI_ALLREDUCE(correc,gcorr,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  !17 June 2018
        correc(1:6)=gcorr(1:6) !17 June 2018
      else
        correc(1:6)=0.d0
      endif
!
!   give readable names to coordinate entries
      xbar=av(1)
      ybar=av(3)
      tbar=av(5)
      pxbar=av(2)
      pybar=av(4)
      ptbar=av(6)
      if(ncent.eq.1)then  !if ncent=1, don't remove centroid from 2nd moment calc
        av(1:6)=0.d0
      endif
!   compute first-order and second-order moments
      vecm(1)= sum((coords(1,1:nraysp)-av(1))*(coords(1,1:nraysp)-av(1)),msk(1:nraysp))*den1-correc(1)*correc(1)
      vecm(2)= sum((coords(1,1:nraysp)-av(1))*(coords(2,1:nraysp)-av(2)),msk(1:nraysp))*den1-correc(1)*correc(2)
      vecm(3)= sum((coords(1,1:nraysp)-av(1))*(coords(3,1:nraysp)-av(3)),msk(1:nraysp))*den1-correc(1)*correc(3)
      vecm(4)= sum((coords(1,1:nraysp)-av(1))*(coords(4,1:nraysp)-av(4)),msk(1:nraysp))*den1-correc(1)*correc(4)
      vecm(5)= sum((coords(1,1:nraysp)-av(1))*(coords(5,1:nraysp)-av(5)),msk(1:nraysp))*den1-correc(1)*correc(5)
      vecm(6)= sum((coords(1,1:nraysp)-av(1))*(coords(6,1:nraysp)-av(6)),msk(1:nraysp))*den1-correc(1)*correc(6)
      vecm(7)= sum((coords(2,1:nraysp)-av(2))*(coords(2,1:nraysp)-av(2)),msk(1:nraysp))*den1-correc(2)*correc(2)
      vecm(8)= sum((coords(2,1:nraysp)-av(2))*(coords(3,1:nraysp)-av(3)),msk(1:nraysp))*den1-correc(2)*correc(3)
      vecm(9)= sum((coords(2,1:nraysp)-av(2))*(coords(4,1:nraysp)-av(4)),msk(1:nraysp))*den1-correc(2)*correc(4)
      vecm(10)=sum((coords(2,1:nraysp)-av(2))*(coords(5,1:nraysp)-av(5)),msk(1:nraysp))*den1-correc(2)*correc(5)
      vecm(11)=sum((coords(2,1:nraysp)-av(2))*(coords(6,1:nraysp)-av(6)),msk(1:nraysp))*den1-correc(2)*correc(6)
      vecm(12)=sum((coords(3,1:nraysp)-av(3))*(coords(3,1:nraysp)-av(3)),msk(1:nraysp))*den1-correc(3)*correc(3)
      vecm(13)=sum((coords(3,1:nraysp)-av(3))*(coords(4,1:nraysp)-av(4)),msk(1:nraysp))*den1-correc(3)*correc(4)
      vecm(14)=sum((coords(3,1:nraysp)-av(3))*(coords(5,1:nraysp)-av(5)),msk(1:nraysp))*den1-correc(3)*correc(5)
      vecm(15)=sum((coords(3,1:nraysp)-av(3))*(coords(6,1:nraysp)-av(6)),msk(1:nraysp))*den1-correc(3)*correc(6)
      vecm(16)=sum((coords(4,1:nraysp)-av(4))*(coords(4,1:nraysp)-av(4)),msk(1:nraysp))*den1-correc(4)*correc(4)
      vecm(17)=sum((coords(4,1:nraysp)-av(4))*(coords(5,1:nraysp)-av(5)),msk(1:nraysp))*den1-correc(4)*correc(5)
      vecm(18)=sum((coords(4,1:nraysp)-av(4))*(coords(6,1:nraysp)-av(6)),msk(1:nraysp))*den1-correc(4)*correc(6)
      vecm(19)=sum((coords(5,1:nraysp)-av(5))*(coords(5,1:nraysp)-av(5)),msk(1:nraysp))*den1-correc(5)*correc(5)
      vecm(20)=sum((coords(5,1:nraysp)-av(5))*(coords(6,1:nraysp)-av(6)),msk(1:nraysp))*den1-correc(5)*correc(6)
      vecm(21)=sum((coords(6,1:nraysp)-av(6))*(coords(6,1:nraysp)-av(6)),msk(1:nraysp))*den1-correc(6)*correc(6)
!!    write(6,*)sum((coords(5,1:nraysp)-av(5))*(coords(5,1:nraysp)-av(5)),msk(1:nraysp))*den1
!!    write(6,*)sum((coords(6,1:nraysp)-av(6))*(coords(6,1:nraysp)-av(6)),msk(1:nraysp))*den1
      call MPI_ALLREDUCE(vecm,gvecm,21,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) !could replace this by MPI_REDUCE
!     gvecm(1:21)=vecm(1:21)
!!!insert code below to compute eigen-emittances
!     wex=0.d0
!     wey=0.d0
!     wet=0.d0
!!!
!give readable names to entries in the uncorrelated case:
      sq1=gvecm(1)
      sq2=gvecm(7)
      sq3=gvecm(12)
      sq4=gvecm(16)
      sq5=gvecm(19)
      sq6=gvecm(21)
      xpx=gvecm(2)
      ypy=gvecm(13)
      tpt=gvecm(20)
!   compute derived quantities
      epsx2=(sq1*sq2-xpx*xpx)
      epsy2=(sq3*sq4-ypy*ypy)
      epst2=(sq5*sq6-tpt*tpt)
      xrms=sqrt(sq1)
      yrms=sqrt(sq3)
      trms=sqrt(sq5)
      pxrms=sqrt(sq2)
      pyrms=sqrt(sq4)
      ptrms=sqrt(sq6)
      epsx=sqrt(max(epsx2,0.d0))
      epsy=sqrt(max(epsy2,0.d0))
      epst=sqrt(max(epst2,0.d0))
      if(includepi.eq.1)then
        pi=2.d0*asin(1.d0)
        epsx=epsx/pi
        epsy=epsy/pi
        epst=epst/pi
      endif
      if(ncross.eq.1)then
        xpxfac=0.d0
        ypyfac=0.d0
        tptfac=0.d0
        if(xrms.ne.0. .and. pxrms.ne.0.)xpxfac=1./(xrms*pxrms)
        if(yrms.ne.0. .and. pyrms.ne.0.)ypyfac=1./(yrms*pyrms)
        if(trms.ne.0. .and. ptrms.ne.0.)tptfac=1./(trms*ptrms)
        xpx=xpx*xpxfac
        ypy=ypy*ypyfac
        tpt=tpt*tptfac
      endif
!
!
! output results
      if(myrank.eq.0)then
!   set format for amount of data and desired precision
        nlength=nprecision+7
        call num2string2(nlength,a2,2)
        myformat(11:12)=a2
        call num2string2(nprecision,a2,2)
        myformat(14:15)=a2
!!!!!   write(nfile0,myformat)tors,av(1:6),gvecm(1:21),rgood
        write(nfile0,myformat)tors,rv(1:6),gvecm(1:21),rgood !"raw" centroid values
        write(nfile1,myformat)tors,xrms,pxrms,xpx,epsx,av(1),av(2),rgood,theta
        write(nfile2,myformat)tors,yrms,pyrms,ypy,epsy,av(3),av(4),rgood,theta
        write(nfile3,myformat)tors,trms,ptrms,tpt,epst,av(5),av(6),rgood,theta
        call flush(nfile0)
        call flush(nfile1)
        call flush(nfile2)
        call flush(nfile3)
      endif
!here is the reverse rotation to restore coords:
!17 June 2018 commented out since now coords is just a copy of the input array coordsin
!     theta=-theta
!     call xzrot(theta,coords,msk,n1,nraysp)
!     coords(1,1:nraysp)=coords(1,1:nraysp)+x00
!     coords(5,1:nraysp)=coords(5,1:nraysp)+z00
      return
      end

      subroutine xzrot(theta,coords,msk,n1,nraysp)
      implicit none
      real*8 :: theta
      integer :: n1,nraysp
      real*8, dimension(n1,*) :: coords
      logical, dimension(*) :: msk
      real*8 :: co,si,xnew,pxnew,znew,pznew
      integer :: n
      if(n1.lt.6)write(6,*)'xzrot error'
      if(n1.lt.6)call myexit
      co=cos(theta)
      si=sin(theta)
      do n=1,nraysp
        if(.not.msk(n))cycle
        xnew=coords(1,n)*co+coords(5,n)*si
        pxnew=coords(2,n)*co+coords(6,n)*si
        znew=-coords(1,n)*si+coords(5,n)*co
        pznew=-coords(2,n)*si+coords(6,n)*co
        coords(1,n)=xnew
        coords(2,n)=pxnew
        coords(5,n)=znew
        coords(6,n)=pznew
      enddo
      return
      end

      subroutine getcentroid(coords,av,n1,nraysp)
      use mpi
      implicit none
      integer :: n1,nraysp
      real*8, dimension(n1,*) :: coords
      real*8, dimension(6) :: avloc,av
      logical, dimension(nraysp) :: msk
      integer :: n,ngoodloc,ngood,ierr
      real*8 :: den1
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
!hack to use particle 1 on PE0
!     if(myrank.eq.0)write(6,*)'subroutine getcentroid() is hacked to return ptcl#1 on PE0'
!     if(myrank.eq.0)av(1:6)=coords(1:6,1)
!     call MPI_Bcast(av,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
 
!     if(myrank.eq.0)write(6,*)'subroutine getcentroid() is returning the global centroid'
      msk(1:nraysp)=.true. !in the future this should be in the argument list
      ngoodloc=count(msk(1:nraysp))
      call MPI_ALLREDUCE(ngoodloc,ngood,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      den1=1.d0/ngood

!!!   do n=1,6
!!!     avloc(n)=sum(coords(n,1:nraysp),msk(1:nraysp))*den1
!!!   enddo
      avloc(1:6)=0.d0
      do n=1,nraysp
        if(.not.msk(n))cycle
        avloc(1:6)=avloc(1:6)+coords(1:6,n)
      enddo
      avloc(1:6)=avloc(1:6)*den1

      call MPI_ALLREDUCE(avloc,av,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      return
      end

      subroutine profile1d(ptcls,rwall,ncol,ngridpts,nunit) !call with rwall<0 to autoscale
      use mpi
      implicit none
      integer :: ncol,ngridpts,nunit,n1,np,n
      real*8, dimension(:,:) :: ptcls
      real*8 :: rwall,xmin,xmax,xval,hx,hxi,ab
      real*8, dimension(ngridpts) :: rho1d,rhogbl
      real*8, dimension(2) :: bndy,rbndy
      integer :: mprocs,myrank,mpierr,ierr,indx,l
      integer ::  mpistat(MPI_STATUS_SIZE)  ! or should this simply be integer :: mpistat ???
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      n1=ubound(ptcls,1)
      np=ubound(ptcls,2)
      if(rwall.gt.0.d0)then
        xmin=-rwall
        xmax=+rwall
      else
!       determine the range from the particle coords:
        bndy(1)= maxval(ptcls(ncol,1:np))
        bndy(2)=-minval(ptcls(ncol,1:np))
        call MPI_ALLREDUCE(bndy,rbndy,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        xmax=rbndy(1)
        xmin=-rbndy(2)
!       a little bigger to make sure nothing falls off the end:
        if(xmin.lt.0.d0)then
          xmin=xmin*(1.d0+1.d-13)
        else
          xmin=xmin*(1.d0-1.d-13)
        endif
        if(xmax.gt.0.d0)then
          xmax=xmax*(1.d0+1.d-13)
        else
          xmax=xmax*(1.d0-1.d-13)
        endif
      endif

      hx=(xmax-xmin)/(ngridpts-1)
      hxi=1.d0/hx
      rho1d=0.d0
      do n=1,np
        indx=(ptcls(ncol,n)-xmin)*hxi + 1
        if(indx.lt.1 .or. indx.gt.ngridpts-1)then
          write(6,*)'particle #',n,' is out of range; indx=',indx
          write(6,*)'ptcls(' ,ncol, ',' ,n, ')=' ,ptcls(ncol,n)
          cycle
        endif
        ab=((xmin-ptcls(ncol,n))+indx*hx)*hxi
        rho1d(indx)=rho1d(indx)+ab
        rho1d(indx+1)=rho1d(indx+1)+(1.d0-ab)
      enddo
!
! now that every processor has its own 1d profile, sum them together on processor 0:
      call MPI_REDUCE(rho1d,rhogbl,ngridpts,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
! write results:
      if(myrank.eq.0)then
        do n=1,ngridpts
          xval=xmin+(n-1)*hx
          write(nunit,'(3(1pe19.12,1x))')xval,rhogbl(n),rhogbl(n)*hxi
        enddo
      endif
      return
      end

      subroutine num2string2(num,a,ndigits)
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

      subroutine num2string5(num,a,ndigits)
! converts an integer "num" with at most "ndigits" digits into a character string "a"
! if ndigits exceeds the actual number of digits in num, then leading zeroes are inserted in the string
      implicit none
      integer num,ndigits
!     character a*(*)
      character*5 :: a
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

      subroutine openfile(fname,nunit,numsuffix,ndigits)
      implicit none
      character*32 :: fname
      integer :: nunit,numsuffix,ndigits
!local variables:
      character*5 aseq
      integer :: j
!
!!!!!!if(idproc.eq.0)then
        if(ndigits.gt.0)then
          call num2string5(numsuffix,aseq,ndigits)
          j=len_trim(fname)
          fname=fname(1:j)//aseq(1:ndigits)
        endif
        open(unit=nunit,file=fname,status='unknown',err=110)
        goto 150
  110   continue
        write(6,*)'something wrong; error opening file ',fname
!!!!!!endif
  150 continue
      return
      end

end module
