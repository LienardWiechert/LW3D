module lienwiech

use, intrinsic :: iso_fortran_env

implicit none

! Fortran 2008
integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64
integer, parameter, private :: qp = REAL128

contains

      subroutine getselfconsistentfield(t,np_lcl,ptcls,ebself)
      use mpi
      use lwprobcons, only : npart_gbl,chrgperbunch,imax,jmax,kmax
      use diagnostics, only : getcentroid
      implicit none
      real*8 :: t
      integer :: np_lcl
      real*8, dimension(:,:) :: ptcls,ebself
      real*8, dimension(3) :: umin,umax,uspread
      real*8, dimension(6) :: lvec,gvec,av
      real*8, dimension(3) :: obspt
      real*8 :: xmin,xmax,hx,hxi,ymin,ymax,hy,hyi,zmin,zmax,hz,hzi,xval,yval,zval,xcent,ycent,zcent
      real*8 :: obstime
      real*8, dimension(6) :: ebvel,ebrad,ebradtot,ebveltot
      real*8, parameter :: clite=299792458.d0
      real*8, parameter :: fpei=clite**2*1.e-7 !use this if will multiply by charge/bunch
      real*8, parameter :: efpei=clite**2*1.e-7*1.602176487e-19 !use if will multiply by #particles/bunchs
      real*8, dimension(12,imax,jmax,kmax) :: ebcomp,ebcompgbl
      real*8, dimension(6,imax,jmax,kmax) :: ebgrid
      real*8 :: tsearch,coneval,tsearchmin,conevalmax
      integer :: ifailconv,ibadnewt
      integer :: i,j,k,lth,n,mlo,mhi,mlomin,mhimax,ierr
      integer :: indx,jndx,kndx,indxp1,jndxp1,kndxp1
      real*8 :: ab,de,gh
      integer :: ihertz,iticks0,iticks1
      real*8 :: elapsed
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      call system_clock(count_rate=ihertz)
      call system_clock(count=iticks0)
      ebcomp(:,:,:,:)=0.d0 !do I need this?
!get the min and max sizes needed to contain the bunch:
      lvec(1)=-minval(ptcls(1,1:np_lcl))
      lvec(2)=-minval(ptcls(3,1:np_lcl))
      lvec(3)=-minval(ptcls(5,1:np_lcl))
      lvec(4)= maxval(ptcls(1,1:np_lcl))
      lvec(5)= maxval(ptcls(3,1:np_lcl))
      lvec(6)= maxval(ptcls(5,1:np_lcl))
      call MPI_ALLREDUCE(lvec,gvec,6,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      umin(1:3)=-gvec(1:3)
      umax(1:3)= gvec(4:6)
      uspread(1:3)=umax(1:3)-umin(1:3)
      umin(1:3)=umin(1:3)-1.d-9*uspread(1:3)
      umax(1:3)=umax(1:3)+1.d-9*uspread(1:3)
      hx=(umax(1)-umin(1))/dble(imax-1)
      hy=(umax(2)-umin(2))/dble(jmax-1)
      hz=(umax(3)-umin(3))/dble(kmax-1)
      hxi=1.d0/hx
      hyi=1.d0/hy
      hzi=1.d0/hz
!diagnostics (not really needed):
      tsearchmin=9999999.
      conevalmax=-9999999.
      mlomin=999999
      mhimax=-999999
!loop over grid points and compute the LW field at the grid points:
      obstime=t
      do i=1,imax
        xval=umin(1)+(i-1)*hx
        do j=1,jmax
        yval=umin(2)+(j-1)*hy
          do k=1,kmax
            zval=umin(3)+(k-1)*hz
            obspt(1)=xval
            obspt(2)=yval
            obspt(3)=zval
            ebradtot(:)=0.d0
            ebveltot(:)=0.d0
            do lth=1,np_lcl
              ebvel(:)=0.d0 !do I need this?
              ebrad(:)=0.d0 !do I need this?
              call getfieldat4dpt(obstime,obspt,ebvel,ebrad,tsearch,coneval,mlo,mhi,lth,ifailconv,ibadnewt,ierr)
              if(ierr.ne.0)cycle
              if(tsearch.lt.tsearchmin)tsearchmin=tsearch
              if(coneval.gt.conevalmax)conevalmax=coneval
              if(mlo.lt.mlomin)mlomin=mlo
              if(mhi.gt.mhimax)mhimax=mhi
              ebveltot(:)=ebveltot(:)+ebvel
              ebradtot(:)=ebradtot(:)+ebrad
            enddo !over particles
            ebcomp(1:3,i,j,k)=ebveltot(1:3)   !electric velocity field
            ebcomp(4:6,i,j,k)=ebradtot(1:3)   !electric radiation field
            ebcomp(7:9,i,j,k)=ebveltot(4:6)   !magnetic velocity field
            ebcomp(10:12,i,j,k)=ebradtot(4:6) !magnetic radiation field
          enddo
        enddo
      enddo
      call MPI_ALLREDUCE(ebcomp,ebcompgbl,12*imax*jmax*kmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      ebcomp(1:12,1:imax,1:jmax,1:kmax)=ebcompgbl(1:12,1:imax,1:jmax,1:kmax)*fpei*chrgperbunch/npart_gbl
      call system_clock(count=iticks1)
      elapsed=(iticks1-iticks0)/(1.d0*ihertz)
!     if(myrank.eq.0)write(6,*)'PE0 diagnostics:',mlomin,mhimax,conevalmax,elapsed
      if(myrank.eq.0)write(6,"('PE0 diagnostics: ',i5,1x,i5,1x,1pe10.3,1x,1pe12.3)")mlomin,mhimax,conevalmax,elapsed
!
!the following could be edited to keep/exclude various x-y-z components, keep/exclude velocity or radiation terms...
      ebgrid(1,:,:,:)=ebcompgbl(1,:,:,:)+ebcompgbl(4,:,:,:) !Ex velocity + radiation
      ebgrid(2,:,:,:)=ebcompgbl(2,:,:,:)+ebcompgbl(5,:,:,:) !Ey velocity + radiation
      ebgrid(3,:,:,:)=ebcompgbl(3,:,:,:)+ebcompgbl(6,:,:,:) !Ez velocity + radiation
      ebgrid(4,:,:,:)=ebcompgbl(7,:,:,:)+ebcompgbl(10,:,:,:) !Bx velocity + radiation
      ebgrid(5,:,:,:)=ebcompgbl(8,:,:,:)+ebcompgbl(11,:,:,:) !By velocity + radiation
      ebgrid(6,:,:,:)=ebcompgbl(9,:,:,:)+ebcompgbl(12,:,:,:) !Bz velocity + radiation
!
!interpolate the field at the particle locations:
      do n=1,np_lcl
        indx=(ptcls(1,n)-umin(1))*hxi + 1 !arrays are 1-based in this routine
        jndx=(ptcls(3,n)-umin(2))*hyi + 1
        kndx=(ptcls(5,n)-umin(3))*hzi + 1
        indxp1=indx+1
        jndxp1=jndx+1
        kndxp1=kndx+1
        ab=((umin(1)-ptcls(1,n))+indx*hx)*hxi
        de=((umin(2)-ptcls(3,n))+jndx*hy)*hyi
        gh=((umin(3)-ptcls(5,n))+kndx*hz)*hzi
!<this would be the place check for out-of-bounds but I'm skipping it to save time>
        ebself(1:6,n)=                                                  &
     &  ebgrid(1:6,indx,jndx,kndx)*ab*de*gh                             &
     & +ebgrid(1:6,indx,jndxp1,kndx)*ab*(1.-de)*gh                      &
     & +ebgrid(1:6,indx,jndxp1,kndxp1)*ab*(1.-de)*(1.-gh)               &
     & +ebgrid(1:6,indx,jndx,kndxp1)*ab*de*(1.-gh)                      &
     & +ebgrid(1:6,indxp1,jndx,kndxp1)*(1.-ab)*de*(1.-gh)               &
     & +ebgrid(1:6,indxp1,jndxp1,kndxp1)*(1.-ab)*(1.-de)*(1.-gh)        &
     & +ebgrid(1:6,indxp1,jndxp1,kndx)*(1.-ab)*(1.-de)*gh               &
     & +ebgrid(1:6,indxp1,jndx,kndx)*(1.-ab)*de*gh
!!!     ebself(1:6,n)=ebself(1:6,n)*0.d0 ! DEBUGGING debugging
      enddo
      return
      end

      subroutine getfield_atallpts(t,av,np_lcl,nsuffix)
      use mpi
      use lwprobcons, only : npart_gbl,chrgperbunch,xoffmin,xoffmax,yoffmin,yoffmax,zoffmin,zoffmax,&
     &                       imax,jmax,kmax,iwritezeroes,irotategrid,relativetocentroid,&
     &                       rowmajor,blankafter2dblock
      use diagnostics, only : openfile
      implicit none
      real*8 :: t
      real*8, dimension(*) :: av
      integer :: np_lcl,nsuffix
      real*8 :: obstime
      real*8 :: xmin,xmax,hx,ymin,ymax,hy,zmin,zmax,hz,theta,xval,yval,zval
      real*8 :: exvel,eyvel,ezvel,exrad,eyrad,ezrad
      real*8 :: bxvel,byvel,bzvel,bxrad,byrad,bzrad
      real*8 :: poyntxrad,poyntyrad,poyntzrad
      real*8, parameter :: rmu0inv=1.d0/(8.d0*asin(1.d0)*1.d-7)
      real*8 :: tsearch,coneval,tsearchmin,conevalmax,gbltsearchmin,gblconevalmax
      real*8, dimension(3) :: obspt,obsptmin,obsptmax
      real*8, dimension(6) :: ebvel,ebrad,ebradtot,ebveltot
      integer :: i,j,k,lth,mlo,mhi,mlomin,mhimax,gblmlomin,gblmhimax,ierr
      integer :: ifailconv,ibadnewt,ifailconvtot,ibadnewttot,ifailconvgbl,ibadnewtgbl
      character*32 :: fname

      real*8, parameter :: clite=299792458.d0
      real*8, parameter :: fpei=clite**2*1.e-7 !use this if will multiply by charge/bunch
      real*8, parameter :: efpei=clite**2*1.e-7*1.602176487e-19 !use if will multiply by #particles/bunchs

      real*8, dimension(12,imax,jmax,kmax) :: ebcomp,ebcompgbl
      integer :: ihertz,iticks0,iticks1
      real*8 :: elapsed_getfieldatallpts,elapsedmax
      integer :: mprocs,myrank,mpierr
      real*8, allocatable, dimension(:) :: elapsedarray,elapsedarraygbl

      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      ifailconvtot=0
      ibadnewttot=0

      ebcomp(:,:,:,:)=0.d0 !do I need this?

      call system_clock(count_rate=ihertz)
      allocate(elapsedarray(0:mprocs-1))
      allocate(elapsedarraygbl(0:mprocs-1))
      elapsedarray(:)=0.d0
      elapsedarraygbl(:)=0.d0

          obstime=t
          if(myrank.eq.0)write(6,"('centroid=',6(1pe15.8,1x))")av(1:6)

          theta=0.d0  !for undulator radiation simulation
          if(irotategrid.eq.1)theta=atan2(av(2),av(6))  !theta=atan("y","x") where "y"=av(2)=betagamma-x and "x"=av(6)=betagamma-z
          if(myrank.eq.0.and.theta.ne.0.d0)write(6,*)'theta=',theta
          if(myrank.eq.0.and.theta.ne.0.d0)write(6,*)'In what follows, obs pts (x,z) are in a frame aligned w/ direction of motion'

          xmin=relativetocentroid*av(1)+xoffmin
          xmax=relativetocentroid*av(1)+xoffmax
          if(imax.eq.1)then
            hx=0.d0          !sets x to xmin
          else
            hx=(xmax-xmin)/dble(imax-1)
          endif

          ymin=relativetocentroid*av(3)+yoffmin
          ymax=relativetocentroid*av(3)+yoffmax
          if(jmax.eq.1)then
            hy=0.d0          !sets y to ymin
          else
            hy=(ymax-ymin)/dble(jmax-1)
          endif

          zmin=relativetocentroid*av(5)+zoffmin
          zmax=relativetocentroid*av(5)+zoffmax
          if(kmax.eq.1)then
            hz=0.d0          !sets z to zmin
          else
            hz=(zmax-zmin)/dble(kmax-1)
          endif
          if(myrank.eq.0)write(6,"('xmin,xmax,zmin,zmax=',4(1pe15.8,1x))")xmin,xmax,zmin,zmax

          tsearchmin=9999999.
          conevalmax=-9999999.
          mlomin=999999
          mhimax=-999999
          call system_clock(count=iticks0) !time the triple-do-loop
          do i=1,imax
          if(myrank.eq.0)write(6,*)'i=',i
          xval=xmin+(i-1)*hx
          do j=1,jmax
          yval=ymin+(j-1)*hy
          do k=1,kmax
          zval=zmin+(k-1)*hz
            if(irotategrid.eq.0)then
              obspt(1)=xval
              obspt(2)=yval
              obspt(3)=zval
            else
              obspt(1)=av(1) +(xval-av(1))*cos(theta)+(zval-av(5))*sin(theta)
              obspt(2)=av(3)
              obspt(3)=av(5) -(xval-av(1))*sin(theta)+(zval-av(5))*cos(theta)
            endif
            ebradtot(:)=0.d0
            ebveltot(:)=0.d0
            do lth=1,np_lcl
              ebvel(:)=0.d0
              ebrad(:)=0.d0
              call getfieldat4dpt(obstime,obspt,ebvel,ebrad,tsearch,coneval,mlo,mhi,lth,ifailconv,ibadnewt,ierr)
              if(ierr.ne.0)cycle
              ifailconvtot=ifailconvtot+ifailconv
              ibadnewttot=ibadnewttot+ibadnewt
              if(tsearch.lt.tsearchmin.and.ierr.eq.0)tsearchmin=tsearch
              if(coneval.gt.conevalmax.and.ierr.eq.0)conevalmax=coneval
              if(mlo.lt.mlomin.and.ierr.eq.0)then
                mlomin=mlo
                obsptmin(1:3)=obspt(1:3)  !the observation point when the returned value of mlo was least
              endif
              if(mhi.gt.mhimax.and.ierr.eq.0)then
                mhimax=mhi
                obsptmax(1:3)=obspt(1:3)  !the observation point when the returned value of mhi was most
              endif
              ebveltot(:)=ebveltot(:)+ebvel
              ebradtot(:)=ebradtot(:)+ebrad
            enddo !over particles
!ebxxxtot(1:3) are electric field; 4:6 are magnetic field
            ebcomp(1:3,i,j,k)=ebveltot(1:3)   !electric velocity field
            ebcomp(4:6,i,j,k)=ebradtot(1:3)   !electric radiation field
            ebcomp(7:9,i,j,k)=ebveltot(4:6)   !magnetic velocity field
            ebcomp(10:12,i,j,k)=ebradtot(4:6) !magnetic radiation field
          enddo !k
          enddo !j
          enddo !i
          call system_clock(count=iticks1)
          elapsedarray(myrank)=(iticks1-iticks0)/(1.d0*ihertz)
          if(myrank.eq.0)write(6,*)'PE0 tsearchmin=',tsearchmin
          if(myrank.eq.0)write(6,*)'PE0 conevalmax=',conevalmax
          if(myrank.eq.0)write(6,*)'PE0 mlomin, mhimax =',mlomin,mhimax
          if(myrank.eq.0)write(6,*)'PE0 obspt at mlomin =',obsptmin
          if(myrank.eq.0)write(6,*)'PE0 obspt at mhimax =',obsptmax
          call MPI_REDUCE(tsearchmin,gbltsearchmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
          call MPI_REDUCE(conevalmax,gblconevalmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
          call MPI_REDUCE(mlomin,gblmlomin,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_WORLD,ierr)
          call MPI_REDUCE(mhimax,gblmhimax,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
          call MPI_REDUCE(elapsedarray,elapsedarraygbl,mprocs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if(myrank.eq.0)write(6,*)'global tsearchmin=',gbltsearchmin
          if(myrank.eq.0)write(6,*)'global conevalmax=',gblconevalmax
          if(myrank.eq.0)write(6,*)'global mlomin,mhimax=',gblmlomin,gblmhimax
!         if(myrank.eq.0)then
!           do i=0,mprocs-1
!           write(999,*)elapsedarraygbl(i)
!           enddo
!           close(999)
!         endif

          call MPI_REDUCE(ebcomp,ebcompgbl,12*imax*jmax*kmax,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if(myrank.eq.0)ebcomp(1:12,1:imax,1:jmax,1:kmax)=ebcompgbl(1:12,1:imax,1:jmax,1:kmax)*fpei*chrgperbunch/npart_gbl
!
!row-major ordering of the output:
          if(rowmajor)then
          if(myrank.eq.0)then
            fname='field'
            call openfile(fname,300,nsuffix,5)
            do i=1,imax
              xval=xmin+(i-1)*hx
              do j=1,jmax
                yval=ymin+(j-1)*hy
                do k=1,kmax
                zval=zmin+(k-1)*hz
                if(irotategrid.eq.0)then
                  obspt(1)=xval
                  obspt(2)=yval
                  obspt(3)=zval
                else
                  obspt(1)=av(1) +(xval-av(1))*cos(theta)+(zval-av(5))*sin(theta)
                  obspt(2)=av(3)
                  obspt(3)=av(5) -(xval-av(1))*sin(theta)+(zval-av(5))*cos(theta)
                endif
                exvel=ebcomp(1,i,j,k)*cos(theta)-ebcomp(3,i,j,k)*sin(theta)
                eyvel=ebcomp(2,i,j,k)
                ezvel=ebcomp(1,i,j,k)*sin(theta)+ebcomp(3,i,j,k)*cos(theta)
                exrad=ebcomp(4,i,j,k)*cos(theta)-ebcomp(6,i,j,k)*sin(theta)
                eyrad=ebcomp(5,i,j,k)
                ezrad=ebcomp(4,i,j,k)*sin(theta)+ebcomp(6,i,j,k)*cos(theta)
                bxvel=ebcomp(7,i,j,k)*cos(theta)-ebcomp(9,i,j,k)*sin(theta)
                byvel=ebcomp(8,i,j,k)
                bzvel=ebcomp(7,i,j,k)*sin(theta)+ebcomp(9,i,j,k)*cos(theta)
                bxrad=ebcomp(10,i,j,k)*cos(theta)-ebcomp(12,i,j,k)*sin(theta)
                byrad=ebcomp(11,i,j,k)
                bzrad=ebcomp(10,i,j,k)*sin(theta)+ebcomp(12,i,j,k)*cos(theta)
                poyntxrad=(eyrad*bzrad-ezrad*byrad)*rmu0inv
                poyntyrad=(ezrad*bxrad-exrad*bzrad)*rmu0inv
                poyntzrad=(exrad*byrad-eyrad*bxrad)*rmu0inv
                if(ezrad.eq.0.d0.and.iwritezeroes.eq.0)cycle !should check other components not just ezrad
                write(300,'(3(1pe20.13,1x),15(1pe15.8,1x),i7,1x,i7,1x,2(1pe20.13,1x))') &
     &                obspt(1),obspt(2),obspt(3),          &
     &                exvel,eyvel,ezvel,exrad,eyrad,ezrad, &
     &                bxvel,byvel,bzvel,bxrad,byrad,bzrad, &
     &                poyntxrad,poyntyrad,poyntzrad,       &
     &                i,k,xval-av(1),zval-av(5)
                enddo !k
              enddo !j
              if(blankafter2dblock)write(300,*)' '
            enddo !i
            close(300)
          endif !myrank.eq.0
          endif !rowmajor
!
!column-major ordering of the output:
          if(.not.rowmajor)then
          if(myrank.eq.0)then
            fname='field'
            call openfile(fname,300,nsuffix,5)
            do k=1,kmax
              zval=zmin+(k-1)*hz
              do j=1,jmax
                yval=ymin+(j-1)*hy
                do i=1,imax
                xval=xmin+(i-1)*hx
                if(irotategrid.eq.0)then
                  obspt(1)=xval
                  obspt(2)=yval
                  obspt(3)=zval
                else
                  obspt(1)=av(1) +(xval-av(1))*cos(theta)+(zval-av(5))*sin(theta)
                  obspt(2)=av(3)
                  obspt(3)=av(5) -(xval-av(1))*sin(theta)+(zval-av(5))*cos(theta)
                endif
                exvel=ebcomp(1,i,j,k)*cos(theta)-ebcomp(3,i,j,k)*sin(theta)
                eyvel=ebcomp(2,i,j,k)
                ezvel=ebcomp(1,i,j,k)*sin(theta)+ebcomp(3,i,j,k)*cos(theta)
                exrad=ebcomp(4,i,j,k)*cos(theta)-ebcomp(6,i,j,k)*sin(theta)
                eyrad=ebcomp(5,i,j,k)
                ezrad=ebcomp(4,i,j,k)*sin(theta)+ebcomp(6,i,j,k)*cos(theta)
                bxvel=ebcomp(7,i,j,k)*cos(theta)-ebcomp(9,i,j,k)*sin(theta)
                byvel=ebcomp(8,i,j,k)
                bzvel=ebcomp(7,i,j,k)*sin(theta)+ebcomp(9,i,j,k)*cos(theta)
                bxrad=ebcomp(10,i,j,k)*cos(theta)-ebcomp(12,i,j,k)*sin(theta)
                byrad=ebcomp(11,i,j,k)
                bzrad=ebcomp(10,i,j,k)*sin(theta)+ebcomp(12,i,j,k)*cos(theta)
                poyntxrad=(eyrad*bzrad-ezrad*byrad)*rmu0inv
                poyntyrad=(ezrad*bxrad-exrad*bzrad)*rmu0inv
                poyntzrad=(exrad*byrad-eyrad*bxrad)*rmu0inv
                if(ezrad.eq.0.d0.and.iwritezeroes.eq.0)cycle !should check other components not just ezrad
                write(300,'(3(1pe20.13,1x),15(1pe15.8,1x),i7,1x,i7,1x,2(1pe20.13,1x))') &
     &                obspt(1),obspt(2),obspt(3),          &
     &                exvel,eyvel,ezvel,exrad,eyrad,ezrad, &
     &                bxvel,byvel,bzvel,bxrad,byrad,bzrad, &
     &                poyntxrad,poyntyrad,poyntzrad,       &
     &                i,k,xval-av(1),zval-av(5)
                enddo !i
              enddo !j
              if(blankafter2dblock)write(300,*)' '
            enddo !k
            close(300)
          endif !myrank.eq.0
          endif !.not.rowmajor
      return
      end

      subroutine getfieldat4dpt(tobs,obspt,ebvel,ebrad,tsrch,fnormx,mlortn,mhirtn,ipart,ifailconv,ibadnewt,ierr)
      use mpi
      use lwprobcons, only : refptclhist,refptcltime,laststoredrefpart,nsteps
      implicit none
      real*8 :: tobs,tsrch,fnormx
      real*8, dimension(3) :: obspt
      real*8, dimension(6) :: ebvel,ebrad
      integer :: mlortn,mhirtn,ipart,ifailconv,ibadnewt,ierr
      integer :: k,mlo,mhi

      real*8 :: xobs,yobs,zobs,valueleftend,valuerightend,valuelo,valuehi,valueat
      real*8, parameter :: clite=299792458.d0
      real*8, parameter :: clite2=299792458.d0**2

      real*8, dimension(6) :: zetalcl,zetagbl,zzzeta1,zzzeta2,zzzeta
      real*8, allocatable, dimension(:,:) :: zetasrch,f,zetanewt,rhsnewt,zetatmp
      real*8 :: tnewt,tnewtreturned,errbrent
      real*8 :: exvel,eyvel,ezvel,bxvel,byvel,bzvel,exrad,eyrad,ezrad,bxrad,byrad,bzrad,aphi1,aphi2,aphi3,aphi4
      real*8 :: tzerox,fnormlo,fnormhi

      real*8 :: gam2,gam,alfx,alfy,alfz,betx,bety,betz,aa,bb,cc,tplus,tminu
      real*8 :: xplus,yplus,zplus,xminu,yminu,zminu,valueplus,valueminu
      real*8, dimension(3) :: betadotval
      integer :: icountbrent,ididzb

      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)

      ifailconv=0  !note if newton fails to converge
      ibadnewt=0   !note if newton produces a garbage result (outside of search bounds)

      tsrch=999999.
      fnormx=-999999.

      if(ipart.eq.-1)then
        return
      endif

      xobs=obspt(1)
      yobs=obspt(2)
      zobs=obspt(3)
      valueleftend=(xobs-refptclhist(1,0,ipart))**2 + (yobs-refptclhist(3,0,ipart))**2+  &
     &         (zobs-refptclhist(5,0,ipart))**2-clite2*(tobs-refptcltime(0))**2
      k=laststoredrefpart
      valuerightend=(xobs-refptclhist(1,k,ipart))**2 + (yobs-refptclhist(3,k,ipart))**2+  &
     &         (zobs-refptclhist(5,k,ipart))**2-clite2*(tobs-refptcltime(k))**2

      if(valueleftend.gt.0.d0)then
        ebvel(1:6)=0.d0
        ebrad(1:6)=0.d0
        ierr=-1
        return
      endif
      if(valuerightend.lt.0.d0)then
        ebvel(1:6)=0.d0
        ebrad(1:6)=0.d0
        ierr=1
        return
      endif
      ierr=0
      
      call bisectionsearch(refptcltime,refptclhist(1,0,ipart),xobs,yobs,zobs,tobs,0,nsteps,k,mlo,mhi)
      mlortn=mlo
      mhirtn=mhi
      if(mhirtn.ne.mlortn+1)then
        ebvel(1:6)=0.d0
        ebrad(1:6)=0.d0
        ierr=9
        return
      endif

!this could be commented out to save time:
      valuelo=(xobs-refptclhist(1,mlo,ipart))**2+(yobs-refptclhist(3,mlo,ipart))**2+(zobs-refptclhist(5,mlo,ipart))**2-&
     &         clite2*(tobs-refptcltime(mlo))**2
      valuehi=(xobs-refptclhist(1,mhi,ipart))**2+(yobs-refptclhist(3,mhi,ipart))**2+(zobs-refptclhist(5,mhi,ipart))**2-&
     &         clite2*(tobs-refptcltime(mhi))**2
      if(valuelo*valuehi.gt.0.d0)then
        write(6,*)'problem with valuelo*valuehi'
        write(6,*)'mlo,mhi=',mlo,mhi
        write(6,*)'valuelo,valuehi=',valuelo,valuehi
        stop
      endif

!Interpolate to estimate the zero-crossing:
      tzerox=refptcltime(mlo)-valuelo*(refptcltime(mhi)-refptcltime(mlo))/(valuehi-valuelo)

      allocate(zetanewt(6,1),rhsnewt(6,1),zetatmp(6,1),zetasrch(6,1),f(6,1))
!Check for before start of history:
      if(mlo.eq.0 .and. mhi.eq.1)then
              gam2=1.d0+refptclhist(2,mlo,ipart)**2+refptclhist(4,mlo,ipart)**2+refptclhist(6,mlo,ipart)**2 !same for mlo,mhi
              gam=sqrt(gam2)
              alfx=clite/gam*refptclhist(2,mlo,ipart)
              alfy=clite/gam*refptclhist(4,mlo,ipart)
              alfz=clite/gam*refptclhist(6,mlo,ipart)
              betx=refptclhist(1,mlo,ipart)-alfx*refptcltime(mlo)-xobs
              bety=refptclhist(3,mlo,ipart)-alfy*refptcltime(mlo)-yobs
              betz=refptclhist(5,mlo,ipart)-alfz*refptcltime(mlo)-zobs
              aa=alfx**2+alfy**2+alfz**2-clite**2
              bb=2.d0*(alfx*betx+alfy*bety+alfz*betz+tobs*clite**2 )
              cc=betx**2+bety**2+betz**2-(tobs*clite)**2
              if(bb.ge.0.d0)then
                tplus=(-bb-sqrt(bb**2-4.d0*aa*cc))/(2.d0*aa)
                tminu=2.d0*cc/(-bb-sqrt(bb**2-4.d0*aa*cc))
              else
                tplus=2.d0*cc/(-bb+sqrt(bb**2-4.d0*aa*cc))
                tminu=(-bb+sqrt(bb**2-4.d0*aa*cc))/(2.d0*aa)
              endif
              if(tplus.ge.refptcltime(mlo).and.tplus.le.refptcltime(mhi))then
                xplus=refptclhist(1,mlo,ipart)+alfx*(tplus-refptcltime(mlo))
                yplus=refptclhist(3,mlo,ipart)+alfy*(tplus-refptcltime(mlo))
                zplus=refptclhist(5,mlo,ipart)+alfz*(tplus-refptcltime(mlo))
                tnewt=tplus
                zetanewt(1,1)=xplus
                zetanewt(2,1)=refptclhist(2,mlo,ipart)
                zetanewt(3,1)=yplus
                zetanewt(4,1)=refptclhist(4,mlo,ipart)
                zetanewt(5,1)=zplus
                zetanewt(6,1)=refptclhist(6,mlo,ipart)
!               valueplus=(xobs-xplus)**2 + (yobs-yplus)**2+ (zobs-zplus)**2-clite2*(tobs-tplus)**2
                goto 1357
              else
                xminu=refptclhist(1,mlo,ipart)+alfx*(tminu-refptcltime(mlo))
                yminu=refptclhist(3,mlo,ipart)+alfy*(tminu-refptcltime(mlo))
                zminu=refptclhist(5,mlo,ipart)+alfz*(tminu-refptcltime(mlo))
                tnewt=tminu
                zetanewt(1,1)=xminu
                zetanewt(2,1)=refptclhist(2,mlo,ipart)
                zetanewt(3,1)=yminu
                zetanewt(4,1)=refptclhist(4,mlo,ipart)
                zetanewt(5,1)=zminu
                zetanewt(6,1)=refptclhist(6,mlo,ipart)
                goto 1357
              endif
      endif

!this version uses a brent search instead of a newton search
      zzzeta1(1:6)=refptclhist(1:6,mlo,ipart)
      zzzeta2(1:6)=refptclhist(1:6,mhi,ipart)
      call zbrent(xobs,yobs,zobs,tobs,zzzeta1,refptcltime(mlo),zzzeta2,refptcltime(mhi),1.d-16,0.d0,0.d0,tnewt,&
     &            errbrent,icountbrent,ididzb)
      if(ididzb.ne.1)write(6,*)'ididzb=',ididzb
        zetanewt(1:6,1)=zzzeta1(1:6)
        call rk4silent(zetanewt,refptcltime(mlo),tnewt,10,evalt1)   !changed 1 to 10

 1357 continue
      tsrch=tnewt
      zetalcl(:)=zetanewt(1:6,1)
      zetasrch(1:6,1)=zetanewt(1:6,1)

      call fourvecnormsq(xobs,yobs,zobs,tobs,zetasrch,tsrch,fnormx)

      zetagbl(:)=zetalcl(:)
      call lwfieldcalc(xobs,yobs,zobs,tsrch,zetalcl,zetagbl,betadotval, &
     &    exvel,eyvel,ezvel,bxvel,byvel,bzvel, &
     &    exrad,eyrad,ezrad,bxrad,byrad,bzrad, &
     &    aphi1,aphi2,aphi3,aphi4,ipart)
      ebvel(1)=exvel; ebvel(2)=eyvel; ebvel(3)=ezvel
      ebvel(4)=bxvel; ebvel(5)=byvel; ebvel(6)=bzvel
      ebrad(1)=exrad; ebrad(2)=eyrad; ebrad(3)=ezrad
      ebrad(4)=bxrad; ebrad(5)=byrad; ebrad(6)=bzrad

      return
      end

      subroutine bisectionsearch(tvec,coords,x,y,z,tval,n1,n2,nlast,mlo,mhi)
      implicit none
      integer :: n1,n2,nlast
      real*8 :: x,y,z,tval
      real*8, dimension(n1:n2) :: tvec
      real*8, dimension(6,n1:n2) :: coords !in the future, consider first index other than 6
      real*8, parameter :: clite2=299792458.d0**2
      real*8 :: trhs0,tvecm,x1,y1,z1
      integer :: mlo,mhi,m
      mlo=n1
      mhi=nlast
    9 if (mhi-mlo.gt.1)then
        m=(mhi+mlo)/2
        x1=coords(1,m)
        y1=coords(3,m)
        z1=coords(5,m)
        tvecm=tvec(m)
        trhs0=(x-x1)**2+(y-y1)**2+(z-z1)**2-clite2*(tval-tvecm)**2
        if(trhs0.gt.0.d0)then
           mhi=m
        else
           mlo=m
        endif
        goto 9
      endif
      return
      end subroutine bisectionsearch
! 
      subroutine fourvecnormsq(x,y,z,tval,zeta1,t1,value)
      implicit none
      real*8 :: x,y,z,tval,t1,value
      real*8, dimension(6) :: zeta1
      real*8 :: clite2
      clite2=299792458.d0**2
      value=(x-zeta1(1))**2+(y-zeta1(3))**2+(z-zeta1(5))**2-clite2*(tval-t1)**2
      return
      end subroutine fourvecnormsq
!
      subroutine lwfieldcalc(x,y,z,t,zetalcl,zetagbl,betadotval, & !hardcoded to use/notuse betadotval
     &    exvel,eyvel,ezvel,bxvel,byvel,bzvel, &
     &    exrad,eyrad,ezrad,bxrad,byrad,bzrad, &
     &    aphi1,aphi2,aphi3,aphi4,ipart)
      use mpi
      use externalfield, only : getexternalfield
      implicit none
      integer :: ipart
      real*8, dimension(6) :: zetalcl,zetagbl
      real*8, dimension(3) :: betadotval
      real*8 :: x,y,z,t,exvel,eyvel,ezvel,bxvel,byvel,bzvel,exrad,eyrad,ezrad,bxrad,byrad,bzrad,aphi1,aphi2,aphi3,aphi4
      real*8 :: xoft,xprimeoft,yoft,yprimeoft,zoft,zprimeoft
      real*8 :: udiff1,udiff2,udiff3,rinv2,rinv,unorm,un1,un2,un3,b1,b2,b3,ub1,ub2,ub3,xi,xi2inv,wh1,wh2,wh3,uwdot
      real*8 :: gamma0,gam2inv,gaminv,betacrossb1,betacrossb2,betacrossb3,gammadot,bdot1,bdot2,bdot3,bnfac
      real*8 :: denominv3,triprodx,triprody,triprodz
      real*8, parameter :: clite=299792458.d0
      real*8, parameter :: cliteinv=1.d0/299792458.d0
      real*8, dimension(3) :: e0,b0
      real*8 :: xmc2,qbymc,qbym
      integer, parameter :: iusebetadotval=0 !this is an old option; is no longer used.
      integer :: ierr
      integer :: mprocs,myrank,mpierr
      call MPI_COMM_SIZE(MPI_COMM_WORLD,mprocs,mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      xmc2=0.510998910d6
      qbymc=1.d0/xmc2*clite
      qbym=1.d0/xmc2*clite**2
      gamma0=sqrt(1.d0+zetalcl(2)**2+zetalcl(4)**2+zetalcl(6)**2)
      xoft=zetalcl(1)
      xprimeoft=zetalcl(2)/gamma0*clite
      yoft=zetalcl(3)
      yprimeoft=zetalcl(4)/gamma0*clite
      zoft=zetalcl(5)
      zprimeoft=zetalcl(6)/gamma0*clite
        udiff1=x-xoft
        udiff2=y-yoft
        udiff3=z-zoft
        rinv2=1.d0/(udiff1**2+udiff2**2+udiff3**2)
        rinv=sqrt(rinv2)
        unorm=1.d0/rinv
        un1=udiff1/unorm
        un2=udiff2/unorm
        un3=udiff3/unorm
        b1=xprimeoft/clite
        b2=yprimeoft/clite
        b3=zprimeoft/clite
        ub1=un1-b1
        ub2=un2-b2
        ub3=un3-b3
        xi=sqrt(ub1**2+ub2**2+ub3**2)
        xi2inv=1.d0/(ub1**2+ub2**2+ub3**2)
        wh1=ub1/xi
        wh2=ub2/xi
        wh3=ub3/xi
        uwdot=(un1*wh1+un2*wh2+un3*wh3)
        denominv3=1.d0/uwdot**3

        gam2inv=1-b1**2-b2**2-b3**2
       if(iusebetadotval.eq.0)then
        call getexternalfield(zetagbl(1),zetagbl(3),zetagbl(5),t,e0,b0,ierr)
        gaminv=sqrt(gam2inv)
        betacrossb1=b2*b0(3)-b3*b0(2)
        betacrossb2=b3*b0(1)-b1*b0(3)
        betacrossb3=b1*b0(2)-b2*b0(1)
        gammadot=qbymc*(b1*e0(1)+b2*e0(2)+b3*e0(3))
        bdot1=gaminv*(qbymc*e0(1)+qbym*betacrossb1-b1*gammadot)
        bdot2=gaminv*(qbymc*e0(2)+qbym*betacrossb2-b2*gammadot)
        bdot3=gaminv*(qbymc*e0(3)+qbym*betacrossb3-b3*gammadot)
       else
        bdot1=betadotval(1)
        bdot2=betadotval(2)
        bdot3=betadotval(3)
       endif
!   
        bnfac=1.d0-(b1*un1+b2*un2+b3*un3)
!       aphi1=rinv/bnfac
!       aphi2=b1*rinv/bnfac
!       aphi3=b2*rinv/bnfac
!       aphi4=b3*rinv/bnfac
        exvel=gam2inv*rinv2*xi2inv*wh1*denominv3
        eyvel=gam2inv*rinv2*xi2inv*wh2*denominv3
        ezvel=gam2inv*rinv2*xi2inv*wh3*denominv3
        bxvel=(un2*ezvel-un3*eyvel)*cliteinv        !uncommented 4/20   and included cliteinv 4/21
        byvel=(un3*exvel-un1*ezvel)*cliteinv        !uncommented 4/20   and included cliteinv 4/21
        bzvel=(un1*eyvel-un2*exvel)*cliteinv        !uncommented 4/20   and included cliteinv 4/21
        triprodx=un2*(wh1*bdot2-wh2*bdot1)-un3*(wh3*bdot1-wh1*bdot3)
        triprody=un3*(wh2*bdot3-wh3*bdot2)-un1*(wh1*bdot2-wh2*bdot1)
        triprodz=un1*(wh3*bdot1-wh1*bdot3)-un2*(wh2*bdot3-wh3*bdot2)
        exrad=rinv*xi2inv*denominv3*triprodx*cliteinv
        eyrad=rinv*xi2inv*denominv3*triprody*cliteinv
        ezrad=rinv*xi2inv*denominv3*triprodz*cliteinv
        bxrad=(un2*ezrad-un3*eyrad)*cliteinv        !uncommented 4/20   and included cliteinv 4/21
        byrad=(un3*exrad-un1*ezrad)*cliteinv        !uncommented 4/20   and included cliteinv 4/21
        bzrad=(un1*eyrad-un2*exrad)*cliteinv        !uncommented 4/20   and included cliteinv 4/21
      return
      end

      subroutine rk4silent(y,tini,tfin,nsteps,eval)
      use mpi
      implicit none
      real*8, dimension(:,:) :: y
      real*8 :: tini,tfin
      integer :: nsteps,ipart
      real*8, allocatable, dimension(:,:) :: yt,f,a,b,c,d
      integer :: n1,np,n,irtn,k
      real*8 :: h,t,tt
      integer :: mprocs,myrank,mpierr

      real*8 :: obstime,zmin,zmax,hz
      real*8, dimension(3) :: obspt
      real*8, dimension(6) :: ebvel,ebrad

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

      subroutine zbrent(x,y,z,tval,zeta1,t1,zeta2,t2,tol,rtol,atol,tret,errbrent,icountbrent,ididzb)
! for a given observation point in space (x,y,z) and and observation time tval,
! find retarded quantities zetaret=(x,vx,y,vy,z,vz)_ret and the retarded time tret for which the 4-vector norm is zero
! the 4-vector norm is (x-zetaret(1))**2+(y-zetaret(3))**2+(z-zetaret(5))**2-clite**2*(tval-tret))**2
! on input, the 4-vector norm is known to change sign between t1 and t2
! adapted from numerical recipes function zbrent
      implicit none
      real*8 :: x,y,z,tval,t1,t2,tret,errbrent
      integer :: icountbrent,ididzb,ididdop
      real*8, dimension(6) :: zeta1,zeta2
      real*8, dimension(6) :: zetaret
      real*8, allocatable, dimension(:,:) :: zetaret2d
      integer, parameter :: itmax=20
      real*8, parameter :: eps=3.d-8
      real*8 :: tol,rtol,atol !tol=zbrent tolerance; rtol,atol=dop853 tolerances
      integer, parameter :: iout=0
      integer :: iter
      real*8 :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      real*8, dimension(1) :: rpar
      integer, dimension(1) :: ipar
      allocate(zetaret2d(6,1))
      a=t1
      b=t2
!!!!! fa=func(a)
!!!!! fb=func(b)
      call fourvecnormsq(x,y,z,tval,zeta1,a,fa)
      call fourvecnormsq(x,y,z,tval,zeta2,b,fb)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))then
!!!!!   write(6,*)'error: root not bracketed'
!!!!!   stop
        ididzb=3
        return
      endif
      c=b
      fc=fb
      do iter=1,itmax
        icountbrent=iter
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*eps*abs(b)+0.5*tol !Convergence check.
        xm=0.5d0*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          tret=b
          errbrent=xm
          ididzb=1 !success
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa !Attempt inverse quadratic interpolation.
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q !Check whether in bounds.
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d   !accept interpolation
            d=p/q
          else
            d=xm !Interpolation failed, use bisection.
            e=d
          endif
        else !Bounds decreasing too slowly, use bisection.
          d=xm
          e=d
        endif
        a=b !Move last best guess to a.
        fa=fb
        if(abs(d) .gt. tol1) then !Evaluate new trial root.
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        zetaret2d(1:6,1)=zeta1(1:6)
        call rk4silent(zetaret2d,t1,b,1,evalt1)
!       call rk4silent(zetaret2d,t1,b,10,evalt1)
        zetaret(1:6)=zetaret2d(1:6,1)
        call fourvecnormsq(x,y,z,tval,zetaret,b,fb)
      enddo
      write(6,*)'zbrent failure: exceeded maximum iterations'
      ididzb=4
      tret=b
      return
      end

      subroutine evalt1(t,y,f) ! 1particle, t is the independent variable
      use mpi
      use externalfield, only : getexternalfield
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
      clite=299792458.d0
      cliteinv=1.d0/clite
      xmc2=0.510998910d6
      qbymcc=1.d0/xmc2
      qbymc=1.d0/xmc2*clite
      qbym=1.d0/xmc2*clite**2
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
      return
      end

end module
