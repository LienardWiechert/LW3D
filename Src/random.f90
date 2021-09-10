      subroutine initrandom(inpseed)
      implicit none
      integer :: inpseed
      integer :: iseedsize,n,i
      integer, dimension(:), allocatable :: iputarray
! determine the seed size and allocate the seed array:
      call random_seed(size=iseedsize)
      allocate(iputarray(iseedsize))
! set a unique seed for each proc:
      n=inpseed
      do i=1,iseedsize
        n=mod(9301*n+49297,233280)
        iputarray(i)=n
      enddo
      call random_seed(put=iputarray)
      return
      end subroutine initrandom
