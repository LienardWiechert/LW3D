      subroutine myexit
      use mpi
      implicit none
      integer :: myrank,mpierr
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,mpierr)
      if(myrank.eq.0)write(6,*)'job finished'
      call MPI_FINALIZE(mpierr)
!     stop
      call exit()
      end
