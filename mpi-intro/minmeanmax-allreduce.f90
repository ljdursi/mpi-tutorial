       program randomdata
       use mpi
       implicit none

       integer,parameter :: nx=1500
       real, allocatable :: dat(:)

       integer :: i,n
       real    :: datamin, datamax, datamean
       real    :: globmin, globmax, globmean
       integer :: rank, comsize, ierr

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,comsize,ierr)
       
!
! random data
!
       allocate(dat(nx))
       call random_seed(size=n)
       call random_seed(put=[(rank,i=1,n)])
       call random_number(dat)
       dat = 2*dat - 1.

!
! find local min/mean/max
! 
       datamin = minval(dat)
       datamax = maxval(dat)
       datamean= (1.*sum(dat))/nx
       deallocate(dat)

       print *,rank,': min/mean/max = ', datamin, datamean, datamax
!
! combine data
!
       call MPI_ALLREDUCE(datamin, globmin, 1, MPI_REAL, MPI_MIN, &
                          MPI_COMM_WORLD, ierr)
!
! to just send to task 0:
!       call MPI_REDUCE(datamin, globmin, 1, MPI_REAL, MPI_MIN,
!     &                    0, MPI_COMM_WORLD, ierr)
!
       call MPI_ALLREDUCE(datamax, globmax, 1, MPI_REAL, MPI_MAX, &
                         MPI_COMM_WORLD, ierr)
       call MPI_ALLREDUCE(datamean, globmean, 1, MPI_REAL, MPI_SUM, &
                         MPI_COMM_WORLD, ierr)
       globmean = globmean/comsize
       if (rank == 0) then
           print *, rank,': Global min/mean/max=',globmin,globmean,globmax 
       endif

       call MPI_FINALIZE(ierr)
       end program randomdata
