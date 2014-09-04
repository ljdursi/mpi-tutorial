       program randomdata
       use mpi
       implicit none
       integer,parameter :: nx=1500
       real, allocatable :: dat(:)

       integer :: i,n
       real    :: datamin, datamax, datamean
       real    :: globmin, globmax, globmean
       integer :: ierr, rank, comsize
       integer :: ourtag=5
       real, dimension(3) :: sendbuffer, recvbuffer
       integer, dimension(MPI_STATUS_SIZE) :: status

       call MPI_INIT(ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, comsize, ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
!
! random data
!
       allocate(dat(nx))
       call random_seed(size=n)
       call random_seed(put=[(rank,i=1,n)])
       call random_number(dat)
       dat = 2*dat - 1.

!
! find min/mean/max
! 
       datamin = minval(dat)
       datamax = maxval(dat)
       datamean= (1.*sum(dat))/nx
       deallocate(dat)

       if (rank /= 0) then
            sendbuffer(1) = datamin
            sendbuffer(2) = datamean
            sendbuffer(3) = datamax
            call MPI_SSEND()  ! stuff goes here
       else
            globmin = datamin
            globmax = datamax
            globmean = datamean
            do i=2,comsize 
                call MPI_RECV() ! stuff goes here
                if (recvbuffer(1) < globmin) globmin=recvbuffer(1)
                if (recvbuffer(3) > globmax) globmax=recvbuffer(3)
                globmean = globmean + recvbuffer(2)
            enddo
            globmean = globmean / comsize
       endif

       print *,rank, ': min/mean/max = ', datamin, datamean, datamax

       if (rank==0) then
           print *, 'global min/mean/max = ', globmin, globmean, globmax
       endif 
  
       call MPI_FINALIZE(ierr)
       end
