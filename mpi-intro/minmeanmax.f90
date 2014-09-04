       program randomdata
       implicit none
       integer,parameter :: nx=1500
       real, allocatable :: dat(:)

       integer :: i,n
       real    :: datamin, datamax, datamean

!
! random data
!
       allocate(dat(nx))
       call random_seed(size=n)
       call random_seed(put=[(i,i=1,n)])
       call random_number(dat)
       dat = 2*dat - 1.

!
! find min/mean/max
! 
       datamin = minval(dat)
       datamax = maxval(dat)
       datamean= (1.*sum(dat))/nx

       deallocate(dat)

       print *,'min/mean/max = ', datamin, datamean, datamax
  
       return
       end


  
