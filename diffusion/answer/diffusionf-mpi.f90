       program diffuse
       use mpi
       implicit none
   
!
! simulation parameters
!
       integer, parameter :: totpoints=1000
       real, parameter    :: xleft=-12., xright=+12.
       real, parameter    :: kappa=1.
       integer, parameter :: nsteps=100000
       integer, parameter :: plotsteps=50

!
! the calculated temperature, and the known correct
! solution from theory
!
       real, allocatable :: x(:)
       real, allocatable, target :: temperature(:,:)
       real, allocatable :: theory(:)
       real, dimension(:), pointer :: old, new, tmp

       integer :: step
       integer :: i
       integer :: red, grey, white
       real :: time
       real :: dt, dx
       real :: error

!
!  parameters of the original temperature distribution
!
       real, parameter :: ao=1., sigmao = 1.
       real a, sigma

!
!  mpi variables
! 
       integer :: ierr, rank, comsize
       integer :: locpoints, startn, endn
       real    :: locxleft, locxright
       integer :: left, right
       integer :: lefttag=1, righttag=2
       integer, dimension(MPI_STATUS_SIZE) :: rstatus

       call MPI_Init(ierr)
       call MPI_Comm_size(MPI_COMM_WORLD,comsize,ierr)
       call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

       locpoints = totpoints/comsize
       startn = rank*locpoints+1
       endn   = startn + locpoints
       if (rank == comsize-1) endn=totpoints
       locpoints = endn-startn+1
   
       left = rank-1
       if (left < 0) left = MPI_PROC_NULL
       right = rank+1
       if (right >= comsize) right = MPI_PROC_NULL
!
! set parameters
!
       dx = (xright-xleft)/(totpoints-1)
       dt = dx**2 * kappa/10.

       locxleft = xleft + dx*(startn-1)
       locxright= xleft + dx*endn

! 
! allocate data, including ghost cells: old and new timestep
! theory doesn't need ghost cells, but we include it for simplicity
!
       allocate(temperature(locpoints+2,2))
       allocate(theory(locpoints+2))
       allocate(x(locpoints+2))
!
!  setup initial conditions
!
       old => temperature(:,1)
       new => temperature(:,2)
       time = 0.
       x = locxleft + [((i-1)*dx,i=1,locpoints+2)]
       old    = ao*exp(-(x)**2 / (2.*sigmao**2))
       theory= ao*exp(-(x)**2 / (2.*sigmao**2))

       call pgbeg(0, "/xwindow", 1, 1)
       call pgask(0)
       call pgenv(locxleft, locxright, 0., 1.5*ao, 0, 0)
       call pglab('x', 'Temperature', 'Diffusion Test')
       red = 2
       call pgscr(red,1.,0.,0.)
       grey = 3
       call pgscr(grey,.2,.2,.2)
       white = 4
       call pgscr(white,1.,1.,1.)
 
 
       call pgsls(1)
       call pgsci(white)
       call pgline(locpoints, x, theory(2:locpoints+1))
       call pgsls(2)
       call pgsci(red)
       call pgline(locpoints, x(2:locpoints+1), old(2:locpoints+1))

!
!  evolve
!
       do step=1, nsteps
!
! boundary conditions: keep endpoint temperatures fixed.
!
           new(1) = old(1)
           new(locpoints+2) = old(locpoints+2)

           call MPI_Sendrecv(old(locpoints+1), 1, MPI_REAL, right, righttag,  &
                     old(1), 1, MPI_REAL, left,  righttag, MPI_COMM_WORLD, rstatus, ierr)

           call MPI_Sendrecv(old(2), 1, MPI_REAL, left, lefttag,  &
                     old(locpoints+2), 1, MPI_REAL, right,  lefttag, MPI_COMM_WORLD, rstatus, ierr)
!
! update solution
!
           forall (i=2:locpoints+1)
               new(i) = old(i) + dt*kappa/(dx**2) * &
                         (old(i+1) - 2*old(i) + old(i-1))
           end forall
           time = time + dt

           if (mod(step, plotsteps) == 0) then
                call pgbbuf()
                call pgeras
                call pgsls(2)   !dashed
                call pgslw(12)  !thick
                call pgsci(red)
                call pgline(locpoints, x(2:locpoints+1), new(2:locpoints+1))
            endif
! 
! update correct solution
!
           sigma = sqrt(2.*kappa*time + sigmao**2)
           a = ao*sigmao/sigma
           theory = a*exp(-(x)**2 / (2.*sigma**2))

           if (mod(step, plotsteps) == 0) then
               call pgsls(1)   
               call pgslw(6)   
               call pgsci(white)
               call pgline(locpoints, x(2:locpoints+1), theory(2:locpoints+1))
               call pgebuf()
           endif
           error = sqrt(sum((theory(1:locpoints+1) - new(1:locpoints+1))**2))
           print *, 'Step = ', step, 'Time = ', time, ' Err = ', error

           tmp => old
           old => new
           new => tmp
       enddo

       call pgend
       deallocate(temperature)
       deallocate(theory)
       deallocate(x)
       call MPI_Finalize(ierr)
       end program diffuse
