!
! F77-ized version of a Fortran90 TVD split MHD code
!
!  Original F90 code  copyright (C) 2001,2003, Ue-Li Pen
!  written November 2001 by Ue-Li Pen, pen@cita.utoronto.ca
!
! see http://arxiv.org/abs/astro-ph/0305088
! or http://www.cita.utoronto.ca/~pen/MHD
!
! This code is licencensed under the GPL
!
program hydro
  use mpi
  use ppm, only : outputppm, outputppm_mpiio
  use solver, only : nvars, timestep, idens, initialconditions, nguard
  use plot, only : openplot, showplot, closeplot
  implicit none

  integer :: n,nx,ny
  real  :: t,dt
  integer :: iter
  real, allocatable, dimension(:,:,:) :: u
  character(len=30) :: arg
  character(len=3) :: rankstr

  integer :: ierr
  integer :: rank, comsize
  integer :: commcart
  integer :: left, right
  integer :: gridcoords(2)
  integer :: startn, endn, locnx

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
  call MPI_Cart_create(MPI_COMM_WORLD, 2, [comsize,1],  &
                       [.true.,.true.], .true., commcart, ierr)
  call MPI_Comm_rank(commcart, rank, ierr)
  call MPI_Cart_shift(commcart, 0, 1, left, right, ierr)
  call MPI_Cart_coords(commcart, rank, 2, gridcoords, ierr)

  write(rankstr,fmt='(I0)') rank
  ! get domain size from command line; if nothing
  ! on command line, prompt user

  if (command_argument_count() < 1) then
      print *,'Enter size of domain, in pixels '
      do 
          read *, n
          if (n < 2 .or. n>500) then
            print *, 'invalid n = ', n, ' try again.'
          else
            exit
      endif
      enddo 
  else
      call get_command_argument(1, arg)
      read( arg,'(I30)'), n
      if (n < 2 .or. n > 500) then
         print *,'invalid n = ', trim(arg), ' using 100'
         n = 100
      endif
  endif

  ! calculate an nx and an ny from n, rank
  locnx = n/comsize
  startn= rank*locnx+1
  endn  = (rank+1)*locnx
  if (gridcoords(1) == comsize-1) endn = n
  locnx = endn - startn + 1

  nx = locnx+2*nguard   
  ny = n+2*nguard
  allocate(u(nvars,nx,ny))

  call initialconditions(u,n,startn)
  call outputppm(u,trim(rankstr)//'-ics.ppm',idens)
  !call outputppm_mpiio(u,n,startn,'ics.ppm',idens)
  call openplot(locnx, ny)
  t=0
  timesteps: do iter=1,n*12
      call timestep(u,dt,commcart,left,right)
      t = t + 2*dt
      if (mod(iter,10) == 1) then 
        if (rank == 0) print *, iter, 'dt = ', dt, ' t = ', t
        call showplot(u, startn)
      endif 
  end do timesteps
  call outputppm(u,trim(rankstr)//'-dens.ppm',idens)
  !call outputppm_mpiio(u,n,startn,'dens.ppm',idens)

  deallocate(u)
  call MPI_Finalize(ierr)
end
