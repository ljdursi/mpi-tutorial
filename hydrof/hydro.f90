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
  use ppm, only : outputppm
  use solver, only : nvars, timestep, idens, initialconditions, nguard
  use plot, only : openplot, showplot, closeplot

  implicit none
  integer :: n,nx,ny
  real  :: t,dt
  integer :: iter
  real, allocatable, dimension(:,:,:) :: u
  character(len=30) :: arg

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
  nx = n+2*nguard   ! boundary condition zones on each side
  ny = n+2*nguard
  allocate(u(nvars,nx,ny))

  call initialconditions(u)
  call outputppm(u,'ics.ppm',idens)
  call openplot(nx, ny)
  t=0
  timesteps: do iter=1,nx*12
      call timestep(u,dt)
      t = t + 2*dt
      if (mod(iter,10) == 1) then 
        print *, iter, 'dt = ', dt, ' t = ', t
        call showplot(u)
      endif 
  end do timesteps
  call outputppm(u,'dens.ppm',idens)

  deallocate(u)
end
