module plot

    integer, parameter, private :: ncontours=5
    public :: openplot, showplot, closeplot

contains
    subroutine openplot(nx, ny)
    use solver, only : nguard, idens, iener, imomx, imomy
    implicit none
    integer, intent(in) :: nx, ny
    real :: xl, xr, yl, yr

    xl = 1.
    yl = 1.
    xr = nx - 2*nguard 
    yr = ny - 2*nguard 

    call pgbeg(0, "/XWINDOW", 1, 1);
    call pgask(0)
    call pgenv(xl, xr, yl, yr, 1, 1)
    end subroutine openplot


    subroutine showplot(u) 
    use solver
    implicit none
    real, dimension(:,:,:), intent(IN) :: u
    integer :: nx, ny

    real, dimension(ncontours) :: denscontours, prescontours
    real, dimension(:,:), allocatable :: dens, pres
    real :: maxd, mind, maxp, minp
    real,dimension(6) :: tr
    integer :: i

    nx = size(u,2)
    ny = size(u,3)

    allocate(dens((nx-2*nguard),(ny-2*nguard)))
    allocate(pres((nx-2*nguard),(ny-2*nguard)))
    dens = u(idens, nguard+1:nx-nguard, nguard+1:ny-nguard)
    pres = u(iener, nguard+1:nx-nguard, nguard+1:ny-nguard) - &
              0.5*((u(imomx, nguard+1:nx-nguard, nguard+1:ny-nguard)**2)/dens + &
                   (u(imomy, nguard+1:nx-nguard, nguard+1:ny-nguard)**2)/dens)
    mind = minval(dens)
    maxd = maxval(dens)
    minp = minval(pres)
    maxp = maxval(pres)

    tr = [0.,1.,0.,0.,0.,1.]

    denscontours = [(mind+(i*(maxd-mind)/(ncontours+1)), i=1,ncontours)]
    prescontours = [(minp+(i*(maxp-minp)/(ncontours+1)), i=1,ncontours)]

    call pgenv(1., nx-2.*nguard , 1., nx-2.*nguard, 1, 1)
    call pgsci(2)
    call pgcont(dens, nx-2*nguard, ny-2*nguard, 1, nx-2*nguard, 1, ny-2*nguard, &
                denscontours, ncontours, tr)
    if (minp /= maxp) then
        call pgsci(3)
        call pgcont(pres, nx-2*nguard, ny-2*nguard, 1, nx-2*nguard, 1, ny-2*nguard, &
                    prescontours, ncontours, tr)
    endif
    call pgsci(1)
    deallocate(dens,pres)
    end subroutine showplot

    subroutine closeplot
    call pgend()
    end subroutine closeplot

end module plot
