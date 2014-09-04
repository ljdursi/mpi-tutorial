program SimulateNBody

    implicit none

    type Nbody
        integer :: id
        double precision, dimension(3) :: x
        double precision, dimension(3) :: vel
        double precision, dimension(3) :: force
        double precision :: mass
        double precision :: potentialE
    end type Nbody    

    real :: display = 1.3
    integer :: npts=1500
    integer :: nsteps=500
    integer :: outevery = 5
    integer :: simulation = 1
    integer :: output=1
    integer :: i
    integer :: seed = 5
    double precision :: dt = 0.01
    double precision :: time = 0.
    double precision :: tote = 0.

    type(Nbody), allocatable :: pdata(:)

    allocate(pdata(npts))

    call random_seed(seed)

    call initialize_particles(pdata, npts, simulation)
    call calculate_forces_fastest(pdata, npts)
    call calculate_energy(pdata, npts, tote)

    do i=1,nsteps
        call nbody_step(pdata, npts, dt)
        call calculate_energy(pdata, npts, tote)
        time = time + dt
        if (output /= 0) then
            print *, i, dt, time, tote
            if (mod(i,outevery) == 0) then 
                call display_particles(pdata, npts, display)
            endif
        endif
    enddo

    contains

    subroutine initialize_particles(pdata, n, simulation)
        implicit none
        integer, intent(IN) :: n
        type(Nbody), intent(INOUT), dimension(n) :: pdata
        integer, intent(IN) :: simulation

        integer :: i, j
        double precision, parameter :: vamp = 0.1
        double precision :: r
        double precision :: rnd
        double precision ::  mass

        mass = 1./(1.*n)

        if (simulation == 1) then
            i = 1
            do while (i<=n)
                r = 0.
                do j=1,3
                    call random_number(rnd)
                    pdata(i) % x(j) = 2.0 * rnd - 1.
                    r = r + pdata(i) % x(j)**2
                enddo
                if (r < 1.) then
                    pdata(i) % mass = mass
                    do j=1,3
                        call random_number(rnd)
                        pdata(i) % vel(j) = vamp * ( 2.0 * rnd - 1.)
                    enddo
                    pdata(i) % id = i
                    i = i + 1
                endif
            enddo
        endif
    return
    end subroutine initialize_particles


    subroutine calculate_forces_fastest(pdata, n)
        implicit none
        integer, intent(IN) :: n
        type(Nbody), intent(INOUT), dimension(n) :: pdata
        integer :: i, j, d
        double precision :: rsq
        double precision, parameter :: EPS=0.1
        double precision, parameter :: gravconst=1.
        double precision, dimension(3) :: dx
        double precision :: ir
        double precision :: forcex

        do i=1,n
            pdata(i) % force = 0.
            pdata(i) % potentialE = 0.
        enddo

        do i=1,n
            do j=i+1,n
                rsq = EPS*EPS
                dx = 0.
                do d=1,3
                    dx(d) = pdata(j)%x(d) - pdata(i)%x(d)
                    rsq = rsq + dx(d)*dx(d)
                enddo
                ir = 1./sqrt(rsq)
                rsq = ir/rsq
                do d=1,3
                    forcex = rsq*dx(d) * pdata(i)%mass * pdata(j)%mass * gravconst
                    pdata(i)%force(d) = pdata(i)%force(d) + forcex
                    pdata(j)%force(d) = pdata(j)%force(d) - forcex
                enddo
                pdata(i)%potentialE = pdata(i)%potentialE -             &
                     gravconst * pdata(i)%mass * pdata(j)%mass * ir  
                pdata(j)%potentialE = pdata(i)%potentialE -             &
                     gravconst * pdata(i)%mass * pdata(j)%mass * ir
            enddo
        enddo

        return
    end subroutine calculate_forces_fastest

    subroutine calculate_energy(pdata, n, energy)
        implicit none
        integer, intent(IN) :: n
        type(NBody), intent(INOUT), dimension(n) :: pdata
        double precision, intent(OUT) :: energy
        
        integer :: i

        energy = 0.
        do i=1,n
           energy = energy + pdata(i)%potentialE  + &
            0.5 * sum(pdata(i)%vel*pdata(i)%vel) * pdata(i)%mass
        enddo

        return
    end subroutine calculate_energy



    subroutine nbody_step(pdata, n, dt)
        implicit none
        integer, intent(IN) :: n
        double precision, intent(IN) :: dt
        type(NBody), intent(INOUT), dimension(n) :: pdata

        integer :: i
        double precision, dimension(3) :: accel

        ! kick
        do i=1,n
            pdata(i)%x = pdata(i)%x + pdata(i)%vel*dt
        enddo

        ! drift
        call calculate_forces_fastest(pdata, n)

        do i=1,n
            accel = pdata(i)%force / pdata(i)%mass
            pdata(i)%vel = pdata(i)%vel + accel*dt
        enddo
       
        return
    end subroutine nbody_step


    subroutine display_particles(pdata, n, cmax)
        implicit none
        integer, intent(IN) :: n
        type(NBody), intent(IN), dimension(n) :: pdata
        real, intent(IN) :: cmax

        integer, parameter :: npix=256
        real, dimension(6) :: tr
        real, dimension(npix,npix) :: plot
        logical, save :: firsttime = .true.

        integer :: i, ii, jj

        if (firsttime) then
            call pgbeg(0, "/xwindow", 1, 1)
            call pgask(0)
            call pgenv(1.,npix*1.,1.,npix*1.,1,2)
            call pglab('x','y','View from top')
            firsttime = .false.
        endif

        tr = (/ 0,1,0,0,0,1 /)
        plot = 0.

        do i=1,n
            ii = npix*int((pdata(i)%x(1) + cmax)/(2.*cmax))
            jj = npix*int((pdata(i)%x(2) + cmax)/(2.*cmax))
            if ( (ii > 0) .and. (ii <= npix) .and. (jj > 0) .and. (jj <= npix) ) then
                plot(ii,jj) = plot(ii,jj) - 1.
            endif
        enddo

        call pgimag(plot,npix,npix,1,npix,1,npix,-5.0,.0,tr)
        return
    end subroutine display_particles  

end program 

