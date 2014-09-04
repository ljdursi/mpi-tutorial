module solver
  use mpi
  implicit none

  integer, parameter :: idens=1, iener=4, imomx=2, imomy=3, nvars=4
  integer, parameter :: nguard = 2
  real, parameter :: gamma=5./3.

  public :: initialconditions, timestep, nguard, gamma
  private:: cfl, xytranspose, xsweep, tvd1d, hydroflux

contains

subroutine periodicBCs(u,xydim)
    implicit none
    real, dimension(:,:,:), intent(INOUT) :: u
    character, intent(IN) :: xydim
    integer :: n

    if (xydim == 'y') then
        n = size(u,3)
        u(:,:,1:2)   = u(:,:,n-3:n-2)
        u(:,:,n-1:n) = u(:,:,3:4)
    else if (xydim == 'x') then
        n = size(u,2)
        u(:,1:2,:)   = u(:,n-3:n-2,:)
        u(:,n-1:n,:) = u(:,3:4,:)
    endif
end subroutine periodicBCs

subroutine bufferGuardcells(u,xydim,commcart,left,right)
    real, dimension(:,:,:), intent(INOUT) :: u
    character, intent(IN) :: xydim
    integer, intent(IN) :: commcart,left,right
    integer :: nx, ny
    integer :: datasize
    real, dimension(:), allocatable:: sendbuffer, rcvbuffer
    integer :: ierr
    integer :: righttag, lefttag
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    nx = size(u,2) 
    ny = size(u,3) 

    lefttag = 1
    righttag = 2

    datasize = (ny-2*nguard)*nguard*nvars
    allocate(sendbuffer(datasize),rcvbuffer(datasize))

    ! rightward
    ! Copy data into buffer
    !
    ! call MPI_Sendrecv(sendbuffer, datasize,  ..?
    !                   rcvbuffer,  datasize,  ..?
    ! copy data out of buffer.

    ! leftward
    ! Copy data into buffer
    !
    ! call MPI_Sendrecv(sendbuffer, datasize,  ..?
    !                   rcvbuffer,  datasize,  ..?
    ! copy data out of buffer.


    deallocate(sendbuffer,rcvbuffer)
end subroutine bufferGuardcells

subroutine vectorGuardcells(u,xydim,commcart,left,right)
    real, dimension(:,:,:), intent(INOUT) :: u
    character, intent(IN) :: xydim
    integer, intent(IN) :: commcart, left, right
    integer :: nx, ny
    integer :: xbctype
    integer :: lefttag, righttag
    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: rstatus

    nx = size(u,2) 
    ny = size(u,3) 

    lefttag = 1
    righttag = 2

    ! call MPI_Type_vector(...?)
    call MPI_Type_commit(xbctype, ierr)

    ! call MPI_Sendrecv(...?
    ! call MPI_Sendrecv(...?

    call MPI_Type_free(xbctype, ierr)
end subroutine vectorGuardcells

subroutine timestep(u,dt,commcart,left,right)
    implicit none
    real, dimension(:,:,:), intent(INOUT) :: u
    real, intent(OUT) :: dt
    integer, intent(in) :: commcart,left,right

    real, dimension(nvars,size(u,3),size(u,2)) :: ut

    dt=0.25*cfl(u)
! the x sweep
    call periodicBCs(u,'x')
    !call vectorGuardcells(u,'x',commcart,left,right)
    call xsweep(u,dt)
! the y sweeps
    call xytranspose(ut,u)
    call periodicBCs(ut,'x')
    call xsweep(ut,dt)
    call periodicBCs(ut,'x')
    call xsweep(ut,dt)
! 2nd x sweep
    call xytranspose(u,ut)
    call periodicBCs(u,'x')
    !call vectorGuardcells(u,'x',commcart,left,right)
    call xsweep(u,dt)
end subroutine timestep

real function cfl(u)
   implicit none
   real, intent(IN) :: u(:,:,:)
   real, dimension(size(u,2),size(u,3)) :: v, p
   real :: c
   real :: globc
   integer :: ierr

   v = abs(u(imomx,:,:)/u(idens,:,:))
   where (v < abs(u(imomy,:,:))/u(idens,:,:))
      v = abs(u(imomy,:,:)/u(idens,:,:))
   endwhere

   p = (u(iener,:,:)-0.5*(u(imomx,:,:)**2+u(imomy,:,:)**2)/u(idens,:,:))*(gamma-1.)
   c = maxval(v+sqrt(abs((gamma*p))/u(idens,:,:))) 
   
   ! find global max of all c's
   call MPI_Allreduce(c, globc, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
   cfl=1./globc
end function cfl

subroutine initialconditions(u,n,startn)
  implicit none
  real, intent(OUT)   :: u(:,:,:)
  integer, intent(IN) :: n, startn
  integer :: nx, ny
   
  integer :: i,j
  real, parameter :: backgrounddens=1., projdens=50., projvel=4., p0=1.
  real :: r(size(u,2), size(u,3))

  nx = size(u,2)
  ny = size(u,3)
 
  forall (i=1:nx,j=1:ny)
      r(i,j) =sqrt((j*1-n/2.)**2. + ((i+startn-1)*1.-n/2.)**2)
  endforall

  where (r < 0.1*sqrt(n*n+n*n*1.))
    u(idens,:,:) =projdens
    u(imomx,:,:) =projdens*projvel
    u(imomy,:,:) =0
    u(iener,:,:) =0.5*(projdens*projvel*projvel)+p0/(gamma-1.)
  elsewhere
    u(idens,:,:) =backgrounddens
    u(imomx,:,:) =0.
    u(imomy,:,:) =0.
    u(iener,:,:) =p0/(gamma-1.)
  endwhere

end subroutine initialconditions

subroutine xsweep(u,dt)
  implicit none
  real, intent(INOUT), dimension(:,:,:) :: u
  real, intent(IN) :: dt
  integer :: j

  do j=1,size(u,3)
     call tvd1d(u(:,:,j),dt)
  enddo
end subroutine xsweep

subroutine tvd1d(u,dt)
  implicit none
  real, intent(IN) :: dt
  real, dimension(:,:), intent(INOUT) :: u

  real, dimension(size(u,1),size(u,2)) :: v, u1, wr, wl
  real, dimension(size(u,1),size(u,2)) :: dfrp, dfrm, dflm, dflp
  real, dimension(size(u,1),size(u,2)) :: dfl, dfr
  real, dimension(size(u,1),size(u,2)) :: flux
  real, dimension(size(u,1),size(u,2)) :: du
  real :: c
  real, parameter :: smallp = 1.e-3
  integer :: i, n

  n = size(u,2)

  call hydroflux(v,c,u)
  wr = u+v
  wl = u-v
  flux=c*(wr - eoshift(wl,1,dim=2))/2.
  u1 = u - (flux - eoshift(flux,-1,dim=2))*dt/2.

  call hydroflux(v,c,u1)
  wr = u1+v
  wl = u1-v

  dfrp = c*(eoshift(wr,1,dim=2)-wr)/2.
  dfrm = c*(wr-eoshift(wr,-1,dim=2))/2.
  where (dfrp*dfrm > 0) 
     dfr = 2*dfrp*dfrm/(dfrp+dfrm)
  elsewhere
     dfr = 0
  endwhere
  
  dflp = c*(eoshift(wl,1,dim=2)-eoshift(wl,2,dim=2))/2.
  dflm = c*(wl-eoshift(wl,1,dim=2))/2.
  where (dflp*dflm > 0) 
     dfl = 2*dflp*dflm/(dflp+dflm)
  elsewhere
     dfl = 0
  endwhere
  
  flux = (c*(wr-eoshift(wl,1,dim=2)) + (dfr-dfl))/2.
  du = (flux-eoshift(flux,-1,dim=2))*dt
  u(:,3:n-2) = u(:,3:n-2) - du(:,3:n-2)

  do i=1,n
    if ( ((u(iener,i)-0.5*(u(imomx,i)**2+u(imomy,i)**2)/u(idens,i))*(gamma-1.)) < smallp ) &
              u(iener,i) = 0.5*(u(imomx,i)**2+u(imomy,i)**2)/u(idens,i) + smallp/(gamma-1.)
  enddo
end subroutine tvd1d

subroutine hydroflux(v,c,u)
  implicit none
  real, intent(OUT)    :: c
  real, dimension(:,:), intent(IN) :: u
  real, dimension(:,:), intent(OUT) :: v
  real, dimension(size(u,2)) :: vx, p

  c  = 0.
  vx = u(imomx,:)/u(idens,:)
  p  = (u(iener,:)-(u(imomx,:)**2+u(imomy,:)**2)/u(idens,:)/2)*(gamma-1)

  v(idens,:)=u(imomx,:)
  v(imomx,:)=u(imomx,:)*vx + p
  v(imomy,:)=u(imomy,:)*vx
  v(iener,:)=(u(iener,:)+p)*vx

  c=maxval(abs(vx)+sqrt(abs(gamma*p/u(idens,:))))
  if (c > 0) v = v/c
end subroutine hydroflux

subroutine xytranspose(ut,u)
  implicit none
  real, dimension(:,:,:), intent(IN) :: u
  real, dimension(:,:,:), intent(OUT) :: ut

  ut(idens,:,:)=transpose(u(idens,:,:))
  ut(imomx,:,:)=transpose(u(imomy,:,:))
  ut(imomy,:,:)=transpose(u(imomx,:,:))
  ut(iener,:,:)=transpose(u(iener,:,:))
end subroutine xytranspose

end module solver
