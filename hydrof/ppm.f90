module ppm
  implicit none

  private :: color
  public :: outputppm

contains

pure function color(val, lo, hi) result(rgb)
  implicit none
  character, dimension(3) :: rgb
  real, intent(IN) :: val, lo, hi
  integer :: tmp

  real :: scaleval 

  scaleval = (val-lo)/(hi-lo) * 255.
  tmp = int(scaleval*3)
  if (tmp > 255) tmp = 255
  rgb(1) = achar(tmp)
  rgb(2) = achar(int(scaleval))
  rgb(3) = achar(0)
end function color

subroutine outputppm(var,filename,varnum)
  implicit none
  character(len=*), intent(IN), optional :: filename
  integer, intent(IN), optional :: varnum
  real, dimension(:,:,:), intent(IN) :: var
  character, dimension(:,:,:), allocatable :: rgb
  real :: lo, hi
  integer :: i,j,c
  integer :: nx, ny, nvars
  character(len=30) :: chosenfilename
  integer :: chosenvarnum

  if (present(varnum)) then
      chosenvarnum = varnum
  else
      chosenvarnum = 1
  endif
  if (present(filename)) then
      chosenfilename = filename
  else
      chosenfilename = 'output.ppm'
  endif

  nvars = size(var,1)
  nx = size(var,2)
  ny = size(var,3)

  if (chosenvarnum > nvars) then
      print *, 'Invalid variable number ', chosenvarnum
  else
      allocate(rgb(3,nx,ny))

      lo = minval(var(chosenvarnum,:,:))
      hi = maxval(var(chosenvarnum,:,:))

      forall (i=1:nx,j=1:ny)
          rgb(1:3,i,j) = color(var(chosenvarnum,i,j), lo, hi)
      end forall
  
      open(unit=44,file=chosenfilename)
      write(44,'(A)') 'P6'
      write(44,'(2(1x,i4))') nx, ny
      write(44,'(i3)'), 255
      write(44,'(9999999A1)',advance='no') (((rgb(c,i,j),c=1,3), i=1,nx),j=ny,1,-1)
     
      deallocate(rgb)
      close(44)
  endif
end subroutine outputppm

end module ppm
