module tmp
  implicit none
  real(kind=8) :: ang
end module
program main
  use tmp
  implicit none
  real(kind=8) :: minvalue
  real(kind=8),external :: PotenPack
  real(kind=8),parameter :: tol = 1.0D-8
  integer,parameter :: Nvar = 1
  real(kind=8),parameter :: maxstep(Nvar) = 1.0D-1
  real(kind=8) :: var(Nvar)
  real(kind=8) :: xmin(Nvar)
  real(kind=8) :: ynewlo
  real(kind=8), parameter :: DegtoRad = 1.74532925199D-2
  real(kind=8), parameter :: RadtoDeg = 5.72957795131D1
  real(kind=8) :: anglist(91) = (/ (real(i), i = 0, 90) /)
  integer :: i
  integer :: konvge = 3
  integer :: kcount = 200
  integer :: icount, numres, ifault
  
  anglist = anglist * 2.0D0 * DegtoRad
  call ReadParameters()
  var(1) = 4.0D0
    
  write(*,*) "Ang   Rad    Potential"
  do i = 1, 91
    ang = anglist(i)
    call nelmin ( PotenPack, Nvar, var, xmin, ynewlo, tol, &
           maxstep, konvge, kcount, icount, numres, ifault )
    if (ifault==0) then
       write(*,"(F4.0,' ',F5.3,' ',F7.2)") ang*RadtoDeg,xmin,ynewlo
    else
      write(*,*) "found error"
      write(*,*) ifault
      stop
    end if
  end do
  
end program

function PotenPack(var)
  use tmp
  use FitFuncJac
  use FitPara
  implicit none
  real(kind=8) :: PotenPack
  real(kind=8) :: var(1)
  
  PotenPack = Potential(var(1),ang,para_g,para_b,para_d,para_C)
  
end function

subroutine ReadParameters()
  use FitPara
  implicit none
  real(kind=8) :: x(40)
  integer :: i
  
  open(12, file="parameters.txt")
  do i = 1, 40
    read(12,*) x(i)
  end do
  close(12)
  
  para_g = reshape(x(1:24),(/4,6/))
  para_d = x(25:30)
  para_b = x(31:36)
  para_C = x(37:40)

end subroutine
