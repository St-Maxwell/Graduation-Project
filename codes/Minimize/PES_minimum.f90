program main
  implicit none
  real(kind=8) :: minvalue
  real(kind=8),external :: praxis, PotenPack
  real(kind=8),parameter :: tol = 5.0D-7
  real(kind=8),parameter :: maxstep = 1.0D-2
  integer,parameter :: Nvar = 2
  integer :: prin = 3
  real(kind=8) :: var(Nvar)
  real(kind=8), parameter :: DegtoRad = 1.74532925199D-2
  real(kind=8), parameter :: RadtoDeg = 5.72957795131D1
  
  var(1) = 4.0D0
  var(2) = 9.0D1 * DegtoRad
  write(*,*) " Initial Guess:"
  write(*,*) " R = 4.0 Ang"
  write(*,*) " Theta = 90 Deg"
  call ReadParameters()
  minvalue = praxis(tol, maxstep, Nvar, prin, var, PotenPack)
  write(*,"(' PES Minimum = ',F7.2,' cm-1')") minvalue
  write(*,"(' R = ',F5.3,' Ang','  Theta = ',F5.2,' Deg')") &
    var(1), var(2) * RadtoDeg
  
end program

function PotenPack(var,n)
  use FitFuncJac
  use FitPara
  implicit none
  real(kind=8) :: PotenPack
  integer :: n
  real(kind=8) :: var(n)
  
  PotenPack = Potential(var(1),var(2),para_g,para_b,para_d,para_C)
  
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
