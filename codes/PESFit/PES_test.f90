program main
  use FitFuncJac, m => nydat
  use lapack95
  implicit none
  integer, parameter :: n = 40 ! 40 parameters 
  real(kind=8) :: fvec(m)
  real(kind=8) :: fjac(m,40)
  real(kind=8) :: x(n)
  real(kind=8) :: para_g(4,6), para_d(6), para_b(6), para_C(4)
  real(kind=8) :: rdat(nrdat) ! r = 2.7, 7.0, 0.1
  real(kind=8) :: adat(nadat) ! theta = 0, 180, 15
  real(kind=8) :: ipiv(n)
  real(kind=8) :: A(n,n), B(n)
  integer :: i, j, info
  
  real(kind=8) :: ang,rad,V
  
  call ReadGuess(x)
  para_g = reshape(x(1:24),(/4,6/))
  para_d = x(25:30)
  para_b = x(31:36)
  para_C = x(37:40)
  ang = 0.0D0
  do i = 1, 5
    rad = 2.7D0
    do j = 1, 21
      V = Potential(rad,ang,para_g,para_b,para_d,para_C)
      write(*,"(F4.0,' ',F4.1,' ',F9.2)") ang,rad,V
      rad = rad + 0.1D0
    end do
    ang = ang + 45.0D0
  end do 
  
  stop
  
  
  !do i = 1, nrdat
  !  rdat(i) = 2.7D0 + 0.1D0 * (i-1)
  !end do
  !do i = 1, nadat
  !  adat(i) = 0.0D0 + 1.5D1 * (i-1)
  !end do
  !call ReadYData()
  !call ReadGuess(x)
  !
  !do i = 1, 10
  !  para_g = reshape(x(1:24),(/4,6/))
  !  para_d = x(25:30)
  !  para_b = x(31:36)
  !  para_C = x(37:40)
  !  
  !  call FuncArray(fvec,para_g,para_d,para_b,para_C,rdat,adat,ydat)
  !  call JacMat(fjac,para_g,para_d,para_b,para_C,rdat,adat)
  !  
  !  A = matmul(transpose(fjac),fjac)
  !  B = -matmul(transpose(fjac),fvec)
  !  call dgesv( n, 1, A, n, ipiv, B, n, info )
  !  x = x + B
  !  do j = 1, 40
  !    write(*,*) x(i)
  !  end do
  !end do
    
end program

subroutine ReadYData()
  use FitFuncJac, only : nrdat, nadat, ydat
  integer :: i, node1, node2
  open(12, file="pesdata.txt")
  do i = 1, nrdat
    node1 = 1 + nadat * (i-1)
    node2 = i * nadat
    read(12,*) ydat(node1:node2)
  end do
  close(12)
end subroutine

subroutine ReadGuess(x)
  implicit none
  real(kind=8) :: x(40)
  integer :: i
  ! initial guess for parameters  
  ! based on the results of 10.1021/jp303573a
  open(12, file="parameters.txt")
  do i = 1, 40
    read(12,*) x(i)
  end do
  close(12)
end subroutine