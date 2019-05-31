program main
  use FitFuncJac, m => nydat
  implicit none
  integer, parameter :: n = 40 ! 40 parameters 
  real(kind=8) :: fvec(m)
  real(kind=8) :: x(n)
  real(kind=8) :: para_g(4,6), para_d(6), para_b(6), para_C(4)
  real(kind=8) :: rdat(nrdat) ! r = 2.7, 7.0, 0.1
  real(kind=8) :: adat(nadat) ! theta = 0, 180, 15
  real(kind=8) :: rmsd
  integer :: i
  
  call ReadYData()
  
  open(12, file="parameters.txt")
  do i = 1, 40
    read(12,*) x(i)
  end do
  close(12)
  
  do i = 1, nrdat
    rdat(i) = 2.7D0 + 0.1D0 * (i-1)
  end do
  do i = 1, nadat
    adat(i) = 0.0D0 + 1.5D1 * (i-1)
  end do

  para_g = reshape(x(1:24),(/4,6/))
  para_d = x(25:30)
  para_b = x(31:36)
  para_C = x(37:40)
  call FuncArray(fvec,para_g,para_d,para_b,para_C,rdat,adat,ydat)
  call r8vec_print ( m, fvec, '  F(X):' )
  fvec = fvec * fvec
  rmsd = sqrt(sum(fvec)/m)
  write(*,*) ' rmsd = ', rmsd

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