program main
  use FitFuncJac, only : m => nydat
  implicit none
  integer, parameter :: n = 40 ! 40 parameters 
  integer, parameter :: ldfjac = m
  external :: FuncJac
  real(kind=8) :: fjac(ldfjac,n)
  real(kind=8) :: fvec(m)
  integer :: iflag
  integer :: info
  real(kind=8), parameter :: tol = 1.0D-14
  real(kind=8) :: x(n)
  integer :: i
  
    call ReadData()
    
    ! initial guess for parameters  
    ! based on the results of 10.1021/jp303573a
    open(12, file="guess.txt")
    do i = 1, 40
      read(12,*) x(i)
    end do
    close(12)
    
    !call r8vec_print ( n, x, '  Initial X:' )
    iflag = 1
    call FuncJac ( m, n, x, fvec, fjac, ldfjac, iflag )
    !call r8vec_print ( m, fvec, '  F(X):' )
    call lmder1 ( FuncJac, m, n, x, fvec, fjac, ldfjac, tol, info )
    
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
    !call r8vec_print ( n, x, '  X:' )
    !call r8vec_print ( m, fvec, '  F(X):' )
    
    call WriteParameters(x)
    call RMSD(x)
  
end program
  
subroutine FuncJac ( m, n, x, fvec, fjac, ldfjac, iflag )
  use FitFuncJac
  implicit none
  integer :: ldfjac
  integer :: m, n
  real(kind=8) :: fjac(ldfjac,n)
  real(kind=8) :: fvec(m)
  integer :: iflag
  real(kind=8) :: x(n)
  real(kind=8) :: para_g(4,6), para_d(6), para_b(6), para_C(4)
  integer :: i
  
  if ( iflag == 0 ) then
    write ( *, '(a)' ) ''
    do i = 1, n
      write ( *, '(g14.6)' ) x(i)
    end do
  else if ( iflag == 1 ) then
    ! x = g, d, b, C
    para_g = reshape(x(1:24),(/4,6/))
    para_d = x(25:30)
    para_b = x(31:36)
    para_C = x(37:40)
    call FuncArray(fvec,para_g,para_d,para_b,para_C)
  else if ( iflag == 2 ) then
    para_g = reshape(x(1:24),(/4,6/))
    para_d = x(25:30)
    para_b = x(31:36)
    para_C = x(37:40)
    call JacMat(fjac,para_g,para_d,para_b,para_C)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LMDER1_F - Fatal error!'
    write ( *, '(a,i6)' ) '  Called with unexpected value of IFLAG = ', iflag
    stop
  end if

end subroutine

subroutine ReadData()
  use FitFuncJac, only : nydat, ydat, rdat, adat, DegtoRad
  integer :: i, node1, node2

  open(12, file="pesdata.txt")
  do i = 1, nydat
    read(12,*) adat(i), rdat(i), ydat(i)
  end do
  close(12)
  adat = DegtoRad * adat
  
end subroutine

subroutine WriteParameters(x)
  implicit none
  real(kind=8) :: x(40)
  integer :: i
  
  open(12, file="parameters.txt")
  do i = 1, 40
    write(12,"(F15.5)") x(i)
  end do
  close(12)

end subroutine

subroutine RMSD(x)
  use FitFuncJac
  implicit none
  real(kind=8) :: x(40)
  real(kind=8) :: para_g(4,6), para_d(6), para_b(6), para_C(4)
  real(kind=8) :: fvec(nydat)
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
  
  call FuncArray(fvec,para_g,para_d,para_b,para_C)
  
  open(12, file="output.txt")
  do i = 1, nydat
    write(12,"(F4.0,' ',F4.1,' ',F9.3,' ',F9.3)") &
      adat(i)*RadtoDeg,rdat(i),ydat(i),fvec(i)
  end do
  
  fvec = fvec * fvec
  write(12,*) "RMSD = ", sqrt(sum(fvec)/nydat)
  close(12)
  
end subroutine