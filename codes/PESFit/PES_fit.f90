program main
  use FitFuncJac, only : m => nydat
  implicit none
  integer,parameter :: n = 40
  external :: Func, Jac, Avv, callback
  real(kind=8) :: x(n)
  real(kind=8) :: fvec(m)
  real(kind=8) :: fjac(m,n)
  integer :: info = 0
  logical :: analytic_jac = .FALSE., analytic_Avv = .FALSE.
  logical :: center_diff = .TRUE.
  real(kind=8) :: h1 = 5.0D-2, h2 = 5.0D-2
  real(kind=8) :: dtd(n,n)
  integer :: damp_mode = 0
  integer :: niteres, nfev, njev, naev
  integer :: maxiter = 500
  integer :: maxfev = 0, maxjev = 0, maxaev = 0
  real(kind=8) :: maxlam = 1.0D16, minlam = -1.0D0
  real(kind=8) :: artol = 1.0D-7, Cgoal = 5.0D-2
  real(kind=8) :: gtol = 1.0D-11
  real(kind=8) :: xtol = 1.0D-11
  real(kind=8) :: xrtol = 1.0D-11
  real(kind=8) :: ftol = 1.0D0-11
  real(kind=8) :: frtol = 1.0D-11
  integer :: converged
  integer :: print_level = 2
  integer :: print_unit = 6
  integer :: imethod = 0
  integer :: iaccel = 1
  integer :: ibold = 4
  integer :: ibroyden = 0
  real(kind=8) :: initialfactor = 1.0D0
  real(kind=8) :: factoraccept = 15.0D0
  real(kind=8) :: factorreject = 2.0D0
  real(kind=8) :: avmax = 0.75D0
  integer :: i
  
  call ReadData()
  
  ! initial guess for parameters  
  ! based on the results of 10.1021/jp303573a
  open(12, file="guess.txt")
  do i = 1, 40
    read(12,*) x(i)
  end do
  close(12)
    
  call geodesicLM(Func, Jac, Avv, x, fvec, fjac, n, m, callback, info, &
          analytic_jac, analytic_Avv, center_diff, h1, h2, &
          dtd, damp_mode, niteres, nfev, njev, naev, &
          maxiter, maxfev, maxjev, maxaev, maxlam, minlam, &
          artol, Cgoal, gtol, xtol, xrtol, ftol, frtol, &
          converged, print_level, print_unit, &
          imethod, iaccel, ibold, ibroyden, &
          initialfactor, factoraccept, factorreject, avmax)

  call WriteParameters(x)
  call RMSD(x)
  
end program
  
subroutine Func( m, n, x, fvec)
  use FitFuncJac
  implicit none
  integer :: m, n
  real(kind=8) :: fvec(m)
  real(kind=8) :: x(n)
  real(kind=8) :: para_g(4,6), para_d(6), para_b(6), para_C(4)

  para_g = reshape(x(1:24),(/4,6/))
  para_d = x(25:30)
  para_b = x(31:36)
  para_C = x(37:40)
  call FuncArray(fvec,para_g,para_d,para_b,para_C)

end subroutine

subroutine Jac ( m, n, x, fjac)
  use FitFuncJac
  implicit none
  integer :: m, n
  real(kind=8) :: fjac(m,n)
  real(kind=8) :: x(n)
  real(kind=8) :: para_g(4,6), para_d(6), para_b(6), para_C(4)

  para_g = reshape(x(1:24),(/4,6/))
  para_d = x(25:30)
  para_b = x(31:36)
  para_C = x(37:40)
  call JacMat(fjac,para_g,para_d,para_b,para_C)

end subroutine

subroutine Avv ( m, n, x, v, acc)
  implicit none
  integer :: m, n
  real(kind=8) :: x(n), v(n), acc(m) 
  ! do nothing
end subroutine

subroutine callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
  implicit none
  integer :: m, n, accepted, info
  real(kind=8) :: x(n), v(n), a(n), fvec(m), fjac(m,n), acc(m), lam, dtd(n,n), fvec_new(m)
  ! do nothing
end subroutine callback

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