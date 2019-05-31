module FitFuncJac
  implicit none
  !integer, parameter :: nrdat = 44
  !integer, parameter :: nadat = 13
  real(kind=8), parameter :: DegtoRad = 1.74532925199D-2
  real(kind=8), parameter :: RadtoDeg = 5.72957795131D1
  !real(kind=8), parameter :: AngtoBohr = 1.8897260D0
  !real(kind=8), parameter :: BohrtoAng = 0.5291772D0
  integer, parameter :: nydat = 459
  real(kind=8) :: ydat(nydat)
  real(kind=8) :: rdat(nydat)
  real(kind=8) :: adat(nydat)
  contains
  function Legendre(x,degree)
  !---------------------------------------
  ! Legendre polynomials, up to degree 5
  ! Inupt:
  ! x       real*8
  ! degree  integer 
  !---------------------------------------
    implicit none
    real(kind=8) :: Legendre, x
    integer :: degree
    select case(degree)
      case(0)
        Legendre = 1.0D0
      case(1)
        Legendre = x
      case(2)
        Legendre = 1.5D0 * x*x - 0.5D0
      case(3)
        Legendre = x * (2.5D0 * x*x - 1.5D0)
      case(4)
        Legendre = x*x * (4.375D0 * x*x - 3.75D0) + 3.75D-1
      case(5)
        Legendre = x * (x*x * (7.875D0 * x*x - 8.75D0) + 1.875D0)
      case default
        write(*,*) " Degree n Out of Range(n=0-5). "
        stop
    end select
  end function

  function Potential(rad,ang,g,b,d,C)
    implicit none
    real(kind=8) :: Potential, rad, ang, g(4,6), b(6), d(6), C(4)
    real(kind=8), parameter :: mHtocm = 219.47463067D0 ! mHartree to cm-1
    Potential = mHtocm * (V_sh(rad,ang,g,b,d) + V_as(rad,ang,C,b))
  end function

  function V_sh(rad,ang,g,b,d)
    implicit none
    real(kind=8) :: V_sh, rad, ang, g(4,6), b(6), d(6)
    V_sh = V_sh_G(rad,ang,g) * exp(V_sh_D(ang,d) + V_sh_B(ang,b) * rad)
  end function

  function V_sh_G(rad,ang,g)
    implicit none
    real(kind=8) :: V_sh_G, rad, ang, g(4,6)
    real(kind=8) :: costheta
    integer :: i
    costheta = cos(ang)
    V_sh_G = 0.0D0
    do i = 1, 6
      V_sh_G = V_sh_G + (g(1,i)+rad*(g(2,i)+rad*(g(3,i)+g(4,i)*rad))) * Legendre(costheta,i-1)
    end do
  end function

  function V_sh_B(ang,b)
    implicit none
    real(kind=8) :: V_sh_B, ang, b(6)
    real(kind=8) :: costheta
    integer :: i
    costheta = cos(ang)
    V_sh_B = 0.0D0
    do i = 1, 6
      V_sh_B = V_sh_B + b(i) * Legendre(costheta,i-1)
    end do
  end function

  function V_sh_D(ang,d)
    implicit none
    real(kind=8) :: V_sh_D, ang, d(6)
    real(kind=8) :: costheta
    integer :: i
    costheta = cos(ang)
    V_sh_D = 0.0D0
    do i = 1, 6
      V_sh_D = V_sh_D + d(i) * Legendre(costheta,i-1)
    end do
  end function

  function TTDampFunc(x,degree)
  !---------------------------------------
  ! Tang-Toennies damping function
  ! doi: 10.1063/1.481043
  ! Inupt:
  ! x       real*8
  ! degree  integer 
  !---------------------------------------
    implicit none
    real(kind=8) :: TTDampFunc, x
    integer :: degree, i
    if (degree<0) then
      write(*,*) " Degree n Is Expected to Be Non-negative. "
      stop
    end if
    TTDampFunc = 0.0D0
    do i = 0, degree
      TTDampFunc = TTDampFunc + x**i / factorial(i)
    end do
    TTDampFunc = 1.0D0 - TTDampFunc * exp(x)
  end function

  function V_as(rad,ang,C,b)
    implicit none
    real(kind=8) :: V_as, rad, ang, C(4), b(6)
    real(kind=8) :: costheta, BR, temp1, temp2
    real(kind=8) :: f6, f7
    costheta = cos(ang)
    BR = V_sh_B(ang,b)*rad
    f6 = TTDampFunc(BR,6)
    f7 = TTDampFunc(BR,7)
    temp1 = f6 * (C(1)*Legendre(costheta,0) + C(2)*Legendre(costheta,2)) / rad**6
    temp2 = f7 * (C(3)*Legendre(costheta,1) + C(4)*Legendre(costheta,3)) / rad**7
    V_as = temp1 + temp2
  end function

  function factorial(n)
  !---------------------------------------
  ! Factorial function
  ! Input:
  ! n  integer
  ! Return:
  ! factorial  integer
  !---------------------------------------
    implicit none
    integer :: factorial, n
    integer :: i
    if (n < 0) then
      write(*,*) "Negative Input in Factorial Funcion. "
      stop
    else if (n == 0) then
      factorial = 1
    else
      factorial = 1
      do i = 1, n
        factorial = factorial * i
      end do
    end if
  end function

  subroutine FuncArray(fvec,para_g,para_d,para_b,para_C)
  !------------------------------------------------------------
  ! fvec: Function Array f_i(x1,x2,...,xi) - y_i
  !------------------------------------------------------------
  implicit none
  real(kind=8) :: fvec(:), para_g(:,:), para_d(:), para_b(:), para_C(:)
  integer :: i

  do i = 1, nydat
    fvec(i) = Potential(rdat(i),adat(i),para_g,para_b,para_d,para_C)
  end do
  fvec = fvec - ydat

  end subroutine

  subroutine JacMat(fjac,para_g,para_d,para_b,para_C)
  !------------------------------------------------------------
  ! fjac: Jacobi Matrix
  !        | ∂V/∂g00(R1,theta1),  ∂V/∂g10(R1,theta1),  ... |
  !        | ∂V/∂g00(R1,theta2),  ∂V/∂g10(R1,theta2),  ... |
  !        | ∂V/∂g00(R1,theta3),  ∂V/∂g10(R1,theta3),  ... |
  !        |        :                    :                 |
  !        | ∂V/∂g00(R1,theta13), ∂V/∂g10(R1,theta13), ... |
  ! fjac = | ∂V/∂g00(R2,theta1),  ∂V/∂g10(R2,theta1),  ... |
  !        | ∂V/∂g00(R2,theta2),  ∂V/∂g10(R2,theta2),  ... |
  !        | ∂V/∂g00(R2,theta3),  ∂V/∂g10(R2,theta3),  ... |
  !        |        :                    :                 |
  !        | ∂V/∂g00(R2,theta13), ∂V/∂g10(R2,theta13), ... |
  !        |        :                    :                 |
  !------------------------------------------------------------
    implicit none
    real(kind=8) :: fjac(:,:), para_g(:,:), para_d(:), para_b(:), para_C(:)
    real(kind=8) :: JacG(nydat,24), JacD(nydat,6), JacB(nydat,6), JacC(nydat,4)

    call JacVsG(JacG,rdat,adat,para_d,para_b)
    fjac(:,1:24) = JacG

    call JacVsD(JacD,rdat,adat,para_g,para_d,para_b)
    fjac(:,25:30) = JacD

    call JacVsB(JacB,rdat,adat,para_g,para_d,para_b,para_C)
    fjac(:,31:36) = JacB

    call JacVsC(JacC,rdat,adat,para_b)
    fjac(:,37:40) = JacC

  end subroutine

  subroutine JacVsG(JacG,rdat,adat,d,b)
    implicit none
    real(kind=8) :: JacG(:,:), rdat(:), adat(:), d(:), b(:)
    real(kind=8) :: costheta, expDB
    integer :: i
    do i = 1, nydat
      costheta = cos(adat(i))
      expDB = exp(V_sh_D(adat(i),d) + V_sh_B(adat(i),b) * rdat(i))
      JacG(i,1) = Legendre(costheta,0) * expDB
      JacG(i,2) = JacG(i,1) * rdat(i)
      JacG(i,3) = JacG(i,2) * rdat(i)
      JacG(i,4) = JacG(i,3) * rdat(i)
      JacG(i,5) = Legendre(costheta,1) * expDB
      JacG(i,6) = JacG(i,5) * rdat(i)
      JacG(i,7) = JacG(i,6) * rdat(i)
      JacG(i,8) = JacG(i,7) * rdat(i)
      JacG(i,9) = Legendre(costheta,2) * expDB
      JacG(i,10) = JacG(i,9) * rdat(i)
      JacG(i,11) = JacG(i,10) * rdat(i)
      JacG(i,12) = JacG(i,11) * rdat(i)
      JacG(i,13) = Legendre(costheta,3) * expDB
      JacG(i,14) = JacG(i,13) * rdat(i)
      JacG(i,15) = JacG(i,14) * rdat(i)
      JacG(i,16) = JacG(i,15) * rdat(i)
      JacG(i,17) = Legendre(costheta,4) * expDB
      JacG(i,18) = JacG(i,17) * rdat(i)
      JacG(i,19) = JacG(i,18) * rdat(i)
      JacG(i,20) = JacG(i,19) * rdat(i)
      JacG(i,21) = Legendre(costheta,5) * expDB
      JacG(i,22) = JacG(i,21) * rdat(i)
      JacG(i,23) = JacG(i,22) * rdat(i)
      JacG(i,24) = JacG(i,23) * rdat(i)
    end do
  end subroutine

  subroutine JacVsD(JacD,rdat,adat,g,d,b)
    implicit none
    real(kind=8) :: JacD(:,:), rdat(:), adat(:), g(:,:), d(:), b(:)
    real(kind=8) :: costheta, Vsh
    integer :: i
    do i = 1, nydat
      costheta = cos(adat(i))
      Vsh = V_sh(rdat(i),adat(i),g,b,d)
      JacD(i,1) = Vsh * Legendre(costheta,0)
      JacD(i,2) = Vsh * Legendre(costheta,1)
      JacD(i,3) = Vsh * Legendre(costheta,2)
      JacD(i,4) = Vsh * Legendre(costheta,3)
      JacD(i,5) = Vsh * Legendre(costheta,4)
      JacD(i,6) = Vsh * Legendre(costheta,5)
    end do
  end subroutine

  subroutine JacVsB(JacB,rdat,adat,g,d,b,C)
    implicit none
    real(kind=8) :: JacB(:,:), rdat(:), adat(:), g(:,:), d(:), b(:), C(:)
    real(kind=8) :: JacBsh(nydat,6), JacBas6(nydat,6), JacBas7(nydat,6)
    call JacVsB_sh(JacBsh,rdat,adat,g,d,b)
    call JacVsB_as6(JacBas6,rdat,adat,b,C)
    call JacVsB_as7(JacBas7,rdat,adat,b,C)
    JacB = JacBsh + JacBas6 + JacBas7
  end subroutine
  
  subroutine JacVsB_sh(JacBsh,rdat,adat,g,d,b)
    implicit none
    real(kind=8) :: JacBsh(:,:), rdat(:), adat(:), g(:,:), d(:), b(:)
    real(kind=8) :: costheta, VshR
    integer :: i
    do i = 1, nydat
      costheta = cos(adat(i))
      VshR = V_sh(rdat(i),adat(i),g,b,d) * rdat(i)
      JacBsh(i,1) = VshR * Legendre(costheta,0)
      JacBsh(i,2) = VshR * Legendre(costheta,1)
      JacBsh(i,3) = VshR * Legendre(costheta,2)
      JacBsh(i,4) = VshR * Legendre(costheta,3)
      JacBsh(i,5) = VshR * Legendre(costheta,4)
      JacBsh(i,6) = VshR * Legendre(costheta,5)
    end do
  end subroutine
  
  subroutine JacVsB_as6(JacBas6,rdat,adat,b,C)
    implicit none
    real(kind=8) :: JacBas6(:,:), rdat(:), adat(:), b(:), C(:)
    real(kind=8) :: costheta, tmp, JacBasTm1(6), JacBasTm2(6)
    integer :: i
    do i = 1, nydat
      costheta = cos(adat(i))
      tmp = (C(1)*Legendre(costheta,0) + C(2)*Legendre(costheta,2)) / rdat(i)**6
      call JacBasTerm1(JacBasTm1,rdat(i),adat(i),6,b)
      call JacBasTerm2(JacBasTm2,rdat(i),adat(i),6,b)
      JacBas6(i,:) = tmp * (JacBasTm1 + JacBasTm2)
    end do
  end subroutine

  subroutine JacVsB_as7(JacBas7,rdat,adat,b,C)
    implicit none
    real(kind=8) :: JacBas7(:,:), rdat(:), adat(:), b(:), C(:)
    real(kind=8) :: costheta, tmp, JacBasTm1(6), JacBasTm2(6)
    integer :: i
    do i = 1, nydat
      costheta = cos(adat(i))
      tmp = (C(3)*Legendre(costheta,1) + C(4)*Legendre(costheta,3)) / rdat(i)**7
      call JacBasTerm1(JacBasTm1,rdat(i),adat(i),7,b)
      call JacBasTerm2(JacBasTm2,rdat(i),adat(i),7,b)
      JacBas7(i,:) = tmp * (JacBasTm1 + JacBasTm2)
    end do
  end subroutine

  subroutine JacBasTerm1(JacBasTm1,rad,ang,n,b)
    implicit none
    real(kind=8) :: JacBasTm1(6), rad, ang, b(:)
    integer :: n
    real(kind=8) :: costheta, tmp, BR
    integer :: i
    
    BR = V_sh_B(ang, b) * rad
    costheta = cos(ang)
    tmp = 0.0D0
    do i = 0, n
      tmp = tmp + BR**i / factorial(i)
    end do
    tmp = tmp * (-exp(BR)) * rad
    do i = 1, 6
      JacBasTm1(1) = tmp * Legendre(costheta,0)
      JacBasTm1(2) = tmp * Legendre(costheta,1)
      JacBasTm1(3) = tmp * Legendre(costheta,2)
      JacBasTm1(4) = tmp * Legendre(costheta,3)
      JacBasTm1(5) = tmp * Legendre(costheta,4)
      JacBasTm1(6) = tmp * Legendre(costheta,5)
    end do

  end subroutine

  subroutine JacBasTerm2(JacBasTm2,rad,ang,n,b)
    implicit none
    real(kind=8) :: JacBasTm2(6), rad, ang, b(:)
    integer :: n
    real(kind=8) :: costheta, tmp, BR
    integer :: i

    costheta = cos(ang)
    BR = V_sh_B(ang, b) * rad
    !if (BR>0.0D0) then
      tmp = rad / factorial(1)
      tmp = tmp + 2.0D0 * BR * rad / factorial(2)
      tmp = tmp + 3.0D0 * BR**2 * rad / factorial(3)
      tmp = tmp + 4.0D0 * BR**3 * rad / factorial(4)
      tmp = tmp + 5.0D0 * BR**4 * rad / factorial(5)
      tmp = tmp + 6.0D0 * BR**5 * rad / factorial(6)
      if (n==7) then
        tmp = tmp + 7.0D0 * BR**6 * rad / factorial(7)
      end if
    !else
    !  tmp = -rad / factorial(1)
    !  tmp = tmp + 2.0D0 * BR * rad / factorial(2)
    !  tmp = tmp - 3.0D0 * BR**2 * rad / factorial(3)
    !  tmp = tmp + 4.0D0 * BR**3 * rad / factorial(4)
    !  tmp = tmp - 5.0D0 * BR**4 * rad / factorial(5)
    !  tmp = tmp + 6.0D0 * BR**5 * rad / factorial(6)
    !  if (n==7) then
    !    tmp = tmp - 7.0D0 * BR**6 * rad / factorial(7)
    !  end if
    !end if

    tmp = tmp * (-exp(BR))
    do i = 1, 6
      JacBasTm2(1) = tmp * Legendre(costheta,0)
      JacBasTm2(2) = tmp * Legendre(costheta,1)
      JacBasTm2(3) = tmp * Legendre(costheta,2)
      JacBasTm2(4) = tmp * Legendre(costheta,3)
      JacBasTm2(5) = tmp * Legendre(costheta,4)
      JacBasTm2(6) = tmp * Legendre(costheta,5)
    end do

  end subroutine

  subroutine JacVsC(JacC,rdat,adat,b)
    implicit none
    real(kind=8) :: JacC(:,:), rdat(:), adat(:), b(:)
    real(kind=8) :: costheta, BR, f6, f7
    integer :: i
    do i = 1, nydat
      costheta = cos(adat(i))
      BR = V_sh_B(adat(i),b) * rdat(i)
      f6 = TTDampFunc(BR,6)
      f7 = TTDampFunc(BR,7)
      JacC(i,1) = f6 * Legendre(costheta,0) / rdat(i)**6
      JacC(i,2) = f6 * Legendre(costheta,2) / rdat(i)**6
      JacC(i,3) = f7 * Legendre(costheta,1) / rdat(i)**7
      JacC(i,4) = f7 * Legendre(costheta,3) / rdat(i)**7
    end do
  end subroutine

end module