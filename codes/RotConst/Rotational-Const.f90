module global
  implicit none
  real(kind=8),parameter :: rh = 1.054571800D0 ! 1.054 571 800 x 10^-34 J s
  real(kind=8),parameter :: mC = 1.20D-2
  real(kind=8),parameter :: mO = 1.599491461956D-2
  real(kind=8),parameter :: mN = 1.40030740048D-2
  real(kind=8),parameter :: mtot = mC + mO + mN
  real(kind=8),parameter :: NA = 6.022140857D0 ! 6.022 140 857 x 10^23 mol-1
  real(kind=8),parameter :: pi4 = 1.256637061D1
  integer,parameter :: input = 11
end module
program main
  use global
  implicit none
  real(kind=8) :: mol(3,3) ! (1,:) : N ; (2,:) : C ; (3,:) : O
  real(kind=8) :: r1, r2, rm ! r1 : N-C ; r2 : C-O ; rm : mass center
  real(kind=8) :: MoInt, RotConst
  character(len=1) :: symbol
  real(kind=8) :: temp(3), temp1, temp2
  integer :: i
  
  open(unit=input,file="NCO.txt")
  do i = 1, 3
    read(input,*) symbol, temp(:)
    call ArrangeCoord(symbol,temp,mol)
  end do
  close(unit=input)
  
  call CalcBondLen(mol,r1,r2)
  write(*,"(' R_NC = ',F8.6,' Ang')") r1
  write(*,"(' R_CO = ',F8.6,' Ang')") r2
  temp1 = mN * r1 ; temp2 = mO * r2
  rm = (temp2 - temp1) / mtot
  write(*,"(' r_Nm = ',F8.6,' Ang')") r1 + rm
  write(*,"(' r_mO = ',F8.6,' Ang')") r2 - rm
  MoInt = temp1 * r1 + temp2 * r2 - (temp1 - temp2)*(temp1 - temp2)/mtot
  RotConst = rh * NA * 1.0D3 / (pi4 * MoInt) ! MHz
  write(*,"(' Rotational Constant = ',F10.4,' MHz')") RotConst
  
end program

subroutine ArrangeCoord(element,temp,mol)
  implicit none
  character(len=1),intent(in) :: element
  real(kind=8),intent(in) :: temp(3)
  real(kind=8) :: mol(3,3)
  
  select case (element)
    case('N')
      mol(1,:) = temp
    case('C')
      mol(2,:) = temp
    case('O')
      mol(3,:) = temp
    case default
      write(*,*) "Unknown Element"
      stop
  end select
  
end subroutine

subroutine CalcBondLen(mol,r1,r2)
  implicit none
  real(kind=8),intent(in) :: mol(3,3)
  real(kind=8),intent(out) :: r1, r2 ! r1 : N-C ; r2 : C-O
  real(kind=8) :: temp(3)
  
  temp = mol(1,:) - mol(2,:)
  r1 = sqrt( sum(temp * temp) )
  temp = mol(2,:) - mol(3,:)
  r2 = sqrt( sum(temp * temp) )
  
end subroutine