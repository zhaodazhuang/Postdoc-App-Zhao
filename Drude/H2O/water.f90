module basis_set
  implicit none
  integer, parameter :: nM = 64, nT = 256
  real(8), parameter :: lX = 12.4171267044140910d0
  real(8), parameter :: lY = 12.4171267044140910d0
  real(8), parameter :: lZ = 12.4171267044140910d0
  real(8), parameter :: mO = 15.9990d0, mH = 1.0080d0
  real(8), parameter :: cgO = -2.050d0, zO = 1.250d0, zH = 0.40d0
  real(8) :: invlX, invlY, invlZ

end module basis_set

module unit_conversion
  implicit none

  !!! 1.0 elementraY charge = 1.602176634d-19 coulomb
  real(8), parameter :: ele_to_clbp19 = 1.602176634d0 !!! C
  real(8), parameter :: epslonp12 = 8.8541878128d0 !!! C/(V*m)
  real(8),parameter:: nan23=6.02214076d0
  real(8), parameter :: kbp23 = 1.38064852d0 !!! Kg*m^2/s^2/k
end module unit_conversion

module short_range
  implicit none

  real(8), parameter :: Rcut = 6.0d0
  real(8), parameter :: Rcut2 = Rcut*Rcut
  real(8), parameter :: Kappa = 209.450d0

  real(8), parameter :: LJA = 39344.980d0, LJB = 42.150d0 !!! Lennard-Jones
  real(8), parameter :: pi = acos(-1.0d0)
  real(8), parameter :: hamA = 396.270d0, rho = 0.250d0, hamC = 10.0d0
  real(8), parameter :: invrho = 1.0d0/rho

  real(8), parameter :: Kijk = 4.199780d0
  real(8), parameter :: theta0 = 1.89699836399263670d0 ! 108.69du

  real(8), parameter :: MorD = 6.2037130d0, Moralpha = 2.220030d0
  real(8), parameter :: Morr0 = 0.92376
  
  real(8), parameter :: LJc3=  3.8815535414125965d-5
  real(8), parameter :: LJc4= -5.8557279375538297d-6
  real(8), parameter :: c3e=  7.6457767218140786d-8
  real(8), parameter :: c4e= -9.2385348286440336d-9
  real(8), parameter :: hamc3=  -1.9830846299845627d-5
  real(8), parameter :: hamc4=   2.2308524701510606d-6
 
end module short_range

module long_range
  implicit none

  real(8), parameter :: alpha = 0.60d0 !!! angstrom^-1
  real(8), parameter :: pi = acos(-1.0d0)
  integer, parameter :: kLmax = 6, kMmax = 6, kNmax = 6
  real(8), parameter :: alpha2 = alpha*alpha

end module long_range


program Water
  use basis_set
  use long_range
  implicit none

  real(8), parameter :: reqTemp = 300.0d0 ! K
  real(8), parameter :: dt = 2.0d-4
  real(8), parameter :: dt2 = dt*dt
  real(8), parameter :: invmO = 1.0d0/mO, invmH = 1.0d0/mH
  real(8) :: r4pie, eang, coitemp, coiv, betapeps
  
  real(8), dimension(1:nT) :: RX, RY, RZ
  real(8), dimension(nM+1:nT) :: VX, VY, VZ
  real(8), dimension(nM+1:nT) :: FRX, FRY, FRZ
  real(8), dimension(nM+1:nT) :: FKX, FKY, FKZ
  real(8), dimension(nM+1:nT) :: FXold, FYold, FZold
  real(8), dimension(nM+1:nT) :: FX, FY, FZ
  real(8), dimension(nM+1:nT) :: invm

  integer :: i
  integer :: step
  real(8) :: nowTemp
  real(8) :: t1, t2, t3

  invlX = 1.0d0/lX
  invlY = 1.0d0/lY
  invlZ = 1.0d0/lZ

  do i = nM+1, nT, 3
    
    invm(i) = 1.0d0/mO; invm(i+1) = 1.0d0/mH; invm(i+2) = 1.0d0/mH 
     
  end do

  call gen_para( reqTemp, r4pie, eang, coitemp, coiv, betapeps )

  call read_R( nT, nM, RX, RY, RZ )

  call read_V( nT, nM, VX, VY, VZ ) 

!  call gen_V( nM, nT, mO, mH, coiv, reqTemp, VX, VY, VZ )

!  open(32,file="temp.dat")
  open(33,file="Tr.xyz")

  !!! RUN !!!

!  call cpu_time(t1)

  call iterate_re( eang, r4pie, RX, RY, RZ )

  call Rspace_core( alpha, eang, r4pie, RX, RY, RZ, FRX, FRY, FRZ )
 
  call Kspace_core( r4pie, RX, RY, RZ, FKX, FKY, FKZ )

  FXold = FRX + FKX; FYold = FRY + FKY; FZold = FRZ + FKZ 

  do step = 1, 50000

    RX(nM+1:nT) = RX(nM+1:nT) + VX*dt + 0.50d0*FXold*invm*dt2
    RY(nM+1:nT) = RY(nM+1:nT) + VY*dt + 0.50d0*FYold*invm*dt2
    RZ(nM+1:nT) = RZ(nM+1:nT) + VZ*dt + 0.50d0*FZold*invm*dt2

    call iterate_re( eang, r4pie, RX, RY, RZ )

    if ( step > 10000 .and. mod(step,5) == 0 ) then

    write(33,*) nT
    write(33,*) lX, lY, lZ
    do i = 1, nM

      write(33,*) "W", RX(i), RY(i), RZ(i)
 
    end do
    do i = nM+1, nT, 3

      write(33,*) "O", RX(i), RY(i), RZ(i)
      write(33,*) "H", RX(i+1), RY(i+1), RZ(i+1)
      write(33,*) "H", RX(i+2), RY(i+2), RZ(i+2)
    end do

    end if

    call Rspace_core( alpha, eang, r4pie, RX, RY, RZ, FRX, FRY, FRZ )

    call Kspace_core( r4pie, RX, RY, RZ, FKX, FKY, FKZ )

    FX = FRX + FKX; FY = FRY + FKY; FZ = FRZ + FKZ

    VX = VX + 0.50d0*(FXold+FX)*invm*dt
    VY = VY + 0.50d0*(FYold+FY)*invm*dt
    VZ = VZ + 0.50d0*(FZold+FZ)*invm*dt

!    if ( step <= 20000 ) then

!      call contral_temp( nM, nT, VX, VY, VZ, invm, coitemp, reqTemp)

!    else if ( step > 20000 .and. mod(step,10) == 0 ) then
!    if ( mod(step,100) == 0 ) then

!      call calcul_temp( nM, nT, VX, VY, VZ, invm, coitemp, reqTemp, nowTemp )                          
!      write(32,*) step, nowTemp

!    end if

    FXold = FX; FYold = FY; FZold = FZ

  end do !!! loop for step

!  call cpu_time(t2)

!  write(*,*) "TOT TIME", (t2-t1)/60.0d0, "min"

!  close(32)
  close(33)

  call output_R_V( RX, RY, RZ, VX, VY, VZ )

end program Water



subroutine gen_para( reqTemp, r4pie, eang, coitemp, coiv, betapeps )
  use unit_conversion
  implicit none

  real(8), parameter :: pi = acos(-1.0d0)
  real(8) :: reqTemp
  real(8) :: r4pie, eang, coitemp, coiv, betapeps

  r4pie = ele_to_clbp19**2*nan23/(4.0d0*pi*epslonp12)*1.0d6 !!! unit: 10J/mol/ang
  eang = ele_to_clbp19*nan23*1.0d3 !!! unit: 10J/mol/ang
  coitemp = 10.0d0/(3.0d0*nan23*kbp23) !!! unit: K
  coiv = nan23*kbp23*1.0d-1 !!! unit: ans^2/ps^2
  betapeps = kbp23*reqTemp*epslonp12*1.0d-7/ele_to_clbp19**2 !!! unit: e^2/ang

return

end subroutine gen_para

subroutine read_R( nT, nM, RX, RY, RZ )
  implicit none

  integer :: nT, nM
  real(8), dimension(nT) ::  RX, RY, RZ

  integer :: i
  character(1) :: a
  
  open(33,file="/DATA/users/yihao/shell_module/H2O/L12/data4/100.xyz")
   
    read(33,*) 
    read(33,*)

    do i = 1, nT
 
      read(33,*) a, RX(i), RY(i), RZ(i)

    end do

  close(33) 

  
return 

end subroutine read_R

subroutine read_V( nT, nM, VX, VY, VZ )
  implicit none

  integer :: nT, nM
  real(8), dimension(nM+1:nT) ::  VX, VY, VZ

  integer :: i
  character(1) :: a

  open(33,file="/DATA/users/yihao/shell_module/H2O/L12/data4/100.v")

    read(33,*)
    read(33,*)

    do i = nM+1, nT

      read(33,*) a, VX(i), VY(i), VZ(i)

    end do

  close(33)


return

end subroutine read_V
 

subroutine gen_V( nM, nT, mO, mH, coiv, reqTemp, VX, VY, VZ )
  implicit none

  integer :: nM, nT
  real(8) :: mO, mH, coiv, reqTemp
  real(8), dimension(nM+1:nT) :: VX, VY, VZ

  integer :: i
  real(8) :: pi, twopi
  real(8), dimension(nM+1:nT) :: x, y, z

  pi = acos(-1.0d0)
  twopi = 2.0d0*pi

  call random_seed()
  call random_number(x)
  call random_number(y)
  z = sqrt( -2.0d0*log(x) )*( cos(twopi*y) )

  do i = nM+1, nT, 3
    
    VX(i) = sqrt(coiv*reqtemp/mO)*z(i)
    VX(i+1) = sqrt(coiv*reqtemp/mH)*z(i+1)
    VX(i+2) = sqrt(coiv*reqtemp/mH)*z(i+2)

  end do

  call random_seed()
  call random_number(x)
  call random_number(y)
  z = sqrt( -2.0d0*log(x) )*( cos(twopi*y) )

  do i = nM+1, nT, 3

    VY(i) = sqrt(coiv*reqtemp/mO)*z(i)
    VY(i+1) = sqrt(coiv*reqtemp/mH)*z(i+1)
    VY(i+2) = sqrt(coiv*reqtemp/mH)*z(i+2)

  end do

  call random_seed()
  call random_number(x)
  call random_number(y)
  z = sqrt( -2.0d0*log(x) )*( cos(twopi*y) )

  do i = nM+1, nT, 3

    VZ(i) = sqrt(coiv*reqtemp/mO)*z(i)
    VZ(i+1) = sqrt(coiv*reqtemp/mH)*z(i+1)
    VZ(i+2) = sqrt(coiv*reqtemp/mH)*z(i+2)

  end do

return

end subroutine gen_V

subroutine iterate_re( eang, r4pie, RX, RY, RZ ) 
  use basis_set
  use long_range
  implicit none

  real(8) :: eang, r4pie
  real(8), dimension(1:nT) :: RX, RY, RZ

  integer :: j, loop
  real(8), dimension(1:nM) :: FRX, FRY, FRZ
  real(8), dimension(1:nM) :: FKX, FKY, FKZ
  real(8) :: diffF, enlarge, narrow
  real(8), dimension(1:nM) :: etaX, etaY, etaZ
  real(8), dimension(1:nM) :: FX, FY, FZ, FFX, FFY, FFZ

  diffF = 1.0d0
  etaX = 1.40d-7
  etaY = 1.40d-7
  etaZ = 1.40d-7

  enlarge = 1.05d0
  narrow  = 0.95d0

  loop = 0

  do while( diffF .gt. 5.0d-6 )

    call Rspace_shell( alpha, eang, r4pie, RX, RY, RZ, FRX, FRY, FRZ)

    call Kspace_shell( r4pie, RX, RY, RZ, FKX, FKY, FKZ)

    FX = FRX+FKX
    FY = FRY+FKY
    FZ = FRZ+FKZ

    do j = 1, nM

      if ( FFX(j)*FX(j) .gt. 0 ) then

        etaX(j) = etaX(j)*enlarge

      else

        etaX(j) = etaX(j)*narrow

      end if

      if ( FFY(j)*FY(j) .gt. 0 ) then

        etaY(j) = etaY(j)*enlarge

      else

        etaY(j) = etaY(j)*narrow

      end if

      if ( FFZ(j)*FZ(j) .gt. 0 ) then

        etaZ(j) = etaZ(j)*enlarge

      else

        etaZ(j) = etaZ(j)*narrow

      end if

    end do
  
    RX(1:nM) = RX(1:nM) + FX*etaX
    RY(1:nM) = RY(1:nM) + FY*etaY
    RZ(1:nM) = RZ(1:nM) + FZ*etaZ

    FFX = FX
    FFY = FY
    FFZ = FZ

    diffF =  sum( abs(FX(1:nM)) ) + sum( abs(FY(1:nM)) ) + sum( abs(FZ(1:nM)) )

    loop = loop + 1

    if ( mod(loop, 50) .eq. 0 ) then

      write(*,*) diffF

    end if

  end do

return

end subroutine iterate_re 

subroutine contral_temp( nM, nT, VX, VY, VZ, invm, coitemp, reqTemp)
  implicit none
  
  integer :: nM, nT
  real(8), dimension(nM+1:nT) :: VX, VY, VZ       
  real(8), dimension(nM+1:nT) :: invm
  real(8) :: coitemp, reqTemp

  real(8) :: nowTemp, beta

  nowtemp = sum( invm*(VX*VX + VY*VY +VZ*VZ ) ) 

  nowTemp = nowTemp*coitemp/real(nT-nM ,8)

  beta = sqrt( reqTemp/nowTemp )

  VX = VX*beta
  VY = VY*beta
  VZ = VZ*beta

return

end subroutine contral_temp

subroutine calcul_temp( nM, nT, VX, VY, VZ, invm, coitemp, reqTemp, nowTemp)
  implicit none

  integer :: nM, nT
  real(8), dimension(nM+1:nT) :: VX, VY, VZ
  real(8), dimension(nM+1:nT) :: invm
  real(8) :: coitemp, reqTemp

  real(8) :: nowTemp

  nowtemp = sum( invm*(VX*VX + VY*VY +VZ*VZ ) )

  nowTemp = nowTemp*coitemp/real(nT-nM ,8)

return

end subroutine calcul_temp


subroutine output_R_V( RX, RY, RZ, VX, VY, VZ )
  use basis_set
  implicit none
   
  real(8), dimension(1:nT) :: RX, RY, RZ
  real(8), dimension(nM+1:nT) :: VX, VY, VZ

  integer :: i

  open(34, file="coor.xyz")

    write(34,*) nT
    write(34,*) lX, lY, lZ
    do i = 1, nM

      write(34,*) "W", RX(i), RY(i), RZ(i)

    end do
    do i = nM+1, nT, 3

      write(34,*) "O", RX(i), RY(i), RZ(i)
      write(34,*) "H", RX(i+1), RY(i+1), RZ(i+1)
      write(34,*) "H", RX(i+2), RY(i+2), RZ(i+2)
    end do

  close(34)  

  open(35, file="V.dat")
 
    write(35,*) nT-nM
    write(35,*) lX, lY, lZ

    do i = nM+1, nT, 3

      write(35,*) "O", VX(i), VY(i), VZ(i)
      write(35,*) "H", VX(i+1), VY(i+1), VZ(i+1)
      write(35,*) "H", VX(i+2), VY(i+2), VZ(i+2)
    end do

  close(35)

return

end subroutine output_R_V

