module basis_set
  implicit none

  integer, parameter :: nion = 216
  real(8), parameter :: lX = 24.0959889334770030d0, lY=24.0959889334770030d0, lZ=24.0959889334770030d0
  real(8), dimension(1:2*nion) :: rX, rY, rZ
  real(8) :: invlX, invlY, invlZ
  integer :: hnion, twonion
end module basis_set  

module unit_conversion
  implicit none

  !!! 1.0 elementraY charge = 1.602176634d-19 coulomb
  real(8), parameter :: ele_to_clbp19 = 1.602176634d0 !!! C
  real(8), parameter :: epslonp12 = 8.8541878128d0 !!! C/(V*m)
  real(8),parameter:: nan23=6.02214076d0
  real(8), parameter :: kbp23 = 1.38064852 !!! Kg*m^2/s^2/k
end module unit_conversion

module short_range
  implicit none

  real(8), parameter :: cutoff = 10.0d0
  real(8), parameter :: kNa = 63.0140d0, kCl = 25.7240d0
  real(8), parameter :: ANaNa = 487.0d0, ANaCl = 145134.0d0, AClCl = 405774.0d0
  real(8), parameter :: rho = 0.237680d0
  real(8), parameter :: CNaNa = 1.050d0, CNaCl = 6.990d0, CClCl = 72.40d0
  real(8), parameter :: DNaNa = 0.50d0, DNaCl = 8.70d0, DClCl = 145.40d0 
  real(8) :: c3NaNa, c3NaCl, c3ClCl, c4NaNa, c4NaCl, c4ClCl
  real(8) :: eang, invrho, cutoff2, c33, c44

end module short_range

module Coulombic
  implicit none
  
  real(8), parameter :: cgNa = -0.50560d0, cgCl = -2.50050d0
  real(8) :: zNa, zCl
  real(8), parameter :: alpha = 1.0d0/4.50d0 !!! angstrom^-1
  real(8), parameter :: pi = acos(-1.0d0)  
  integer, parameter :: kLmax = 8, kMmax = 8, kNmax = 8
  real(8) :: c3e, c4e
  real(8) :: r4pie, alpha2

end module Coulombic

module update
  use basis_set
  implicit none

  real(8), parameter :: mNa = 22.989769d0, mCl = 35.4530d0
  real(8), parameter :: dt = 2.0d-4
  real(8), dimension(nion+1:2*nion) :: vX, vY, vZ
  real(8), dimension(nion+1:2*nion) :: hvX, hvY, hvZ
  real(8) :: coiv

end module update

module RDfmod
 implicit none

 integer :: maxbin
 integer, parameter :: kN = 2000
 real(8), parameter :: dx = 1.0d-2
 integer, allocatable :: histNaNa(:), histNaCl(:), histClCl(:)
 real(8), dimension(1:kN) :: Sk 

end module RDfmod
  
program main_NaCl
  use basis_set
  use unit_conversion
  use short_range
  use Coulombic
  use update
  use RDfmod
  implicit none

  integer :: step, nrecord
  real(8) :: invmNa, invmCl, hrhoB
  real(8) :: coitemp, nowtemp, reqtemp, betapeps
  real(8), dimension(nion+1:2*nion) :: FRX, FRY, FRZ
  real(8), dimension(nion+1:2*nion) :: FKX, FKY, FKZ
  real(8), dimension(nion+1:2*nion) :: FX, FY, FZ
  real(8) :: t1, t2

  reqtemp = 1427.0d0

  zNa = 1.0d0-cgNa
  zCl = -1.0d0-cgCl
  invrho = 1.0d0/rho

  nrecord = 0
  maxbin = aint( cutoff/dx + dx/100.0d0 ) 
  allocate( histNaNa(0:maxbin) ) ; histNaNa = 0
  allocate( histNaCl(0:maxbin) ) ; histNaCl = 0
  allocate( histClCl(0:maxbin) ) ; histClCl = 0  
                                   Sk = 0.0d0 

  invlX = 1.0d0/lX
  invlY = 1.0d0/lY
  invlZ = 1.0d0/lZ
 
  hnion = nion/2
  twonion = 2*nion 

  hrhoB = real(hnion,8)*invlX*invlY*invlZ
 
  r4pie = ele_to_clbp19**2*nan23/(4.0d0*pi*epslonp12)*1.0d6 !!! unit: 10J/mol/ang 
  eang = ele_to_clbp19*nan23*1.0d3 !!! unit: 10J/mol/ang
  coitemp = 10.0d0/(3.0d0*nan23*kbp23) !!! unit: K
  coiv = nan23*kbp23*1.0d-1 !!! unit: ans^2/ps^2  
  betapeps = kbp23*reqtemp*epslonp12*1.0d-7/ele_to_clbp19**2 !!! unit: e^2/ang

  !!! now cutoff = 10.0d0, alpha = 1.0d0/4.50d0
  c3NaNa = -2.1159999835524911d-8 
  c4NaNa = 1.4284999879380471d-9
  c3NaCl = -1.4258395098372078d-7
  c4NaCl = 9.6278964053398648d-9
  c3ClCl = -1.4945278629574622d-6
  c4ClCl = 1.0093878994984208d-7
  
  c3e = 9.1906709115610602d-6
  c4e = -6.4018719925921016d-7

  cutoff2 = cutoff*cutoff
  alpha2 = alpha*alpha

  invmNa = 1.0d0/mNa
  invmCl = 1.0d0/mCl

  call initialize_pos(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ)

  call gauss_velocity(nion, hnion, twonion, hvX, hvY, hvZ, mNa, mCl, coiv, reqtemp)

!  call read_v_coor_hv(nion, twonion, rX, rY, rZ, vX, vY, vZ, hvX, hvY, hvZ)
 
  call cpu_time(t1)

  do step = 1, 500

    call iterate_re(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ, &
       & invlX, invlY, invlZ)

    call Rspace_core(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ, &
       & invlX, invlY, invlZ, FRX, FRY, FRZ, kNa, kCl, cutoff2, eang)

    call Kspace_core(nion, hnion, twonion, rX, rY, rZ, &
       & invlX, invlY, invlZ, FKX, FKY, FKZ) 

    FX = FRX + FKX ; FY = FRY + FKY ; FZ = FRZ + FKZ
    
   !if ( mod(step,1) .eq. 0 ) then

     call contral_temp(nion, hnion, twonion, hvX, hvY, hvZ, mNa, mCl,&
                     & coitemp, nowtemp, reqtemp)

      write(*,*) step, nowtemp

     call record_hist(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ,&
                    & invlX, invlY, invlZ, cutoff2 )

     call record_Sk(nion, hnion, twonion, kN, invlZ,&
        & pi, cgNa, cgCl, zNa, zCl, rZ, Sk)

     nrecord = nrecord + 1

   !end if

    call update_rc(nion, hnion, twonion, invmNa, invmCl, rX, rY, rZ, &
       & lX, lY, lZ, invlX, invlY, invlZ, FX, FY, FZ, vX, vY, vZ, &
       & hvX, hvY, hvZ, dt)

  end do 
  
  call cpu_time(t2)

  write(*,*) (t2-t1)/60.0d0, 'min'

  call outputcoor(twonion, rX, rY, rZ)
  call outputv(nion, twonion, vX, vY, vZ, hvX, hvY, hvZ)

  call output_RDF(hnion, nrecord, pi, hrhoB)
  call output_Sk(hnion, kN, nrecord, invlX, invlY, invlZ,&
     & cgNa, cgCl, zNa, zCl, pi, Sk, betapeps)


  call output_finalcoor(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ)


end


subroutine iterate_re(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ, &
         & invlX, invlY, invlZ)
  use Coulombic
  implicit none

  integer :: nion, hnion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8) :: lX, lY, lZ, invlX, invlY, invlZ

  integer :: j, loop
  real(8), dimension(1:nion) :: FRX, FRY, FRZ
  real(8), dimension(1:nion) :: FKX, FKY, FKZ
  real(8) :: diffF, enlarge, narrow
  real(8), dimension(1:nion) :: etaX, etaY, etaZ
  real(8), dimension(1:nion) :: FX, FY, FZ, FFX, FFY, FFZ

  diffF = 1.0d0
  etaX = 1.40d-6
  etaY = 1.40d-6
  etaZ = 1.40d-6

  enlarge = 1.05d0
  narrow  = 0.95d0

  loop = 0

  do while( diffF .gt. 2.0d-5 )

    call Rspace_shell(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ, &
       & invlX, invlY, invlZ, FRX, FRY, FRZ)

    call Kspace_shell(nion, hnion, twonion, rX, rY, rZ, &
       & invlX, invlY, invlZ, FKX, FKY, FKZ)

    do j = 1, nion

      FX(j) = FRX(j)+FKX(j)
      FY(j) = FRY(j)+FKY(j)
      FZ(j) = FRZ(j)+FKZ(j)

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

      rX(j) = rX(j) + FX(j)*etaX(j)
      rY(j) = rY(j) + FY(j)*etaY(j)
      rZ(j) = rZ(j) + FZ(j)*etaZ(j)

      rX(j) = rX(j) - floor(rX(j)*invlX)*lX
      rY(j) = rY(j) - floor(rY(j)*invlY)*lY
      rZ(j) = rZ(j) - floor(rZ(j)*invlZ)*lZ

      FFX(j) = FX(j)
      FFY(j) = FY(j)
      FFZ(j) = FZ(j)

    end do

      diffF =  sum( abs(FX(1:nion)) ) + sum( abs(FY(1:nion)) ) + sum( abs(FZ(1:nion)) )

      loop = loop + 1

    if ( mod(loop, 500) .eq. 0 ) then

      write(*,*) diffF

    end if

  end do

return

end subroutine iterate_re




