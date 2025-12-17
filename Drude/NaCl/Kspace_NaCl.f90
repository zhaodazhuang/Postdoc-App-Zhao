subroutine Kspace_shell(nion, hnion, twonion, rX, rY, rZ, &
         & invlX, invlY, invlZ, FKX, FKY, FKZ)
  use Coulombic
  implicit none

  integer :: nion, hnion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8) :: invlX, invlY, invlZ
  real(8), dimension(1:nion) :: FKX, FKY, FKZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! 2*pi         e^{-k^2/(4*alpha^2)}
!!! ———— sum_{k} ———————————————————— sum_{i} q_i*e^{i*k*r_i} sum_{j} q_j*e^{-i*k*r_j}
!!!  V                  k^2
!!!                      |                   |                   |
!!!                      |                   |                   |
!!!                      |     sumqcos(k*r)+i*sumqsin(k*r)       |
!!!                     akk                                      |
!!!                                               sumqcos(k*r)-i*sumqsin(k*r)
!!!
!!!                                 ( sumqcos(kr) )^2+ ( sumqsin(kr) )^2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! k is a three-dimensional vector k(rkL, rkM, rkN)
!!!
!!!       2*pi*kL
!!! rkL = ———————  ...   ...      kL, kM, kN is integer
!!!         LX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            4*pi         \vec{k}
!!! F_j = -q_j ———— sum_{k} ——————— e^{-k^2/(4*alpha^2)} sum_{i} q_i*sin[k*(r_i-r_j)]
!!!             V             k^2
!!!
!!!           4*pi         \vec{k}
!!!     = q_j ———— sum_{k} ——————— e^{-k^2/(4*alpha^2)} sum_{i} q_i*sin[k*(r_j-r_i)]
!!!            V             k^2
!!!           q_j*sum_{i} q_i*sin[k*(r_j-r_i)]
!!!                       |
!!!                       |
!!!  q_j*sin(k*r_j)*sumqcos(k*r) - q_j*cos(k*r_j)*sumqsin(k*r)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Why use only kx symmetry then X2 instead of kx, ky, kz sysmetry then X8?
!!!
!!! https://physics.stackexchange.com/questions/205794/symmetry-in-program-for-ewald-summation
!!! And to add to this answer: sin^2(x) is also an even function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: j
  integer :: kL, kM, kN, LL, MM, NN, kMmin, kNmin
  real(8) :: ssX, ssY, ssZ, rkX, rkY, rkZ, rksq, invrksq
  real(8) :: p4alp2, twopi, akk, sumqcoskr, sumqsinkr, cs
  real(8) :: F, Fcoeffi, rcpcut, rcpct2
  real(8), dimension(1:twonion,0:kLmax) :: eLc, eLs
  real(8), dimension(1:twonion,0:kMmax) :: eMc, eMs
  real(8), dimension(1:twonion,0:kNmax) :: eNc, eNs 
  real(8), dimension(1:twonion) :: clm, slm, coskr, sinkr 
   
  p4alp2 = -0.250d0/alpha2 
  twopi = 2.0d0*pi
  
  do j = 1, nion 

    FKX(j) = 0.0d0
    FKY(j) = 0.0d0
    FKZ(j) = 0.0d0

  end do

  ! set up reciprocal space cutoff, note: rcpcut is twice of rcpcut in pme/gse
  rcpcut = min(real(kLmax,8)*invlX,real(kMmax,8)*invlY,real(kNmax,8)*invlZ)
  rcpcut = rcpcut*1.050d0*twopi
  rcpct2 = rcpcut**2 

  do j = 1, twonion
 
    eLc(j,0) = 1.0d0;  eLs(j,0) = 0.0d0
    eMc(j,0) = 1.0d0;  eMs(j,0) = 0.0d0
    eNc(j,0) = 1.0d0;  eNs(j,0) = 0.0d0

    ! e^{i*k*r} = e^{i*(rkL*rX + rkM*rY + rkL*rz)}
    !           = e^{i*rkL*rX}*e^{...}*e^{...}
    !           = e^{i*2*pi*kL*rX/Lx}*...*...
    ssX = rX(j)*invlX
    ssY = rY(j)*invlY
    ssZ = rZ(j)*invlZ

    eLc(j,1) = cos(twopi*ssX); eLs(j,1) = sin(twopi*ssX)
    eMc(j,1) = cos(twopi*ssY); eMs(j,1) = sin(twopi*ssY)
    eNc(j,1) = cos(twopi*ssZ); eNs(j,1) = sin(twopi*ssZ)

 end do

 do kM = 2, kMMax

   do j = 1, twonion
     ! cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
     eMc(j,kM) = eMc(j,kM-1)*eMc(j,1) - eMs(j,kM-1)*eMs(j,1)
     ! sin(a+b) = sin(a)cos(b) + cos(a)sin(b)
     eMs(j,kM) = eMs(j,kM-1)*eMc(j,1) + eMc(j,kM-1)*eMs(j,1)

   end do

 end do

 do kN = 2, kNMax

   do j = 1, twonion

     eNc(j,kN) = eNc(j,kN-1)*eNc(j,1) - eNs(j,kN-1)*eNs(j,1)

     eNs(j,kN) = eNs(j,kN-1)*eNc(j,1) + eNc(j,kN-1)*eNs(j,1)

   end do

 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 kMmin = 0
 kNmin = 1

 do LL = 0, kLmax
 
   kL = LL
   rkX = twopi*real(kL,8)*invlX 
  ! put eLc(j, kL) to eLc(j,0) 

   if ( kL .eq. 1 ) then

     do j = 1, twonion

       eLc(j,0) = eLc(j,1); eLs(j,0) = eLs(j,1)

     end do

   else if ( kL .gt. 1) then

     do j = 1, twonion
       ! now, the stored value of eLc(j,0) is eLc(j,kL-1)
       cs = eLc(j,0)
       eLc(j,0) = cs      *eLc(j,1) - eLs(j,0)*eLs(j,1)
       eLs(j,0) = eLs(j,0)*eLc(j,1) + cs      *eLs(j,1)
       ! now, the stored value of eLc(j,0) is eLc(j,kL)
     end do

   end if

   do MM = kMmin, kMmax

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     kM = iabs(MM)
     rkY = twopi*real(MM,8)*invlY

     if ( MM .ge. 0) then

       ! e^[i ( rkL*rX + rkM*rY )]  !!! e^ikr = e^[i ( rkL*rX + rkM*rY + rkN*rZ )]
       do j = 1, twonion

         clm(j) = eLc(j,0)*eMc(j,kM) - eLs(j,0)*eMs(j,kM)
         slm(j) = eLs(j,0)*eMc(j,kM) + eLc(j,0)*eMs(j,kM)

       end do

     else
       ! now MM is negative but kM is positive so cos(a-b)
       do j = 1, twonion

         clm(j) = eLc(j,0)*eMc(j,kM) + eLs(j,0)*eMs(j,kM)
         slm(j) = eLs(j,0)*eMc(j,kM) - eLc(j,0)*eMs(j,kM)

       end do

     end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do NN = kNmin, kNmax

         kN = iabs(NN)
        rkZ = twopi*real(NN,8)*invlZ
       rksq = rkX*rkX + rkY*rkY + rkZ*rkZ
 
       if ( rksq .le. rcpct2 ) then
  
         if ( NN .ge. 0 ) then

           do j = 1, twonion

             coskr(j) = clm(j)*enc(j,kN) - slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) + clm(j)*ens(j,kN)

           end do

         else

           do j = 1, twonion

             coskr(j) = clm(j)*enc(j,kN) + slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) - clm(j)*ens(j,kN)

           end do

         end if

         ! sum_{i} q_i*e^{i*k*r_i}
         sumqcoskr = 0.0d0
         sumqsinkr = 0.0d0
 
         do j = 1, hnion

           sumqcoskr = sumqcoskr +cgNa*coskr(j)
           sumqsinkr = sumqsinkr +cgNa*sinkr(j)

         end do

         do j = hnion+1, nion

           sumqcoskr = sumqcoskr + cgCl*coskr(j)
           sumqsinkr = sumqsinkr + cgCl*sinkr(j)

         end do

         do j = nion+1, nion+hnion

           sumqcoskr = sumqcoskr + zNa*coskr(j)
           sumqsinkr = sumqsinkr + zNa*sinkr(j)

         end do

         do j = nion+hnion+1, twonion

           sumqcoskr = sumqcoskr + zCl*coskr(j)
           sumqsinkr = sumqsinkr + zCl*sinkr(j)

         end do
         
         invrksq = 1.0d0/rksq
         akk = exp(p4alp2*rksq)*invrksq

         !!! calculate force between shell and shell 

         do j = 1, nion

           if ( j .le. hnion ) then
             !!! there is no potential between core and own shell 
             F = akk*cgNa*( sinkr(j)*( sumqcoskr-zNa*coskr(j+nion) ) & 
                        & - coskr(j)*( sumqsinkr-zNa*sinkr(j+nion) ) )
 
           else if ( j .ge. hnion + 1) then

             F = akk*cgCl*( sinkr(j)*( sumqcoskr-zCl*coskr(j+nion) ) & 
                        & - coskr(j)*( sumqsinkr-zCl*sinkr(j+nion) ) )
            
           end if

           FKX(j) = FKX(j) + rkX*F
           FKY(j) = FKY(j) + rkY*F
           FKZ(j) = FKZ(j) + rkZ*F

         end do

       end if
  
     end do
 
     kNmin = -kNmax

   end do

   kMmin = -kMmax
   
 end do

 Fcoeffi = twopi*r4pie*invlX*invlY*invlZ

 Fcoeffi = 4.0d0*Fcoeffi 

 do j = 1, nion
   
   FKX(j) = FKX(j)*Fcoeffi
   FKY(j) = FKY(j)*Fcoeffi
   FKZ(j) = FKZ(j)*Fcoeffi
   
 end do
  
return

end subroutine Kspace_shell


























subroutine Kspace_core(nion, hnion, twonion, rX, rY, rZ, &
         & invlX, invlY, invlZ, FKX, FKY, FKZ)
  use Coulombic
  implicit none

  integer :: nion, hnion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8) :: invlX, invlY, invlZ
  real(8), dimension(nion+1:twonion) :: FKX, FKY, FKZ

  integer :: j
  integer :: kL, kM, kN, LL, MM, NN, kMmin, kNmin
  real(8) :: ssX, ssY, ssZ, rkX, rkY, rkZ, rksq, invrksq
  real(8) :: p4alp2, twopi, akk, sumqcoskr, sumqsinkr, cs
  real(8) :: F, Fcoeffi, rcpcut, rcpct2
  real(8), dimension(1:twonion,0:kLmax) :: eLc, eLs
  real(8), dimension(1:twonion,0:kMmax) :: eMc, eMs
  real(8), dimension(1:twonion,0:kNmax) :: eNc, eNs
  real(8), dimension(1:twonion) :: clm, slm, coskr, sinkr

  p4alp2 = -0.250d0/alpha2
  twopi = 2.0d0*pi

  do j = nion+1, twonion

    FKX(j) = 0.0d0
    FKY(j) = 0.0d0
    FKZ(j) = 0.0d0

  end do

  ! set up reciprocal space cutoff, note: rcpcut is twice of rcpcut in pme/gse
  rcpcut = min(real(kLmax,8)*invlX,real(kMmax,8)*invlY,real(kNmax,8)*invlZ)
  rcpcut = rcpcut*1.050d0*twopi
  rcpct2 = rcpcut**2

  do j = 1, twonion

    eLc(j,0) = 1.0d0;  eLs(j,0) = 0.0d0
    eMc(j,0) = 1.0d0;  eMs(j,0) = 0.0d0
    eNc(j,0) = 1.0d0;  eNs(j,0) = 0.0d0

    ssX = rX(j)*invlX
    ssY = rY(j)*invlY
    ssZ = rZ(j)*invlZ

    eLc(j,1) = cos(twopi*ssX); eLs(j,1) = sin(twopi*ssX)
    eMc(j,1) = cos(twopi*ssY); eMs(j,1) = sin(twopi*ssY)
    eNc(j,1) = cos(twopi*ssZ); eNs(j,1) = sin(twopi*ssZ)

 end do

 do kM = 2, kMMax

   do j = 1, twonion
     ! cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
     eMc(j,kM) = eMc(j,kM-1)*eMc(j,1) - eMs(j,kM-1)*eMs(j,1)
     ! sin(a+b) = sin(a)cos(b) + cos(a)sin(b)
     eMs(j,kM) = eMs(j,kM-1)*eMc(j,1) + eMc(j,kM-1)*eMs(j,1)

   end do

 end do

 do kN = 2, kNMax

   do j = 1, twonion

     eNc(j,kN) = eNc(j,kN-1)*eNc(j,1) - eNs(j,kN-1)*eNs(j,1)

     eNs(j,kN) = eNs(j,kN-1)*eNc(j,1) + eNc(j,kN-1)*eNs(j,1)

   end do

 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 kMmin = 0
 kNmin = 1

 do LL = 0, kLmax

   kL = LL
   rkX = twopi*real(kL,8)*invlX
  ! put eLc(j, kL) to eLc(j,0)

   if ( kL .eq. 1 ) then

     do j = 1, twonion

       eLc(j,0) = eLc(j,1); eLs(j,0) = eLs(j,1)

     end do

   else if ( kL .gt. 1) then

     do j = 1, twonion
       ! now, the stored value of eLc(j,0) is eLc(j,kL-1)
       cs = eLc(j,0)
       eLc(j,0) = cs      *eLc(j,1) - eLs(j,0)*eLs(j,1)
       eLs(j,0) = eLs(j,0)*eLc(j,1) + cs      *eLs(j,1)
       ! now, the stored value of eLc(j,0) is eLc(j,kL)
     end do

   end if

   do MM = kMmin, kMmax

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     kM = iabs(MM)
     rkY = twopi*real(MM,8)*invlY

     if ( MM .ge. 0) then

       ! e^[i ( rkL*rX + rkM*rY )]  !!! e^ikr = e^[i ( rkL*rX + rkM*rY + rkN*rZ )]
       do j = 1, twonion

         clm(j) = eLc(j,0)*eMc(j,kM) - eLs(j,0)*eMs(j,kM)
         slm(j) = eLs(j,0)*eMc(j,kM) + eLc(j,0)*eMs(j,kM)

       end do

     else
       ! now MM is negative but kM is positive so cos(a-b)
       do j = 1, twonion

         clm(j) = eLc(j,0)*eMc(j,kM) + eLs(j,0)*eMs(j,kM)
         slm(j) = eLs(j,0)*eMc(j,kM) - eLc(j,0)*eMs(j,kM)

       end do

     end if

     do NN = kNmin, kNmax

         kN = iabs(NN)
        rkZ = twopi*real(NN,8)*invlZ
       rksq = rkX*rkX + rkY*rkY + rkZ*rkZ

       if ( rksq .le. rcpct2 ) then

         if ( NN .ge. 0 ) then

           do j = 1, twonion

             coskr(j) = clm(j)*enc(j,kN) - slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) + clm(j)*ens(j,kN)

           end do

         else

           do j = 1, twonion

             coskr(j) = clm(j)*enc(j,kN) + slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) - clm(j)*ens(j,kN)

           end do

         end if

         ! sum_{i} q_i*e^{i*k*r_i}
         sumqcoskr = 0.0d0
         sumqsinkr = 0.0d0

         do j = 1, hnion

           sumqcoskr = sumqcoskr +cgNa*coskr(j)
           sumqsinkr = sumqsinkr +cgNa*sinkr(j)

         end do

         do j = hnion+1, nion

           sumqcoskr = sumqcoskr + cgCl*coskr(j)
           sumqsinkr = sumqsinkr + cgCl*sinkr(j)

         end do

         do j = nion+1, nion+hnion

           sumqcoskr = sumqcoskr + zNa*coskr(j)
           sumqsinkr = sumqsinkr + zNa*sinkr(j)

         end do

         do j = nion+hnion+1, twonion

           sumqcoskr = sumqcoskr + zCl*coskr(j)
           sumqsinkr = sumqsinkr + zCl*sinkr(j)

         end do

         invrksq = 1.0d0/rksq
         akk = exp(p4alp2*rksq)*invrksq

         !!! calculate force between shell and shell

         do j = nion+1, twonion

           if ( j .le. nion+hnion ) then

             F = akk*zNa*( sinkr(j)*( sumqcoskr - cgNa*coskr(j-nion) ) &
                       & - coskr(j)*( sumqsinkr - cgNa*sinkr(j-nion) ) ) 

           else if ( j .ge. nion+hnion + 1) then

             F = akk*zCl*( sinkr(j)*( sumqcoskr - cgCl*coskr(j-nion) ) & 
                       & - coskr(j)*( sumqsinkr - cgCl*sinkr(j-nion) ) )

           end if

           FKX(j) = FKX(j) + rkX*F
           FKY(j) = FKY(j) + rkY*F
           FKZ(j) = FKZ(j) + rkZ*F

         end do

       end if

     end do

     kNmin = -kNmax

   end do

   kMmin = -kMmax

 end do

 Fcoeffi = twopi*r4pie*invlX*invlY*invlZ

 Fcoeffi = 4.0d0*Fcoeffi

 do j = nion+1, twonion

   FKX(j) = FKX(j)*Fcoeffi
   FKY(j) = FKY(j)*Fcoeffi
   FKZ(j) = FKZ(j)*Fcoeffi

 end do

return

end subroutine Kspace_core

subroutine record_Sk(nion, hnion, twonion, kN, invlZ,&
         & pi, cgNa, cgCl, zNa, zCl, rZ, Sk)
  implicit none

  integer :: nion, hnion, twonion, kN
  real(8) :: invlZ, pi, cgNa, cgCl, zNa, zCl
  real(8), dimension(1:twonion) :: rZ
  real(8), dimension(1:kN) :: Sk

  integer :: j, NN
  real(8) :: twopi, kZ, kr, sumqcoskr, sumqsinkr

  twopi = 2.0d0*pi

  do NN = 1, kN

    kZ = twopi*real(NN,8)*invlZ
    sumqcoskr = 0.0d0
    sumqsinkr = 0.0d0

    do j = 1, hnion

      kr = kZ*rZ(j)
      sumqcoskr = sumqcoskr + cgNa*cos(kr)
      sumqsinkr = sumqsinkr + cgNa*sin(kr)

    end do

    do j = hnion+1, nion

      kr = kZ*rZ(j)
      sumqcoskr = sumqcoskr + cgCl*cos(kr)
      sumqsinkr = sumqsinkr + cgCl*sin(kr)

    end do

    do j = nion+1, nion+hnion

      kr = kZ*rZ(j)
      sumqcoskr = sumqcoskr + zNa*cos(kr)
      sumqsinkr = sumqsinkr + zNa*sin(kr)

    end do

    do j = nion+hnion+1, twonion

      kr = kZ*rZ(j)
      sumqcoskr = sumqcoskr + zCl*cos(kr)
      sumqsinkr = sumqsinkr + zCl*sin(kr)

    end do

    Sk(NN) = Sk(NN) + sumqcoskr*sumqcoskr + sumqsinkr*sumqsinkr

  end do

return

end subroutine record_Sk

subroutine output_Sk(hnion, kN, nrecord, invlX, invlY, invlZ,&
     & cgNa, cgCl, zNa, zCl, pi, Sk, betapeps)
  implicit none

  integer :: hnion, kN, nrecord
  real(8) :: invlX, invlY, invlZ
  real(8) :: cgNa, cgCl, zNa, zCl, pi, betapeps
  real(8), dimension(1:kN) :: Sk

  integer :: NN
  real(8) :: twopi, kZ, invV, I

  twopi = 2.0d0*pi
  invV = invlX*invlY*invlZ

  I = hnion*(cgNa*cgNa + cgCl*cgCl + zNa*zNa + zCl*zCl)*invV

  Sk = Sk*invV/real(nrecord,8)

  open(30,file='Sk_I.dat')
  open(31,file='Sk_ksq.dat')

  do NN = 1, kN

    kZ = twopi*real(NN,8)*invlZ

    write(30,*) kZ, Sk(NN)/I
    write(31,*) kZ, Sk(NN)/(kZ*kZ)/betapeps

  end do

  close(30)
  close(31)

return

end subroutine output_Sk

