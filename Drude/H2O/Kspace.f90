subroutine Kspace_shell( r4pie, RX, RY, RZ, FX, FY, FZ)
  use basis_set
  use long_range
  implicit none

  real(8) :: r4pie
  real(8), dimension(1:nT) :: RX, RY, RZ
  real(8), dimension(1:nM) :: FX, FY, FZ

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

  integer :: j, k
  integer :: kL, kM, kN, LL, MM, NN, kMmin, kNmin
  real(8) :: ssX, ssY, ssZ, rkX, rkY, rkZ, rksq, invrksq
  real(8) :: p4alp2, twopi, akk, sumqcoskr, sumqsinkr, cs
  real(8) :: F, Fcoeffi, rcpcut, rcpct2
  real(8), dimension(1:nT,0:kLmax) :: eLc, eLs
  real(8), dimension(1:nT,0:kMmax) :: eMc, eMs
  real(8), dimension(1:nT,0:kNmax) :: eNc, eNs 
  real(8), dimension(1:nT) :: clm, slm, coskr, sinkr 
   
  p4alp2 = -0.250d0/alpha2 
  twopi = 2.0d0*pi
  

  FX = 0.0d0
  FY = 0.0d0
  FZ = 0.0d0

  ! set up reciprocal space cutoff, note: rcpcut is twice of rcpcut in pme/gse
  rcpcut = min(real(kLmax,8)*invlX,real(kMmax,8)*invlY,real(kNmax,8)*invlZ)
  rcpcut = rcpcut*1.050d0*twopi
  rcpct2 = rcpcut**2 

  do j = 1, nT
 
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

   do j = 1, nT
     ! cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
     eMc(j,kM) = eMc(j,kM-1)*eMc(j,1) - eMs(j,kM-1)*eMs(j,1)
     ! sin(a+b) = sin(a)cos(b) + cos(a)sin(b)
     eMs(j,kM) = eMs(j,kM-1)*eMc(j,1) + eMc(j,kM-1)*eMs(j,1)

   end do

 end do

 do kN = 2, kNMax

   do j = 1, nT

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

     do j = 1, nT

       eLc(j,0) = eLc(j,1); eLs(j,0) = eLs(j,1)

     end do

   else if ( kL .gt. 1) then

     do j = 1, nT
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
       do j = 1, nT

         clm(j) = eLc(j,0)*eMc(j,kM) - eLs(j,0)*eMs(j,kM)
         slm(j) = eLs(j,0)*eMc(j,kM) + eLc(j,0)*eMs(j,kM)

       end do

     else
       ! now MM is negative but kM is positive so cos(a-b)
       do j = 1, nT

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

           do j = 1, nT

             coskr(j) = clm(j)*enc(j,kN) - slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) + clm(j)*ens(j,kN)

           end do

         else

           do j = 1, nT

             coskr(j) = clm(j)*enc(j,kN) + slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) - clm(j)*ens(j,kN)

           end do

         end if

         ! sum_{i} q_i*e^{i*k*r_i}
         sumqcoskr = 0.0d0
         sumqsinkr = 0.0d0
 
        do j = 1, nM 

           sumqcoskr = sumqcoskr +cgO*coskr(j)
           sumqsinkr = sumqsinkr +cgO*sinkr(j)

        end do 

        do j = nM+1, nT, 3

           sumqcoskr = sumqcoskr +zO*coskr(j)
           sumqsinkr = sumqsinkr +zO*sinkr(j)

           sumqcoskr = sumqcoskr +zH*coskr(j+1)
           sumqsinkr = sumqsinkr +zH*sinkr(j+1)

           sumqcoskr = sumqcoskr +zH*coskr(j+2)
           sumqsinkr = sumqsinkr +zH*sinkr(j+2)

         end do

         invrksq = 1.0d0/rksq
         akk = exp(p4alp2*rksq)*invrksq

         !!! calculate force between shell and shell 

         do j = 1, nM

           k = nM + 3*(j-1)           
             !!! No electrostatic interactions within the molecule 
           F = akk*cgO*( sinkr(j)*( sumqcoskr-zO*coskr(k+1)   & 
                                           & -zH*coskr(k+2)   &
                                           & -zH*coskr(k+3) ) & 
                     & - coskr(j)*( sumqsinkr-zO*sinkr(k+1)   &
                                           & -zH*sinkr(k+2)   &  
                                           & -zH*sinkr(k+3) ) )   

           FX(j) = FX(j) + rkX*F
           FY(j) = FY(j) + rkY*F
           FZ(j) = FZ(j) + rkZ*F

         end do

       end if
  
     end do
 
     kNmin = -kNmax

   end do

   kMmin = -kMmax
   
 end do

 Fcoeffi = twopi*r4pie*invlX*invlY*invlZ

 Fcoeffi = 4.0d0*Fcoeffi 

 do j = 1, nM
   
   FX(j) = FX(j)*Fcoeffi
   FY(j) = FY(j)*Fcoeffi
   FZ(j) = FZ(j)*Fcoeffi
   
 end do
  
return

end subroutine Kspace_shell



subroutine Kspace_core( r4pie, RX, RY, RZ, FX, FY, FZ)
  use basis_set
  use long_range
  implicit none

  real(8) :: r4pie
  real(8), dimension(1:nT) :: RX, RY, RZ
  real(8), dimension(nM+1:nT) :: FX, FY, FZ

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

  integer :: j, k
  integer :: kL, kM, kN, LL, MM, NN, kMmin, kNmin
  real(8) :: ssX, ssY, ssZ, rkX, rkY, rkZ, rksq, invrksq
  real(8) :: p4alp2, twopi, akk, sumqcoskr, sumqsinkr, cs
  real(8) :: F, Fcoeffi, rcpcut, rcpct2
  real(8), dimension(1:nT,0:kLmax) :: eLc, eLs
  real(8), dimension(1:nT,0:kMmax) :: eMc, eMs
  real(8), dimension(1:nT,0:kNmax) :: eNc, eNs 
  real(8), dimension(1:nT) :: clm, slm, coskr, sinkr 
   
  p4alp2 = -0.250d0/alpha2 
  twopi = 2.0d0*pi
  

  FX = 0.0d0
  FY = 0.0d0
  FZ = 0.0d0

  ! set up reciprocal space cutoff, note: rcpcut is twice of rcpcut in pme/gse
  rcpcut = min(real(kLmax,8)*invlX,real(kMmax,8)*invlY,real(kNmax,8)*invlZ)
  rcpcut = rcpcut*1.050d0*twopi
  rcpct2 = rcpcut**2 

  do j = 1, nT
 
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

   do j = 1, nT
     ! cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
     eMc(j,kM) = eMc(j,kM-1)*eMc(j,1) - eMs(j,kM-1)*eMs(j,1)
     ! sin(a+b) = sin(a)cos(b) + cos(a)sin(b)
     eMs(j,kM) = eMs(j,kM-1)*eMc(j,1) + eMc(j,kM-1)*eMs(j,1)

   end do

 end do

 do kN = 2, kNMax

   do j = 1, nT

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

     do j = 1, nT

       eLc(j,0) = eLc(j,1); eLs(j,0) = eLs(j,1)

     end do

   else if ( kL .gt. 1) then

     do j = 1, nT
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
       do j = 1, nT

         clm(j) = eLc(j,0)*eMc(j,kM) - eLs(j,0)*eMs(j,kM)
         slm(j) = eLs(j,0)*eMc(j,kM) + eLc(j,0)*eMs(j,kM)

       end do

     else
       ! now MM is negative but kM is positive so cos(a-b)
       do j = 1, nT

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

           do j = 1, nT

             coskr(j) = clm(j)*enc(j,kN) - slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) + clm(j)*ens(j,kN)

           end do

         else

           do j = 1, nT

             coskr(j) = clm(j)*enc(j,kN) + slm(j)*ens(j,kN)
             sinkr(j) = slm(j)*enc(j,kN) - clm(j)*ens(j,kN)

           end do

         end if

         ! sum_{i} q_i*e^{i*k*r_i}
         sumqcoskr = 0.0d0
         sumqsinkr = 0.0d0
 
        do j = 1, nM 

           sumqcoskr = sumqcoskr +cgO*coskr(j)
           sumqsinkr = sumqsinkr +cgO*sinkr(j)

        end do 

        do j = nM+1, nT, 3

           sumqcoskr = sumqcoskr +zO*coskr(j)
           sumqsinkr = sumqsinkr +zO*sinkr(j)

           sumqcoskr = sumqcoskr +zH*coskr(j+1)
           sumqsinkr = sumqsinkr +zH*sinkr(j+1)

           sumqcoskr = sumqcoskr +zH*coskr(j+2)
           sumqsinkr = sumqsinkr +zH*sinkr(j+2)

         end do

         invrksq = 1.0d0/rksq
         akk = exp(p4alp2*rksq)*invrksq

         !!! calculate force between shell and shell 

         do j = nM+1, nT
           !!! No electrostatic interactions within the molecule
           k = floor( real(j-nM-1,8)/3.0d0 + 1.0d-10 ) + 1       

           if ( mod((j-nM), 3) == 1 ) then !!! O core
 
             F = akk*zO*( sinkr(j)*( sumqcoskr-cgO*coskr(k)   &
                                            & -zH*coskr(j+1)   &
                                            & -zH*coskr(j+2) ) &
                      & - coskr(j)*( sumqsinkr-cgO*sinkr(k)   &
                                            & -zH*sinkr(j+1)   &
                                            & -zH*sinkr(j+2) ) )

           else if (  mod((j-nM), 3) == 2 ) then !!! H1 

             F = akk*zH*( sinkr(j)*( sumqcoskr-cgO*coskr(k)   &
                                            & -zO*coskr(j-1)   &
                                            & -zH*coskr(j+1) ) &
                      & - coskr(j)*( sumqsinkr-cgO*sinkr(k)   &
                                            & -zO*sinkr(j-1)   &
                                            & -zH*sinkr(j+1) ) )

           else !!! mod((j-nM), 3) == 0 !!! H2

             F = akk*zH*( sinkr(j)*( sumqcoskr-cgO*coskr(k)   &
                                            & -zO*coskr(j-2)   &
                                            & -zH*coskr(j-1) ) &
                      & - coskr(j)*( sumqsinkr-cgO*sinkr(k)   &
                                            & -zO*sinkr(j-2)   &
                                            & -zH*sinkr(j-1) ) )

           end if

           FX(j) = FX(j) + rkX*F
           FY(j) = FY(j) + rkY*F
           FZ(j) = FZ(j) + rkZ*F

         end do

       end if
  
     end do
 
     kNmin = -kNmax

   end do

   kMmin = -kMmax
   
 end do

 Fcoeffi = twopi*r4pie*invlX*invlY*invlZ

 Fcoeffi = 4.0d0*Fcoeffi 

 do j = nM+1, nT
   
   FX(j) = FX(j)*Fcoeffi
   FY(j) = FY(j)*Fcoeffi
   FZ(j) = FZ(j)*Fcoeffi
   
 end do
  
return

end subroutine Kspace_core
























