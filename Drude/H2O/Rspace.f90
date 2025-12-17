subroutine Rspace_shell( alpha, eang, r4pie, RX, RY, RZ, FX, FY, FZ)
  use basis_set
  use short_range
  implicit none

  real(8) ::  alpha, eang, r4pie
  real(8), dimension(1:nT) :: RX, RY, RZ
  real(8), dimension(1:nM) :: FX, FY, FZ

  integer :: i, j, k
  real(8) :: alpha2, coeffi
  real(8) :: rXj, rYj, rZj, rXij, rYij, rZij
  real(8) :: rij, rij2, invrij, invrij2
  real(8) :: XH1j, YH1j, ZH1j, rH1j, ll
  real(8) :: XH2j, YH2j, ZH2j, rH2j, theta
  real(8) :: f, ff, fxx, fyy, fzz

  FX = 0.0d0
  FY = 0.0d0
  FZ = 0.0d0
  
  alpha2 = alpha*alpha

  coeffi = 2.0d0*alpha/sqrt(pi)

  do j = 1, nM

    rXj = RX(j)
    rYj = RY(j)
    rZj = RZ(j)
  
    k = nM + 3*(j-1)

    do i = 1, nT

      if ( i == j ) then

        XH1j = rX(k+2)-rXj ; XH1j = XH1j - anint(XH1j*invlX)*lX
        YH1j = rY(k+2)-rYj ; YH1j = YH1j - anint(YH1j*invlY)*lY
        ZH1j = rZ(k+2)-rZj ; ZH1j = ZH1j - anint(ZH1j*invlZ)*lZ
        rH1j = sqrt(XH1j*XH1j + YH1j*YH1j + ZH1j*ZH1j)

        XH2j = rX(k+3)-rXj ; XH2j = XH2j - anint(XH2j*invlX)*lX
        YH2j = rY(k+3)-rYj ; YH2j = YH2j - anint(YH2j*invlY)*lY
        ZH2j = rZ(k+3)-rZj ; ZH2j = ZH2j - anint(ZH2j*invlZ)*lZ
        rH2j = sqrt(XH2j*XH2j + YH2j*YH2j + ZH2j*ZH2j)

        ll = (XH1j*XH2j + YH1j*YH2j + ZH1j*ZH2j)/(rH1j*rH2j)
        
        theta = acos(ll)

        f = -eang*Kijk*(theta-theta0)*( -1.0d0/(sqrt(1-ll*ll)) )
   
        fxx = f*( 1.0d0/(rH1j*rH2j)*(-XH1j-XH2j) + &
                 ll/(rH1j*rH2j)*(rH2j/rH1j*XH1j + rH1j/rH2j*XH2j) )
        fyy = f*( 1.0d0/(rH1j*rH2j)*(-YH1j-YH2j) + &
                 ll/(rH1j*rH2j)*(rH2j/rH1j*YH1j + rH1j/rH2j*YH2j) )
        fzz = f*( 1.0d0/(rH1j*rH2j)*(-ZH1j-ZH2j) + &
                 ll/(rH1j*rH2j)*(rH2j/rH1j*ZH1j + rH1j/rH2j*ZH2j) )

        FX(j) = FX(j) + fxx
        FY(j) = FY(j) + fyy
        FZ(j) = FZ(j) + fzz
        
      else 

        rXij = rX(i)-rXj ; rXij = rXij - anint(rXij*invlX)*lX
        rYij = rY(i)-rYj ; rYij = rYij - anint(rYij*invlY)*lY
        rZij = rZ(i)-rZj ; rZij = rZij - anint(rZij*invlZ)*lZ
       
        rij2 = rXij*rXij + rYij*rYij + rZij*rZij

        if ( rij2 .le. Rcut2 ) then

          invrij2 = 1.0d0/rij2
          rij = sqrt(rij2)
          invrij = 1.0d0/rij

          if ( i .eq. k + 1 ) then

            f = Kappa*eang

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else if ( i == k + 2 .or. i == k + 3 ) then

            !!! Morse
            f = 2.0d0*Moralpha*MorD*( 1.0d0 - exp(-Moralpha*(rij-Morr0)) )&
               *exp(-Moralpha*(rij-Morr0))/rij          

            f = f*eang

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij
           
          else if ( i <= nM ) then !!! shell --  shell

            !!! erfc + Lennard-Jones 
            f = - (12.0d0*LJA*invrij2**7 - 6.0d0*LJB*invrij2**4 &
                  - 3.0d0*LJc3*rij - 4.0d0*LJc4*rij2 ) 

            ff= -cgO*cgO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )

            f = f*eang + ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else  !!! shell - core

            if ( mod(i-nM, 3) == 1 ) then  !!! O core -- shell

              ff= -cgO*zO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
              f = ff*r4pie

            else  !!! H core -- shell
               
              f = -( hamA*invrho*exp(-rij*invrho)*invrij - 6.0d0*hamC*invrij2*invrij2*invrij2*invrij2&
                    -3.0d0*hamc3*rij - 4.0d0*hamc4*rij2 )

              ff = -cgO*zH*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
 
              f = f*eang + ff*r4pie

            end if 

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          end if !!! j <= nM 

        end if !!! i .eq. j + nM

      end if !!! i == j

    end do !!! loop for i

  end do !!! loop for j
    
return

end subroutine Rspace_shell


subroutine Rspace_core( alpha, eang, r4pie, RX, RY, RZ, FX, FY, FZ )
  use basis_set
  use short_range
  implicit none

  real(8) ::  alpha, eang, r4pie
  real(8), dimension(1:nT) :: RX, RY, RZ
  real(8), dimension(nM+1:nT) :: FX, FY, FZ

  integer :: i, j, k
  real(8) :: alpha2, coeffi
  real(8) :: rXj, rYj, rZj, rXij, rYij, rZij
  real(8) :: rij, rij2, invrij, invrij2
  real(8) :: Xjs, Yjs, Zjs, rjs, ll
  real(8) :: XHs, YHs, ZHs, rHs, theta
  real(8) :: f, ff, fxx, fyy, fzz

  FX = 0.0d0
  FY = 0.0d0
  FZ = 0.0d0

  alpha2 = alpha*alpha

  coeffi = 2.0d0*alpha/sqrt(pi)

  do j = nM+1, nT

    rXj = RX(j)
    rYj = RY(j)
    rZj = RZ(j)

    k = floor( real(j-nM-1,8)/3.0d0 + 1.0d-10 ) + 1
  
    do i = 1, nT

    if ( i == j ) then

      if ( mod((j-nM),3) == 2 ) then !!! H1

        Xjs = rXj - rX(k) ; Xjs = Xjs - anint(Xjs*invlX)*lX
        Yjs = rYj - rY(k) ; Yjs = Yjs - anint(Yjs*invlY)*lY
        Zjs = rZj - rZ(k) ; Zjs = Zjs - anint(Zjs*invlZ)*lZ
        rjs = sqrt(Xjs*Xjs + Yjs*Yjs + Zjs*Zjs)

        XHs = rX(j+1) - rX(k) ; XHs = XHs - anint(XHs*invlX)*lX
        YHs = rY(j+1) - rY(k) ; YHs = YHs - anint(YHs*invlY)*lY
        ZHs = rZ(j+1) - rZ(k) ; ZHs = ZHs - anint(ZHs*invlZ)*lZ
        rHs = sqrt(XHs*XHs + YHs*YHs + ZHs*ZHs)

        ll = (Xjs*XHs + Yjs*YHs + Zjs*ZHs)/(rjs*rHs)

        theta = acos(ll)

        f = -eang*Kijk*(theta-theta0)*( -1.0d0/(sqrt(1-ll*ll)) )

        fxx = f*( 1.0d0/(rjs*rHs)*XHs - &
                 ll/(rjs*rHs)*(rHs/rjs)*Xjs )
        fyy = f*( 1.0d0/(rjs*rHs)*YHs - &
                 ll/(rjs*rHs)*(rHs/rjs)*Yjs )
        fzz = f*( 1.0d0/(rjs*rHs)*ZHs - &
                 ll/(rjs*rHs)*(rHs/rjs)*Zjs )

        FX(j) = FX(j) + fxx
        FY(j) = FY(j) + fyy
        FZ(j) = FZ(j) + fzz
        
      else if ( mod((j-nM),3) == 0 ) then !!! H2
 
        Xjs = rXj - rX(k) ; Xjs = Xjs - anint(Xjs*invlX)*lX
        Yjs = rYj - rY(k) ; Yjs = Yjs - anint(Yjs*invlY)*lY
        Zjs = rZj - rZ(k) ; Zjs = Zjs - anint(Zjs*invlZ)*lZ
        rjs = sqrt(Xjs*Xjs + Yjs*Yjs + Zjs*Zjs)

        XHs = rX(j-1) - rX(k) ; XHs = XHs - anint(XHs*invlX)*lX
        YHs = rY(j-1) - rY(k) ; YHs = YHs - anint(YHs*invlY)*lY
        ZHs = rZ(j-1) - rZ(k) ; ZHs = ZHs - anint(ZHs*invlZ)*lZ
        rHs = sqrt(XHs*XHs + YHs*YHs + ZHs*ZHs)

        ll = (Xjs*XHs + Yjs*YHs + Zjs*ZHs)/(rjs*rHs)

        theta = acos(ll)

        f = -eang*Kijk*(theta-theta0)*( -1.0d0/(sqrt(1-ll*ll)) )

        fxx = f*( 1.0d0/(rjs*rHs)*XHs - &
                 ll/(rjs*rHs)*(rHs/rjs)*Xjs )
        fyy = f*( 1.0d0/(rjs*rHs)*YHs - &
                 ll/(rjs*rHs)*(rHs/rjs)*Yjs )
        fzz = f*( 1.0d0/(rjs*rHs)*ZHs - &
                 ll/(rjs*rHs)*(rHs/rjs)*Zjs )

        FX(j) = FX(j) + fxx
        FY(j) = FY(j) + fyy
        FZ(j) = FZ(j) + fzz

      end if  

      !!! O core no thing 

    else 

      rXij = rX(i)-rXj ; rXij = rXij - anint(rXij*invlX)*lX
      rYij = rY(i)-rYj ; rYij = rYij - anint(rYij*invlY)*lY
      rZij = rZ(i)-rZj ; rZij = rZij - anint(rZij*invlZ)*lZ

      rij2 = rXij*rXij + rYij*rYij + rZij*rZij

      if ( rij2 <= Rcut2 ) then  

          invrij2 = 1.0d0/rij2
          rij = sqrt(rij2)
          invrij = 1.0d0/rij

        if ( mod((j-nM),3) == 1 ) then !!! O core
 
          if ( i == k ) then

            f = Kappa*eang

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij
        
          else if ( i <= nM ) then

            ff = -zO*cgO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else if ( mod((i-nM),3) == 1 ) then

            ff = -zO*zO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else 

            if ( i .ne. j + 1 .and. i .ne. j + 2 ) then

              ff = -zO*zH*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
              f = ff*r4pie

              FX(j) = FX(j) + f*rXij
              FY(j) = FY(j) + f*rYij
              FZ(j) = FZ(j) + f*rZij

            end if !!! i .ne. j + 1 .and. i .ne. j + 2
 
          end if !!! i == k 

        else if ( mod((j-nM),3) == 2 ) then !!! H1

          if ( i == k ) then 

            f = 2.0d0*Moralpha*MorD*( 1.0d0 - exp(-Moralpha*(rij-Morr0)) )&
                *exp(-Moralpha*(rij-Morr0))/rij

            f = f*eang

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else if ( i <= nM )  then  !!! O shell

            f = -( hamA*invrho*exp(-rij*invrho)*invrij - 6.0d0*hamC*invrij2*invrij2*invrij2*invrij2&
                  -3.0d0*hamc3*rij - 4.0d0*hamc4*rij2 )

            ff = -zH*cgO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )

            f = f*eang + ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else if ( mod((i-nM),3) == 1 ) then !!! O core

            if ( i .ne. j - 1 ) then !!! no erfc/r between H and own O core

            ff = -zH*zO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij
               
            end if

          else if ( mod((i-nM),3) == 0 ) then !!! H2 

            if ( i .ne. j + 1 ) then !!! no erfc/r between H1 and own O core and H2

            ff = -zH*zH*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

            end if

          else !!! H1

            ff = -zH*zH*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          end if !!! i == k 

        else !!! H2

          if ( i == k ) then

            f = 2.0d0*Moralpha*MorD*( 1.0d0 - exp(-Moralpha*(rij-Morr0)) )&
                *exp(-Moralpha*(rij-Morr0))/rij

            f = f*eang

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else if ( i <= nM )  then  !!! O shell

            f = -( hamA*invrho*exp(-rij*invrho)*invrij - 6.0d0*hamC*invrij2*invrij2*invrij2*invrij2&
                  -3.0d0*hamc3*rij - 4.0d0*hamc4*rij2 )

            ff = -zH*cgO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = f*eang + ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          else if ( mod((i-nM),3) == 1 ) then !!! O core

            if ( i .ne. j - 2 ) then !!! no erfc/r between H and own O core

            ff = -zH*zO*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

            end if
            
          else if ( mod((i-nM),3) == 2 ) then !!! H1

            if ( i .ne. j - 1 ) then !!! no erfc/r between H1 and own O core and H2

            ff = -zH*zH*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

            end if

          else !!! H2

            ff = -zH*zH*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            f = ff*r4pie

            FX(j) = FX(j) + f*rXij
            FY(j) = FY(j) + f*rYij
            FZ(j) = FZ(j) + f*rZij

          end if !!! i == k

        end if !!! O core

      end if !!!  rij2 <= Rcut2 

    end if !!! i == j

    end do !!! loop for i
    
  end do !!! loop for j
   
return

end subroutine Rspace_core
