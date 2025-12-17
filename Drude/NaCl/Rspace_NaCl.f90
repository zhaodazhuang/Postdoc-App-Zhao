subroutine Rspace_shell(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ, &
     & invlX, invlY, invlZ, FRX, FRY, FRZ) 
  use short_range
  use Coulombic
  implicit none

  integer :: nion, hnion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8) :: lX, lY, lZ, invlX, invlY, invlZ
  real(8), dimension(1:nion) :: FRX, FRY, FRZ
  
  integer :: i, j
  real(8) :: F, FF, coeffi
  real(8) :: rXj, rYj, rZj, rijX, rijY, rijZ
  real(8) :: rij, rij2, invrij, invrij2

  coeffi = 2.0d0*alpha/sqrt(pi)

  FRX = 0.0d0
  FRY = 0.0d0
  FRZ = 0.0d0

  do j =1, nion

    rXj = rX(j); rYj = rY(j); rZj = rZ(j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Na shell
    if ( j .le. hnion ) then

      do i = 1, twonion

        rijX = rX(i)-rXj ; rijX = rijX - anint(rijX*invlX)*lX
        rijY = rY(i)-rYj ; rijY = rijY - anint(rijY*invlY)*lY
        rijZ = rZ(i)-rZj ; rijZ = rijZ - anint(rijZ*invlZ)*lZ

        rij2 = rijX*rijX + rijY*rijY + rijZ*rijZ

        if ( i .eq. j + nion ) then
           
          F = kNa*eang

          FRX(j) = FRX(j) + F*rijX
          FRY(j) = FRY(j) + F*rijY
          FRZ(j) = FRZ(j) + F*rijZ

        end if

        !!!----------------- shell -- other shell
        if ( rij2 .le. cutoff2 ) then
 
            invrij2 = 1.0d0/rij2
            rij = sqrt(rij2)
            invrij = 1.0d0/rij
            
          if ( i .le. hnion ) then
 
            if ( i .ne. j  ) then
         
            F = -ANaNa*exp(-rij*invrho)*invrho*invrij + 6.0d0*CNaNa*invrij2*invrij2*invrij2*invrij2&
              & + 8.0d0*DNaNa*invrij2*invrij2*invrij2*invrij2*invrij2&
              & + 3.0d0*c3NaNa*rij + 4.0d0*c4NaNa*rij2  
            F = F*eang
           
            FF= -cgNa*cgNa*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            FF= FF*r4pie

            F = F + FF

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

            end if
           
          else if ( i .gt. hnion .and. i .le. nion ) then

            F = -ANaCl*exp(-rij*invrho)*invrho*invrij + 6.0d0*CNaCl*invrij2*invrij2*invrij2*invrij2&
              & + 8.0d0*DNaCl*invrij2*invrij2*invrij2*invrij2*invrij2&
              & + 3.0d0*c3NaCl*rij + 4.0d0*c4NaCl*rij2
            F = F*eang

            FF= -cgNa*cgCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            FF= FF*r4pie

            F = F + FF

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

          !!!-------------------- shell -- core          
          else if ( i .gt. nion .and. i .le. nion+hnion ) then
  
            if ( i .ne. j + nion ) then

            F = -zNa*cgNa*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                         & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ
 
            end if
 
          else if ( i .gt. nion+hnion .and. i .le. twonion ) then
  
            F = -zCl*cgNa*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                         & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie
  
            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ
   
          else

            write(*,*) i
            write(*,*) 'error: index beyond 2*nion at Rspace_NaCl'
            write(*,*) 'stop at one'
            stop
  
          end if

        end if !!! if ( rij2 .le. cutoff2 )

      end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Cl shell
    else if ( j .gt. hnion ) then 

      do i = 1, twonion

        rijX = rX(i)-rXj ; rijX = rijX - anint(rijX*invlX)*lX
        rijY = rY(i)-rYj ; rijY = rijY - anint(rijY*invlY)*lY
        rijZ = rZ(i)-rZj ; rijZ = rijZ - anint(rijZ*invlZ)*lZ

        rij2 = rijX*rijX + rijY*rijY + rijZ*rijZ

        if ( i .eq. j + nion ) then

          F = kCl*eang

          FRX(j) = FRX(j) + F*rijX
          FRY(j) = FRY(j) + F*rijY
          FRZ(j) = FRZ(j) + F*rijZ

        end if

        !!!----------------- shell -- other shell
        if ( rij2 .le. cutoff2 ) then

            invrij2 = 1.0d0/rij2
            rij = sqrt(rij2)
            invrij = 1.0d0/rij

          if ( i .le. hnion ) then

            F = -ANaCl*exp(-rij*invrho)*invrho*invrij + 6.0d0*CNaCl*invrij2*invrij2*invrij2*invrij2&
              & + 8.0d0*DNaCl*invrij2*invrij2*invrij2*invrij2*invrij2&
              & + 3.0d0*c3NaCl*rij + 4.0d0*c4NaCl*rij2
            F = F*eang

            FF= -cgNa*cgCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            FF= FF*r4pie

            F = F + FF

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ
 
          else if ( i .gt. hnion .and. i .le. nion ) then

            if ( i .ne. j  ) then

            F = -AClCl*exp(-rij*invrho)*invrho*invrij + 6.0d0*CClCl*invrij2*invrij2*invrij2*invrij2&
              & + 8.0d0*DClCl*invrij2*invrij2*invrij2*invrij2*invrij2&
              & + 3.0d0*c3ClCl*rij + 4.0d0*c4ClCl*rij2
            F = F*eang

            FF= -cgCl*cgCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij& 
                          & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            FF= FF*r4pie

            F =  F + FF

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

            end if

          else if ( i .gt. nion .and. i .le. nion+hnion ) then

            F = -zNa*cgCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                         & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

          else if ( i .gt. nion+hnion .and. i .le. twonion ) then

            if ( i .ne. j + nion ) then

            F = -zCl*cgCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                         & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

            end if

          else

            write(*,*) 'error: index beyond 2*nion at Rspace_NaCl'
            write(*,*) 'stop at two'
            stop

          end if

        end if !!! if ( rij2 .le. cutoff2 )
     
      end do

    end if !!! if ( j .le. hnion )

  end do

return        

end subroutine Rspace_shell 






subroutine Rspace_core(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ, &
     & invlX, invlY, invlZ, FRX, FRY, FRZ, kNa, kCl, cutoff2, eang)
  use Coulombic
  implicit none

  integer :: nion, hnion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8) :: lX, lY, lZ, invlX, invlY, invlZ
  real(8), dimension(nion+1:twonion) :: FRX, FRY, FRZ
  real(8) :: kNa, kCl, cutoff2, eang

  integer :: i, j
  real(8) :: F, FF, coeffi
  real(8) :: rXj, rYj, rZj, rijX, rijY, rijZ
  real(8) :: rij, rij2, invrij, invrij2

  coeffi = 2.0d0*alpha/sqrt(pi)

  FRX = 0.0d0
  FRY = 0.0d0
  FRZ = 0.0d0

  do j = nion+1, twonion

    rXj = rX(j); rYj = rY(j); rZj = rZ(j)

    if ( j .le. nion+hnion ) then

      do i = 1, twonion

        rijX = rX(i)-rXj ; rijX = rijX - anint(rijX*invlX)*lX
        rijY = rY(i)-rYj ; rijY = rijY - anint(rijY*invlY)*lY
        rijZ = rZ(i)-rZj ; rijZ = rijZ - anint(rijZ*invlZ)*lZ

        rij2 = rijX*rijX + rijY*rijY + rijZ*rijZ

        if ( i .eq. j - nion ) then

          F = kNa*eang

          FRX(j) = FRX(j) + F*rijX
          FRY(j) = FRY(j) + F*rijY
          FRZ(j) = FRZ(j) + F*rijZ

        end if

        !!!----------------- core -- other shell
        if ( rij2 .le. cutoff2 ) then

            invrij2 = 1.0d0/rij2
            rij = sqrt(rij2)
            invrij = 1.0d0/rij

          if ( i .le. hnion ) then

            if ( i .ne. j - nion ) then 

            F = -cgNa*zNa*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                         & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

            end if

          else if ( i .gt. hnion .and. i .le. nion ) then

            F = -cgCl*zNa*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                         & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

          !!!-------------------- core -- core
          else if ( i .gt. nion .and. i .le. nion+hnion ) then

            if ( i .ne. j ) then

            F = -zNa*zNa*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                        & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

            end if

          else if ( i .gt. nion+hnion .and. i .le. twonion ) then

            F = -zCl*zNa*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                        & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

          else

            write(*,*) i
            write(*,*) 'error: index beyond 2*nion at Rspace_NaCl'
            write(*,*) 'stop at one'
            stop

          end if

        end if !!! if ( rij2 .le. cutoff2 )

      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Cl core
    else if ( j .gt. nion+hnion ) then

      do i = 1, twonion

        rijX = rX(i)-rXj ; rijX = rijX - anint(rijX*invlX)*lX
        rijY = rY(i)-rYj ; rijY = rijY - anint(rijY*invlY)*lY
        rijZ = rZ(i)-rZj ; rijZ = rijZ - anint(rijZ*invlZ)*lZ

        rij2 = rijX*rijX + rijY*rijY + rijZ*rijZ

        if ( i .eq. j - nion ) then

          F = kCl*eang

          FRX(j) = FRX(j) + F*rijX
          FRY(j) = FRY(j) + F*rijY
          FRZ(j) = FRZ(j) + F*rijZ

        end if

        !!!----------------- shell -- other shell
        if ( rij2 .le. cutoff2 ) then

            invrij2 = 1.0d0/rij2
            rij = sqrt(rij2)
            invrij = 1.0d0/rij

          if ( i .le. hnion ) then

            F= -cgNa*zCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                        & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F= F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

          else if ( i .gt. hnion .and. i .le. nion ) then

            if ( i .ne. j - nion ) then

            F= -cgCl*zCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                        & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F= F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

            end if

          else if ( i .gt. nion .and. i .le. nion+hnion ) then

            F = -zNa*zCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                        & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

          else if ( i .gt. nion+hnion .and. i .le. twonion ) then

            if ( i .ne. j ) then

            F = -zCl*zCl*(coeffi*exp(-alpha2*rij2)*invrij2 + erfc(alpha*rij)*invrij2*invrij&
                        & - 3.0d0*c3e*rij - 4.0d0*c4e*rij2 )
            F = F*r4pie

            FRX(j) = FRX(j) + F*rijX
            FRY(j) = FRY(j) + F*rijY
            FRZ(j) = FRZ(j) + F*rijZ

            end if

          else

            write(*,*) 'error: index beyond 2*nion at Rspace_NaCl'
            write(*,*) 'stop at two'
            stop

          end if

        end if !!! if ( rij2 .le. cutoff2 )

      end do

    end if !!! if ( j .le. hnion )

  end do

return

end subroutine Rspace_core












   
