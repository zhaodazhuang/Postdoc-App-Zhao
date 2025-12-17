subroutine initialize_pos(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ) 
   implicit none

  integer :: nion, hnion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ  
  real(8) :: lX, lY, lZ

  integer :: i, j, k, mark, c, d
  real(8) :: x, y, z, hlX, hlY, hlZ, dlX, dlY, dlZ
  character(2), dimension(1:nion) :: aa, aaa
  real(8), dimension(1:twonion) :: XX, YY, ZZ

  c = 1
  d = 0

  dlX = lX/6.0d0
  dlY = lY/6.0d0
  dlZ = lZ/6.0d0

  hlX = 0.50d0*dlX
  hlY = 0.50d0*dlY
  hlZ = 0.50d0*dlZ
 
  mark = 1

  do i = 0, 5
     x = real(i,8)*dlX
    do j = 0, 5
       y = real(j,8)*dlY
      do k = 0, 5
         z = real(k,8)*dlZ

         rX(mark) = x; rY(mark) = y; rZ(mark) = z;

         if ( mod(mark,2) .eq. c ) then

           aa(mark) = 'Cl'

         else if ( mod(mark,2) .eq. d ) then

           aa(mark) = 'Na'

         end if

         if ( mod(mark,6) .eq. 0 ) then

           d = c + d; c = d - c; d = d - c

         end if

         if ( mod(mark,36) .eq. 0 ) then

           d = c + d; c = d - c; d = d - c

         end if

         mark = mark + 1

      end do
    end do
  end do

  rX = rX + hlX
  rY = rY + hlY
  rZ = rZ + hlZ

  j = 1
  k = hnion+1
  
  XX = rX(1:nion)
  YY = rY(1:nion)
  ZZ = rZ(1:nion)
  
  do i = 1, nion

    if ( aa(i) .eq. 'Na' ) then
     
      aaa(j) = aa(i)  
      rX(j) = XX(i); rY(j) = YY(i); rZ(j) = ZZ(i)
      j = j + 1

    else if ( aa(i) .eq. 'Cl' ) then

      aaa(k) = aa(i)
      rX(k) = XX(i); rY(k) = YY(i); rZ(k) = ZZ(i)
      k = k +1
  
    end if

  end do

  rX(nion+1:twonion) = rX(1:nion)
  rY(nion+1:twonion) = rY(1:nion)
  rZ(nion+1:twonion) = rZ(1:nion)

  open(25,file='initial_core.xyz')

  write(25,*) nion
  write(25,*) lX, lY, lZ

  do i = 1, nion

    write(25,*) aaa(i), rX(i), rY(i), rZ(i)

  end do

  close(25)

return

end subroutine initialize_pos

subroutine gauss_velocity(nion, hnion, twonion, hvX, hvY, hvZ, mNa, mCl, coiv, reqtemp)
  implicit none

  integer :: nion, hnion, twonion
  real(8) :: mNa, mCl, coiv, reqtemp
  real(8), dimension(nion+1:twonion) :: hvX, hvY, hvZ

  integer :: j
  real(8) :: pi, twopi
  real(8), dimension(nion+1:twonion) :: x, y, z

  pi = acos(-1.0d0) 
  twopi = 2.0d0*pi

  call random_seed()
  call random_number(x)
  call random_number(y)
  z = sqrt( -2.0d0*log(x) )*( cos(twopi*y) ) 

  hvX(nion+1:nion+hnion) = sqrt(coiv*reqtemp/mNa)*z(nion+1:nion+hnion)
  hvX(nion+hnion+1:twonion) = sqrt(coiv*reqtemp/mCl)*z(nion+hnion+1:twonion)

  call random_seed()
  call random_number(x)
  call random_number(y)
  z = sqrt( -2.0d0*log(x) )*( cos(twopi*y) )

  hvY(nion+1:nion+hnion) = sqrt(coiv*reqtemp/mNa)*z(nion+1:nion+hnion)
  hvY(nion+hnion+1:twonion) = sqrt(coiv*reqtemp/mCl)*z(nion+hnion+1:twonion)

  call random_seed()
  call random_number(x)
  call random_number(y)
  z = sqrt( -2.0d0*log(x) )*( cos(twopi*y) )

  hvZ(nion+1:nion+hnion) = sqrt(coiv*reqtemp/mNa)*z(nion+1:nion+hnion)
  hvZ(nion+hnion+1:twonion) = sqrt(coiv*reqtemp/mCl)*z(nion+hnion+1:twonion)

return

end subroutine gauss_velocity

subroutine read_v_coor_hv(nion, twonion, rX, rY, rZ, vX, vY, vZ, hvX, hvY, hvZ)
  implicit none

  integer :: nion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8), dimension(nion+1:twonion) :: vX, vY, vZ, hvX, hvY, hvZ

  integer :: j

  open(33, file='coor.xyz')

  do j = 1, twonion

    read(33,*) rX(j), rY(j), rZ(j)

  end do

  close(33)

  open(34, file='v.dat')

  do j = nion+1, twonion

    read(34,*) vX(j), vY(j), vZ(j)

  end do

  close(34)

  open(35, file='hv.dat')

  do j = nion+1, twonion

    read(35,*) hvX(j), hvY(j), hvZ(j)

  end do

  close(35)

return

end subroutine read_v_coor_hv


subroutine update_rc(nion, hnion, twonion, invmNa, invmCl, rX, rY, rZ, &
         & lX, lY, lZ, invlX, invlY, invlZ, FX, FY, FZ, vX, vY, vZ, &
         & hvX, hvY, hvZ, dt)

  implicit none

  integer :: nion, hnion, twonion
  real(8) :: invmNa, invmCl, lX, lY, lZ, invlX, invlY, invlZ, dt
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8), dimension(nion+1:twonion) :: FX, FY, FZ
  real(8), dimension(nion+1:twonion) :: vX, vY, vZ
  real(8), dimension(nion+1:twonion) :: hvX, hvY, hvZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  use standard leapforg
!!!
!!!  x_n+1 = x_n + dT*v_n+1/2
!!!  v_n+1/2 = v_n-1/2 + dT*a(x_n)
!!!  v_n = 1/2*( v_n+1/2 + v_n-1/2 )
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: j
  real(8) :: aX, aY, aZ, rXj, rYj, rZj
  real(8) :: hdt

  hdt = 0.50d0*dt

  do j = nion+1, twonion

    rXj = rX(j) ; rYj = rY(j) ; rZj = rZ(j)

    if ( j .le. nion+hnion ) then

    aX = FX(j)*invmNa
    aY = FY(j)*invmNa
    aZ = FZ(j)*invmNa

    else if ( j .gt. nion+hnion ) then

    aX = FX(j)*invmCl
    aY = FY(j)*invmCl
    aZ = FZ(j)*invmCl

    end if

    hvX(j) = hvX(j) + dt*aX
    hvY(j) = hvY(j) + dt*aY
    hvZ(j) = hvZ(j) + dt*aZ

    rX(j) = rX(j) + dt*hvX(j)
    rY(j) = rY(j) + dt*hvY(j)
    rZ(j) = rZ(j) + dt*hvZ(j)

    rX(j) = rX(j) - floor(rX(j)*invlX)*lX
    rY(j) = rY(j) - floor(rY(j)*invlY)*lY
    rZ(j) = rZ(j) - floor(rZ(j)*invlZ)*lZ

    vX(j) = vX(j) + hdt*aX
    vY(j) = vY(j) + hdt*aY
    vZ(j) = vZ(j) + hdt*aZ

  end do

return

end subroutine update_rc


subroutine contral_temp(nion, hnion, twonion, hvX, hvY, hvZ, mNa, mCl,&
                      & coitemp, nowtemp, reqtemp)
  implicit none

  integer :: nion, hnion, twonion
  real(8) :: mNa, mCl, coitemp, nowtemp, reqtemp
  real(8), dimension(nion+1:twonion) :: hvX, hvY, hvZ

  real(8) :: beta

  nowtemp = mNa*sum( hvX(nion+1:nion+hnion)*hvX(nion+1:nion+hnion) &
                 & + hvY(nion+1:nion+hnion)*hvY(nion+1:nion+hnion) &
                 & + hvZ(nion+1:nion+hnion)*hvZ(nion+1:nion+hnion) )&
          & + mCl*sum( hvX(nion+hnion+1:twonion)*hvX(nion+hnion+1:twonion) &
                   & + hvY(nion+hnion+1:twonion)*hvY(nion+hnion+1:twonion) &
                   & + hvZ(nion+hnion+1:twonion)*hvZ(nion+hnion+1:twonion) )
  
  nowtemp = nowtemp*coitemp/nion

!  beta = sqrt( reqtemp/nowtemp )
!
!  hvX = hvX*beta 
!  hvY = hvY*beta 
!  hvZ = hvZ*beta

return

end subroutine contral_temp


subroutine record_hist(nion, hnion, twonion, rX, rY, rZ, lX, lY, lZ,&
         & invlX, invlY, invlZ, cutoff2)
  use RDFmod
  implicit none

  integer :: nion, hnion, twonion
  real(8) :: lX, lY, lZ, invlX, invlY, invlZ, cutoff2
  real(8), dimension(1:twonion) :: rX, rY, rZ

  integer :: i, j, bin
  real(8) :: rXj, rYj, rZj, rijX, rijY, rijZ, rij, rij2

  do j = nion+1, nion+hnion-1

      rXj = rX(j)
      rYj = rY(j)
      rZj = rZ(j)

    do i = j+1, nion+hnion

      rijX = rX(i) - rXj
      rijY = rY(i) - rYj
      rijZ = rZ(i) - rZj

      rijX = rijX - anint(rijX*invlX)*lX
      rijY = rijY - anint(rijY*invlY)*lY
      rijZ = rijZ - anint(rijZ*invlZ)*lZ

      rij2 = rijX*rijX + rijY*rijY + rijZ*rijZ

      if ( rij2 .lt. cutoff2 ) then

        rij = sqrt(rij2)

        bin = aint(rij/dx)

        histNaNa(bin) = histNaNa(bin) + 2

      end if

    end do

  end do
!!!-----------------------------------------
  do j = nion+1, nion+hnion

      rXj = rX(j)
      rYj = rY(j)
      rZj = rZ(j)

    do i = nion+hnion+1, twonion

      rijX = rX(i) - rXj
      rijY = rY(i) - rYj
      rijZ = rZ(i) - rZj

      rijX = rijX - anint(rijX*invlX)*lX
      rijY = rijY - anint(rijY*invlY)*lY
      rijZ = rijZ - anint(rijZ*invlZ)*lZ

      rij2 = rijX*rijX + rijY*rijY + rijZ*rijZ

      if ( rij2 .lt. cutoff2 ) then

        rij = sqrt(rij2)

        bin = aint(rij/dx)

        histNaCl(bin) = histNaCl(bin) + 1

      end if

    end do

  end do
!!!---------------------------------------------
  do j = nion+hnion+1, twonion-1

      rXj = rX(j)
      rYj = rY(j)
      rZj = rZ(j)

    do i = j+1, twonion

      rijX = rX(i) - rXj
      rijY = rY(i) - rYj
      rijZ = rZ(i) - rZj

      rijX = rijX - anint(rijX*invlX)*lX
      rijY = rijY - anint(rijY*invlY)*lY
      rijZ = rijZ - anint(rijZ*invlZ)*lZ

      rij2 = rijX*rijX + rijY*rijY + rijZ*rijZ

      if ( rij2 .lt. cutoff2 ) then

        rij = sqrt(rij2)

        bin = aint(rij/dx)

        histClCl(bin) = histClCl(bin) + 2

      end if

    end do

  end do

return

end subroutine record_hist


subroutine output_RDF(hion, nrecord, pi, hrhoB)
  use RDFmod
  implicit none

  integer :: hion, nrecord
  real(8) :: pi, hrhoB

  integer :: bin
  real(8) :: rhead, rtail, cst1, cst2, dv
  real(8), dimension(0:maxbin) :: grNaNa, grNaCl, grClCl

  cst1 = real(nrecord,8)*real(hion-1,8)*4.0d0/3.0d0*pi*hrhoB
  cst2 = real(nrecord,8)*real(hion,8)*4.0d0/3.0d0*pi*hrhoB

  open(27, file = 'grNaNa.dat')
  open(28, file = 'grNaCl.dat')
  open(29, file = 'grClCl.dat')

  do bin = 0, maxbin-1

    rhead = real(bin,8)*dx
    rtail = rhead + dx

    dv = rtail**3 - rhead**3

    grNaNa(bin) = real(histNaNa(bin),8)/( cst1*dv )

    write(27,*) rhead + 0.50d0*dx, grNaNa(bin)

    grNaCl(bin) = real(histNaCl(bin),8)/( cst2*dv )

    write(28,*) rhead + 0.50d0*dx, grNaCl(bin)

    grClCl(bin) = real(histClCl(bin),8)/( cst1*dv )

    write(29,*) rhead + 0.50d0*dx, grClCl(bin)

  end do

  close(27)
  close(28)
  close(29)

return

end subroutine output_RDF


subroutine outputcoor(twonion, rX, rY, rZ)
  implicit none

  integer :: twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ

  integer :: j

  open(30,file = 'coor.xyz' )

    do j = 1, twonion

      write(30,*) rX(j), rY(j), rZ(j)

    end do

  close(30)

return

end subroutine outputcoor


subroutine outputv(nion, twonion, vX, vY, vZ, hvX, hvY, hvZ)
  implicit none

  integer :: nion, twonion
  real(8), dimension(nion+1:twonion) :: vX, vY, vZ, hvX, hvY, hvZ

  integer :: j

  open(34, file='v.dat')

  do j = nion+1, twonion

    write(34,*) vX(j), vY(j), vZ(j)

  end do

  close(34)

  open(35, file='hv.dat')

  do j = nion+1, twonion

    write(35,*) hvX(j), hvY(j), hvZ(j)

  end do

  close(35)

return

end subroutine outputv


subroutine outputforce(nion, FX, FY, FZ)
  implicit none

  integer :: nion
  real(8), dimension(1:nion) :: FX, FY, FZ

  integer :: j

  open(29,file = 'force_NaCl.dat')

    do j = 1, nion

      write(29,*) FX(j), FY(j), FZ(j)

    end do

  close(29)

return

end subroutine outputforce


subroutine outputforceRK(nion, FRX, FRY, FRZ, FKX, FKY, FKZ)
  implicit none

  integer :: nion
  real(8), dimension(1:nion) :: FRX, FRY, FRZ, FKX, FKY, FKZ

  integer :: j

  open(27,file = 'forceR_NaCl.dat')

    do j = 1, nion

      write(27,*) FRX(j), FRY(j), FRZ(j)

    end do

  close(27)

  open(28,file = 'forceK_NaCl.dat')

    do j = 1, nion

      write(28,*) FKX(j), FKY(j), FKZ(j)

    end do

  close(28)

return

end subroutine outputforceRK


subroutine output_finalcoor(nion, hnion, twonion, rX, rY, rZ, &
                    & lX, lY, lZ)
  implicit none

  integer :: nion, hnion, twonion
  real(8), dimension(1:twonion) :: rX, rY, rZ
  real(8) :: lX, lY, lZ

  integer :: j

  open(30,file = 'final_pos.xyz' )

      write(30,*) nion
      write(30,*) lX, lY, lZ

    do j = nion+1, nion+hnion

      write(30,*) 'Na', rX(j), rY(j), rZ(j)

    end do

    do j = nion+hnion+1, twonion

      write(30,*) 'Cl', rX(j), rY(j), rZ(j)

    end do

  close(30)

return

end subroutine output_finalcoor


