module Box
  implicit none
  
  integer, parameter :: natm = 12288
  real(8), parameter :: lX = 49.6685068176563650d0
  real(8), parameter :: lY = lX
  real(8), parameter :: lZ = lX

end module Box

module unit_conversion
  implicit none

  !!! 1.0 elementraY charge = 1.602176634d-19 coulomb
  real(8), parameter :: ele_to_clbp19 = 1.602176634d0 !!! C
  real(8), parameter :: epslonp12 = 8.8541878128d0 !!! C/(V*m)
  real(8), parameter :: nan23 = 6.02214076d0
  real(8), parameter :: kbp23 = 1.38064852 !!! Kg*m^2/s^2/k
end module unit_conversion

module setPME
  implicit none

  real(8), parameter :: alpha = 0.55d0
  integer, parameter :: nospl = 6
  integer, parameter :: kmax1 = 32
  integer, parameter :: kmax2 = 32
  integer, parameter :: kmax3 = 32

end module setPME

program main
    use Box
    use unit_conversion
    use setPME
    implicit none
  
    real(8), dimension(natm) :: q
    real(8), dimension(natm) :: rx, ry, rz
    real(8), dimension(natm) :: fx, fy, fz
    real(8) :: invlx, invly, invlz
  
    real(8) :: r4pie
    real(8) :: eng_io
    integer,dimension(1:kmax1):: key1
    integer,dimension(1:kmax2):: key2
    integer,dimension(1:kmax3):: key3
    complex(8),dimension(1:kmax1):: ww1
    complex(8),dimension(1:kmax2):: ww2
    complex(8),dimension(1:kmax3):: ww3
    complex(8),dimension(1:kmax1):: bscx
    complex(8),dimension(1:kmax2):: bscy
    complex(8),dimension(1:kmax3):: bscz
    real(8), dimension(1:kmax1,1:kmax2,1:kmax3) :: pentical_mesh
    real(8), parameter :: pi = acos(-1.0d0)
  
    integer, parameter :: KKmax1 = 20
    integer, parameter :: KKmax2 = 20
    integer, parameter :: KKmax3 = 20
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: Epartself
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: FXpartself
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: FYpartself
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: FZpartself
    integer :: kk1, kk2, kk3
    real(8) :: delx, dely, delz, dx, dy, dz, ee
    real(8), dimension(1) :: X, Y, Z, qq, FFX, FFY, FFZ

    integer :: i, j, k
    real(8) :: t1, t2, l
    character(1) :: ele

    invlx = 1.d0/lx
    invly = 1.d0/ly
    invlz = 1.d0/lz

    open(12,file="H2O.xyz")
        read(12,*)
        read(12,*)
        do i = 1, natm
            read(12,*) ele, rX(i), rY(i), rZ(i)
        end do
    close(12)

    do i = 1, natm, 3
        q(i) = -0.8476d0; q(i+1) = 0.4238d0; q(i+2) = 0.4238d0;
    end do

    r4pie = ele_to_clbp19**2*nan23/(4.0d0*pi*epslonp12)*1.0d6 !!! unit: 10J/mol/ang

    call fft_efb_3(kmax1,kmax2,kmax3,key1,key2,key3,ww1,ww2,ww3)

    call bsp_c_3(nospl,kmax1,kmax2,kmax3,ww1,ww2,ww3,bscx,bscy,bscz)

    call calc_pentical_mesh(invlX, invlY, invlZ, key1,key2,key3,ww1,ww2,ww3,&
                            pentical_mesh ) 

!!!!!!!!!!!!
    delx  = lX/real(Kmax1,8)
    dely  = lY/real(Kmax2,8)
    delz  = lZ/real(Kmax3,8)

    dx  = lX/real(Kmax1*KKmax1,8)
    dy  = lY/real(Kmax2*KKmax2,8)
    dz  = lZ/real(Kmax3*KKmax3,8)

    qq(1) = 1.d0
    do kk1 = 0, KKmax1
        X(1) = real(kk1,8)*dx
    do kk2 = 0, KKmax2
        Y(1) = real(kk2,8)*dy
    do kk3 = 0, KKmax3
        Z(1) = real(kk3,8)*dz

        call fr_pmel_p3(1, invlx, invly, invlz, r4pie, X, Y, Z, qq, &
               key1,key2,key3,ww1,ww2,ww3, bscx,bscy,bscz, pentical_mesh,&
               eng_io, FFX, FFY, FFZ )
        Epartself(kk1, kk2, kk3) = eng_io
        FXpartself(kk1, kk2, kk3) = FFX(1)
        FYpartself(kk1, kk2, kk3) = FFY(1)
        FZpartself(kk1, kk2, kk3) = FFZ(1)
    end do
    end do
    end do
!!!!!!!!!!!!

    call fr_pmel_p3(natm, invlx, invly, invlz, r4pie, rx, ry, rz, q, &
         key1,key2,key3,ww1,ww2,ww3, bscx,bscy,bscz, pentical_mesh,&
         eng_io, fx, fy, fz )
 
    call deleta_partself( natm, invlX, invlY, invlZ, rX, rY, rZ, q, KKmax1, KKmax2, KKmax3,&
                        & delx, dely, delz, dx, dy, dz, Epartself, eng_io,&
                        & FXpartself, FYpartself, FZpartself, fx, fy, fz)

    open(34,file="Energy.dat")
        write(34,*) eng_io, "10J/mol"
    close(34)

    open(33,file="Force.dat")
         write(33,*) "10J/mol/angstron"
         do i = 1, natm
             write(33,*) fx(i), fy(i), fz(i)
         end do
    close(33)
    

end program main


subroutine deleta_partself( nion, invlX, invlY, invlZ, rX, rY, rZ, q, KKmax1, KKmax2, KKmax3,&
                          & delx, dely, delz, dx, dy, dz, Epartself, recE, &
                          & FXpartself, FYpartself, FZpartself, fx, fy, fz )
    implicit none

    integer :: nion
    integer :: KKmax1, KKmax2, KKmax3
    real(8) :: invlX, invlY, invlZ
    real(8), dimension(1:nion) :: rX, rY, rZ, q
    real(8), dimension(1:nion) :: fX, fY, fZ
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: Epartself
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: FXpartself
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: FYpartself
    real(8), dimension(0:KKmax1,0:KKmax2,0:KKmax3) :: FZpartself
    real(8) :: delx, dely, delz, dx, dy, dz
    real(8) :: recE

    integer :: i
    integer :: kk1, kk2, kk3
    real(8) :: X, Y, Z

    do i = 1, nion

        X = rX(i) - anint(rX(i)*invlX)/invlX + 0.5d0/invlX
        Y = rY(i) - anint(rY(i)*invlY)/invlY + 0.5d0/invlY
        Z = rZ(i) - anint(rZ(i)*invlZ)/invlZ + 0.5d0/invlZ

        X = X - floor(X/delX)*delX
        Y = Y - floor(Y/delY)*delY
        Z = Z - floor(Z/delZ)*delZ

        kk1 =  floor(X/dx)
        kk2 =  floor(Y/dy)
        kk3 =  floor(Z/dz)

        recE = recE - q(i)**2*Epartself(kk1, kk2, kk3)
        fX(i)= fX(i) - q(i)**2*FXpartself(kk1, kk2, kk3)
        fY(i)= fY(i) - q(i)**2*FYpartself(kk1, kk2, kk3)
        fZ(i)= fZ(i) - q(i)**2*FZpartself(kk1, kk2, kk3) 

    end do

    return;
end subroutine deleta_partself



subroutine fft_efb_3(ndiv1,ndiv2,ndiv3,key1,key2,key3,ww1,ww2,ww3)
  implicit none
  integer:: ndiv1,ndiv2,ndiv3

  integer,dimension(1:ndiv1):: key1
  integer,dimension(1:ndiv2):: key2
  integer,dimension(1:ndiv3):: key3
  complex(8),dimension(1:ndiv1):: ww1
  complex(8),dimension(1:ndiv2):: ww2
  complex(8),dimension(1:ndiv3):: ww3
  !------------------------------------------------------------------------------
  !
  ! fft_efb_3():
  !
  !  fast fourier transfrom complex exponential factors and bit address in 3 dim
  !  f    f       t_                e           f           b
  !
  !------------------------------------------------------------------------------
  real(8),parameter:: tpi = 2.d0*acos(-1.0d0)
  integer::  nu1,nu2,nu3,iii,jjj,i,j,jj2, kkk
  logical:: lgndiv
  real(8):: arg
  !
  ! check array dimensions
  nu1 = int(log(dble(ndiv1)+1.0d-10)/log(2.0d0))
  nu2 = int(log(dble(ndiv2)+1.0d-10)/log(2.0d0))
  nu3 = int(log(dble(ndiv3)+1.0d-10)/log(2.0d0))
  lgndiv = (ndiv1.ne.2**nu1) .or.  (ndiv2.ne.2**nu2) .or. (ndiv3.ne.2**nu3)
  if (lgndiv) stop 'fft array not 2**N; stop fft_efb_3'
  !
  ! set reverse bit address arrays
  do kkk=1,ndiv1
    iii=0; jjj=kkk-1
    do j=1,nu1
      jj2=jjj/2; iii=2*(iii-jj2)+jjj; jjj=jj2
    end do
    key1(kkk)=iii+1
  end do
  do kkk=1,ndiv2
    iii=0; jjj=kkk-1
    do j=1,nu2
      jj2=jjj/2; iii=2*(iii-jj2)+jjj; jjj=jj2
    end do
    key2(kkk)=iii+1
  end do
  do kkk=1,ndiv3
    iii=0; jjj=kkk-1
    do j=1,nu3
      jj2=jjj/2; iii=2*(iii-jj2)+jjj; jjj=jj2
    end do
    key3(kkk)=iii+1
  end do
  !
  ! initialise complex exponential factors
  ww1(1)=(1.d0,0.d0); ww1(ndiv1/2+1)=(-1.d0,0.d0)
  do i=1,ndiv1/2-1
    arg=(tpi/dble(ndiv1))*dble(i)
    ww1(i+1)=cmplx(cos(arg),sin(arg),kind=8); ww1(ndiv1+1-i)=conjg(ww1(i+1))
  end do
  ww2(1)=(1.d0,0.d0); ww2(ndiv2/2+1)=(-1.d0,0.d0)
  do i=1,ndiv2/2-1
    arg=(tpi/dble(ndiv2))*dble(i)
    ww2(i+1)=cmplx(cos(arg),sin(arg),kind=8); ww2(ndiv2+1-i)=conjg(ww2(i+1))
  end do
  ww3(1)=(1.d0,0.d0); ww3(ndiv3/2+1)=(-1.d0,0.d0)
  do i=1,ndiv3/2-1
    arg=(tpi/dble(ndiv3))*dble(i)
    ww3(i+1)=cmplx(cos(arg),sin(arg),kind=8); ww3(ndiv3+1-i)=conjg(ww3(i+1))
  end do
  return
end subroutine fft_efb_3

subroutine bsp_c_3(nospl,kmax1,kmax2,kmax3,ww1,ww2,ww3,bscx,bscy,bscz)
  implicit none
  integer:: nospl,kmax1,kmax2,kmax3
  complex(8),dimension(1:kmax1):: ww1
  complex(8),dimension(1:kmax2):: ww2
  complex(8),dimension(1:kmax3):: ww3

  complex(8),dimension(1:kmax1):: bscx
  complex(8),dimension(1:kmax2):: bscy
  complex(8),dimension(1:kmax3):: bscz
  !------------------------------------------------------------------------------
  !  bsp_c_3()
  !
  !  b-spline coefficient for Euler exponential splines in 3 dim
  !  b sp_    c_                                           3
  !
  !  see Essmann et.al. A smooth particle mesh Ewald method, 1995  JCP, 8577
  !
  !  see Eq. (4.1) and Eq. (4.4) in  page 8581
  !
  !  second part see: bsp_d_3(nospl,natm,tx,ty,tz,bspx,bspy,bspz,bsdx,bsdy,bsdz)
  !
  !------------------------------------------------------------------------------
  real(8),dimension(1:nospl):: csp
  integer:: i,j,k
  complex(8):: ccc
  !
  ! calculate B-splines at knots
  !
  ! This part computes    M_n( u )  =  csp( u-1 )  using Eq. (4.1)
  !
  ! for integer u = 0, 1, 2, ..., n-1
  !
  ! where  n = nospl
  !
  ! Note that to keep value of csp at the order k-1:  M_{k-1} (u-1)
  ! we must carry out the computation in a reverse order for j = k, 2, -1
  ! otherwise the value of csp(u-1) at the old order k-1 is replaced by
  ! its value at the new order k
  !
  csp(1)=0.d0; csp(2)=1.d0
  do k=3,nospl
    csp(k)=0.d0
    do j=k,2,-1
      csp(j)=(dble(j-1)*csp(j)+dble(k-j+1)*csp(j-1))/dble(k-1)
    end do
  end do
  !
  ! calculate B-spline coefficients
  !
  ! see Eq. (4.4)
  !
  ! bscx(i+1) = b_1( i )  of  Eq. (4.4)
  ! bscy(i+1) = b_2( i )  of  Eq. (4.4)
  ! bscz(i+1) = b_3( i )  of  Eq. (4.4)
  !
  ! ww1(), ww2(), ww3() are complex exponential factors which
  !                     were defined in the first part of dlpfft3()
  !
  ! csp(k+2) = M_n(k+1)
  !
  do i=0,kmax1-1
    ccc=(0.d0,0.d0)
    do k=0,nospl-2
      ccc=ccc+csp(k+2)*ww1(mod(i*k,kmax1)+1)
    end do
    bscx(i+1)=ww1(mod(i*(nospl-1),kmax1)+1)/ccc
  end do
  do i=0,kmax2-1
    ccc=(0.d0,0.d0)
    do k=0,nospl-2
      ccc=ccc+csp(k+2)*ww2(mod(i*k,kmax2)+1)
    end do
    bscy(i+1)=ww2(mod(i*(nospl-1),kmax2)+1)/ccc
  end do
  do i=0,kmax3-1
    ccc=(0.d0,0.d0)
    do k=0,nospl-2
      ccc=ccc+csp(k+2)*ww3(mod(i*k,kmax3)+1)
    end do
    bscz(i+1)=ww3(mod(i*(nospl-1),kmax3)+1)/ccc
  end do
  return
end subroutine bsp_c_3

subroutine  calc_pentical_mesh(invlX, invlY, invlZ, key1,key2,key3,ww1,ww2,ww3, pentical_mesh )
    use setPME
    implicit none

    real(8) :: invlx, invly, invlz
    integer,dimension(1:kmax1):: key1
    integer,dimension(1:kmax2):: key2
    integer,dimension(1:kmax3):: key3
    complex(8),dimension(1:kmax1):: ww1
    complex(8),dimension(1:kmax2):: ww2
    complex(8),dimension(1:kmax3):: ww3   
    real(8), dimension(1:kmax1,1:kmax2,1:kmax3) :: pentical_mesh
    complex(8), dimension(1:kmax1,1:kmax2,1:kmax3) :: ZZU

    integer, parameter :: pkmax1 = 128
    integer, parameter :: pkmax2 = 128
    integer, parameter :: pkmax3 = 128
    integer,dimension(1:pkmax1) :: pkey1
    integer,dimension(1:pkmax2) :: pkey2
    integer,dimension(1:pkmax3) :: pkey3
    complex(8),dimension(1:pkmax1) :: pww1
    complex(8),dimension(1:pkmax2) :: pww2
    complex(8),dimension(1:pkmax3) :: pww3
    complex(8), dimension(1:pkmax1,1:pkmax2,1:pkmax3) :: SDU
    real(8),parameter:: tpi = 2.d0*acos(-1.0d0)
    real(8), parameter :: r4alpsq = -0.25d0/(alpha*alpha)

    integer :: j, k, l, jj, kk, ll
    real(8) :: rk3sq, rk2sq, rksq
    integer :: multi1, multi2, multi3


    call fft_efb_3(pkmax1,pkmax2,pkmax3,pkey1,pkey2,pkey3,pww1,pww2,pww3) 

    do l=1, pkmax3
        ll=l-1; if(l.gt.pkmax3/2)ll=l-pkmax3-1; 
        rk3sq = (tpi*dble(ll)*invlz)**2

        do k=1, pkmax2
            kk=k-1; if(k.gt.pkmax2/2)kk=k-pkmax2-1; 
            rk2sq = rk3sq + (tpi*dble(kk)*invly)**2

            do j=1, pkmax1
                jj=j-1; if(j.gt.pkmax1/2)jj=j-pkmax1-1; 

                rksq = rk2sq + (tpi*dble(jj)*invlx)**2
 
                if (rksq.gt.1.d-6) then
                    SDU(j,k,l)= complex(exp(r4alpsq*rksq)/rksq, 0.d0)
                else
                    SDU(j,k,l)=(0.d0,0.d0)
                end if
            end do
        end do
    end do

    call fft_acx_lg3(pkmax1,pkmax2,pkmax3,pkey1,pkey2,pkey3,pww1,pww2,pww3,SDU,0)

    multi1 = pkmax1/kmax1; multi2 = pkmax2/kmax2; multi3 = pkmax3/kmax3;

    do l = 1, kmax3
        ll = multi3*(l-1)+1
    do k = 1, kmax2
        kk = multi2*(k-1)+1
    do j = 1, kmax1
        jj = multi1*(j-1)+1
         
        ZZU(j,k,l) = SDU(jj,kk,ll)
    end do
    end do
    end do

    call fft_acx_lg3(kmax1,kmax2,kmax3,key1,key2,key3,ww1,ww2,ww3,ZZU,1)


    do l = 1, kmax3
    do k = 1, kmax2
    do j = 1, kmax1
        pentical_mesh(j,k,l) = real(ZZU(j,k,l))
    end do
    end do
    end do

    return;
end subroutine  calc_pentical_mesh

subroutine fr_pmel_p3(natm, invlx, invly, invlz, r4pie, rx, ry, rz, q,&
           key1,key2,key3,ww1,ww2,ww3,bscx,bscy,bscz, pentical_mesh, &
           eng_io, fx, fy, fz )
    use setPME
    implicit none

    integer :: natm
    real(8) :: invlx, invly, invlz
    real(8) :: r4pie
    real(8), dimension(natm) :: q
    real(8), dimension(natm) :: rx, ry, rz
    integer,dimension(1:kmax1):: key1
    integer,dimension(1:kmax2):: key2
    integer,dimension(1:kmax3):: key3
    complex(8),dimension(1:kmax1):: ww1
    complex(8),dimension(1:kmax2):: ww2
    complex(8),dimension(1:kmax3):: ww3
    complex(8),dimension(1:kmax1):: bscx
    complex(8),dimension(1:kmax2):: bscy
    complex(8),dimension(1:kmax3):: bscz
    real(8), dimension(1:kmax1,1:kmax2,1:kmax3) :: pentical_mesh
    real(8) :: eng_io, eng_non_neu
    real(8), dimension(natm) :: fx, fy, fz
  
    real(8), dimension(natm) :: tx, ty, tz
    real(8), dimension(1:natm,1:nospl) :: bspx,bspy,bspz
    real(8), dimension(1:natm,1:nospl) :: bsdx,bsdy,bsdz
    real(8), dimension(1:kmax1,1:kmax2,1:kmax3) :: qqc;
    complex(8), dimension(1:kmax1,1:kmax2,1:kmax3) :: qqq
   
    integer :: i, j, k, l, jj, kk, ll
    real(8) :: tmp
    real(8),parameter:: tpi = 2.d0*acos(-1.0d0)
    real(8):: bb1,bb2,bb3, r4alpsq,engcoeffi
    real(8):: rcpcut,rcpct2, rkx,rky,rkz,rksq
    complex(8):: cpetot
    real(8):: fcoeffi,qsum,bdx,bdy,bdz
    real(8),dimension(1:3):: fff
    real(8) :: qtotsq
  
    qtotsq = ( sum(q) )**2
  
    eng_io = 0.d0
    fx = 0.d0; fy = 0.d0; fz = 0.d0
  
    do i=1, natm
      tmp=rx(i) - anint(rx(i)*invlx)/invlX; tx(i)=dble(kmax1)*(invlx*tmp+0.5d0)
      tmp=ry(i) - anint(ry(i)*invly)/invlY; ty(i)=dble(kmax2)*(invly*tmp+0.5d0)
      tmp=rz(i) - anint(rz(i)*invlz)/invlZ; tz(i)=dble(kmax3)*(invlz*tmp+0.5d0)
    end do
 
  call bsp_d_3(natm,tx,ty,tz,bspx,bspy,bspz,bsdx,bsdy,bsdz)


  qqc(:,:,:) = 0.0d0
 
  do i=1, natm
    do l=1, nospl
      ll=int(tz(i))-l+2; if(ll.gt.kmax3)ll=1; if(ll.lt.1)ll=ll+kmax3
      do k=1, nospl
        kk=int(ty(i))-k+2; if(kk.gt.kmax2)kk=1; if(kk.lt.1)kk=kk+kmax2
        do j=1, nospl
          jj=int(tx(i))-j+2; if(jj.gt.kmax1)jj=1; if(jj.lt.1)jj=jj+kmax1
          qqc(jj,kk,ll)=qqc(jj,kk,ll)+ q(i)*bspx(i,j)*bspy(i,k)*bspz(i,l)
        end do
      end do
    end do

  end do

  qqq(:,:,:) = cmplx(qqc(:,:,:),0.0d0,kind=8)


  call fft_acx_lg3(kmax1,kmax2,kmax3,key1,key2,key3,ww1,ww2,ww3,qqq,0)


  rcpcut=0.5d0*min(dble(kmax1)*invlx,dble(kmax2)*invly,dble(kmax3)*invlz)
  rcpcut=rcpcut*1.05d0*tpi;  rcpct2=rcpcut**2;  r4alpsq=-0.25d0/alpha**2

  do l=1, kmax3
    bb3=real(bscz(l)*conjg(bscz(l)))
     ll=l-1; if(l.gt.kmax3/2)ll=l-kmax3-1; rkz=tpi*dble(ll)*invlz
    do k=1, kmax2
      bb2=bb3*real(bscy(k)*conjg(bscy(k)))
      kk=k-1; if(k.gt.kmax2/2)kk=k-kmax2-1; rky=tpi*dble(kk)*invly
      do j=1, kmax1
        bb1=bb2*real(bscx(j)*conjg(bscx(j)))    ! see Eq. (4.8), page 8581
        jj=j-1; if(j.gt.kmax1/2)jj=j-kmax1-1; rkx=tpi*dble(jj)*invlx
        rksq=rkx*rkx+rky*rky+rkz*rkz

          qqq(j,k,l)=bb1*pentical_mesh(j,k,l)*qqq(j,k,l)
      end do
    end do
  end do


  call fft_acx_lg3(kmax1,kmax2,kmax3,key1,key2,key3,ww1,ww2,ww3,qqq,1)

  cpetot = sum(qqq(:,:,:)*qqc(:,:,:))
  engcoeffi = tpi*(invlx*invly*invlz)*r4pie/real(Kmax1*Kmax2*Kmax3, 8); 
  eng_io = real(cpetot)*engcoeffi

  cpetot = sum(pentical_mesh(:,:,:))*qtotsq
  eng_non_neu = real(cpetot)*engcoeffi

  eng_io = eng_io - eng_non_neu

  ! STEP 8): calculate forces
  fcoeffi=-2.0d0*engcoeffi
  do i=1, natm
    do l=1, nospl
      ll=int(tz(i))-l+2; if(ll.gt.kmax3)ll=1; if(ll.lt.1)ll=ll+kmax3
      do k=1, nospl
        kk=int(ty(i))-k+2; if(kk.gt.kmax2)kk=1; if(kk.lt.1)kk=kk+kmax2
        do j=1, nospl
          jj=int(tx(i))-j+2; if(jj.gt.kmax1)jj=1; if(jj.lt.1)jj=jj+kmax1
          qsum=real(qqq(jj,kk,ll))
          bdx=qsum*bsdx(i,j)*bspy(i,k)*bspz(i,l)*dble(kmax1)
          bdy=qsum*bspx(i,j)*bsdy(i,k)*bspz(i,l)*dble(kmax2)
          bdz=qsum*bspx(i,j)*bspy(i,k)*bsdz(i,l)*dble(kmax3)
          fx(i)=fx(i)+fcoeffi*q(i)*(bdx*invlx)
          fy(i)=fy(i)+fcoeffi*q(i)*(bdy*invly)
          fz(i)=fz(i)+fcoeffi*q(i)*(bdz*invlz)
        end do
      end do
    end do
  end do
  !
  !   remove COM drift arising from SPME approximations

!    if ( natm .ne. 1 ) then 
!        fff(1)=0.d0; fff(2)=0.d0; fff(3)=0.d0
!        do i=1, natm
!          fff(1)=fff(1)+fx(i);  fff(2)=fff(2)+fy(i);  fff(3)=fff(3)+fz(i)
!        end do
!        fff(1)=fff(1)/dble(natm); fff(2)=fff(2)/dble(natm); fff(3)=fff(3)/dble(natm)
!        do i=1, natm
!          fx(i)=fx(i)-fff(1); fy(i)=fy(i)-fff(2); fz(i)=fz(i)-fff(3)
!        end do
!    end if

    return

end subroutine fr_pmel_p3

subroutine fft_acx_lg3(ndiv1,ndiv2,ndiv3,key1,key2,key3,ww1,ww2,ww3,f_io,lgi_back)
  implicit none
  integer:: ndiv1,ndiv2,ndiv3
  integer,dimension(1:ndiv1):: key1
  integer,dimension(1:ndiv2):: key2
  integer,dimension(1:ndiv3):: key3
  complex(8),dimension(1:ndiv1):: ww1
  complex(8),dimension(1:ndiv2):: ww2
  complex(8),dimension(1:ndiv3):: ww3
  integer:: lgi_back

  complex(8),dimension(1:ndiv1,1:ndiv2,1:ndiv3):: f_io

  !------------------------------------------------------------------------------
  !
  ! fft_efb_3():
  !
  ! fast fourier transfrom array of complex with logic forward/backward in 3 dim
  ! f    f       t_        a        c     x_     l g                       3
  !
  ! f_io:  input/output array of complex
  ! lgi_back: 0 for forward, 1 for backward without dividing N=ndiv1*ndiv2*ndiv3
  !
  ! key1(1:ndiv1), key2(1:ndiv2), key3(1:ndiv3): input bit address
  ! ww1(1:ndiv1),ww2(1:ndiv2),ww3(1:ndiv3): input exponential factors
  !
  ! this algorithm is supposed to be O( N * log_2(N) ) complexity
  !------------------------------------------------------------------------------
  integer:: nu1,nu2,nu3
  integer:: i,j,k,l, iii,jjj,kkk, num,kk1,k12
  complex(8):: cdy ! complex dummy variable
  logical:: lgndiv
  !
  ! check array dimensions
  nu1 = int(log(dble(ndiv1)+1.0d-10)/log(2.0d0))
  nu2 = int(log(dble(ndiv2)+1.0d-10)/log(2.0d0))
  nu3 = int(log(dble(ndiv3)+1.0d-10)/log(2.0d0))
  lgndiv = (ndiv1.ne.2**nu1) .or.  (ndiv2.ne.2**nu2) .or. (ndiv3.ne.2**nu3)
  if (lgndiv) stop 'fft array not 2**N; stop fft_efb_3'
  !
  ! take conjugate of exponentials if required
  if(lgi_back.eq.1) then
    do i=1,ndiv1
      ww1(i)=conjg(ww1(i))
    end do
    do i=1,ndiv2
      ww2(i)=conjg(ww2(i))
    end do
    do i=1,ndiv3
      ww3(i)=conjg(ww3(i))
    end do
  end if
  !
  ! perform fourier transform in X direction
  kkk=0; num=ndiv1/2
  do l=1,nu1
    do while(kkk.lt.ndiv1)
      do i=1,num
        iii=key1(kkk/num+1); kk1=kkk+1; k12=kk1+num
        do j=1,ndiv2
          do k=1,ndiv3
            cdy=f_io(k12,j,k)*ww1(iii)
            f_io(k12,j,k)=f_io(kk1,j,k)-cdy; f_io(kk1,j,k)=f_io(kk1,j,k)+cdy
          end do
        end do
        kkk=kkk+1
      end do
      kkk=kkk+num
    end do
    kkk=0; num=num/2
  end do
  !
  ! unscramble the fft using bit address array
  do kkk=1,ndiv1
    iii=key1(kkk)
    if(iii.gt.kkk)then
      do j=1,ndiv2
        do k=1,ndiv3
          cdy=f_io(kkk,j,k)
          f_io(kkk,j,k)=f_io(iii,j,k); f_io(iii,j,k)=cdy
        end do
      end do
    end if
  end do
  !
  ! perform fourier transform in Y direction
  kkk=0; num=ndiv2/2
  do l=1,nu2
    do while(kkk.lt.ndiv2)
      do i=1,num
        iii=key2(kkk/num+1); kk1=kkk+1; k12=kk1+num
        do j=1,ndiv1
          do k=1,ndiv3
            cdy=f_io(j,k12,k)*ww2(iii)
            f_io(j,k12,k)=f_io(j,kk1,k)-cdy; f_io(j,kk1,k)=f_io(j,kk1,k)+cdy
          end do
        end do
        kkk=kkk+1
      end do
      kkk=kkk+num
    end do
    kkk=0; num=num/2
  end do
  !
  ! unscramble the fft using bit address array
  do kkk=1,ndiv2
    iii=key2(kkk)
    if(iii.gt.kkk)then
      do j=1,ndiv1
        do k=1,ndiv3
          cdy=f_io(j,kkk,k)
          f_io(j,kkk,k)=f_io(j,iii,k); f_io(j,iii,k)=cdy
        end do
      end do
    end if
  end do
  !
  ! perform fourier transform in Z direction
  kkk=0; num=ndiv3/2
  do l=1,nu3
    do while(kkk.lt.ndiv3)
      do i=1,num
        iii=key3(kkk/num+1); kk1=kkk+1; k12=kk1+num
        do j=1,ndiv1
          do k=1,ndiv2
            cdy=f_io(j,k,k12)*ww3(iii)
            f_io(j,k,k12)=f_io(j,k,kk1)-cdy; f_io(j,k,kk1)=f_io(j,k,kk1)+cdy
          end do
        end do
        kkk=kkk+1
      end do
      kkk=kkk+num
    end do
    kkk=0; num=num/2
  end do
  !
  ! unscramble the fft using bit address array
  do kkk=1,ndiv3
    iii=key3(kkk)
    if(iii.gt.kkk)then
      do j=1,ndiv1
        do k=1,ndiv2
          cdy=f_io(j,k,kkk)
          f_io(j,k,kkk)=f_io(j,k,iii); f_io(j,k,iii)=cdy
        end do
      end do
    end if
  end do
  !
  ! restore exponentials to unconjugated values if necessary
  if(lgi_back.eq.1)then
    do i=1,ndiv1
      ww1(i)=conjg(ww1(i))
    end do
    do i=1,ndiv2
      ww2(i)=conjg(ww2(i))
    end do
    do i=1,ndiv3
      ww3(i)=conjg(ww3(i))
    end do
  end if
  return
end subroutine fft_acx_lg3


subroutine bsp_d_3(natm,tx,ty,tz,bspx,bspy,bspz,bsdx,bsdy,bsdz)
  use setPME
  implicit none

  integer :: natm
  real(8), dimension(natm) :: tx, ty, tz
  real(8), dimension(1:natm,1:nospl) :: bspx,bspy,bspz
  real(8), dimension(1:natm,1:nospl) :: bsdx,bsdy,bsdz

  integer :: i, j, k
  real(8) :: aaa, bbb, ccc

  !  bsdx(iatm,n) == dM_n(u_i)
  !  bspx(iatm,n) == M_n(u_i)

  do i=1,natm
    bsdx(i,1)= 1.d0; bsdy(i,1)= 1.d0; bsdz(i,1)= 1.d0 ! dM_1(u_i) = 1
    bsdx(i,2)=-1.d0; bsdy(i,2)=-1.d0; bsdz(i,2)=-1.d0 ! dM_2(u_i) = 1
    bspx(i,1)=tx(i)-int(tx(i)) ! M_1(u_i) = ux(i) - ux^int(i)
    bspy(i,1)=ty(i)-int(ty(i))  
    bspz(i,1)=tz(i)-int(tz(i))
    bspx(i,2)=1.d0-tx(i)+int(tx(i)) ! M_2(u_i) =  1 - (ux(i) - ux^int(i))
    bspy(i,2)=1.d0-ty(i)+int(ty(i))
    bspz(i,2)=1.d0-tz(i)+int(tz(i))
  end do

  do k=3,nospl

    do i=1,natm  !  M_k(u_i) = 0
      bspx(i,k)=0.d0; bspy(i,k)=0.d0; bspz(i,k)=0.d0 !  
    end do

    do j=k,2,-1

      if(k.eq.nospl)then
        do i=1,natm ! ls
          bsdx(i,j)=bspx(i,j)-bspx(i,j-1) 
          bsdy(i,j)=bspy(i,j)-bspy(i,j-1)
          bsdz(i,j)=bspz(i,j)-bspz(i,j-1)
        end do
      end if
                   !              uu             k-uu
      do i=1,natm  ! M_j(u_i) = ---- M_j((u_i) + ----- M_{j-1}(u_i)
                   !             k-1             k-1
        aaa=tx(i)+dble(j-1)-int(tx(i))
        bbb=ty(i)+dble(j-1)-int(ty(i))
        ccc=tz(i)+dble(j-1)-int(tz(i))
        bspx(i,j)=(aaa*bspx(i,j)+(dble(k)-aaa)*bspx(i,j-1))/dble(k-1)
        bspy(i,j)=(bbb*bspy(i,j)+(dble(k)-bbb)*bspy(i,j-1))/dble(k-1)
        bspz(i,j)=(ccc*bspz(i,j)+(dble(k)-ccc)*bspz(i,j-1))/dble(k-1)
      end do

    end do

    if(k.eq.nospl)then
      do i=1,natm
        bsdx(i,1)=bspx(i,1); bsdy(i,1)=bspy(i,1); bsdz(i,1)=bspz(i,1)
      end do
    end if

    do i=1,natm
      bspx(i,1)=(tx(i)-int(tx(i)))*bspx(i,1)/dble(k-1)
      bspy(i,1)=(ty(i)-int(ty(i)))*bspy(i,1)/dble(k-1)
      bspz(i,1)=(tz(i)-int(tz(i)))*bspz(i,1)/dble(k-1)
    end do
     
  end do

  
return

end subroutine bsp_d_3


 





