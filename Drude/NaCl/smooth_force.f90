program test_force
  implicit none
  
  real(8), parameter :: cutoff = 10.0d0
  real(8), parameter :: kNa = 63.0140d0, kCl = 25.7240d0
  real(8), parameter :: ANaNa = 487.0d0, ANaCl = 145134.0d0, AClCl = 405774.0d0
  real(8), parameter :: rho = 0.237680d0
  real(8), parameter :: CNaNa = 1.050d0, CNaCl = 6.990d0, CClCl = 72.40d0
  real(8), parameter :: DNaNa = 0.50d0, DNaCl = 8.70d0, DClCl = 145.40d0
 
  real(8), parameter :: cgNa = -0.50560d0, cgCl = -2.50050d0
  real(8) :: zNa, zCl
  real(8), parameter :: alpha = 1.0d0/4.50d0 !!! angstrom^-1
  real(8), parameter :: pi = acos(-1.0d0)

  integer :: i
  real(8), parameter :: dx = 1.0d-3
  real(8) :: r4pie, eang
  real(8) :: A, C, D, c0, c3, c4
  real(8) :: r, r2, x, u, f
  
  open(23, file='xi_shu.dat')

  write(23,*) 'cutoff=', cutoff
  write(23,*) 'alpha=', alpha

  r4pie = 138935.457644382060d0
  eang  = 9648.53321233100360d0 

  r = cutoff
  r2 = cutoff*cutoff

  A = ANaNa; C = CNaNa; D = DNaNa

  c3 = 1.0d0/r2*( (A/rho**2*exp(-r/rho) - 42.0d0*C/r2**4 - 72.0d0*D/r2**5)*r/3.0d0 +&
               & -(-A/rho*exp(-r/rho) + 6.0d0*C/(r2**3*r) + 8.0d0*D/(r2**4*r) )   )

  c4 =-1.0d0/(2.0d0*r2*r)*( (A/rho**2*exp(-r/rho) - 42.0d0*C/r2**4 - 72.0d0*D/r2**5)*r/2.0d0 +&
               & -(-A/rho*exp(-r/rho) + 6.0d0*C/(r2**3*r) + 8.0d0*D/(r2**4*r) )   )

  c0 =-( A*exp(-r/rho) - C/r2**3 - D/r2**4 ) - c3*r**3 - c4*r**4 

  write(23,*) 'NaNa'
  write(23,*) 'c3=',  c3, 'c4=', c4
  write(23,*)

  A = ANaCl; C = CNaCl; D = DNaCl

  c3 = 1.0d0/r2*( (A/rho**2*exp(-r/rho) - 42.0d0*C/r2**4 - 72.0d0*D/r2**5)*r/3.0d0 +&
               & -(-A/rho*exp(-r/rho) + 6.0d0*C/(r2**3*r) + 8.0d0*D/(r2**4*r) )   )

  c4 =-1.0d0/(2.0d0*r2*r)*( (A/rho**2*exp(-r/rho) - 42.0d0*C/r2**4 - 72.0d0*D/r2**5)*r/2.0d0 +&
               & -(-A/rho*exp(-r/rho) + 6.0d0*C/(r2**3*r) + 8.0d0*D/(r2**4*r) )   )

  write(23,*) 'NaCl'
  write(23,*) 'c3=',  c3, 'c4=', c4
  write(23,*)

  A = AClCl; C = CClCl; D = DClCl

  c3 = 1.0d0/r2*( (A/rho**2*exp(-r/rho) - 42.0d0*C/r2**4 - 72.0d0*D/r2**5)*r/3.0d0 +&
               & -(-A/rho*exp(-r/rho) + 6.0d0*C/(r2**3*r) + 8.0d0*D/(r2**4*r) )   )

  c4 =-1.0d0/(2.0d0*r2*r)*( (A/rho**2*exp(-r/rho) - 42.0d0*C/r2**4 - 72.0d0*D/r2**5)*r/2.0d0 +&
               & -(-A/rho*exp(-r/rho) + 6.0d0*C/(r2**3*r) + 8.0d0*D/(r2**4*r) )   )

  write(23,*) 'ClCl'
  write(23,*) 'c3=',  c3, 'c4=', c4
  write(23,*)


  c3 = 1.0d0/r2*( &
(4.0d0*alpha**3/sqrt(pi)*exp(-alpha**2*r2) + 4.0d0*alpha/sqrt(pi)*exp(-alpha**2*r2)/r2 &
 + 2.0d0*erfc(alpha*r)/(r2*r) )*r/3.0d0      -              &
( (-2.0d0*alpha/sqrt(pi)*exp(-alpha**2*r2)*r - erfc(alpha*r))/r2 )  )

  c4 = -1.0d0/(2.0d0*r2*r)*( &
(4.0d0*alpha**3/sqrt(pi)*exp(-alpha**2*r2) + 4.0d0*alpha/sqrt(pi)*exp(-alpha**2*r2)/r2 &
 + 2.0d0*erfc(alpha*r)/(r2*r) )*r/2.0d0 - &
( (-2.0d0*alpha/sqrt(pi)*exp(-alpha**2*r2)*r - erfc(alpha*r))/r2 )         )

  write(23,*) 'jingdian'
  write(23,*) 'c3=', c3, 'c4=', c4

   
return

end  
