! This is where the action lives
! Units are cgs, except mass
! Follow the recipe: https://arxiv.org/abs/2111.06895


subroutine spawn(x,y,z,vx,vy,vz)
  use star
  implicit none
  double precision :: x,y,z,vx,vy,vz,T
  logical accept
  double precision :: r, N,pdf,a,b,ctheta,phi
  double precision random_normal

  ! print*,"Spawining a walker"

  !Basic rejection method to get r
  accept = .false.
  do while (accept .eqv. .false.)
  call random_number(a)
  a = a*Rsun
  pdf = a**2*exp(-(a/rchi)**2)
  call random_number(b)
  b= b*rchi**2*exp(-1.d0)
  ! print*," Rsun ", rsun, " rchi ", rchi,  " a ", a, ' b ', b
  if (b .lt. pdf) then
    accept = .true.
  end if
  end do
  r = a
  ! get theta, recycle a and b
  call random_number(a)
  ctheta = 2.*a-1.
  a = acos(ctheta)
  call random_number(b)
  phi = 2.*pi*b
  x = r*sin(a)*sin(phi)
  y = r*sin(a)*cos(phi)
  z = r*ctheta
  T = temperature(r)

  ! print*,"Spawining velocities"
  vx = Sqrt(kB*T/mdm)*random_normal()
  vy = Sqrt(kB*T/mdm)*random_normal()
  vz = Sqrt(kB*T/mdm)*random_normal()

! print*,'A random number ', random_normal()

end subroutine spawn

subroutine propagate(xin,vin,xout,vout,tout)
use star
implicit none
double precision :: xin(3),vin(3),xout(3),vout(3)
double precision :: a, tau
double precision :: x,y,z,vx,vy,vz,T
double precision :: phase_i(3), amplitude_i(3) !initial conditions

!optical depth that we travel before collision
!1) draw an optical depth
call random_number(a)
tau = -log(1.d0-a)

! 2) Integrate the time it takes to reach that optical depth
!first get the initial conditions for this step
phase_i = atan(-vin,OmegaSHO*xin)
amplitude_i = xin/cos(phase_i)





end subroutine propagate



subroutine collide()
  use star
  implicit none


end subroutine collide

  !this goes into the RK solver
  !y in this subroutine is the function we're integrating (i.e. the optical depth, tau)
subroutine step(t,y,yprime,phase_i,amplitude_i)
  implicit none
  double precision :: t, y, yprime
  double precision :: ri(3), vi(3)
  double precision :: x(3),vx(3) !the integrator does not need to know these
  x = amplitude_i*cos(OmegaSHO*t+phase_i)
  vx = -amplitude_i*OmegaSHO*sin(OmegaSHO*t+phase_i)

!this needs to be a loop if you have multiple species
  call omega(x,vx,phase_i,amplitude_i,y,yprime)

end subroutine step

subroutine omega(xin,vin,phase_i,amplitude_i,omega,omegaprime)
  !compute omega and its derivative given the position and velocity of a particle
  !only scattering with a single species (hydrogen)
  use star
  implicit none
  double precision, intent(in) :: xin(3),vin(3),phase_i(3),amplitude_i(3)
  double precision :: vT,r,v2,y,omega,omegaprime,yprime,accel(3),wprefactor
  r = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
  vT = sqrt(2.*kB*temperature(r)/mdm)
  v2 = vin(1)**2 + vin(2)**2 + vin(3)**2
  y = sqrt(v2/mdm)
  !niso = 1 = hydrogen hardcoded
  wprefactor = 2.*sigSD*ndensity(r,1)*vT*sqrt(mu)
  omega = wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))
  accel = -OmegaSHO**2*xin

  yprime = 2.d0*sum(accel*vin)/mdm !y^2'
  yprime = .5/y
  omegaprime = dndr(r,1)/ndensity(r,1)*omega + yprime*wprefactor*(erf(y)*(1.-y**-2) + exp(-y**2)/sqrt(pi)/y)

end subroutine omega


! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
