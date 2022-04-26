! This is where the action lives
! Units are cgs, except mass


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

subroutine step(t,y,yprime,phase_i,amplitude_i)
  implicit none
  double precision :: t, y, yprime
  double precision :: ri(3), vi(3)
  double precision :: x(3),vx(3),accel(3) !the integrator does not need to know these


  !this goes into the RK solver


end subroutine step

subroutine omega(xin,vin,phase_i,amplitude_i,omega_out,omegaprime_out)
  use star
  implicit none
double precision, intent(in) :: xin(3),vin(3),phase_i(3),amplitude_i(3)
double precision :: vT,r,v,y
r = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
vT = sqrt(2.*kB*temperature(r)/mdm)



!compute omega and its derivative given the position and velocity of a particle







end subroutine omega


! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
