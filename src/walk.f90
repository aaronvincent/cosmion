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

subroutine propagate()
use star
implicit none
double precision :: a, tau

!optical depth that we travel before collision
call random_number(a)
tau = -log(1.d0-a)



!1) draw an optical depth



end subroutine propagate



subroutine collide()
  use star
  implicit none


end subroutine collide

! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
