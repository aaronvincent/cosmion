! This is where the action lives
! Units are cgs, except mass


subroutine spawn(x,y,z,vx,vy,vz)
  use star
double precision :: x,y,z,vx,vy,vz
logical accept
double precision :: r, N,pdf,a,b,ctheta,phi

!Basic rejection method to get r
accept = .false.
do while (accept == .false.)
call random_number(a)
a = a*Rsun
pdf = a**2*exp(-(a/rchi)**2)
call random_number(b)
b= b*rchi**2*exp(-1.d0)
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







end subroutine spawn

subroutine propagate()
use star


end subroutine propagate



subroutine collide()
  use star

end subroutine collide

! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
