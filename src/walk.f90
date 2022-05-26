! This is where the action lives
! Units are cgs, except mass
! Follow the recipe: https://arxiv.org/abs/2111.06895


subroutine spawn(x,v)
  use star
  implicit none
  double precision :: x(3),v(3),T
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
  x(1) = r*sin(a)*sin(phi)
  x(2) = r*sin(a)*cos(phi)
  x(3) = r*ctheta
  T = temperature(r)

  ! print*,"Spawining velocities"
  v(1) = Sqrt(kB*T/mdm)*random_normal()
  v(2) = Sqrt(kB*T/mdm)*random_normal()
  v(3) = Sqrt(kB*T/mdm)*random_normal()

  ! Throws the particle out of the star in one step, for testing
  !x = (/37d9,37d9,37d9/)
  !v = (/0d7,54d7,32d7/)

! print*,'A random number ', random_normal()

end subroutine spawn

subroutine propagate(xin,vin,xout,vout,tout)
use star
use init_conds
implicit none
interface
  subroutine step(t,y,yp)
    !should be rk kind but whatever
    double precision, intent(in) :: t,y
    double precision, intent(out) :: yp
  end subroutine
  subroutine pets(y,t,yp)
    !should be rk kind but whatever
    double precision, intent(in) :: t,y
    double precision, intent(out) :: yp
  end subroutine
  subroutine keplerian(xin,vin,xout,vout,tout)
    double precision, intent(in) :: xin(3),vin(3)
    double precision, intent(out) :: xout(3),vout(3),tout
  end subroutine
  end interface

double precision, intent(in) :: xin(3),vin(3)
double precision, intent(out) :: xout(3),vout(3)
double precision :: a, tau,r,vx
double precision :: T,n,mfp,dt
! double precision :: phase_i(3), amplitude_i(3) !initial conditions now in module
!for the rkf45
integer :: neqn, flag,counter
double precision :: y, yp, time, tout, relerr, abserr

time = 0.d0
tout = 0.d0
relerr = 1.d-10
abserr = 1.d-10
flag = 1
counter = 0

!optical depth that we travel before collision
!1) draw an optical depth
call random_number(a)
tau = -log(1.d0-a)

! 2) Integrate the time it takes to reach that optical depth
!first get the initial conditions for this step
phase_i = atan(-vin,OmegaSHO*xin)
amplitude_i = xin/cos(phase_i)

!this is where we call RKF45, passin the step sub

!Replaced the external block with an interface
! Instead of integrating until we hit the right optical depth, we'll guess then adjust
!a little janky
! 1) Guess a time to get to tau. Hardcoded to H only

r  = sqrt(sum(xin**2))
vx = sqrt(sum(vin**2))

n = ndensity(R,1)
!initial guess for how long we should go
mfp = 1./(n*sigSD*vx) !a mean free path in seconds
! print*,"mean free path in cm ", mfp*vx, "and in s ", mfp
tout = tau*mfp
! tout = tout/4. !this speeds things up
!tout is no longer used
!old iterative method:
! do while ((abs(y-tau)/tau .gt. 1.d-6 ) .and. (counter .lt. 10000))
! counter = counter + 1
! time = 0.d0
! flag = 1
! y = 0.d0 !initial column density
! call rkf45 (step, y, yp, time, tout, relerr, abserr, flag )
! ! print*, "flag is", flag
!  ! Newton's method to find the stopping point-- this does very well for constant density/small cross section,
!  ! not well at all if n(r) varies a lot over the trajectory
! ! print*," tau = ", tau, " y = ", y, ' tol ', abs(y-tau)/tau, ' tout ', tout
! tout = tout - (y-tau)/yp
!
! ! I'm gonna do what you might call a pro-gamer move
! !this integrates t until we reach optical depth tau
!
! ! print*,'time ', tout
! end do
! print*," tau = ", tau, " y = ", y, ' tol ', abs(y-tau)/tau, ' tout ', time
t = 0.d0
time = 0.d0
flag = 1
y = 0.d0
call rkf45 (pets, t, yp, y,tau, relerr, abserr, flag )
tout = t
xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))
r = sqrt(sum(xout**2))
! This checks if t is NaN or if r is greater than Rsun, which means the particle has left the star.
if (t /= t .or. r>Rsun*0.9995) then  ! Use 0.9995 times the radius because of interp1 subroutine
  ! print*,"Elvis has left the building"
  ! Here we include code for the Keplerian orbit of the particle.
  ! First, we determine the path from its initial position to the solar radius.
  ! Then, we use the orbit to find the re-entry parameters.

  ! The following lines calculate the time that it would take the particle to leave the star.
  ! This didn't work
  !vr = sum(xin*vin) / r
  !c = sqrt(3./(4.*pi*GN*rhoSHO))
  !b = sqrt(4.*pi*GN*rhoSHO/(3.*vr**2.+4.*pi*GN*rhoSHO*r**2.))
  !tout = c*(asin(Rsun*0.9995*b)-asin(r*b)) ! Use 0.9995 times the radius because of interp1 subroutine
  ! This didn't work either
  !vr = sum(xin*vin) / r
  !tout = 2/OmegaSHO * atan((vr-sqrt(vr**2+omegaSHO**2*(r**2-Rsun**2)))/(omegaSHO*(r+Rsun)))
  ! Third time's the charm (just loop and keep checking if it's at the boundary yet)

  ! The Keplerian orbit places the particle at the boundary of the star.
  ! So, if the particle is already outside the star, this moves inside the boudary.
  ! Should only take one step, which is a negligible time.
  tout = 0.d0
  dt = 1.d-4 ! Determines precision
  r = sqrt(sum(xin**2))
  if (r>=Rsun*0.99949) then
    do while (r>=Rsun*0.99949) ! Use 0.9995 times the radius because of interp1 subroutine
      tout = tout+dt
      xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
      xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
      xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))
      r = sqrt(sum(xout**2))
    end do
    if (tout/dt > 1.) then
      print*,"Too many steps for particle to re-enter star."
      print*,"steps=",tout/dt
    end if
  end if
  do while (r<Rsun*0.99949) ! Use 0.9995 times the radius because of interp1 subroutine
    tout = tout+dt
    xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
    xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
    xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))
    r = sqrt(sum(xout**2))
  end do
end if

! print*, "took ", counter, ' tries ', ' guess ', tau*mfp/4., 'actual ', tout
! print*,"using pro move: "
! print*," tau = ", tau, " y = ", y, ' tol ', abs(y-tau)/tau, ' tout ', t
! print*," Did the loop, tout = ", tout


!xout(1) =  amplitude_i(1)*cos(OmegaSHO*tout+phase_i(1))
!xout(2) =  amplitude_i(2)*cos(OmegaSHO*tout+phase_i(2))
!xout(3) =  amplitude_i(3)*cos(OmegaSHO*tout+phase_i(3))

vout(1) = -amplitude_i(1)*OmegaSHO*sin(OmegaSHO*tout+phase_i(1))
vout(2) = -amplitude_i(2)*OmegaSHO*sin(OmegaSHO*tout+phase_i(2))
vout(3) = -amplitude_i(3)*OmegaSHO*sin(OmegaSHO*tout+phase_i(3))
! print*,'xout ', xout

end subroutine propagate



subroutine collide(x,vin,vout)
  !turns old velocity into new velocity
  use star
  implicit none
  integer niso
  double precision, intent(in) :: x(3),vin(3)
  double precision, intent(out) :: vout(3)
  double precision :: v(3),vnuc(3),unuc,s(3),T,r,vcm,a,b
  double precision :: ctheta, phi,random_normal !outgoing angles in COM frame
!) select a species to collide with. Hardcoded for now
  niso = 1
  !a little different from Hannah's method: we draw 3 nuclear velocities from a local MB distribution
  v = vin
  r = sqrt(sum(x**2))
  T = temperature(r)

  vnuc(1) = Sqrt(kB*T/(AtomicNumber(niso)*mnucg))*random_normal()
  vnuc(2) = Sqrt(kB*T/(AtomicNumber(niso)*mnucg))*random_normal()
  vnuc(3) = Sqrt(kB*T/(AtomicNumber(niso)*mnucg))*random_normal()

! 2) Boost to CM frame
s = (mdm*v + AtomicNumber(niso)*mnucg*vnuc)/(mdm+AtomicNumber(niso)*mnucg)
vcm = sqrt(sum((v-s)**2)) !velocity in CM frame

!Scattering is isotropic, so the new angle does not depend on the old one
call random_number(a)
ctheta = 2.*a-1.
a = acos(ctheta)
call random_number(b)
phi = 2.*pi*b
vout(1) = vcm*sin(a)*sin(phi)
vout(2) = vcm*sin(a)*cos(phi)
vout(3) = vcm*ctheta
!boost back to star frame

vout = vout + s


end subroutine collide

  !this goes into the RK solver
  !y in this subroutine is the function we're integrating (i.e. the optical depth, tau)
subroutine step(t,y,yprime)
  use init_conds
  use star
  implicit none
  double precision, intent(in) :: t,y
  double precision, intent(out) ::  yprime
  double precision :: ri(3), vi(3)
  double precision :: x(3),vx(3) !the integrator does not need to know these
  x = amplitude_i*cos(OmegaSHO*t+phase_i)
  vx = -amplitude_i*OmegaSHO*sin(OmegaSHO*t+phase_i)
  ! print*,'calling step, x = ', x
!this needs to be a loop if you have multiple species
  !y is not used
  call omega(x,vx,phase_i,amplitude_i,yprime)

end subroutine step

!inverse of steps(). dt/dtau = yprimeinv = 1/w
subroutine pets(y,t,yprimeinv)
  use init_conds
  use star
  implicit none
  double precision, intent(in) :: t,y
  double precision, intent(out) ::  yprimeinv
  double precision :: ri(3), vi(3),yprime
  double precision :: x(3),vx(3) !the integrator does not need to know these
  x = amplitude_i*cos(OmegaSHO*t+phase_i)
  vx = -amplitude_i*OmegaSHO*sin(OmegaSHO*t+phase_i)
  ! print*,'calling step, x = ', x
!this needs to be a loop if you have multiple species
  !y is not used
  call omega(x,vx,phase_i,amplitude_i,yprime)
  yprimeinv = 1.d0/yprime

end subroutine pets


subroutine omega(xin,vin,phase_i,amplitude_i,omega_out) !,omegaprime)
  !compute omega and its derivative given the position and velocity of a particle
  !only scattering with a single species (hydrogen)
  use star
  implicit none
  double precision, intent(in) :: xin(3),vin(3),phase_i(3),amplitude_i(3)
  double precision :: vT,r,v2,y,omega_out,omegaprime,yprime,accel(3),wprefactor
  ! The following line and conditional check if the particle is inside the star.
  ! If it left the star, it returns NaN, which makes the integrator return NaN for the time t.
  ! It's sketchy, but we can then easily identify if the particle left after integration.
  if (sqrt(sum(xin**2)) .ge. Rsun*0.9995) then ! Use 0.9995 times the radius because of interp1 subroutine
    omega_out = 0.0
    omega_out = omega_out/omega_out
    return
  end if
  r = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
  vT = sqrt(2.*kB*temperature(r)/mdm)
  v2 = vin(1)**2 + vin(2)**2 + vin(3)**2
  y = sqrt(v2/mu/vT**2)

! print*, "Omega: R ", r, " vT ", vT, " v ", sqrt(v2), " y ", y

  !niso = 1 = hydrogen hardcoded
  wprefactor = 2.*sigSD*ndensity(r,1)*vT*sqrt(mu)
  omega_out = wprefactor*((y+.5/y)*erf(y)+1./sqrt(pi)*exp(-y**2))
  !the next bit is me not understanding wtf is going on. I made derivatives yay. Ignore
!   accel = -OmegaSHO**2*xin
!
!   yprime = 2.d0*sum(accel*vin)/mdm !y^2'
!   yprime = .5/y
!   omegaprime = dndr(r,1)/ndensity(r,1)*omega_out + yprime*wprefactor*(erf(y)*(1.-1./y**2) + exp(-y**2)/sqrt(pi)/y)
! print*,'y ', y, 'yprime ', yprime
end subroutine omega

subroutine cross(x,y,z)
  implicit none
  double precision, intent(in) :: x(3),y(3)
  double precision, intent(out) :: z(3)
  z(1) = x(2)*y(3)-x(3)*y(2)
  z(2) = -x(1)*y(3)+x(3)*y(1)
  z(3) = x(1)*y(2)-x(2)*y(1)
end subroutine

subroutine keplerian(xin,vin,xout,vout,tout)
  use star
  implicit none
  interface
    subroutine cross(x,y,z)
      double precision, intent(in) :: x(3),y(3)
      double precision, intent(out) :: z(3)
    end subroutine
  end interface
  double precision, intent(in) :: xin(3),vin(3)
  double precision, intent(out) :: xout(3),vout(3),tout
  double precision :: r,vr,vtot,vesc
  double precision :: h(3),vh(3),smaj,e,theta,norm(3),plvec(3),ang
  double precision :: area,areatot,period,c
  ! Here we must include code for the Keplerian orbit of the particle
  ! First, we integrate the path from its initial position to the solar radius.
  ! Then, we use the orbit to find the re-entry parameters.
  ! Do Keplerian stuff
  r = sqrt(sum(xin**2))
  vtot = sqrt(sum(vin**2))
  vr = sum(xin*vin) / r
  vesc = 2.*Rsun*sqrt(2.*pi*GN*rhoSHO/3.)
  ! Check if the particle exceeds the escape velocity
  if (vtot >= vesc) then
    print*,"The particle has escaped!"
    ! We'll have to stop and respawn the particle in this case
  else
    !print*,"It's coming back"
    ! Do Keplerian stuff
    smaj = 4*pi*GN*Rsun**3*rhoSHO*r / (8*pi*GN*Rsun**3*rhoSHO-3*vtot**2*r)
    call cross(xin,vin,h)
    e = sqrt(1-sum(h**2)/(smaj*4*pi*GN*Rsun**3*rhoSHO/3))
    theta = acos((smaj-e**2*smaj-r)/(e*r)) ! fixed?
    call cross(xin,vin,norm)
    norm = norm/sqrt(sum(norm**2))
    call cross(norm,xin,plvec)
    xout = cos(2*pi-theta*2)*xin+sin(2*pi-theta*2)*plvec
    vr = sum(xin*vin) / r
    ang = acos(vr/vtot)
    call cross(norm,xout,plvec)
    vout = cos(pi-ang)*xout+sin(pi-ang)*plvec
    vout = (vout/r)*vtot
    period = sqrt(3*pi*smaj**3/(GN*Rsun**3*rhoSHO))
    areatot = pi*smaj**2*sqrt(1-e**2)
    c = 2*atanh((e-1)*tan(theta/2)/sqrt(cmplx(e**2-1)))/sqrt(cmplx(e**2-1))
    area = smaj**2*(e**2-1)*(e*sin(theta)/(e*cos(theta)+1)-c)
    tout = (1-area/areatot)*period
    print*,"a=",smaj
    print*,"e=",e
    print*,"theta=",theta
    print*,"norm=",norm
    print*,"xin=",xin
    print*,"vin=",vin
    print*,"xout=",xout
    print*,"vout=",vout
    !print*,"area",area
    !print*,"areatot",areatot
    !print*,"tout=",tout
    print*,"tfrac=",tout/period
    print*,""
  end if
end subroutine

! function spawnCDF(r,r0)
!   double precision :: r,r0, CDF,N
!   N = sqrt(pi)*r0**3/4.
!   CDF =
!   return CDF
! end function spawnCDF
